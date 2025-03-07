#include <vector>
#include <limits>
#include <cstdint>
#include <cassert>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <unordered_map>
#include <ranges>

#include "kphf.hpp"
#include "sux/bits/EliasFano.hpp"
#include "consensus.hpp"
#include "util.hpp"

#include <ips2ra.hpp>

#if 0
/* highly illegal hack */
namespace tlx {
	template<>
	inline unsigned clz<unsigned __int128>(unsigned __int128 i) {
		return std::countl_zero(i);
	}

	template<>
	inline unsigned ctz<unsigned __int128>(unsigned __int128 i) {
		return std::countr_zero(i);
	}
};
#endif

namespace kphf {
namespace ThresholdBasedBumpingConsensus {

template<uint64_t n_thresholds>
constexpr std::array<uint64_t, n_thresholds>
  compute_thresholds(uint64_t _k, double bucket_size) {
	double lo = 0, hi = _k;

	for (int i = 0; i < 100; i++) {
		double mid = (lo+hi)/2;
		double prev = 0, cur = mid;
		for (uint64_t j = 0; j < n_thresholds-1; j++) {
			double x = cur / _k - prev / _k * std::pow(prev / cur, _k-1);
			prev = cur;
			if (x >= 1) {
				cur = std::numeric_limits<double>::infinity();
				break;
			} else {
				cur -= std::log1p(-x);
			}
		}
		if (cur > bucket_size
		  || cur / _k - prev / _k * std::pow(prev / cur, _k-1) > 1) {
			hi = mid;
		} else {
			lo = mid;
		}
	}

	auto convert = [](double x) -> uint64_t {
		x = std::round(std::ldexp(x, 64));
		if (x >= std::ldexp(1, 64)) {
			return std::numeric_limits<uint64_t>::max();
		} else {
			return uint64_t(x);
		}
	};

	std::array<uint64_t, n_thresholds> res;
	double prev = 0, cur = lo;
	for (uint64_t idx = 0; idx < n_thresholds; idx++) {
		res[idx] = convert(cur / bucket_size);
		double x = cur / _k - prev / _k * std::pow(prev / cur, _k-1);
		prev = cur;
		if (x == 1.0) cur = std::numeric_limits<double>::infinity();
		else cur -= std::log1p(-x);
	}

	return res;
}

#if 0
template<uint64_t n_thresholds>
constexpr auto exp_success(uint64_t n, uint64_t k, const std::array<uint64_t, n_thresholds> &thresholds) -> double {
	double total = 0;
	for (uint64_t ti: thresholds) {
		double t = ti / (double)uint64_t(-1);
		total += exp(lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1) - k * log(t) - (n-k) * log(1-t));
	}
	return total;
}

template<uint64_t n_thresholds>
constexpr auto choose_error(uint64_t k, uint64_t n, const std::array<uint64_t, n_thresholds> &thresholds) -> std::pair<uint64_t, uint64_t> {
	uint64_t e = std::min(n, k);
	double need = 1.05 - exp_succ(n, e, thresholds);
	if (need <= 0.0) return {1,0};
	for (uint64_t x = 1; x <= e; x++) {
		double s = exp_succ(n, e-x, thresholds);
		if (s >= need) {
			return {e, (uint64_t)(need/s * (double)uint64_t(-1))};
		}
		need -= s;
	}
	return {e+1,0};
}

template<uint64_t n_thresholds, uint64_t limit>
constexpr auto choose_errors(uint64_t k, const std::array<uint64_t, n_thresholds> &thresholds) -> std::array<std::pair<uint64_t, uint64_t>, limit> {
	std::array<std::pair<uint64_t, uint64_t>, limit> res;
	for (uint64_t n = 0; n < limit; n++) {
		res[n] = choose_error(k, n, thresholds);
	}
	return res;
}
#endif

#ifdef STATS
uint64_t perfect_thresholds = 0, total_thresholds = 0, extra_bumped = 0;
uint64_t bumped_keys = 0, total_keys = 0;
uint64_t filled_buckets = 0, overfull_buckets = 0;

void reset_stats() {
	perfect_thresholds = total_thresholds = extra_bumped = 0;
	bumped_keys = total_keys = 0;
	filled_buckets = overfull_buckets = 0;
}

void dump_stats() {
	std::cout << "Perfect thresholds: "
		<< 100.0 * perfect_thresholds / total_thresholds << "%\n";
	std::cout << "Avg. extra bumped: "
		<< (double) extra_bumped / total_thresholds << "\n";
	std::cout << "Bumped keys: " << 100.0 * bumped_keys / total_keys << "%\n";
	std::cout << "Filled buckets: " << 100.0 * filled_buckets / total_thresholds << "%\n";
	std::cout << "Overfull buckets: " << 100.0 * overfull_buckets / total_thresholds << "%\n";
}
#endif

struct Key {
	uint64_t bucket, fingerprint;
	Hash128 hash;
};

namespace {
	template<typename R>
	void sort_buckets(R &&r) {
		ips2ra::sort(r.begin(), r.end(),
		  [](const Key &key) -> uint64_t { return key.bucket; });
	}

	template<typename R>
	void sort_fingerprints(R &&r) {
		ips2ra::sort(r.begin(), r.end(),
		  [](const Key &key) -> uint64_t { return key.fingerprint; });
	}
}

template<uint64_t K, double OVERLOAD, int THRESHOLD_SIZE, typename PHF>
class ThresholdBasedBumpingConsensus {
	static_assert(OVERLOAD > 1.0);
	static_assert(THRESHOLD_SIZE >= 1);
	static_assert(K > 0);

private:
	static constexpr uint64_t _k = K;
	static constexpr double overload = OVERLOAD;
	static constexpr uint64_t threshold_size = THRESHOLD_SIZE;
	static constexpr uint64_t n_regions = uint64_t(1) << threshold_size;
	static constexpr uint64_t n_thresholds = n_regions;
	static constexpr double overload_bucket_size = _k * overload;

	//static constexpr std::array<uint64_t, n_thresholds> avail_thresholds = compute_thresholds<n_thresholds>(_k, overload_bucket_size);
	static constexpr std::array<uint64_t, n_thresholds> avail_thresholds = {0u,7905747460161224704u,8359656723032226816u,11209380810981552128u,11209564803788365824u,11221594445794713600u,11276579379956709376u,11288430531065952256u,11332799341920528384u,11342481969361567744u,11344680738781210624u,11359181991566505984u,11363474910187046912u,11568734981661212672u,11832245536448208896u,11885686365579614208u,11980505356276629504u,11995954856412207104u,12017535134456393728u,12364569529833445376u,12390947090743441408u,12783335344582895616u,12870216783903963136u,12920039626888509440u,12923028927110004736u,13050390665884866560u,13295816408610588672u,13303062202347474944u,13323680574124595200u,13421228870775138304u,13448223989524359168u,13490547365284313088u,14150172056635172864u,14230181064797011968u,14343485657561276416u,14375294195601678336u,14474162065318088704u,14486947284777529344u,14521922665052733440u,14564431634759616512u,14899691169639270400u,14923491076626616320u,15177323683791149056u,15476934380018302976u,15606531015835721728u,15630731142253086720u,15689191871037216768u,15993127398723125248u,16030189648959875072u,16106253426251331584u,16204964460477093888u,16641624385250715648u,16672703444383512576u,16804830772426848256u,16875549673689362432u,17085179239725324288u,17224996427921891328u,17503828862107926528u,17567973965200832512u,17731697708159930368u,17807922867090583552u,18053269036571664384u,18249541562612258816u,18346382559634016256u};
	//static constexpr std::array<std::pair<uint64_t, uint64_t>, 261> errors = choose_errors<n_thresholds, 261>(_k, avail_thresholds);
	//static constexpr std::array<std::pair<uint64_t, uint64_t>, 261> errors = {std::pair{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,392013137091346240u},{1,1139708386398813056u},{1,1931748125539233792u},{1,2771391464666411520u},{1,3661975118988153344u},{1,4606940686312779264u},{1,5609858848618733568u},{1,6674451318415364096u},{1,7804611168095846400u},{1,9004422044777951232u},{1,10278176671441104896u},{1,11630394958486155264u},{1,13065841991678441472u},{1,14589546118251909120u},{1,16206817319257651200u},{1,17923266030667259904u},{2,627446048028978304u},{2,1546914389409404672u},{2,2504899856496407040u},{2,3503443477781552128u},{2,4544670611517535744u},{2,5630798930477588480u},{2,6764146162196584448u},{2,7947137683050243072u},{2,9182314050706025472u},{2,10472338548376496128u},{2,11820004805351882752u},{2,13228244551086372864u},{2,14700135554358161408u},{2,16238909794461071360u},{2,17847961907784026112u},{3,623869206330489984u},{3,1625819611784296192u},{3,2660300363526443520u},{3,3728750723680119296u},{3,4832672973067215872u},{3,5973636153106227200u},{3,7153279818781456384u},{3,8373317828112084992u},{3,9635542191308679168u},{3,10941827001005129728u},{3,12294132463401400320u},{3,13694509048972201984u},{3,15145101780368353280u},{3,16648154674374500352u},{3,18206015354196637696u},{4,878346910536907520u},{4,1938473135499285760u},{4,3026773638033993216u},{4,4144344642927478784u},{4,5292330077827712000u},{4,6471923741970545664u},{4,7684371514899043328u},{4,8930973614808476672u},{4,10213086915677685760u},{4,11532127331984629760u},{4,12889572279505745920u},{4,14286963220373712896u},{4,15725908300534827008u},{4,17208085087435808768u},{5,198813346228043488u},{5,1278018029045627136u},{5,2381213555462550016u},{5,3509238604367326208u},{5,4662967581805503488u},{5,5843312013873160192u},{5,7051221969962784768u},{5,8287687521225723904u},{5,9553740238918518784u},{5,10850454737277255680u},{5,12178950265450330112u},{5,13540392352961196032u},{5,14935994513133283328u},{5,16367020008921571328u},{5,17834783685533145088u},{6,648396475952692352u},{5,1843188744037168640u},{4,3511431356373527552u},{3,6495438616357370880u},{2,12021723011582734336u},{2,2262425024359645952u},{1,13192067935048235008u},{1,7967138465736371200u},{1,4392539093295746560u},{1,1800763234010587904u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,676964084583448832u},{1,2103755056022692352u},{1,3723705974823207936u},{1,5568142508722549760u},{1,7674336853159241728u},{1,10086732677233711104u},{1,12858451936832993280u},{1,16053154864476989440u},{2,1521759670419965696u},{2,6589055588438309888u},{2,12576426520225544192u},{3,1490997079016282112u},{3,11863801964940926976u},{4,7546618532547097600u},{5,11133712664063131648u},{8,748915998701075456u},{100,1168542735169222912u},{100,3710919974667152384u},{100,5933947739751833600u},{100,7868361324895917056u},{100,9543418539512940544u},{100,10986747911477641216u},{100,12224231877257705472u},{100,13279925644925581312u},{100,14176011133730926592u},{100,14932784247153338368u},{100,15568672774975367168u},{100,16100281475098978304u},{100,16542460370839117824u},{100,16908392011801534464u},{100,17209693370000580608u},{100,17456528151523598336u},{100,17657725565051314176u},{100,17820901965613733888u},{100,17952582248626886656u},{100,18058318370802210816u},{100,18142802890202781696u},{100,18209975921631270912u},{100,18263124375378302976u},{100,18304972772491368448u},{100,18337765298824517632u},{100,18363339068626933760u},{100,18383188815788087296u},{100,18398523419601901568u},{100,18410314806933297152u},{100,18419339860269045760u},{100,18426216008482527232u},{100,18431431191558330368u},{100,18435368879159482368u},{100,18438328792346236928u},{100,18440543933805074432u},{100,18442194479583827968u},{100,18443419028639946752u},{100,18444323648699865088u},{100,18444989100421781504u},{100,18445476568382019584u},{100,18445832178084583424u},{100,18446090533685401600u},{100,18446277471677945856u},{100,18446412191420618752u},{100,18446508893839005696u},{100,18446578034601334784u},{100,18446627277080598528u},{100,18446662213034809344u},{100,18446686904685664256u},{100,18446704290306168832u},{100,18446716486107502592u},{100,18446725009803548672u},{100,18446730945342969856u},{100,18446735063701231616u},{100,18446737911022589952u},{100,18446739872638539776u},{100,18446741219345143808u},{100,18446742140706592768u},{100,18446742768910309376u},{100,18446743195785861120u},{100,18446743484885080064u},{100,18446743680028219392u},{100,18446743811319896064u},{100,18446743899365939200u},{100,18446743958221985792u},{100,18446743997440563200u},{100,18446744023491909632u},{100,18446744040742819840u},{100,18446744052131403776u},{100,18446744059626625024u},{100,18446744064544970752u},{100,18446744067762526208u},{100,18446744069861513216u},{100,18446744071226497024u},{100,18446744072111757312u},{100,18446744072684279808u},{100,18446744073053378560u},{100,18446744073290618880u},{100,18446744073442924544u},{100,18446744073540179968u},{100,18446744073602308096u},{100,18446744073641891840u},{100,18446744073666795520u},{100,18446744073682786304u},{100,18446744073692747776u},{100,18446744073699039232u},{100,18446744073702971392u},{100,18446744073705330688u},{100,18446744073706903552u},{100,18446744073707952128u},{100,18446744073708476416u},{100,18446744073709000704u},{100,18446744073709262848u},{100,18446744073709262848u}};
	static constexpr std::array<std::pair<uint64_t, uint64_t>, 261> errors = {std::pair{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,0u},{1,263853258438344224u},{1,630213768999269632u},{1,1114913473115578496u},{1,1736756047603053312u},{1,2517236559532976640u},{1,3481310413746565120u},{1,4658304974004731904u},{1,6083008485716885504u},{1,7796982269833041920u},{1,9850154757564135424u},{1,12302770320472934400u},{1,15227782328302909440u},{2,317433873964596800u},{2,5316620689067447296u},{2,11402897223015800832u},{3,450160397806681472u},{3,11853638089916563456u},{4,9747176399764981760u},{6,7300020384972789u},{9,1581441414762297088u},{12,8888444466983422976u},{14,6567507240073812992u},{15,8868272148931039232u},{16,4032418609023805952u},{16,12991601024551958528u},{17,421336675963633024u},{17,3733338357908564480u},{17,5189354776670561280u},{17,5123900944900285440u},{17,3807838056869652992u},{17,1459700502175545344u},{16,16708187575450542080u},{16,12847948671773816832u},{16,8457963724120313856u},{16,3637082399395802112u},{15,16936469420650373120u},{15,11619779611335278592u},{15,6150533058443402240u},{15,568979469967588544u},{14,13508823142073249792u},{14,8048181796669641728u},{14,2642242829528986112u},{13,15852258114589992960u},{13,10856337540001972224u},{13,6042207868327910400u},{13,1426528636310419456u},{12,15616121362558298112u},{12,11706134884888887296u},{12,8105944117313065984u},{12,4840467983689901056u},{12,1940362589180556544u},{11,17926047346485008384u},{11,16027281315980945408u},{11,14616724434264338432u},{11,13741929695996194816u},{11,13459249753899776000u},{11,13834918531652726784u},{11,14946295096093386752u},{11,16883297849879207936u},{12,1473971343152216064u},{12,5943400880275251200u},{12,11834944185665208320u},{13,1071420869088410112u},{13,12183269904847058944u},{14,9141478694905138176u},{15,14128403900375248896u},{18,3592172323928075264u},{100,216402355847167680u},{100,1025827800192907776u},{100,1833925970542982656u},{100,2638312747812516352u},{100,3436670143719533568u},{100,4226764624008570880u},{100,5006463785422545920u},{100,5773751234189582336u},{100,6526739550222146560u},{100,7263681257571866624u},{100,7982977757841183744u},{100,8683186217580482560u},{100,9363024433434576896u},{100,10021373728450922496u},{100,10657279959917903872u},{100,11269952742580535296u},{100,11858763010915942400u},{100,12423239060443871232u},{100,12963061220630312960u},{100,13478055320865454080u},{100,13968185116480745472u},{100,14433543844156911616u},{100,14874345075077001216u},{100,15290913030928574464u},{100,15683672521827919872u},{100,16053138657360193536u},{100,16399906472203624448u},{100,16724640596768038912u},{100,17028065091106349056u},{100,17310953547515273216u},{100,17574119553922473984u},{100,17818407596659546112u},{100,18044684467864836096u},{100,18253831229657815040u},{100,18446735774666692608u}};

	uint64_t n;
	std::vector<std::pair<uint64_t, Consensus>> layers;
	PHF phf;
	mutable sux::bits::EliasFano<> gaps;

	ThresholdBasedBumpingConsensus(uint64_t n,
	  std::vector<std::pair<uint64_t, Consensus>> &&layers,
	  PHF &&phf,
	  sux::bits::EliasFano<> &&gaps): n(n),
	  layers(std::move(layers)),
	  phf(std::move(phf)), gaps(std::move(gaps)) {}

public:
	template<typename F>
	ThresholdBasedBumpingConsensus(const std::vector<Hash128> &keys, F &&build_phf):
	  ThresholdBasedBumpingConsensus(build(keys, std::forward<F>(build_phf))) {}

	uint64_t operator()(Hash128 key) const {
		uint64_t offset = 0;
		for (uint64_t i = 0; i < layers.size(); i++) {
			auto &[cur_buckets, consensus] = layers[i];
			uint64_t h = remix(key.hi + i);
			uint64_t b = rescale(h, cur_buckets);
			auto [seed,tidx] =
			  consensus.get(b * threshold_size, threshold_size);
			tidx = decrypt(seed, tidx);
			uint64_t f = remix(key.lo + seed + i);

			if (f < avail_thresholds[tidx]) return offset + b;

			offset += cur_buckets;
		}

		return gaps.select(phf(key));
	}

	size_t count_bits() const {
		size_t s = 0;
		for (auto &[nb,c]: layers) s += c.count_bits() - sizeof(c) * 8;
		return sizeof(*this) * 8
			+ layers.capacity() * sizeof(layers[0])
		    + s
			+ phf.count_bits() - sizeof(phf) * 8
			+ gaps.bitCount() - sizeof(gaps) * 8
		;
	}

private:
	static uint64_t encrypt(uint64_t key, uint64_t tidx) {
		return (tidx ^ remix(key)) & ((uint64_t(1) << threshold_size) - 1);
	}

	static uint64_t decrypt(uint64_t key, uint64_t tidx) {
		return (tidx ^ remix(key)) & ((uint64_t(1) << threshold_size) - 1);
	}

	struct Builder {
		std::vector<Key> keys;
		std::ranges::subrange<std::vector<Key>::iterator> current;
		uint64_t cur_bucket, total_buckets, layer, offset;
		std::vector<uint64_t> *spots;
		std::vector<Hash128> *bumped;

		Builder(std::vector<Hash128> &hashes, uint64_t layer,
		  uint64_t buckets, uint64_t offset,
		  std::vector<uint64_t> *spots):
		  keys(std::move(prepare(hashes, layer, buckets))),
		  current(keys.begin(), keys.begin()),
		  cur_bucket(0), total_buckets(buckets), layer(layer), offset(offset),
		  spots(spots), bumped(&hashes) {
			bumped->clear();
		}

		static std::vector<Key> prepare(const std::vector<Hash128> &hashes,
		  uint64_t seed, uint64_t n_buckets) {
			auto r = std::views::transform(hashes,
			  [seed, n_buckets](Hash128 key) {
				uint64_t b = rescale(remix(key.hi + seed), n_buckets);
				return Key { b, 0, key };
			});
			std::vector<Key> keys(r.begin(), r.end());
			sort_buckets(keys);
			return keys;
		}

		std::optional<uint64_t> advance(uint64_t seed) {
			if (cur_bucket == total_buckets) return {};

			auto it = current.end();
			while (it != keys.end() && it->bucket == cur_bucket) {
				it->fingerprint = remix(it->hash.lo + seed + layer);
				++it;
			}

			current = std::ranges::subrange(current.end(), it);
			sort_fingerprints(current);

			cur_bucket++;

			return threshold_size;
		}

		uint64_t backtrack() {
			cur_bucket--;

			if (cur_bucket == 0) {
				current = std::ranges::subrange(keys.begin(), keys.begin());
				return 0;
			}

			auto it = current.begin();
			while (it != keys.begin() && prev(it)->bucket == cur_bucket-1) {
				--it;
			}
			current = std::ranges::subrange(it, current.begin());

			uint64_t cnt = 0;
			while (!spots->empty() &&
			  spots->back() == offset + cur_bucket - 1) {
				spots->pop_back();
				cnt++;
			}
			bumped->resize(bumped->size() - (current.size() + cnt - _k));

			return threshold_size;
		}

		std::optional<uint64_t> find_first(uint64_t seed) {
			if (current.size() > _k) {
				uint64_t fp = current[_k].fingerprint;
				auto it = std::ranges::upper_bound(avail_thresholds, fp);
				return find(seed, it - begin(avail_thresholds));
			} else {
				return find(seed, n_thresholds);
			}
		}

		std::optional<uint64_t> find_next(uint64_t seed, uint64_t prev) {
			return find(seed, decrypt(seed, prev));
		}

		std::optional<uint64_t> find(uint64_t seed, uint64_t prev) {
			uint64_t goal = std::min(_k, current.size());
			uint64_t error_bound, error_prob;
			if (current.size() < errors.size()) {
				std::tie(error_bound, error_prob) = errors[current.size()];
			} else {
				error_bound = _k+1, error_prob = 0;
			}
			for (uint64_t tidx = prev; tidx > 0;) {
				tidx--;
				auto it = std::ranges::lower_bound(current, avail_thresholds[tidx], {}, [](Key k) -> std::uint64_t { return k.fingerprint; });
				uint64_t idx = it - current.begin();
				uint64_t error = goal - idx;
				if (error > error_bound) break;
				if (error < error_bound || remix(seed + tidx) < error_prob) {
					for (Key &k: current | std::views::drop(idx)) {
						bumped->push_back(k.hash);
					}
					for (uint64_t x = idx; x < _k; x++) {
						spots->push_back(offset + cur_bucket - 1);
					}
					return encrypt(seed, tidx);
				}
			}
			return {};
		}
	};

	template<typename F>
	static ThresholdBasedBumpingConsensus
	  build(std::vector<Hash128> keys, F &&build_phf) {
		uint64_t n = keys.size();
		uint64_t total_buckets = (n + _k - 1) / _k;

		std::vector<std::pair<uint64_t, Consensus>> layers;

#ifdef STATS
		total_keys += n;
#endif

		uint64_t offset = 0;
		std::vector<uint64_t> spots;
		for (uint64_t i = 0; offset != total_buckets; i++) {
			uint64_t remaining = total_buckets - offset;
			uint64_t cur_buckets =
			  std::ceil(keys.size() / overload_bucket_size);
			if (cur_buckets >= remaining) cur_buckets = remaining;
			else if ((double)(keys.size() - _k * cur_buckets) / (remaining - cur_buckets) > overload_bucket_size) cur_buckets = remaining;

			Builder builder(keys, i, cur_buckets, offset, &spots);
			Consensus consensus(builder);

			layers.emplace_back(cur_buckets, std::move(consensus));

			offset += cur_buckets;
		}
		layers.shrink_to_fit();

#ifdef STATS
		bumped_keys += keys.size();
#endif

		PHF phf = std::forward<F>(build_phf)(std::as_const(keys));
		std::vector<uint64_t> actual_spots;
		for (Hash128 k: keys) {
			uint64_t h = phf(k);
			if (h >= actual_spots.size()) actual_spots.resize(h+1);
			actual_spots[h] = 1;
		}
		auto it = spots.begin();
		for (uint64_t &x: actual_spots) {
			uint64_t v = *it;
			if (x) ++it;
			x = v;
		}

		return ThresholdBasedBumpingConsensus(n, std::move(layers),
			std::move(phf),
			sux::bits::EliasFano<>(actual_spots, total_buckets)
		);
	}
};

}
}
