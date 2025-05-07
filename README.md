# Engineering Minimal k-Perfect Hash Functions

A k-perfect hash function is a hash function that maps a set of keys to a set of integers.
Each output value might collide with up to k other keys.

This repository contains a selection of k-perfect hash function implementations.
It also includes benchmarks comparing them to other approaches from the literature.

### Benchmark framework

The benchmark framework for k-perfect hash functions can be found in the `contenders` and `benchmarklib` folders.
To build the benchmarks, clone this repository (with submodules) and run the following commands.

```
cmake -B ./build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

You can then generate the plots using `./build/Comparison`.
Use `./build/Comparison --help` to view available command line arguments.

### Library Usage

Clone this repository (with submodules) and add the following to your `CMakeLists.txt`.

```
add_subdirectory(path/to/kPerfectHashing)
target_link_libraries(YourTarget PRIVATE Kphf::<hash function to use>)
```

You can link with one of the following targets.
Each target is a different k-perfect hash function with its own set of dependencies.

- `Kphf::PaCHash`
- `Kphf::RecSplit`
- `Kphf::HashDisplace`
- `Kphf::ThresholdBasedBumping`

### License

This code is licensed under the [GPLv3](/LICENSE).

If you use this work in an academic context or publication, please cite our paper:

```
@article{hermann2025engineering,
  author = {Stefan Hermann and Sebastian Kirmayer and Hans-Peter Lehmann and Peter Sanders and Stefan Walzer},
  title = {Engineering Minimal k-Perfect Hash Functions},
  journal = {CoRR},
  volume = {abs/2504.20001},
  year = {2025},
  doi = {10.48550/ARXIV.2504.20001}
}
```
