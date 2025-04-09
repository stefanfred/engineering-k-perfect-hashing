#!/bin/bash
#SBATCH --array=0-17

params="--numKeys 10M --numQueries 10M"
jobsBase=(
    "$params --chd"
    "$params --hashDisplace"
    "$params --pachash"
    "$params --thresholdBased"
    "$params --thresholdBasedConsensus"
    "$params --kRecSplit"
)
jobsAll=()
for i in "${!jobsBase[@]}"; do
    jobsAll+=("${jobsBase[$i]} -k 10")
    jobsAll+=("${jobsBase[$i]} -k 100")
    jobsAll+=("${jobsBase[$i]} -k 1000")
done

if [ "$#" -ne 0 ]; then
    # Use "--list" or "--jobs" or actually anything
    # Print available jobs
    for i in "${!jobsAll[@]}"; do
        echo -e "$i \t ${jobsAll[$i]}"
    done
    exit 0
fi

echo "Host: $(hostname)"
echo "Git commit: $(git rev-parse HEAD)"
git status
echo "Gcc version: $(gcc --version)"
echo "-----"

make Comparison
strings Comparison | grep fPIC

if [[ -z "$SLURM_ARRAY_TASK_COUNT" ]]; then
    echo "This is a simple batch execution. Running all tasks."
    for i in "${!jobsAll[@]}"; do
        echo "./Comparison ${jobsAll[$i]}"
        ./Comparison ${jobsAll[$i]}
    done
else
    echo "This is a slurm array job. Running only the requested task."
    sleep 5
    echo "./Comparison ${jobsAll[$SLURM_ARRAY_TASK_ID]}"
    ./Comparison ${jobsAll[$SLURM_ARRAY_TASK_ID]}
fi

