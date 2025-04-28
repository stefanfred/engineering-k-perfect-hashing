#!/bin/bash

echo "Host: $(hostname)"
echo "Git commit: $(git rev-parse HEAD)"
git status
echo "Gcc version: $(gcc --version)"
echo "-----"

make Table
strings Table | grep fPIC

./Table -n 100m -q 100m

