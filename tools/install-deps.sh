#!/usr/bin/env bash

set -euo pipefail


curl -L https://github.com/symengine/symengine/archive/v0.14.0.tar.gz | tar -xz 
rm -f v0.14.0.tar.gz 
cmake -S"symengine-0.14.0" -B"symengine-0.14.0/build" -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF -GNinja 
cmake --build "symengine-0.14.0/build" --target install 
rm -rf symengine-0.14.0 