#!/bin/bash

# Check if an argument was provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <tree_string>"
    exit 1
fi

tree="$1"
n=${#tree}
X=$((2 * n - 2))

mpic++ --prefix /usr/local/share/OpenMPI -o vuv vuv.cpp

mpirun --prefix /usr/local/share/OpenMPI -np "$X" ./vuv "$tree"

rm -f vuv
