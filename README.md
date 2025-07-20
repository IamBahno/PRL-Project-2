# MPI Tree Node Level Computation

**Subject:** PRL – Parallel and Distributed Algorithms  
**Project:** Level Computation of Tree Nodes using MPI  
**School Year:** 2024/25  
**Author:** Ondřej Bahounek

## Project Overview

This project implements a parallel algorithm using **MPI (Message Passing Interface)** in C++ to compute the levels of nodes in a binary tree. Each node of the binary tree is assigned to a process, and the algorithm determines how deep each node lies in the tree, starting from the root.

The tree is represented as a sequence of characters. For example:
```
ABC
```

is interpreted as a binary tree with three nodes: A (root), B (left child), and C (right child).

### Example Input and Output

For the input `ABC`, the output will be:
```A:0,B:1,C:1```


Each output line contains the name of a node and its level (distance from the root), separated by a colon, and entries are comma-separated.

## Files

- `vuv.cpp` – Well-documented C++ source file implementing the algorithm using MPI.
- `test.sh` – Bash script that compiles the code and runs it with the appropriate number of processes. Takes one argument – the node definition string (e.g., `ABC`).

## Usage and Testing

To test the program on a system with MPI (e.g., the **merlin** server), run:

```bash
./test.sh ABC
```

