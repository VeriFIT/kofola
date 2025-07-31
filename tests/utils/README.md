# Test Utils

This directory contains utility test programs that are separate from the main test suite.

## test_complement_equivalence

A standalone program for testing the equivalence between `kofola::complement_tela` and Spot's complement algorithm.

### Usage

```bash
./test_complement_equivalence <hoa_file>
```

### Building

The utility is built automatically as part of the main build process:

```bash
cd build
cmake ..
make
```

The executable will be located at `build/tests/utils/test_complement_equivalence`.

### Description

This program:
1. Loads an automaton from a HOA format file
2. Computes its complement using both kofola's `complement_tela` and Spot's `complement` functions
3. Checks if the resulting automata are language-equivalent
4. Reports success or failure, including information about distinguishing words if the complements differ

The program is useful for validating that kofola's complement implementation produces correct results compared to the reference implementation in Spot.
