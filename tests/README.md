# Kofola Tests

This directory contains unit tests for the Kofola project using the Catch2 testing framework.

## Running Tests

To build and run the tests:

```bash
cd build
cmake .. -DBUILD_TESTS=ON
make kofola_tests
./tests/kofola_tests
```

## Test Structure

- `test_main.cpp`: Main test runner (uses Catch2WithMain)
- `tela/`: Tests for TELA (Transition-based Emerson-Lei Automata) operations
  - Acceptance code DNF conversion, complementation, transition colors, safe models, partitions
- `modular/`: Tests for modular complementation algorithms
  - Modular complementation, SCC type detection
- `inclusion/`: Language inclusion tests
- `utils/`: Utility tests and complement equivalence checks
- `test_data/`: Test automata files

## Dependencies

The tests require:
- Catch2 v3.x (automatically fetched if not found)
- SPOT library
- BDDX library

## Test Coverage

The current tests cover:
- Basic functionality of `acc_code_dnf` with various acceptance conditions
- Büchi, co-Büchi, and generalized Büchi acceptance conditions
- Complex acceptance conditions (Rabin, Streett, mixed)
- Edge cases (true/false conditions, large number of acceptance sets)
- DNF structure validation
- Deterministic behavior of the computation
- Multiple acceptance sets handling
