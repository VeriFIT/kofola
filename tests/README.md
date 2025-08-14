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
- `test_acc_code_dnf.cpp`: Tests for the `acc_code_dnf` method

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
