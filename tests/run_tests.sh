#!/bin/bash

# Script to build and run Kofola tests

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$SCRIPT_DIR/../build"

echo "Building Kofola tests..."
cd "$BUILD_DIR"

# Configure with tests enabled
cmake .. -DBUILD_TESTS=ON

# Build the test executable
make kofola_tests

echo "Running tests..."
./tests/kofola_tests "$@"

echo "All tests completed successfully!"

# Example usage:
# ./run_tests.sh                    # Run all tests
# ./run_tests.sh "[acc_code_dnf]"   # Run only acc_code_dnf tests
