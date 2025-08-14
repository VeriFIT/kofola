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

## generate_random_tela.sh

A script to generate random semi-deterministic automata with general TELA accepting conditions using spot's autfilt utility.

### Usage

```bash
./generate_random_tela.sh [OPTIONS]
```

### Options

- `-n, --num-automata NUM`: Number of automata to generate (default: 10)
- `-o, --output-dir DIR`: Output directory for generated automata (default: random_tela_automata)
- `--min-states NUM`: Minimum number of states (default: 3)
- `--max-states NUM`: Maximum number of states (default: 15)
- `--min-aps NUM`: Minimum number of atomic propositions (default: 1)
- `--max-aps NUM`: Maximum number of atomic propositions (default: 3)
- `--density FLOAT`: Transition density [0.0-1.0] (default: 0.3)
- `--seed NUM`: Random seed for reproducibility
- `-h, --help`: Display help message

### Examples

```bash
# Generate 5 automata with default settings
./generate_random_tela.sh -n 5

# Generate 20 automata with larger state space
./generate_random_tela.sh -n 20 --min-states 5 --max-states 25

# Generate automata with specific seed for reproducibility
./generate_random_tela.sh --seed 12345 -n 10
```

### Requirements

- spot library with `autfilt` and `randaut` utilities must be installed and available in PATH

## batch_test_tela.sh

A script to batch test generated TELA automata with kofola complement algorithms.

### Usage

```bash
./batch_test_tela.sh [OPTIONS]
```

### Options

- `-d, --automata-dir DIR`: Directory containing .hoa automata files (default: random_tela_automata)
- `-k, --kofola-binary PATH`: Path to kofola binary (default: auto-detect from build)
- `-o, --output-dir DIR`: Output directory for test results (default: test_results)
- `-a, --algorithm ALG`: Complement algorithm to test (default: ncsb)
- `-t, --timeout SECONDS`: Timeout for each test in seconds (default: 30)
- `--stats-only`: Only collect statistics, don't run full complement
- `-h, --help`: Display help message

### Examples

```bash
# Test all automata in default directory with NCSB algorithm
./batch_test_tela.sh

# Test with specific algorithm and timeout
./batch_test_tela.sh -a rank -t 60

# Just collect statistics without running complement
./batch_test_tela.sh --stats-only
```

### Workflow Example

Complete workflow for generating and testing random automata:

```bash
# 1. Generate random TELA automata
./generate_random_tela.sh -n 50 --min-states 5 --max-states 20

# 2. Build kofola (if not already built)
cd ../../build && make

# 3. Test with NCSB algorithm
cd ../tests/utils
./batch_test_tela.sh -a ncsb

# 4. Test with different algorithm
./batch_test_tela.sh -a rank -o test_results_rank

# 5. Compare results from test_results/ and test_results_rank/ directories
```

## run_tela_experiment.sh

A convenience script that combines generation and testing of random TELA automata with multiple complement algorithms.

### Usage

```bash
./run_tela_experiment.sh [OPTIONS]
```

### Options

- `-n, --num-automata NUM`: Number of automata to generate (default: 10)
- `-a, --algorithms ALGS`: Comma-separated list of algorithms to test (default: ncsb,rank)
- `-o, --output-base DIR`: Base name for output directories (default: experiment_TIMESTAMP)
- `--min-states NUM`: Minimum number of states (passed to generation)
- `--max-states NUM`: Maximum number of states (passed to generation)  
- `--cleanup`: Remove automata files after testing
- `--seed NUM`: Random seed for reproducibility
- `-h, --help`: Display help message

### Examples

```bash
# Quick experiment with 5 automata
./run_tela_experiment.sh -n 5

# Comprehensive experiment with multiple algorithms
./run_tela_experiment.sh -n 20 -a ncsb,rank,safra --min-states 5 --max-states 25

# Reproducible experiment with cleanup
./run_tela_experiment.sh --seed 42 -n 15 --cleanup
```

### Output

The experiment script creates:
- `<output_base>_automata/`: Directory with generated automata
- `<output_base>_results_<algorithm>/`: Results for each algorithm tested
- `<output_base>_comparison.txt`: Summary comparison of all algorithms
