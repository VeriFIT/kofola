# Random TELA Generation Quick Reference

This directory contains utilities for generating and testing random semi-deterministic automata with TELA accepting conditions.

## Quick Start

```bash
# Generate 5 automata and test with NCSB and Rank algorithms
./run_tela_experiment.sh -n 5 -a ncsb,rank

# Generate more complex automata with larger state space
./run_tela_experiment.sh -n 10 --min-states 10 --max-states 30 -a ncsb,ncsb-delay,rank
```

## Available Scripts

| Script | Purpose |
|--------|---------|
| `generate_random_tela.sh` | Generate random TELA automata using spot |
| `batch_test_tela.sh` | Test generated automata with kofola algorithms |
| `run_tela_experiment.sh` | Combined generation and testing workflow |

## Common Use Cases

### 1. Performance Testing
```bash
# Generate many small automata for quick performance tests
./generate_random_tela.sh -n 100 --max-states 8
./batch_test_tela.sh -a ncsb -t 10
```

### 2. Stress Testing
```bash
# Generate large complex automata for stress testing
./generate_random_tela.sh -n 20 --min-states 15 --max-states 50 --density 0.5
./batch_test_tela.sh -a rank -t 300
```

### 3. Algorithm Comparison
```bash
# Compare multiple algorithms on the same set of automata
./run_tela_experiment.sh -n 25 -a ncsb,ncsb-delay,rank,safra
```

### 4. Reproducible Testing
```bash
# Use fixed seed for reproducible results
./run_tela_experiment.sh --seed 12345 -n 20
```

## Output Structure

```
experiment_TIMESTAMP_automata/           # Generated automata
├── random_tela_001.hoa
├── random_tela_002.hoa
├── ...
└── generation_summary.txt

experiment_TIMESTAMP_results_ncsb/       # NCSB algorithm results
├── batch_test_log.txt
├── batch_test_summary.txt
└── random_tela_001_result.txt

experiment_TIMESTAMP_results_rank/       # Rank algorithm results
├── batch_test_log.txt
├── batch_test_summary.txt
└── random_tela_001_result.txt

experiment_TIMESTAMP_comparison.txt      # Algorithm comparison
```

## Requirements

- Spot library with `autfilt` and `randaut` utilities
- Built kofola binary (run `make` in the build directory)

## Tips

- Start with small automata (`-n 5 --max-states 10`) to verify setup
- Use `--stats-only` for quick validation without running full complement
- Monitor disk space when generating many large automata
- Use timeouts (`-t`) appropriate for your automata complexity
- The `--cleanup` option removes automata files after testing to save space
