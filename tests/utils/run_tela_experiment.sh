#!/bin/bash

# Convenience script that generates random TELA automata and runs batch tests
# Combines generate_random_tela.sh and batch_test_tela.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default parameters
NUM_AUTOMATA=10
TEST_ALGORITHMS=("ncsb" "rank")
OUTPUT_BASE="experiment_$(date +%Y%m%d_%H%M%S)"
CLEAN_UP=false

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Generate random TELA automata and test them with multiple complement algorithms.

OPTIONS:
    -n, --num-automata NUM      Number of automata to generate (default: $NUM_AUTOMATA)
    -a, --algorithms ALGS       Comma-separated list of algorithms to test (default: ncsb,rank)
                               Available: ncsb, ncsb-delay, rank, safra, mh
    -o, --output-base DIR       Base name for output directories (default: experiment_TIMESTAMP)
    --min-states NUM            Minimum number of states (default: 3)
    --max-states NUM            Maximum number of states (default: 15)
    --cleanup                   Remove automata files after testing
    --seed NUM                  Random seed for reproducibility
    -h, --help                  Display this help message

EXAMPLES:
    # Quick test with 5 automata using NCSB and Rank algorithms
    $0 -n 5

    # Comprehensive test with multiple algorithms
    $0 -n 20 -a ncsb,rank,safra --min-states 5 --max-states 25

    # Reproducible experiment
    $0 --seed 42 -n 15 -o my_experiment

EOF
}

# Parse arguments
EXTRA_ARGS=()
while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--num-automata)
            NUM_AUTOMATA="$2"
            EXTRA_ARGS+=("$1" "$2")
            shift 2
            ;;
        -a|--algorithms)
            IFS=',' read -ra TEST_ALGORITHMS <<< "$2"
            shift 2
            ;;
        -o|--output-base)
            OUTPUT_BASE="$2"
            shift 2
            ;;
        --cleanup)
            CLEAN_UP=true
            shift
            ;;
        --min-states|--max-states|--seed)
            EXTRA_ARGS+=("$1" "$2")
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

AUTOMATA_DIR="${OUTPUT_BASE}_automata"
RESULTS_BASE="${OUTPUT_BASE}_results"

echo "================================================"
echo "Kofola Random TELA Experiment"
echo "================================================"
echo "Number of automata: $NUM_AUTOMATA"
echo "Algorithms to test: ${TEST_ALGORITHMS[*]}"
echo "Automata directory: $AUTOMATA_DIR"
echo "Results base: $RESULTS_BASE"
echo ""

# Step 1: Generate random automata
echo "Step 1: Generating random TELA automata..."
echo "----------------------------------------"
if ! "$SCRIPT_DIR/generate_random_tela.sh" -o "$AUTOMATA_DIR" -n "$NUM_AUTOMATA" "${EXTRA_ARGS[@]}"; then
    echo "Error: Failed to generate automata"
    exit 1
fi

echo ""

# Step 2: Test with each algorithm
echo "Step 2: Testing with complement algorithms..."
echo "---------------------------------------------"

all_success=true
for algorithm in "${TEST_ALGORITHMS[@]}"; do
    echo "Testing with algorithm: $algorithm"
    results_dir="${RESULTS_BASE}_${algorithm}"
    
    if "$SCRIPT_DIR/batch_test_tela.sh" -d "$AUTOMATA_DIR" -o "$results_dir" -a "$algorithm"; then
        echo "✓ $algorithm tests completed"
    else
        echo "✗ $algorithm tests failed"
        all_success=false
    fi
    echo ""
done

# Step 3: Generate comparison report
echo "Step 3: Generating comparison report..."
echo "--------------------------------------"

comparison_file="${OUTPUT_BASE}_comparison.txt"
cat > "$comparison_file" << EOF
Kofola Random TELA Experiment Comparison
========================================

Generated: $(date)
Number of automata: $NUM_AUTOMATA
Algorithms tested: ${TEST_ALGORITHMS[*]}

EOF

# Extract results from each algorithm's summary
for algorithm in "${TEST_ALGORITHMS[@]}"; do
    results_dir="${RESULTS_BASE}_${algorithm}"
    summary_file="$results_dir/batch_test_summary.txt"
    
    if [[ -f "$summary_file" ]]; then
        echo "Algorithm: $algorithm" >> "$comparison_file"
        echo "$(printf '=%.s' {1..20})" >> "$comparison_file"
        grep -E "(Successful|Failed|Timeouts|Success rate)" "$summary_file" >> "$comparison_file"
        echo "" >> "$comparison_file"
    else
        echo "Algorithm: $algorithm - NO RESULTS" >> "$comparison_file"
        echo "" >> "$comparison_file"
    fi
done

# Add automata statistics
echo "Automata Statistics:" >> "$comparison_file"
echo "===================" >> "$comparison_file"
if command -v autfilt &> /dev/null && [[ -f "$AUTOMATA_DIR/generation_summary.txt" ]]; then
    tail -n +10 "$AUTOMATA_DIR/generation_summary.txt" >> "$comparison_file"
else
    echo "Statistics not available" >> "$comparison_file"
fi

echo "Comparison report saved to: $comparison_file"

# Display summary
echo ""
echo "==============================================="
echo "EXPERIMENT SUMMARY"
echo "==============================================="
cat "$comparison_file"

# Cleanup if requested
if [[ "$CLEAN_UP" == "true" ]]; then
    echo ""
    echo "Cleaning up automata files..."
    rm -rf "$AUTOMATA_DIR"
    echo "Automata directory removed: $AUTOMATA_DIR"
fi

if [[ "$all_success" == "true" ]]; then
    echo ""
    echo "✓ All tests completed successfully!"
    exit 0
else
    echo ""
    echo "✗ Some tests failed. Check individual result directories for details."
    exit 1
fi
