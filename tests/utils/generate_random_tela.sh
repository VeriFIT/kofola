#!/bin/bash

# Script to generate random semi-deterministic automata with general TELA accepting condition
# Uses spot's autfilt utility for generation
# Note: Generates TGBA (Transition-based Generalized Büchi Automata) which is equivalent to TELA

set -e

# Default parameters
NUM_AUTOMATA=10
OUTPUT_DIR="random_tela_automata"
NUM_STATES_MIN=6
NUM_STATES_MAX=15
NUM_APS_MIN=1
NUM_APS_MAX=3
DENSITY=0.1
SEED=""

# Function to display usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Generate random semi-deterministic automata with general TELA accepting conditions using spot's autfilt.

Note: This script generates TGBA (Transition-based Generalized Büchi Automata) which is
equivalent to TELA (Transition-based Emerson-Lei Automata) for practical purposes.

OPTIONS:
    -n, --num-automata NUM      Number of automata to generate (default: $NUM_AUTOMATA)
    -o, --output-dir DIR        Output directory for generated automata (default: $OUTPUT_DIR)
    --min-states NUM            Minimum number of states (default: $NUM_STATES_MIN)
    --max-states NUM            Maximum number of states (default: $NUM_STATES_MAX)
    --min-aps NUM               Minimum number of atomic propositions (default: $NUM_APS_MIN)
    --max-aps NUM               Maximum number of atomic propositions (default: $NUM_APS_MAX)
    --density FLOAT             Transition density [0.0-1.0] (default: $DENSITY)
    --seed NUM                  Random seed (default: random)
    -h, --help                  Display this help message

EXAMPLES:
    # Generate 5 automata with default settings
    $0 -n 5

    # Generate 20 automata with larger state space
    $0 -n 20 --min-states 5 --max-states 25

    # Generate automata with specific seed for reproducibility
    $0 --seed 12345 -n 10

    # Generate automata to specific directory
    $0 -o my_test_automata -n 15

REQUIREMENTS:
    - spot library with autfilt utility must be installed and available in PATH
    - randaut utility from spot (for random automaton generation)

EOF
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--num-automata)
            NUM_AUTOMATA="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --min-states)
            NUM_STATES_MIN="$2"
            shift 2
            ;;
        --max-states)
            NUM_STATES_MAX="$2"
            shift 2
            ;;
        --min-aps)
            NUM_APS_MIN="$2"
            shift 2
            ;;
        --max-aps)
            NUM_APS_MAX="$2"
            shift 2
            ;;
        --density)
            DENSITY="$2"
            shift 2
            ;;
        --seed)
            SEED="$2"
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

# Validate parameters
if ! command -v autfilt &> /dev/null; then
    echo "Error: autfilt command not found. Please install spot library."
    exit 1
fi

if ! command -v randaut &> /dev/null; then
    echo "Error: randaut command not found. Please install spot library."
    exit 1
fi

if [[ $NUM_AUTOMATA -lt 1 ]]; then
    echo "Error: Number of automata must be at least 1"
    exit 1
fi

if [[ $NUM_STATES_MIN -gt $NUM_STATES_MAX ]]; then
    echo "Error: Minimum states ($NUM_STATES_MIN) cannot be greater than maximum states ($NUM_STATES_MAX)"
    exit 1
fi

if [[ $NUM_APS_MIN -gt $NUM_APS_MAX ]]; then
    echo "Error: Minimum APs ($NUM_APS_MIN) cannot be greater than maximum APs ($NUM_APS_MAX)"
    exit 1
fi

# Set random seed if provided
if [[ -n "$SEED" ]]; then
    export SPOT_SRAND="$SEED"
    echo "Using random seed: $SEED"
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Generating $NUM_AUTOMATA random semi-deterministic TELA automata..."
echo "Output directory: $OUTPUT_DIR"
echo "States range: $NUM_STATES_MIN-$NUM_STATES_MAX"
echo "APs range: $NUM_APS_MIN-$NUM_APS_MAX"
echo "Transition density: $DENSITY"
echo ""

# Function to generate a single automaton
generate_automaton() {
    local index=$1
    local output_file="$OUTPUT_DIR/random_tela_$(printf "%03d" $index).hoa"
    local max_attempts=100
    local attempt=0
    local found=0
    while [[ $attempt -lt $max_attempts ]]; do
        local num_states=$((RANDOM % (NUM_STATES_MAX - NUM_STATES_MIN + 1) + NUM_STATES_MIN))
        local num_aps=$((RANDOM % (NUM_APS_MAX - NUM_APS_MIN + 1) + NUM_APS_MIN))
        # Generate and filter automaton
        randaut -A "Streett 1..3" -Q "$num_states" -e "$DENSITY" "2..3" \
            | autfilt --is-semi-deterministic --name="Random Streett automaton $index (${num_states}s, ${num_aps}ap)" \
            > "$output_file"
        # Check if automaton satisfies constraints (autfilt output is non-empty and valid)
        if [[ -s "$output_file" ]] && autfilt --stats="%f,%s,%e,%a" "$output_file" > /dev/null 2>&1; then
            echo "Generated: $output_file (${num_states} states, ${num_aps} APs)"
            found=1
            break
        else
            rm -f "$output_file"
        fi
        ((attempt++))
    done
    if [[ $found -eq 1 ]]; then
        return 0
    else
        echo "Error: Could not generate valid automaton $index after $max_attempts attempts."
        return 1
    fi
}

# Generate automata
successful=0
failed=0

for i in $(seq 1 $NUM_AUTOMATA); do
    if generate_automaton $i; then
        ((successful++))
    else
        ((failed++))
    fi
done

echo ""
echo "Generation complete!"
echo "Successfully generated: $successful automata"
if [[ $failed -gt 0 ]]; then
    echo "Failed to generate: $failed automata"
fi

# Generate summary report
summary_file="$OUTPUT_DIR/generation_summary.txt"
cat > "$summary_file" << EOF
Random TELA Automata Generation Summary
======================================

Generated on: $(date)
Number of automata requested: $NUM_AUTOMATA
Successfully generated: $successful
Failed: $failed

Parameters used:
- States range: $NUM_STATES_MIN-$NUM_STATES_MAX
- APs range: $NUM_APS_MIN-$NUM_APS_MAX
- Transition density: $DENSITY
- Random seed: ${SEED:-"random"}
- Output directory: $OUTPUT_DIR

Automaton details:
EOF

# Add details for each generated automaton
for file in "$OUTPUT_DIR"/*.hoa; do
    if [[ -f "$file" ]]; then
        filename=$(basename "$file")
        stats=$(autfilt --stats="%s states, %e edges, %a acceptance sets, %g acc-name" "$file" 2>/dev/null || echo "N/A")
        echo "- $filename: $stats" >> "$summary_file"
    fi
done

echo "Summary written to: $summary_file"

# Optional: Display first few lines of a sample automaton
if [[ $successful -gt 0 ]]; then
    sample_file=$(find "$OUTPUT_DIR" -name "*.hoa" | head -1)
    echo ""
    echo "Sample automaton (first 15 lines of $(basename "$sample_file")):"
    echo "----------------------------------------"
    head -15 "$sample_file"
    echo "... (rest truncated)"
fi

echo ""
echo "To use these automata in tests, you can iterate over the generated files:"
echo "  for file in $OUTPUT_DIR/*.hoa; do"
echo "    echo \"Testing \$file\""
echo "    # Run your test with \$file"
echo "  done"
