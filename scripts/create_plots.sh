#!/usr/bin/env bash

set -euo pipefail

# Check input arguments
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 INPUT_DIR"
    exit 1
fi

INPUT_DIR="$1"

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: INPUT_DIR does not exist or is not a directory"
    exit 1
fi

# ------------------------------------------------------------------
# Loop over subdirectories
# ------------------------------------------------------------------
for METHOD_DIR in "$INPUT_DIR"/*; do
    # Only process directories
    [[ -d "$METHOD_DIR" ]] || continue

    METHOD_NAME="$(basename "$METHOD_DIR")"

    # Skip directories starting with 'plots'
    if [[ "$METHOD_NAME" == plots*  ]]; then
      echo "Skipping $METHOD_NAME"
      continue
    fi

    echo "Processing method: $METHOD_NAME"

    fccanalysis plots src/plotting/debugging_plots.py -- \
        --inputDir "$METHOD_DIR"

    python src/plotting/resolution_plots.py \
        --inputDir "$METHOD_DIR"

    python src/plotting/mass_plots.py \
        --inputDir "$METHOD_DIR"

    echo "Finished $METHOD_NAME"
    echo "----------------------------------------"
done

# ------------------------------------------------------------------
# Produce the summary matrix plots (comparison of different methods with respect to different metrics)
# ------------------------------------------------------------------
python src/plotting/joint_plots.py --inputDir $INPUT_DIR

# ------------------------------------------------------------------
# Run final statistics command
# ------------------------------------------------------------------
echo "Running basic statistics on full input directory..."
python src/plotting/print_basic_stats.py --inputDir "$INPUT_DIR"

