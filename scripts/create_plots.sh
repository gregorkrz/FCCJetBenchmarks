#!/usr/bin/env bash

set -euo pipefail

# Parse arguments: an optional --html-only flag, then the positional INPUT_DIR.
HTML_ONLY=false
ARGS=()
for arg in "$@"; do
    if [[ "$arg" == "--html-only" ]]; then
        HTML_ONLY=true
    else
        ARGS+=("$arg")
    fi
done

if [[ ${#ARGS[@]} -ne 1 ]]; then
    echo "Usage: $0 [--html-only] INPUT_DIR"
    echo "  --html-only  Skip the ROOT/PDF-producing steps and only rebuild the interactive"
    echo "               dashboard (dashboard_data.json + dashboard.html) from the pickles"
    echo "               already on disk from a prior full run. Fast iteration on the"
    echo "               dashboard HTML/JS itself."
    exit 1
fi

INPUT_DIR="${ARGS[0]}"

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: INPUT_DIR does not exist or is not a directory"
    exit 1
fi

if [[ "$HTML_ONLY" == false ]]; then
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

        #fccanalysis plots src/plotting/debugging_plots.py -- \
        #    --inputDir "$METHOD_DIR"

        python src/plotting/extract_resolution_data.py \
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
    python src/plotting/joint_plots.py --inputDir $INPUT_DIR --AK-comparison
    python src/plotting/joint_plots.py --inputDir $INPUT_DIR --AK-comparison --energy-recovery

    # ------------------------------------------------------------------
    # Run final statistics command
    # ------------------------------------------------------------------
    python src/plotting/print_basic_stats.py --inputDir "$INPUT_DIR" --all-folders
else
    echo "--html-only set: skipping ROOT/PDF steps, rebuilding the dashboard from existing pickles."
fi

# ------------------------------------------------------------------
# Build the interactive resolution dashboard (packages all methods'
# resolution_dashboard_data.pkl / mass_dashboard_data.pkl and
# basic_stats_summary.json into one JSON + self-contained HTML page)
# ------------------------------------------------------------------
python src/plotting/build_dashboard_data.py --inputDir "$INPUT_DIR"
python src/plotting/make_interactive_dashboard.py --data "$INPUT_DIR/plots/dashboard_data.json"
