#!/usr/bin/env bash
# End-to-end event displays for the mH "Physics" curve (h_mH_gen).
#
#   source env.sh
#   bash scripts/make_event_displays.sh                       # both stages
#   bash scripts/make_event_displays.sh --dump-only
#   bash scripts/make_event_displays.sh --draw-only
#   N_PER_WINDOW=5 MAX_FILES=1 bash scripts/make_event_displays.sh
#
# Any further arguments are passed through to src/event_displays.py, e.g.
#   bash scripts/make_event_displays.sh --process p8_ee_ZH_qqqq_ecm240
#
# Stage 1 (dump) needs ROOT + FCCAnalyses, so it runs inside the container: the
# clustering has to be redone from the source EDM4hep files, because
# $PATH_TO_HISTOGRAMS holds only TH1Ds. Stage 2 (draw) is plain matplotlib and
# runs directly, which is what makes iterating on the figure cheap.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"

SIF="${SIF:-/sdf/scratch/atlas/gregork/apptainer_tmp/alma_v0.sif}"
KEY4HEP="${KEY4HEP:-/cvmfs/sw.hsf.org/key4hep/setup.sh}"
KEY4HEP_RELEASE="${KEY4HEP_RELEASE:-2025-05-29}"

PROCESS="${PROCESS:-p8_ee_ZH_6jet_LF_ecm240}"
N_PER_WINDOW="${N_PER_WINDOW:-20}"
MAX_FILES="${MAX_FILES:-2}"

DO_DUMP=true
DO_DRAW=true
EXTRA=()
for arg in "$@"; do
    case "$arg" in
        --dump-only) DO_DRAW=false ;;
        --draw-only) DO_DUMP=false ;;
        *) EXTRA+=("$arg") ;;
    esac
done

: "${PATH_TO_HISTOGRAMS:?source env.sh first}"
OUT_DIR="$PATH_TO_HISTOGRAMS/plots/event_displays"
PAYLOAD="$OUT_DIR/${PROCESS}_payload.pkl"
PDF="$OUT_DIR/${PROCESS}_event_displays.pdf"
mkdir -p "$OUT_DIR"

if [[ "$DO_DUMP" == true ]]; then
    echo "== stage 1: selecting events and dumping the payload (in container) =="
    singularity exec -B /sdf -B /cvmfs -B /fs --nv "$SIF" /bin/bash -c "
        . '$KEY4HEP' -r '$KEY4HEP_RELEASE' >/dev/null 2>&1 &&
        cd '$REPO' &&
        python src/event_displays.py \
            --process '$PROCESS' \
            --n-per-window '$N_PER_WINDOW' \
            --max-files '$MAX_FILES' \
            --output '$PAYLOAD' \
            ${EXTRA[*]+${EXTRA[*]}}"
fi

if [[ "$DO_DRAW" == true ]]; then
    echo "== stage 2: drawing the multi-page PDF =="
    python src/plotting/event_display_plots.py --payload "$PAYLOAD" --output "$PDF"
    echo
    echo "Open: $PDF"
fi
