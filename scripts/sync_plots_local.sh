#!/usr/bin/env bash
# Download the current figures from S3DF.
#
# Run this on your own machine, not on S3DF: it pulls from $HOST to the local
# directory. Running it on S3DF would copy files sideways on the same filesystem.
#
#   bash scripts/sync_plots_local.sh          # into ./fcc_plots
#   bash scripts/sync_plots_local.sh -n       # dry run, list what would transfer
#   DEST=~/talks bash scripts/sync_plots_local.sh
#
# Override HOST or REMOTE if your login node or histogram tree differs.
set -euo pipefail

HOST=${HOST:-gregork@sdfiana005}
REMOTE=${REMOTE:-/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks_histograms/IDEA_20260902_seeded}
REPO=${REPO:-/fs/ddn/sdf/group/atlas/d/gregork/fastsim/FCCJetBenchmarks}
DEST=${DEST:-./fcc_plots}

DRY=""
[[ "${1:-}" == "-n" || "${1:-}" == "--dry-run" ]] && DRY="--dry-run"
RS=(rsync -av --protect-args)
[[ -n "$DRY" ]] && RS+=("$DRY")

mkdir -p "$DEST/mh_grids" "$DEST/decomposition" "$DEST/event_displays"

# Higgs-mass grids. The per-radius figures overlay all three exponents
# (p = -1, 0, +1) plus Durham at one radius; the per-exponent figures overlay
# one exponent's six radii. Both come in a plain and an _Erecovery version.
echo "== Higgs-mass grids =="
"${RS[@]}" "$HOST:$REMOTE/plots/mh_grids/*.pdf" "$DEST/mh_grids/"

# Physics | Detector | Detector+Physics. One figure per algorithm, so the
# method name has to go into the filename or the copies would collide.
# Energy-recovery variants only, at R=0.8 (all processes) and R=1.4 (2- and
# 4-jet only: the 6-jet events do not survive the filter at that radius).
echo "== Physics | Detector | Detector+Physics =="
for m in EEAntiKtR08 EEAntiKtR14 EECambridgeR08 EECambridgeR14 EEKtR08 EEKtR14; do
  for v in "" "_by_jets" "_wide_log" "_by_jets_wide_log"; do
    src="$REMOTE/plots/mh_decomposition_PF_E_recovery_${m}/pflow/mH_decomposition_genphys_common${v}.pdf"
    "${RS[@]}" "$HOST:$src" "$DEST/decomposition/PDP_${m}${v}.pdf" 2>/dev/null \
      || echo "   missing: $m$v"
  done
done

# Event displays: one multi-page PDF per process.
echo "== event displays =="
"${RS[@]}" "$HOST:$REMOTE/plots/event_displays/*.pdf" "$DEST/event_displays/" 2>/dev/null \
  || echo "   none yet"

# The two reference documents.
echo "== documents =="
"${RS[@]}" "$HOST:$REPO/doc/*.pdf" "$DEST/"

cat <<EOF

Done. Under $DEST:
  mh_grids/mH_grid_radius_R{04..14}[_Erecovery].pdf       three exponents at fixed R
  mh_grids/mH_grid_exponent_p{m1,0,pp1}[_Erecovery].pdf   six radii at fixed exponent
  decomposition/PDP_<algo>R<rr>[_by_jets][_wide_log].pdf  Physics | Detector | Detector+Physics
  event_displays/<process>_event_displays.pdf             per-event eta-phi displays
  jet_algorithms.pdf, clustering_algorithms_slides.pdf
EOF
