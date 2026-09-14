#!/usr/bin/env bash
# Rename the mislabelled "AntiKt" method directories to the algorithm actually run:
# e+e- Cambridge/Aachen.
#
# Until 2026-09, src/histmaker_tools/jets.py built the radius scan with
#
#     JetClustering::clustering_ee_genkt(R, 0, 0, 0)
#
# but that constructor is (radius, exclusive, cut, sorted, recombination=0,
# exponent=0.). Only four arguments were passed, so the exponent fell back to 0,
# which is Cambridge/Aachen - not anti-kT (-1). The histograms are correct and
# unchanged; only the name was wrong. No histmaker re-run is needed, just the
# plotting stages that key off directory names.
#
# Per-method pickles live *inside* each method directory and embed no method name
# or path, so a plain mv carries them along.
#
#   DRY_RUN=1 bash scripts/rename_ca_method_dirs.sh "$PATH_TO_HISTOGRAMS"   # show
#             bash scripts/rename_ca_method_dirs.sh "$PATH_TO_HISTOGRAMS"   # do it
#
# Only $PATH_TO_HISTOGRAMS and $PATH_TO_HISTOGRAMS_OLD have these directories;
# the _nofilter and _nofilter_common trees hold Durham only.
set -euo pipefail

TREE="${1:?usage: $0 PATH_TO_HISTOGRAMS}"
TREE="${TREE%/}"
[[ -d "$TREE" ]] || { echo "no such tree: $TREE" >&2; exit 1; }

RADII=(04 06 08 10 12 14)
renamed=0
skipped=0

for r in "${RADII[@]}"; do
    for pair in "PF_AntiKtR$r:PF_EECambridgeR$r" \
                "PF_E_recovery_AntiKtR$r:PF_E_recovery_EECambridgeR$r"; do
        src="$TREE/${pair%%:*}"
        dst="$TREE/${pair##*:}"
        [[ -d "$src" ]] || continue
        if [[ -e "$dst" ]]; then
            echo "SKIP (destination exists): $dst"
            skipped=$((skipped + 1))
            continue
        fi
        if [[ "${DRY_RUN:-0}" == "1" ]]; then
            echo "would mv $src -> $dst"
        else
            mv -v "$src" "$dst"
        fi
        renamed=$((renamed + 1))
    done
done

# The SLURM logs are named after the job (AK06_<process>.stderr etc.). Harmless
# either way, so they are only reported, not moved.
if [[ -d "$TREE/logs" ]]; then
    n_logs=$(find "$TREE/logs" -maxdepth 1 -name 'AK*' -o -maxdepth 1 -name 'e_recovery_AK*' 2>/dev/null | wc -l)
    [[ "$n_logs" -gt 0 ]] && echo "note: $n_logs log file(s) under $TREE/logs still carry AK* job names (harmless)"
fi

echo
if [[ "${DRY_RUN:-0}" == "1" ]]; then
    echo "DRY RUN: $renamed directory/ies would be renamed, $skipped skipped."
    exit 0
fi
echo "Renamed $renamed directory/ies, skipped $skipped."
cat <<'EOF'

Next: regenerate the derived artefacts (no histmaker, no batch jobs). Stale
output carrying the old names has to go first, or it will linger next to the new:

  rm -rf   "$TREE"/plots/full_hist \
           "$TREE"/plots/comparison_AK "$TREE"/plots/comparison_AK_energy_recovery
  rm -f    "$TREE"/plots/JER_points_only/JER_grid_points_Durham_vs_AntiKt*.pdf

  python src/plotting/print_basic_stats.py --inputDir "$TREE" --all-folders
  python src/plotting/joint_plots.py --inputDir "$TREE"
  for F in ee-ca ee-ca-er ee-kt; do
      python src/plotting/joint_plots.py --inputDir "$TREE" --family "$F"
  done
  python src/plotting/build_dashboard_data.py --inputDir "$TREE"
  python src/plotting/make_interactive_dashboard.py --data "$TREE"/plots/dashboard_data.json
  PATH_TO_HISTOGRAMS="$TREE" python src/plotting/presentation_jer_plots.py

Publishing: publish_plots.sh syncs without --delete by default, so the
old-named objects would remain on the live site. Use PRUNE=1 to add --delete.
EOF
