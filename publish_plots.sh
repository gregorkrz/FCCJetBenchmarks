#!/usr/bin/env bash
# Publish a plots directory to the public site (S3 + CloudFront invalidation).
#
#   source env.sh && bash publish_plots.sh              # $PATH_TO_HISTOGRAMS/plots
#   bash publish_plots.sh /path/to/other/plots          # explicit directory
#   DRY_RUN=1 bash publish_plots.sh                     # show what would sync
#   PRUNE=1 bash publish_plots.sh                       # also delete remote-only objects
#
# This overwrites what is live at the CloudFront URL in the README, so it prints
# what it is about to do and waits for confirmation (skip with FORCE=1).
#
# By default the sync only adds and overwrites: objects that exist remotely but
# not locally are left alone. That is the safe default, but it means a *rename*
# leaves the old names live alongside the new ones - e.g. the comparison_AK/ and
# *_vs_AntiKt*.pdf paths from before the e+e- Cambridge/Aachen relabelling. Use
# PRUNE=1 to add --delete and clear them out.
set -euo pipefail

PLOTS_DIR="${1:-${PATH_TO_HISTOGRAMS:?set PATH_TO_HISTOGRAMS or pass a directory}/plots}"
BUCKET="s3://gregor-research/FCCJetBenchmarks/"
DISTRIBUTION_ID="E268ANYYX3SEJZ"

if [[ ! -d "$PLOTS_DIR" ]]; then
    echo "No such directory: $PLOTS_DIR" >&2
    exit 1
fi
if [[ ! -f "$PLOTS_DIR/dashboard.html" ]]; then
    echo "Warning: $PLOTS_DIR has no dashboard.html - run scripts/create_plots.sh first?" >&2
fi

echo "About to publish:"
echo "  from: $PLOTS_DIR"
echo "        ($(find "$PLOTS_DIR" -type f | wc -l) files, $(du -sh "$PLOTS_DIR" | cut -f1); "
echo "         dashboard.html $(date -r "$PLOTS_DIR/dashboard.html" '+%Y-%m-%d %H:%M' 2>/dev/null || echo 'missing'))"
echo "  to:   $BUCKET  (CloudFront $DISTRIBUTION_ID)"

SYNC_ARGS=()
if [[ "${PRUNE:-0}" == "1" ]]; then
    SYNC_ARGS+=(--delete)
    echo "  PRUNE=1: objects that exist only remotely will be DELETED"
fi

if [[ "${DRY_RUN:-0}" == "1" ]]; then
    aws s3 sync --dryrun "${SYNC_ARGS[@]}" "$PLOTS_DIR" "$BUCKET"
    exit 0
fi
if [[ "${FORCE:-0}" != "1" ]]; then
    read -r -p "This replaces the live site. Continue? [y/N] " reply
    [[ "$reply" == "y" || "$reply" == "Y" ]] || { echo "Aborted."; exit 1; }
fi

aws s3 sync "${SYNC_ARGS[@]}" "$PLOTS_DIR" "$BUCKET"
aws cloudfront create-invalidation --distribution-id "$DISTRIBUTION_ID" --paths "/*"
