#!/usr/bin/env bash
# Publish a plots directory to the NERSC web portal (project m3246, user gregork).
#
#   source env.sh && bash publish_plots.sh              # $PATH_TO_HISTOGRAMS/plots
#   bash publish_plots.sh /path/to/other/plots          # explicit directory
#   DRY_RUN=1 bash publish_plots.sh                     # show what would sync
#
# Served from the project CFS www/ dir:
#   dtn:/global/cfs/cdirs/m3246/www/gregork/fccjetbenchmarks/
#   https://portal.nersc.gov/cfs/m3246/gregork/fccjetbenchmarks/dashboard.html
#
# Only the plots tree is synced -- raw histograms are not published.
# The sync does not --delete: files removed locally remain on NERSC.
#
# Auth: SSH key at ~/.ssh/nersc (sshproxy cert, valid 24 h). Refresh with
# `sshproxy -u gregork` (interactive). Expired key shows as "Permission denied
# (publickey)" or "Too many authentication failures".
#
# Portal only serves world-readable files in world-executable dirs, so
# --chmod=Do+rx,Fo+r is set and -a/-p are deliberately NOT used.
set -euo pipefail

PLOTS_DIR="${1:-${PATH_TO_HISTOGRAMS:?set PATH_TO_HISTOGRAMS or pass a directory}/plots}"
DEST="${NERSC_WWW_DEST:-dtn:/global/cfs/cdirs/m3246/www/gregork}"
REMOTE_SUBDIR="${NERSC_REMOTE_SUBDIR:-fccjetbenchmarks}"
URL_BASE="https://portal.nersc.gov/cfs/m3246/gregork/${REMOTE_SUBDIR}"

if [[ ! -d "$PLOTS_DIR" ]]; then
    echo "No such directory: $PLOTS_DIR" >&2
    exit 1
fi
if [[ ! -f "$PLOTS_DIR/dashboard.html" ]]; then
    echo "Warning: $PLOTS_DIR has no dashboard.html - run scripts/create_plots.sh first?" >&2
fi

if [[ -f "$PLOTS_DIR/dashboard.html" ]]; then
    # Pre-gzip for the Apache content negotiation on NERSC (see .htaccess).
    gzip -k9 -f "$PLOTS_DIR/dashboard.html"
fi
if [[ ! -f "$PLOTS_DIR/.htaccess" ]]; then
    cat > "$PLOTS_DIR/.htaccess" <<EOF
<IfModule mod_rewrite.c>
  RewriteEngine On
  RewriteBase /cfs/m3246/gregork/$REMOTE_SUBDIR/
  RewriteCond %{HTTP:Accept-Encoding} gzip
  RewriteCond %{REQUEST_FILENAME}.gz -s
  RewriteRule ^(.+)\.html\$ \$1.html.gz [L]
</IfModule>
<Files "*.html.gz">
  AddType "text/html" .gz
  AddEncoding gzip .gz
</Files>
EOF
fi

echo "About to publish:"
echo "  from: $PLOTS_DIR"
echo "        ($(find "$PLOTS_DIR" -type f | wc -l) files, $(du -sh "$PLOTS_DIR" | cut -f1); "
echo "         dashboard.html $(date -r "$PLOTS_DIR/dashboard.html" '+%Y-%m-%d %H:%M' 2>/dev/null || echo 'missing'))"
echo "  to:   $DEST/$REMOTE_SUBDIR/"
echo "  URL:  $URL_BASE/dashboard.html"

RSYNC_ARGS=(-rtm --chmod=Do+rx,Fo+r -e "ssh -o BatchMode=yes"
            --include="*/"
            --include=".htaccess"
            --include="*.html" --include="*.html.gz"
            --include="*.json" --include="*.js" --include="*.css"
            --include="*.pdf" --include="*.png" --include="*.jpg" --include="*.jpeg"
            --include="*.gif" --include="*.svg"
            --exclude="*")

if [[ "${DRY_RUN:-0}" == "1" ]]; then
    rsync --dry-run "${RSYNC_ARGS[@]}" "$PLOTS_DIR/" "$DEST/$REMOTE_SUBDIR/"
    exit 0
fi
if [[ "${FORCE:-0}" != "1" ]]; then
    read -r -p "This replaces the live site. Continue? [y/N] " reply
    [[ "$reply" == "y" || "$reply" == "Y" ]] || { echo "Aborted."; exit 1; }
fi

rsync "${RSYNC_ARGS[@]}" "$PLOTS_DIR/" "$DEST/$REMOTE_SUBDIR/"

echo
echo "Published:"
echo "  Dashboard: $URL_BASE/dashboard.html"
echo "  Verify:    curl -s -o /dev/null -w '%{http_code}' $URL_BASE/dashboard.html"
