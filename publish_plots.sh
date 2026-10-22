#!/usr/bin/env bash
# Sync a plots tree to any rsync-able web host.
#
# Required env vars:
#   PUBLISH_DEST   rsync destination, e.g.  user@host:/var/www/project
#   PUBLISH_URL    base URL it is served from, e.g.  https://host/project
#
# Optional:
#   PUBLISH_HTACCESS=1   write an Apache .htaccess serving pre-gzipped .html.gz
#   DRY_RUN=1            show what rsync would do, no changes
#   FORCE=1              skip the "replace live site?" prompt
#
# Usage:
#   source env.sh
#   PUBLISH_DEST=user@host:/path PUBLISH_URL=https://host/path bash publish_plots.sh
#   bash publish_plots.sh /path/to/other/plots
#
# Only the plots tree is synced (not the raw histograms). rsync does NOT --delete,
# so files removed locally remain on the remote.
set -euo pipefail

PLOTS_DIR="${1:-${PATH_TO_HISTOGRAMS:?set PATH_TO_HISTOGRAMS or pass a directory}/plots}"
DEST="${PUBLISH_DEST:?set PUBLISH_DEST (rsync destination, e.g. user@host:/var/www/project)}"
URL_BASE="${PUBLISH_URL:?set PUBLISH_URL (base URL, e.g. https://host/project)}"
URL_BASE="${URL_BASE%/}"

if [[ ! -d "$PLOTS_DIR" ]]; then
    echo "No such directory: $PLOTS_DIR" >&2
    exit 1
fi
if [[ ! -f "$PLOTS_DIR/dashboard.html" ]]; then
    echo "Warning: $PLOTS_DIR has no dashboard.html - run scripts/create_plots.sh first?" >&2
fi

if [[ -f "$PLOTS_DIR/dashboard.html" ]]; then
    gzip -k9 -f "$PLOTS_DIR/dashboard.html"
fi

if [[ "${PUBLISH_HTACCESS:-0}" == "1" && ! -f "$PLOTS_DIR/.htaccess" ]]; then
    rewrite_base=$(printf '%s' "$URL_BASE" | sed -E 's|^https?://[^/]+||')/
    cat > "$PLOTS_DIR/.htaccess" <<EOF
<IfModule mod_rewrite.c>
  RewriteEngine On
  RewriteBase $rewrite_base
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
echo "        ($(find "$PLOTS_DIR" -type f | wc -l) files, $(du -sh "$PLOTS_DIR" | cut -f1))"
echo "  to:   $DEST"
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
    rsync --dry-run "${RSYNC_ARGS[@]}" "$PLOTS_DIR/" "$DEST/"
    exit 0
fi
if [[ "${FORCE:-0}" != "1" ]]; then
    read -r -p "This replaces the live site. Continue? [y/N] " reply
    [[ "$reply" == "y" || "$reply" == "Y" ]] || { echo "Aborted."; exit 1; }
fi

rsync "${RSYNC_ARGS[@]}" "$PLOTS_DIR/" "$DEST/"

echo
echo "Published: $URL_BASE/dashboard.html"
