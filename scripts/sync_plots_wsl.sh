#!/usr/bin/env bash
# Download the current figures from S3DF. WSL version.
#
# WSL has its own ~/.ssh, separate from C:\Users\<you>\.ssh, so a Host alias
# that works in PowerShell does not resolve here. This script sidesteps that by
# calling the Windows OpenSSH client (scp.exe), which reads the Windows config.
# Nothing needs to be copied or configured.
#
#   bash scripts/sync_plots_wsl.sh            # into ./fcc_plots
#   bash scripts/sync_plots_wsl.sh -n         # list what would transfer
#   DEST=/mnt/c/Users/Gregor/talks bash scripts/sync_plots_wsl.sh
set -uo pipefail

HOST=${HOST:-sdfiana005}
USER_AT=${USER_AT:-gregork}
REMOTE=${REMOTE:-/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks_histograms/IDEA_20260902_seeded}
REPO=${REPO:-/fs/ddn/sdf/group/atlas/d/gregork/fastsim/FCCJetBenchmarks}
DEST=${DEST:-./fcc_plots}
DRY=0
[[ "${1:-}" == "-n" || "${1:-}" == "--dry-run" ]] && DRY=1

TARGET="$USER_AT@$HOST"

# Prefer the Windows scp so the Windows ssh config (and its ProxyJump) applies.
WINSCP=/mnt/c/Windows/System32/OpenSSH/scp.exe
if [[ -x "$WINSCP" ]]; then
    SCP="$WINSCP"
    WINPATH=1
    echo "using Windows OpenSSH: $WINSCP"
elif command -v scp >/dev/null 2>&1; then
    SCP=scp
    WINPATH=0
    echo "using WSL scp (needs $HOST in ~/.ssh/config or DNS)"
else
    echo "no scp found. Install openssh-client: sudo apt install openssh-client" >&2
    exit 1
fi

mkdir -p "$DEST/mh_grids" "$DEST/decomposition" "$DEST/event_displays"

# scp.exe needs a Windows-style destination path.
topath() { if [[ $WINPATH -eq 1 ]]; then wslpath -w "$1"; else printf '%s' "$1"; fi; }

ok=0; fail=0
get() {  # get <remote path> <local path>
    if [[ $DRY -eq 1 ]]; then echo "  would fetch: $1"; return; fi
    if "$SCP" -q "$TARGET:$1" "$(topath "$2")" 2>/dev/null; then
        ok=$((ok+1))
    else
        fail=$((fail+1)); echo "  failed: $1"
    fi
}

# Higgs-mass grids. Per-radius figures overlay all three exponents
# (p = -1, 0, +1) plus Durham at one radius; per-exponent figures overlay one
# exponent's six radii. Both come plain and as _Erecovery.
echo "== Higgs-mass grids =="
get "$REMOTE/plots/mh_grids/*.pdf" "$DEST/mh_grids"

# Physics | Detector | Detector+Physics. Same filename in every method
# directory, so the method is folded into the local name.
echo "== Physics | Detector | Detector+Physics =="
for m in EEAntiKtR08 EEAntiKtR14 EECambridgeR08 EECambridgeR14 EEKtR08 EEKtR14; do
    for v in "" "_by_jets" "_wide_log" "_by_jets_wide_log"; do
        get "$REMOTE/plots/mh_decomposition_PF_E_recovery_$m/pflow/mH_decomposition_genphys_common$v.pdf" \
            "$DEST/decomposition/PDP_$m$v.pdf"
    done
done

echo "== event displays =="
get "$REMOTE/plots/event_displays/*.pdf" "$DEST/event_displays"

echo "== documents =="
get "$REPO/doc/*.pdf" "$DEST"

echo
if [[ $DRY -eq 1 ]]; then echo "dry run, nothing transferred"; exit 0; fi
echo "$ok transfers ok, $fail failed. Files under $DEST:"
find "$DEST" -name '*.pdf' | sed 's|^|  |' | head -20
n=$(find "$DEST" -name '*.pdf' | wc -l)
[[ $n -gt 20 ]] && echo "  ... and $((n-20)) more"
