#!/usr/bin/env bash
# Download the current figures from S3DF into $DEST (default ./fcc_plots).
#   bash sync_plots_local.sh              # into ./fcc_plots
#   DEST=~/talks bash sync_plots_local.sh
set -euo pipefail

HOST=${HOST:-gregork@s3dflogin.slac.stanford.edu}
REMOTE=${REMOTE:-/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks_histograms/IDEA_20260902_seeded_v2}
REPO=${REPO:-/fs/ddn/sdf/group/atlas/d/gregork/fastsim/FCCJetBenchmarks}
DEST=${DEST:-./fcc_plots}

export MSYS_NO_PATHCONV=1 MSYS2_ARG_CONV_EXCL='*'

mkdir -p "$DEST"

# ssh offers every key it can find before falling back to a password, and the
# server drops the connection after about six attempts ("Too many
# authentication failures"). So pin exactly one method: the S3DF key if it is
# there, otherwise password only with key offering switched off.
if [[ -z "${SSH_OPTS:-}" ]]; then
    if [[ -f "$HOME/.ssh/s3df/id_ed25519" ]]; then
        SSH_OPTS="-i $HOME/.ssh/s3df/id_ed25519 -o IdentitiesOnly=yes"
    else
        SSH_OPTS="-o PubkeyAuthentication=no -o PreferredAuthentications=keyboard-interactive,password"
    fi
fi

# Quoted heredoc: the remote's own $D and $m must survive. The two placeholders
# are filled in locally, just below.
remote_script=$(cat <<'REMOTE_SCRIPT'
set -eu
D=$(mktemp -d)
trap 'rm -rf "$D"' EXIT
mkdir -p "$D/mh_grids" "$D/decomposition" "$D/event_displays"

cp __REMOTE__/plots/mh_grids/*.pdf       "$D/mh_grids/"       2>/dev/null || true
cp __REMOTE__/plots/event_displays/*.pdf "$D/event_displays/" 2>/dev/null || true
cp __REPO__/doc/*.pdf                    "$D/"                2>/dev/null || true

# The decomposition figure has the same filename in every method directory, so
# rename each copy to carry the method before taring.
for m in EEAntiKtR08 EEAntiKtR14 EECambridgeR08 EECambridgeR14 EEKtR08 EEKtR14; do
  for v in "" _by_jets _wide_log _by_jets_wide_log; do
    src="__REMOTE__/plots/mh_decomposition_PF_E_recovery_$m/pflow/mH_decomposition_genphys_common$v.pdf"
    [ -f "$src" ] && cp "$src" "$D/decomposition/PDP_$m$v.pdf"
  done
done

echo "remote staged $(find "$D" -name '*.pdf' | wc -l) PDFs" >&2
tar cf - -C "$D" . | base64
REMOTE_SCRIPT
)
remote_script=${remote_script//__REMOTE__/$REMOTE}
remote_script=${remote_script//__REPO__/$REPO}

tmp=$(mktemp)
trap 'rm -f "$tmp"' EXIT

printf '%s\n' "$remote_script" | ssh $SSH_OPTS "$HOST" "bash -s" > "$tmp"

# A login banner on stdout would sit in front of the payload.
grep -E '^[A-Za-z0-9+/=]+$' "$tmp" | base64 -d | tar xf - -C "$DEST"

n=$(find "$DEST" -name '*.pdf' | wc -l)
if [[ "$n" -eq 0 ]]; then
    echo "The remote staged nothing. Check that these paths exist on $HOST:" >&2
    echo "  $REMOTE/plots/mh_grids" >&2
    echo "  $REPO/doc" >&2
    exit 1
fi

echo
echo "Downloaded $n PDFs into $DEST:"
echo "  mh_grids/mH_grid_radius_R{04..14}[_Erecovery].pdf       three exponents at fixed R"
echo "  mh_grids/mH_grid_exponent_p{m1,0,pp1}[_Erecovery].pdf   six radii at fixed exponent"
echo "  decomposition/PDP_<algo>R<rr>[_by_jets][_wide_log].pdf  Physics | Detector | Detector+Physics"
echo "  event_displays/<process>_event_displays.pdf             per-event eta-phi displays"
echo "  event_displays/<process>_event_displays_simple.pdf      the same, particles coloured by jet"
echo "  jet_algorithms.pdf, clustering_algorithms_slides.pdf"
