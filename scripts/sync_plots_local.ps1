# Download the current figures from S3DF. PowerShell version.
#
# Run on your own machine, from PowerShell. Asks for your password once.
#
#   .\sync_plots_local.ps1
#   .\sync_plots_local.ps1 -Dest C:\Users\Gregor\talks
#
# If PowerShell refuses to run it:
#   powershell -ExecutionPolicy Bypass -File .\sync_plots_local.ps1
#
# Uses s3dflogin.slac.stanford.edu directly, so no Host alias or ProxyCommand is
# involved. An sdfiana* name would not work here: those resolve only inside
# SLAC and reach you through a ProxyCommand that only your Windows ssh config
# knows about. The login node mounts the same Lustre filesystem, so it sees the
# same files.
#
# Everything comes over a single ssh connection: the remote side stages the
# files under the names we want and streams one tar, which Windows' built-in
# tar.exe unpacks. Otherwise password auth would prompt once per file.

param(
    [string]$SshHost = "s3dflogin.slac.stanford.edu",
    [string]$User    = "gregork",
    [string]$Dest    = ".\fcc_plots",
    [string]$Remote  = "/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks_histograms/IDEA_20260902_seeded",
    [string]$Repo    = "/fs/ddn/sdf/group/atlas/d/gregork/fastsim/FCCJetBenchmarks"
)

$ErrorActionPreference = "Stop"
New-Item -ItemType Directory -Force -Path $Dest | Out-Null
$DestFull = (Resolve-Path $Dest).Path

# Single-quoted here-string: PowerShell must not expand $D, $m and friends.
# Only $Remote and $Repo are substituted, on the two lines below.
$remoteScript = @'
set -eu
D=$(mktemp -d)
trap 'rm -rf "$D"' EXIT
mkdir -p "$D/mh_grids" "$D/decomposition" "$D/event_displays"
cp __REMOTE__/plots/mh_grids/*.pdf       "$D/mh_grids/"       2>/dev/null || true
cp __REMOTE__/plots/event_displays/*.pdf "$D/event_displays/" 2>/dev/null || true
cp __REPO__/doc/*.pdf                    "$D/"                2>/dev/null || true
for m in EEAntiKtR08 EEAntiKtR14 EECambridgeR08 EECambridgeR14 EEKtR08 EEKtR14; do
  for v in "" _by_jets _wide_log _by_jets_wide_log; do
    src="__REMOTE__/plots/mh_decomposition_PF_E_recovery_$m/pflow/mH_decomposition_genphys_common$v.pdf"
    [ -f "$src" ] && cp "$src" "$D/decomposition/PDP_$m$v.pdf"
  done
done
tar cf - -C "$D" .
'@ -replace '__REMOTE__', $Remote -replace '__REPO__', $Repo

$tarPath = Join-Path $env:TEMP "fcc_plots.tar"
Write-Host "Connecting to $SshHost (one password prompt)..." -ForegroundColor Cyan

# Write the tar to a file rather than piping: PowerShell pipelines mangle binary.
$remoteScript | & ssh "$User@$SshHost" "bash -s" | Set-Content -Path $tarPath -Encoding Byte
if ($LASTEXITCODE -ne 0) { throw "ssh failed with exit code $LASTEXITCODE" }

& tar -xf $tarPath -C $DestFull
Remove-Item $tarPath -Force

$n = (Get-ChildItem -Path $DestFull -Recurse -Filter *.pdf).Count
Write-Host ""
Write-Host "Downloaded $n PDFs into $DestFull :" -ForegroundColor Green
Write-Host "  mh_grids\mH_grid_radius_R{04..14}[_Erecovery].pdf       three exponents at fixed R"
Write-Host "  mh_grids\mH_grid_exponent_p{m1,0,pp1}[_Erecovery].pdf   six radii at fixed exponent"
Write-Host "  decomposition\PDP_<algo>R<rr>[_by_jets][_wide_log].pdf  Physics | Detector | Detector+Physics"
Write-Host "  event_displays\<process>_event_displays.pdf             per-event eta-phi displays"
Write-Host "  jet_algorithms.pdf, clustering_algorithms_slides.pdf"
