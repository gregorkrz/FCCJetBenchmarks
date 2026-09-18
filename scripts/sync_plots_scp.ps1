# Download the current figures from S3DF using scp. PowerShell version.
#
# Run on your own machine, from PowerShell:
#
#   .\sync_plots_scp.ps1
#   .\sync_plots_scp.ps1 -SshHost s3df -Dest C:\Users\Gregor\talks
#
# If PowerShell refuses to run it:
#   powershell -ExecutionPolicy Bypass -File .\sync_plots_scp.ps1
#
# Unlike sync_plots_local.ps1 this goes through whatever Host alias your
# ~/.ssh/config defines, so an internal name such as sdfiana005 works: Windows'
# own ssh.exe reads the config and applies the ProxyCommand. Nothing is piped:
# the remote builds one tar, scp fetches that single file, tar.exe unpacks it.
# Two password prompts, one for ssh and one for scp.
#
# The remote paths are substituted into the script text, which travels over
# stdin. They must not become ssh arguments: under Git Bash MSYS would rewrite
# a leading slash into a Windows path, and every remote cp would fail silently.

param(
    [string]$SshHost   = "sdfiana005",
    [string]$Dest      = ".\fcc_plots",
    [string]$RemoteTar = "/tmp/fcc_plots_sync.tar",
    [string]$Remote    = "/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks_histograms/IDEA_20260902_seeded",
    [string]$Repo      = "/fs/ddn/sdf/group/atlas/d/gregork/fastsim/FCCJetBenchmarks"
)

$ErrorActionPreference = "Stop"
New-Item -ItemType Directory -Force -Path $Dest | Out-Null
$DestFull = (Resolve-Path $Dest).Path

# Single-quoted here-string: PowerShell must not expand $D, $m and friends.
$remoteScript = @'
set -eu
rm -f __TAR__
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

n=$(find "$D" -name '*.pdf' | wc -l)
[ "$n" -gt 0 ] || { echo "staged nothing, check the remote paths" >&2; exit 1; }
tar cf __TAR__ -C "$D" .
echo "staged $n PDFs"
'@ -replace '__REMOTE__', $Remote -replace '__REPO__', $Repo -replace '__TAR__', $RemoteTar

Write-Host "1/3 building the archive on $SshHost ..." -ForegroundColor Cyan
$remoteScript | & ssh $SshHost "bash -s"
if ($LASTEXITCODE -ne 0) { throw "remote staging failed with exit code $LASTEXITCODE" }

$tarPath = Join-Path $env:TEMP "fcc_plots.tar"
Write-Host "2/3 copying it down ..." -ForegroundColor Cyan
& scp "${SshHost}:${RemoteTar}" $tarPath
if ($LASTEXITCODE -ne 0) { throw "scp failed with exit code $LASTEXITCODE" }

Write-Host "3/3 unpacking ..." -ForegroundColor Cyan
& tar -xf $tarPath -C $DestFull
if ($LASTEXITCODE -ne 0) { throw "tar failed with exit code $LASTEXITCODE" }
Remove-Item $tarPath -Force

$n = (Get-ChildItem -Path $DestFull -Recurse -Filter *.pdf).Count
if ($n -eq 0) { throw "archive unpacked but contains no PDFs" }
Write-Host ""
Write-Host "Downloaded $n PDFs into $DestFull :" -ForegroundColor Green
Write-Host "  mh_grids\mH_grid_radius_R{04..14}[_Erecovery].pdf       three exponents at fixed R"
Write-Host "  mh_grids\mH_grid_exponent_p{m1,0,pp1}[_Erecovery].pdf   six radii at fixed exponent"
Write-Host "  decomposition\PDP_<algo>R<rr>[_by_jets][_wide_log].pdf  Physics | Detector | Detector+Physics"
Write-Host "  event_displays\<process>_event_displays.pdf             per-event eta-phi displays"
Write-Host "  jet_algorithms.pdf, clustering_algorithms_slides.pdf"
Write-Host ""
Write-Host "The remote copy stays at ${SshHost}:${RemoteTar} and is overwritten on the next run."
