# Download the current figures from S3DF. PowerShell version.
#
# Run on your own machine, from PowerShell. Uses the Windows OpenSSH client, so
# it picks up the Host aliases in C:\Users\<you>\.ssh\config - the same config
# that makes `ssh sdfiana005` work. Nothing else needs configuring.
#
#   .\sync_plots_local.ps1
#   .\sync_plots_local.ps1 -Dest C:\Users\Gregor\talks
#
# If PowerShell refuses to run it:
#   powershell -ExecutionPolicy Bypass -File .\sync_plots_local.ps1

param(
    [string]$SshHost = "sdfiana005",
    [string]$User    = "gregork",
    [string]$Dest    = ".\fcc_plots",
    [string]$Remote  = "/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks_histograms/IDEA_20260902_seeded",
    [string]$Repo    = "/fs/ddn/sdf/group/atlas/d/gregork/fastsim/FCCJetBenchmarks"
)

$target = "$User@$SshHost"
$grids  = Join-Path $Dest "mh_grids"
$decomp = Join-Path $Dest "decomposition"
$events = Join-Path $Dest "event_displays"
foreach ($d in @($Dest, $grids, $decomp, $events)) {
    New-Item -ItemType Directory -Force -Path $d | Out-Null
}

function Get-Remote($remotePath, $localPath) {
    # Single-quoted remote path: the remote shell expands any glob, not PowerShell.
    & scp -q "${target}:$remotePath" $localPath
    if ($LASTEXITCODE -ne 0) { Write-Host "   failed: $remotePath" -ForegroundColor DarkYellow }
}

# Higgs-mass grids. The per-radius figures overlay all three exponents
# (p = -1, 0, +1) plus Durham at one radius; the per-exponent figures overlay
# one exponent's six radii. Both come plain and as _Erecovery.
Write-Host "== Higgs-mass grids ==" -ForegroundColor Cyan
Get-Remote "$Remote/plots/mh_grids/*.pdf" $grids

# Physics | Detector | Detector+Physics. The figure has the same name in every
# method directory, so the method goes into the local filename.
Write-Host "== Physics | Detector | Detector+Physics ==" -ForegroundColor Cyan
foreach ($m in @("EEAntiKtR08","EEAntiKtR14","EECambridgeR08",
                 "EECambridgeR14","EEKtR08","EEKtR14")) {
    foreach ($v in @("","_by_jets","_wide_log","_by_jets_wide_log")) {
        $src = "$Remote/plots/mh_decomposition_PF_E_recovery_$m/pflow/mH_decomposition_genphys_common$v.pdf"
        Get-Remote $src (Join-Path $decomp "PDP_$m$v.pdf")
    }
}

Write-Host "== event displays ==" -ForegroundColor Cyan
Get-Remote "$Remote/plots/event_displays/*.pdf" $events

Write-Host "== documents ==" -ForegroundColor Cyan
Get-Remote "$Repo/doc/*.pdf" $Dest

Write-Host ""
Write-Host "Done. Under $Dest :" -ForegroundColor Green
Write-Host "  mh_grids\mH_grid_radius_R{04..14}[_Erecovery].pdf       three exponents at fixed R"
Write-Host "  mh_grids\mH_grid_exponent_p{m1,0,pp1}[_Erecovery].pdf   six radii at fixed exponent"
Write-Host "  decomposition\PDP_<algo>R<rr>[_by_jets][_wide_log].pdf  Physics | Detector | Detector+Physics"
Write-Host "  event_displays\<process>_event_displays.pdf             per-event eta-phi displays"
Write-Host "  jet_algorithms.pdf, clustering_algorithms_slides.pdf"
