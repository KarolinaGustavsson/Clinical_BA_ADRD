# =====================================================================
# Windows / PowerShell driver: one fresh R process per outcome.
# =====================================================================
#
# The VDI is Windows, so this replaces run_fg_project1.sh (which needs bash).
#
# Each outcome gets its own R process, so memory is fully released between
# outcomes and a failure on one outcome does not stop the rest. Results
# accumulate in the same CSVs and finished models are skipped on restart.
#
#   .\run_fg_project1.ps1                 # everything
#   .\run_fg_project1.ps1 -Which main     # main analysis only
#   .\run_fg_project1.ps1 -Which svedem   # SveDem analysis only
#   .\run_fg_project1.ps1 -Which test     # just the test suite
#
# If Rscript.exe is not found automatically, pass it explicitly:
#   .\run_fg_project1.ps1 -RscriptPath "C:\Program Files\R\R-4.4.1\bin\Rscript.exe"
#
# If PowerShell refuses to run this file ("running scripts is disabled"),
# either unblock it for this session:
#   Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
# or run the commands by hand - see README_fg_project1.md.

param(
    [ValidateSet("all", "main", "svedem", "test")]
    [string]$Which = "all",
    [string]$RscriptPath = ""
)

$ErrorActionPreference = "Continue"

# --- locate Rscript.exe ------------------------------------------------
function Find-Rscript {
    # Already on PATH?
    $onPath = Get-Command Rscript.exe -ErrorAction SilentlyContinue
    if ($onPath) { return $onPath.Source }

    # Usual install locations, newest version first.
    $candidates = @(
        "C:\Program Files\R\*\bin\x64\Rscript.exe",
        "C:\Program Files\R\*\bin\Rscript.exe",
        "C:\Program Files (x86)\R\*\bin\Rscript.exe",
        "$env:LOCALAPPDATA\Programs\R\*\bin\x64\Rscript.exe",
        "$env:LOCALAPPDATA\Programs\R\*\bin\Rscript.exe"
    )
    foreach ($pattern in $candidates) {
        $hit = Get-ChildItem $pattern -ErrorAction SilentlyContinue |
               Sort-Object FullName -Descending | Select-Object -First 1
        if ($hit) { return $hit.FullName }
    }
    return $null
}

if ($RscriptPath -ne "") {
    $R = $RscriptPath
} else {
    $R = Find-Rscript
}

if (-not $R -or -not (Test-Path $R)) {
    Write-Host "Could not find Rscript.exe." -ForegroundColor Red
    Write-Host "Locate it with:"
    Write-Host '  Get-ChildItem C:\ -Filter Rscript.exe -Recurse -ErrorAction SilentlyContinue | Select-Object -First 5 FullName'
    Write-Host "then re-run:"
    Write-Host '  .\run_fg_project1.ps1 -RscriptPath "C:\path\to\Rscript.exe"'
    exit 1
}

Write-Host "Using R: $R" -ForegroundColor Cyan
Write-Host "Folder : $(Get-Location)" -ForegroundColor Cyan

$MainOutcomes   = @("ADRD", "Alzheimers", "Vascular", "Other")
$SveDemOutcomes = @("S_ADRD", "S_e_alz", "S_l_alz", "S_alz_vas", "S_vas",
                    "S_lew", "S_fro", "S_par", "S_uns", "S_oth")

$failed = @()

# $ScriptArgs is the full argument list passed to Rscript, e.g.
# @("fg_project1.R", "main", "ADRD") or @("test_fg_project1.R").
function Invoke-Step {
    param([string]$Label, [string[]]$ScriptArgs)

    Write-Host ""
    Write-Host "==============================================" -ForegroundColor DarkGray
    Write-Host ">>> $Label   ($(Get-Date -Format 'yyyy-MM-dd HH:mm:ss'))"
    Write-Host "==============================================" -ForegroundColor DarkGray

    & $R @ScriptArgs

    if ($LASTEXITCODE -ne 0) {
        Write-Host "!!! FAILED: $Label" -ForegroundColor Red
        $script:failed += $Label
    }
}

if ($Which -eq "test") {
    Invoke-Step "test_fg_project1.R" @("test_fg_project1.R")
} else {
    if ($Which -eq "all" -or $Which -eq "main") {
        foreach ($o in $MainOutcomes) {
            Invoke-Step "fg_project1.R main $o" @("fg_project1.R", "main", $o)
        }
    }
    if ($Which -eq "all" -or $Which -eq "svedem") {
        foreach ($o in $SveDemOutcomes) {
            Invoke-Step "fg_project1.R svedem $o" @("fg_project1.R", "svedem", $o)
        }
    }
}

Write-Host ""
Write-Host "==============================================" -ForegroundColor DarkGray
if ($failed.Count -gt 0) {
    Write-Host "Finished with failures:" -ForegroundColor Yellow
    $failed | ForEach-Object { Write-Host "  $_" -ForegroundColor Yellow }
    Write-Host "Re-run this script to retry them; completed models are skipped."
    exit 1
}
Write-Host "All outcomes completed." -ForegroundColor Green
