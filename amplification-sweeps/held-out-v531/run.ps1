# Held-out exact-tag replay (boys v5.3.1), end to end.
#
# Rebuilds the held-out sweep against the exact v5.3.1 tag tree and reruns it
# on the deposited 1300-argument grid, then compares every output with the
# v5.1.1 record in ../held-out.  The reference is the deposited one: the
# reference chain is byte-identical between the two kernels, so it is copied
# and its digest checked rather than regenerated.  Run from a Visual Studio
# developer shell's cmd through build.bat; the PowerShell steps only need
# python on PATH.
#
# This sweep measures no duration, so it has no canary; its figures are
# deterministic arithmetic from a fixed argument list.
#
# The kernel revision is proved by content, not by git: kernel_record_v531.py
# compares all 76 evaluated files against the v5.3.1 digest column of the
# certified-bound v531-check record (the column the exact-tag replay proved
# against the tag tree at harness v1.0.35), and the three reference-chain
# files against the second sweep's record.
#
#   powershell -ExecutionPolicy Bypass -File run.ps1

#Requires -Version 5.1
[CmdletBinding()]
param(
    [string]$BoysTree = $env:BOYS_TREE,
    [string]$Python = "python",
    [string]$Digests,
    [string]$Previous,
    [string]$ReferenceSource,
    [string]$Deposited
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

Set-Location $PSScriptRoot

# $PSScriptRoot is not set while a param default is evaluated in Windows
# PowerShell 5.1, so the paths are resolved here.
if (-not $BoysTree) { throw "set BOYS_TREE to the v5.3.1 tag checkout, or pass -BoysTree" }
if (-not $Digests) { $Digests = Join-Path $PSScriptRoot "..\certified-bound\v531-check\digests.txt" }
if (-not $Previous) { $Previous = Join-Path $PSScriptRoot "..\amplification-form\results\provenance.json" }
if (-not $ReferenceSource) { $ReferenceSource = Join-Path $PSScriptRoot "..\held-out\data\reference.csv" }
if (-not $Deposited) { $Deposited = Join-Path $PSScriptRoot "..\held-out\results" }

$env:BOYS_TREE = $BoysTree

if (-not (Test-Path (Join-Path $BoysTree "include/boys/version.hpp"))) {
    throw "not a boys tree: $BoysTree"
}

New-Item -ItemType Directory -Force -Path build, data, results | Out-Null

& $Python (Join-Path $PSScriptRoot "kernel_record_v531.py") `
    --tree $BoysTree --v531-digests $Digests --previous $Previous `
    --out (Join-Path $PSScriptRoot "results/kernel_record_v531.json")
if ($LASTEXITCODE -ne 0) { throw "kernel record failed" }

& $Python (Join-Path $PSScriptRoot "grid.py") `
    --boys-tree $BoysTree --out (Join-Path $PSScriptRoot "data/grid.csv")
if ($LASTEXITCODE -ne 0) { throw "grid failed" }

Copy-Item $ReferenceSource (Join-Path $PSScriptRoot "data/reference.csv") -Force
Copy-Item (Join-Path (Split-Path $ReferenceSource) "reference_validation.txt") `
    (Join-Path $PSScriptRoot "data/reference_validation.txt") -Force
$refHash = (Get-FileHash (Join-Path $PSScriptRoot "data/reference.csv") -Algorithm SHA256).Hash
if ($refHash -ne (Get-FileHash $ReferenceSource -Algorithm SHA256).Hash) {
    throw "the copied reference does not match the deposited one"
}
Write-Host "reference: $($refHash.ToLower())"

$Include = Join-Path $BoysTree "include"
$Flags = "/nologo /std:c++latest /O2 /EHsc /fp:precise /arch:AVX2 /I $Include sweep_heldout.cpp build\sweep_heldout.obj build\sweep_heldout.exe"
& $Python (Join-Path $PSScriptRoot "provenance.py") `
    --boys-tree $BoysTree --kernel-record (Join-Path $PSScriptRoot "results/kernel_record_v531.json") `
    --harness sweep_heldout --harness-dir $PSScriptRoot --flags $Flags `
    --out (Join-Path $PSScriptRoot "results/provenance.json")
if ($LASTEXITCODE -ne 0) { throw "provenance failed" }

& cmd /c "`"$(Join-Path $PSScriptRoot 'build.bat')`""
if ($LASTEXITCODE -ne 0) { throw "compile failed" }

& (Join-Path $PSScriptRoot "build/sweep_heldout.exe") `
    --grid (Join-Path $PSScriptRoot "data/grid.csv") `
    --cells (Join-Path $PSScriptRoot "results/cells.csv") `
    --admission (Join-Path $PSScriptRoot "results/admission.csv")
if ($LASTEXITCODE -ne 0) { throw "sweep_heldout failed" }

& $Python (Join-Path $PSScriptRoot "analyze.py") `
    --tree $PSScriptRoot --provenance (Join-Path $PSScriptRoot "results/provenance.json")
if ($LASTEXITCODE -ne 0) { throw "analysis failed" }

# Compare every value-bearing output with the v5.1.1 deposit.  summary.txt
# carries the run's own paths on two lines; every other line must agree.
$lines = @()
foreach ($name in "cells.csv", "admission.csv", "prediction_heldout.csv",
                  "band_summary.csv", "per_order.csv") {
    $a = Get-FileHash (Join-Path $PSScriptRoot "results/$name") -Algorithm SHA256
    $b = Get-FileHash (Join-Path $Deposited $name) -Algorithm SHA256
    $state = if ($a.Hash -eq $b.Hash) { "identical" } else { "DIFFERS" }
    $lines += "$name`t$state`t$($a.Hash.ToLower())"
}
$mine = Get-Content (Join-Path $PSScriptRoot "results/summary.txt")
$theirs = Get-Content (Join-Path $Deposited "summary.txt")
$diff = Compare-Object $theirs $mine
$summaryState = if (-not $diff) { "identical" } else { "$($diff.Count) differing lines" }
$lines += "summary.txt`t$summaryState"
if ($diff) {
    $diff | ForEach-Object { $lines += "  $($_.SideIndicator) $($_.InputObject)" }
}
$lines | Set-Content (Join-Path $PSScriptRoot "results/comparison.txt") -Encoding utf8
$lines | ForEach-Object { Write-Host $_ }
Write-Host ""
Write-Host "held-out exact-tag replay complete"
