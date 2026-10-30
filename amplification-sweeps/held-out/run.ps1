# Held-out sweep, end to end.
#
# Builds one harness translation unit against the pinned boys v5.1.1 tree,
# evaluates both entries on the held-out grid (rule H), and scores the paper's
# frozen prediction rule on the result.  Run from a Visual Studio developer
# shell, without python -I (mpmath lives in the user site).
#
# This sweep measures no duration, so it has no canary; its figures are
# deterministic arithmetic from a fixed argument list.
#
# The pipeline does not run git.  The kernel revision is the one the second
# sweep's provenance record was taken from (v5.1.1, 834bdbb); provenance.py
# proves the tree is that revision by matching every kernel file digest
# against that record, and fails if any differs.
#
#   powershell -ExecutionPolicy Bypass -File run.ps1

#Requires -Version 5.1
[CmdletBinding()]
param(
    [string]$BoysTree = $(if ($env:BOYS_TREE) { $env:BOYS_TREE } else { "C:\Users\myaml\dev\_boys-v511" }),
    [string]$Python = "python"
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$Here = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $Here

if (-not (Test-Path (Join-Path $BoysTree "include/boys/version.hpp"))) {
    throw "not a boys tree: $BoysTree"
}

$BuildDir = Join-Path $Here "build"
$ResultsDir = Join-Path $Here "results"
$DataDir = Join-Path $Here "data"
New-Item -ItemType Directory -Force -Path $BuildDir, $ResultsDir, $DataDir | Out-Null

# The held-out argument list.  grid.py asserts exact-double disjointness from
# both earlier sweeps and fails the run on any collision.
& $Python (Join-Path $Here "grid.py") `
    --boys-tree $BoysTree --out (Join-Path $DataDir "grid.csv")
if ($LASTEXITCODE -ne 0) { throw "grid failed" }

$Include = Join-Path $BoysTree "include"
$CompileFlags = @(
    "/nologo", "/std:c++latest", "/O2", "/EHsc", "/fp:precise", "/arch:AVX2",
    "/I", $Include
)
$Source = Join-Path $Here "sweep_heldout.cpp"
$Object = Join-Path $BuildDir "sweep_heldout.obj"
$Exe = Join-Path $BuildDir "sweep_heldout.exe"

$FlagLine = (@($CompileFlags) + @($Source) + @($Object, $Exe)) -join " "
& $Python (Join-Path $Here "provenance.py") `
    --boys-tree $BoysTree --harness sweep_heldout --harness-dir $Here `
    --flags $FlagLine --out (Join-Path $ResultsDir "provenance.json")
if ($LASTEXITCODE -ne 0) { throw "provenance failed" }

& $Python (Join-Path $Here "../common/reference.py") `
    --boys-tree $BoysTree --grid (Join-Path $DataDir "grid.csv") `
    --out (Join-Path $DataDir "reference.csv") `
    --validation (Join-Path $DataDir "reference_validation.txt") `
    --dps 70 --digits 60
if ($LASTEXITCODE -ne 0) { throw "reference failed" }

& cl @CompileFlags "/Fo$Object" "/Fe$Exe" $Source
if ($LASTEXITCODE -ne 0) { throw "compile failed" }

& $Exe --grid (Join-Path $DataDir "grid.csv") `
    --cells (Join-Path $ResultsDir "cells.csv") `
    --admission (Join-Path $ResultsDir "admission.csv")
if ($LASTEXITCODE -ne 0) { throw "sweep_heldout failed" }

& $Python (Join-Path $Here "analyze.py") `
    --tree $Here --provenance (Join-Path $ResultsDir "provenance.json")
if ($LASTEXITCODE -ne 0) { throw "analysis failed" }

Write-Host ""
Write-Host "held-out sweep complete"
Write-Host "  cells:    $(Join-Path $ResultsDir 'cells.csv')"
Write-Host "  evidence: $(Join-Path $ResultsDir 'summary.txt')"
Write-Host "  table:    $(Join-Path $ResultsDir 'prediction_heldout.csv')"
