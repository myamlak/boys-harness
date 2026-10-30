# Seed-spliced sweep, end to end.
#
# Builds one harness translation unit against the pinned boys v5.1.1 tree and
# runs the spliced-ladder comparison.  Run from a Visual Studio developer
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

$Include = Join-Path $BoysTree "include"
$CompileFlags = @(
    "/nologo", "/std:c++latest", "/O2", "/EHsc", "/fp:precise", "/arch:AVX2",
    "/I", $Include
)
$Source = Join-Path $Here "sweep_seed.cpp"
$Object = Join-Path $BuildDir "sweep_seed.obj"
$Exe = Join-Path $BuildDir "sweep_seed.exe"

$FlagLine = (@($CompileFlags) + @($Source) + @($Object, $Exe)) -join " "
& $Python (Join-Path $Here "provenance.py") `
    --boys-tree $BoysTree --harness sweep_seed --harness-dir $Here `
    --flags $FlagLine --out (Join-Path $ResultsDir "provenance.json")
if ($LASTEXITCODE -ne 0) { throw "provenance failed" }

& $Python (Join-Path $Here "../common/reference.py") `
    --boys-tree $BoysTree --grid (Join-Path $DataDir "grid.csv") `
    --out (Join-Path $DataDir "reference.csv") `
    --validation (Join-Path $DataDir "reference_validation.txt") `
    --dps 70 --digits 60
if ($LASTEXITCODE -ne 0) { throw "reference failed" }

# The same grid and generator as the second sweep: the reference must be the
# same bytes, so the two records are read against one reference.
$SourceReference = Join-Path $Here "../amplification_form/data/reference.csv"
$NewHash = (Get-FileHash -Algorithm SHA256 (Join-Path $DataDir "reference.csv")).Hash
$OldHash = (Get-FileHash -Algorithm SHA256 $SourceReference).Hash
Write-Host "reference sha256: $NewHash"
Write-Host "second sweep hash: $OldHash"
if ($NewHash -ne $OldHash) {
    throw "the regenerated reference is not byte-identical to the second sweep's"
}

& cl @CompileFlags "/Fo$Object" "/Fe$Exe" $Source
if ($LASTEXITCODE -ne 0) { throw "compile failed" }

& $Exe --grid (Join-Path $DataDir "grid.csv") `
    --cells (Join-Path $ResultsDir "cells.csv") `
    --admission (Join-Path $ResultsDir "admission.csv") `
    --steps (Join-Path $ResultsDir "replica_steps.csv")
if ($LASTEXITCODE -ne 0) { throw "sweep_seed failed" }

& $Python (Join-Path $Here "analyze.py") `
    --tree $Here --provenance (Join-Path $ResultsDir "provenance.json")
if ($LASTEXITCODE -ne 0) { throw "analysis failed" }

Write-Host ""
Write-Host "seed-spliced sweep complete"
Write-Host "  cells:    $(Join-Path $ResultsDir 'cells.csv')"
Write-Host "  evidence: $(Join-Path $ResultsDir 'summary.txt')"
Write-Host "  table:    $(Join-Path $ResultsDir 'per_order.csv')"
