# The v5.2.0 execution check, from the runbooks' own commands.
#
# Compiles sweep1.cpp and sweep2.cpp with the compile line documented in
# sweeps/product_amplification/SWEEP.md and sweeps/amplification_form/SWEEP.md,
# substituting the v5.2.0 include tree for the v5.1.1 one (as the v5.1.1-check
# did in the other direction), then runs them with the argument lists in the
# runbooks, then each analyze.py with its runbook arguments, then the five
# table scripts.  The grids, fits and 60-digit references are the deposited
# ones, copied into this scratch tree: the check re-runs the evaluation, not
# the reference generation.
#
# The boys identity was established by verify_pub_tree.py and is recorded in
# pub-tree-identity.txt; no git command is run here.

$ErrorActionPreference = 'Stop'
$Here = Split-Path -Parent $MyInvocation.MyCommand.Path
$Boys = 'C:\Users\myaml\dev\_boys-v520-pub'
$Log = Join-Path $Here 'logs\execution-check.log'
$Logs = Join-Path $Here 'logs'
$DevShell = 'C:\Program Files\Microsoft Visual Studio\18\Community\Common7\Tools\Launch-VsDevShell.ps1'
$S1Args = @('/nologo', '/std:c++latest', '/O2', '/EHsc', '/fp:precise',
            '/arch:AVX2', '/I', "$Boys\include",
            '/Fobuild\sweep1.obj', '/Febuild\sweep1.exe', 'sweep1.cpp')
$S2Args = @('/nologo', '/std:c++latest', '/O2', '/EHsc', '/fp:precise',
            '/arch:AVX2', '/I', "$Boys\include",
            '/Fobuild\sweep2.obj', '/Febuild\sweep2.exe', 'sweep2.cpp')

$S1 = Join-Path $Here 'sweeps\product_amplification'
$S2 = Join-Path $Here 'sweeps\amplification_form'
$L1 = Join-Path $Here 'logs\sweep1-run.log'
$L1a = Join-Path $Here 'logs\sweep1-analyze.log'
$L2 = Join-Path $Here 'logs\sweep2-run.log'
$L2a = Join-Path $Here 'logs\sweep2-analyze.log'

function Step([string]$name, [scriptblock]$block) {
    $sw = [Diagnostics.Stopwatch]::StartNew()
    & $block
    $code = $LASTEXITCODE
    $line = "{0,-20} EXIT={1} SECONDS={2}" -f $name, $code, [math]::Round($sw.Elapsed.TotalSeconds, 1)
    Add-Content -Path $Log -Value $line
    Write-Output $line
    if ($code -ne 0) { throw "$name exited $code" }
}

$identity = Get-Content (Join-Path $Here 'pub-tree-identity.txt') -Raw
$TagName = ([regex]::Match($identity, 'refs/tags/(v[0-9.]+)')).Groups[1].Value
$Commit = ([regex]::Match($identity, 'tag target commit\s+([0-9a-f]{40})')).Groups[1].Value
if (-not $TagName -or -not $Commit) {
    throw 'pub-tree-identity.txt does not carry the tag and its commit'
}
$env:BOYS_TAG_AT_HEAD = $TagName
$env:BOYS_DESCRIBE = $TagName
$env:BOYS_COMMIT = $Commit
$env:BOYS_DIRTY = '0'

function Enter-DevShell {
    if (-not (Get-Command cl.exe -ErrorAction SilentlyContinue)) {
        & $DevShell -Arch amd64 -HostArch amd64 -SkipAutomaticLocation | Out-Null
    }
    (Get-Command cl.exe -ErrorAction Stop).Source |
        Set-Content -Path (Join-Path $Logs 'compiler-path.log') -Encoding utf8
    cmd /c 'cl 2>&1' | Select-Object -First 2 |
        Set-Content -Path (Join-Path $Logs 'compiler-version.log') -Encoding utf8
}

# The dev shell comes first: the provenance recorder reads the compiler's
# banner and environment, and must see the same compiler the builds use.
Enter-DevShell

Set-Content -Path $Log -Value "v5.2.0 execution check, runbook commands" -Encoding utf8

Step 'sweep1-provenance' {
    Set-Location $S1
    & python (Join-Path $S1 '..\common\boysprovenance.py') --boys-tree $Boys `
        --harness sweep1 --harness-dir $S1 --flags ($S1Args -join ' ') `
        --out (Join-Path $S1 'results\provenance.json') *> (Join-Path $Logs 'sweep1-provenance.log')
}

Step 'sweep1-build' {
    Enter-DevShell
    Set-Location $S1
    & cl $S1Args *> (Join-Path $Logs 'sweep1-build.log')
}

Step 'sweep1-run' {
    Set-Location $S1
    & .\build\sweep1.exe --grid (Join-Path $S1 'data\grid.csv') `
        --cells (Join-Path $S1 'results\cells.csv') `
        --admission (Join-Path $S1 'results\admission.csv') `
        --gain (Join-Path $S1 'results\replica_gain.csv') --batch full *> $L1
}

Step 'sweep1-analyze' {
    Set-Location $S1
    & python analyze.py --tree $S1 --boys-tree $Boys `
        --provenance (Join-Path $S1 'results\provenance.json') `
        --validation (Join-Path $S1 'data\reference_validation.txt') *> $L1a
}

Step 'sweep2-provenance' {
    Set-Location $S2
    & python (Join-Path $S2 '..\common\boysprovenance.py') --boys-tree $Boys `
        --harness sweep2 --harness-dir $S2 --flags ($S2Args -join ' ') `
        --out (Join-Path $S2 'results\provenance.json') *> (Join-Path $Logs 'sweep2-provenance.log')
}

Step 'sweep2-build' {
    Enter-DevShell
    Set-Location $S2
    & cl $S2Args *> (Join-Path $Logs 'sweep2-build.log')
}

Step 'sweep2-run' {
    Set-Location $S2
    & .\build\sweep2.exe --grid (Join-Path $S2 'data\grid.csv') `
        --fits (Join-Path $S2 'data\fits.csv') `
        --paths (Join-Path $S2 'results\paths.csv') `
        --cells (Join-Path $S2 'results\cells.csv') `
        --admission (Join-Path $S2 'results\admission.csv') *> $L2
}

Step 'sweep2-analyze' {
    Set-Location $S2
    & python analyze.py --tree $S2 --boys-tree $Boys `
        --provenance (Join-Path $S2 'results\provenance.json') `
        --validation (Join-Path $S2 'data\reference_validation.txt') *> $L2a
}

Step 'tables' {
    Set-Location $Here
    $env:BOYS_SWEEP_TREE = $S1
    & python make_tables.py *> (Join-Path $Here 'logs\tables-make_tables.log')
}

Step 'merge-per-order' {
    Set-Location $Here
    & python merge_per_order.py *> (Join-Path $Here 'logs\tables-merge_per_order.log')
}

Step 'endtoend' {
    Set-Location $Here
    & python make_endtoend.py *> (Join-Path $Here 'logs\tables-make_endtoend.log')
}

Step 'frontier' {
    Set-Location $Here
    & python make_frontier.py *> (Join-Path $Here 'logs\tables-make_frontier.log')
}

Step 'prediction-record' {
    Set-Location $Here
    & python make_prediction_record.py *> (Join-Path $Here 'logs\tables-make_prediction_record.log')
}

Write-Output "FINISH $(Get-Date -Format o)"
