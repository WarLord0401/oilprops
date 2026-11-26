<#
.SYNOPSIS
    Bump version, run tests, recompile wheel/sdist, tag, and push to GitHub.
    GitHub Actions (.github/workflows/publish.yml) will automatically publish the tagged release to PyPI.

.EXAMPLE
    .\scripts\release.ps1
    .\scripts\release.ps1 -Bump minor
    .\scripts\release.ps1 -Version 0.2.0
#>
param(
    [ValidateSet("patch", "minor", "major")]
    [string]$Bump = "patch",
    [string]$Version = ""
)

$ErrorActionPreference = "Stop"

$Python = if (Test-Path ".\venv\Scripts\python.exe") { ".\venv\Scripts\python.exe" } else { "python" }

& $Python -c "import setuptools_scm" 2>$null
if ($LASTEXITCODE -ne 0) {
    Write-Host "Installing setuptools-scm into environment..." -ForegroundColor Yellow
    & $Python -m pip install -q setuptools-scm
}

Write-Host "Running pytest..." -ForegroundColor Cyan
& $Python -m pytest -v
if ($LASTEXITCODE -ne 0) { throw "Tests failed. Release aborted." }

$Status = git status --porcelain
if ($Status) {
    throw "Working tree has uncommitted changes. Please commit your changes before releasing."
}

if (-not $Version) {
    $Tags = git tag -l "v*" | Where-Object { $_ -match '^v(\d+)\.(\d+)\.(\d+)$' }
    if (-not $Tags) {
        $Version = "0.1.0"
    } else {
        $Latest = $Tags | Sort-Object {
            $m = [regex]::Match($_, '^v(\d+)\.(\d+)\.(\d+)$')
            [version]::new([int]$m.Groups[1].Value, [int]$m.Groups[2].Value, [int]$m.Groups[3].Value)
        } | Select-Object -Last 1
        $m = [regex]::Match($Latest, '^v(\d+)\.(\d+)\.(\d+)$')
        $maj = [int]$m.Groups[1].Value
        $min = [int]$m.Groups[2].Value
        $pat = [int]$m.Groups[3].Value
        switch ($Bump) {
            "major" { $Version = "$($maj + 1).0.0" }
            "minor" { $Version = "$maj.$($min + 1).0" }
            "patch" { $Version = "$maj.$min.$($pat + 1)" }
        }
    }
}

$Tag = "v$Version"
Write-Host "Preparing release $Tag ..." -ForegroundColor Green

git tag -a $Tag -m "Release $Tag"
if ($LASTEXITCODE -ne 0) { throw "Failed to create git tag $Tag" }

Remove-Item -Recurse -Force -ErrorAction SilentlyContinue dist, build, src\oilprops.egg-info
$env:SETUPTOOLS_SCM_PRETEND_VERSION = $Version
& $Python -m pip install -q --no-build-isolation --no-deps -e .
& $Python -m build --no-isolation
Remove-Item Env:\SETUPTOOLS_SCM_PRETEND_VERSION -ErrorAction SilentlyContinue
if ($LASTEXITCODE -ne 0) { throw "Build failed." }

Write-Host "Built distributions:" -ForegroundColor Cyan
Get-ChildItem dist\ | Select-Object Name, Length

git push origin HEAD $Tag
Write-Host "Pushed commit and $Tag to GitHub. PyPI Trusted Publishing workflow triggered!" -ForegroundColor Green
