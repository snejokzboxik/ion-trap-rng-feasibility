param([switch]$DryRun)
$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$py = Get-Command py -ErrorAction SilentlyContinue
if (-not $py) { throw 'Install Python 3.12 x64 and ensure py.exe is on PATH, then rerun.' }
& py -3.12 --version
if ($LASTEXITCODE -ne 0) { throw 'Python 3.12 not detected.' }
$venv = Join-Path $repo '.venv'
if ($DryRun) { Write-Host "DRY RUN: would create/use $venv and install requirements-lock.txt"; exit 0 }
if (-not (Test-Path (Join-Path $venv 'Scripts/python.exe'))) { & py -3.12 -m venv $venv }
& (Join-Path $venv 'Scripts/python.exe') -m pip install -r (Join-Path $repo 'requirements-lock.txt')
if ($LASTEXITCODE -ne 0) { throw 'Locked dependency installation failed' }
# Tests use tests/conftest.py to expose the src-layout package; avoiding editable
# installation avoids unpinned build-backend downloads.
& (Join-Path $venv 'Scripts/python.exe') -m pytest -q
if ($LASTEXITCODE -ne 0) { throw 'Regression suite failed; no FEM should be launched.' }
