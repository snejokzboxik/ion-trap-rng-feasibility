param(
  [string]$ResultsDirectory,
  [switch]$DryRun
)
$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$python = Join-Path $repo '.venv/Scripts/python.exe'
$validator = Join-Path $repo 'scripts/verify_remote_outputs.py'
if (-not (Test-Path -LiteralPath $python)) { Write-Output 'Verification failed: project Python environment is missing.'; exit 2 }
if (-not (Test-Path -LiteralPath $validator)) { Write-Output 'Verification failed: validator script is missing.'; exit 2 }
$pyArgs = @($validator, '--repo', $repo)
if ($ResultsDirectory) { $pyArgs += @('--output-dir', [IO.Path]::GetFullPath($ResultsDirectory)) }
$output = & $python @pyArgs 2>&1
$code = $LASTEXITCODE
$output | ForEach-Object { Write-Output $_ }
if ($code -ne 0) { Write-Output 'Output verification failed. File integrity is separate from numerical convergence.'; exit 2 }
Write-Output 'Output files and P2 dimensions verified; numerical force convergence is not assessed.'
exit 0
