param([Parameter(Mandatory=$true)][string]$Destination, [switch]$DryRun)
$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
& (Join-Path $PSScriptRoot 'verify_results.ps1')
$dest = [IO.Path]::GetFullPath($Destination)
if ($DryRun) { Write-Host "DRY RUN: would copy solver outputs and audit JSON into $dest"; exit 0 }
New-Item -ItemType Directory -Force -Path $dest | Out-Null
Copy-Item -LiteralPath (Join-Path $repo 'results/m3a3_octant_fields_plateau014.npz') -Destination $dest
Copy-Item -LiteralPath (Join-Path $repo 'results/m3a3_octant_boundary_plateau014.json') -Destination $dest
Get-ChildItem -File $dest | Get-FileHash -Algorithm SHA256 | Select-Object Path,Hash | ConvertTo-Json | Set-Content (Join-Path $dest 'SHA256SUMS.json')
Write-Host "Collected outputs with SHA-256 manifest in $dest"
