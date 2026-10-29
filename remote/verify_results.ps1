param([switch]$DryRun)
$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$required = @('results/m3a3_octant_fields_plateau014.npz','results/m3a3_octant_boundary_plateau014.json')
foreach ($rel in $required) {
  $path = Join-Path $repo $rel
  if (-not (Test-Path $path)) { throw "Required output missing: $rel" }
  if ((Get-Item $path).Length -le 0) { throw "Empty output: $rel" }
}
$audit = Get-Content (Join-Path $repo $required[1]) -Raw | ConvertFrom-Json
if ($audit.mesh_tag -ne 'plateau014' -or -not $audit.p2_dofs -or -not $audit.factor_nnz_LU) {
  throw 'Output audit lacks expected mesh/solver fields.'
}
if ($DryRun) { Write-Host 'DRY RUN: required result set appears complete.' }
else { Write-Host 'Required result set validated; this does not validate physics/convergence.' }
