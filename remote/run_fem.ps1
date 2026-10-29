param([switch]$DryRun, [switch]$Execute)
$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
& (Join-Path $PSScriptRoot 'preflight.ps1') -DryRun
if ($LASTEXITCODE -ne 0) { throw 'Preflight did not pass.' }
$venvPy = Join-Path $repo '.venv/Scripts/python.exe'
if (-not (Test-Path $venvPy)) { throw 'Run remote/setup_windows.ps1 first.' }
$tag = 'plateau014'
$args = @('scripts/m3a3_octant_p2.py', '--tag', $tag, '--solve')
$log = Join-Path $repo ('results/r0_fem_' + (Get-Date -Format 'yyyyMMdd_HHmmss') + '.log')
if ($DryRun -or -not $Execute) {
  Write-Host "DRY RUN: $venvPy $($args -join ' ')"
  Write-Host "Log would be: $log"
  Write-Host 'No solver launched. Explicit -Execute and prior project-stage approval are required.'
  exit 0
}
# Guardrail only for this single bounded solve. It can stop its own Python process
# if observed working set or system free RAM crosses configured bounds. It never
# terminates another job or changes process priorities. SuperLU is not checkpointable mid-LU.
$cfg = Get-Content (Join-Path $repo 'configs/r0_remote_fem.json') -Raw | ConvertFrom-Json
$os = Get-CimInstance Win32_OperatingSystem
$freeGiB = $os.FreePhysicalMemory / 1MB
if ($freeGiB -lt $cfg.minimum_free_memory_gib) { throw 'Available RAM fell below configured preflight threshold; solve not started.' }
Push-Location $repo
try {
  $stdout = "$log.stdout"
  $stderr = "$log.stderr"
  $proc = Start-Process -FilePath $venvPy -ArgumentList $args -WorkingDirectory $repo `
    -PassThru -RedirectStandardOutput $stdout -RedirectStandardError $stderr
  $guardTripped = $false
  while (-not $proc.HasExited) {
    Start-Sleep -Seconds 2
    $proc.Refresh()
    $osNow = Get-CimInstance Win32_OperatingSystem
    $freeNow = $osNow.FreePhysicalMemory / 1MB
    $rssNow = $proc.WorkingSet64 / 1GB
    if ($rssNow -gt $cfg.max_worker_memory_gib -or $freeNow -lt $cfg.reserved_memory_gib) {
      $guardTripped = $true
      Stop-Process -Id $proc.Id -Force
      $proc.WaitForExit()
      break
    }
  }
  Get-Content $stdout -ErrorAction SilentlyContinue | Tee-Object -FilePath $log
  Get-Content $stderr -ErrorAction SilentlyContinue | Add-Content -LiteralPath $log
  if ($guardTripped) { throw 'Memory guard stopped the FEM process; inspect preserved logs and partial files.' }
  if ($proc.ExitCode -ne 0) { throw "FEM solve exited with $($proc.ExitCode); preserve log and partial outputs." }
} finally { Pop-Location }
& (Join-Path $PSScriptRoot 'verify_results.ps1')
