param(
  [switch]$DryRun,
  [switch]$Execute,
  [switch]$TestMode,
  [ValidateSet('pass','failed-solver','incomplete-output')]
  [string]$TestScenario = 'pass'
)
$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path

if ($TestMode) {
  if ($TestScenario -eq 'pass') { Write-Output 'TEST RUNNER PASS; simulated dry-run did not start a solver.'; exit 0 }
  if ($TestScenario -eq 'failed-solver') { Write-Output 'TEST RUNNER FAIL: simulated solver exit code was nonzero; artifacts preserved.'; exit 3 }
  Write-Output 'TEST RUNNER FAIL: simulated output set was incomplete; verification rejected it.'; exit 4
}
if ($DryRun -and $Execute) { Write-Output 'Choose either -DryRun or -Execute, not both.'; exit 2 }
$doExecute = [bool]$Execute
$cfgPath = Join-Path $repo 'configs/r0_remote_fem.json'
$preflight = Join-Path $PSScriptRoot 'preflight.ps1'
$verify = Join-Path $PSScriptRoot 'verify_results.ps1'
$python = Join-Path $repo '.venv/Scripts/python.exe'
$solver = Join-Path $repo 'scripts/m3a3_octant_p2.py'
$meshTag = 'plateau014'
$fieldOutput = Join-Path $repo "results/m3a3_octant_fields_$meshTag.npz"
$auditOutput = Join-Path $repo "results/m3a3_octant_solve_boundary_$meshTag.json"
$runId = (Get-Date -Format 'yyyyMMdd_HHmmss') + '_' + [guid]::NewGuid().ToString('N').Substring(0,8)
$log = Join-Path $repo "results/r0_fem_$runId.log"
$stdout = "$log.stdout"
$stderr = "$log.stderr"
$mutex = $null
$mutexHeld = $false
$proc = $null
$workerGuardTripped = $false
$oldEnv = @{}

try {
  $mutex = New-Object System.Threading.Mutex($false, 'Global\IonTrapRngM3a3P2Solve')
  try { $mutexHeld = $mutex.WaitOne(0) } catch [System.Threading.AbandonedMutexException] { $mutexHeld = $true }
  if (-not $mutexHeld) { throw 'Another bounded ion-trap FEM run holds the single-run lock.' }

  if (-not (Test-Path -LiteralPath $preflight)) { throw 'Required preflight.ps1 is missing; refusing to launch.' }
  $psExe = (Get-Process -Id $PID).Path
  $preflightArgs = @('-NoProfile','-ExecutionPolicy','Bypass','-File',"`"$preflight`"",'-LockHeldByCaller')
  $pfOut = "$log.preflight.stdout"; $pfErr = "$log.preflight.stderr"
  $pf = Start-Process -FilePath $psExe -ArgumentList $preflightArgs -WorkingDirectory $repo -PassThru -Wait -NoNewWindow -RedirectStandardOutput $pfOut -RedirectStandardError $pfErr
  Get-Content $pfOut -ErrorAction SilentlyContinue | Add-Content -LiteralPath $log
  Get-Content $pfErr -ErrorAction SilentlyContinue | Add-Content -LiteralPath $log
  if ($pf.ExitCode -ne 0) { throw "Preflight failed with exit code $($pf.ExitCode); no FEM solver was launched. See $log" }

  if (-not $doExecute) {
    Write-Output "DRY RUN PASSED: validated host, dependencies, configuration, mesh, resources, workload, output permissions and single-run lock."
    Write-Output "Would run: `"$python`" scripts/m3a3_octant_p2.py --tag plateau014 --solve"
    Write-Output 'No FEM solver launched.'
    exit 0
  }

  if (-not (Test-Path -LiteralPath $python) -or -not (Test-Path -LiteralPath $solver)) { throw 'Python environment or solver script missing; no solve launched.' }
  if ((Test-Path -LiteralPath $fieldOutput) -or (Test-Path -LiteralPath $auditOutput)) {
    throw 'A prior plateau014 result exists. Refusing to overwrite; archive/inspect it explicitly before a new run.'
  }
  $cfg = Get-Content -LiteralPath $cfgPath -Raw | ConvertFrom-Json
  $configSha = (Get-FileHash -LiteralPath $cfgPath -Algorithm SHA256).Hash.ToLowerInvariant()
  $envNames = @('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','IONTRAP_RUN_CONFIG_SHA256')
  foreach ($n in $envNames) { $oldEnv[$n] = [Environment]::GetEnvironmentVariable($n, 'Process') }
  $threadCount = [string]$cfg.threads
  $env:OMP_NUM_THREADS = $threadCount; $env:OPENBLAS_NUM_THREADS = $threadCount
  $env:MKL_NUM_THREADS = $threadCount; $env:NUMEXPR_NUM_THREADS = $threadCount
  $env:IONTRAP_RUN_CONFIG_SHA256 = $configSha
  $proc = Start-Process -FilePath $python -ArgumentList @(('"' + $solver + '"'),'--tag',$meshTag,'--solve') -WorkingDirectory $repo -PassThru -RedirectStandardOutput $stdout -RedirectStandardError $stderr
  while (-not $proc.HasExited) {
    Start-Sleep -Seconds 2
    $proc.Refresh()
    if ($proc.HasExited) { break }
    try {
      $osNow = Get-CimInstance Win32_OperatingSystem
      $freeNow = ($osNow.FreePhysicalMemory * 1KB) / 1GB
      $rssNow = $proc.WorkingSet64 / 1GB
      $unsafe = ($freeNow -lt $cfg.reserved_memory_gib) -or ($rssNow -gt $cfg.max_worker_memory_gib)
    } catch { $unsafe = $true }
    if ($unsafe) {
      $workerGuardTripped = $true
      if ($proc -and -not $proc.HasExited) {
        # Kill through the Process handle returned by Start-Process: this can
        # target only this child, never another user's workload or PID search.
        $proc.Kill()
        $proc.WaitForExit()
      }
      break
    }
  }

  Get-Content $stdout -ErrorAction SilentlyContinue | Add-Content -LiteralPath $log
  Get-Content $stderr -ErrorAction SilentlyContinue | Add-Content -LiteralPath $log
  if ($workerGuardTripped) { throw 'Memory guard stopped only the verified FEM worker; logs and partial files were preserved.' }
  if ($proc.ExitCode -ne 0) { throw "FEM solver failed with exit code $($proc.ExitCode); preserve logs and any partial output." }
  & $verify
  if ($LASTEXITCODE -ne 0) { throw 'Solver returned success but output-integrity verification failed; preserve files for inspection.' }
  Write-Output "P2 solve and output-integrity verification completed. Log: $log"
  exit 0
} catch {
  if ($proc -and -not $proc.HasExited) {
    try { $proc.Kill(); $proc.WaitForExit() }
    catch { Write-Output 'Unable to stop the child FEM worker through its process handle; no other process was touched.' }
  }
  Write-Output ("FEM run refused or failed: " + $_.Exception.Message)
  if (Test-Path -LiteralPath $stdout) { Get-Content $stdout -ErrorAction SilentlyContinue | Add-Content -LiteralPath $log }
  if (Test-Path -LiteralPath $stderr) { Get-Content $stderr -ErrorAction SilentlyContinue | Add-Content -LiteralPath $log }
  Write-Output "Logs and any partial outputs preserved. Log: $log"
  exit 3
} finally {
  foreach ($n in $oldEnv.Keys) { [Environment]::SetEnvironmentVariable($n, $oldEnv[$n], 'Process') }
  if ($mutexHeld -and $mutex) { try { $mutex.ReleaseMutex() } catch {} }
  if ($mutex) { $mutex.Dispose() }
}
