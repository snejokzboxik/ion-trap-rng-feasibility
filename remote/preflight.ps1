param(
  [switch]$DryRun,
  [switch]$TestMode,
  [ValidateSet('pass','missing-python','missing-dependency','low-memory','low-disk','missing-mesh','invalid-config','high-cpu','solver-running','unwritable-results')]
  [string]$TestScenario = 'pass',
  [string]$ConfigPath,
  [switch]$LockHeldByCaller
)
$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
if (-not $ConfigPath) { $ConfigPath = Join-Path $repo 'configs/r0_remote_fem.json' }

function Get-Sha256([string]$Path) {
  return (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant()
}

function Get-PreflightErrors($Snapshot, $Config) {
  $errors = New-Object System.Collections.Generic.List[string]
  if (-not $Snapshot.Windows) { $errors.Add('Host is not Windows.') }
  if (-not $Snapshot.PythonValid) { $errors.Add('Python 3.12 x64 virtual environment or required packages are invalid.') }
  if (-not $Snapshot.ConfigValid) { $errors.Add('Configuration JSON, checksum, schema, or immutable limits are invalid.') }
  if (-not $Snapshot.MeshValid) { $errors.Add('plateau014 mesh identity, SHA-256, or geometry statistics are invalid.') }
  if ($Snapshot.FreeMemoryGiB -lt $Config.minimum_free_memory_gib) { $errors.Add("Free RAM $($Snapshot.FreeMemoryGiB) GiB is below $($Config.minimum_free_memory_gib) GiB.") }
  if ($Snapshot.FreeDiskGiB -lt $Config.minimum_free_disk_gib) { $errors.Add("Free disk $($Snapshot.FreeDiskGiB) GiB is below $($Config.minimum_free_disk_gib) GiB.") }
  if ($Snapshot.CpuUtilizationPct -gt $Config.maximum_preflight_cpu_utilization_percent) { $errors.Add("CPU utilization $($Snapshot.CpuUtilizationPct)% exceeds the configured workload limit.") }
  if ($Snapshot.ConflictingSolver) { $errors.Add('Another m3a3_octant_p2 solve is already running.') }
  if (-not $Snapshot.ResultsWritable) { $errors.Add('Results directory is not writable.') }
  if (-not $Snapshot.RunLockAvailable) { $errors.Add('Another ion-trap FEM run holds the single-run lock.') }
  return ,$errors
}

function Test-ConfigPolicy($Config) {
  if ($Config.stage -notmatch '^R0 preparation') { return $false }
  if ($Config.reference_stage -ne 'M3a3 symmetry-reduced octant P2 electrostatics') { return $false }
  if ($Config.mesh_tag -ne 'plateau014') { return $false }
  if ($Config.existing_mesh_required -ne 'results/m3a3_octant_mesh_plateau014.npz') { return $false }
  if ($Config.expected_mesh_vertices -ne 9335 -or $Config.expected_mesh_tetrahedra -ne 39615 -or $Config.expected_p2_dofs -ne 63759) { return $false }
  if ($Config.python -ne '3.12' -or $Config.threads -ne 1) { return $false }
  if ($Config.minimum_free_memory_gib -lt 13 -or $Config.minimum_free_disk_gib -lt 10) { return $false }
  if ($Config.reserved_memory_gib -lt 5 -or $Config.max_worker_memory_gib -gt 8) { return $false }
  if ($Config.estimated_peak_memory_gib -le 0 -or $Config.max_worker_memory_gib -le 0 -or $Config.max_worker_memory_gib -gt $Config.estimated_peak_memory_gib) { return $false }
  if ($Config.minimum_free_memory_gib -lt ($Config.reserved_memory_gib + $Config.estimated_peak_memory_gib)) { return $false }
  if ($Config.threads -ne 1 -or $Config.scope -notmatch 'single sequential P2 solve') { return $false }
  if ($Config.maximum_preflight_cpu_utilization_percent -le 0 -or $Config.maximum_preflight_cpu_utilization_percent -gt 100) { return $false }
  if ($Config.outputs.Count -ne 2) { return $false }
  if ($Config.outputs -notcontains 'results/m3a3_octant_fields_plateau014.npz') { return $false }
  if ($Config.outputs -notcontains 'results/m3a3_octant_solve_boundary_plateau014.json') { return $false }
  if ($Config.mesh_sha256 -notmatch '^[0-9a-f]{64}$' -or $Config.mesh_audit_sha256 -notmatch '^[0-9a-f]{64}$') { return $false }
  return $true
}

if ($TestMode) {
  $cfg = Get-Content -LiteralPath $ConfigPath -Raw | ConvertFrom-Json
  $s = [ordered]@{Windows=$true;PythonValid=$true;ConfigValid=(Test-ConfigPolicy $cfg);MeshValid=$true;FreeMemoryGiB=16.0;FreeDiskGiB=20.0;CpuUtilizationPct=5.0;ConflictingSolver=$false;ResultsWritable=$true;RunLockAvailable=$true}
  switch ($TestScenario) {
    'missing-python' { $s.PythonValid=$false }
    'missing-dependency' { $s.PythonValid=$false }
    'low-memory' { $s.FreeMemoryGiB=12.99 }
    'low-disk' { $s.FreeDiskGiB=9.99 }
    'missing-mesh' { $s.MeshValid=$false }
    'invalid-config' { $s.ConfigValid=$false }
    'high-cpu' { $s.CpuUtilizationPct=90.0 }
    'solver-running' { $s.ConflictingSolver=$true }
    'unwritable-results' { $s.ResultsWritable=$false }
  }
  $testErrors = Get-PreflightErrors $s $cfg
  if ($testErrors.Count -eq 0) { Write-Output 'TEST PREFLIGHT PASS; no host inspection or FEM execution.'; exit 0 }
  Write-Output ('TEST PREFLIGHT FAIL: ' + ($testErrors -join ' ')); exit 2
}

try {
  if ([Environment]::OSVersion.Platform -ne [PlatformID]::Win32NT) { throw 'Remote FEM setup requires Windows.' }
  $cfg = Get-Content -LiteralPath $ConfigPath -Raw | ConvertFrom-Json
  $signaturePath = [IO.Path]::ChangeExtension($ConfigPath, '.sha256')
  $configValid = (Test-Path -LiteralPath $signaturePath) -and (Test-ConfigPolicy $cfg)
  if ($configValid) {
    $expectedConfigHash = ((Get-Content -LiteralPath $signaturePath -Raw).Trim() -split '\s+')[0].ToLowerInvariant()
    $configValid = ($expectedConfigHash -match '^[0-9a-f]{64}$') -and ((Get-Sha256 $ConfigPath) -eq $expectedConfigHash)
  }

  $venvPy = Join-Path $repo '.venv/Scripts/python.exe'
  $pythonValid = Test-Path -LiteralPath $venvPy
  $pythonProbe = $null
  if ($pythonValid) {
    $probeCode = @'
import importlib.metadata as md, json, struct, sys
names={'numpy':'numpy','scipy':'scipy','PyYAML':'yaml','scikit-fem':'skfem','gmsh':'gmsh','matplotlib':'matplotlib','pytest':'pytest'}
mods={}; errors=[]
for dist,mod in names.items():
 try:
  __import__(mod); mods[dist]=md.version(dist)
 except Exception as exc: errors.append(f'{dist}: {exc}')
print(json.dumps({'version':list(sys.version_info[:3]),'bits':struct.calcsize('P')*8,'modules':mods,'errors':errors}))
'@
    $probeRaw = & $venvPy -c $probeCode 2>&1
    if ($LASTEXITCODE -eq 0) {
      try {
        $pythonProbe = ($probeRaw -join "`n") | ConvertFrom-Json
        $pythonValid = ($pythonProbe.version[0] -eq 3 -and $pythonProbe.version[1] -eq 12 -and $pythonProbe.bits -eq 64 -and $pythonProbe.errors.Count -eq 0)
        $lockLines = Get-Content -LiteralPath (Join-Path $repo 'requirements-lock.txt') | Where-Object { $_ -match '^([A-Za-z0-9_.-]+)==([^#\s]+)$' }
        foreach ($line in $lockLines) {
          if ($line -match '^([A-Za-z0-9_.-]+)==([^#\s]+)$') {
            $distName=$Matches[1];$want=$Matches[2];$got=$pythonProbe.modules.$distName
            if ($got -ne $want) { $pythonValid=$false }
          }
        }
      } catch { $pythonValid = $false }
    } else { $pythonValid = $false }
  }

  $meshPath = Join-Path $repo $cfg.existing_mesh_required
  $meshAuditPath = Join-Path $repo $cfg.mesh_audit_required
  $meshValid = (Test-Path -LiteralPath $meshPath) -and (Test-Path -LiteralPath $meshAuditPath)
  if ($meshValid) {
    $meshValid = ((Get-Sha256 $meshPath) -eq $cfg.mesh_sha256) -and ((Get-Sha256 $meshAuditPath) -eq $cfg.mesh_audit_sha256)
    try {
      $meshAudit = Get-Content -LiteralPath $meshAuditPath -Raw | ConvertFrom-Json
      $meshProbeCode = @'
import json,numpy as np,sys
p=sys.argv[1];d=np.load(p,allow_pickle=False);v=d['points'];t=d['tetrahedra']
edges=np.concatenate([t[:,[0,1]],t[:,[0,2]],t[:,[0,3]],t[:,[1,2]],t[:,[1,3]],t[:,[2,3]]]);nedge=len(np.unique(np.sort(edges,axis=1),axis=0))
print(json.dumps({'point_shape':list(v.shape),'tet_shape':list(t.shape),'finite':bool(np.isfinite(v).all()),'unique_edges':nedge,'p2_dofs':int(len(v)+nedge)}))
'@
      $meshRaw=& $venvPy -c $meshProbeCode $meshPath 2>&1
      if ($LASTEXITCODE -ne 0) { $meshValid=$false }
      else {
        $meshProbe=($meshRaw -join "`n") | ConvertFrom-Json
        $meshValid = $meshValid -and ($meshAudit.mesh_targets_r0.hcenter -eq 0.014) -and ($meshAudit.vertices -eq $cfg.expected_mesh_vertices) -and ($meshAudit.tetrahedra -eq $cfg.expected_mesh_tetrahedra) -and ($meshAudit.estimated_P2_dofs -eq $cfg.expected_p2_dofs) -and ($meshProbe.point_shape -join ',') -eq '9335,3' -and ($meshProbe.tet_shape -join ',') -eq '39615,4' -and $meshProbe.finite -and ($meshProbe.p2_dofs -eq $cfg.expected_p2_dofs)
      }
    } catch { $meshValid=$false }
  }

  $os = Get-CimInstance Win32_OperatingSystem
  $freeGiB = [math]::Round(($os.FreePhysicalMemory * 1KB) / 1GB, 3)
  $driveName = ([IO.Path]::GetPathRoot($repo).Substring(0,1))
  $disk = Get-PSDrive -Name $driveName
  $freeDiskGiB = [math]::Round($disk.Free / 1GB, 3)
  $cpuRow = Get-CimInstance Win32_PerfFormattedData_PerfOS_Processor -Filter "Name='_Total'"
  if (-not $cpuRow) { throw 'Could not read system CPU utilization counter; refusing to launch.' }
  $cpuPct = [double]$cpuRow.PercentProcessorTime
  $processes = @(Get-Process | Select-Object Id,ProcessName,CPU,WorkingSet64 | Sort-Object WorkingSet64 -Descending | Select-Object -First 8)
  $wmiProcesses = @(Get-CimInstance Win32_Process | Where-Object { $_.CommandLine -match 'm3a3_octant_p2\.py.*--solve' })
  $conflictingSolver = ($wmiProcesses.Count -gt 0)
  $resultsPath = Join-Path $repo 'results'
  $resultsWritable = Test-Path -LiteralPath $resultsPath
  $probePath = Join-Path $resultsPath ('.preflight-' + [guid]::NewGuid().ToString('N') + '.tmp')
  if ($resultsWritable) {
    try {
      $stream=[IO.File]::Open($probePath,[IO.FileMode]::CreateNew,[IO.FileAccess]::ReadWrite,[IO.FileShare]::None)
      $stream.WriteByte(0);$stream.Flush();$stream.Dispose()
      $probeBytes=[IO.File]::ReadAllBytes($probePath)
      if ($probeBytes.Length -ne 1) { $resultsWritable=$false }
    }
    catch { $resultsWritable=$false }
    finally { if (Test-Path -LiteralPath $probePath) { Remove-Item -LiteralPath $probePath -Force } }
  }
  $runLock = New-Object System.Threading.Mutex($false, 'Global\IonTrapRngM3a3P2Solve')
  $runLockAvailable=$true
  if (-not $LockHeldByCaller) {
    try { $runLockAvailable=$runLock.WaitOne(0) }
    catch [System.Threading.AbandonedMutexException] { $runLockAvailable=$true }
    if ($runLockAvailable) { $runLock.ReleaseMutex() }
  }
  $snapshot=[ordered]@{Windows=$true;PythonValid=$pythonValid;ConfigValid=$configValid;MeshValid=$meshValid;FreeMemoryGiB=$freeGiB;FreeDiskGiB=$freeDiskGiB;CpuUtilizationPct=$cpuPct;ConflictingSolver=$conflictingSolver;ResultsWritable=$resultsWritable;RunLockAvailable=$runLockAvailable}
  $errors=Get-PreflightErrors $snapshot $cfg
  $reportedMeshSha = $null
  if ($meshValid) { $reportedMeshSha = $cfg.mesh_sha256 }
  $report=[ordered]@{passed=($errors.Count -eq 0);dryRun=[bool]$DryRun;repository=$repo;config_sha256=(Get-Sha256 $ConfigPath);mesh_sha256=$reportedMeshSha;python=$pythonProbe;free_memory_gib=$freeGiB;free_disk_gib=$freeDiskGiB;cpu_utilization_percent=$cpuPct;cpu_limit_percent=$cfg.maximum_preflight_cpu_utilization_percent;reserved_memory_gib=$cfg.reserved_memory_gib;worker_memory_limit_gib=$cfg.max_worker_memory_gib;threads=$cfg.threads;process_memory_snapshot=@($processes | ForEach-Object {[ordered]@{pid=$_.Id;name=$_.ProcessName;working_set_gib=[math]::Round($_.WorkingSet64/1GB,2);cpu_seconds=$_.CPU}});errors=@($errors)}
  $report | ConvertTo-Json -Depth 8
  if ($errors.Count -gt 0) { exit 2 }
  if ($DryRun) { Write-Output 'Preflight passed in dry-run mode; no FEM solver was started.' }
  else { Write-Output 'Preflight passed; resource and input gates are currently satisfied.' }
} catch {
  Write-Output ('Preflight failed closed: ' + $_.Exception.Message)
  exit 2
}
