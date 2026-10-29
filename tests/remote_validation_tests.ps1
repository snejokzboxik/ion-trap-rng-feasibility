param()
$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$preflight = Join-Path $repo 'remote/preflight.ps1'
$runner = Join-Path $repo 'remote/run_fem.ps1'
$python = Join-Path $repo '.venv/Scripts/python.exe'
if (-not (Test-Path $python) -and $env:IONTRAP_TEST_PYTHON) { $python = $env:IONTRAP_TEST_PYTHON }
$psExe = (Get-Process -Id $PID).Path
if (-not (Test-Path $python)) { throw 'Run remote/setup_windows.ps1 before these tests.' }

function Invoke-ChildScript([string]$Path, [string[]]$Arguments) {
  $quoted = @('-NoProfile','-ExecutionPolicy','Bypass','-File',('"' + $Path + '"')) + $Arguments
  $child = Start-Process -FilePath $psExe -ArgumentList $quoted -PassThru -Wait -NoNewWindow
  return $child.ExitCode
}

$checks = @(
  @{name='normal dry-run';path=$preflight;args=@('-TestMode','-DryRun','-TestScenario','pass');expected=0},
  @{name='missing Python';path=$preflight;args=@('-TestMode','-TestScenario','missing-python');expected=2},
  @{name='missing dependency';path=$preflight;args=@('-TestMode','-TestScenario','missing-dependency');expected=2},
  @{name='insufficient RAM';path=$preflight;args=@('-TestMode','-TestScenario','low-memory');expected=2},
  @{name='insufficient disk';path=$preflight;args=@('-TestMode','-TestScenario','low-disk');expected=2},
  @{name='missing mesh';path=$preflight;args=@('-TestMode','-TestScenario','missing-mesh');expected=2},
  @{name='invalid configuration';path=$preflight;args=@('-TestMode','-TestScenario','invalid-config');expected=2},
  @{name='busy CPU';path=$preflight;args=@('-TestMode','-TestScenario','high-cpu');expected=2},
  @{name='existing solver';path=$preflight;args=@('-TestMode','-TestScenario','solver-running');expected=2},
  @{name='unwritable output directory';path=$preflight;args=@('-TestMode','-TestScenario','unwritable-results');expected=2},
  @{name='runner dry-run';path=$runner;args=@('-TestMode','-TestScenario','pass');expected=0},
  @{name='failed solver';path=$runner;args=@('-TestMode','-TestScenario','failed-solver');expected=3},
  @{name='incomplete output';path=$runner;args=@('-TestMode','-TestScenario','incomplete-output');expected=4}
)
$failures = New-Object System.Collections.Generic.List[string]
foreach ($check in $checks) {
  $actual = Invoke-ChildScript $check.path $check.args
  if ($actual -ne $check.expected) { $failures.Add("$($check.name): expected $($check.expected), got $actual") }
  else { Write-Output "PASS: $($check.name) (exit $actual)" }
}
& $python -m pytest -q --basetemp (Join-Path $repo 'results/pytest_remote_validation') (Join-Path $repo 'tests/test_remote_output_validation.py')
if ($LASTEXITCODE -ne 0) { $failures.Add('Python output-integrity fixture tests failed.') }
if ($failures.Count) { $failures | ForEach-Object { Write-Error $_ }; exit 1 }
Write-Output 'All PowerShell workflow and result-integrity validation tests passed. No FEM solver was started.'
