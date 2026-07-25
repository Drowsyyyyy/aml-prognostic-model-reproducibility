$ErrorActionPreference = "Stop"
$ProjectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$Python = Join-Path $ProjectRoot ".venv\Scripts\python.exe"

& $Python (Join-Path $ProjectRoot "scripts\01_prepare_data.py")
& $Python (Join-Path $ProjectRoot "scripts\02_validate_models.py")
& $Python (Join-Path $ProjectRoot "scripts\03_summarize_results.py")
& $Python (Join-Path $ProjectRoot "scripts\99_quality_checks.py")

