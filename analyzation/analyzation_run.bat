@echo off
setlocal

set "ROOT=%~dp0"
set "PS1=%ROOT%merge_translucency_summaries.ps1"

if not exist "%PS1%" (
  echo [ERROR] Script not found: %PS1%
  pause
  exit /b 1
)

powershell -NoProfile -ExecutionPolicy Bypass -File "%PS1%"
set "CODE=%ERRORLEVEL%"

if not "%CODE%"=="0" (
  echo.
  echo [ERROR] merge script failed with code %CODE%
  pause
  exit /b %CODE%
)

echo.
echo [OK] Done. Check folder: %ROOT%analyzation
pause
exit /b 0

