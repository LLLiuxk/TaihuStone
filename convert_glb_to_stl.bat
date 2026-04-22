@echo off
setlocal
chcp 65001 >nul

set "SCRIPT=%~dp0glb_to_stl.py"

where python >nul 2>nul
if errorlevel 1 (
    echo Error: Python was not found in PATH.
    echo Please install Python or add it to PATH, then try again.
    pause
    exit /b 1
)

if not exist "%SCRIPT%" (
    echo Error: missing converter script:
    echo "%SCRIPT%"
    pause
    exit /b 1
)

if "%~1"=="" (
    echo No file was dropped. Converting all .glb files in:
    echo %~dp0
    python "%SCRIPT%" "%~dp0"
) else (
    python "%SCRIPT%" %*
)

echo.
pause
