@echo off
pushd "%~dp0"
echo ===================================================
echo STL Mesh Simplification Tool (Limit: 8MB)
echo ===================================================

:: Check if python is in PATH
where python >nul 2>nul
if %errorlevel% neq 0 (
    echo Error: python is not found in PATH!
    echo Please install Python and add it to your system PATH.
    pause
    exit /b 1
)

:: Check if the Python script exists
if not exist "simplify_stls.py" (
    echo Error: simplify_stls.py not found in the same directory!
    pause
    exit /b 1
)

echo Scanning and processing STL files...
python simplify_stls.py "%~dp0"

echo.
echo ===================================================
echo Processing Finished!
echo ===================================================
pause
popd
