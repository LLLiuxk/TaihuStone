@echo off
setlocal ENABLEEXTENSIONS

rem ====== 候选 Visual Studio 路径 ======
set VS1=C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\devenv.exe
set VS2=D:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\devenv.exe

set DEVENV=

if exist "%VS1%" (
    set DEVENV=%VS1%
) else if exist "%VS2%" (
    set DEVENV=%VS2%
)

if not defined DEVENV (
    echo [ERROR] Cannot find Visual Studio devenv.exe
    pause
    exit /b 1
)

start "Compare files" /B /MIN "%DEVENV%" /diff %2 %1 First:"%2" Second:"%1"
