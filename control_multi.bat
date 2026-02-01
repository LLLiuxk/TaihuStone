@echo off
set EXE_PATH=D:\VSprojects\TaihuStone\x64\Release\TaihuStone.exe

:: 检查 exe 是否存在
if not exist %EXE_PATH% (
    echo Error: %EXE_PATH% not found!
    pause
    exit /b
)

echo ==========================================
echo Starting Batch Processing...
echo ==========================================

:: 运行 para1
echo [Task 1] Running with para file...
%EXE_PATH% D:\VSprojects\TaihuStone\para\RockSetEr_rotated60_0.5_0.88.txt
echo Task 1 Finished.
echo.

:: 运行 para2 (假设你已经创建了)
echo [Task 2] Running with para2.txt...
%EXE_PATH% D:\VSprojects\TaihuStone\para\RockSetAr_rotated60_0.5_0.88.txt
echo Task 2 Finished.
echo.

echo [Task 3] Running with para3.txt...
%EXE_PATH% D:\VSprojects\TaihuStone\para\RockSetCr_rotated60_0.5_0.88.txt
echo Task 3 Finished.
echo.

echo [Task 4] Running with para3.txt...
%EXE_PATH% D:\VSprojects\TaihuStone\para\RockSetDr_rotated60_0.5_0.88.txt
echo Task 4 Finished.
echo.

echo ==========================================
echo All Tasks Completed!
pause