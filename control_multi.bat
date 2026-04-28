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
::%EXE_PATH% D:\VSprojects\TaihuStone\para\comp\namaqualand_128_60_loc0.txt
::%EXE_PATH% D:\VSprojects\TaihuStone\para\asteroid_smooth_128_50_70_0.7.txt
echo Task 1 Finished.
echo.

:: 运行 para2 (假设你已经创建了)
echo [Task 2] Running with para2.txt...
::%EXE_PATH% D:\VSprojects\TaihuStone\para\comp\namaqualand_128_60_len0.txt
::%EXE_PATH% D:\VSprojects\TaihuStone\para\comp\namaqualand_128_60_default_thres_0.25.txt
echo Task 2 Finished.
echo.

echo [Task 3] Running with para3.txt...
::%EXE_PATH% D:\VSprojects\TaihuStone\para\comp\namaqualand_128_60_dir0.txt
%EXE_PATH% D:\VSprojects\TaihuStone\para\asteroid_smooth_128_50_70_0.8.txt
echo Task 3 Finished.
echo.

echo [Task 4] Running with para4.txt...
::%EXE_PATH% D:\VSprojects\TaihuStone\para\comp\namaqualand_128_60_ang0.txt
%EXE_PATH% D:\VSprojects\TaihuStone\para\asteroid_smooth_128_50_70_0.9.txt
echo Task 4 Finished.
echo.


echo ==========================================
echo All Tasks Completed!
pause