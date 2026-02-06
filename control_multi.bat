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
%EXE_PATH% D:\VSprojects\TaihuStone\para\high_rock_rotated_128_50_0.6_3_25_1.7_0.7_0.89_0.35_11011.txt
echo Task 1 Finished.
echo.

:: 运行 para2 (假设你已经创建了)
echo [Task 2] Running with para2.txt...
%EXE_PATH% D:\VSprojects\TaihuStone\para\high_rock_rotated_128_60_0.6_3_25_1.7_0.7_0.89_0.35_11011.txt
echo Task 2 Finished.
echo.

echo [Task 3] Running with para3.txt...
%EXE_PATH% D:\VSprojects\TaihuStone\para\high_rock_rotated_128_70_0.6_3_25_1.7_0.7_0.89_0.35_11011.txt
echo Task 3 Finished.
echo.

echo [Task 4] Running with para3.txt...
%EXE_PATH% D:\VSprojects\TaihuStone\para\high_rock_rotated_128_80_0.6_3_25_1.7_0.7_0.89_0.35_11011.txt
echo Task 4 Finished.
echo.

echo [Task 5] Running with para3.txt...
%EXE_PATH% D:\VSprojects\TaihuStone\para\high_rock_rotated_128_90_0.6_3_25_1.7_0.7_0.89_0.35_11011.txt
echo Task 5 Finished.
echo.

echo [Task 6] Running with para3.txt...
%EXE_PATH% D:\VSprojects\TaihuStone\para\high_rock_rotated_128_60_0.6_3_25_1.7_0.7_0.86_0.35_01011.txt
echo Task 6 Finished.
echo.

echo ==========================================
echo All Tasks Completed!
pause