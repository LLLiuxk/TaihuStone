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
%EXE_PATH% D:\VSprojects\TaihuStone\para\comp5\moonrock_128_60_default.txt
echo Task 1 Finished.
echo.


:: 运行 para2 (假设你已经创建了)
echo [Task 2] Running with para2.txt...
%EXE_PATH% D:\VSprojects\TaihuStone\para\comp5\moonrock_128_60_ang_highhigh.txt
echo Task 2 Finished.
echo.

:: 运行 para2 (假设你已经创建了)
echo [Task 2] Running with para2.txt...
%EXE_PATH% D:\VSprojects\TaihuStone\para\comp5\moonrock_128_60_ang_lowlow.txt
echo Task 2 Finished.
echo.


echo [Task 3] Running with para3.txt...
%EXE_PATH% D:\VSprojects\TaihuStone\para\comp5\moonrock_128_60_len_low.txt
echo Task 3 Finished.
echo.

echo [Task 3] Running with para4.txt...
%EXE_PATH% D:\VSprojects\TaihuStone\para\comp5\moonrock_128_60_len_high.txt
echo Task 3 Finished.
echo.

echo [Task 4] Running with para5.txt...
%EXE_PATH% D:\VSprojects\TaihuStone\para\comp5\moonrock_128_60_loc_low.txt
echo Task 4 Finished.
echo.

echo [Task 4] Running with para6.txt...
%EXE_PATH% D:\VSprojects\TaihuStone\para\comp5\moonrock_128_60_loc_high.txt
echo Task 4 Finished.
echo.

echo [Task 5] Running with para6.txt...
%EXE_PATH% D:\VSprojects\TaihuStone\para\comp5\moonrock_128_60_dir_high.txt
echo Task 5 Finished.
echo.

echo [Task 5] Running with para6.txt...
%EXE_PATH% D:\VSprojects\TaihuStone\para\comp5\moonrock_128_60_dir_low.txt
echo Task 5 Finished.
echo.

echo [Task 6] Running with para6.txt...
%EXE_PATH% D:\VSprojects\TaihuStone\para\comp5\moonrock_128_60+a-d.txt
echo Task 6 Finished.
echo.

echo [Task 6] Running with para6.txt...
%EXE_PATH% D:\VSprojects\TaihuStone\para\comp5\moonrock_128_60+d-a.txt
echo Task 6 Finished.
echo.

echo [Task 6] Running with para6.txt...
%EXE_PATH% D:\VSprojects\TaihuStone\para\comp5\moonrock_128_60+d-ll.txt
echo Task 6 Finished.
echo.

echo [Task 6] Running with para6.txt...
%EXE_PATH% D:\VSprojects\TaihuStone\para\comp5\moonrock_128_60+ll-a.txt
echo Task 6 Finished.
echo.

echo [Task 7] Running with para7.txt...
%EXE_PATH% D:\VSprojects\TaihuStone\para\comp5\moonrock_128_60_ave.txt
echo Task 5 Finished.
echo.

echo ==========================================
echo All Tasks Completed!
pause