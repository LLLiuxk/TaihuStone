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
%EXE_PATH% D:\VSprojects\TaihuStone\para\comp\namaqualand_128_60_0.6_3_25_1.7_0.7_0.89_0.35_11011.txt
:: namaqualand60_0.5_0.88.txt
:: namaqualand_128_60_0.6_4_30_1.7_0.7_0.89_0.3_01011.txt
:: %EXE_PATH% D:\VSprojects\TaihuStone\para\Block_Island80_0.5_0.9.txt
echo Task 1 Finished.
echo.

:: 运行 para2 (假设你已经创建了)
:: echo [Task 2] Running with para2.txt...
:: %EXE_PATH% para2.txt
:: echo Task 2 Finished.
:: echo.


echo ==========================================
echo All Tasks Completed!
pause