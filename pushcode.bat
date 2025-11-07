@echo off
@title bat execute git auto commit
echo DOCS PUSH BAT

:: 获取当前日期（月日格式）
for /f "tokens=1,2 delims=/" %%a in ('date /t') do (
    set mm=%%a
    set dd=%%b
)

:: 组合commit消息
set commit_msg=%mm%%dd%
echo Commit消息: %commit_msg%

echo 1. Start submitting code to the local repository
git add .
echo=

echo 2. Commit the changes to the local repository
git commit -m "%commit_msg%"
echo=

echo 4. Push the changes to the remote git server
git push origin master
 
echo=
echo Batch execution complete!
pause