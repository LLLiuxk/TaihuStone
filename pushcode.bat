@echo off
@title bat execute git auto commit
echo DOCS PUSH BAT

for /f "tokens=2-4 delims=/.- " %%a in ('date /t') do (
    set month=%%a
    set day=%%b
)
set commit_msg=%month%%day%

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