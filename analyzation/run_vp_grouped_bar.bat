@echo off
cd /d "%~dp0"
echo Running draw_vp_weight_grouped_bar.py with taihu_env...
"%~dp0..\limitstl\taihu_env\Scripts\python.exe" draw_vp_weight_grouped_bar.py
echo.
echo Process finished.
pause
