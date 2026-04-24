@echo off
cd /d "%~dp0"
echo Running plot_translucency_final_nama.py with taihu_env...
"%~dp0..\limitstl\taihu_env\Scripts\python.exe" plot_translucency_final_nama.py --input translucency_final_nama.xlsx --output translucency_final_nama_lines.png
echo.
echo Process finished.
pause
