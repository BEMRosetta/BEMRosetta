@echo off
if exist clang_path.bat call clang_path.bat
echo Compiling test.c

clang -Wall test.c .test\\libbemrosetta.lib -o test.exe
if errorlevel 1 (
	echo Compilation failed
	exit /b 1
)
echo Compilation successful

test.exe
pause