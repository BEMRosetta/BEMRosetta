@echo off
if exist clang_path.bat call clang_path.bat

echo Compiling test.cpp using libbemrosetta.lib. libbemrosetta.dll has to be with the .exe or in the PATH

clang++ -Wall test.cpp .test\\libbemrosetta.lib -o .test\\test.exe
if errorlevel 1 (
	echo Compilation failed
	exit /b 1 
)
echo Compilation successful

.test\\test.exe

del .test\\test.exe

echo Compiling test.cpp using directly libbemrosetta.dll

clang++ -Wall -DBEMROSETTA_DYNAMIC test.cpp -o .test\\test.exe
if errorlevel 1 (
	echo Compilation failed
	exit /b 1 
)
echo Compilation successful

.test\\test.exe

del .test\\test.exe

pause