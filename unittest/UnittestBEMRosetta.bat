@echo off
mkdir .\test

umk BEMRosetta BEMRosetta    CLANGX64 -r +GUI test\BEMRosettaCLANG.exe
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
umk BEMRosetta BEMRosetta_cl CLANGX64 -r      test\BEMRosetta_clCLANG.exe
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
test\BEMRosetta_clCLANG.exe
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1

umk BEMRosetta BEMRosetta    MSBT19x64 -r +GUI test\BEMRosetta.exe
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
umk BEMRosetta BEMRosetta_cl MSBT19x64 -r      test\BEMRosetta_cl.exe
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
test\BEMRosetta_cl.exe
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1

umk BEMRosetta BEMRosetta_cl MSBT19x64 -r +DLL		test\libbemrosetta.dll
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
umk BEMRosetta BEMRosetta_cl MSBT19x64 -r +TEST_DLL	test\testdll_bemrosetta.exe
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
test\testdll_bemrosetta.exe  test ..

copy test\BEMRosetta.exe ..\bin
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
copy test\BEMRosetta_cl.exe ..\bin
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
copy test\libbemrosetta.dll ..\bin
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
copy test\libbemrosetta.txt ..\bin
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
