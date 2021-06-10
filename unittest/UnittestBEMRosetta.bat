mkdir .\.test
del .\.test\*.* /q

umk BEMRosetta BEMRosetta    %1 -r +GUI .\.test\BEMRosetta.exe
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
umk BEMRosetta BEMRosetta_cl %1 -r      .\.test\BEMRosetta_cl.exe
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
.\.test\BEMRosetta_cl.exe
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1

umk BEMRosetta BEMRosetta_cl %1 -r +DLL		.\.test\libbemrosetta.dll
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
umk BEMRosetta BEMRosetta_cl %1 -r +TEST_DLL	.\.test\testdll_bemrosetta.exe
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
.\.test\testdll_bemrosetta.exe  .\.test ..

copy .\.test\BEMRosetta.exe ..\bin
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
copy .\.test\BEMRosetta.exe ..\other\test\BEMRosetta_experimental.exe
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
copy .\.test\BEMRosetta_cl.exe ..\bin
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
copy .\.test\libbemrosetta.dll ..\bin
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
copy .\.test\libbemrosetta.txt ..\bin
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
copy .\.test\libbemrosetta.py ..\bin
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1