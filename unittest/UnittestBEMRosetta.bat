mkdir .\.test
del .\.test\*.* /q

title BEMRosetta_cl
umk BEMRosetta BEMRosetta_cl %1 %2 +BEMR_CL    		.\.test\BEMRosetta_cl.exe
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
.\.test\BEMRosetta_cl.exe
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
title BEMRosetta
umk BEMRosetta BEMRosetta    %1 %2 +GUI 			.\.test\BEMRosetta.exe
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1

title BEMRosetta_cl DLL
umk BEMRosetta BEMRosetta_cl %1 %2 +BEMR_DLL,DLL	.\.test\libbemrosetta.dll
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
title BEMRosetta_cl TEST_DLL
umk BEMRosetta BEMRosetta_cl %1 %2 +BEMR_TEST_DLL	.\.test\testdll_bemrosetta.exe
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
.\.test\testdll_bemrosetta.exe  .\.test ..

copy .\.test\BEMRosetta.exe ..\bin
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
copy .\.test\BEMRosetta_cl.exe ..\bin
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
copy .\.test\libbemrosetta.dll ..\bin
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
copy .\.test\libbemrosetta.txt ..\bin
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
copy .\.test\libbemrosetta.py ..\bin
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1