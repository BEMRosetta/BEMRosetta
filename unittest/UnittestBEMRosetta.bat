@mkdir .\.test
@del .\.test\*.* /q

@title Compiling BEMRosetta_cl %1
@echo Compiling BEMRosetta_cl %1
umk BEMRosetta BEMRosetta_cl %1 %2 +BEMR_CL  		.\.test\BEMRosetta_cl.exe  		
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
@.\.test\BEMRosetta_cl -paramfile TestBEMRosetta_CL-time.txt
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1

@title Compiling BEMRosetta GUI %1
@echo Compiling BEMRosetta GUI %1
umk BEMRosetta BEMRosetta    %1 %2 +GUI 			.\.test\BEMRosetta.exe
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1

@title Compiling BEMRosetta_cl DLL %1
@echo Compiling BEMRosetta_cl DLL %1
umk BEMRosetta BEMRosetta_cl %1 %2 +BEMR_DLL,DLL	.\.test\libbemrosetta.dll
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1

rem @.\.test\BEMRosetta_cl -paramfile TestBEMRosetta_CL.txt
rem @IF %ERRORLEVEL% NEQ 0 EXIT /B 1

@title Compiling BEMRosetta_cl TEST_DLL %1
@echo Compiling BEMRosetta_cl TEST_DLL %1
umk BEMRosetta BEMRosetta_cl %1 %2 +BEMR_TEST_DLL	.\.test\testdll_bemrosetta.exe
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
.\.test\testdll_bemrosetta.exe  .\.test ..
 
@title Copying BEMRosetta %1
@echo Copying BEMRosetta %1
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