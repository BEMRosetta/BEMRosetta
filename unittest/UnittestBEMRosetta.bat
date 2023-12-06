@mkdir .\.test
@del .\.test\*.* /q

@title Compiling BEMRosetta_cl %1
@echo Compiling BEMRosetta_cl %1
umk BEMRosetta BEMRosetta_cl %1 %2 +BEMR_CL  		.\.test\BEMRosetta_cl.exe  		
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
@title Testing BEMRosetta_cl %1
@echo Testing BEMRosetta_cl %1
@.\.test\BEMRosetta_cl -paramfile TestBEMRosetta_CL-mesh.txt
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
@.\.test\BEMRosetta_cl -paramfile TestBEMRosetta_CL-bem.txt
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

@title Compiling BEMRosetta_cl TEST_DLL %1
@echo Compiling BEMRosetta_cl TEST_DLL %1
umk BEMRosetta BEMRosetta_cl %1 %2 +BEMR_TEST_DLL	.\.test\testdll_bemrosetta.exe
@IF %ERRORLEVEL% NEQ 0 EXIT /B 1
.\.test\testdll_bemrosetta.exe  .\.test ..
 
rem @title Copying BEMRosetta %1
rem @echo Copying BEMRosetta %1
rem copy .\.test\BEMRosetta.exe ..\_bin
rem @IF %ERRORLEVEL% NEQ 0 EXIT /B 1
rem copy .\.test\BEMRosetta_cl.exe ..\_bin
rem @IF %ERRORLEVEL% NEQ 0 EXIT /B 1
rem copy .\.test\libbemrosetta.dll ..\_bin
rem @IF %ERRORLEVEL% NEQ 0 EXIT /B 1
rem copy .\.test\libbemrosetta.txt ..\_bin
rem @IF %ERRORLEVEL% NEQ 0 EXIT /B 1
rem copy .\.test\libbemrosetta.py ..\_bin
rem @IF %ERRORLEVEL% NEQ 0 EXIT /B 1