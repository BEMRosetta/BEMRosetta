@if not exist .\.test mkdir .\.test
@del .\.test\*.* /q

@title Compiling BEMRosetta_cl %1
@echo Compiling BEMRosetta_cl %1
umk BEMRosetta BEMRosetta_cl %1 %2 +BEMR_CL  -r	.\.test\BEMRosetta_cl.exe  		
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error compiling BEMRosetta"
@title Testing BEMRosetta_cl %1
@echo Testing BEMRosetta_cl %1
@.\.test\BEMRosetta_cl -paramfile TestBEMRosetta_CL-mesh.txt
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error testing BEMRosetta"
@.\.test\BEMRosetta_cl -paramfile TestBEMRosetta_CL-bem.txt
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error testing BEMRosetta"
@.\.test\BEMRosetta_cl -paramfile TestBEMRosetta_CL-time.txt
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error testing BEMRosetta"
@.\.test\BEMRosetta_cl -paramfile TestBEMRosetta_CL-wind.txt
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error testing BEMRosetta"
@.\.test\BEMRosetta_cl -paramfile TestBEMRosetta_CL-orca.txt
rem @IF %ERRORLEVEL% NEQ 0 PAUSE "Error testing BEMRosetta"

@title Compiling BEMRosetta_cl DLL %1
@echo Compiling BEMRosetta_cl DLL %1
umk BEMRosetta BEMRosetta_cl %1 %2 +BEMR_DLL,DLL -r	.\.test\libbemrosetta.dll
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error compiling and testing BEMRosetta"

@title Compiling BEMRosetta_cl TEST_DLL %1
@echo Compiling BEMRosetta_cl TEST_DLL %1
umk BEMRosetta BEMRosetta_cl %1 %2 +BEMR_TEST_DLL -r	.\.test\testdll_bemrosetta.exe
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error compiling BEMRosetta"
.\.test\testdll_bemrosetta.exe  .
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error testing BEMRosetta"
del /Q /F .\.test\testdll_bemrosetta.exe

python test.py 
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error testing BEMRosetta"

@title Compiling BEMRosetta %1
@echo Compiling BEMRosetta %1
umk BEMRosetta BEMRosetta %1 %2 +GUI  -r	.\.test\BEMRosetta.exe  		
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error compiling BEMRosetta"

@title Copying BEMRosetta %1
@echo Copying BEMRosetta %1
copy .\.test\BEMRosetta.exe ..\_bin
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error copying BEMRosetta"
copy .\.test\BEMRosetta_cl.exe ..\_bin
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error copying BEMRosetta"
copy .\.test\libbemrosetta.dll ..\_bin
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error copying BEMRosetta"
copy .\.test\libbemrosetta.txt ..\_bin
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error copying BEMRosetta"
copy .\.test\libbemrosetta.py ..\_bin
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error copying BEMRosetta"

del /Q /F .\.test\TurbSim2.bts
del /Q /F .\.test\hello.*
del /Q /F .\.test\*.log
del /Q /F .\.test\*.txt
del /Q /F .\.test\*.csv
del /Q /F .\.test\__pycache__\*.*
rd  .\.test\__pycache__

