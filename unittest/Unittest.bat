@title "Testing BEMRosetta for CLANGX64"
call UnittestBEMRosetta.bat CLANGX64 -r
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error testing and deploying BEMRosetta CLANGX64"
@echo "Testing BEMRosetta for MSVS22x64"

call IsMSVCInstalled
if "%MSVC_INSTALLED%"=="1" (
	@title "Testing BEMRosetta for MSVS22x64"
	call UnittestBEMRosetta.bat MSVS22x64 -r
	IF %ERRORLEVEL% NEQ 0 PAUSE "Error testing and deploying BEMRosetta MSVS22x64"	
)

@title "Testing BEMRosetta for Linux
@start bash UnittestBEMRosetta.sh
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error testing and deploying BEMRosetta Linux"