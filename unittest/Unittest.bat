@echo "Testing BEMRosetta for CLANGX64"
call UnittestBEMRosetta.bat CLANGX64 -r
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error testing and deploying BEMRosetta CLANGX64"
@echo "Testing BEMRosetta for MSBT19x64"
call UnittestBEMRosetta.bat MSBT19x64 -r
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error testing and deploying BEMRosetta MSBT19x64"
cd ..\install
@echo "Making BEMRosetta installer"
call make_installer.bat
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error making the installer"
