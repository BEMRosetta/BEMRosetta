call UnittestBEMRosetta.bat CLANGX64
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error testing and deploying BEMRosetta CLANGX64"
call UnittestBEMRosetta.bat MSBT19x64
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error testing and deploying BEMRosetta MSBT19x64"
cd ..\install
call make_installer.bat
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error making the installer"
