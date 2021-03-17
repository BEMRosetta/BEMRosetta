call UnittestBEMRosetta.bat CLANGX64
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error testing and deploying BEMRosetta CLANGX64"
call UnittestBEMRosetta.bat MSBT19x64
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error testing and deploying BEMRosetta MSBT19x64"