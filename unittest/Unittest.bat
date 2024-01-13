@title "Testing BEMRosetta for CLANGX64"
call UnittestBEMRosetta.bat CLANGX64 -r
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error testing and deploying BEMRosetta CLANGX64"
@echo "Testing BEMRosetta for MSVS22x64"
@echo call UnittestBEMRosetta.bat MSVS22x64 -r
@echo IF %ERRORLEVEL% NEQ 0 PAUSE "Error testing and deploying BEMRosetta MSVS22x64"
@title "Testing BEMRosetta for Linux
@start call bash UnittestBEMRosetta.sh
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error testing and deploying BEMRosetta Linux"