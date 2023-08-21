cd ..\install
@echo "Making BEMRosetta installer"
title "Making BEMRosetta installer"
call make_installer.bat
@IF %ERRORLEVEL% NEQ 0 PAUSE "Error making the installer"