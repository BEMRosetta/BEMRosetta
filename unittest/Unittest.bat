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


setlocal enabledelayedexpansion

set OUTFILE=SHA256SUMS.txt

if exist "%OUTFILE%" del "%OUTFILE%"

for %%F in (".\.test\*.exe" ".\.test\*.dll") do (
    echo Processing %%F
    set HASH=
    for /f "tokens=1" %%H in ('
        certutil -hashfile "%%F" SHA256 ^| findstr /r /v /c:"hash" /c:"CertUtil"
    ') do (
        if not defined HASH set HASH=%%H
    )
    echo !HASH!  %%~nxF>>".\.test\%OUTFILE%"
)