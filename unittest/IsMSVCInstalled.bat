@echo off
setlocal
set "MSVC_INSTALLED=0"

set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
if not exist "%VSWHERE%" (
    goto :end
)

for /f "delims=" %%I in ('
    "%VSWHERE%" -latest ^
               -products * ^
               -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 ^
               -property installationPath ^
               -format value
') do set "VSPATH=%%I"

if not defined VSPATH goto :end

set "CLPATH=%VSPATH%\VC\Tools\MSVC"
if not exist "%CLPATH%" goto :end

set "MSVC_INSTALLED=1"

:end
endlocal & set "MSVC_INSTALLED=%MSVC_INSTALLED%"
exit /b
