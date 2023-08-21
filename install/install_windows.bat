..\..\..\umk examples-bazaar BEMRosetta/src/BEMRosetta    CLANGx64 -r +GUI ..\unittest\.test\BEMRosetta.exe
if %errorlevel% neq 0 ( 
 	pause "Error compiling BEMRosetta.exe"
 	exit
)
..\..\..\umk examples-bazaar BEMRosetta/src/BEMRosetta_cl CLANGx64 -r      ..\unittest\.test\BEMRosetta_cl.exe
if %errorlevel% neq 0 ( 
 	pause "Error compiling BEMRosetta_cl.exe"
 	exit
)