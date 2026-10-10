topic "Help";
[ $$0,0#00000000000000000000000000000000:Default]
[{_} 
[s0; [*+184 API Help]&]
[s0;+92 &]
[s0; [*+92 Introduction]&]
[s0;+92 &]
[s0; [+92 The BEMRosetta API provides a simple interface to many of 
the library’s functions from C, C`+`+ and Python. It allows 
applications and scripts to work with meshes, hydrodynamic models 
and simulation data without needing to understand BEMRosetta’s 
internal implementation.]&]
[s0;+92 &]
[s0; [+92 The C and C`+`+ interfaces use a small set of straightforward 
data structures, while the Python interface uses Python data 
structures. The three interfaces provide similar functionality, 
adapted to the conventions of each language.]&]
[s0;+92 &]
[s0; [*+92 Simple integration]&]
[s0;+92 &]
[s0; [+92 The API separates the application from BEMRosetta’s internal 
code and provides a more stable integration interface. Applications 
can use the exposed functions without depending directly on the 
library’s internal classes or implementation details.]&]
[s0;+92 &]
[s0; [+92 For C and C`+`+, there is no need to incorporate BEMRosetta’s 
source code into the application’s build system, whether it 
uses CMake or another compilation strategy. The library can be 
accessed in either of two ways:]&]
[s0;+92 &]
[s0;i150;O0; [+92 Linking with libbemrosetta.lib.]&]
[s0;i150;O0; [+92 Loading the dynamic library at runtime (DLL on Windows)]&]
[s0;+92 &]
[s0; [+92 Python provides access to similar operations through its 
own interface, using Python data structures rather than the C/C`+`+ 
structures.]&]
[s0;+92 &]
[s0; [*+92 Available functionality]&]
[s0;+92 &]
[s0; [+92 The API groups related operations into several functional 
areas:]&]
[s0;+92 &]
[s0;i150;O0; [*+92 Library management:][+92  initialise the library, retrieve 
its build date and time, enable or disable printed messages, 
and retrieve or clear the last error.]&]
[s0;i150;O0; [*+92 Mesh operations:][+92  load, save, select, duplicate 
and transform meshes; set body properties such as mass, centre 
of gravity, inertia, damping and mooring stiffness; retrieve 
geometric and hydrostatic properties; and extract or generate 
hull meshes, waterplane lids and control surfaces.]&]
[s0;i150;O0; [*+92 Hydrodynamic models: ][+92 create, load, save, select 
and duplicate BEM cases; define water depth, gravity, density, 
wave origin, frequencies and headings; associate meshes with 
bodies; and export cases for supported solvers.]&]
[s0;i150;O0; [*+92 Hydrodynamic data processing: ][+92 calculate Froude–Krylov 
forces, combine them with diffraction forces to obtain excitation 
forces, reset or transfer force and coefficient data between 
cases, and map potentials onto another mesh.]&]
[s0;i150;O0; [*+92 Solver settings: ][+92 configure selected Wamit and 
AQWA options and query solver capabilities.]&]
[s0;i150;O0; [*+92 OpenFAST file processing:][+92  read output files, 
access time histories and parameter information, calculate statistics, 
and read or modify parameters in .dat and .fst input files.]&]
[s0;+92 &]
[s0; [*+92 Working with the API]&]
[s0;+92 &]
[s0; [+92 Mesh, BEM and OpenFAST operations use an active mesh or model, 
selected through its identifier, rather than a pointer or a complex 
structure. Applications use these IDs to identify the objects 
they want to work with, while BEMRosetta manages their internal 
representation. This allows several meshes or cases to be loaded 
and the required one to be selected for subsequent operations.]&]
[s0;+92 &]
[s0; [+92 All features available in the GUI may be also available in 
the API, and vice versa. If you notice that a feature is missing 
from one that you have seen in the other, please let us know, 
since both use the same library and such features are usually 
very easy to add.]&]
[s0;+92 &]
[s0; [*+92 How to start in C`+`+]&]
[s0;+92 &]
[s0; [+92 Include ][C@(160.80.0)+92 `"libbemrosetta.hpp`"][+92  at the 
beginning.]&]
[s0;+92 &]
[s0; [+92 And create a variable with the path to the dynamic library 
to the constructor:]&]
[s0;+92 &]
[s0; [C+92 BEMRosetta mybemr(][C@(163.21.21)+92 `"PATH TO/libbemrosetta.dll`"][C+92 );]&]
[s0;+92 &]
[s0; [*+92 How to start in C]&]
[s0;+92 &]
[s0; [+92 Include ][C@(160.80.0)+92 `"libbemrosetta.h`"][+92  at the beginning.]&]
[s0;+92 &]
[s0; [+92 And create a variable with the path to the dynamic library 
to the constructor:]&]
[s0;+92 &]
[s0; [C+92 BEMRosetta mybemr `= BEMRosetta`_Init(][C@(163.21.21)+92 `"PATH`_TO/libbemrose
tta.dll`"][C+92 );]&]
[s0;+92 &]
[s0; [*+92 How to start in Python]&]
[s0;*+92 &]
[s0; [+92 Include ][C@(160.80.0)+92 from][C+92  libbemrosetta ][C@(160.80.0)+92 import][C+92  
BEMRosetta][+92  at the beginning.]&]
[s0;+92 &]
[s0; [+92 And create a variable with the path to the dynamic library 
to the constructor:]&]
[s0;+92 &]
[s0; [C+92 mybemr `= BEMRosetta(][C@(163.21.21)+92 `"PATH`_TO/libbemrosetta.dll`"][C+92 )]&]
[s0;+92 &]
[s0; [*+92 API function descriptions and examples]&]
[s0;+92 &]
[s0;#i150;b17;a17;O0; [^topic`:`/`/BEMRosetta`/main`/Python`_en`-us^+92 Python 
API][+92  and ][^topic`:`/`/BEMRosetta`/main`/Python`_Example`_en`-us^+92 Python 
Example]&]
[s0;#i150;b17;a17;O0; [^topic`:`/`/BEMRosetta`/main`/Cpp`_en`-us^+92 C`+`+ 
API][+92  and ][^topic`:`/`/BEMRosetta`/main`/Cpp`_Example`_en`-us^+92 C`+`+ 
Example][+92 .]&]
[s0;#i150;b17;a17;O0; [^topic`:`/`/BEMRosetta`/main`/C`_en`-us^+92 C 
API][+92  and ][^topic`:`/`/BEMRosetta`/main`/C`_Example`_en`-us^+92 C 
Example][+92 .]&]
[s0; ]]