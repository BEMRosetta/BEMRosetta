#!/bin/bash
mkdir .\test
set -x
set -e

cwd=$(pwd)

umk BEMRosetta BEMRosetta    CLANG 	-bd +GUI,SHARED 	~/bemrosetta
umk BEMRosetta BEMRosetta    GCC 	-bd +GUI,SHARED 	~/bemrosetta

umk BEMRosetta BEMRosetta_cl CLANG 	-bd +BEMR_CL,SHARED ~/bemrosetta_cl
~/bemrosetta_cl

umk BEMRosetta BEMRosetta_cl GCC 	-bd +BEMR_CL,SHARED ~/bemrosetta_cl
~/bemrosetta_cl
echo -----------All done. NO error!
