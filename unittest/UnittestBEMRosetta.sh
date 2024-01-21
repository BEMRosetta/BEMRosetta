#!/bin/bash
mkdir ./.test
set -x
set -e
export PATH=$PATH:/home/aupa/upp

cwd=$(pwd)

umk BEMRosetta BEMRosetta_cl CLANG_17 -bd +BEMR_CL,SHARED ./.test/bemrosetta_cl
umk BEMRosetta BEMRosetta_cl GCC 	  -bd +BEMR_CL,SHARED ./.test/bemrosetta_cl

umk BEMRosetta BEMRosetta    CLANG_17 -bd +GUI,SHARED 	  ./.test/bemrosetta
umk BEMRosetta BEMRosetta    GCC 	  -bd +GUI,SHARED     ./.test/bemrosetta

echo -----------All done. NO error!
