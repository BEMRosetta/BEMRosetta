#!/bin/bash
mkdir ./.test
set -x
set -e
export PATH=$PATH:/home/aupa/upp

cwd=$(pwd)

echo "Compiling BEMRosetta_cl"
umk BEMRosetta BEMRosetta_cl CLANG_17 -bd +BEMR_CL,SHARED ./.test/bemrosetta_cl
umk BEMRosetta BEMRosetta_cl GCC 	  -r  +BEMR_CL,SHARED ./.test/bemrosetta_cl

echo "Testing BEMRosetta_cl"
./.test/bemrosetta_cl -paramfile TestBEMRosetta_CL-mesh.txt
./.test/bemrosetta_cl -paramfile TestBEMRosetta_CL-bem.txt
./.test/bemrosetta_cl -paramfile TestBEMRosetta_CL-time.txt

echo "Compiling BEMRosetta"
umk BEMRosetta BEMRosetta    CLANG_17 -bd +GUI,SHARED 	  ./.test/bemrosetta
umk BEMRosetta BEMRosetta    GCC 	  -bd +GUI,SHARED     ./.test/bemrosetta

echo -----------All done. NO error!
