#!/bin/bash
echo -----------Compiling CLI files
set -x 
set -e
rm -f Makefile		
cp ./Makefile_cli ./Makefile
make
cp ./BEMRosetta/BEMRosetta_cl.out ../BEMRosetta_cl