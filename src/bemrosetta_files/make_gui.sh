#!/bin/bash
echo -----------Compiling GUI files
set -x 
set -e
rm -f Makefile		
cp ./Makefile_gui ./Makefile
make
cp ./BEMRosetta/BEMRosetta.out ../BEMRosetta