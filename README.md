<!--
  Title: BEMRosetta
  Description: Hydrodynamic coefficients viewer and converter for Boundary Element Method solver formats.
  Authors: Iñaki Zabala.
  -->

# BEMRosetta
**Hydrodynamic coefficients viewer and converter for Boundary Element Method solver formats.**

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/platforms-windows_linux-blue.svg" alt="Platforms">
<img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/build-passed-success.svg" alt="Build status">
<img src="https://img.shields.io/github/last-commit/izabala123/bemrosetta.svg" alt="Last commit">

<p align="center">
  <img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/deepcwind.jpg" width="45%" title="DeepCWind mesh in Windows">
  <img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/deepcwind_linux.JPG" width="45%" title="DeepCWind mesh in Linux">
</p>

[Boundary Element Methods](https://en.wikipedia.org/wiki/Boundary_element_method) are extensively used to model hydrodynamic forces in offshore devices like ships, offshore wind platforms and wave energy converters. These solvers use device geometry mesh to get some hydrodynamics coefficients as radiation damping, added mass, wave diffraction force, and wave excitation force. All these data is saved in file formats incompatible between them. These may avoid to use the coefficients between programs. 

BEMRosetta allows to load the hydrodynamic coefficients from a format saving it in another. In addition it allows to compare the results obtained between programs, the results between similar geometries and the same geometry with different discretization levels.

Moreover, BEMRosetta allows to view and visually compare the meshes from different programs.

BEMRosetta runs on Windows and Linux, **no install is necessary in Windows** [(see Install)](https://github.com/izabala123/BEMRosetta/tree/master/install), and it includes a GUI, [a command line version](https://github.com/izabala123/BEMRosetta/blob/master/other/test), a library (DLL), and glue code for Python. 

## Features
### - Supported file formats

* BEM coefficients
  * Load-View
    * [Wamit](https://www.wamit.com/): .out, .3sc, 3fk, .1, .3, .4, .hst, .7, .8, .9, .12s, .12d
    * [HAMS](https://github.com/YingyiLiu/HAMS): ControlFile.in
    * [Nemoh](https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-presentation-192863.kjsp) and [Capytaine](https://github.com/mancellin/capytaine): Nemoh.cal, Mesh/Hydrostatics*.dat, Mesh/KH*.dat, RadiationCoefficients.tec, ExcitationForce.tec, DiffractionForce.tec, FKForce.tec, IRF.tec
    * [Diodore](https://www.principia-group.com/blog/product/diodore/): .hdb
    * [Hydrostar](https://marine-offshore.bureauveritas.com/hydrostar-software-powerful-hydrodynamic): .out
    * [OrcaWave](https://www.orcina.com/orcaflex/): .yml
    * [OpenFAST-Wamit](https://www.nrel.gov/wind/nwtc/openfast.html): HydroDyn.dat
    * [SeaFEM/TDyn-Nemoh](http://www.compassis.com/compass/en/Productos/SeaFEM): .flavia.inf, RadiationCoefficients.tec, ExcitationForce.tec, DiffractionForce.tec, FKForce.tec
    * [Ansys AQWA](https://www.ansys.com/products/structures/ansys-aqwa): .lis, .ah1, .qtf
    * [FOAMM](http://www.eeng.nuim.ie/coer/downloads/): .mat
    
  * Save
    * [Wamit](https://www.wamit.com/): .out, .1, .3, .hst, .4, .7, .8, .9, .12s, .12d
    * [HAMS](https://github.com/YingyiLiu/HAMS): ControlFile.in and all the folder structure.
	* [Diodore](https://www.principia-group.com/blog/product/diodore/): .hdb
    * [Ansys AQWA](https://www.ansys.com/products/structures/ansys-aqwa): .qtf
    * [OpenFAST-Wamit](https://nwtc.nrel.gov/FAST): HydroDyn.dat

* Case files
    * Load-View
      * [HAMS](https://github.com/YingyiLiu/HAMS): ControlFile.in
      * [Ansys AQWA](https://www.ansys.com/products/structures/ansys-aqwa): .dat
      * [Nemoh](https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-presentation-192863.kjsp) and [Capytaine](https://github.com/mancellin/capytaine): Nemoh.cal and all the folder structure.
    * Save
      * [HAMS](https://github.com/YingyiLiu/HAMS): ControlFile.in and all the folder structure.
      * [Nemoh](https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-presentation-192863.kjsp) and [Capytaine](https://github.com/mancellin/capytaine): Nemoh.cal and all the folder structure.
      
* Mesh files
  * Load-View
    * [Wamit](https://www.wamit.com/): .gdf, pan.dat
    * [HAMS](https://github.com/YingyiLiu/HAMS): .pnl
    * [Nemoh](https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-mesh-192932.kjsp?RH=1489593406974) and [Capytaine](https://github.com/mancellin/capytaine): .dat
    * [Ansys AQWA](https://www.ansys.com/products/structures/ansys-aqwa): .dat
    * [Hydrostar](https://marine-offshore.bureauveritas.com/hydrostar-software-powerful-hydrodynamic): .hst
    * [Salome](https://www.salome-platform.org/): .dat
    * [STL format](https://en.wikipedia.org/wiki/STL_(file_format)): .stl (binary and text)
    * [SeaFEM/TDyn](http://www.compassis.com/compass/en/Productos/SeaFEM): .msh
  * Save
    * [Wamit](https://www.wamit.com/): .gdf
    * [HAMS](https://github.com/YingyiLiu/HAMS): HullMesh.pnl, WaterplaneMesh.pnl
    * [Ansys AQWA](https://www.ansys.com/products/structures/ansys-aqwa): .dat
    * [STL format](https://en.wikipedia.org/wiki/STL_(file_format)): .stl (binary and text)

* Time domain simulations
  * Load-View
    * [OpenFAST](https://www.nrel.gov/wind/nwtc/openfast.html): .out, .outb
    * [Deeplines Wind](https://www.principia-group.com/blog/product/deeplines-wind/): .db
    * [Ansys AQWA Naut](https://www.ansys.com/products/structures/ansys-aqwa): .lis
    * CSV: .csv
  * Save
    * [OpenFAST](https://www.nrel.gov/wind/nwtc/openfast.html): .out
    * CSV: .csv


### - Load the hydrodynamic coefficients from one format and save them in another

The goal is to have a good robustness in the handling of files


### - Compare the hydrodynamic coefficients for the same geometry from different software

- Damping for the same geometry got from different solvers
  
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/2%20solvers%20B.jpg" width="300" title="Damping for the same geometry got from different solvers"></p>

- Excitation force for the same geometry got from different solvers_
  
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/2%20solvers%20exc.jpg" width="300" title="Excitation force for the same geometry got from different solvers"></p>

### - Forces handling

It simmetrizes the available forces in all directions, averaging them when they are available on both possitive and negative headings. Some examples cases:
* Only the forces on positive headings from 0 to 180º have been processed: Symmetrize duplicates them to the negative heading values from 0 to -180º
* Both positive and negative headings forces have been processed: Symmetrize averages them

### - Compare the hydrodynamic coefficients for the same geometry for different discretization levels
### - Compare the hydrodynamic coefficients for different geometries

- Damping for different offshore wind floating platforms_
  
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/offshore%20wind%20platforms%20B.png" width="300" title="Damping for different offshore wind floating platforms"></p>

- Excitation force for different offshore wind floating platforms_
  
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/offshore%20wind%20platforms%20exc.jpg" width="300" title="Excitation force for different offshore wind floating platforms"></p>

### - FOAMM connection

[Finite Order Approximation by Moment-Matching (FOAMM)](http://www.eeng.nuim.ie/coer/wp-content/uploads/2019/02/FOAMM-Manual.pdf) is an application developed by N. Faedo, Y. Peña-Sanchez and J. V. Ringwood in the [Maynooth University](https://www.maynoothuniversity.ie/)'s [Centre for Ocean Energy Research (COER)](http://www.eeng.nuim.ie/coer/), that implements the moment-matching based frequency-domain identification algorithm.

BEMRosetta allows an interactive and seamless FOAMM connection to get state space coefficients.

### - Mesh loading, combining them for visual comparison 

Several meshes can be loaded in this basic viewer, allowing a visual comparison of geometries.

<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/deepcwind.jpg" width="300" title="Mesh loading"></p>


### - Mesh handling

- Interactive mesh rotation and translation around user defined center
- Automatic free surface, underwater surface, center of buoyancy, hydrostatic stiffness matrix, and other parameters calculation
- Improved viewer including dropdown menu in viewer screen
- Hydrostatic stiffness matrix viewer
- Mesh healing option
    
### - Case launcher, Nemoh & HAMS

Added Nemoh and [HAMS](https://github.com/YingyiLiu/HAMS) launcher. It can loadexisting files from HAMS, Nemoh or ANSYS AQWA, it lets you editing it, and creates the set of files to launch Nemoh and HAMS from a .bat file (it replaces the classic Nemoh MATLAB launcher)

### - Time domain simulations

BEMRosetta includes a time domain simulations viewer supporting OpenFAST, Deeplines Wind, Ansys AQWA Naut and csv formats, designed to be very easy to use.
Files may be opened by drag and drop, and parameters are filtered by name or units.

<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/FAST_Reader.png" width="800" title="OpenFAST .out/.outb reader"></p>

### - OrcaFlex command line

If you have an [OrcaFlex](https://www.orcina.com/orcaflex/) licence, the command line version allows you to perform operations not available directly in OrcaWave/OrcaFlex, like:
- Calculating hydrodynamic coefficients with OrcaWave.
- Performing time domain simulations with OrcaWave.
- Save the results of time domain simulations to .csv files.

### - Other

All files, mesh, case or BEM files, can be loaded by Drag and Drop or Copy and Paste from file explorer in Windows and Linux.

<p align="center">
  <img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/Underwater.png" width="45%" title="Underwater mesh and waterline">
  <img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/Mesh.png" width="45%" title="All mesh ans waterline">
</p>

## Acknowledgments

J. C. C. Portillo, J. C. C. Henriques, M. J. Sanchez-Lara, J. Galvan, A. Otter, M. Alonso, A. Aristondo.<br/>
Some file parsing strategies taken from the [BEMIO project](https://wec-sim.github.io/bemio/).<br/>
Done with the [U++ multiplatform library](https://www.ultimatepp.org/).

## License

Copyright © 2019-2021 Iñaki Zabala, Markel Peñalba, Yerai Peña-Sanchez, Thomas Kelly.

BEMRosetta is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.\
BEMRosetta is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for details. You should have received a copy of the GNU General Public License along with BEMRosetta. If not, see http://www.gnu.org/licenses/.<br/>
<br/>
<br/>
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/logos/BEMRosetta.png" width="200" title="BEMRosetta"></p>

