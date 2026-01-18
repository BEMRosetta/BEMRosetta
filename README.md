<!--
  Title: BEMRosetta
  Description: Hydrodynamic coefficients viewer and converter for Boundary Element Method solver formats.
  Authors: Iñaki Zabala.
  -->

# BEMRosetta
**Hydrodynamic solvers viewer and converter.**

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/platforms-windows_linux-blue.svg" alt="Platforms">
<img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/build-passed-success.svg" alt="Build status">
<img src="https://img.shields.io/github/last-commit/izabala123/bemrosetta.svg" alt="Last commit">

<p align="center">
  <img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/DeepCWind2.png" width="100%" title="DeepCWind mesh in Windows and Linux">
</p>

[Boundary Element Methods](https://en.wikipedia.org/wiki/Boundary_element_method) are extensively used to model hydrodynamic forces in offshore devices like ships, offshore wind platforms and wave energy converters. These solvers use device geometry mesh to get some hydrodynamics coefficients as radiation damping, added mass, wave diffraction force, and wave excitation force. All these data is saved in file formats incompatible between them. These may avoid to use the coefficients between programs. 

BEMRosetta allows to load the hydrodynamic coefficients from a format saving it in another. In addition it allows to compare the results obtained between programs, the results between similar geometries and the same geometry with different discretization levels.

Moreover, BEMRosetta allows to view and visually compare the meshes from different programs.

BEMRosetta runs on Windows and Linux, **no install is necessary in Windows** [(see Install)](https://github.com/izabala123/BEMRosetta/tree/master/install), and it includes a GUI, [a command line version](https://github.com/izabala123/BEMRosetta/blob/master/other/test), a library (DLL), and glue code for Python. 

## Features
### - Supported file formats

* BEM coefficients
  * Load-View
    * [Wamit](https://www.wamit.com/): .out, .3sc, 3fk, .1, .3, .4, .hst, .5p, .7, .8, .9, .12s, .12d
    * [Ansys AQWA](https://www.ansys.com/products/structures/ansys-aqwa): .lis, .ah1, .qtf
	* [OrcaFlex](https://www.orcina.com/orcaflex/): .yml
    * [OrcaWave](https://www.orcina.com/orcaflex/): .owr (requires OrcaWave installed)
    * [HAMS](https://github.com/YingyiLiu/HAMS): ControlFile.in
    * [Nemoh](https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-presentation-192863.kjsp) and [Capytaine](https://github.com/mancellin/capytaine): Nemoh.cal, Mesh/Hydrostatics*.dat, Mesh/KH*.dat, RadiationCoefficients.tec, ExcitationForce.tec, DiffractionForce.tec, FKForce.tec, IRF.tec
    * [Capytaine](https://github.com/mancellin/capytaine): .nc
    * [Bemio](https://wec-sim.github.io/bemio/api.html#writing-the-data-to-the-standard-bemio-format): .h5
	* [Matlab]: .mat
    * [Diodore](https://www.principia-group.com/blog/product/diodore/): .hdb
    * [Hydrostar](https://marine-offshore.bureauveritas.com/hydrostar-software-powerful-hydrodynamic): .out
    * [OpenFAST-Wamit](https://www.nrel.gov/wind/nwtc/openfast.html): HydroDyn.dat
    * [SeaFEM/TDyn-Nemoh](http://www.compassis.com/compass/en/Productos/SeaFEM): .flavia.inf, RadiationCoefficients.tec, ExcitationForce.tec, DiffractionForce.tec, FKForce.tec
    * [FOAMM](http://www.eeng.nuim.ie/coer/downloads/): .mat
    
  * Save
    * [Wamit](https://www.wamit.com/): .out, .1, .3, .hst, .4, .7, .8, .9, .12s, .12d
    * [HAMS](https://github.com/YingyiLiu/HAMS): ControlFile.in and all the folder structure.
    * [Bemio](https://wec-sim.github.io/bemio/api.html#writing-the-data-to-the-standard-bemio-format): .h5
    * [Matlab]: .mat
	* [Diodore](https://www.principia-group.com/blog/product/diodore/): .hdb
    * [Ansys AQWA](https://www.ansys.com/products/structures/ansys-aqwa): .qtf
    * [OpenFAST-Wamit](https://nwtc.nrel.gov/FAST): HydroDyn.dat

* Case files
    * Load-View
      * [HAMS](https://github.com/YingyiLiu/HAMS) and [HAMS-MREL](https://research.tudelft.nl/en/datasets/hydrodynamic-analysis-of-marine-structures-marine-renewable-energ/): ControlFile.in
      * [Nemoh](https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-presentation-192863.kjsp) and [Capytaine](https://github.com/mancellin/capytaine): Nemoh.cal and all the folder structure.
      * [Ansys AQWA](https://www.ansys.com/products/structures/ansys-aqwa): .dat
	  * [Wamit](https://www.wamit.com/): .pot, .frc, .cfg
	  * [OrcaWave](https://www.orcina.com/orcaflex/): .yml
    * Save
      * [HAMS](https://github.com/YingyiLiu/HAMS) and [HAMS-MREL](https://research.tudelft.nl/en/datasets/hydrodynamic-analysis-of-marine-structures-marine-renewable-energ/): ControlFile.in and all the folder structure.
      * [Nemoh](https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-presentation-192863.kjsp) and [Capytaine](https://github.com/mancellin/capytaine): Nemoh.cal and all the folder structure.
	  * [Ansys AQWA](https://www.ansys.com/products/structures/ansys-aqwa): .dat
	  * [Wamit](https://www.wamit.com/): .pot, .frc, .cfg 
	  * [OrcaWave](https://www.orcina.com/orcaflex/): .yml
	  * [Hydrostar](https://marine-offshore.bureauveritas.com/hydrostar-software-powerful-hydrodynamic): .hsg, .mcn, ,rao, ,qtf, .rdf, .dft
      
* Mesh files
  * Load-View
    * [Wamit](https://www.wamit.com/): .gdf, .idf, pan.dat
    * [HAMS](https://github.com/YingyiLiu/HAMS): .pnl
    * [Nemoh](https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-mesh-192932.kjsp?RH=1489593406974) and [Capytaine](https://github.com/mancellin/capytaine): .dat
    * [Ansys AQWA](https://www.ansys.com/products/structures/ansys-aqwa): .dat
    * [Hydrostar](https://marine-offshore.bureauveritas.com/hydrostar-software-powerful-hydrodynamic): .hst
    * [Salome](https://www.salome-platform.org/): .dat
    * [STL format](https://en.wikipedia.org/wiki/STL_(file_format)): .stl (binary and text)
    * [SeaFEM/TDyn](http://www.compassis.com/compass/en/Productos/SeaFEM): .msh
    * [GMSH](https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format): .msh
    * [OrcaFlex](https://www.orcina.com/orcaflex/): .yml
    * [OrcaWave](https://www.orcina.com/orcaflex/): .owr (requires OrcaWave installed)
	* [MIKE21](https://www.mikepoweredbydhi.com/products/mike-21-3): [.grd](https://www.mikepoweredbydhi.com/products/-/media/B1DE5EDAD83C4D51B36B380882771D2A.ashx)
	* [GeomView off](http://www.geomview.org/docs/html/OFF.html#OFF): .off
  * Save
    * [Wamit](https://www.wamit.com/): .gdf
    * [HAMS](https://github.com/YingyiLiu/HAMS): HullMesh.pnl, WaterplaneMesh.pnl
    * [Ansys AQWA](https://www.ansys.com/products/structures/ansys-aqwa): .dat
    * [Hydrostar](https://marine-offshore.bureauveritas.com/hydrostar-software-powerful-hydrodynamic): .hst
    * [STL format](https://en.wikipedia.org/wiki/STL_(file_format)): .stl (binary and text)
	* [MIKE21](https://www.mikepoweredbydhi.com/products/mike-21-3): [.grd](https://www.mikepoweredbydhi.com/products/-/media/B1DE5EDAD83C4D51B36B380882771D2A.ashx)
	* [GeomView off](http://www.geomview.org/docs/html/OFF.html#OFF): .off
	
* Time domain simulations
  * Load-View
    * [OpenFAST](https://www.nrel.gov/wind/nwtc/openfast.html): .out, .outb
    * [Deeplines Wind](https://www.principia-group.com/blog/product/deeplines-wind/): .db
    * [Ansys AQWA Naut](https://www.ansys.com/products/structures/ansys-aqwa): .lis
    * CSV: .csv
  * Save
    * [OpenFAST](https://www.nrel.gov/wind/nwtc/openfast.html): .out
    * CSV: .csv

* Turbulent wind
  * Load
    * [OpenFAST](https://www.nrel.gov/wind/nwtc/openfast.html): .bts, .sim
    * [DNV Bladed](https://www.dnv.com/services/wind-turbine-design-software-bladed-3775): .wnd

  * Save
    * [OpenFAST](https://www.nrel.gov/wind/nwtc/openfast.html): .bts

### - Load the hydrodynamic coefficients from one format and save them in another

The goal is to have a good robustness in the handling of files


### - Compare the hydrodynamic coefficients for the same geometry from different software

- Damping for the same geometry got from different solvers
  
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/2%20solvers%20B.jpg" width="300" title="Damping for the same geometry got from different solvers"></p>

- Excitation force for the same geometry got from different solvers_
  
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/2%20solvers%20exc.jpg" width="300" title="Excitation force for the same geometry got from different solvers"></p>

### - Forces handling

It symmetrizes the available forces in all directions, averaging them when they are available on both possitive and negative headings. Some examples cases:
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

## Citations

**Processing of potentials and hydrodynamic coefficients at panel level, and remapping according to a defined device structure.**<br/>
[Dynamic Modelling of HiveWind Floating Wind Substructure in OpenFAST](https://www.osti.gov/biblio/2482263). R. Bergua, I. Zabala, A. Gomez, L. Wang, O. Pena, J. Jonkman, J. Peña (2024). NAWEA/WindTech 2024, New Brunswick. USA..

**A method for translating the hydrodynamic coefficients to another reference frame. A convergence acceleration method applied to the hydrodynamic coefficients.**<br/>
[Post-processing techniques to improve the results of hydrodynamic boundary element method solvers](https://www.sciencedirect.com/science/article/pii/S0029801824002506). I. Zabala, J.C.C. Henriques, T.E. Kelly, P.P. Ricci, J.M. Blanco. Ocean Engineering 295, 116913.

**Post-processing for irregular frequency removal.**<br/>
[A post-processing technique for removing ‘irregular frequencies’ and other issues in the results from BEM solvers](https://mural.maynoothuniversity.ie/id/eprint/16316/). T. Kelly, I. Zabala, Y. Peña-Sanchez, J. Ringwood, J. Henriques, J.M. Blanco (2022). International Journal of Marine Energy 5 (1), 123-131

**Review of the main features.**<br/>
[BEMRosetta: An open-source hydrodynamic coefficients converter and viewer integrated with Nemoh and FOAMM](https://mural.maynoothuniversity.ie/id/eprint/16257/). I. Zabala, Y. Peña-Sanchez, T.E. Kelly, J.C.C. Henriques, M. Peñalba, N. Faedo, J. Ringwood, J.M. Blanco (2021). 14th European Wave and Tidal Energy Conference, 5-9th Sept 2021, Plymouth, UK.

## Short Videos

[A series of short videos](https://www.youtube.com/@BEMRosetta) has been created to help you use the software's features in a simple way. The following table shows their characteristics.

| Videos | Introduction | Vessel Mesh | BEM Solvers | Hydrodynamic Coefficients | Mooring |
| --- |:---:|:---:|:---:|:---:|:---:|
| [General Presentation](https://www.youtube.com/watch?v=sEjo7yI6rl0) | X | X | X | X | |
| [Windows Install](https://www.youtube.com/watch?v=VMt0QTr-3SI) | X | | | |
| [GNU/Linux Install](https://www.youtube.com/watch?v=VPGCHyuXHa4) | X | | | | |
| [Command line Presentation](https://www.youtube.com/watch?v=vY7_t6-_f-o) | X | X | | X | |
| [Animation](https://www.youtube.com/watch?v=lEvHlYaOB6k) | X | X | | | |
| [Mooring Editor](https://www.youtube.com/watch?v=qy9UWf7wj_U) | X | | | | X |
| [Mooring Animation](https://www.youtube.com/watch?v=-5bqJii5ZG4) | X | X | | | X |
| [From mesh to case in less than 4 mins](https://www.youtube.com/watch?v=ACKvyERhNYI) | | X | X | X | |
| [Revolution Mesh](https://www.youtube.com/watch?v=45blMIfCnYc) | | X | | | |
| [Get Lid and Hull](https://www.youtube.com/watch?v=mOueU1Hnh4U) | | X | | | |
| [Mesh Extrusion](https://youtu.be/T-To-PC9lFo) | | X | | | |
| [AQWA QTF and OpenFAST](https://www.youtube.com/watch?v=jCG9eQxWajY) | | | | X | |

If you are missing any video, please request it.
	
## Acknowledgments

J. C. Portillo, J. C. C. Henriques, J. M. Blanco, M. J. Sanchez-Lara, M. Alonso, A. Aristondo, P. P. Ricci, A. Otter, J. Galvan, K. Ruehl, S. Husain, S. Zheng, L. Garcia..<br/>
Some file parsing strategies taken from the [BEMIO project](https://wec-sim.github.io/bemio/).<br/>
Done with the [U++ multiplatform library](https://www.ultimatepp.org/).

## License

Copyright © 2019-2025 Iñaki Zabala, Markel Peñalba, Yerai Peña-Sanchez, Thomas Kelly.

BEMRosetta is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.\
BEMRosetta is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for details. You should have received a copy of the GNU General Public License along with BEMRosetta. If not, see http://www.gnu.org/licenses/.<br/>
<br/>
<br/>
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/logos/BEMRosetta.png" width="200" title="BEMRosetta"></p>

