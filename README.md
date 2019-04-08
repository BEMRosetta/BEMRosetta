# BEMRosetta
Hydrodynamic coefficients viewer and converter for Boundary Element Method solver formats

<p align="center">
  <img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/deepcwind.jpg" width="40%" title="DeepCWind in Windows">
  <img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/deepcwind_linux.JPG" width="40%" title="DeepCWind in Linux">
</p>

Boundary Element Methods are extensively used to model hydrodynamic forces in offshore devices like ships, offshore wind platforms and wave energy converters. These solvers use device geometry mesh to get some hydrodynamics coefficients as radiation damping, added mass, wave diffraction force, and wave excitation force. All these data is saved in file formats incompatible between them. These may avoid to use the coefficients between programs. 

BEMRosetta allows to load the hydrodynamic coefficients from a format saving it in another. In addition it allows to compare the results obtained between programs, the results between similar geometries and the same geometry with different discretization levels.

In addition, BEMRosetta allows to view and visually compare the meshes from different programs.

BEMRosetta runs in Windows and Linux and is done in C++, so it does not depend on any scripting language.


## Features
### Supported file formats

* BEM coefficients
  * Load
    * [Wamit](https://www.wamit.com/): .out, .3sc, 3fk, .1, .3, .4, .hst
    * [Nemoh](https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-presentation-192863.kjsp): Nemoh.cal, Mesh/Hydrostatics*.dat, Mesh/KH*.dat, RadiationCoefficients.tec, ExcitationForce.tec, DiffractionForce.tec, FKForce.tec, IRF.tec
    * [FAST-Wamit](https://nwtc.nrel.gov/FAST): HydroDyn.dat, .1, .3, .hst
    * [SeaFEM-Nemoh](http://www.compassis.com/compass/en/Productos/SeaFEM): .flavia.inf, RadiationCoefficients.tec, ExcitationForce.tec, DiffractionForce.tec, FKForce.tec
  * Save
    * [Wamit](https://www.wamit.com/): .1, .3, .hst
    * [FAST-Wamit](https://nwtc.nrel.gov/FAST): HydroDyn.dat, .1, .3, .hst
    * Other solvers that may use these files: [Bladed](https://www.dnvgl.com/services/bladed-3775), [Orcaflex](https://www.orcina.com/SoftwareProducts/OrcaFlex/)    
* Mesh view
  * [Wamit](https://www.wamit.com/): .gdf
  * [Nemoh](https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-mesh-192932.kjsp?RH=1489593406974): .dat	

### Load the hydrodynamic coefficients from one format and save them in another

### Load the hydrodynamic coefficients for the same geometry from different softwares and compare the results

_Damping for the same geometry got from different solvers_
  
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/2%20solvers%20B.jpg" width="300" title="Damping for the same geometry got from different solvers"></p>

_Excitation force for the same geometry got from different solvers_
  
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/2%20solvers%20exc.jpg" width="300" title="Excitation force for the same geometry got from different solvers"></p>

### Load the hydrodynamic coefficients for the same geometry for different discretization levels and compare the results
### Load the hydrodynamic coefficients for different geometries to compare them

_Damping for different offshore wind floating platforms_
  
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/offshore%20wind%20platforms%20B.png" width="300" title="Damping for different offshore wind floating platforms"></p>

_Excitation force for different offshore wind floating platforms_
  
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/offshore%20wind%20platforms%20exc.jpg" width="300" title="Excitation force for different offshore wind floating platforms"></p>

### Mesh loading, combining them for visual comparison 

Several meshes can be loaded in this basic viewer, allowing a visual comparison of geometries.

<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/deepcwind.jpg" width="300" title="Mesh loading"></p>


## Acknowledgments

Markel Peñalba, Juan C. C. Portillo, João C. C. Henriques, Yerai Peña.<br/>
Some file parsing strategies taken from the [BEMIO project](https://wec-sim.github.io/bemio/).<br/>
Based on the [U++ multiplatform library](https://www.ultimatepp.org/).


## New features in process

* Command line version
* Loading files in [Ansys AQWA](https://www.ansys.com/products/structures/ansys-aqwa) format


## License

Copyright (c) 2019 Iñaki Zabala

BEMRosetta is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.\
BEMRosetta is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for details. You should have received a copy of the GNU General Public License along with BEMRosetta. If not, see http://www.gnu.org/licenses/.<br/>
Rosetta stone image from https://www.iconspng.com/. The pictures are free for personal and even for commercial use. You can modify, copy and distribute the vectors on The Rosetta Stone in iconspng.com. All without asking for permission or setting a link to the source. So, attribution is not required.<br/>
<br/><br/><br/>
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/Rosetta.png" width="100" title="BEMRosetta"></p>

