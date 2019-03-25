# BEMRosetta
Hydrodynamic coefficients viewer and converter for Boundary Element Method solver formats

<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/deepcwind.jpg" width="500" title="Github Logo"></p>

Boundary Element Methods are extensively used to model hydrodynamic forces in offshore devices like ships, offshore wind platforms and wave energy converters. These solvers use device geometry mesh to get some hydrodynamics coefficients as radiation damping, added mass, wave diffraction force, and wave excitation force. All these data is saved in file formats incompatible between them. These may avoid to use the coefficients between programs. 

BEMRosetta allows to load the hydrodynamic coefficients from a format saving it in another. In addition it allows to compare the results obtained between programs, the results between similar geometries and the same geometry with different discretization levels.

As an additional advantage, BEMRosetta allows to view and visually compare the meshes from different programs.


## Features

* Supported file formats
  * BEM coefficients
    * Load
      * Wamit
      * Nemoh
    * Save	
	
  * Mesh load
    * Wamit .gdf
    * Nemoh .dat	

* Load the hydrodynamic coefficients from one format and save them in another.

* Load the hydrodynamic coefficients for the same geometry from different softwares and compare the results

_Damping for the same geometry got from different solvers_
  
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/2%20solvers%20B.jpg" width="300" title="Damping for the same geometry got from different solvers"></p>

_Excitation force for the same geometry got from different solvers_
  
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/2%20solvers%20exc.jpg" width="300" title="Excitation force for the same geometry got from different solvers"></p>

* Load the hydrodynamic coefficients for the same geometry for different discretization levels and compare the results
* Load the hydrodynamic coefficients for different geometries to compare them

_Damping for different offshore wind floating platforms_
  
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/offshore%20wind%20platforms%20B.png" width="300" title="Damping for different offshore wind floating platforms"></p>

_Excitation force for different offshore wind floating platforms_
  
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/offshore%20wind%20platforms%20exc.jpg" width="300" title="Excitation force for different offshore wind floating platforms"></p>

* Mesh loading, combining them for visual comparison 

<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/deepcwind.jpg" width="300" title="Mesh loading"></p>


## Acknowledgments

Markel Peñalba, João C. C. Henriques, Yerai Peña, Juan C. C. Portillo.<br/>
Some file parsing strategies taken from the [BEMIO project](https://wec-sim.github.io/bemio/).<br/>
Based on the [U++ multiplatform library](https://www.ultimatepp.org/).


## License

Copyright 2019 Iñaki Zabala

BEMRosetta is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.\
BEMRosetta is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for details. You should have received a copy of the GNU General Public License along with BEMRosetta. If not, see http://www.gnu.org/licenses/.<br/>
Rosetta stone image from https://www.iconspng.com/. The pictures are free for personal and even for commercial use. You can modify, copy and distribute the vectors on The Rosetta Stone in iconspng.com. All without asking for permission or setting a link to the source. So, attribution is not required.<br/>
<br/><br/><br/>
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/BEMRosetta/icon256x256.png" width="100" title="BEMRosetta"></p>

