# BEMRosetta command line utility test

| It has been a long time since any features were added to the command line application. If you want it to be improved, please pull your requests |
| ----------------------------------------------------------------------------------------------------------------------------------------------- |

BEMRosetta_cl includes the features of BEMRosetta that do not need a graphical user interface.

To get the options, run 

```
bemrosetta_cl -help"
```

This demo just loads a Nemoh model (-i --input option), prints some data about the model (-r --report option), and exports it to Wamit format and to FAST format (-c --convert option). The command is:

```
bemrosetta_cl -i "..\..\examples\nemoh\ellipsoid\Nemoh.cal" -r -c ".\testWamit\ellip.1" -c ".\testFAST\ellip.dat"
```

Parameters are done in sequence: if a physical parameter is changed after export, saved files will not include the change. In addition, they can be repeated as desired.

The complete test output is shown below:

```
BEMRosetta\other\test>md ".\testWamit"
BEMRosetta\other\test>md ".\testFAST"
BEMRosetta\other\test>..\..\install\bemrosetta_cl -i "..\..\examples\nemoh\ellipsoid\Nemoh.cal" -r -c ".\testWamit\ellip.1" -c ".\testFAST\ellip.dat"
BEMRosetta Copyright (c) 2019 Iñaki Zabala
Hydrodynamic coefficients converter for Boundary Element Method solver formats
Version beta 2019042317, release, 64 bits

Loading '..\..\examples\nemoh\ellipsoid\Nemoh.cal'
- Hydrostatics file(s) 'Mesh/Hydrostatics*.dat': **Not found**
- KH file(s) 'Mesh/KH*.dat': **Not found**
- Radiation file 'RadiationCoefficients.tec'
- Excitation force file 'ExcitationForce.tec'
- Diffraction force file 'DiffractionForce.tec'
- Froude Krylov file 'FKForce.tec'
- IRF file(s) 'IRF.tec': **Not found**
File '..\..\examples\nemoh\ellipsoid\Nemoh.cal' loaded
Nemoh file '..\..\examples\nemoh\ellipsoid\Nemoh.cal'
g [m/s2]: 9.810, h [m]: INFINITY, rho [kg/m3]: 1000.000 length scale [m]: 1.0
#freqs: 308 (0.030 to 9.240 steps 0.030 [rad/s]), #headings: 1 (0.0 [º])
#bodies: 1
1. 'ellipsoid' dof: 6
- Hydrodynamic coefficients A and B file 'ellip.1'
- Diffraction exciting file 'ellip.3'
- Hydrostatic restoring file 'ellip.hst'
File '.\testWamit\ellip.1' converted
- Hydrodynamic coefficients A and B file 'ellipsoid.1'
- Diffraction exciting file 'ellipsoid.3'
- Hydrostatic restoring file 'ellipsoid.hst'
File '.\testFAST\ellip.dat' converted
BEMRosetta\other\test>pause
Press any key to continue . . .
```
