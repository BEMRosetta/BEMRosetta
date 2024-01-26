// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2024, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"

#include <iostream>
#include <fstream>
#include <unsupported/Eigen/CXX11/Tensor>

using namespace Eigen;

// Based on https://github.com/OpenFAST/matlab-toolbox/blob/main/Utilities/readfile_BTS.m
// Bonnie Jonkman, National Renewable Energy Laboratory
void BTSWind::LoadBTSHeader(FileInBinary &file, VectorXf &Vslope, VectorXf &Voffset) {
	turbSimFormat = file.Read<int16>();  

    nz    = file.Read<int32>();      // the number of grid points vertically, INT(4)
    ny    = file.Read<int32>();      // the number of grid points laterally, INT(4)
    ntwr  = file.Read<int32>();      // the number of tower points, INT(4)
    nt    = file.Read<int32>();      // the number of time steps, INT(4)

    dz    = file.Read<float>();      // grid spacing in vertical direction, REAL(4), in m
    dy    = file.Read<float>();      // grid spacing in lateral direction, REAL(4), in m
    dt    = file.Read<float>();      // grid spacing in delta time, REAL(4), in m/s
    mffws = file.Read<float>();      // the mean wind speed at hub height, REAL(4), in m/s
    zHub  = file.Read<float>();      // height of the hub, REAL(4), in m
    zGrid = file.Read<float>();      // height of the bottom of the grid, REAL(4), in m

    Vslope(0)  = file.Read<float>(); // the U-component slope for scaling, REAL(4)
    Voffset(0) = file.Read<float>(); // the U-component offset for scaling, REAL(4)
    Vslope(1)  = file.Read<float>(); // the V-component slope for scaling, REAL(4)
    Voffset(1) = file.Read<float>(); // the V-component offset for scaling, REAL(4)
    Vslope(2)  = file.Read<float>(); // the W-component slope for scaling, REAL(4)
    Voffset(2) = file.Read<float>(); // the W-component offset for scaling, REAL(4)
        
    int nchar = file.Read<int32>();  // the number of characters in the description string, max 200, INT(4)
    StringBuffer str(nchar);
    file.Read(str, nchar); 			 // the ASCII integer representation of the character string
    description = str;
}
	              
String BTSWind::LoadBTS(String fileName) {
	FileInBinary file(ForceExt(fileName, ".bts"));
	if (!file.IsOpen())
		return t_("Impossible to load file");

	this->fileName = fileName;

	try {
	    VectorXf Vslope(3), Voffset(3);
	 	LoadBTSHeader(file, Vslope, Voffset);
	
		String ret;
		int64 pos = file.GetPos();
		if (!IsEmpty(ret = LoadBTSBody<float>(file, Vslope, Voffset))) {// Try with float
			file.Seek(pos);
			ret = LoadBTSBody<int16>(file, Vslope, Voffset);			// If not float, try with int16
		}
		return ret;	
	} catch(Exc e) {
		return e;
	}
}

void BTSWind::SaveBTSHeader(FileOutBinary &file, VectorXf &Vslope, VectorXf &Voffset, int fmtSz) const {
    file.Write(int16(turbSimFormat));  

    file.Write(int32(nz));      // the number of grid points vertically, INT(4)
    file.Write(int32(ny));      // the number of grid points laterally, INT(4)
    file.Write(int32(ntwr));    // the number of tower points, INT(4)
    file.Write(int32(nt));      // the number of time steps, INT(4)

    file.Write(float(dz));      // grid spacing in vertical direction, REAL(4), in m
    file.Write(float(dy));      // grid spacing in lateral direction, REAL(4), in m
    file.Write(float(dt));      // grid spacing in delta time, REAL(4), in m/s
    file.Write(float(mffws));   // the mean wind speed at hub height, REAL(4), in m/s
    file.Write(float(zHub));    // height of the hub, REAL(4), in m
    file.Write(float(zGrid));   // height of the bottom of the grid, REAL(4), in m

    VectorXd mn = VectorXd::Constant(3, std::numeric_limits<double>::max());
    VectorXd mx = VectorXd::Constant(3, std::numeric_limits<double>::lowest());
    for (int it = 0; it < nt; ++it) {
        for (int iz = 0; iz < nz; ++iz)
            for (int iy = 0; iy < ny; ++iy)
                for (int k = 0; k < 3; ++k) {
                    mx(k) = std::max(velocity(it,k,iy,iz), mx(k));
                    mn(k) = std::min(velocity(it,k,iy,iz), mn(k));
                }
			
		if (ntwr > 0) {
            for (int k = 0; k < 3; ++k)      // scale the data
                for (int itw = 0; itw < ntwr; ++itw) {
                	mx(k) = std::max(twrVelocity(it,k,itw), mx(k));     
                	mn(k) = std::min(twrVelocity(it,k,itw), mn(k));     
                }
		}
    }
    double maxType, minType;
    if (fmtSz == 2) {
        maxType = std::numeric_limits<int16>::max();
        minType = std::numeric_limits<int16>::lowest();
    } else {
        maxType = std::numeric_limits<float>::max();
        minType = std::numeric_limits<float>::lowest();
    }
    double rangeType = maxType - minType;
    for (int k = 0; k < 3; ++k) {
        Vslope(k) = rangeType/(mx(k) - mn(k));
       	Voffset(k) = maxType - mx(k)*Vslope(k);
    }	
    file.Write(float(Vslope(0)));  	// the U-component slope for scaling, REAL(4)
    file.Write(float(Voffset(0))); 	// the U-component offset for scaling, REAL(4)
    file.Write(float(Vslope(1))); 	// the U-component offset for scaling, REAL(4)
    file.Write(float(Voffset(1)));  // the U-component slope for scaling, REAL(4)
    file.Write(float(Vslope(2)));  	// the U-component slope for scaling, REAL(4)
    file.Write(float(Voffset(2))); 	// the U-component offset for scaling, REAL(4)
    
    int nchar = description.GetLength();
    file.Write(int32(nchar));
    file.Write(description.begin(), nchar);
}

String BTSWind::SaveBTS(String fileName, int fmtSz) const {
    if (fmtSz < 0) 
        fmtSz = 2;
    
	FileOutBinary file(ForceExt(fileName, ".bts"));
	if (!file.IsOpen())
		return t_("Impossible to open file");

    VectorXf Vslope(3);
    VectorXf Voffset(3);
    
    SaveBTSHeader(file, Vslope, Voffset, fmtSz);
    
    if (fmtSz == 2)
    	SaveBTSBody<int16>(file, Vslope, Voffset);
    else
    	SaveBTSBody<float>(file, Vslope, Voffset);
    return "";
}
