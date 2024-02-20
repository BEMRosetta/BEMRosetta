// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2024, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"

#include <iostream>
#include <fstream>
#include <unsupported/Eigen/CXX11/Tensor>

using namespace Eigen;


// Based on https://github.com/OpenFAST/matlab-toolbox/blob/main/Utilities/readfile_WND.m
// Bonnie Jonkman, National Renewable Energy Laboratory
bool WNDWind::LoadSum(String fileName, const UVector<String> &varNames, UVector<float> &varValues, float &ZGoffset) {
	FileInLine file(ForceExt(fileName, ".sum"));
	if (!file.IsOpen())
		return false;

	while (std::any_of(varValues.begin(), varValues.end(), [](float d) {return d == 0;})) {		//MFFWS and the TIs should not be zero
		if (file.IsEof())
	    	return false; 		// Reached the end of summary file without all necessary data.');            
	    	
	    String line = ToUpper(file.GetLine());
	
	    int findx = line.FindAfter("=");  // first index
	    if (findx < 0)
	        findx = 0;
	    int lindx = line.GetCount();      // last index
	
		for (int i = 0; i < varValues.size(); ++i) {
			if (varValues[i] == 0) {
				int k = line.Find(varNames[i]);
				if (k >= 0) {
	                k = line.Find("%");
	                if (k >= 0)
	                    lindx = std::max(findx, k-1);
	                
	                String tmp = Trim(line.Mid(findx, lindx - findx + 1));
	                if (tmp[0] == 'T')
	                    varValues[i] = 1;
	                else if (tmp[0] == 'F')
	                    varValues[i] = -1;
	                else
	                    varValues[i] = float(ScanDouble(tmp));
				}
			}
		}
	}
	while (true) {
		if (file.IsEof())
			break;
		
	    String line = ToUpper(file.GetLine());
	
	    int findx = line.Find("HEIGHT OFFSET");	// used to allow grids that are not centered vertically on the turbine hub-height. 
		if (findx >= 0) {
			int lindx = line.GetCount();        
	        int findx = line.FindAfter("=");
	        if (findx < 0)
	            findx = 0;
	        ZGoffset = float(ScanDouble(line.Mid(findx, lindx-findx+1)));	 //z grid offset
	        break;
		}
	}
	
	LHR = false; 		// Default value for TurbSim
	while (true) {
		if (file.IsEof())
			break;
		
	    String line = ToUpper(file.GetLine());
	
	    int findx = line.Find("BLADED LEFT-HAND RULE");	// left-hand rule (to flip sign for v wind component)
		if (findx >= 0)	{	
	        LHR = true;
	        break;
	    }
	}
	return true;
}

String WNDWind::LoadWND(String fileName, double _zHub) {
	FileInBinary file(ForceExt(fileName, ".wnd"));
	if (!file.IsOpen())
		return t_("Impossible to load file");
	
	float ZGoffset = 0;	// used to allow grids that are not centered vertically on the turbine hub-height. 
	ntwr = 0;
	this->fileName = fileName;
	
	try {
		double ConvFact = 1.0; // results in meters and seconds
	
		UVector<String> varNames = {"HUB HEIGHT","CLOCKWISE","UBAR","TI(U)","TI(V)","TI(W)"};  // MUST be in UPPER case
		UVector<float> varValues(varNames.size(), 0);
		
		// READ THE HEADER OF THE BINARY FILE 
	
		float dx;
		
		int nffc  = file.Read<int16>();  	           // number of components
	
		if (nffc != -99) {  // AN OLD-STYLE AERODYN WIND FILE
		    dz    = file.Read<int16>();               // delta z in mm
		    dy    = file.Read<int16>();               // delta y in mm
		    dx    = file.Read<int16>();               // delta x (actually t in this case) in mm
		    nt    = file.Read<int16>();               // half number of time steps
		    mffws = file.Read<int16>();               // 10 times mean FF wind speed, should be equal to MWS
		              file.SeekCur(5*sizeof(int16));  // unnecessary lines
		    nz    = file.Read<int16>();               // 1000 times number of points in vertical direction, max 32
		    ny    = file.Read<int16>();               // 1000 times the number of points in horizontal direction, max 32
		    file.SeekCur(3*(-nffc-1)*sizeof(int16)); 
		
		    // convert the integers to real numbers 
		    nffc  = -nffc;
		    dz    = float(0.001*ConvFact*dz);
		    dy    = float(0.001*ConvFact*dy);
		    dx    = float(0.001*ConvFact*dx);
		    mffws = float(0.1*ConvFact*mffws);
	
		    nz    = static_cast<int>(std::round(std::fmod(nz, std::pow(2, 16))/1000.)); // the mod 2^16 is a work around for somewhat larger grids
		    ny    = static_cast<int>(std::round(std::fmod(ny, std::pow(2, 16))/1000.)); // the mod 2^16 is a work around for somewhat larger grids
		} else {
	    	fc = file.Read<int16>();// should be 4 to allow turbulence intensity to be stored in the header
									// 1 = 1-component von Karman
									// 2 = 1-component Kaimal
									// 3 = 3-component von Karman
									// 4 = improved von Karman
									// 5 = IEC-2 Kaimal
									// 6 = (not supported)
									// 7 = General Kaimal
									// 8 = Mann model
			float TI_U, TI_V, TI_W;
			
	    	if (fc == 4) {
		        nffc     = file.Read<int32>();        	// number of components (should be 3)
		        double lat      = file.Read<float>();	// latitude (deg)
		        double z0       = file.Read<float>();   // Roughness length (m)
		        double zOffset  = file.Read<float>();   // Reference height (m) = Z(1) + GridHeight / 2.0
		        TI_U 		    = file.Read<float>();   // Turbulence Intensity of u component (%)
		        TI_V    		= file.Read<float>();   // Turbulence Intensity of v component (%)
		        TI_W     		= file.Read<float>();   // Turbulence Intensity of w component (%)
		    } else {
		        if (fc > 2)
		            nffc = 3;
		        else
		            nffc = 1;
		        
		        TI_U  = 1;
		        TI_V  = 1;
		        TI_W  = 1;
		        
		        if (fc == 8  //MANN model 
		         || fc == 7) { // General Kaimal      
		            int HeadRec = file.Read<int32>();	// Number of bytes in the header (Nheader)
		            nffc    = file.Read<int32>();  		// Number of turbulence components (1, 2 or 3)
		        }
		    } 
		    dz    = file.Read<float>();  // Grid point spacing z in m 
		    dy    = file.Read<float>();  // Grid point spacing y in m
		    dx    = file.Read<float>();  // Grid point spacing x in m           
		    nt    = file.Read<int32>();  // half the number of time steps
		    mffws = file.Read<float>();  // mean full-field wind speed
		
		    file.SeekCur(3*sizeof(float) + 	// zLu, yLu, xLu: unused variables (for BLADED)
		   				 2*sizeof(int32));	// unused variables (for BLADED) [unused integer, random seed]
		    nz    = file.Read<int32>();  // number of points in vertical direction
		    ny    = file.Read<int32>();  // number of points in horizontal direction
		    if (nffc == 3)
		    	file.SeekCur(2*nffc*sizeof(int32));     // other length scales: unused variables (for BLADED)                
		    
		    if (fc == 7) {
		        float CohDec = file.Read<float>();		// Coherence decay constant
		        float CohLc  = file.Read<float>();		// Coherence scale parameter in m
		    } else if (fc == 8) {        
		        float gamma  = file.Read<float>();      // MANN model shear parameter
		        float Scale  = file.Read<float>();      // MANN model scale length
		                 file.SeekCur(4*sizeof(float) +
		                 			  3*sizeof(int32) +
		                 			  2*sizeof(float) +
		                 			  3*sizeof(int32) +
		                 			  2*sizeof(float));
		    }
	    	varValues[3] = TI_U;
	    	varValues[4] = TI_V;
	    	varValues[5] = TI_W;
		}	
	
		varValues[2] = mffws;
		
		nt = std::max(nt*2, 1);
		dt = dx/mffws;
	                
		// READ THE SUMMARY FILE FOR SCALING FACTORS
		
		if (!LoadSum(fileName, varNames, varValues, ZGoffset)) {	// default values for Bladed
			LHR = true;		
			ZGoffset = 0;
			if (IsNull(_zHub))
				varValues[0] = dz*(nz-1)/2;		// z1 == 0
			else 
				varValues[0] = float(_zHub);
			varValues[1] = 1;		//clockwise rotation
		}
		
		zHub = varValues[0];
		zGrid   = zHub - ZGoffset - dz*(nz-1)/2;  //this is the bottom of the grid
		
		//READ THE GRID DATA FROM THE BINARY FILE
	
		int nv       = nffc*ny*nz;               // the size of one time step
		double factor = 0.00001*mffws;
		UVector<double> Scale = {factor*varValues[3], factor*varValues[4], factor*varValues[5]};
		UVector<double> Offset = {mffws, 0, 0};
		if (LHR) 		// Bladed defined the v-component opposite of the right-hand rule...
		    Scale[1] = -Scale[1];
		
		velocity = Tensor<double, 4>(nt,nffc,ny,nz);
		velocity.setZero();
		
		UVector<int> y_ix;
		if (varValues[1] > 0) //clockwise rotation
		    //flip the y direction....
		    //let's change the dimension of velocity so that it's 4-d instead of 3-d   
		    Arange(y_ix, ny-1, 0, -1);
		else
			Arange(y_ix, 0, ny-1, 1);
		
		Buffer<int16> v(nv);
		for (int it = 0; it < nt; ++it) {
		    file.Read(v, nv*sizeof(int16));
		    
		    int cnt2 = 0;
		    for (int iz = 0; iz < nz; ++iz)
		        for (int iy : y_ix)
		            for (int k = 0; k < nffc; ++k)
		                velocity(it,k,iy,iz) = v[cnt2++]*Scale[k] + Offset[k];
		}
	} catch(Exc e) {
		return e;
	}
	
	return "";
}
