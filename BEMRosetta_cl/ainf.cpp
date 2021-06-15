//// Infinite added-mass (Ainf) calculator
//
//   This function calculates the impulse response function and added-mass radiation coefficient 
//	 at the infinite-frequency (Ainf) using the Ogilvie's formula [A], which is 
//   necessary for time-domain simulations of semi-submerged or totally submerged bodies.
//   
//   	    [A] Ogilvie, F. T. Recent Progress Towards the Understanding and Prediction of Ship
//   		5th Symposium on Naval Hydrodynamics, vol. 112, Washington DC, USA, 1964.
//
// Markel Penalba
// Fluid-Mechanics Research Group, Mondragon University

#include "BEMRosetta.h"
#include <STEM4U/Integral.h>
#include "functions.h"

using namespace Eigen;

double Hydro::GetK_IRF_MaxT() {
	if (w.size() < 2)
		return -1;
	double delta = 0;
	int num = 0;
	for (int iw = 1; iw < w.size(); ++iw)
		if (w[iw] != w[iw-1]) {
			delta += w[iw] - w[iw-1];
			num++;
		}
	delta = delta/num;
		
	return M_PI/delta;		// (2*M_PI/delta)/2;
}

void Hydro::K_IRF(double maxT, int numT) {
    Tirf = VectorXd::LinSpaced(numT, 0, maxT);
    
    Kirf.SetCount(Nb*6); 			
    for (int i = 0; i < Nb*6; ++i) {
    	Kirf[i].SetCount(Nb*6); 			 
   		for (int j = 0; j < Nb*6; ++j)
			Kirf[i][j].setConstant(numT, Null);
    }
	
	if (B.IsEmpty())
		return;
	
	Vector<double> y(Nf);
  	for (int i = 0; i < Nb*6; ++i)
    	for (int j = 0; j < Nb*6; ++j) 
			if (!IsNull(B[i][j][0])) 
				for (int it = 0; it < numT; ++it) {
					const VectorXd &b = B[i][j];
					VectorXd &kirf = Kirf[i][j];
					for (int iw = 0; iw < Nf; ++iw)
						y[iw] = b(iw)*cos(w[iw]*Tirf(it));
					kirf(it) = Integral(w, y, SIMPSON_1_3)*2/M_PI;
				}
}  

	
void Hydro::Ainf() {
	if (Nf == 0 || A.size() < Nb*6)
		return;	Awinf.setConstant(Nb*6, Nb*6, Null);
	int numT = int(Tirf.size());
	double dt = Tirf[1] - Tirf[0];
	
    for (int i = 0; i < Nb*6; ++i)
        for (int j = 0; j < Nb*6; ++j)
		    Awinf(i, j) = GetAinf(Kirf[i][j], Tirf, Get_w(), A[i][j], dt, Tirf[Tirf.size()-1]);
}
