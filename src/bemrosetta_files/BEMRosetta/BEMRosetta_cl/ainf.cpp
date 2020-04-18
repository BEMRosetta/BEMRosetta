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

void Hydro::K_IRF(double maxT, int numT) {
	LinSpaced(Tirf, numT, 0., maxT);
    
    Kirf.SetCount(numT); 			
	for (int it = 0; it < numT; ++it) 
		Kirf[it].setConstant(Nb*6, Nb*6, Null);
	
	Buffer<double> y(Nf);
    for (int it = 0; it < numT; ++it) 
        for (int i = 0; i < Nb*6; ++i)
        	for (int j = 0; j < Nb*6; ++j) {
        		if (!B.IsEmpty() && !IsNull(B[0](i, j))) {
        			for (int iw = 0; iw < Nf; ++iw)
        				y[iw] = B[iw](i, j)*cos(w[iw]*Tirf[it]);
        			double kirf = 0;
        			for (int iw = 1; iw < Nf; ++iw)
            			kirf += Avg(y[iw-1], y[iw])*(w[iw] - w[iw-1]);
        			Kirf[it](i, j) = kirf*2/M_PI;
        		}
    		}
}
	
void Hydro::Ainf() {
	if (Nf == 0)
		return;
	Awinf.setConstant(Nb*6, Nb*6, Null);
	int numT = Tirf.GetCount();
	double dt = Tirf[1] - Tirf[0];
	
	Buffer<double> y(numT);
    for (int i = 0; i < Nb*6; ++i)
        for (int j = 0; j < Nb*6; ++j) {
		    double awinf = 0;
		    bool isnull = false;
		    for (int iw = 0; iw < Nf; ++iw) {
		        for (int it = 0; it < numT; ++it) {
		            if (IsNull(Kirf[it](i, j))) {
		                isnull = true;
		                break;
		            }
		            y[it] = Kirf[it](i, j)*sin(w[iw]*Tirf[it]);
		        }
		        if (isnull)
		            break;
				double kint = 0;
		        for (int it = 1; it < numT; ++it) 
		            kint += Avg(y[it-1], y[it])*dt;
		        // Ogilvie's formula
		        awinf += A[iw](i, j) + kint/w[iw];
			}
			if (isnull)
				Awinf(i, j) = Null;
			else
		    	Awinf(i, j) = awinf/Nf;
        }
}