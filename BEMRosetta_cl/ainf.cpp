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

void Hydro::GetK_IRF(double maxT, int numT) {
	if (Nf == 0 || B.IsEmpty())
		return;
	
    Kirf.SetCount(Nb*6); 			
    for (int i = 0; i < Nb*6; ++i) {
    	Kirf[i].SetCount(Nb*6); 			 
   		for (int j = 0; j < Nb*6; ++j)
			Kirf[i][j].setConstant(numT, Null);
    }
		
	double dt = maxT/numT;
	
	GetTirf(Tirf, dt, maxT);
	
	Vector<double> y(Nf);
  	for (int i = 0; i < Nb*6; ++i)
    	for (int j = 0; j < Nb*6; ++j) 
			if (!IsNull(B[i][j][0]))
				GetKirf(Kirf[i][j], Tirf, Get_w(), B[i][j], dt, maxT); 
}  

void Hydro::GetAinf() {
	if (Nf == 0 || A.size() < Nb*6 || !IsLoadedKirf())
		return;	
	
	Awinf.setConstant(Nb*6, Nb*6, Null);
	int numT = int(Tirf.size());
	double dt = Tirf[1] - Tirf[0];
	
    for (int i = 0; i < Nb*6; ++i)
        for (int j = 0; j < Nb*6; ++j)
		    Awinf(i, j) = ::GetAinf(Kirf[i][j], Tirf, Get_w(), A[i][j], dt, Tirf[Tirf.size()-1]);
}

void Hydro::GetAinfw() {
	if (Nf == 0 || A.size() < Nb*6 || !IsLoadedKirf())
		return;	
	
	Ainfw.SetCount(Nb*6); 			
    for (int i = 0; i < Nb*6; ++i) {
    	Ainfw[i].SetCount(Nb*6); 			 
   		for (int j = 0; j < Nb*6; ++j)
			Ainfw[i][j].setConstant(Nf, Null);
    }
    
	double dt = Tirf[1] - Tirf[0];
	
    for (int i = 0; i < Nb*6; ++i)
        for (int j = 0; j < Nb*6; ++j)
		    ::GetAinfw(Ainfw[i][j], Kirf[i][j], Tirf, Get_w(), A[i][j], dt, Tirf[Tirf.size()-1]);
}