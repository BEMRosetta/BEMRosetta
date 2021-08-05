#include "BEMRosetta.h"
#include <STEM4U/Integral.h>
#include "functions.h"
#include "heal.h"


using namespace Eigen;

double Hydro::GetK_IRF_MaxT(const Vector<double> &w) {
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

double Hydro::GetK_IRF_MaxT() const {
	return GetK_IRF_MaxT(w);
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
		
	GetTirf(Tirf, numT, maxT);
	
	Vector<double> y(Nf);
  	for (int idf = 0; idf < Nb*6; ++idf) {
    	for (int jdf = 0; jdf < Nb*6; ++jdf) { 
			if (B[idf][jdf].size() == 0 || IsNull(B[idf][jdf][0])) 
				continue;
			if (dimen)
				GetKirf(Kirf[idf][jdf], Tirf, Get_w(), B[idf][jdf]);
			else {
				GetKirf(Kirf[idf][jdf], Tirf, Get_w(), B_dim(idf, jdf));
				Kirf[idf][jdf] /= g_rho_dim();
			}
    	}
  	}
}  

void Hydro::GetAinf() {
	if (Nf == 0 || A.size() < Nb*6 || !IsLoadedKirf())
		return;	
	
	Awinf.setConstant(Nb*6, Nb*6, Null);
	int numT = int(Tirf.size());
	
    for (int i = 0; i < Nb*6; ++i)
        for (int j = 0; j < Nb*6; ++j)
		    Awinf(i, j) = ::GetAinf(Kirf[i][j], Tirf, Get_w(), A[i][j]);
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
    
    for (int idf = 0; idf < Nb*6; ++idf)
        for (int jdf = 0; jdf < Nb*6; ++jdf) {
            if (B[idf][jdf].size() == 0 || IsNull(B[idf][jdf][0])) 
                continue;
            if (dimen)
		    	::GetAinfw(Ainfw[idf][jdf], Kirf[idf][jdf], Tirf, Get_w(), A[idf][jdf]);
            else {
                ::GetAinfw(Ainfw[idf][jdf], Kirf[idf][jdf]*g_rho_dim(), Tirf, Get_w(), A_dim(idf, jdf));
                Ainfw[idf][jdf] *= (1/(rho_dim()*pow(len, GetK_AB(idf, jdf))));
            }
        }
}

void Hydro::GetOgilvieCompliance(bool zremoval, bool thinremoval, bool decayingTail) {
	if (Nf == 0 || A.size() < Nb*6)
		return;	
	
	HealBEM data;
	
	if (Ainfw.size() == 0) {
		Ainfw.SetCount(Nb*6); 			
	    for (int idf = 0; idf < Nb*6; ++idf) {
    		Ainfw[idf].SetCount(Nb*6); 			 
   			for (int jdf = 0; jdf < Nb*6; ++jdf)
				Ainfw[idf][jdf].setConstant(Nf, Null);
	    }
    }
    double maxT = min(bem->maxTimeA, Hydro::GetK_IRF_MaxT(w));
    int numT = bem->numValsA;
    
    if (Kirf.size() == 0) {
        Kirf.SetCount(Nb*6); 			
	    for (int idf = 0; idf < Nb*6; ++idf) {
	    	Kirf[idf].SetCount(Nb*6); 			 
	   		for (int jdf = 0; jdf < Nb*6; ++jdf)
				Kirf[idf][jdf].setConstant(numT, Null);
	    }
	}
		
    for (int idf = 0; idf < Nb*6; ++idf) {
        for (int jdf = 0; jdf < Nb*6; ++jdf) {
            if (B[idf][jdf].size() == 0 || IsNull(B[idf][jdf][0])) 
                continue;
    		if (data.Load(Get_w(), A_dim(idf, jdf), B_dim(idf, jdf), numT, maxT)) {
				data.Heal(zremoval, thinremoval, decayingTail);
            	data.Save(Get_w(), A[idf][jdf], Ainfw[idf][jdf], Awinf(idf, jdf), B[idf][jdf], Tirf, Kirf[idf][jdf]); 
    		} else
    			data.Reset(Get_w(), A[idf][jdf], Ainfw[idf][jdf], Awinf(idf, jdf), B[idf][jdf], Tirf, Kirf[idf][jdf]);
    		if (dimen) {
    			dimen = false;
    			A[idf][jdf] = A_ndim(idf, jdf);
    			Ainfw[idf][jdf] *= (rho_ndim()/rho_dim());
    			Awinf(idf, jdf) *= (rho_ndim()/rho_dim());
    			B[idf][jdf] = B_ndim(idf, jdf);
    			Kirf[idf][jdf] = Kirf_ndim(idf, jdf);
    			dimen = true;
    		} else {
    			dimen = true;
    			A[idf][jdf] = A_ndim(idf, jdf);
    			Ainfw[idf][jdf] *= (1/(rho_ndim()*pow(len, GetK_AB(idf, jdf))));
    			Awinf(idf, jdf) *= (1/(rho_ndim()*pow(len, GetK_AB(idf, jdf))));
    			B[idf][jdf] = B_ndim(idf, jdf);
    			Kirf[idf][jdf] = Kirf_ndim(idf, jdf);
    			dimen = false;
    		}
        }
    }
    rao.Reset();	// Previous RAO is now invalid
}

void Heal();
void Load(const VectorXd &w, const VectorXd &A, const VectorXd &B, double maxT, int num);
void Save(const VectorXd &w, VectorXd &A, VectorXd &Ainfw, double &ainf, VectorXd &B, 
			VectorXd &Tirf, VectorXd &Kinf);
			   				

	
	