// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2021, the BEMRosetta author and contributors
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
	
	Ainf.setConstant(Nb*6, Nb*6, Null);
	int numT = int(Tirf.size());
	
    for (int i = 0; i < Nb*6; ++i)
        for (int j = 0; j < Nb*6; ++j)
		    Ainf(i, j) = ::GetAinf(Kirf[i][j], Tirf, Get_w(), A[i][j]);
}

void Hydro::InitAinf_w() {
	Ainf_w.SetCount(Nb*6); 			
    for (int i = 0; i < Nb*6; ++i) {
    	Ainf_w[i].SetCount(Nb*6); 			 
   		for (int j = 0; j < Nb*6; ++j)
			Ainf_w[i][j].setConstant(Nf, Null);
    }
}

void Hydro::GetAinf_w() {
	if (Nf == 0 || A.size() < Nb*6 || !IsLoadedKirf())
		return;	
	
	InitAinf_w();
    
    for (int idf = 0; idf < Nb*6; ++idf)
        for (int jdf = 0; jdf < Nb*6; ++jdf) {
            if (B[idf][jdf].size() == 0 || IsNull(B[idf][jdf][0])) 
                continue;
            if (dimen)
		    	::GetAinf_w(Ainf_w[idf][jdf], Kirf[idf][jdf], Tirf, Get_w(), A[idf][jdf]);
            else {
                ::GetAinf_w(Ainf_w[idf][jdf], Kirf[idf][jdf]*g_rho_dim(), Tirf, Get_w(), A_dim(idf, jdf));
                Ainf_w[idf][jdf] *= (1/(rho_dim()*pow(len, GetK_AB(idf, jdf))));
            }
        }
}

void Hydro::GetOgilvieCompliance(bool zremoval, bool thinremoval, bool decayingTail) {
	if (Nf == 0 || A.size() < Nb*6)
		return;	
	
	HealBEM data;
	
	if (Ainf_w.size() == 0) {
		Ainf_w.SetCount(Nb*6); 			
	    for (int idf = 0; idf < Nb*6; ++idf) {
    		Ainf_w[idf].SetCount(Nb*6); 			 
   			for (int jdf = 0; jdf < Nb*6; ++jdf)
				Ainf_w[idf][jdf].setConstant(Nf, Null);
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
            	data.Save(Get_w(), A[idf][jdf], Ainf_w[idf][jdf], Ainf(idf, jdf), B[idf][jdf], Tirf, Kirf[idf][jdf]); 
    		} else
    			data.Reset(Get_w(), A[idf][jdf], Ainf_w[idf][jdf], Ainf(idf, jdf), B[idf][jdf], Tirf, Kirf[idf][jdf]);
    		if (dimen) {
    			dimen = false;
    			A[idf][jdf] = A_ndim(idf, jdf);
    			Ainf_w[idf][jdf] *= (rho_ndim()/rho_dim());
    			Ainf(idf, jdf)   *= (rho_ndim()/rho_dim());
    			B[idf][jdf] = B_ndim(idf, jdf);
    			Kirf[idf][jdf] = Kirf_ndim(idf, jdf);
    			dimen = true;
    		} else {
    			dimen = true;
    			A[idf][jdf] = A_ndim(idf, jdf);
    			Ainf_w[idf][jdf] *= (1/(rho_ndim()*pow(len, GetK_AB(idf, jdf))));
    			Ainf(idf, jdf)   *= (1/(rho_ndim()*pow(len, GetK_AB(idf, jdf))));
    			B[idf][jdf] = B_ndim(idf, jdf);
    			Kirf[idf][jdf] = Kirf_ndim(idf, jdf);
    			dimen = false;
    		}
        }
    }
    rao.Clear();	// Previous RAO is now invalid
}

void Hydro::GetTranslationTo(double xto, double yto, double zto) {
	double xg = xto - c0(0);
	double yg = yto - c0(1);
	double zg = zto - c0(2);
	
	auto CalcAB = [&](auto &A) {
        auto An = clone(A);
	
		for (int ib = 0; ib < Nb; ++ib) {
			int ib6 = ib*6;
			for (int jb = 0; jb < Nb; ++jb) {
				int jb6 = jb*6;
				
				for (int iif = 0; iif < Nf; ++iif) {
					An[ib6 + 0][jb6 + 3][iif] += - yg*A[ib6 + 0][jb6 + 2][iif] + zg*A[ib6 + 0][jb6 + 1][iif];
					An[ib6 + 1][jb6 + 3][iif] += - yg*A[ib6 + 1][jb6 + 2][iif] + zg*A[ib6 + 1][jb6 + 1][iif];
					An[ib6 + 2][jb6 + 3][iif] += - yg*A[ib6 + 2][jb6 + 2][iif] + zg*A[ib6 + 2][jb6 + 1][iif];

					An[ib6 + 0][jb6 + 4][iif] += - zg*A[ib6 + 0][jb6 + 0][iif] + xg*A[ib6 + 0][jb6 + 2][iif];
					An[ib6 + 1][jb6 + 4][iif] += - zg*A[ib6 + 1][jb6 + 0][iif] + xg*A[ib6 + 1][jb6 + 2][iif];
					An[ib6 + 2][jb6 + 4][iif] += - zg*A[ib6 + 2][jb6 + 0][iif] + xg*A[ib6 + 2][jb6 + 2][iif];

					An[ib6 + 0][jb6 + 5][iif] += - xg*A[ib6 + 0][jb6 + 1][iif] + yg*A[ib6 + 0][jb6 + 0][iif];
					An[ib6 + 1][jb6 + 5][iif] += - xg*A[ib6 + 1][jb6 + 1][iif] + yg*A[ib6 + 1][jb6 + 0][iif];
					An[ib6 + 2][jb6 + 5][iif] += - xg*A[ib6 + 2][jb6 + 1][iif] + yg*A[ib6 + 2][jb6 + 0][iif];

					An[ib6 + 3][jb6 + 0][iif] += - yg*A[ib6 + 2][jb6 + 0][iif] + zg*A[ib6 + 1][jb6 + 0][iif];	
					An[ib6 + 3][jb6 + 1][iif] += - yg*A[ib6 + 2][jb6 + 1][iif] + zg*A[ib6 + 1][jb6 + 1][iif];
					An[ib6 + 3][jb6 + 2][iif] += - yg*A[ib6 + 2][jb6 + 2][iif] + zg*A[ib6 + 1][jb6 + 2][iif];

					An[ib6 + 4][jb6 + 0][iif] += - zg*A[ib6 + 0][jb6 + 0][iif] + xg*A[ib6 + 2][jb6 + 0][iif];	
					An[ib6 + 4][jb6 + 1][iif] += - zg*A[ib6 + 0][jb6 + 1][iif] + xg*A[ib6 + 2][jb6 + 1][iif];
					An[ib6 + 4][jb6 + 2][iif] += - zg*A[ib6 + 0][jb6 + 2][iif] + xg*A[ib6 + 2][jb6 + 2][iif];

					An[ib6 + 5][jb6 + 0][iif] += - xg*A[ib6 + 1][jb6 + 0][iif] + yg*A[ib6 + 0][jb6 + 0][iif];	
					An[ib6 + 5][jb6 + 1][iif] += - xg*A[ib6 + 1][jb6 + 1][iif] + yg*A[ib6 + 0][jb6 + 1][iif];
					An[ib6 + 5][jb6 + 2][iif] += - xg*A[ib6 + 1][jb6 + 2][iif] + yg*A[ib6 + 0][jb6 + 2][iif];

					An[ib6 + 3][jb6 + 3][iif] += -    2*yg*A[ib6 + 2][jb6 + 3][iif] +  2*zg*A[ib6 + 1][jb6 + 3][iif]
												 +   yg*yg*A[ib6 + 2][jb6 + 2][iif] + zg*zg*A[ib6 + 1][jb6 + 1][iif]
												 - 2*yg*zg*A[ib6 + 1][jb6 + 2][iif];   

		    		An[ib6 + 4][jb6 + 4][iif] += -    2*zg*A[ib6 + 0][jb6 + 4][iif] +  2*xg*A[ib6 + 2][jb6 + 4][iif]
												 +   zg*zg*A[ib6 + 0][jb6 + 0][iif] + xg*xg*A[ib6 + 2][jb6 + 2][iif]
												 - 2*zg*xg*A[ib6 + 0][jb6 + 2][iif];

					An[ib6 + 5][jb6 + 5][iif] += -    2*xg*A[ib6 + 1][jb6 + 5][iif] +  2*yg*A[ib6 + 0][jb6 + 5][iif]
												 +   xg*xg*A[ib6 + 1][jb6 + 1][iif] + yg*yg*A[ib6 + 0][jb6 + 0][iif]
												 - 2*xg*yg*A[ib6 + 0][jb6 + 1][iif];
				}
				An[ib6 + 3][jb6 + 0] = clone(A[ib6 + 0][jb6 + 3]);
				An[ib6 + 3][jb6 + 1] = clone(A[ib6 + 1][jb6 + 3]);
				An[ib6 + 3][jb6 + 2] = clone(A[ib6 + 2][jb6 + 3]);

				An[ib6 + 4][jb6 + 0] = clone(A[ib6 + 0][jb6 + 4]);
				An[ib6 + 4][jb6 + 1] = clone(A[ib6 + 1][jb6 + 4]);
				An[ib6 + 4][jb6 + 2] = clone(A[ib6 + 2][jb6 + 4]);
			}
		}
		A = pick(An);
    };
	
	if (IsLoadedA())
		CalcAB(A);
	if (IsLoadedAinf_w())
		CalcAB(Ainf_w);
	if (IsLoadedB())
		CalcAB(B);
    
    auto CalcA = [&](auto &A) {
        auto An = clone(A);
	
		for (int ib = 0; ib < Nb; ++ib) {
			int ib6 = ib*6;
			for (int jb = 0; jb < Nb; ++jb) {
				int jb6 = jb*6;
				
				An(ib6 + 0, jb6 + 3) += - yg*A(ib6 + 0, jb6 + 2) + zg*A(ib6 + 0, jb6 + 1);
				An(ib6 + 1, jb6 + 3) += - yg*A(ib6 + 1, jb6 + 2) + zg*A(ib6 + 1, jb6 + 1);
				An(ib6 + 2, jb6 + 3) += - yg*A(ib6 + 2, jb6 + 2) + zg*A(ib6 + 2, jb6 + 1);
	
				An(ib6 + 0, jb6 + 4) += - zg*A(ib6 + 0, jb6 + 0) + xg*A(ib6 + 0, jb6 + 2);
				An(ib6 + 1, jb6 + 4) += - zg*A(ib6 + 1, jb6 + 0) + xg*A(ib6 + 1, jb6 + 2);
				An(ib6 + 2, jb6 + 4) += - zg*A(ib6 + 2, jb6 + 0) + xg*A(ib6 + 2, jb6 + 2);
	
				An(ib6 + 0, jb6 + 5) += - xg*A(ib6 + 0, jb6 + 1) + yg*A(ib6 + 0, jb6 + 0);
				An(ib6 + 1, jb6 + 5) += - xg*A(ib6 + 1, jb6 + 1) + yg*A(ib6 + 1, jb6 + 0);
				An(ib6 + 2, jb6 + 5) += - xg*A(ib6 + 2, jb6 + 1) + yg*A(ib6 + 2, jb6 + 0);
	
				An(ib6 + 3, jb6 + 0) += - yg*A(ib6 + 2, jb6 + 0) + zg*A(ib6 + 1, jb6 + 0);	
				An(ib6 + 3, jb6 + 1) += - yg*A(ib6 + 2, jb6 + 1) + zg*A(ib6 + 1, jb6 + 1);
				An(ib6 + 3, jb6 + 2) += - yg*A(ib6 + 2, jb6 + 2) + zg*A(ib6 + 1, jb6 + 2);
	
				An(ib6 + 4, jb6 + 0) += - zg*A(ib6 + 0, jb6 + 0) + xg*A(ib6 + 2, jb6 + 0);	
				An(ib6 + 4, jb6 + 1) += - zg*A(ib6 + 0, jb6 + 1) + xg*A(ib6 + 2, jb6 + 1);
				An(ib6 + 4, jb6 + 2) += - zg*A(ib6 + 0, jb6 + 2) + xg*A(ib6 + 2, jb6 + 2);
	
				An(ib6 + 5, jb6 + 0) += - xg*A(ib6 + 1, jb6 + 0) + yg*A(ib6 + 0, jb6 + 0);	
				An(ib6 + 5, jb6 + 1) += - xg*A(ib6 + 1, jb6 + 1) + yg*A(ib6 + 0, jb6 + 1);
				An(ib6 + 5, jb6 + 2) += - xg*A(ib6 + 1, jb6 + 2) + yg*A(ib6 + 0, jb6 + 2);
				
				An(ib6 + 3, jb6 + 3) += -    2*yg*A(ib6 + 2, jb6 + 3) +  2*zg*A(ib6 + 1, jb6 + 3)
									 	+   yg*yg*A(ib6 + 2, jb6 + 2) + zg*zg*A(ib6 + 1, jb6 + 1)
										- 2*yg*zg*A(ib6 + 1, jb6 + 2);   
	
	    		An(ib6 + 4, jb6 + 4) += -    2*zg*A(ib6 + 0, jb6 + 4) +  2*xg*A(ib6 + 2, jb6 + 4)
										+   zg*zg*A(ib6 + 0, jb6 + 0) + xg*xg*A(ib6 + 2, jb6 + 2)
										- 2*zg*xg*A(ib6 + 0, jb6 + 2);
	
				An(ib6 + 5, jb6 + 5) += -    2*xg*A(ib6 + 1, jb6 + 5) +  2*yg*A(ib6 + 0, jb6 + 5)
										+   xg*xg*A(ib6 + 1, jb6 + 1) + yg*yg*A(ib6 + 0, jb6 + 0)
										- 2*xg*yg*A(ib6 + 0, jb6 + 1);
			}
		}
		A = pick(An);
    };
    
    if (IsLoadedA0())
		CalcA(A0);
    if (IsLoadedAinf())
		CalcA(Ainf);
	if (IsLoadedDlin())
		CalcA(Dlin);
		    
    auto CalcF = [&](auto &ex) {
    	auto exnre = clone(ex.re);
    	auto exnim = clone(ex.im);
    	
	    for (int ih = 0; ih < Nh; ++ih) {
	    	for (int ib = 0; ib < Nb; ++ib) {
	    		int ib6 = ib*6;
				for (int ifr = 0; ifr < Nf; ++ifr) {
					exnre[ih](ifr, 3 + ib6) += -yg*ex.re[ih](ifr, 2 + ib6) + zg*ex.re[ih](ifr, 1 + ib6);
					exnim[ih](ifr, 3 + ib6) += -yg*ex.im[ih](ifr, 2 + ib6) + zg*ex.im[ih](ifr, 1 + ib6);
	    			exnre[ih](ifr, 4 + ib6) += -zg*ex.re[ih](ifr, 0 + ib6) + xg*ex.re[ih](ifr, 2 + ib6);
	    			exnim[ih](ifr, 4 + ib6) += -zg*ex.im[ih](ifr, 0 + ib6) + xg*ex.im[ih](ifr, 2 + ib6);
	    			exnre[ih](ifr, 5 + ib6) += -xg*ex.re[ih](ifr, 1 + ib6) + yg*ex.re[ih](ifr, 0 + ib6);
					exnim[ih](ifr, 5 + ib6) += -xg*ex.im[ih](ifr, 1 + ib6) + yg*ex.im[ih](ifr, 0 + ib6);
				}
	    	}
	    }
		ex.re = pick(exnre);
		ex.im = pick(exnim);
		GetMaPh(ex);
    };
    
    if (IsLoadedFex())
    	CalcF(ex);
	if (IsLoadedFsc())
		CalcF(sc);
	if (IsLoadedFfk())
		CalcF(fk);


    auto CalcQTF = [&](auto &qtf) {
	    for (int i = 0; i < qtf.size(); ++i) {
	        QTF &q = qtf[i];
			q.fre[3] += -yg*q.fre[2] + zg*q.fre[1];
			q.fim[3] += -yg*q.fim[2] + zg*q.fim[1];
			q.fma[3] = sqrt(sqr(q.fre[3]) + sqr(q.fim[3]));
			q.fph[3] = atan2(q.fim[3], q.fre[3]);
			
			q.fre[4] += -yg*q.fre[0] + zg*q.fre[2];
			q.fim[4] += -yg*q.fim[0] + zg*q.fim[2];
			q.fma[4] = sqrt(sqr(q.fre[4]) + sqr(q.fim[4]));
			q.fph[4] = atan2(q.fim[4], q.fre[4]);
			
			q.fre[5] += -yg*q.fre[1] + zg*q.fre[0];
			q.fim[5] += -yg*q.fim[1] + zg*q.fim[0];
			q.fma[5] = sqrt(sqr(q.fre[5]) + sqr(q.fim[5]));
			q.fph[5] = atan2(q.fim[5], q.fre[5]);
	    }
    };

	if (IsLoadedQTF()) {
		CalcQTF(qtfsum);		
		CalcQTF(qtfdif);
	}

	
	c0(0) = xto;
	c0(1) = yto;
	c0(2) = zto;
	
	// Some previous data are now invalid
	Kirf.Clear();
	rao.Clear();	
	C.Clear();
	
	
	if (!AfterLoad()) {
		String error = GetLastError();
		throw Exc(Format(t_("Problem translating model: '%s'\n%s"), error));	
	}
}

	
void Heal();
void Load(const VectorXd &w, const VectorXd &A, const VectorXd &B, double maxT, int num);
void Save(const VectorXd &w, VectorXd &A, VectorXd &Ainfw, double &ainf, VectorXd &B, 
			VectorXd &Tirf, VectorXd &Kinf);
			   				

	
	