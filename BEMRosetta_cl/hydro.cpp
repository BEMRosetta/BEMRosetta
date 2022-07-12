// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include <STEM4U/Integral.h>
#include <STEM4U/Utility.h>
#include "functions.h"
#include "heal.h"


using namespace Eigen;

double Hydro::GetK_IRF_MaxT(const UVector<double> &w) {
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
	
	UVector<double> y(Nf);
  	for (int idf = 0; idf < Nb*6; ++idf) {
    	for (int jdf = 0; jdf < Nb*6; ++jdf) { 
			if (B[idf][jdf].size() == 0 || !IsNum(B[idf][jdf][0])) 
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
	
	Ainf.setConstant(Nb*6, Nb*6, 0);
	//int numT = int(Tirf.size());
	
    for (int i = 0; i < Nb*6; ++i) 
        for (int j = 0; j < Nb*6; ++j) 
            if (IsNum(Kirf[i][j][0]))
		    	Ainf(i, j) = ::GetAinf(Kirf[i][j], Tirf, Get_w(), A[i][j]);
}

void Hydro::GetRAO() {
	if (Nf == 0 || A.size() < Nb*6)
		return;	

	Initialize_Forces(rao);

	MatrixXd D = MatrixXd::Zero(6, 6);
	MatrixXd D2 = MatrixXd::Zero(6, 6);
	
	for (int ib = 0; ib < Nb; ++ib) {
		MatrixXd C = C_(!dimen, ib);
		const MatrixXd &M_ = M[ib];
		for (int ih = 0; ih < Nh; ++ih) {	
			for (int ifr = 0; ifr < Nf; ++ifr) {
				VectorXcd RAO = GetRAO(w[ifr], A_(!dimen, ifr, ib), B_(!dimen, ifr, ib), 
								F_(!dimen, ex, ih, ifr), C, M_, D, D2);
				for (int idf = 0; idf < 6; ++idf)
					rao.Set(ib, ih, ifr, idf, RAO[idf]);
			}
		}
	}
}

VectorXcd Hydro::GetRAO(double w, const MatrixXd &Aw, const MatrixXd &Bw, const VectorXcd &Fwh, 
		const MatrixXd &C, const MatrixXd &M, const MatrixXd &D, const MatrixXd &D2) {
	const std::complex<double> j = std::complex<double>(0, 1);
	VectorXcd RAO = (-sqr(w)*(M + Aw) - j*w*(Bw + D) + C).inverse()*Fwh;
	return RAO;
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
            if (!IsLoadedB(idf, jdf)) 
                continue;
            if (dimen)
		    	::GetAinf_w(Ainf_w[idf][jdf], Kirf[idf][jdf], Tirf, Get_w(), A[idf][jdf]);
            else {
                ::GetAinf_w(Ainf_w[idf][jdf], Kirf[idf][jdf]*g_rho_dim(), Tirf, Get_w(), A_dim(idf, jdf));
                Ainf_w[idf][jdf] *= (1/(rho_dim()*pow(len, GetK_AB(idf, jdf))));
            }
        }
}

void Hydro::GetOgilvieCompliance(bool zremoval, bool thinremoval, bool decayingTail, bool haskind, UVector<int> &vidof, UVector<int> &vjdof) {
	vidof.Clear();
	vjdof.Clear();
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
        MatrixXd ex_hf(Nh, Nf);
        
        for (int jdf = 0; jdf < Nb*6; ++jdf) {
            if (B[idf][jdf].size() == 0 || !IsNum(B[idf][jdf][0])) 
                ;
            else {
	    		if (data.Load(Get_w(), A_dim(idf, jdf), B_dim(idf, jdf), numT, maxT, ex_hf) &&
					data.Heal(zremoval, thinremoval, decayingTail, haskind && idf == jdf)) {
	            	data.Save(Get_w(), A[idf][jdf], Ainf_w[idf][jdf], Ainf(idf, jdf), B[idf][jdf], Tirf, Kirf[idf][jdf]); 
	            	vidof << idf;
	            	vjdof << jdf;
	    		} else
	    			data.Reset(Get_w(), A[idf][jdf], Ainf_w[idf][jdf], Ainf(idf, jdf), B[idf][jdf], Tirf, Kirf[idf][jdf]);
	    		if (dimen) {
	    			dimen = false;
	    			A[idf][jdf] 	= pick(A_ndim(idf, jdf));
	    			Ainf_w[idf][jdf] *= (rho_ndim()/rho_dim());
	    			Ainf(idf, jdf)   *= (rho_ndim()/rho_dim());
	    			B[idf][jdf] 	= pick(B_ndim(idf, jdf));
	    			Kirf[idf][jdf]  = pick(Kirf_ndim(idf, jdf));
	    			dimen = true;
	    		} else {
	    			dimen = true;
	    			A[idf][jdf] 	= pick(A_ndim(idf, jdf));
	    			Ainf_w[idf][jdf] *= (1/(rho_ndim()*pow(len, GetK_AB(idf, jdf))));
	    			Ainf(idf, jdf)   *= (1/(rho_ndim()*pow(len, GetK_AB(idf, jdf))));
	    			B[idf][jdf] 	= pick(B_ndim(idf, jdf));
	    			Kirf[idf][jdf] 	= pick(Kirf_ndim(idf, jdf));
	    			dimen = false;
	    		}
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
				
				for (int idof = 0; idof < 6; ++idof) {		// All dof are available?
					for (int jdof = 0; jdof < 6; ++jdof)
						if (!IsNum(A[ib6 + idof][jb6 + jdof][0])) 
							throw Exc("Coefficient translations requires all DOFs to be available");
				}
				
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
				
				for (int idof = 0; idof < 6; ++idof) {
					for (int jdof = 0; jdof < 6; ++jdof)
						if (!IsNum(A(ib6 + idof, jb6 + jdof))) 
							throw Exc("Coefficient translations requires all DOFs to be available");
				}
				
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
	if (IsLoadedLinearDamping())
		CalcA(linearDamping);
		    
    auto CalcF = [&](auto &ex) {
    	auto exforce = clone(ex.force);
    	
	    for (int ih = 0; ih < Nh; ++ih) {
	    	for (int ib = 0; ib < Nb; ++ib) {
	    		int ib6 = ib*6;
				for (int ifr = 0; ifr < Nf; ++ifr) {
					exforce[ih](ifr, 3 + ib6) += -yg*ex.force[ih](ifr, 2 + ib6) + zg*ex.force[ih](ifr, 1 + ib6);
	    			exforce[ih](ifr, 4 + ib6) += -zg*ex.force[ih](ifr, 0 + ib6) + xg*ex.force[ih](ifr, 2 + ib6);
	    			exforce[ih](ifr, 5 + ib6) += -xg*ex.force[ih](ifr, 1 + ib6) + yg*ex.force[ih](ifr, 0 + ib6);
				}
	    	}
	    }
		ex.force = pick(exforce);
		//GetMaPh(ex);
    };
    
	if (IsLoadedFex())
		CalcF(ex);
	if (IsLoadedFsc())
		CalcF(sc);
	if (IsLoadedFfk())
		CalcF(fk);


    auto CalcQTF = [&](auto &qtf) {
		for (int ib = 0; ib < Nb; ++ib)
		        for (int ih = 0; ih < qh.size(); ++ih) 
					for (int ifr1 = 0; ifr1 < qw.size(); ++ifr1) 
						for (int ifr2 = 0; ifr2 < qw.size(); ++ifr2) {
							std::complex<double> &v0 = qtfdif[ib][ih][0](ifr1, ifr2),
												 &v1 = qtfdif[ib][ih][1](ifr1, ifr2),
												 &v2 = qtfdif[ib][ih][2](ifr1, ifr2),
												 &v3 = qtfdif[ib][ih][3](ifr1, ifr2),
												 &v4 = qtfdif[ib][ih][4](ifr1, ifr2),
												 &v5 = qtfdif[ib][ih][5](ifr1, ifr2);
							v3 += -yg*v2 + zg*v1;
							v4 += -yg*v0 + zg*v2;
							v5 += -yg*v1 + zg*v0;
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

void Hydro::ResetForces(Hydro::FORCE force) {
	if (force == Hydro::FK) {
		if (!IsLoadedFfk())
			return;
		if (!IsLoadedFsc() && !IsLoadedFex())
			return;
		
		if (IsLoadedFsc()) 
			ex = clone(sc);
		else {
			for (int ih = 0; ih < Nh; ++ih) {
				for (int ifr = 0; ifr < Nf; ++ifr) 
					for (int i = 0; i < Nb*6; ++i) 
						if (!IsNull(sc.force[ih](ifr, i))) 
							ex.force[ih](ifr, i) = ex.force[ih](ifr, i) - fk.force[ih](ifr, i);
			}		
		}
		fk.Clear();
	} else if (force == Hydro::SCATTERING) {
		if (!IsLoadedFsc())
			return;
		if (!IsLoadedFfk() && !IsLoadedFex())
			return;
		
		if (IsLoadedFfk()) 
			ex = clone(fk);
		else {
			for (int ih = 0; ih < Nh; ++ih) 
				for (int ifr = 0; ifr < Nf; ++ifr) 
					for (int i = 0; i < Nb*6; ++i) 
						if (!IsNull(sc.force[ih](ifr, i))) 
							ex.force[ih](ifr, i) = ex.force[ih](ifr, i) - sc.force[ih](ifr, i);
		}
		sc.Clear();		
	} else {
		ex.Clear();		
		sc.Clear();		
		fk.Clear();		
	}
}

void Hydro::DeleteFrequencies(const UVector<int> &idFreq) {
	if (idFreq.size() > 0) {
		auto DeleteAB = [&](UArray<UArray<VectorXd>> &A) {
	        UArray<UArray<VectorXd>> An;
		
			An.SetCount(6*Nb);
			for (int ib = 0; ib < 6*Nb; ++ib) {
				An[ib].SetCount(6*Nb);
				for (int jb = 0; jb < 6*Nb; ++jb) {
					An[ib][jb].resize(Nf - idFreq.size());	
					int i = 0, j = 0;
					for (int iif = 0; iif < Nf; ++iif) {
						if (j >= idFreq.size() || iif != idFreq[j])
							An[ib][jb][i++] = A[ib][jb][iif];		
						else 
							j++;
					}
				}
			}
			A = pick(An);
	    };
		
		if (IsLoadedA())
			DeleteAB(A);
		if (IsLoadedAinf_w())
			DeleteAB(Ainf_w);
		if (IsLoadedB())
			DeleteAB(B);
	
		auto DeleteF = [&](Forces &ex) {
	        Forces _ex;
		
			_ex.force.SetCount(Nh);
		    for (int ih = 0; ih < Nh; ++ih) {
		        _ex.force[ih].resize(Nf-idFreq.size(), 6*Nb);
		    	for (int ib = 0; ib < 6*Nb; ++ib) {
					int i = 0, j = 0;
					for (int iif = 0; iif < Nf; ++iif) {
						if (j >= idFreq.size() || iif != idFreq[j]) {
							_ex.force[ih](i, ib) = ex.force[ih](iif, ib);
							i++;
						} else 
							j++;
					}
		    	}
		    }
		    ex = pick(_ex);
	    };	
	
		if (IsLoadedFex())
			DeleteF(ex);
		if (IsLoadedFsc())
			DeleteF(sc);
		if (IsLoadedFfk())
			DeleteF(fk);	
		if (IsLoadedRAO())
			DeleteF(rao);
	
		int j = idFreq.size()-1;	
		for (int i = w.size()-1; i >= 0 && j >= 0; --i) {
			if (i == idFreq[j]) {	
				w.Remove(i);
				T.Remove(i);
				j--;
			}
		}
		Nf = w.size();
	}
}

void Hydro::DeleteFrequenciesQTF(const UVector<int> &idFreqQTF) {
	if (idFreqQTF.size() > 0) {
		UVector<int> vids;
		LinSpaced(vids, int(qw.size()), 0, int(qw.size())-1);
		for (int i = idFreqQTF.size()-1; i >= 0; --i) 
			vids.Remove(idFreqQTF[i]);
		VectorXi ids;
		::Copy(vids, ids);
		qw = VectorXd(qw(ids));
		
		auto DeleteSumDif = [&](UArray<UArray<UArray<MatrixXcd>>> &qtf) {
			for (int ib = 0; ib < Nb; ++ib)
		        for (int ih = 0; ih < qh.size(); ++ih) 
					for (int idf = 0; idf < 6; ++idf) {
						MatrixXcd &m = qtf[ib][ih][idf];
						m = MatrixXcd(m(all, ids));
						m = MatrixXcd(m(ids, all));
					}
		};
		if (IsLoadedQTF()) {
			DeleteSumDif(qtfsum);
			DeleteSumDif(qtfdif);
		}
	}
}

void Hydro::DeleteHeadings(const UVector<int> &idHead) {
	if (idHead.size() > 0) {
		auto DeleteF = [&](Forces &ex) {
			int j = idHead.size()-1;	
			for (int i = head.size()-1; i >= 0 && j >= 0; --i) {
				if (i == idHead[j]) {	
					ex.force.Remove(i);
					j--;
				}
			}
	    };	
	
		if (IsLoadedFex())
			DeleteF(ex);
		if (IsLoadedFsc())
			DeleteF(sc);
		if (IsLoadedFfk())
			DeleteF(fk);	
		if (IsLoadedRAO())
			DeleteF(rao);
	
		int j = idHead.size()-1;	
		for (int i = head.size()-1; i >= 0 && j >= 0; --i) {
			if (i == idHead[j]) {	
				head.Remove(i);
				j--;
			}
		}
		Nh = head.size();
	}
}

void Hydro::DeleteHeadingsQTF(const UVector<int> &idHeadQTF) {
	if (idHeadQTF.size() > 0) {
		UVector<int> vids;
		LinSpaced(vids, int(qh.size()), 0, int(qh.size())-1);
		for (int i = idHeadQTF.size()-1; i >= 0; --i) 
			vids.Remove(idHeadQTF[i]);
		VectorXi ids;
		::Copy(vids, ids);
		qh = VectorXcd(qh(ids));
			
		auto DeleteSumDif = [&](UArray<UArray<UArray<MatrixXcd>>> &qtf) {
			for (int ib = 0; ib < Nb; ++ib)
		        for (int ih = idHeadQTF.size()-1; ih >= 0; --ih)
		        	qtf[ib].Remove(idHeadQTF[ih]);
		};
		if (IsLoadedQTF()) {
			DeleteSumDif(qtfsum);
			DeleteSumDif(qtfdif);
		}
	}
}

void Hydro::FillFrequencyGapsABForces() {
	if (w.size() == 0)
		return;
	
	double dw = std::numeric_limits<double>::max();
	for (int i = 1; i < w.size(); ++i) {
		if (w[i] - w[i-1] < dw)
			dw = w[i] - w[i-1];
	}
	VectorXd ww = Eigen::Map<VectorXd>(w, w.size());
	
	UVector<double> neww, newT;
	neww << w[0];
	newT << T[0];
	for (int i = 1; w[0] + i*dw <= w[w.size()-1]; ++i) {
		neww << w[0] + i*dw;
		newT << T[0] + i*dw/2/M_PI;
	}
	Nf = neww.size();
	
	auto FillAB = [&](UArray<UArray<VectorXd>> &A) {
        UArray<UArray<VectorXd>> An;
	
		An.SetCount(6*Nb);
		for (int ib = 0; ib < 6*Nb; ++ib) {
			An[ib].SetCount(6*Nb);
			for (int jb = 0; jb < 6*Nb; ++jb) {
				An[ib][jb].resize(Nf);	
				for (int i = 0; i < Nf; ++i) 
					An[ib][jb][i] = LinearInterpolate(neww[i], ww, A[ib][jb]);
			}
		}
		A = pick(An);
    };
		
	if (IsLoadedA())
		FillAB(A);
	if (IsLoadedAinf_w())
		FillAB(Ainf_w);
	if (IsLoadedB())
		FillAB(B);
	
	auto FillF = [&](Forces &ex) {
        Forces _ex;
	
		_ex.force.SetCount(Nh);
	    for (int ih = 0; ih < Nh; ++ih) {
	        _ex.force[ih].resize(Nf, 6*Nb);
	    	for (int ib = 0; ib < 6*Nb; ++ib) {
				for (int i = 0; i < Nf; ++i) {
					const VectorXcd &ref = ex.force[ih].col(ib);
					_ex.force[ih](i, ib) = LinearInterpolateC(neww[i], ww, ref);
				}
	    	}
	    }
	    ex = pick(_ex);
    };	

	if (IsLoadedFex())
		FillF(ex);
	if (IsLoadedFsc())
		FillF(sc);
	if (IsLoadedFfk())
		FillF(fk);	
	if (IsLoadedRAO())
		FillF(rao);	
	
	w = pick(neww);
	T = pick(newT);
}

void Hydro::FillFrequencyGapsQTF() {
	if (qw.size() == 0)
		return;
	
	double dw = std::numeric_limits<double>::max();
	for (int i = 1; i < qw.size(); ++i) {
		if (qw[i] - qw[i-1] < dw)
			dw = qw[i] - qw[i-1];
	}	
	UVector<double> neww;
	neww << qw[0];
	for (int i = 1; qw[0] + i*dw <= qw[qw.size()-1]; ++i) 
		neww << qw[0] + i*dw;
	
	auto FillSumDif = [&](UArray<UArray<UArray<MatrixXcd>>> &qtf) {
		for (int ib = 0; ib < Nb; ++ib)
		        for (int ih = 0; ih < qh.size(); ++ih) 
					for (int idof = 0; idof < 6; ++idof) {
						MatrixXcd nm, &m = qtf[ib][ih][idof];
						VectorXd nw;
						Resample(qw, qw, m, nw, nw, nm, dw, dw);
						m = pick(nm);
					}
	};
	if (IsLoadedQTF()) {
		FillSumDif(qtfsum);
		FillSumDif(qtfdif);
	}
	
	::Copy(neww, qw);
}

void Heal();
void Load(const VectorXd &w, const VectorXd &A, const VectorXd &B, double maxT, int num);
void Save(const VectorXd &w, VectorXd &A, VectorXd &Ainfw, double &ainf, VectorXd &B, 
			VectorXd &Tirf, VectorXd &Kinf);
			   				

	
	