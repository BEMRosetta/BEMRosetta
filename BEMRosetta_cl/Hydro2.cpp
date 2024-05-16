// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2024, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#include <ScatterDraw/DataSource.h>
#include <ScatterDraw/Equation.h>
#include "functions.h"
#include <SysInfo/Crash.h>
#include <STEM4U/SeaWaves.h>
#include <MatIO/matio.h>

using namespace Upp;
using namespace Eigen;


const char *Hydro::strDataToPlot[] = {t_("A(ω)"), t_("A∞"), t_("A₀"), t_("B(ω)"), t_("A∞(ω)"), t_("Kirf"),
				t_("|Fsc|"), t_("arg(Fsc)"), t_("|Ffk|"), t_("arg(Ffk)"), t_("|Fex|"), t_("arg(Fex)"),
				t_("|RAO|"), t_("arg(RAO)"), t_("|Z|"), t_("arg(Z)"), t_("|Kr|"), t_("arg(Kr)"), 
				t_("|TFS|"), t_("arg(TFS)")};

// enum BEM_FMT 					  {WAMIT, 		  WAMIT_1_3, 					 FAST_WAMIT, 				 	  HAMS Wamit    HAMS,   WADAM_WAMIT,   NEMOH,     NEMOHv115, 	NEMOHv3, 	SEAFEM_NEMOH,  AQWA,   					  FOAMM,   DIODORE,	  BEMROSETTA, 	     ORCAWAVE,   CSV_MAT,    CSV_TABLE,    BEMIOH5,		CAPYTAINE, HYDROSTAR,   CAPY.nc, 		 ORCAWAVE.owr,		UNKNOWN, NUMBEM};
const char *Hydro::bemStr[]         = {"Wamit .out", "Wamit .1.2.3.hst.7.8.9.ss.12", "FAST .dat.1.2.3.hst.789.ss.12", "HAMS Wamit", "HAMS", "Wadam Wamit","Nemoh v2", "Nemoh v115", "Nemoh v3", "SeaFEM Nemoh","AQWA .lis .ah1 .qtf", 	  "FOAMM", "Diodore", "BEMRosetta .bem", "OrcaWave", ".csv mat", ".csv table", "BEMIO .h5",	".dat",    ".out", 		"Capytaine .nc", 
#ifdef PLATFORM_WIN32	
"OrcaWave .owr", 	
#endif
				"By extension"};
const bool Hydro::bemCanSave[] 		= {true, 	      true,	     				     true,		 			 	 	  false,		false,  false,		   false,     false,	 	false,	   	false, 		   true,  					  false,   true,	  true,			     false,	     true, 	     true, 		   true,		false, 	   false, 		false,			 
#ifdef PLATFORM_WIN32	
false,
#endif
				true};       
const char *Hydro::bemExt[]	   		= {"*.out", 	  "*.1",	     				 "*.1",		 			 	      "",		   	"",	    "",		       "",        "", 		   	"",			"",			   "*.qtf", 				  "",      "*.hdb",	  "*.bem",		   	 "*.yml",	 "*.csv",    "*.csv", 	   "*.h5",		"*.dat",   "*.out", 	"*.nc", 		 
#ifdef PLATFORM_WIN32	
"*.owr",	
#endif		
				"*.*"};       
	
const bool Hydro::caseCanSave[]     = {false, 	      false,	     				 false,		 			 	 	  false,		true,	false,		   true,      true,	 		true,	   	false, 		   false,  					  false,   false,	  false,			 false,	     false, 	 false, 	   false,		false, 	   false, 		false,			 
#ifdef PLATFORM_WIN32	
false,
#endif
				true};

int Hydro::idCount = 0;	

void Hydro::Initialize_Sts() {
	dt.sts.SetCount(6*dt.Nb);
	for (int ib = 0; ib < 6*dt.Nb; ++ib) 
		dt.sts[ib].SetCount(6*dt.Nb);
}

void Hydro::GetFexFromFscFfk() {
	Initialize_Forces(dt.ex);
	for (int ih = 0; ih < dt.Nh; ++ih) {
		for (int ifr = 0; ifr < dt.Nf; ++ifr) {
			for (int i = 0; i < dt.Nb*6; ++i) {
				if (IsNum(dt.sc.force[ih](ifr, i)) && IsNum(dt.fk.force[ih](ifr, i))) 
					dt.ex.force[ih](ifr, i) = dt.sc.force[ih](ifr, i) + dt.fk.force[ih](ifr, i);
			}
		}
	}
}

void Hydro::GetFscFromFexFfk() {
	Initialize_Forces(dt.sc);
	for (int ih = 0; ih < dt.Nh; ++ih) {
		for (int ifr = 0; ifr < dt.Nf; ++ifr) {
			for (int i = 0; i < dt.Nb*6; ++i) {
				if (IsNum(dt.ex.force[ih](ifr, i)) && IsNum(dt.fk.force[ih](ifr, i))) 
					dt.sc.force[ih](ifr, i) = dt.ex.force[ih](ifr, i) - dt.fk.force[ih](ifr, i);
			}
		}
	}	
}

void Hydro::GetFfkFromFexFsc() {
	Initialize_Forces(dt.fk);
	for (int ih = 0; ih < dt.Nh; ++ih) {
		for (int ifr = 0; ifr < dt.Nf; ++ifr) {
			for (int i = 0; i < dt.Nb*6; ++i) {
				if (IsNum(dt.ex.force[ih](ifr, i)) && IsNum(dt.sc.force[ih](ifr, i))) 
					dt.fk.force[ih](ifr, i) = dt.ex.force[ih](ifr, i) - dt.sc.force[ih](ifr, i);
			}
		}
	}	
}

void Hydro::Normalize() {
	if (IsLoadedC()) {
		for (int ib = 0; ib < dt.Nb; ++ib) {
			for (int idf = 0; idf < 6; ++idf) 
				for (int jdf = 0; jdf < 6; ++jdf) 
					dt.msh[ib].dt.C(idf, jdf) = C_ndim(ib, idf, jdf);
		}
	}
	if (IsLoadedA() && IsLoadedB()) {
		for (int ifr = 0; ifr < dt.Nf; ++ifr) {
			for (int idf = 0; idf < 6*dt.Nb; ++idf) {
				for (int jdf = 0; jdf < 6*dt.Nb; ++jdf) {	
					dt.A[idf][jdf][ifr] = A_ndim(ifr, idf, jdf);
					dt.B[idf][jdf][ifr] = B_ndim(ifr, idf, jdf);
				}
			}
		}
	}
	if (IsLoadedAinf()) {
		for (int i = 0; i < 6*dt.Nb; ++i) 
			for (int j = 0; j < 6*dt.Nb; ++j) 
				dt.Ainf(i, j) = Ainf_ndim(i, j);
	}
	if (IsLoadedA0()) {
		for (int i = 0; i < 6*dt.Nb; ++i) 
			for (int j = 0; j < 6*dt.Nb; ++j) 
				dt.A0(i, j) = A0_ndim(i, j);
	}
	if (IsLoadedFex())
    	Normalize_Forces(dt.ex);
	if (IsLoadedFsc())
		Normalize_Forces(dt.sc);
	if (IsLoadedFfk())
		Normalize_Forces(dt.fk);
	if (IsLoadedRAO()) 
		Normalize_RAO(dt.rao);
}

void Hydro::Dimensionalize() {
	if (IsLoadedC()) {
		for (int ib = 0; ib < dt.Nb; ++ib) {
			for (int idf = 0; idf < 6; ++idf) 
				for (int jdf = 0; jdf < 6; ++jdf) 
					dt.msh[ib].dt.C(idf, jdf) = C_dim(ib, idf, jdf);
		}
	}
	if (IsLoadedA() && IsLoadedB()) {
		for (int ifr = 0; ifr < dt.Nf; ++ifr) {
			for (int idf = 0; idf < 6*dt.Nb; ++idf) {
				for (int jdf = 0; jdf < 6*dt.Nb; ++jdf) {	
					dt.A[idf][jdf][ifr] = A_dim(ifr, idf, jdf);
					dt.B[idf][jdf][ifr] = B_dim(ifr, idf, jdf);
				}
			}
		}
	}
	if (IsLoadedAinf()) {	
		for (int i = 0; i < 6*dt.Nb; ++i) 
			for (int j = 0; j < 6*dt.Nb; ++j) 
				dt.Ainf(i, j) = Ainf_dim(i, j);
	}
	if (IsLoadedA0()) {	
		for (int i = 0; i < 6*dt.Nb; ++i) 
			for (int j = 0; j < 6*dt.Nb; ++j) 
				dt.A0(i, j) = A0_ndim(i, j);
	}
	if (IsLoadedFex())
    	Dimensionalize_Forces(dt.ex);
	if (IsLoadedFsc())
		Dimensionalize_Forces(dt.sc);
	if (IsLoadedFfk())
		Dimensionalize_Forces(dt.fk);
	if (IsLoadedRAO()) 
		Dimensionalize_RAO(dt.rao);
}

// 		-180	0	180
//		180		0	180
//		180		0
//		1		0

// 		0	180		360
//		0	180		0
//		0	180
//		0	1

void Hydro::SortHeadings() {
	UVector<double> nhead;
	for (int ih = 0; ih < dt.head.size(); ++ih)
		FindAdd(nhead, FixHeading_0_360(dt.head[ih]));
	
	UVector<int> indices_h = GetSortOrderX(nhead);
	SetSortOrder(nhead, indices_h);
	dt.head = pick(nhead);
	dt.Nh = indices_h.size();

	UArray<std::complex<double>> nmhead;
	for (int ih = 0; ih < dt.mdhead.size(); ++ih)
		FindAdd(nmhead, FixHeading_0_360(dt.mdhead[ih]));
	
	UVector<int> indices_md = GetSortOrderX(nmhead, SortComplex);
	
	dt.mdhead.resize(nmhead.size());		// Cannot pick()
	for (int ih = 0; ih < dt.mdhead.size(); ++ih)
		dt.mdhead[ih] = nmhead[indices_md[ih]];
		
	auto SortF = [&](Forces &F) {
		Forces f = clone(F);
		for (int ih = 0; ih < dt.Nh; ++ih) 
			for (int ifr = 0; ifr < dt.Nf; ++ifr) 
				for (int idf = 0; idf < 6*dt.Nb; ++idf) 
					F.force[ih](ifr, idf) = f.force[indices_h[ih]](ifr, idf);
		F.force.SetCount(dt.Nh);
	};

	auto SortMD = [&](UArray<UArray<UArray<VectorXd>>> &MD) {
		UArray<UArray<UArray<VectorXd>>> _md = clone(MD);
		for (int ib = 0; ib < dt.Nb; ++ib) {
			MD[ib].SetCount(int(dt.mdhead.size()));
			for (int ih = 0; ih < dt.mdhead.size(); ++ih)
				for (int idf = 0; idf < 6; ++idf) 
					for (int ifr = 0; ifr < dt.Nf; ++ifr) 		
						MD[ib][ih][idf][ifr] = _md[ib][indices_md[ih]][idf][ifr];
		}
	};
			
	if (IsLoadedFex())
		SortF(dt.ex);
	if (IsLoadedFsc())
		SortF(dt.sc);
	if (IsLoadedFfk())
		SortF(dt.fk);
	if (IsLoadedRAO()) 
		SortF(dt.rao);	

	if(IsLoadedMD())
		SortMD(dt.md);

	UArray<std::complex<double>> nqh;
	for (int ih = 0; ih < dt.qh.size(); ++ih)
		FindAdd(nqh, FixHeading_0_360(dt.qh[ih]));
	
	UVector<int> indices_qw = GetSortOrderX(nqh, SortComplex);
	
	dt.qh.resize(nqh.size());
	for (int ih = 0; ih < dt.qh.size(); ++ih)
		dt.qh[ih] = nqh[ih];
	
	auto SortQTF = [&](UArray<UArray<UArray<MatrixXcd>>> &QTF) {
		UArray<UArray<UArray<MatrixXcd>>> qtf = clone(QTF);
		for (int ib = 0; ib < dt.Nb; ++ib) {
			QTF[ib].SetCount(int(dt.qh.size()));
	        for (int ih = 0; ih < dt.qh.size(); ++ih) 
	        	for (int idf = 0; idf < 6; ++idf) 
					for (int ifr1 = 0; ifr1 < dt.qw.size(); ++ifr1) 
						for (int ifr2 = 0; ifr2 < dt.qw.size(); ++ifr2) 
							QTF[ib][ih][idf](ifr1, ifr2) = qtf[ib][indices_qw[ih]][idf](ifr1, ifr2);
		}
	};
		
	if (IsLoadedQTF(true))
		SortQTF(dt.qtfsum);
	if (IsLoadedQTF(false))
		SortQTF(dt.qtfdif);
}

void Hydro::SortFrequencies() {
	if (!IsSorted(dt.w)) {
		UVector<int> indices = GetSortOrderX(dt.w);
		dt.w = ApplyIndex(dt.w, indices);
		//T = ApplyIndex(T, indices);
	
		auto SortAB = [&](UArray<UArray<VectorXd>> &_A) {
			UArray<UArray<VectorXd>> a = clone(_A);
			for (int ifr = 0; ifr < dt.Nf; ++ifr) 
				for (int idf = 0; idf < 6*dt.Nb; ++idf) 
					for (int jdf = 0; jdf < 6*dt.Nb; ++jdf) 	
						_A[idf][jdf][ifr] = a[idf][jdf][indices[ifr]];
		};
	
		auto SortF = [&](Forces &F) {
			Forces f = clone(F);
			for (int ih = 0; ih < dt.Nh; ++ih) 
				for (int ifr = 0; ifr < dt.Nf; ++ifr) 
					for (int idf = 0; idf < 6*dt.Nb; ++idf) 
						F.force[ih](ifr, idf) = f.force[ih](indices[ifr], idf);
		};
		
		auto SortMD = [&](UArray<UArray<UArray<VectorXd>>> &MD) {
			UArray<UArray<UArray<VectorXd>>> _md = clone(MD);
			for (int ib = 0; ib < dt.Nb; ++ib) 
				for (int ih = 0; ih < dt.mdhead.size(); ++ih)
					for (int idf = 0; idf < 6; ++idf) 
						for (int ifr = 0; ifr < dt.Nf; ++ifr) 		
							MD[ib][ih][idf][ifr] = _md[ib][ih][idf][indices[ifr]];
		};
		
		if (IsLoadedA()) 
			SortAB(dt.A);
		if (IsLoadedB()) 
			SortAB(dt.B);
		
		if (IsLoadedFex())
			SortF(dt.ex);
		if (IsLoadedFsc())
			SortF(dt.sc);
		if (IsLoadedFfk())
			SortF(dt.fk);
		if (IsLoadedRAO()) 
			SortF(dt.rao);
		
		if(IsLoadedMD())
			SortMD(dt.md);
	}
	if (!IsSorted(dt.qw)) {
		UVector<int> indices = GetSortOrderX(dt.qw);
		dt.qw = ApplyIndex(dt.qw, indices);
		
		auto SortQTF = [&](UArray<UArray<UArray<MatrixXcd>>> &QTF) {
			UArray<UArray<UArray<MatrixXcd>>> qtf = clone(QTF);
			for (int ib = 0; ib < dt.Nb; ++ib) 
		        for (int ih = 0; ih < dt.qh.size(); ++ih) 
		        	for (int idf = 0; idf < 6; ++idf) 
						for (int ifr1 = 0; ifr1 < dt.qw.size(); ++ifr1) 
							for (int ifr2 = 0; ifr2 < dt.qw.size(); ++ifr2) 
								QTF[ib][ih][idf](ifr1, ifr2) = qtf[ib][ih][idf](indices[ifr1], indices[ifr2]);
		};
			
		if (IsLoadedQTF(true))
			SortQTF(dt.qtfsum);
		if (IsLoadedQTF(false))
			SortQTF(dt.qtfdif);
	}
}

void Hydro::Initialize_AB(UArray<UArray<VectorXd>> &a, double val) {
	a.SetCount(6*dt.Nb);
	for (int i = 0; i < 6*dt.Nb; ++i) {
		a[i].SetCount(6*dt.Nb);
		for (int j = 0; j < 6*dt.Nb; ++j) 
			a[i][j].setConstant(dt.Nf, val);
	}
}

void Hydro::Initialize_ABpan(UArray<UArray<UArray<UArray<UArray<double>>>>> &a, double val) {
	a.SetCount(dt.Nb);
	for (int ib = 0; ib < dt.Nb; ++ib) {
		a[ib].SetCount(dt.pots_rad[ib].size());
		for (int idp = 0; idp < dt.pots_rad[ib].size(); ++idp) {
			a[ib][idp].SetCount(6);
			for (int idf1 = 0; idf1 < 6; ++idf1) {
				a[ib][idp][idf1].SetCount(6);
				for (int idf2 = 0; idf2 < 6; ++idf2) 
					a[ib][idp][idf1][idf2].SetCount(dt.Nf, val);
			}
		}
	}
}

void Hydro::Initialize_PotsRad() {
	dt.pots_rad.SetCount(dt.Nb);
	for (int ib = 0; ib < dt.Nb; ++ib) {
		if (dt.pots_rad[ib].IsEmpty())
			dt.pots_rad[ib].SetCount(dt.msh[ib].dt.mesh.GetNumPanels());
		for (int ipot = 0; ipot < dt.pots_rad[ib].size(); ++ipot) {
			dt.pots_rad[ib][ipot].SetCount(6);
			for (int idf = 0; idf < 6; ++idf)
				dt.pots_rad[ib][ipot][idf].SetCount(dt.Nf, 0);
		}
	}
}

void Hydro::Initialize_PotsIncDiff(UArray<UArray<UArray<UArray<std::complex<double>>>>> &pots) {
	pots.SetCount(dt.Nb);
	for (int ib = 0; ib < dt.Nb; ++ib) {
		if (pots[ib].IsEmpty())
			pots[ib].SetCount(dt.msh[ib].dt.mesh.GetNumPanels());
		for (int ipot = 0; ipot < pots[ib].size(); ++ipot) {
			pots[ib][ipot].SetCount(dt.Nh);
			for (int ih = 0; ih < dt.Nh; ++ih)
				pots[ib][ipot][ih].SetCount(dt.Nf, 0);
		}
	}
}

void Hydro::Initialize_Forces() {
	Initialize_Forces(dt.ex);
	Initialize_Forces(dt.sc);
	Initialize_Forces(dt.fk);
}

void Hydro::Initialize_Forces(Forces &f, int _Nh, double val) {
	if (_Nh == -1)
		_Nh = dt.Nh;
	f.force.SetCount(_Nh);
	for (int ih = 0; ih < _Nh; ++ih) 
		f.force[ih].setConstant(dt.Nf, dt.Nb*6, val);
}
	
void Hydro::Normalize_Forces(Forces &f) {
	for (int ih = 0; ih < dt.Nh; ++ih) 
		for (int ifr = 0; ifr < dt.Nf; ++ifr) 
			for (int idf = 0; idf < 6*dt.Nb; ++idf) 
				f.force[ih](ifr, idf) = F_dim(f, ih, ifr, idf);
}

void Hydro::Normalize_RAO(RAO &f) {
	for (int ih = 0; ih < dt.Nh; ++ih) 
		for (int ifr = 0; ifr < dt.Nf; ++ifr) 
			for (int idf = 0; idf < 6*dt.Nb; ++idf) 
				f.force[ih](ifr, idf) = RAO_dim(f, ih, ifr, idf);
}

void Hydro::Dimensionalize_Forces(Forces &f) {
	for (int ih = 0; ih < dt.Nh; ++ih) 
		for (int ifr = 0; ifr < dt.Nf; ++ifr) 
			for (int idf = 0; idf < 6*dt.Nb; ++idf) 
				f.force[ih](ifr, idf) = F_dim(f, ih, ifr, idf);
}

void Hydro::Dimensionalize_RAO(RAO &f) {
	for (int ih = 0; ih < dt.Nh; ++ih) 
		for (int ifr = 0; ifr < dt.Nf; ++ifr) 
			for (int idf = 0; idf < 6*dt.Nb; ++idf) 
				f.force[ih](ifr, idf) = RAO_dim(f, ih, ifr, idf);
}

void Hydro::Add_Forces(Forces &to, const Hydro &hy, const Forces &from) {
	if (hy.IsLoadedForce(from)) {
		for (int ihhy = 0; ihhy < hy.dt.Nh; ++ihhy) {
			int ih = FindClosest(dt.head, hy.dt.head[ihhy]);
			for (int ifrhy = 0; ifrhy < hy.dt.Nf; ++ifrhy) {
				int ifr = FindClosest(dt.w, hy.dt.w[ifrhy]);
				for (int idf = 0; idf < 6*dt.Nb; ++idf) 	 
					if (IsNum(from.force[ihhy](ifrhy, idf))) 
						to.force[ih](ifr, idf) = hy.F_ndim(from, ihhy, ifrhy, idf);
			}
		} 
	}
}

void Hydro::Add_RAO(RAO &to, const Hydro &hy, const RAO &from) {
	if (hy.IsLoadedForce(from)) {
		for (int ihhy = 0; ihhy < hy.dt.Nh; ++ihhy) {
			int ih = FindClosest(dt.head, hy.dt.head[ihhy]);
			for (int ifrhy = 0; ifrhy < hy.dt.Nf; ++ifrhy) {
				int ifr = FindClosest(dt.w, hy.dt.w[ifrhy]);
				for (int idf = 0; idf < 6*dt.Nb; ++idf) 	 
					if (IsNum(from.force[ihhy](ifrhy, idf))) 
						to.force[ih](ifr, idf) = hy.RAO_ndim(from, ihhy, ifrhy, idf);
			}
		} 
	}
}

static double MirrorHead(double head, bool xAxis) {
	if (xAxis)
		return -head;
	else {
		head = FixHeading_180(head);
		head += 90;
		head = FixHeading_180(-head);
		head -= 90;
		return head;		
	}
}

static std::complex<double> MirrorHead(const std::complex<double> &head, bool xAxis) {
	double h1 = head.real(), h2 = head.imag();
	if (xAxis) {
		h1 = -h1;
		h2 = -h2;	
	} else {
		auto Fix = [](double h) {
			h = FixHeading_180(h);
			h += 90;
			h = FixHeading_180(-h);
			h -= 90;
		};
		Fix(h1);
		Fix(h2);
	}
	return std::complex<double>(h1, h2);
}

void Hydro::Copy(const Hydro &hyd) {
	dt.file = hyd.dt.file;
	dt.name = hyd.dt.name;
	dt.g = hyd.dt.g;
    dt.h = hyd.dt.h;
    dt.rho = hyd.dt.rho;
    dt.len = hyd.dt.len;
    dt.dimen = hyd.dt.dimen;
    dt.Nb = hyd.dt.Nb;
    dt.Nf = hyd.dt.Nf;
    dt.Nh = hyd.dt.Nh;

    dt.A = clone(hyd.dt.A);
	dt.Ainf_w = clone(hyd.dt.Ainf_w);
    dt.Ainf = clone(hyd.dt.Ainf);
    dt.A0 = clone(hyd.dt.A0);
	
	//Dlin = clone(Dlin);
	//Dquad = clone(Dquad);
	//Cmoor = clone(Cmoor);
	
    dt.B = clone(hyd.dt.B);
    
    dt.msh = clone(hyd.dt.msh);
    
    dt.head = clone(hyd.dt.head);
    //names = clone(hyd.names);
    //C = clone(hyd.C);
    //M = clone(hyd.M);
    //cb = clone(hyd.cb);
    //cg = clone(hyd.cg);
    //c0 = clone(hyd.c0);
    dt.solver = hyd.dt.solver;     
    //dof = clone(hyd.dof); 
    
    dt.Kirf = clone(hyd.dt.Kirf);
    dt.Tirf = clone(hyd.dt.Tirf);
    
    dt.ex = clone(hyd.dt.ex);
    dt.sc = clone(hyd.dt.sc);
    dt.fk = clone(hyd.dt.fk);
    dt.rao = clone(hyd.dt.rao);
    
    dt.description = hyd.dt.description;

    dt.sts = clone(hyd.dt.sts);
    dt.dimenSTS = hyd.dt.dimenSTS;
    dt.stsProcessor = hyd.dt.stsProcessor;
    
    dt.qtfsum = clone(hyd.dt.qtfsum);
    dt.qtfdif = clone(hyd.dt.qtfdif);
    dt.qw = clone(hyd.dt.qw);
    dt.qh = clone(hyd.dt.qh);
    dt.qtfdataFromW = hyd.dt.qtfdataFromW;
    dt.qtftype = hyd.dt.qtftype;
    
    dt.mdhead = clone(hyd.dt.mdhead);
	dt.md = clone(hyd.dt.md);
	dt.mdtype = hyd.dt.mdtype;
	    
    //T = clone(hyd.T);
    dt.w = clone(hyd.dt.w);
    dt.dataFromW = hyd.dt.dataFromW;
    //Vo = clone(hyd.Vo); 
    
    dt.SetId(hyd.dt.GetId());
}

void AvgB(Eigen::MatrixXd &ret, const UArray<const Eigen::MatrixXd*> &d) {
	int numT = d.size();
	if (numT == 0) 
		return;
	
	int num = int(d[0]->size());
	for (int it = 1; it < numT; ++it) 
		if (d[it]->size() != num)
			throw Exc(t_("Avg() has to have same number of values"));
	
	for (int i = 0; i < num; ++i) {
		Eigen::VectorXd r(numT);
		for (int it = 0; it < numT; ++it) 
			r[it] = (*d[it])(i);
		ret(i) = r.mean();
	}
}

using UA = UArray<UArray<VectorXd>>;
void AvgB(UA &ret, const UArray<const UA*> &d) {
	int numT = d.size();
	if (numT == 0) 
		return;

	for (int i = 0; i < ret.size(); ++i) {
		for (int j = 0; j < ret[i].size(); ++j) {
			for (int k = 0; k < ret[i][j].size(); ++k) {
				Eigen::VectorXd r(numT);
				for (int it = 0; it < numT; ++it) 
					r[it] = (*d[it])[i][j][k];
				ret[i][j][k] = r.mean();
			}
		}
	}
}

using UAM = UArray<MatrixXd>;
void AvgB(UAM &ret, const UArray<const UAM*> &d) {
	int numT = d.size();
	if (numT == 0) 
		return;

	for (int i = 0; i < ret.size(); ++i) {
		for (int j = 0; j < ret[i].size(); ++j) {
			Eigen::VectorXd r(numT);
			for (int it = 0; it < numT; ++it) 
				r[it] = (*d[it])[i](j);
			ret[i](j) = r.mean();
		}
	}
}

void AvgB(Hydro::Forces &ret, const UArray<const Hydro::Forces*> &d) {
	int numT = d.size();
	if (numT == 0) 
		return;

	for (int i = 0; i < ret.force.size(); ++i) {
		for (int j = 0; j < ret.force[i].cols(); ++j) {
			for (int k = 0; k < ret.force[i].rows(); ++k) {
				Eigen::VectorXcd r(numT);
				for (int it = 0; it < numT; ++it) 
					r[it] = (*d[it]).force[i](k, j);
				ret.force[i](k, j) = r.mean();
			}
		}
	}
}

using UMD = UArray<UArray<UArray<VectorXd>>>;
void AvgB(UMD &ret, const UArray<const UMD*> &d) {
	int numT = d.size();
	if (numT == 0) 
		return;

	for (int i = 0; i < ret.size(); ++i) {
		for (int j = 0; j < ret[i].size(); ++j) {
			for (int k = 0; k < ret[i][j].size(); ++k) {
				for (int l = 0; l < ret[i][j][k].size(); ++l) {
					Eigen::VectorXd r(numT);
					for (int it = 0; it < numT; ++it) 
						r[it] = (*d[it])[i][j][k][l];
					ret[i][j][k][l] = r.mean();
				}
			}
		}
	}
}

using UMQ = UArray<UArray<UArray<MatrixXcd>>>;
void AvgB(UMQ &ret, const UArray<const UMQ*> &d) {
	int numT = d.size();
	if (numT == 0) 
		return;

	for (int i = 0; i < ret.size(); ++i) {
		for (int j = 0; j < ret[i].size(); ++j) {
			for (int k = 0; k < ret[i][j].size(); ++k) {
				for (int l = 0; l < ret[i][j][k].rows(); ++l) {
					for (int m = 0; l < ret[i][j][k].cols(); ++m) {
						Eigen::VectorXcd r(numT);
						for (int it = 0; it < numT; ++it) 
							r[it] = (*d[it])[i][j][k](l, m);
						ret[i][j][k](l, m) = r.mean();
					}
				}
			}
		}
	}
}

void Hydro::Average(const UArray<Hydro> &hydros, const UVector<int> &ids) {
	const Hydro &h0 = hydros[0];
	dt.Nb = h0.dt.Nb;
    dt.Nf = h0.dt.Nf;
    dt.Nh = h0.dt.Nh;
	dt.h = h0.dt.h;
	dt.rho = h0.dt.rho;
	dt.g = h0.dt.g;

	//T = clone(h0.T);
    dt.w = clone(h0.dt.w);
    dt.head = clone(h0.dt.head);
	dt.qw = clone(h0.dt.qw);
    dt.qh = clone(h0.dt.qh);
	dt.mdhead = clone(h0.dt.mdhead);
	dt.msh = clone(h0.dt.msh);
	//c0 = clone(h0.c0);
	//dof = clone(h0.dof);
	
	for (int i = 1; i < ids.size(); ++i) {
		const Hydro &hy = hydros[i];
		if (dt.Nb != hy.dt.Nb)
			throw Exc(t_("All models have to have the same number of bodies"));
		if (dt.Nf != hy.dt.Nf)
			throw Exc(t_("All models have to have the same number of frequencies"));
		if (dt.Nh != hy.dt.Nh)
			throw Exc(t_("All models have to have the same number of headings"));
		if (dt.rho != hy.dt.rho)
			throw Exc(t_("All models have to have the same density"));
		if (dt.g != hy.dt.g)
			throw Exc(t_("All models have to have the same gravity"));
		if (dt.h != hy.dt.h) 
			dt.h = -1;
		//if (!CompareDecimals(T, hy.T, 2))
		//	throw Exc(t_("All models have to have the same periods"));
		if (!CompareDecimals(dt.w, hy.dt.w, 2))
			throw Exc(t_("All models have to have the same frequencies"));
		if (!CompareDecimals(dt.head, hy.dt.head, 2))
			throw Exc(t_("All models have to have the same headings"));
		if (!CompareDecimals(dt.qw, hy.dt.qw, 2))
			throw Exc(t_("All models have to have the same frequencies in QTF"));
		if (!Compare(dt.qh, hy.dt.qh))
			throw Exc(t_("All models have to have the same headings in QTF"));
		if (!Compare(dt.mdhead, hy.dt.mdhead))
			throw Exc(t_("All models have to have the same headings in mean drift"));
		for (int ib = 0; ib < dt.Nb; ++ib) {
			if (dt.msh[ib].dt.c0 != hy.dt.msh[ib].dt.c0)
				throw Exc(t_("All models and bodies have to have the same centre of reference"));
		}
		//if (c0 != hy.c0)
		//	throw Exc(t_("All models have to have the same centre of reference"));
		//if (dof != hy.dof)
			//throw Exc(t_("All models have to have the same number of dof"));
	}
	
	dt.file = "Average";
	dt.name = "Average";
	dt.description = "Average";
	dt.solver = BEMROSETTA;
	
    dt.len = 1;
    dt.dimen = true;
    
    dt.dataFromW = dt.qtfdataFromW = true;
    
    //names = clone(h0.names);
    
    dt.qtftype = h0.dt.qtftype;
    dt.mdtype = h0.dt.mdtype;
    
    //Vo.SetCount(Nb, NaNDouble);
	//cg.setConstant(3, Nb, NaNDouble);
	//cb.setConstant(3, Nb, NaNDouble);
	dt.Ainf.setConstant(dt.Nb*6, dt.Nb*6, NaNDouble);
	dt.A0.setConstant(dt.Nb*6, dt.Nb*6, NaNDouble);
	
	Initialize_AB(dt.A);
	Initialize_AB(dt.B);
	
	//C.SetCount(Nb);
	for (int ib = 0; ib < dt.Nb; ++ib) 
		dt.msh[ib].dt.C.setConstant(6, 6, 0);
	//M.SetCount(Nb);
	for (int ib = 0; ib < dt.Nb; ++ib) 
		dt.msh[ib].dt.M.setConstant(6, 6, 0);
	
	Initialize_Forces(dt.ex);
	Initialize_Forces(dt.sc);
	Initialize_Forces(dt.fk);
	Initialize_Forces(dt.rao);
	
	Hydro::Initialize_MD(dt.md, dt.Nb, int(dt.mdhead.size()), dt.Nf);
		
	Hydro::Initialize_QTF(dt.qtfsum, dt.Nb, int(dt.qh.size()), int(dt.qw.size()));
	Hydro::Initialize_QTF(dt.qtfdif, dt.Nb, int(dt.qh.size()), int(dt.qw.size()));
			
	//UArray<const UVector<double>*> Vos;
	UArray<const MatrixXd*> Ainfs, A0s;
	UArray<const UArray<UArray<VectorXd>>*> As, Bs;
	UArray<const Forces*> exs, scs, fks, raos;
	UArray<const UArray<UArray<UArray<VectorXd>>>*> mds;
	UArray<const UArray<UArray<UArray<MatrixXcd>>>*> qtfsums, qtfdifs;
	UArray<UArray<const MatrixXd*>> Cs(dt.Nb), Ms(dt.Nb);
	UArray<UArray<const Point3D*>> cgs(dt.Nb), cbs(dt.Nb);
	VectorXd Vos(dt.Nb);
	
	for (int i = 0; i < ids.size(); ++i) {
		const Hydro &hy = hydros[i];
		
		for (int ib = 0; ib < dt.Nb; ++ib) {
			if (!IsNull(hy.dt.msh[ib].dt.cg))
				cgs[ib].Add(&(hy.dt.msh[ib].dt.cg));
			if (!IsNull(hy.dt.msh[ib].dt.cb))
				cbs[ib].Add(&(hy.dt.msh[ib].dt.cb));
			if (!IsNull(hy.dt.msh[ib].dt.Vo))
				Vos[ib] = hy.dt.msh[ib].dt.Vo;
			Cs[ib].Add(&(hy.dt.msh[ib].dt.C));
			Ms[ib].Add(&(hy.dt.msh[ib].dt.M));
		}
		//if (hy.Vo.size() > i)
        	//Vos.Add(&(hy.Vo));
		//if (hy.cg.size() > 0)
        //	cgs.Add(&(hydros[i].cg));
		//if (hy.cb.size() > 0)
        //	cbs.Add(&(hydros[i].cb));
		if (hy.dt.Ainf.size() > 0)
        	Ainfs.Add(&hy.dt.Ainf);
		if (hy.dt.A0.size() > 0)
        	A0s.Add(&hy.dt.A0);
        As.Add(&hy.dt.A);
        Bs.Add(&hy.dt.B);
        //Cs.Add(&hy.C);
        //Ms.Add(&hy.M);
        exs.Add(&hy.dt.ex);
        scs.Add(&hy.dt.sc);
        fks.Add(&hy.dt.fk);
        raos.Add(&hy.dt.rao);
        mds.Add(&hy.dt.md);
        qtfsums.Add(&hy.dt.qtfsum);
        qtfdifs.Add(&hy.dt.qtfdif);
    }
    
    for (int ib = 0; ib < dt.Nb; ++ib) {
    	AvgB(dt.msh[ib].dt.cg, cgs[ib]);
    	AvgB(dt.msh[ib].dt.cb, cbs[ib]);
    	AvgB(dt.msh[ib].dt.C, Cs[ib]);
    	AvgB(dt.msh[ib].dt.M, Ms[ib]);
    	
    	dt.msh[ib].dt.Vo = Vos.mean();
    }
    
	//AvgB(Vo, Vos);   
	
	AvgB(dt.Ainf, Ainfs); 
	AvgB(dt.A0, A0s); 
	
	AvgB(dt.A, As); 
	AvgB(dt.B, Bs);
	 
	//AvgB(C, Cs);
	//AvgB(M, Ms);
	
    AvgB(dt.ex, exs);
    AvgB(dt.sc, scs);
    AvgB(dt.fk, fks);
    AvgB(dt.rao, raos);
    
    AvgB(dt.md, mds);
    
	AvgB(dt.qtfsum, qtfsums);
    AvgB(dt.qtfdif, qtfdifs);

	
	/*
	Ainf_w
	
	Dlin = clone(Dlin);
	Cmoor
   
    Kirf = clone(hyd.Kirf);
    Tirf = clone(hyd.Tirf);
    */
}

bool Hydro::SymmetryRule(int idf6, bool xAxis) {
	const bool symmetryRulesXZ[] = {false, true,  false, true,  false, true};
	const bool symmetryRulesYZ[] = {true,  false, false, false, true,  true};
	
	ASSERT(idf6 >= 0);
	idf6 = idf6%6;
	if (xAxis)
		return symmetryRulesXZ[idf6];
	else
		return symmetryRulesYZ[idf6];
}
	
void Hydro::Symmetrize_Forces(bool xAxis) {
	if (!IsLoadedFex() && !IsLoadedFsc() && !IsLoadedFfk() && !IsLoadedRAO())
		return;
	
	UVector<double> newHead;
	for (int ih = 0; ih < dt.Nh; ++ih) {
		FindAddDelta(newHead, FixHeading_180(dt.head[ih]), 0.01);
		FindAddDelta(newHead, FixHeading_180(MirrorHead(dt.head[ih], xAxis)), 0.01);
	}
	Sort(newHead);
	
	auto Symmetrize_ForcesEach = [&](const Forces &f, Forces &newf) {
		auto Symmetrize_Forces_Each0 = [&](const Forces &f, Forces &newf, double hh, int ih, int idf, bool applysym) {
			int nih = FindClosest(newHead, hh);
			bool avg = IsNum(newf.force[nih](0, idf));
			for (int ifr = 0; ifr < dt.Nf; ++ifr) {
				std::complex<double> force = f.force[ih](ifr, idf);
				if (applysym)
					force = std::polar(abs(force), arg(force) + M_PI);
				std::complex<double> &nf = newf.force[nih](ifr, idf);
				if (avg)
					nf = Avg(nf, force);
				else 
					nf = force;
			}
		};
	
		Initialize_Forces(newf, newHead.size());
		
		for (int idf = 0; idf < 6*dt.Nb; ++idf) {
			for (int ih = 0; ih < dt.Nh; ++ih) {
				Symmetrize_Forces_Each0(f, newf, FixHeading_180(dt.head[ih]), ih, idf, false);
				double he = FixHeading_180(MirrorHead(dt.head[ih], xAxis));
				bool applysym = false;
				if (SymmetryRule(idf, xAxis)) {
					if (xAxis) {
						if (he != 0 && he != 180)
							applysym = true;
					} else {
						if (he != 90 && he != -90)
							applysym = true;
					}
				}
				Symmetrize_Forces_Each0(f, newf, he, ih, idf, applysym);
			}
		}
	};

	if (IsLoadedFex()) {
		Forces newex;
		Symmetrize_ForcesEach(dt.ex, newex);
		dt.ex = pick(newex);
	}
	if (IsLoadedFsc()) {
		Forces newsc;
		Symmetrize_ForcesEach(dt.sc, newsc);
		dt.sc = pick(newsc);
	}
	if (IsLoadedFfk()) {
		Forces newfk;
		Symmetrize_ForcesEach(dt.fk, newfk);
		dt.fk = pick(newfk);
	}
	if (IsLoadedRAO()) {
		RAO newrao;
		Symmetrize_ForcesEach(dt.rao, newrao);
		dt.rao = pick(newrao);
	}
	dt.Nh = newHead.size();
	dt.head = pick(newHead);		// New headings are set between -180 and 180
}

void Hydro::Symmetrize_QTF(bool xAxis) {
	if (!IsLoadedQTF(true) && !IsLoadedQTF(false))
		return;
	
	UArray<std::complex<double>> newHead;
	for (int ih = 0; ih < dt.qh.size(); ++ih) {
		FindAddDelta(newHead, FixHeading_180(dt.qh[ih]), 0.01);
		FindAddDelta(newHead, FixHeading_180(MirrorHead(dt.qh[ih], xAxis)), 0.01);
	}
	Sort(newHead, SortComplex);

	auto Symmetrize_ForcesEach = [&](const UArray<UArray<UArray<MatrixXcd>>> &f, UArray<UArray<UArray<MatrixXcd>>> &newf) {
		auto Symmetrize_Forces_Each0 = [&](const UArray<UArray<UArray<MatrixXcd>>> &f, UArray<UArray<UArray<MatrixXcd>>> &newf, const std::complex<double> &hh, int ib, int ih, int idf, bool applysym) {
			int nih = FindClosest(newHead, hh);
			bool avg = IsNum(newf[ib][nih][idf](0, 0));
			for (int ifr1 = 0; ifr1 < dt.qw.size(); ++ifr1) {
				for (int ifr2 = 0; ifr2 < dt.qw.size(); ++ifr2) {
					std::complex<double> force = f[ib][ih][idf](ifr1, ifr2);
					if (applysym) 
						force = std::polar(abs(force), arg(force) + M_PI);
					std::complex<double> &nf = newf[ib][nih][idf](ifr1, ifr2);
					if (avg) 
						nf = Avg(nf, force);
					else 
						nf = force;
				}
			}
		};
	
		Initialize_QTF(newf, dt.Nb, newHead.size(), int(dt.qw.size()));
		
		for (int ib = 0; ib < dt.Nb; ++ib) {
			for (int ih = 0; ih < dt.qh.size(); ++ih) {
				for (int idf = 0; idf < 6; ++idf) {
					Symmetrize_Forces_Each0(f, newf, FixHeading_180(dt.qh[ih]), ib, ih, idf, false);
					std::complex<double> he = FixHeading_180(MirrorHead(dt.qh[ih], xAxis));
					bool applysym = false;
					if (SymmetryRule(idf, xAxis)) {
						if (xAxis) {
							if (he != std::complex<double>(0, 0) && he != std::complex<double>(180, 180))
								applysym = true;
						} else {
							if (he != std::complex<double>(90, 90) && he != std::complex<double>(-90, -90))
								applysym = true;
						}
					}
					Symmetrize_Forces_Each0(f, newf, he, ib, ih, idf, applysym);
				}
			}
		}
	};	
	
	if (IsLoadedQTF(true)) {
		UArray<UArray<UArray<MatrixXcd>>> newqtf;
		Symmetrize_ForcesEach(dt.qtfsum, newqtf);
		dt.qtfsum = pick(newqtf);
	}
	if (IsLoadedQTF(false)) {
		UArray<UArray<UArray<MatrixXcd>>> newqtf;
		Symmetrize_ForcesEach(dt.qtfdif, newqtf);
		dt.qtfdif = pick(newqtf);
	}
	
	::Copy(newHead, dt.qh);	
}

void Hydro::Symmetrize_MD(bool xAxis) {
	if (!IsLoadedMD())
		return;
	
	UArray<std::complex<double>> newHead;
	for (int ih = 0; ih < dt.mdhead.size(); ++ih) {
		FindAddDelta(newHead, FixHeading_180(dt.mdhead[ih]), 0.01);
		FindAddDelta(newHead, FixHeading_180(MirrorHead(dt.mdhead[ih], xAxis)), 0.01);
	}
	Sort(newHead, SortComplex);
	
	auto Symmetrize_ForcesEach = [&](const UArray<UArray<UArray<VectorXd>>> &f, UArray<UArray<UArray<VectorXd>>> &newf) {
		auto Symmetrize_Forces_Each0 = [&](const UArray<UArray<UArray<VectorXd>>> &f, UArray<UArray<UArray<VectorXd>>> &newf, const std::complex<double> &hh, int ib, int ih, int idf, bool applysym) {
			int nih = FindClosest(newHead, hh);
			bool avg = IsNum(newf[ib][nih][idf](0));
			for (int ifr = 0; ifr < dt.w.size(); ++ifr) {
				double force = f[ib][ih][idf](ifr);
				if (applysym) 
					force = -force;
				double& nf = newf[ib][nih][idf](ifr);
				if (avg) 
					nf = Avg(nf, force);
				else 
					nf = force;
			}
		};
	
		Initialize_MD(newf, dt.Nb, newHead.size(), dt.w.size());
		
		for (int ib = 0; ib < dt.Nb; ++ib) {
			for (int ih = 0; ih < dt.mdhead.size(); ++ih) {
				for (int idf = 0; idf < 6; ++idf) {
					Symmetrize_Forces_Each0(f, newf, FixHeading_180(dt.mdhead[ih]), ib, ih, idf, false);
					std::complex<double> he = FixHeading_180(MirrorHead(dt.mdhead[ih], xAxis));
					bool applysym = false;
					if (SymmetryRule(idf, xAxis)) {
						if (xAxis) {
							if (he != std::complex<double>(0, 0) && he != std::complex<double>(180, 180))
								applysym = true;
						} else {
							if (he != std::complex<double>(90, 90) && he != std::complex<double>(-90, -90))
								applysym = true;
						}
					}
					Symmetrize_Forces_Each0(f, newf, he, ib, ih, idf, applysym);
				}
			}
		}
	};	
	
	UArray<UArray<UArray<VectorXd>>> newmd;
	Symmetrize_ForcesEach(dt.md, newmd);
	dt.md = pick(newmd);
	
	::Copy(newHead, dt.mdhead);	
}

void Hydro::RemoveThresDOF_A(double thres) {
	if (!IsLoadedA())
		return;
	for (int idf = 0; idf < 6*dt.Nb; ++idf) {
		for (int jdf = 0; jdf < 6*dt.Nb; ++jdf) {
			double mx = -DBL_MAX, mn = DBL_MAX;
			for (int ifr = 0; ifr < dt.Nf; ifr++) {
				double val = A_ndim(ifr, idf, jdf);
				mx = max(mx, val);
				mn = min(mn, val);
			}
			double delta = mx - mn;
			if (IsNum(mx) && IsNum(mn)) {
				double res = 0;
				for (int ifr = 1; ifr < dt.Nf; ifr++) 
					res += abs(A_ndim(ifr, idf, jdf) - A_ndim(ifr-1, idf, jdf));
				res /= delta*(dt.Nf - 1);
				if (res > thres) {
					for (int ifr = 0; ifr < dt.Nf; ifr++) 
						dt.A[idf][jdf][ifr] = Null;
					dt.Ainf(idf, jdf) = Null;		
				}
			}
		}
	}
}

void Hydro::RemoveThresDOF_B(double thres) {
	if (!IsLoadedB())
		return;
	for (int idf = 0; idf < 6*dt.Nb; ++idf) {
		for (int jdf = 0; jdf < 6*dt.Nb; ++jdf) {
			double mx = -DBL_MAX, mn = DBL_MAX;
			for (int ifr = 0; ifr < dt.Nf; ifr++) {
				double val = B_ndim(ifr, idf, jdf);
				mx = max(mx, val);
				mn = min(mn, val);
			}
			double delta = mx - mn;
			if (IsNum(mx) && IsNum(mn)) {
				double res = 0;
				for (int ifr = 1; ifr < dt.Nf; ifr++) 
					res += abs(B_ndim(ifr, idf, jdf) - B_ndim(ifr-1, idf, jdf));
				res /= delta*(dt.Nf - 1);
				if (res > thres) {
					for (int ifr = 0; ifr < dt.Nf; ifr++) 
						dt.B[idf][jdf][ifr] = Null;
				}
			}
		}
	}
}

void Hydro::RemoveThresDOF_Force(Forces &f, double thres) {
	if (!IsLoadedForce(f))
		return;
	for (int ih = 0; ih < dt.Nh; ++ih) {
		for (int i = 0; i < 6*dt.Nb; ++i) {
			double mx = -DBL_MAX, mn = DBL_MAX;
			for (int ifr = 0; ifr < dt.Nf; ifr++) {
				double val = abs(F_ndim(f, ih, ifr, i));
				mx = max(mx, val);
				mn = min(mn, val);
			}
			if (mx != -DBL_MAX && mn != DBL_MAX) {
				double delta = mx - mn;
				double res = 0;
				for (int ifr = 1; ifr < dt.Nf; ifr++) 
					res += abs(F_ndim(f, ih, ifr, i) - F_ndim(f, ih, ifr-1, i));
				res /= delta*(dt.Nf - 1);
				if (res > thres) {
					for (int ifr = 0; ifr < dt.Nf; ifr++) 
						f.force[ih](ifr, i) = Null;
				}
			}
		}
	}
}

void Hydro::Compare_rho(Hydro &a) {
	if (a.dt.rho != dt.rho)
		throw Exc(Format(t_("%s is not the same %f<>%f"), t_("Density rho"), a.dt.rho, dt.rho));
}

void Hydro::Compare_g(Hydro &a) {
	if (a.dt.g != dt.g)
		throw Exc(Format(t_("%s is not the same %f<>%f"), t_("Gravity g"), a.dt.g, dt.g));
}

void Hydro::Compare_h(Hydro &a) {
	if (a.dt.h != dt.h)
		throw Exc(Format(t_("%s is not the same %f<>%f"), t_("Water depth h"), a.dt.h, dt.h));
}

void Hydro::Compare_Nb(Hydro &a) {
	if (a.dt.Nb != dt.Nb)
		throw Exc(Format(t_("%s is not the same %d<>%d"), t_("Number of bodies"), a.dt.Nb, dt.Nb));
}

void Hydro::Compare_w(Hydro &a) {
	if (a.dt.Nf != dt.Nf)	
		throw Exc(Format(t_("%s is not the same %f<>%f"), t_("Number of frequencies"), a.dt.Nf, dt.Nf));
	for (int i = 0; i < a.dt.Nf; ++i) {
		if (!EqualRatio(a.dt.w[i], dt.w[i], 0.0001))
			throw Exc(Format(t_("%s is not the same %f<>%f"), 
							Format(t_("#%d %s"), i+1, t_("frequency")), a.dt.w[i], dt.w[i]));
	}
}

void Hydro::Compare_head(Hydro &a) {
	if (a.dt.Nh != dt.Nh)	
		throw Exc(Format(t_("%s is not the same %f<>%f"), t_("Number of headings"), a.dt.Nh, dt.Nh));
	for (int i = 0; i < a.dt.Nh; ++i) {
		if (a.dt.head[i] != dt.head[i])
			throw Exc(Format(t_("%s is not the same %f<>%f"), 
							Format(t_("#%d %s"), i+1, t_("frequency")), a.dt.w[i], dt.w[i]));
	}
}

void Hydro::Compare_A(const UArray<UArray<VectorXd>> &a) {
	for (int ifr = 0; ifr < dt.Nf; ifr++) {
		for (int idf = 0; idf < 6*dt.Nb; ++idf) {
			for (int jdf = 0; jdf < 6*dt.Nb; ++jdf) {
				double Aa = a[idf][jdf][ifr];
				double Ab = dt.A[idf][jdf][ifr];
				if (IsNum(Aa) && IsNum(Ab) && Aa != Ab)
					throw Exc(Format(t_("%s is not the same %f<>%f"), 
							Format(t_("%s[%d](%d, %d)"), t_("A"), ifr+1, idf+1, jdf+1), 
							Aa, Ab));
			}
		}
	}
}

void Hydro::Compare_B(const UArray<UArray<VectorXd>> &b) {
	for (int ifr = 0; ifr < dt.Nf; ifr++) {
		for (int idf = 0; idf < 6*dt.Nb; ++idf) {
			for (int jdf = 0; jdf < 6*dt.Nb; ++jdf) {
				double Ba = b[idf][jdf][ifr];
				double Bb = dt.B[idf][jdf][ifr];
				if (IsNum(Ba) && IsNum(Bb) && Ba != Bb)
					throw Exc(Format(t_("%s is not the same %f<>%f"), 
							Format(t_("%s[%d](%d, %d)"), t_("B"), ifr+1, idf+1, jdf+1), 
							Ba, Bb));
			}
		}
	}
}

void Hydro::Compare_C(Hydro &a) {
	for (int ib = 0; ib < a.dt.Nb; ib++) {
		for (int idf = 0; idf < 6; ++idf) {
			for (int jdf = 0; jdf < 6; ++jdf) {
				double Ca = a.dt.msh[ib].dt.C(idf, jdf);
				double Cb = dt.msh[ib].dt.C(idf, jdf);
				if (IsNum(Ca) && IsNum(Cb) && !EqualRatio(Ca, Cb, 0.0001))
					throw Exc(Format(t_("%s is not the same %f<>%f"), 
							Format(t_("%s[%d](%d, %d)"), t_("C"), ib+1, idf+1, jdf+1), 
							Ca, Cb));
			}
		}
	}
}

void Hydro::Compare_cg(Hydro &a) {
	for (int i = 0; i < 3; i++) {
		for (int ib = 0; ib < a.dt.Nb; ib++) {
			if (a.dt.msh[ib].dt.cg[i] != dt.msh[ib].dt.cg[i])
				throw Exc(Format(t_("%s is not the same %f<>%f"), 
						Format(t_("%s(%d, %d)"), t_("cg"), i+1, ib+1), 
							a.dt.msh[ib].dt.cg[i], dt.msh[ib].dt.cg[i]));
		}
	}
}

void Hydro::Compare_F(const Forces &a, const Forces &b, String type) {
	for (int ih = 0; ih < dt.Nh; ++ih) 	
		for (int ifr = 0; ifr < dt.Nf; ++ifr)
			for (int idf = 0; idf < 6*dt.Nb; ++idf) {
				std::complex<double> Fa = a.force[ih](ifr, idf);
				std::complex<double> Fb = b.force[ih](ifr, idf);
				if (IsLoadedForce(a, idf, ih) && Fa != Fb)
					throw Exc(Format(t_("%s is not the same %f:%f<>%f:%f"), 
							Format(t_("%s[%d](%d, %d)"), type, ifr+1, idf+1, ih+1), 
							Fa.real(), Fa.imag(), Fb.real(), Fb.imag()));
			}
}

void Hydro::SaveAs(String fileName, Function <bool(String, int)> Status, BEM_FMT type, int qtfHeading) {
	int realNh = dt.Nh;
	int realNf = dt.Nf;
	
	if (type == UNKNOWN) {
		String ext = ToLower(GetFileExt(fileName));
		
		if (ext == ".1" || ext == ".2" || ext == ".3" || ext == ".3sc" || ext == ".3fk" || 
			ext == ".hst" || ext == ".4" || ext == ".12s" || ext == ".12d") 
			type = WAMIT_1_3;
		else if (ext == ".out")
			type = WAMIT;
		else if (ext == ".dat")
			type = FAST_WAMIT;	
		else if (ext == ".bem")
			type = BEMROSETTA;
		else if (ext == ".csv")
			type = CSV_TABLE;
		else if (ext == ".hdb")
			type = DIODORE;
		else if (ext == ".h5")
			type = BEMIOH5;
		else
			throw Exc(Format(t_("Conversion to file type '%s' not supported"), fileName));
	}
	if (type == WAMIT) {
		Wamit data;
		data.Save_out(fileName);			
	} else if (type == WAMIT_1_3) {
		Wamit data;
		data.Save(fileName, Status, true, qtfHeading);	
	} else if (type == FAST_WAMIT) {
		Fast data;
		data.Save(fileName, Status, qtfHeading);		
	} else if (type == BEMROSETTA) {
		Hydro data;
		data.SaveSerialization(fileName);		
	} else if (type == AQWA) {
		Aqwa data;
		data.Save(fileName, Status);		
	} else if (type == CSV_MAT) {
		Hydro data;
		data.SaveCSVMat(fileName);		
	} else if (type == CSV_TABLE) {
		Hydro data;
		data.SaveCSVTable(fileName);		
	} else if (type == DIODORE) {
		Hydro data;
		data.SaveDiodoreHDB(fileName);		
	} else if (type == BEMIOH5) {
		BemioH5 data;
		data.Save(fileName);
	} else
		throw Exc(Format(t_("Conversion to file type '%s' not supported"), fileName));
	
	//solver = type;
	dt.Nh = realNh;
	dt.Nf = realNf;
}

void Hydro::Join(const UVector<Hydro *> &hydrosp) {
	dt.name = t_("Joined files");
	for (int ihy = 0; ihy < hydrosp.size(); ++ihy) {
		const Hydro &hy = *hydrosp[ihy];
		if (hy.dt.name.Find("Nemoh_Part") >= 0) {
			dt.name = GetFileTitle(GetFileFolder(GetFileFolder(hy.dt.file)));
			break;
		}
	}
	
	dt.dimen = false;
	dt.g = Bem().g;
	dt.rho = Bem().rho;
	dt.len = 1;

	dt.h = Null;
	for (int ihy = 0; ihy < hydrosp.size(); ++ihy) {
		const Hydro &hy = *hydrosp[ihy];
		if (IsNum(hy.dt.h))
			dt.h = hy.dt.h;
		else if (dt.h != hy.dt.h)
			throw Exc(Format(t_("Water depth does not match between '%s'(%d) and '%s'(%d)"), 
					hydrosp[0]->dt.name, hydrosp[0]->dt.h, hy.dt.name, hy.dt.h));			
	}	
	if (!IsNum(dt.h))
		throw Exc(t_("No water depth found in models"));
			
	dt.Nb = Null;
	for (int ihy = 0; ihy < hydrosp.size(); ++ihy) {
		const Hydro &hy = *hydrosp[ihy];
		if (!IsNum(dt.Nb))
			dt.Nb = hy.dt.Nb;
		else if (dt.Nb != hy.dt.Nb)
			throw Exc(Format(t_("Number of bodies does not match between '%s'(%d) and '%s'(%d)"), 
					hydrosp[0]->dt.name, hydrosp[0]->dt.Nb, hy.dt.name, hy.dt.Nb));			
	}
	if (!IsNum(dt.Nb))
		throw Exc(t_("No body found in models"));
		
	dt.head.Clear();
	for (int ihy = 0; ihy < hydrosp.size(); ++ihy) {
		const Hydro &hy = *hydrosp[ihy];
		if (hy.IsLoadedFex() || hy.IsLoadedFsc() || hy.IsLoadedFfk()) {
			for (int ih = 0; ih < hy.dt.head.size(); ih++) {
				double head_v = hy.dt.head[ih];
				FindAddRatio(dt.head, head_v, 0.001);
			}
		}
	}
	dt.Nh = dt.head.size();
	if (dt.Nh == 0)
		throw Exc(t_("No head found in models"));
	
	dt.w.Clear();
	for (int ihy = 0; ihy < hydrosp.size(); ++ihy) {
		const Hydro &hy = *hydrosp[ihy];
		if (hy.IsLoadedA() && hy.IsLoadedB()) {
			for (int ifr = 0; ifr < hy.dt.w.size(); ifr++) {
				double w_v = hy.dt.w[ifr];
				FindAddRatio(dt.w, w_v, 0.001);
			}
		}
	}
	Sort(dt.w);
	dt.Nf = dt.w.size();
	/*T.Clear();
	for (int i = 0; i < Nf; ++i)
		T << 2*M_PI/w[i];*/
	
	if (dt.Nf == 0)
		throw Exc(t_("No frequency found in models"));
	
	dt.msh.SetCount(dt.Nb);
//	dof.SetCount(Nb);
//	cg.setConstant(3, Nb, NaNDouble);
//	c0.setConstant(3, Nb, NaNDouble);
//	cb.setConstant(3, Nb, NaNDouble);
//	Vo.SetCount(Nb, NaNDouble);
	
	Initialize_AB(dt.A);
	Initialize_AB(dt.B);
	
	//C.SetCount(Nb);
	for (int ib = 0; ib < dt.Nb; ++ib) 
		dt.msh[ib].dt.C.setConstant(6, 6, NaNDouble);
	
	Initialize_Forces(dt.ex);
	Initialize_Forces(dt.sc);
	Initialize_Forces(dt.fk);
		
	for (int ihy = 0; ihy < hydrosp.size(); ++ihy) {
		const Hydro &hy = *hydrosp[ihy];
		
		// All this block should have to be the same. Now it is not tested
		
		dt.solver = hy.dt.solver;
		
		for (int ib = 0; ib < dt.Nb; ++ib) {
			if (!hy.dt.msh[ib].dt.name.IsEmpty())
				dt.msh[ib].dt.name = hy.dt.msh[ib].dt.name;
			//if (IsNum(hy.Vo[ib]))
			dt.msh[ib].dt.Vo = hy.dt.msh[ib].dt.Vo;
			dt.msh[ib].dt.cg = clone(hy.dt.msh[ib].dt.cg);
			dt.msh[ib].dt.c0 = clone(hy.dt.msh[ib].dt.c0);
			for (int i = 0; i < 3; ++i) {
				/*if (IsNum(hy.cg(i, ib)))
					msh[ib].cg[i] = hy.cg(i, ib);*/
				if (IsNum(hy.dt.msh[ib].dt.cb[i]))
					dt.msh[ib].dt.cb[i] = hy.dt.msh[ib].dt.cb[i];
				/*if (IsNum(hy.c0(i, ib)))
					c0(i, ib) = hy.c0(i, ib);*/
			}
			//dof[ib] = hy.dof[ib];
		}
		
		if (/*IsLoadedC() && */ hy.IsLoadedC()) {
			for (int ib = 0; ib < dt.Nb; ++ib) {
				for (int idf = 0; idf < 6; ++idf) 
					for (int jdf = 0; jdf < 6; ++jdf) 
						dt.msh[ib].dt.C(idf, jdf) = hy.C_ndim(ib, idf, jdf);
			}
		}
		///////////////////////////////////////////////////////////////////
		
		if (/*IsLoadedA() && IsLoadedB() && */hy.IsLoadedA() && hy.IsLoadedB()) {
			for (int ifrhy = 0; ifrhy < hy.dt.Nf; ++ifrhy) {
				int ifr = FindClosest(dt.w, hy.dt.w[ifrhy]);
				for (int idf = 0; idf < 6*dt.Nb; ++idf) {
					for (int jdf = 0; jdf < 6*dt.Nb; ++jdf) {	
						if (IsNum(hy.dt.A[idf][jdf][ifrhy]))
							dt.A[idf][jdf][ifr] = hy.A_ndim(ifrhy, idf, jdf);
						if (IsNum(hy.dt.B[idf][jdf][ifrhy]))
							dt.B[idf][jdf][ifr] = hy.B_ndim(ifrhy, idf, jdf);
					}
				}
			}
		}	
		Add_Forces(dt.ex, hy, hy.dt.ex);
		Add_Forces(dt.sc, hy, hy.dt.sc);
		Add_Forces(dt.fk, hy, hy.dt.fk);
		
		Add_RAO(dt.rao, hy, hy.dt.rao);
	}

	
	// Aw0 has to be recalculated
	/*
	// A0 is set from the lower frequency data set
	
	int ihminw = -1;
	double minw = DBL_MAX;
	for (int ihy = 0; ihy < hydrosp.size(); ++ihy) {
		const Hydro &hy = *hydrosp[ihy];
		if (hydro.IsLoadedAw0()) {
			double minwhy = Min(hydro.w);
			if (minw > minwhy) {
				ihminw = ihy;
				minw = minwhy;
			}
		}
	}
	if (ihminw >= 0) {
		for (int i = 0; i < 6*Nb; ++i) 
			for (int j = 0; j < 6*Nb; ++j) 
				Aw0(i, j) = hydrosp[ihminw]->Aw0_ndim(i, j);	
	}*/
	
	//bem->calcAinf = true;
	//bem->calcAinf_w = true;
	
	// sts has to be recalculated	
}

void Hydro::Report() const {
	BEM::Print("\n" + Format(t_("%s file '%s'"), GetCodeStr(), dt.file));
	BEM::Print("\n" + Format(t_("g [m/s2]: %s, h [m]: %s, rho [kg/m3]: %s, length scale [m]: %s"), 
								S_g(), S_h(), S_rho(), S_len()));
	String freqs;
	if (dt.w.IsEmpty()) 
		freqs = t_("NONE");
	else if (dt.w.size() > 1) {
		String strDeltaH;
		if (GetIrregularFreq() < 0) 
			strDeltaH = Format(t_("delta %s [rad/s]"), FDS(dt.w[1] - dt.w[0], 8, false));
		else {
			String strHead;
			for (int i = 0; i < dt.w.size(); ++i) {
				if (i > 0)
					strHead << ", ";
				strHead << dt.w[i];
			}
			strDeltaH = Format(t_("Non constant delta (%s)"), strHead); 
		}
	 	freqs = Format(t_("%s to %s %s"), FDS(dt.w[0], 8, false), 
	 									  FDS(dt.w[dt.w.size()-1], 8, false), strDeltaH);	
	} else
		freqs = Format(t_("%s [rad/s]"), FDS(dt.w[0], 8, false));
	
	String heads;
	if (dt.head.IsEmpty())
		heads = t_("NONE");
	else if (dt.head.size() > 1) {
		String strDeltaH;
		if (GetIrregularHead() < 0) 
			strDeltaH = Format(t_("delta %.1f [º]"), dt.head[1] - dt.head[0]);
		else {
			String strHead;
			for (int i = 0; i < dt.head.size(); ++i) {
				if (i > 0)
					strHead << ", ";
				strHead << dt.head[i];
			}
			strDeltaH = Format(t_("Non constant delta (%s)"), strHead); 
		}
	 	heads = Format(t_("%.1f to %.1f %s"), dt.head[0], dt.head[dt.head.size()-1], strDeltaH);	
	} else
		heads = Format(t_("%.1f [º]"), dt.head[0]);
	
	BEM::Print("\n" + Format(t_("#freqs: %d (%s)"), dt.Nf, freqs)); 
	BEM::Print("\n" + Format(t_("#1st order headings: %d (%s)"), dt.Nh, heads)); 
	BEM::Print("\n" + Format(t_("#bodies: %d"), dt.Nb));
	for (int ib = 0; ib < dt.Nb; ++ib) {
		String str = Format("\n%d.", ib+1);
		str += " '" + dt.msh[ib].dt.name + "'";
		//if (dof.size() > ib)
		//	str += S(" ") + t_("dof") + ": " + FormatInt(dof[ib]);
		if (/*Vo.size() > ib && */IsNum(dt.msh[ib].dt.Vo))
			str += S(" ") + t_("vol [m3]") + ": " + FDS(dt.msh[ib].dt.Vo, 8, false);
		if (IsNum(dt.msh[ib].dt.cg))
			str += " " + Format("Cg(%.3f, %.3f, %.3f)[m]", dt.msh[ib].dt.cg.x, dt.msh[ib].dt.cg.y, dt.msh[ib].dt.cg.z);
		if (IsNum(dt.msh[ib].dt.cb))
			str += " " + Format("Cb(%.3f, %.3f, %.3f)[m]", dt.msh[ib].dt.cb.x, dt.msh[ib].dt.cb.y, dt.msh[ib].dt.cb.z);
		if (IsNum(dt.msh[ib].dt.c0))
			str += " " + Format("C0(%.3f, %.3f, %.3f)[m]", dt.msh[ib].dt.c0.x, dt.msh[ib].dt.c0.y, dt.msh[ib].dt.c0.z);
		
		BEM::Print(str);
	}
}
/*
void Hydro::GetBodyDOF() {
	dof.Clear();	 dof.SetCount(Nb, 0);
	for (int ib = 0; ib < Nb; ++ib)
		for (int idf = 0; idf < 6; ++idf)
			if (IsAvailableDOF(ib, idf))
				dof[ib]++;
}
*/
String Hydro::AfterLoad(Function <bool(String, int)> Status) {
	SortFrequencies();
	SortHeadings();
	
	if (!IsLoadedA0())  
		GetA0();
	
	if ((!IsLoadedAinf() || !IsLoadedKirf()) && Bem().calcAinf) {
		if (!IsNum(Bem().maxTimeA) || Bem().maxTimeA == 0) 
			return t_("Incorrect time for A∞ calculation. Please review it in Options");
		if (!IsNum(Bem().numValsA) || Bem().numValsA < 10) 
			return t_("Incorrect number of time values for A∞ calculation. Please review it in Options");
		if (!IsLoadedKirf()) {
			if (Status && !Status(t_("Obtaining the Impulse Response Function"), 40)) 
				return t_("Cancelled by the user");
			GetK_IRF(min(Bem().maxTimeA, GetK_IRF_MaxT()), Bem().numValsA);
		}
		if (!IsLoadedAinf()) {
			if (Status && !Status(t_("Obtaining the infinite-frequency added mass (A∞)"), 70)) 
				return t_("Cancelled by the user");
			GetAinf();
		}
	}
	if (Bem().calcAinf_w) {
		if (!IsLoadedKirf())
			GetK_IRF(min(Bem().maxTimeA, GetK_IRF_MaxT()), Bem().numValsA);
		if (!IsLoadedAinf())
			GetAinf();
		if (Status && !Status(t_("Obtaining the frequency-dependent infinite-frequency added mass (A∞(ω))"), 90)) 
			return t_("Cancelled by the user");
		GetAinf_w();
	}
	
	CompleteForces1st();
	
//	if (Ainf_w.size() == 0)
//		InitAinf_w();
//	if (Ainf.size() == 0)
//		Ainf.setConstant(Nb*6, Nb*6, 0);
	
	/*try {
		CheckNaN();
	} catch (Exc e) {
		lastError = e;
		return false;
	}*/
	
	// Fill the other side of the diagonal. If Null, fill with zero
	auto FillNullQTF = [&](UArray<UArray<UArray<MatrixXcd>>> &qtf, bool isSum) {
		for (int ib = 0; ib < dt.Nb; ++ib) {
	        for (int ih = 0; ih < dt.qh.size(); ++ih) {
				for (int idf = 0; idf < 6; ++idf) { 
					MatrixXcd &c = qtf[ib][ih][idf];
					Eigen::Index rows = c.rows();
					for (int iw = 0; iw < rows; ++iw) {
						for (int jw = iw+1; jw < rows; ++jw) {
							std::complex<double> &cij = c(iw, jw), &cji = c(jw, iw);
							if (IsNull(cij)) {
								if (IsNull(cji))
									cij = cji = 0;
								else {
									if (isSum) 
										cij = cji;
									else
										cij = std::complex<double>(cji.real(), -cji.imag());;
								}
							} else {
								if (IsNull(cji)) {
									if (isSum) 
										cji = cij;
									else
										cji = std::complex<double>(cij.real(), -cij.imag());
								}
							}
						}
					}
				}
	        }
		}
	};
	if (IsLoadedQTF(true)) 
		FillNullQTF(dt.qtfsum, true);
	if (IsLoadedQTF(false))
		FillNullQTF(dt.qtfdif, false);
	
	for (int ib = 0; ib < dt.msh.size(); ++ib) {
		Mesh &m = dt.msh[ib];
		String ret;
		if (!m.dt.mesh.IsEmpty())
			ret = m.dt.mesh.CheckErrors();
		if (!ret.IsEmpty())
			return ret;
		
		m.dt.fileName = dt.file;
		if (m.dt.name.IsEmpty())
			m.dt.name = InitCaps(GetFileTitle(dt.file));

		//m.AfterLoad(dt.rho, dt.g, false, true, true);

		// Symmetrize radiation potentials
		if (dt.symX) {
			int npan = m.dt.mesh.panels.size();
			m.dt.mesh.DeployXSymmetry();
			if (IsLoadedPotsRad()) {
				dt.pots_rad[ib].SetCount(2*npan);
				for (int ipan = 0; ipan < npan; ++ipan) {
					dt.pots_rad[ib][ipan + npan].SetCount(6);
					for (int idf = 0; idf < 6; ++idf) 
						dt.pots_rad[ib][ipan + npan][idf].SetCount(dt.Nf);
					for (int ifr = 0; ifr < dt.Nf; ++ifr) {
						dt.pots_rad[ib][ipan + npan][0][ifr] = -dt.pots_rad[ib][ipan][0][ifr];
						dt.pots_rad[ib][ipan + npan][1][ifr] =  dt.pots_rad[ib][ipan][1][ifr];
						dt.pots_rad[ib][ipan + npan][2][ifr] =  dt.pots_rad[ib][ipan][2][ifr];
						dt.pots_rad[ib][ipan + npan][3][ifr] =  dt.pots_rad[ib][ipan][3][ifr];
						dt.pots_rad[ib][ipan + npan][4][ifr] = -dt.pots_rad[ib][ipan][4][ifr];
						dt.pots_rad[ib][ipan + npan][5][ifr] = -dt.pots_rad[ib][ipan][5][ifr];
					}
				}
			}
		}
		if (dt.symY) {
			int npan = m.dt.mesh.panels.size();
			m.dt.mesh.DeployYSymmetry();
			if (IsLoadedPotsRad()) {
				dt.pots_rad[ib].SetCount(2*npan);
				for (int ipan = 0; ipan < npan; ++ipan) {
					dt.pots_rad[ib][ipan + npan].SetCount(6);
					for (int idf = 0; idf < 6; ++idf) 
						dt.pots_rad[ib][ipan + npan][idf].SetCount(dt.Nf);
					for (int ifr = 0; ifr < dt.Nf; ++ifr) {
						dt.pots_rad[ib][ipan + npan][0][ifr] =  dt.pots_rad[ib][ipan][0][ifr];
						dt.pots_rad[ib][ipan + npan][1][ifr] = -dt.pots_rad[ib][ipan][1][ifr];
						dt.pots_rad[ib][ipan + npan][2][ifr] =  dt.pots_rad[ib][ipan][2][ifr];
						dt.pots_rad[ib][ipan + npan][3][ifr] = -dt.pots_rad[ib][ipan][3][ifr];
						dt.pots_rad[ib][ipan + npan][4][ifr] =  dt.pots_rad[ib][ipan][4][ifr];
						dt.pots_rad[ib][ipan + npan][5][ifr] = -dt.pots_rad[ib][ipan][5][ifr];
					}
				}
			}
		}
		m.AfterLoad(dt.rho, dt.g, false, true);
	}
	if (IsLoadedPotsRad()) {
		Initialize_AB(dt.A_P, 0);
		Initialize_AB(dt.B_P, 0);
		
		dt.Apan = Tensor<double, 5>(dt.Nb, dt.pots_rad[0].size(), 6, 6, dt.Nf);
		dt.Bpan = Tensor<double, 5>(dt.Nb, dt.pots_rad[0].size(), 6, 6, dt.Nf);
		
		UVector<double> n(6);
		for (int ib = 0; ib < dt.Nb; ++ib) {
			const Value3D &c0_ = dt.msh[ib].dt.c0;
			for (int ip = 0; ip < dt.pots_rad[ib].size(); ++ip) {
				const Panel &pan = dt.msh[ib].dt.mesh.panels[ip];
				double s = pan.surface0 + pan.surface1;
				const Value3D &n1 = pan.normalPaint;
				n[0] = n1[0];	n[1] = n1[1];	n[2] = n1[2];
				Value3D r = pan.centroidPaint - c0_;
				Value3D n2 = r%n1;
				n[3] = n2[0];	n[4] = n2[1];	n[5] = n2[2];
				for (int ifr = 0; ifr < dt.Nf; ++ifr) {
					double rho_w = dt.rho/dt.w[ifr];
						
					for (int idf2 = 0; idf2 < 6; ++idf2) {
						const std::complex<double> &comp = dt.pots_rad[ib][ip][idf2][ifr];
						for (int idf1 = 0; idf1 < 6; ++idf1) {
							dt.A_P[idf1 + ib*6][idf2 + ib*6][ifr] += (dt.Apan(ib, ip, idf1, idf2, ifr) = rho_w*comp.imag()*n[idf1]*s);
							dt.B_P[idf1 + ib*6][idf2 + ib*6][ifr] -= (dt.Bpan(ib, ip, idf1, idf2, ifr) = dt.rho*comp.real()*n[idf1]*s);
						}
					}
				}
			}
		}
	}
	if (!dt.msh.IsEmpty()) {
		Initialize_PotsIncDiff(dt.pots_inc);

    	for (int ifr = 0; ifr < dt.Nf; ++ifr) {		
    		double k = SeaWaves::WaveNumber_w(dt.w[ifr], dt.h, g_dim());
    		double g_w = g_dim()/dt.w[ifr];
			for (int ib = 0; ib < dt.Nb; ++ib) {
				Mesh &m = dt.msh[ib];
				int npan = m.dt.mesh.panels.size();
				for (int ip = 0; ip < npan; ++ip) {
					const Panel &pan = dt.msh[ib].dt.mesh.panels[ip];
					const Point3D &p = pan.centroidPaint;
					for (int ih = 0; ih < dt.Nh; ++ih) {	
						if (p.z >= 0) {
							double th = ToRad(dt.head[ih]);
							double cs;
							if (dt.h > 0 && k*dt.h < 700) 
							 	cs = cosh(k*(p.z + dt.h))/cosh(k*dt.h);
							else
								cs = exp(k*p.z);
							double ex = k*(p.x*cos(th) + p.y*sin(th));
							
							dt.pots_inc[ib][ip][ih][ifr] = g_w*cs*std::complex<double>(sin(ex), -cos(ex));
						}
					}
				}
			}
    	}
    	
    	Initialize_Forces(dt.fk_pot, -1, 0);
    	
    	const std::complex<double> i = std::complex<double>(0, 1);
    	
    	UVector<double> n(6);
		for (int ih = 0; ih < dt.Nh; ++ih) {
			for (int ifr = 0; ifr < dt.Nf; ++ifr) {
				double rho_w = dt.rho/dt.w[ifr];
				for (int ib = 0; ib < dt.Nb; ++ib) {
					const Value3D &c0 = dt.msh[ib].dt.c0;
					for (int ip = 0; ip < dt.pots_inc[ib].size(); ++ip) {
						const Panel &pan = dt.msh[ib].dt.mesh.panels[ip];
						double s = pan.surface0 + pan.surface1;
						const Value3D &n13 = pan.normalPaint;
						n[0] = n13[0];	n[1] = n13[1];	n[2] = n13[2];
						Value3D r = pan.centroidPaint - c0;
						Value3D n26 = r%n13;
						n[3] = n26[0];	n[4] = n26[1];	n[5] = n26[2];	
												
						for (int idf = 0; idf < 6; ++idf) 
							dt.fk_pot.force[ih](ifr, idf + 6*ib) += rho_w*dt.pots_inc[ib][ip][ih][ifr]*n[idf]*s*i;
					}
				}
			}
		}
	}

	return String();
}

int Hydro::GetW0() {
	for (int i = 0; i < dt.w.size(); ++i) {
		if (dt.w[i] < 0.0001)
			return i;
	}
	return Null;
}

void Hydro::Get3W0(int &id1, int &id2, int &id3) {
	UVector<double> ww = clone(dt.w);
	
	Sort(ww);
	id1 = FindAdd(dt.w, ww[0]);
	id2 = FindAdd(dt.w, ww[1]); 
	id3 = FindAdd(dt.w, ww[2]); 
}

void Hydro::GetA0() {
	if (!IsLoadedA())
		return;
	
	int iw0 = GetW0();
	if (IsNum(iw0)) {
		dt.A0.setConstant(dt.Nb*6, dt.Nb*6, Null);
		for (int i = 0; i < dt.Nb*6; ++i)
	        for (int j = 0; j < dt.Nb*6; ++j)
				dt.A0(i, j) = dt.A[i][j][iw0];
	} else if (dt.w.size() < 3)
		return;
	else {
		int iw1, iw2, iw3;
		Get3W0(iw1, iw2, iw3);
		double wiw1 = dt.w[iw1];
		double wiw2 = dt.w[iw2];
		double wiw3 = dt.w[iw3];
		if (wiw1 > 0.5 && wiw1 > 3*(wiw2 - wiw1))	// Too high to guess A[0]
			return;
		dt.A0.setConstant(dt.Nb*6, dt.Nb*6, Null);
		for (int i = 0; i < dt.Nb*6; ++i)
	        for (int j = 0; j < dt.Nb*6; ++j) {
				if (!IsNum(dt.A[i][j][iw1]) || !IsNum(dt.A[i][j][iw2]) || !IsNum(dt.A[i][j][iw3]))
	                dt.A0(i, j) = Null;
	            else {
	                double val = QuadraticInterpolate<double>(0, wiw1, wiw2, wiw3, dt.A[i][j][iw1], dt.A[i][j][iw2], dt.A[i][j][iw3]);
	        		if (abs(val) < 0.0000001) {
			    		if (abs(dt.A[i][j][iw1]) < 0.0000001) 
	                        dt.A0(i, j) = val;
	                    else
	                    	dt.A0(i, j) = Null;
	                } else if (abs(1 - val/dt.A[i][j][iw1]) > 0.3)
	                    dt.A0(i, j) = Null;	// Too far to be good
	                else
						dt.A0(i, j) = val;
	            }
	        }
	}
}

String Hydro::C_units_base(int i, int j) {
	if (i == 2 && j == 2)
		return "N/m";
	else if ((i == 2 && j == 3) || (i == 2 && j == 4))
		return "N/rad";
	else if (i >= 3 && j >= 3)
		return "Nm/rad";
	else
		return "";
}

String Hydro::C_units(int i, int j) {
	String ret = C_units_base(i, j);
	if (ret.IsEmpty())
		return C_units_base(j, i);
	return ret;
}

void Hydro::SetC(int ib, const MatrixXd &K) {		// K is supposed to be dimensionalized
	//if (C.IsEmpty())
	//	C.SetCount(Nb);
	if (dt.msh[ib].dt.C.size() == 0)
		dt.msh[ib].dt.C.setConstant(6, 6, Null); 
	for (int idf = 0; idf < 6; ++idf) {
		for (int jdf = 0; jdf < 6; ++jdf) {
			double k = dt.dimen ? g_rho_dim()/g_rho_ndim() : g_rho_ndim()*pow(dt.len, GetK_C(idf, jdf));
	      	dt.msh[ib].dt.C(idf, jdf) = K(idf, jdf)/k;
		}
	}
}

int Hydro::GetIrregularHead() const {
	if (dt.Nh <= 2)
		return -1;
	double delta0 = dt.head[1] - dt.head[0];
	for (int i = 1; i < dt.Nh - 1; ++i) {
		double delta = dt.head[i+1] - dt.head[i];
		if (!EqualRatio(delta, delta0, 0.001))
			return i;
	}
	return -1;
}

int Hydro::GetIrregularFreq() const {
	if (dt.Nf <= 2)
		return -1;
	double delta0 = dt.w[1] - dt.w[0];
	for (int i = 1; i < dt.Nf - 1; ++i) {
		double delta = dt.w[i+1] - dt.w[i];
		if (!EqualRatio(delta, delta0, delta0/10))
			return i;
	}
	return -1;
}

double Hydro::g_dim() 		const {return Bem().g;}							// Dimensionalize only with system data
double Hydro::g_ndim()		const {return IsNum(dt.g) ? dt.g : Bem().g;}	// Nondimensionalize with model data, if possible
double Hydro::rho_dim() 	const {return Bem().rho;}		
double Hydro::rho_ndim()	const {return IsNum(dt.rho) ? dt.rho : Bem().rho;}
double Hydro::g_rho_dim() 	const {return Bem().rho*Bem().g;}
double Hydro::g_rho_ndim()	const {return g_ndim()*rho_ndim();}

void Hydro::StateSpace::GetTFS(const UVector<double> &ww) {
	Eigen::Index sz = A_ss.rows();
	TFS.SetCount(ww.size());
	for (int ifr = 0; ifr < ww.size(); ++ifr) {
		std::complex<double> wi = std::complex<double>(0, ww[ifr]);
		
		MatrixXcd Iwi_A = MatrixXd::Identity(sz, sz)*wi - A_ss;
		
		if (FullPivLU<MatrixXcd>(Iwi_A).isInvertible())
			TFS[ifr] = (C_ss.transpose()*Iwi_A.inverse()*B_ss)(0);	// C_ss*inv(I*w*i-A_ss)*B_ss
		else
			TFS[ifr] = Null;
	}		
}

int Hydro::GetHeadId(double hd) const {
	hd = FixHeading_180(hd);
	for (int i = 0; i < dt.head.size(); ++i) {
		if (EqualRatio(FixHeading_180(dt.head[i]), hd, 0.01))
			return i;
	}
	return -1;
}

int Hydro::GetHeadIdMD(const std::complex<double> &hd) const {
	std::complex<double> hd_ = FixHeading_180(hd);
	for (int i = 0; i < dt.mdhead.size(); ++i) {
		if (EqualRatio(FixHeading_180(dt.mdhead[i]), hd_, 0.01))
			return i;
	}
	return -1;
}

VectorXd Hydro::B_dim(int idf, int jdf) const {
	if (dt.dimen)
		return dt.B[idf][jdf]*(rho_dim()/rho_ndim());
	else {
		VectorXd ret = dt.B[idf][jdf]*(rho_dim()*pow(dt.len, GetK_AB(idf, jdf)));
		VectorXd ww = Get_w();
		return ret.array()*ww.array();
	}
}

VectorXd Hydro::B_ndim(int idf, int jdf) const {
	if (dt.B[idf][jdf].size() == 0)
		return VectorXd();
	if (!dt.dimen)
		return dt.B[idf][jdf]*(rho_ndim()/rho_dim());
	else {
		VectorXd ret = dt.B[idf][jdf]/(rho_ndim()*pow(dt.len, GetK_AB(idf, jdf)));
		VectorXd ww = Get_w();
		return ret.array()/ww.array();
	}
}


void Hydro::GetOldAB(const UArray<MatrixXd> &oldAB, UArray<UArray<VectorXd>> &AB) {
	AB.Clear();
	int Nf = oldAB.size();
	int Nb = 0;
	if (Nf > 0)
		Nb = int(oldAB[0].rows())/6;
	AB.SetCount(6*Nb);
	for (int i = 0; i < 6*Nb; ++i) {
		AB[i].SetCount(6*Nb);
		for (int j = 0; j < 6*Nb; ++j) {
			AB[i][j].resize(Nf);	
			for (int idf = 0; idf < Nf; ++idf) 
				AB[i][j][idf] = oldAB[idf](i, j);	
		}
	}
}

void Hydro::SetOldAB(UArray<MatrixXd> &oldAB, const UArray<UArray<VectorXd>> &AB) {
	oldAB.Clear();
	int Nb = AB.size()/6;
	int Nf = 0;
	if (Nb > 0)
		Nf = int(AB[0][0].size());
	oldAB.SetCount(Nf);
	for (int idf = 0; idf < Nf; ++idf) {
		oldAB[idf].resize(6*Nb, 6*Nb);
		for (int i = 0; i < 6*Nb; ++i) 
			for (int j = 0; j < 6*Nb; ++j) 
				oldAB[idf](i, j) = AB[i][j][idf];
	}
}

MatrixXd Hydro::A_mat(bool ndim, int ifr, int ib1, int ib2) const {
	MatrixXd ret;
	if (!IsLoadedA())
		return ret;
	ret.resize(6, 6);
	for (int idf = 0; idf < 6; ++idf) 	
		for (int jdf = 0; jdf < 6; ++jdf) 
			ret(idf, jdf) = A_(ndim, ifr, idf + 6*ib1, jdf + 6*ib2);
	return ret;
}

MatrixXd Hydro::Ainf_mat(bool ndim, int ib1, int ib2) const {
	MatrixXd ret;
	if (!IsLoadedAinf())
		return ret;
	ret.resize(6, 6);
	for (int idf = 0; idf < 6; ++idf) 	
		for (int jdf = 0; jdf < 6; ++jdf) 
			ret(idf, jdf) = Ainf_(ndim, idf + 6*ib1, jdf + 6*ib2);
	return ret;
}

MatrixXd Hydro::B_mat(bool ndim, int ifr, int ib1, int ib2) const {
	MatrixXd ret;
	if (!IsLoadedA())
		return ret;
	ret.resize(6, 6);
	for (int idf = 0; idf < 6; ++idf) 	
		for (int jdf = 0; jdf < 6; ++jdf) 
			ret(idf, jdf) = B_(ndim, ifr, idf + 6*ib1, jdf + 6*ib2);
	return ret;
}

MatrixXd Hydro::C_(bool ndim, int ib) const {
	MatrixXd ret;
	if (dt.msh[ib].dt.C.size() == 0)
		return ret;
	ret.resize(6, 6);
	for (int idf = 0; idf < 6; ++idf) 	
		for (int jdf = 0; jdf < 6; ++jdf) 
			ret(idf, jdf) = C_(ndim, ib, idf, jdf);
	return ret;
}

MatrixXd Hydro::CMoor_(bool ndim, int ib) const {
	MatrixXd ret;
	if (dt.msh[ib].dt.Cmoor.size() == 0)
		return ret;
	ret.resize(6, 6);
	for (int idf = 0; idf < 6; ++idf) 	
		for (int jdf = 0; jdf < 6; ++jdf) 
			ret(idf, jdf) = CMoor_(ndim, ib, idf, jdf);
	return ret;
}

MatrixXd Hydro::Dlin_dim(int ib) const {
	MatrixXd ret;
	if (!IsLoadedDlin())
		return ret;
	ret.resize(6, 6);
	for (int idf = 0; idf < 6; ++idf) 	
		for (int jdf = 0; jdf < 6; ++jdf) 
			ret(idf, jdf) = dt.msh[ib].dt.Dlin(idf, jdf);
	return ret;
}

MatrixXd Hydro::Dquad_dim(int ib) const {
	MatrixXd ret;
	if (!IsLoadedDquad())
		return ret;
	ret.resize(6, 6);
	for (int idf = 0; idf < 6; ++idf) 	
		for (int jdf = 0; jdf < 6; ++jdf) 
			ret(idf, jdf) = dt.msh[ib].dt.Dquad(idf, jdf);
	return ret;
}

void Hydro::C_dim() {
	//if (C.IsEmpty())
	//	return;
	for (int ib = 0; ib < dt.Nb; ++ib) 
		for (int idf = 0; idf < 6; ++idf) 	
			for (int jdf = 0; jdf < 6; ++jdf) 
				dt.msh[ib].dt.C(idf, jdf) = C_dim(ib, idf, jdf);
}

void Hydro::F_dim(Forces &f) {
	if (f.force.IsEmpty())
		return;
	for (int ih = 0; ih < dt.Nh; ++ih) 	
		for (int ifr = 0; ifr < dt.Nf; ++ifr)
			for (int idf = 0; idf < 6*dt.Nb; ++idf) 
				f.force[ih](ifr, idf) = F_dim(f, ih, ifr, idf);
}

VectorXcd Hydro::F_(bool ndim, const Forces &f, int _h, int ifr) const {
	VectorXcd ret;
	if (f.force.IsEmpty())
		return ret;
	ret.resize(6);
	for (int idf = 0; idf < 6; ++idf) 
		ret[idf] = F_(ndim, f, _h, ifr, idf);
	return ret;
}

VectorXcd Hydro::F_dof(bool ndim, const Forces &f, int _h, int idf) const {
	VectorXcd ret;
	if (f.force.IsEmpty())
		return ret;
	ret.resize(dt.Nf);
	for (int ifr = 0; ifr < dt.Nf; ++ifr) 
		ret[ifr] = F_(ndim, f, _h, ifr, idf);
	return ret;
}

void Hydro::RAO_dim(RAO &f) {
	if (f.force.IsEmpty())
		return;
	for (int ih = 0; ih < dt.Nh; ++ih) 	
		for (int ifr = 0; ifr < dt.Nf; ++ifr)
			for (int idf = 0; idf < 6*dt.Nb; ++idf) 
				f.force[ih](ifr, idf) = RAO_dim(f, ih, ifr, idf);
}

VectorXcd Hydro::RAO_(bool ndim, const RAO &f, int _h, int ifr) const {
	VectorXcd ret;
	if (f.force.IsEmpty())
		return ret;
	ret.resize(6);
	for (int idf = 0; idf < 6; ++idf) 
		ret[idf] = RAO_(ndim, f, _h, ifr, idf);
	return ret;
}

VectorXcd Hydro::RAO_dof(bool ndim, int _h, int idf) const {
	return RAO_dof(ndim, dt.rao, _h, idf);
}

VectorXcd Hydro::RAO_dof(bool ndim, const RAO &f, int _h, int idf) const {
	VectorXcd ret;
	if (f.force.IsEmpty())
		return ret;
	ret.resize(dt.Nf);
	for (int ifr = 0; ifr < dt.Nf; ++ifr) 
		ret[ifr] = RAO_(ndim, f, _h, ifr, idf);
	return ret;
}

VectorXd Hydro::Md_dof(bool ndim, int _h, int idf) const {
	VectorXd ret;
	if (dt.md.IsEmpty())
		return ret;
	ret.resize(dt.Nf);
	for (int ifr = 0; ifr < dt.Nf; ++ifr) 
		ret[ifr] = Md_(ndim, idf, _h, ifr);
	return ret;
}

MatrixXcd Hydro::QTF_dof(bool ndim, bool isSum, int _h, int idf, int ib) const {
	const UArray<UArray<UArray<MatrixXcd>>> &qtf = isSum ? dt.qtfsum : dt.qtfdif;
	
	MatrixXcd ret;
	if (qtf.IsEmpty())
		return ret;
	ret.resize(dt.qw.size(), dt.qw.size());
	for (int ifr1 = 0; ifr1 < dt.qw.size(); ++ifr1) 
		for (int ifr2 = 0; ifr2 < dt.qw.size(); ++ifr2) 
			ret(ifr1, ifr2) = F_(ndim, qtf[ib][_h][idf](ifr1, ifr2), idf);
	return ret;
}


void Hydro::CheckNaN() {
	if (!IsNum(dt.A))
		throw Exc("Error loading A. NaN found");
	if (!IsNum(dt.Ainf_w))
		throw Exc("Error loading Ainfw. NaN found");
	if (!IsNum(dt.Ainf))
		throw Exc("Error loading Awinf. NaN found");
	if (!IsNum(dt.A0))
		throw Exc("Error loading A_0. NaN found");
	if (!IsNum(dt.B))
		throw Exc("Error loading B. NaN found");
	if (!IsNum(dt.head))
		throw Exc("Error loading head. NaN found");
//	if (!IsNum(M))
//		throw Exc("Error loading M. NaN found");
//	if (!IsNum(C))
//		throw Exc("Error loading C. NaN found");
	//if (!IsNum(dof))
		//throw Exc("Error loading dof. NaN found");
	if (!IsNum(dt.Kirf))
		throw Exc("Error loading Kirf. NaN found");
	if (!IsNum(dt.Tirf))
		throw Exc("Error loading Tirf. NaN found");
	if (!IsNum(dt.ex))
		throw Exc("Error loading ex. NaN found");
	if (!IsNum(dt.sc))
		throw Exc("Error loading sc. NaN found");
	if (!IsNum(dt.fk))
		throw Exc("Error loading fk. NaN found");
	if (!IsNum(dt.rao))
		throw Exc("Error loading rao. NaN found");	
	
	for (int ib = 0; ib < dt.Nb; ++ib) {
		if (!IsNum(dt.msh[ib].dt.cb))
			throw Exc("Error loading cb. NaN found");
		if (!IsNum(dt.msh[ib].dt.cg))
			throw Exc("Error loading cg. NaN found");
		if (!IsNum(dt.msh[ib].dt.c0))
			throw Exc("Error loading c0. NaN found");
		
		if (!IsNum(dt.msh[ib].dt.M))
			throw Exc("Error loading M. NaN found");
		if (!IsNum(dt.msh[ib].dt.C))
			throw Exc("Error loading C. NaN found");
	}
}

double Hydro::Tdof(int ib, int idf) const {
	if (!IsLoadedAinf(idf+6*ib, idf+6*ib) || !IsLoadedM(ib, idf, idf) || !IsLoadedC(ib, idf, idf))
		return Null;
	double c = C_dim(ib, idf, idf);
	if (c == 0)
		return Null;
	double m = dt.msh[ib].dt.M(idf, idf);
	double a = Ainf_dim(ib*6+idf, ib*6+idf);
	double d = (m + a)/c;
	if (d < 0)
		return Null;
	return 2*M_PI*sqrt(d); 
}

double Hydro::Tdofw(int ib, int idf) const {
	if (!IsLoadedA(idf+6*ib, idf+6*ib) || !IsLoadedM(ib, idf, idf) || !IsLoadedC(ib, idf, idf))
		return Null;
	double c = C_dim(ib, idf, idf);
	if (c == 0)
		return Null;
	double m = dt.msh[ib].dt.M(idf, idf);
	VectorXd nw(dt.w.size());
	for (int ifr = 0; ifr < dt.w.size(); ++ifr) {
		double a = A_dim(ifr, ib*6+idf, ib*6+idf);
		double d = (m + a)/c;
		if (d < 0)
			return Null;
		nw[ifr] = 1/sqrt(d); 
	}
	VectorXd delta = Get_w() - nw;
	
	UVector<double> zeros; 
	ZeroCrossing(Get_w(), delta, true, true, zeros);
	if (zeros.IsEmpty())
		return Null;
	return 2*M_PI/First(zeros); 
}

double Hydro::Theave(int ib)  const {return Tdof(ib, 2);}
double Hydro::Theavew(int ib) const {return Tdofw(ib, 2);}
double Hydro::Troll(int ib)   const {return Tdof(ib, 3);}
double Hydro::Trollw(int ib)  const {return Tdofw(ib, 3);}
double Hydro::Tpitch(int ib)  const {return Tdof(ib, 4);}
double Hydro::Tpitchw(int ib) const {return Tdofw(ib, 4);}

double Hydro::GM(int ib, int idf) const {
	if (/*Vo.size() <= ib*/ !IsNum(dt.msh[ib].dt.Vo) || dt.msh[ib].dt.Vo == 0)
		return Null;
	double den = rho_dim()*g_dim()*dt.msh[ib].dt.Vo;
	if (den == 0)
		return Null;
	if (IsLoadedC(ib))
		return C_dim(ib, idf, idf)/den;
	else
		return Null;
}

double Hydro::GMroll(int ib) const {
	return GM(ib, 3);
}

double Hydro::GMpitch(int ib) const {
	return GM(ib, 4);
}

void Hydro::Jsonize(JsonIO &json) {
	int icode;
	UArray<MatrixXd> oldA, oldB, oldKirf;
	if (json.IsStoring()) {
		icode = dt.solver;
		SetOldAB(oldA, dt.A);
		SetOldAB(oldB, dt.B);
		SetOldAB(oldKirf, dt.Kirf);
	}
	json
		("file", dt.file)
		("name", dt.name)
		("g", dt.g)
		("h", dt.h)
		("rho", dt.rho)
		("len", dt.len)
		("dimen", dt.dimen)
		("Nb", dt.Nb)
		("Nf", dt.Nf)
		("Nh", dt.Nh)
		("A", oldA)
		("Awinf", dt.Ainf)
		("Aw0", dt.A0)
		("B", oldB)
		("head", dt.head)
		//("names", names)
		//("C", C)
		//("cb", cb)
		//("cg", cg)
		//("c0", c0)
		("code", icode)
		//("dof", dof)
//		("dofOrder", dofOrder)
		("Kirf", oldKirf)
		("Tirf", dt.Tirf)
		("ex", dt.ex)
		("sc", dt.sc)
		("fk", dt.fk)
		("rao", dt.rao)
		("sts", dt.sts)
		//("T", T)
		("w", dt.w)
		("dataFromW", dt.dataFromW)
		//("Vo", Vo)
		("stsProcessor", dt.stsProcessor)
		("dimenSTS", dt.dimenSTS)
		("description", dt.description)
		("qtfsum", dt.qtfsum)
		("qtfdif", dt.qtfdif)
		//("Dlin", Dlin)
		("msh", dt.msh)
	;
	if(json.IsLoading()) {
		dt.solver = static_cast<Hydro::BEM_FMT>(icode);
		GetOldAB(oldA, dt.A);
		GetOldAB(oldB, dt.B);
		GetOldAB(oldKirf, dt.Kirf);
	}
}