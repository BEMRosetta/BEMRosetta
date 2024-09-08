// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2024, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#include <ScatterDraw/DataSource.h>
#include <ScatterDraw/Equation.h>
#include "functions.h"
//#include <SysInfo/Crash.h>
#include <STEM4U/SeaWaves.h>
#include <MatIO/matio.h>

using namespace Upp;
using namespace Eigen;


const char *Hydro::strDataToPlot[] = {t_("A(ω)"), t_("A∞"), t_("A₀"), t_("B(ω)"), t_("A∞(ω)"), t_("Kirf"),
				t_("|Fsc|"), t_("arg(Fsc)"), t_("|Ffk|"), t_("arg(Ffk)"), t_("|Fex|"), t_("arg(Fex)"),
				t_("|RAO|"), t_("arg(RAO)"), t_("|Z|"), t_("arg(Z)"), t_("|Kr|"), t_("arg(Kr)"), 
				t_("|TFS|"), t_("arg(TFS)")};

// enum BEM_FMT 					  {WAMIT, 		  WAMIT_1_3, 					FAST_WAMIT, 				 	   HAMS_WAMIT,  HAMS,   WADAM_WAMIT,   NEMOH,     NEMOHv115,    NEMOHv3,    SEAFEM_NEMOH,   AQWA,   			    AQWA_QTF,	  AQWA_DAT, 	FOAMM,   DIODORE,			BEMROSETTA, 	   ORCAFLEX_YML,   		CSV_MAT,    CSV_TABLE,    BEMIOH5,		CAPYTAINE, 			HYDROSTAR_OUT, 	CAPYNC, 		ORCAWAVE_YML, 		CAPYTAINE_PY, 	 BEMROSETTA_H5,		AKSELOS_NPZ,	ORCAWAVE.owr,		UNKNOWN, NUMBEM};
const char *Hydro::bemStr[]         = {"Wamit .out", "Wamit .1.2.3.hst.7.8.9.ss.12", "FAST .dat.1.2.3.hst.789.ss.12", "HAMS Wamit", "HAMS", "Wadam Wamit","Nemoh v2", "Nemoh v115", "Nemoh v3", "SeaFEM Nemoh","AQWA .lis .ah1 .qtf", 	"AQWA .qtf",  "AQWA .dat", "FOAMM", "Diodore .hdb", 	"BEMRosetta .bem", "OrcaFlex .yml", 	".csv mat", ".csv table", "BEMIO .h5",	"Capytaine .cal",    ".out", 		"Capytaine .nc", "OrcaWave .yml",	"Capytaine .py", "BEMRosetta .h5",	"Akselos .npz",
#ifdef PLATFORM_WIN32	
"OrcaWave .owr", 	
#endif
				"By extension"};
const bool Hydro::bemCanSave[] 		= {true, 	      true,	     				     true,		 			 	 	  false,		false,  false,		   false,     false,	 	false,	   	false, 		   false,  					 true, 	 	  false,		false,   true,	  			true,			     false,	     		true, 	     true, 		   true,		false, 	   			false, 			false,			 false,				false,		 	  false,			true,
#ifdef PLATFORM_WIN32	
false,
#endif
				true};       
const char *Hydro::bemExt[]	   		= {"*.out", 	  "*.1",	     				 "*.1",		 			 	      "",		   	"",	    "",		       "",        "", 		   	"",			"",			   "", 				  		 "*.qtf",	  ".dat",		"",      "*.hdb",	  		"*.bem",		   	 "*.yml",	 		"*.csv",    "*.csv", 	   "*.h5",		"",        			"*.out", 		"*.nc", 		 "*.yml",			"*.py",	 		  "*.h5",			"*.npz",
#ifdef PLATFORM_WIN32	
"*.owr",	
#endif		
				"*.*"};       
	
const bool Hydro::caseCanSave[]     = {false, 	      false,	     				 false,		 			 	 	  false,		true,	false,		   true,      true,	 		true,	   	false, 		   false,  					 false, 	 true,			false,   false,	  			false,			 	false,	     		false, 	 	false, 	   		false,		true, 	   			false, 			false,			 true,				true,			  true,				false,
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
			for (int ib = 0; ib < dt.Nb; ++ib) {
				for (int i = 0; i < 6; ++i) {
					if (IsNum(dt.sc[ib][ih](ifr, i)) && IsNum(dt.fk[ib][ih](ifr, i))) 
						dt.ex[ib][ih](ifr, i) = dt.sc[ib][ih](ifr, i) + dt.fk[ib][ih](ifr, i);
				}
			}
		}
	}
}

void Hydro::GetFscFromFexFfk() {
	Initialize_Forces(dt.sc);
	for (int ih = 0; ih < dt.Nh; ++ih) {
		for (int ifr = 0; ifr < dt.Nf; ++ifr) {
			for (int ib = 0; ib < dt.Nb; ++ib) {
				for (int i = 0; i < 6; ++i) {
					if (IsNum(dt.ex[ib][ih](ifr, i)) && IsNum(dt.fk[ib][ih](ifr, i))) 
						dt.sc[ib][ih](ifr, i) = dt.ex[ib][ih](ifr, i) - dt.fk[ib][ih](ifr, i);
				}
			}
		}
	}	
}

void Hydro::GetFfkFromFexFsc() {
	Initialize_Forces(dt.fk);
	for (int ih = 0; ih < dt.Nh; ++ih) {
		for (int ifr = 0; ifr < dt.Nf; ++ifr) {
			for (int ib = 0; ib < dt.Nb; ++ib) {
				for (int i = 0; i < 6; ++i) {
					if (IsNum(dt.ex[ib][ih](ifr, i)) && IsNum(dt.sc[ib][ih](ifr, i))) 
						dt.fk[ib][ih](ifr, i) = dt.ex[ib][ih](ifr, i) - dt.sc[ib][ih](ifr, i);
				}
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

void Hydro::SortHeadings(BasicBEM::HeadingType range, BasicBEM::HeadingType rangeMD, BasicBEM::HeadingType rangeQTF) {
	{	// Forces
		UVector<double> nhead;
		UVector<int> indices;
		for (int ih = 0; ih < dt.head.size(); ++ih) {
			double h = FixHeading(dt.head[ih], range);
			int id = Find(nhead, h);
			if (id < 0) {	// Duplicated headings are discarded
				nhead << h;
				indices << ih;
			}
		}
		if (!Compare(dt.head, nhead)) {
			UVector<int> order = GetSortOrderX(nhead);
			SetSortOrder(nhead, order);
			SetSortOrder(indices, order);
			
			dt.head = pick(nhead);
			dt.Nh = dt.head.size();
				
			auto SortF = [&](Forces &F) {
				Forces f(dt.Nb);
				for (int ib = 0; ib < dt.Nb; ++ib) {
					f[ib].SetCount(dt.Nh);
					for (int ih = 0; ih < dt.Nh; ++ih) {
						f[ib][ih].resize(dt.Nf, 6);
						for (int ifr = 0; ifr < dt.Nf; ++ifr) 
							for (int idf = 0; idf < 6; ++idf) 
								f[ib][ih](ifr, idf) = F[ib][indices[ih]](ifr, idf);
					}
				}
				F = pick(f);
			};
			
			auto SortPotsDifInc = [&](UArray<UArray<UArray<UArray<std::complex<double>>>>> &pots) {
				UArray<UArray<UArray<UArray<std::complex<double>>>>> _pots(dt.Nb);
				for (int ib = 0; ib < dt.Nb; ++ib) {
					_pots[ib].SetCount(pots[ib].size());
					for (int ip = 0; ip < pots[ib].size(); ++ip) {
						_pots[ib][ip].SetCount(dt.Nh);	
						for (int ih = 0; ih < dt.Nh; ++ih) {
							_pots[ib][ip][ih].SetCount(dt.Nf);		
							for (int ifr = 0; ifr < dt.Nf; ++ifr) 
								_pots[ib][ip][ih][ifr] = pots[ib][ip][indices[ih]][ifr];
						}
					}
				}
				pots = pick(_pots);
			};
				
			if (IsLoadedFex())
				SortF(dt.ex);
			if (IsLoadedFsc())
				SortF(dt.sc);
			if (IsLoadedFfk())
				SortF(dt.fk);
			if (IsLoadedRAO()) 
				SortF(dt.rao);	
			
			if (IsLoadedFsc_pot())
				SortF(dt.sc_pot);
			if (IsLoadedFfk_pot())
				SortF(dt.fk_pot);
			if (IsLoadedFfk_pot_bmr())
				SortF(dt.fk_pot_bmr);
			
			if (IsLoadedPotsDif()) 
				SortPotsDifInc(dt.pots_dif);
			if (IsLoadedPotsInc()) 
				SortPotsDifInc(dt.pots_inc);
			if (IsLoadedPotsIncB()) 
				SortPotsDifInc(dt.pots_inc_bmr);
		}
	}
	{	// Mean Drift
		UArray<std::complex<double>> nhead;
		UVector<int> indices;
		for (int ih = 0; ih < dt.mdhead.size(); ++ih) {
			std::complex<double> h = FixHeading(dt.mdhead[ih], rangeMD);
			int id = Find(nhead, h);
			if (id < 0) {	// Duplicated headings are discarded
				nhead << h;
				indices << ih;
			}
		}
		if (!Compare(dt.mdhead, nhead)) {
			UVector<int> order = GetSortOrderX(nhead, SortComplex);
			SetSortOrder(nhead, order);
			SetSortOrder(indices, order);
			
			dt.mdhead.resize(nhead.size());		// Cannot pick()
			for (int ih = 0; ih < dt.mdhead.size(); ++ih)
				dt.mdhead[ih] = nhead[ih];
				
			auto SortMD = [&](UArray<UArray<UArray<VectorXd>>> &MD) {
				UArray<UArray<UArray<VectorXd>>> md(dt.Nb);
				for (int ib = 0; ib < dt.Nb; ++ib) {
					md[ib].SetCount(int(dt.mdhead.size()));
					for (int ih = 0; ih < dt.mdhead.size(); ++ih) {
						md[ib][ih].SetCount(6);
						for (int idf = 0; idf < 6; ++idf) {
							md[ib][ih][idf].resize(dt.Nf);
							for (int ifr = 0; ifr < dt.Nf; ++ifr) 		
								md[ib][ih][idf][ifr] = MD[ib][indices[ih]][idf][ifr];
						}
					}
				}
				MD = pick(md);
			};
			
			if(IsLoadedMD())
				SortMD(dt.md);
		}
	}
	{	// Full QTF
		UArray<std::complex<double>> nhead;
		UVector<int> indices;
		for (int ih = 0; ih < dt.qhead.size(); ++ih) {
			std::complex<double> h = FixHeading(dt.qhead[ih], rangeQTF);
			int id = Find(nhead, h);
			if (id < 0) {	// Duplicated headings are discarded
				nhead << h;
				indices << ih;
			}
		}
		if (!Compare(dt.qhead, nhead)) {
			UVector<int> order = GetSortOrderX(nhead, SortComplex);
			SetSortOrder(nhead, order);
			SetSortOrder(indices, order);		
			
			dt.qhead.resize(nhead.size());		// Cannot pick()
			for (int ih = 0; ih < dt.qhead.size(); ++ih)
				dt.qhead[ih] = nhead[ih];		
			
			auto SortQTF = [&](UArray<UArray<UArray<MatrixXcd>>> &QTF) {
				UArray<UArray<UArray<MatrixXcd>>> qtf(dt.Nb);
				for (int ib = 0; ib < dt.Nb; ++ib) {
					qtf[ib].SetCount(int(dt.qhead.size()));
			        for (int ih = 0; ih < dt.qhead.size(); ++ih) {
			            qtf[ib][ih].SetCount(6);
			        	for (int idf = 0; idf < 6; ++idf) {
			        		qtf[ib][ih][idf].resize(dt.qw.size(), dt.qw.size());
							for (int ifr1 = 0; ifr1 < dt.qw.size(); ++ifr1) 
								for (int ifr2 = 0; ifr2 < dt.qw.size(); ++ifr2) 
									qtf[ib][ih][idf](ifr1, ifr2) = QTF[ib][indices[ih]][idf](ifr1, ifr2);
			        	}
			        }
				}
				QTF = pick(qtf);
			};
				
			if (IsLoadedQTF(true))
				SortQTF(dt.qtfsum);
			if (IsLoadedQTF(false))
				SortQTF(dt.qtfdif);
		}
	}
}

BasicBEM::HeadingType Hydro::ShortestHeadingRange(const UVector<double> &head) {
	if (head.size() == 0)
		return BasicBEM::HEAD_0_360;

	UVector<double> head180, head360;
	for (int ih = 0; ih < head.size(); ++ih) {
		head180 << FixHeading_180(head[ih]);
		head360 << FixHeading_0_360(head[ih]);
	}
	auto MaxAbs = [](const UVector<double> &hd)->double {
		double ret = First(hd);
		for (int i = 0; i < hd.size(); ++i) {
			if (abs(hd[i]) > ret)
				ret = abs(hd[i]);
		}
		return ret;
	};
	double mx180 = MaxAbs(head180);
	double mx360 = MaxAbs(head360);
	if (mx360 <= mx180)
		return BasicBEM::HEAD_0_360;
	else
		return BasicBEM::HEAD_180_180;
}

BasicBEM::HeadingType Hydro::ShortestHeadingRange(const VectorXcd &head) {
	if (head.size() == 0)
		return BasicBEM::HEAD_0_360;
	
	VectorXcd head180(head.size()), head360(head.size());
	for (int ih = 0; ih < head.size(); ++ih) {
		head180[ih] = FixHeading_180(head[ih]);
		head360[ih] = FixHeading_0_360(head[ih]);
	}
	auto MaxAbs = [](const VectorXcd &hd)->std::complex<double> {
		std::complex<double> ret = First(hd);
		for (int i = 0; i < hd.size(); ++i) {
			if (abs(hd[i].real()) > ret.real())
				ret.real(abs(hd[i].real()));
			if (abs(hd[i].imag()) > ret.imag())
				ret.imag(abs(hd[i].imag()));
		}
		return ret;
	};	
	std::complex<double> range180 = MaxAbs(head180);
	std::complex<double> range360 = MaxAbs(head360);
	if ((range360.real() + range360.imag()) <= (range180.real() + range180.imag()))
		return BasicBEM::HEAD_0_360;
	else
		return BasicBEM::HEAD_180_180;
}

void Hydro::SortFrequencies() {
	if (!IsSorted(dt.w)) {
		UVector<int> indices = GetSortOrderX(dt.w);
		dt.w = ApplyIndex(dt.w, indices);
	
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
					for (int ib = 0; ib < dt.Nb; ++ib) 
						for (int idf = 0; idf < 6; ++idf) 
							F[ib][ih](ifr, idf) = f[ib][ih](indices[ifr], idf);
		};
		
		auto SortMD = [&](UArray<UArray<UArray<VectorXd>>> &MD) {
			UArray<UArray<UArray<VectorXd>>> _md = clone(MD);
			for (int ib = 0; ib < dt.Nb; ++ib) 
				for (int ih = 0; ih < dt.mdhead.size(); ++ih)
					for (int idf = 0; idf < 6; ++idf) 
						for (int ifr = 0; ifr < dt.Nf; ++ifr) 		
							MD[ib][ih][idf][ifr] = _md[ib][ih][idf][indices[ifr]];
		};
		
		auto SortPotsRad = [&]() {
			UArray<UArray<UArray<UArray<std::complex<double>>>>> _pots = clone(dt.pots_rad);
			for (int ib = 0; ib < dt.Nb; ++ib) 
				for (int ip = 0; ip < dt.pots_rad[ib].size(); ++ip)
					for (int idf = 0; idf < 6; ++idf) 
						for (int ifr = 0; ifr < dt.Nf; ++ifr)
							dt.pots_rad[ib][ip][idf][ifr] = _pots[ib][ip][idf][indices[ifr]];
		};
		
		auto SortPotsDifInc = [&](UArray<UArray<UArray<UArray<std::complex<double>>>>> &pots) {
			UArray<UArray<UArray<UArray<std::complex<double>>>>> _pots = clone(pots);
			for (int ib = 0; ib < dt.Nb; ++ib) 
				for (int ip = 0; ip < pots[ib].size(); ++ip)
					for (int ih = 0; ih < dt.Nh; ++ih) 
						for (int ifr = 0; ifr < dt.Nf; ++ifr)
							pots[ib][ip][ih][ifr] = _pots[ib][ip][ih][indices[ifr]];		
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
		
		if (IsLoadedFsc_pot())
			SortF(dt.sc_pot);
		if (IsLoadedFfk_pot())
			SortF(dt.fk_pot);
		if (IsLoadedFfk_pot_bmr())
			SortF(dt.fk_pot_bmr);
		
		if(IsLoadedMD())
			SortMD(dt.md);
		
		if (IsLoadedPotsRad())
			SortPotsRad();
		if (IsLoadedPotsDif()) 
			SortPotsDifInc(dt.pots_dif);
		if (IsLoadedPotsInc()) 
			SortPotsDifInc(dt.pots_inc);
		if (IsLoadedPotsIncB()) 
			SortPotsDifInc(dt.pots_inc_bmr);
	}
	if (!IsSorted(dt.qw)) {
		UVector<int> indices = GetSortOrderX(dt.qw);
		dt.qw = ApplyIndex(dt.qw, indices);
		
		auto SortQTF = [&](UArray<UArray<UArray<MatrixXcd>>> &QTF) {
			UArray<UArray<UArray<MatrixXcd>>> qtf = clone(QTF);
			for (int ib = 0; ib < dt.Nb; ++ib) 
		        for (int ih = 0; ih < dt.qhead.size(); ++ih) 
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

int Hydro::Data::FindClosestHead(double hd) const {
	return FindClosest(head, FixHeading_0_360(hd));
}

int Hydro::Data::FindClosestHead(const VectorXcd &list, const std::complex<double> &hd) {
	return FindClosest(list, FixHeading_0_360(hd));
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
		for (int ip = 0; ip < dt.pots_rad[ib].size(); ++ip) {
			a[ib][ip].SetCount(6);
			for (int idf1 = 0; idf1 < 6; ++idf1) {
				a[ib][ip][idf1].SetCount(6);
				for (int idf2 = 0; idf2 < 6; ++idf2) 
					a[ib][ip][idf1][idf2].SetCount(dt.Nf, val);
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
	f.SetCount(dt.Nb);
	for (int ib = 0; ib < dt.Nb; ++ib) {
		f[ib].SetCount(_Nh);
		for (int ih = 0; ih < _Nh; ++ih) 
			f[ib][ih].setConstant(dt.Nf, 6, val);
	}
}
	
void Hydro::Normalize_Forces(Forces &f) {
	for (int ih = 0; ih < dt.Nh; ++ih) 
		for (int ifr = 0; ifr < dt.Nf; ++ifr) 
			for (int ib = 0; ib < dt.Nb; ++ib) 
				for (int idf = 0; idf < 6; ++idf) 
					f[ib][ih](ifr, idf) = F_dim(f, ih, ifr, idf, ib);
}

void Hydro::Normalize_RAO(RAO &f) {
	for (int ih = 0; ih < dt.Nh; ++ih) 
		for (int ifr = 0; ifr < dt.Nf; ++ifr) 
			for (int ib = 0; ib < dt.Nb; ++ib) 
				for (int idf = 0; idf < 6; ++idf) 
					f[ib][ih](ifr, idf) = RAO_dim(f, ih, ifr, idf, ib);
}

void Hydro::Dimensionalize_Forces(Forces &f) {
	for (int ih = 0; ih < dt.Nh; ++ih) 
		for (int ifr = 0; ifr < dt.Nf; ++ifr) 
			for (int ib = 0; ib < dt.Nb; ++ib) 
				for (int idf = 0; idf < 6; ++idf) 
					f[ib][ih](ifr, idf) = F_dim(f, ih, ifr, idf, ib);
}

void Hydro::Dimensionalize_RAO(RAO &f) {
	for (int ih = 0; ih < dt.Nh; ++ih) 
		for (int ifr = 0; ifr < dt.Nf; ++ifr) 
			for (int ib = 0; ib < dt.Nb; ++ib) 
				for (int idf = 0; idf < 6; ++idf) 
					f[ib][ih](ifr, idf) = RAO_dim(f, ih, ifr, idf, ib);
}

void Hydro::Add_Forces(Forces &to, const Hydro &hy, const Forces &from) {
	if (hy.IsLoadedForce(from)) {
		for (int ib = 0; ib < dt.Nb; ++ib) 
			for (int ihhy = 0; ihhy < hy.dt.Nh; ++ihhy) {
				int ih = dt.FindClosestHead(hy.dt.head[ihhy]);
				for (int ifrhy = 0; ifrhy < hy.dt.Nf; ++ifrhy) {
					int ifr = FindClosest(dt.w, hy.dt.w[ifrhy]);
					for (int idf = 0; idf < 6; ++idf) 	 
						if (IsNum(from[ib][ihhy](ifrhy, idf))) 
							to[ib][ih](ifr, idf) = hy.F_ndim(from, ihhy, ifrhy, idf, ib);
				}
			} 
	}
}

void Hydro::Add_RAO(RAO &to, const Hydro &hy, const RAO &from) {
	if (hy.IsLoadedForce(from)) {
		for (int ib = 0; ib < dt.Nb; ++ib) 
			for (int ihhy = 0; ihhy < hy.dt.Nh; ++ihhy) {
				int ih = dt.FindClosestHead(hy.dt.head[ihhy]);
				for (int ifrhy = 0; ifrhy < hy.dt.Nf; ++ifrhy) {
					int ifr = FindClosest(dt.w, hy.dt.w[ifrhy]);
					for (int idf = 0; idf < 6; ++idf) 	 
						if (IsNum(from[ib][ihhy](ifrhy, idf))) 
							to[ib][ih](ifr, idf) = hy.RAO_ndim(from, ihhy, ifrhy, idf, ib);
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
	dt.Copy(hyd.dt);
}

void Hydro::Data::Copy(const Hydro::Data &hyd) {
	file = hyd.file;
	name = hyd.name;
	g = hyd.g;
    h = hyd.h;
    rho = hyd.rho;
    len = hyd.len;
    dimen = hyd.dimen;
    Nb = hyd.Nb;
    Nf = hyd.Nf;
    Nh = hyd.Nh;

    A = clone(hyd.A);
	Ainf_w = clone(hyd.Ainf_w);
    Ainf = clone(hyd.Ainf);
    A0 = clone(hyd.A0);
    A_P = clone(hyd.A_P);
	
    B = clone(hyd.B);
    B_H = clone(hyd.B_H);
    B_P = clone(hyd.B_P);
    
    head = clone(hyd.head);
	
	x_w = hyd.x_w;
	y_w = hyd.y_w;

    solver = hyd.solver;     
    
    Kirf = clone(hyd.Kirf);
    Tirf = clone(hyd.Tirf);
    
    ex = clone(hyd.ex);
    sc = clone(hyd.sc);
    fk = clone(hyd.fk);
    sc_pot = clone(hyd.sc_pot);
    fk_pot = clone(hyd.fk_pot);
    fk_pot_bmr = clone(hyd.fk_pot_bmr);
    
    rao = clone(hyd.rao);
    
    description = hyd.description;

    sts = clone(hyd.sts);
    dimenSTS = hyd.dimenSTS;
    stsProcessor = hyd.stsProcessor;
    
    qtfsum = clone(hyd.qtfsum);
    qtfdif = clone(hyd.qtfdif);
    qw = clone(hyd.qw);
    qhead = clone(hyd.qhead);
    qtfdataFromW = hyd.qtfdataFromW;
    qtftype = hyd.qtftype;
    
    mdhead = clone(hyd.mdhead);
	md = clone(hyd.md);
	mdtype = hyd.mdtype;
	    
    //T = clone(hyd.T);
    w = clone(hyd.w);
    //dataFromW = hyd.dataFromW;
    //Vo = clone(hyd.Vo); 
    
    msh = clone(hyd.msh);
    
    pots_rad = clone(hyd.pots_rad);
    pots_dif = clone(hyd.pots_dif);
    pots_inc = clone(hyd.pots_inc);
    
    Apan = clone(hyd.Apan);
    
    symX = hyd.symX;
    symY = hyd.symY;
    
    SetId(hyd.GetId());
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

	for (int ib = 0; ib < ret.size(); ++ib) 
		for (int i = 0; i < ret[ib].size(); ++i) 
			for (int j = 0; j < ret[ib][i].cols(); ++j) 
				for (int k = 0; k < ret[ib][i].rows(); ++k) {
					Eigen::VectorXcd r(numT);
					for (int it = 0; it < numT; ++it) 
						r[it] = (*d[it])[ib][i](k, j);
					ret[ib][i](k, j) = r.mean();
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

    dt.w = clone(h0.dt.w);
    dt.head = clone(h0.dt.head);
	dt.qw = clone(h0.dt.qw);
    dt.qhead = clone(h0.dt.qhead);
	dt.mdhead = clone(h0.dt.mdhead);
	dt.msh = clone(h0.dt.msh);
	
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
		if (!CompareDecimals(dt.w, hy.dt.w, 2))
			throw Exc(t_("All models have to have the same frequencies"));
		if (!CompareDecimals(dt.head, hy.dt.head, 2))
			throw Exc(t_("All models have to have the same headings"));
		if (!CompareDecimals(dt.qw, hy.dt.qw, 2))
			throw Exc(t_("All models have to have the same frequencies in QTF"));
		if (!Compare(dt.qhead, hy.dt.qhead))
			throw Exc(t_("All models have to have the same headings in QTF"));
		if (!Compare(dt.mdhead, hy.dt.mdhead))
			throw Exc(t_("All models have to have the same headings in mean drift"));
		for (int ib = 0; ib < dt.Nb; ++ib) {
			if (dt.msh[ib].dt.c0 != hy.dt.msh[ib].dt.c0)
				throw Exc(t_("All models and bodies have to have the same centre of reference"));
		}
	}
	
	dt.file = "Average";
	dt.name = "Average";
	dt.description = "Average";
	dt.solver = BEMROSETTA;
	
    dt.len = 1;
    dt.dimen = true;
    
    /*dt.dataFromW = */dt.qtfdataFromW = true;
    
    dt.qtftype = h0.dt.qtftype;
    dt.mdtype = h0.dt.mdtype;
    
	dt.Ainf.setConstant(dt.Nb*6, dt.Nb*6, NaNDouble);
	dt.A0.setConstant(dt.Nb*6, dt.Nb*6, NaNDouble);
	
	Initialize_AB(dt.A);
	Initialize_AB(dt.B);
	
	for (int ib = 0; ib < dt.Nb; ++ib) 
		dt.msh[ib].dt.C.setConstant(6, 6, 0);
	for (int ib = 0; ib < dt.Nb; ++ib) 
		dt.msh[ib].dt.M.setConstant(6, 6, 0);
	
	Initialize_Forces(dt.ex);
	Initialize_Forces(dt.sc);
	Initialize_Forces(dt.fk);
	Initialize_Forces(dt.rao);
	
	Hydro::Initialize_MD(dt.md, dt.Nb, int(dt.mdhead.size()), dt.Nf);
		
	Hydro::Initialize_QTF(dt.qtfsum, dt.Nb, int(dt.qhead.size()), int(dt.qw.size()));
	Hydro::Initialize_QTF(dt.qtfdif, dt.Nb, int(dt.qhead.size()), int(dt.qw.size()));
			
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
		if (hy.dt.Ainf.size() > 0)
        	Ainfs.Add(&hy.dt.Ainf);
		if (hy.dt.A0.size() > 0)
        	A0s.Add(&hy.dt.A0);
        As.Add(&hy.dt.A);
        Bs.Add(&hy.dt.B);
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
    
	AvgB(dt.Ainf, Ainfs); 
	AvgB(dt.A0, A0s); 
	
	AvgB(dt.A, As); 
	AvgB(dt.B, Bs);
	 
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
		auto Symmetrize_Forces_Each0 = [&](const Forces &f, Forces &newf, double hh, int ih, int idf, int ib, bool applysym) {
			int nih = FindClosest(newHead, hh);
			bool avg = IsNum(newf[ib][nih](0, idf));
			for (int ifr = 0; ifr < dt.Nf; ++ifr) {
				std::complex<double> force = f[ib][ih](ifr, idf);
				if (applysym)
					force = std::polar(abs(force), arg(force) + M_PI);
				std::complex<double> &nf = newf[ib][nih](ifr, idf);
				if (avg)
					nf = Avg(nf, force);
				else 
					nf = force;
			}
		};
	
		Initialize_Forces(newf, newHead.size());
		
		for (int ib = 0; ib < dt.Nb; ++ib) {
			for (int idf = 0; idf < 6; ++idf) {
				for (int ih = 0; ih < dt.Nh; ++ih) {
					Symmetrize_Forces_Each0(f, newf, FixHeading_180(dt.head[ih]), ih, idf, ib, false);
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
					Symmetrize_Forces_Each0(f, newf, he, ih, idf, ib, applysym);
				}
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
	for (int ih = 0; ih < dt.qhead.size(); ++ih) {
		FindAddDelta(newHead, FixHeading_180(dt.qhead[ih]), 0.01);
		FindAddDelta(newHead, FixHeading_180(MirrorHead(dt.qhead[ih], xAxis)), 0.01);
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
			for (int ih = 0; ih < dt.qhead.size(); ++ih) {
				for (int idf = 0; idf < 6; ++idf) {
					Symmetrize_Forces_Each0(f, newf, FixHeading_180(dt.qhead[ih]), ib, ih, idf, false);
					std::complex<double> he = FixHeading_180(MirrorHead(dt.qhead[ih], xAxis));
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
	
	::Copy(newHead, dt.qhead);	
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
		for (int ib = 0; ib < dt.Nb; ++ib) {
			for (int i = 0; i < 6; ++i) {
				double mx = -DBL_MAX, mn = DBL_MAX;
				for (int ifr = 0; ifr < dt.Nf; ifr++) {
					double val = abs(F_ndim(f, ih, ifr, i, ib));
					mx = max(mx, val);
					mn = min(mn, val);
				}
				if (mx != -DBL_MAX && mn != DBL_MAX) {
					double delta = mx - mn;
					double res = 0;
					for (int ifr = 1; ifr < dt.Nf; ifr++) 
						res += abs(F_ndim(f, ih, ifr, i, ib) - F_ndim(f, ih, ifr-1, i, ib));
					res /= delta*(dt.Nf - 1);
					if (res > thres) {
						for (int ifr = 0; ifr < dt.Nf; ifr++) 
							f[ib][ih](ifr, i) = Null;
					}
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
	for (int ib = 0; ib < dt.Nb; ++ib)
		for (int ih = 0; ih < dt.Nh; ++ih) 	
			for (int ifr = 0; ifr < dt.Nf; ++ifr)
				for (int idf = 0; idf < 6; ++idf) {
					std::complex<double> Fa = a[ib][ih](ifr, idf);
					std::complex<double> Fb = b[ib][ih](ifr, idf);
					if (IsLoadedForce(a, idf, ih) && Fa != Fb)
						throw Exc(Format(t_("%s is not the same %f:%f<>%f:%f"), 
								Format(t_("%s[%d][%d](%d, %d)"), type, ib+1, ih+1, ifr+1, idf+1), 
								Fa.real(), Fa.imag(), Fb.real(), Fb.imag()));
				}
}

void Hydro::SaveAs(String fileName, Function <bool(String, int)> Status, BEM_FMT type, int qtfHeading, int ib) {
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
			type = BEMIO_H5;
		else if (ext == ".npz")
			type = AKSELOS_NPZ;
		else
			throw Exc(Format(t_("Conversion to file type '%s' not supported"), fileName));
	}
	BasicBEM::HeadingType htp = Hydro::ShortestHeadingRange(dt.head);
	BasicBEM::HeadingType mtp = Hydro::ShortestHeadingRange(dt.mdhead);
	BasicBEM::HeadingType qtp = Hydro::ShortestHeadingRange(dt.qhead);
	Hydro tosave, *save;
	if (htp != BasicBEM::HEAD_0_360 || mtp != BasicBEM::HEAD_0_360 || qtp != BasicBEM::HEAD_0_360) {
		tosave = clone(*this);
		tosave.SortHeadings(htp, mtp, qtp);
		save = &tosave;
	} else
		save = this;
	
	if (type == WAMIT)
		static_cast<Wamit&>(*save).Save_out(fileName);			
	else if (type == WAMIT_1_3)
		static_cast<Wamit&>(*save).Save(fileName, Status, true, qtfHeading);	
	else if (type == FAST_WAMIT)
		static_cast<Fast&>(*save).Save(fileName, Status, qtfHeading);		
	else if (type == BEMROSETTA)
		save->SaveSerialization(fileName);		
	else if (type == AQWA)
		static_cast<Aqwa&>(*save).Save(fileName, Status);		
	else if (type == CSV_MAT)
		save->SaveCSVMat(fileName);		
	else if (type == CSV_TABLE)
		save->SaveCSVTable(fileName);		
	else if (type == DIODORE)
		save->SaveDiodoreHDB(fileName);		
	else if (type == BEMIO_H5)
		static_cast<BemioH5&>(*save).Save(fileName);
	else if (type == AKSELOS_NPZ)
		SaveAkselos(ib, fileName);
	else
		throw Exc(Format(t_("Conversion to file type '%s' not supported"), fileName));
	
	//solver = type;
	//dt.Nh = realNh;
	//dt.Nf = realNf;
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
	
	if (dt.Nf == 0)
		throw Exc(t_("No frequency found in models"));
	
	dt.msh.SetCount(dt.Nb);
	
	Initialize_AB(dt.A);
	Initialize_AB(dt.B);
	
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
			dt.msh[ib].dt.Vo = hy.dt.msh[ib].dt.Vo;
			dt.msh[ib].dt.cg = clone(hy.dt.msh[ib].dt.cg);
			dt.msh[ib].dt.c0 = clone(hy.dt.msh[ib].dt.c0);
			for (int i = 0; i < 3; ++i) {
				if (IsNum(hy.dt.msh[ib].dt.cb[i]))
					dt.msh[ib].dt.cb[i] = hy.dt.msh[ib].dt.cb[i];
			}
		}
		
		if (hy.IsLoadedC()) {
			for (int ib = 0; ib < dt.Nb; ++ib) {
				for (int idf = 0; idf < 6; ++idf) 
					for (int jdf = 0; jdf < 6; ++jdf) 
						dt.msh[ib].dt.C(idf, jdf) = hy.C_ndim(ib, idf, jdf);
			}
		}
		///////////////////////////////////////////////////////////////////
		
		if (hy.IsLoadedA() && hy.IsLoadedB()) {
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
}

void Hydro::Report() const {
	BEM::Print("\n" + Format(t_("%s file '%s'"), GetCodeStr(), dt.file));
	BEM::Print("\n" + Format(t_("g [m/s2]: %s, h [m]: %s, rho [kg/m³]: %s, length scale [m]: %s"), 
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
			str += S(" ") + t_("vol [m³]") + ": " + FDS(dt.msh[ib].dt.Vo, 8, false);
		if (IsNum(dt.msh[ib].dt.cg))
			str += " " + Format("Cg(%.3f, %.3f, %.3f)[m]", dt.msh[ib].dt.cg.x, dt.msh[ib].dt.cg.y, dt.msh[ib].dt.cg.z);
		if (IsNum(dt.msh[ib].dt.cb))
			str += " " + Format("Cb(%.3f, %.3f, %.3f)[m]", dt.msh[ib].dt.cb.x, dt.msh[ib].dt.cb.y, dt.msh[ib].dt.cb.z);
		if (IsNum(dt.msh[ib].dt.c0))
			str += " " + Format("C0(%.3f, %.3f, %.3f)[m]", dt.msh[ib].dt.c0.x, dt.msh[ib].dt.c0.y, dt.msh[ib].dt.c0.z);
		
		BEM::Print(str);
	}
}

String Hydro::AfterLoad(Function <bool(String, int)> Status) {
	Status(t_("Sorting frequencies and headings"), -1);
	SortFrequencies();
	SortHeadings(BasicBEM::HEAD_0_360, BasicBEM::HEAD_0_360, BasicBEM::HEAD_0_360);
		
	if ((!IsLoadedAinf() || !IsLoadedKirf()) && Bem().calcAinf) {
		Status(t_("Obtaining Ainf, Kirf, and A0"), -1);
		
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
		if (!IsLoadedA0())  
			GetA0();
	}
		
	if (Bem().calcAinf_w) {
		Status(t_("Obtaining Ainf_w"), -1);
		
		if (!IsLoadedKirf())
			GetK_IRF(min(Bem().maxTimeA, GetK_IRF_MaxT()), Bem().numValsA);
		if (!IsLoadedAinf())
			GetAinf();
		if (Status && !Status(t_("Obtaining the frequency-dependent infinite-frequency added mass (A∞(ω))"), 90)) 
			return t_("Cancelled by the user");
		GetAinf_w();
	}
	
	// Fill the other side of the diagonal. If Null, fill with zero
	auto FillNullQTF = [&](UArray<UArray<UArray<MatrixXcd>>> &qtf, bool isSum) {
		for (int ib = 0; ib < dt.Nb; ++ib) {
	        for (int ih = 0; ih < dt.qhead.size(); ++ih) {
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
										cij = std::complex<double>(cji.real(), -cji.imag());
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
	Status(t_("Filling QTF matrices"), -1);
	if (IsLoadedQTF(true)) 
		FillNullQTF(dt.qtfsum, true);
	if (IsLoadedQTF(false))
		FillNullQTF(dt.qtfdif, false);
	
	Status(t_("Postprocessing meshes and symmetrizing potentials"), -1);
	for (int ib = 0; ib < dt.msh.size(); ++ib) {
		Body &m = dt.msh[ib];
		String ret;
		if (!m.dt.mesh.IsEmpty())
			ret = m.dt.mesh.CheckErrors();
		if (!ret.IsEmpty())
			return ret;
		
		//m.dt.fileName = dt.file;
		if (m.dt.name.IsEmpty())
			m.dt.name = InitCaps(GetFileTitle(dt.file));

		// Symmetrize potentials
		auto DeploySymRadiation = [&](int npan, const UVector<double> &signs) {
			if (IsLoadedPotsRad(ib)) {
				dt.pots_rad[ib].SetCount(2*npan);
				for (int ipan = 0; ipan < npan; ++ipan) {
					dt.pots_rad[ib][ipan + npan].SetCount(6);
					for (int idf = 0; idf < 6; ++idf) {
						dt.pots_rad[ib][ipan + npan][idf].SetCount(dt.Nf);
						for (int ifr = 0; ifr < dt.Nf; ++ifr) 
							dt.pots_rad[ib][ipan + npan][idf][ifr] = dt.pots_rad[ib][ipan][idf][ifr]*signs[idf];
					}
				}
			}
		};
		/*auto DeploySymIncDif = [&](int npan, const UVector<double> &signs, UArray<UArray<UArray<UArray<std::complex<double>>>>> &p) {
			if (IsLoadedPotsIncDif(ib, p)) {
				dt.pots_rad[ib].SetCount(2*npan);
				for (int ipan = 0; ipan < npan; ++ipan) {
					dt.pots_rad[ib][ipan + npan].SetCount(dt.Nh);
					for (int ih = 0; ih < dt.Nh; ++ih) {
						dt.pots_rad[ib][ipan + npan][ih].SetCount(dt.Nf);
						for (int ifr = 0; ifr < dt.Nf; ++ifr) 
							p[ib][ipan + npan][ih][ifr] = p[ib][ipan][ih][ifr]*signs[idf];
					}
				}
			}
		};*/
		if (dt.symX) {
			const UVector<double> signs = {-1, 1, 1, 1, -1, -1};
			int npan = m.dt.mesh.panels.size();
			m.dt.mesh.DeployXSymmetry();
			DeploySymRadiation(npan, signs);
			//DeploySymIncDif(npan, signs, dt.pots_inc);
			//DeploySymIncDif(npan, signs, dt.pots_dif);
		}
		if (dt.symY) {
			const UVector<double> signs = {1, -1, 1, -1, 1, -1};
			int npan = m.dt.mesh.panels.size();
			m.dt.mesh.DeployYSymmetry();
			DeploySymRadiation(npan, signs);
			//DeploySymIncDif(npan, signs, dt.pots_inc);
			//DeploySymIncDif(npan, signs, dt.pots_dif);
		}
		if (m.dt.mesh0.IsEmpty()) {
			Point3D cb = m.dt.cb;
			MatrixXd C = m.dt.C;
			m.AfterLoad(dt.rho, dt.g, false, true);
			if (!IsNull(cb))		// restores values from original hydrodynamic file
				m.dt.cb = cb;
			if (C.size() != 0)
				m.dt.C = C;
		}
	}
	if (IsLoadedPotsRad()) {
		Status(t_("Obtaining A and B from potentials"), -1);	
		GetABFromPotentials();
	}  
    if (IsLoadedPotsInc()) {	
        Status(t_("Obtaining Ffk from potentials"), -1);
        GetForcesFromPotentials(dt.pots_inc, dt.fk_pot);
	}
	if (!dt.msh.IsEmpty() && !IsLoadedPotsIncB()) {
		Status(t_("Obtaining incident potentials from mesh"), -1);
		GetPotentialsIncident();
	}	
	if (IsLoadedPotsIncB()) {	
		Status(t_("Obtaining Ffk from bmr potentials"), -1);
		GetForcesFromPotentials(dt.pots_inc_bmr, dt.fk_pot_bmr);
	}
	if (IsLoadedPotsDif()) {
		Status(t_("Obtaining Fsc from potentials"), -1);		
		GetForcesFromPotentials(dt.pots_dif, dt.sc_pot);
	}
	
	CompleteForces1st();

	return String();
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
	
	//int iw0 = GetW0();
	if (dt.w[0] < 0.001) {
		int iw0 = 0;
		dt.A0.setConstant(dt.Nb*6, dt.Nb*6, Null);
		for (int i = 0; i < dt.Nb*6; ++i)
	        for (int j = 0; j < dt.Nb*6; ++j)
				dt.A0(i, j) = dt.A[i][j][iw0];
/*	} else {
		if (dt.Nf == 0 || !IsLoadedKirf())
			return;	
	
		dt.A0.setConstant(dt.Nb*6, dt.Nb*6, NaNDouble);
		
	    for (int i = 0; i < dt.Nb*6; ++i) 
	        for (int j = 0; j < dt.Nb*6; ++j) 
	    		if (IsNum(dt.Kirf[i][j][0]))
			    	dt.A0(i, j) = ::GetA0(dt.Kirf[i][j], dt.Tirf, dt.Ainf(i, j));
	}*/	
	
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
	i %= 6;
	j %= 6;
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
	for (int ib = 0; ib < dt.Nb; ++ib) 
		for (int idf = 0; idf < 6; ++idf) 	
			for (int jdf = 0; jdf < 6; ++jdf) 
				dt.msh[ib].dt.C(idf, jdf) = C_dim(ib, idf, jdf);
}

void Hydro::F_dim(Forces &f) {
	if (f.IsEmpty())
		return;
	for (int ib = 0; ib < dt.Nb; ++ib)
		for (int ih = 0; ih < dt.Nh; ++ih) 	
			for (int ifr = 0; ifr < dt.Nf; ++ifr)
				for (int idf = 0; idf < 6; ++idf) 
					f[ib][ih](ifr, idf) = F_dim(f, ih, ifr, idf, ib);
}

VectorXcd Hydro::F_(bool ndim, const Forces &f, int _h, int ifr, int ib) const {
	VectorXcd ret;
	if (f.IsEmpty())
		return ret;
	ret.resize(6);
	for (int idf = 0; idf < 6; ++idf) 
		ret[idf] = F_(ndim, f, _h, ifr, idf, ib);
	return ret;
}

VectorXcd Hydro::F_dof(bool ndim, const Forces &f, int _h, int idf, int ib) const {
	VectorXcd ret;
	if (f.IsEmpty())
		return ret;
	ret.resize(dt.Nf);
	for (int ifr = 0; ifr < dt.Nf; ++ifr) 
		ret[ifr] = F_(ndim, f, _h, ifr, idf, ib);
	return ret;
}

void Hydro::RAO_dim(RAO &f) {
	if (f.IsEmpty())
		return;
	for (int ib = 0; ib < dt.Nb; ++ib)
		for (int ih = 0; ih < dt.Nh; ++ih) 	
			for (int ifr = 0; ifr < dt.Nf; ++ifr)
				for (int idf = 0; idf < 6; ++idf) 
					f[ib][ih](ifr, idf) = RAO_dim(f, ih, ifr, idf, ib);
}

VectorXcd Hydro::RAO_(bool ndim, const RAO &f, int _h, int ifr, int ib) const {
	VectorXcd ret;
	if (f.IsEmpty())
		return ret;
	ret.resize(6);
	for (int idf = 0; idf < 6; ++idf) 
		ret[idf] = RAO_(ndim, f, _h, ifr, idf, ib);
	return ret;
}

VectorXcd Hydro::RAO_dof(bool ndim, int _h, int idf, int ib) const {
	return RAO_dof(ndim, dt.rao, _h, idf, ib);
}

VectorXcd Hydro::RAO_dof(bool ndim, const RAO &f, int _h, int idf, int ib) const {
	VectorXcd ret;
	if (f.IsEmpty())
		return ret;
	ret.resize(dt.Nf);
	for (int ifr = 0; ifr < dt.Nf; ++ifr) 
		ret[ifr] = RAO_(ndim, f, _h, ifr, idf, ib);
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
	if (!IsNum(dt.msh[ib].dt.Vo) || dt.msh[ib].dt.Vo == 0)
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

void Hydro::GetPotentialsIncident() {
	Initialize_PotsIncDiff(dt.pots_inc_bmr);

	for (int ifr = 0; ifr < dt.Nf; ++ifr) {		
		double k = SeaWaves::WaveNumber_w(dt.w[ifr], dt.h, g_dim());
		double g_w = g_dim()/dt.w[ifr];
		for (int ib = 0; ib < dt.Nb; ++ib) {
			Body &m = dt.msh[ib];
			int npan = m.dt.mesh.panels.size();
			for (int ip = 0; ip < npan; ++ip) {
				const Panel &pan = m.dt.mesh.panels[ip];
				const Point3D &p = pan.centroidPaint;
				for (int ih = 0; ih < dt.Nh; ++ih) {	
					if (p.z < 0) {
						double th = ToRad(dt.head[ih]);
						double cs;
						if (dt.h > 0 && k*dt.h < 700) 
						 	cs = cosh(k*(p.z + dt.h))/cosh(k*dt.h);
						else
							cs = exp(k*p.z);
						double ex = k*(p.x*cos(th) + p.y*sin(th));
													
						dt.pots_inc_bmr[ib][ip][ih][ifr] = g_w*cs*std::complex<double>(sin(ex), cos(ex));// Φ = i g/ω cs (cos(ex)-isin(ex))	Wamit criterion	
					}
				}
			}
		}
	}
}

void Hydro::GetABFromPotentials() {
	Initialize_AB(dt.A_P, 0);
	Initialize_AB(dt.B_P, 0);
	
	dt.Apan = Tensor<double, 5>(dt.Nb, dt.pots_rad[0].size(), 6, 6, dt.Nf);
	
	for (int ib = 0; ib < dt.Nb; ++ib)  {
		const Point3D &c0 = dt.msh[ib].dt.c0;
		for (int ip = 0; ip < dt.pots_rad[ib].size(); ++ip) {
			Value6D n = dt.msh[ib].dt.mesh.panels[ip].NormalExt(c0);	
			for (int ifr = 0; ifr < dt.Nf; ++ifr) {
				for (int idf2 = 0; idf2 < 6; ++idf2) {
					for (int idf1 = 0; idf1 < 6; ++idf1) {
						double A = dt.Apan(ib, ip, idf1, idf2, ifr) = A_pan(ib, ip, idf1, idf2, ifr, n);
						dt.A_P[idf1 + ib*6][idf2 + ib*6][ifr] += A_fromDimFactor(idf1, idf2)*A;
						dt.B_P[idf1 + ib*6][idf2 + ib*6][ifr] += B_fromDimFactor(ifr, idf1, idf2)*B_pan(ib, ip, idf1, idf2, ifr, n);
					}
				}
			}
		}
	}
}

void Hydro::GetForcesFromPotentials(const UArray<UArray<UArray<UArray<std::complex<double>>>>> &pot, Forces &f) {
	Initialize_Forces(f, -1, 0);
	
	for (int ih = 0; ih < dt.Nh; ++ih) {
		for (int ib = 0; ib < dt.Nb; ++ib) {
			const Point3D &c0 = dt.msh[ib].dt.c0;
			for (int ip = 0; ip < pot[ib].size(); ++ip) {
				Value6D n = dt.msh[ib].dt.mesh.panels[ip].NormalExt(c0);	
				for (int idf = 0; idf < 6; ++idf) 
					for (int ifr = 0; ifr < dt.Nf; ++ifr) 
						f[ib][ih](ifr, idf) += F_fromDimFactor(idf)*F_pan(pot, ib, ip, ih, idf, ifr, n);
			}
		}
	}
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
		("code", icode)
		("Kirf", oldKirf)
		("Tirf", dt.Tirf)
		("ex", dt.ex)
		("sc", dt.sc)
		("fk", dt.fk)
		("rao", dt.rao)
		("sts", dt.sts)
		("w", dt.w)
		//("dataFromW", dt.dataFromW)
		("stsProcessor", dt.stsProcessor)
		("dimenSTS", dt.dimenSTS)
		("description", dt.description)
		("qtfsum", dt.qtfsum)
		("qtfdif", dt.qtfdif)
		("msh", dt.msh)
	;
	if(json.IsLoading()) {
		dt.solver = static_cast<Hydro::BEM_FMT>(icode);
		GetOldAB(oldA, dt.A);
		GetOldAB(oldB, dt.B);
		GetOldAB(oldKirf, dt.Kirf);
	}
}