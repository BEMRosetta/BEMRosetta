// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include <ScatterDraw/DataSource.h>
#include <ScatterDraw/Equation.h>
#include "functions.h"

#include <plugin/matio/lib/matio.h>

using namespace Upp;

Function <void(String)> BEM::Print 		  = [](String s) {Cout() << s;};
Function <void(String)> BEM::PrintWarning = [](String s) {Cout() << s;};
Function <void(String)> BEM::PrintError   = [](String s) {Cout() << s;};

const char *BEM::strDOFtext[] 	 = {t_("surge"), t_("sway"), t_("heave"), t_("roll"), t_("pitch"), t_("yaw")};
const char *BEM::strDOFtextAbrev[] = {t_("s"), t_("w"), t_("h"), t_("r"), t_("p"), t_("y")};
const char *BEM::strDOFnum[] 	 = {t_("1"), t_("2"), t_("3"), t_("4"), t_("5"), t_("6")};
const char *BEM::strDOFxyz[] 	 = {t_("x"), t_("y"), t_("z"), t_("rx"), t_("ry"), t_("rz")};

const char *Hydro::strDataToPlot[] = {t_("A(ω)"), t_("A∞"), t_("A0"), t_("B(ω)"), t_("A∞(ω)"), t_("Kirf"),
				t_("Fsc_ma"), t_("Fsc_ph"), t_("Ffk_ma"), t_("Ffk_ph"), t_("Fex_ma"), t_("Fex_ph"),
				t_("RAO_ma"), t_("RAO_ph"), t_("Z_ma"), t_("Z_ph"), t_("Kr_ma"), t_("Kr_ph"), 
				t_("TFS_ma"), t_("TFS_ph")};

const char *BEM::strDOFType[] = {t_("1,2,3,4,5,6"), t_("surge,sway,"), t_("x,y,z,rx,ry,rz"), ""};
BEM::DOFType BEM::dofType = BEM::DOFSurgeSway;

const char *BEM::strHeadingType[] = {t_("-180->180º"), t_("0->360º"), ""};
BEM::HeadingType BEM::headingType = BEM::HEAD_180_180;
	
const char *BEMCase::solverStr[] = {t_("Nemoh"), t_("Nemoh v115"), t_("Capytaine"), t_("HAMS"), t_("AQWA")};

int Hydro::idCount = 0;	

bool PrintStatus(String s, int) {
	Cout() << "\n" << RemoveAccents(s);
	return true;
};

void Hydro::InitializeSts() {
	sts.SetCount(6*Nb);
	for (int ib = 0; ib < 6*Nb; ++ib) 
		sts[ib].SetCount(6*Nb);
}

void Hydro::GetFexFromFscFfk() {
	for (int ih = 0; ih < Nh; ++ih) {
		for (int ifr = 0; ifr < Nf; ++ifr) {
			for (int i = 0; i < Nb*6; ++i) {
				if (!IsNull(sc.ma[ih](ifr, i))) {
					double exre = sc.re[ih](ifr, i) + fk.re[ih](ifr, i);
					double exim = sc.im[ih](ifr, i) + fk.im[ih](ifr, i);
					ex.re[ih](ifr, i) = exre;
					ex.im[ih](ifr, i) = exim;
					ex.ma[ih](ifr, i) = sqrt(exre*exre + exim*exim);
					ex.ph[ih](ifr, i) = atan2(exim, exre);
				}
			}
		}
	}
}

void Hydro::Initialize_RAO() {
	Initialize_Forces(rao);
}

void Hydro::Normalize() {
	if (IsLoadedC()) {
		for (int ib = 0; ib < Nb; ++ib) {
			for (int idf = 0; idf < 6; ++idf) 
				for (int jdf = 0; jdf < 6; ++jdf) 
					C[ib](idf, jdf) = C_ndim(ib, idf, jdf);
		}
	}
	if (IsLoadedA() && IsLoadedB()) {
		for (int ifr = 0; ifr < Nf; ++ifr) {
			for (int idf = 0; idf < 6*Nb; ++idf) {
				for (int jdf = 0; jdf < 6*Nb; ++jdf) {	
					A[idf][jdf][ifr] = A_ndim(ifr, idf, jdf);
					B[idf][jdf][ifr] = B_ndim(ifr, idf, jdf);
				}
			}
		}
	}
	if (IsLoadedAinf()) {
		for (int i = 0; i < 6*Nb; ++i) 
			for (int j = 0; j < 6*Nb; ++j) 
				Ainf(i, j) = Ainf_ndim(i, j);
	}
	if (IsLoadedA0()) {
		for (int i = 0; i < 6*Nb; ++i) 
			for (int j = 0; j < 6*Nb; ++j) 
				A0(i, j) = A0_ndim(i, j);
	}
	if (IsLoadedFex())
    	Normalize_Forces(ex);
	if (IsLoadedFsc())
		Normalize_Forces(sc);
	if (IsLoadedFfk())
		Normalize_Forces(fk);
	if (IsLoadedRAO()) {
		for (int h = 0; h < Nh; ++h) {
			for (int ifr = 0; ifr < Nf; ++ifr) {
				for (int i = 0; i < 6*Nb; ++i) {	 
					rao.ma[h](ifr, i) = R_ma_ndim(rao, h, ifr, i);
					rao.re[h](ifr, i) = R_re_ndim(rao, h, ifr, i);
					rao.im[h](ifr, i) = R_im_ndim(rao, h, ifr, i);
				}
			}
		}
	}
}

void Hydro::Dimensionalize() {
	if (IsLoadedC()) {
		for (int ib = 0; ib < Nb; ++ib) {
			for (int idf = 0; idf < 6; ++idf) 
				for (int jdf = 0; jdf < 6; ++jdf) 
					C[ib](idf, jdf) = C_dim(ib, idf, jdf);
		}
	}
	if (IsLoadedA() && IsLoadedB()) {
		for (int ifr = 0; ifr < Nf; ++ifr) {
			for (int idf = 0; idf < 6*Nb; ++idf) {
				for (int jdf = 0; jdf < 6*Nb; ++jdf) {	
					A[idf][jdf][ifr] = A_dim(ifr, idf, jdf);
					B[idf][jdf][ifr] = B_dim(ifr, idf, jdf);
				}
			}
		}
	}
	if (IsLoadedAinf()) {	
		for (int i = 0; i < 6*Nb; ++i) 
			for (int j = 0; j < 6*Nb; ++j) 
				Ainf(i, j) = Ainf_dim(i, j);
	}
	if (IsLoadedA0()) {	
		for (int i = 0; i < 6*Nb; ++i) 
			for (int j = 0; j < 6*Nb; ++j) 
				A0(i, j) = A0_ndim(i, j);
	}
	if (IsLoadedFex())
    	Dimensionalize_Forces(ex);
	if (IsLoadedFsc())
		Dimensionalize_Forces(sc);
	if (IsLoadedFfk())
		Dimensionalize_Forces(fk);
	if (IsLoadedRAO()) {
		for (int h = 0; h < Nh; ++h) {
			for (int ifr = 0; ifr < Nf; ++ifr) {
				for (int i = 0; i < 6*Nb; ++i) {	 
					rao.ma[h](ifr, i) = R_ma_dim(rao, h, ifr, i);
					rao.re[h](ifr, i) = R_re_dim(rao, h, ifr, i);
					rao.im[h](ifr, i) = R_im_dim(rao, h, ifr, i);
				}
			}
		}
	}
}

void Hydro::Initialize_Forces() {
	Initialize_Forces(ex);
	Initialize_Forces(sc);
	Initialize_Forces(fk);
}

void Hydro::Initialize_Forces(Forces &f, int _Nh) {
	if (_Nh == -1)
		_Nh = Nh;
	f.ma.SetCount(_Nh);
	f.ph.SetCount(_Nh);
	f.re.SetCount(_Nh);
	f.im.SetCount(_Nh);
	for (int ih = 0; ih < _Nh; ++ih) {
		f.ma[ih].setConstant(Nf, Nb*6, Null);
		f.ph[ih].setConstant(Nf, Nb*6, Null);
		f.re[ih].setConstant(Nf, Nb*6, Null);
		f.im[ih].setConstant(Nf, Nb*6, Null);
	}
}

void Hydro::GetMaPh(Forces &f) {
	for (int ih = 0; ih < Nh; ++ih) {
		for (int ifr = 0; ifr < Nf; ++ifr) {
			for (int idf = 0; idf < 6*Nb; ++idf) {	 
				f.ma[ih](ifr, idf) = sqrt(sqr(f.re[ih](ifr, idf)) + sqr(f.im[ih](ifr, idf)));
				f.ph[ih](ifr, idf) = atan2(f.im[ih](ifr, idf), f.re[ih](ifr, idf));
			}
		}
	} 
}
    	
void Hydro::Normalize_Forces(Forces &f) {
	for (int ih = 0; ih < Nh; ++ih) {
		for (int ifr = 0; ifr < Nf; ++ifr) {
			for (int idf = 0; idf < 6*Nb; ++idf) {	 
				f.ma[ih](ifr, idf) = F_ma_dim(f, ih, ifr, idf);
				f.re[ih](ifr, idf) = F_re_dim(f, ih, ifr, idf);
				f.im[ih](ifr, idf) = F_im_dim(f, ih, ifr, idf);
			}
		}
	} 
}

void Hydro::Dimensionalize_Forces(Forces &f) {
	for (int ih = 0; ih < Nh; ++ih) {
		for (int ifr = 0; ifr < Nf; ++ifr) {
			for (int idf = 0; idf < 6*Nb; ++idf) {	 
				f.ma[ih](ifr, idf) = F_ma_dim(f, ih, ifr, idf);
				f.re[ih](ifr, idf) = F_re_dim(f, ih, ifr, idf);
				f.im[ih](ifr, idf) = F_im_dim(f, ih, ifr, idf);
			}
		}
	}
}

void Hydro::Add_Forces(Forces &to, const Hydro &hydro, const Forces &from) {
	if (hydro.IsLoadedForce(from)) {
		for (int ihhy = 0; ihhy < hydro.Nh; ++ihhy) {
			int ih = FindClosest(head, hydro.head[ihhy]);
			for (int ifrhy = 0; ifrhy < hydro.Nf; ++ifrhy) {
				int ifr = FindClosest(w, hydro.w[ifrhy]);
				for (int idf = 0; idf < 6*Nb; ++idf) {	 
					if (!IsNull(from.ma[ihhy](ifrhy, idf))) {
						to.ma[ih](ifr, idf) = hydro.F_ma_ndim(from, ihhy, ifrhy, idf);
						to.ph[ih](ifr, idf) = from.ph[ihhy](ifrhy, idf); 
						to.re[ih](ifr, idf) = hydro.F_re_ndim(from, ihhy, ifrhy, idf);
						to.im[ih](ifr, idf) = hydro.F_im_ndim(from, ihhy, ifrhy, idf);
					}
				}
			}
		} 
	}
}

void Hydro::Symmetrize_Forces_Each0(const Forces &f, Forces &newf, const Upp::Vector<double> &newHead, double h, int ih, int idb) {
	int nih  = FindClosest(newHead, h);
	bool avg  = !IsNull(newf.re[nih](0, idb));
	for (int ifr = 0; ifr < Nf; ++ifr) {
		if (avg) {
			double re = newf.re[nih](ifr, idb) = Avg(newf.re[nih](ifr, idb), f.re[ih](ifr, idb));
			double im = newf.im[nih](ifr, idb) = Avg(newf.im[nih](ifr, idb), f.im[ih](ifr, idb));			
			newf.ma[nih](ifr, idb) = sqrt(re*re + im*im);
			newf.ph[nih](ifr, idb) = atan2(im, re);
		} else {
			newf.ma[nih](ifr, idb) = f.ma[ih](ifr, idb);
			newf.ph[nih](ifr, idb) = f.ph[ih](ifr, idb); 
			newf.re[nih](ifr, idb) = f.re[ih](ifr, idb);
			newf.im[nih](ifr, idb) = f.im[ih](ifr, idb);
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

void Hydro::Symmetrize_ForcesEach(const Forces &f, Forces &newf, const Upp::Vector<double> &newHead, int newNh, bool xAxis) {
	Initialize_Forces(newf, newNh);
	
	for (int idb = 0; idb < 6*Nb; ++idb) {
		for (int ih = 0; ih < Nh; ++ih) {
			Symmetrize_Forces_Each0(f, newf, newHead, FixHeading_180(head[ih]), ih, idb);
			Symmetrize_Forces_Each0(f, newf, newHead, FixHeading_180(MirrorHead(head[ih], xAxis)), ih, idb);
		}
	}
}

void Hydro::Copy(const Hydro &hyd) {
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
	
	Dlin = clone(Dlin);
	
    B = clone(hyd.B);
    
    head = clone(hyd.head);
    names = clone(hyd.names);
    C = clone(hyd.C);
    M = clone(hyd.M);
    cb = clone(hyd.cb);
    cg = clone(hyd.cg);
    c0 = clone(hyd.c0);
    code = hyd.code;     
    dof = clone(hyd.dof); 
    dofOrder = clone(hyd.dofOrder);
    
    Kirf = clone(hyd.Kirf);
    Tirf = clone(hyd.Tirf);
    
    ex = clone(hyd.ex);
    sc = clone(hyd.sc);
    fk = clone(hyd.fk);
    rao = clone(hyd.rao);
    
    description = hyd.description;

    sts = clone(hyd.sts);
    dimenSTS = hyd.dimenSTS;
    stsProcessor = hyd.stsProcessor;
    
    qtfsum = clone(hyd.qtfsum);
    qtfdif = clone(hyd.qtfdif);
    qtfw = clone(hyd.qtfw);
    qtfT = clone(hyd.qtfT);
    qtfhead = clone(hyd.qtfhead);
    qtfdataFromW = hyd.qtfdataFromW;
    
    qtfCases = clone(hyd.qtfCases);
    
    T = clone(hyd.T);
    w = clone(hyd.w);
    dataFromW = hyd.dataFromW;
    Vo = clone(hyd.Vo); 
    
    bem = hyd.bem;
		
	//id = hyd.id;
}

void Hydro::Symmetrize_Forces(bool xAxis) {
	if (!IsLoadedFex() && !IsLoadedFsc() && !IsLoadedFfk() && !IsLoadedRAO())
		return;
	
	Upp::Vector<double> newHead;
	for (int ih = 0; ih < Nh; ++ih) {
		FindAddRatio(newHead, FixHeading_180(head[ih]), 0.001);
		FindAddRatio(newHead, FixHeading_180(MirrorHead(head[ih], xAxis)), 0.001);
	}
	Sort(newHead);
	int newNh = newHead.size();
	
	Forces newex, newsc, newfk;
	RAO newrao;
	
	if (IsLoadedFex()) {
		Symmetrize_ForcesEach(ex, newex, newHead, newNh, xAxis);
		ex = pick(newex);
	}
	if (IsLoadedFsc()) {
		Symmetrize_ForcesEach(sc, newsc, newHead, newNh, xAxis);
		sc = pick(newsc);
	}
	if (IsLoadedFfk()) {
		Symmetrize_ForcesEach(fk, newfk, newHead, newNh, xAxis);
		fk = pick(newfk);
	}
	if (IsLoadedRAO()) {
		Symmetrize_ForcesEach(rao, newrao, newHead, newNh, xAxis);
		rao = pick(newrao);
	}
	head = pick(newHead);		// New headings are set between -180 and 180
	Nh = newNh;
}

void Hydro::RemoveThresDOF_A(double thres) {
	if (!IsLoadedA())
		return;
	for (int idf = 0; idf < 6*Nb; ++idf) {
		for (int jdf = 0; jdf < 6*Nb; ++jdf) {
			double mx = -DBL_MAX, mn = DBL_MAX;
			for (int ifr = 0; ifr < Nf; ifr++) {
				double val = A_ndim(ifr, idf, jdf);
				mx = max(mx, val);
				mn = min(mn, val);
			}
			double delta = mx - mn;
			if (!IsNull(mx) && !IsNull(mn)) {
				double res = 0;
				for (int ifr = 1; ifr < Nf; ifr++) 
					res += abs(A_ndim(ifr, idf, jdf) - A_ndim(ifr-1, idf, jdf));
				res /= delta*(Nf - 1);
				if (res > thres) {
					for (int ifr = 0; ifr < Nf; ifr++) 
						A[idf][jdf][ifr] = Null;
					Ainf(idf, jdf) = Null;		
				}
			}
		}
	}
}

void Hydro::RemoveThresDOF_B(double thres) {
	if (!IsLoadedB())
		return;
	for (int idf = 0; idf < 6*Nb; ++idf) {
		for (int jdf = 0; jdf < 6*Nb; ++jdf) {
			double mx = -DBL_MAX, mn = DBL_MAX;
			for (int ifr = 0; ifr < Nf; ifr++) {
				double val = B_ndim(ifr, idf, jdf);
				mx = max(mx, val);
				mn = min(mn, val);
			}
			double delta = mx - mn;
			if (!IsNull(mx) && !IsNull(mn)) {
				double res = 0;
				for (int ifr = 1; ifr < Nf; ifr++) 
					res += abs(B_ndim(ifr, idf, jdf) - B_ndim(ifr-1, idf, jdf));
				res /= delta*(Nf - 1);
				if (res > thres) {
					for (int ifr = 0; ifr < Nf; ifr++) 
						B[idf][jdf][ifr] = Null;
				}
			}
		}
	}
}

void Hydro::RemoveThresDOF_Force(Forces &f, double thres) {
	if (!IsLoadedForce(f))
		return;
	for (int h = 0; h < Nh; ++h) {
		for (int i = 0; i < 6*Nb; ++i) {
			double mx = -DBL_MAX, mn = DBL_MAX;
			for (int ifr = 0; ifr < Nf; ifr++) {
				double val = F_ma_ndim(f, h, ifr, i);
				mx = max(mx, val);
				mn = min(mn, val);
			}
			double delta = mx - mn;
			if (!IsNull(mx) && !IsNull(mn)) {
				double res = 0;
				for (int ifr = 1; ifr < Nf; ifr++) 
					res += abs(F_ma_ndim(f, h, ifr, i) - F_ma_ndim(f, h, ifr-1, i));
				res /= delta*(Nf - 1);
				if (res > thres) {
					for (int ifr = 0; ifr < Nf; ifr++) 
						f.ma[h](ifr, i) = Null;
				}
			}
		}
	}
}

void Hydro::Compare_rho(Hydro &a) {
	if (a.rho != rho)
		throw Exc(Format(t_("%s is not the same %f<>%f"), t_("Density rho"), a.rho, rho));
}

void Hydro::Compare_g(Hydro &a) {
	if (a.g != g)
		throw Exc(Format(t_("%s is not the same %f<>%f"), t_("Gravity g"), a.g, g));
}

void Hydro::Compare_h(Hydro &a) {
	if (a.h != h)
		throw Exc(Format(t_("%s is not the same %f<>%f"), t_("Water depth h"), a.h, h));
}

void Hydro::Compare_Nb(Hydro &a) {
	if (a.Nb != Nb)
		throw Exc(Format(t_("%s is not the same %d<>%d"), t_("Number of bodies"), a.Nb, Nb));
}

void Hydro::Compare_w(Hydro &a) {
	if (a.Nf != Nf)	
		throw Exc(Format(t_("%s is not the same %f<>%f"), t_("Number of frequencies"), a.Nf, Nf));
	for (int i = 0; i < a.Nf; ++i) {
		if (!EqualRatio(a.w[i], w[i], 0.0001))
			throw Exc(Format(t_("%s is not the same %f<>%f"), 
							Format(t_("#%d %s"), i+1, t_("frequency")), a.w[i], w[i]));
	}
}

void Hydro::Compare_head(Hydro &a) {
	if (a.Nh != Nh)	
		throw Exc(Format(t_("%s is not the same %f<>%f"), t_("Number of headings"), a.Nh, Nh));
	for (int i = 0; i < a.Nh; ++i) {
		if (a.head[i] != head[i])
			throw Exc(Format(t_("%s is not the same %f<>%f"), 
							Format(t_("#%d %s"), i+1, t_("frequency")), a.w[i], w[i]));
	}
}

void Hydro::Compare_A(Hydro &a) {
	for (int ifr = 0; ifr < a.Nf; ifr++) {
		for (int idf = 0; idf < 6*a.Nb; ++idf) {
			for (int jdf = 0; jdf < 6*a.Nb; ++jdf) {
				double Aa = a.A[idf][jdf][ifr];
				double Ab = A[idf][jdf][ifr];
				if (!IsNull(Aa) && !IsNull(Ab) && Aa != Ab)
					throw Exc(Format(t_("%s is not the same %f<>%f"), 
							Format(t_("%s[%d](%d, %d)"), t_("A"), ifr+1, idf+1, jdf+1), 
							Aa, Ab));
			}
		}
	}
}

void Hydro::Compare_B(Hydro &a) {
	for (int ifr = 0; ifr < a.Nf; ifr++) {
		for (int idf = 0; idf < 6*a.Nb; ++idf) {
			for (int jdf = 0; jdf < 6*a.Nb; ++jdf) {
				double Ba = a.B[idf][jdf][ifr];
				double Bb = B[idf][jdf][ifr];
				if (!IsNull(Ba) && !IsNull(Bb) && Ba != Bb)
					throw Exc(Format(t_("%s is not the same %f<>%f"), 
							Format(t_("%s[%d](%d, %d)"), t_("B"), ifr+1, idf+1, jdf+1), 
							Ba, Bb));
			}
		}
	}
}

void Hydro::Compare_C(Hydro &a) {
	for (int ib = 0; ib < a.Nb; ib++) {
		for (int idf = 0; idf < 6; ++idf) {
			for (int jdf = 0; jdf < 6; ++jdf) {
				double Ca = a.C[ib](idf, jdf);
				double Cb = C[ib](idf, jdf);
				if (!IsNull(Ca) && !IsNull(Cb) && !EqualRatio(Ca, Cb, 0.0001))
					throw Exc(Format(t_("%s is not the same %f<>%f"), 
							Format(t_("%s[%d](%d, %d)"), t_("C"), ib+1, idf+1, jdf+1), 
							Ca, Cb));
			}
		}
	}
}

void Hydro::Compare_cg(Hydro &a) {
	for (int i = 0; i < 3; i++) {
		for (int ib = 0; ib < a.Nb; ib++) {
			if (a.cg(i, ib) != cg(i, ib))
				throw Exc(Format(t_("%s is not the same %f<>%f"), 
						Format(t_("%s(%d, %d)"), t_("cg"), i+1, ib+1), 
							a.cg(i, ib), cg(i, ib)));
		}
	}
}

bool Hydro::SaveAs(String file, Function <bool(String, int)> Status, BEM_SOFT type, int qtfHeading) {
	int realNh = Nh;
	int realNf = Nf;
	
	if (type == UNKNOWN) {
		String ext = ToLower(GetFileExt(file));
		
		if (ext == ".1" || ext == ".2" || ext == ".3" || ext == ".hst" || ext == ".4" || ext == ".12s" || ext == ".12d") 
			type = Hydro::WAMIT_1_3;
		else if (ext == ".dat")
			type = Hydro::FAST_WAMIT;	
		else if (ext == ".bem")
			type = Hydro::BEMROSETTA;
		else
			throw Exc(Format(t_("Conversion to file type '%s' not supported"), file));
	}
	bool ret = false;
	if (type == WAMIT) {
		Wamit data(*bem, this);
		ret = data.Save_out(file, bem->g, bem->rho);			
	} else if (type == WAMIT_1_3) {
		Wamit data(*bem, this);
		ret = data.Save(file, Status, true, qtfHeading);	
	} else if (type == FAST_WAMIT) {
		Fast data(*bem, this);
		ret = data.Save(file, Status, qtfHeading);		
	} else if (type == BEMROSETTA) {
		HydroClass data(*bem, this);
		ret = data.Save(file);		
	}
	code = type;
	Nh = realNh;
	Nf = realNf;
	
	return ret;
}

void Hydro::Join(const Upp::Vector<Hydro *> &hydrosp) {
	name = t_("Joined files");
	for (int ihy = 0; ihy < hydrosp.size(); ++ihy) {
		const Hydro &hydro = *hydrosp[ihy];
		if (hydro.name.Find("Nemoh_Part") >= 0) {
			name = GetFileTitle(GetFileFolder(GetFileFolder(hydro.file)));
			break;
		}
	}
	
	dimen = false;
	g = bem->g;
	rho = bem->rho;
	len = 1;

	h = Null;
	for (int ihy = 0; ihy < hydrosp.size(); ++ihy) {
		const Hydro &hydro = *hydrosp[ihy];
		if (!IsNull(hydro.h))
			h = hydro.h;
		else if (h != hydro.h)
			throw Exc(Format(t_("Water depth does not match between '%s'(%d) and '%s'(%d)"), 
					hydrosp[0]->name, hydrosp[0]->h, hydro.name, hydro.h));			
	}	
	if (IsNull(h))
		throw Exc(t_("No water depth found in models"));
			
	Nb = Null;
	for (int ihy = 0; ihy < hydrosp.size(); ++ihy) {
		const Hydro &hydro = *hydrosp[ihy];
		if (IsNull(Nb))
			Nb = hydro.Nb;
		else if (Nb != hydro.Nb)
			throw Exc(Format(t_("Number of bodies does not match between '%s'(%d) and '%s'(%d)"), 
					hydrosp[0]->name, hydrosp[0]->Nb, hydro.name, hydro.Nb));			
	}
	if (IsNull(Nb))
		throw Exc(t_("No body found in models"));
		
	head.Clear();
	for (int ihy = 0; ihy < hydrosp.size(); ++ihy) {
		const Hydro &hydro = *hydrosp[ihy];
		if (hydro.IsLoadedFex() || hydro.IsLoadedFsc() || hydro.IsLoadedFfk()) {
			for (int ih = 0; ih < hydro.head.size(); ih++) {
				double head_v = hydro.head[ih];
				FindAddRatio(head, head_v, 0.001);
			}
		}
	}
	Nh = head.size();
	if (Nh == 0)
		throw Exc(t_("No head found in models"));
	
	w.Clear();
	for (int ihy = 0; ihy < hydrosp.size(); ++ihy) {
		const Hydro &hydro = *hydrosp[ihy];
		if (hydro.IsLoadedA() && hydro.IsLoadedB()) {
			for (int ifr = 0; ifr < hydro.w.size(); ifr++) {
				double w_v = hydro.w[ifr];
				FindAddRatio(w, w_v, 0.001);
			}
		}
	}
	Sort(w);
	Nf = w.size();
	T.Clear();
	for (int i = 0; i < Nf; ++i)
		T << 2*M_PI/w[i];
	
	if (Nf == 0)
		throw Exc(t_("No frequency found in models"));
	
	names.SetCount(Nb);
	dof.SetCount(Nb);
	cg.setConstant(3, Nb, Null);
	c0.setConstant(3, Nb, Null);
	cb.setConstant(3, Nb, Null);
	Vo.SetCount(Nb, Null);
	
	A.SetCount(6*Nb);
	B.SetCount(6*Nb);
	for (int i = 0; i < 6*Nb; ++i) {
		A[i].SetCount(6*Nb);
		B[i].SetCount(6*Nb);
		for (int j = 0; j < 6*Nb; ++j) {
			A[i][j].setConstant(Nf, Null);	
			B[i][j].setConstant(Nf, Null);	
		}
	}
	
	C.SetCount(Nb);
	for (int ib = 0; ib < Nb; ++ib) 
		C[ib].setConstant(6, 6, Null);
	
	Initialize_Forces(ex);
	Initialize_Forces(sc);
	Initialize_Forces(fk);
		
	for (int ihy = 0; ihy < hydrosp.size(); ++ihy) {
		const Hydro &hydro = *hydrosp[ihy];
		
		// All this block should have to be the same. Now it is not tested
		
		code = hydro.code;
		
		for (int ib = 0; ib < Nb; ++ib) {
			if (!hydro.names[ib].IsEmpty())
				names[ib] = hydro.names[ib];
			if (!IsNull(hydro.Vo[ib]))
				Vo[ib] = hydro.Vo[ib];
			for (int i = 0; i < 3; ++i) {
				if (!IsNull(hydro.cg(i, ib)))
					cg(i, ib) = hydro.cg(i, ib);
				if (!IsNull(hydro.cb(i, ib)))
					cb(i, ib) = hydro.cb(i, ib);
				if (!IsNull(hydro.c0(i, ib)))
					c0(i, ib) = hydro.c0(i, ib);
			}
			dof[ib] = hydro.dof[ib];
		}
		
		if (IsLoadedC() && hydro.IsLoadedC()) {
			for (int ib = 0; ib < Nb; ++ib) {
				for (int idf = 0; idf < 6; ++idf) 
					for (int jdf = 0; jdf < 6; ++jdf) 
						C[ib](idf, jdf) = hydro.C_ndim(ib, idf, jdf);
			}
		}
		///////////////////////////////////////////////////////////////////
		
		if (IsLoadedA() && IsLoadedB() && hydro.IsLoadedA() && hydro.IsLoadedB()) {
			for (int ifrhy = 0; ifrhy < hydro.Nf; ++ifrhy) {
				int ifr = FindClosest(w, hydro.w[ifrhy]);
				for (int idf = 0; idf < 6*Nb; ++idf) {
					for (int jdf = 0; jdf < 6*Nb; ++jdf) {	
						if (!IsNull(hydro.A[idf][jdf][ifrhy]))
							A[idf][jdf][ifr] = hydro.A_ndim(ifrhy, idf, jdf);
						if (!IsNull(hydro.B[idf][jdf][ifrhy]))
							B[idf][jdf][ifr] = hydro.B_ndim(ifrhy, idf, jdf);
					}
				}
			}
		}	
		Add_Forces(ex, hydro, hydro.ex);
		Add_Forces(sc, hydro, hydro.sc);
		Add_Forces(fk, hydro, hydro.fk);
		
		if (IsLoadedRAO() && hydro.IsLoadedRAO()) {
			for (int ihhy = 0; ihhy < Nh; ++ihhy) {
				int ih = FindClosest(head, hydro.head[ihhy]);
				for (int ifrhy = 0; ifrhy < Nf; ++ifrhy) {
					int ifr = FindClosest(w, hydro.w[ifrhy]);
					for (int idf = 0; idf < 6*Nb; ++idf) {	 
						if (!IsNull(rao.ma[ihhy](ifrhy, idf))) {
							rao.ma[ih](ifr, idf) = hydro.R_ma_ndim(rao, ihhy, ifrhy, idf);
							rao.ph[ih](ifr, idf) = rao.ph[ihhy](ifrhy, idf); 
							rao.re[ih](ifr, idf) = hydro.R_re_ndim(rao, ihhy, ifrhy, idf);
							rao.im[ih](ifr, idf) = hydro.R_im_ndim(rao, ihhy, ifrhy, idf);
						}
					}
				}
			}
		}
	}
	
	// Aw0 has to be recalculated
	/*
	// A0 is set from the lower frequency data set
	
	int ihminw = -1;
	double minw = DBL_MAX;
	for (int ihy = 0; ihy < hydrosp.size(); ++ihy) {
		const Hydro &hydro = *hydrosp[ihy];
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
	
	bem->calcAinf = true;
	bem->calcAinf_w = true;
	
	// sts has to be recalculated	
}

void Hydro::Report() const {
	BEM::Print("\n" + Format(t_("%s file '%s'"), GetCodeStr(), file));
	BEM::Print("\n" + Format(t_("g [m/s2]: %s, h [m]: %s, rho [kg/m3]: %s, length scale [m]: %s"), 
								S_g(), S_h(), S_rho(), S_len()));
	String freqs;
	if (w.IsEmpty()) 
		freqs = t_("NONE");
	else if (w.size() > 1) {
		String strDeltaH;
		if (GetIrregularFreq() < 0) 
			strDeltaH = Format(t_("delta %s [rad/s]"), FormatDoubleSize(w[1] - w[0], 8, false));
		else {
			String strHead;
			for (int i = 0; i < w.size(); ++i) {
				if (i > 0)
					strHead << ", ";
				strHead << w[i];
			}
			strDeltaH = Format(t_("Non constant delta (%s)"), strHead); 
		}
	 	freqs = Format(t_("%s to %s %s"), FormatDoubleSize(w[0], 8, false), 
	 									  FormatDoubleSize(w[w.size()-1], 8, false), strDeltaH);	
	} else
		freqs = Format(t_("%s [rad/s]"), FormatDoubleSize(w[0], 8, false));
	
	String heads;
	if (head.IsEmpty())
		heads = t_("NONE");
	else if (head.size() > 1) {
		String strDeltaH;
		if (GetIrregularHead() < 0) 
			strDeltaH = Format(t_("delta %.1f [º]"), head[1] - head[0]);
		else {
			String strHead;
			for (int i = 0; i < head.size(); ++i) {
				if (i > 0)
					strHead << ", ";
				strHead << head[i];
			}
			strDeltaH = Format(t_("Non constant delta (%s)"), strHead); 
		}
	 	heads = Format(t_("%.1f to %.1f %s"), head[0], head[head.size()-1], strDeltaH);	
	} else
		heads = Format(t_("%.1f [º]"), head[0]);
	
	BEM::Print("\n" + Format(t_("#freqs: %d (%s)"), Nf, freqs)); 
	BEM::Print("\n" + Format(t_("#headings: %d (%s)"), Nh, heads)); 
	BEM::Print("\n" + Format(t_("#bodies: %d"), Nb));
	for (int ib = 0; ib < Nb; ++ib) {
		String str = Format("\n%d.", ib+1);
		if (names.size() > ib)
			str += " '" + names[ib] + "'";
		if (dof.size() > ib)
			str += S(" ") + t_("dof") + ": " + FormatInt(dof[ib]);
		if (Vo.size() > ib && !IsNull(Vo[ib]))
			str += S(" ") + t_("vol [m3]") + ": " + FormatDoubleSize(Vo[ib], 8, false);
		if (cg.size() > 3*ib && !IsNull(cg(0, ib)))
			str += " " + Format("Cg(%.3f, %.3f, %.3f)[m]", cg(0, ib), cg(1, ib), cg(2, ib));
		if (cb.size() > 3*ib && !IsNull(cb(0, ib)))
			str += " " + Format("Cb(%.3f, %.3f, %.3f)[m]", cb(0, ib), cb(1, ib), cb(2, ib));
		if (c0.size() > 3*ib && !IsNull(c0(0, ib)))
			str += " " + Format("C0(%.3f, %.3f, %.3f)[m]", c0(0, ib), c0(1, ib), c0(2, ib));
		
		BEM::Print(str);
	}
}

void Hydro::GetBodyDOF() {
	dof.Clear();	 dof.SetCount(Nb, 0);
	for (int ib = 0; ib < Nb; ++ib)
		for (int idf = 0; idf < 6; ++idf)
			if (IsAvailableDOF(ib, idf))
				dof[ib]++;
}

bool Hydro::AfterLoad(Function <bool(String, int)> Status) {
	if (dofOrder.IsEmpty()) {
		dofOrder.SetCount(6*Nb);
		for (int i = 0, order = 0; i < 6*Nb; ++i, ++order) 
			dofOrder[i] = order;
	}
	
	if (!IsLoadedA0())  
		GetA0();
	
	if ((!IsLoadedAinf() || !IsLoadedKirf()) && bem->calcAinf) {
		if (IsNull(bem->maxTimeA) || bem->maxTimeA == 0) {
			lastError = t_("Incorrect time for A∞ calculation. Please review it in Options");
			return false;
		}
		if (IsNull(bem->numValsA) || bem->numValsA < 10) {
			lastError = t_("Incorrect number of time values for A∞ calculation. Please review it in Options");
			return false;
		}
		if (!IsLoadedKirf()) {
			if (Status && !Status(t_("Obtaining the Impulse Response Function"), 40)) {
				lastError = t_("Cancelled by the user");
				return false;
			}
			GetK_IRF(min(bem->maxTimeA, GetK_IRF_MaxT()), bem->numValsA);
		}
		if (!IsLoadedAinf()) {
			if (Status && !Status(t_("Obtaining the infinite-frequency added mass (A∞)"), 70)) {
				lastError = t_("Cancelled by the user");
				return false;
			}
			GetAinf();
		}
	}
	if (bem->calcAinf_w) {
		if (!IsLoadedKirf())
			GetK_IRF(min(bem->maxTimeA, GetK_IRF_MaxT()), bem->numValsA);
		if (!IsLoadedAinf())
			GetAinf();
		if (Status && !Status(t_("Obtaining the frequency-dependent infinite-frequency added mass (A∞(ω))"), 90)) {
			lastError = t_("Cancelled by the user");
			return false;
		}
		GetAinf_w();
	}
	if (Ainf_w.size() == 0)
		InitAinf_w();
	if (Ainf.size() == 0)
		Ainf.setConstant(Nb*6, Nb*6, 0);
	
	/*try {
		CheckNaN();
	} catch (Exc e) {
		lastError = e;
		return false;
	}*/
	
	return true;
}

int Hydro::GetW0() {
	for (int i = 0; i < w.size(); ++i) {	
		if (w[i] < 0.0001)
			return i;
	}
	return Null;
}

void Hydro::Get3W0(int &id1, int &id2, int &id3) {
	Upp::Vector<double> ww = clone(w);
	
	Sort(ww);
	id1 = FindAdd(w, ww[0]); 	
	id2 = FindAdd(w, ww[1]); 
	id3 = FindAdd(w, ww[2]); 
}

void Hydro::GetA0() {
	if (!IsLoadedA())
		return;
	
	int iw0 = GetW0();
	if (!IsNull(iw0)) {
		A0.setConstant(Nb*6, Nb*6, Null);
		for (int i = 0; i < Nb*6; ++i)
	        for (int j = 0; j < Nb*6; ++j) 
				A0(i, j) = A[i][j][iw0];
	} else if (w.size() < 3)
		return;
	else { 
		int iw1, iw2, iw3;
		Get3W0(iw1, iw2, iw3);
		double wiw1 = w[iw1];
		double wiw2 = w[iw2];
		double wiw3 = w[iw3];
		if (wiw1 > 1. && wiw1 > 3*(wiw2 - wiw1))	// Too high to guess A[0]
			return;
		A0.setConstant(Nb*6, Nb*6, Null);
		for (int i = 0; i < Nb*6; ++i)
	        for (int j = 0; j < Nb*6; ++j) {
	            if (IsNull(A[i][j][iw1]) || IsNull(A[i][j][iw2]) || IsNull(A[i][j][iw3]))
	                A0(i, j) = Null;
	            else
					A0(i, j) = QuadraticInterpolate<double>(0, wiw1, wiw2, wiw3, A[i][j][iw1], A[i][j][iw2], A[i][j][iw3]);
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

void Hydro::SetC(int ib, const Eigen::MatrixXd &K) {		// K is supposed to be dimensionalized
	if (C.IsEmpty())
		C.SetCount(Nb);
	if (C[ib].size() == 0)
		C[ib].setConstant(6, 6, Null); 
	for (int idf = 0; idf < 6; ++idf) {
		for (int jdf = 0; jdf < 6; ++jdf) {
			double k = dimen ? g_rho_dim()/g_rho_ndim() : g_rho_ndim()*pow(len, GetK_C(idf, jdf));
	      	C[ib](idf, jdf) = K(idf, jdf)/k;
		}
	}
}

int Hydro::GetIrregularHead() const {
	if (Nh <= 2)
		return -1;
	double delta0 = head[1] - head[0];
	for (int i = 1; i < Nh - 1; ++i) {
		double delta = head[i+1] - head[i];
		if (!EqualRatio(delta, delta0, 0.001))
			return i;
	}
	return -1;
}

int Hydro::GetIrregularFreq() const {
	if (Nf <= 2)
		return -1;
	double delta0 = w[1] - w[0];
	for (int i = 1; i < Nf - 1; ++i) {
		double delta = w[i+1] - w[i];
		if (!EqualRatio(delta, delta0, delta0/10))
			return i;
	}
	return -1;
}

double Hydro::g_dim() 		const {return bem->g;}					// Dimensionalize only with system data
double Hydro::g_ndim()		const {return !IsNull(g) ? g : bem->g;}	// Nondimensionalize with model data, if possible
double Hydro::rho_dim() 	const {return bem->rho;}		
double Hydro::rho_ndim()	const {return !IsNull(rho) ? rho : bem->rho;}
double Hydro::g_rho_dim() 	const {return bem->rho*bem->g;}
double Hydro::g_rho_ndim()	const {return g_ndim()*rho_ndim();}

void Hydro::StateSpace::GetTFS(const Upp::Vector<double> &w) {
	Eigen::Index sz = A_ss.rows();
	TFS.SetCount(w.size());
	for (int ifr = 0; ifr < w.size(); ++ifr) {
		std::complex<double> wi = std::complex<double>(0, w[ifr]);
		TFS[ifr] = (C_ss.transpose()*(Eigen::MatrixXd::Identity(sz, sz)*wi - A_ss).inverse()*B_ss)(0);	// C_ss*inv(I*w*i-A_ss)*B_ss
	}		
}

int Hydro::GetHeadId(double hd) const {
	hd = FixHeading_180(hd);
	for (int i = 0; i < head.size(); ++i) {
		if (EqualRatio(head[i], hd, 0.01))
			return i;
	}
	return -1;
}
	
int Hydro::GetQTFHeadId(double hd) const {
	for (int i = 0; i < qtfhead.size(); ++i) {
		if (EqualRatio(qtfhead[i], hd, 0.01))
			return i;
	}
	return -1;
}
	
int Hydro::GetQTFId(int lastid, const Upp::Array<Hydro::QTF> &qtfList, 
			const QTFCases &qtfCases, int ib, int ih1, int ih2, int ifr1, int ifr2) {
	if (qtfCases.ib.size() > 0) {
		bool found = false;
		for (int i = 0; i < qtfCases.ib.size(); ++i) {
			if (qtfCases.ib[i] == ib && qtfCases.ih1[i] == ih1 && qtfCases.ih2[i] == ih2) {
				found = true;
				break;
			}
		}
		if (!found)
			return -1;
	}
	if (lastid < 0)
		lastid = 0;
	for (int i = lastid; i < qtfList.size(); ++i) {
		const QTF &qtf = qtfList[i];
		if (qtf.ib == ib && qtf.ih1 == ih1 && qtf.ih2 == ih2 && qtf.ifr1 == ifr1 && qtf.ifr2 == ifr2) 
			return i;
	}
	return -1;
}

void Hydro::GetQTFList(const Upp::Array<Hydro::QTF> &qtfList, QTFCases &qtfCases, const Vector<double> &headings) {
	qtfCases.Clear();
	for (int i = 0; i < qtfList.size(); ++i) {
		const QTF &qtf = qtfList[i];
		bool found = false;
		for (int j = 0; j < qtfCases.ib.size(); ++j) {
			if (qtf.ib == qtfCases.ib[j] && qtf.ih1 == qtfCases.ih1[j] && qtf.ih2 == qtfCases.ih2[j]) {
				found = true;
				break;
			}
		}
		if (!found) {
			qtfCases.ib << qtf.ib;
			qtfCases.ih1 << qtf.ih1;
			qtfCases.ih2 << qtf.ih2;	
		}
	}
	qtfCases.Sort(headings);
}

Eigen::VectorXd Hydro::B_dim(int idf, int jdf) const {
	if (dimen)
		return B[idf][jdf]*(rho_dim()/rho_ndim());
	else {
		Eigen::VectorXd ret = B[idf][jdf]*(rho_dim()*pow(len, GetK_AB(idf, jdf)));
		Eigen::VectorXd ww = Eigen::Map<Eigen::VectorXd>((double *)w.begin(), w.size());
		return ret.array()*ww.array();
	}
}

Eigen::VectorXd Hydro::B_ndim(int idf, int jdf) const {
	if (!dimen)
		return B[idf][jdf]*(rho_ndim()/rho_dim());
	else {
		Eigen::VectorXd ret = B[idf][jdf]/(rho_ndim()*pow(len, GetK_AB(idf, jdf)));
		Eigen::VectorXd ww = Eigen::Map<Eigen::VectorXd>((double *)w.begin(), w.size());
		return ret.array()/ww.array();
	}
}


void Hydro::GetOldAB(const Upp::Array<Eigen::MatrixXd> &oldAB, Upp::Array<Upp::Array<Eigen::VectorXd>> &AB) {
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

void Hydro::SetOldAB(Upp::Array<Eigen::MatrixXd> &oldAB, const Upp::Array<Upp::Array<Eigen::VectorXd>> &AB) {
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

Eigen::MatrixXd Hydro::Ainf_(int ib, bool ndim) const {
	Eigen::MatrixXd ret;
	if (!IsLoadedAinf())
		return ret;
	ret.resize(6, 6);
	for (int idf = 0; idf < 6; ++idf) 	
		for (int jdf = 0; jdf < 6; ++jdf) 
			ret(idf, jdf) = Ainf_(ndim, idf + 6*ib, jdf + 6*ib);	// It doesn't return added mass between bodies...
	return ret;
}

Eigen::MatrixXd Hydro::C_(int ib, bool ndim) const {
	Eigen::MatrixXd ret;
	if (C.IsEmpty())
		return ret;
	ret.resize(6, 6);
	for (int idf = 0; idf < 6; ++idf) 	
		for (int jdf = 0; jdf < 6; ++jdf) 
			ret(idf, jdf) = C_(ndim, ib, idf, jdf);
	return ret;
}

Eigen::MatrixXd Hydro::Dlin_dim(int ib) const {
	Eigen::MatrixXd ret;
	if (!IsLoadedDlin())
		return ret;
	ret.resize(6, 6);
	for (int idf = 0; idf < 6; ++idf) 	
		for (int jdf = 0; jdf < 6; ++jdf) 
			ret(idf, jdf) = Dlin(ib*6 + idf, ib*6 + jdf);
	return ret;
}

void Hydro::C_dim() {
	if (C.IsEmpty())
		return;
	for (int ib = 0; ib < Nb; ++ib) 
		for (int idf = 0; idf < 6; ++idf) 	
			for (int jdf = 0; jdf < 6; ++jdf) 
				C[ib](idf, jdf) = C_dim(ib, idf, jdf);
}

void Hydro::F_dim(Forces &f) {
	if (f.ma.IsEmpty())
		return;
	for (int ih = 0; ih < Nh; ++ih) 	
		for (int ifr = 0; ifr < Nf; ++ifr)
			for (int idf = 0; idf < 6*Nb; ++idf) {
				f.ma[ih](ifr, idf) = F_ma_dim(f, ih, ifr, idf);
				f.re[ih](ifr, idf) = F_re_dim(f, ih, ifr, idf);
				f.im[ih](ifr, idf) = F_im_dim(f, ih, ifr, idf);
			}
}

void Hydro::CheckNaN() {
	if (!IsNum(A))
		throw Exc("Error loading A. NaN found");
	if (!IsNum(Ainf_w))
		throw Exc("Error loading Ainfw. NaN found");
	if (!IsNum(Ainf))
		throw Exc("Error loading Awinf. NaN found");
	if (!IsNum(A0))
		throw Exc("Error loading A_0. NaN found");
	if (!IsNum(B))
		throw Exc("Error loading B. NaN found");
	if (!IsNum(head))
		throw Exc("Error loading head. NaN found");
	if (!IsNum(M))
		throw Exc("Error loading M. NaN found");
	if (!IsNum(C))
		throw Exc("Error loading C. NaN found");
	if (!IsNum(cb))
		throw Exc("Error loading cb. NaN found");
	if (!IsNum(cg))
		throw Exc("Error loading cg. NaN found");
	if (!IsNum(c0))
		throw Exc("Error loading c0. NaN found");
	if (!IsNum(dof))
		throw Exc("Error loading dof. NaN found");
	if (!IsNum(dofOrder))
		throw Exc("Error loading dofOrder. NaN found");
	if (!IsNum(Kirf))
		throw Exc("Error loading Kirf. NaN found");
	if (!IsNum(Tirf))
		throw Exc("Error loading Tirf. NaN found");
	if (!IsNum(ex))
		throw Exc("Error loading ex. NaN found");
	if (!IsNum(sc))
		throw Exc("Error loading sc. NaN found");
	if (!IsNum(fk))
		throw Exc("Error loading fk. NaN found");
	if (!IsNum(rao))
		throw Exc("Error loading rao. NaN found");	
}

double Hydro::Theave(int ib) const {
	return 2*M_PI*sqrt((M[ib](0, 0) + Ainf_dim(ib*6+2, ib*6+2))/C_dim(ib, 2, 2)); 
}

double Hydro::Troll(int ib) const {
	return 2*M_PI*sqrt((M[ib](3, 3) + Ainf_dim(ib*6+3, ib*6+3))/C_dim(ib, 3, 3)); 
}

double Hydro::Tpitch(int ib) const {
	return 2*M_PI*sqrt((M[ib](4, 4) + Ainf_dim(ib*6+4, ib*6+4))/C_dim(ib, 4, 4)); 
}

double Hydro::GMroll(int ib) const {
	return C_dim(ib, 3, 3)/(rho*g*Vo[ib]);
}

double Hydro::GMpitch(int ib) const {
	return C_dim(ib, 4, 4)/(rho*g*Vo[ib]);
}

void Hydro::Jsonize(JsonIO &json) {
	int icode;
	Upp::Array<Eigen::MatrixXd> oldA, oldB, oldKirf;
	if (json.IsStoring()) {
		icode = code;
		SetOldAB(oldA, A);
		SetOldAB(oldB, B);
		SetOldAB(oldKirf, Kirf);
	}
	json
		("file", file)
		("name", name)
		("g", g)
		("h", h)
		("rho", rho)
		("len", len)
		("dimen", dimen)
		("Nb", Nb)
		("Nf", Nf)
		("Nh", Nh)
		("A", oldA)
		("Awinf", Ainf)
		("Aw0", A0)
		("B", oldB)
		("head", head)
		("names", names)
		("C", C)
		("cb", cb)
		("cg", cg)
		("c0", c0)
		("code", icode)
		("dof", dof)
		("dofOrder", dofOrder)
		("Kirf", oldKirf)
		("Tirf", Tirf)
		("ex", ex)
		("sc", sc)
		("fk", fk)
		("rao", rao)
		("sts", sts)
		("T", T)
		("w", w)
		("dataFromW", dataFromW)
		("Vo", Vo)
		("stsProcessor", stsProcessor)
		("dimenSTS", dimenSTS)
		("description", description)
		("qtfsum", qtfsum)
		("qtfdif", qtfdif)
		("Dlin", Dlin)
	;
	if(json.IsLoading()) {
		code = static_cast<Hydro::BEM_SOFT>(icode);
		GetOldAB(oldA, A);
		GetOldAB(oldB, B);
		GetOldAB(oldKirf, Kirf);
	}
}
	
BEM::BEM() {
	bemFilesAst = clone(bemFilesExt);
	bemFilesAst.Replace(".", "*.");
	experimental = ToLower(GetExeTitle()).Find("experimental") >= 0;
}

void BEM::LoadBEM(String file, Function <bool(String, int)> Status, bool checkDuplicated) {
	Status(t_("Loading files"), 10);
	if (checkDuplicated) {
		for (int i = 0; i < hydros.size(); ++i) {
			if (hydros[i].hd().file == file) 
				throw Exc(Format(t_("Model '%s' is already loaded"), file));
		}
	}
	String ext = ToLower(GetFileExt(file));
	if (ext == ".cal" || ext == ".tec") {
		Nemoh &data = hydros.Create<Nemoh>(*this);
		if (!data.Load(file)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.size()-1);
			throw Exc(Format(t_("Problem loading '%s'\n%s"), file, error));	
		}
	} else if (ext == ".inf") {
		Nemoh &data = hydros.Create<Nemoh>(*this);
		if (!data.Load(file)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.size()-1);
			throw Exc(Format(t_("Problem loading '%s'\n%s"), file, error));	
		}
	} else if (ext == ".out") {
		Wamit &data = hydros.Create<Wamit>(*this);
		if (!data.Load(file, Status)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.size()-1);
			throw Exc(Format(t_("Problem loading '%s'\n%s"), file, error));
		}
	} else if (ext == ".in") {
		HAMS &data = hydros.Create<HAMS>(*this);
		if (!data.Load(file, Status)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.size()-1);
			throw Exc(Format(t_("Problem loading '%s'\n%s"), file, error));
		}
	} else if (ext == ".dat" || ext == ".fst") {
		Fast &data = hydros.Create<Fast>(*this);
		if (!data.Load(file, Status)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.size()-1);
			throw Exc(Format(t_("Problem loading '%s'\n%s"), file, error));		
		}
	} else if (ext == ".1" || ext == ".2" || ext == ".3" || ext == ".hst" || ext == ".4" || ext == ".12s" || ext == ".12d") {
		Wamit &data = hydros.Create<Wamit>(*this);
		if (!data.Load(file, Status)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.size()-1);
			throw Exc(Format(t_("Problem loading '%s'\n%s"), file, error));		
		}
	} else if (ext == ".ah1" || ext == ".lis" || ext == ".qtf") {
		Aqwa &data = hydros.Create<Aqwa>(*this);
		if (!data.Load(file)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.size()-1);
			throw Exc(Format(t_("Problem loading '%s'\n%s"), file, error));	
		}
	} else if (ext == ".mat") {
		Foamm &data = hydros.Create<Foamm>(*this);
		if (!data.Load(file)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.size()-1);
			throw Exc(Format(t_("Problem loading '%s'\n%s"), file, error));	
		}
	} else if (ext == ".bem") {
		HydroClass &data = hydros.Create<HydroClass>(*this);
		if (!data.Load(file)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.size()-1);
			throw Exc(Format(t_("Problem loading '%s'\n%s"), file, error));	
		}
	} else 
		throw Exc(Format(t_("Unknown BEM file extension in '%s'"), file));
	
	Hydro &justLoaded = hydros.Top().hd();

	if (!justLoaded.AfterLoad(Status)) {
		String error = RemoveAccents(justLoaded.GetLastError());
		hydros.SetCount(hydros.size()-1);
		throw Exc(Format(t_("Problem processing '%s'\n%s"), file, error));	
	}
	
	/*if (discardNegDOF) {
		if (Status && !Status(t_("Discarding negligible DOF"), 90)) {
			hydros.SetCount(hydros.size()-1);	
			throw Exc(t_("Cancelled by user"));
		}
		justLoaded.RemoveThresDOF_A(thres);
		justLoaded.RemoveThresDOF_B(thres);
		justLoaded.RemoveThresDOF_Force(justLoaded.ex, thres);
		justLoaded.RemoveThresDOF_Force(justLoaded.sc, thres);
		justLoaded.RemoveThresDOF_Force(justLoaded.fk, thres);
		justLoaded.RemoveThresDOF_Force(justLoaded.rao, thres);
	}*/
	if (hydros.size() == 1)
		Nb = justLoaded.Nb;
	else {
		if (justLoaded.Nb > Nb) {
			//int justLoaded_Nb = justLoaded.Nb;	Error removed for now
			//hydros.SetCount(hydros.size()-1);
			//throw Exc(Format(t_("Model has more bodies (%d) than previously loaded (%d)"), justLoaded_Nb, Nb));
			Nb = justLoaded.Nb;
		}
	}
	UpdateHeadAll();
}

HydroClass &BEM::Join(Upp::Vector<int> &ids, Function <bool(String, int)> Status) {
	Upp::Vector<Hydro *>hydrosp;
	
	hydrosp.SetCount(ids.size());
	for (int i = 0; i < ids.size(); ++i) 
		hydrosp[i] = &hydros[ids[i]].hd(); 
	
	HydroClass &data = hydros.Create<HydroClass>(*this);
	data.hd().Join(hydrosp);
	if (!data.hd().AfterLoad(Status)) {
		String error = data.hd().GetLastError();
		throw Exc(Format(t_("Problem joining models: '%s'\n%s"), error));	
	}
	Sort(ids, StdLess<int>());
	for (int i = ids.size()-1; i >= 0; --i)
		hydros.Remove(ids[i]);
	return data;
}

HydroClass &BEM::Duplicate(int id) {
	HydroClass &data = hydros.Create<HydroClass>(*this);
	data.hd().Copy(hydros[id].hd());
	
	return data;
}

void BEM::Symmetrize(int id, bool xAxis) {
	hydros[id].hd().Symmetrize_Forces(xAxis);
	
	UpdateHeadAll();
}

void BEM::UpdateHeadAll() {
	headAll.Clear();
	orderHeadAll.Clear();
				
	for (int id = 0; id < hydros.size(); ++id) {
		for (int ih = 0; ih < hydros[id].hd().head.size(); ++ih) 
			FindAddDelta(headAll, FixHeading(hydros[id].hd().head[ih], headingType), 0.1);
	}
	orderHeadAll = GetSortOrder(headAll);
}

void BEM::A0(int id) {
	hydros[id].hd().GetA0();
}

void BEM::Kirf(int id, double maxT) {
	hydros[id].hd().GetK_IRF(maxT, numValsA);
}

void BEM::Ainf(int id) {
	hydros[id].hd().GetAinf();
}

void BEM::Ainf_w(int id) {
	hydros[id].hd().GetAinf_w();
}

void BEM::OgilvieCompliance(int id, bool zremoval, bool thinremoval, bool decayingTail, bool haskind) {
	hydros[id].hd().GetOgilvieCompliance(zremoval, thinremoval, decayingTail, haskind);
}

void BEM::DeleteHeadingsFrequencies(int id, const Vector<int> &idFreq, const Vector<int> &idFreqQTF, const Vector<int> &idHead, const Vector<int> &idHeadQTF) {
	hydros[id].hd().DeleteFrequencies(idFreq);
	hydros[id].hd().DeleteFrequenciesQTF(idFreqQTF);
	hydros[id].hd().DeleteHeadings(idHead);
	hydros[id].hd().DeleteHeadingsQTF(idHeadQTF);
}

void BEM::TranslationTo(int id, double xto, double yto, double zto) {
	hydros[id].hd().GetTranslationTo(xto, yto, zto);
}


void BEM::LoadMesh(String fileName, Function <bool(String, int pos)> Status, bool cleanPanels, bool checkDuplicated) {
	Status(Format(t_("Loading mesh '%s'"), fileName), 10);
	
	if (checkDuplicated) {
		for (int i = 0; i < surfs.size(); ++i) {
			if (surfs[i].fileName == fileName) {
				BEM::Print(S("\n") + t_("Model is already loaded"));
				throw Exc(t_("Model is already loaded"));
			}
		}
	}
	Mesh &mesh = surfs.Add();
	String error = mesh.Load(fileName, rho, g, cleanPanels);
	if (!error.IsEmpty()) {
		BEM::Print("\n" + Format(t_("Problem loading '%s'") + S("\n%s"), fileName, error));
		RemoveMesh(surfs.size()-1);
		throw Exc(Format(t_("Problem loading '%s'") + S("\n%s"), fileName, error));
	}
}

void BEM::HealingMesh(int id, bool basic, Function <bool(String, int)> Status) {
	Status(Format(t_("Healing mesh '%s'"), surfs[id].fileName), 10);
	Print(S("\n\n") + Format(t_("Healing mesh '%s'"), surfs[id].fileName));
	
	String ret;
	try {
		ret = surfs[id].Heal(basic, rho, g, Status);
	} catch (Exc e) {
		surfs.SetCount(surfs.size()-1);
		Print("\n" + Format(t_("Problem healing '%s': %s") + S("\n%s"), e));
		throw std::move(e);
	}
	if (!ret.IsEmpty()) {
		ret.Replace("\n", "\n- ");
		Print(ret);
	} else
		Print(S(". ") + t_("The mesh is in good condition"));
}

void BEM::OrientSurface(int id, Function <bool(String, int)> Status) {
	Status(Format(t_("Orienting surface mesh '%s'"), surfs[id].fileName), 10);
	Print(S("\n\n") + Format(t_("Orienting surface mesh '%s'"), surfs[id].fileName));
	
	try {
		surfs[id].Orient();
	} catch (Exc e) {
		surfs.SetCount(surfs.size()-1);
		Print("\n" + Format(t_("Problem orienting surface '%s': %s") + S("\n%s"), e));
		throw std::move(e);
	}
}

void BEM::UnderwaterMesh(int id, Function <bool(String, int pos)> Status) {
	Status(Format(t_("Getting underwater mesh '%s'"), surfs[id].fileName), 10);
	
	Mesh &mesh = surfs.Add();
	Mesh &orig = surfs[id];
	mesh.fileName = orig.fileName;
	
	try {
		mesh.mesh.CutZ(orig.mesh, -1);
	} catch (Exc e) {
		surfs.SetCount(surfs.size()-1);
		Print("\n" + Format(t_("Problem loading '%s': %s") + S("\n%s"), e));
		throw std::move(e);
	}
}

void BEM::RemoveMesh(int id) {
	surfs.Remove(id);
}

void BEM::JoinMesh(int idDest, int idOrig) {
	const Mesh &orig = surfs[idOrig];
	Mesh &dest = surfs[idDest];
	dest.fileName << "/" << orig.fileName;
	
	try {
		dest.Join(orig.mesh, rho, g);
		RemoveMesh(idOrig);
	} catch (Exc e) {
		surfs.SetCount(surfs.size()-1);
		Print("\n" + Format(t_("Problem loading '%s': %s") + S("\n%s"), e));
		throw std::move(e);
	}
}

Upp::Vector<int> BEM::SplitMesh(int id, Function <bool(String, int pos)> Status) {
	Status(Format(t_("Splitting mesh '%s'"), surfs[id].fileName), 0);
	Mesh &orig = surfs[id];
	
	Upp::Vector<int> ret;
	try {
		Upp::Vector<Upp::Vector<int>> sets = orig.mesh.GetPanelSets(Status);
		if (sets.size() == 1)
			return ret;
		for (int i = 0; i < sets.size(); ++i) {		
			Mesh &surf = surfs.Add();
			ret << surfs.size()-1-1;		// One more as id is later removed
			for (int ii = 0; ii < sets[i].size(); ++ii) 
				surf.mesh.panels << clone(orig.mesh.panels[sets[i][ii]]);	
			
			surf.mesh.nodes = clone(orig.mesh.nodes);
			surf.SetCode(orig.GetCode());
			surf.AfterLoad(rho, g, false, true);
		}
		RemoveMesh(id);
	} catch (Exc e) {
		surfs.SetCount(surfs.size() - ret.size());
		Print("\n" + Format(t_("Problem loading '%s': %s") + S("\n%s"), e));
		throw std::move(e);
	}
	return ret;
}

void BEM::AddFlatPanel(double x, double y, double z, double size, double panWidth, double panHeight) {
	try {
		Mesh &surf = surfs.Add();

		surf.SetCode(Mesh::EDIT);
		surf.mesh.AddFlatPanel(panWidth, panHeight, size); 
		surf.mesh.Translate(x, y, z);
	} catch (Exc e) {
		surfs.SetCount(surfs.size() - 1);
		Print("\n" + Format(t_("Problem adding flat panel: %s"), e));
		throw std::move(e);
	}	
}

void BEM::AddRevolution(double x, double y, double z, double size, Upp::Vector<Pointf> &vals) {
	try {
		Mesh &surf = surfs.Add();

		surf.SetCode(Mesh::EDIT);
		surf.mesh.AddRevolution(vals, size); 
		surf.mesh.Translate(x, y, z);
	} catch (Exc e) {
		surfs.SetCount(surfs.size() - 1);
		Print("\n" + Format(t_("Problem adding revolution surface: %s"), e));
		throw std::move(e);
	}	
}

void BEM::AddPolygonalPanel(double x, double y, double z, double size, Upp::Vector<Pointf> &vals) {
	try {
		Mesh &surf = surfs.Add();

		surf.SetCode(Mesh::EDIT);
		surf.mesh.AddPolygonalPanel(vals, size, true); 
		surf.mesh.Translate(x, y, z);
	} catch (Exc e) {
		surfs.SetCount(surfs.size() - 1);
		Print("\n" + Format(t_("Problem adding revolution surface: %s"), e));
		throw std::move(e);
	}	
}

void BEM::AddWaterSurface(int id, char c) {
	try {
		Mesh &surf = surfs.Add();

		surf.SetCode(Mesh::EDIT);
		surf.mesh.AddWaterSurface(surfs[id].mesh, surfs[id].under, c); 
		
		if (c == 'r')
			surf.name = t_("Water surface removed");
		else if (c == 'f')
			surf.name = t_("Water surface");
		else if (c == 'e')
			surf.name = t_("Water surface extracted");
		surf.name = surfs[id].name + " " + surf.name;
		surf.fileName =  "";
		
		surf.AfterLoad(rho, g, false, false);
		
		//surf.Report(rho);
	} catch (Exc e) {
		surfs.SetCount(surfs.size() - 1);
		Print("\n" + Format(t_("Problem adding revolution surface: %s"), e));
		throw std::move(e);
	}	
}
			
bool BEM::LoadSerializeJson() {
	bool ret;
	String folder = AppendFileNameX(GetAppDataFolder(), "BEMRosetta");
	if (!DirectoryCreateX(folder))
		ret = false;
	else {
		String fileName = AppendFileNameX(folder, "configdata.cf");
		if (!FileExists(fileName)) 
			ret = false;
		else {
			String jsonText = LoadFile(fileName);
			if (jsonText.IsEmpty())
				ret = false;
			else {
				if (!LoadFromJson(*this, jsonText))
					ret = false;
				else
					ret = true;
			}
		}
	}
	if (!ret || IsNull(g)) 
		g = 9.81;
	if (!ret || IsNull(depth)) 
		depth = 100;
	if (!ret || IsNull(rho)) 
		rho = 1000;
	if (!ret || IsNull(len)) 
		len = 1;
	//if (!ret || IsNull(discardNegDOF))
	//	discardNegDOF = false;
	//if (!ret || IsNull(thres)) 
	//	thres = 0.01;
	if (!ret || IsNull(calcAinf))
		calcAinf = true;
	if (!ret || IsNull(calcAinf_w))
		calcAinf_w = true;
	if (!ret || IsNull(maxTimeA))
		maxTimeA = 120;
	if (!ret || IsNull(numValsA))
		numValsA = 1000;
	if (!ret || IsNull(onlyDiagonal))
		onlyDiagonal = false;
	if (!ret || IsNull(volWarning))	
		volWarning = 1;
	if (!ret || IsNull(volError))
		volError = 10;
			
	return ret;
}

bool BEM::ClearTempFiles() {
	String folder = GetTempFilesFolder();
	DeleteFolderDeepWildcardsX(folder, "*.*");	Sleep(100);
	return DirectoryCreateX(folder);
}
	
bool BEM::StoreSerializeJson() {
	String folder = AppendFileNameX(GetAppDataFolder(), "BEMRosetta");
	if (!DirectoryCreateX(folder))
		return 0;
	String fileName = AppendFileNameX(folder, "configdata.cf");
	return StoreAsJsonFile(*this, fileName, true);
}


bool HydroClass::Load(String file) {
	BEM::Print("\n\n" + Format(t_("Loading '%s'"), file));
	
	if (!LoadFromJsonFile(hd(), file)) {
		BEM::PrintError("\n" + Format(t_("Error loading '%s'"), file));
		hd().lastError = "\n" + Format(t_("Error loading '%s'"), file);
		return false;
	}
	hd().file = file;
	return true;
}
	
bool HydroClass::Save(String file) {
	BEM::Print("\n\n" + Format(t_("Saving '%s'"), file));
	if (!StoreAsJsonFile(hd(), file, true)) {
		BEM::PrintError("\n" + Format(t_("Error saving '%s'"), file));
		hd().lastError = "\n" + Format(t_("Error saving '%s'"), file);
		return false;
	}
	return true;
}

int IsTabSpace(int c) {
	if (c == '\t' || c == ' ' || c == '!')
		return true;
	return false;
}

Upp::Vector<int> NumSets(int num, int numsets) {
	ASSERT(numsets > 0);
	Upp::Vector<int> ret;
	ret.SetCount(numsets);
	
	for (int i = 0; numsets > 0; ++i) {
		int delta = int(num/numsets);
		ret[i] = delta;
		num -= delta;
		numsets--;
	}
	return ret;
}

String FormatWam(double d) {
	if (IsNull(d))
		return "0.0";
	return (d >= 0 ? " " : "-") + Format("%12E", abs(d));
}

void FieldSplitWamit::LoadWamitJoinedFields(String _line) {	
	line = _line;
	fields.Clear();
	Upp::Vector<String> prefields = Split(line, IsTabSpace, true);
	for (int id = 0; id < prefields.size(); ++id) {
		String s = prefields[id];
		String ns;
		for (int i = 0; i < s.GetCount(); ++i) {	
			int c = s[i];
			if (c == '-') {
				if (i == 0)
					ns.Cat(c);
				else if (s[i-1] == 'E')
					ns.Cat(c);
				else {
					fields << ns;
					ns.Clear();
					ns.Cat(c);
				}
			} else
				ns.Cat(c);
		}
		fields << ns;
	}
}
	
BEMBody::BEMBody() {
	dof.SetCount(6, false);	
	cg = Eigen::Vector3d::Zero();
	c0 = Eigen::Vector3d::Zero();
	mass.setConstant(6, 6, 0);
	linearDamping.setConstant(6, 6, 0);
	quadraticDamping.setConstant(6, 6, 0);
	hydrostaticRestoring.setConstant(6, 6, 0);
	externalRestoring.setConstant(6, 6, 0);
}
	
int BEMBody::GetNDOF() const {
	int ret = 0;
	for (auto &d : dof)
		ret += d;
	return ret;
}

void BEMCase::Load(String file, const BEM &bem) {
	if (ToLower(GetFileName(file)) == "nemoh.cal") {
		if (!static_cast<NemohCase&>(*this).Load(file)) 
			throw Exc(Format(t_("Problem loading '%s' file"), file));
	} else if (ToLower(GetFileExt(file)) == ".in") {
		if (!static_cast<HamsCase&>(*this).Load(file)) 
			throw Exc(Format(t_("Problem loading '%s' file"), file));
	} else if (ToLower(GetFileExt(file)) == ".dat") {
		if (!static_cast<AQWACase&>(*this).Load(file)) 
			throw Exc(Format(t_("Problem loading '%s' file"), file));
	} else
		throw Exc(t_("Unknown BEM input format"));
	
	if (IsNull(rho))
		rho = bem.rho;
	if (IsNull(g))
		g = bem.g;	
}

void BEMCase::SaveFolder(String folder, bool bin, int numCases, int numThreads, const BEM &bem, int solver) const {
	if (solver <= CAPYTAINE)
		static_cast<const NemohCase &>(*this).SaveFolder(folder, bin, numCases, numThreads, bem, solver);
	else if (solver == HAMS)
		static_cast<const HamsCase &>(*this).SaveFolder(folder, bin, numCases, numThreads, bem, solver);
	else if (solver == NEMOH)
		static_cast<const NemohCase &>(*this).SaveFolder(folder, bin, numCases, Null, bem, solver);
	else if (solver == NEMOHv115)
		static_cast<const NemohCase &>(*this).SaveFolder(folder, bin, numCases, Null, bem, solver);
	else
		static_cast<const AQWACase &>(*this).SaveFolder(folder, bin, numCases, numThreads, bem, solver);
}

void BEMCase::BeforeSave(String folderBase, int numCases, bool deleteFolder) const {
	if (numCases < 1)
		throw Exc(Format(t_("Number cases must be higher than 1 (%d)"), numCases));
	
	if (numCases > Nf)
		throw Exc(Format(t_("Number of cases %d must not be higher than number of frequencies %d"), numCases, Nf));
	
	if (deleteFolder) {		// If called from GUI, user has been warned
		if (!DeleteFileDeepWildcardsX(folderBase))
			throw Exc(Format(t_("Impossible to clean folder '%s'. Maybe it is in use"), folderBase));
		Sleep(100);
	}
	if (!DirectoryCreateX(folderBase))
		throw Exc(Format(t_("Problem creating '%s' folder"), folderBase));
}

Vector<String> BEMCase::Check(int solver) const {
	Vector<String> ret;
	
	bool isNemoh = solver != BEMCase::HAMS;
	bool oneBody = solver == BEMCase::HAMS;
	
	if (isNemoh) 
		ret.Append(static_cast<const NemohCase &>(*this).Check());
	else if (solver == BEMCase::HAMS)
		ret.Append(static_cast<const HamsCase &>(*this).Check());
	
	if (IsNull(h) || h < -1 || h > 100000)
		ret << Format(t_("Incorrect depth %s"), FormatDoubleEmpty(h));

	if (IsNull(Nf) || Nf < 1 || Nf > 1000)
		ret << Format(t_("Incorrect number of frequencies %s"), FormatIntEmpty(Nf));
	if (IsNull(minF) || minF < 0)
		ret << Format(t_("Incorrect min frequency %s"), FormatDoubleEmpty(minF));
	if (IsNull(maxF) || maxF < minF)
		ret << Format(t_("Minimum frequency %s has to be lower than maximum frequency %s"), FormatDoubleEmpty(minF), FormatDoubleEmpty(maxF));	
	
	if (IsNull(Nh) || Nh < 1 || Nh > 1000)
		ret << Format(t_("Incorrect number of headings %s"), FormatIntEmpty(Nh));
	if (IsNull(minH) || minH < -180)
		ret << Format(t_("Incorrect min heading %s"), FormatDoubleEmpty(minH));
	if (IsNull(maxH) || maxH > 360)
		ret << Format(t_("Incorrect max heading %s"), FormatDoubleEmpty(maxH));
	if (maxH < minH)
		ret << Format(t_("Minimum heading %s has to be lower than maximum heading %s"), FormatDoubleEmpty(minH), FormatDoubleEmpty(maxH));	

	if(bodies.size() < 1)
		ret << t_("The case has to include at least one body");
	if (oneBody && bodies.size() > 1)
		ret << t_("This solver just processes one body");
	
	return ret;
}

String FormatDoubleEmpty(double val) {
	if (IsNull(val))
		return t_("'empty'");
	else
		return FormatDoubleSize(val, 10, false);
}

String FormatIntEmpty(int val) {
	if (IsNull(val))
		return t_("'empty'");
	else
		return FormatInt(val);
}

bool IsNum(const Hydro::Forces &f) {
	return IsNum(f.ma) && IsNum(f.ph) && IsNum(f.re) && IsNum(f.im);
}

