#include "BEMRosetta.h"
#include <ScatterDraw/DataSource.h>

#include <plugin/matio/matio.h>

Function <void(String)> BEMData::Print 		  = [](String s) {Cout() << s;};
Function <void(String)> BEMData::PrintWarning = [](String s) {Cout() << s;};
Function <void(String)> BEMData::PrintError   = [](String s) {Cout() << s;};

const char *Hydro::strDOF[] 	 = {t_("surge"), t_("sway"), t_("heave"), t_("roll"), t_("pitch"), t_("yaw")};
const char *Hydro::strDOFAbrev[] = {t_("s"), t_("w"), t_("h"), t_("r"), t_("p"), t_("y")};
int Hydro::idCount = 0;	

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
	if (IsLoadedAwinf()) {
		for (int i = 0; i < 6*Nb; ++i) 
			for (int j = 0; j < 6*Nb; ++j) 
				Awinf(i, j) = Awinf_ndim(i, j);
	}
	if (IsLoadedAw0()) {
		for (int i = 0; i < 6*Nb; ++i) 
			for (int j = 0; j < 6*Nb; ++j) 
				Aw0(i, j) = Aw0_ndim(i, j);
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
	if (IsLoadedAwinf()) {	
		for (int i = 0; i < 6*Nb; ++i) 
			for (int j = 0; j < 6*Nb; ++j) 
				Awinf(i, j) = Awinf_dim(i, j);
	}
	if (IsLoadedAw0()) {	
		for (int i = 0; i < 6*Nb; ++i) 
			for (int j = 0; j < 6*Nb; ++j) 
				Aw0(i, j) = Aw0_ndim(i, j);
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

void Hydro::Symmetrize_ForcesEach(const Forces &f, Forces &newf, const Upp::Vector<double> &newHead, int newNh) {
	Initialize_Forces(newf, newNh);
	
	for (int idb = 0; idb < 6*Nb; ++idb) {
		for (int ih = 0; ih < Nh; ++ih) {
			Symmetrize_Forces_Each0(f, newf, newHead, head[ih], ih, idb);
			Symmetrize_Forces_Each0(f, newf, newHead, -head[ih], ih, idb);
		}
	}
}

void Hydro::Symmetrize_Forces() {
	if (!IsLoadedFex() && !IsLoadedFsc() && !IsLoadedFfk() && !IsLoadedRAO())
		return;
	
	Upp::Vector<double> newHead;
	for (int ih = 0; ih < Nh; ++ih) {
		FindAddRatio(newHead, head[ih], 0.001);
		FindAddRatio(newHead, -head[ih], 0.001);
	}
	Sort(newHead);
	int newNh = newHead.size();
	
	Forces newex, newsc, newfk;
	RAO newrao;
	
	if (IsLoadedFex()) {
		Symmetrize_ForcesEach(ex, newex, newHead, newNh);
		ex = pick(newex);
	}
	if (IsLoadedFsc()) {
		Symmetrize_ForcesEach(sc, newsc, newHead, newNh);
		sc = pick(newsc);
	}
	if (IsLoadedFfk()) {
		Symmetrize_ForcesEach(fk, newfk, newHead, newNh);
		fk = pick(newfk);
	}
	if (IsLoadedRAO()) {
		Symmetrize_ForcesEach(rao, newrao, newHead, newNh);
		rao = pick(newrao);
	}
	head = pick(newHead);
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
					Awinf(idf, jdf) = Null;		
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

void Hydro::SaveAs(String file, BEM_SOFT type, int qtfHeading) {
	int realNh = Nh;
	int realNf = Nf;
	
	if (type == UNKNOWN) {
		String ext = ToLower(GetFileExt(file));
		
		if (ext == ".1" || ext == ".3" || ext == ".hst" || ext == ".4" || ext == ".12s" || ext == ".12d") 
			type = Hydro::WAMIT_1_3;
		else if (ext == ".dat")
			type = Hydro::FAST_WAMIT;	
		else if (ext == ".bem")
			type = Hydro::BEMROSETTA;
		else
			throw Exc(Format(t_("Conversion to type of file '%s' not supported"), file));
	}
	if (type == WAMIT) {
		Wamit data(*bem, this);
		data.Save_out(file);			
	} else if (type == WAMIT_1_3) {
		Wamit data(*bem, this);
		data.Save(file, true, qtfHeading);	
	} else if (type == FAST_WAMIT) {
		Fast data(*bem, this);
		data.Save(file, qtfHeading);		
	} else if (type == BEMROSETTA) {
		HydroClass data(*bem, this);
		data.Save(file);		
	}
	code = type;
	Nh = realNh;
	Nf = realNf;
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
			}
			dof[ib] = hydro.dof[ib];
		}
		
		if (hydro.IsLoadedC()) {
			for (int ib = 0; ib < Nb; ++ib) {
				for (int idf = 0; idf < 6; ++idf) 
					for (int jdf = 0; jdf < 6; ++jdf) 
						C[ib](idf, jdf) = hydro.C_ndim(ib, idf, jdf);
			}
		}
		///////////////////////////////////////////////////////////////////
		
		if (hydro.IsLoadedA() && hydro.IsLoadedB()) {
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
		
		if (hydro.IsLoadedRAO()) {
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
	
	bem->calcAwinf = true;
	
	// sts has to be recalculated	
}

void Hydro::Report() {
	BEMData::Print("\n" + Format(t_("%s file '%s'"), GetCodeStr(), file));
	BEMData::Print("\n" + Format(t_("g [m/s2]: %s, h [m]: %s, rho [kg/m3]: %s, length scale [m]: %s"), 
								S_g(), S_h(), S_rho(), S_len()));
	String freqs;
	if (w.IsEmpty()) 
		freqs = t_("NONE");
	else if (w.size() > 1) {
		String strDeltaH;
		if (GetIrregularFreq() < 0) 
			strDeltaH = Format(t_("delta %s [rad/s]"), FormatDouble(w[1] - w[0], 5, FD_EXP));
		else {
			String strHead;
			for (int i = 0; i < w.size(); ++i) {
				if (i > 0)
					strHead << ", ";
				strHead << w[i];
			}
			strDeltaH = Format(t_("Non constant delta (%s)"), strHead); 
		}
	 	freqs = Format(t_("%s to %s %s"), FormatDouble(w[0], 3, FD_EXP), 
	 									  FormatDouble(w[w.size()-1], 3, FD_EXP), strDeltaH);	
	} else
		freqs = Format(t_("%s [rad/s]"), FormatDouble(w[0], 3, FD_EXP));
	
	String heads;
	if (head.IsEmpty())
		heads = t_("NONE");
	else if (head.size() > 1) {
		String strDeltaH;
		if (GetIrregularHead() < 0) 
			strDeltaH = Format(t_("delta %.1f [deg]"), head[1] - head[0]);
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
		heads = Format(t_("%.1f [deg]"), head[0]);
	
	BEMData::Print("\n" + Format(t_("#freqs: %d (%s)"), Nf, freqs)); 
	BEMData::Print("\n" + Format(t_("#headings: %d (%s)"), Nh, heads)); 
	BEMData::Print("\n" + Format(t_("#bodies: %d"), Nb));
	for (int ib = 0; ib < Nb; ++ib) {
		String str = Format("\n%d.", ib+1);
		if (names.size() > ib)
			str += " '" + names[ib] + "'";
		if (dof.size() > ib)
			str += S(" ") + t_("dof") + ": " + FormatInt(dof[ib]);
		if (Vo.size() > ib && !IsNull(Vo[ib]))
			str += S(" ") + t_("vol [m3]") + ": " + FormatDouble(Vo[ib]);
		if (cg.size() > 3*ib && !IsNull(cg(0, ib)))
			str += " " + Format("Cg(%.3f, %.3f, %.3f)[m]", cg(0, ib), cg(1, ib), cg(2, ib));
		if (cb.size() > 3*ib && !IsNull(cb(0, ib)))
			str += " " + Format("Cb(%.3f, %.3f, %.3f)[m]", cb(0, ib), cb(1, ib), cb(2, ib));
		
		BEMData::Print(str);
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
	
	if (!IsLoadedAw0())  
		A0();
	
	if (!IsLoadedAwinf() && bem->calcAwinf) {
		if (IsNull(bem->maxTimeA) || bem->maxTimeA == 0) {
			lastError = t_("Incorrect time for Ainf calculation. Please review it in Options");
			return false;
		}
		if (IsNull(bem->numValsA) || bem->numValsA < 10) {
			lastError = t_("Incorrect number of time values for Ainf calculation. Please review it in Options");
			return false;
		}
		if (!Status(t_("Obtaining Impulse Response Function"), 40)) {
			lastError = t_("Cancelled by user");
			return false;
		}
		K_IRF(min(bem->maxTimeA, GetK_IRF_MaxT()), bem->numValsA);
		if (!Status(t_("Obtaining Infinite-Frequency Added Mass (A_inf)"), 70)) {
			lastError = t_("Cancelled by user");
			return false;
		}
		Ainf();
	}
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

void Hydro::A0() {
	if (!IsLoadedA())
		return;
	
	int iw0 = GetW0();
	if (!IsNull(iw0)) {
		Aw0.setConstant(Nb*6, Nb*6, Null);
		for (int i = 0; i < Nb*6; ++i)
	        for (int j = 0; j < Nb*6; ++j) 
				Aw0(i, j) = A[i][j][iw0];
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
		Aw0.setConstant(Nb*6, Nb*6, Null);
		for (int i = 0; i < Nb*6; ++i)
	        for (int j = 0; j < Nb*6; ++j) {
	            if (IsNull(A[i][j][iw1]) || IsNull(A[i][j][iw2]) || IsNull(A[i][j][iw3]))
	                Aw0(i, j) = Null;
	            else
					Aw0(i, j) = QuadraticInterpolate<double>(0, wiw1, wiw2, wiw3, A[i][j][iw1], A[i][j][iw2], A[i][j][iw3]);
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
		if (!EqualRatio(delta, delta0, 0.001))
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
	
int Hydro::GetQTFId(const Upp::Array<Hydro::QTF> &qtfList, int ib, int ih1, int ih2, int ifr1, int ifr2) {
	for (int i = 0; i < qtfList.size(); ++i) {
		const QTF &qtf = qtfList[i];
		if (qtf.ib == ib && qtf.ih1 == ih1 && qtf.ih2 == ih2 && qtf.ifr1 == ifr1 && qtf.ifr2 == ifr2)
			return i;
	}
	return -1;
}

void Hydro::GetQTFList(const Upp::Array<Hydro::QTF> &qtfList, Upp::Vector<int> &ibL, Upp::Vector<int> &ih1L, Upp::Vector<int> &ih2L) {
	ibL.Clear();
	ih1L.Clear();
	ih2L.Clear();
	for (int i = 0; i < qtfList.size(); ++i) {
		const QTF &qtf = qtfList[i];
		bool found = false;
		for (int j = 0; j < ibL.size(); ++j) {
			if (qtf.ib == ibL[j] && qtf.ih1 == ih1L[j] && qtf.ih2 == ih2L[j]) {
				found = true;
				break;
			}
		}
		if (!found) {
			ibL << qtf.ib;
			ih1L << qtf.ih1;
			ih2L << qtf.ih2;	
		}
	}
}

Eigen::VectorXd Hydro::B_dim(int idf, int jdf) const {
	if (dimen)
		return B[idf][jdf]*(rho_dim()/rho_ndim());
	else {
		Eigen::VectorXd ret = B[idf][jdf]*(rho_dim()*pow(len, GetK_AB(idf, jdf)));
		return ret.cwiseProduct(Eigen::Map<Eigen::VectorXd>((double *)w.begin(), w.size()));
	}
}

Eigen::VectorXd Hydro::B_ndim(int idf, int jdf) const {
	if (!dimen)
		return B[idf][jdf]*(rho_ndim()/rho_dim());
	else {
		Eigen::VectorXd ret = B[idf][jdf]/(rho_ndim()*pow(len, GetK_AB(idf, jdf)));
		return ret.cwiseProduct(Eigen::Map<Eigen::VectorXd>((double *)w.begin(), w.size()).inverse());
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
		("Awinf", Awinf)
		("Aw0", Aw0)
		("B", oldB)
		("head", head)
		("names", names)
		("C", C)
		("cb", cb)
		("cg", cg)
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
	;
	if(json.IsLoading()) {
		code = static_cast<Hydro::BEM_SOFT>(icode);
		GetOldAB(oldA, A);
		GetOldAB(oldB, B);
		GetOldAB(oldKirf, Kirf);
	}
}
	
BEMData::BEMData() {
	bemFilesAst = clone(bemFilesExt);
	bemFilesAst.Replace(".", "*.");
	experimental = ToLower(GetExeTitle()).Find("experimental") >= 0;
}

void BEMData::Load(String file, Function <bool(String, int)> Status, bool checkDuplicated) {
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
		if (!data.Load(file)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.size()-1);
			throw Exc(Format(t_("Problem loading '%s'\n%s"), file, error));
		}
	} else if (ext == ".dat") {
		Fast &data = hydros.Create<Fast>(*this);
		if (!data.Load(file)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.size()-1);
			throw Exc(Format(t_("Problem loading '%s'\n%s"), file, error));		
		}
	} else if (ext == ".1" || ext == ".3" || ext == ".hst" || ext == ".4" || ext == ".12s" || ext == ".12d") {
		Wamit &data = hydros.Create<Wamit>(*this);
		if (!data.Load(file)) {
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
		String error = justLoaded.GetLastError();
		hydros.SetCount(hydros.size()-1);
		throw Exc(Format(t_("Problem processing '%s'\n%s"), file, error));	
	}
	
	if (discardNegDOF) {
		if (!Status(t_("Discarding negligible DOF"), 90)) {
			hydros.SetCount(hydros.size()-1);	
			throw Exc(t_("Cancelled by user"));
		}
		justLoaded.RemoveThresDOF_A(thres);
		justLoaded.RemoveThresDOF_B(thres);
		justLoaded.RemoveThresDOF_Force(justLoaded.ex, thres);
		justLoaded.RemoveThresDOF_Force(justLoaded.sc, thres);
		justLoaded.RemoveThresDOF_Force(justLoaded.fk, thres);
		justLoaded.RemoveThresDOF_Force(justLoaded.rao, thres);
	}
	if (hydros.size() == 1)
		Nb = justLoaded.Nb;
	else {
		if (justLoaded.Nb > Nb) {
			int justLoaded_Nb = justLoaded.Nb;
			hydros.SetCount(hydros.size()-1);
			throw Exc(Format(t_("Model has more bodies (%d) than previously loaded (%d)"), justLoaded_Nb, Nb));
		}
	}
	for (int i = 0; i < justLoaded.head.size(); ++i) 
		FindAddRatio(headAll, justLoaded.head[i], 0.01);
	Sort(headAll);
}

HydroClass &BEMData::Join(Upp::Vector<int> &ids, Function <bool(String, int)> Status) {
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

void BEMData::Symmetrize(int id) {
	hydros[id].hd().Symmetrize_Forces();
	
	for (int i = 0; i < hydros[id].hd().head.size(); ++i) 
		FindAddRatio(headAll, hydros[id].hd().head[i], 0.01);
	Sort(headAll);
}

void BEMData::A0(int id) {
	hydros[id].hd().A0();
}

void BEMData::Ainf(int id, double maxT) {
	hydros[id].hd().K_IRF(maxT, numValsA);
	hydros[id].hd().Ainf();
}

void BEMData::LoadMesh(String fileName, Function <void(String, int pos)> Status, bool cleanPanels, bool checkDuplicated) {
	Status(Format(t_("Loaded mesh '%s'"), fileName), 10);
	
	if (checkDuplicated) {
		for (int i = 0; i < surfs.size(); ++i) {
			if (surfs[i].fileName == fileName) {
				BEMData::Print(S("\n") + t_("Model is already loaded"));
				throw Exc(t_("Model is already loaded"));
			}
		}
	}
	MeshData &mesh = surfs.Add();
	String error = mesh.Load(fileName, rho, g, cleanPanels);
	if (!error.IsEmpty()) {
		BEMData::Print("\n" + Format(t_("Problem loading '%s'") + S("\n%s"), fileName, error));
		RemoveMesh(surfs.size()-1);
		throw Exc(Format(t_("Problem loading '%s'") + S("\n%s"), fileName, error));
	}
}

void BEMData::HealingMesh(int id, bool basic, Function <void(String, int)> Status) {
	Status(Format(t_("Healing mesh '%s'"), surfs[id].fileName), 10);
	Print(S("\n\n") + Format(t_("Healing mesh '%s'"), surfs[id].fileName));
	
	String ret;
	try {
		ret = surfs[id].Heal(basic, Status);
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

void BEMData::OrientSurface(int id, Function <void(String, int)> Status) {
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

void BEMData::UnderwaterMesh(int id, Function <void(String, int pos)> Status) {
	Status(Format(t_("Getting underwater mesh '%s'"), surfs[id].fileName), 10);
	
	MeshData &mesh = surfs.Add();
	MeshData &orig = surfs[id];
	mesh.fileName = orig.fileName;
	
	try {
		mesh.mesh.CutZ(orig.mesh, -1);
	} catch (Exc e) {
		surfs.SetCount(surfs.size()-1);
		Print("\n" + Format(t_("Problem loading '%s': %s") + S("\n%s"), e));
		throw std::move(e);
	}
}

void BEMData::RemoveMesh(int id) {
	surfs.Remove(id);
}

void BEMData::JoinMesh(int idDest, int idOrig) {
	const MeshData &orig = surfs[idOrig];
	MeshData &dest = surfs[idDest];
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

Upp::Vector<int> BEMData::SplitMesh(int id, Function <void(String, int pos)> Status) {
	Status(Format(t_("Splitting mesh '%s'"), surfs[id].fileName), 0);
	MeshData &orig = surfs[id];
	
	Upp::Vector<int> ret;
	try {
		Upp::Vector<Upp::Vector<int>> sets = orig.mesh.GetPanelSets(Status);
		if (sets.size() == 1)
			return ret;
		for (int i = 0; i < sets.size(); ++i) {		
			MeshData &surf = surfs.Add();
			ret << surfs.size()-1-1;		// One more as id is later removed
			for (int ii = 0; ii < sets[i].size(); ++ii) 
				surf.mesh.panels << clone(orig.mesh.panels[sets[i][ii]]);	
			
			surf.mesh.nodes = clone(orig.mesh.nodes);
			surf.SetCode(orig.GetCode());
			surf.AfterLoad(rho, g, false);
		}
		RemoveMesh(id);
	} catch (Exc e) {
		surfs.SetCount(surfs.size() - ret.size());
		Print("\n" + Format(t_("Problem loading '%s': %s") + S("\n%s"), e));
		throw std::move(e);
	}
	return ret;
}

void BEMData::AddFlatPanel(double x, double y, double z, double size, double panWidth, double panHeight) {
	try {
		MeshData &surf = surfs.Add();

		surf.SetCode(MeshData::EDIT);
		surf.mesh.AddFlatPanel(panWidth, panHeight, size); 
		surf.mesh.Translate(x, y, z);
	} catch (Exc e) {
		surfs.SetCount(surfs.size() - 1);
		Print("\n" + Format(t_("Problem adding flat panel: %s"), e));
		throw std::move(e);
	}	
}

void BEMData::AddRevolution(double x, double y, double z, double size, Upp::Vector<Pointf> &vals) {
	try {
		MeshData &surf = surfs.Add();

		surf.SetCode(MeshData::EDIT);
		surf.mesh.AddRevolution(vals, size); 
		surf.mesh.Translate(x, y, z);
	} catch (Exc e) {
		surfs.SetCount(surfs.size() - 1);
		Print("\n" + Format(t_("Problem adding revolution surface: %s"), e));
		throw std::move(e);
	}	
}

void BEMData::AddPolygonalPanel(double x, double y, double z, double size, Upp::Vector<Pointf> &vals) {
	try {
		MeshData &surf = surfs.Add();

		surf.SetCode(MeshData::EDIT);
		surf.mesh.AddPolygonalPanel(vals, size); 
		surf.mesh.Translate(x, y, z);
	} catch (Exc e) {
		surfs.SetCount(surfs.size() - 1);
		Print("\n" + Format(t_("Problem adding revolution surface: %s"), e));
		throw std::move(e);
	}	
}

			
bool BEMData::LoadSerializeJson(bool &firstTime) {
	bool ret;
	String folder = AppendFileName(GetAppDataFolder(), "BEMRosetta");
	DirectoryCreate(folder);
	if (!DirectoryExists(folder))
		ret = false;
	else {
		String fileName = AppendFileName(folder, "configdata.cf");
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
	if (!ret || IsNull(discardNegDOF))
		discardNegDOF = false;
	if (!ret || IsNull(thres)) 
		thres = 0.01;
	if (!ret || IsNull(calcAwinf))
		calcAwinf = true;
	if (!ret || IsNull(maxTimeA))
		maxTimeA = 120;
	if (!ret || IsNull(numValsA))
		numValsA = 1000;
	if (!ret || IsNull(onlyDiagonal))
		onlyDiagonal = false;
	
	firstTime = !ret;
	return true;
}

bool BEMData::ClearTempFiles() {
	String folder = GetTempFilesFolder();
	DeleteFolderDeepWildcardsX(folder, "*.*");	Sleep(100);
	DirectoryCreate(folder);
	return DirectoryExists(folder);
}
	
bool BEMData::StoreSerializeJson() {
	String folder = AppendFileName(GetAppDataFolder(), "BEMRosetta");
	DirectoryCreate(folder);
	if (!DirectoryExists(folder))
		return 0;
	String fileName = AppendFileName(folder, "configdata.cf");
	return StoreAsJsonFile(*this, fileName, true);
}


bool HydroClass::Load(String file) {
	BEMData::Print("\n\n" + Format(t_("Loading '%s'"), file));
	
	if (!LoadFromJsonFile(hd(), file)) {
		BEMData::PrintError("\n" + Format(t_("Error loading '%s'"), file));
		hd().lastError = "\n" + Format(t_("Error loading '%s'"), file);
		return false;
	}
	hd().file = file;
	return true;
}
	
bool HydroClass::Save(String file) {
	BEMData::Print("\n\n" + Format(t_("Saving '%s'"), file));
	if (!StoreAsJsonFile(hd(), file, true)) {
		BEMData::PrintError("\n" + Format(t_("Error saving '%s'"), file));
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
	return (d >= 0 ? " " : "-") + Format("%12E", abs(d));
}

void ShowHelp(BEMData &md) {
	Cout() << "\n" << t_("Usage: bemrosetta_cl [options] [-i infile]... [-e outfile]");
	Cout() << "\n";
	Cout() << "\n" << t_("Options:");
	Cout() << "\n" << t_("-h  --help     -- print options");
	Cout() << "\n" << t_("-p  --params   -- set physical parameters:");
	Cout() << "\n" << t_("                 parameter description   units  default value");
	Cout() << "\n" << t_("                    g      gravity       [m/s2]    ") << md.g;
	Cout() << "\n" << t_("                    length length scale  []        ") << md.len;
	Cout() << "\n" << t_("                    rho    water density [Kg/m3]   ") << md.rho;
	Cout() << "\n" << t_("                    depth  water depth   [m]       ") << md.depth;
	//Cout() << "\n" << t_("                    thres  threshold to discard DOF") << md.thres;
	Cout() << "\n" << t_("-i  --input    -- load model");
	Cout() << "\n" << t_("-e  --export   -- export from input file to output file");
	Cout() << "\n" << t_("-c  --compare  -- compare input files");
	Cout() << "\n" << t_("-r  --report   -- output last loaded model data");
	Cout() << "\n" << t_("-cl --clear    -- clear loaded model");
	Cout() << "\n";
	Cout() << "\n" << t_("Actions");
	Cout() << "\n" << t_("- are done in sequence: if a physical parameter is changed after export, saved files will not include the change");
	Cout() << "\n" << t_("- can be repeated as desired");
}

void CheckNumArgs(const Upp::Vector<String>& command, int i, String param) {
	if (i >= command.size())	
		throw Exc(Format(t_("Missing parameters when reading '%s'"), param));
}

void ConsoleMain(const Upp::Vector<String>& command, bool gui) {	
	String str = t_("BEMRosetta Copyright (c) 2019 Iñaki Zabala\nHydrodynamic coefficients converter for Boundary Element Method solver formats\nVersion beta BUILDINFO");
	SetBuildInfo(str);
	Cout() << str;
	
	BEMData md;
	
	bool firstTime = false;
	if (!md.LoadSerializeJson(firstTime))
		Cout() << "\n" << t_("BEM configuration data are not loaded. Defaults are set");
	
	String errorStr;
	try {
		if (command.IsEmpty()) {
			Cout() << "\n" << t_("Command argument list is empty");
			ShowHelp(md);
		} else {
			for (int i = 0; i < command.size(); i++) {
				if (command[i] == "-h" || command[i] == "--help") {
					ShowHelp(md);
					break;
				} else if (command[i] == "-i" || command[i] == "--input") {
					i++;
					CheckNumArgs(command, i, "--input");
					
					String file = command[i];
					if (!FileExists(file)) 
						throw Exc(Format(t_("File '%s' not found"), file)); 
					
					md.Load(file, [&](String str, int) {Cout() << str; return true;}, true);
					Cout() << "\n" << Format(t_("File '%s' loaded"), file);
				} else if (command[i] == "-r" || command[i] == "--report") {
					if (md.hydros.IsEmpty()) 
						throw Exc(t_("No file loaded"));
					int lastId = md.hydros.size() - 1;
					md.hydros[lastId].hd().Report();
				} else if (command[i] == "-cl" || command[i] == "--clear") {
					md.hydros.Clear();
					Cout() << "\n" << t_("Series cleared");
				} else if (command[i] == "-c" || command[i] == "--convert") {
					if (md.hydros.IsEmpty()) 
						throw Exc(t_("No file loaded"));
					i++;
					CheckNumArgs(command, i, "--convert");
					
					String file = command[i];
					
					md.hydros[0].hd().SaveAs(file);
					Cout() << "\n" << Format(t_("File '%s' converted"), file);
				} else if (command[i] == "-p" || command[i] == "--params") {
					i++;
					CheckNumArgs(command, i, "--params");
					
					while (i < command.size()) {
						if (command[i] == "g") {
							i++;
							CheckNumArgs(command, i, "-p g");
							double g = ScanDouble(command[i]);
							if (IsNull(g))
								throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
							md.g = g;
						} else if (command[i] == "length") {
							i++;
							CheckNumArgs(command, i, "-p length");
							double len = ScanDouble(command[i]);
							if (IsNull(len))
								throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
							md.len = len;
						} else if (command[i] == "rho") {
							i++;
							CheckNumArgs(command, i, "-p rho");
							double rho = ScanDouble(command[i]);
							if (IsNull(rho))
								throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
							md.rho = rho;
						} else if (command[i] == "depth") {
							i++;
							CheckNumArgs(command, i, "-p depth");
							double depth = ScanDouble(command[i]);
							if (IsNull(depth))
								throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
							md.depth = depth;
						} /*else if (command[i] == "thres") {
							i++;
							CheckNumArgs(command, i, "-p thres");
							double thres = ScanDouble(command[i]);
							if (IsNull(thres))
								throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
							md.thres = thres;
						} */else 
							throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
					}
				} else 
					throw Exc(Format(t_("Unknown argument '%s'"), command[i]));
			}
		}
	} catch (Exc e) {
		errorStr = e;
	} catch(const char *cad) {
		errorStr = cad;
	} catch(const std::string &e) {
		errorStr = e.c_str();	
	} catch (const std::exception &e) {
		errorStr = e.what();
	} catch(...) {
		errorStr = t_("Unknown error");
	}	
	if (!errorStr.IsEmpty()) {
		Cerr() << Format("\n%s: %s", t_("Error"), errorStr);
		Cerr() << S("\n\n") + t_("In case of doubt try option -h or --help");
		if (gui)
			Cerr() << S("\n") + t_("or just call command line without arguments to open GUI window");
	}
	Cout() << "\n";
}

void SetBuildInfo(String &str) {
	String name, mode;
	Time date;
	int version, bits;
	GetCompilerInfo(name, version, date, mode, bits);
	str.Replace("BUILDINFO", Format("%4d%02d%02d%02d, %s, %d bits", 
				date.year, date.month, date.day, date.hour, mode, bits)); 
}