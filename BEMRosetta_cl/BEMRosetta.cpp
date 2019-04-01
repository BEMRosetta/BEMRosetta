#include "BEMRosetta.h"


Function <void(String)> Hydro::Print = [](String s) {Cout() << s;};
Function <void(String)> Hydro::PrintError = [](String s) {Cout() << s;};

const char *Hydro::strDOF[] 	 = {"surge", "sway", "heave", "roll", "pitch", "yaw"};
const char *Hydro::strDOFAbrev[] = {"s", "w", "h", "r", "p", "y"};
	
void Hydro::Initialize_Forces() {
	Initialize_Forces(ex);
	Initialize_Forces(sc);
	Initialize_Forces(fk);
}

void Hydro::Initialize_Forces(Forces &f) {
	f.ma.SetCount(Nh);
	f.ph.SetCount(Nh);
	f.re.SetCount(Nh);
	f.im.SetCount(Nh);
	for (int ih = 0; ih < Nh; ++ih) {
		f.ma[ih].setConstant(Nf, Nb*6, nan(""));
		f.ph[ih].setConstant(Nf, Nb*6, nan(""));
		f.re[ih].setConstant(Nf, Nb*6, nan(""));
		f.im[ih].setConstant(Nf, Nb*6, nan(""));
	}
}

void Hydro::Initialize_RAO() {
	Initialize_Forces(rao);
}

void Hydro::Normalize() {
	if (IsLoadedC()) {
		for (auto &c : C) 
			for (int i = 0; i < 6; ++i) 
				for (int j = 0; j < 6; ++j) 
					c(i, j) /= (g*rho*pow(len, GetK_C(i, j)));
	}
	if (IsLoadedA() && IsLoadedB()) {
		for (int f = 0; f < Nf; ++f) {
			MatrixXd &a = A[f];
			MatrixXd &b = B[f];
			for (int i = 0; i < 6*Nb; ++i) {
				for (int j = 0; j < 6*Nb; ++j) {	
					int k = GetK_AB(i, j);
					a(i, j) /= (rho*pow(len, k));
					b(i, j) /= (rho*pow(len, k)*w[f]);
				}
			}
		}
	}
	if (IsLoadedAwinf()) {
		for (int i = 0; i < 6*Nb; ++i) 
			for (int j = 0; j < 6*Nb; ++j) 
				Awinf(i, j) /= (rho*pow(len, GetK_AB(i, j)));
	}
	if (IsLoadedAw0()) {
		for (int i = 0; i < 6*Nb; ++i) 
			for (int j = 0; j < 6*Nb; ++j) 
				Aw0(i, j) /= (rho*pow(len, GetK_AB(i, j)));
	}
	if (IsLoadedFex())
    	Normalize_Forces(ex);
	if (IsLoadedFsc())
		Normalize_Forces(sc);
	if (IsLoadedFfk())
		Normalize_Forces(fk);
	double A = 1;
	if (IsLoadedRAO()) {
		for (int h = 0; h < Nh; ++h) {
			for (int ifr = 0; ifr < Nf; ++ifr) {
				for (int i = 0; i < 6*Nb; ++i) {	 
					rao.ma[h](ifr, i) /= A/pow(len, GetK_RAO(i));
					rao.re[h](ifr, i) /= A/pow(len, GetK_RAO(i));
					rao.im[h](ifr, i) /= A/pow(len, GetK_RAO(i));
				}
			}
		}
	}
}

void Hydro::Dimensionalize() {
	if (IsLoadedC()) {
		for (auto &c : C) 
			for (int i = 0; i < 6; ++i) 
				for (int j = 0; j < 6; ++j) 
					c(i, j) *= (g*rho*pow(len, GetK_C(i, j)));
	}
	if (IsLoadedA() && IsLoadedB()) {
		for (int f = 0; f < Nf; ++f) {
			MatrixXd &a = A[f];
			MatrixXd &b = B[f];
			for (int i = 0; i < 6*Nb; ++i) {
				for (int j = 0; j < 6*Nb; ++j) {	
					int k = GetK_AB(i, j);
					a(i, j) *= (rho*pow(len, k));
					b(i, j) *= (rho*pow(len, k)*w[f]);
				}
			}
		}
	}
	if (IsLoadedAwinf()) {	
		for (int i = 0; i < 6*Nb; ++i) 
			for (int j = 0; j < 6*Nb; ++j) 
				Awinf(i, j) *= (rho*pow(len, GetK_AB(i, j)));
	}
	if (IsLoadedAw0()) {	
		for (int i = 0; i < 6*Nb; ++i) 
			for (int j = 0; j < 6*Nb; ++j) 
				Aw0(i, j) *= (rho*pow(len, GetK_AB(i, j)));
	}
	if (IsLoadedFex())
    	Dimensionalize_Forces(ex);
	if (IsLoadedFsc())
		Dimensionalize_Forces(sc);
	if (IsLoadedFfk())
		Dimensionalize_Forces(fk);
	double A = 1;
	if (IsLoadedRAO()) {
		for (int h = 0; h < Nh; ++h) {
			for (int ifr = 0; ifr < Nf; ++ifr) {
				for (int i = 0; i < 6*Nb; ++i) {	 
					rao.ma[h](ifr, i) *= A/pow(len, GetK_RAO(i));
					rao.re[h](ifr, i) *= A/pow(len, GetK_RAO(i));
					rao.im[h](ifr, i) *= A/pow(len, GetK_RAO(i));
				}
			}
		}
	}
}

void Hydro::Normalize_Forces(Forces &f) {
	for (int h = 0; h < Nh; ++h) {
		for (int ifr = 0; ifr < Nf; ++ifr) {
			for (int i = 0; i < 6*Nb; ++i) {	 
				f.ma[h](ifr, i) /= g*rho*pow(len, GetK_F(i));
				f.re[h](ifr, i) /= g*rho*pow(len, GetK_F(i));
				f.im[h](ifr, i) /= g*rho*pow(len, GetK_F(i));
			}
		}
	} 
}

void Hydro::Dimensionalize_Forces(Forces &f) {
	for (int h = 0; h < Nh; ++h) {
		for (int ifr = 0; ifr < Nf; ++ifr) {
			for (int i = 0; i < 6*Nb; ++i) {	 
				f.ma[h](ifr, i) *= g*rho*pow(len, GetK_F(i));
				f.re[h](ifr, i) *= g*rho*pow(len, GetK_F(i));
				f.im[h](ifr, i) *= g*rho*pow(len, GetK_F(i));
			}
		}
	}
}

void Hydro::SaveAs(String file, BEM_SOFT type) {
	if (type == WAMIT_1_3) {
		Fast data(this);
		data.Save(file, false);	
	} else if (type == FAST_WAMIT) {
		Fast data(this);
		data.Save(file, true);		
	}
}

void Hydro::Report() {
	Print(Format("\n%s file '%s'", GetCodeStr(), file));
	Print(Format("\ng [m/s2]: %.3f, h [m]: %.3f, rho [kg/m3]: %.3f length scale [m]: %.1f", g, h, rho, len));
	String freqs;
	if (w.IsEmpty()) 
		freqs = "NONE";
	else
		freqs = Format("%.3f to %.3f steps %.3f [rad/s]", w[0], w[w.GetCount()-1], w[1]-w[0]);
	String heads;
	if (head.IsEmpty())
		heads = "NONE";
	else if (head.GetCount() > 1)
	 	heads = Format("%.1f to %.1f steps %.1f [º]", head[0], head[head.GetCount()-1], head[1] - head[0]);	
	else
		heads = Format("%.1f [º]", head[0]);
		
	Print(Format("\n#freqs: %d (%s), #headings: %d (%s)", Nf, freqs, Nh, heads)); 
	Print(Format("\n#bodies: %d", Nb));
	for (int ib = 0; ib < Nb; ++ib) {
		String str = Format("\n%d.", ib+1);
		if (names.GetCount() > ib)
			str += " '" + names[ib] + "'";
		if (dof.GetCount() > ib)
			str += " dof: " + dof[ib];
		if (Vo.size() > ib)
			str += " vol [m3]: " + FormatDouble(Vo[ib]);
		if (cg.size() > 3*ib)
			str += Format(" cg(%.3f, %.3f, %.3f)[m]", cg(0, ib), cg(1, ib), cg(2, ib));
		if (cb.size() > 3*ib)
			str += Format(" cb(%.3f, %.3f, %.3f)[m]", cb(0, ib), cb(1, ib), cb(2, ib));
		
		Print(str);
	}
}

bool HydroClass::MatchCoeffStructure(Upp::Array<HydroClass> &hydro, String &strError) {
	strError.Clear();
	if (hydro.IsEmpty()) {
		strError = "No data loaded";
		return false;
	}
	int Nb = hydro[0].hd().Nb;
	int Nh = hydro[0].hd().Nh;
	for (int i = 1; i < hydro.GetCount(); ++i) {
		if (hydro[i].hd().Nb != Nb) {
			strError = "Different number of bodies";
			return false;
		} else if (hydro[i].hd().Nh != Nh) {
			strError = "Different number of wave headings";
			return false;
		} else {
			for (int ih = 0; ih < Nh; ++ih) {
				if (hydro[i].hd().head[ih] != hydro[0].hd().head[ih]) {
					strError = "Wave headings do not match";
					return false;
				}
			}
		}
	}
	return true;
}

void Hydro::GetBodyDOF() {
	dof.Clear();	 dof.SetCount(Nb, 0);
	for (int ib = 0; ib < Nb; ++ib)
		for (int idof = 0; idof < 6; ++idof)
			if (IsAvailableDOF(ib, idof))
				dof[ib]++;
}

void Hydro::AfterLoad() {
	dofOrder.SetCount(6*Nb);
	for (int i = 0, order = 0; i < 6*Nb; ++i, ++order) {
		//if (order >= 6)
		//	order = 0;
		dofOrder[i] = order;
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

int IsTabSpace(int c) {
	if (c == '\t' || c == ' ' || c == '!')
		return true;
	return false;
}

String FormatWam(double d) {
	return (d >= 0 ? " " : "-") + Format("%12E", abs(d));
}
