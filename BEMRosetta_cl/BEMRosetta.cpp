#include "BEMRosetta.h"


Function <void(String)> Hydro::Print 		= [](String s) {Cout() << s;};
Function <void(String)> Hydro::PrintWarning = [](String s) {Cout() << s;};
Function <void(String)> Hydro::PrintError 	= [](String s) {Cout() << s;};

const char *Hydro::strDOF[] 	 = {t_("surge"), t_("sway"), t_("heave"), t_("roll"), t_("pitch"), t_("yaw")};
const char *Hydro::strDOFAbrev[] = {t_("s"), t_("w"), t_("h"), t_("r"), t_("p"), t_("y")};
	
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
		f.ma[ih].setConstant(Nf, Nb*6, Null);
		f.ph[ih].setConstant(Nf, Nb*6, Null);
		f.re[ih].setConstant(Nf, Nb*6, Null);
		f.im[ih].setConstant(Nf, Nb*6, Null);
	}
}

void Hydro::GetFexFromFscFfk() {
	for (int ih = 0; ih < Nh; ++ih) {
		for (int ifr = 0; ifr < Nf; ++ifr) {
			for (int i = 0; i < Nb*6; ++i) {
				if (!IsNull(sc.ma[ih](ifr, i))) {
					double scm = sc.ma[ih](ifr, i);
					double scp = sc.ph[ih](ifr, i);
					double fkm = fk.ma[ih](ifr, i);
					double fkp = fk.ph[ih](ifr, i);
					double exre = scm*cos(scp) + fkm*cos(fkp);
					double exim = scm*sin(scp) + fkm*sin(fkp);
					ex.re[ih](ifr, i) = exre;
					ex.im[ih](ifr, i) = exim;
					ex.ma[ih](ifr, i) = sqrt(exre*exre + exim*exim);
					ex.ph[ih](ifr, i) = atan2(exre, exim);
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

void Hydro::RemoveThresDOF_A(double thres) {
	if (!IsLoadedA())
		return;
	for (int idof = 0; idof < 6*Nb; ++idof) {
		for (int jdof = 0; jdof < 6*Nb; ++jdof) {
			double mx = 0;
			for (int ifr = 0; ifr < Nf; ifr++) 
				mx = max(mx, abs(A[ifr](idof, jdof)));
			if (mx < thres) {
				for (int ifr = 0; ifr < Nf; ifr++) 
					A[ifr](idof, jdof) = Null;	
			}
		}
	}
}

void Hydro::RemoveThresDOF_B(double thres) {
	if (!IsLoadedB())
		return;
	for (int idof = 0; idof < 6*Nb; ++idof) {
		for (int jdof = 0; jdof < 6*Nb; ++jdof) {
			double mx = 0;
			for (int ifr = 0; ifr < Nf; ifr++) 
				mx = max(mx, abs(B[ifr](idof, jdof)));
			if (mx < thres) {
				for (int ifr = 0; ifr < Nf; ifr++) 
					B[ifr](idof, jdof) = Null;	
			}
		}
	}
}

void Hydro::RemoveThresDOF_Force(Forces &f, double thres) {
	if (!IsLoadedForce(f))
		return;
	for (int h = 0; h < Nh; ++h) {
		for (int i = 0; i < 6*Nb; ++i) {	 
			double mx = 0;
			for (int ifr = 0; ifr < Nf; ++ifr) 
				mx = max(mx, f.ma[h](ifr, i));
			if (mx < thres) {
				for (int ifr = 0; ifr < Nf; ifr++) 
					f.ma[h](ifr, i) = Null;	
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
		if (abs((a.w[i] - w[i])/w[i]) > 0.0001)
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
		for (int idof = 0; idof < 6*a.Nb; ++idof) {
			for (int jdof = 0; jdof < 6*a.Nb; ++jdof) {
				double Aa = a.A[ifr](idof, jdof);
				double Ab = A[ifr](idof, jdof);
				if (!IsNull(Aa) && !IsNull(Ab) && Aa != Ab)
					throw Exc(Format(t_("%s is not the same %f<>%f"), 
							Format(t_("%s[%d](%d, %d)"), t_("A"), ifr+1, idof+1, jdof+1), 
							Aa, Ab));
			}
		}
	}
}

void Hydro::Compare_B(Hydro &a) {
	for (int ifr = 0; ifr < a.Nf; ifr++) {
		for (int idof = 0; idof < 6*a.Nb; ++idof) {
			for (int jdof = 0; jdof < 6*a.Nb; ++jdof) {
				double Ba = a.B[ifr](idof, jdof);
				double Bb = B[ifr](idof, jdof);
				if (!IsNull(Ba) && !IsNull(Bb) && Ba != Bb)
					throw Exc(Format(t_("%s is not the same %f<>%f"), 
							Format(t_("%s[%d](%d, %d)"), t_("B"), ifr+1, idof+1, jdof+1), 
							Ba, Bb));
			}
		}
	}
}

void Hydro::Compare_C(Hydro &a) {
	for (int ib = 0; ib < a.Nb; ib++) {
		for (int idof = 0; idof < 6; ++idof) {
			for (int jdof = 0; jdof < 6; ++jdof) {
				double Ca = a.C[ib](idof, jdof);
				double Cb = C[ib](idof, jdof);
				if (!IsNull(Ca) && !IsNull(Cb) && abs((Ca-Cb)/Cb) > 0.0001 )
					throw Exc(Format(t_("%s is not the same %f<>%f"), 
							Format(t_("%s[%d](%d, %d)"), t_("C"), ib+1, idof+1, jdof+1), 
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

void Hydro::SaveAs(String file, BEM_SOFT type) {
	int realNh = Nh;
	int realNf = Nf;
	
//	int lastHead = GetIrregularHead();
//	if (lastHead >= 0) {
//		Nh = lastHead+1;
//		Print("\n" + t_("Head spread is not homogeneous, so it has been reduced for Wamit"));
//	}
	
//	int lastFreq = GetIrregularFreq();
//	if (lastFreq >= 0) {
//		Nf = lastFreq+1;
//		Print("\n" + t_("Frequency spread is not homogeneous, so it has been reduced for Wamit"));
//	}
	
	if (type == UNKNOWN) {
		String ext = ToLower(GetFileExt(file));
		
		if (ext == ".1" || ext == ".3" || ext == ".hst")
			type = Hydro::WAMIT_1_3;
		else if (ext == ".dat")
			type = Hydro::FAST_WAMIT;	
		else
			throw Exc(Format(t_("Conversion to type of file '%s' not supported"), file));
	}
	if (type == WAMIT_1_3) {
		Wamit data(this);
		data.Save(file);	
	} else if (type == FAST_WAMIT) {
		Fast data(this);
		data.Save(file);		
	}
	Nh = realNh;
	Nf = realNf;
}

void Hydro::Report() {
	Print("\n" + Format(t_("%s file '%s'"), GetCodeStr(), file));
	Print("\n" + Format(t_("g [m/s2]: %.3f, h [m]: %s, rho [kg/m3]: %.3f length scale [m]: %.1f"), g, h < 0 ? x_(t_("INFINITY")) : FormatDouble(h), rho, len));
	String freqs;
	if (w.IsEmpty()) 
		freqs = t_("NONE");
	else if (w.GetCount() > 1) {
		String strDeltaH;
		if (GetIrregularFreq() < 0) 
			strDeltaH = Format(t_("delta %.1f [rad/s]"), w[1] - w[0]);
		else {
			String strHead;
			for (int i = 0; i < w.GetCount(); ++i) {
				if (i > 0)
					strHead << ", ";
				strHead << w[i];
			}
			strDeltaH = Format(t_("non constant delta (%s)"), strHead); 
		}
	 	freqs = Format(t_("%.1f to %.1f %s"), w[0], w[w.GetCount()-1], strDeltaH);	
	} else
		freqs = Format(t_("%.1f [rad/s]"), w[0]);
	
	String heads;
	if (head.IsEmpty())
		heads = t_("NONE");
	else if (head.GetCount() > 1) {
		String strDeltaH;
		if (GetIrregularHead() < 0) 
			strDeltaH = Format(t_("delta %.1f [º]"), head[1] - head[0]);
		else {
			String strHead;
			for (int i = 0; i < head.GetCount(); ++i) {
				if (i > 0)
					strHead << ", ";
				strHead << head[i];
			}
			strDeltaH = Format(t_("non constant delta (%s)"), strHead); 
		}
	 	heads = Format(t_("%.1f to %.1f %s"), head[0], head[head.GetCount()-1], strDeltaH);	
	} else
		heads = Format(t_("%.1f [º]"), head[0]);
	
	Print("\n" + Format(t_("#freqs: %d (%s)"), Nf, freqs)); 
	Print("\n" + Format(t_("#headings: %d (%s)"), Nh, heads)); 
	Print("\n" + Format(t_("#bodies: %d"), Nb));
	for (int ib = 0; ib < Nb; ++ib) {
		String str = Format("\n%d.", ib+1);
		if (names.GetCount() > ib)
			str += " '" + names[ib] + "'";
		if (dof.GetCount() > ib)
			str += x_(" ") + t_("dof") + ": " + FormatInt(dof[ib]);
		if (Vo.size() > ib && !IsNull(Vo[ib]))
			str += x_(" ") + t_("vol [m3]") + ": " + FormatDouble(Vo[ib]);
		if (cg.size() > 3*ib && !IsNull(cg(0, ib)))
			str += " " + Format("Cg(%.3f, %.3f, %.3f)[m]", cg(0, ib), cg(1, ib), cg(2, ib));
		if (cb.size() > 3*ib && !IsNull(cb(0, ib)))
			str += " " + Format("Cb(%.3f, %.3f, %.3f)[m]", cb(0, ib), cb(1, ib), cb(2, ib));
		
		Print(str);
	}
}

bool HydroClass::MatchCoeffStructure(Upp::Array<HydroClass> &hydro, String &strError) {
	strError.Clear();
	if (hydro.IsEmpty()) {
		strError = t_("No data loaded");
		return false;
	}
	int Nb = hydro[0].hd().Nb;
	int Nh = hydro[0].hd().Nh;
	for (int i = 1; i < hydro.GetCount(); ++i) {
		if (hydro[i].hd().Nb != Nb) {
			strError = t_("Different number of bodies");
			return false;
		} else if (hydro[i].hd().Nh != Nh) {
			strError = t_("Different number of wave headings");
			return false;
		} else {
			for (int ih = 0; ih < Nh; ++ih) {
				if (hydro[i].hd().head[ih] != hydro[0].hd().head[ih]) {
					strError = t_("Wave headings do not match");
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

int Hydro::GetIrregularHead() {
	if (head.GetCount() <= 2)
		return -1;
	double delta0 = head[1] - head[0];
	for (int i = 1; i < head.GetCount() - 1; ++i) {
		double delta = head[i+1] - head[i];
		if (abs((delta - delta0)/delta0) > 0.001)
			return i;
	}
	return -1;
}

int Hydro::GetIrregularFreq() {
	if (w.GetCount() <= 2)
		return -1;
	double delta0 = w[1] - w[0];
	for (int i = 1; i < w.GetCount() - 1; ++i) {
		double delta = w[i+1] - w[i];
		if (abs((delta - delta0)/delta0) > 0.001)
			return i;
	}
	return -1;
}


void BEMData::Load(String file, Function <void(BEMData &, HydroClass&)> AdditionalData) {
	for (int i = 0; i < hydros.GetCount(); ++i) {
		if (hydros[i].hd().file == file) 
			throw Exc(Format(t_("Model '%s' already loaded"), file));
	}
	String ext = ToLower(GetFileExt(file));
	if (ext == ".cal") {
		Nemoh &data = hydros.Create<Nemoh>();
		if (!data.Load(file)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.GetCount()-1);
			throw Exc(Format(t_("Problem loading '%s'\n%s"), file, error));	
		}
	} else if (ext == ".inf") {
		Nemoh &data = hydros.Create<Nemoh>();
		if (!data.Load(file)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.GetCount()-1);
			throw Exc(Format(t_("Problem loading '%s'\n%s"), file, error));	
		}
	} else if (ext == ".out") {
		Wamit &data = hydros.Create<Wamit>();
		if (!data.Load(file, Null)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.GetCount()-1);
			throw Exc(Format(t_("Problem loading '%s'") + x_("\n%s"), file, error));
		}
		AdditionalData(*this, data);
		
		data.hd().Dimensionalize();
	} else if (ext == ".dat") {
		Fast &data = hydros.Create<Fast>();
		if (!data.Load(file, Null)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.GetCount()-1);
			throw Exc(Format(t_("Problem loading '%s'") + x_("\n%s"), file, error));		
		}
		AdditionalData(*this, data);
				
		data.hd().Dimensionalize();
	} else if (ext == ".1" || ext == ".3" || ext == ".hst" || ext == ".4") {
		Wamit &data = hydros.Create<Wamit>();
		if (!data.Load(file, Null)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.GetCount()-1);
			throw Exc(Format(t_("Problem loading '%s'") + x_("\n%s"), file, error));		
		}
		AdditionalData(*this, data);

		data.hd().Dimensionalize();		
	} else if (ext == ".ah1" || ext == ".lis") {
		Aqwa &data = hydros.Create<Aqwa>();
		if (!data.Load(file)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.GetCount()-1);
			throw Exc(Format(t_("Problem loading '%s'\n%s"), file, error));	
		}
	} else if (ext == ".mat") {
		Foamm &data = hydros.Create<Foamm>();
		if (!data.Load(file)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.GetCount()-1);
			throw Exc(Format(t_("Problem loading '%s'\n%s"), file, error));	
		}
		AdditionalData(*this, data);
	} else 
		throw Exc(Format(t_("Unknown file extension in '%s'"), file));
	
	/*Hydro &justLoaded = hydros[hydros.GetCount()-1].hd();
	justLoaded.RemoveThresDOF_A(thres);
	justLoaded.RemoveThresDOF_B(thres);
	justLoaded.RemoveThresDOF_Force(justLoaded.ex, thres);
	justLoaded.RemoveThresDOF_Force(justLoaded.sc, thres);
	justLoaded.RemoveThresDOF_Force(justLoaded.fk, thres);
	justLoaded.RemoveThresDOF_Force(justLoaded.rao, thres);*/
}

int IsTabSpace(int c) {
	if (c == '\t' || c == ' ' || c == '!')
		return true;
	return false;
}

String FormatWam(double d) {
	return (d >= 0 ? " " : "-") + Format("%12E", abs(d));
}
