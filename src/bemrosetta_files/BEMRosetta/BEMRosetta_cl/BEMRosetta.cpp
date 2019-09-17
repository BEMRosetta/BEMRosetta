#include "BEMRosetta.h"

#include <plugin/matio/matio.h>

Function <void(String)> BEMData::Print 		  = [](String s) {Cout() << s;};
Function <void(String)> BEMData::PrintWarning = [](String s) {Cout() << s;};
Function <void(String)> BEMData::PrintError   = [](String s) {Cout() << s;};

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
					A[ifr](idf, jdf) = A_ndim(ifr, idf, jdf);
					B[ifr](idf, jdf) = B_ndim(ifr, idf, jdf);
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
					rao.ma[h](ifr, i) = F_ma_ndim(rao, h, ifr, i);
					rao.re[h](ifr, i) = F_re_ndim(rao, h, ifr, i);
					rao.im[h](ifr, i) = F_im_ndim(rao, h, ifr, i);
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
					A[ifr](idf, jdf) = A_dim(ifr, idf, jdf);
					B[ifr](idf, jdf) = B_dim(ifr, idf, jdf);
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
					rao.ma[h](ifr, i) = F_ma_dim(rao, h, ifr, i);
					rao.re[h](ifr, i) = F_re_dim(rao, h, ifr, i);
					rao.im[h](ifr, i) = F_im_dim(rao, h, ifr, i);
				}
			}
		}
	}
}

void Hydro::Normalize_Forces(Forces &f) {
	for (int h = 0; h < Nh; ++h) {
		for (int ifr = 0; ifr < Nf; ++ifr) {
			for (int idf = 0; idf < 6*Nb; ++idf) {	 
				f.ma[h](ifr, idf) = F_ma_dim(f, h, ifr, idf);
				f.re[h](ifr, idf) = F_re_dim(f, h, ifr, idf);
				f.im[h](ifr, idf) = F_im_dim(f, h, ifr, idf);
			}
		}
	} 
}

void Hydro::Dimensionalize_Forces(Forces &f) {
	for (int h = 0; h < Nh; ++h) {
		for (int ifr = 0; ifr < Nf; ++ifr) {
			for (int idf = 0; idf < 6*Nb; ++idf) {	 
				f.ma[h](ifr, idf) = F_ma_dim(f, h, ifr, idf);
				f.re[h](ifr, idf) = F_re_dim(f, h, ifr, idf);
				f.im[h](ifr, idf) = F_im_dim(f, h, ifr, idf);
			}
		}
	}
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
						A[ifr](idf, jdf) = Null;
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
						B[ifr](idf, jdf) = Null;
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
				double Aa = a.A[ifr](idf, jdf);
				double Ab = A[ifr](idf, jdf);
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
				double Ba = a.B[ifr](idf, jdf);
				double Bb = B[ifr](idf, jdf);
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
		Wamit data(*bem, this);
		data.Save(file);	
	} else if (type == FAST_WAMIT) {
		Fast data(*bem, this);
		data.Save(file);		
	}
	Nh = realNh;
	Nf = realNf;
}

void Hydro::GetFOAMM(String file, Function <bool(String)> Running) {
	MatFile mat;
	
	if (!mat.OpenWrite("C:\\Desarrollo\\Aplicaciones\\BEMRosetta\\FOAMM_Win\\pru\\temp_file.mat")) 
		throw Exc(Format(t_("Problem launching FOAMM from '%s'"), file));
	
	mat.VarDelete("Options");
	
	Vector<String> optionsVars;
	optionsVars << "Mode" << "Method" << "FreqRangeChoice" << "FreqChoice";
	MatVar options("Options", 1, 1, optionsVars);
	
	options.VarWriteStruct<double>("Mode", 0);
	options.VarWriteStruct<double>("Method", 0);
	options.VarWriteStruct("FreqRangeChoice", "G");
	options.VarWriteStruct("FreqChoice", "G");
	
	Vector<String> optimVars;
	optimVars << "InitCond" << "Tol" << "maxEval" << "maxIter" << "StepTol" << "ThresRel" << "ThresAbs";
	MatVar optim("Optim", 1, 1, optimVars);
	optim.VarWriteStruct<double>("InitCond", 50);
	optim.VarWriteStruct<double>("Tol", 1E-5);
	optim.VarWriteStruct<double>("maxEval", 1E3);
	optim.VarWriteStruct<double>("maxIter", 200);
	optim.VarWriteStruct<double>("StepTol", 1E-6);
	optim.VarWriteStruct<double>("ThresRel", 0.03);
	optim.VarWriteStruct<double>("ThresAbs", 0.1);
	
	options.VarWriteStruct("Optim", optim);
	mat.VarWrite(options, true);
	
	mat.Close();
	
	
	
	
	
	
	
	
	
	LocalProcess process;
	if (!process.Start(bem->foammPath, NULL, "C:\\Desarrollo\\Aplicaciones\\BEMRosetta\\FOAMM_Win\\pru"))
		throw Exc(Format(t_("Problem launching FOAMM from '%s'"), file));

	String msg, reso, rese;
	bool endProcess = false;
	while (process.IsRunning()) {
		if (!endProcess) {
			msg.Clear();
			if (process.Read2(reso, rese)) {
				if (!reso.IsEmpty())
					msg << reso;
				if (!rese.IsEmpty()) {
					if (!msg.IsEmpty())
						msg << ";";
					msg << rese;
				}
				/*if (!msg.IsEmpty()) {
					process.Kill();
					Running(msg);
					return;
				}*/
			}
			endProcess = Running(msg); 
			if (endProcess)
				process.Kill();
		}
		Sleep(200);
	}
	
	if (!mat.OpenRead("C:\\Desarrollo\\Aplicaciones\\BEMRosetta\\FOAMM_Win\\pru\\temp_file.mat")) 
		throw Exc(Format(t_("Problem launching FOAMM from '%s'"), file));
	Cout() << "\nData saved" << "\n";
	Cout() << "Version: " << mat.GetVersionName() << "\n";
	for (int i = 0; i < mat.GetVarCount(); ++i) {
		String name = mat.GetVarName(i);
		MatVar var = mat.GetVar(name);
		Cout() << name << " (" << var.GetTypeString() << ")";
		for (int d = 0; d < var.GetDimCount(); ++d)
			Cout() << "(" << var.GetDimCount(d) << ")";
		Cout() << "\n";
		if (var.GetType() == MAT_C_STRUCT) {
			for (int s = 0; s < var.GetFieldCount(); ++s) {
				Cout() << "\t" << var.GetFieldName(s) << "\n";
			}
		} 
		
	}
}

void Hydro::Report() {
	BEMData::Print("\n" + Format(t_("%s file '%s'"), GetCodeStr(), file));
	String sg   = IsNull(g)   ? x_("unknown") : Format("%.3f", g);
	String srho = IsNull(rho) ? x_("unknown") : Format("%.3f", rho);
	String slen = IsNull(len) ? x_("unknown") : Format("%.1f", len);
	BEMData::Print("\n" + Format(t_("g [m/s2]: %s, h [m]: %s, rho [kg/m3]: %s, length scale [m]: %s"), 
								sg, h < 0 ? x_(t_("INFINITY")) : FormatDouble(h), srho, slen));
	String freqs;
	if (w.IsEmpty()) 
		freqs = t_("NONE");
	else if (w.GetCount() > 1) {
		String strDeltaH;
		if (GetIrregularFreq() < 0) 
			strDeltaH = Format(t_("delta %s [rad/s]"), FormatDouble(w[1] - w[0], 5, FD_EXP));
		else {
			String strHead;
			for (int i = 0; i < w.GetCount(); ++i) {
				if (i > 0)
					strHead << ", ";
				strHead << w[i];
			}
			strDeltaH = Format(t_("non constant delta (%s)"), strHead); 
		}
	 	freqs = Format(t_("%s to %s %s"), FormatDouble(w[0], 3, FD_EXP), 
	 									  FormatDouble(w[w.GetCount()-1], 3, FD_EXP), strDeltaH);	
	} else
		freqs = Format(t_("%s [rad/s]"), FormatDouble(w[0], 3, FD_EXP));
	
	String heads;
	if (head.IsEmpty())
		heads = t_("NONE");
	else if (head.GetCount() > 1) {
		String strDeltaH;
		if (GetIrregularHead() < 0) 
			strDeltaH = Format(t_("delta %.1f [deg]"), head[1] - head[0]);
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
		heads = Format(t_("%.1f [deg]"), head[0]);
	
	BEMData::Print("\n" + Format(t_("#freqs: %d (%s)"), Nf, freqs)); 
	BEMData::Print("\n" + Format(t_("#headings: %d (%s)"), Nh, heads)); 
	BEMData::Print("\n" + Format(t_("#bodies: %d"), Nb));
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

void Hydro::AfterLoad(Function <void(String, int)> Status) {
	dofOrder.SetCount(6*Nb);
	for (int i = 0, order = 0; i < 6*Nb; ++i, ++order) 
		dofOrder[i] = order;
	
	if (!IsLoadedAw0())  
		A0();
	
	if (!IsLoadedAwinf() && bem->calcAwinf) {
		if (IsNull(bem->maxTimeA) || bem->maxTimeA == 0)
			throw Exc(t_("Incorrect time for Ainf calculation. Please review it in Options"));
		if (IsNull(bem->numValsA) || bem->numValsA < 10)
			throw Exc(t_("Incorrect number of time values for Ainf calculation. Please review it in Options"));

		Status(t_("Obtaining Impulse Response Function"), 40);
		K_IRF(bem->maxTimeA, bem->numValsA);
		Status(t_("Obtaining Infinite-Frequency Added Mass (A_inf)"), 70);
		Ainf();
	}
}

int Hydro::GetW0() {
	for (int i = 0; i < w.GetCount(); ++i) {	
		if (w[i] < 0.0001)
			return i;
	}
	return Null;
}

void Hydro::A0() {
	int iw0 = GetW0();
	if (IsNull(iw0)) 
		return;
	
	Aw0.setConstant(Nb*6, Nb*6, Null);
	for (int i = 0; i < Nb*6; ++i)
        for (int j = 0; j < Nb*6; ++j) 
			Aw0(i, j) = A[iw0](i, j);
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

void BEMData::Load(String file, Function <void(String, int pos)> Status) {
	Status(t_("Loading files"), 10);
	for (int i = 0; i < hydros.GetCount(); ++i) {
		if (hydros[i].hd().file == file) 
			throw Exc(Format(t_("Model '%s' already loaded"), file));
	}
	String ext = ToLower(GetFileExt(file));
	if (ext == ".cal") {
		Nemoh &data = hydros.Create<Nemoh>(*this);
		if (!data.Load(file)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.GetCount()-1);
			throw Exc(Format(t_("Problem loading '%s'\n%s"), file, error));	
		}
	} else if (ext == ".inf") {
		Nemoh &data = hydros.Create<Nemoh>(*this);
		if (!data.Load(file)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.GetCount()-1);
			throw Exc(Format(t_("Problem loading '%s'\n%s"), file, error));	
		}
	} else if (ext == ".out") {
		Wamit &data = hydros.Create<Wamit>(*this);
		if (!data.Load(file)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.GetCount()-1);
			throw Exc(Format(t_("Problem loading '%s'") + x_("\n%s"), file, error));
		}
	} else if (ext == ".dat") {
		Fast &data = hydros.Create<Fast>(*this);
		if (!data.Load(file)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.GetCount()-1);
			throw Exc(Format(t_("Problem loading '%s'") + x_("\n%s"), file, error));		
		}
	} else if (ext == ".1" || ext == ".3" || ext == ".hst" || ext == ".4") {
		Wamit &data = hydros.Create<Wamit>(*this);
		if (!data.Load(file)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.GetCount()-1);
			throw Exc(Format(t_("Problem loading '%s'") + x_("\n%s"), file, error));		
		}
	} else if (ext == ".ah1" || ext == ".lis") {
		Aqwa &data = hydros.Create<Aqwa>(*this);
		if (!data.Load(file)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.GetCount()-1);
			throw Exc(Format(t_("Problem loading '%s'\n%s"), file, error));	
		}
	} else if (ext == ".mat") {
		Foamm &data = hydros.Create<Foamm>(*this);
		if (!data.Load(file)) {
			String error = data.hd().GetLastError();
			hydros.SetCount(hydros.GetCount()-1);
			throw Exc(Format(t_("Problem loading '%s'\n%s"), file, error));	
		}
	} else 
		throw Exc(Format(t_("Unknown file extension in '%s'"), file));
	
	Hydro &justLoaded = hydros[hydros.GetCount()-1].hd();
	
	justLoaded.AfterLoad(Status);
	if (discardNegDOF) {
		Status(t_("Discarding negligible DOF"), 90);
		justLoaded.RemoveThresDOF_A(thres);
		justLoaded.RemoveThresDOF_B(thres);
		justLoaded.RemoveThresDOF_Force(justLoaded.ex, thres);
		justLoaded.RemoveThresDOF_Force(justLoaded.sc, thres);
		justLoaded.RemoveThresDOF_Force(justLoaded.fk, thres);
		justLoaded.RemoveThresDOF_Force(justLoaded.rao, thres);
	}
	if (hydros.GetCount() == 1)
		Nb = justLoaded.Nb;
	else {
		if (justLoaded.Nb > Nb)
			throw Exc(Format(t_("Model has more bodies (%d) than previously loaded (%d)"), justLoaded.Nb, Nb));
	}
	for (int i = 0; i < justLoaded.head.GetCount(); ++i) 
		FindAddRatio(head, justLoaded.head[i], 0.01);
	Sort(head);
}

void BEMData::LoadMesh(String file, Function <void(String, int pos)> Status) {
	Status(Format(t_("Loaded mesh '%s'"), file), 10);
	
	for (int i = 0; i < surfs.GetCount(); ++i) {
		if (surfs[i].file == file) {
			BEMData::Print(x_("\n") + t_("Model already loaded"));
			throw Exc(t_("Model already loaded"));
		}
	}

	MeshData &mesh = surfs.Add();
	String error = mesh.Load(file, rho, g);
	if (!error.IsEmpty()) {
		BEMData::Print("\n" + Format(t_("Problem loading '%s'") + x_("\n%s"), file, error));
		surfs.Remove(surfs.GetCount()-1);
		throw Exc(Format(t_("Problem loading '%s'") + x_("\n%s"), file, error));
	}
}

void BEMData::HealingMesh(int id, Function <void(String, int pos)> Status) {
	Status(Format(t_("Healing mesh '%s'"), surfs[id].file), 10);
	Print(x_("\n\n") + Format(t_("Healing mesh '%s'"), surfs[id].file));
	
	String ret;
	try {
		ret = surfs[id].Heal(Status);
	} catch (Exc e) {
		surfs.SetCount(surfs.GetCount()-1);
		Print("\n" + Format(t_("Problem loading '%s': %s") + x_("\n%s"), e));
		throw e;
	}
	if (!ret.IsEmpty()) {
		ret.Replace("\n", "\n- ");
		Print(ret);
	} else
		Print(x_(". ") + t_("The mesh is in good condition"));
}

void BEMData::UnderwaterMesh(int id, Function <void(String, int pos)> Status) {
	Status(Format(t_("Getting underwater mesh '%s'"), surfs[id].file), 10);
	
	MeshData &mesh = surfs.Add();
	MeshData &orig = surfs[id];
	mesh.file = orig.file;
	
	String ret;
	try {
		mesh.mesh.Underwater(orig.mesh);
	} catch (Exc e) {
		surfs.SetCount(surfs.GetCount()-1);
		Print("\n" + Format(t_("Problem loading '%s': %s") + x_("\n%s"), e));
		throw e;
	}
}

int IsTabSpace(int c) {
	if (c == '\t' || c == ' ' || c == '!')
		return true;
	return false;
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
	Cout() << "\n" << t_("                    length length scale  []        ") << md.length;
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

void CheckNumArgs(const Vector<String>& command, int i, String param) {
	if (i >= command.GetCount())	
		throw Exc(Format(t_("Missing parameters when reading '%s'"), param));
}

void ConsoleMain(const Vector<String>& command, bool gui) {	
	String str = t_("BEMRosetta Copyright (c) 2019 Iñaki Zabala\nHydrodynamic coefficients converter for Boundary Element Method solver formats\nVersion beta BUILDINFO");
	Hydro::SetBuildInfo(str);
	Cout() << str;
	
	BEMData md;
	
	if (!md.LoadSerializeJson())
		Cout() << "\n" << t_("BEM configuration data are not loaded. Defaults are set");
	
	String errorStr;
	try {
		if (command.IsEmpty()) {
			Cout() << "\n" << t_("Command argument list is empty");
			ShowHelp(md);
		} else {
			for (int i = 0; i < command.GetCount(); i++) {
				if (command[i] == "-h" || command[i] == "--help") {
					ShowHelp(md);
					break;
				} else if (command[i] == "-i" || command[i] == "--input") {
					i++;
					CheckNumArgs(command, i, "--input");
					
					String file = command[i];
					if (!FileExists(file)) 
						throw Exc(Format(t_("File '%s' not found"), file)); 
					
					md.Load(file, [&](String str, int ) {Cout() << str;});
					Cout() << "\n" << Format(t_("File '%s' loaded"), file);
				} else if (command[i] == "-r" || command[i] == "--report") {
					if (md.hydros.IsEmpty()) 
						throw Exc(t_("No file loaded"));
					int lastId = md.hydros.GetCount() - 1;
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
					
					while (i < command.GetCount()) {
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
							double length = ScanDouble(command[i]);
							if (IsNull(length))
								throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
							md.length = length;
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
		Cerr() << x_("\n\n") + t_("In case of doubt try option -h or --help");
		if (gui)
			Cerr() << x_("\n") + t_("or just call command line without arguments to open GUI window");
	}
	Cout() << "\n";
}

bool OUTB(int id, int total) {
	if (id < 0	|| id >= total)
		return true;
	return false;
}