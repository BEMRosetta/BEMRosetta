// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"


bool Diodore::Load(String file, double) {
	hd().file = file;
	hd().name = GetFileTitle(file);
	hd().dimen = true;
	hd().len = 1;
	hd().code = Hydro::DIODORE;
	hd().Nb = Null;
	
	try {
		BEM::Print("\n\n" + Format(t_("Loading '%s'"), file));

		BEM::Print("\n- " + S(t_("HDB file")));
		if (!Load_HDB()) 
			BEM::PrintWarning(S(": ** HDB file ") + t_("Not found") + "**");
		
		if (IsNull(hd().Nb))
			return false;
	
		hd().dof.Clear();	hd().dof.SetCount(hd().Nb, 0);
		for (int i = 0; i < hd().Nb; ++i)
			hd().dof[i] = 6;
	} catch (Exc e) {
		BEM::PrintError(Format("\n%s: %s", t_("Error"), e));
		hd().lastError = e;
		return false;
	}
	
	return true;
}

bool Diodore::Load_HDB() {
	String fileName = ForceExt(hd().file, ".hdb");
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		return false;
	
	String line; 
	FieldSplit f(in);
	f.IsSeparator = IsTabSpace;
	
	hd().dataFromW = false;
	hd().dimen = true;
	
	hd().Nb = hd().Nf = hd().Nh = Null;
	
	int nqw, ib = Null, id;
	
	auto ReadForce = [&](Hydro::Forces &force, String text) {
		bool ismod;
		for (int ih = 0; ih < 2*hd().Nh; ++ih) {
			f.GetLine_discard_empty();
			if ((id = f.GetText().FindAfter(text)) > 0) {
				String str = f.GetText().Mid(id);
				int idh;
				if (str.StartsWith("MOD_")) {
					ismod = true;
					idh = ScanInt(str.Mid(4));
				} else if (str.StartsWith("PH_")) {
					ismod = false;
					idh = ScanInt(str.Mid(3));
				} else
					throw Exc(in.Str() + "\n"  + Format(t_("Wrong text in %s"), ib, f.GetText()));
					
				if (IsNull(idh) || idh < 1 || idh > hd().Nh)
					throw Exc(in.Str() + "\n"  + Format(t_("Wrong force id in %s"), ib, f.GetText()));
		
				idh--;
				
				for (int iw = 0; iw < hd().Nf; ++iw) {
					f.GetLine_discard_empty();
					for (int idof = 0; idof < 6; ++idof) 
						if (ismod)
							force.force[idh](iw, idof + 6*ib) = std::complex<double>(f.GetDouble(1+idof), 0);
						else
							SetPhaseToMag(force.force[idh](iw, idof + 6*ib), f.GetDouble(1+idof));
				}
			}
		}
	};
	
	while(!f.IsEof()) {
		f.GetLine_discard_empty();
		if (f.size() == 0)
			;
		else if (f.GetText(0) == "[PERIODS_NUMBER]") {
			hd().Nf = f.GetInt(1);
			hd().T.SetCount(hd().Nf);
			hd().w.SetCount(hd().Nf);
			hd().qw.resize(hd().Nf);
		} else if (f.GetText(0) == "[QTF_PERIODS_NUMBER]") {
			if (f.GetInt(1) > 0) {
				nqw = f.GetInt(1);
				hd().qw.resize(nqw);
			}
		} else if (f.GetText(0) == "[HEADINGS_NUMBER]") {
			hd().Nh = f.GetInt(1);
			hd().head.SetCount(hd().Nh);
			hd().mdhead.resize(hd().Nh);	
		} else if (f.GetText(0) == "[List_calculated_periods]") {
			for (int iw = 0; iw < hd().Nf; ++iw) {
				f.GetLine_discard_empty();
				hd().T[iw] = f.GetDouble(0);
				hd().w[iw] = hd().qw[iw] = 2*M_PI/hd().T[iw];
			}
		} else if (f.GetText(0) == "[List_calculated_headings]") {
			for (int ih = 0; ih < hd().Nh; ++ih) {
				f.GetLine_discard_empty();
				hd().head[ih] = f.GetDouble(0);
				hd().mdhead[ih] = std::complex<double>(f.GetDouble(0), f.GetDouble(0));
			}
		} else if (f.GetText(0) == "[STRUCTURES_NUMBER]") {
			hd().Nb = f.GetInt(1);
			hd().names.SetCount(hd().Nb);
			hd().Vo.SetCount(hd().Nb);
			hd().cb.resize(3, hd().Nb);
			hd().cg.resize(3, hd().Nb);
			hd().c0.resize(3, hd().Nb);
			hd().C.SetCount(hd().Nb);
			for (int ib = 0; ib < hd().Nb; ++ib) 
				hd().C[ib].setConstant(6, 6, 0);
			hd().M.SetCount(hd().Nb);
			for (int ib = 0; ib < hd().Nb; ++ib) 
				hd().M[ib].setConstant(6, 6, 0);
		    hd().moor.SetCount(hd().Nb);
			for (int ib = 0; ib < hd().Nb; ++ib) 
				hd().moor[ib].setConstant(6, 6, 0);
		    hd().A.SetCount(6*hd().Nb);
			hd().B.SetCount(6*hd().Nb);
			for (int i = 0; i < 6*hd().Nb; ++i) {
				hd().A[i].SetCount(6*hd().Nb);
				hd().B[i].SetCount(6*hd().Nb);
				for (int j = 0; j < 6*hd().Nb; ++j) {
					hd().A[i][j].setConstant(hd().Nf, NaNDouble);	
					hd().B[i][j].setConstant(hd().Nf, NaNDouble);	
				}
			}
			hd().linearDamping = Eigen::MatrixXd::Zero(6*hd().Nb, 6*hd().Nb);
			Hydro::InitMD(hd().md, hd().Nb, hd().Nh, hd().Nf);
		} else if ((id = f.GetText().FindAfter("[STRUCTURE_")) > 0) {
			ib = ScanInt(f.GetText().Mid(id));
			if (ib < 1 || ib > hd().Nb)
				throw Exc(in.Str() + "\n"  + Format(t_("Wrong body %d (%s)"), ib, f.GetText()));
			ib--;
		} else if (f.GetText(0) == "[UNDERWATER_VOLUME]") 
			hd().Vo[ib] = f.GetDouble(1);
		else if (f.GetText(0) == "[CENTER_OF_BUOYANCY]") {
			hd().cb(0, ib) = f.GetDouble(1);
			hd().cb(1, ib) = f.GetDouble(2);
			hd().cb(2, ib) = f.GetDouble(3);
		} else if (f.GetText(0) == "[CENTER_OF_GRAVITY]") {
			hd().cg(0, ib) = f.GetDouble(1);
			hd().cg(1, ib) = f.GetDouble(2);
			hd().cg(2, ib) = f.GetDouble(3);
			hd().c0 = clone(hd().cg);
		} else if (f.GetText(0) == "[Mass_Inertia_matrix]") {
			for (int r = 0; r < 6; ++r) {
				f.GetLine_discard_empty();
				for (int c = 0; c < 6; ++c) 
					hd().M[ib](r, c) = f.GetDouble(c);
			}
		} else if (f.GetText(0) == "[Hydrostatic_matrix]") {
			for (int r = 0; r < 6; ++r) {
				f.GetLine_discard_empty();
				for (int c = 0; c < 6; ++c) 
					hd().C[ib](r, c) = f.GetDouble(c);
			}
		} else if (f.GetText(0) == "[Stiffness_matrix_of_the_mooring_system]") {
			for (int r = 0; r < 6; ++r) {
				f.GetLine_discard_empty();
				for (int c = 0; c < 6; ++c) 
					hd().moor[ib](r, c) = f.GetDouble(c);
			}
		} else if (f.GetText(0) == "[EXCITATION_FORCES_AND_MOMENTS]") {
			hd().Initialize_Forces(hd().ex);
			ReadForce(hd().ex, "[INCIDENCE_EFM_");
		} else if (f.GetText(0) == "[FROUDEKRYLOV_FORCES_AND_MOMENTS]") {
			hd().Initialize_Forces(hd().fk);
			ReadForce(hd().fk, "[INCIDENCE_EFM_FFK_");
		} else if (f.GetText(0) == "[DIFFRACTION_FORCES_AND_MOMENTS]") {
			hd().Initialize_Forces(hd().sc);
			ReadForce(hd().sc, "[INCIDENCE_EFM_DIFF_");
		} else if (f.GetText(0) == "[RAO]" && !IsNull(hd().Nf)) {		// There are two [RAO] fields!
			hd().Initialize_Forces(hd().rao);
			ReadForce(hd().rao, "[INCIDENCE_RAO_");
		// What is 'INTER_RAO'?
		} else if (f.GetText(0) == "[DRIFT_FORCES_AND_MOMENTS]") {
			for (int ih = 0; ih < hd().Nh; ++ih) {
				f.GetLine_discard_empty();
				if ((id = f.GetText().FindAfter("[INCIDENCE_DFM_")) > 0) {
					String str = f.GetText().Mid(id);
					int idh = ScanInt(str);
	
					if (IsNull(idh) || idh < 1 || idh > hd().Nh)
						throw Exc(in.Str() + "\n"  + Format(t_("Wrong force id in %s"), ib, f.GetText()));
			
					idh--;
				
					for (int iw = 0; iw < hd().Nf; ++iw) {
						f.GetLine_discard_empty();
						for (int idof = 0; idof < 6; ++idof) 
							hd().md[ib][idh][idof](iw) = f.GetDouble(1+idof);
					}
				}
			}
		} else if (f.GetText(0) == "[Added_mass_Radiation_Damping]") {
			for (int i = 0; i < 6*2; ++i) {
				f.GetLine_discard_empty();
				if ((id = f.GetText().FindAfter("[ADDED_MASS_LINE_")) > 0) {
					String str = f.GetText().Mid(id);
					int idrow = ScanInt(str);
			
					if (IsNull(idrow) || idrow < 1 || idrow > 6)
						throw Exc(in.Str() + "\n"  + Format(t_("Wrong row id in %s"), ib, f.GetText()));
			
					idrow--;
			
					for (int iw = 0; iw < hd().Nf; ++iw) {
						f.GetLine_discard_empty();
						for (int idof = 0; idof < 6; ++idof) 
							hd().A[idrow + 6*ib][idof + 6*ib][iw] = f.GetDouble(1+idof);
					}
				} else if ((id = f.GetText().FindAfter("[DAMPING_TERM_")) > 0) {
					String str = f.GetText().Mid(id);
					int idrow = ScanInt(str);
			
					if (IsNull(idrow) || idrow < 1 || idrow > 6)
						throw Exc(in.Str() + "\n"  + Format(t_("Wrong row id in %s"), ib, f.GetText()));
			
					idrow--;
			
					for (int iw = 0; iw < hd().Nf; ++iw) {
						f.GetLine_discard_empty();
						for (int idof = 0; idof < 6; ++idof) 
							hd().B[idrow + 6*ib][idof + 6*ib][iw] = f.GetDouble(1+idof);
					}
				}
			}
		}
	}
	hd().mdtype = 9;
	
	return true;
}

