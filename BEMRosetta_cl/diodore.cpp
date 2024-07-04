// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"


String Diodore::Load(String file, double) {
	dt.file = file;
	dt.name = GetFileTitle(file);
	dt.dimen = true;
	dt.len = 1;
	dt.solver = Hydro::DIODORE;
	dt.Nb = Null;
	
	try {
		BEM::Print("\n\n" + Format(t_("Loading '%s'"), file));

		BEM::Print("\n- " + S(t_("HDB file")));
		
		Load_HDB();
		
		if (IsNull(dt.Nb))
			return t_("No data found");
	
		/*dt.dof.Clear();	dt.dof.SetCount(dt.Nb, 0);
		for (int i = 0; i < dt.Nb; ++i)
			dt.dof[i] = 6;*/
	} catch (Exc e) {
		//BEM::PrintError(Format("\n%s: %s", t_("Error"), e));
		//dt.lastError = e;
		return e;
	}
	
	return String();
}

void Diodore::Load_HDB() {
	String fileName = ForceExtSafer(dt.file, ".hdb");
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		throw Exc(in.Str() + "\n" + t_("File not found or blocked"));
	
	String line; 
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
	
	//dt.dataFromW = false;
	dt.dimen = true;
	
	dt.Nb = dt.Nf = dt.Nh = Null;
	
	int nqw, ib = Null, id;
	
	auto ReadForce = [&](Hydro::Forces &force, String text) {
		bool ismod;
		for (int ih = 0; ih < 2*dt.Nh; ++ih) {
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
					throw Exc(in.Str() + "\n"  + Format(t_("Wrong text in %s"), f.GetText()));
					
				if (IsNull(idh) || idh < 1 || idh > dt.Nh)
					throw Exc(in.Str() + "\n"  + Format(t_("Wrong force id in %s"), f.GetText()));
		
				idh--;
				
				for (int iw = 0; iw < dt.Nf; ++iw) {
					f.GetLine_discard_empty();
					for (int idf = 0; idf < 6; ++idf) 
						if (ismod)
							force[ib][idh](iw, idf) = std::complex<double>(f.GetDouble(1+idf), 0);
						else
							SetPhaseToMag(force[ib][idh](iw, idf), f.GetDouble(1+idf));
				}
			}
		}
	};
	
	while(!f.IsEof()) {
		f.GetLine_discard_empty();
		if (f.size() == 0)
			;
		else if (f.GetText(0) == "[PERIODS_NUMBER]") {
			dt.Nf = f.GetInt(1);
			//dt.T.SetCount(dt.Nf);
			dt.w.SetCount(dt.Nf);
		} else if (f.GetText(0) == "[QTF_PERIODS_NUMBER]") {
			if (f.GetInt(1) > 0) {
				nqw = f.GetInt(1);
				dt.qw.resize(nqw);
			}
		} else if (f.GetText(0) == "[HEADINGS_NUMBER]") {
			dt.Nh = f.GetInt(1);
			dt.head.SetCount(dt.Nh);
			dt.mdhead.resize(dt.Nh);	
		} else if (f.GetText(0) == "[List_calculated_periods]") {
			for (int iw = 0; iw < dt.Nf; ++iw) {
				f.GetLine_discard_empty();
				double T = f.GetDouble(0);
				dt.w[iw] = 2*M_PI/T;
			}
		} else if (f.GetText(0) == "[List_calculated_headings]") {
			for (int ih = 0; ih < dt.Nh; ++ih) {
				f.GetLine_discard_empty();
				dt.head[ih] = f.GetDouble(0);
				dt.mdhead[ih] = std::complex<double>(f.GetDouble(0), f.GetDouble(0));
			}
		} else if (f.GetText(0) == "[STRUCTURES_NUMBER]") {
			dt.Nb = f.GetInt(1);
			dt.msh.SetCount(dt.Nb);
			//dt.names.SetCount(dt.Nb);
			//dt.Vo.SetCount(dt.Nb);
			//dt.cb.resize(3, dt.Nb);
			//dt.cg.resize(3, dt.Nb);
			//dt.c0.resize(3, dt.Nb);
			//dt.C.SetCount(dt.Nb);
			for (int iib = 0; iib < dt.Nb; ++iib) 
				dt.msh[iib].dt.C.setConstant(6, 6, 0);
			//dt.M.SetCount(dt.Nb);
			for (int iib = 0; iib < dt.Nb; ++iib) 
				dt.msh[iib].dt.M.setConstant(6, 6, 0);
		    //dt.Cmoor.SetCount(dt.Nb);
			for (int iib = 0; iib < dt.Nb; ++iib) 
				dt.msh[iib].dt.Cmoor.setConstant(6, 6, 0);
		    Initialize_AB(dt.A);
			Initialize_AB(dt.B);
			//dt.Dlin = Eigen::MatrixXd::Zero(6*dt.Nb, 6*dt.Nb);
			Hydro::Initialize_MD(dt.md, dt.Nb, dt.Nh, dt.Nf);
		} else if ((id = f.GetText().FindAfter("[STRUCTURE_")) > 0) {
			ib = ScanInt(f.GetText().Mid(id));
			if (ib < 1 || ib > dt.Nb)
				throw Exc(in.Str() + "\n"  + Format(t_("Wrong body %d (%s)"), ib, f.GetText()));
			ib--;
			dt.msh[ib].dt.name = f.GetText(1);
		} else if (f.GetText(0) == "[UNDERWATER_VOLUME]") 
			dt.msh[ib].dt.Vo = f.GetDouble(1);
		else if (f.GetText(0) == "[CENTER_OF_BUOYANCY]") {
			dt.msh[ib].dt.cb.x = f.GetDouble(1);
			dt.msh[ib].dt.cb.y = f.GetDouble(2);
			dt.msh[ib].dt.cb.z = f.GetDouble(3);
		} else if (f.GetText(0) == "[CENTER_OF_GRAVITY]") {
			dt.msh[ib].dt.cg.x = f.GetDouble(1);
			dt.msh[ib].dt.cg.y = f.GetDouble(2);
			dt.msh[ib].dt.cg.z = f.GetDouble(3);
			dt.msh[ib].dt.c0 = clone(dt.msh[ib].dt.cg);
		} else if (f.GetText(0) == "[Mass_Inertia_matrix]") {
			for (int r = 0; r < 6; ++r) {
				f.GetLine_discard_empty();
				for (int c = 0; c < 6; ++c) 
					dt.msh[ib].dt.M(r, c) = f.GetDouble(c);
			}
		} else if (f.GetText(0) == "[Hydrostatic_matrix]") {
			for (int r = 0; r < 6; ++r) {
				f.GetLine_discard_empty();
				for (int c = 0; c < 6; ++c) 
					dt.msh[ib].dt.C(r, c) = f.GetDouble(c);
			}
		} else if (f.GetText(0) == "[Stiffness_matrix_of_the_mooring_system]") {
			for (int r = 0; r < 6; ++r) {
				f.GetLine_discard_empty();
				for (int c = 0; c < 6; ++c) 
					dt.msh[ib].dt.Cmoor(r, c) = f.GetDouble(c);
			}
		} else if (f.GetText(0) == "[EXCITATION_FORCES_AND_MOMENTS]") {
			Initialize_Forces(dt.ex);
			ReadForce(dt.ex, "[INCIDENCE_EFM_");
		} else if (f.GetText(0) == "[FROUDEKRYLOV_FORCES_AND_MOMENTS]") {
			Initialize_Forces(dt.fk);
			ReadForce(dt.fk, "[INCIDENCE_EFM_FFK_");
		} else if (f.GetText(0) == "[DIFFRACTION_FORCES_AND_MOMENTS]") {
			Initialize_Forces(dt.sc);
			ReadForce(dt.sc, "[INCIDENCE_EFM_DIFF_");
		} else if (f.GetText(0) == "[RAO]" && !IsNull(dt.Nf)) {		// There are two [RAO] fields!
			Initialize_Forces(dt.rao);
			ReadForce(dt.rao, "[INCIDENCE_RAO_");
		// What is 'INTER_RAO'?
		} else if (f.GetText(0) == "[DRIFT_FORCES_AND_MOMENTS]") {
			for (int ih = 0; ih < dt.Nh; ++ih) {
				f.GetLine_discard_empty();
				if ((id = f.GetText().FindAfter("[INCIDENCE_DFM_")) > 0) {
					String str = f.GetText().Mid(id);
					int idh = ScanInt(str);
	
					if (IsNull(idh) || idh < 1 || idh > dt.Nh)
						throw Exc(in.Str() + "\n"  + Format(t_("Wrong force id in %s"), f.GetText()));
			
					idh--;
				
					for (int iw = 0; iw < dt.Nf; ++iw) {
						f.GetLine_discard_empty();
						for (int idof = 0; idof < 6; ++idof) 
							dt.md[ib][idh][idof](iw) = f.GetDouble(1+idof);
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
						throw Exc(in.Str() + "\n"  + Format(t_("Wrong row id in %s"), f.GetText()));
			
					idrow--;
			
					for (int iw = 0; iw < dt.Nf; ++iw) {
						f.GetLine_discard_empty();
						for (int idof = 0; idof < 6; ++idof) 
							dt.A[idrow + 6*ib][idof + 6*ib][iw] = f.GetDouble(1+idof);
					}
				} else if ((id = f.GetText().FindAfter("[DAMPING_TERM_")) > 0) {
					String str = f.GetText().Mid(id);
					int idrow = ScanInt(str);
			
					if (IsNull(idrow) || idrow < 1 || idrow > 6)
						throw Exc(in.Str() + "\n"  + Format(t_("Wrong row id in %s"), f.GetText()));
			
					idrow--;
			
					for (int iw = 0; iw < dt.Nf; ++iw) {
						f.GetLine_discard_empty();
						for (int idof = 0; idof < 6; ++idof) 
							dt.B[idrow + 6*ib][idof + 6*ib][iw] = f.GetDouble(1+idof);
					}
				}
			}
		}
	}
	dt.qtftype = dt.mdtype = 9;
}

void Hydro::SaveDiodoreHDB(String file) const {
	BEM::Print("\n\n" + Format(t_("Saving '%s'"), file));
	
	String folder = GetFileFolder(file);
	String name = GetFileTitle(file);
	Time t = GetSysTime();
	
	FileOut out(file);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to save '%s'. File already used."), file));
	
	out << "[SOFT]          DIODORE" << "\n" 	
		<< "[VERSION]       5.5.2" << "\n"
		<< "[Date]          " << Format("%02d:%02d:%02d    %Mon %02d %d\n", t.hour, t.minute, t.second, t.month, t.day, t.year)
		<< "[INPUT_FILE]    " << name << "\n"
		<< "[Locally_At]    " << folder << "\n"
		<< "[UNIT]" << "\n"          
		<< "[FORWARD_SPEED]        0.00" << "\n"
		<< "[HYDRO_PARA]        Y" << "\n"
		<< "[DETAIL_EXC]        Y" << "\n"
		<< "[RAO]               " << (IsLoadedRAO() ? 'Y' : 'N') << "\n"
		<< "[DRIFT_FORCE]       " << (IsLoadedMD() ? 'Y' : 'N') << "\n"			
		<< "[QTF]               " << (IsLoadedQTF(false) ? 'Y' : 'N') << "\n"
		<< "[PERIODS_NUMBER]    " << dt.Nf << "\n"
		<< "[INTER_PERIODS_NB]  " << 0 << "\n"
		<< "[QTF_PERIODS_NUMBER]    " << dt.qw.size() << "\n"							
		<< "[QTF_SWEEPING_NUMBER]   " << 0 << "\n"	
		<< "[QTF_SWEEPING_STEP]          0.000\n"		
		<< "[HEADINGS_NUMBER]   " << dt.Nh << "\n"
		<< Format("[LOWEST_HEADING]    %.2f\n", (dt.Nh > 0 ? First(dt.head) : 0))
		<< Format("[HIGHEST_HEADING]   %.2f\n", (dt.Nh > 0 ? Last(dt.head) : 0));	

	out << "[List_calculated_periods]   \n";
	
	VectorXd TT	= Get_T();							
	for (auto T : TT)
		out << Format("%10.3f\n", T);							
								
	out << "[List_calculated_headings]  \n";								
	for (auto h : dt.head)
		out << Format("%10.3f\n", h);	
										
	out << "[STRUCTURES_NUMBER]        " << dt.Nb << "\n";							


	auto WriteForce = [&](const Hydro::Forces &force, int ib, String text) {
		for (int ih = 0; ih < dt.Nh; ++ih) {
			out << Format("[%s_MOD_%03d]   %10.3f\n", text, ih+1, dt.head[ih]);
			for (int iw = 0; iw < dt.Nf; ++iw) {
				out << Format("%7.2f", TT[iw]);
				for (int idof = 0; idof < 6; ++idof) 
					out << Format(" % .7E", abs(F_dim(force, ih, iw, idof, ib)));	
				out << "\n";
			}
		}
		for (int ih = 0; ih < dt.Nh; ++ih) {
			out << Format("[%s_PH_%03d]   %10.3f\n", text, ih+1, dt.head[ih]);
			for (int iw = 0; iw < dt.Nf; ++iw) {
				out << Format("%7.2f", TT[iw]);
				for (int idof = 0; idof < 6; ++idof) 
					out << Format(" % .7E", arg(F_dim(force, ih, iw, idof, ib)));	
				out << "\n";
			}
		}
	};
		
	for (int ib = 0; ib < dt.Nb; ++ib) {
		out << Format("[STRUCTURE_%02d] %s\n", ib+1, dt.msh[ib].dt.name);
		out << Format("[UNDERWATER_VOLUME]        %.3f\n", dt.msh[ib].dt.Vo);
		out << Format("[CENTER_OF_BUOYANCY]     %.4f   %.4f   %.4f\n", dt.msh[ib].dt.cb.x, dt.msh[ib].dt.cb.x, dt.msh[ib].dt.cb.x);
		out << Format("[CENTER_OF_GRAVITY]      %.4f   %.4f   %.4f\n", dt.msh[ib].dt.cg.x, dt.msh[ib].dt.cg.y, dt.msh[ib].dt.cg.z);
		
		out << "[Mass_Inertia_matrix]   \n";		
		for (int r = 0; r < 6; ++r) {
			for (int c = 0; c < 6; ++c) {
				double d = IsLoadedM() ? dt.msh[ib].dt.M(r, c) : 0;
				out << Format("% .2E   ", d);
			}
			out << "\n";
		}
		out << "[Hydrostatic_matrix]    \n";		
		for (int r = 0; r < 6; ++r) {
			for (int c = 0; c < 6; ++c) {
				double d = IsLoadedC() ? C_dim(ib, r, c) : 0;
				out << Format("% .7E ", d);
			}
			out << "\n";
		}
		out << "[Stiffness_matrix_of_the_mooring_system]\n";		
		for (int r = 0; r < 6; ++r) {
			for (int c = 0; c < 6; ++c) {
				double d = IsLoadedCMoor() ? CMoor_dim(ib, r, c) : 0;
				out << Format("% .7E ", d);
			}
			out << "\n";
		}
		out << "\n";
		
		if (IsLoadedFex()) {
			out << "[EXCITATION_FORCES_AND_MOMENTS]  \n";
			WriteForce(dt.ex, ib, "INCIDENCE_EFM");
		}
		if (IsLoadedFfk()) {
			out << "[FROUDEKRYLOV_FORCES_AND_MOMENTS]\n";
			WriteForce(dt.fk, ib, "INCIDENCE_EFM_FFK");
		}
		if (IsLoadedFsc()) {
			out << "[DIFFRACTION_FORCES_AND_MOMENTS]\n";
			WriteForce(dt.sc, ib, "INCIDENCE_EFM_DIFF");
		}
		if (IsLoadedRAO()) {
			out << "[RAO]\n";
			WriteForce(dt.rao, ib, "INCIDENCE_RAO");
		}
		if (IsLoadedMD()) {
			out << "[DRIFT_FORCES_AND_MOMENTS]\n";
			for (int ih = 0; ih < dt.mdhead.size(); ++ih) {
				out << Format("[INCIDENCE_DFM_%03d]   %10.3f\n", ih+1, real(dt.mdhead[ih]));
				for (int iw = 0; iw < dt.Nf; ++iw) {
					out << Format("%7.2f", TT[iw]);
					for (int idf = 0; idf < 6; ++idf) 
						out << Format(" % .7E", Md_dim(idf, ih, iw));
					out << "\n";
				}
			}
		}
		if (IsLoadedA() || IsLoadedB()) {
			out << "[Added_mass_Radiation_Damping]\n";
			if (IsLoadedA()) {
				for (int idrow = 0; idrow < 6; ++idrow) {
					out << Format("[ADDED_MASS_LINE_%d]\n", idrow+1);
					for (int iw = 0; iw < dt.Nf; ++iw) {
						out << Format("%7.2f", TT[iw]);
						for (int idof = 0; idof < 6; ++idof) 
							out << Format(" % .7E", A_dim(iw, idrow + 6*ib, idof + 6*ib));
						out << "\n";
					}
				}
			}
			if (IsLoadedB()) {
				for (int idrow = 0; idrow < 6; ++idrow) {
					out << Format("[DAMPING_TERM_%d]\n", idrow+1);
					for (int iw = 0; iw < dt.Nf; ++iw) {
						out << Format("%7.2f", TT[iw]);
						for (int idof = 0; idof < 6; ++idof) 
							out << Format(" % .7E", B_dim(iw, idrow + 6*ib, idof + 6*ib));
						out << "\n";
					}
				}
			}
		}					
	}
}