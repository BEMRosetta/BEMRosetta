// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"


String Diodore::Load(String file, double) {
	hd().file = file;
	hd().name = GetFileTitle(file);
	hd().dimen = true;
	hd().len = 1;
	hd().solver = Hydro::DIODORE;
	hd().Nb = Null;
	
	try {
		BEM::Print("\n\n" + Format(t_("Loading '%s'"), file));

		BEM::Print("\n- " + S(t_("HDB file")));
		
		Load_HDB();
		
		if (IsNull(hd().Nb))
			return t_("No data found");
	
		/*hd().dof.Clear();	hd().dof.SetCount(hd().Nb, 0);
		for (int i = 0; i < hd().Nb; ++i)
			hd().dof[i] = 6;*/
	} catch (Exc e) {
		//BEM::PrintError(Format("\n%s: %s", t_("Error"), e));
		//hd().lastError = e;
		return e;
	}
	
	return String();
}

void Diodore::Load_HDB() {
	String fileName = ForceExtSafer(hd().file, ".hdb");
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		throw Exc(in.Str() + "\n" + t_("File not found or blocked"));
	
	String line; 
	LineParser f(in);
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
					throw Exc(in.Str() + "\n"  + Format(t_("Wrong text in %s"), f.GetText()));
					
				if (IsNull(idh) || idh < 1 || idh > hd().Nh)
					throw Exc(in.Str() + "\n"  + Format(t_("Wrong force id in %s"), f.GetText()));
		
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
			//hd().T.SetCount(hd().Nf);
			hd().w.SetCount(hd().Nf);
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
				double T = f.GetDouble(0);
				hd().w[iw] = 2*M_PI/T;
			}
		} else if (f.GetText(0) == "[List_calculated_headings]") {
			for (int ih = 0; ih < hd().Nh; ++ih) {
				f.GetLine_discard_empty();
				hd().head[ih] = f.GetDouble(0);
				hd().mdhead[ih] = std::complex<double>(f.GetDouble(0), f.GetDouble(0));
			}
		} else if (f.GetText(0) == "[STRUCTURES_NUMBER]") {
			hd().Nb = f.GetInt(1);
			hd().msh.SetCount(hd().Nb);
			//hd().names.SetCount(hd().Nb);
			//hd().Vo.SetCount(hd().Nb);
			//hd().cb.resize(3, hd().Nb);
			//hd().cg.resize(3, hd().Nb);
			//hd().c0.resize(3, hd().Nb);
			//hd().C.SetCount(hd().Nb);
			for (int ib = 0; ib < hd().Nb; ++ib) 
				hd().msh[ib].C.setConstant(6, 6, 0);
			//hd().M.SetCount(hd().Nb);
			for (int ib = 0; ib < hd().Nb; ++ib) 
				hd().msh[ib].M.setConstant(6, 6, 0);
		    //hd().Cmoor.SetCount(hd().Nb);
			for (int ib = 0; ib < hd().Nb; ++ib) 
				hd().msh[ib].Cmoor.setConstant(6, 6, 0);
		    hd().Initialize_AB(hd().A);
			hd().Initialize_AB(hd().B);
			//hd().Dlin = Eigen::MatrixXd::Zero(6*hd().Nb, 6*hd().Nb);
			Hydro::Initialize_MD(hd().md, hd().Nb, hd().Nh, hd().Nf);
		} else if ((id = f.GetText().FindAfter("[STRUCTURE_")) > 0) {
			ib = ScanInt(f.GetText().Mid(id));
			if (ib < 1 || ib > hd().Nb)
				throw Exc(in.Str() + "\n"  + Format(t_("Wrong body %d (%s)"), ib, f.GetText()));
			ib--;
			hd().msh[ib].name = f.GetText(1);
		} else if (f.GetText(0) == "[UNDERWATER_VOLUME]") 
			hd().msh[ib].Vo = f.GetDouble(1);
		else if (f.GetText(0) == "[CENTER_OF_BUOYANCY]") {
			hd().msh[ib].cb.x = f.GetDouble(1);
			hd().msh[ib].cb.y = f.GetDouble(2);
			hd().msh[ib].cb.z = f.GetDouble(3);
		} else if (f.GetText(0) == "[CENTER_OF_GRAVITY]") {
			hd().msh[ib].cg.x = f.GetDouble(1);
			hd().msh[ib].cg.y = f.GetDouble(2);
			hd().msh[ib].cg.z = f.GetDouble(3);
			hd().msh[ib].c0 = clone(hd().msh[ib].cg);
		} else if (f.GetText(0) == "[Mass_Inertia_matrix]") {
			for (int r = 0; r < 6; ++r) {
				f.GetLine_discard_empty();
				for (int c = 0; c < 6; ++c) 
					hd().msh[ib].M(r, c) = f.GetDouble(c);
			}
		} else if (f.GetText(0) == "[Hydrostatic_matrix]") {
			for (int r = 0; r < 6; ++r) {
				f.GetLine_discard_empty();
				for (int c = 0; c < 6; ++c) 
					hd().msh[ib].C(r, c) = f.GetDouble(c);
			}
		} else if (f.GetText(0) == "[Stiffness_matrix_of_the_mooring_system]") {
			for (int r = 0; r < 6; ++r) {
				f.GetLine_discard_empty();
				for (int c = 0; c < 6; ++c) 
					hd().msh[ib].Cmoor(r, c) = f.GetDouble(c);
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
						throw Exc(in.Str() + "\n"  + Format(t_("Wrong force id in %s"), f.GetText()));
			
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
						throw Exc(in.Str() + "\n"  + Format(t_("Wrong row id in %s"), f.GetText()));
			
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
						throw Exc(in.Str() + "\n"  + Format(t_("Wrong row id in %s"), f.GetText()));
			
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
	hd().qtftype = hd().mdtype = 9;
}

void Hydro::SaveDiodoreHDB(String file) {
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
		<< "[RAO]               " << (hd().IsLoadedRAO() ? 'Y' : 'N') << "\n"
		<< "[DRIFT_FORCE]       " << (hd().IsLoadedMD() ? 'Y' : 'N') << "\n"			
		<< "[QTF]               " << (hd().IsLoadedQTF(false) ? 'Y' : 'N') << "\n"
		<< "[PERIODS_NUMBER]    " << hd().Nf << "\n"
		<< "[INTER_PERIODS_NB]  " << 0 << "\n"
		<< "[QTF_PERIODS_NUMBER]    " << hd().qw.size() << "\n"							
		<< "[QTF_SWEEPING_NUMBER]   " << 0 << "\n"	
		<< "[QTF_SWEEPING_STEP]          0.000\n"		
		<< "[HEADINGS_NUMBER]   " << hd().Nh << "\n"
		<< Format("[LOWEST_HEADING]    %.2f\n", (hd().Nh > 0 ? First(hd().head) : 0))
		<< Format("[HIGHEST_HEADING]   %.2f\n", (hd().Nh > 0 ? Last(hd().head) : 0));	

	out << "[List_calculated_periods]   \n";
	
	VectorXd TT	= hd().Get_T();							
	for (auto T : TT)
		out << Format("%10.3f\n", T);							
								
	out << "[List_calculated_headings]  \n";								
	for (auto h : hd().head)
		out << Format("%10.3f\n", h);	
										
	out << "[STRUCTURES_NUMBER]        " << hd().Nb << "\n";							


	auto WriteForce = [&](Hydro::Forces &force, String text) {
		for (int ih = 0; ih < hd().Nh; ++ih) {
			out << Format("[%s_MOD_%03d]   %10.3f\n", text, ih+1, hd().head[ih]);
			for (int iw = 0; iw < hd().Nf; ++iw) {
				out << Format("%7.2f", TT[iw]);
				for (int idof = 0; idof < 6; ++idof) 
					out << Format(" % .7E", abs(hd().F_dim(force, ih, iw, idof)));	
				out << "\n";
			}
		}
		for (int ih = 0; ih < hd().Nh; ++ih) {
			out << Format("[%s_PH_%03d]   %10.3f\n", text, ih+1, hd().head[ih]);
			for (int iw = 0; iw < hd().Nf; ++iw) {
				out << Format("%7.2f", TT[iw]);
				for (int idof = 0; idof < 6; ++idof) 
					out << Format(" % .7E", arg(hd().F_dim(force, ih, iw, idof)));	
				out << "\n";
			}
		}
	};
		
	for (int ib = 0; ib < hd().Nb; ++ib) {
		out << Format("[STRUCTURE_%02d] %s\n", ib+1, hd().msh[ib].name);
		out << Format("[UNDERWATER_VOLUME]        %.3f\n", hd().msh[ib].Vo);
		out << Format("[CENTER_OF_BUOYANCY]     %.4f   %.4f   %.4f\n", hd().msh[ib].cb.x, hd().msh[ib].cb.x, hd().msh[ib].cb.x);
		out << Format("[CENTER_OF_GRAVITY]      %.4f   %.4f   %.4f\n", hd().msh[ib].cg.x, hd().msh[ib].cg.y, hd().msh[ib].cg.z);
		
		out << "[Mass_Inertia_matrix]   \n";		
		for (int r = 0; r < 6; ++r) {
			for (int c = 0; c < 6; ++c) {
				double d = hd().IsLoadedM() ? hd().msh[ib].M(r, c) : 0;
				out << Format("% .2E   ", d);
			}
			out << "\n";
		}
		out << "[Hydrostatic_matrix]    \n";		
		for (int r = 0; r < 6; ++r) {
			for (int c = 0; c < 6; ++c) {
				double d = hd().IsLoadedC() ? hd().C_dim(ib, r, c) : 0;
				out << Format("% .7E ", d);
			}
			out << "\n";
		}
		out << "[Stiffness_matrix_of_the_mooring_system]\n";		
		for (int r = 0; r < 6; ++r) {
			for (int c = 0; c < 6; ++c) {
				double d = hd().IsLoadedCMoor() ? hd().CMoor_dim(ib, r, c) : 0;
				out << Format("% .7E ", d);
			}
			out << "\n";
		}
		out << "\n";
		
		if (hd().IsLoadedFex()) {
			out << "[EXCITATION_FORCES_AND_MOMENTS]  \n";
			WriteForce(hd().ex, "INCIDENCE_EFM");
		}
		if (hd().IsLoadedFfk()) {
			out << "[FROUDEKRYLOV_FORCES_AND_MOMENTS]\n";
			WriteForce(hd().fk, "INCIDENCE_EFM_FFK");
		}
		if (hd().IsLoadedFsc()) {
			out << "[DIFFRACTION_FORCES_AND_MOMENTS]\n";
			WriteForce(hd().sc, "INCIDENCE_EFM_DIFF");
		}
		if (hd().IsLoadedRAO()) {
			out << "[RAO]\n";
			WriteForce(hd().rao, "INCIDENCE_RAO");
		}
		if (hd().IsLoadedMD()) {
			out << "[DRIFT_FORCES_AND_MOMENTS]\n";
			for (int ih = 0; ih < hd().mdhead.size(); ++ih) {
				out << Format("[INCIDENCE_DFM_%03d]   %10.3f\n", ih+1, real(hd().mdhead[ih]));
				for (int iw = 0; iw < hd().Nf; ++iw) {
					out << Format("%7.2f", TT[iw]);
					for (int idf = 0; idf < 6; ++idf) 
						out << Format(" % .7E", hd().Md_dim(idf, ih, iw));
					out << "\n";
				}
			}
		}
		if (hd().IsLoadedA() || hd().IsLoadedB()) {
			out << "[Added_mass_Radiation_Damping]\n";
			if (hd().IsLoadedA()) {
				for (int idrow = 0; idrow < 6; ++idrow) {
					out << Format("[ADDED_MASS_LINE_%d]\n", idrow+1);
					for (int iw = 0; iw < hd().Nf; ++iw) {
						out << Format("%7.2f", TT[iw]);
						for (int idof = 0; idof < 6; ++idof) 
							out << Format(" % .7E", hd().A_dim(iw, idrow + 6*ib, idof + 6*ib));
						out << "\n";
					}
				}
			}
			if (hd().IsLoadedB()) {
				for (int idrow = 0; idrow < 6; ++idrow) {
					out << Format("[DAMPING_TERM_%d]\n", idrow+1);
					for (int iw = 0; iw < hd().Nf; ++iw) {
						out << Format("%7.2f", TT[iw]);
						for (int idof = 0; idof < 6; ++idof) 
							out << Format(" % .7E", hd().B_dim(iw, idrow + 6*ib, idof + 6*ib));
						out << "\n";
					}
				}
			}
		}					
	}
}