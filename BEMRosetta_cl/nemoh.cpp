// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#include <STEM4U/Utility.h>

bool Nemoh::Load(String file, Function <bool(String, int)> Status, double) {
	try {
		String ext = GetFileExt(file); 

		if (ext == ".tec") {
			String folder = GetFileFolder(file);
			String folderTitle = GetFileName(folder);
			if (ToLower(folderTitle) != "results") 
				throw Exc(Format(t_(".tec file '%s' should have to be in 'results' folder"), file));
			bool found = false;
			String upperFolder = GetUpperFolder(folder);
			for (FindFile ff(AppendFileNameX(upperFolder, "*.*")); ff; ++ff) {
				if (ff.IsFile()) {
					if (ToLower(ff.GetName()) == "nemoh.cal") {
						file = ff.GetPath();
						found = true;
						break;
					}
				}
			}
			if (!found)
				throw Exc(Format(t_("nemoh.cal file not found in '%s' folder"), upperFolder));
		}
	
		if (ext == ".cal" || ext == ".tec")
			hd().code = Hydro::NEMOH;
		else
			hd().code = Hydro::SEAFEM_NEMOH;
	
		hd().file = file;
		hd().name = GetFileTitle(GetFileFolder(file));
		folder = GetFileFolder(file);
		hd().len = 1;
		hd().dimen = true;
		hd().Nb = Null;
	
		String fileCal;
		BEM::Print("\n\n" + Format(t_("Loading '%s'"), file));
		if (hd().code == Hydro::NEMOH) 
			fileCal = file;
		else 
			fileCal = AppendFileNameX(folder, "Nemoh_output/Nemoh.cal");
		if (!Load_Cal(fileCal)) 
			throw Exc(Format(t_("File '%s' not found"), fileCal));
		
		String fileRad, folderForces;
		if (hd().code == Hydro::NEMOH) {
			BEM::Print(S("\n- ") + t_("Hydrostatics file(s) 'Mesh/Hydrostatics*.dat'"));
			if (!Load_Hydrostatics("Mesh"))
				BEM::PrintWarning(S(": ** Mesh/Hydrostatics*.dat ") + t_("Not found") + "**");
			BEM::Print(S("\n- ") + t_("KH file(s) 'Mesh/KH*.dat'"));
			if (!Load_KH("Mesh"))
				BEM::PrintWarning(S(": ** Mesh/KH ") + t_("Not found") + "**");
			
			dynamic_cast<Wamit *>(this)->Load_frc2(ForceExt(fileCal, ".frc"));
			
			fileRad = AppendFileNameX(folder, AppendFileNameX("Results", "RadiationCoefficients.tec"));
			folderForces = folder;
		} else {
			if (!Load_Inf(file)) 
				throw Exc(Format(t_("File '%s' not found"), file));

			fileRad = AppendFileNameX(folder, "Nemoh_output/Results", "RadiationCoefficients.tec");
			folderForces = AppendFileNameX(folder, "Nemoh_output");
		} 
		BEM::Print(S("\n- ") + t_("Radiation file 'RadiationCoefficients.tec'"));
		if (!Load_Radiation(fileRad))
			BEM::PrintWarning(S(": ** RadiationCoefficients.tec ") + t_("Not found") + "**");

		BEM::Print(S("\n- ") + t_("Excitation force file 'ExcitationForce.tec'"));
		if (!Load_Excitation(folderForces))
			BEM::PrintWarning(S(": ** ExcitationForce.tec ") + t_("Not found") + "**");
		
		BEM::Print(S("\n- ") + t_("Diffraction force file 'DiffractionForce.tec'"));
		if (!Load_Diffraction(folderForces))
			BEM::PrintWarning(S(": ** DiffractionForce.tec ") + t_("Not found") + "**");
		BEM::Print(S("\n- ") + t_("Froude Krylov file 'FKForce.tec'"));
		if (!Load_FroudeKrylov(folderForces))
			BEM::PrintWarning(S(": ** FKForce.tec ") + t_("Not found") + "**");
		
		if (dcase.solver == BEMCase::NEMOHv3) {
			Load_Inertia("mechanics");
			Load_LinearDamping("mechanics");
			Load_QTF(AppendFileNameX("results", "qtf"), Status);
		}
		
		UVector<int> idsRemove;
		for (int ih = 0; ih < hd().Nh; ++ih) {
			int id = FindDelta(hd().head, hd().head[ih], 0.001, ih+1);
			if (id > 0)
				idsRemove << ih;
		}
		hd().DeleteHeadings(idsRemove);
		
		if (hd().code == Hydro::NEMOH) {
			if (!hd().dof.IsEmpty()) {
				BEM::Print(S("\n- ") + t_("IRF file(s) 'IRF.tec'"));
				if (!Load_IRF(AppendFileNameX(folder, "Results", "IRF.tec")))
					BEM::PrintWarning(S(": ** IRF.tec ") + t_("Not found") + "**");
			}
		}
		if (IsNull(hd().Nb))
			return false;
	} catch (Exc e) {
		BEM::PrintError(Format("\n%s: %s", t_("Error"), e));
		//hd().lastError = e;
		return false;
	}
	return true;
}

bool Nemoh::Load_Cal(String fileName) {	
	if (!static_cast<NemohCase&>(dcase).Load(fileName))
		return false;

	hd().rho = dcase.rho;
	hd().g = dcase.g;
	hd().h = dcase.h;
	if (hd().h == 0)
		hd().h = -1;
	hd().Nb = dcase.bodies.size();
	for (int i = 0; i < hd().Nb; ++i) 	
		hd().names << GetFileTitle(dcase.bodies[i].meshFile);
	hd().dof.SetCount(hd().Nb, 6);
	//for (int i = 0; i < hd().Nb; ++i)
	//	hd().dof[i] = data.bodies[i].ndof;
	hd().Nf = dcase.Nf;
	LinSpaced(hd().w, hd().Nf, dcase.minF, dcase.maxF); 
 	hd().T.SetCount(hd().Nf);
    for (int i = 0; i < hd().Nf; ++i) 
		hd().T[i] = 2*M_PI/hd().w[i];  
   
	hd().Nh = dcase.Nh;  						
    LinSpaced(hd().head, hd().Nh, dcase.minH, dcase.maxH); 		
    for (int ih = 0; ih < hd().head.size(); ++ih)
		hd().head[ih] = FixHeading_180(hd().head[ih]);

	hd().dataFromW = true;
	
	if (dcase.irf) {
		hd().Tirf.resize(int(dcase.irfDuration/dcase.irfStep));
		for (int i = 0; i < hd().Tirf.size(); ++i) 
			hd().Tirf[i] = i*dcase.irfStep;
	}
	
	hd().c0.resize(3, hd().Nb);
	for (int i = 0; i < hd().Nb; ++i) 
		hd().c0.col(i) = dcase.bodies[i].c0.transpose();
		
	return true;
}

int NemohCase::GetNumArgs(const LineParser &f) {
	for (int i = 0; i < f.size(); ++i) {
		double num = ScanDouble(f.GetText(i));
		if (IsNull(num))
			return i;
	}
	return f.size();
}
	
void NemohCase::LoadFreeSurface(const FileInLine &in, const LineParser &f) {
	nFreeX = f.GetInt(0);		nFreeY = f.GetInt(1);	
	domainX = f.GetDouble(2);	domainY = f.GetDouble(3);
	if (nFreeX < 0)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect number of points in x direction %s"), f.GetText(0)));
	if (nFreeX > 0 && (nFreeY <= 0 || domainX < 0 || domainY < 0))
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect free surface elevation %s"), f.GetText()));	
}


void NemohCase::LoadKochin(const FileInLine &in, const LineParser &f) {
	nKochin = f.GetInt(0);	minK = f.GetDouble(1);	maxK = f.GetDouble(2);
	if (nKochin < 0)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect number of Kochin function directions %s"), f.GetText(0)));
	if (nKochin > 0) {
		if (minK < -360)
			throw Exc(in.Str() + "\n"  + Format(t_("Incorrect Kochin direction %s"), f.GetText(1)));
		if (maxK > 360)
			throw Exc(in.Str() + "\n"  + Format(t_("Incorrect Kochin direction %s"), f.GetText(2)));
		if (maxK <= minK)
			throw Exc(in.Str() + "\n"  + Format(t_("Minimum Kochin direction %s has to be lower than maximum direction %s"), f.GetText(1), f.GetText(2)));	
	}
}

bool NemohCase::Load(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	solver = NEMOH;
	
	String line;
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
	
	in.GetLine();
	f.GetLine();	
	rho = f.GetDouble(0);
	if (rho < 0 || rho > 10000)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect rho %s"), f.GetText(0)));
	
	f.GetLine();	
	g = f.GetDouble(0);
	if (g < 0 || g > 100)
		throw Exc(in.Str() + "\n" + Format(t_("Incorrect g %s"), f.GetText(0)));
	
	f.GetLine();	
	String sh = ToLower(f.GetText(0));
	if (sh == "inf")
		h = -1;
	else {
		h = ScanDouble(sh);
		if (h < 0 || h > 100000)
			throw Exc(in.Str() + "\n"  + Format(t_("Incorrect depth %s"), f.GetText(0)));
		else if (h == 0)
			h = -1;
	}
	f.GetLine();	xeff = f.GetDouble(0);	yeff = f.GetDouble(1);
	in.GetLine();
	f.GetLine();	int Nb = f.GetInt(0);
	if (Nb < 1 || Nb > 100)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect number of bodies %s"), f.GetText(0)));
	bodies.SetCount(Nb);
	for (int ib = 0; ib < Nb; ++ib) {
		int npoints, npanels;
		
		BEMBody &body = bodies[ib];
		in.GetLine();
		f.GetLine();	
		body.meshFile = f.GetText(0);
		f.GetLine();	
		npoints = f.GetInt(0);		
		npanels = f.GetInt(1);
		if (!FileExists(body.meshFile)) {
			body.meshFile = AppendFileNameX(GetFileFolder(fileName), body.meshFile);
			if (!FileExists(body.meshFile)) 
				BEM::PrintWarning(in.Str() + "\n"  + Format(t_("Mesh file '%s ' not found"), body.meshFile));
				//throw Exc(in.Str() + "\n"  + Format(t_("Mesh file '%s ' not found"), file));
		}
		if (npoints < 1 || npoints > 100000000)
			throw Exc(in.Str() + "\n"  + Format(t_("Incorrect number of points %s"), f.GetText(0)));
		if (npanels < 1 || npanels > 100000000)
			throw Exc(in.Str() + "\n"  + Format(t_("Incorrect number of panels %s"), f.GetText(1)));	
		f.GetLine();	
		body.ndof = f.GetInt(0);
		if (body.ndof < 0 || body.ndof > 6)
			throw Exc(in.Str() + "\n"  + Format(t_("Incorrect DOF number %s in body %d"), f.GetText(0), ib+1));
		for (int idf = 0; idf < body.ndof; ++idf) {
			f.GetLine();
			int type = f.GetInt(0);
			bool x = f.GetDouble(1) > 0;		
			bool y = f.GetDouble(2) > 0;
			bool z = f.GetDouble(3) > 0;
			double cx = f.GetDouble(4);
			double cy = f.GetDouble(5);
			double cz = f.GetDouble(6);
			if (type == 1) {
				if (x) 
					body.dof[BEM::SURGE] = true;
				else if (y)
					body.dof[BEM::SWAY] = true;
				else if (z)
					body.dof[BEM::HEAVE] = true;
			} else if (type == 2) {
				body.c0[0] = cx;
				body.c0[1] = cy;
				body.c0[2] = cz;
				if (x) 
					body.dof[BEM::ROLL] = true;
				else if (y)
					body.dof[BEM::PITCH] = true;
				else if (z)
					body.dof[BEM::YAW] = true;
			} else
				throw Exc(in.Str() + "\n"  + Format(t_("Incorrect DOF type %d set in body %d"), f.GetText(0), ib+1));
		}
		f.GetLine();	int nforces = f.GetInt(0);
		in.GetLine(nforces);	// Discarded
		f.GetLine();	int nadditional = f.GetInt(0);	
		in.GetLine(nadditional);// Discarded	
	}
	in.GetLine();
	f.GetLine();
	char type;
	int pos;
	if (f.size() > 3 && f.IsDouble(3)) {
		solver = NEMOHv3;
		pos = 1;
		switch (f.GetInt(0)) {
		case 1:		type = 'r';		break;
		case 2:		type = 'h';		break;
		case 3:		type = 't';		break;
		default:	throw Exc(in.Str() + "\n"  + Format(t_("Incorrect frequency type %d"), f.GetText(0)));
		}
	} else {	
		pos = 0;
		type = 'r';
	} 
	
	Nf = f.GetInt(pos + 0);	minF = f.GetDouble(pos + 1);	maxF = f.GetDouble(pos + 2);
	if (Nf < 1 || Nf > 1000)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect number of frequencies %s"), f.GetText(0)));
	if (minF < 0)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect frequency %s"), f.GetText(1)));
	if (maxF < minF)
		throw Exc(in.Str() + "\n"  + Format(t_("Minimum frequency %s has to be lower than maximum frequency %s"), f.GetText(1), f.GetText(2)));	
	
	if (type == 'h') {
		minF *= 2*M_PI;
		maxF *= 2*M_PI;
	} else if (type == 't') {
		minF = 2*M_PI/maxF;
		maxF = 2*M_PI/minF;
	}
	
	f.GetLine();	Nh = f.GetInt(0);	minH = f.GetDouble(1);	maxH = f.GetDouble(2);
	if (Nh < 1 || Nh > 1000)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect number of headings %s"), f.GetText(0)));
	if (minH < -360)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect direction %s"), f.GetText(1)));
	if (maxH > 360)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect direction %s"), f.GetText(2)));
	if (maxH < minH)
		throw Exc(in.Str() + "\n"  + Format(t_("Minimum direction %s has to be lower than maximum direction %s"), f.GetText(1), f.GetText(2)));	
	
	in.GetLine();
	f.GetLine();	irf = f.GetInt(0) > 0;	irfStep = f.GetDouble(1);	irfDuration = f.GetDouble(2);
	if (irf && irfStep <= 0)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect IRF step %s"), f.GetText(1)));
	if (irf && irfDuration <= irfStep)
		throw Exc(in.Str() + "\n"  + Format(t_("IRF step %s has to be lower than duration %s"), f.GetText(1), f.GetText(2)));	
	f.GetLine();	showPressure = f.GetInt(0) > 0;
	
	bool loadedFree = false, loadedKochin = false;
	f.GetLine();	
	if (GetNumArgs(f) == 4) {
		LoadFreeSurface(in, f);
		loadedFree = true;
	} else if (GetNumArgs(f) == 3) {
		LoadKochin(in, f);
		loadedKochin = true;
	} else
		throw Exc(in.Str() + "\n"  + Format(t_("Unexpected data %s"), f.GetText()));			
	
	f.GetLine();	
	if (GetNumArgs(f) == 4) {
		if (!loadedFree)
			LoadFreeSurface(in, f);
		else
			throw Exc(in.Str() + "\n"  + Format(t_("Free surface data is already loaded %s"), f.GetText()));				
	} else if (GetNumArgs(f) == 3) {
		if (!loadedKochin)
			LoadKochin(in, f);
		else
			throw Exc(in.Str() + "\n"  + Format(t_("Kochin data is already loaded %s"), f.GetText()));				
	} else if (GetNumArgs(f) > 0)
		throw Exc(in.Str() + "\n"  + Format(t_("Unexpected data %s"), f.GetText()));			
	else {
		if (!loadedFree) {
			nFreeX = nFreeY = 0;	
			domainX = domainY = 0;
		}
		if (!loadedKochin) {	
			nKochin = 0;	
			minK = maxK = 0;
		}
	}
	return true;
}

UVector<String> NemohCase::Check() const {
	UVector<String> ret;
	
	if (IsNull(rho) || rho < 0 || rho > 10000)
		 ret << Format(t_("Incorrect rho %s"), FormatDoubleEmpty(rho));
	if (IsNull(g) || g < 0 || g > 100)
		ret << Format(t_("Incorrect g %s"), FormatDoubleEmpty(g));

	if (irf) {
		if (IsNull(irfStep) || irfStep <= 0)
			ret << Format(t_("Incorrect IRF step %s"), FormatDoubleEmpty(irfStep));
		if (IsNull(irfDuration) || irfDuration <= irfStep)
			ret << Format(t_("IRF step %s has to be lower than duration %s"), FormatDoubleEmpty(irfStep), FormatDoubleEmpty(irfDuration));	
	}
	
	if (IsNull(nFreeX) || nFreeX < 0)
		ret << Format(t_("Incorrect number of points in x direction %s (0 for no free surface calculation)"), FormatIntEmpty(nFreeX));
	if (nFreeX > 0) {
		if (IsNull(nFreeY) || nFreeY <= 0)
			ret << Format(t_("Incorrect number of points in x direction %s"), FormatIntEmpty(nFreeY));
		if (IsNull(domainX) || domainX <= 0)
			ret << Format(t_("Incorrect free surface domain X %s"), FormatDoubleEmpty(domainX));
		if (IsNull(domainY) || domainY <= 0)
			ret << Format(t_("Incorrect free surface domain Y %s"), FormatDoubleEmpty(domainY));
	}
	
	if (IsNull(nKochin) || nKochin < 0)
		ret << Format(t_("Incorrect number of Kochin function directions %s"), FormatIntEmpty(nKochin));
	if (nKochin > 0) {
		if (IsNull(minK) || minK < -180)
			ret << Format(t_("Incorrect Kochin direction %s"), FormatDoubleEmpty(minK));
		if (IsNull(maxK) || maxK < 180)
			ret << Format(t_("Incorrect Kochin direction %s"),FormatDoubleEmpty(minK));
		if (maxK <= minK)
			Format(t_("Minimum Kochin direction %s has to be lower than maximum direction %s"), FormatDoubleEmpty(minK), FormatDoubleEmpty(maxK));	
	}	
	
	return ret;
}

void NemohCase::Save_Id(String folder) const {
	String fileName = AppendFileNameX(folder, "ID.dat");
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to create '%s'"), fileName));
	out << "1\n.";
}

void NemohCase::Save_Mesh_bat(String folder, String caseFolder, const UVector<String> &meshes, String meshName, bool bin) const {
	String fileName = AppendFileNameX(folder, "Mesh_cal.bat");
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to create '%s'"), fileName));

	String strBin;
	if (bin)
		strBin = AppendFileNameX(caseFolder.IsEmpty() ? "." : "..", "bin");
	
	if (meshes.size() == 1) {
		out << "copy Mesh_0.cal Mesh.cal\n";
		out << "\"" << AppendFileNameX(strBin, meshName) << "\"";
	} else {
		for (int i = 0; i < meshes.size(); ++i) {
			out << Format("copy Mesh_%d.cal Mesh.cal\n", i);
			out << "\"" << AppendFileNameX(strBin, meshName) << "\"\n";
			out << Format("ren mesh\\KH.dat KH_%d.dat\n", i);
		}
	}
}

void NemohCase::Save_Bat(String folder, String batname, String caseFolder, bool bin, 
		String preName, String hydroName, String solvName, String postName) const {
	String fileName = AppendFileNameX(folder, batname);
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to create '%s'"), fileName));
	
	out << Format("title %s in '%s'\n", solvName, caseFolder);
	
	if (!IsNull(caseFolder))
		out << "cd \"" << caseFolder << "\"\n";
	String strBin;
	if (bin)
		strBin = AppendFileNameX(caseFolder.IsEmpty() ? "." : "..", "bin");
	//out << "call Mesh_cal.bat\n";
	if (!preName.IsEmpty()) 
		out << "\"" << AppendFileNameX(strBin, preName) << "\"\n";
	if (!hydroName.IsEmpty()) 
		out << "\"" << AppendFileNameX(strBin, hydroName) << "\"\n";
	if (!solvName.IsEmpty()) {
		if (solvName == "capytaine")
			out << "\"" << solvName << "\"\n";
		else if (preName.IsEmpty()) 
			out << "\"" << AppendFileNameX(strBin, solvName) << "\" -all\n";
		else
			out << "\"" << AppendFileNameX(strBin, solvName) << "\"\n";
	}
	if (!postName.IsEmpty()) 
		out << "\"" << AppendFileNameX(strBin, postName) << "\"\n";
}

void NemohCase::Save_Input(String folder, int solver) const {
	if (solver == BEMCase::NEMOHv3) {
		String fileName = AppendFileNameX(folder, "input_solver.txt");
		FileOut out(fileName);
		if (!out.IsOpen())
			throw Exc(Format(t_("Impossible to create '%s'"), fileName));
		
		out << "2				! Gauss quadrature (GQ) surface integration, N^2 GQ Nodes, specify N(1,4)" << "\n";
		out << "0.001			! eps_zmin for determine minimum z of flow and source points of panel, zmin=eps_zmin*body_diameter" << "\n";
		out << "1 				! 0 GAUSS ELIM.; 1 LU DECOMP.: 2 GMRES	!Linear system solver" << "\n";
		out << "10 1e-5 1000  	! Restart parameter, Relative Tolerance, max iter -> additional input for GMRES";
	} else {
		String fileName = AppendFileNameX(folder, "Input.txt");
		FileOut out(fileName);
		if (!out.IsOpen())
			throw Exc(Format(t_("Impossible to create '%s'"), fileName));
		
		out << "--- Calculation parameters ------------------------------------------------------------------------------------" << "\n";
		out << "0          ! Indiq_solver   ! - ! Solver (0) Direct Gauss (1) GMRES (2) GMRES with FMM acceleration (2 not implemented yet)" << "\n";
		out << "20         ! IRES           ! - ! Restart parameter for GMRES" << "\n";
		out << "5.E-07     ! TOL_GMRES      ! - ! Stopping criterion for GMRES" << "\n";
		out << "100        ! MAXIT          ! - ! Maximum iterations for GMRES" << "\n";
		out << "1          ! Sav_potential  ! - ! Save potential for visualization";
	} 
}

String NemohHeader(String str) {
	String ret = "--- " << str << " ";
	ret << String('-', 130 - ret.GetCount());
	return ret;
}

String NemohField(String str, int length) {
	String ret = str;
	if (length > ret.GetCount())
		ret << String(' ', length - ret.GetCount());
	return ret + " ";
}

void NemohCase::SaveFolder(String folderBase, bool bin, int numCases, int, const BEM &bem, int solver) const {
	SaveFolder0(folderBase, bin, 1, bem, true, solver);
	if (numCases > 1)
		SaveFolder0(folderBase, bin, numCases, bem, false, solver);
}

void NemohCase::SaveFolder0(String folderBase, bool bin, int numCases, const BEM &bem, 
							bool deleteFolder, int solver) const {
	BeforeSave(folderBase, numCases, deleteFolder);
	
	#define MIN_F_NEMOH 0.01
	
	double fixminF = minF;
	if (fixminF < MIN_F_NEMOH)
		fixminF = MIN_F_NEMOH;
	
	UVector<int> valsf;
	int _nf;
	double _minf, _maxf;
	int ifr = 0;
	UVector<double> freqs;
	if (numCases > 1) { 
		LinSpaced(freqs, Nf, fixminF, maxF);
		valsf = NumSets(Nf, numCases);
	}
	
	String binResults = AppendFileNameX(folderBase, "bin");
	if (!DirectoryCreateX(binResults))
		throw Exc(Format(t_("Problem creating '%s' folder"), binResults));
		
	//String meshName = GetFileName(bem.nemohPathMesh);
	//String destMesh = AppendFileNameX(binResults, meshName);
	//if (!FileCopy(bem.nemohPathMesh, destMesh)) 
	//	throw Exc(Format(t_("Problem copying mesh binary from '%s'"), bem.nemohPathMesh));
	
	String preName, hydroName, solvName, postName, batName = "Nemoh";
	if (solver == BEMCase::CAPYTAINE) {
		solvName = "capytaine";
		batName = "Capytaine_bat";
	} else if (bin) {
		if (solver == BEMCase::NEMOHv115) {
			solvName = "nemoh.exe";
			if (!FileCopy(AppendFileNameX(bem.nemoh115Path, solvName), 
						  AppendFileNameX(binResults, solvName))) 
				throw Exc(Format(t_("Problem copying solver binary from '%s'"), bem.nemoh115Path));	
		} else if (solver == BEMCase::NEMOH) {
			preName = "preprocessor.exe";
			if (!FileCopy(AppendFileNameX(bem.nemohPath, preName), 
						  AppendFileNameX(binResults, preName))) 
				throw Exc(Format(t_("Problem copying preprocessor binary from '%s'"), bem.nemohPath));	
			solvName = "solver.exe";
			if (!FileCopy(AppendFileNameX(bem.nemohPath, solvName), 
						  AppendFileNameX(binResults, solvName))) 
				throw Exc(Format(t_("Problem copying solver binary from '%s'"), bem.nemohPath));
			postName = "postprocessor.exe";
			if (!FileCopy(AppendFileNameX(bem.nemohPath, postName), 
						  AppendFileNameX(binResults, postName))) 
				throw Exc(Format(t_("Problem copying postprocessor binary from '%s'"), bem.nemohPath));
		} else if (solver == BEMCase::NEMOHv3) {
			preName = "preproc.exe";
			if (!FileCopy(AppendFileNameX(bem.nemoh3Path, preName), 
						  AppendFileNameX(binResults, preName))) 
				throw Exc(Format(t_("Problem copying preproc binary from '%s'"), bem.nemoh3Path));	
			hydroName = "hydroscal.exe";
			if (!FileCopy(AppendFileNameX(bem.nemoh3Path, hydroName), 
						  AppendFileNameX(binResults, hydroName))) 
				throw Exc(Format(t_("Problem copying hydroscal binary from '%s'"), bem.nemoh3Path));
			solvName = "solver.exe";
			if (!FileCopy(AppendFileNameX(bem.nemoh3Path, solvName), 
						  AppendFileNameX(binResults, solvName))) 
				throw Exc(Format(t_("Problem copying solver binary from '%s'"), bem.nemoh3Path));
			postName = "postproc.exe";
			if (!FileCopy(AppendFileNameX(bem.nemoh3Path, postName), 
						  AppendFileNameX(binResults, postName))) 
				throw Exc(Format(t_("Problem copying postproc binary from '%s'"), bem.nemoh3Path));
		} else if (solver == BEMCase::CAPYTAINE) 
			solvName = "capytaine";

	} else {
		if (solver == BEMCase::NEMOHv115) 
			solvName = "nemoh.exe";
		else if (solver == BEMCase::NEMOH) {
			preName = "preprocessor.exe";
			solvName = "solver.exe";
			postName = "postprocessor.exe";
		} else if (solver == BEMCase::NEMOHv3) {
			preName = "preproc.exe";
			hydroName = "hydroscal.exe";
			solvName = "solver.exe";
			postName = "postproc.exe";
		}
	}
		
	//String sumcases;
	for (int i = 0; i < numCases; ++i) {
		String folder;
		if (numCases > 1) {
			folder = AppendFileNameX(folderBase, Format("%s_Part_%d", batName, i+1));
			if (!DirectoryCreateX(folder))
				throw Exc(Format(t_("Problem creating '%s' folder"), folder));
			//sumcases << " " << AppendFileNameX(folder, "Nemoh.cal");
			_minf = freqs[ifr];
			int deltaf = valsf[i];
			_maxf = freqs[ifr + deltaf - 1];
			_nf = deltaf;
			ifr += deltaf;
		} else {
			folder = folderBase;
			_nf = Nf;
			_minf = fixminF;
			_maxf = maxF;
		}
		
		String folderMech;
		if (solver == BEMCase::NEMOHv3) {
			folderMech = AppendFileNameX(folder, "mechanics");
			if (!DirectoryCreateX(folderMech))
				throw Exc(Format(t_("Problem creating '%s' folder"), folderMech));
		}
		
		if (solver != BEMCase::NEMOHv3)
			Save_Id(folder);
		Save_Input(folder, solver);
		String folderMesh = AppendFileNameX(folder, "mesh");
		if (!DirectoryCreateX(folderMesh))
			throw Exc(Format(t_("Problem creating '%s' folder"), folderMesh));
	
		UVector<String> meshes(bodies.size());
		UVector<int> nodes(bodies.size()), panels(bodies.size());
		for (int ib = 0; ib < bodies.size(); ++ib) {
			String name = GetFileName(bodies[ib].meshFile);
			name = RemoveAccents(name);
			name.Replace(" ", "_");
			String dest = AppendFileNameX(folderMesh, name);
			
			bool y0z = false, x0z = false;
			Mesh mesh;
			String err = Mesh::Load(mesh, bodies[ib].meshFile, rho, g, false, y0z, x0z);
			if (!err.IsEmpty()) {
				err = Mesh::Load(mesh, AppendFileNameX(folderMesh, GetFileName(bodies[ib].meshFile)), rho, g, false, y0z, x0z);
				if (!err.IsEmpty()) {
					throw Exc(err);
				}
			}
			mesh.c0 = clone(bodies[ib].c0);
			mesh.cg = clone(bodies[ib].cg);
			mesh.AfterLoad(rho, g, true, false);
			Mesh::SaveAs(mesh, dest, Mesh::NEMOH_DAT, Mesh::UNDERWATER, rho, g, false, x0z, nodes[ib], panels[ib]);
			
			if (solver == BEMCase::NEMOHv3) {
				Save_Mesh_cal(folder, bodies.size() == 1 ? -1 : ib, 
						bodies[ib].meshFile, mesh, panels[ib], x0z, bodies[ib].cg, rho, g);
				String inertiaName;
				if (bodies.size() == 1)
					inertiaName = "inertia.dat";
				else
					inertiaName = Format("inertia_%d.dat", i);
				Nemoh::Save_6x6(bodies[ib].mass, AppendFileNameX(folderMech, inertiaName));
				meshes[ib] = GetFileTitle(bodies[ib].meshFile);
			} else {
				String khName;
				if (bodies.size() == 1)
					khName = "KH.dat";
				else
					khName = Format("KH_%d.dat", i);
				static_cast<NemohMesh&>(mesh).SaveKH(AppendFileNameX(folderMesh, khName));			
			}
		}
		Save_Cal(folder, _nf, _minf, _maxf, nodes, panels, solver);
				
		String folderResults = AppendFileNameX(folder, "results");
		if (!DirectoryCreateX(folderResults))
			throw Exc(Format(t_("Problem creating '%s' folder"), folderResults));
		
		if (bin && solver != BEMCase::NEMOHv3 && !GetFileName(bem.nemohPathGREN).IsEmpty()) {
			String destGREN = AppendFileNameX(folder, GetFileName(bem.nemohPathGREN));
			if (!FileCopy(bem.nemohPathGREN, destGREN)) 
				throw Exc(Format(t_("Problem copying gren file '%s'"), bem.nemohPathGREN));
		}
		
		if (numCases > 1) {
			String caseFolder = Format("%s_Part_%d", batName, i+1);
			//Save_Mesh_bat(folder, caseFolder, meshes, meshName, bin || BEMCase::CAPYTAINE); 
			Save_Bat(folderBase, Format("%s_Part_%d.bat", batName, i+1), caseFolder, bin, preName, hydroName, solvName, postName);
		} else {
			//Save_Mesh_bat(folder, Null, meshes, meshName, bin || BEMCase::CAPYTAINE);
			Save_Bat(folder, Format("%s.bat", batName), Null, bin, preName, hydroName, solvName, postName);
		}
	}
}

void NemohCase::Save_Mesh_cal(String folder, int ib, String meshFile, Mesh &mesh, int npanels, bool x0z, Eigen::Vector3d cg, double rho, double g) const {
	String title = ForceExt(GetFileTitle(meshFile), ".pmsh");
	title = RemoveAccents(title);
	title.Replace(" ", "_");
	Mesh::SaveAs(mesh, AppendFileNameX(folder, "Mesh", title), 
				Mesh::NEMOH_PRE, Mesh::UNDERWATER, rho, g, false, x0z);
	
	String fileName;
	if (ib < 0)
		fileName = AppendFileNameX(folder, "Mesh.cal");
	else
		fileName = AppendFileNameX(folder, Format("Mesh_%d.cal", ib));
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'"), fileName));
	
	out << Format("%20s ! Mesh file\n", "\"" + title + "\"");
	out << (x0z ? "1" : "0") << "                    ! 1 if a symmetry about (x0z) is used. 0 otherwise\n";	
	out << "0    0 	             ! Possible translation about x axis (first number) and y axis (second number)\n";
	out << Format("%s %s %s ! Coordinates of gravity centre\n", FDS(cg[0], 6, true), FDS(cg[1], 6, true), 
																FDS(cg[2], 6, true));
	out << Format("%6d               ! Target for the number of panels in refined mesh\n", npanels);
	out << "2\n";
	out << "0                    ! Mesh z-axis translation from origin\n";
	out << "1                    ! Mesh scale\n";
	out << Format("%s               ! Water density (kg/m3)\n", FDS(rho, 6, true));
	out << Format("%s               ! Gravity (m/s2)", FDS(g, 6, true));
	
}
	
void NemohCase::Save_Cal(String folder, int _nf, double _minf, double _maxf, const UVector<int> &nodes, 
		const UVector<int> &panels, int solver) const {
	String fileName = AppendFileNameX(folder, "Nemoh.cal");
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'"), fileName));
	
	int cp = 28;	
	out << NemohHeader("Environment - Created with BEMRosetta") << "\n";
	out << NemohField(Format("%f", rho), cp) 		   << "! RHO             ! KG/M**3   ! Fluid specific volume" << "\n";
	out << NemohField(Format("%f", g), cp)   		   << "! G               ! M/S**2    ! Gravity " << "\n";
	double seth = h;
	if (h < 0)
		seth = 0;
	out << NemohField(Format("%f", seth), cp)   	   << "! DEPTH           ! M         ! Water depth" << "\n";
	out << NemohField(Format("%f %f", xeff, yeff), cp) << "! XEFF YEFF       ! M         ! Wave measurement point" << "\n";
	
	out << NemohHeader("Description of floating bodies") << "\n";
	out << NemohField(Format("%d", bodies.size()), cp) << "! Number of bodies" << "\n";
	
	for (int i = 0; i < bodies.size(); ++i) {
		const BEMBody &b = bodies[i];
		out << NemohHeader(Format("Body %d", i+1)) << "\n";	
		String name = GetFileName(b.meshFile);
		name = RemoveAccents(name);
		name.Replace(" ", "_");
		String file = AppendFileNameX("mesh", name);
		
		out << NemohField(Format("%s", file), cp) << "! Name of mesh file" << "\n";
		out << NemohField(Format("%d %d", nodes[i], panels[i]), cp) << "! Number of points and number of panels" << "\n";	
		out << NemohField(Format("%d", b.GetNDOF()), cp) << "! Number of degrees of freedom" << "\n";	
		if (b.dof[BEM::SURGE])
			out << NemohField("1 1. 0. 0. 0. 0. 0.", cp) << "! Surge" << "\n";	
		if (b.dof[BEM::SWAY])
			out << NemohField("1 0. 1. 0. 0. 0. 0.", cp) << "! Sway" << "\n";	
		if (b.dof[BEM::HEAVE])
			out << NemohField("1 0. 0. 1. 0. 0. 0.", cp) << "! Heave" << "\n";	
		if (b.dof[BEM::ROLL])
			out << NemohField(Format("2 1. 0. 0. %.2f %.2f %.2f", b.c0[0], b.c0[1], b.c0[2]), cp) << "! Roll about a point" << "\n";	
		if (b.dof[BEM::PITCH])
			out << NemohField(Format("2 0. 1. 0. %.2f %.2f %.2f", b.c0[0], b.c0[1], b.c0[2]), cp) << "! Pitch about a point" << "\n";	
		if (b.dof[BEM::YAW])		
			out << NemohField(Format("2 0. 0. 1. %.2f %.2f %.2f", b.c0[0], b.c0[1], b.c0[2]), cp) << "! Yaw about a point" << "\n";	
		out << NemohField(Format("%d", b.GetNDOF()), cp) << "! Number of resulting generalised forces" << "\n";	
		if (b.dof[BEM::SURGE])
			out << NemohField("1 1. 0. 0. 0. 0. 0.", cp) << "! Force in x direction" << "\n";	
		if (b.dof[BEM::SWAY])
			out << NemohField("1 0. 1. 0. 0. 0. 0.", cp) << "! Force in y direction" << "\n";	
		if (b.dof[BEM::HEAVE])
			out << NemohField("1 0. 0. 1. 0. 0. 0.", cp) << "! Force in z direction" << "\n";	
		if (b.dof[BEM::ROLL])
			out << NemohField(Format("2 1. 0. 0. %.2f %.2f %.2f", b.c0[0], b.c0[1], b.c0[2]), cp) << "! Moment force in x direction about a point" << "\n";	
		if (b.dof[BEM::PITCH])
			out << NemohField(Format("2 0. 1. 0. %.2f %.2f %.2f", b.c0[0], b.c0[1], b.c0[2]), cp) << "! Moment force in y direction about a point" << "\n";	
		if (b.dof[BEM::YAW])		
			out << NemohField(Format("2 0. 0. 1. %.2f %.2f %.2f", b.c0[0], b.c0[1], b.c0[2]), cp) << "! Moment force in z direction about a point" << "\n";	
		out << NemohField("0", cp) << "! Number of lines of additional information" << "\n";
	}
	out << NemohHeader("Load cases to be solved") << "\n";
	if (solver == NEMOHv3)
		out << NemohField(Format("1 %d %f %f", _nf, _minf, _maxf), cp) << "! Freq type 1,2,3=[rad/s,Hz,s], Number of wave frequencies/periods, Min, and Max" << "\n";
	else	
		out << NemohField(Format("%d %f %f", _nf, _minf, _maxf), cp)   << "! Number of wave frequencies, Min, and Max (rad/s)" << "\n";
	
	//double _minH = !isCapy ? minH : ToRad(minH);	// 29/12/2021 Capytaine issue
	//double _maxH = !isCapy ? maxH : ToRad(maxH);
	out << NemohField(Format("%d %f %f", Nh, minH, maxH), cp) << "! Number of wave directions, Min and Max (degrees)" << "\n";
	
	out << NemohHeader("Post processing") << "\n";
	out << NemohField(Format("%4<d %.2f %.2f", irf ? 1 : 0, irfStep, irfDuration), cp) << "! IRF                    ! IRF calculation (0 for no calculation), time step and duration" << "\n";
	out << NemohField(Format("%d", showPressure ? 1 : 0), cp) << "! Show pressure" << "\n";	
	out << NemohField(Format("%4<d %.2f %.2f", nKochin, minK, maxK), cp) << "! Kochin function        ! Number of directions of calculation (0 for no calculations), Min and Max (degrees)" << "\n";
	out << NemohField(Format("%4<d %4<d %.2f %.2f", nFreeX, nFreeY, domainX, domainY), cp) << "! Free surface elevation ! Number of points in x direction (0 for no calculations) and y direction and dimensions of domain in x and y direction" << "\n";
	
	if (solver == NEMOHv3) {
		out << NemohField("0						! Response Amplitude Operator (RAO), 0 no calculation, 1 calculated", cp) << "\n";	
		out << NemohField("1						! output freq type, 1,2,3=[rad/s,Hz,s]", cp) << "\n";	 
		out << NemohHeader("QTF") << "\n";	
		out << NemohField("0", cp) << "! QTF flag, 1 is calculated" << "\n";	
	}
	out << "---";
}

bool Nemoh::Load_Inf(String fileName) {
	if (hd().Nb != 1)
		throw Exc(Format(t_("SeaFEM_Nemoh only allows one body, found %d"), hd().Nb));

	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
				
	hd().cg.setConstant(3, 1, NaNDouble);
	hd().cb.setConstant(3, 1, NaNDouble);
	//hd().c0.setConstant(3, 1, NaNDouble);
	hd().Vo.SetCount(1, NaNDouble);
	hd().C.SetCount(1);
	hd().C[0].setConstant(6, 6, NaNDouble);   
	
	double minimumDirectionAngle = 0;
	
	String line;
	while(!in.IsEof()) {
		line = in.GetLine();
		int pos;
		if ((pos = line.FindAfter("XG [m]=")) >= 0) 
			hd().cg(0, 0) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("YG [m]=")) >= 0) 
			hd().cg(1, 0) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("ZG [m]=")) >= 0) 
			hd().cg(2, 0) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("XC [m]=")) >= 0) 
			hd().cb(0, 0) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("YC [m]=")) >= 0) 
			hd().cb(1, 0) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("ZC [m]=")) >= 0) 
			hd().cb(2, 0) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("Displacement [m3]=")) >= 0) 
			hd().Vo[0] = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [3][3] [N/m]=")) >= 0) 
			hd().C[0](2, 2) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [3][4] [N/rad]=")) >= 0) 
			hd().C[0](2, 3) = hd().C[0](3, 2) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [3][5] [N/rad]=")) >= 0) 
			hd().C[0](2, 4) = hd().C[0](4, 2) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [4][4] [Nm/rad]=")) >= 0) 
			hd().C[0](3, 3) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [4][5] [Nm/rad]=")) >= 0) 
			hd().C[0](3, 4) = hd().C[0](4, 3) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [4][6] [Nm/rad]=")) >= 0) 
			hd().C[0](3, 5) = hd().C[0](5, 3) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [5][5] [Nm/rad]=")) >= 0) 
			hd().C[0](4, 4) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [5][6] [Nm/rad]=")) >= 0) 
			hd().C[0](4, 5) = hd().C[0](5, 4) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("Minimum direction angle =")) >= 0)
			minimumDirectionAngle = ScanDouble(line.Mid(pos))*180/M_PI;
	}
	
	for (int ih = 0; ih < hd().head.size(); ++ih)
		hd().head[ih] -= minimumDirectionAngle;

	return true;	
}

bool Nemoh::Load_Hydrostatics(String subfolder) {
	return Load_Hydrostatics_static(AppendFileNameX(folder, subfolder), hd().Nb, hd().cg, hd().cb, hd().Vo);	
}
	
bool Nemoh::Load_Hydrostatics_static(String subfolder, int Nb, MatrixXd &cg, MatrixXd &cb, UVector<double> &Vo) {
	cg.setConstant(3, Nb, NaNDouble);
	cb.setConstant(3, Nb, NaNDouble);
	Vo.SetCount(Nb, NaNDouble);
	
	for (int ib = 0; ib < Nb; ++ib) {
	    String fileHydro;
	    if (Nb == 1)
	        fileHydro = AppendFileNameX(subfolder, "Hydrostatics.dat");
	    else
	        fileHydro = AppendFileNameX(subfolder, Format("Hydrostatics_%d.dat", ib));
	    
	    FileInLine in(fileHydro);
	    if (!in.IsOpen())
	        return false;
	    
	    LineParser f(in);
	    f.IsSeparator = IsTabSpace;
	    for (int i = 0; i < 3 && !in.IsEof(); ++i) {
			f.Load(in.GetLine());
			cg(i, ib) = f.GetDouble(6);
			cb(i, ib) = f.GetDouble(2);
	    }
	    if (!f.IsEof() && f.GetCount() > 0) {
			f.Load(in.GetLine());
		    Vo[ib] = f.GetDouble(2); 		
	    }
	}
	return true;
}

void Nemoh::Save_Hydrostatics(String subfolder) const {
	Save_Hydrostatics_static(subfolder, hd().Nb, hd().cg, hd().cb, hd().Vo);	
}

void Nemoh::Save_Hydrostatics_static(String folder, int Nb, const MatrixXd &cg, const MatrixXd &cb, const UVector<double> &Vo) {
	for (int ib = 0; ib < Nb; ++ib) {
	    String fileHydro;
	    if (Nb == 1)
	        fileHydro = AppendFileNameX(folder, "Hydrostatics.dat");
	    else
	        fileHydro = AppendFileNameX(folder, Format("Hydrostatics_%d.dat", ib));
	    
	    FileOut out(fileHydro);
		if (!out.IsOpen())
			throw Exc(Format(t_("Impossible to create '%s'"), fileHydro));
	    
		out << Format(" XF =   %.3f - XG =   %.3f\n", cb(0, ib), cg(0, ib));
		out << Format(" YF =   %.3f - YG =   %.3f\n", cb(1, ib), cg(1, ib));
		out << Format(" ZF =   %.3f - ZG =   %.3f\n", cb(2, ib), cg(2, ib));
		if (Vo.size() == 3) 
			out << Format(" Displacement =  %.7G\n", Vo[ib]);
	}
}

bool Nemoh::Load_KH(String subfolder) {
	hd().C.SetCount(hd().Nb);
	for (int ib = 0; ib < hd().Nb; ++ib) {
	    String fileKH;
		if (hd().Nb == 1) 
			fileKH = AppendFileNameX(folder, subfolder, "KH.dat");
		else 
			fileKH = AppendFileNameX(folder, subfolder, Format("KH_%d.dat", ib));
	    
	    if (!Load_6x6(hd().C[ib], fileKH))
	        return false;	        
	}
	return true;
}

bool Nemoh::Load_Inertia(String subfolder) {
	hd().M.SetCount(hd().Nb);
	for (int ib = 0; ib < hd().Nb; ++ib) {
	    String file;
		if (hd().Nb == 1) 
			file = AppendFileNameX(folder, subfolder, "inertia.dat");
		else 
			file = AppendFileNameX(folder, subfolder, Format("inertia_%d.dat", ib));
	    
	    if (!Load_6x6(hd().M[ib], file))
	        return false;	        
	}
	return true;
}

bool Nemoh::Load_LinearDamping(String subfolder) {
	hd().Dlin = MatrixXd::Zero(6*hd().Nb, 6*hd().Nb);
	
	for (int ib = 0; ib < hd().Nb; ++ib) {
	    String file;
		if (hd().Nb == 1) 
			file = AppendFileNameX(folder, subfolder, "inertia.dat");
		else 
			file = AppendFileNameX(folder, subfolder, Format("inertia_%d.dat", ib));
	    
	    MatrixXd m;
	    if (!Load_6x6(m, file))
	        return false;
	    hd().Dlin.block<6,6>(ib*6, ib*6) = m;	        
	}
	return true;
}

bool Nemoh::Load_12(String fileName, bool isSum, Function <bool(String, int)> Status) {
	hd().dimen = true;
	if (IsNull(hd().len))
		hd().len = 1;
	
	String ext    = isSum ? "summation" : "difference";
	double phmult = isSum ? 1 : -1;		// Difference is conjugate-symmetric
	
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	Status(Format("Loading %s base data", ext), 0);
	
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
	
	FileInLine::Pos fpos = in.GetPos();
	f.GetLine();		
	if (!IsNull(f.GetDouble_nothrow(0)))
		in.SeekPos(fpos);		// No header, rewind
	else
		fpos = in.GetPos();		// Avoid header
	
	UArray<UArray<UArray<MatrixXcd>>> &qtf = isSum ? hd().qtfsum : hd().qtfdif;		
	
	qtf.Clear();
	
	UVector<double> w;
    UArray<std::complex<double>> head;
	
	int Nb = 0;
	int nline = 0;
	while (!in.IsEof()) {
		f.GetLine();
		nline++;
		double freq = f.GetDouble(0);
		FindAdd(w, freq);
		double freq2 = f.GetDouble(1);
		FindAdd(w, freq2);
		double hd1 = f.GetDouble(2);
		double hd2 = f.GetDouble(3);
		FindAdd(head, std::complex<double>(hd1, hd2));	
		Nb = max(Nb, 1 + (f.GetInt(4)-1)/6);
	}
	if (IsNull(hd().Nb))
		hd().Nb = Nb;
	else {
		if (hd().Nb < Nb)
			throw Exc(Format(t_("Number of bodies loaded is lower than previous (%d != %d)"), hd().Nb, Nb));
	}
	
	if (hd().names.IsEmpty())
		hd().names.SetCount(hd().Nb);
		
	int Nf = w.size();
	int Nh = head.size();
	if (Nh == 0)
		throw Exc(Format(t_("Wrong format in Wamit file '%s'. No headings found"), hd().file));
	
	Hydro::Initialize_QTF(qtf, Nb, Nh, Nf);
	
	Status(Format("Loading %s data", ext), 20);
	
	in.SeekPos(fpos);
	int iline = 0;
	while (!in.IsEof()) {
		f.GetLine();
		iline++;
		
		if (Status && !(iline%(nline/10)) && !Status(Format("Loading %s", ext), 20 + (80*iline)/nline))
			throw Exc(t_("Stop by user"));

		int ifr1 = Find(w, f.GetDouble(0));
		int ifr2 = Find(w, f.GetDouble(1));
		int ih = Find(head, std::complex<double>(f.GetDouble(2), f.GetDouble(3)));
		int idf = f.GetInt(4)-1;
		int ib = int(idf/6);
		idf -= 6*ib;
		double ma = f.GetDouble(5)*hd().rho*hd().g;
		double ph = -ToRad(f.GetDouble(6));			//-Phase to follow Wamit
		qtf[ib][ih][idf](ifr1, ifr2) = std::polar(ma, ph);
		qtf[ib][ih][idf](ifr2, ifr1) = std::polar(ma, ph*phmult);
	}
	
	Copy(w, hd().qw);
	Copy(head, hd().qh);
	
	hd().qtfdataFromW = !(w[0] > w[1]);
	
	if (!hd().qtfdataFromW) 
		for (int i = 0; i < hd().qw.size(); ++i)
   			hd().qw(i) = 2*M_PI/hd().qw(i);
	
	for (int i = 0; i < hd().qh.size(); ++i)
   		hd().qh(i) = std::complex<double>(hd().qh(i).real(), hd().qh(i).imag()); //FixHeading_180(hd().qh(i).real()), FixHeading_180(hd().qh(i).imag()));
	
	return true;
}	

bool Nemoh::Load_QTF(String subfolder, Function <bool(String, int)> Status) {
	String file;
	
	file = AppendFileNameX(folder, subfolder, "OUT_QTFM_N.dat");
	if (!Load_12(file, false, Status))
		return false;
	file = AppendFileNameX(folder, subfolder, "OUT_QTFP_N.dat");
	if (!Load_12(file, true, Status))
		return false;
	
	return true;
}
	
bool Nemoh::Load_6x6(Eigen::MatrixXd &C, String file) {
    FileInLine in(file);
	if (!in.IsOpen()) 
        return false;

	C.setConstant(6, 6, 0);
    
    LineParser f(in);
    f.IsSeparator = IsTabSpace;
	for (int i = 0; i < 6 && !in.IsEof(); ++i) {
		f.Load(in.GetLine());
		for (int ifr = 0; ifr < 6; ++ifr)
			C(i, ifr) = f.GetDouble(ifr);
	}	
	return true;
}

bool Nemoh::Save_KH(String folder) const {
	for (int ib = 0; ib < hd().Nb; ++ib) {
	    String file;
		if (hd().Nb == 1) 
			file = AppendFileNameX(folder, "Mesh", "KH.dat");
		else 
			file = AppendFileNameX(folder, "Mesh", Format("KH_%d.dat", ib));
	    
	    if (!Save_6x6(hd().C[ib], file))
	        return false;
	}
	return true;
}

bool Nemoh::Save_Inertia(String folder) const {
	for (int ib = 0; ib < hd().Nb; ++ib) {
	    String file;
		if (hd().Nb == 1) 
			file = AppendFileNameX(folder, "mechanics", "inertia.dat");
		else 
			file = AppendFileNameX(folder, "mechanics", Format("inertia_%d.dat", ib));
	    
	    if (!Save_6x6(hd().M[ib], file))
	        return false;
	}
	return true;
}

bool Nemoh::Save_6x6(const Eigen::MatrixXd &C, String file) {
    FileOut out(file);
	if (!out.IsOpen()) 
        return false;

	for (int i = 0; i < 6; ++i) {
		for (int ifr = 0; ifr < 6; ++ifr)
			out << "  " << FormatE(C(i, ifr), 13);
		out << "\n";
	}
	return true;
}

bool Nemoh::Load_Radiation(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	String line;
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
	in.GetLine();
	while(!in.IsEof()) {
		line = in.GetLine();
	    if (line.Find("Motion of body") >= 0 || line.Find("dof_") >= 0 || 
	    	line.Find("Surge") >= 0 || line.Find("Sway") >= 0 || line.Find("Heave") >= 0 || 
	    	line.Find("Roll") >= 0 || line.Find("Pitch") >= 0 || line.Find("Yaw") >= 0)
	        break;
	}
	hd().Initialize_AB(hd().A);
	hd().Initialize_AB(hd().B);
	
	for (int ibody = 0; ibody < hd().Nb; ++ibody) {
		for (int idof = 0; idof < 6; ++idof) {
			if (dcase.IsDof(ibody, idof)) {
				for (int ifr = 0; ifr < hd().Nf; ++ifr) {	
					f.Load(in.GetLine());
					int col = 1;
					for (int idof2 = 0; idof2 < 6; ++idof2) {			
						if (dcase.IsDof(ibody, idof2)) {
							hd().A[ibody*6 + idof][ibody*6 + idof2][ifr] = f.GetDouble(col++);
		        			hd().B[ibody*6 + idof][ibody*6 + idof2][ifr] = f.GetDouble(col++);
			        	}
					}
				}
		    	if (idof < 6)
		    		in.GetLine();
	    	}
		}
	}
	return true;
}

bool Nemoh::Load_Excitation(String folder) {	
	return Load_Forces(hd().ex, folder, "ExcitationForce.tec");
}

bool Nemoh::Load_Diffraction(String folder) {
	return Load_Forces(hd().sc, folder, "DiffractionForce.tec");
}

bool Nemoh::Load_FroudeKrylov(String folder) {
	return Load_Forces(hd().fk, folder, "FKForce.tec");
}

bool Nemoh::Load_Forces(Hydro::Forces &fc, String nfolder, String fileName) {
	FileInLine in(AppendFileNameX(nfolder, "Results", fileName));
	if (!in.IsOpen())
		return false;
	String line;
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
	in.GetLine();
	UVector<UVector<int>> dof;
	dof.SetCount(hd().Nb);
	while(!in.IsEof()) {
		line = in.GetLine();
		if (line.StartsWith("Zone") || line.StartsWith("angle"))
	        break;
	}
	hd().Initialize_Forces(fc);
	for (int ih = 0; ih < hd().Nh; ++ih) {
		int ifr = 0;
		while(!in.IsEof()) {
			line = in.GetLine();
			if (line.StartsWith("Zone") || line.StartsWith("angle"))
				break;
			f.Load(line);
			int il = 0;
			for (int ib = 0; ib < hd().Nb; ++ib) {
				for (int ibdof = 0; ibdof < 6; ++ibdof) {
					if (dcase.IsDof(ib, ibdof)) {
						if (ifr >= hd().Nf)
							throw Exc(in.Str() + "\n"  + t_("Number of frequencies higher than the defined in Nemoh.cal file"));		
						double ma = f.GetDouble(1 + 2*il);	
						double ph = -f.GetDouble(1 + 2*il + 1); //-Phase to follow Wamit
						fc.force[ih](ifr, ibdof+6*ib) = std::polar(ma, ph); 
						il++; 
					}
				}
			}
			ifr++;
		}
	}
	return true;
}

bool Nemoh::Load_IRF(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	String line;
	LineParser f(in);	
	f.IsSeparator = IsTabSpace;
	hd().Ainf.setConstant(hd().Nb*6, hd().Nb*6, 0);
	//int ibodydof = 0;
	hd().Kirf.SetCount(hd().Nb*6); 	// Initialize Kirf		
    for (int i = 0; i < hd().Nb*6; ++i) {
    	hd().Kirf[i].SetCount(hd().Nb*6); 			 
   		for (int j = 0; j < hd().Nb*6; ++j)
			hd().Kirf[i][j].setConstant(hd().Tirf.size(), NaNDouble);
    }
    while(!in.IsEof()) {
		line = in.GetLine();	
		if (line.Find("Zone t=") >= 0) 
			break;
	}
	for (int ib = 0; ib < hd().Nb; ++ib) {
		for (int ibdof = 0; ibdof < 6; ++ibdof) {
			if (dcase.IsDof(ib, ibdof)) {
				for (int iNt = 0; iNt < hd().Tirf.size(); ++iNt) {
					f.Load(in.GetLine());
					int col = 1;
					for (int idf = 0; idf < 6; ++idf) {
						if (dcase.IsDof(ib, idf)) {
							hd().Ainf(ibdof, idf) = f.GetDouble(col++);
							hd().Kirf[ib*6+ibdof][ib*6+idf][iNt] = f.GetDouble(col++);
						}
					}
				}
				in.GetLine();
			}
		}
	}
	return true;
}

void Nemoh::Save(String ) {
	throw Exc("Option not implemented");
}		

