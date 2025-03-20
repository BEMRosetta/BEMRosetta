// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#include <STEM4U/Utility.h>

String Nemoh::Load(String file, Function <bool(String, int)> Status, double) {
	try {
		String ext = GetFileExt(file); 
		String folder = GetFileFolder(file);
		
		if (ext == ".tec") {
			String folderTitle = GetFileName(folder);
			if (ToLower(folderTitle) != "results") 
				return Format(t_(".tec file '%s' should have to be in 'results' folder"), file);
			bool found = false;
			String upperFolder = GetUpperFolder(folder);
			for (FindFile ff(AFX(upperFolder, "*.*")); ff; ++ff) {
				if (ff.IsFile()) {
					if (ToLower(ff.GetName()) == "nemoh.cal") {
						file = ff.GetPath();
						found = true;
						break;
					}
				}
			}
			if (!found)
				return Format(t_("nemoh.cal file not found in '%s' folder"), upperFolder);
		}
	
		if (ext == ".cal" || ext == ".tec")
			dt.solver = Hydro::NEMOH;
		else
			dt.solver = Hydro::SEAFEM_NEMOH;
	
		dt.file = file;
		dt.name = GetFileTitle(GetFileFolder(file));
		dt.len = 1;
		dt.dimen = true;
		dt.Nb = Null;
	
		String fileCal;
		BEM::Print("\n\n" + Format(t_("Loading '%s'"), file));
		if (dt.solver == Hydro::NEMOH) 
			fileCal = file;
		else 
			fileCal = AFX(folder, "Nemoh_output/Nemoh.cal");
		if (!Load_Cal(fileCal)) 
			return Format(t_("File '%s' not found"), fileCal);
		
		String fileRad, folderForces;
		if (dt.solver != Hydro::SEAFEM_NEMOH) {
			BEM::Print(S("\n- ") + t_("Hydrostatics file(s) 'Mesh/Hydrostatics*.dat'"));
			if (!Load_Hydrostatics(folder, "Mesh"))
				BEM::PrintWarning(S(": ** Mesh/Hydrostatics*.dat ") + t_("Not found") + "**");
			BEM::Print(S("\n- ") + t_("KH file(s) 'Mesh/KH*.dat'"));
			if (!Load_KH(folder, "Mesh"))
				BEM::PrintWarning(S(": ** Mesh/KH ") + t_("Not found") + "**");
			
			dynamic_cast<Wamit *>(this)->Load_frc2(ForceExtSafer(fileCal, ".frc"));
			
			fileRad = AFX(folder, "Results", "RadiationCoefficients.tec");
			folderForces = folder;
		} else {
			if (!Load_Inf(file)) 
				throw Exc(Format(t_("File '%s' not found"), file));

			fileRad = AFX(folder, "Nemoh_output/Results", "RadiationCoefficients.tec");
			folderForces = AFX(folder, "Nemoh_output");
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
		
		if (dt.solver == Hydro::NEMOHv3) {
			Load_Inertia(folder, "mechanics");
			Load_LinearDamping(folder, "mechanics");
			Load_QTF(folder, AFX("results", "qtf"), Status);
		}
		
		UVector<int> idsRemove;
		for (int ih = 0; ih < dt.Nh; ++ih) {
			int id = FindDelta(dt.head, dt.head[ih], 0.001, ih+1);
			if (id > 0)
				idsRemove << ih;
		}
		DeleteHeadings(idsRemove);
		
		if (dt.solver == Hydro::NEMOH) {
			BEM::Print(S("\n- ") + t_("IRF file(s) 'IRF.tec'"));
			if (!Load_IRF(AFX(folder, "Results", "IRF.tec")))
				BEM::PrintWarning(S(": ** IRF.tec ") + t_("Not found") + "**");
		}
		if (IsNull(dt.Nb))
			return t_("No body found");
	} catch (Exc e) {
		BEM::PrintError(Format("\n%s: %s", t_("Error"), e));
		//dt.lastError = e;
		return e;
	}
	return String();
}

bool Nemoh::Load_Cal(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	if (dt.solver == UNKNOWN)
		dt.solver = Hydro::NEMOH;
	//dt.dataFromW = true;
	
	String line;
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
	
	in.GetLine();
	f.GetLine();	
	dt.rho = f.GetDouble(0);
	if (dt.rho < 0 || dt.rho > 10000)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect rho %s"), f.GetText(0)));
	
	f.GetLine();	
	dt.g = f.GetDouble(0);
	if (dt.g < 0 || dt.g > 100)
		throw Exc(in.Str() + "\n" + Format(t_("Incorrect g %s"), f.GetText(0)));
	
	f.GetLine();	
	String sh = ToLower(f.GetText(0));
	if (sh == "inf")
		dt.h = -1;
	else {
		dt.h = ScanDouble(sh);
		if (dt.h > 100000)
			throw Exc(in.Str() + "\n"  + Format(t_("Incorrect depth %s"), f.GetText(0)));
		else if (dt.h <= 0)
			dt.h = -1;
	}
	f.GetLine();	dt.x_w = f.GetDouble(0);	dt.y_w = f.GetDouble(1);
	in.GetLine();
	f.GetLine();	dt.Nb = f.GetInt(0);
	if (dt.Nb < 1 || dt.Nb > 100)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect number of bodies %s"), f.GetText(0)));
	dt.msh.SetCount(dt.Nb);
	
	UVector<String> names;
	for (int ib = 0; ib < dt.Nb; ++ib) {
		int npoints, npanels;
		
		Body &body = dt.msh[ib];
		in.GetLine();
		f.GetLine();	
		body.dt.fileName = f.GetText(0);
		f.GetLine();	
		npoints = f.GetInt(0);		
		npanels = f.GetInt(1);
		if (!FileExists(body.dt.fileName)) {
			body.dt.fileName = AFX(GetFileFolder(fileName), body.dt.fileName);
			if (!FileExists(body.dt.fileName)) 
				BEM::PrintWarning(in.Str() + "\n"  + Format(t_("Mesh file '%s ' not found"), body.dt.fileName));
				//throw Exc(in.Str() + "\n"  + Format(t_("Mesh file '%s ' not found"), file));
		}
		body.dt.name = GetFileTitle(body.dt.fileName);
		if (Find(names, body.dt.name) >= 0) {
			body.dt.name << "_" << (ib+1);
			names << body.dt.name;
		}
		if (npoints < 1 || npoints > 1e8)
			throw Exc(in.Str() + "\n"  + Format(t_("Incorrect number of points %s"), f.GetText(0)));
		if (npanels < 1 || npanels > 1e8)
			throw Exc(in.Str() + "\n"  + Format(t_("Incorrect number of panels %s"), f.GetText(1)));	
		f.GetLine();	
		/*body.*/ int ndof = f.GetInt(0);
		if (/*body.*/ndof < 0 || /*body.*/ndof > 6)
			throw Exc(in.Str() + "\n"  + Format(t_("Incorrect DOF number %s in body %d"), f.GetText(0), ib+1));
		for (int idf = 0; idf < /*body.*/ndof; ++idf) {
			f.GetLine();
			int type = f.GetInt(0);
			//bool x = f.GetDouble(1) > 0;		
			//bool y = f.GetDouble(2) > 0;
			//bool z = f.GetDouble(3) > 0;
			double cx = f.GetDouble(4);
			double cy = f.GetDouble(5);
			double cz = f.GetDouble(6);
			if (type == 1) {
			//	if (x) 
			//		body.dof[BEM::SURGE] = true;
			//	else if (y)
			//		body.dof[BEM::SWAY] = true;
			//	else if (z)
			;//		body.dof[BEM::HEAVE] = true;
			} else if (type == 2) {
				body.dt.c0[0] = cx;
				body.dt.c0[1] = cy;
				body.dt.c0[2] = cz;
			//	if (x) 
			//		body.dof[BEM::ROLL] = true;
			//	else if (y)
			//		body.dof[BEM::PITCH] = true;
			//	else if (z)
			//		body.dof[BEM::YAW] = true;
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
		dt.solver = Hydro::NEMOHv3;
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
	
	dt.Nf = f.GetInt(pos + 0);	
	double minF = f.GetDouble(pos + 1);	
	double maxF = f.GetDouble(pos + 2);
	if (dt.Nf < 1 || dt.Nf > 1000)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect number of frequencies %s"), f.GetText(0)));
	if (minF < 0)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect frequency %s"), f.GetText(1)));
	//if (maxF == minF)
	//	throw Exc(in.Str() + "\n"  + Format(t_("Minimum frequency %s has to be different than maximum frequency %s"), f.GetText(1), f.GetText(2)));	
	//if (maxF < minF)
	//	throw Exc(in.Str() + "\n"  + Format(t_("Minimum frequency %s has to be lower than maximum frequency %s"), f.GetText(1), f.GetText(2)));	
	
	if (type == 'h') {
		minF *= 2*M_PI;
		maxF *= 2*M_PI;
	} else if (type == 't') {
		minF = 2*M_PI/maxF;
		maxF = 2*M_PI/minF;
	}
	
	if (minF > maxF)
		Swap(minF, maxF);
	
	LinSpaced(dt.w, dt.Nf, minF, maxF); 
 	/*dt.T.SetCount(dt.Nf);
    for (int i = 0; i < dt.Nf; ++i) 
		dt.T[i] = 2*M_PI/dt.w[i];  */
	
	f.GetLine();	
	dt.Nh = f.GetInt(0);	
	double minH = f.GetDouble(1);	
	double maxH = f.GetDouble(2);
	if (dt.Nh < 1 || dt.Nh > 1000)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect number of headings %s"), f.GetText(0)));
	if (minH < -360)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect direction %s"), f.GetText(1)));
	if (maxH > 360)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect direction %s"), f.GetText(2)));
	if (maxH < minH)
		throw Exc(in.Str() + "\n"  + Format(t_("Minimum direction %s has to be lower than maximum direction %s"), f.GetText(1), f.GetText(2)));	
	
    LinSpaced(dt.head, dt.Nh, minH, maxH); 		
    //for (int ih = 0; ih < dt.head.size(); ++ih)
	//	dt.head[ih] = FixHeading_180(dt.head[ih]);	
	
	in.GetLine();
	f.GetLine();	
	int irf = f.GetInt(0) > 0;		
	if (irf) {
		double irfStep = f.GetDouble(1);	
		double irfDuration = f.GetDouble(2);
		if (irfStep <= 0)
			throw Exc(in.Str() + "\n"  + Format(t_("Incorrect IRF step %s"), f.GetText(1)));
		if (irfDuration <= irfStep)
			throw Exc(in.Str() + "\n"  + Format(t_("IRF step %s has to be lower than duration %s"), f.GetText(1), f.GetText(2)));	
	
		dt.Tirf.resize(int(irfDuration/irfStep));
		for (int i = 0; i < dt.Tirf.size(); ++i) 
			dt.Tirf[i] = i*irfStep;
	}
	return true;
}
	

void Nemoh::Save_Id(String folder) const {
	String fileName = AFX(folder, "ID.dat");
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to create '%s'"), fileName));
	out << "1\n.";
}

void Nemoh::Save_Body_bat(String folder, String caseFolder, const UVector<String> &meshes, String meshName, bool bin) const {
	String fileName = AFX(folder, "Mesh_cal.bat");
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to create '%s'"), fileName));

	String strBin;
	if (bin)
		strBin = AFX(caseFolder.IsEmpty() ? "." : "..", "bin");
	
	if (meshes.size() == 1) {
		out << "copy Mesh_0.cal Mesh.cal\n";
		out << "\"" << AFX(strBin, meshName) << "\"";
	} else {
		for (int i = 0; i < meshes.size(); ++i) {
			out << Format("copy Mesh_%d.cal Mesh.cal\n", i);
			out << "\"" << AFX(strBin, meshName) << "\"\n";
			out << Format("ren mesh\\KH.dat KH_%d.dat\n", i);
		}
	}
}

void Nemoh::Save_Bat(String folder, String batname, String caseFolder, bool bin, 
		String preName, String hydroName, String solvName, String postName, int numThreads) const {
	String fileName = AFX(folder, batname);
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to create '%s'"), fileName));
	
	if (!IsEmpty(caseFolder))
		out << Format("title \"%s in '%s'\"\n", solvName, caseFolder);
	else
		out << Format("title %s\n", solvName);
	
	out << "\necho Start: \%date\% \%time\% > time.txt\n";
	
	if (!IsNull(caseFolder))
		out << "cd \"" << caseFolder << "\"\n";
	String strBin;
	if (bin)
		strBin = AFX(caseFolder.IsEmpty() ? "." : "..", "bin");
	//out << "call Mesh_cal.bat\n";
	if (!preName.IsEmpty()) 
		out << "\"" << AFX(strBin, preName) << "\"\n";
	if (!hydroName.IsEmpty()) 
		out << "\"" << AFX(strBin, hydroName) << "\"\n";
	if (!solvName.IsEmpty()) {
		if (solvName == "capytaine") {
			if (!IsNull(numThreads) && numThreads > 0) 
				out << "set OMP_NUM_THREADS=" << numThreads << "\n"
					<< "set MKL_NUM_THREADS=" << numThreads << "\n";
			if (!IsEmpty(Bem().pythonEnv)) 
				out << Format("call activate %s\n", Bem().pythonEnv); 
			out << "\"" << solvName << "\"\n";
		} else if (preName.IsEmpty()) 
			out << "\"" << AFX(strBin, solvName) << "\" -all\n";
		else
			out << "\"" << AFX(strBin, solvName) << "\"\n";
	}
	if (!postName.IsEmpty()) 
		out << "\"" << AFX(strBin, postName) << "\"\n";
	
	out << "\necho End:   \%date\% \%time\% >> time.txt\n";
}

void Nemoh::Save_Input(String folder, int solver) const {
	if (solver == Hydro::NEMOHv3) {
		String fileName = AFX(folder, "input_solver.txt");
		FileOut out(fileName);
		if (!out.IsOpen())
			throw Exc(Format(t_("Impossible to create '%s'"), fileName));
		
		out << "2				! Gauss quadrature (GQ) surface integration, N^2 GQ Nodes, specify N(1,4)" << "\n";
		out << "0.001			! eps_zmin for determine minimum z of flow and source points of panel, zmin=eps_zmin*body_diameter" << "\n";
		out << "1 				! 0 GAUSS ELIM.; 1 LU DECOMP.: 2 GMRES	!Linear system solver" << "\n";
		out << "10 1e-5 1000  	! Restart parameter, Relative Tolerance, max iter -> additional input for GMRES";
	} else if (solver != Hydro::CAPYTAINE) {
		String fileName = AFX(folder, "Input.txt");
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

void Nemoh::SaveCase(String folderBase, bool bin, int numCases, int solver, int numThreads, bool x0z, bool y0z, const UArray<Body> &lids) const {
	SaveFolder0(folderBase, bin, 1, true, solver, numThreads, x0z, y0z, lids);
	if (numCases > 1)
		SaveFolder0(folderBase, bin, numCases, false, solver, numThreads, x0z, y0z, lids);
}

void Nemoh::SaveFolder0(String folderBase, bool bin, int numCases, bool deleteFolder, int solver, int numThreads, bool x0z, bool y0z, const UArray<Body> &lids) const {
	BeforeSaveCase(folderBase, numCases, deleteFolder);

	UVector<int> valsf;
	int ifr = 0;
	if (numCases > 1) 
		valsf = NumSets(dt.Nf, numCases);
	
	String binResults = AFX(folderBase, "bin");
	if (!DirectoryCreateX(binResults))
		throw Exc(Format(t_("Problem creating '%s' folder"), binResults));
	
	String preName, hydroName, solvName, postName, batName = "Nemoh";
	if (solver == Hydro::CAPYTAINE) {
		solvName = "capytaine";
		batName = GetFileTitle(folderBase);
		bin = false;
	} else if (bin) {
		if (solver == Hydro::NEMOHv115) {
			solvName = "nemoh.exe";
			if (!FileCopy(AFX(Bem().nemoh115Path, solvName), 
						  AFX(binResults, solvName))) 
				throw Exc(Format(t_("Problem copying solver binary from '%s'"), Bem().nemoh115Path));	
		} else if (solver == Hydro::NEMOH) {
			preName = "preprocessor.exe";
			if (!FileCopy(AFX(Bem().nemohPath, preName), 
						  AFX(binResults, preName))) 
				throw Exc(Format(t_("Problem copying preprocessor binary from '%s'"), Bem().nemohPath));	
			solvName = "solver.exe";
			if (!FileCopy(AFX(Bem().nemohPath, solvName), 
						  AFX(binResults, solvName))) 
				throw Exc(Format(t_("Problem copying solver binary from '%s'"), Bem().nemohPath));
			postName = "postprocessor.exe";
			if (!FileCopy(AFX(Bem().nemohPath, postName), 
						  AFX(binResults, postName))) 
				throw Exc(Format(t_("Problem copying postprocessor binary from '%s'"), Bem().nemohPath));
		} else if (solver == Hydro::NEMOHv3) {
			preName = "preproc.exe";
			if (!FileCopy(AFX(Bem().nemoh3Path, preName), 
						  AFX(binResults, preName))) 
				throw Exc(Format(t_("Problem copying preproc binary from '%s'"), Bem().nemoh3Path));	
			hydroName = "hydroscal.exe";
			if (!FileCopy(AFX(Bem().nemoh3Path, hydroName), 
						  AFX(binResults, hydroName))) 
				throw Exc(Format(t_("Problem copying hydroscal binary from '%s'"), Bem().nemoh3Path));
			solvName = "solver.exe";
			if (!FileCopy(AFX(Bem().nemoh3Path, solvName), 
						  AFX(binResults, solvName))) 
				throw Exc(Format(t_("Problem copying solver binary from '%s'"), Bem().nemoh3Path));
			postName = "postproc.exe";
			if (!FileCopy(AFX(Bem().nemoh3Path, postName), 
						  AFX(binResults, postName))) 
				throw Exc(Format(t_("Problem copying postproc binary from '%s'"), Bem().nemoh3Path));
		} else if (solver == Hydro::CAPYTAINE) 
			solvName = "capytaine";

	} else {
		if (solver == Hydro::NEMOHv115) 
			solvName = "nemoh.exe";
		else if (solver == Hydro::NEMOH) {
			preName = "preprocessor.exe";
			solvName = "solver.exe";
			postName = "postprocessor.exe";
		} else if (solver == Hydro::NEMOHv3) {
			preName = "preproc.exe";
			hydroName = "hydroscal.exe";
			solvName = "solver.exe";
			postName = "postproc.exe";
		}
	}
		
	for (int i = 0; i < numCases; ++i) {
		String folder;
		UVector<double> freqs;
		if (numCases > 1) {
			folder = AFX(folderBase, Format("%s_Part_%d", batName, i+1));
			if (!DirectoryCreateX(folder))
				throw Exc(Format(t_("Problem creating '%s' folder"), folder));
			Upp::Block(dt.w, freqs, ifr, valsf[i]);
			ifr += valsf[i];
		} else {
			folder = folderBase;
			freqs = clone(dt.w);
		}
		
		String folderMech;
		if (solver == Hydro::NEMOHv3) {
			folderMech = AFX(folder, "mechanics");
			if (!DirectoryCreateX(folderMech))
				throw Exc(Format(t_("Problem creating '%s' folder"), folderMech));
		}
		
		if (solver != Hydro::NEMOHv3)
			Save_Id(folder);
		Save_Input(folder, solver);
		String folderMesh = AFX(folder, "mesh");
		if (!DirectoryCreateX(folderMesh))
			throw Exc(Format(t_("Problem creating '%s' folder"), folderMesh));
	
		UVector<int> numNodes(dt.Nb), numPanels(dt.Nb);
		for (int ib = 0; ib < dt.msh.size(); ++ib) {
			String dest = AFX(folderMesh, Format(t_("Body_%d.dat"), ib+1));
			
			if (solver == Hydro::NEMOHv3) {
				Save_Body_cal(folder, dt.msh.size() == 1 ? -1 : ib, 
						dest, dt.msh[ib], dt.symY, dt.msh[ib].dt.cg, dt.rho, dt.g);
				String inertiaName;
				if (dt.msh.size() == 1)
					inertiaName = "inertia.dat";
				else
					inertiaName = Format("inertia_%d.dat", i);
				Nemoh::Save_6x6(dt.msh[ib].dt.M, AFX(folderMech, inertiaName));
			} else {
				String khName;
				if (dt.msh.size() == 1)
					khName = "KH.dat";
				else
					khName = Format("KH_%d.dat", i);
				static_cast<const NemohBody&>(dt.msh[ib]).SaveKH(AFX(folderMesh, khName));			
			}
			numNodes[ib] = dt.msh[ib].dt.under.nodes.size();
			numPanels[ib] = dt.msh[ib].dt.under.panels.size();
		}
		Save_Cal(folder, freqs, numNodes, numPanels, solver, y0z, x0z, lids);
				
		String folderResults = AFX(folder, "results");
		if (!DirectoryCreateX(folderResults))
			throw Exc(Format(t_("Problem creating '%s' folder"), folderResults));
		
		if (bin && solver != Hydro::NEMOHv3 && !GetFileName(Bem().nemohPathGREN).IsEmpty()) {
			String destGREN = AFX(folder, GetFileName(Bem().nemohPathGREN));
			if (!FileCopy(Bem().nemohPathGREN, destGREN)) 
				throw Exc(Format(t_("Problem copying gren file '%s'"), Bem().nemohPathGREN));
		}
		
		if (numCases > 1) {
			String caseFolder = Format("%s_Part_%d", batName, i+1);
			Save_Bat(folderBase, Format("%s_Part_%d.bat", batName, i+1), caseFolder, bin, preName, hydroName, solvName, postName, numThreads);
		} else 
			Save_Bat(folder, Format("%s.bat", batName), Null, bin, preName, hydroName, solvName, postName, numThreads);
	}
}

void Nemoh::Save_Body_cal(String folder, int ib, String meshFile, const Body &mesh, bool x0z, const Point3D &cg, double rho, double g) const {
	String title = ForceExt(GetFileTitle(meshFile), ".pmsh");
	title = RemoveAccents(title);
	title.Replace(" ", "_");
	Body::SaveAs(mesh, AFX(folder, "Mesh", title), 
				Body::NEMOH_PRE, Body::ALL, rho, g, false, x0z);
	
	int npanels = mesh.dt.mesh.panels.size();
	String fileName;
	if (ib < 0)
		fileName = AFX(folder, "Mesh.cal");
	else
		fileName = AFX(folder, Format("Mesh_%d.cal", ib));
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
	out << Format("%s               ! Water density (kg/m³)\n", FDS(rho, 6, true));
	out << Format("%s               ! Gravity (m/s2)", FDS(g, 6, true));
}
	
void Nemoh::Save_Cal(String folder, const UVector<double> &freqs, const UVector<int> &nodes, 
		const UVector<int> &panels, int solver, bool y0z, bool x0z, const UArray<Body> &lids) const {
	String fileName = AFX(folder, "Nemoh.cal");
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'"), fileName));
	
	int Nf = freqs.size();
	double minF = First(freqs);
	double maxF = Last(freqs);
	
	int cp = 28;	
	out << NemohHeader("Environment - Created with BEMRosetta") << "\n";
	out << NemohField(Format("%f", dt.rho), cp) 		   << "! RHO             ! KG/M**3   ! Fluid specific volume" << "\n";
	out << NemohField(Format("%f", dt.g), cp)   		   << "! G               ! M/S**2    ! Gravity " << "\n";
	double seth = dt.h;
	if (seth < 0)
		seth = 0;
	out << NemohField(Format("%f", seth), cp)   	   << "! DEPTH           ! M         ! Water depth" << "\n";
	out << NemohField(Format("%f %f", dt.x_w, dt.y_w), cp) << "! XEFF YEFF       ! M         ! Wave measurement point" << "\n";
	
	out << NemohHeader("Description of floating bodies") << "\n";
	out << NemohField(Format("%d", dt.msh.size()), cp) << "! Number of bodies" << "\n";
	
	String folderMesh = AFX(folder, "mesh");
	DirectoryCreate(folderMesh);
	
	for (int ib = 0; ib < dt.msh.size(); ++ib) {
		const Body &b = dt.msh[ib];
		String name = Format("Body_%d.dat", ib+1);
		
		Body under = clone(dt.msh[ib]);
		bool isLid = solver == Hydro::NEMOHv3 && lids.size() > ib && !lids[ib].dt.mesh.panels.IsEmpty();
		if (isLid) 
			under.Append(lids[ib].dt.mesh, dt.rho, dt.g);
			
		int nNodes, nPanels;
		Body::SaveAs(under, AFX(folderMesh, name), Body::NEMOH_DAT, Body::ALL, dt.rho, dt.g, y0z, x0z, nNodes, nPanels);
				
		out << NemohHeader(name) << "\n";	
		
		String file = AFX("mesh", name);
		
		out << NemohField(Format("%s", file), cp) << "! Name of mesh file" << "\n";
		out << NemohField(Format("%d %d", nNodes, nPanels), cp) << "! Number of points and number of panels" << "\n";	
		out << NemohField(Format("%d", 6/*b.GetNDOF()*/), cp) << "! Number of degrees of freedom" << "\n";	
		//if (b.dof[BEM::SURGE])
			out << NemohField("1 1. 0. 0. 0. 0. 0.", cp) << "! Surge" << "\n";	
		//if (b.dof[BEM::SWAY])
			out << NemohField("1 0. 1. 0. 0. 0. 0.", cp) << "! Sway" << "\n";	
		//if (b.dof[BEM::HEAVE])
			out << NemohField("1 0. 0. 1. 0. 0. 0.", cp) << "! Heave" << "\n";	
		//if (b.dof[BEM::ROLL])
			out << NemohField(Format("2 1. 0. 0. %.2f %.2f %.2f", b.dt.c0[0], b.dt.c0[1], b.dt.c0[2]), cp) << "! Roll about a point" << "\n";	
		//if (b.dof[BEM::PITCH])
			out << NemohField(Format("2 0. 1. 0. %.2f %.2f %.2f", b.dt.c0[0], b.dt.c0[1], b.dt.c0[2]), cp) << "! Pitch about a point" << "\n";	
		//if (b.dof[BEM::YAW])		
			out << NemohField(Format("2 0. 0. 1. %.2f %.2f %.2f", b.dt.c0[0], b.dt.c0[1], b.dt.c0[2]), cp) << "! Yaw about a point" << "\n";	
		out << NemohField(Format("%d", 6/*b.GetNDOF()*/), cp) << "! Number of resulting generalised forces" << "\n";	
		//if (b.dof[BEM::SURGE])
			out << NemohField("1 1. 0. 0. 0. 0. 0.", cp) << "! Force in x direction" << "\n";	
		//if (b.dof[BEM::SWAY])
			out << NemohField("1 0. 1. 0. 0. 0. 0.", cp) << "! Force in y direction" << "\n";	
		//if (b.dof[BEM::HEAVE])
			out << NemohField("1 0. 0. 1. 0. 0. 0.", cp) << "! Force in z direction" << "\n";	
		//if (b.dof[BEM::ROLL])
			out << NemohField(Format("2 1. 0. 0. %.2f %.2f %.2f", b.dt.c0[0], b.dt.c0[1], b.dt.c0[2]), cp) << "! Moment force in x direction about a point" << "\n";	
		//if (b.dof[BEM::PITCH])
			out << NemohField(Format("2 0. 1. 0. %.2f %.2f %.2f", b.dt.c0[0], b.dt.c0[1], b.dt.c0[2]), cp) << "! Moment force in y direction about a point" << "\n";	
		//if (b.dof[BEM::YAW])		
			out << NemohField(Format("2 0. 0. 1. %.2f %.2f %.2f", b.dt.c0[0], b.dt.c0[1], b.dt.c0[2]), cp) << "! Moment force in z direction about a point" << "\n";	
		out << NemohField("0", cp) << "! Number of lines of additional information" << "\n";
	}
	out << NemohHeader("Load cases to be solved") << "\n";
	if (solver == Hydro::NEMOHv3)
		out << NemohField(Format("1 %d %f %f", Nf, minF, maxF), cp) << "! Freq type 1,2,3=[rad/s,Hz,s], Number of wave frequencies/periods, Min, and Max" << "\n";
	else	
		out << NemohField(Format("%d %f %f", Nf, minF, maxF), cp)   << "! Number of wave frequencies, Min, and Max (rad/s)" << "\n";
	
	out << NemohField(Format("%d %f %f", dt.Nh, First(dt.head), Last(dt.head)), cp) << "! Number of wave directions, Min and Max (degrees)" << "\n";
	
	out << NemohHeader("Post processing") << "\n";
	if (dt.Tirf.size() > 0)
		out << NemohField(Format("%4<d %.2f %.2f", 1, dt.Tirf[1]-dt.Tirf[0], Last(dt.Tirf)), cp) << "! IRF                    ! IRF calculation (0 for no calculation), time step and duration" << "\n";
	else
		out << NemohField(Format("%4<d %.2f %.2f", 0, 0, 0), cp) << "! IRF                    ! IRF calculation (0 for no calculation), time step and duration" << "\n";
	out << NemohField(Format("%d", 0), cp) << "! Show pressure" << "\n";	
	out << NemohField(Format("%4<d %.2f %.2f", 0, 0, 0), cp) << "! Kochin function        ! Number of directions of calculation (0 for no calculations), Min and Max (degrees)" << "\n";
	out << NemohField(Format("%4<d %4<d %.2f %.2f", 0, 0, 0, 0), cp) << "! Free surface elevation ! Number of points in x direction (0 for no calculations) and y direction and dimensions of domain in x and y direction" << "\n";
	
	if (solver == Hydro::NEMOHv3) {
		out << NemohField("0						! Response Amplitude Operator (RAO), 0 no calculation, 1 calculated", cp) << "\n";	
		out << NemohField("1						! output freq type, 1,2,3=[rad/s,Hz,s]", cp) << "\n";	 
		out << NemohHeader("QTF") << "\n";	
		out << NemohField("0", cp) << "! QTF flag, 1 is calculated" << "\n";	
	}
	out << "---";
}

bool Nemoh::Load_Inf(String fileName) {
	if (dt.Nb != 1)
		throw Exc(Format(t_("SeaFEM_Nemoh only allows one body, found %d"), dt.Nb));

	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	dt.msh.SetCount(dt.Nb);	
	dt.msh[0].dt.C.setConstant(6, 6, NaNDouble);   
	
	double minimumDirectionAngle = 0;
	
	String line;
	while(!in.IsEof()) {
		line = in.GetLine();
		int pos;
		if ((pos = line.FindAfter("XG [m]=")) >= 0) 
			dt.msh[0].dt.cg.x = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("YG [m]=")) >= 0) 
			dt.msh[0].dt.cg.y = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("ZG [m]=")) >= 0) 
			dt.msh[0].dt.cg.z = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("XC [m]=")) >= 0) 
			dt.msh[0].dt.cb.x = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("YC [m]=")) >= 0) 
			dt.msh[0].dt.cb.y = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("ZC [m]=")) >= 0) 
			dt.msh[0].dt.cb.z = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("Displacement [m³]=")) >= 0) 
			dt.msh[0].dt.Vo = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [3][3] [N/m]=")) >= 0) 
			dt.msh[0].dt.C(2, 2) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [3][4] [N/rad]=")) >= 0) 
			dt.msh[0].dt.C(2, 3) = dt.msh[0].dt.C(3, 2) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [3][5] [N/rad]=")) >= 0) 
			dt.msh[0].dt.C(2, 4) = dt.msh[0].dt.C(4, 2) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [4][4] [Nm/rad]=")) >= 0) 
			dt.msh[0].dt.C(3, 3) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [4][5] [Nm/rad]=")) >= 0) 
			dt.msh[0].dt.C(3, 4) = dt.msh[0].dt.C(4, 3) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [4][6] [Nm/rad]=")) >= 0) 
			dt.msh[0].dt.C(3, 5) = dt.msh[0].dt.C(5, 3) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [5][5] [Nm/rad]=")) >= 0) 
			dt.msh[0].dt.C(4, 4) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [5][6] [Nm/rad]=")) >= 0) 
			dt.msh[0].dt.C(4, 5) = dt.msh[0].dt.C(5, 4) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("Minimum direction angle =")) >= 0)
			minimumDirectionAngle = ScanDouble(line.Mid(pos))*180/M_PI;
	}
	
	for (int ih = 0; ih < dt.head.size(); ++ih)
		dt.head[ih] -= minimumDirectionAngle;

	return true;	
}

bool Nemoh::Load_Hydrostatics(String folder, String subfolder) {
	return Load_Hydrostatics_static(AFX(folder, subfolder), dt.Nb, dt.msh);	
}
	
bool Nemoh::Load_Hydrostatics_static(String subfolder, int Nb, UArray<Body> &msh) {
	for (int ib = 0; ib < Nb; ++ib) {
	    String fileHydro;
	    if (Nb == 1)
	        fileHydro = AFX(subfolder, "Hydrostatics.dat");
	    else
	        fileHydro = AFX(subfolder, Format("Hydrostatics_%d.dat", ib));
	    
	    FileInLine in(fileHydro);
	    if (!in.IsOpen())
	        return false;
	    
	    LineParser f(in);
	    f.IsSeparator = IsTabSpace;
	    for (int i = 0; i < 3 && !in.IsEof(); ++i) {
			String line = f.GetLine();
			line.Replace(" ", "");
			line = ToUpper(line);
			int ibu = line.FindAfter("=");
			if (ibu > 0)
				msh[ib].dt.cb[i] = ScanDouble(line.Mid(ibu));
			int ig = line.FindAfter("=", ibu); 
			if (ig > 0)
				msh[ib].dt.cg[i] = ScanDouble(line.Mid(ig));
	    }
	    if (!f.IsEof() && f.GetCount() > 0) {
			f.Load(in.GetLine());
		    msh[ib].dt.Vo = f.GetDouble_nothrow(2); 		
	    }
	}
	return true;
}

void Nemoh::Save_Hydrostatics(String subfolder) const {
	Save_Hydrostatics_static(subfolder, dt.Nb, dt.msh);	
}

void Nemoh::Save_Hydrostatics_static(String folder, int Nb, const UArray<Body> &msh) {
	for (int ib = 0; ib < Nb; ++ib) {
		const Body &m = msh[ib];
		if (m.dt.Vo > 0) {
		    String fileHydro;
		    if (Nb == 1)
		        fileHydro = AFX(folder, "Hydrostatics.dat");
		    else
		        fileHydro = AFX(folder, Format("Hydrostatics_%d.dat", ib));
		    
		    FileOut out(fileHydro);
			if (!out.IsOpen())
				throw Exc(Format(t_("Impossible to create '%s'"), fileHydro));
		   
			out << Format(" XF =   %.3f - XG =   %.3f\n", m.dt.cb.x, m.dt.cg.x);
			out << Format(" YF =   %.3f - YG =   %.3f\n", m.dt.cb.y, m.dt.cg.y);
			out << Format(" ZF =   %.3f - ZG =   %.3f\n", m.dt.cb.z, m.dt.cg.z);
			out << Format(" Displacement =  %.7G\n", m.dt.Vo);
		}
	}
}

bool Nemoh::Load_KH(String folder, String subfolder) {
	for (int ib = 0; ib < dt.Nb; ++ib) {
	    String fileKH;
		if (dt.Nb == 1) 
			fileKH = AFX(folder, subfolder, "KH.dat");
		else 
			fileKH = AFX(folder, subfolder, Format("KH_%d.dat", ib));
	    
	    if (!Load_6x6(dt.msh[ib].dt.C, fileKH))
	        return false;	        
	}
	return true;
}

bool Nemoh::Load_Inertia(String folder, String subfolder) {
	for (int ib = 0; ib < dt.Nb; ++ib) {
	    String file;
		if (dt.Nb == 1) 
			file = AFX(folder, subfolder, "inertia.dat");
		else 
			file = AFX(folder, subfolder, Format("inertia_%d.dat", ib));
	    
	    if (!Load_6x6(dt.msh[ib].dt.M, file))
	        return false;	        
	}
	return true;
}

bool Nemoh::Load_LinearDamping(String folder, String subfolder) {
	for (int ib = 0; ib < dt.Nb; ++ib) {
	    String file;
		if (dt.Nb == 1) 
			file = AFX(folder, subfolder, "Badd.dat");
		else 
			file = AFX(folder, subfolder, Format("Badd_%d.dat", ib));
	    
	    MatrixXd m;
	    if (!Load_6x6(m, file))
	        return false;
	    dt.msh[ib].dt.Dlin = m;	        
	}
	return true;
}

bool Nemoh::Load_12(String fileName, bool isSum, Function <bool(String, int)> Status) {
	dt.dimen = true;
	if (IsNull(dt.len))
		dt.len = 1;
	
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
	
	UArray<UArray<UArray<MatrixXcd>>> &qtf = isSum ? dt.qtfsum : dt.qtfdif;		
	
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
	
	if (nline == 0)		// Only the header, but no data
		return false;
	
	if (IsNull(dt.Nb) || Nb == 0)
		dt.Nb = Nb;
	else {
		if (dt.Nb < Nb)
			throw Exc(Format(t_("Number of bodies loaded is lower than previous (%d != %d)"), dt.Nb, Nb));
	}
	
	if (dt.msh.IsEmpty())
		dt.msh.SetCount(dt.Nb);
		
	int Nf = w.size();
	int Nh = head.size();
	if (Nh == 0)
		throw Exc(Format(t_("Wrong format in QTF file '%s'. No headings found"), dt.file));
	
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
		double ma = f.GetDouble(5)*dt.rho*dt.g;
		double ph = -ToRad(f.GetDouble(6));			//-Phase to follow Wamit
		qtf[ib][ih][idf](ifr1, ifr2) = std::polar(ma, ph);
		qtf[ib][ih][idf](ifr2, ifr1) = std::polar(ma, ph*phmult);
	}
	
	::Copy(w, dt.qw);
	::Copy(head, dt.qhead);
	
	dt.qtfdataFromW = !(w[0] > w[1]);
	
	if (!dt.qtfdataFromW) 
		for (int i = 0; i < dt.qw.size(); ++i)
   			dt.qw(i) = 2*M_PI/dt.qw(i);
	
	for (int i = 0; i < dt.qhead.size(); ++i)
   		dt.qhead(i) = std::complex<double>(dt.qhead(i).real(), dt.qhead(i).imag()); //FixHeading_180(dt.qh(i).real()), FixHeading_180(dt.qh(i).imag()));
	
	return true;
}	

bool Nemoh::Load_QTF(String folder, String subfolder, Function <bool(String, int)> Status) {
	String file;
	
	file = AFX(folder, subfolder, "OUT_QTFM_N.dat");
	if (!Load_12(file, false, Status))
		return false;
	file = AFX(folder, subfolder, "OUT_QTFP_N.dat");
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
	for (int ib = 0; ib < dt.Nb; ++ib) {
	    String file;
		if (dt.Nb == 1) 
			file = AFX(folder, "Mesh", "KH.dat");
		else 
			file = AFX(folder, "Mesh", Format("KH_%d.dat", ib));
	    
	    if (!Save_6x6(dt.msh[ib].dt.C, file))
	        return false;
	}
	return true;
}

bool Nemoh::Save_Inertia(String folder) const {
	for (int ib = 0; ib < dt.Nb; ++ib) {
	    String file;
		if (dt.Nb == 1) 
			file = AFX(folder, "mechanics", "inertia.dat");
		else 
			file = AFX(folder, "mechanics", Format("inertia_%d.dat", ib));
	    
	    if (!Save_6x6(dt.msh[ib].dt.M, file))
	        return false;
	}
	return true;
}

bool Nemoh::Save_6x6(const Eigen::MatrixXd &C, String file) {
	if (C.size() < 36)
		return false;
	
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
	Initialize_AB(dt.A);
	Initialize_AB(dt.B);
	
	for (int ibody = 0; ibody < dt.Nb; ++ibody) {
		for (int idof = 0; idof < 6; ++idof) {
			//if (dcase.IsDof(ibody, idof)) {
				for (int ifr = 0; ifr < dt.Nf; ++ifr) {	
					f.Load(in.GetLine());
					if (IsNull(f.GetDouble_nothrow(0)))
						throw Exc(in.Str() + "\n"  + t_("Frecuency not found. May be radiation file doesn't match with Nemoh.cal file"));		
					int col = 1;
					for (int idof2 = 0; idof2 < 6; ++idof2) {			
						//if (dcase.IsDof(ibody, idof2)) {
							dt.A[ibody*6 + idof][ibody*6 + idof2][ifr] = f.GetDouble(col++);
		        			dt.B[ibody*6 + idof][ibody*6 + idof2][ifr] = f.GetDouble(col++);
			        	//}
					}
				}
		    	if (idof < 6)
		    		in.GetLine();
	    	//}
		}
	}
	return true;
}

bool Nemoh::Load_Excitation(String folder) {	
	return Load_Forces(dt.ex, folder, "ExcitationForce.tec");
}

bool Nemoh::Load_Diffraction(String folder) {
	return Load_Forces(dt.sc, folder, "DiffractionForce.tec");
}

bool Nemoh::Load_FroudeKrylov(String folder) {
	return Load_Forces(dt.fk, folder, "FKForce.tec");
}

bool Nemoh::Load_Forces(Hydro::Forces &fc, String nfolder, String fileName) {
	FileInLine in(AFX(nfolder, "Results", fileName));
	if (!in.IsOpen())
		return false;
	String line;
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
	in.GetLine();
	UVector<UVector<int>> dof;
	dof.SetCount(dt.Nb);
	while(!in.IsEof()) {
		line = in.GetLine();
		if (line.StartsWith("Zone") || line.StartsWith("angle"))
	        break;
	}
	Initialize_Forces(fc);
	for (int ih = 0; ih < dt.Nh; ++ih) {
		int ifr = 0;
		while(!in.IsEof()) {
			line = in.GetLine();
			if (line.StartsWith("Zone") || line.StartsWith("angle"))
				break;
			f.Load(line);
			if (IsNull(f.GetDouble_nothrow(0)))
				throw Exc(in.Str() + "\n"  + t_("Frecuency not found. May be excitation file doesn't match with Nemoh.cal file"));		
			int il = 0;
			for (int ib = 0; ib < dt.Nb; ++ib) {
				for (int ibdof = 0; ibdof < 6; ++ibdof) {
					//if (dcase.IsDof(ib, ibdof)) {
						if (ifr >= dt.Nf)
							throw Exc(in.Str() + "\n"  + t_("Number of frequencies higher than the defined in Nemoh.cal file"));		
						double ma = f.GetDouble(1 + 2*il);	
						double ph = -f.GetDouble(1 + 2*il + 1); //-Phase to follow Wamit
						fc[ib][ih](ifr, ibdof) = std::polar(ma, ph); 
						il++; 
					//}
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
	dt.Ainf.setConstant(dt.Nb*6, dt.Nb*6, 0);
	dt.Kirf.SetCount(dt.Nb*6); 	// Initialize Kirf		
    for (int i = 0; i < dt.Nb*6; ++i) {
    	dt.Kirf[i].SetCount(dt.Nb*6); 			 
   		for (int j = 0; j < dt.Nb*6; ++j)
			dt.Kirf[i][j].setConstant(dt.Tirf.size(), NaNDouble);
    }
    while(!in.IsEof()) {
		line = in.GetLine();	
		if (line.Find("Zone t=") >= 0) 
			break;
	}
	for (int ib = 0; ib < dt.Nb; ++ib) {
		for (int ibdof = 0; ibdof < 6; ++ibdof) {
			for (int iNt = 0; iNt < dt.Tirf.size(); ++iNt) {
				f.Load(in.GetLine());
				int col = 1;
				for (int idf = 0; idf < 6; ++idf) {
					dt.Ainf(ibdof, idf) = f.GetDouble(col++);
					dt.Kirf[ib*6+ibdof][ib*6+idf][iNt] = f.GetDouble(col++);
				}
			}
			in.GetLine();
		}
	}
	return true;
}

void Nemoh::Save(String ) {
	throw Exc("Option not implemented");
}		


