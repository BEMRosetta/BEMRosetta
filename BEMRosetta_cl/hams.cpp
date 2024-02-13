// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#include <STEM4U/Utility.h>


bool HAMS::Load(String file, Function <bool(String, int)> Status) {
	hd().file = file;	
	hd().name = GetFileTitle(file);
	
	//hd().g = g;
	
	hd().x_w = hd().y_w = 0;
	
	String baseFolder = GetFileFolder(GetFileFolder(file));
	
	try {
		if (GetFileExt(file) != ".in") 
			throw Exc("\n" + Format(t_("File '%s' is not of HAMS type"), file));
			
		BEM::Print("\n\n" + Format(t_("Loading '%s'"), file));

		hd().code = Hydro::HAMS;
		
		String wamitFile = AFX(baseFolder, "Output", "Wamit_format", "Buoy.1");
		
		if (!Wamit::Load(wamitFile, Status)) 
			return false;
		
		if (IsNull(hd().Nb))
			return false;

	} catch (Exc e) {
		BEM::PrintError("\nError: " + e);
		hd().lastError = e;
		return false;
	}
	try {		
		double rhog = hd().g_rho_dim();
		String settingsFile = AFX(baseFolder, "Settings.ctrl");
		if (FileExists(settingsFile) && !Load_Settings(settingsFile))
			throw Exc("\n" + Format(t_("Problem loading Settings.ctrl file '%s'"), settingsFile));
		
		String controlFile = AFX(baseFolder, "Input", "ControlFile.in");
		if (FileExists(controlFile)) {
			HamsCase cas;
			cas.Load(controlFile);
			hd().h = cas.h;
			if (cas.bodies.size() > 0) {
				hd().c0.setConstant(3, hd().Nb, Null);
				hd().c0(0, 0) = cas.bodies[0].c0(0);
				hd().c0(1, 0) = cas.bodies[0].c0(1);
				hd().c0(2, 0) = cas.bodies[0].c0(2);
				hd().cg.setConstant(3, hd().Nb, Null);
				hd().cg(0, 0) = cas.bodies[0].cg(0);
				hd().cg(1, 0) = cas.bodies[0].cg(1);
				hd().cg(2, 0) = cas.bodies[0].cg(2);
			}
		}
			
		String hydrostaticFile = AFX(baseFolder, /*"Input", */"Hydrostatic.in");
		if (FileExists(hydrostaticFile)) {
			bool iszero = true;
			if (hd().IsLoadedC()) {
				for (int i = 0; i < hd().C[0].size(); ++i) {
					if (abs(hd().C[0].array()(i)) > 1E-8) {
						iszero = false;
						break;
					}
				}
			}
			if (iszero && !Load_HydrostaticMesh(hydrostaticFile, rhog))
				throw Exc("\n" + Format(t_("Problem loading Hydrostatic file '%s'"), hydrostaticFile));
		}
	} catch (Exc e) {
		BEM::PrintWarning("\nWarning: " + e);
	}
	
	return true;
}

bool HAMS::Load_Settings(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
	
	f.Load(in.GetLine());
	hd().g = f.GetDouble(0);
	f.Load(in.GetLine());
	hd().rho = f.GetDouble(0);
	
	in.GetLine();
	
	hd().c0.setConstant(3, hd().Nb, Null);
	f.Load(in.GetLine());	
	hd().c0(0, 0) = f.GetDouble(0);
	hd().c0(1, 0) = f.GetDouble(1);
	hd().c0(2, 0) = f.GetDouble(2);
	
	hd().cg.setConstant(3, hd().Nb, Null);
	f.Load(in.GetLine());	
	hd().cg(0, 0) = f.GetDouble(0);
	hd().cg(1, 0) = f.GetDouble(1);
	hd().cg(2, 0) = f.GetDouble(2);	
	
	return true;
}

// Load Hydrostatic.in obtained from WAMIT_MeshTran.exe
// Warning: This file is not normally the one at Input folder
bool HAMS::Load_HydrostaticMesh(String fileName, double rhog) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
 
	hd().Nb = 1;
	
	hd().c0.setConstant(3, hd().Nb, Null);
	
	hd().C.SetCount(hd().Nb);
	for (int ib = 0; ib < hd().Nb; ++ib)
		hd().C[ib].setConstant(6, 6, 0);

	in.GetLine();
	f.GetLine();	
	hd().c0(0, 0) = f.GetDouble(0);
	hd().c0(1, 0) = f.GetDouble(1);
	hd().c0(2, 0) = f.GetDouble(2);
		
	in.GetLine(8);
	for (int r = 0; r < 6; ++r) {
		f.GetLine();	
		for (int c = 0; c < 6; ++c) 
			hd().C[0](r, c) = f.GetDouble(c)/rhog;
	}
	return true;
}
	

UVector<String> HamsCase::Check() const {
	UVector<String> ret;
	
	if (bodies.size() != 1)
		ret << t_("HAMS just allows one body");

	return ret;
}

bool HamsCase::Load(String fileName) {
	fileName = AFX(GetFileFolder(fileName), "ControlFile.in");
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	solver = HAMS;
	
	String line;
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
	
	int input_frequency_type = 0, output_frequency_type = 0;
	
	while (!f.IsEof()) {
		f.GetLine();
		
		if (f.GetText().Find("End HAMS Control file") >= 0)
			break;
		if (f.size() < 2)
			continue;
		String var = Trim(f.GetText(0));
		if (var.StartsWith("#"))
			continue;
		if (var == "Waterdepth")
			h = f.GetDouble(1);
		else if (var == "Input_frequency_type") {
			input_frequency_type = f.GetInt(1); 
			if (input_frequency_type != 3 && input_frequency_type != 4)
				throw Exc(t_("HAMS loader just allows loading input_frequency_type = 3 wave frequency or 4 wave period"));
		} else if (var == "Output_frequency_type") {
			output_frequency_type = f.GetInt(1);
			if (output_frequency_type != 3 && output_frequency_type != 4)
				throw Exc(t_("HAMS loader just allows loading output_frequency_type = 3 wave frequency or 4 wave period"));
		} else if (var == "Number_of_frequencies") {
			Nf = f.GetInt(1);
			if (Nf < 0) {
				Nf = -Nf;
				f.GetLine();
				if (f.size() < 2 || Trim(f.GetText(0)) != "Minimum_frequency_Wmin")
					throw Exc(t_("Minimum_frequency_Wmin not found"));
				minF = f.GetDouble(1);
				f.GetLine();
				if (f.size() < 2 || Trim(f.GetText(0)) != "Frequency_step")
					throw Exc(t_("Frequency_step not found"));
				maxF = minF + f.GetDouble(1)*Nf;
			} else {
				f.GetLine();
				if (f.size() < 2)
					throw Exc(t_("Frequencies or periods not found"));
				UVector<double> data;
				double delta;
				for (int i = 0; i < f.size(); ++i) {
					if (input_frequency_type == 3)
						data << f.GetDouble(i);
					else
						data.At(0, 2*M_PI/f.GetDouble(i));
					if (i == 1)
						delta = data[i] - data[i-1];
					else if (i > 1 && !EqualDecimals(delta, data[i] - data[i-1], 2)) {
						BEM::PrintWarning(t_("HAMS loader just allows equidistant frequencies"));
						continue;
					}
				}
				minF = data[0];
				maxF = data[data.size()-1];
			}
		} else if (var == "Number_of_headings") {
			Nh = f.GetInt(1);
			if (Nh < 0) {
				Nh = -Nh;
				f.GetLine();
				if (f.size() < 2 || Trim(f.GetText(0)) != "Minimum_heading")
					throw Exc(t_("Minimum_heading not found"));
				minH = f.GetDouble(1);
				f.GetLine();
				if (f.size() < 2 || Trim(f.GetText(0)) != "Heading_step")
					throw Exc(t_("Heading_step not found"));
				maxH = minH + f.GetDouble(1)*Nh;
			} else {
				f.GetLine();
				if (f.size() < 1)
					throw Exc(t_("Headings not found"));
				UVector<double> data;
				double delta = 0;
				for (int i = 0; i < f.size(); ++i) {
					data << f.GetDouble(i);
					if (i == 1)
						delta = data[i] - data[i-1];
					else if (i > 1 && !EqualDecimals(delta, data[i] - data[i-1], 2)) {
						BEM::PrintWarning(t_("HAMS loader just allows equidistant headings"));
						continue;
					}
				}
				minH = data[0];
				maxH = data[data.size()-1];
			}
		} else if (var == "Reference_body_centre") {
			if (f.size() < 4)
				throw Exc(t_("Lack of data in Reference_body_center"));
			bodies.SetCount(1);
			BEMBody &body = bodies[0];
			body.c0[0] = f.GetDouble(1);
			body.c0[1] = f.GetDouble(2);
			body.c0[2] = f.GetDouble(3);
			String meshFile = AFX(GetFileFolder(fileName), "HullMesh.pnl");
			if (FileExists(meshFile))
				body.meshFile = meshFile;
			String lidFile = AFX(GetFileFolder(fileName), "WaterplaneMesh.pnl");
			if (FileExists(lidFile))
				body.lidFile = lidFile;
		}
	}
	return LoadHydrostatic(AFX(GetFileFolder(fileName), "Hydrostatic.in"));
}

bool HamsCase::LoadHydrostatic(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	bodies.SetCount(1);
	BEMBody &body = bodies[0];
	
	String line;
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
	
	while (!f.IsEof()) {
		f.GetLine();
		
		if (f.size() < 1)
			continue;	
		
		String line = Trim(f.GetText());
	
		if (line == "Centre of Gravity:") {
			f.GetLine();
			if (f.size() < 3)
				throw Exc(t_("Centre of Gravity data is not complete"));
			body.cg[0] = f.GetDouble(0);
			body.cg[1] = f.GetDouble(1);
			body.cg[2] = f.GetDouble(2);
		} else if (line == "Body Mass Matrix:") 
			InMatrix(f, body.M);
		else if (line == "External Linear Damping Matrix:") 
			InMatrix(f, body.Dlin);
		else if (line == "External Quadratic Damping Matrix:") 
			InMatrix(f, body.Dquad);	
		else if (line == "Hydrostatic Restoring Matrix:") 
			InMatrix(f, body.C);	
		else if (line == "External Restoring Matrix:") 
			InMatrix(f, body.Cext);	
	}
	return true;
}

void HamsCase::SaveFolder(String folderBase, bool bin, int numCases, int numThreads, const BEM &bem, int) const {
	SaveFolder0(folderBase, bin, 1, bem, true, numThreads);
	if (numCases > 1)
		SaveFolder0(folderBase, bin, numCases, bem, false, numThreads);
}

void HamsCase::SaveFolder0(String folderBase, bool bin, int numCases, const BEM &bem, bool deleteFolder, int numThreads) const {
	BeforeSave(folderBase, numCases, deleteFolder);
	
	#define MIN_F_HAMS 0.01
	
	double fixminF = minF;
	if (fixminF < MIN_F_HAMS)
		fixminF = MIN_F_HAMS;
	
	UVector<int> valsf;
	int _nf;
	double _minf, _maxf;
	int ifr = 0;
	UVector<double> freqs;
	if (numCases > 1) { 
		LinSpaced(freqs, Nf, fixminF, maxF);
		valsf = NumSets(Nf, numCases);
	}
	
	String solvName = "HAMS_x64.exe";
	if (bin) {
		String source = bem.hamsPath;
		solvName = GetFileName(source);
		String destNew = AFX(folderBase, solvName);
		if (!FileCopy(source, destNew)) 
			throw Exc(Format(t_("Problem copying Hams exe file from '%s'"), bem.hamsPath));
		source = AFX(GetFileFolder(bem.hamsPath), "libiomp5md.dll");		
		destNew = AFX(folderBase, "libiomp5md.dll");		
		if (!FileCopy(source, destNew)) 
			throw Exc(Format(t_("Problem copying Hams dll file from '%s'"), source));					
	} 
	String meshName = "WAMIT_MeshTran.exe";
	if (bin) {
		String source = bem.hamsMeshPath;
		meshName = GetFileName(source);
		String destNew = AFX(folderBase, meshName);
		if (!FileCopy(source, destNew)) 
			throw Exc(Format(t_("Problem copying Hams mesh exe file from '%s'"), bem.hamsMeshPath));				
	} 
		
	//String sumcases;
	for (int i = 0; i < numCases; ++i) {
		String folder;
		if (numCases > 1) {
			folder = AFX(folderBase, Format("HAMS_Part_%d", i+1));
			if (!DirectoryCreateX(folder))
				throw Exc(Format(t_("Problem creating '%s' folder"), folder));
			//sumcases << " " << AFX(folder, "Nemoh.cal");
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
		String folderInput = AFX(folder, "Input");
		if (!DirectoryCreateX(folderInput))
			throw Exc(Format(t_("Problem creating '%s' folder"), folderInput));
		
		Save_ControlFile(folderInput, _nf, _minf, _maxf, numThreads);
		Save_Hydrostatic(folderInput);
	
		bool y0zmesh = false, x0zmesh = false;
		UArray<Mesh> msh;
		int ib = 0;		// Just one file
		
		String err = Mesh::Load(msh, bodies[ib].meshFile, rho, g, false, y0zmesh, x0zmesh);
		if (!err.IsEmpty())
			throw Exc(err);
		
		Mesh &mesh = First(msh);
		
		if (y0zmesh == true && x0zmesh == true) 
			y0zmesh = false;

		String dest = AFX(folderInput, "HullMesh.pnl");
		Mesh::SaveAs(mesh, dest, Mesh::HAMS_PNL, Mesh::UNDERWATER, rho, g, y0zmesh, x0zmesh);
		
		bool y0zlid = false, x0zlid = false;	// Hull symmetries rules over lid ones
		if (!bodies[ib].lidFile.IsEmpty()) {
			String err = Mesh::Load(msh, bodies[ib].lidFile, rho, g, false, y0zlid, x0zlid);
			if (!err.IsEmpty())
				throw Exc(err);
			
			Mesh &mesh = First(msh);
				
			String dest = AFX(folderInput, "WaterplaneMesh.pnl");
			Mesh::SaveAs(mesh, dest, Mesh::HAMS_PNL, Mesh::ALL, rho, g, y0zmesh, x0zmesh);
		}
		
		Save_Settings(folder, !bodies[ib].lidFile.IsEmpty(), bem);
		
		String folderOutput = AFX(folder, "Output");
		if (!DirectoryCreateX(folderOutput))
			throw Exc(Format(t_("Problem creating '%s' folder"), folderOutput));
		if (!DirectoryCreateX(AFX(folderOutput, "Hams_format")))
			throw Exc(Format(t_("Problem creating '%s' folder"), AFX(folderOutput, "Hams_format")));
		if (!DirectoryCreateX(AFX(folderOutput, "Hydrostar_format")))
			throw Exc(Format(t_("Problem creating '%s' folder"), AFX(folderOutput, "Hydrostar_format")));
		if (!DirectoryCreateX(AFX(folderOutput, "Wamit_format")))
			throw Exc(Format(t_("Problem creating '%s' folder"), AFX(folderOutput, "Wamit_format")));
		
		if (numCases > 1) 
			Save_Bat(folderBase, Format("HAMS_Part_%d.bat", i+1), Format("HAMS_Part_%d", i+1), bin, AFX("..", solvName), AFX("..", meshName));
		else
			Save_Bat(folder, "HAMS.bat", Null, bin, solvName, meshName);
	}
}

void HamsCase::Save_Bat(String folder, String batname, String caseFolder, bool bin, String solvName, String meshName) const {
	String fileName = AFX(folder, batname);
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to create '%s'"), fileName));
	
	out << Format("title %s in '%s'\n", solvName, caseFolder);
	
	if (!IsNull(caseFolder))
		out << "cd \"" << caseFolder << "\"\n";
	
	out << "\"" << meshName << "\"\n";
	out << "\"" << solvName << "\"\n";
}

void HamsCase::OutMatrix(FileOut &out, String header, const Eigen::MatrixXd &mat) {
	out << "\n " << header << ":";
	for (int y = 0; y < 6; ++y) {
		out << "\n";
		for (int x = 0; x < 6; ++x)
			out << "   " << FDS(mat(x, y), 11, true);
	}
}

void HamsCase::InMatrix(LineParser &f, Eigen::MatrixXd &mat) {
	for (int y = 0; y < 6; ++y) {
		f.GetLine();
		for (int x = 0; x < 6; ++x)
			mat(x, y) = f.GetDouble(x);
	}
}
	
void HamsCase::Save_Hydrostatic(String folderInput) const {
	String fileName = AFX(folderInput, "Hydrostatic.in");
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to create '%s'"), fileName));
	out << " Centre of Gravity:";
	
	if (bodies.IsEmpty())
		throw Exc(t_("No bodies found"));
	
	const BEMBody &b = bodies[0];
	out << Format("\n%s %s %s", FDS(b.cg[0], 15, true), 
								FDS(b.cg[1], 15, true), 
								FDS(b.cg[2], 15, true));

	OutMatrix(out, "Body Mass Matrix", b.M);
	OutMatrix(out, "External Linear Damping Matrix", b.Dlin);
	OutMatrix(out, "External Quadratic Damping Matrix", b.Dquad);
	OutMatrix(out, "Hydrostatic Restoring Matrix", b.C);
	OutMatrix(out, "External Restoring Matrix", b.Cext);
}


void HamsCase::Save_Settings(String folderInput, bool thereIsLid, const BEM &bem) const {
	String fileName = AFX(folderInput, "Settings.ctrl");
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to create '%s'"), fileName));
	
	Mesh mesh;
	String res = Mesh::Load(mesh, AFX(folderInput, "Input", "HullMesh.pnl"), rho, g, Null, Null, false);
	if (!res.IsEmpty())
		throw Exc(res);
	
	if (thereIsLid) {
		Mesh lid;
		lid.mesh.AddWaterSurface(mesh.mesh, mesh.under, 'f', bem.roundVal, bem.roundEps); 
		lid.AfterLoad(rho, g, false, false);
		
		mesh.Append(lid.mesh, rho, g);
	}
	Mesh::SaveAs(mesh, AFX(folderInput, "Input", "mesh.gdf"), Mesh::WAMIT_GDF, Mesh::ALL, rho, g, false, false);	
	
	out << g << "\n";
	out << rho << "\n";
	out << ".\\Input\\mesh.gdf" << "\n";
	out << bodies[0].c0[0] << "   " << bodies[0].c0[1] << "   " << bodies[0].c0[2] << "\n";
	out << bodies[0].cg[0] << "   " << bodies[0].cg[1] << "   " << bodies[0].cg[2] << "\n";
	out << "\n"
		<< "The format is as below:\n"
		<< "gravity acceleration\n"
		<< "sea water density\n"
		<< "mesh file name (gdf format)\n"
		<< "Centre of rotation\n"
		<< "Centre of gravity";
}

void HamsCase::Save_ControlFile(String folderInput, int _nf, double _minf, double _maxf,
					int numThreads) const {
	String fileName = AFX(folderInput, "ControlFile.in");
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to create '%s'"), fileName));	
	
	out << "   --------------HAMS Control file---------------"
		   "\n";
	double depth = h;
	if (depth >= 400)
		depth = -1;
	out << "\n   Waterdepth  " << Format("%.4f", depth) << "D0";
	out << "\n"
		   "\n   #Start Definition of Wave Frequencies"
		   "\n    0_inf_frequency_limits      1"
		   "\n    Input_frequency_type        3"
		   "\n    Output_frequency_type       3";
	out << "\n    Number_of_frequencies      " << _nf;
	UVector<double> freqs;
	LinSpaced(freqs, _nf, _minf, _maxf);
	out << "\n    ";
	for (const auto &freq : freqs)
		out << Format("%.4f ", freq);	
	out << "\n   #End Definition of Wave Frequencies"
		   "\n"
		   "\n   #Start Definition of Wave Headings";
	out << "\n    Number_of_headings         " << Nh;
	UVector<double> headings;
	LinSpaced(headings, Nh, minH, maxH);
	out << "\n    ";
	for (const auto &heading : headings)
		out << Format("%.4f ", heading);	
	out << "\n   #End Definition of Wave Headings"
		   "\n";
	out << "\n    Reference_body_centre " << Format("%11.3f %11.3f %11.3f", 
								bodies[0].c0[0], bodies[0].c0[1], bodies[0].c0[2]) << 
		   "\n    Reference_body_length   1.D0"
		   "\n    Wave_diffrac_solution    2";
	bool remove_irr_freq = !bodies[0].lidFile.IsEmpty();
	out << "\n    If_remove_irr_freq      " << (remove_irr_freq ? 1 : 0);
	out << "\n    Number of threads       " << numThreads;
	out << "\n"
		   "\n   #Start Definition of Pressure and/or Elevation (PE)"
		   "\n    Number_of_field_points     1                           # number of field points where to calculate PE"
		   "\n 0.000000    0.000000    0.000000    Global_coords_point_1"
		   "\n   #End Definition of Pressure and/or Elevation"
   		   "\n";
	out << "\n    ----------End HAMS Control file---------------"
		   "\n   Input_frequency_type options:"
		   "\n   1--deepwater wave number; 2--finite-depth wave number; 3--wave frequency; 4--wave period; 5--wave length"
   		   "\n   Output_frequency_type options: same as Input_frequency_type options";
}