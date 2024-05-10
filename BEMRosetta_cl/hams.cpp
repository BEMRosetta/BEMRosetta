// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#include <STEM4U/Utility.h>


String Hams::Load(String file, Function <bool(String, int)> Status) {
	hd().file = file;	
	hd().name = GetFileTitle(file);
	
	//hd().g = g;
	
	hd().x_w = hd().y_w = 0;
	
	String baseFolder = GetFileFolder(GetFileFolder(file));
	
	try {
		if (GetFileExt(file) != ".in") 
			return Format(t_("File '%s' is not of HAMS type"), file);
			
		BEM::Print("\n\n" + Format(t_("Loading '%s'"), file));

		hd().solver = Hydro::HAMS_WAMIT;
		
		String wamitFile = AFX(baseFolder, "Output", "Wamit_format", "Buoy.1");
		
		String ret = Wamit::Load(wamitFile, Status);
		if (!ret.IsEmpty()) 
			return ret;
		
		if (IsNull(hd().Nb))
			return t_("No body found");

	} catch (Exc e) {
		BEM::PrintError("\nError: " + e);
		//hd().lastError = e;
		return e;
	}
	try {		
		double rhog = hd().g_rho_dim();
		String settingsFile = AFX(baseFolder, "Settings.ctrl");
		if (FileExists(settingsFile) && !Load_Settings(settingsFile))
			throw Exc("\n" + Format(t_("Problem loading Settings.ctrl file '%s'"), settingsFile));
		
		String controlFile = AFX(baseFolder, "Input", "ControlFile.in");
		if (FileExists(controlFile)) 
			Load_In(controlFile);
			
		String hydrostaticFile = AFX(baseFolder, /*"Input", */"Hydrostatic.in");
		if (FileExists(hydrostaticFile)) {
			bool iszero = true;
			if (hd().IsLoadedC()) {
				for (int i = 0; i < 36; ++i) {
					if (abs(hd().msh[0].C.array()(i)) > 1E-8) {
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
	
	return String();
}

bool Hams::Load_Settings(String fileName) {
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
	
	//hd().c0.setConstant(3, hd().Nb, Null);
	f.Load(in.GetLine());	
	hd().msh[0].c0.x = f.GetDouble(0);
	hd().msh[0].c0.y = f.GetDouble(1);
	hd().msh[0].c0.z = f.GetDouble(2);
	
	//hd().cg.setConstant(3, hd().Nb, Null);
	f.Load(in.GetLine());	
	hd().msh[0].cg.x = f.GetDouble(0);
	hd().msh[0].cg.y = f.GetDouble(1);
	hd().msh[0].cg.z = f.GetDouble(2);	
	
	return true;
}

// Load Hydrostatic.in obtained from WAMIT_MeshTran.exe
// Warning: This file is not normally the one at Input folder
bool Hams::Load_HydrostaticMesh(String fileName, double rhog) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
 
	hd().Nb = 1;
	hd().msh.SetCount(1);
	
	//hd().c0.setConstant(3, hd().Nb, Null);
	
	//hd().C.SetCount(hd().Nb);
	for (int ib = 0; ib < hd().Nb; ++ib)
		hd().msh[ib].C.setConstant(6, 6, 0);

	in.GetLine();
	f.GetLine();	
	hd().msh[0].c0.x = f.GetDouble(0);
	hd().msh[0].c0.y = f.GetDouble(1);
	hd().msh[0].c0.z = f.GetDouble(2);
		
	in.GetLine(8);
	for (int r = 0; r < 6; ++r) {
		f.GetLine();	
		for (int c = 0; c < 6; ++c) 
			hd().msh[0].C(r, c) = f.GetDouble(c)/rhog;
	}
	return true;
}
	

UVector<String> Hams::Check() const {
	UVector<String> ret;
	
	if (hd().msh.size() != 1)
		ret << t_("HAMS just allows one body");

	return ret;
}

bool Hams::Load_In(String fileName) {
	fileName = AFX(GetFileFolder(fileName), "ControlFile.in");
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	hd().solver = Hydro::HAMS_WAMIT;
	
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
			hd().h = f.GetDouble(1);
		else if (var == "Input_frequency_type") {
			input_frequency_type = f.GetInt(1); 
			if (input_frequency_type != 3 && input_frequency_type != 4)
				throw Exc(t_("HAMS loader just allows loading input_frequency_type = 3 wave frequency or 4 wave period"));
		} else if (var == "Output_frequency_type") {
			output_frequency_type = f.GetInt(1);
			if (output_frequency_type != 3 && output_frequency_type != 4)
				throw Exc(t_("HAMS loader just allows loading output_frequency_type = 3 wave frequency or 4 wave period"));
		} else if (var == "Number_of_frequencies") {
			int Nf = f.GetInt(1);
			if (!IsNull(hd().Nf)) {
				if (hd().Nf != abs(Nf))
					throw Exc(Format(t_("Number of frequencies %d does not match with previously loaded %d"), Nf, hd().Nf));
			} else 
				hd().Nf = Nf;
			UVector<double> data;
			if (Nf < 0) {
				Nf = -Nf;
				f.GetLine();
				if (f.size() < 2 || Trim(f.GetText(0)) != "Minimum_frequency_Wmin")
					throw Exc(t_("Minimum_frequency_Wmin not found"));
				double minF = f.GetDouble(1);
				f.GetLine();
				if (f.size() < 2 || Trim(f.GetText(0)) != "Frequency_step")
					throw Exc(t_("Frequency_step not found"));
				double maxF = minF + f.GetDouble(1)*Nf;
				
				LinSpaced(data, Nf, minF, maxF);
			} else {
				f.GetLine();
				if (f.size() < 2)
					throw Exc(t_("Frequencies or periods not found"));

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
			}
			if (hd().w.IsEmpty()) {
				hd().w = pick(data);
				/*for (int ifr = 0; ifr < hd().Nf; ++ifr) 
					hd().T[ifr] = 2*M_PI/hd().w[ifr];*/
			} else if (!CompareDecimals(hd().w, data, 3))
				throw Exc(t_("List of frequencies does not match with previously loaded"));

		} else if (var == "Number_of_headings") {
			int Nh = f.GetInt(1);
			if (!IsNull(hd().Nh)) {
				if (hd().Nh != abs(Nh))
					throw Exc(Format(t_("Number of headings %d does not match with previously loaded %d"), Nh, hd().Nh));
			} else 
				hd().Nh = Nh;
			UVector<double> data;
			if (Nh < 0) {
				Nh = -Nh;
				f.GetLine();
				if (f.size() < 2 || Trim(f.GetText(0)) != "Minimum_heading")
					throw Exc(t_("Minimum_heading not found"));
				double minH = f.GetDouble(1);
				f.GetLine();
				if (f.size() < 2 || Trim(f.GetText(0)) != "Heading_step")
					throw Exc(t_("Heading_step not found"));
				double maxH = minH + f.GetDouble(1)*Nh;
				
				LinSpaced(data, Nh, minH, maxH);
			} else {
				f.GetLine();
				if (f.size() < 1)
					throw Exc(t_("Headings not found"));
				
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
			}
			if (hd().head.IsEmpty()) 
				hd().head = pick(data);
			else if (!CompareDecimals(hd().head, data, 2))
				throw Exc(t_("List of headings does not match with previously loaded"));

		} else if (var == "Reference_body_centre") {
			if (f.size() < 4)
				throw Exc(t_("Lack of data in Reference_body_center"));
			hd().msh.SetCount(1);
			Mesh &body = hd().msh[0];
			body.c0[0] = f.GetDouble(1);
			body.c0[1] = f.GetDouble(2);
			body.c0[2] = f.GetDouble(3);
			String meshFile = AFX(GetFileFolder(fileName), "HullMesh.pnl");
			if (FileExists(meshFile))
				body.fileName = meshFile;
			String lidFile = AFX(GetFileFolder(fileName), "WaterplaneMesh.pnl");
			if (FileExists(lidFile))
				body.lidFile = lidFile;
		}
	}
	return LoadHydrostatic(AFX(GetFileFolder(fileName), "Hydrostatic.in"));
}

bool Hams::LoadHydrostatic(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	hd().msh.SetCount(1);
	Mesh &body = hd().msh[0];
	
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
			InMatrix(f, body.Cadd);	
	}
	return true;
}

void Hams::SaveFolder(String folderBase, bool bin, int numCases, int numThreads, int) const {
	SaveFolder0(folderBase, bin, 1, true, numThreads);
	if (numCases > 1)
		SaveFolder0(folderBase, bin, numCases, false, numThreads);
}

void Hams::SaveFolder0(String folderBase, bool bin, int numCases, bool deleteFolder, int numThreads) const {
	hd().BeforeSaveCase(folderBase, numCases, deleteFolder);
	
	#define MIN_F_HAMS 0.01
	
	double fixminF = max(MIN_F_HAMS, First(hd().w));
	
	UVector<int> valsf;
	int _nf;
	double _minf, _maxf;
	int ifr = 0;
	UVector<double> freqs;
	if (numCases > 1)  
		valsf = NumSets(hd().Nf, numCases);
	
	String solvName = "HAMS_x64.exe";
	if (bin) {
		String source = Bem().hamsPath;
		solvName = GetFileName(source);
		String destNew = AFX(folderBase, solvName);
		if (!FileCopy(source, destNew)) 
			throw Exc(Format(t_("Problem copying Hams exe file from '%s'"), Bem().hamsPath));
		source = AFX(GetFileFolder(Bem().hamsPath), "libiomp5md.dll");		
		destNew = AFX(folderBase, "libiomp5md.dll");		
		if (!FileCopy(source, destNew)) 
			throw Exc(Format(t_("Problem copying Hams dll file from '%s'"), source));					
	} 
	String meshName = "WAMIT_MeshTran.exe";
	if (bin) {
		String source = Bem().hamsMeshPath;
		meshName = GetFileName(source);
		String destNew = AFX(folderBase, meshName);
		if (!FileCopy(source, destNew)) 
			throw Exc(Format(t_("Problem copying Hams mesh exe file from '%s'"), Bem().hamsMeshPath));				
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
			_maxf = hd().w[ifr + deltaf - 1];
			_nf = deltaf;
			ifr += deltaf;
		} else {
			folder = folderBase;
			_nf = hd().Nf;
			_minf = fixminF;
			_maxf = Last(hd().w);
		}
		String folderInput = AFX(folder, "Input");
		if (!DirectoryCreateX(folderInput))
			throw Exc(Format(t_("Problem creating '%s' folder"), folderInput));
		
		Save_ControlFile(folderInput, _nf, _minf, _maxf, numThreads);
		Save_Hydrostatic(folderInput);
	
		bool y0zmesh = false, x0zmesh = false;
		UArray<Mesh> msh;
		int ib = 0;		// Just one file
		
		String err = Mesh::Load(msh, hd().msh[ib].fileName, hd().rho, hd().g, false, y0zmesh, x0zmesh);
		if (!err.IsEmpty())
			throw Exc(err);
		
		Mesh &mesh = First(msh);
		
		if (y0zmesh == true && x0zmesh == true) 
			y0zmesh = false;

		String dest = AFX(folderInput, "HullMesh.pnl");
		Mesh::SaveAs(mesh, dest, Mesh::HAMS_PNL, Mesh::UNDERWATER, hd().rho, hd().g, y0zmesh, x0zmesh);
		
		bool y0zlid = false, x0zlid = false;	// Hull symmetries rules over lid ones
		if (!hd().msh[ib].lidFile.IsEmpty()) {
			String err = Mesh::Load(msh, hd().msh[ib].lidFile, hd().rho, hd().g, false, y0zlid, x0zlid);
			if (!err.IsEmpty())
				throw Exc(err);
			
			Mesh &mesh = First(msh);
				
			String dest = AFX(folderInput, "WaterplaneMesh.pnl");
			Mesh::SaveAs(mesh, dest, Mesh::HAMS_PNL, Mesh::ALL, hd().rho, hd().g, y0zmesh, x0zmesh);
		}
		
		Save_Settings(folder, !hd().msh[ib].lidFile.IsEmpty());
		
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

void Hams::Save_Bat(String folder, String batname, String caseFolder, bool bin, String solvName, String meshName) const {
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

void Hams::OutMatrix(FileOut &out, String header, const Eigen::MatrixXd &mat) {
	out << "\n " << header << ":";
	for (int y = 0; y < 6; ++y) {
		out << "\n";
		for (int x = 0; x < 6; ++x)
			out << "   " << FDS(mat(x, y), 11, true);
	}
}

void Hams::InMatrix(LineParser &f, Eigen::MatrixXd &mat) {
	for (int y = 0; y < 6; ++y) {
		f.GetLine();
		for (int x = 0; x < 6; ++x)
			mat(x, y) = f.GetDouble(x);
	}
}
	
void Hams::Save_Hydrostatic(String folderInput) const {
	String fileName = AFX(folderInput, "Hydrostatic.in");
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to create '%s'"), fileName));
	out << " Centre of Gravity:";
	
	if (hd().msh.IsEmpty())
		throw Exc(t_("No bodies found"));
	
	const Mesh &b = hd().msh[0];
	out << Format("\n%s %s %s", FDS(b.cg[0], 15, true), 
								FDS(b.cg[1], 15, true), 
								FDS(b.cg[2], 15, true));

	OutMatrix(out, "Body Mass Matrix", b.M);
	OutMatrix(out, "External Linear Damping Matrix", b.Dlin);
	OutMatrix(out, "External Quadratic Damping Matrix", b.Dquad);
	OutMatrix(out, "Hydrostatic Restoring Matrix", b.C);
	OutMatrix(out, "External Restoring Matrix", b.Cadd);
}


void Hams::Save_Settings(String folderInput, bool thereIsLid) const {
	String fileName = AFX(folderInput, "Settings.ctrl");
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to create '%s'"), fileName));
	
	Mesh mesh;
	String res = Mesh::Load(mesh, AFX(folderInput, "Input", "HullMesh.pnl"), hd().rho, hd().g, Null, Null, false);
	if (!res.IsEmpty())
		throw Exc(res);
	
	if (thereIsLid) {
		Mesh lid;
		lid.mesh.AddWaterSurface(mesh.mesh, mesh.under, 'f', Bem().roundVal, Bem().roundEps); 
		lid.AfterLoad(hd().rho, hd().g, false, false);
		
		mesh.Append(lid.mesh, hd().rho, hd().g);
	}
	Mesh::SaveAs(mesh, AFX(folderInput, "Input", "mesh.gdf"), Mesh::WAMIT_GDF, Mesh::ALL, hd().rho, hd().g, false, false);	
	
	out << hd().g << "\n";
	out << hd().rho << "\n";
	out << ".\\Input\\mesh.gdf" << "\n";
	out << hd().msh[0].c0[0] << "   " << hd().msh[0].c0[1] << "   " << hd().msh[0].c0[2] << "\n";
	out << hd().msh[0].cg[0] << "   " << hd().msh[0].cg[1] << "   " << hd().msh[0].cg[2] << "\n";
	out << "\n"
		<< "The format is as below:\n"
		<< "gravity acceleration\n"
		<< "sea water density\n"
		<< "mesh file name (gdf format)\n"
		<< "Centre of rotation\n"
		<< "Centre of gravity";
}

void Hams::Save_ControlFile(String folderInput, int _nf, double _minf, double _maxf,
					int numThreads) const {
	String fileName = AFX(folderInput, "ControlFile.in");
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to create '%s'"), fileName));	
	
	out << "   --------------HAMS Control file---------------"
		   "\n";
	double depth = hd().h;
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
	out << "\n    Number_of_headings         " << hd().Nh;
	UVector<double> headings;
	LinSpaced(headings, hd().Nh, First(hd().head), Last(hd().head));
	out << "\n    ";
	for (const auto &heading : headings)
		out << Format("%.4f ", heading);	
	out << "\n   #End Definition of Wave Headings"
		   "\n";
	out << "\n    Reference_body_centre " << Format("%11.3f %11.3f %11.3f", 
								hd().msh[0].c0[0], hd().msh[0].c0[1], hd().msh[0].c0[2]) << 
		   "\n    Reference_body_length   1.D0"
		   "\n    Wave_diffrac_solution    2";
	bool remove_irr_freq = !hd().msh[0].lidFile.IsEmpty();
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