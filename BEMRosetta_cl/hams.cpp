// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#include <STEM4U/Utility.h>


String Hams::Load(String file, Function <bool(String, int)> Status) {
	dt.file = file;	
	dt.name = GetFileTitle(file);
	
	dt.x_w = dt.y_w = 0;
	
	String baseFolder = GetFileFolder(GetFileFolder(file));
	
	try {
		if (GetFileExt(file) != ".in") 
			return Format(t_("File '%s' is not of HAMS type"), file);
			
		BEM::Print("\n\n" + Format(t_("Loading '%s'"), file));

		dt.solver = Hydro::HAMS_WAMIT;
		
		String wamitFile = AFX(baseFolder, "Output", "Wamit_format", "Buoy.1");
		
		String ret = Wamit::Load(wamitFile, Status);
		if (!ret.IsEmpty()) 
			return ret;
		
		if (IsNull(dt.Nb))
			return t_("No body found");

	} catch (Exc e) {
		BEM::PrintError("\nError: " + e);
		//dt.lastError = e;
		return e;
	}
	try {		
		double rhog = g_rho_dim();
		String settingsFile = AFX(baseFolder, "Settings.ctrl");
		if (FileExists(settingsFile) && !Load_Settings(settingsFile))
			throw Exc("\n" + Format(t_("Problem loading Settings.ctrl file '%s'"), settingsFile));
		
		String controlFile = AFX(baseFolder, "Input", "ControlFile.in");
		if (FileExists(controlFile)) 
			Load_ControlFile(controlFile);
			
		String hydrostaticFile = AFX(baseFolder, /*"Input", */"Hydrostatic.in");
		if (FileExists(hydrostaticFile)) {
			bool iszero = true;
			if (IsLoadedC()) {
				for (int i = 0; i < 36; ++i) {
					if (abs(dt.msh[0].dt.C.array()(i)) > 1E-8) {
						iszero = false;
						break;
					}
				}
			}
			if (iszero && !Load_HydrostaticBody(hydrostaticFile, rhog))
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
	dt.g = f.GetDouble(0);
	f.Load(in.GetLine());
	dt.rho = f.GetDouble(0);
	
	in.GetLine();
	
	//dt.c0.setConstant(3, dt.Nb, Null);
	f.Load(in.GetLine());	
	dt.msh[0].dt.c0.x = f.GetDouble(0);
	dt.msh[0].dt.c0.y = f.GetDouble(1);
	dt.msh[0].dt.c0.z = f.GetDouble(2);
	
	//dt.cg.setConstant(3, dt.Nb, Null);
	f.Load(in.GetLine());	
	dt.msh[0].dt.cg.x = f.GetDouble(0);
	dt.msh[0].dt.cg.y = f.GetDouble(1);
	dt.msh[0].dt.cg.z = f.GetDouble(2);	
	
	return true;
}

// Load Hydrostatic.in obtained from WAMIT_BodyTran.exe
// Warning: This file is not normally the one at Input folder
bool Hams::Load_HydrostaticBody(String fileName, double rhog) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
 
	dt.Nb = 1;
	dt.msh.SetCount(1);
	
	//dt.c0.setConstant(3, dt.Nb, Null);
	
	//dt.C.SetCount(dt.Nb);
	for (int ib = 0; ib < dt.Nb; ++ib)
		dt.msh[ib].dt.C.setConstant(6, 6, 0);

	in.GetLine();
	f.GetLine();	
	dt.msh[0].dt.c0.x = f.GetDouble(0);
	dt.msh[0].dt.c0.y = f.GetDouble(1);
	dt.msh[0].dt.c0.z = f.GetDouble(2);
		
	in.GetLine(8);
	for (int r = 0; r < 6; ++r) {
		f.GetLine();	
		for (int c = 0; c < 6; ++c) 
			dt.msh[0].dt.C(r, c) = f.GetDouble(c)/rhog;
	}
	return true;
}
	

UVector<String> Hams::Check() const {
	UVector<String> ret;
	
	if (dt.msh.size() != 1)
		ret << t_("HAMS just allows one body");

	return ret;
}

bool Hams::Load_ControlFile(String fileName) {
	fileName = AFX(GetFileFolder(fileName), "ControlFile.in");
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	dt.solver = Hydro::HAMS_WAMIT;
	
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
			dt.h = f.GetDouble(1);
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
			if (!IsNull(dt.Nf)) {
				if (dt.Nf != abs(Nf))
					throw Exc(Format(t_("Number of frequencies %d does not match with previously loaded %d"), Nf, dt.Nf));
			} else 
				dt.Nf = abs(Nf);
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
				double maxF = minF + f.GetDouble(1)*(Nf-1);
				
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
			if (dt.w.IsEmpty()) {
				dt.w = pick(data);
				/*for (int ifr = 0; ifr < dt.Nf; ++ifr) 
					dt.T[ifr] = 2*M_PI/dt.w[ifr];*/
			} else if (!CompareDecimals(dt.w, data, 3))
				throw Exc(t_("List of frequencies does not match with previously loaded"));

		} else if (var == "Number_of_headings") {
			int Nh = f.GetInt(1);
			if (!IsNull(dt.Nh)) {
				if (dt.Nh != abs(Nh))
					throw Exc(Format(t_("Number of headings %d does not match with previously loaded %d"), Nh, dt.Nh));
			} else 
				dt.Nh = abs(Nh);
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
				double maxH = minH + f.GetDouble(1)*(Nh-1);
				
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
			if (dt.head.IsEmpty()) 
				dt.head = pick(data);
			else if (!CompareDecimals(dt.head, data, 2))
				throw Exc(t_("List of headings does not match with previously loaded"));

		} else if (var == "Reference_body_centre" || var == "Reference_body_center") {
			if (f.size() < 4)
				throw Exc(t_("Lack of data in Reference_body_center"));
			dt.msh.SetCount(1);
			Body &body = dt.msh[0];
			body.dt.c0[0] = f.GetDouble(1);
			body.dt.c0[1] = f.GetDouble(2);
			body.dt.c0[2] = f.GetDouble(3);
			String meshFile = AFX(GetFileFolder(fileName), "HullMesh.pnl");
			if (FileExists(meshFile))
				body.dt.fileName = meshFile;
			String lidFile = AFX(GetFileFolder(fileName), "WaterplaneMesh.pnl");
			if (FileExists(lidFile))
				body.dt.lidFile = lidFile;
		}
	}
	return LoadHydrostatic(AFX(GetFileFolder(fileName), "Hydrostatic.in"));
}

bool Hams::LoadHydrostatic(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	dt.msh.SetCount(1);
	Body &body = dt.msh[0];
	
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
	
	while (!f.IsEof()) {
		f.GetLine();
		
		if (f.size() < 1)
			continue;	
		
		String line = Trim(f.GetText());
	
		if (line == "Centre of Gravity:" || line == "Center of Gravity:") {
			f.GetLine();
			if (f.size() < 3)
				throw Exc(t_("Centre of Gravity data is not complete"));
			body.dt.cg[0] = f.GetDouble(0);
			body.dt.cg[1] = f.GetDouble(1);
			body.dt.cg[2] = f.GetDouble(2);
		} else if (line == "Body Mass Matrix:") 
			InMatrix(f, body.dt.M);
		else if (line == "External Linear Damping Matrix:") 
			InMatrix(f, body.dt.Dlin);
		else if (line == "External Quadratic Damping Matrix:") 
			InMatrix(f, body.dt.Dquad);	
		else if (line == "Hydrostatic Restoring Matrix:") 
			InMatrix(f, body.dt.C);	
		else if (line == "External Restoring Matrix:") 
			InMatrix(f, body.dt.Cadd);	
	}
	return true;
}

void Hams::SaveCase(String folderBase, bool bin, int numCases, int numThreads, bool x0z, bool y0z, const UArray<Body> &lids) const {
	SaveFolder0(folderBase, bin, 1, true, numThreads,  x0z, y0z, lids);
	if (numCases > 1)
		SaveFolder0(folderBase, bin, numCases, false, numThreads, x0z, y0z, lids);
}

void Hams::SaveFolder0(String folderBase, bool bin, int numCases, bool deleteFolder, int numThreads, bool x0z, bool y0z, const UArray<Body> &lids) const {
	BeforeSaveCase(folderBase, numCases, deleteFolder);
	
	UVector<int> valsf;
	int ifr = 0;
	if (numCases > 1)  
		valsf = NumSets(dt.Nf, numCases);
	
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
		String source = Bem().hamsBodyPath;
		meshName = GetFileName(source);
		String destNew = AFX(folderBase, meshName);
		if (!FileCopy(source, destNew)) 
			throw Exc(Format(t_("Problem copying Hams WAMIT_MeshTran exe file from '%s'"), Bem().hamsBodyPath));				
	} 
		
	//String sumcases;
	for (int i = 0; i < numCases; ++i) {
		String folder;
		UVector<double> freqs;
		if (numCases > 1) {
			folder = AFX(folderBase, Format("HAMS_Part_%d", i+1));
			if (!DirectoryCreateX(folder))
				throw Exc(Format(t_("Problem creating '%s' folder"), folder));
			//sumcases << " " << AFX(folder, "Nemoh.cal");
			Upp::Block(dt.w, freqs, ifr, valsf[i]);
			ifr += valsf[i];
		} else {
			folder = folderBase;
			freqs = clone(dt.w);
		}
		String folderInput = AFX(folder, "Input");
		if (!DirectoryCreateX(folderInput))
			throw Exc(Format(t_("Problem creating '%s' folder"), folderInput));
		
		bool irrRemoval = !lids.IsEmpty() && !lids[0].IsEmpty();
		
		if (IsNull(numThreads) || numThreads <= 0)
			numThreads = 8;
		Save_ControlFile(folderInput, freqs, numThreads, irrRemoval);
		Save_Hydrostatic(folderInput);
	
		int ib = 0;			// Just one file
		
		if (y0z == true && x0z == true) 
			y0z = false;

		String dest = AFX(folderInput, "HullMesh.pnl");
		Body::SaveAs(dt.msh[ib], dest, Body::HAMS_PNL, Body::UNDERWATER, dt.rho, dt.g, y0z, x0z);
		
		if (irrRemoval) {
			dest = AFX(folderInput, "WaterplaneMesh.pnl");
			Body::SaveAs(lids[0], dest, Body::HAMS_PNL, Body::ALL, dt.rho, dt.g, y0z, x0z);
		}
		
		Save_Settings(folder, lids);
		
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
	
	out << "echo Start: \%date\% \%time\% > time.txt\n";
	
	if (!IsNull(caseFolder))
		out << "cd \"" << caseFolder << "\"\n";
	
	out << "\"" << meshName << "\"\n";
	out << "\"" << solvName << "\"\n";
	
	out << "echo End:   \%date\% \%time\% >> time.txt\n";
}

void Hams::OutMatrix(FileOut &out, String header, const Eigen::MatrixXd &mat) {
	bool isvoid = mat.size() != 36;

	out << "\n " << header << ":";
	for (int y = 0; y < 6; ++y) {
		out << "\n";
		for (int x = 0; x < 6; ++x)
			out << "   " << Format("%.5E", isvoid ? 0. : mat(x, y));
	}
}

void Hams::InMatrix(LineParser &f, Eigen::MatrixXd &mat) {
	if (mat.size() != 36)
		mat.resize(6, 6);
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
	out << " Center of Gravity:";
	
	if (dt.msh.IsEmpty())
		throw Exc(t_("No bodies found"));
	
	const Body &b = dt.msh[0];
	out << Format("\n  %.15E  %.15E  %.15E", b.dt.cg[0], b.dt.cg[1], b.dt.cg[2]);

	OutMatrix(out, "Body Mass Matrix", b.dt.M);
	OutMatrix(out, "External Linear Damping Matrix", b.dt.Dlin);
	OutMatrix(out, "External Quadratic Damping Matrix", b.dt.Dquad);
	OutMatrix(out, "Hydrostatic Restoring Matrix", b.dt.C);
	OutMatrix(out, "External Restoring Matrix", b.dt.Cadd);
	out << "\n";
}


void Hams::Save_Settings(String folderInput, const UArray<Body> &lids) const {
	String fileName = AFX(folderInput, "Settings.ctrl");
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to create '%s'"), fileName));
	
	Body mesh;
	String res = Body::Load(mesh, AFX(folderInput, "Input", "HullMesh.pnl"), dt.rho, dt.g, Null, Null, false);
	if (!res.IsEmpty())
		throw Exc(res);
	
	if (!lids.IsEmpty() && !lids[0].IsEmpty()) {
		//Body lid;
		//lid.dt.mesh.AddWaterSurface(mesh.dt.mesh, mesh.dt.under, 'f', Bem().roundVal, Bem().roundEps); 
		//lid.AfterLoad(dt.rho, dt.g, false, false);
		
		mesh.Append(lids[0].dt.mesh, dt.rho, dt.g);
	}
	Body::SaveAs(mesh, AFX(folderInput, "Input", "mesh.gdf"), Body::WAMIT_GDF, Body::ALL, dt.rho, dt.g, false, false);	
	
	out << dt.g << "\n";
	out << dt.rho << "\n";
	out << ".\\Input\\mesh.gdf" << "\n";
	out << dt.msh[0].dt.c0[0] << "   " << dt.msh[0].dt.c0[1] << "   " << dt.msh[0].dt.c0[2] << "\n";
	out << dt.msh[0].dt.cg[0] << "   " << dt.msh[0].dt.cg[1] << "   " << dt.msh[0].dt.cg[2] << "\n";
	out << "\n"
		<< "The format is as below:\n"
		<< "gravity acceleration\n"
		<< "sea water density\n"
		<< "mesh file name (gdf format)\n"
		<< "Center of rotation\n"
		<< "Center of gravity";
}

void Hams::Save_ControlFile(String folderInput, const UVector<double> &freqs,
					int numThreads, bool remove_irr_freq) const {
	String fileName = AFX(folderInput, "ControlFile.in");
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to create '%s'"), fileName));	
	
	out << "   --------------HAMS Control file---------------"
		   "\n";
	double depth = dt.h;
	if (depth >= 400)
		depth = -1;
	out << "\n   Waterdepth  " << Format("%.4f", depth) << "D0";
	out << "\n"
		   "\n   #Start Definition of Wave Frequencies"
		   "\n    0_inf_frequency_limits      1"
		   "\n    Input_frequency_type        3"
		   "\n    Output_frequency_type       3";
	out << "\n    Number_of_frequencies      " << freqs.size();
	
	out << "\n    ";
	for (const auto &freq : freqs)
		out << Format("%.4f ", freq);	
	out << "\n   #End Definition of Wave Frequencies"
		   "\n"
		   "\n   #Start Definition of Wave Headings";
	out << "\n    Number_of_headings         " << dt.Nh;
	UVector<double> headings;
	LinSpaced(headings, dt.Nh, First(dt.head), Last(dt.head));
	out << "\n    ";
	for (const auto &heading : headings)
		out << Format("%.4f ", heading);	
	out << "\n   #End Definition of Wave Headings"
		   "\n";
	out << "\n    Reference_body_center " << Format("%11.3f %11.3f %11.3f", 
								dt.msh[0].dt.c0[0], dt.msh[0].dt.c0[1], dt.msh[0].dt.c0[2]) << 
		   "\n    Reference_body_length   1.D0"
		   "\n    Wave_diffrac_solution    2";
	
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