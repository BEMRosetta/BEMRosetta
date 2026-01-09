// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#include <STEM4U/Utility.h>
#include <STEM4U/SeaWaves.h>


String Hams::Load(String filein, bool onlycase, Function <bool(String, int)> Status) {
	dt.file = filein;	
	dt.name = GetFileTitle(filein);
	
	dt.dimen = true;
	dt.x_w = dt.y_w = 0;
	
	String baseFolder = GetFileFolder(GetFileFolder(filein));
	
	try {
		if (GetFileExt(filein) != ".in") 
			return Format(t_("File '%s' is not of HAMS type"), filein);
			
		BEM::Print("\n\n" + Format(t_("Loading '%s'"), filein));

		dt.solver = Hydro::HAMS_WAMIT;
		
		int output_frequency_type = 0;
		String controlFile = AFX(baseFolder, "Input", "ControlFile.in");
		if (FileExists(controlFile)) 
			output_frequency_type = Load_ControlFile(controlFile);
		
		int iperout;
		switch (output_frequency_type) {// HAMS						-> Wamit
		case 1:	iperout = 3;	break;	// deepwater wave number    -> Infinite-depth wavenumber	KL = ω2L∕g
		case 2:	iperout = 4;	break;	// finite-depth wave number	-> Finite-depth wavenumber
		case 3:	iperout = 2;	break;	// wave frequency			-> Frequency 	ω = 2π∕T
		case 4:	iperout = 1;	break;	// wave period				-> Period in seconds 	T
		case 5:	iperout = 5;	break;	// 							   wave length
		}

		if (!onlycase) {
			String wamitFile = AFX(baseFolder, "Output", "Wamit_format", "Buoy.1");
			String ret = Wamit::Load(wamitFile, true, iperout, Status);
			if (!ret.IsEmpty()) 
				return ret;
		}
		
		if (dt.Nb == 1) {
			String hydroFile = AFX(baseFolder, "Input", "Hydrostatic.in");
			if (FileExists(hydroFile)) 
				LoadHydrostatic(hydroFile, 0);
		} else {
			for (int ib = 0; ib < dt.Nb; ++ib) {
				String hydroFile = AFX(baseFolder, "Input", Format("Hydrostatic_%d.in", ib+1));
				if (FileExists(hydroFile)) 
					LoadHydrostatic(hydroFile, ib);
			}
		}
		String settingsFile = AFX(baseFolder, "Settings.ctrl");
		if (FileExists(settingsFile))
			Load_Settings(settingsFile);
		
		// If Hydrostatic.in exists and the .hst file is filled with zeroes
		String hydrostaticFile = AFX(baseFolder, "Hydrostatic.in");
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
			if (iszero && !onlycase && !Load_HydrostaticBody(hydrostaticFile, dt.rho*dt.g))
				return Format(t_("Problem loading Hydrostatic file '%s'"), hydrostaticFile);
		}
			
		if (IsNull(dt.Nb))
			return t_("No body found");

	} catch (Exc e) {
		BEM::PrintError("\nError: " + e);
		return e;
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
	
	f.Load(in.GetLine());	
	dt.msh[0].dt.c0.x = f.GetDouble(0);
	dt.msh[0].dt.c0.y = f.GetDouble(1);
	dt.msh[0].dt.c0.z = f.GetDouble(2);
	
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
		ret << t_("HAMS issue: Just first body will be saved");

	return ret;
}

int Hams::Load_ControlFile(String fileName) {
	dt.g = 9.80665;		// Is constant

	String folder = GetFileFolder(fileName);
	fileName = AFX(folder, "ControlFile.in");
	FileInLine in(fileName);
	if (!in.IsOpen())
		return 0;
	
	dt.solver = Hydro::HAMS_WAMIT;	
	dt.Nb = 1;				// If it is HAMS
	dt.msh.SetCount(1);	
	dt.msh[0].dt.fileName = AFX(folder, "HullMesh.pnl");
	dt.msh[0].dt.name = "Body";
		
	String line;
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
	
	bool ismrel = false;
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
			if (input_frequency_type < 1 || input_frequency_type > 5)
				throw Exc(t_("Wrong input_frequency_type"));
		} else if (var == "Output_frequency_type") {
			output_frequency_type = f.GetInt(1); 
			if (output_frequency_type < 1 || output_frequency_type > 5)
				throw Exc(t_("Wrong output_frequency_type"));
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

				for (int i = 0; i < f.size(); ++i)
					data << f.GetDouble(i);
				if (Nf != data.size())
					throw Exc(Format(t_("Number of frequencies mismatch (%d != %d)"), Nf, data.size()));
			}
			for (double &dat : data) {
				switch (input_frequency_type) {
				case 1:	dat = sqrt(dt.g*dat);										break;
				case 2:	dat = SeaWaves::FrequencyFromWaveNumber(dat, dt.h, dt.g);	break;
				case 3:	break;
				case 4:	dat = 2*M_PI/dat;											break;
				case 5:	dat = SeaWaves::FrequencyFromWaveLength(dat, dt.h, dt.g);	break;
				}
			}
			if (dt.w.IsEmpty())
				dt.w = pick(data);
			else if (!CompareDecimals(dt.w, data, 3))
				throw Exc(t_("List of frequencies does not match with previously loaded"));
		} else if (var == "Number_of_bodies") {
			ismrel = true;
			dt.Nb = f.GetInt(1);
			dt.msh.SetCount(dt.Nb);
			if (dt.Nb > 1) {
				for (int ib = 0; ib < dt.Nb; ++ib) {
					dt.msh[ib].dt.fileName = AFX(folder, Format("HullMesh_%d.pnl", ib+1));
					dt.msh[ib].dt.name = Format("Body_%d", ib+1);
				}
			}
			UVector<Point3D> lcs(dt.Nb);
			for (int ib = 0; ib < dt.Nb; ++ib) {
				f.GetLine();
				lcs[ib].x = f.GetDouble(1);
				lcs[ib].y = f.GetDouble(2);
				lcs[ib].z = f.GetDouble(3);
				if (f.GetDouble(4) != 0)
					throw Exc(t_("BEMRosetta does not support rotations in local coordinate system LCS"));
				
				Body &b = dt.msh[ib];
				String ret = Body::Load(b, b.dt.fileName, dt.rho, Bem().g, Null, Null, false);
				if (!IsEmpty(ret))
					BEM::PrintWarning(Format(t_("Problem loading mesh '%s': %s"), b.dt.fileName, ret));
				b.dt.mesh.Translate(lcs[ib].x, lcs[ib].y, lcs[ib].z);
				b.AfterLoad(dt.rho, Bem().g, false, false, false, false);
			}
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
				
				for (int i = 0; i < f.size(); ++i)
					data << f.GetDouble(i);
			}
			if (dt.head.IsEmpty())
				dt.head = pick(data);
			else if (!CompareDecimals(dt.head, data, 2))
				throw Exc(t_("List of headings does not match with previously loaded"));

		} else if (var.StartsWith("Reference_body_centre") || var.StartsWith("Reference_body_center")) {
			int id = var.FindAfter("Reference_body_center_");
			int ib;
			if (id < 0)
				ib = 0;
			else {
				ib = ScanInt(var.Mid(id));
				if (ib > dt.Nb || ib < 1)
					throw Exc(t_("Reference_body_center of wrong body"));
				ib--;
			}
			if (f.size() < 4)
				throw Exc(t_("Lack of data in Reference_body_center"));
			
			Body &body = dt.msh[ib];
			body.dt.c0[0] = f.GetDouble(1);
			body.dt.c0[1] = f.GetDouble(2);
			body.dt.c0[2] = f.GetDouble(3);
		} else if (var == "Reference_body_length")
			dt.len = f.GetDouble(1);
		else if (var == "Number_of_field_points") {
			int numPoints = f.GetInt(1);
			listPointsTemp.SetCount(numPoints);
			for (int r = 0; r < numPoints; ++r) {
				f.GetLine();
				if (f.IsEof())
					throw Exc(Format(t_("Field points list is incomplete. read %d from %d"), r, numPoints));
				listPointsTemp[r].x = f.GetDouble(0);
				listPointsTemp[r].y = f.GetDouble(1);
				listPointsTemp[r].z = f.GetDouble(2);
			}
		} else if (var == "If_remove_irr_freq") {
			int remove_irr_freq = f.GetInt(1);	
			if (remove_irr_freq > 0) {
				if (dt.Nb == 1)	
					dt.msh[0].dt.lidFile = AFX(folder, "WaterplaneMesh.pnl");
				else {
					for (int ib = 0; ib < dt.Nb; ++ib)
						dt.msh[ib].dt.lidFile = AFX(folder, Format("WaterplaneMesh_%d.pnl", ib+1));
				}
			}
		}
	}
	if (!ismrel) {
		Body &b = dt.msh[0];
		String ret = Body::Load(b, b.dt.fileName, dt.rho, Bem().g, Null, Null, false);
		if (!IsEmpty(ret))
			BEM::PrintWarning(Format(t_("Problem loading mesh '%s': %s"), b.dt.fileName, ret));
		b.AfterLoad(dt.rho, Bem().g, false, false, false, false);
	}
	return output_frequency_type;
}

bool Hams::LoadHydrostatic(String fileName, int ib) {
	ASSERT(dt.Nb > ib);
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	if (dt.msh.size() < dt.Nb)
		dt.msh.SetCount(dt.Nb);
	Body &body = dt.msh[ib];
	
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
	
	while (!f.IsEof()) {
		f.GetLine();
		
		if (f.size() < 1)
			continue;	
		
		String line = Trim(f.GetText());
	
		if (line.StartsWith("Centre of Gravity:") || line.StartsWith("Center of Gravity:")) {
			f.GetLine();
			if (f.size() < 3)
				throw Exc(t_("Centre of Gravity data is not complete"));
			body.dt.cg[0] = f.GetDouble(0);
			body.dt.cg[1] = f.GetDouble(1);
			body.dt.cg[2] = f.GetDouble(2);
		} else if (line.StartsWith("Body Mass Matrix:"))
			InMatrix(f, body.dt.M);
		else if (line.StartsWith("External Linear Damping Matrix:")) 
			InMatrix(f, body.dt.Dlin);
		else if (line.StartsWith("External Quadratic Damping Matrix:")) 
			InMatrix(f, body.dt.Dquad);	
		else if (line.StartsWith("Hydrostatic Restoring Matrix:")) 
			InMatrix(f, body.dt.C);	
		else if (line.StartsWith("External Restoring Matrix:")) 
			InMatrix(f, body.dt.Cadd);	
	}
	return true;
}

void Hams::SaveCase(String folderBase, bool bin, int numCases, int numThreads, bool x0z, bool y0z, 
					const UArray<Body> &lids,const UVector<Point3D> &listPoints, bool ismrel) const {
	SaveFolder0(folderBase, bin, 1, true, numThreads,  x0z, y0z, lids, listPoints, ismrel);
	if (numCases > 1)
		SaveFolder0(folderBase, bin, numCases, false, numThreads, x0z, y0z, lids, listPoints, ismrel);
}

void Hams::SaveFolder0(String folderBase, bool bin, int numCases, bool deleteFolder, int numThreads, bool x0z, bool y0z, 
		const UArray<Body> &lids, const UVector<Point3D> &listPoints, bool ismrel) const {
	BeforeSaveCase(folderBase, numCases, deleteFolder);
	
	UVector<int> valsf;
	int ifr = 0;
	if (numCases > 1)  
		valsf = NumSets(dt.Nf, numCases);
	
	String solvName = "HAMS_x64.exe";
	if (bin) {
		String source;
		if (!ismrel)
			source = Bem().hamsPath;
		else
			source = Bem().hamsmrelPath;
		solvName = GetFileName(source);
		String destNew = AFX(folderBase, solvName);
		if (!FileCopy(source, destNew)) 
			throw Exc(Format(t_("Problem copying exe file from '%s'"), source));
		source = AFX(GetFileFolder(source), "libiomp5md.dll");		
		destNew = AFX(folderBase, "libiomp5md.dll");		
		if (!FileCopy(source, destNew)) 
			throw Exc(Format(t_("Problem copying dll file from '%s'"), source));					
	}
	 
	String meshName;
	if (bin && !ismrel) {
		String source = Bem().hamsBodyPath;
		meshName = GetFileName(source);
		String destNew = AFX(folderBase, meshName);
		if (!FileCopy(source, destNew)) 
			throw Exc(Format(t_("Problem copying Hams WAMIT_MeshTran exe file from '%s'"), Bem().hamsBodyPath));				
	} 
		
	for (int i = 0; i < numCases; ++i) {
		String folder;
		UVector<double> freqs;
		if (numCases > 1) {
			folder = AFX(folderBase, Format("HAMS_Part_%d", i+1));
			if (!DirectoryCreateX(folder))
				throw Exc(Format(t_("Problem creating '%s' folder"), folder));
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
		Save_ControlFile(folderInput, freqs, numThreads, irrRemoval, listPoints, ismrel);
		Save_Hydrostatic(folderInput);
		
		if (y0z == true && x0z == true) 
			y0z = false;

		if (!ismrel) {
			int ib = 0;
			
			String dest = AFX(folderInput, "HullMesh.pnl");
			Body::SaveAs(dt.msh[ib], dest, Body::HAMS_PNL, Body::UNDERWATER, dt.rho, dt.g, y0z, x0z);
			
			if (irrRemoval) {
				dest = AFX(folderInput, "WaterplaneMesh.pnl");
				Body::SaveAs(lids[ib], dest, Body::HAMS_PNL, Body::ALL, dt.rho, dt.g, y0z, x0z);
			}
		} else {
			for (int ib = 0; ib < dt.Nb; ++ib) {
				String dest = AFX(folderInput, Format("HullMesh_%d.pnl", ib+1));
				Body::SaveAs(dt.msh[ib], dest, Body::HAMS_PNL, Body::UNDERWATER, dt.rho, dt.g, y0z, x0z);
				
				if (irrRemoval) {
					dest = AFX(folderInput, Format("WaterplaneMesh_%d.pnl", ib+1));
					Body::SaveAs(lids[ib], dest, Body::HAMS_PNL, Body::ALL, dt.rho, dt.g, y0z, x0z);
				}
			}
		}
		
		if (!ismrel)
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
	
	if (!meshName.IsEmpty())
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
			out << Format("  %12.5E", isvoid ? 0. : mat(x, y));
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
					int numThreads, bool remove_irr_freq, const UVector<Point3D> &listPoints, bool ismrel) const {
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
	for (const double &freq : freqs)
		out << Format("%.4f ", freq);	
	out << "\n   #End Definition of Wave Frequencies"
		   "\n";
	
	if (ismrel) { 
   		out << "\n   # Start of definition of bodies for multi-body interaction";
    	out << "\n    Number_of_bodies          " << FormatInt(dt.Nb);
    	for (int ib = 0; ib < dt.Nb; ++ib) {
    		out << "\n    LCS_" << FormatInt(ib+1) << "                     0.000        0.000        0.000        0.000";
    		if (ib == 0)
    			out << "             # The first three entries refer to the x, y and z coordinates of the origin of LCS w.r.t origin of GCS and the fourth entry refers to the rotation of LCS x-axis w.r.t GCS x-axis";
    	}
    	out << "\n   # End of definition of bodies for multi-body interaction";
		out << "\n";		   
	}
	
	out << "\n   #Start Definition of Wave Headings";
	out << "\n    Number_of_headings         " << dt.Nh;
	UVector<double> headings;
	LinSpaced(headings, dt.Nh, First(dt.head), Last(dt.head));
	out << "\n    ";
	for (const double &heading : headings)
		out << Format("%.4f ", heading);	
	out << "\n   #End Definition of Wave Headings"
		   "\n";
	
	if (!ismrel || dt.Nb == 1)
		out << "\n    Reference_body_center " << Format("%11.3f %11.3f %11.3f", 
									dt.msh[0].dt.c0[0], dt.msh[0].dt.c0[1], dt.msh[0].dt.c0[2]);
	else {
		for (int ib = 0; ib < dt.Nb; ++ib)
			out << "\n    Reference_body_center_" << FormatInt(ib+1) << " " << Format("%11.3f %11.3f %11.3f", 
									dt.msh[0].dt.c0[0], dt.msh[0].dt.c0[1], dt.msh[0].dt.c0[2]);
	}
	
	out << "\n    Reference_body_length   1.D0"
		   "\n    Wave_diffrac_solution    2";
	
	out << "\n    If_remove_irr_freq      " << (remove_irr_freq ? 1 : 0);
	out << "\n    Number of threads       " << numThreads;
	
	out << "\n"
		   "\n   #Start Definition of Pressure and/or Elevation (PE)";
	if (listPoints.IsEmpty())
		out << "\n    Number_of_field_points     1                           # number of field points where to calculate PE"
			   "\n    0.000000    0.000000    0.000000    Global_coords_point_1";
	else {
		out << Format("\n    Number_of_field_points     %d                           # number of field points where to calculate PE", listPoints.size());
		for (int i = 0; i < listPoints.size(); ++i) {
			const Point3D &p = listPoints[i];
			out << Format("\n    %s    %s    %s    Global_coords_point_%d", FDS(p.x, 10), FDS(p.y, 10), FDS(p.z, 10), i+1);
		}
	}
	out << "\n   #End Definition of Pressure and/or Elevation"
   		   "\n";	   
   		   
	out << "\n    ----------End HAMS Control file---------------"
		   "\n   Input_frequency_type options:"
		   "\n   1--deepwater wave number; 2--finite-depth wave number; 3--wave frequency; 4--wave period; 5--wave length"
   		   "\n   Output_frequency_type options: same as Input_frequency_type options";
}