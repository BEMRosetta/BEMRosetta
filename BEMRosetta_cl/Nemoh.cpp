#include "BEMRosetta.h"


bool Nemoh::Load(String file, double) {
	this->file = file;
	this->name = GetFileTitle(GetFileFolder(file));
	folder = GetFileFolder(file);
	len = 1;
	String ext = GetFileExt(file);
	if (ext == ".cal")
		code = NEMOH;
	else
		code = SEAFEM_NEMOH;
	
	try {
		String fileCal;
		Print("\n\n" + Format("Loading '%s'", file));
		if (code == NEMOH) 
			fileCal = file;
		else 
			fileCal = AppendFileName(folder, "Nemoh_output/Nemoh.cal");
		if (!Load_Cal(fileCal)) {
			PrintError("\n" + Format("File '%s' not found", fileCal));
			return false;
		}
		
		String fileRad, folderForces;
		if (code == NEMOH) {
			Print("\n- Hydrostatics file(s) '/Mesh/Hydrostatics*.dat'");
			if (!Load_Hydrostatics())
				PrintError(": **Not found**");
			Print("\n- KH file(s) '/Mesh/KH*.dat'");
			if (!Load_KH())
				PrintError(": **Not found**");
			fileRad = AppendFileName(folder, AppendFileName("Results", "RadiationCoefficients.tec"));
			folderForces = folder;
		} else {
			if (!Load_Inf(file)) {
				PrintError("\n" + Format("File '%s' not found", file));
				return false;
			}
			fileRad = AppendFileName(folder, AppendFileName("Nemoh_output/Results", "RadiationCoefficients.tec"));
			folderForces = AppendFileName(folder, "Nemoh_output");
		}
		
		Print("\n- Radiation file 'RadiationCoefficients.tec'");
		if (!Load_Radiation(fileRad))
			PrintError(": **Not found**");
		
		Print("\n- Excitation force file 'ExcitationForce.tec'");
		if (!Load_Excitation(folderForces))
			PrintError(": **Not found**");
		
		if (!dof.IsEmpty()) {
			Print("\n- Diffraction force file 'DiffractionForce.tec'");
			if (!Load_Diffraction(folderForces))
				PrintError(": **Not found**");
			Print("\n- Froude Krylov file 'FKForce.tec'");
			if (!Load_FroudeKrylov(folderForces))
				PrintError(": **Not found**");
		}
		if (code == NEMOH) {
			Print("\n- IRF file(s) 'IRF.tec'");
			if (!Load_IRF(AppendFileName(folder, AppendFileName("Results", "IRF.tec"))))
				PrintError(": **Not found**");
		}
		AfterLoad();
	} catch (Exc e) {
		PrintError("\nError: " + e);
		lastError = e;
		return false;
	}
	
	return true;
}

bool Nemoh::Load_Cal(String fileName) {
	FileIn in(fileName);
	if (!in.IsOpen())
		return false;
	
	rho = g = h = Null;
	
	String line;
	FieldSplit f;
	while(!in.IsEof()) {
		line = in.GetLine();
		
		if (line.Find("Fluid specific volume") >= 0) 
			rho = ScanDouble(line);
		else if (line.Find("Gravity") >= 0) 
			g = ScanDouble(line);
		else if (line.Find("Water depth") >= 0) {
			h = ScanDouble(line);
			if (h == 0)
				h = INFINITY;
		} else if (line.Find("Number of bodies") >= 0) 
			Nb = ScanInt(line);
		else if (line.Find("Name of mesh file") >= 0) { 
			f.Load(line);
			names << GetFileTitle(f.GetText(0));
		} else if (line.Find("Number of wave frequencies") >= 0) {
			f.Load(line);
			Nf = f.GetInt(0);  						// Number of wave frequencies
			double minF = f.GetDouble(1);
			double maxF = f.GetDouble(2);
        	LinSpaced(w, Nf, minF, maxF); 			// Wave frequencies
        	T.SetCount(Nf);
        	for (int i = 0; i < Nf; ++i)
        		T[i] = 2*M_PI/w[i];  						// Wave periods
		} else if (line.Find("Number of wave directions") >= 0) {
			f.Load(line);
			Nh = f.GetInt(0);  						// Number of wave headings
			double minD = f.GetDouble(1);
			double maxD = f.GetDouble(2);
        	LinSpaced(head, Nh, minD, maxD); 				// Wave frequencies
		}
	}
	if (Nb == 0 || Nf == 0 || Nh == 0 || IsNull(rho) || IsNull(g) || IsNull(h))
		throw Exc(Format("Wrong format in Nemoh file '%s'", fileName));
		
	return true;
}

bool Nemoh::Load_Inf(String fileName) {
	if (Nb != 1)
		throw Exc(Format("SeaFEM_Nemoh only allows one body, found %d", Nb));
		
	cg.setConstant(3, 1, nan(""));
	cb.setConstant(3, 1, nan(""));
	Vo.SetCount(1, nan(""));
	C.SetCount(1);
	C[0].setConstant(6, 6, nan(""));   
	
	FileIn in(fileName);
	if (!in.IsOpen())
		return false;
	
	String line;
	while(!in.IsEof()) {
		line = in.GetLine();
		int pos;
		if ((pos = line.FindAfter("XG [m]=")) >= 0) 
			cg(0, 0) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("YG [m]=")) >= 0) 
			cg(1, 0) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("ZG [m]=")) >= 0) 
			cg(2, 0) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("XC [m]=")) >= 0) 
			cb(0, 0) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("YC [m]=")) >= 0) 
			cb(1, 0) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("ZC [m]=")) >= 0) 
			cb(2, 0) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("Displacement [m3]=")) >= 0) 
			Vo[0] = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [3][3] [N/m]=")) >= 0) 
			C[0](2, 2) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [3][4] [N/rad]=")) >= 0) 
			C[0](2, 3) = C[0](3, 2) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [3][5] [N/rad]=")) >= 0) 
			C[0](2, 4) = C[0](4, 2) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [4][4] [Nm/rad]=")) >= 0) 
			C[0](3, 3) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [4][5] [Nm/rad]=")) >= 0) 
			C[0](3, 4) = C[0](4, 3) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [4][6] [Nm/rad]=")) >= 0) 
			C[0](3, 5) = C[0](5, 3) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [5][5] [Nm/rad]=")) >= 0) 
			C[0](4, 4) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [5][6] [Nm/rad]=")) >= 0) 
			C[0](4, 5) = C[0](5, 4) = ScanDouble(line.Mid(pos));
	}
	return true;	
}
	
bool Nemoh::Load_Hydrostatics() {
	cg.setConstant(3, Nb, nan(""));
	cb.setConstant(3, Nb, nan(""));
	Vo.SetCount(Nb, nan(""));
	String line;
	FieldSplit f;
	for (int b = 0; b < Nb; ++b) {
	    String fileHydro;
	    if (Nb == 1)
	        fileHydro = AppendFileName(folder, AppendFileName("Mesh", "Hydrostatics.dat"));
	    else
	        fileHydro = AppendFileName(folder, AppendFileName("Mesh", Format("Hydrostatics_%d.dat", b)));
	    
	    FileIn in(fileHydro);
	    if (!in.IsOpen())
	        return false;
	    for (int i = 0; i < 3 && !in.IsEof(); ++i) {
			f.Load(in.GetLine());
			cg(i, b) = f.GetDouble(6);
			cb(i, b) = f.GetDouble(2);
	    }
		f.Load(in.GetLine());
	    Vo[b] = f.GetDouble(2); 		// Displacement volume
	}
	return true;
}

bool Nemoh::Load_KH() {
	C.SetCount(Nb);
	FieldSplit f;
	for (int ib = 0; ib < Nb; ++ib) {
	    String fileKH;
	    if (Nb == 1)
	        fileKH = AppendFileName(folder, AppendFileName("Mesh", "KH.dat"));
	    else
	        fileKH = AppendFileName(folder, AppendFileName("Mesh", Format("KH_%d.dat", ib)));
	
		C[ib].setConstant(6, 6, nan(""));    
	    FileIn in(fileKH);
	    if (!in.IsOpen())
	        return false;
		for (int i = 0; i < 6 && !in.IsEof(); ++i) {
			f.Load(in.GetLine());
			for (int ifr = 0; ifr < 6; ++ifr)
				C[ib](i, ifr) = f.GetDouble(ifr);
		}
	}
	return true;
}

bool Nemoh::Load_Radiation(String fileName) {
	FileIn in(fileName);
	if (!in.IsOpen())
		return false;
	String line;
	FieldSplit f;
	dof.Clear();	 dof.SetCount(Nb, 0);
	in.GetLine();
	while(!in.IsEof()) {
		line = in.GetLine();
	    if (line.Find("Motion of body") >= 0)
	        break;
	    f.Load(line);
	    int ibody = f.GetInt(1) - 1;
	    int ndof = f.GetInt(2);
		dof[ibody] = ndof;    
	}
	A.SetCount(Nf);
	B.SetCount(Nf);
	for (int k = 0; k < Nf; ++k) {
		A[k].setConstant(Nb*6, Nb*6, nan(""));
		B[k].setConstant(Nb*6, Nb*6, nan(""));
	}
	int ibodydof = 0;
	for (int ibody = 0; ibody < Nb; ++ibody) {
		for (int idof = 0; idof < dof[ibody]; ++idof) {
			for (int k = 0; k < Nf; ++k) {	
				f.Load(in.GetLine());
				for (int df = 0; df < dof[ibody]; ++df) {		
					A[k](ibodydof, df) = f.GetDouble(1 + 2*df);
	        		B[k](ibodydof, df) = f.GetDouble(2 + 2*df);
				}
	        }
	        ++ibodydof;
	        in.GetLine();
	    }
	}
	return true;
}

bool Nemoh::Load_Excitation(String folder) {	
	return Load_Forces(ex, folder, "ExcitationForce.tec", "Diffraction force");
}

bool Nemoh::Load_Diffraction(String folder) {
	return Load_Forces(sc, folder, "DiffractionForce.tec", "Diffraction force");
}

bool Nemoh::Load_FroudeKrylov(String folder) {
	return Load_Forces(fk, folder, "FKForce.tec", "FKforce");
}

bool Nemoh::Load_Forces(Forces &fc, String nfolder, String fileName, String textDelim) {
	Initialize_Forces(fc);
	FileIn in(AppendFileName(nfolder, AppendFileName("Results", fileName)));
	if (!in.IsOpen())
		return false;
	String line;
	FieldSplit f;
	while(!in.IsEof()) {
		line = in.GetLine();
		if (line.Find(textDelim) >= 0)
	        break;
	}
	for (int h = 0; h < Nh; ++h) {
		int ifr = 0;
		while(!in.IsEof()) {
			line = in.GetLine();
			if (line.Find(textDelim) >= 0)
				break;
			f.Load(line);
			int ib = 0, idof = 0, ibdof = 0;
			for (int i = 0; i < Nb*6; ++i) {
				double ma = fc.ma[h](ifr, ibdof) = f.GetDouble(1 + 2*i);	// Magnitude of exciting force
				double ph = fc.ph[h](ifr, ibdof) = -f.GetDouble(1 + 2*i + 1);//	Phase of exciting force (-ph, since NEMOH's x-dir is flipped)
				fc.re[h](ifr, ibdof) = ma*cos(ph);  // Real part of exciting force
				fc.im[h](ifr, ibdof) = ma*sin(ph);  // Imaginary part of exciting force
				idof++;
				ibdof++;
				if (idof >= dof[ib]) {
					idof = 0;
					ib++;
					ibdof = 6*ib;
					if (ib >= Nb)
						break;
				}
			}
			ifr++;
		}
	}
	return true;
}

bool Nemoh::Load_IRF(String fileName) {
	FileIn in(fileName);
	if (!in.IsOpen())
		return false;
	String line;
	FieldSplit f;	
	Awinf.setConstant(Nb*6, Nb*6, nan(""));
	int ibodydof = 0;
	for (int ibody = 0; ibody < Nb; ++ibody) {
		for (int idof = 0; idof < dof[ibody]; ++idof) {
			while(!in.IsEof()) {
				line = in.GetLine();	
				if (line.Find("Zone t=") >= 0) 
					break;
			}
			line = in.GetLine();	
			f.Load(line);
			for (int df = 0; df < dof[ibody]; ++df) 
				Awinf(ibodydof, df) = f.GetDouble(1 + 2*df);
			
			++ibodydof;
		}
	}
	return true;
}

void Nemoh::Save(String file) {
	throw Exc("Option not implemented");
}		