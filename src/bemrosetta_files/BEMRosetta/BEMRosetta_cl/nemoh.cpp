#include "BEMRosetta.h"

bool Nemoh::Load(String file, double) {
	hd().file = file;
	hd().name = GetFileTitle(GetFileFolder(file));
	folder = GetFileFolder(file);
	hd().len = 1;
	hd().dimen = true;
	hd().Nb = Null;
	
	String ext = GetFileExt(file); 
	if (ext == ".cal")
		hd().code = Hydro::NEMOH;
	else
		hd().code = Hydro::SEAFEM_NEMOH;
	
	try {
		String fileCal;
		hd().Print("\n\n" + Format(t_("Loading '%s'"), file));
		if (hd().code == Hydro::NEMOH) 
			fileCal = file;
		else 
			fileCal = AppendFileName(folder, "Nemoh_output/Nemoh.cal");
		if (!Load_Cal(fileCal)) 
			throw Exc("\n" + Format(t_("File '%s' not found"), fileCal));
		
		String fileRad, folderForces;
		if (hd().code == Hydro::NEMOH) {
			hd().Print(x_("\n- ") + t_("Hydrostatics file(s) 'Mesh/Hydrostatics*.dat'"));
			if (!Load_Hydrostatics())
				hd().PrintWarning(x_(": **") + t_("Not found") + "**");
			hd().Print(x_("\n- ") + t_("KH file(s) 'Mesh/KH*.dat'"));
			if (!Load_KH())
				hd().PrintWarning(x_(": **") + t_("Not found") + "**");
			fileRad = AppendFileName(folder, AppendFileName("Results", "RadiationCoefficients.tec"));
			folderForces = folder;
		} else {
			if (!Load_Inf(file)) 
				throw Exc("\n" + Format(t_("File '%s' not found"), file));

			fileRad = AppendFileName(folder, AppendFileName("Nemoh_output/Results", "RadiationCoefficients.tec"));
			folderForces = AppendFileName(folder, "Nemoh_output");
		} 
		
		hd().Print("\n- Radiation file 'RadiationCoefficients.tec'");
		if (!Load_Radiation(fileRad))
			hd().PrintWarning(x_(": **") + t_("Not found") + "**");
		
		hd().Print("\n- Excitation force file 'ExcitationForce.tec'");
		if (!Load_Excitation(folderForces))
			hd().PrintWarning(x_(": **") + t_("Not found") + "**");
		
		if (!hd().dof.IsEmpty()) {
			hd().Print(x_("\n- ") + t_("Diffraction force file 'DiffractionForce.tec'"));
			if (!Load_Diffraction(folderForces))
				hd().PrintWarning(x_(": **") + t_("Not found") + "**");
			hd().Print(x_("\n- ") + t_("Froude Krylov file 'FKForce.tec'"));
			if (!Load_FroudeKrylov(folderForces))
				hd().PrintWarning(x_(": **") + t_("Not found") + "**");
		}
		if (hd().code == Hydro::NEMOH) {
			hd().Print(x_("\n- ") + t_("IRF file(s) 'IRF.tec'"));
			if (!Load_IRF(AppendFileName(folder, AppendFileName("Results", "IRF.tec"))))
				hd().PrintWarning(x_(": **") + t_("Not found") + "**");
		}
		if (IsNull(hd().Nb))
			return false;
	} catch (Exc e) {
		hd().PrintError(Format("\n%s: %s", t_("Error"), e));
		hd().lastError = e;
		return false;
	}
	
	return true;
}

bool Nemoh::Load_Cal(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	hd().rho = hd().g = hd().h = Null;
	
	hd().dataFromW = true;
	
	String line;
	FieldSplit f(in);
	while(!in.IsEof()) {
		line = in.GetLine();
		
		if (line.Find("Fluid specific volume") >= 0) 
			hd().rho = ScanDouble(line);
		else if (line.Find("Gravity") >= 0) 
			hd().g = ScanDouble(line);
		else if (line.Find("Water depth") >= 0) {
			hd().h = ScanDouble(line);
			if (hd().h < 0)
				throw(Exc(t_("Water depth has to be positive")));
			if (hd().h == 0)
				hd().h = -1;
		} else if (line.Find("Number of bodies") >= 0) 
			hd().Nb = ScanInt(line);
		else if (line.Find("Name of mesh file") >= 0) { 
			f.Load(line);
			hd().names << GetFileTitle(f.GetText(0));
		} else if (line.Find("Number of wave frequencies") >= 0) {
			f.Load(line);
			hd().Nf = f.GetInt(0);  						
			double minF = f.GetDouble(1);
			double maxF = f.GetDouble(2);
        	LinSpaced(hd().w, hd().Nf, minF, maxF); 
        	hd().T.SetCount(hd().Nf);
        	for (int i = 0; i < hd().Nf; ++i)
        		hd().T[i] = 2*M_PI/hd().w[i];  				
		} else if (line.Find("Number of wave directions") >= 0) {
			f.Load(line);
			hd().Nh = f.GetInt(0);  						
			double minD = f.GetDouble(1);
			double maxD = f.GetDouble(2);
        	LinSpaced(hd().head, hd().Nh, minD, maxD); 		
		}
	}
	if (hd().Nb == 0 || hd().Nf == 0 || hd().Nh == 0 || IsNull(hd().rho) || IsNull(hd().g) || IsNull(hd().h))
		throw Exc(Format(t_("Wrong format in Nemoh file '%s'"), fileName));
		
	return true;
}

bool Nemoh::Load_Inf(String fileName) {
	if (hd().Nb != 1)
		throw Exc(Format(t_("SeaFEM_Nemoh only allows one body, found %d"), hd().Nb));
		
	hd().cg.setConstant(3, 1, Null);
	hd().cb.setConstant(3, 1, Null);
	hd().Vo.SetCount(1, Null);
	hd().C.SetCount(1);
	hd().C[0].setConstant(6, 6, Null);   
	
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
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
	}
	return true;	
}
	
bool Nemoh::Load_Hydrostatics() {
	hd().cg.setConstant(3, hd().Nb, Null);
	hd().cb.setConstant(3, hd().Nb, Null);
	hd().Vo.SetCount(hd().Nb, Null);
	String line;
	
	for (int b = 0; b < hd().Nb; ++b) {
	    String fileHydro;
	    if (hd().Nb == 1)
	        fileHydro = AppendFileName(folder, AppendFileName("Mesh", "Hydrostatics.dat"));
	    else
	        fileHydro = AppendFileName(folder, AppendFileName("Mesh", Format("Hydrostatics_%d.dat", b)));
	    
	    FileInLine in(fileHydro);
	    if (!in.IsOpen())
	        return false;
	    
	    FieldSplit f(in);
	    for (int i = 0; i < 3 && !in.IsEof(); ++i) {
			f.Load(in.GetLine());
			hd().cg(i, b) = f.GetDouble(6);
			hd().cb(i, b) = f.GetDouble(2);
	    }
		f.Load(in.GetLine());
	    hd().Vo[b] = f.GetDouble(2); 		
	}
	return true;
}

bool Nemoh::Load_KH() {
	hd().C.SetCount(hd().Nb);
	for (int ib = 0; ib < hd().Nb; ++ib) {
	    String fileKH;
	    if (hd().Nb == 1) {
	        fileKH = AppendFileName(folder, AppendFileName("Mesh", "KH.dat"));
	    } else {
	        fileKH = AppendFileName(folder, AppendFileName("Mesh", Format("KH_%d.dat", ib)));
	    }
	    
		hd().C[ib].setConstant(6, 6, Null);    
	    FileInLine in(fileKH);
	    if (!in.IsOpen()) 
	        return false;
	    
	    FieldSplit f(in);
		for (int i = 0; i < 6 && !in.IsEof(); ++i) {
			f.Load(in.GetLine());
			for (int ifr = 0; ifr < 6; ++ifr)
				hd().C[ib](i, ifr) = f.GetDouble(ifr);
		}
	}
	return true;
}

bool Nemoh::Load_Radiation(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	String line;
	FieldSplit f(in);
	hd().dof.Clear();	 hd().dof.SetCount(hd().Nb, 0);
	in.GetLine();
	while(!in.IsEof()) {
		line = in.GetLine();
	    if (line.Find("Motion of body") >= 0)
	        break;
	    f.Load(line);
	    int ibody = f.GetInt(1) - 1;
	    int ndof = f.GetInt(2);
		hd().dof[ibody] = ndof;    
	}
	hd().A.SetCount(hd().Nf);
	hd().B.SetCount(hd().Nf);
	for (int k = 0; k < hd().Nf; ++k) {
		hd().A[k].setConstant(hd().Nb*6, hd().Nb*6, Null);
		hd().B[k].setConstant(hd().Nb*6, hd().Nb*6, Null);
	}
	int ibodydof = 0;
	for (int ibody = 0; ibody < hd().Nb; ++ibody) {
		for (int idof = 0; idof < hd().dof[ibody]; ++idof) {
			for (int k = 0; k < hd().Nf; ++k) {	
				f.Load(in.GetLine());
				for (int df = 0; df < hd().dof[ibody]; ++df) {		
					hd().A[k](ibodydof, df) = f.GetDouble(1 + 2*df);
	        		hd().B[k](ibodydof, df) = f.GetDouble(2 + 2*df);
				}
	        }
	        ++ibodydof;
	        in.GetLine();
	    }
	}
	return true;
}

bool Nemoh::Load_Excitation(String folder) {	
	return Load_Forces(hd().ex, folder, "ExcitationForce.tec", "Diffraction force");
}

bool Nemoh::Load_Diffraction(String folder) {
	return Load_Forces(hd().sc, folder, "DiffractionForce.tec", "Diffraction force");
}

bool Nemoh::Load_FroudeKrylov(String folder) {
	return Load_Forces(hd().fk, folder, "FKForce.tec", "FKforce");
}

bool Nemoh::Load_Forces(Hydro::Forces &fc, String nfolder, String fileName, String textDelim) {
	hd().Initialize_Forces(fc);
	FileInLine in(AppendFileName(nfolder, AppendFileName("Results", fileName)));
	if (!in.IsOpen())
		return false;
	String line;
	FieldSplit f(in);
	while(!in.IsEof()) {
		line = in.GetLine();
		if (line.Find(textDelim) >= 0)
	        break;
	}
	for (int h = 0; h < hd().Nh; ++h) {
		int ifr = 0;
		while(!in.IsEof()) {
			line = in.GetLine();
			if (line.Find(textDelim) >= 0)
				break;
			f.Load(line);
			int ib = 0, idof = 0, ibdof = 0;
			for (int i = 0; i < hd().Nb*6; ++i) {
				double ma = fc.ma[h](ifr, ibdof) = f.GetDouble(1 + 2*i);	
				double ph = fc.ph[h](ifr, ibdof) = -f.GetDouble(1 + 2*i + 1); //-Phase to follow Wamit
				fc.re[h](ifr, ibdof) = ma*cos(ph); 
				fc.im[h](ifr, ibdof) = ma*sin(ph); 
				idof++;
				ibdof++;
				if (idof >= hd().dof[ib]) {
					idof = 0;
					ib++;
					ibdof = 6*ib;
					if (ib >= hd().Nb)
						break;
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
	FieldSplit f(in);	
	hd().Awinf.setConstant(hd().Nb*6, hd().Nb*6, Null);
	int ibodydof = 0;
	for (int ibody = 0; ibody < hd().Nb; ++ibody) {
		for (int idof = 0; idof < hd().dof[ibody]; ++idof) {
			while(!in.IsEof()) {
				line = in.GetLine();	
				if (line.Find("Zone t=") >= 0) 
					break;
			}
			line = in.GetLine();	
			f.Load(line);
			for (int df = 0; df < hd().dof[ibody]; ++df) 
				hd().Awinf(ibodydof, df) = f.GetDouble(1 + 2*df);
			
			++ibodydof;
		}
	}
	return true;
}

void Nemoh::Save(String file) {
	throw Exc("Option not implemented");
}		

bool Nemoh::LoadDatMesh(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen()) {
		hd().PrintError("\n" + Format(t_("Impossible to open file '%s'"), fileName));
		mh().lastError = Format(t_("Impossible to open file '%s'"), fileName);
		return false;
	}
	mh().file = fileName;
	mh.SetCode(MeshData::NEMOH_DAT);
	
	String line;
	FieldSplit f(in);	
		
	try {
		line = in.GetLine();	
		f.Load(line);
	
		if (f.GetInt(0) != 2)
			throw Exc(t_("Format error in Nemoh .dat mesh file"));
		if (f.GetInt(1) == 1)
			mh().x0z = true;
		
		mh().nodes.Clear();
		mh().panels.Clear();
		
		while(!in.IsEof()) {
			line = in.GetLine();	
			f.Load(line);
			int id = f.GetInt(0);	
			if (id == 0)
				break;
			Point3D &node = mh().nodes.Add();
			node.x = f.GetDouble(1);
			node.y = f.GetDouble(2);
			node.z = f.GetDouble(3);
		}
		while(!in.IsEof()) {
			line = in.GetLine();	
			f.Load(line);
			int id0 = f.GetInt(0);	
			if (id0 == 0)
				break;
			Panel &panel = mh().panels.Add();
			panel.id[0] = id0-1;
			panel.id[1] = f.GetInt(1)-1;	
			panel.id[2] = f.GetInt(2)-1;	
			panel.id[3] = f.GetInt(3)-1;	
		}	
		//if (mh().Check())
		//	throw Exc(t_("Wrong nodes found in Nemoh .dat mesh file"));
		mh().GetLimits();
	} catch (Exc e) {
		hd().PrintError(Format("\n%s: %s", t_("Error"), e));
		mh().lastError = e;
		return false;
	}
	
	return true;
}

