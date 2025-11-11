	// Copyright 2020 - 2025, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#include "functions.h"
#include <STEM4U/SeaWaves.h>


String Wamit::Load(String file, bool isHams, int iperout, Function <bool(String, int)> Status) {
	dt.name = GetFileTitle(file);
	dt.file = file;	
	
	try {
		String folder = GetFileFolder(file);
		String ext = ToLower(GetFileExt(file));
		int iperin = Null;
		String filepot, filefrc, filecfg;
		
		if (!isHams) {
			String config_wam = AFX(folder, "config.wam");
			if (FileExists(config_wam))
				Load_cfg(config_wam, iperin, iperout);
			
			String fnames_wam = AFX(folder, "fnames.wam");
			BEM::Print("\n- " + Format(t_("Wamit file .wam file '%s'"), GetFileName(fnames_wam)));
			if (!Load_wam(fnames_wam, filepot, filefrc, filecfg))
				BEM::Print(S(": ** wam ") + t_("Not found") + "**");
			
			if (!filecfg.IsEmpty()) {
				if (!FileExists(filecfg))
					filecfg = AFX(folder, filecfg);
				if (!Load_cfg(filecfg, iperin, iperout))
					BEM::Print(S(": ** cfg ") + t_("Not found") + "**");
			}
			if (!filefrc.IsEmpty())
				if (!FileExists(filefrc))
					filefrc = AFX(folder, filefrc);			
			if (!filepot.IsEmpty())
				if (!FileExists(filepot))
					filepot = AFX(folder, filepot);				
		}	
		if (ext == ".out") {
			String fileout = ForceExtSafer(file, ".out");
			BEM::Print("\n\n" + Format(t_("Output file '%s'"), GetFileName(fileout)));
			if (!Load_out(fileout, Status)) 
				BEM::Print(S(": ** out ") + t_("Not found") + "**");
			
			if (FileExists(filefrc)) {
				BEM::Print("\n- " + Format(t_("Force Control file .frc file '%s'"), GetFileName(filefrc)));
				try {
					if (!Load_frc(filefrc))
						BEM::Print(S(": ** frc ") + t_("Not found") + "**");
				} catch  (Exc e) {
					BEM::Print(S(": ** frc ") + e + "**");
				}
			}
		} else if (S(".mcn.hdf").Find(ext) >= 0) {
			FindFile ff(AFX(folder, "*.out"));
			if (!ff)
				BEM::Print(S(": ** out associated to input file ") + t_("Not found") + "**");
			else if (!Load_out(ff.GetPath(), Status)) 
				BEM::Print(S(": ** out ") + t_("Not found") + "**");
				
		} else if (S(".1.2.3.3sc.3fk.hst.4.7.8.9.12d.12s.cfg.frc.pot.mmx.wam.6p").Find(ext) >= 0) {
			if (isHams)
				dt.solver = Hydro::HAMS_WAMIT;
			else if (dt.name == "WAMIT_5S")
				dt.solver = Hydro::WADAM_WAMIT;
			else
				dt.solver = Hydro::WAMIT;
	
			if (!isHams) {
				if (filepot.IsEmpty())
					filepot = ForceExtSafer(file, ".pot");
				BEM::Print("\n- " + Format(t_("Potential Control file .pot file '%s'"), GetFileName(filepot)));
				if (!Load_pot(filepot))
					BEM::Print(S(": ** pot ") + t_("Not found") + "**");
		
				if (filefrc.IsEmpty())
					filefrc = ForceExtSafer(file, ".frc");
				BEM::Print("\n- " + Format(t_("Force Control file .frc file '%s'"), GetFileName(filefrc)));
				try {
					if (!Load_frc(filefrc))
						BEM::Print(S(": ** frc ") + t_("Not found") + "**");
				} catch  (Exc e) {
					BEM::Print(S(": ** frc ") + e + "**");
				}
				
				if (dt.msh.size() > 0) {
					String filegdf = AFX(folder, dt.msh[0].dt.fileName);		// To get LEN and g
					BEM::Print("\n- " + Format(t_("Mesh file .gdf file '%s'"), GetFileName(filegdf)));
					if (!Load_gdf(filegdf))
						BEM::Print(S(": ** gdf ") + t_("Not found") + "**");
				}
			}
			if (iperout == 0) {
				if (Bem().opT == 0) {
					BEM::PrintWarning(t_("\nThere is not information about Wamit files format."
								         "\nThe data is assumed to follow Wamit's default criterion, i.e. it includes periods [s]."
								         "\nPlease review if the results are correct"));
					iperout = 1;
				} else
					BEM::PrintWarning(t_("\nThere is not information about Wamit files format."
								         "\nAn attempt will be made to guess whether the data is included with periods, frequencies, or wavelengths."
								         "\nPlease review if the results are correct"));
			}
			String file1;
			if (FileExists(filefrc))
				file1 = ForceExtSafer(filefrc, ".1");
			else
				file1 = ForceExtSafer(file, ".1");
			
			BEM::Print("\n- " + Format(t_("Hydrodynamic coefficients A and B .1 file '%s'"), GetFileName(file1)));
			if (!Load_1(file1, iperout, isHams))
				BEM::Print(S(": ** .1 ") + t_("Not found or empty") + "**");
			
			String file2 = ForceExtSafer(file1, ".2"),
				   file3 = ForceExtSafer(file1, ".3");
				   
			if (ext == ".2")
				;
			else {
				file = file3;
				if (!FileExists(file3) && FileExists(file2))
					file = file2;
			}
			BEM::Print("\n- " + Format(t_("Diffraction exciting %s file '%s'"), GetFileExt(file), GetFileName(file)));
			if (!Load_3(file, iperout))
				BEM::Print(S(": ** .3 ") + t_("Not found or empty") + "**");
			
			String fileHST = ForceExtSafer(file, ".hst");
			BEM::Print("\n- " + Format(t_("Hydrostatic restoring file '%s'"), GetFileName(fileHST)));
			if (!Load_hst(fileHST))
				BEM::Print(S(": ** .hst ") + t_("Not found or empty") + "**");
		
			String fileRAO = ForceExtSafer(file, ".4");
			BEM::Print("\n- " + Format(t_("RAO file '%s'"), GetFileName(fileRAO)));
			if (!Load_4(fileRAO, iperout))
				BEM::Print(S(": ** .4 ") + t_("Not found or empty") + "**");
			
			BEM::Print("\n- " + Format(t_("Mean drift file '%s.7/.8/.9'"), GetFileTitle(file)));
			if (!Load_789(file, iperout))
				BEM::Print(S(": ** .7.8.9 ") + t_("Not found or empty") + "**");
			
			String file12s = ForceExtSafer(file, ".12s");
			BEM::Print("\n- " + Format(t_("Second order sum coefficients .12s file '%s'"), GetFileName(file12s)));
			if (!Load_12(file12s, true, iperout, Status))
				BEM::Print(S(": ** .12s ") + t_("Not found") + "**");
			
			String file12d = ForceExtSafer(file, ".12d");
			BEM::Print("\n- " + Format(t_("Second order mean drift coefficients .12d file '%s'"), GetFileName(file12d)));
			if (!Load_12(file12d, false, iperout, Status))
				BEM::Print(S(": ** .12d ") + t_("Not found") + "**");
		}
		
		String filemmx = ForceExtSafer(file, ".mmx");
		BEM::Print("\n- " + Format(t_("Mesh file .mmx file '%s'"), GetFileName(filemmx)));
		if (!Load_mmx(filemmx))
			BEM::Print(S(": ** mmx ") + t_("Not found") + "**");
					
		String fileSC = ForceExtSafer(file, ".3sc");
		BEM::Print("\n- " + Format(t_("Scattering file '%s'"), GetFileName(fileSC)));
		if (!Load_Scattering(fileSC, iperout))
			BEM::Print(S(": ** 3sc ") + t_("Not found") + "**");
		String fileFK = ForceExtSafer(file, ".3fk");
		BEM::Print("\n- " + Format(t_("Froude-Krylov file '%s'"), GetFileName(fileFK)));
		if (!Load_FK(fileFK, iperout))
			BEM::Print(S(": ** 3fk ") + t_("Not found") + "**");
	#ifdef flagDEBUG
		String file5p = ForceExtSafer(file, ".5p");
		BEM::Print("\n- " + Format(t_("Pressures file .5p file '%s'"), GetFileName(file5p)));
		if (!Load_5p(file5p, iperout, Status))
			BEM::Print(S(": ** 6p ") + t_("Not found") + "**");
	#endif	
		/*String file6p = ForceExtSafer(file, ".6p");
		BEM::Print("\n- " + Format(t_("Pressures file .6p file '%s'"), GetFileName(file6p)));
		if (!Load_6p(file6p, iperout, Status))
			BEM::Print(S(": ** 6p ") + t_("Not found") + "**");*/
		
		if (IsNull(dt.Nh))
			dt.Nh = 0;
		if (IsNull(dt.Nf))
			dt.Nf = 0;

		if (IsNull(dt.Nb)/* || IsNull(dt.Nh) || IsNull(dt.Nf) || dt.Nh == 0 || dt.Nf == 0*/) 
			throw Exc(t_("No data found"));
		
	} catch (Exc e) {
		Status("", -1);
		return e;
	}
	
	if (IsNum(dt.msh[0].dt.c0)) {		// In BEMRosetta, cg and cb are referred to global axis, not to XBODY
		for (int ib = 0; ib < dt.Nb; ++ib) {
			if (!IsNull(dt.msh[ib].dt.cg)) 	
				dt.msh[ib].dt.cg += dt.msh[ib].dt.c0;
			if (!IsNull(dt.msh[ib].dt.cb)) 	
				dt.msh[ib].dt.cb += dt.msh[ib].dt.c0;
		}
	}

	dt.x_w = dt.y_w = 0;
	
	Status("", -1);
	return String();
}

void Wamit::Save(String file, Function <bool(String, int)> Status, bool force_T, int qtfHeading, double heading) const {
	String fileext;
	
	UVector<Point3D> listPoints = GetListPointsTemp(IsLoadedPotsRad() || IsLoadedPotsDif());
	
	if (!IsNull(dt.msh[0].dt.cg)) {
		BEM::Print("\n- " + Format(t_("Force Control file '%s'"), GetFileName(fileext = ForceExt(file, ".frc"))));
		Save_frc2(fileext, false, IsLoadedQTF(true) || IsLoadedQTF(false), listPoints);
	}
	BEM::Print("\n- " + Format(t_("Configurarion file '%s'"), GetFileName(fileext = ForceExt(file, ".cfg"))));
	Save_cfg(fileext, IsLoadedQTF(true) || IsLoadedQTF(false), false, force_T, !listPoints.IsEmpty());
	
	if (dt.Nh > 0 && dt.Nf > 0) {
		BEM::Print("\n- " + Format(t_("Potential Control file '%s'"), GetFileName(fileext = ForceExt(file, ".pot"))));
		Save_pot(fileext, false, false, false, UArray<Body>());
	}
	if (IsLoadedA() && IsLoadedB()) {
		BEM::Print("\n- " + Format(t_("Hydrodynamic coefficients A and B file '%s'"), GetFileName(fileext = ForceExt(file, ".1"))));
		Save_1(fileext, force_T);
	}
	if (IsLoadedFex()) {
		BEM::Print("\n- " + Format(t_("Diffraction exciting file '%s'"), GetFileName(fileext = ForceExt(file, ".3"))));
		Save_3(fileext, force_T);
	}
	if (IsLoadedC()) {
		BEM::Print("\n- " + Format(t_("Hydrostatic restoring file '%s'"), GetFileName(fileext = ForceExt(file, ".hst"))));
		Save_hst(fileext);
	}
	if (IsLoadedRAO()) {
		BEM::Print("\n- " + Format(t_("RAO file '%s'"), GetFileName(fileext = ForceExt(file, ".4"))));
		Save_4(fileext, force_T);
	}
	
	if (IsLoadedMD()) {
		BEM::Print("\n- " + Format(t_("Mean drift file '%s'"), GetFileName(fileext = ForceExt(file, Format(".%d", dt.mdtype)))));
		Save_789(fileext, force_T, true);
	}

	if (IsLoadedQTF(true)) {
		BEM::Print("\n- " + Format(t_("QTF file '%s'"), GetFileName(fileext = ForceExt(file, ".12s"))));
		Save_12(fileext, true, Status, force_T, true, qtfHeading, heading);
	}
	if (IsLoadedQTF(false)) {
		BEM::Print("\n- " + Format(t_("QTF file '%s'"), GetFileName(fileext = ForceExt(file, ".12d"))));
		Save_12(fileext, false, Status, force_T, true, qtfHeading, heading);
	}
}

bool Wamit::Load_out(String fileName, Function <bool(String, int)> Status) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	Status(Format("Loading .out file %s", fileName), 0);
			
	dt.Nf = dt.Nh = 0;
	
	dt.solver = Hydro::WAMIT;
	
	int pos;
	int ibody = -1;
	dt.dimen = false;
	
	String line;
	LineParserWamit f(in);
	f.IsSeparator = IsTabSpace;
	
	bool isHydrostar = false;
	
	UArray<UArray<UArray<VectorXd>>> md7, md8, md9;
	
	auto LoadDrift = [&](UArray<UArray<UArray<VectorXd>>> &md, int ifr) {
		int ih = 0;
		while (!in.IsEof()) {		
			line = in.GetLine();
			if (line.Find("Wave Heading (deg) :") >= 0) {
				f.Load(line);
				double hd1 = f.GetDouble(4);
				double hd2 = f.GetDouble(5);
				int iih = Find(dt.mdhead, std::complex<double>(hd1, hd2));
				if (iih >= 0) {
					in.GetLine(3); 
					while (!TrimBoth(line = in.GetLine()).IsEmpty()) {
						f.Load(line);
						int idf = abs(f.GetInt(0)) - 1;
						int ib = int(idf/6);
						idf -= ib*6;
						if (OUTB(ih, dt.Nh) || OUTB(ifr, dt.Nf) || OUTB(ib, dt.Nb*6))
							throw Exc(in.Str() + "\n"  + Format(t_("Index [%d](%d, %d) out of bounds"), ih, ifr, ib));
						double ma = f.GetDouble(1);
						double ph = ToRad(f.GetDouble(2));
						md[ib][iih][idf][ifr] = (std::polar<double>(ma, ph)).real();	// To get the sign properly
					}
				}
				ih++;
				if (ih >= dt.Nh)
					break;
			}
		}
	};
	
	UArray<std::complex<double>> mdhead;
	
	while(!in.IsEof()) {
		line = in.GetLine();
		f.Load(line);
		
		if (line.Find("HydroStar") > 0 || line.Find("HSrao") > 0) {
			dt.description = "Created with HydroStar";
			isHydrostar = true;
		}
		if ((pos = line.FindAfter("N=")) >= 0) {
			ibody = ScanInt(line.Mid(pos, 10));
			if (IsNull(ibody))
				throw Exc(in.Str() + "\n" +  t_("Wrong body index"));
			dt.Nb = max(dt.Nb, ibody);
			if (ibody > dt.msh.size()) 
				dt.msh.SetCount(dt.Nb);
			ibody--;
			if (dt.msh[ibody].dt.name.IsEmpty() && !isHydrostar)
				dt.msh[ibody].dt.name = GetFileTitle(f.GetText(2));
		} else if ((pos = line.FindAfter("Input from Geometric Data File:")) >= 0) {
			ibody = 0;
			dt.Nb = 1;
			dt.msh.SetCount(dt.Nb);
			dt.msh[ibody].dt.name = GetFileTitle(TrimBoth(line.Mid(pos)));		// 1 body
		} else if (line.Find("Gravity:") >= 0) {
			dt.g = f.GetDouble(1);
			dt.len = f.GetDouble(4);
		} else if (line.Find("Water depth:") >= 0) {
			if (ToLower(f.GetText(2)) == "infinite")
				dt.h = -1;
			else {
				dt.h = f.GetDouble(2);
				if (dt.h < 0)
					throw Exc(in.Str() + "\n" +  t_("Water depth has to be positive"));
			}
			if (line.Find("Water density:") >= 0) 
				dt.rho = f.GetDouble(5);			
		} else if (line.Find("Water density:") >= 0) 
			dt.rho = f.GetDouble(2);			
		else if (line.Find("XBODY =") >= 0) {
			if (IsNull(dt.Nb)) {
				dt.Nb = 1;
				dt.msh.SetCount(dt.Nb);
				ibody = 0;
			}
			if (ibody >= dt.Nb)
				throw Exc(in.Str() + "\n"  + Format(t_("Found additional bodies over %d"), dt.Nb));
			dt.msh[ibody].dt.c0.x = f.GetDouble(2);
			dt.msh[ibody].dt.c0.y = f.GetDouble(5);
			dt.msh[ibody].dt.c0.z = f.GetDouble(8);
		} else if ((pos = line.FindAfter("Volumes (VOLX,VOLY,VOLZ):")) >= 0) {		
			dt.msh[ibody].dt.Vo = ScanDouble(line.Mid(pos));
		} else if (line.Find("Center of Gravity  (Xg,Yg,Zg):") >= 0) {
			dt.msh[ibody].dt.cg.x = f.GetDouble(4);
			dt.msh[ibody].dt.cg.y = f.GetDouble(5);
			dt.msh[ibody].dt.cg.z = f.GetDouble(6);
		} else if (line.Find("Center of Buoyancy (Xb,Yb,Zb):") >= 0) {
			dt.msh[ibody].dt.cb.x = f.GetDouble(4);
			dt.msh[ibody].dt.cb.y = f.GetDouble(5);
			dt.msh[ibody].dt.cb.z = f.GetDouble(6);
		} else if (line.Find("Radii of gyration:") >= 0) {
			if (IsNum(dt.rho) && IsNum(dt.msh[ibody].dt.Vo)) {
				double mass = dt.rho*dt.msh[ibody].dt.Vo;
				if (dt.msh[ibody].dt.M.size() > 0)
					throw Exc(in.Str() + "\n"  + t_("Problem in M matrix. Please report"));
				dt.msh[ibody].dt.M.setConstant(6, 6, 0);
				Eigen::MatrixXd &inertia = dt.msh[ibody].dt.M;
				inertia(0, 0) = inertia(1, 1) = inertia(2, 2) = 1;
				inertia(3, 3) = sqr(f.GetDouble(3));
				inertia(3, 4) = sqr(f.GetDouble(4));
				inertia(3, 5) = sqr(f.GetDouble(5));
				f.GetLine();
				inertia(4, 3) = sqr(f.GetDouble(0));
				inertia(4, 4) = sqr(f.GetDouble(1));
				inertia(4, 5) = sqr(f.GetDouble(2));
				f.GetLine();
				inertia(5, 3) = sqr(f.GetDouble(0));
				inertia(5, 4) = sqr(f.GetDouble(1));
				inertia(5, 5) = sqr(f.GetDouble(2));
				double cx = dt.msh[ibody].dt.cg.x;
				double cy = dt.msh[ibody].dt.cg.y;
				double cz = dt.msh[ibody].dt.cg.z;
				inertia(1, 5) = inertia(5, 1) =  cx;
				inertia(2, 4) = inertia(4, 2) = -cx;
				inertia(2, 3) = inertia(3, 2) =  cy;
				inertia(0, 5) = inertia(5, 0) = -cy;
				inertia(0, 4) = inertia(4, 0) =  cz;
				inertia(1, 3) = inertia(3, 1) = -cz;
				
				inertia *= mass;
			}
		} else if (line.Find("Global body and external mass matrix:") >= 0) {
			if (dt.msh[ibody].dt.M.size() > 0)
				throw Exc(in.Str() + "\n"  + t_("Problem in M matrix. Please report"));
			dt.msh[ibody].dt.M.setConstant(6, 6, NaNDouble);
			for (int r = 0; r < 6; ++r) {
				f.GetLine();
				for (int c = 0; c < 6; ++c)
					dt.msh[ibody].dt.M(r, c) = f.GetDouble(c);
			}		
		} else if (line.Find("Hydrostatic and gravitational") >= 0) {
			//if (dt.C.size() < dt.Nb)
			// 	throw Exc(in.Str() + "\n"  + t_("C matrix is not dimensioned"));
			dt.msh[ibody].dt.C.setConstant(6, 6, 0);
			f.LoadWamitJoinedFields(in.GetLine());
			dt.msh[ibody].dt.C(2, 2) = f.GetDouble(1);
			dt.msh[ibody].dt.C(2, 3) = dt.msh[ibody].dt.C(3, 2) = f.GetDouble(2);
			dt.msh[ibody].dt.C(2, 4) = dt.msh[ibody].dt.C(4, 2) = f.GetDouble(3);
			f.LoadWamitJoinedFields(in.GetLine());
			dt.msh[ibody].dt.C(3, 3) = f.GetDouble(1);
			dt.msh[ibody].dt.C(3, 4) = dt.msh[ibody].dt.C(4, 3) = f.GetDouble(2);
			dt.msh[ibody].dt.C(3, 5) = f.GetDouble(3);
			f.LoadWamitJoinedFields(in.GetLine());
			dt.msh[ibody].dt.C(4, 4) = f.GetDouble(1);
			dt.msh[ibody].dt.C(4, 5) = f.GetDouble(2);
			
		} else if (line.Find("Output from") >= 0) {
			dt.head.Clear();
			FileInLine::Pos fpos = in.GetPos();
			
			bool foundNh = false, found2ndorder = false, found7 = false, found8 = false, found9 = false;
			while (!in.IsEof()) {
				line = in.GetLine();
				if (line.Find("Wave period (sec)") >= 0) {
					++dt.Nf;
					if (dt.head.size() > 0 && !foundNh)
						foundNh = true;
				} else if (!foundNh) {
					if (dt.head.size() > 0 && (line.Find("*********************") >= 0))// ||
								   				 //line.Find("FORCES AND MOMENTS") >= 0)) 
						foundNh = true;
					else if (line.Find("Wave Heading (deg) :") >= 0) {
						f.Load(line);
						if (f.size() == 5)
							FindAddDelta(dt.head, f.GetDouble(4)/*FixHeading_180(f.GetDouble(4))*/, 0.001);
						else if (f.size() == 6)
							FindAdd(mdhead, std::complex<double>(f.GetDouble(4), f.GetDouble(5)));
					}
				} else if (line.Find("2nd-order") >= 0) {
					found2ndorder = true;
					break;
				} else if (line.Find("(Momentum Conservation)") >= 0) 
					found8 = true;
				else if (line.Find("(Pressure Integration)") >= 0) 
					found9 = true;
				else if (line.Find("(Control Surface)") >= 0) 
					found7 = true;
			}
			dt.Nh = dt.head.size();
			if (dt.Nb == 0)
				throw Exc(Format(t_("No bodies found in Wamit file '%s'"), dt.file));
			if (dt.Nf == 0)
				throw Exc(Format(t_("No frequencies found in Wamit file '%s'"), dt.file));
			
			dt.w.SetCount(dt.Nf);

			Initialize_AB(dt.A);
			Initialize_AB(dt.B);
			
			dt.mdhead.resize(mdhead.size());
			::Copy(mdhead, dt.mdhead);
			
			if (found7)
				Hydro::Initialize_MD(md7, dt.Nb, int(dt.mdhead.size()), dt.Nf);
			if (found8)
				Hydro::Initialize_MD(md8, dt.Nb, int(dt.mdhead.size()), dt.Nf);
			if (found9)
				Hydro::Initialize_MD(md9, dt.Nb, int(dt.mdhead.size()), dt.Nf);
						
			int qtfNh = 0, qtfNf;
			UVector<double> qw;
			UArray<std::complex<double>> head;
			if (found2ndorder) {
				bool foundSum = false, foundDif = false;
				while (!in.IsEof()) {
					f.GetLine();
					if (f.IsInLine("Period indices:")) {
						FindAdd(qw, f.GetDouble(f.size()-1));
						FindAdd(qw, f.GetDouble(f.size()-2));
					} else if (f.IsInLine("Headings (deg):")) {
						double hd1 = f.GetDouble(6);
						double hd2 = f.GetDouble(7);
						FindAdd(head, std::complex<double>(hd1, hd2));	
					} else if (!foundSum && f.IsInLine("SUM-FREQUENCY")) 
						foundSum = true;
					else if (!foundDif && f.IsInLine("DIFFERENCE-FREQUENCY")) 
						foundDif = true;
				}
				dt.qhead.resize(qtfNh = head.size());
				for (int i = 0; i < qtfNh; ++i)
					dt.qhead[i] = head[i];
				dt.qw.resize(qtfNf = qw.size());
				for (int i = 0; i < qtfNf; ++i)
					dt.qw[i] = 2*M_PI/qw[i];
					
				if (foundSum)
					Hydro::Initialize_QTF(dt.qtfsum, dt.Nb, qtfNh, qtfNf);
				if (foundDif)
					Hydro::Initialize_QTF(dt.qtfdif, dt.Nb, qtfNh, qtfNf);
			}
			
			in.SeekPos(fpos);
			while (in.GetLine().Find("Wave period = infinite") < 0 && !in.IsEof())
				; 
			if (!in.IsEof()) {
				dt.A0.setConstant(dt.Nb*6, dt.Nb*6, NaNDouble);
				Load_A(in, dt.A0);
			}
			in.SeekPos(fpos);
			while (in.GetLine().Find("Wave period = zero") < 0 && !in.IsEof())
				; 
			if (!in.IsEof()) {
				dt.Ainf.setConstant(dt.Nb*6, dt.Nb*6, NaNDouble);
				Load_A(in, dt.Ainf);
			}
			
			in.SeekPos(fpos);
			
			int ifr = -1;
			while (!in.IsEof()) {
				line = in.GetLine();
				if (line.Find("2nd-order") >= 0) 
					break;
				else if (line.Find("Wave period (sec)") < 0)
					continue;
				
				f.Load(line);
				
				ifr++;
				if (OUTB(ifr, dt.Nf))
					throw Exc(in.Str() + "\n" + Format(t_("Found additional frequencies over %d"), dt.Nf));
				
				int idT = 3;
				if (f.GetText(3) == "=")	// Hydrostar
					idT = 4;
	            dt.w[ifr] = 2*M_PI/f.GetDouble(idT);
	            
	            bool nextFreq = false;
	            while (!in.IsEof() && !nextFreq) {
	            	line = in.GetLine();
	            	if (line.Find("ADDED-MASS AND DAMPING COEFFICIENTS") >= 0) {
						in.GetLine(2);
		            
			            while (!in.IsEof()) {
							line = TrimBoth(in.GetLine());
							if (line.IsEmpty())
			                	break;
							f.Load(line);
							int i = f.GetInt(0) - 1;
							int j = f.GetInt(1) - 1;
							double Aij = f.GetDouble(2);
							double Bij = f.GetDouble(3);
							if (OUTB(i, dt.Nb*6) || OUTB(j, dt.Nb*6))
								throw Exc(in.Str() + "\n"  + Format(t_("Index (%d, %d) out of bounds"), i, j));
							dt.A[i][j][ifr] = Aij;
							dt.B[i][j][ifr] = Bij;
						}
	            	} else if (line.Find("DIFFRACTION EXCITING FORCES AND MOMENTS") >= 0) {
						if (dt.ex.IsEmpty()) 
							Initialize_Forces(dt.ex);
						
						int ih = 0;
						while (!in.IsEof()) {		
							line = in.GetLine();
							if (line.Find("Wave Heading (deg) :") >= 0) {
								f.Load(line);
								double hd = f.GetDouble(4); //FixHeading_180(f.GetDouble(4)); 
								int iih = FindDelta(dt.head, hd, 0.001);
								if (iih >= 0) {
									in.GetLine(3); 
									while (!TrimBoth(line = in.GetLine()).IsEmpty()) {
										f.Load(line);
										double ma = f.GetDouble(1);
										double ph = ToRad(f.GetDouble(2));
										int idf = abs(f.GetInt(0)) - 1;
										if (OUTB(ih, dt.Nh) || OUTB(ifr, dt.Nf) || OUTB(idf, dt.Nb*6))
											throw Exc(in.Str() + "\n"  + Format(t_("Index [%d](%d, %d) out of bounds"), ih, ifr, idf));
										int ib = idf/6;
										idf %= 6;
										dt.ex[ib][ih](ifr, idf) = std::polar(ma, ph);	
									}
								}
								ih++;
								if (ih >= dt.Nh)
									break;
							}
						} 
					} else if (line.Find("RESPONSE AMPLITUDE OPERATORS") >= 0) {
						if (dt.rao.IsEmpty()) 
							Initialize_Forces(dt.rao);
						
						int ih = 0;
						while (!in.IsEof()) {		
							line = in.GetLine();
							if (line.Find("Wave Heading (deg) :") >= 0) {
								f.Load(line);
								double hd = f.GetDouble(4); 	//FixHeading_180(f.GetDouble(4)); 
								int iih = FindDelta(dt.head, hd, 0.001);
								if (iih >= 0) {
									in.GetLine(3); 
									while (!TrimBoth(line = in.GetLine()).IsEmpty()) {
										f.Load(line);
										double ma = f.GetDouble(1);
										double ph = ToRad(f.GetDouble(2));
										int idf = abs(f.GetInt(0)) - 1;
										if (OUTB(ih, dt.Nh) || OUTB(ifr, dt.Nf) || OUTB(idf, dt.Nb*6))
											throw Exc(in.Str() + "\n"  + Format(t_("Index [%d](%d, %d) out of bounds"), ih, ifr, idf));
										int ib = idf/6;
										idf %= 6;
										dt.rao[ib][iih](ifr, idf) = std::polar(ma, ph);	
									}
								}
								ih++;
								if (ih >= dt.Nh)
									break;
							}
						}
					} else if (line.Find("SURGE, SWAY & YAW DRIFT FORCES (Momentum Conservation)") >= 0) 
						LoadDrift(md8, ifr);
					else if (line.Find("SURGE, SWAY, HEAVE, ROLL, PITCH & YAW DRIFT FORCES (Pressure Integration)") >= 0)
						LoadDrift(md9, ifr);
					else if (line.Find("SURGE, SWAY, HEAVE, ROLL, PITCH & YAW DRIFT FORCES (Control Surface)") >= 0)
						LoadDrift(md7, ifr);					
					else if (line.Find("VELOCITY VECTOR IN FLUID DOMAIN") >= 0 ||
							   line.Find("HYDRODYNAMIC PRESSURE IN FLUID DOMAIN") >= 0 ||
							   line.Find("*************************************") >= 0) {
						nextFreq = true;
						break;
					} 
				}
			}
			
			if (Hydro::IsLoadedMD(md7) && md7[0][0][0][md7[0][0][0].size()/2] > 0.) {
				dt.md = pick(md7);
				dt.mdtype = dt.qtftype = 7;
			} else if (Hydro::IsLoadedMD(md9) && md9[0][0][0][md9[0][0][0].size()/2] > 0.) {
				dt.md = pick(md9);
				dt.mdtype = dt.qtftype = 9;
			} else if (Hydro::IsLoadedMD(md8) && md8[0][0][0][md8[0][0][0].size()/2] > 0.) {
				dt.md = pick(md8);
				dt.mdtype = dt.qtftype = 8;
			}
			
			if (found2ndorder) {
				int ifr1, ifr2, ih; 
				UArray<UArray<UArray<MatrixXcd>>> *qtf = nullptr;
				while (!in.IsEof()) {
					f.GetLine();
					if (f.IsInLine("Period indices:")) {
						ifr1 = Find(dt.qw, 2*M_PI/f.GetDouble(f.size()-2));
						ifr2 = Find(dt.qw, 2*M_PI/f.GetDouble(f.size()-1));
						if (ifr1 < 0 || ifr2 < 0)
							throw Exc(in.Str() + "\n"  + t_("Periods not found"));
					} else if (f.IsInLine("SUM-FREQUENCY")) 
						qtf = &dt.qtfsum;
					else if (f.IsInLine("DIFFERENCE-FREQUENCY")) 
						qtf = &dt.qtfdif;
					else if (f.IsInLine("Heading indices:")) {
						double hd1 = f.GetDouble(6);
						double hd2 = f.GetDouble(7);
						ih = Find(head, std::complex<double>(hd1, hd2));
					} else if (f.IsInLine("Mod[F2(I)]")) {					
						bool start = true;									// Because of Hydrostar has not an empty line
						while (!Trim(f.GetLine()).IsEmpty() || start) {	
							if (f.size() == 0) 					
								;
							else {
								start = false;
								int idof = f.GetInt(0)-1;
								int ib = idof/6;
								idof -= ib*6;
								double ma = f.GetDouble(1);
								double ph = ToRad(f.GetDouble(2));
								if (OUTB(ih, qtfNh) || OUTB(ifr1, qtfNf) || OUTB(ifr2, qtfNf) || OUTB(ib, dt.Nb))
									throw Exc(in.Str() + "\n"  + Format(t_("Index [%d][%d](%d, %d) out of bounds"), ib, ih, idof, ifr1, ifr2));
								(*qtf)[ib][ih][idof](ifr1, ifr2) = std::polar(ma, ph);
							}
						}
					}
				}
			}
		}
	}
	
	if (dt.Nb == 0)
		throw Exc(t_("Incorrect .out format"));
		
	if (isHydrostar) {		// HydroStar corrections due to mismatch with Wamit rules
		Status("Processing HydroStar data", 0);
		
		dt.solver = Hydro::HYDROSTAR_OUT;
		
		UVector<Point3D> refPoint;
		UVector<Pointf> refWave;
		if (Wamit::Load_mcn(fileName, dt.Nb, refPoint, refWave)) {
			for (int ib = 0; ib < dt.Nb; ++ib)
				dt.msh[ib].dt.c0 = clone(refPoint[ib]);
			if (IsNull(refWave[0])) {
				for (int ib = 0; ib < dt.Nb; ++ib)
					refWave[ib] = Pointf(dt.msh[ib].dt.cb.x + dt.msh[ib].dt.c0.x, 		// Real cb = cb + c0
										 dt.msh[ib].dt.cb.y + dt.msh[ib].dt.c0.y);
			}
		} else {
			BEM::PrintError(t_("HydroStar .mcn file not found.\nIt is advisable to include it to avoid misunderstandings regarding the global and body axis.\nThe default criteria defined by HydroStar have been considered."));
			refWave.SetCount(dt.Nb);
			for (int ib = 0; ib < dt.Nb; ++ib)
				refWave[ib] = Pointf(dt.msh[ib].dt.cb.x, dt.msh[ib].dt.cb.y);
		}
		if (Load_HDF(Status)) {
			MatrixXd delta(3, dt.Nb);
			for (int ib = 0; ib < dt.Nb; ++ib)
				delta.col(ib) = Vector3d(dt.msh[ib].dt.cb);
			TranslateRadiationPotentials(-delta);// Radiation are referred to the cb, so they are translated to c0 (remember that here cb is not absolute but referred to c0)
		}
			
		for (int ib = 0; ib < dt.Nb; ++ib)	// Translates all bodies phase to 0,0, by translating -refWave
			AddWave(ib, -refWave[ib].x, -refWave[ib].y, dt.g);	
	}
	
	return true;
}

// This is only for HydroStar
bool Wamit::Load_mcn(String fileName, int nb, UVector<Point3D> &refPoint, UVector<Pointf> &refWave) {
	if (nb == 0)
		return false;
	
	String folder = GetFileFolder(fileName);
	FindFile mcn(AFX(folder, "*.mcn"));
	if (!mcn)
		return false;
	
	String mcnFile = mcn.GetPath();
	
	FileInLine in(mcnFile);
	if (!in.IsOpen())
		return false;

	LineParserWamit f(in);
	f.IsSeparator = IsTabSpace;
	
	UVector<Point3D> cog(nb, Null);
	refPoint.SetCount(nb, Null);
	refWave.SetCount(nb, Null);
	
	while(!in.IsEof()) {
		f.GetLine_discard_empty();
		if (f.IsEof())
			break;

		if (f.GetText(0) == "COGPOINT_BODY") {
			int ib = f.GetInt(1);
			if (ib < 1 || ib > nb)
				throw Exc(in.Str() + "\n"  + t_("Wrong body in COGPOINT_BODY"));
			cog[ib-1] = Point3D(f.GetDouble(2), f.GetDouble(3), f.GetDouble(4));
		} else if (f.GetText(0).StartsWith("REFPOINT")) {
			int ib = f.GetInt(1);
			if (ib < 1 || ib > nb)
				throw Exc(in.Str() + "\n"  + t_("Wrong body in REFPOINT"));
			refPoint[ib-1] = Point3D(f.GetDouble(2), f.GetDouble(3), f.GetDouble(4));
		} else if (f.GetText(0).StartsWith("REFWAVE")) 
			refWave.Set(0, Pointf(f.GetDouble(1), f.GetDouble(2)), nb);
	}
	if (IsNull(cog[0]))
		return false;
	
	if (IsNull(refPoint[0]))
		refPoint = clone(cog);
	
	return true;
}
			
void Wamit::Save_A(FileOut &out, Function <double(int, int)> fun, const Eigen::MatrixXd &base, String wavePeriod) const {
	out << 	" ************************************************************************\n\n"
			" Wave period = " << wavePeriod << "\n"
			" ------------------------------------------------------------------------\n\n\n"
			"    ADDED-MASS COEFFICIENTS\n"
			"     I     J         A(I,J)\n\n";
	for (int r = 0; r < dt.Nb*6; ++r) 
		for (int c = 0; c < dt.Nb*6; ++c) {
			double a;
			if (base.rows() > r && base.cols() > c && IsNum(base(r, c))) 
				a = fun(r, c);
			else
				a = 0;
			out << Format("%6>d%6>d  % E\n", r+1, c+1, a);
		}
	out << "\n\n";
}

void Wamit::Save_AB(FileOut &out, int ifr) const {
	out <<	"    ADDED-MASS AND DAMPING COEFFICIENTS\n"
			"     I     J         A(I,J)         B(I,J)\n\n";
	for (int r = 0; r < dt.Nb*6; ++r) 
		for (int c = 0; c < dt.Nb*6; ++c) {
			double a, b;
			if (IsLoadedA(r, c) && IsLoadedB(r, c) && IsNum(dt.A[r][c][ifr]) && IsNum(dt.B[r][c][ifr])) {
				a = A_ndim(ifr, r, c);
				b = B_ndim(ifr, r, c);
			} else
				a = b = 0;
			out << Format("%6>d%6>d  % E  % E\n", r+1, c+1, a, b);
		}
	out << "\n\n\n\n";
}

void Wamit::Save_Forces(FileOut &out, int ifr) const {
	out <<	"    DIFFRACTION EXCITING FORCES AND MOMENTS\n\n";
	for (int ih = 0; ih < dt.Nh; ++ih) {
		out << "  Wave Heading (deg) :      " << dt.head[ih] << "\n\n"
			<< "     I     Mod[Xh(I)]     Pha[Xh(I)]\n\n";
		for (int ib = 0; ib < dt.Nb; ++ib) 
			for (int idof = 0; idof < 6; ++idof) {
				std::complex<double> c;
				if (IsLoadedFex(idof, ih, ib) && IsNum(dt.ex[ib][ih](ifr, idof)))
					c = F_ndim(dt.ex, ih, ifr, idof, ib);
				else
					c = std::complex<double>(0, 0);
				out << Format("%6>d   %E         %6>d\n", idof + 6*ib + 1, abs(c), round(ToDeg(arg(c))));
			}
		out << "\n\n";
	}
}

void Wamit::Save_RAO(FileOut &out, int ifr) const {
	out <<	"    RESPONSE AMPLITUDE OPERATORS\n\n";
	for (int ih = 0; ih < dt.Nh; ++ih) {
		out << "  Wave Heading (deg) :      " << dt.head[ih] << "\n\n"
			<< "     I     Mod[Xh(I)]     Pha[Xh(I)]\n\n";
		for (int ib = 0; ib < dt.Nb; ++ib) 
			for (int idof = 0; idof < 6; ++idof) {
				std::complex<double> c;
				if (IsLoadedRAO(idof, ih, ib) && IsNum(dt.rao[ib][ih](ifr, idof)))
					c = RAO_ndim(dt.rao, ih, ifr, idof, ib);
				else
					c = std::complex<double>(0, 0);
				out << Format(" %7>d   %E   %f\n", idof + 6*ib + 1, abs(c), ToDeg(arg(c)));
			}
		out << "\n\n\n\n";
	}
}

void Wamit::Save_MD(FileOut &out, int ifr) const {
	if (dt.mdtype == 7)
		out << " SURGE, SWAY, HEAVE, ROLL, PITCH & YAW DRIFT FORCES (Control Surface)";
	else if (dt.mdtype == 8)
		out << " SURGE, SWAY & YAW DRIFT FORCES (Momentum Conservation)";
	else if (dt.mdtype == 9)
		out << " SURGE, SWAY, HEAVE, ROLL, PITCH & YAW DRIFT FORCES (Pressure Integration)  PRE";
	else
		throw Exc("Unknown drift type");
	
	out << "\n\n";
	
	for (int ih = 0; ih < dt.Nh; ++ih) {
		out << Format("  Wave Heading (deg) : %s%s\n\n", FDS(dt.head[ih], 12, true), FDS(dt.head[ih], 12, true))
			<< "     I      Mod[F(I)]    Pha[F(I)]\n\n";
		for (int ib = 0; ib < dt.Nb; ++ib) {
			for (int idf = 0; idf < 6; ++idf) {
				if (IsLoadedMD(ib, ih) && IsNum(dt.md[ib][ih][idf][ifr])) {
					const std::complex<double> &c = Md_ndim(idf+ib*6, ih, ifr);
					out << Format(" %7>d   %E         %6>d\n", idf+1+6*ib, abs(c), round(ToDeg(arg(c))));
				}
			}
		}
		out << "\n\n\n\n";
	}
}

void Wamit::Save_out(String file) const {
	String filename = GetFileTitle(file);
	FileOut out(file);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to save '%s'. File already used."), file));

	out << Format(" %s\n\n", String('-', 71))
		<< "                        WAMIT  Version 6.1/6.107S \n\n"
	    << "                        BEMRosetta generated .out format\n\n\n\n"
	    << Format(" %s\n\n\n", String('-', 72))
	    << Format(" %s\n\n\n", String('-', 72))
		<< " Low-order panel method  (ILOWHI=0)\n\n";
	if (dt.Nb == 1)
		out << " Input from Geometric Data File:         " << filename << ".gdf\n"
			<< " Unknown gdf file source\n\n";
	else {
		out << " Input from Geometric Data Files:\n";
		for (int ib = 0; ib < dt.Nb; ++ib)
			out << Format("                               N=  %d     %s%d.gdf\n"
						  " Unknown gdf file source\n\n", ib+1, filename, ib);
	}
	out	<< " Input from Potential Control File:      " << filename << ".pot\n"
		<< " " << filename << ".pot -- file type .gdf, ILOWHI=0, IRR=1\n\n\n"
		<< " POTEN run date and starting time:        01-Jan-2000  --  00:00:00\n"
		<< "   Period       Time           RAD      DIFF  (max iterations)\n";
	if (IsLoadedA0())
		out << "   -1.0000    00:00:00          -1\n";
	if (IsLoadedAinf())
		out << "    0.0000    00:00:00          -1\n";
	VectorXd T = Get_T();
	for (int it = 0; it < T.size(); ++it)
		out << " " << Format("%9.4f", T[it]) << "    00:00:00          -1      -1\n";
	out << "\n"
	   	<< " Gravity:     " << Bem().g
	    << "                Length scale:        " << dt.len << "\n"
		<< " Water depth:        " << (dt.h < 0 ? "infinite" : FDS(dt.h, 9)) << "    "
		<< " Water density:      " << (Bem().rho) << "\n"
		<< " Logarithmic singularity index:              ILOG =     1\n"
		<< " Source formulation index:                   ISOR =     0\n"
		<< " Diffraction/scattering formulation index: ISCATT =     0\n"
		<< " Number of blocks used in linear system:   ISOLVE =     1\n"
		<< " Number of unknowns in linear system:        NEQN =  111\n\n";      
	
	out << " BODY PARAMETERS:\n\n";
	
	for (int ibody = 0; ibody < dt.Nb; ++ibody) {
		if (dt.Nb > 1)
			out << " Body number: N= " << ibody+1 << "   ";
		out	<< " Total panels:  1111    Waterline panels:   11      Symmetries: none\n";
		out	<< " Irregular frequency index: IRR =1\n"; 
		out	<< " Free surface panels:     111\n\n";
		out	<< Format(" XBODY =    %s YBODY =    %s ZBODY =    %s PHIBODY =   0.0\n", 
						FormatWam(dt.msh[ibody].dt.c0.x), 
						FormatWam(dt.msh[ibody].dt.c0.y), 
						FormatWam(dt.msh[ibody].dt.c0.z));
		double Vo = /*dt.Vo.size() > ibody ?*/ dt.msh[ibody].dt.Vo;// : 0;
		out	<< Format(" Volumes (VOLX,VOLY,VOLZ):      %s %s %s\n", 
					FormatWam(Vo), FormatWam(Vo), FormatWam(Vo));
		
		if (!IsNull(dt.msh[ibody].dt.cb)) {
			double cbx = dt.msh[ibody].dt.cb.x - dt.msh[ibody].dt.c0.x;
			double cby = dt.msh[ibody].dt.cb.y - dt.msh[ibody].dt.c0.y;
			double cbz = dt.msh[ibody].dt.cb.z - dt.msh[ibody].dt.c0.z;
			
			out	<< Format(" Center of Buoyancy (Xb,Yb,Zb): %s %s %s\n", 
						FormatWam(cbx), FormatWam(cby), FormatWam(cbz));
		}
		
		if (IsLoadedC()) {
			out	<< " Hydrostatic and gravitational restoring coefficients:\n"; 
			out	<< " C(3,3),C(3,4),C(3,5): " << Format("%s %s %s\n", 
						FormatWam(C_ndim(ibody, 2, 2)), 
						FormatWam(C_ndim(ibody, 2, 3)), 
						FormatWam(C_ndim(ibody, 2, 4)));
			out	<< " C(4,4),C(4,5),C(4,6):               " << Format("%s %s %s\n", 
						FormatWam(C_ndim(ibody, 3, 3)), 
						FormatWam(C_ndim(ibody, 3, 4)), 
						FormatWam(C_ndim(ibody, 3, 5)));
			out	<< "        C(5,5),C(5,6):                             " << Format("%s %s\n", 
						FormatWam(C_ndim(ibody, 4, 4)), 
						FormatWam(C_ndim(ibody, 4, 5)));
		}
		double cgx = 0, cgy = 0, cgz = 0;
		if (!IsNull(dt.msh[ibody].dt.cg)) {
			cgx = dt.msh[ibody].dt.cg.x - dt.msh[ibody].dt.c0.x;
			cgy = dt.msh[ibody].dt.cg.y - dt.msh[ibody].dt.c0.y;
			cgz = dt.msh[ibody].dt.cg.z - dt.msh[ibody].dt.c0.z;
		}
		out	<< Format(" Center of Gravity  (Xg,Yg,Zg): %s %s %s\n\n", 
					FormatWam(cgx), FormatWam(cgy), FormatWam(cgz));
		if (IsLoadedM()) {
			out	<< " Global body and external mass matrix:        ";
			for (int r = 0; r < 6; ++r) {
				out << "\n ";
				for (int c = 0; c < 6; ++c)
					out << Format(" % 10.4E", dt.msh[ibody].dt.M(r, c));
			}
			out	<< "\n";
		}
		out	<< "\n";
	}
	out << "\n\n\n";
	out << " ------------------------------------------------------------------------\n"
    	<< "                            Output from  WAMIT\n"
		<< " ------------------------------------------------------------------------\n"
		<< " FORCE run date and starting time:                22-Dec-2021 -- 11:54:03\n"
		<< " ------------------------------------------------------------------------\n"
		<< " I/O Files:         " << filename << ".frc       " << filename << ".p2f       " << filename << ".out\n"
		<< "  " << filename << ".frc -- file type .gdf, ILOWHI=0, IRR=1\n \n\n"
	;

	if (IsLoadedA0())
		Save_A(out, [&](int idf, int jdf)->double {return A0_ndim(idf, jdf);},   dt.A0,   "infinite");
	if (IsLoadedAinf())
		Save_A(out, [&](int idf, int jdf)->double {return Ainf_ndim(idf, jdf);}, dt.Ainf, "zero");
	
	for (int ifr = 0; ifr < T.size(); ++ifr) {
		out << 	" ************************************************************************\n\n"
				" Wave period (sec) = " << FormatWam(T[ifr]) << "\n"
				" ------------------------------------------------------------------------\n\n\n";
		if (IsLoadedA() && IsLoadedB()) 
			Save_AB(out, ifr);
		if (IsLoadedFex())
			Save_Forces(out, ifr);
		if (IsLoadedRAO())
			Save_RAO(out, ifr);
		if (IsLoadedMD())
			Save_MD(out, ifr);
	}
	UVector<double> heads;
	UVector<int> id1(int(dt.qhead.size())), id2(int(dt.qhead.size()));
	if (IsLoadedQTF(true) || IsLoadedQTF(false)) {
		for (int ih = 0; ih < dt.qhead.size(); ++ih) {
			id1[ih] = 1+FindAdd(heads, dt.qhead[ih].real());
			id2[ih] = 1+FindAdd(heads, dt.qhead[ih].imag());
		}
	}
		
	if (IsLoadedQTF(true)) {
		for (int ih = 0; ih < dt.qhead.size(); ++ih) {
			for (int ifr1 = 0; ifr1 < dt.qw.size(); ++ifr1) {
				for (int ifr2 = ifr1; ifr2 < dt.qw.size(); ++ifr2) {
					double TT = 2*M_PI/(dt.qw[ifr1] + dt.qw[ifr2]);
					double T1 = 2*M_PI/(dt.qw[ifr1]);
					double T2 = 2*M_PI/(dt.qw[ifr2]);
					out << " ************************************************************************\n\n"
						   " 2nd-order period (sec) =  " << Format("%12E", TT) << "\n\n"
						   " Period indices:     " << Format("%2d   %2d", ifr1+1, ifr2+1) << "          Periods:  " 
						   		<< Format("%12E %12E", T1, 	T2) << "\n"
						   " ------------------------------------------------------------------------\n\n\n"
						   "SUM-FREQUENCY EXCITING FORCES AND MOMENTS-DIRECT METHOD\n\n"
						   "  Heading indices:    " << Format("%2d   %2d", id1[ih], id2[ih]) << "  Headings (deg):       " 
						   		<< Format("%4.1f", dt.qhead[ih].real()) << "      " 
						   		<< Format("%4.1f", dt.qhead[ih].imag()) << "\n\n\n"
						   "      I     Mod[F2(I)]     Pha[F2(I)]\n\n";			// Hydrostar only one \n
					for (int ib = 0; ib < dt.Nb; ++ib) {
						for (int idf = 0; idf < 6; ++idf) {
							static int idf12[] = {1, 3, 5, 2, 4, 6};
			        		int iidf = idf12[idf]-1;
							out << Format("    %2d   %12.6E     %10d\n", 1+iidf + 6*ib, 
									F_ndim(abs(dt.qtfsum[ib][ih][iidf](ifr1, ifr2)), iidf),
 									int(ToDeg(arg(dt.qtfsum[ib][ih][iidf](ifr1, ifr2)))));
						}
					}
					out << "\n\n";
				}
			}
		}
	}
	if (IsLoadedQTF(false)) {
		for (int ih = 0; ih < dt.qhead.size(); ++ih) {
			for (int ifr2 = 0; ifr2 < dt.qw.size(); ++ifr2) {
				for (int ifr1 = ifr2; ifr1 < dt.qw.size(); ++ifr1) {
					String sT, units;
					if (ifr1 == ifr2)
						sT = "infinite";
					else {
						units = "(sec) ";
						sT = Format(" %12E", 2*M_PI/(dt.qw[ifr1] - dt.qw[ifr2]));
					}
					double T1 = 2*M_PI/(dt.qw[ifr1]);
					double T2 = 2*M_PI/(dt.qw[ifr2]);
					out << " ************************************************************************\n\n"
						   " 2nd-order period " << units << "= " << sT << "\n\n"
						   " Period indices:     " << Format("%2d   %2d", ifr1+1, ifr2+1) << "          Periods:  " 
						   		<< Format("%12E %12E", T1, T2) << "\n"
						   " ------------------------------------------------------------------------\n\n\n"
						   "DIFFERENCE-FREQUENCY EXCITING FORCES AND MOMENTS-DIRECT METHOD\n\n"
						   "  Heading indices:    " << Format("%2d   %2d", id1[ih], id2[ih]) << "  Headings (deg):       "
						   		<< Format("%4.1f", dt.qhead[ih].real()) << "      " 
						   		<< Format("%4.1f", dt.qhead[ih].imag()) << "\n\n\n"
						   "      I     Mod[F2(I)]     Pha[F2(I)]\n\n";			// Hydrostar only one \n
					for (int ib = 0; ib < dt.Nb; ++ib) {
						for (int idf = 0; idf < 6; ++idf) {
							static int idf12[] = {1, 3, 5, 2, 4, 6};
			        		int iidf = idf12[idf]-1;
							out << Format("    %2d   %12.6E     %10d\n", 1+iidf + 6*ib, 
									F_ndim(abs(dt.qtfdif[ib][ih][iidf](ifr1, ifr2)), iidf),
 									int(ToDeg(arg(dt.qtfdif[ib][ih][iidf](ifr1, ifr2)))));
						}
					}
					out << "\n\n";
				}
			}
		}
	}     									     				
}

void Wamit::Load_A(FileInLine &in, Eigen::MatrixXd &A) {
	in.GetLine(6);
	while (!in.IsEof()) {
		String line = TrimBoth(in.GetLine());
		if (line.IsEmpty())
           	break;
		LineParser f(in);
		f.IsSeparator = IsTabSpace;
		f.Load(line);
		int i = f.GetInt(0) - 1;
		int j = f.GetInt(1) - 1;
		double Aij = f.GetDouble(2);
		if (OUTB(i, A.rows()) || OUTB(j, A.cols()))
			throw Exc(in.Str() + "\n"  + Format(t_("Index (%d, %d) out of bounds"), i, j));
		A(i, j) = Aij;
	}
}

bool Wamit::Load_cfg(String fileName, int &iperin, int &iperout) {
	iperout = 0;	// If cfg exists, but no iperout, then it is 1 by default
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	LineParser f(in);
 	f.IsSeparator = [](int c)->int {
		if (c == '\t' || c == ' ' || c == '=')
			return true;
		return false;
	};
 	iperin = iperout = 1;
 	while (!in.IsEof()) {
		f.Load(in.GetLine());
		if (!f.IsEmpty()) {
			if (f.GetText(0) == "IPERIN") 
				iperin = f.GetInt(1);
			else if (f.GetText(0) == "IPEROUT") 
				iperout = f.GetInt(1);
			else if (f.GetText(0) == "IPERIO") 
				iperin = iperout = f.GetInt(1);
		}
 	}
 	return true;
}

bool Wamit::Load_wam(String fileName, String &filepot, String &filefrc, String &filecfg) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	String folder = GetFileFolder(fileName);
	LineParser f(in);
	
	while (!f.IsEof()) {
		f.GetLine();
		if (f.size() < 1)
			continue;
		String file = f.GetText(0);
		if (!FileExists(file))
			file = AFX(folder, file);
		String ext = GetFileExt(file);
		if (ext == ".pot")
			filepot = file;
		else if (ext == ".frc" || ext == ".gfrc")
			filefrc = file;
		else if (ext == ".cfg")
			filecfg = file;
	}
	return true;
}

bool Wamit::Load_pot(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	LineParser f(in);
 	f.IsSeparator = IsTabSpace;
 	
 	in.GetLine();
 	f.GetLine();
 	dt.h = f.GetDouble(0);
	if (dt.h <= 0)
		dt.h = -1;
	
	in.GetLine();
 	f.GetLine();
 	
 	// Commented, as there is no relationship between the units of input and output time/frequency. 

	int Nf = f.GetInt(0);		
	
 	if (abs(Nf) > 1000)
 		throw Exc(in.Str() + "\n" + Format(t_("Wrong number of periods %s"), Nf));
 	
 	if (Nf != 0) {
 		int nf = Nf;
 		if (nf > 0) {
	 		while (nf > 0 && !f.IsEof()) {
		 		f.GetLine();
		 		nf -= f.GetCount();
	 		}
 		} else {
 			f.GetLine();
 			nf = abs(nf);
 		}
 	}
 	
 	f.GetLine();
 	dt.Nh = f.GetInt(0);
 	if (abs(dt.Nh) > 1000)
 		throw Exc(in.Str() + "\n" + Format(t_("Wrong number of headings %s"), dt.Nh));
 	
 	if (dt.Nh != 0) {
 		f.GetLine();
 	 	if (dt.Nh > 0) {
	 		dt.head.SetCount(dt.Nh);
		 	if (dt.Nh > f.GetCount())
		 		throw Exc(in.Str() + "\n" + Format(t_("Wrong number of headings %d. Found %d"), dt.Nh, f.GetCount()));
		 	for (int i = 0; i < dt.Nh; ++i) {
				dt.head[i] = f.GetDouble(i);
				if (i > 0 && dt.head[i] <= dt.head[i-1])
					throw Exc(in.Str() + "\n" + Format(t_("Wrong heading %f, it should be higher than previous one"), dt.head[i]));	
		 	}
	 	} else {
	 		dt.Nh = -dt.Nh;
	 		double init = f.GetDouble(0);
	 		dt.head.SetCount(dt.Nh);
			dt.head[0] = init;
	 		double delta = f.GetDouble(1);
		 	for (int i = 1; i < dt.Nh; ++i) 
		 		dt.head[i] = init + i*delta;
	 	}
 	}
 	
 	f.GetLine();
 	dt.Nb = f.GetInt(0);
 	if (dt.Nb < 1 || dt.Nb > 100)
 		throw Exc(in.Str() + "\n" + Format(t_("Wrong number of bodies %s"), f.GetText(0)));
	if (dt.msh.IsEmpty())
		dt.msh.SetCount(dt.Nb);
	
	UVector<int> idxs;
	Upp::Index<String> names;
	for (int ib = 0; ib < dt.Nb; ++ib) {
		Body &b = dt.msh[ib];
		
		f.GetLine();
		
		b.dt.fileName = f.GetText(0);
		if (!FileExists(b.dt.fileName))
			b.dt.fileName = AFX(GetFileFolder(fileName), b.dt.fileName);
		
		String ret = Body::Load(b, b.dt.fileName, dt.rho, Bem().g, Null, Null, false, idxs);
		if (!IsEmpty(ret))
			throw Exc(ret);

		b.dt.name = GetFileTitle(f.GetText(0));
		if (names.Find(b.dt.name) >= 0)
			b.dt.name << "_" << FormatInt(ib+1); 
		names << b.dt.name;
				
		f.GetLine();
		b.dt.c0.x = f.GetDouble(0);
		b.dt.c0.y = f.GetDouble(1);
		b.dt.c0.z = f.GetDouble(2);
		if (f.GetCount() >= 4 && f.GetDouble(3) != 0)
			BEM::PrintWarning("\n" + Format(t_("XBODY angle %f cannot be handled by BEMRosetta"), f.GetDouble(3)));
		
		b.dt.mesh.Translate(b.dt.c0.x, b.dt.c0.y, b.dt.c0.z);
		b.AfterLoad(dt.rho, Bem().g, false, false, false, false);
		
		b.dt.cb -= b.dt.c0;		// This is corrected later
		
		in.GetLine();
	}
 	return true;
}

bool Wamit::Load_frc(String fileName) {
	if (ToLower(GetFileExt(fileName) == ".gfrc"))
		Load_frc3(fileName);
	else {
		try {
			if (!Load_frc1(fileName))
				return false;
		} catch(...) {
			if (!Load_frc2(fileName))
				return false;
		}
	}
	return true;
}

bool Wamit::Load_frc1(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	LineParser f(in);
 	f.IsSeparator = IsTabSpace;
 	
 	dt.rho = Bem().rho;		// There is no other source
 	
	in.GetLine(2);
	
	for (int ib = 0; ib < dt.Nb; ++ib) {
		UVector<int> idxs;
		Body &b = dt.msh[ib];
		
		if (b.IsEmpty()) {
			Point3D c0 = b.dt.c0;
			String ret = Body::Load(b, b.dt.fileName, dt.rho, Bem().g, Null, Null, false, idxs);
			if (!IsEmpty(ret))
				throw Exc(ret);
			b.dt.c0 = c0;
			b.dt.mesh.Translate(c0.x, c0.y, c0.z);
			b.AfterLoad(dt.rho, Bem().g, false, false, false, false);
			
			b.dt.cb -= b.dt.c0;		// This is corrected later
		}
		
		f.GetLine();
		double vcg = f.GetDouble(0);
		b.dt.cg.x = b.dt.cb.x;
		b.dt.cg.y = b.dt.cb.y;
		b.dt.cg.z = vcg;
		
		Matrix3d inertia3;
		for (int r = 0; r < 3; ++r) {
			f.GetLine();
			for (int c = 0; c < 3; ++c) 
				inertia3(r, c) = f.GetDouble(c);
		}
		if (inertia3.maxCoeff() > EPS_LEN) {
			Surface::GetInertia66(b.dt.M, inertia3, b.dt.cg, b.dt.c0, false);
			double m = dt.rho*b.dt.Vo;
			b.dt.M *= m;
		} else
			Clear(b.dt.M);
	}
	
	f.GetLine();		// NBETAH
	int nbeta = f.GetInt(0);
	if (nbeta < 0)
		throw Exc(in.Str() + "\n" + Format(t_("NBETAH '%s' has to be higher or equal to zero"), f.GetText(0)));
	if (nbeta > 0)
		f.GetLine();	// BETAH(1) BETAH(2) ... BETAH(NBETAH)
	
	f.GetLine();		// NFIELD
	int numPoints = f.GetInt(0);
	if (numPoints < 0)
		throw Exc(in.Str() + "\n" + Format(t_("NFIELD '%s' has to be higher or equal to zero"), f.GetText(0)));
	
	listPointsTemp.SetCount(numPoints);
	for (int r = 0; r < numPoints; ++r) {
		f.GetLine();
		if (f.size() < 3 && f.IsEof())
			throw Exc(Format(t_("Field points list is incomplete. Read %d from %d"), r, numPoints));	
		listPointsTemp[r].x = f.GetDouble(0);	// XFIELD(1,1) XFIELD(2,1) XFIELD(3,1)
		listPointsTemp[r].y = f.GetDouble(1);
		listPointsTemp[r].z = f.GetDouble(2);
	}
	
	return true;
}

bool Wamit::Load_frc2(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	LineParser f(in);
 	f.IsSeparator = IsTabSpace;

	in.GetLine(2);
	f.GetLine_discard_empty();
	double rho = f.GetDouble(0);
	if (rho <= 0 || rho > 5000)
		throw Exc(in.Str() + "\n" + Format(t_("Wrong density %s"), f.GetText(0)));
	dt.rho = rho;
	
	f.GetLine_discard_empty();
	int mxNb = f.size()/3;
	int Nb;
	for (Nb = 0; Nb < mxNb; ++Nb) {
		if (IsNull(f.GetDouble_nothrow(3*Nb + 0)) || 
			IsNull(f.GetDouble_nothrow(3*Nb + 1)) || 
			IsNull(f.GetDouble_nothrow(3*Nb + 2)))
			break;
	}
		
	if (!IsNull(dt.Nb) && dt.Nb != Nb) 
		throw Exc(in.Str() + "\n" + Format(t_("Wrong number of bodies %d. They should be %d"), Nb, dt.Nb));
	
	dt.Nb = Nb;
	
	if (dt.msh.IsEmpty())
		dt.msh.SetCount(dt.Nb);
		
	for (int ib = 0; ib < Nb; ++ib) {
		dt.msh[ib].dt.cg.x = f.GetDouble(0 + ib*Nb);
		dt.msh[ib].dt.cg.y = f.GetDouble(1 + ib*Nb);
		dt.msh[ib].dt.cg.z = f.GetDouble(2 + ib*Nb);
	}
	
	f.GetLine_discard_empty();
	int imass = f.GetInt(0);
	if (imass < 0 || imass > 1)
		throw Exc(in.Str() + "\n" + Format(t_("Wrong IMASS %d"), imass));
	if (imass == 1) {
		for (int ib = 0; ib < Nb; ++ib) {
			dt.msh[ib].dt.M.resize(6, 6);
			for (int r = 0; r < 6; ++r) {
				f.GetLine_discard_empty();
				for (int c = 0; c < 6; ++c) 	// Discards crossed mass elements
					dt.msh[ib].dt.M(r, c) = f.GetDouble(c + 6*ib);
			}
		}
	}
				
	f.GetLine_discard_empty();
	int idamp = f.GetInt(0);
	if (idamp < 0 || idamp > 1)
		throw Exc(in.Str() + "\n" + Format(t_("Wrong IDAMP %d"), idamp));
	if (idamp == 1) {
		for (int ib = 0; ib < Nb; ++ib) {
			dt.msh[ib].dt.Dlin.resize(6, 6);
			for (int r = 0; r < 6; ++r) {
				f.GetLine_discard_empty();
				for (int c = 0; c < 6; ++c) 	// Discards crossed mass elements
					dt.msh[ib].dt.Dlin(r, c) = f.GetDouble(c + 6*ib);
			}
		}
	}

	f.GetLine_discard_empty();
	int imoor = f.GetInt(0);
	if (imoor < 0 || imoor > 1)
		throw Exc(in.Str() + "\n" + Format(t_("Wrong ISTIF %d"), imoor));
	if (imoor == 1) {
		for (int ib = 0; ib < Nb; ++ib) {
			dt.msh[ib].dt.Cmoor.resize(6, 6);
			for (int r = 0; r < 6; ++r) {
				f.GetLine_discard_empty();
				for (int c = 0; c < 6; ++c) 	// Discards crossed mass elements
					dt.msh[ib].dt.Cmoor(r, c) = f.GetDouble(c + 6*ib);
			}
		}
	}
	
	f.GetLine_discard_empty();	// NBETAH
	int nbeta = f.GetInt(0);
	if (nbeta > 0)
		f.GetLine_discard_empty();	// BETAH(1) BETAH(2) ... BETAH(NBETAH)
	
	f.GetLine_discard_empty();	// NFIELD
	int numPoints = f.GetInt(0);
	if (numPoints < 0)
		throw Exc(in.Str() + "\n" + Format(t_("NFIELD '%s' has to be higher or equal to zero"), f.GetText(0)));
	
	listPointsTemp.SetCount(numPoints);
	for (int r = 0; r < numPoints; ++r) {
		f.GetLine_discard_empty();
		if (f.size() < 3 && f.IsEof())
			throw Exc(Format(t_("Field points list is incomplete. Read %d from %d"), r, numPoints));	
		listPointsTemp[r].x = f.GetDouble(0);	// XFIELD(1,1) XFIELD(2,1) XFIELD(3,1)
		listPointsTemp[r].y = f.GetDouble(1);
		listPointsTemp[r].z = f.GetDouble(2);
	}
	
	return true;
}

bool Wamit::Load_frc3(String fileName) {
	
	
	///// To be implemented

	return true;	
}

bool Wamit::Load_gdf(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	LineParser f(in);
 	f.IsSeparator = IsTabSpace;
 	
 	in.GetLine();
 	f.GetLine();
 	
 	dt.len = f.GetDouble(0);
	dt.g = f.GetDouble(1);

 	return true;
}

bool Wamit::Load_mmx(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	LineParser f(in);
 	f.IsSeparator = IsTabSpace;
 	
 	int ib = Null;
 	int ialtfrc = 0;
 	
 	while (!in.IsEof()) {
		f.GetLine();
		
		if (f.IsEmpty())
			;
		else if (f.GetText(0) == "Gravity:") {
			dt.g = f.GetDouble(1);
			if (f.GetText(2) != "Length" || f.GetText(3) != "scale:") 
				throw Exc(t_("Expected to find 'Length scale:' text in .mmx file"));
			dt.len = f.GetDouble(4);
		} else if (f.GetText(0) == "NBODY") {
			dt.Nb = f.GetInt(2);
			ialtfrc = f.GetInt(5);

			dt.msh.SetCount(dt.Nb);

			if (!IsLoadedM()) {
				for (int i = 0; i < dt.Nb; ++i)
					dt.msh[i].dt.M = MatrixXd::Zero(6, 6);
			}
			if (!IsLoadedCMoor()) {
				for (int i = 0; i < dt.Nb; ++i)
					dt.msh[i].dt.Cmoor = MatrixXd::Zero(6, 6);
			}
			if (!IsLoadedDlin()) {
				for (int i = 0; i < dt.Nb; ++i)
					dt.msh[i].dt.Dlin = MatrixXd::Zero(6, 6);
			}
		} else if (f.GetText(0) == "WAMIT" && (f.GetText(1) == "Ouputs" || f.GetText(1) == "Outputs")) {
			ib = f.GetInt(6) - 1; 
			if (ib >= dt.Nb)
				throw Exc(Format(t_("Unexpected body %d found"), ib+1));
 		} else if (f.GetText(0) == "Volumes") 
 			dt.msh[ib].dt.Vo = Avg(f.GetDouble(2), f.GetDouble(3), f.GetDouble(4));
 		else if (f.GetText(0) == "Center" && f.GetText(2).StartsWith("Buoyancy")) {
 			dt.msh[ib].dt.cb.x = f.GetDouble(4);
 			dt.msh[ib].dt.cb.y = f.GetDouble(5);
 			dt.msh[ib].dt.cb.z = f.GetDouble(6);
 		} else if (f.GetText(0) == "Center" && f.GetText(2).StartsWith("Gravity")) {
 			dt.msh[ib].dt.cg.x = f.GetDouble(4);
 			dt.msh[ib].dt.cg.y = f.GetDouble(5);
 			dt.msh[ib].dt.cg.z = f.GetDouble(6);
 		} else if (f.IsInt(0) && f.IsInt(1)) {
 			int i = f.GetInt(0) - 1;
 			int j = f.GetInt(1) - 1;
 			
 			i = i%6;
 			j = j%6;
 			
 			if (ialtfrc == 1 && !IsNull(dt.rho)) 
 				dt.msh[ib].dt.M(i, j) = f.GetDouble(2)*dt.rho;
 			else if (ialtfrc == 2) {
	 			dt.msh[ib].dt.M(i, j) = f.GetDouble(2);
	 			if (f.size() > 3)
	 				dt.msh[ib].dt.Cmoor(i, j) = f.GetDouble(3);
	 			if (f.size() > 4)
	 				dt.msh[ib].dt.Dlin(i, j) = f.GetDouble(4);
 			}
 		}
 	}
 	return true;
}

static double w_iperout3(double KL, double g, double len) {
	return sqrt(KL*g/len);
}

static double w_iperout4(double nuL, double g, double len, double h) {
	if (h < 0)
		return w_iperout3(nuL, g, len);
	
	double nu = nuL/len;
	return SeaWaves::FrequencyFromWaveNumber(nu, h, g);
}

double RatioSameDelta(const UVector<double> &w) {
	ASSERT(w.size() > 0);
	
	UVector<double> deltas;
	UVector<int> nums;
	for (int i = 1; i < w.size(); ++i) {
		double delta = w[i] - w[i-1];
		int id = FindRatio(deltas, delta, 0.001);
		if (id < 0) {
			deltas << delta;
			nums << 1;
		} else
			nums[id]++;
	}
	return Max(nums)/double(w.size()-1);
}

/*
	1 	Period in seconds 			T
	2 	Frequency 					ω = 2π∕T
	3 	Infinite-depth wavenumber	KL = ω2L∕g
	4 	Finite-depth wavenumber 	νL (νL tanh νH = ω2L∕g)
*/
int Wamit::GuessIperin(const UVector<double> &w) {
	ASSERT(w.size() > 0);
	
	int iperin = 0, iperout = 0;
	double ratio = 0;
	
	double g = Nvl2(dt.g, Bem().g);
	double len = Nvl2(dt.len, Bem().len);
	for (int ipot = 1; ipot <= 4; ++ipot) {
		UVector<double> wot(w.size());
		for (int i = 0; i < w.size(); ++i) {
			if (ipot == 1)
				wot[i] = 2*M_PI/w[i];
			else if (ipot == 2)
				wot[i] = w[i];
			else if (ipot == 3)
				wot[i] = w_iperout3(w[i], g, len);
			else if (ipot == 4)
				wot[i] = w_iperout4(w[i], g, len, dt.h);
		}
		for (int ipin = 1; ipin <= 4; ++ipin) {
			UVector<double> win(w.size());	
			for (int i = 0; i < w.size(); ++i) {
				if (ipin == 1)
					win[i] = 2*M_PI/wot[i];
				else if (ipin == 2)
					win[i] = wot[i];
				else if (ipin == 3)
					win[i] = sqr(wot[i])*len/g;
				else if (ipin == 4)
					win[i] = SeaWaves::WaveNumber_w(wot[i], dt.h, g/len);
			}
			double r = RatioSameDelta(win);
			if (r > ratio) {
				ratio = r;
				iperin = ipin;
				iperout = ipot;
			}
		}
	}
	if (ratio > 0.8)
		return iperout;
	else
		return 1;
}

double Wamit::InputToFreq(double input, int iper, double g, double h, double len) {
	if (iper == 1)
		return 2*M_PI/input;
	else if (iper == 2)
		return input;
	else if (iper == 3) 		// 	Infinite-depth wavenumber
		return w_iperout3(input, g, len);
	else if (iper == 4)
		return w_iperout4(input, g, len, h);
	else if (iper == 5)
		return SeaWaves::FrequencyFromWaveLength(input/len, h, g);
	
	NEVER();
	return Null;
}

void Wamit::ProcessFirstColumn1_3(UVector<double> &w, int iperout) {	
	if (IsNull(iperout) || iperout < 1)
		iperout = 1;

	double len = Nvl2(dt.len, Bem().len);
	if (iperout == 1) {
		for (double &ww : w)
			ww = 2*M_PI/ww;
	} else if (iperout == 2)
		;
	else if (iperout == 3) {		// 	Infinite-depth wavenumber
		double g = Nvl2(dt.g, Bem().g);
		for (double &ww : w)
			ww = w_iperout3(ww, g, len);
	} else if (iperout == 4) {
		double g = Nvl2(dt.g, Bem().g);
		if (IsNull(dt.h))
			throw Exc(t_("Wamit .1 file with finite water depth wavenumber requires .pot file"));
		for (auto &ww : w) 
			ww = w_iperout4(ww, g, len, dt.h);
	} else if (iperout == 5) {
		for (auto &ww : w)
			ww = SeaWaves::FrequencyFromWaveLength(ww/len, dt.h, dt.g);
	} else
		NEVER();
}

void Wamit::ProcessFirstColumnPot(UVector<double> &w, int iperin) {	
	if (IsNull(iperin) || iperin < 1)
		iperin = 1;
	
	double len = Nvl2(dt.len, Bem().len);
	if (iperin == 1) {
		for (auto &ww : w) 
			ww = 2*M_PI/ww;
	} else if (iperin == 2) 
		;
	else if (iperin == 3) {
		double g = Nvl2(dt.g, Bem().g);
		for (auto &ww : w)
			ww = w_iperout3(ww, g, len);
	} else if (iperin == 4) {
		double g = Nvl2(dt.g, Bem().g);
		for (auto &ww : w) 
			ww = w_iperout4(ww, g, len, dt.h);
	} else
		NEVER();
}

bool Wamit::Load_1(String fileName, int &iperout, bool isHams) {
	dt.dimen = false;
	
	if (IsNull(dt.len))
		dt.len = 1;
	
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	LineParser f(in);
 	f.IsSeparator = IsTabSpace;
 	
 	FileInLine::Pos fpos;
 	while (IsNull(ScanDouble(in.GetLine())) && !in.IsEof())
 		fpos = in.GetPos();
	if (in.IsEof())
		return false;
	
	UVector<double> w; 	
    
	in.SeekPos(fpos);
	
	int maxDof = 0;
	bool thereIsA0 = false, thereIsAinf = false; 
	while (!in.IsEof()) {
		f.GetLine();
		if (IsNull(f.GetDouble_nothrow(3))) {
			BEM::PrintWarning(t_("Wrong data found before file end"));
			break;
		}
		
		double freq = f.GetDouble(0);
		if (freq < 0)
			thereIsA0 = true;
		else if (freq == 0)
			thereIsAinf = true;
		else
			FindAdd(w, freq);
		
		int dof = f.GetInt(1);
		if (dof > maxDof)
			maxDof = dof-1;
	}
	
	int Nb = 1 + int(maxDof/6);
	if (!IsNull(dt.Nb) && dt.Nb < Nb)
		throw Exc(in.Str() + "\n"  + Format(t_("The files read have different number of bodies.\nIn the previous is %d, in the .1 is %d"), dt.Nb, Nb));
	dt.Nb = Nb;
	if (dt.msh.IsEmpty())
		dt.msh.SetCount(dt.Nb);	
	
	int Nf = w.size();
	if (!IsNull(dt.Nf) && dt.Nf != Nf)
		throw Exc(in.Str() + "\n"  + Format(t_("The files read have different number of frequencies.\nIn the previous is %d, in the .1 is %d"), dt.Nf, Nf));
	dt.Nf = Nf;
	
	if (dt.Nb == 0)// || dt.Nf < 2)
		throw Exc(in.Str() + "\n"  + Format(t_("Wrong format in Wamit file '%s'"), dt.file));
	
	UVector<double> src = clone(w);
	
	if (iperout == 0)
		iperout = GuessIperin(w);
	
	ProcessFirstColumn1_3(w, iperout);
	
	if (thereIsA0)
		dt.A0.setConstant(dt.Nb*6, dt.Nb*6, NaNDouble);
	if (thereIsAinf)
		dt.Ainf.setConstant(dt.Nb*6, dt.Nb*6, NaNDouble);

	if (Nf > 0) {
		Initialize_AB(dt.A);
		Initialize_AB(dt.B);
	}
	
	if (!dt.w.IsEmpty()) {
		UVector<double> ww = clone(w);		Sort(ww);
		UVector<double> dw = clone(dt.w);	Sort(dw);
		if (!CompareRatio(ww, dw, 0.001))
			throw Exc(in.Str() + "\n"  + Format(t_("The files read have different number of frequencies.\nIn the previous is %s, in the .1 is %s"), ToString(dt.w), ToString(w)));
	}
	dt.w = pick(w);
		
	in.SeekPos(fpos);
	
	while (!in.IsEof()) {
		f.GetLine();
		if (IsNull(f.GetDouble_nothrow(3)))
			break;
				
		double freq = f.GetDouble(0);
 		int i = f.GetInt(1) - 1;
 		int j = f.GetInt(2) - 1;
 		if (i >= Nb*6 || i < 0 || j >= Nb*6 || j < 0)
			throw Exc(in.Str() + "\n"  + Format(t_("DOF # does not match (%d, %d)"), i+1, j+1));
 		
 		double Aij = f.GetDouble(3);
 		
 		if (freq < 0) {
 			if (!isHams)
				dt.A0(i, j) = Aij;
 			else {
 				if (iperout == 1 || iperout == 5)
 					dt.A0(i, j) = Aij;	
 				else
 					dt.Ainf(i, j) = Aij;	
 			}
		} else if (freq == 0) {
			if (!isHams)
				dt.Ainf(i, j) = Aij;
			else {
 				if (iperout == 1 || iperout == 5)
 					dt.Ainf(i, j) = Aij;	
 				else
 					dt.A0(i, j) = Aij;	
 			}	
		} else {
			int ifr = FindRatio(src, freq, 0.001);
			if (ifr < 0)
				throw Exc(in.Str() + "\n"  + Format(t_("Period/Frequency %f is unknown"), freq));

		  	dt.A[i][j][ifr] = Aij;    
		  	dt.B[i][j][ifr] = f.GetDouble(4);   	
		}
	}
	
	return true;	
}

bool Wamit::Load_hst(String fileName) {
	dt.dimen = false;
	
	if (IsNull(dt.len))
		dt.len = 1;
	
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
 
 	FileInLine::Pos fpos;
 	while (IsNull(ScanDouble(in.GetLine())) && !in.IsEof())
 		fpos = in.GetPos();
	if (in.IsEof())
		return false;
	
	in.SeekPos(fpos);
	
	int maxDof = 0;	
	while (!in.IsEof()) {
		f.Load(in.GetLine());

		int dof = f.GetInt(0);
		if (dof > maxDof)
			maxDof = dof-1;
	}
	
	in.SeekPos(fpos);
	
	int Nb = 1 + int(maxDof/6);
	if (!IsNull(dt.Nb) && dt.Nb < Nb)
		throw Exc(in.Str() + "\n"  + Format(t_("The files read have different number of bodies.\nIn the previous is %d, in the .hst is %d"), dt.Nb, Nb));
	dt.Nb = Nb;
	if (dt.msh.IsEmpty())
		dt.msh.SetCount(dt.Nb);
	
	for(int ib = 0; ib < dt.Nb; ++ib)
		dt.msh[ib].dt.C.setConstant(6, 6, 0);

	while (!in.IsEof()) {
		f.Load(in.GetLine());	
		int i = f.GetInt(0) - 1;
		int ib_i = i/6;
		i = i - ib_i*6;
		int j = f.GetInt(1) - 1;
		int ib_j = j/6;
		j = j - ib_j*6;
		if (ib_i == ib_j) 
			dt.msh[ib_i].dt.C(i, j) = f.GetDouble(2);
	}
		
	return true;
}

bool Wamit::Load_3(String fileName, int &iperout) {
	return Load_Forces(fileName, dt.ex, iperout);
}

bool Wamit::Load_Scattering(String fileName, int &iperout) {
	return Load_Forces(fileName, dt.sc, iperout);
}
		
bool Wamit::Load_FK(String fileName, int &iperout) {
	return Load_Forces(fileName, dt.fk, iperout);
}

bool Wamit::Load_4(String fileName, int &iperout) {
	return Load_Forces(fileName, dt.rao, iperout, dt.solver == Hydro::HAMS_WAMIT);
}

bool Wamit::Load_Forces(String fileName, Hydro::Forces &force, int &iperout, bool israohams) {
	dt.dimen = false;
	
	if (IsNull(dt.len))
		dt.len = 1;
	
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
 
 	FileInLine::Pos fpos;
 	while (IsNull(ScanDouble(in.GetLine())) && !in.IsEof())
 		fpos = in.GetPos();
	if (in.IsEof()) 
		return false;
	
	UVector<double> w, T; 	
    
	in.SeekPos(fpos);
	
	int maxDof = 0;
	dt.head.Clear();
	while (!in.IsEof()) {
		f.GetLine();
		if (IsNull(f.GetDouble_nothrow(3))) {
			BEM::PrintWarning(t_("Wrong data found before file end"));
			break;
		}
		
		double freq = f.GetDouble(0);
		double head = f.GetDouble(1);
		FindAdd(w, freq);
		FindAdd(dt.head, head);
		
		int dof = f.GetInt(2);
		if (dof > maxDof)
			maxDof = dof-1;
	}
	
	if (dt.head.size() == 0)
		throw Exc(in.Str() + "\n"  + Format(t_("Wrong format in Wamit file '%s'"), dt.file));
	
	if (!IsNull(dt.Nh) && dt.Nh != dt.head.size())
		throw Exc(in.Str() + "\n"  + Format(t_("The files read have different number of headings.\nIn the previous is %d, in this one is %d"), dt.Nh, dt.head.size()));
	dt.Nh = dt.head.size();
	
	int Nb = 1 + int(maxDof/6);
	if (!IsNull(dt.Nb) && dt.Nb < Nb)
		throw Exc(in.Str() + "\n"  + Format(t_("The files read have different number of bodies.\nIn the previous is %d, in this one is %d"), dt.Nb, Nb));
	dt.Nb = Nb;
	if (dt.msh.IsEmpty())
		dt.msh.SetCount(dt.Nb);
		
	int Nf = w.size();
	if (!IsNull(dt.Nf) && dt.Nf != Nf)
		throw Exc(in.Str() + "\n"  + Format(t_("The files read have different number of frequencies.\nIn the previous is %d, in this one is %d"), dt.Nf, Nf));
	dt.Nf = Nf;
	
	if (dt.Nb == 0 || dt.Nf < 1)
		throw Exc(in.Str() + "\n"  + Format(t_("Wrong format in Wamit file '%s'"), dt.file));
	
	UVector<double> src = clone(w);
	
	if (iperout == 0)
		iperout = GuessIperin(w);
		
	ProcessFirstColumn1_3(w, iperout);
	bool swap = false;
	if (!dt.w.IsEmpty()) {
		UVector<double> ww = clone(w);		Sort(ww);
		UVector<double> dw = clone(dt.w);	Sort(dw);
		if (!CompareRatio(ww, dw, 0.001))
			throw Exc(in.Str() + "\n"  + Format(t_("The files read have different frequencies.\nIn the previous has %s,\nin this one has %s"), ToString(dt.w), ToString(w)));
		
		if (!EqualRatio(dt.w[0], w[0], 0.001))
			swap = true;	// Swap forces to the same order of original w
	} else
		dt.w = pick(w);
	
	in.SeekPos(fpos);
	
	Initialize_Forces(force);
		
	while (!in.IsEof()) {
		f.GetLine();
		if (IsNull(f.GetDouble_nothrow(3))) 
			break;
		
		double freq = f.GetDouble(0);
		int ifr = FindRatio(src, freq, 0.001);
		if (ifr < 0)
			throw Exc(in.Str() + "\n"  + Format(t_("Period/Frequency %f is unknown"), freq));

		double head = f.GetDouble(1);	
		int ih = FindRatio(dt.head, head, 0.001);
		if (ih < 0)
			throw Exc(in.Str() + "\n"  + Format(t_("Heading %f is unknown"), head));
			
		int idof = f.GetInt(2) - 1;		
		int ib = idof/6;
		idof %= 6;
		double re = f.GetDouble(5);
		double im = f.GetDouble(6);
        if (israohams)
            force[ib][ih](ifr, idof) = std::complex<double>(im, -re)*g_rho_ndim()*pow(dt.len, GetK_F(idof));
        else
            force[ib][ih](ifr, idof) = std::complex<double>(re, im);
	}
	if (swap)
		SwapForcesFreqOrder(force);
	
	return true;
}

bool Wamit::Load_12(String fileName, bool isSum, int iperout, Function <bool(String, int)> Status) {
	dt.dimen = false;
	
	if (IsNull(dt.len))
		dt.len = 1;
	
	String ext = isSum ? ".12s" : ".12d";
	
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
	if (IsNull(dt.Nb))
		dt.Nb = Nb;
	else {
		if (dt.Nb < Nb)
			throw Exc(Format(t_("The files read have different number of bodies.\nIn the previous is %d, in the 12 is %d"), dt.Nb, Nb));
	}
	
	if (dt.msh.IsEmpty())
		dt.msh.SetCount(dt.Nb);
		
	int Nf = w.size();
	int Nh = head.size();
	if (Nh == 0)
		throw Exc(Format(t_("Wrong format in Wamit file '%s'. No headings found"), dt.file));
	
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
		qtf[ib][ih][idf](ifr1, ifr2).real(f.GetDouble(7));
		qtf[ib][ih][idf](ifr1, ifr2).imag(f.GetDouble(8));
	}
	
	if (iperout == 0)
		iperout = GuessIperin(w);
		
	ProcessFirstColumn1_3(w, iperout);
	
	::Copy(w, dt.qw);
	::Copy(head, dt.qhead);
	
	return true;
}

bool Wamit::Load_789(String fileName, int &iperout) {
	UArray<UArray<UArray<VectorXd>>> md;
	if (Load_789_0(fileName, 7, md, iperout)) {
		if (Hydro::IsLoadedMD(md) && md[0][0][0](md[0][0][0].size()/2) != 0.) {
			dt.md = pick(md);
			return true;
		}
	}
	if (Load_789_0(fileName, 9, md, iperout)) {
		if (Hydro::IsLoadedMD(md) && md[0][0][0](md[0][0][0].size()/2) != 0.) {
			dt.md = pick(md);
			return true;
		}
	}
	if (Load_789_0(fileName, 8, md, iperout)) {
		if (Hydro::IsLoadedMD(md) && md[0][0][0](md[0][0][0].size()/2) != 0.) {
			dt.md = pick(md);
			return true;
		}
	}
	return false;
}

bool Wamit::Load_789_0(String fileName, int type, UArray<UArray<UArray<VectorXd>>> &mmd, int &iperout) {
	dt.dimen = false;
	
	if (IsNull(dt.len))
		dt.len = 1;
	
	fileName = ForceExtSafer(fileName, Format(".%d", type));
	
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	dt.mdtype = dt.qtftype = type;
	
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
	
	FileInLine::Pos fpos = in.GetPos();
	f.GetLine();		
	if (!IsNull(f.GetDouble_nothrow(0)))
		in.SeekPos(fpos);		// No header, rewind
	else
		fpos = in.GetPos();		// Avoid header
	
	UVector<double> w, T;
    UArray<std::complex<double>> head;
	
	int Nb = 0;
	while (!in.IsEof()) {
		f.GetLine();
		double freq = f.GetDouble(0);
		FindAdd(w, freq);
		double hd1 = f.GetDouble(1);
		double hd2 = f.GetDouble(2);
		FindAdd(head, std::complex<double>(hd1, hd2));	
		Nb = max(Nb, 1 + (f.GetInt(3)-1)/6);
	}
	
	if (IsNull(dt.Nb))
		dt.Nb = Nb;
	else {
		if (dt.Nb < Nb)
			throw Exc(Format(t_("The files read have different number of bodies.\nIn the previous is %d, in this one is %d"), dt.Nb, Nb));
	}
	
	if (dt.msh.IsEmpty())
		dt.msh.SetCount(dt.Nb);
		
	int Nf = w.size();
	if (!IsNull(dt.Nf) && dt.Nf != Nf)
		throw Exc(in.Str() + "\n"  + Format(t_("The files read have different number of frequencies.\nIn the previous is %d, in this one is %d"), dt.Nf, Nf));
	dt.Nf = Nf;
	
	int Nh = head.size();
	if (Nh == 0)
		throw Exc(Format(t_("Wrong format in Wamit file '%s'. No headings found"), dt.file));
	
	UVector<double> src = clone(w);
	
	if (iperout == 0)
		iperout = GuessIperin(w);
		
	ProcessFirstColumn1_3(w, iperout);
	
	Hydro::Initialize_MD(mmd, Nb, Nh, Nf);
	
	if (!dt.w.IsEmpty()) {
		UVector<double> ww = clone(w);		Sort(ww);
		UVector<double> dw = clone(dt.w);	Sort(dw);
		if (!CompareRatio(w, dw, 0.001))
			throw Exc(in.Str() + "\n"  + Format(t_("The files read have different number of frequencies.\nIn the previous is %s,\nin this one is %s"), ToString(dt.w), ToString(w)));
	}
	dt.w = pick(w);
	
	in.SeekPos(fpos);
	
	while (!in.IsEof()) {
		f.GetLine();
		
		double freq = f.GetDouble(0);
		int ifr1 = FindRatio(src, freq, 0.001);
		if (ifr1 < 0)
			throw Exc(in.Str() + "\n"  + Format(t_("Period/Frequency %f is unknown"), freq));

		int ih = Find(head, std::complex<double>(f.GetDouble(1), f.GetDouble(2)));
		int idf = f.GetInt(3)-1;
		if (idf >= 0) {			// idf may be negative!
			int ib = 0;
			while (idf > 5) {
				idf -= 6;
				ib++;
			}
			mmd[ib][ih][idf][ifr1] = f.GetDouble(6);
		}
	}	

	::Copy(head, dt.mdhead);
	
	return true;
}

bool Wamit::Load_5p(String fileName, int &iperout, Function <bool(String, int)> Status) {
	dt.dimen = false;
	
	if (IsNull(dt.len))
		dt.len = 1;
	
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	Status("Loading .5p base data", 0);
	
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
	
	FileInLine::Pos fpos = in.GetPos();
	f.GetLine();		
	if (!IsNull(f.GetDouble_nothrow(0)))
		in.SeekPos(fpos);		// No header, rewind
	else
		fpos = in.GetPos();		// Avoid header

	int np = listPointsTemp.size();
	UVector<bool> idPanels;
	MatrixXi idPanelsM;
	UVector<int> idFs, idRest;
	LoadListPointsTemp(idPanels, idPanelsM, idFs, idRest);
	
	if (!idPanels.IsEmpty() && Max(idPanels) > 0) {
		Initialize_PotsRad();
		Initialize_PotsIncDiff(dt.pots_dif);
	}
	if (!idFs.IsEmpty() && Max(idFs) >= 0) {
		dt.fsPoints.pots_rad.resize(dt.fsPoints.points.size(), 6, dt.Nf);			dt.fsPoints.pots_rad.setConstant(0);
		dt.fsPoints.pots_dif.resize(dt.fsPoints.points.size(), dt.Nh, dt.Nf);		dt.fsPoints.pots_dif.setConstant(0);
	}
	if (!idRest.IsEmpty() && Max(idRest) >= 0) {
		dt.freePoints.pots_rad.resize(dt.fsPoints.points.size(), 6, dt.Nf);			dt.freePoints.pots_rad.setConstant(0);
		dt.freePoints.pots_dif.resize(dt.fsPoints.points.size(), dt.Nh, dt.Nf);		dt.freePoints.pots_dif.setConstant(0);
	}

	double len = Nvl2(dt.len, Bem().len);
	double g = Nvl2(dt.g, Bem().g);
	
	int lifr = -1;
	
	bool isdiff = false, israd = false;
	double rho_g = g_ndim()*rho_ndim();
	
	while (!f.IsEof()) {
		f.GetLine_discard_empty();
		double freq = Wamit::InputToFreq(f.GetDouble(0), iperout, g, dt.h, len);
		int ifr = FindDelta(dt.w, freq, 0.001);
		
		if (ifr < 0 || ifr == 0)	// Inf and zero frec. Not processed for now
			continue;

		if (Status && ifr > lifr && !Status("Loading .5p base data", (100.*ifr)/dt.Nf))
			throw Exc(t_("Stop by user"));

		double rho_w = dt.w[ifr]*rho_ndim();
		
		lifr = ifr;
		
		if (f.size() == 6) {	// pots_diff
			isdiff = true;
			double head = f.GetDouble(1);	
			int ih = FindDelta(dt.head, head, 0.01);
			if (ih < 0)
				throw Exc(in.Str() + "\n"  + Format(t_("Heading %f is unknown"), head));
	   	
	   		int M = f.GetInt(2);		// Quadrant of symmetry. Unused
	   		
	   		int ip = f.GetInt(3);
			if (ip < 1 || ip > np)
				throw Exc(in.Str() + "\n"  + Format(t_("Point number %d is unknown"), ip));
			ip--;
			
			double re, im;
			re = f.GetDouble(4);
			im = f.GetDouble(5);
			
			std::complex<double> val(-im, re); 	// p = -iρωΦ ; Φ = [Im(p) - iRe(p)]/ρω
			val *= g_ndim()/dt.w[ifr];
			
			if (idPanels[ip]) {
				for (int ibb = 0; ibb < dt.Nb; ++ibb) {
					UArray<UArray<UArray<std::complex<double>>>> &pib = dt.pots_dif[ibb];
					int ipp = idPanelsM(ip, ibb);
					if (ipp >= 0) {
						pib[ipp][ih][ifr] = val;
						break;
					}
				}
			} else if (idFs[ip] >= 0)
				dt.fsPoints.pots_dif(idFs[ip], ih, ifr) = val;
			else
				dt.freePoints.pots_dif(idRest[ip], ih, ifr) = val;
		} else {								// pots_rad
			
			int M = f.GetInt(2);		// Quadrant of symmetry. Unused
			
			israd = true;
			int ip = f.GetInt(1);
			if (ip < 1 || ip > np)
				throw Exc(in.Str() + "\n"  + Format(t_("Point number %d is unknown"), ip));
			ip--;
			
			if (idPanels[ip]) {
				for (int ibb = 0; ibb < dt.Nb; ++ibb) {
					UArray<UArray<UArray<std::complex<double>>>> &pib = dt.pots_rad[ibb];
					int ipp = idPanelsM(ip, ibb);
					if (ipp >= 0) {
						for (int ib = 0; ib < dt.Nb; ++ib) {		
							int col0 = 0; 
							if (ib == 0)
								col0 = 1;
							
							for (int idof = 0; idof < 6; ++idof)
								pib[ipp][idof][ifr] += std::complex<double>(f.GetDouble(col0 + idof*2), f.GetDouble(col0 + idof*2 + 1));
							
							if (ib < dt.Nb-1)
								f.GetLine_discard_empty();
						}
						break;
					}
				}
			} else if (idFs[ip] >= 0) {
				for (int ib = 0; ib < dt.Nb; ++ib) {		
					int col0 = 0; 
					if (ib == 0)
						col0 = 1;
					
					for (int idof = 0; idof < 6; ++idof)
						dt.fsPoints.pots_rad(idFs[ip], idof, ifr) += std::complex<double>(f.GetDouble(col0 + idof*2), f.GetDouble(col0 + idof*2 + 1));
					
					if (ib < dt.Nb-1)
						f.GetLine_discard_empty();
				}
			} else {
				for (int ib = 0; ib < dt.Nb; ++ib) {		
					int col0 = 0; 
					if (ib == 0)
						col0 = 1;
					
					for (int idof = 0; idof < 6; ++idof)
						dt.freePoints.pots_rad(idRest[ip], idof, ifr) += std::complex<double>(f.GetDouble(col0 + idof*2), f.GetDouble(col0 + idof*2 + 1));
					
					if (ib < dt.Nb-1)
						f.GetLine_discard_empty();
				}
			}
		} 
	}
	if (!israd) {
		dt.pots_rad.Clear();
		dt.fsPoints.pots_rad.resize(0, 0, 0);
		dt.freePoints.pots_rad.resize(0, 0, 0);
	}
	if (!isdiff) {
		dt.pots_dif.Clear();
		dt.fsPoints.pots_dif.resize(0, 0, 0);
		dt.freePoints.pots_dif.resize(0, 0, 0);
	}

	return true;
}

bool Wamit::Load_6p(String fileName, int &iperout, Function <bool(String, int)> Status) {
	dt.dimen = false;
	
	if (IsNull(dt.len))
		dt.len = 1;
	
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	Status("Loading .6p base data", 0);
	
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
	
	FileInLine::Pos fpos = in.GetPos();
	f.GetLine();		
	if (!IsNull(f.GetDouble_nothrow(0)))
		in.SeekPos(fpos);		// No header, rewind
	else
		fpos = in.GetPos();		// Avoid header

	int np = listPointsTemp.size();
	UVector<bool> idPanels;
	MatrixXi idPanelsM;
	UVector<int> idFs, idRest;
	LoadListPointsTemp(idPanels, idPanelsM, idFs, idRest);
	
	if (!idPanels.IsEmpty() && Max(idPanels) > 0) {
		Initialize_PotsRad();
		Initialize_PotsIncDiff(dt.pots_dif);
	}
	if (!idFs.IsEmpty() && Max(idFs) >= 0) {
		dt.fsPoints.pots_rad.resize(dt.fsPoints.points.size(), 6, dt.Nf);			dt.fsPoints.pots_rad.setConstant(0);
		dt.fsPoints.pots_dif.resize(dt.fsPoints.points.size(), dt.Nh, dt.Nf);		dt.fsPoints.pots_dif.setConstant(0);
	}
	if (!idRest.IsEmpty() && Max(idRest) >= 0) {
		dt.freePoints.pots_rad.resize(dt.fsPoints.points.size(), 6, dt.Nf);			dt.freePoints.pots_rad.setConstant(0);
		dt.freePoints.pots_dif.resize(dt.fsPoints.points.size(), dt.Nh, dt.Nf);		dt.freePoints.pots_dif.setConstant(0);
	}

	double len = Nvl2(dt.len, Bem().len);
	double g = Nvl2(dt.g, Bem().g);
	
	int lifr = -1;
	
	bool isdiff = false, israd = false;
	double rho_g = g_ndim()*rho_ndim();
	
	while (!f.IsEof()) {
		f.GetLine_discard_empty();
		double freq = Wamit::InputToFreq(f.GetDouble(0), iperout, g, dt.h, len);
		int ifr = FindDelta(dt.w, freq, 0.001);
		
		if (ifr < 0 || ifr == 0)	// Inf and zero frec. Not processed for now
			continue;

		if (Status && ifr > lifr && !Status("Loading .6p base data", (100.*ifr)/dt.Nf))
			throw Exc(t_("Stop by user"));

		double rho_w = dt.w[ifr]*rho_ndim();
		
		lifr = ifr;
		
		if (f.size() == 5 || f.size() == 7) {	// pots_diff
			isdiff = true;
			double head = f.GetDouble(1);	
			int ih = FindDelta(dt.head, head, 0.01);
			if (ih < 0)
				throw Exc(in.Str() + "\n"  + Format(t_("Heading %f is unknown"), head));
	   	
	   		int ip = f.GetInt(2);
			if (ip < 1 || ip > np)
				throw Exc(in.Str() + "\n"  + Format(t_("Point number %d is unknown"), ip));
			ip--;
			
			double re, im;
			if (f.size() == 5) {
				re = f.GetDouble(3);
				im = f.GetDouble(4);
			} else {
				re = f.GetDouble(5);
				im = f.GetDouble(6);
			}
			
			std::complex<double> val(-im*70/dt.w[ifr]/dt.w[ifr], re*110/dt.w[ifr]); 	// p = -iρωΦ ; Φ = [Im(p) - iRe(p)]/ρω
			
			if (idPanels[ip]) {
				for (int ibb = 0; ibb < dt.Nb; ++ibb) {
					UArray<UArray<UArray<std::complex<double>>>> &pib = dt.pots_dif[ibb];
					int ipp = idPanelsM(ip, ibb);
					if (ipp >= 0) {
						pib[ipp][ih][ifr] = val;
						break;
					}
				}
			} else if (idFs[ip] >= 0)
				dt.fsPoints.pots_dif(idFs[ip], ih, ifr) = val;
			else
				dt.freePoints.pots_dif(idRest[ip], ih, ifr) = val;
		} else {								// pots_rad
			israd = true;
			int ip = f.GetInt(1);
			if (ip < 1 || ip > np)
				throw Exc(in.Str() + "\n"  + Format(t_("Point number %d is unknown"), ip));
			ip--;
			
			if (idPanels[ip]) {
				for (int ibb = 0; ibb < dt.Nb; ++ibb) {
					UArray<UArray<UArray<std::complex<double>>>> &pib = dt.pots_rad[ibb];
					int ipp = idPanelsM(ip, ibb);
					if (ipp >= 0) {
						for (int ib = 0; ib < dt.Nb; ++ib) {		
							int col0 = 0; 
							if (ib == 0)
								col0 = 1;
							
							for (int idof = 0; idof < 6; ++idof)
								pib[ipp][idof][ifr] += std::complex<double>(f.GetDouble(col0 + idof*2), f.GetDouble(col0 + idof*2 + 1));
							
							if (ib < dt.Nb-1)
								f.GetLine_discard_empty();
						}
						break;
					}
				}
			} else if (idFs[ip] >= 0) {
				for (int ib = 0; ib < dt.Nb; ++ib) {		
					int col0 = 0; 
					if (ib == 0)
						col0 = 1;
					
					for (int idof = 0; idof < 6; ++idof)
						dt.fsPoints.pots_rad(idFs[ip], idof, ifr) += std::complex<double>(f.GetDouble(col0 + idof*2), f.GetDouble(col0 + idof*2 + 1));
					
					if (ib < dt.Nb-1)
						f.GetLine_discard_empty();
				}
			} else {
				for (int ib = 0; ib < dt.Nb; ++ib) {		
					int col0 = 0; 
					if (ib == 0)
						col0 = 1;
					
					for (int idof = 0; idof < 6; ++idof)
						dt.freePoints.pots_rad(idRest[ip], idof, ifr) += std::complex<double>(f.GetDouble(col0 + idof*2), f.GetDouble(col0 + idof*2 + 1));
					
					if (ib < dt.Nb-1)
						f.GetLine_discard_empty();
				}
			}
		} 
	}
	if (!israd) {
		dt.pots_rad.Clear();
		dt.fsPoints.pots_rad.resize(0, 0, 0);
		dt.freePoints.pots_rad.resize(0, 0, 0);
	}
	if (!isdiff) {
		dt.pots_dif.Clear();
		dt.fsPoints.pots_dif.resize(0, 0, 0);
		dt.freePoints.pots_dif.resize(0, 0, 0);
	}

	return true;
}

void Wamit::Save_1(String fileName, bool force_T) const {
	if (!(IsLoadedA() && IsLoadedB())) 
		return; 
		
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to save '%s'. File already used."), fileName));
	
	if (IsLoadedA0()) {
		for (int i = 0; i < dt.Nb*6; ++i)  
			for (int j = 0; j < dt.Nb*6; ++j) 
				out << Format(" %s %5d %5d %s\n", FormatWam(-1), i+1, j+1,		// Always -1 is A0
								FormatWam(Nvl2(dt.A0(i, j), A0_ndim(i, j), 0.)));
	}
	if (IsLoadedAinf()) {
		for (int i = 0; i < dt.Nb*6; ++i)  
			for (int j = 0; j < dt.Nb*6; ++j)
				out << Format(" %s %5d %5d %s\n", FormatWam(0), i+1, j+1,		// Always 0 is Ainf
								FormatWam(Nvl2(dt.Ainf(i, j), Ainf_ndim(i, j), 0.)));
	}
	
	if (dt.Nf < 2)
		throw Exc(t_("No enough data to save (at least 2 frequencies)"));
	
	VectorXd data = force_T ? Get_T() : Get_w();
	
	for (int ifr = 0; ifr < dt.Nf; ++ifr)
		for (int i = 0; i < dt.Nb*6; ++i)  
			for (int j = 0; j < dt.Nb*6; ++j) 
				out << Format(" %s %5d %5d %s %s\n", FormatWam(data[ifr]), i+1, j+1,
							 FormatWam(Nvl2(dt.A[i][j][ifr], A_ndim(ifr, i, j), 0.)), 
							 FormatWam(Nvl2(dt.B[i][j][ifr], B_ndim(ifr, i, j), 0.)));
}

void Wamit::Save_3(String fileName, bool force_T) const {
	if (!IsLoadedFex()) 
		return;
	
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to save '%s'. File already used."), fileName));

	if (dt.Nf < 2)
		throw Exc(t_("No enough data to save (at least 2 frequencies)"));

	VectorXd data = force_T ? Get_T() : Get_w();
	
	// Insert heading -180 if 180 is available
	UVector<int> headids;
	Arange(headids, 0, dt.Nh-1, 1);
	UVector<double> head = clone(dt.head);	
	if (abs(Last(dt.head) - 180) < 0.1 && abs(First(dt.head) + 180) > 0.1) {
		headids.Insert(0, dt.Nh-1);
		head.Insert(0, -180);
	}
	for (int ifr = 0; ifr < dt.Nf; ++ifr)
		for (int ih = 0; ih < head.size(); ++ih)
			for (int ib = 0; ib < dt.Nb; ++ib) 
				for (int idf = 0; idf < 6; ++idf) {
					const std::complex<double> &f = dt.ex[ib][headids[ih]](ifr, idf);
					std::complex<double> fn = F_ndim(dt.ex, headids[ih], ifr, idf, ib);
					out << Format(" %s %s %5d %s %s %s %s\n", 
						FormatWam(data[ifr]), FormatWam(head[ih]), idf+1 + 6*ib,
						FormatWam(Nvl2(abs(f), abs(fn), 0.)), 
						FormatWam(Nvl2(arg(f), arg(fn)*180/M_PI, 0.)),
						FormatWam(Nvl2(f.real(), fn.real(), 0.)), 
						FormatWam(Nvl2(f.imag(), fn.imag(), 0.)));
				}
}

void Wamit::Save_4(String fileName, bool force_T) const {
	if (!IsLoadedRAO()) 
		return;
		
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to save '%s'. File already used."), fileName));
	
	if (dt.Nf < 2)
		throw Exc(t_("No enough data to save (at least 2 frequencies)"));
		
	VectorXd data = force_T ? Get_T() : Get_w();
	
	// Insert heading -180 if 180 is available
	UVector<int> headids;
	Arange(headids, 0, dt.Nh-1, 1);
	UVector<double> head = clone(dt.head);	
	if (abs(Last(dt.head) - 180) < 0.1 && abs(First(dt.head) + 180) > 0.1) {
		headids.Insert(0, dt.Nh-1);
		head.Insert(0, -180);
	}
	for (int ifr = 0; ifr < dt.Nf; ++ifr)
		for (int ih = 0; ih < head.size(); ++ih)
			for (int ib = 0; ib < dt.Nb; ++ib) 
				for (int idf = 0; idf < 6; ++idf) {
					const std::complex<double> &f = dt.rao[ib][headids[ih]](ifr, idf);	
					std::complex<double> fn = RAO_ndim(dt.rao, headids[ih], ifr, idf, ib);
					out << Format(" %s %s %5d %s %s %s %s\n", 
						FormatWam(data[ifr]), FormatWam(head[ih]), idf+1 + 6*ib,
						FormatWam(Nvl2(abs(f), abs(fn), 0.)), 
						FormatWam(Nvl2(arg(f), arg(fn)*180/M_PI, 0.)),
						FormatWam(Nvl2(f.real(), fn.real(), 0.)), 
						FormatWam(Nvl2(f.imag(), fn.imag(), 0.)));
				}
}

String WamitField(String str, int length) {
	String ret = str;
	if (length > ret.GetCount())
		ret << String(' ', length - ret.GetCount());
	return ret + " ";
}
	
void Wamit::Save_hst(String fileName) const {
	if (!IsLoadedC()) 
		return;
		
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to save '%s'. File already used."), fileName));

	for (int i = 0; i < 6*dt.Nb; ++i)  
		for (int j = 0; j < 6*dt.Nb; ++j) {
			int ib_i = i/6;
			int ii = i - ib_i*6;
			int ib_j = j/6;
			int jj = j - ib_j*6;
			
			out << Format(" %5d %5d  %s\n", i+1, j+1, 
					FormatWam(Nvl2(dt.msh[ib_i].dt.C(ii, jj), C_ndim(ib_i, ii, jj), 0.)));
		}
}

// K is supposed to be dimensionalized
void Wamit::Save_hst_static(const Eigen::MatrixXd &C, String fileName, double rho, double g) {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to save '%s'. File already used."), fileName));

	for (int i = 0; i < 6; ++i)  
		for (int j = 0; j < 6; ++j) 
			out << Format(" %5d %5d  %s\n", i+1, j+1, 
					FormatWam(Nvl2(C(i, j), C(i, j)/rho/g, 0.)));
}

void Wamit::Save_12(String fileName, bool isSum, Function <bool(String, int)> Status,
					bool force_T, bool force_Deg, int qtfHeadingId, double heading) const {
	if (!IsLoadedQTF(isSum)) 
		return;
	
	String ext = isSum ? ".12s" : ".12d";
	fileName = ForceExt(fileName, ext);
	
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to save '%s'. File already used."), fileName));
	
	int Nf = int(dt.qw.size());
	int Nh = int(dt.qhead.size()); 
	
	if (Nf < 2)
		throw Exc(t_("Not enough data to save (at least 2 frequencies)"));
	
	VectorXd data = force_T ? Get_qT() : Get_qw();
	
	// Insert heading -180 if 180 is available
	UVector<std::complex<double>> head;		// New heading pairs
	UVector<int> headids;					// Where are they
	
	if (abs(Last(dt.qhead) - std::complex<double>(180, 180)) < 0.1 && abs(First(dt.qhead) - std::complex<double>(-180, -180)) > 0.1) {
		UVector<double> allHeads;
		
		std::complex<double> qtfHeading;
		if (!IsNull(qtfHeadingId) && qtfHeadingId >= 0) {
		 	qtfHeading = dt.qhead[qtfHeadingId];
		 	allHeads << qtfHeading.real();
		} else {
			allHeads << -180;
			for (const std::complex<double> &h : dt.qhead)
				FindAddDelta(allHeads, h.real(), 0.01);
		}
		head.SetCount(pow2(allHeads.size()));
		headids.SetCount(pow2(allHeads.size()));
		int ic = 0;
		for (int ih1 = 0; ih1 < allHeads.size(); ++ih1) {
			for (int ih2 = 0; ih2 < allHeads.size(); ++ih2) {
				std::complex<double> h(allHeads[ih1], allHeads[ih2]);
				head[ic] = h;
				if (h.real() == -180)
					h.real(180);
				if (h.imag() == -180)
					h.imag(180);
				headids[ic] = FindRatio(dt.qhead, h, 0.01);		// Can be -1 if crossed headings do not exist
				ic++;
			}
		}
		if (!IsNull(qtfHeadingId) && qtfHeadingId >= 0)
			qtfHeadingId = FindClosest(head, qtfHeading);
	} else {
		 ::Copy(dt.qhead, head);
		 Arange(headids, 0, Nh-1, 1);	
	}
	
	const UArray<UArray<UArray<MatrixXcd>>> &qtf = isSum ? dt.qtfsum : dt.qtfdif;		
			
	out << " WAMIT Numeric Output -- Filename  " << Format("%20<s", GetFileName(fileName)) << "  " << Format("%", GetSysTime()) << "\n";

	int num = int(pow2(Nf)*head.size());
	int inum = 0;
	for (int ifr1 = 0; ifr1 < Nf; ++ifr1) 
		for (int ifr2 = 0; ifr2 < Nf; ++ifr2) 
			for (int ih = 0; ih < head.size(); ++ih) {
				inum++;
				if (Status && num >= 20 && !(inum%(num/20)) && !Status(Format("Saving %s", fileName), (100*inum)/num))
					throw Exc(t_("Stop by user"));
				
				double h1 = head[ih].real();
				double h2 = head[ih].imag();
				
				if (IsNull(qtfHeadingId) ||
					(qtfHeadingId == -1 && abs(h1 - h2) < 0.01) ||
					(ih == qtfHeadingId)) {
					if (ih == qtfHeadingId)
						h1 = h2 = heading;
					if (headids[ih] >= 0) {
						for (int ib = 0; ib < dt.Nb; ++ib)
		        			for (int idf = 0; idf < 6; ++idf) {
		        				static int idf12[] = {1, 3, 5, 2, 4, 6};
		        				int iidf = idf12[idf]-1;
		        				out << Format("   % 8.6E", data[ifr1]);
		        				out << Format("   % 8.6E", data[ifr2]);
		        				out << Format("   % 8.6E", h1);
		        				out << Format("   % 8.6E", h2);
		        				out << Format("   %2d", ib*6 + iidf+1);
		        				out << Format("   % 8.6E", F_ndim(abs(qtf[ib][headids[ih]][iidf](ifr1, ifr2)), iidf));
		        				out << Format("   % 8.6E", !force_Deg ? arg(qtf[ib][headids[ih]][iidf](ifr1, ifr2)) : ToDeg(arg(qtf[ib][headids[ih]][iidf](ifr1, ifr2))));
		        				out << Format("   % 8.6E", F_ndim(qtf[ib][headids[ih]][iidf](ifr1, ifr2).real(), iidf));
		        				out << Format("   % 8.6E", F_ndim(qtf[ib][headids[ih]][iidf](ifr1, ifr2).imag(), iidf));
		        				out << "\n";
		        			}
					}
				}
			}
}

void Wamit::Save_789(String fileName, bool force_T, bool force_Deg) const {
	if (!IsLoadedMD()) 
		return;
	
	String ext = Format(".%d", dt.mdtype);
	fileName = ForceExt(fileName, ext);
	
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to save '%s'. File already used."), fileName));
	
	int Nf = int(dt.w.size());
	int Nh = int(dt.mdhead.size()); 
	
	if (Nf < 2)
		throw Exc(t_("Not enough data to save (at least 2 frequencies)"));
	
	VectorXd data = force_T ? Get_T() : Get_w();
	
	// Insert heading -180 if 180 is available
	UVector<int> headids;
	Arange(headids, 0, Nh-1, 1);
	VectorXcd head;	
	if (abs(Last(dt.mdhead) - std::complex<double>(180, 180)) < 0.1 && abs(First(dt.mdhead) - std::complex<double>(-180, -180)) > 0.1) {
		headids.Insert(0, dt.Nh-1);
		head.resize(dt.mdhead.size()+1);
		head(0) = std::complex<double>(-180, -180);
		head.tail(dt.mdhead.size()) = dt.mdhead;
	} else 
		 head = dt.mdhead;
		
	const UArray<UArray<UArray<VectorXd>>> &md = dt.md;		
			
	for (int ifr = 0; ifr < dt.Nf; ++ifr)
		for (int ih = 0; ih < head.size(); ++ih) {
			for (int ib = 0; ib < dt.Nb; ++ib)
    			for (int idf = 0; idf < 6; ++idf) {
    				static int idf12[] = {1, 2, 3, 4, 5, 6};
    				int idof = idf12[idf]-1;
    				double re = F_ndim(md[ib][headids[ih]][idof](ifr), idof);
    
    				out << Format("   % 8.6E", data[ifr]);
    				out << Format("   % 8.6E", head[ih].real());
    				out << Format("   % 8.6E", head[ih].imag());
    				out << Format("   %2d", ib*6 + idof+1);
    				out << Format("   % 8.6E", abs(re));
    				out << Format("   % 8.6E", re > 0 ? 0. : 180.);
    				out << Format("   % 8.6E", re);
    				out << Format("   % 8.6E", 0.);
    				out << "\n";
    			}
		}
}

void Wamit::Save_frc2(String fileName, bool force1st, bool withQTF, UVector<Point3D> &listPoints) const {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to save '%s'. File already used."), fileName));

	out << "% BEMRosetta generated .frc file\n";
	out << WamitField(Format("%d %d %d %d 0 %d 0 0 0", force1st || (IsLoadedA() && IsLoadedB()) ? 1 : 0,
										 			   force1st || IsLoadedFex() ? 1 : 0,
										 			   force1st || IsLoadedFex() ? 2 : 0,
										 			   force1st || IsLoadedRAO() ? 1 : 0, 
										 			   !listPoints.IsEmpty() ? 1 : 0), 22);
	if (withQTF)
		out << " 0 0 1 0 0 0 0";
		
	out << "% 9 digits for Wamit .1 ... ." << (!withQTF ? 9 : 16) << " files included\n";
	out << WamitField(Format("%.1f", Bem().rho), 22) << "% RHO\n";
	
	for (int ib = 0; ib < dt.Nb; ++ib)
		out << Format("%.2f %.2f %.2f ", dt.msh[ib].dt.cg.x - dt.msh[ib].dt.c0.x, 
												   dt.msh[ib].dt.cg.y - dt.msh[ib].dt.c0.y, 
												   dt.msh[ib].dt.cg.z - dt.msh[ib].dt.c0.z);
	out << "% XCG, YCG, ZCG\n";
	
	if (IsLoadedM()) {
		out << "1  % IMASS\n";
		for (int ib = 0; ib < dt.Nb; ++ib) {
			for (int r = 0; r < 6; ++r) {
				int cc = 0;
				for (int c = 0; c < 6*dt.Nb; ++c) {
					if (c > 0)
						out << "\t";
					if (c < ib*6 || c >= (ib+1)*6)
						out << 0;
					else
						out << Format("%.2f", dt.msh[ib].dt.M(r, cc++));
				}
				out << "\n";
			}
		}
	} else 
		out << "0 % IMASS\n";	
	
	if (IsLoadedDlin()) {
		out << "1 % IDAMP\n";
		for (int ib = 0; ib < dt.Nb; ++ib) {
			for (int r = 0; r < 6; ++r) {
				int cc = 0;
				for (int c = 0; c < 6*dt.Nb; ++c) {
					if (c > 0)
						out << "\t";
					if (c < ib*6 || c >= (ib+1)*6)
						out << 0;
					else
						out << Format("%.2f", dt.msh[ib].dt.Dlin(r, cc++));
				}
				out << "\n";
			}
		}
	} else 
		out << "0 % IDAMP\n";
	
	
	if (IsLoadedCMoor()) {
		out << "1 % ISTIF\n";
		for (int ib = 0; ib < dt.Nb; ++ib) {
			for (int r = 0; r < 6; ++r) {
				int cc = 0;
				for (int c = 0; c < 6*dt.Nb; ++c) {
					if (c > 0)
						out << " ";
					if (c < ib*6 || c >= (ib+1)*6)
						out << 0;
					else
						out << Format("%.2f", dt.msh[ib].dt.Cmoor(r, cc++));
				}
				out << "\n";
			}
		}
	} else 
		out << "0 % ISTIF\n";
	
	out << "0 % NBETA\n";
	/*out << WamitField(Format("%d", dt.Nh), 12) << "% NBETA\n";
	UVector<double> head = clone(dt.head);
	Sort(head);
	for (int ih = 0; ih < dt.Nh; ++ih) 
		out << Format("%.4f ", head[ih]);
	out << "% BETA\n";*/
	
	out << listPoints.size() << " % NFIELD\n";
	for (const Point3D &p : listPoints)
		out << "   " << p.x << "\t" << p.y << "\t" << p.z << "\n";
}

void Wamit::Save_pot(String fileName, bool withMesh, bool x0z, bool y0z, const UArray<Body> &lids) const {
	String folder = GetFileFolder(fileName);
	
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to save '%s'. File already used."), fileName));

	out << "% BEMRosetta generated .pot file\n";
	if (!IsNum(dt.h) || dt.h < 0 || dt.h > 248)
		out << "-1           ";
	else
		out << WamitField(Format("%.2f", dt.h), 12);
	out << "% HBOT\n";
	out << WamitField("1 1", 12) << "% IRAD, IDIFF\n";
	out << WamitField(Format("%d", dt.Nf + 2), 12) << "% NPER\n";
	VectorXd TT	= Get_T();	
	Sort(TT);
	out << Format("%.4f ", -1);
	out << Format("%.4f ", 0);
	for (int iT = 0; iT < dt.Nf; ++iT) 
		out << Format("%.4f ", TT[iT]);
	out << "% PER(1), increment\n";
	
	// Insert heading -180 if 180 is available
	UVector<double> head = clone(dt.head);	
	if (abs(Last(dt.head) - 180) < 0.1 && abs(First(dt.head) + 180) > 0.1)
		head.Insert(0, -180);
	out << WamitField(Format("%d", head.size()), 12) << "% NBETA\n";
	for (int ih = 0; ih < head.size(); ++ih) 
		out << Format("%.4f ", head[ih]);
	out << "% BETA\n";
	
	out << WamitField(Format("%d ", dt.Nb), 12) << "% NBODY";
	UVector<String> names;
	for (int ib = 0; ib < dt.Nb; ++ib) {
		String name = Format("Body_%d", ib+1);
		name = ForceExt(name, ".gdf");
		out << "\n" << name << "\n";
		out << Format("%.3f %.3f %.3f 0.0 ", dt.msh[ib].dt.c0.x, dt.msh[ib].dt.c0.y, dt.msh[ib].dt.c0.z) << "% XBODY(1-4)\n";
		out << "1 1 1 1 1 1  % MODE(1:6)";
		names << AFX(folder, name);
	}
	if (withMesh) {
		UArray<Body> msh = clone(dt.msh);
		for (int ib = 0; ib < dt.Nb; ++ib) {
			if (lids.size() > ib) 
				msh[ib].dt.under.Append(lids[ib].dt.mesh);
			msh[ib].dt.under.Translate(-msh[ib].dt.c0.x, -msh[ib].dt.c0.y, -msh[ib].dt.c0.z);	
		}
		int nNodes, nPanels;
		Body::SaveAs(msh, names, Body::WAMIT_GDF, Body::UNDERWATER, Bem().rho, Bem().g, y0z, x0z, nNodes, nPanels,
			dt.w, dt.head, false, false, dt.h, Null);
	}
}

void Wamit::Save_Config(String folder, int numThreads) const {
	String fileBat = AFX(folder, "config.wam");
	FileOut file(fileBat);
	if (!file)
		throw Exc(Format(t_("Problem creating '%s' file"), fileBat));
	
	String wamitFolder = GetFileFolder(Bem().wamitPath);
	
	if (IsNull(numThreads) || numThreads <= 0)
		numThreads = 8;
			
	file << "! Generic configuration file:  config.wam\n"
 		<< " RAMGBMAX=4\n"
 		<< " NCPU=" << numThreads << "\n"
  		<< " USERID_PATH=" << wamitFolder << "   (directory for *.exe, *.dll, and userid.wam)\n"
  		<< " LICENSE_PATH=" << AFX(wamitFolder, "license") << "\n"
  	;
}

void Wamit::Save_cfg(String fileName, bool withQTF, bool lid, bool force_T, bool is6p) const {
	FileOut out(fileName);
	if (!out)
		throw Exc(Format(t_("Problem creating '%s' file"), fileName));
	
	bool isUnderWater = true;
	for (int ib = 0; ib < dt.Nb && isUnderWater; ++ib) {
		const UVector<Point3D> &nodes = dt.msh[ib].dt.under.nodes;
		for (const Point3D &n : nodes) {
			if (n.z >= -EPS_LEN) {
				isUnderWater = false;
				break;
			}
		}
	}
	
	out << "! BEMRosetta generated .cfg file\n"
		<< " IPERIN = 1       (Input period)\n";
	if (force_T)
 		out << " IPEROUT = 1      (Output period [s])\n";
 	else
 		out << " IPEROUT = 2      (Output frequency [rad/s])\n";
 	out << " MAXITT = 50\n";
 	if (withQTF)
		out << " I2ND = 1\n";
 	if (!isUnderWater) {
 		if (lid)
 			out << " IRR = 1\n";
 		else
 			out << " IRR = 3\n";
 	}
 	if (is6p)
 		out << " INUMOPT6 = 1\n";
 	out << " ILOG = 1\n"
 		<< " ISOLVE = 1\n"
 		<< " IALTFRC = 2\n"
		<< " ipltdat = 5\n"
		<< " ILOWHI = 0\n"
		<< " NUMHDR = 0\n"
	;
}

void Wamit::Save_Fnames(String folder) const {
	String fileBat = AFX(folder, "fnames.wam");
	FileOut file(fileBat);
	if (!file)
		throw Exc(Format(t_("Problem creating '%s' file"), fileBat));
	
	String folderName = GetFileTitle(folder);
	file << folderName << ".cfg\n"
		 << folderName << ".pot\n"
		 << folderName << ".frc";
}

void Wamit::SaveCase(String folder, int numThreads, bool withQTF, bool x0z, bool y0z, const UArray<Body> &lids, UVector<Point3D> &listPoints) const {
	if (!DirectoryCreateX(folder))
		throw Exc(Format(t_("Problem creating '%s' folder"), folder));

	Save_Fnames(folder);
	Save_Config(folder, numThreads);
	
	String folderName = GetFileTitle(folder);
	
	Save_pot(AFX(folder, folderName + ".pot"), true, x0z, y0z, lids);
	Save_frc2(AFX(folder, folderName + ".frc"), true, withQTF, listPoints);
	Save_cfg(AFX(folder, folderName + ".cfg"), withQTF, !lids.IsEmpty(), false, !listPoints.IsEmpty());
	
	String fileBat = AFX(folder, "Wamit_bat.bat");		
	FileOut bat(fileBat);
	if (!bat)
		throw Exc(Format(t_("Problem creating '%s' file"), fileBat));
	
	bat << "echo Start: \%date\% \%time\% >  time.txt\n";
	bat << "call \"" << Bem().wamitPath << "\" fnames.wam";
	bat << "\necho End:   \%date\% \%time\% >> time.txt\n";
}
