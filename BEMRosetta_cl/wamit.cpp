// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#include "functions.h"


bool Wamit::Load(String file, Function <bool(String, int)> Status) {
	hd().name = GetFileTitle(file);
	hd().file = file;	
	
	try {
		String ext = GetFileExt(file);
		if (ext == ".out") {
			hd().code = Hydro::WAMIT;
			BEM::Print("\n\n" + Format(t_("Loading out file '%s'"), file));
			if (!Load_out()) {
				BEM::PrintWarning("\n" + Format(t_("File '%s' not found"), file));
				return false;
			}
		} else if (S(".1.2.3.3sc.3fk.hst.4.12d.12s").Find(ext) >= 0) {
			if (GetFileName(GetFileFolder(file)) == "Wamit_format")
				hd().code = Hydro::HAMS_WAMIT;
			else if (hd().name == "WAMIT_5S")
				hd().code = Hydro::WADAM_WAMIT;
			else
				hd().code = Hydro::WAMIT_1_3;

			String filecfg = ForceExt(file, ".cfg");
			BEM::Print("\n- " + Format(t_("Configuration file .cfg file '%s'"), GetFileName(filecfg)));
			if (!Load_cfg(filecfg))
				BEM::Print(S(": ** cfg ") + t_("Not found") + "**");

			String filepot = ForceExt(file, ".pot");
			BEM::Print("\n- " + Format(t_("Configuration file .pot file '%s'"), GetFileName(filepot)));
			if (!Load_pot(filepot))
				BEM::Print(S(": ** pot ") + t_("Not found") + "**");
				
			String filegdf = ForceExt(file, ".gdf");
			BEM::Print("\n- " + Format(t_("Mesh file .gdf file '%s'"), GetFileName(filegdf)));
			if (!Load_gdf(filegdf))
				BEM::Print(S(": ** gdf ") + t_("Not found") + "**");
							
			String file1 = ForceExt(file, ".1");
			BEM::Print("\n- " + Format(t_("Hydrodynamic coefficients A and B .1 file '%s'"), GetFileName(file1)));
			if (!Load_1(file1))
				BEM::PrintWarning(S(": ** .1 ") + t_("Not found or empty") + "**");
			
			String file2 = ForceExt(file, ".2"),
				   file3 = ForceExt(file, ".3");
				   
			if (ext == ".2")
				;
			else {
				file = file3;
				if (!FileExists(file3) && FileExists(file2))
					file = file2;
			}
			BEM::Print("\n- " + Format(t_("Diffraction exciting %s file '%s'"), GetFileExt(file), GetFileName(file)));
			if (!Load_3(file))
				BEM::PrintWarning(S(": ** .3 ") + t_("Not found or empty") + "**");
			
			String fileHST = ForceExt(file, ".hst");
			BEM::Print("\n- " + Format(t_("Hydrostatic restoring file '%s'"), GetFileName(fileHST)));
			if (!Load_hst(fileHST))
				BEM::PrintWarning(S(": ** .hst ") + t_("Not found or empty") + "**");
		
			String fileRAO = ForceExt(file, ".4");
			BEM::Print("\n- " + Format(t_("RAO file '%s'"), GetFileName(fileRAO)));
			if (!Load_4(fileRAO))
				BEM::Print(S(": ** .4 ") + t_("Not found or empty") + "**");
			
			String file12s = ForceExt(file, ".12s");
			BEM::Print("\n- " + Format(t_("Second order sum coefficients .12s file '%s'"), GetFileName(file12s)));
			if (!Load_12(file12s, true, Status))
				BEM::Print(S(": ** .12s ") + t_("Not found") + "**");
			
			String file12d = ForceExt(file, ".12d");
			BEM::Print("\n- " + Format(t_("Second order mean drift coefficients .12d file '%s'"), GetFileName(file12d)));
			if (!Load_12(file12d, false, Status))
				BEM::Print(S(": ** .12d ") + t_("Not found") + "**");
		}
		
		String fileSC = ForceExt(file, ".3sc");
		BEM::Print("\n- " + Format(t_("Scattering file '%s'"), GetFileName(fileSC)));
		if (!Load_Scattering(fileSC))
			BEM::Print(S(": ** 3sc ") + t_("Not found") + "**");
		String fileFK = ForceExt(file, ".3fk");
		BEM::Print("\n- " + Format(t_("Froude-Krylov file '%s'"), GetFileName(fileFK)));
		if (!Load_FK(fileFK))
			BEM::Print(S(": ** 3fk ") + t_("Not found") + "**");
		
		if (IsNull(hd().Nh))
			hd().Nh = 0;
		if (IsNull(hd().Nf))
			hd().Nf = 0;

	
		if (IsNull(hd().Nb)/* || IsNull(hd().Nh) || IsNull(hd().Nf) || hd().Nh == 0 || hd().Nf == 0*/)
			return false;
		
		hd().dof.Clear();	hd().dof.SetCount(hd().Nb, 6);

	} catch (Exc e) {
		Status("", -1);
		BEM::PrintError(Format("\n%s: %s", t_("Error"), e));
		hd().lastError = e;
		return false;
	}
	Status("", -1);
	return true;
}

bool Wamit::Save(String file, Function <bool(String, int)> Status, bool force_T, int qtfHeading) {
	try {
		if (hd().IsLoadedA() && hd().IsLoadedB()) {
			String file1 = ForceExt(file, ".1");
			BEM::Print("\n- " + Format(t_("Hydrodynamic coefficients A and B file '%s'"), GetFileName(file1)));
			Save_1(file1, force_T);
		}
		if (hd().IsLoadedFex()) {
			String file3 = ForceExt(file, ".3");
			BEM::Print("\n- " + Format(t_("Diffraction exciting file '%s'"), GetFileName(file3)));
			Save_3(file3, force_T);
		}
		if (hd().IsLoadedC()) {
			String fileHST = ForceExt(file, ".hst");
			BEM::Print("\n- " + Format(t_("Hydrostatic restoring file '%s'"), GetFileName(fileHST)));
			Save_hst(fileHST);
		}
		if (hd().IsLoadedRAO()) {
			String fileRAO = ForceExt(file, ".4");
			BEM::Print("\n- " + Format(t_("RAO file '%s'"), GetFileName(fileRAO)));
			Save_4(fileRAO, force_T);
		}
		if (hd().IsLoadedQTF(true)) {
			String fileQTFs = ForceExt(file, ".12s");
			BEM::Print("\n- " + Format(t_("QTF file '%s'"), GetFileName(fileQTFs)));
			Save_12(fileQTFs, true, Status, force_T, true, qtfHeading);
		}
		if (hd().IsLoadedQTF(false)) {
			String fileQTFd = ForceExt(file, ".12d");
			BEM::Print("\n- " + Format(t_("QTF file '%s'"), GetFileName(fileQTFd)));
			Save_12(fileQTFd, false, Status, force_T, true, qtfHeading);
		}
	} catch (Exc e) {
		BEM::PrintError(Format("\n%s: %s", t_("Error"), e));
		hd().lastError = e;
		return false;
	}
	return true;
}

bool Wamit::Load_out() {
	hd().Nb = 0;
	hd().Nf = 0;
	hd().Nh = 0;
	int pos;
	int ibody = -1;
	hd().dimen = false;
	
	FileInLine in(hd().file);
	if (!in.IsOpen())
		return false;
	String line;
	FieldSplitWamit f(in);
	f.IsSeparator = IsTabSpace;
	
	hd().names.Clear();
	while(!in.IsEof()) {
		line = in.GetLine();
		f.Load(line);
		if (line.Find("N=") >= 0 && line.Find("Body number:") < 0) {
			hd().Nb++;
			hd().names << GetFileTitle(f.GetText(2));
		} else if ((pos = line.FindAfter("Input from Geometric Data File:")) >= 0) {
			hd().Nb = 1;
			hd().names << GetFileTitle(TrimBoth(line.Mid(pos)));
		} else if (line.Find("POTEN run date and starting time:") >= 0) {
			hd().cg.setConstant(3, hd().Nb, NaNDouble);
			hd().cb.setConstant(3, hd().Nb, NaNDouble);
			hd().c0.setConstant(3, hd().Nb, 0);
			hd().Vo.SetCount(hd().Nb, NaNDouble);
			hd().C.SetCount(hd().Nb);
			hd().M.SetCount(hd().Nb);
		} else if (line.Find("Gravity:") >= 0) {
			hd().g = f.GetDouble(1);
			hd().len = f.GetDouble(4);
		} else if (line.Find("Water depth:") >= 0) {
			if (ToLower(f.GetText(2)) == "infinite")
				hd().h = -1;
			else {
				hd().h = f.GetDouble(2);
				if (hd().h < 0)
					throw Exc(in.Str() + "\n" +  t_("Water depth has to be positive"));
			}
			if (line.Find("Water density:") >= 0) 
				hd().rho = f.GetDouble(5);			
		} else if (line.Find("XBODY =") >= 0) {
			ibody++;
			if (ibody >= hd().Nb)
				throw Exc(in.Str() + "\n"  + Format(t_("Found additional bodies over %d"), hd().Nb));
			if (hd().c0.rows() < 3 || hd().c0.cols() < hd().Nb)
			 	throw Exc(in.Str() + "\n"  + t_("c0 matrix is not dimensioned"));
			hd().c0(0, ibody) = f.GetDouble(2);
			hd().c0(1, ibody) = f.GetDouble(5);
			hd().c0(2, ibody) = f.GetDouble(8);
		} else if ((pos = line.FindAfter("Volumes (VOLX,VOLY,VOLZ):")) >= 0) {
			if (hd().Vo.size() < hd().Nb)
			 	throw Exc(in.Str() + "\n"  + t_("Vo matrix is not dimensioned"));		
			hd().Vo[ibody] = ScanDouble(line.Mid(pos));
		} else if (line.Find("Center of Gravity  (Xg,Yg,Zg):") >= 0) {
			if (hd().cg.rows() < 3 || hd().cg.cols() < hd().Nb)
			 	throw Exc(in.Str() + "\n"  + t_("cg matrix is not dimensioned"));
			hd().cg(0, ibody) = f.GetDouble(4);
			hd().cg(1, ibody) = f.GetDouble(5);
			hd().cg(2, ibody) = f.GetDouble(6);
		} else if (line.Find("Center of Buoyancy (Xb,Yb,Zb):") >= 0) {
			if (hd().cb.rows() < 3 || hd().cb.cols() < hd().Nb)
			 	throw Exc(in.Str() + "\n"  + t_("cb matrix is not dimensioned"));
			hd().cb(0, ibody) = f.GetDouble(4);
			hd().cb(1, ibody) = f.GetDouble(5);
			hd().cb(2, ibody) = f.GetDouble(6);
		} else if (line.Find("Radii of gyration:") >= 0) {
			if (IsNull(hd().rho))
				;
			else {
				double mass = hd().rho*hd().Vo[ibody];
				if (hd().M.size() < hd().Nb)
				 	throw Exc(in.Str() + "\n"  + t_("M matrix is not dimensioned"));
				if (hd().M[ibody].size() > 0)
					throw Exc(in.Str() + "\n"  + t_("Problem in M matrix. Please report"));
				hd().M[ibody].setConstant(6, 6, 0);
				Eigen::MatrixXd &inertia = hd().M[ibody];
				inertia(0, 0) = inertia(1, 1) = inertia(2, 2) = mass;
				inertia(3, 3) = f.GetDouble(3)*mass;
				inertia(3, 4) = f.GetDouble(4)*mass;
				inertia(3, 5) = f.GetDouble(5)*mass;
				f.GetLine();
				inertia(4, 3) = f.GetDouble(0)*mass;
				inertia(4, 4) = f.GetDouble(1)*mass;
				inertia(4, 5) = f.GetDouble(2)*mass;
				f.GetLine();
				inertia(5, 3) = f.GetDouble(0)*mass;
				inertia(5, 4) = f.GetDouble(1)*mass;
				inertia(5, 5) = f.GetDouble(2)*mass;
				double cx = mass*hd().cg(0, ibody);
				double cy = mass*hd().cg(1, ibody);
				double cz = mass*hd().cg(2, ibody);
				inertia(1, 5) = inertia(5, 1) =  cx;
				inertia(2, 4) = inertia(4, 2) = -cx;
				inertia(0, 5) = inertia(5, 0) = -cy;
				inertia(2, 3) = inertia(3, 2) =  cy;
				inertia(0, 4) = inertia(4, 0) =  cz;
				inertia(1, 3) = inertia(3, 1) = -cz;
			}
		} else if (line.Find("Global body and external mass matrix:") >= 0) {
			if (hd().M.size() < hd().Nb)
			 	throw Exc(in.Str() + "\n"  + t_("M matrix is not dimensioned"));
			if (hd().M[ibody].size() > 0)
				throw Exc(in.Str() + "\n"  + t_("Problem in M matrix. Please report"));
			hd().M[ibody].setConstant(6, 6, NaNDouble);
			for (int r = 0; r < 6; ++r) {
				f.GetLine();
				for (int c = 0; c < 6; ++c)
					hd().M[ibody](r, c) = f.GetDouble(c);
			}		
		} else if (line.Find("Hydrostatic and gravitational") >= 0) {
			if (hd().C.size() < hd().Nb)
			 	throw Exc(in.Str() + "\n"  + t_("C matrix is not dimensioned"));
			hd().C[ibody].setConstant(6, 6, 0);
			f.LoadWamitJoinedFields(in.GetLine());
			hd().C[ibody](2, 2) = f.GetDouble(1);
			hd().C[ibody](2, 3) = hd().C[ibody](3, 2) = f.GetDouble(2);
			hd().C[ibody](2, 4) = hd().C[ibody](4, 2) = f.GetDouble(3);
			f.LoadWamitJoinedFields(in.GetLine());
			hd().C[ibody](3, 3) = f.GetDouble(1);
			hd().C[ibody](3, 4) = hd().C[ibody](4, 3) = f.GetDouble(2);
			hd().C[ibody](3, 5) = hd().C[ibody](5, 3) = f.GetDouble(3);
			f.LoadWamitJoinedFields(in.GetLine());
			hd().C[ibody](4, 4) = f.GetDouble(1);
			hd().C[ibody](4, 5) = hd().C[ibody](5, 4) = f.GetDouble(2);
		} else if (line.Find("Output from  WAMIT") >= 0) {
			hd().head.Clear();
			FileInLine::Pos fpos = in.GetPos();
			
			bool foundNh = false, found2ndorder = false;
			while (!in.IsEof()) {
				line = in.GetLine();
				if (line.Find("Wave period (sec)") >= 0) {
					++hd().Nf;
					if (hd().head.size() > 0 && !foundNh)
						foundNh = true;
				} else if (!foundNh) {
					if (hd().head.size() > 0 && (line.Find("*********************") >= 0 ||
								   line.Find("FORCES AND MOMENTS") >= 0)) 
						foundNh = true;
					else if (line.Find("Wave Heading (deg) :") >= 0) {
						f.Load(line);
						FindAddDelta(hd().head, f.GetDouble(4)/*FixHeading_180(f.GetDouble(4))*/, 0.001);
					}
				} else if (line.Find("2nd-order") >= 0) {
					found2ndorder = true;
					break;
				}
			}
			Sort(hd().head);
			hd().Nh = hd().head.size();
			if (hd().Nb == 0 || hd().Nh == 0 || hd().Nf == 0)
				throw Exc(in.Str() + "\n"  + Format(t_("Wrong format in Wamit file '%s'"), hd().file));
		
			hd().T.SetCount(hd().Nf);
			hd().w.SetCount(hd().Nf);

			hd().A.SetCount(6*hd().Nb);
			hd().B.SetCount(6*hd().Nb);
			for (int i = 0; i < 6*hd().Nb; ++i) {
				hd().A[i].SetCount(6*hd().Nb);
				hd().B[i].SetCount(6*hd().Nb);
				for (int j = 0; j < 6*hd().Nb; ++j) {
					hd().A[i][j].setConstant(hd().Nf, NaNDouble);// In Wamit, unloaded DOFs are considered negligible	
					hd().B[i][j].setConstant(hd().Nf, NaNDouble);	
				}
			}
						
			int qtfNh = 0, qtfNf = 0;
			UArray<std::complex<double>> head;
			if (found2ndorder) {
				bool foundSum = false, foundDif = false;
				while (!in.IsEof()) {
					f.GetLine();
					if (f.IsInLine("Period indices:")) {
						int if1 = f.GetInt(2);
						int if2 = f.GetInt(3);
						if (if1 > qtfNf)
							qtfNf = if1;
						if (if2 > qtfNf)
							qtfNf = if2;
					} else if (f.IsInLine("Headings (deg):")) {
						double hd1 = f.GetDouble(6);
						double hd2 = f.GetDouble(7);
						FindAdd(head, std::complex<double>(hd1, hd2));	
					} else if (!foundSum && f.IsInLine("SUM-FREQUENCY")) 
						foundSum = true;
					else if (!foundDif && f.IsInLine("DIFFERENCE-FREQUENCY")) 
						foundDif = true;
				}
				hd().qh.resize(qtfNh = head.size());
				for (int i = 0; i < head.size(); ++i)
					hd().qh[i] = head[i];
				hd().qw.resize(qtfNf);
				
				if (foundSum)
					hd().InitQTF(hd().qtfsum, hd().Nb, qtfNh, qtfNf);
				if (foundDif)
					hd().InitQTF(hd().qtfdif, hd().Nb, qtfNh, qtfNf);
			}
			
			in.SeekPos(fpos);
			while (in.GetLine().Find("Wave period = infinite") < 0 && !in.IsEof())
				; 
			if (!in.IsEof()) {
				hd().A0.setConstant(hd().Nb*6, hd().Nb*6, NaNDouble);
				Load_A(in, hd().A0);
			}
			in.SeekPos(fpos);
			while (in.GetLine().Find("Wave period = zero") < 0 && !in.IsEof())
				; 
			if (!in.IsEof()) {
				hd().Ainf.setConstant(hd().Nb*6, hd().Nb*6, NaNDouble);
				Load_A(in, hd().Ainf);
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
				if (OUTB(ifr, hd().Nf))
					throw Exc(in.Str() + "\n" + Format(t_("Found additional frequencies over %d"), hd().Nf));
				
	            hd().T[ifr] = f.GetDouble(4);  			
	            hd().w[ifr] = 2*M_PI/hd().T[ifr];//fround(2*M_PI/hd().T[ifr], 8);
	            hd().dataFromW = false;
	            
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
							if (OUTB(i, hd().Nb*6) || OUTB(j, hd().Nb*6))
								throw Exc(in.Str() + "\n"  + Format(t_("Index (%d, %d) out of bounds"), i, j));
							hd().A[i][j][ifr] = Aij;
							hd().B[i][j][ifr] = Bij;
						}
						hd().GetBodyDOF();
	            	} else if (line.Find("DIFFRACTION EXCITING FORCES AND MOMENTS") >= 0) {
						if (hd().ex.force.IsEmpty()) 
							hd().Initialize_Forces(hd().ex);
						
						int ih = 0;
						while (!in.IsEof()) {		
							line = in.GetLine();
							if (line.Find("Wave Heading (deg) :") >= 0) {
								f.Load(line);
								double head = f.GetDouble(4); //FixHeading_180(f.GetDouble(4)); 
								int iih = FindDelta(hd().head, head, 0.001);
								if (iih >= 0) {
									in.GetLine(3); 
									while (!TrimBoth(line = in.GetLine()).IsEmpty()) {
										f.Load(line);
										double ma = f.GetDouble(1);
										double ph = ToRad(f.GetDouble(2));
										int ib = abs(f.GetInt(0)) - 1;
										if (OUTB(ih, hd().Nh) || OUTB(ifr, hd().Nf) || OUTB(ib, hd().Nb*6))
											throw Exc(in.Str() + "\n"  + Format(t_("Index [%d](%d, %d) out of bounds"), ih, ifr, ib));
										hd().ex.force[ih](ifr, ib) = std::polar(ma, ph);	
									}
								}
								ih++;
								if (ih >= hd().Nh)
									break;
							}
						} 
					} else if (line.Find("RESPONSE AMPLITUDE OPERATORS") >= 0) {
						if (hd().rao.force.IsEmpty()) 
							hd().Initialize_Forces(hd().rao);
						
						int ih = 0;
						while (!in.IsEof()) {		
							line = in.GetLine();
							if (line.Find("Wave Heading (deg) :") >= 0) {
								f.Load(line);
								double head = f.GetDouble(4); 	//FixHeading_180(f.GetDouble(4)); 
								int iih = FindDelta(hd().head, head, 0.001);
								if (iih >= 0) {
									in.GetLine(3); 
									while (!TrimBoth(line = in.GetLine()).IsEmpty()) {
										f.Load(line);
										double ma = f.GetDouble(1);
										double ph = ToRad(f.GetDouble(2));
										int i = abs(f.GetInt(0)) - 1;
										if (OUTB(ih, hd().Nh) || OUTB(ifr, hd().Nf) || OUTB(i, hd().Nb*6))
											throw Exc(in.Str() + "\n"  + Format(t_("Index [%d](%d, %d) out of bounds"), ih, ifr, i));
										hd().rao.force[ih](ifr, i) = std::polar(ma, ph);	
									}
								}
								ih++;
								if (ih >= hd().Nh)
									break;
							}
						}
					} else if (line.Find("SURGE, SWAY & YAW DRIFT FORCES") >= 0 ||
							   line.Find("SURGE, SWAY, HEAVE, ROLL, PITCH & YAW DRIFT FORCES") >= 0 ||
							   line.Find("VELOCITY VECTOR IN FLUID DOMAIN") >= 0 ||
							   line.Find("HYDRODYNAMIC PRESSURE IN FLUID DOMAIN") >= 0 ||
							   line.Find("*************************************") >= 0) {
						nextFreq = true;
						break;
					} 
				}
			}
			if (found2ndorder) {
				hd().qtfdataFromW = false;
				int ifr1, ifr2, ih; 
				UArray<UArray<UArray<MatrixXcd>>> *qtf = nullptr;
				while (!in.IsEof()) {
					f.GetLine();
					if (f.IsInLine("Period indices:")) {
						ifr1 = f.GetInt(2)-1; 
						ifr2 = f.GetInt(3)-1;
						hd().qw[ifr1] = 2*M_PI/f.GetDouble(5);
						hd().qw[ifr2] = 2*M_PI/f.GetDouble(6);
					} else if (f.IsInLine("SUM-FREQUENCY")) 
						qtf = &hd().qtfsum;
					else if (f.IsInLine("DIFFERENCE-FREQUENCY")) 
						qtf = &hd().qtfdif;
					else if (f.IsInLine("Heading indices:")) {
						double hd1 = f.GetDouble(6);
						double hd2 = f.GetDouble(7);
						ih = Find(head, std::complex<double>(hd1, hd2));
					} else if (f.IsInLine("Mod[F2(I)]")) {
						f.GetLine();
						while (!Trim(f.GetLine()).IsEmpty()) {						
							int idof = f.GetInt(0)-1;
							int ib = idof/6;
							double ma = f.GetDouble(1);
							double ph = ToRad(f.GetDouble(2));
							if (OUTB(ih, qtfNh) || OUTB(ifr1, qtfNf) || OUTB(ifr2, qtfNf) || OUTB(ib, hd().Nb))
								throw Exc(in.Str() + "\n"  + Format(t_("Index [%d][%d](%d, %d) out of bounds"), ib, ih, idof, ifr1, ifr2));
							(*qtf)[ib][ih][idof](ifr1, ifr2) = std::polar(ma, ph);
							(*qtf)[ib][ih][idof](ifr2, ifr1) = std::polar(ma, ph);	
						}
					}
				}
			}
		}
	}
	if (hd().Nb == 0)
		throw Exc(t_("Incorrect .out format"));
	
	return true;
}

void Wamit::Save_A(FileOut &out, Function <double(int, int)> fun, const Eigen::MatrixXd &base, String wavePeriod) {
	out << 	" ************************************************************************\n\n"
			" Wave period = " << wavePeriod << "\n"
			" ------------------------------------------------------------------------\n\n\n"
			"    ADDED-MASS COEFFICIENTS\n"
			"     I     J         A(I,J)\n\n";
	for (int r = 0; r < hd().Nb*6; ++r) 
		for (int c = 0; c < hd().Nb*6; ++c) 
			if (IsNum(base(r, c))) 
				out << Format("%6>d%6>d  % E\n", r+1, c+1, fun(r, c));
	out << "\n\n";
}

void Wamit::Save_AB(FileOut &out, int ifr) {
	out <<	"    ADDED-MASS AND DAMPING COEFFICIENTS\n"
			"     I     J         A(I,J)         B(I,J)\n\n";
	for (int r = 0; r < hd().Nb*6; ++r) 
		for (int c = 0; c < hd().Nb*6; ++c) 
			if (IsNum(hd().A[r][c][ifr]) && !IsNull(hd().B[r][c][ifr]))
				out << Format("%6>d%6>d  % E  % E\n", r+1, c+1, hd().A_ndim(ifr, r, c), hd().B_ndim(ifr, r, c));
	out << "\n\n\n\n";
}

void Wamit::Save_Forces(FileOut &out, int ifr) {
	out <<	"    DIFFRACTION EXCITING FORCES AND MOMENTS\n\n";
	for (int ih = 0; ih < hd().Nh; ++ih) {
		out << "  Wave Heading (deg) :      " << hd().head[ih] << "\n\n"
			<< "     I     Mod[Xh(I)]     Pha[Xh(I)]\n\n";
		for (int i = 0; i < hd().ex.force[ih].cols(); ++i)
			if (IsNum(hd().ex.force[ih](ifr, i))) {
				std::complex<double> c = hd().F_ndim(hd().ex, ih, ifr, i);
				out << Format("%6>d   %E         %6>d\n", i+1, abs(c), round(ToDeg(arg(c))));
			}
		out << "\n\n";
	}
}

void Wamit::Save_RAO(FileOut &out, int ifr) {
	out <<	"    RESPONSE AMPLITUDE OPERATORS\n\n";
	for (int ih = 0; ih < hd().Nh; ++ih) {
		out << "  Wave Heading (deg) :      " << hd().head[ih] << "\n\n"
			<< "     I     Mod[Xh(I)]     Pha[Xh(I)]\n\n";
		for (int i = 0; i < hd().rao.force[ih].cols(); ++i)
			if (IsNum(hd().rao.force[ih](ifr, i))) {
				std::complex<double> c = hd().F_ndim(hd().rao, ih, ifr, i);
				out << Format(" %7>d   %E   %f\n", i+1, abs(c), ToDeg(arg(c)));
			}
		out << "\n\n\n\n";
	}
}

bool Wamit::Save_out(String file, double g, double rho) {
	try {
		String filename = GetFileTitle(file);
		FileOut out(file);
		if (!out.IsOpen())
			throw Exc(Format(t_("Impossible to open '%s'"), file));

		out << Format(" %s\n\n", String('-', 71))
			<< "                        WAMIT  Version 6.1/6.107S \n\n"
		    << "                        BEMRosetta generated .out format\n\n\n\n"
		    << Format(" %s\n\n\n", String('-', 72))
		    << Format(" %s\n\n\n", String('-', 72))
			<< " Low-order panel method  (ILOWHI=0)\n\n";
		if (hd().Nb == 1)
			out << " Input from Geometric Data File:         " << filename << ".gdf\n"
				<< " Unknown gdf file source\n\n";
		else {
			out << " Input from Geometric Data Files:\n";
			for (int ib = 0; ib < hd().Nb; ++ib)
				out << Format("                               N=  %d     %s%d.gdf\n"
							  " Unknown gdf file source\n\n", ib+1, filename, ib);
		}
		out	<< " Input from Potential Control File:      " << filename << ".pot\n"
			<< " " << filename << ".pot -- file type .gdf, ILOWHI=0, IRR=1\n\n\n"
			<< " POTEN run date and starting time:        01-Jan-2000  --  00:00:00\n"
			<< "   Period       Time           RAD      DIFF  (max iterations)\n";
		if (hd().IsLoadedA0())
			out << "   -1.0000    00:00:00          -1\n";
		if (hd().IsLoadedAinf())
    		out << "    0.0000    00:00:00          -1\n";
		for (int it = 0; it < hd().T.size(); ++it)
			out << " " << Format("%9.4f", hd().T[it]) << "    00:00:00          -1      -1\n";
		out << "\n"
 		   	<< " Gravity:     " << (IsNull(hd().g) ? g : hd().g)
 		    << "                Length scale:        " << hd().len << "\n"
 			<< " Water depth:        " << (hd().h < 0 ? "infinite" : FDS(hd().h, 9)) << "    "
 			<< " Water density:      " << (IsNull(hd().rho) ? rho : hd().rho) << "\n"
 			<< " Logarithmic singularity index:              ILOG =     1\n"
 			<< " Source formulation index:                   ISOR =     0\n"
 			<< " Diffraction/scattering formulation index: ISCATT =     0\n"
 			<< " Number of blocks used in linear system:   ISOLVE =     1\n"
 			<< " Number of unknowns in linear system:        NEQN =  111\n\n";      
 		
		out << " BODY PARAMETERS:\n\n";
		
		for (int ibody = 0; ibody < hd().Nb; ++ibody) {
			if (hd().Nb > 1)
				out << " Body number: N= " << ibody+1 << "   ";
			out	<< " Total panels:  1111    Waterline panels:   11      Symmetries: none\n";
			out	<< " Irregular frequency index: IRR =1\n"; 
			out	<< " Free surface panels:     111\n\n";
			out	<< Format(" XBODY =    %s YBODY =    %s ZBODY =    %s PHIBODY =   0.0\n", 
							FormatWam(hd().c0(0, ibody)), 
							FormatWam(hd().c0(1, ibody)), 
							FormatWam(hd().c0(2, ibody)));
			double Vo = hd().Vo.size() > ibody ? hd().Vo[ibody] : 0;
			out	<< Format(" Volumes (VOLX,VOLY,VOLZ):      %s %s %s\n", 
						FormatWam(Vo), FormatWam(Vo), FormatWam(Vo));
			double cbx = 0, cby = 0, cbz = 0;
			if (hd().cb.size() > 0) {
				cbx = hd().cb(0, ibody);
				cby = hd().cb(1, ibody);
				cbz = hd().cb(2, ibody);
			}
			out	<< Format(" Center of Buoyancy (Xb,Yb,Zb): %s %s %s\n", 
							FormatWam(cbx), FormatWam(cby), FormatWam(cbz));
			if (hd().IsLoadedC()) {
				out	<< " Hydrostatic and gravitational restoring coefficients:\n"; 
				out	<< " C(3,3),C(3,4),C(3,5): " << Format("%s %s %s\n", 
							FormatWam(hd().C_ndim(ibody, 2, 2)), 
							FormatWam(hd().C_ndim(ibody, 2, 3)), 
							FormatWam(hd().C_ndim(ibody, 2, 4)));
				out	<< " C(4,4),C(4,5),C(4,6):               " << Format("%s %s %s\n", 
							FormatWam(hd().C_ndim(ibody, 3, 3)), 
							FormatWam(hd().C_ndim(ibody, 3, 4)), 
							FormatWam(hd().C_ndim(ibody, 3, 5)));
				out	<< "        C(5,5),C(5,6):                             " << Format("%s %s\n", 
							FormatWam(hd().C_ndim(ibody, 4, 4)), 
							FormatWam(hd().C_ndim(ibody, 4, 5)));
			}
			double cgx = 0, cgy = 0, cgz = 0;
			if (hd().cb.size() > 0) {
				cgx = hd().cg(0, ibody);
				cgy = hd().cg(1, ibody);
				cgz = hd().cg(2, ibody);
			}
			out	<< Format(" Center of Gravity  (Xg,Yg,Zg): %s %s %s\n\n", 
						FormatWam(cgx), FormatWam(cgy), FormatWam(cgz));
			if (hd().IsLoadedM()) {
				out	<< " Global body and external mass matrix:        ";
				for (int r = 0; r < 6; ++r) {
					out << "\n ";
					for (int c = 0; c < 6; ++c)
						out << Format(" % 10.4E", hd().M[ibody](r, c));
				}
			}
		}
		out << "\n\n\n";
		out << " ------------------------------------------------------------------------\n"
        	<< "                            Output from  WAMIT\n"
			<< " ------------------------------------------------------------------------\n"
			<< " FORCE run date and starting time:                22-Dec-2016 -- 11:54:03\n"
			<< " ------------------------------------------------------------------------\n"
			<< " I/O Files:         " << filename << ".frc       " << filename << ".p2f       " << filename << ".out\n"
			<< "  " << filename << ".frc -- file type .gdf, ILOWHI=0, IRR=1\n \n\n"
		;
 
		if (hd().IsLoadedA0())
			Save_A(out, [&](int idf, int jdf)->double {return hd().A0_ndim(idf, jdf);},   hd().A0,   "infinite");
		if (hd().IsLoadedAinf())
			Save_A(out, [&](int idf, int jdf)->double {return hd().Ainf_ndim(idf, jdf);}, hd().Ainf, "zero");
		
		for (int ifr = 0; ifr < hd().T.size(); ++ifr) {
			out << 	" ************************************************************************\n\n"
					" Wave period (sec) = " << FormatWam(hd().T[ifr]) << "\n"
					" ------------------------------------------------------------------------\n\n\n";
			if (hd().IsLoadedA() && hd().IsLoadedB()) 
				Save_AB(out, ifr);
			if (hd().IsLoadedFex())
				Save_Forces(out, ifr);
			if (hd().IsLoadedRAO())
				Save_RAO(out, ifr);
		}
		UVector<double> heads;
		UVector<int> id1(int(hd().qh.size())), id2(int(hd().qh.size()));
		if (hd().IsLoadedQTF(true) || hd().IsLoadedQTF(false)) {
			for (int ih = 0; ih < hd().qh.size(); ++ih) {
				id1[ih] = 1+FindAdd(heads, hd().qh[ih].real());
				id2[ih] = 1+FindAdd(heads, hd().qh[ih].imag());
			}
		}
			
		if (hd().IsLoadedQTF(true)) {
			for (int ih = 0; ih < hd().qh.size(); ++ih) {
				for (int ifr1 = 0; ifr1 < hd().qw.size(); ++ifr1) {
					for (int ifr2 = ifr1; ifr2 < hd().qw.size(); ++ifr2) {
						double T = 2*M_PI/(hd().qw[ifr1] + hd().qw[ifr2]);
						double T1 = 2*M_PI/(hd().qw[ifr1]);
						double T2 = 2*M_PI/(hd().qw[ifr2]);
						out << " ************************************************************************\n\n"
							   " 2nd-order period (sec) =  " << Format("%12E", T) << "\n\n"
							   " Period indices:     " << Format("%2d   %2d", ifr1+1, ifr2+1) << "          Periods:  " 
							   		<< Format("%12E %12E", T1, T2) << "\n"
							   " ------------------------------------------------------------------------\n\n\n"
							   "SUM-FREQUENCY EXCITING FORCES AND MOMENTS-DIRECT METHOD\n\n"
							   "  Heading indices:    " << Format("%2d   %2d", id1[ih], id2[ih]) << "  Headings (deg):       " 
							   		<< Format("%4.1f", hd().qh[ih].real()) << "      " 
							   		<< Format("%4.1f", hd().qh[ih].imag()) << "\n\n\n"
							   "      I     Mod[F2(I)]     Pha[F2(I)]\n\n";
						for (int ib = 0; ib < hd().Nb; ++ib) {
							for (int idf = 0; idf < 6; ++idf) {
								static int idf12[] = {1, 3, 5, 2, 4, 6};
				        		int iidf = idf12[idf]-1;
								out << Format("    %2d   %12.6E     %10d\n", 1+iidf + 6*ib, 
										hd().F_ndim(abs(hd().qtfsum[ib][ih][iidf](ifr1, ifr2)), iidf),
     									int(ToDeg(arg(hd().qtfsum[ib][ih][iidf](ifr1, ifr2)))));
							}
						}
						out << "\n\n";
					}
				}
			}
		}
		if (hd().IsLoadedQTF(false)) {
			for (int ih = 0; ih < hd().qh.size(); ++ih) {
				for (int ifr2 = 0; ifr2 < hd().qw.size(); ++ifr2) {
					for (int ifr1 = ifr2; ifr1 < hd().qw.size(); ++ifr1) {
						String sT, units;
						if (ifr1 == ifr2)
							sT = "infinite";
						else {
							units = "(sec) ";
							sT = Format(" %12E", 2*M_PI/(hd().qw[ifr1] - hd().qw[ifr2]));
						}
						double T1 = 2*M_PI/(hd().qw[ifr1]);
						double T2 = 2*M_PI/(hd().qw[ifr2]);
						out << " ************************************************************************\n\n"
							   " 2nd-order period " << units << "= " << sT << "\n\n"
							   " Period indices:     " << Format("%2d   %2d", ifr1+1, ifr2+1) << "          Periods:  " 
							   		<< Format("%12E %12E", T1, T2) << "\n"
							   " ------------------------------------------------------------------------\n\n\n"
							   "DIFFERENCE-FREQUENCY EXCITING FORCES AND MOMENTS-DIRECT METHOD\n\n"
							   "  Heading indices:    " << Format("%2d   %2d", id1[ih], id2[ih]) << "  Headings (deg):       "
							   		<< Format("%4.1f", hd().qh[ih].real()) << "      " 
							   		<< Format("%4.1f", hd().qh[ih].imag()) << "\n\n\n"
							   "      I     Mod[F2(I)]     Pha[F2(I)]\n\n";
						for (int ib = 0; ib < hd().Nb; ++ib) {
							for (int idf = 0; idf < 6; ++idf) {
								static int idf12[] = {1, 3, 5, 2, 4, 6};
				        		int iidf = idf12[idf]-1;
								out << Format("    %2d   %12.6E     %10d\n", 1+iidf + 6*ib, 
										hd().F_ndim(abs(hd().qtfdif[ib][ih][iidf](ifr1, ifr2)), iidf),
     									int(ToDeg(arg(hd().qtfdif[ib][ih][iidf](ifr1, ifr2)))));
							}
						}
						out << "\n\n";
					}
				}
			}
		}     									     				
	} catch (Exc e) {
		BEM::PrintError(Format("\n%s: %s", t_("Error"), e));
		hd().lastError = e;
		return false;
	}
	return true;
}

void Wamit::Load_A(FileInLine &in, Eigen::MatrixXd &A) {
	in.GetLine(6);
	while (!in.IsEof()) {
		String line = TrimBoth(in.GetLine());
		if (line.IsEmpty())
           	break;
		FieldSplit f(in);
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

bool Wamit::Load_cfg(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	FieldSplit f(in);
 	f.IsSeparator = IsTabSpace;
 	
 	hd().description = in.GetLine();
 	
 	while (!in.IsEof()) {
		f.Load(in.GetLine());
		if (!f.IsEmpty() && f.GetText(0) == "IPEROUT") 
			iperout = f.GetInt(2);
 	}
 	return true;
}

bool Wamit::Load_pot(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	FieldSplit f(in);
 	f.IsSeparator = IsTabSpace;
 	
 	in.GetLine();
 	f.Load(in.GetLine());
 	
 	hd().h = f.GetDouble(0);
	if (hd().h < 0)
		hd().h = -1;

 	return true;
}

bool Wamit::Load_gdf(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	FieldSplit f(in);
 	f.IsSeparator = IsTabSpace;
 	
 	in.GetLine();
 	f.Load(in.GetLine());
 	
 	hd().len = f.GetDouble(0);
	hd().g = f.GetDouble(1);

 	return true;
}

static double w_iperout3(double KL, double g, double len) {
	return sqrt(KL*g/len);
}

static double w_iperout4(double nuL, double g, double len, double h) {
	Eigen::VectorXd x(1);
	x[0] = 1;
	if (!SolveNonLinearEquations(x, [&](const Eigen::VectorXd &x, Eigen::VectorXd &residual)->int {
		double w = x[0];
		double nu = nuL/len;
		residual[0] = nu*tanh(nu*h) - w*w/g;
		return 0;
	}))
		throw Exc(t_("Impossible to convert finite-depth wave number into frequency"));
	return x[0];
}
			
void Wamit::ProcessFirstColumn(UVector<double> &w, UVector<double> &T) {
	if (w.size() < 2)
		return;	
	if (IsNull(iperout)) {
		if (w[0] > w[1]) {
			hd().dataFromW = false;
			T = pick(w);
			w.SetCount(hd().Nf);	
		} else {
			hd().dataFromW = true;
			T.SetCount(hd().Nf);
		}
	} else {
		if (iperout == 1) {
			hd().dataFromW = false;
			T = pick(w);
			w.SetCount(hd().Nf);
		} else if (iperout == 2) {
			hd().dataFromW = true;
			T.SetCount(hd().Nf);
		} else if (iperout == 3) {
			hd().dataFromW = true;
			T.SetCount(hd().Nf);
			double g = Nvl(hd().g, hd().GetBEM().g);
			double len = Nvl(hd().len, hd().GetBEM().len);
			for (auto &ww : w)
				ww = w_iperout3(ww, g, len);
		} else {
			hd().dataFromW = true;
			T.SetCount(hd().Nf);
			double g = Nvl(hd().g, hd().GetBEM().g);
			double len = Nvl(hd().len, hd().GetBEM().len);
			if (IsNull(hd().h))
				throw Exc(t_("Wamit .1 file with finite water depth wavenumber requires .pot file"));
			for (auto &ww : w) 
				ww = w_iperout4(ww, g, len, hd().h);
		}
	}
	for (int ifr = 0; ifr < hd().Nf; ++ifr) {
		if (hd().dataFromW)
			T[ifr] = 2*M_PI/w[ifr];//fround(2*M_PI/w[ifr], 8);
		else
			w[ifr] = 2*M_PI/T[ifr];//fround(2*M_PI/T[ifr], 8);
	}
}

bool Wamit::Load_1(String fileName) {
	hd().dimen = false;
	hd().len = 1;
	
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	FieldSplit f(in);
 	f.IsSeparator = IsTabSpace;
 	
 	FileInLine::Pos fpos;
 	while (IsNull(ScanDouble(in.GetLine())) && !in.IsEof())
 		fpos = in.GetPos();
	if (in.IsEof())
		return false;
	
	UVector<double> w, T; 	
    
	in.SeekPos(fpos);
	
	int maxDof = 0;
	bool thereIsA0 = false, thereIsAinf = false; 
	while (!in.IsEof()) {
		f.GetLine();
		if (IsNull(f.GetDouble_nothrow(3))) {
			BEM::PrintWarning(S("\nWarning: ") + t_("Wrong data found before file end"));
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
	if (!IsNull(hd().Nb) && hd().Nb < Nb)
		throw Exc(in.Str() + "\n"  + Format(t_("Number of bodies loaded is lower than previous (%d != %d)"), hd().Nb, Nb));
	hd().Nb = Nb;
	if (hd().names.IsEmpty())
		hd().names.SetCount(hd().Nb);	
	
	int Nf = w.size();
	if (!IsNull(hd().Nf) && hd().Nf != Nf)
		throw Exc(in.Str() + "\n"  + Format(t_("Number of frequencies loaded is different than previous (%d != %d)"), hd().Nf, Nf));
	hd().Nf = Nf;
	
	if (hd().Nb == 0)// || hd().Nf < 2)
		throw Exc(in.Str() + "\n"  + Format(t_("Wrong format in Wamit file '%s'"), hd().file));
	
	UVector<double> sourcew = clone(w);
	
	ProcessFirstColumn(w, T);
	
	if (thereIsA0)
		hd().A0.setConstant(hd().Nb*6, hd().Nb*6, NaNDouble);
	if (thereIsAinf)
		hd().Ainf.setConstant(hd().Nb*6, hd().Nb*6, NaNDouble);

	if (Nf > 0) {
		hd().A.SetCount(6*hd().Nb);
		hd().B.SetCount(6*hd().Nb);
		for (int i = 0; i < 6*hd().Nb; ++i) {
			hd().A[i].SetCount(6*hd().Nb);
			hd().B[i].SetCount(6*hd().Nb);
			for (int j = 0; j < 6*hd().Nb; ++j) {
				hd().A[i][j].setConstant(hd().Nf, NaNDouble);	
				hd().B[i][j].setConstant(hd().Nf, NaNDouble);	
			}
		}
	}
	
	hd().names.SetCount(Nb);
	
	if (hd().w.IsEmpty()) {
		hd().w = pick(w);
		hd().T = pick(T);
	} else if (!CompareRatio(hd().w, w, 0.001))
		throw Exc(in.Str() + "\n"  + Format(t_("Frequencies loaded are different than previous\nPrevious: %s\nSeries:   %s"), ToString(hd().w), ToString(w)));
	else if (!CompareRatio(hd().T, T, 0.001))
		throw Exc(in.Str() + "\n"  + Format(t_("Periods loaded are different than previous\nPrevious: %s\nSeries:   %s"), ToString(hd().T), ToString(T)));
				
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
 		
 		if ((freq < 0)) {
 			if (!thereIsA0)
				throw Exc(in.Str() + "\n"  + t_("A[w=inf] is not expected"));
			hd().A0(i, j) = Aij;
		} else if (freq == 0) {
			if (!thereIsAinf)
				throw Exc(in.Str() + "\n"  + t_("A[w=0] is not expected"));				
			hd().Ainf(i, j) = Aij;
		} else {
			int ifr = FindRatio(sourcew, freq, 0.001);
			if (ifr < 0) {
				if (hd().dataFromW)
					throw Exc(in.Str() + "\n"  + Format(t_("Frequency %f is unknown"), freq));
				else 
					throw Exc(in.Str() + "\n"  + Format(t_("Period %f is unknown"), freq));
			}
		  	hd().A[i][j][ifr] = Aij;    
		  	hd().B[i][j][ifr] = f.GetDouble(4);   	
		}
	}
	hd().c0.setConstant(3, hd().Nb, 0);
	
	return true;	
}

bool Wamit::Load_hst(String fileName) {
	hd().dimen = false;
	if (IsNull(hd().len))
		hd().len = 1;
	
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	FieldSplit f(in);
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
	if (!IsNull(hd().Nb) && hd().Nb < Nb)
		throw Exc(in.Str() + "\n"  + Format(t_("Number of bodies loaded is lower than previous (%d != %d)"), hd().Nb, Nb));
	hd().Nb = Nb;
	if (hd().names.IsEmpty())
		hd().names.SetCount(hd().Nb);
	
	hd().C.SetCount(hd().Nb);
	for(int ib = 0; ib < hd().Nb; ++ib)
		hd().C[ib].setConstant(6, 6, 0);

	while (!in.IsEof()) {
		f.Load(in.GetLine());	
		int i = f.GetInt(0) - 1;
		int ib_i = i/6;
		i = i - ib_i*6;
		int j = f.GetInt(1) - 1;
		int ib_j = j/6;
		j = j - ib_j*6;
		if (ib_i == ib_j) 
			hd().C[ib_i](i, j) = f.GetDouble(2);
	}
		
	return true;
}

bool Wamit::Load_3(String fileName) {
	return Load_Forces(fileName, hd().ex);
}

bool Wamit::Load_Scattering(String fileName) {
	return Load_Forces(fileName, hd().sc);
}
		
bool Wamit::Load_FK(String fileName) {
	return Load_Forces(fileName, hd().fk);
}

bool Wamit::Load_4(String fileName) {
	return Load_Forces(fileName, hd().rao);
}

bool Wamit::Load_Forces(String fileName, Hydro::Forces &force) {
	hd().dimen = false;
	if (IsNull(hd().len))
		hd().len = 1;
	
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	FieldSplit f(in);
	f.IsSeparator = IsTabSpace;
 
 	FileInLine::Pos fpos;
 	while (IsNull(ScanDouble(in.GetLine())) && !in.IsEof())
 		fpos = in.GetPos();
	if (in.IsEof()) 
		return false;
	
	UVector<double> w, T; 	
    
	in.SeekPos(fpos);
	
	int maxDof = 0;
	hd().head.Clear();
	while (!in.IsEof()) {
		f.GetLine();
		if (IsNull(f.GetDouble_nothrow(3))) {
			BEM::PrintWarning(S("\nWarning: ") + t_("Wrong data found before file end"));
			break;
		}
		
		double freq = f.GetDouble(0);
		double head = f.GetDouble(1);	//FixHeading_180(f.GetDouble(1));
		FindAdd(w, freq);
		FindAdd(hd().head, head);
		
		int dof = f.GetInt(2);
		if (dof > maxDof)
			maxDof = dof-1;
	}
	Sort(hd().head);
	
	if (hd().head.size() == 0)
		throw Exc(in.Str() + "\n"  + Format(t_("Wrong format in Wamit file '%s'"), hd().file));
	
	if (!IsNull(hd().Nh) && hd().Nh != hd().head.size())
		throw Exc(in.Str() + "\n"  + Format(t_("Number of headings loaded is different than previous (%d != %d)"), hd().Nh, hd().head.size()));
	hd().Nh = hd().head.size();
	
	int Nb = 1 + int(maxDof/6);
	if (!IsNull(hd().Nb) && hd().Nb < Nb)
		throw Exc(in.Str() + "\n"  + Format(t_("Number of bodies loaded is lower than previous (%d != %d)"), hd().Nb, Nb));
	hd().Nb = Nb;
	if (hd().names.IsEmpty())
		hd().names.SetCount(hd().Nb);
		
	int Nf = w.size();
	if (!IsNull(hd().Nf) && hd().Nf != Nf)
		throw Exc(in.Str() + "\n"  + Format(t_("Number of frequencies loaded is different than previous (%d != %d)"), hd().Nf, Nf));
	hd().Nf = Nf;
	
	if (hd().Nb == 0 || hd().Nf < 2)
		throw Exc(in.Str() + "\n"  + Format(t_("Wrong format in Wamit file '%s'"), hd().file));
	
	
	UVector<double> sourcew = clone(w);
		
	ProcessFirstColumn(w, T);
	
	if (hd().w.IsEmpty()) {
		hd().w = pick(w);
		hd().T = pick(T);
	} else if (!CompareRatio(hd().w, w, 0.01))
		throw Exc(in.Str() + "\n"  + Format(t_("Frequencies loaded are different than previous\nPrevious: %s\nSeries:   %s"), ToString(hd().w), ToString(w)));
	else if (!CompareRatio(hd().T, T, 0.01))
		throw Exc(Format(t_("[%s] Periods loaded are different than previous\nPrevious: %s\nSeries:   %s"), ToString(hd().T), ToString(T)));
	
	if (hd().names.IsEmpty())
		hd().names.SetCount(hd().Nb);
	
	in.SeekPos(fpos);
	
	hd().Initialize_Forces(force);
		
	while (!in.IsEof()) {
		f.GetLine();
		if (IsNull(f.GetDouble_nothrow(3))) 
			break;
		
		double freq = f.GetDouble(0);
		int ifr = FindRatio(sourcew, freq, 0.001);
		if (ifr < 0) {
			if (hd().dataFromW)
				throw Exc(in.Str() + "\n"  + Format(t_("Frequency %f is unknown"), freq));
			else 
				throw Exc(in.Str() + "\n"  + Format(t_("Period %f is unknown"), freq));
		}		
		double head = f.GetDouble(1);	//FixHeading_180(f.GetDouble(1));
		int ih = FindRatio(hd().head, head, 0.001);
		if (ih < 0)
			throw Exc(in.Str() + "\n"  + Format(t_("Heading %f is unknown"), head));
			
		int i = f.GetInt(2) - 1;		
		
        force.force[ih](ifr, i) = std::complex<double>(f.GetDouble(5), f.GetDouble(6));
	}
	hd().c0.setConstant(3, hd().Nb, 0);
	
	return true;
}

bool Wamit::Load_12(String fileName, bool isSum, Function <bool(String, int)> Status) {
	hd().dimen = false;
	if (IsNull(hd().len))
		hd().len = 1;
	
	String ext = isSum ? ".12s" : ".12d";
	
	Status(Format("Loading %s base data", ext), 0);
	
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	FieldSplit f(in);
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
	
	hd().InitQTF(qtf, Nb, Nh, Nf);
	
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
		int ib = 0;
		while (idf > 5) {
			idf -= 6;
			ib++;
		}
		qtf[ib][ih][idf](ifr1, ifr2).real(f.GetDouble(7));
		qtf[ib][ih][idf](ifr1, ifr2).imag(f.GetDouble(8));
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

void Wamit::Save_1(String fileName, bool force_T) {
	if (!(hd().IsLoadedA() && hd().IsLoadedB())) 
		return;
		
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'"), fileName));
	
	if (hd().IsLoadedA0()) {
		for (int i = 0; i < hd().Nb*6; ++i)  
			for (int j = 0; j < hd().Nb*6; ++j) 
				out << Format(" %s %5d %5d %s\n", FormatWam(-1), i+1, j+1,
								FormatWam(Nvl2(hd().A0(i, j), hd().A0_ndim(i, j), 0.)));
	}
	if (hd().IsLoadedAinf()) {
		for (int i = 0; i < hd().Nb*6; ++i)  
			for (int j = 0; j < hd().Nb*6; ++j)
				out << Format(" %s %5d %5d %s\n", FormatWam(0), i+1, j+1,
								FormatWam(Nvl2(hd().Ainf(i, j), hd().Ainf_ndim(i, j), 0.)));
	}
	
	if (hd().Nf < 2)
		throw Exc(t_("No enough data to save (at least 2 frequencies)"));
		
	UVector<double> *pdata;
	if (force_T)
		pdata = &hd().T;
	else if (hd().dataFromW) 
		pdata = &hd().w;
	else
		pdata = &hd().T;
	UVector<double> &data = *pdata;
	
	int ifr0, ifrEnd, ifrDelta;
	bool growing = data[1] > data[0];
	if ((growing && (hd().dataFromW && !force_T)) || (!growing && !(hd().dataFromW && !force_T))) {
		ifr0 = 0;
		ifrEnd = hd().Nf;
		ifrDelta = 1;
	} else {
		ifr0 = hd().Nf - 1;
		ifrEnd = -1;
		ifrDelta = -1;
	}
	
	for (int ifr = ifr0; ifr != ifrEnd; ifr += ifrDelta)
		for (int i = 0; i < hd().Nb*6; ++i)  
			for (int j = 0; j < hd().Nb*6; ++j)
				out << Format(" %s %5d %5d %s %s\n", FormatWam(data[ifr]), i+1, j+1,
							 FormatWam(Nvl2(hd().A[i][j][ifr], hd().A_ndim(ifr, i, j), 0.)), 
							 FormatWam(Nvl2(hd().B[i][j][ifr], hd().B_ndim(ifr, i, j), 0.)));
}

void Wamit::Save_3(String fileName, bool force_T) {
	if (!hd().IsLoadedFex()) 
		return;
	
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'"), fileName));

	if (hd().Nf < 2)
		throw Exc(t_("No enough data to save (at least 2 frequencies)"));

	UVector<double> *pdata;
	int ifr0, ifrEnd, ifrDelta;
	if (hd().dataFromW && !force_T) 
		pdata = &hd().w;
	else
		pdata = &hd().T;
	UVector<double> &data = *pdata;
	
	if (((data[1] > data[0]) && (hd().dataFromW && !force_T)) || ((data[1] < data[0]) && !(hd().dataFromW && !force_T))) {
		ifr0 = 0;
		ifrEnd = hd().Nf;
		ifrDelta = 1;
	} else {
		ifr0 = hd().Nf - 1;
		ifrEnd = -1;
		ifrDelta = -1;
	}
	
	for (int ifr = ifr0; ifr != ifrEnd; ifr += ifrDelta)
		for (int ih = 0; ih < hd().Nh; ++ih)
			for (int i = 0; i < hd().Nb*6; ++i) {
				std::complex<double> &f = hd().ex.force[ih](ifr, i);
				std::complex<double> fn = hd().F_ndim(hd().ex, ih, ifr, i);
				out << Format(" %s %s %5d %s %s %s %s\n", 
					FormatWam(data[ifr]), FormatWam(hd().head[ih]), i+1,
					FormatWam(Nvl2(abs(f), abs(fn), 0.)), 
					FormatWam(Nvl2(arg(f), arg(fn)*180/M_PI, 0.)),
					FormatWam(Nvl2(f.real(), fn.real(), 0.)), 
					FormatWam(Nvl2(f.imag(), fn.imag(), 0.)));
			}
}

void Wamit::Save_hst(String fileName) {
	if (!hd().IsLoadedC()) 
		return;
		
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'"), fileName));

	for (int i = 0; i < 6*hd().Nb; ++i)  
		for (int j = 0; j < 6*hd().Nb; ++j) {
			int ib_i = i/6;
			int ii = i - ib_i*6;
			int ib_j = j/6;
			int jj = j - ib_j*6;
			out << Format(" %5d %5d  %s\n", i+1, j+1, 
					FormatWam(Nvl2(hd().C[ib_i](ii, jj), hd().C_ndim(ib_i, ii, jj), 0.)));
		}
}

// K is supposed to be dimensionalized
void Wamit::Save_hst_static(const Eigen::MatrixXd &C, String fileName, double rho, double g) {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'"), fileName));

	for (int i = 0; i < 6; ++i)  
		for (int j = 0; j < 6; ++j) 
			out << Format(" %5d %5d  %s\n", i+1, j+1, 
					FormatWam(Nvl2(C(i, j), C(i, j)/rho/g, 0.)));
}

void Wamit::Save_4(String fileName, bool force_T) {
	if (!hd().IsLoadedRAO()) 
		return;
		
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'"), fileName));
	
	if (hd().Nf < 2)
		throw Exc(t_("No enough data to save (at least 2 frequencies)"));
		
	UVector<double> *pdata;
	int ifr0, ifrEnd, ifrDelta;
	if (hd().dataFromW && !force_T) 
		pdata = &hd().w;
	else
		pdata = &hd().T;
	UVector<double> &data = *pdata;
	
	if (((data[1] > data[0]) && (hd().dataFromW && !force_T)) || ((data[1] < data[0]) && !(hd().dataFromW && !force_T))) {
		ifr0 = 0;
		ifrEnd = hd().Nf;
		ifrDelta = 1;
	} else {
		ifr0 = hd().Nf - 1;
		ifrEnd = -1;
		ifrDelta = -1;
	}
	
	for (int ifr = ifr0; ifr != ifrEnd; ifr += ifrDelta)
		for (int ih = 0; ih < hd().Nh; ++ih)
			for (int i = 0; i < hd().Nb*6; ++i) {
				std::complex<double> &f = hd().rao.force[ih](ifr, i);	
				std::complex<double> fn = hd().R_ndim(hd().rao, ih, ifr, i);
				out << Format(" %s %s %5d %s %s %s %s\n", 
					FormatWam(data[ifr]), FormatWam(hd().head[ih]), i+1,
					FormatWam(Nvl2(abs(f), abs(fn), 0.)), 
					FormatWam(Nvl2(arg(f), arg(fn), 0.)),
					FormatWam(Nvl2(f.real(), fn.real(), 0.)), 
					FormatWam(Nvl2(f.imag(), fn.imag(), 0.)));
			}
}
	
void Wamit::Save_12(String fileName, bool isSum, Function <bool(String, int)> Status,
					bool force_T, bool force_Deg, int qtfHeading) {
	if (!hd().IsLoadedQTF(isSum)) 
		return;
	
	String ext = isSum ? ".12s" : ".12d";
	fileName = ForceExt(fileName, ext);
	
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'"), fileName));
	
	int Nf = int(hd().qw.size());
	int Nh = int(hd().qh.size()); 
	
	if (Nf < 2)
		throw Exc(t_("Not enough data to save (at least 2 frequencies)"));
			
	UArray<UArray<UArray<MatrixXcd>>> &qtf = isSum ? hd().qtfsum : hd().qtfdif;		
			
	out << " WAMIT Numeric Output -- Filename  " << Format("%20<s", GetFileName(fileName)) << "  " << Format("%", GetSysTime()) << "\n";

	int num = int(pow2(hd().qw.size())*hd().qh.size());
	int inum = 0;
	for (int ifr1 = 0; ifr1 < hd().qw.size(); ++ifr1) 
		for (int ifr2 = 0; ifr2 < hd().qw.size(); ++ifr2) 
			for (int ih = 0; ih < hd().qh.size(); ++ih) {
				inum++;
				if (Status && num >= 20 && !(inum%(num/20)) && !Status(Format("Saving %s", fileName), (100*inum)/num))
					throw Exc(t_("Stop by user"));
				
				for (int ib = 0; ib < hd().Nb; ++ib)
        			for (int idf = 0; idf < 6; ++idf) {
        				static int idf12[] = {1, 3, 5, 2, 4, 6};
        				int iidf = idf12[idf]-1;
        				out << Format("   % 8.6E", 2*M_PI/hd().qw[ifr1]);
        				out << Format("   % 8.6E", 2*M_PI/hd().qw[ifr2]);
        				out << Format("   % 8.6E", hd().qh[ih].real());
        				out << Format("   % 8.6E", hd().qh[ih].imag());
        				out << Format("   %2d", ib*6 + iidf+1);
        				out << Format("   % 8.6E", hd().F_ndim(abs(qtf[ib][ih][iidf](ifr1, ifr2)), iidf));
        				out << Format("   % 8.6E", !force_Deg ? arg(qtf[ib][ih][iidf](ifr1, ifr2)) : ToDeg(arg(qtf[ib][ih][iidf](ifr1, ifr2))));
        				out << Format("   % 8.6E", hd().F_ndim(qtf[ib][ih][iidf](ifr1, ifr2).real(), iidf));
        				out << Format("   % 8.6E", hd().F_ndim(qtf[ib][ih][iidf](ifr1, ifr2).imag(), iidf));
        				out << "\n";
        			}
			}
}
