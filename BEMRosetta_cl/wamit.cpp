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
			String fileSC = ForceExt(file, ".3sc");
			BEM::Print("\n- " + Format(t_("Scattering file '%s'"), GetFileName(fileSC)));
			if (!Load_Scattering(fileSC))
				BEM::Print(S(": ** 3sc ") + t_("Not found") + "**");
			String fileFK = ForceExt(file, ".3fk");
			BEM::Print("\n- " + Format(t_("Froude-Krylov file '%s'"), GetFileName(fileFK)));
			if (!Load_FK(fileFK))
				BEM::Print(S(": ** 3fk ") + t_("Not found") + "**");
		} else if (S(".1.2.3.hst.4.12d.12s").Find(ext) >= 0) {
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

			if (IsNull(hd().Nh))
				hd().Nh = 0;
			if (IsNull(hd().Nf))
				hd().Nf = 0;

			String file12s = ForceExt(file, ".12s");
			BEM::Print("\n- " + Format(t_("Second order sum coefficients .12s file '%s'"), GetFileName(file12s)));
			if (!Load_12(file12s, true, Status))
				BEM::Print(S(": ** .12s ") + t_("Not found") + "**");
			
			String file12d = ForceExt(file, ".12d");
			BEM::Print("\n- " + Format(t_("Second order mean drift coefficients .12d file '%s'"), GetFileName(file12d)));
			if (!Load_12(file12d, false, Status))
				BEM::Print(S(": ** .12d ") + t_("Not found") + "**");
		}
		
		if (IsNull(hd().Nb)/* || IsNull(hd().Nh) || IsNull(hd().Nf) || hd().Nh == 0 || hd().Nf == 0*/)
			return false;
		
		hd().dof.Clear();	hd().dof.SetCount(hd().Nb, 6);
		//for (int i = 0; i < hd().Nb; ++i)
		//	hd().dof[i] = 6;
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
		if (hd().IsLoadedQTF()) {
			String fileQTFs = ForceExt(file, ".12s");
			BEM::Print("\n- " + Format(t_("QTF file '%s'"), GetFileName(fileQTFs)));
			Save_12(fileQTFs, true, Status, force_T, true, qtfHeading);
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
	
	//double xbody = 0, ybody = 0, zbody = 0;
	
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
			hd().cg.setConstant(3, hd().Nb, Null);
			hd().cb.setConstant(3, hd().Nb, Null);
			hd().c0.setConstant(3, hd().Nb, 0);
			hd().Vo.SetCount(hd().Nb, Null);
			hd().C.SetCount(hd().Nb);
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
			
			bool foundNh = false;
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
						FindAddDelta(hd().head, FixHeading_180(f.GetDouble(4)), 0.001);
					}
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
					hd().A[i][j].setConstant(hd().Nf, 0);		// In Wamit, unloaded DOFs are considered negligible	
					hd().B[i][j].setConstant(hd().Nf, 0);	
				}
			}			
			
			in.SeekPos(fpos);
			while (in.GetLine().Find("Wave period = infinite") < 0 && !in.IsEof())
				; 
			if (!in.IsEof()) {
				hd().A0.setConstant(hd().Nb*6, hd().Nb*6, 0);
				Load_A(in, hd().A0);
			}
			in.SeekPos(fpos);
			while (in.GetLine().Find("Wave period = zero") < 0 && !in.IsEof())
				; 
			if (!in.IsEof()) {
				hd().Ainf.setConstant(hd().Nb*6, hd().Nb*6, 0);
				Load_A(in, hd().Ainf);
			}
			
			in.SeekPos(fpos);
			
			int ifr = -1;
			while (!in.IsEof()) {
				line = in.GetLine();
				while ((line = in.GetLine()).Find("Wave period (sec)") < 0 && !in.IsEof())
					; 
				if (in.IsEof())
					return true;
				
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
	            		//if (hd().A.IsEmpty()) {
						//	hd().A.SetCount(hd().Nf);
						//	hd().B.SetCount(hd().Nf);
						//}
						in.GetLine(2);
						//if (hd().A.size() < hd().Nf)
			 			//	throw Exc(in.Str() + "\n"  + t_("A matrix is not dimensioned"));
						//if (hd().B.size() < hd().Nf)
			 			//	throw Exc(in.Str() + "\n"  + t_("B matrix is not dimensioned"));
		            
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
						if (hd().ex.ma.IsEmpty()) 
							hd().Initialize_Forces(hd().ex);
						
						int ih = 0;
						while (!in.IsEof()) {		
							line = in.GetLine();
							if (line.Find("Wave Heading (deg) :") >= 0) {
								f.Load(line);
								double head = FixHeading_180(f.GetDouble(4)); 
								int iih = FindDelta(hd().head, head, 0.001);
								if (iih >= 0) {
									in.GetLine(3); 
									while (!TrimBoth(line = in.GetLine()).IsEmpty()) {
										f.Load(line);
										double ma = f.GetDouble(1);
										double ph = ToRad(f.GetDouble(2));
										double re = ma*cos(ph);
										double im = ma*sin(ph);
										int i = abs(f.GetInt(0)) - 1;
										if (OUTB(ih, hd().Nh) || OUTB(ifr, hd().Nf) || OUTB(i, hd().Nb*6))
											throw Exc(in.Str() + "\n"  + Format(t_("Index [%d](%d, %d) out of bounds"), ih, ifr, i));
										hd().ex.ma[ih](ifr, i) = ma;	
										hd().ex.ph[ih](ifr, i) = ph;	
										hd().ex.re[ih](ifr, i) = re;	
										hd().ex.im[ih](ifr, i) = im;	
									}
								}
								ih++;
								if (ih >= hd().Nh)
									break;
							}
						} 
					} else if (line.Find("RESPONSE AMPLITUDE OPERATORS") >= 0) {
						if (hd().rao.ma.IsEmpty()) 
							hd().Initialize_Forces(hd().rao);
						
						int ih = 0;
						while (!in.IsEof()) {		
							line = in.GetLine();
							if (line.Find("Wave Heading (deg) :") >= 0) {
								f.Load(line);
								double head = FixHeading_180(f.GetDouble(4)); 
								int iih = FindDelta(hd().head, head, 0.001);
								if (iih >= 0) {
									in.GetLine(3); 
									while (!TrimBoth(line = in.GetLine()).IsEmpty()) {
										f.Load(line);
										double ma = f.GetDouble(1);
										double ph = ToRad(f.GetDouble(2));
										double re = ma*cos(ph);
										double im = ma*sin(ph);
										int i = abs(f.GetInt(0)) - 1;
										if (OUTB(ih, hd().Nh) || OUTB(ifr, hd().Nf) || OUTB(i, hd().Nb*6))
											throw Exc(in.Str() + "\n"  + Format(t_("Index [%d](%d, %d) out of bounds"), ih, ifr, i));
										hd().rao.ma[ih](ifr, i) = ma;	
										hd().rao.ph[ih](ifr, i) = ph;	
										hd().rao.re[ih](ifr, i) = re;	
										hd().rao.im[ih](ifr, i) = im;	
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
		}
	}

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
			if (!IsNull(base(r, c))) 
				out << Format("%6>d%6>d %E\n", r+1, c+1, fun(r, c));
	out << "\n\n";
}

void Wamit::Save_AB(FileOut &out, int ifr) {
	out <<	"    ADDED-MASS AND DAMPING COEFFICIENTS\n"
			"     I     J         A(I,J)         B(I,J)\n\n";
	for (int r = 0; r < hd().Nb*6; ++r) 
		for (int c = 0; c < hd().Nb*6; ++c) 
			if (!IsNull(hd().A[r][c][ifr]) && !IsNull(hd().B[r][c][ifr]))
				out << Format("%6>d%6>d %E %E\n", r+1, c+1, hd().A_ndim(ifr, r, c), hd().B_ndim(ifr, r, c));
	out << "\n\n\n\n";
}

void Wamit::Save_Forces(FileOut &out, int ifr) {
	out <<	"    DIFFRACTION EXCITING FORCES AND MOMENTS\n\n";
	for (int ih = 0; ih < hd().Nh; ++ih) {
		out << "  Wave Heading (deg) :      " << hd().head[ih] << "\n\n"
			<< "     I     Mod[Xh(I)]     Pha[Xh(I)]\n\n";
		for (int i = 0; i < hd().ex.ma[ih].cols(); ++i)
			if (!IsNull(hd().ex.ma[ih](ifr, i))) 
				out << Format(" %7>d   %E   %f\n", i+1, hd().F_ma_ndim(hd().ex, ih, ifr, i), ToDeg(hd().ex.ph[ih](ifr, i)));
		out << "\n\n\n\n";
	}
}

void Wamit::Save_RAO(FileOut &out, int ifr) {
	out <<	"    RESPONSE AMPLITUDE OPERATORS\n\n";
	for (int ih = 0; ih < hd().Nh; ++ih) {
		out << "  Wave Heading (deg) :      " << hd().head[ih] << "\n\n"
			<< "     I     Mod[Xh(I)]     Pha[Xh(I)]\n\n";
		for (int i = 0; i < hd().rao.ma[ih].cols(); ++i)
			if (!IsNull(hd().rao.ma[ih](ifr, i))) 
				out << Format(" %7>d   %E   %f\n", i+1, hd().F_ma_ndim(hd().rao, ih, ifr, i), ToDeg(hd().rao.ph[ih](ifr, i)));
		out << "\n\n\n\n";
	}
}

bool Wamit::Save_out(String file, double g, double rho) {
	try {
		FileOut out(file);
		if (!out.IsOpen())
			throw Exc(Format(t_("Impossible to open '%s'"), file));

		out << Format(" %s\n\n", String('-', 72))
		    <<        "                   BEMRosetta generated .out format\n\n"
		    << Format(" %s\n\n\n", String('-', 72))
			<< " Low-order panel method  (ILOWHI=0)\n\n";
		if (hd().Nb == 1)
			out << " Input from Geometric Data File:         unknown.gdf\n"
				<< " Unknown gdf file source\n\n";
		else {
			out << " Input from Geometric Data Files:\n";
			for (int ib = 0; ib < hd().Nb; ++ib)
				out << Format("                               N=  %d     unknown%d.gdf\n"
							  " Unknown gdf file source\n\n", ib+1, ib);
		}
		out	<< " Input from Potential Control File:      unknown.pot\n"
			<< " unknown.pot -- file type .gdf, ILOWHI=0, IRR=1\n\n\n"
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
			out	<< " XBODY =    0.0000 YBODY =    0.0000 ZBODY =    0.0000 PHIBODY =   0.0\n";
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
			out	<< " Hydrostatic and gravitational restoring coefficients:\n"; 
			out	<< " C(3,3),C(3,4),C(3,5): " << Format("%s %s %s\n", 
						FormatWam(hd().C[ibody](2, 2)), FormatWam(hd().C[ibody](2, 3)), FormatWam(hd().C[ibody](2, 4)));
			out	<< " C(4,4),C(4,5),C(4,6):               " << Format("%s %s %s\n", 
						FormatWam(hd().C[ibody](3, 3)), FormatWam(hd().C[ibody](3, 4)), FormatWam(hd().C[ibody](3, 5)));
			out	<< "        C(5,5),C(5,6):                             " << Format("%s %s\n", 
						FormatWam(hd().C[ibody](4, 4)), FormatWam(hd().C[ibody](4, 5)));
			double cgx = 0, cgy = 0, cgz = 0;
			if (hd().cb.size() > 0) {
				cgx = hd().cg(0, ibody);
				cgy = hd().cg(1, ibody);
				cgz = hd().cg(2, ibody);
			}
			out	<< Format(" Center of Gravity  (Xg,Yg,Zg): %s %s %s\n", 
						FormatWam(cgx), FormatWam(cgy), FormatWam(cgz));
			out	<< " Radii of gyration:     0.000000     0.000000     0.000000\n"
      			   "                        0.000000     0.000000     0.000000\n"
                   "                        0.000000     0.000000     0.000000\n\n";
		}
		out << " ------------------------------------------------------------------------\n"
        	<< "                    Output from  WAMIT\n"
			<< " ------------------------------------------------------------------------\n\n";
		
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

bool Wamit::Load_Scattering(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	FieldSplit f(in);
	f.IsSeparator = IsTabSpace;
	
	hd().Initialize_Forces(hd().sc);
	
	in.GetLine();
    for (int ifr = 0; ifr < hd().Nf; ++ifr) {
        for (int ih = 0; ih < hd().Nh; ++ih) {
            for (int k = 0; k < hd().Nb*6; ++k) {		// Number of non-zero dof
        		f.Load(in.GetLine());
        		
        		int i = f.GetInt(2) - 1;
        		if (OUTB(i, hd().Nb*6))
        			throw Exc(in.Str() + "\n"  + Format(t_("Index (%d) out of bounds"), i));
                hd().sc.ma[ih](ifr, i) = f.GetDouble(3);
                hd().sc.ph[ih](ifr, i) = f.GetDouble(4)*M_PI/180;
                hd().sc.re[ih](ifr, i) = f.GetDouble(5);
                hd().sc.im[ih](ifr, i) = f.GetDouble(6);
            }
        }
    }
    return true;
}
		
bool Wamit::Load_FK(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	FieldSplit f(in);
	f.IsSeparator = IsTabSpace;
	
	hd().Initialize_Forces(hd().fk);
	
	in.GetLine();	
    for (int ifr = 0; ifr < hd().Nf; ++ifr) {
        for (int ih = 0; ih < hd().Nh; ++ih) {
            for (int k = 0; k < hd().Nb*6; ++k) {		// Number of non-zero dof
                f.Load(in.GetLine());
        		
        		int i = f.GetInt(2) - 1;
        		if (OUTB(i, hd().Nb*6))
        			throw Exc(in.Str() + "\n"  + Format(t_("Index (%d) out of bounds"), i));
                hd().fk.ma[ih](ifr, i) = f.GetDouble(3);
                hd().fk.ph[ih](ifr, i) = f.GetDouble(4)*M_PI/180;
                hd().fk.re[ih](ifr, i) = f.GetDouble(5);
                hd().fk.im[ih](ifr, i) = f.GetDouble(6);
            }
        }
    }	
    return true;
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
		hd().A0.setConstant(hd().Nb*6, hd().Nb*6, 0);
	if (thereIsAinf)
		hd().Ainf.setConstant(hd().Nb*6, hd().Nb*6, 0);

	if (Nf > 0) {
		hd().A.SetCount(6*hd().Nb);
		hd().B.SetCount(6*hd().Nb);
		for (int i = 0; i < 6*hd().Nb; ++i) {
			hd().A[i].SetCount(6*hd().Nb);
			hd().B[i].SetCount(6*hd().Nb);
			for (int j = 0; j < 6*hd().Nb; ++j) {
				hd().A[i][j].setConstant(hd().Nf, 0);	
				hd().B[i][j].setConstant(hd().Nf, 0);	
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

bool Wamit::Load_3(String fileName) {
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
		double head = FixHeading_180(f.GetDouble(1));
		FindAdd(w, freq);
		FindAdd(hd().head, head);
		
		int dof = f.GetInt(2);
		if (dof > maxDof)
			maxDof = dof-1;
	}
	Sort(hd().head);
	
	if (hd().head.size() == 0)
		throw Exc(in.Str() + "\n" + Format(t_("Wrong format in Wamit file '%s'"), hd().file));
		
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
	} else if (!CompareRatio(hd().w, w, 0.001))
		throw Exc(in.Str() + "\n"  + Format(t_("Frequencies loaded are different than previous\nPrevious: %s\nSeries:   %s"), ToString(hd().w), ToString(w)));
	else if (!CompareRatio(hd().T, T, 0.001))
		throw Exc(in.Str() + "\n"  + Format(t_("Periods loaded are different than previous\nPrevious: %s\nSeries:   %s"), ToString(hd().T), ToString(T)));
	
	in.SeekPos(fpos);
	
	hd().Initialize_Forces(hd().ex);
	
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
		double head = FixHeading_180(f.GetDouble(1));
		int ih = FindRatio(hd().head, head, 0.001);
		if (ih < 0)
			throw Exc(in.Str() + "\n"  + Format(t_("Heading %f is unknown"), head));
			
		int i = f.GetInt(2) - 1;		
		
       	hd().ex.ma[ih](ifr, i) = f.GetDouble(3);
     	hd().ex.ph[ih](ifr, i) = f.GetDouble(4)*M_PI/180;
        hd().ex.re[ih](ifr, i) = f.GetDouble(5);
        hd().ex.im[ih](ifr, i) = f.GetDouble(6);
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

bool Wamit::Load_4(String fileName) {
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
		double head = FixHeading_180(f.GetDouble(1));
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
	
	hd().Initialize_RAO();
		
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
		double head = FixHeading_180(f.GetDouble(1));
		int ih = FindRatio(hd().head, head, 0.001);
		if (ih < 0)
			throw Exc(in.Str() + "\n"  + Format(t_("Heading %f is unknown"), head));
			
		int i = f.GetInt(2) - 1;		
		
       	hd().rao.ma[ih](ifr, i) = f.GetDouble(3);
     	hd().rao.ph[ih](ifr, i) = f.GetDouble(4)*M_PI/180;
        hd().rao.re[ih](ifr, i) = f.GetDouble(5);
        hd().rao.im[ih](ifr, i) = f.GetDouble(6);
	}
	hd().c0.setConstant(3, hd().Nb, 0);
	
	return true;
}

bool Wamit::Load_12(String fileName, bool isSum, Function <bool(String, int)> Status) {
	hd().dimen = false;
	if (IsNull(hd().len))
		hd().len = 1;
	
	String ext = isSum ? ".12s" : ".12d";
	
	Status(Format("Loading %s base data", ext), -1);
	
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	int nrows = hd().Nh*hd().Nf*hd().Nf;
	UArray<Hydro::QTF> &qtfList = isSum ? hd().qtfsum : hd().qtfdif;
	
	qtfList.Clear();
	qtfList.Reserve(hd().Nb*nrows);
		
	FieldSplit f(in);
	f.IsSeparator = IsTabSpace;
	
	FileInLine::Pos fpos = in.GetPos();
	f.GetLine();		
	if (!IsNull(f.GetDouble_nothrow(0)))
		in.SeekPos(fpos);		// No header, rewind
	else
		fpos = in.GetPos();		// Avoid header
	
    UVector<double> w, T;
    UVector<double> head;
    
    bool isRad = true;
	int maxDof = 0;	
	while (!in.IsEof()) {
		f.GetLine();
		double freq = f.GetDouble(0);
		double hd = FixHeading_180(f.GetDouble(2));
		FindAdd(w, freq);
		FindAdd(head, hd);
		
		int dof = f.GetInt(4);
		if (dof > maxDof)
			maxDof = dof-1;
		
		if (abs(f.GetDouble(6)) > M_PI + 0.01)
			isRad = false;
	}
	Sort(head);
	
	int Nb = 1 + int(maxDof/6);
	if (!IsNull(hd().Nb) && hd().Nb < Nb)
		throw Exc(Format(t_("Number of bodies loaded is lower than previous (%d != %d)"), hd().Nb, Nb));
	hd().Nb = Nb;	
	if (hd().names.IsEmpty())
		hd().names.SetCount(hd().Nb);
		
	int Nh = head.size();
	if (Nh == 0)
		throw Exc(Format(t_("Wrong format in Wamit file '%s'. No headings found"), hd().file));
	if (hd().qtfhead.IsEmpty()) 
		hd().qtfhead = pick(head);
	else {
		for(int i = 0; i < head.size(); i++) {
			bool found = false;
			for(int j = 0; j < hd().qtfhead.size(); j++) {
				if (EqualRatio(head[i], hd().qtfhead[j], 0.001, 0.001)) {
					found = true;
					break;
				}
			}
			if (!found)
				throw Exc(Format(t_("Heading %f not found in previous %s"), head[i], ToString(hd().head)));
		}
	}

	int Nf = w.size();
	
	if (w[0] > w[1]) {
		hd().qtfdataFromW = false;
		T = pick(w);
		w.SetCount(Nf);	
	} else {
		hd().qtfdataFromW = true;
		T.SetCount(Nf);
	}
	
	for (int ifr = 0; ifr < Nf; ++ifr) {
		if (hd().qtfdataFromW)
			T[ifr] = 2*M_PI/w[ifr];//fround(2*M_PI/w[ifr], 8);
		else
			w[ifr] = 2*M_PI/T[ifr];//fround(2*M_PI/T[ifr], 8);
	}
	
	if (hd().qtfw.IsEmpty()) {
		hd().qtfw = pick(w);
		hd().qtfT = pick(T);
	} else if (!CompareRatio(hd().qtfw, w, 0.01))
		throw Exc(in.Str() + "\n"  + Format(t_("Frequencies loaded are different than previous\nPrevious: %s\nSeries:   %s"), ToString(hd().w), ToString(w)));
	else if (!CompareRatio(hd().qtfT, T, 0.01))
		throw Exc(Format(t_("[%s] Periods loaded are different than previous\nPrevious: %s\nSeries:   %s"), ToString(hd().T), ToString(T)));
	
	in.SeekPos(fpos);
	
	int total = Nb*6*pow2(hd().qtfw.size())*pow2(hd().qtfhead.size());
	int it = 0;
	while (!in.IsEof()) {
		it++;
		f.Load(in.GetLine());
	
		double wT1 = f.GetDouble(0);
		int iwT1;
		if (hd().qtfdataFromW) {
			iwT1 = FindClosest(hd().qtfw, wT1);
			if (!EqualRatio(hd().qtfw[iwT1], wT1, 0.01, 0.001))
				throw Exc(in.Str() + "\n"  + Format(t_("Frequency 1 %f not found"), wT1));
		} else {
			iwT1 = FindClosest(hd().qtfT, wT1);
			if (!EqualRatio(hd().qtfT[iwT1], wT1, 0.01, 0.001))
				throw Exc(in.Str() + "\n"  + Format(t_("Period 1 %f not found"), wT1));
		}
		double wT2 = f.GetDouble(1);
		int iwT2;
		if (hd().qtfdataFromW) {
			iwT2 = FindClosest(hd().qtfw, wT2);
			if (!EqualRatio(hd().qtfw[iwT2], wT2, 0.01, 0.001))
				throw Exc(in.Str() + "\n"  + Format(t_("Frequency 2 %f not found"), wT2));
		} else {
			iwT2 = FindClosest(hd().qtfT, wT2);
			if (!EqualRatio(hd().qtfT[iwT2], wT2, 0.01, 0.001))
				throw Exc(in.Str() + "\n"  + Format(t_("Period 2 %f not found"), wT2));
		}
		double h1 = FixHeading_180(f.GetDouble(2));
		int ih1 = FindClosest(hd().qtfhead, h1);
		if (!EqualRatio(hd().qtfhead[ih1], h1, 0.01, 0.001))
			throw Exc(in.Str() + "\n"  + Format(t_("Heading 1 %f not found"), h1));
		double h2 = FixHeading_180(f.GetDouble(3));
		int ih2 = FindClosest(hd().qtfhead, h2);
		if (!EqualRatio(hd().qtfhead[ih2], h2, 0.01, 0.001))
			throw Exc(in.Str() + "\n"  + Format(t_("Heading 2 %f not found"), h2));

		int idof = f.GetInt(4)-1;
		if (idof < 0 || idof > hd().Nb*6 -1)
			throw Exc(in.Str() + "\n"  + Format(t_("Wrong DOF id %d"), idof+1));
	    
		int ib = int(idof/6);
		idof = idof - ib*6;
		
		Hydro::QTF *pqtf;
		int id = 0;
		if ((id = hd().GetQTFId(id, qtfList, hd().qtfCases, ib, ih1, ih2, iwT1, iwT2)) >= 0)
			pqtf = &qtfList[id];
		else {
			pqtf = &qtfList.Add();
			pqtf->Set(ib, ih1, ih2, iwT1, iwT2);
		}
	    
	    Hydro::QTF &qtf = *pqtf;
	    
	    qtf.fma[idof] = f.GetDouble(5);
	    double ph = f.GetDouble(6);
	    qtf.fph[idof] = isRad ? ph : ToRad(ph);
	    qtf.fre[idof] = f.GetDouble(7);
	    qtf.fim[idof] = f.GetDouble(8);    
	        
		double fre = qtf.fma[idof]*cos(qtf.fph[idof]);
		double fim = qtf.fma[idof]*sin(qtf.fph[idof]);
		
		if (abs(fre - qtf.fre[idof]) > 0.1)  
			throw Exc(in.Str() + "\n"  + Format(t_("Real force %f does not match with magnitude %f and phase %f (%f)"), 
										qtf.fre[idof], qtf.fma[idof], ph, fre));
		if (abs(fim - qtf.fim[idof]) > 0.1)  
			throw Exc(in.Str() + "\n"  + Format(t_("Imaginary force %f does not match with magnitude %f and phase %f (%f)"), 
										qtf.fim[idof], qtf.fma[idof], ph, fim));
		if (Status && !(it%500) && !Status(Format("Loading %s %d/%d", ext, it, total), 100*it/total))
			throw Exc(t_("Stop by user"));
	}
	hd().GetQTFList(qtfList, hd().qtfCases, hd().qtfhead);
	
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
			for (int i = 0; i < hd().Nb*6; ++i)
				out << Format(" %s %s %5d %s %s %s %s\n", 
						FormatWam(data[ifr]), FormatWam(hd().head[ih]), i+1,
			FormatWam(Nvl2(hd().ex.ma[ih](ifr, i), hd().F_ma_ndim(hd().ex, ih, ifr, i), 0.)), 
			FormatWam(Nvl2(hd().ex.ph[ih](ifr, i), hd().ex.ph[ih](ifr, i)*180/M_PI, 0.)),
			FormatWam(Nvl2(hd().ex.re[ih](ifr, i), hd().F_re_ndim(hd().ex, ih, ifr, i), 0.)), 
			FormatWam(Nvl2(hd().ex.im[ih](ifr, i), hd().F_im_ndim(hd().ex, ih, ifr, i), 0.)));
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
			for (int i = 0; i < hd().Nb*6; ++i)
				out << Format(" %s %s %5d %s %s %s %s\n", 
					FormatWam(data[ifr]), FormatWam(hd().head[ih]), i+1,
					FormatWam(Nvl2(hd().rao.ma[ih](ifr, i), hd().R_ma_ndim(hd().rao, ih, ifr, i), 0.)), 
					FormatWam(Nvl2(hd().rao.ph[ih](ifr, i), hd().rao.ph[ih](ifr, i)*180/M_PI, 0.)),
					FormatWam(Nvl2(hd().rao.re[ih](ifr, i), hd().R_re_ndim(hd().rao, ih, ifr, i), 0.)), 
					FormatWam(Nvl2(hd().rao.im[ih](ifr, i), hd().R_im_ndim(hd().rao, ih, ifr, i), 0.)));
}
	
void Wamit::Save_12(String fileName, bool isSum, Function <bool(String, int)> Status,
					bool force_T, bool force_Deg, int qtfHeading) {
	if (!hd().IsLoadedQTF()) 
		return;
	
	String ext = isSum ? ".12s" : ".12d";
	fileName = ForceExt(fileName, ext);
	
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'"), fileName));
	
	int Nf = hd().qtfw.size();
	int Nh = hd().qtfhead.size(); 
	
	if (Nf < 2)
		throw Exc(t_("Not enough data to save (at least 2 frequencies)"));
			
	UArray<Hydro::QTF> &qtfList = isSum ? hd().qtfsum : hd().qtfdif;
		
	out << " WAMIT Numeric Output -- Filename  " << Format("%20<s", GetFileName(fileName)) << "  " << Format("%", GetSysTime()) << "\n";
	
	UVector<double> *pdata;
	int ifr0, ifrEnd, ifrDelta;
	if (hd().qtfdataFromW && !force_T) 
		pdata = &hd().qtfw;
	else
		pdata = &hd().qtfT;
	UVector<double> &data = *pdata;
	
	if (((data[1] > data[0]) && (hd().dataFromW && !force_T)) || ((data[1] < data[0]) && !(hd().dataFromW && !force_T))) {
		ifr0 = 0;
		ifrEnd = hd().qtfw.size();
		ifrDelta = 1;
	} else {
		ifr0 = hd().qtfw.size() - 1;
		ifrEnd = -1;
		ifrDelta = -1;
	}
	
	static int idf12[] = {1, 3, 5, 2, 4, 6};
	
	int num = hd().qtfw.size();
	int it = 0;
	for (int ifr1 = ifr0; ifr1 != ifrEnd; ifr1 += ifrDelta) {
		for (int ifr2 = ifr0; ifr2 != ifrEnd; ifr2 += ifrDelta) {
			if (Status && !(it%500) && !Status(Format("Saving %s %d/%d", ext, it, num*num), 100*it/(num*num)))
				throw Exc(t_("Stop by user"));
			it++;
			int id = 0;
			for (int ih1 = 0; ih1 < Nh; ++ih1) {
				if (qtfHeading >= 0 && ih1 != qtfHeading)
					continue;
				for (int ih2 = 0; ih2 < Nh; ++ih2) {
					if (qtfHeading >= 0 && ih1 != qtfHeading)
						continue;	 
					for (int ib = 0; ib < hd().Nb; ++ib) {
						id = hd().GetQTFId(id, qtfList, hd().qtfCases, ib, ih1, ih2, ifr1, ifr2);
						if (id >= 0 || qtfHeading == -2) {
							double qtfhead1, qtfhead2;
							if (qtfHeading >= 0) 
								qtfhead1 = qtfhead2 = 0;		// Just 0º
							else {
								qtfhead1 = hd().qtfhead[ih1];
								qtfhead2 = hd().qtfhead[ih2];
							}							
							for (int _idf = 0; _idf < 6; ++_idf) {
								int idf = idf12[_idf]-1;
								out << Format("   %s   %s   %s   %s   %2d", 
										FormatWam(data[ifr1]),
										FormatWam(data[ifr2]), 
										FormatWam(qtfhead1),
										FormatWam(qtfhead2),
										ib*6 + idf + 1);
								if ((qtfHeading >= 0 || id >= 0) /*&& idf < 2*/) {	/*// To remove dof */
									Hydro::QTF &qtf = qtfList[id];
									out << Format("   %s   %s   %s   %s\n",		
											FormatWam(hd().F_ndim(qtf.fma[idf], idf)), 
											FormatWam(!force_Deg ? qtf.fph[idf] : ToDeg(qtf.fph[idf])),
											FormatWam(hd().F_ndim(qtf.fre[idf], idf)), 
											FormatWam(hd().F_ndim(qtf.fim[idf], idf)));
								} else // qtfHeading == -2
									out << Format("   %s   %s   %s   %s\n",		
											FormatWam(0), 
											FormatWam(0),
											FormatWam(0), 
											FormatWam(0));
							}
						}
					}
				}
			}
		}
	}
}
