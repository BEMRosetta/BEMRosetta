// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2021, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"

const char *textDOF[] = {"X", "Y", "Z", "RX", "RY", "RZ"};

bool Aqwa::Load(String file, double) {
	hd().file = file;
	hd().name = GetFileTitle(file);
	hd().dimen = true;
	hd().len = 1;
	hd().code = Hydro::AQWA;
	hd().Nb = Null;
	
	try {
		BEM::Print("\n\n" + Format(t_("Loading '%s'"), file));

		BEM::Print("\n- " + S(t_("LIS file")));
		if (!Load_LIS()) {
			BEM::PrintWarning(S(": ** LIS file ") + t_("Not found") + "**");
			BEM::Print("\n- " + S(t_("AH1 file")));
			if (!Load_AH1()) {
				BEM::PrintWarning(S(": ** AH1 file ") + t_("Not found") + "**");
				hd().Nh = hd().Nf = 0;
			//	throw Exc(t_("No .AH1 or .LIS file found"));
			}
		}
		//if (IsNull(hd().Nb))
		//	return false;
		
		BEM::Print("\n- " + S(t_("QTF file")));
		if (!Load_QTF()) 
			BEM::Print(S(": ** QTF file ") + t_("Not found") + "**");
		
		if (IsNull(hd().Nb))
			return false;
	
		hd().dof.Clear();	hd().dof.SetCount(hd().Nb, 0);
		for (int i = 0; i < hd().Nb; ++i)
			hd().dof[i] = 6;
	} catch (Exc e) {
		BEM::PrintError(Format("\n%s: %s", t_("Error"), e));
		hd().lastError = e;
		return false;
	}
	
	return true;
}

bool Aqwa::Load_AH1() {
	String fileName = ForceExt(hd().file, ".AH1");
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;

	hd().Nb = hd().Nh = hd().Nf = Null;
	hd().head.Clear();
	hd().w.Clear();
	hd().rho = hd().g = hd().h = Null;
	hd().dataFromW = true;
	
	String line;
	FieldSplit f(in);
	f.IsSeparator = IsTabSpace;
	while(!in.IsEof()) {
		line = in.GetLine();
		
		if (TrimBoth(line).StartsWith("*"))
			break;
	}
	f.Load(in.GetLine());
	hd().Nb = f.GetInt(0);
	hd().Nh = f.GetInt(1);
	hd().Nf = f.GetInt(2);
	
	if (hd().Nh != f.size() - 3)
		throw Exc(in.Str() + "\n"  + Format(t_("Number of headings do not match %d<>%d"), hd().Nh, f.size() - 3));
	for (int i = 3; i < f.size(); ++i)
		hd().head << f.GetDouble(i);
	
	hd().names.SetCount(hd().Nb);
	hd().cg.setConstant(3, hd().Nb, Null);
	hd().c0.setConstant(3, hd().Nb, Null);
	hd().C.SetCount(hd().Nb);
	for (int ib = 0; ib < hd().Nb; ++ib) 
		hd().C[ib].setConstant(6, 6, Null); 
	hd().A.SetCount(6*hd().Nb);
	hd().B.SetCount(6*hd().Nb);
	for (int i = 0; i < 6*hd().Nb; ++i) {
		hd().A[i].SetCount(6*hd().Nb);
		hd().B[i].SetCount(6*hd().Nb);
		for (int j = 0; j < 6*hd().Nb; ++j) {
			hd().A[i][j].setConstant(hd().Nf, Null);	
			hd().B[i][j].setConstant(hd().Nf, Null);	
		}
	}
	hd().Initialize_Forces(hd().ex);
	
	while(!in.IsEof() && hd().w.size() < hd().Nf) {
		line = TrimBoth(in.GetLine());
		f.Load(line);
		
		for (int i = 0; i < f.size(); ++i) {
			double w = f.GetDouble(i);
			hd().w << w;
			hd().T << 2*M_PI/w;
		}
	}	
	if (hd().Nf != hd().w.size())
		throw Exc(in.Str() + "\n"  + Format(t_("Number of frequencies do not match %d<>%d"), hd().Nf, hd().w.size()));

	while(!in.IsEof()) {
		line = TrimBoth(in.GetLine());
		if (line.StartsWith("GENERAL")) {
			f.Load(in.GetLine());
			hd().h   = f.GetDouble(0);
			hd().rho = f.GetDouble(1);
			hd().g   = f.GetDouble(2);
		} else if (line.StartsWith("COG")) {
			for (int i = 0; i < hd().Nb; ++i) {
				f.Load(in.GetLine());
				int ib = f.GetInt(0) - 1;
				if (ib < 0 || ib >= hd().Nb)
					throw Exc(in.Str() + "\n"  + Format(t_("Unknown body found in COG %d (%d)"), ib+1, hd().Nb));
				
				hd().cg(0, ib) = f.GetDouble(1);
				hd().cg(1, ib) = f.GetDouble(2);
				hd().cg(2, ib) = f.GetDouble(3);
				hd().c0 = clone(hd().cg);
			}
		} else if (line.StartsWith("HYDSTIFFNESS")) {
			for (int ib = 0; ib < hd().Nb; ++ib) {
				hd().C[ib].setConstant(6, 6, 0);
	            for (int idf = 0; idf < 6; ++idf) {
	                f.Load(in.GetLine());
	                int did = 0;
	                if (idf == 0) {
	                    did = 1;
	                    int itb = f.GetInt(0);
	                    if (itb - 1 != ib)
	                        throw Exc(in.Str() + "\n"  + Format(t_("Body # does not match in 'HYDSTIFFNESS' %d<>%d"), itb, ib+1));
	                }
	                for (int jdf = 0; jdf < 6; ++jdf) 
	                    hd().C[ib](idf, jdf) = f.GetDouble(jdf + did);
	            }
			}
		} else if (line.StartsWith("ADDEDMASS") || line.StartsWith("DAMPING")) {
			bool am = line.StartsWith("ADDEDMASS");
			String sarea = am ? "ADDEDMASS" : "DAMPING";
			for (int ib0 = 0; ib0 < hd().Nb; ++ib0) {
				for (int ib1 = 0; ib1 < hd().Nb; ++ib1) {
					for (int ifr = 0; ifr < hd().Nf; ++ifr) {
			            for (int idf = 0; idf < 6; ++idf) {
			                f.Load(in.GetLine());
			                int did = 0;
			                if (idf == 0) {
			                    did = 3;
			                    int itb0 = f.GetInt(0);
			                    if (itb0 - 1 != ib0)
			                        throw Exc(in.Str() + "\n"  + Format(t_("Body # does not match in '%s' %d<>%d"), sarea, itb0, ib0+1));
			                    int itb1 = f.GetInt(1);
			                    if (itb1 - 1 != ib1)
			                        throw Exc(in.Str() + "\n"  + Format(t_("Body # does not match in '%s' %d<>%d"), sarea, itb1, ib1+1));
			                    int itfr = f.GetInt(2);
			                    if (itfr - 1 != ifr)
			                        throw Exc(in.Str() + "\n"  + Format(t_("Frequency # does not match in '%s' %d<>%d"), sarea, itfr, ifr+1));
			                } else
			                    did = 0;
			                for (int jdf = 0; jdf < 6; ++jdf) {
			                    if (am)
			                    	hd().A[6*ib0 + idf][6*ib1 + jdf][ifr] = f.GetDouble(jdf + did);
			                    else
			                        hd().B[6*ib0 + idf][6*ib1 + jdf][ifr] = f.GetDouble(jdf + did);
			                }
			            }
					}
				}
			}
		} else if (line.StartsWith("FORCERAO")) {
	        for (int ib = 0; ib < hd().Nb; ++ib) {
	            for (int ih = 0; ih < hd().Nh; ++ih) {
	                for (int ifr = 0; ifr < hd().Nf; ++ifr) {
	                    f.Load(in.GetLine());
	                  	for (int idf = 0; idf < 6; ++idf) {
		                    int itb = f.GetInt(0);
		                    if (itb - 1 != ib)
		                        throw Exc(in.Str() + "\n"  + Format(t_("Body # does not match in 'FORCERAO' %d<>%d"), itb, ib+1));
		                    int ith = f.GetInt(1);
		                    if (ith - 1 != ih)
		                        throw Exc(in.Str() + "\n"  + Format(t_("Heading # does not match in 'FORCERAO' %d<>%d"), ith, ih+1));
		                    int itfr = f.GetInt(2);
							if (itfr - 1 != ifr)
								throw Exc(in.Str() + "\n"  + Format(t_("Frequency # does not match in 'FORCERAO' %d<>%d"), itfr, ifr+1));
			                hd().ex.ma[ih](ifr, idf + 6*ib) = f.GetDouble(idf + 3);
	                  	}
	                  	f.Load(in.GetLine());
	                  	for (int idf = 0; idf < 6; ++idf) 
	                       	hd().ex.ph[ih](ifr, idf + 6*ib) = -f.GetDouble(idf)*M_PI/180;
	                    for (int idf = 0; idf < 6; ++idf) {   	
		                    hd().ex.re[ih](ifr, idf + 6*ib) = hd().ex.ma[ih](ifr, idf + 6*ib)*cos(hd().ex.ph[ih](ifr, idf*ib));
		       				hd().ex.im[ih](ifr, idf + 6*ib) = hd().ex.ma[ih](ifr, idf + 6*ib)*sin(hd().ex.ph[ih](ifr, idf*ib));
	                    }
	                }
	            }
	        }
		}
	}
	return true;
}

bool Aqwa::Load_LIS() {
	String fileName = ForceExt(hd().file, ".LIS");
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;

	hd().Nb = hd().Nh = hd().Nf = Null;
	hd().head.Clear();
	hd().w.Clear();
	hd().T.Clear();
	hd().rho = hd().g = hd().h = Null;
	hd().dataFromW = true;
	
	String line; 
	FieldSplit f(in);
	f.IsSeparator = IsTabSpace;
	int pos;
	
	FileInLine::Pos fpos = in.GetPos();
	
	hd().Nb = 0;
	while(!in.IsEof()) {
		line = Trim(in.GetLine());
		
		if ((pos = line.FindAfter("F O R   S T R U C T U R E")) >= 0) {
			int ib = ScanInt(line.Mid(pos)); 
			hd().Nb = max(ib, hd().Nb);
		}
	}
	if (hd().Nb == 0)
		throw Exc(t_("Number of bodies not found"));
					
	in.SeekPos(fpos);			
				
	while(!in.IsEof()) {
		line = TrimBoth(in.GetLine());
		f.Load(line);
		
		if (line.StartsWith("WATER  DEPTH")) 
			hd().h = f.GetDouble(f.LAST);
		else if (line.StartsWith("DENSITY  OF  WATER")) 
			hd().rho = f.GetDouble(f.LAST);
		else if (line.StartsWith("ACCELERATION  DUE  TO GRAVITY")) {
			hd().g = f.GetDouble(f.LAST);	
			break;
		}
	}
	
	hd().names.SetCount(hd().Nb);
	hd().Vo.SetCount(hd().Nb, Null);
	hd().cg.setConstant(3, hd().Nb, Null);
	hd().c0.setConstant(3, hd().Nb, Null);
	hd().cb.setConstant(3, hd().Nb, Null);
	hd().C.SetCount(hd().Nb);
	for (int ib = 0; ib < hd().Nb; ++ib) 
		hd().C[ib].setConstant(6, 6, 0);
	hd().M.SetCount(hd().Nb);
	for (int ib = 0; ib < hd().Nb; ++ib) 
		hd().M[ib].setConstant(6, 6, 0);
	
	int ib;	
	
	while(!in.IsEof()) {
		line = TrimBoth(in.GetLine());
		
		if ((pos = line.FindAfter("S T R U C T U R E")) >= 0) {
			ib = ScanInt(line.Mid(pos)) - 1; 
			if (ib >= hd().Nb)
				throw Exc(in.Str() + "\n"  + Format(t_("Wrong body %d"), ib));
		} else if (line.StartsWith("PMAS")) {
			f.Load(line);
			double mass = f.GetDouble(2);
			hd().M[ib](0, 0) = hd().M[ib](1, 1) = hd().M[ib](2, 2) = mass;
		} else if (line.StartsWith("CENTRE OF GRAVITY")) {
			f.Load(line);
			hd().cg(0, ib) = f.GetDouble(3);
			hd().cg(1, ib) = f.GetDouble(4);
			hd().cg(2, ib) = f.GetDouble(5);
			hd().c0 = clone(hd().cg);
		} else if (line.StartsWith("INERTIA MATRIX")) {
			Eigen::MatrixXd &inertia = hd().M[ib];
			f.Load(line);
			inertia(3, 3) = f.GetDouble(2);
			inertia(3, 4) = f.GetDouble(3);
			inertia(3, 5) = f.GetDouble(4);
			in.GetLine();
			f.Load(TrimBoth(in.GetLine()));
			inertia(4, 3) = f.GetDouble(0);
			inertia(4, 4) = f.GetDouble(1);
			inertia(4, 5) = f.GetDouble(2);
			in.GetLine();
			f.Load(TrimBoth(in.GetLine()));
			inertia(5, 3) = f.GetDouble(0);
			inertia(5, 4) = f.GetDouble(1);
			inertia(5, 5) = f.GetDouble(2);
			double mass = inertia(0, 0);
			double cx = mass*hd().cg(0, ib);
			double cy = mass*hd().cg(1, ib);
			double cz = mass*hd().cg(2, ib);
			inertia(1, 5) = inertia(5, 1) =  cx;
			inertia(2, 4) = inertia(4, 2) = -cx;
			inertia(0, 5) = inertia(5, 0) = -cy;
			inertia(2, 3) = inertia(3, 2) =  cy;
			inertia(0, 4) = inertia(4, 0) =  cz;
			inertia(1, 3) = inertia(3, 1) = -cz;
		} else if (IsNull(hd().Nf) && line.Find("W A V E   F R E Q U E N C I E S / P E R I O D S   A N D   D I R E C T I O N S") >= 0) {
			in.GetLine(5);
			hd().Nf = 0;
			bool newparagraph = false;
			while (!in.IsEof()) {
				line = in.GetLine();
				if (line[0] == '1') {
					in.GetLine(6);
					line = in.GetLine();
					newparagraph = true;
				}
				if (TrimBoth(line).StartsWith("--------"))
					break;
				f.Load(line);
				double w;
				if (hd().w.size() == 0 || newparagraph) {
					w = f.GetDouble(2);
					newparagraph = false;
				} else
					w = f.GetDouble(1);
				hd().w << w;
				hd().T << 2*M_PI/w;
			}
			hd().Nf = hd().w.size();
		} else if (IsNull(hd().Nh) && line.Find("DIRECTIONS") >= 0) {
			int idini = (line.Find("STRUCTURE") >= 0) ? 1 : 0;
			while (!in.IsEof()) {
				line = in.GetLine();
				if (TrimBoth(line).StartsWith("--------"))
					break;
			}
			while (!in.IsEof()) {
				line = in.GetLine();
				if (TrimBoth(line).StartsWith("--------"))
					break;
				f.Load(line);
				for (int i = idini; i < f.size(); ++i)
					hd().head << f.GetDouble(i);
				idini = 0;
			}
			hd().Nh = hd().head.size();
			break;
		}
	}
	if (IsNull(hd().Nf))
		throw Exc(t_("Number of frequencies not found"));
	if (IsNull(hd().Nh))
		throw Exc(t_("Number of headings not found"));
					
	hd().A.SetCount(6*hd().Nb);
	hd().B.SetCount(6*hd().Nb);
	for (int i = 0; i < 6*hd().Nb; ++i) {
		hd().A[i].SetCount(6*hd().Nb);
		hd().B[i].SetCount(6*hd().Nb);
		for (int j = 0; j < 6*hd().Nb; ++j) {
			hd().A[i][j].setConstant(hd().Nf, Null);	
			hd().B[i][j].setConstant(hd().Nf, Null);	
		}
	}
	
	hd().Initialize_Forces(hd().ex);
	hd().Initialize_Forces(hd().fk);
	hd().Initialize_Forces(hd().sc);
	hd().Initialize_RAO();
	
	while(!in.IsEof()) {
		line = TrimBoth(in.GetLine());
		
		if ((pos = line.FindAfter("S T R U C T U R E")) >= 0) {
			ib = ScanInt(line.Mid(pos));
			if (!IsNull(ib)) {
				ib -= 1; 
				if (ib >= hd().Nb)
					throw Exc(in.Str() + "\n"  + Format(t_("Wrong body %d"), ib));
			}
		} else if (line.Find("STIFFNESS MATRIX AT THE CENTRE OF GRAVITY") >= 0 ||
				   line.Find("TOTAL HYDROSTATIC STIFFNESS") >= 0) {
			in.GetLine(3);
			line = in.GetLine();
			if (Trim(line).StartsWith("Z"))
				line = in.GetLine();
			pos = line.FindAfter("=");
			if (pos < 0)
				throw Exc(in.Str() + "\n"  + t_("Format error, '=' not found"));
			f.Load(line.Mid(pos));
			hd().C[ib](2, 2) = f.GetDouble(0);
			hd().C[ib](2, 3) = f.GetDouble(1);
			hd().C[ib](2, 4) = f.GetDouble(2);
			if (f.size() > 3)
				hd().C[ib](2, 5) = f.GetDouble(3);	
			line = in.GetLine();
			pos = line.FindAfter("=");
			if (pos < 0)
				throw Exc(in.Str() + "\n"  + t_("Format error, '=' not found"));
			f.Load(line.Mid(pos));
			hd().C[ib](3, 2) = f.GetDouble(0);
			hd().C[ib](3, 3) = f.GetDouble(1);
			hd().C[ib](3, 4) = f.GetDouble(2);
			if (f.size() > 3)
				hd().C[ib](3, 5) = f.GetDouble(3);	
			line = in.GetLine();
			pos = line.FindAfter("=");
			if (pos < 0)
				throw Exc(in.Str() + "\n"  + t_("Format error, '=' not found"));
			f.Load(line.Mid(pos));
			hd().C[ib](4, 2) = f.GetDouble(0);
			hd().C[ib](4, 3) = f.GetDouble(1);
			hd().C[ib](4, 4) = f.GetDouble(2);
			if (f.size() > 3)
	       		hd().C[ib](4, 5) = f.GetDouble(3);	
		} else if (line.Find("STIFFNESS MATRIX") >= 0 && IsNull(hd().C[ib](0,0))) {	// 2nd option to get stiffness matrix
			hd().C[ib].setConstant(6, 6, 0);
			in.GetLine(6);
			for (int r = 0; r < 6; ++r) {
				f.Load(in.GetLine());
				for (int c = 0; c < 6; ++c) 
					hd().C[ib](r, c) = f.GetDouble(c + 1);
				in.GetLine();
			}
		} else if (line.StartsWith("MESH BASED DISPLACEMENT")) {
			pos = line.FindAfter("=");
			if (pos < 0)
				throw Exc(in.Str() + "\n"  + t_("= not found"));
			hd().Vo[ib] = ScanDouble(line.Mid(pos));
		} else if (line.StartsWith("POSITION OF THE CENTRE OF BUOYANCY")) {
			f.Load(line);
			hd().cb(0, ib) = f.GetDouble(8);	
			hd().cb(1, ib) = f.Load(in.GetLine()).GetDouble(2);
			hd().cb(2, ib) = f.Load(in.GetLine()).GetDouble(2);
		} else if (line.StartsWith("FROUDE KRYLOV + DIFFRACTION FORCES-VARIATION WITH WAVE PERIOD/FREQUENCY") ||
				   line.StartsWith(				 "FROUDE KRYLOV FORCES-VARIATION WITH WAVE PERIOD/FREQUENCY") ||
				   line.StartsWith(  			   "DIFFRACTION FORCES-VARIATION WITH WAVE PERIOD/FREQUENCY") ||
				   line.StartsWith(             			  "R.A.O.S-VARIATION WITH WAVE PERIOD/FREQUENCY")) {
			Hydro::Forces *pfrc;
			if (line.StartsWith("FROUDE KRYLOV + DIFFRACTION FORCES-VARIATION WITH WAVE PERIOD/FREQUENCY"))
				pfrc = &hd().ex;
			else if (line.StartsWith("FROUDE KRYLOV FORCES-VARIATION WITH WAVE PERIOD/FREQUENCY")) 
				pfrc = &hd().fk;
			else if (line.StartsWith("DIFFRACTION FORCES-VARIATION WITH WAVE PERIOD/FREQUENCY")) 
				pfrc = &hd().sc;
			else //if (line.StartsWith("R.A.O.S-VARIATION WITH WAVE PERIOD/FREQUENCY")) 
				pfrc = &hd().rao;
			Hydro::Forces &frc = *pfrc; 
			
			in.GetLine(5);
			line = in.GetLine();
			while(!in.IsEof()) {
				if (line[0] == '1')
					break;
				
				static const Vector<int> separatorsh = {8,16,26,36,44,54,62,72,80,90,98,108,116,126};
				f.Load(in.GetLine(), separatorsh);

				double heading = f.GetDouble(2);
				int idh = FindClosest(hd().head, heading);
				if (idh < 0)
					throw Exc(in.Str() + "\n"  + Format(t_("Heading %f is unknown"), heading));
				int dd = 1;
				for (int ifr = 0; ifr < hd().Nf; ++ifr) {
					double freq = f.GetDouble(1);
					int ifrr = FindClosest(hd().w, freq);
					if (ifrr < 0)
						throw Exc(in.Str() + "\n"  + Format(t_("Frequency %f is unknown"), freq));
					for (int idf = 0; idf < 6; ++idf) {
						frc.ma[idh](ifr, idf + 6*ib) = f.GetDouble(2 + dd + idf*2);
						frc.ph[idh](ifr, idf + 6*ib) = -f.GetDouble(2 + dd + idf*2 + 1)*M_PI/180; // Negative to follow Wamit
						frc.re[idh](ifr, idf + 6*ib) = frc.ma[idh](ifr, idf + 6*ib)*cos(frc.ph[idh](ifr, idf + 6*ib));
			       		frc.im[idh](ifr, idf + 6*ib) = frc.ma[idh](ifr, idf + 6*ib)*sin(frc.ph[idh](ifr, idf + 6*ib));
					}
					dd = 0;
					line = in.GetLine();
					static const Vector<int> separators = {8,16,36,44,54,62,72,80,90,98,108,116,126};
					f.Load(line, separators);
				}
			}
		} else if (line.Find("WAVE PERIOD") >= 0 && line.Find("WAVE FREQUENCY") >= 0) {
			int ieq = line.FindAfter("="); ieq = line.FindAfter("=", ieq);
			f.Load(line);
			double freq = ScanDouble(line.Mid(ieq));
			if (IsNull(freq))
				throw Exc(in.Str() + "\n"  + t_("Problem loading frequency"));
			int ifr = FindClosest(hd().w, freq);
			if (ifr < 0)
				throw Exc(in.Str() + "\n"  + Format(t_("Frequency %f is unknown"), freq));
			
			in.GetLine(2);
			if (TrimBoth(in.GetLine()) == "ADDED  MASS") {
				in.GetLine(5);
			
				for (int idf = 0; idf < 6; ++idf) {
					in.GetLine();
					f.Load(in.GetLine());
					if (f.GetText(0) != textDOF[idf])
						throw Exc(in.Str() + "\n"  + Format(t_("Expected %s data, found '%s'"), textDOF[idf], f.GetText())); 
					for (int jdf = 0; jdf < 6; ++jdf) 
						hd().A[6*ib + idf][6*ib + jdf][ifr] = f.GetDouble(1 + jdf);
				}
				in.GetLine(8);
				for (int idf = 0; idf < 6; ++idf) {
					in.GetLine();
					f.Load(in.GetLine());
					if (f.GetText(0) != textDOF[idf])
						throw Exc(in.Str() + "\n"  + Format(t_("Expected %s data, found '%s'"), textDOF[idf], f.GetText())); 
					for (int jdf = 0; jdf < 6; ++jdf) 
						hd().B[6*ib + idf][6*ib + jdf][ifr] = f.GetDouble(1 + jdf);
				}
			}
		} else if (line.Find("H Y D R O D Y N A M I C   P A R A M E T E R S   A T   L O W   &   H I G H") >= 0) {
			int ib = -1;
			double freq = -1;
			while(!in.IsEof() && !(line[0] == '1')) {
				line = in.GetLine();
				int id = line.FindAfter("F O R   S T R U C T U R E");
				if (id >= 0) {
					ib = ScanInt(line.Mid(id)) - 1;
					for (int ii = 0; ii < 2; ++ii) {
						while(!in.IsEof() && !(line[0] == '1')) {
							line = in.GetLine();
							int idr = line.FindAfter("WAVE FREQUENCY");
							if (idr >= 0) {
								freq = ScanDouble(line.Mid(idr + 3));
								if (freq < 0.1 && hd().A0.size() == 0)
									hd().A0.resize(6*hd().Nb, 6*hd().Nb);
								else if (freq > 90 && hd().Ainf.size() == 0)
									hd().Ainf.resize(6*hd().Nb, 6*hd().Nb);
								break;
							}
						}
						while(!in.IsEof() && !(line[0] == '1')) {
							line = in.GetLine();
							if (line.Find("ADDED MASS") >= 0) {
								while(!in.IsEof() && !(line[0] == '1')) {
									f.LoadLine();
									if (f.size() == 7 && f.GetText(0) == "X") {
										for (int idf = 0; idf < 6; ++idf) {
											for (int jdf = 0; jdf < 6; ++jdf) {
												if (freq < 0.1)
													hd().A0(6*ib + idf, 6*ib + jdf) = f.GetDouble(jdf + 1);
												else if (freq > 90)
													hd().Ainf(6*ib + idf, 6*ib + jdf) = f.GetDouble(jdf + 1);
											}
											f.LoadLine(); 
											f.LoadLine();
										}
										break;
									}
								}
								freq = -1;
								break;
							}
						}
					}
				}
			}
		} else if (line.Find("FREQUENCY INDEPENDENT DAMPING") >= 0) {
			if (hd().Dlin.size() == 0)
				hd().Dlin = Eigen::MatrixXd::Zero(6*hd().Nb, 6*hd().Nb);
			in.GetLine(6);
			for (int idf = 0; idf < 6; ++idf) {
				f.LoadLine();
				for (int jdf = 0; jdf < 6; ++jdf) 
					hd().Dlin(6*ib + idf, 6*ib + jdf) = f.GetDouble(jdf + 1);
				f.LoadLine(); 
			}
		}
	}
		
	return true;
}

bool Aqwa::Load_QTF() {
	String fileName = ForceExt(hd().file, ".QTF");
	FileInLine in(fileName);
	if (!in.IsOpen()) {
		fileName = AppendFileNameX(GetFileFolder(fileName), "analysis.qtf"); 
		in.Open(fileName);
		if (!in.IsOpen()) 
			return false;
	}
	
	String line; 
	FieldSplit f(in);
	f.IsSeparator = IsTabSpace;
	
	in.GetLine();	// AQWA version	
	
	hd().qtfsum.Clear();
	hd().qtfdif.Clear();
	hd().qtfCases.Clear();
	
	int nrows = hd().Nh*hd().Nf*hd().Nf;
	hd().qtfdif.Reserve(hd().Nb*nrows);
	hd().qtfsum.Reserve(hd().Nb*nrows);
	
	while (!in.IsEof()) {
		f.Load(in.GetLine());
		int ib = f.GetInt(0);
		if (IsNull(hd().Nb)) {
			hd().Nb = ib;
			if (hd().names.IsEmpty())
				hd().names.SetCount(hd().Nb);
		}
		else if (ib > hd().Nb)
			throw Exc(in.Str() + "\n"  + Format(t_("#%d body found when max are %d"), ib+1, hd().Nb));
		int Nh = f.GetInt(1);
		int Nf = f.GetInt(2);
		
		int col = 3;
		int ih = 0;
		while (!in.IsEof()) {		// Check headings
			while (col < f.size() && ih < Nh) {
				double head = f.GetDouble(col++);
				FindAddRatio(hd().qtfhead, head, 0.001);
				ih++;
			}
			if (ih >= Nh)
				break;
			f.Load(in.GetLine());
			col = 0;
		}
		f.Load(in.GetLine());
		col = 0;
		int ifr = 0;
		while (!in.IsEof()) {		// Check frequencies
			while (col < f.size() && ifr < Nf) {
				double w = f.GetDouble(col++);
				FindAddRatio(hd().qtfw, w, 0.001);
				ifr++;
			}
			if (ifr >= Nf)
				break;
			f.Load(in.GetLine());
			col = 0;
		}
		hd().qtfT.Clear();
		for (int ifr = 0; ifr < hd().qtfw.size(); ++ifr)
			hd().qtfT << 2*M_PI/hd().qtfw[ifr];
		
		nrows = Nh*Nf*Nf;
		for (int i = 0; i < nrows; ++i) {
			f.Load(in.GetLine());
			int ib = f.GetInt(0)-1;
			if (ib >= hd().Nb)
				throw Exc(in.Str() + "\n"  + Format(t_("Body id %d higher than number of bodies"), ib+1, hd().Nb));
			int ih = f.GetInt(1)-1;
			if (ih >= Nh)
				throw Exc(in.Str() + "\n"  + Format(t_("Heading id %d higher than number of headings"), ih+1, Nh));
			int ifr1 = f.GetInt(2)-1;
			if (ifr1 >= Nf)
				throw Exc(in.Str() + "\n"  + Format(t_("Frequency id %d higher than number of frequencies"), ifr1+1, Nf));
			int ifr2 = f.GetInt(3)-1;
			if (ifr2 >= Nf)
				throw Exc(in.Str() + "\n"  + Format(t_("Frequency id %d higher than number of frequencies"), ifr1+2, Nf));
							
			Hydro::QTF &qtfdif = hd().qtfdif.Add();
			qtfdif.Set(ib, ih, ih, ifr1, ifr2);
	        
			for (int idof = 0; idof < 6; ++idof) 
				qtfdif.fre[idof] = f.GetDouble(4 + idof);
				
	        f.Load(in.GetLine());
	        for (int idof = 0; idof < 6; ++idof)
				qtfdif.fim[idof] = f.GetDouble(idof);
		
			for (int idof = 0; idof < 6; ++idof) {
				qtfdif.fma[idof] = sqrt(sqr(qtfdif.fre[idof]) + sqr(qtfdif.fim[idof])); 
				qtfdif.fph[idof] = atan2(qtfdif.fim[idof], qtfdif.fre[idof]);
			}
	
			Hydro::QTF &qtfsum = hd().qtfsum.Add();
			qtfsum.Set(ib, ih, ih, ifr1, ifr2);
					
			f.Load(in.GetLine());
			for (int idof = 0; idof < 6; ++idof) 
				qtfsum.fre[idof] = f.GetDouble(idof);
				
	        f.Load(in.GetLine());
	        for (int idof = 0; idof < 6; ++idof)
				qtfsum.fim[idof] = -f.GetDouble(idof);		// Negative to follow Wamit. Just to the sum term
		
			for (int idof = 0; idof < 6; ++idof) {
				qtfsum.fma[idof] = sqrt(sqr(qtfsum.fre[idof]) + sqr(qtfsum.fim[idof])); 
				qtfsum.fph[idof] = atan2(qtfsum.fim[idof], qtfsum.fre[idof]);
			}
		}
	}
	hd().GetQTFList(hd().qtfsum, hd().qtfCases);
	
	return true;
}

void Aqwa::Save(String ) {
	throw Exc("Option not implemented");
}		

bool AQWACase::Load(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	solver = AQWA;
	
	String line;
	FieldSplit f(in);
	f.IsSeparator = IsTabSpace;
	
	bodies.SetCount(1);
	BEMBody &body = bodies[0];
	
	body.meshFile = fileName;
	body.ndof = 6;
	Resize(body.dof, 6, true);
	
	
	Vector<double> hrtz, head;
	
	while (!in.IsEof()) {
		f.Load(in.GetLine());
		
		if (f.size() == 2) {
			if (f.GetText(0) == "DPTH")
				h = f.GetDouble(1);
			else if (f.GetText(0) == "DENS")
				rho = f.GetDouble(1);
			else if (f.GetText(0) == "ACCG")
				g = f.GetDouble(1);	
		}
		
		if (f.size() == 4) {
			if (f.GetText(0) == "1HRTZ")
				hrtz << f.GetDouble(3);
			else if (f.GetText(0) == "1DIRN") 
				head << f.GetDouble(3);
			else if (f.GetInt_nothrow(0) == 198000) {
				body.cg[0] = f.GetDouble(1);
				body.cg[1] = f.GetDouble(2);
				body.cg[2] = f.GetDouble(3);
				body.c0 = clone(body.cg);
			}
		}
	}
	Nf = hrtz.size();
	if (Nf <= 0)
		throw Exc("Number of frequencies should have to be higher than zero");
	Nh = head.size();
	if (Nh <= 0)
		throw Exc("Number of headings should have to be higher than zero");
	
	minF = fround(hrtz[0]*2*M_PI, 5);				// 5 decimals is enough to filter conversion
	maxF = fround(hrtz[hrtz.size()-1]*2*M_PI, 5);
	
	minH = head[0];
	maxH = head[head.size()-1];
	
	return true;
}