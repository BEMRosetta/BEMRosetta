#include "BEMRosetta.h"

bool Aqwa::Load(String file, double) {
	hd().file = file;
	hd().name = GetFileTitle(GetFileFolder(file));
	hd().len = 1;
	hd().code = Hydro::AQWA;
	
	try {
		hd().Print("\n\n" + Format(t_("Loading '%s'"), file));

		hd().Print("\n- " + x_(t_("LIS file")));
		if (!Load_LIS()) {
			hd().PrintWarning(x_(": **") + t_("Not found") + "**");
			hd().Print("\n- " + x_(t_("AH1 file")));
			if (!Load_AH1()) {
				hd().PrintWarning(x_(": **") + t_("Not found") + "**");
				return false;
			}
		} else {
			hd().Print("\n- " + x_(t_("AH1 file")));
			Aqwa ah1;
			ah1.hd().file = file;
			if (!ah1.Load_AH1()) 
				hd().PrintWarning(x_(": **") + t_("Not found") + "**");
			else {
				hd().Print("\n" + Format(t_("Comparing LIS and AH1 files..."), file));
				hd().Compare_rho(ah1.hd());
				hd().Compare_g(ah1.hd());
				hd().Compare_h(ah1.hd());
				hd().Compare_w(ah1.hd());
				hd().Compare_head(ah1.hd());
				hd().Compare_Nb(ah1.hd());
				hd().Compare_A(ah1.hd());
				hd().Compare_B(ah1.hd());
				hd().Compare_C(ah1.hd());
				hd().Compare_cg(ah1.hd());
				
				hd().ex = pick(ah1.hd().ex);
				hd().Print(t_("Comparison is OK"));
			}
		}
		hd().dof.Clear();	hd().dof.SetCount(hd().Nb, 0);
		for (int i = 0; i < hd().Nb; ++i)
			hd().dof[i] = 6;
	
		hd().AfterLoad();
	} catch (Exc e) {
		hd().PrintError(Format("\n%s: %s", t_("Error"), e));
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
	hd().T.Clear();	
	hd().rho = hd().g = hd().h = Null;
	
	String line;
	FieldSplit f(in);
	while(!in.IsEof()) {
		line = in.GetLine();
		
		if (TrimBoth(line).StartsWith("*"))
			break;
	}
	f.Load(in.GetLine());
	hd().Nb = f.GetInt(0);
	hd().Nh = f.GetInt(1);
	hd().Nf = f.GetInt(2);
	
	if (hd().Nh != f.GetCount() - 3)
		throw Exc(Format(t_("[%d] Number of headings do not match %d<>%d"), in.GetLineNumber(), hd().Nh, f.GetCount() - 3));
	for (int i = 3; i < f.GetCount(); ++i)
		hd().head << f.GetDouble(i);
	
	hd().cg.setConstant(3, hd().Nb, Null);
	hd().C.SetCount(hd().Nb);
	for (int ib = 0; ib < hd().Nb; ++ib) 
		hd().C[ib].setConstant(6, 6, Null); 
	hd().A.SetCount(hd().Nf);
	hd().B.SetCount(hd().Nf);
	for (int ifr = 0; ifr < hd().Nf; ++ifr) {
		hd().A[ifr].setConstant(hd().Nb*6, hd().Nb*6, Null);
		hd().B[ifr].setConstant(hd().Nb*6, hd().Nb*6, Null);
	}
	hd().Initialize_Forces(hd().ex);
	
	while(!in.IsEof() && hd().w.GetCount() < hd().Nf) {
		line = TrimBoth(in.GetLine());
		f.Load(line);
		
		for (int i = 0; i < f.GetCount(); ++i) {
			double w = f.GetDouble(i);
			hd().w << w;
			hd().T << 2*M_PI/w;
		}
	}	
	if (hd().Nf != hd().w.GetCount())
		throw Exc(Format(t_("[%d] Number of frequencies do not match %d<>%d"), in.GetLineNumber(), hd().Nf, hd().w.GetCount()));

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
					throw Exc(Format(t_("[%d] Unknown body found in COG %d (%d)"), in.GetLineNumber(), ib+1, hd().Nb));
				
				hd().cg(0, ib) = f.GetDouble(1);
				hd().cg(1, ib) = f.GetDouble(2);
				hd().cg(2, ib) = f.GetDouble(3);
			}
		} else if (line.StartsWith("HYDSTIFFNESS")) {
			for (int ib = 0; ib < hd().Nb; ++ib) {
	            for (int idof = 0; idof < 6; ++idof) {
	                f.Load(in.GetLine());
	                int did = 0;
	                if (idof == 0) {
	                    did = 1;
	                    int itb = f.GetInt(0);
	                    if (itb - 1 != ib)
	                        throw Exc(Format(t_("[%d] Body # does not match in 'HYDSTIFFNESS' %d<>%d"), in.GetLineNumber(), itb, ib+1));
	                }
	                for (int jdof = 0; jdof < 6; ++jdof) 
	                    hd().C[ib](idof, jdof) = f.GetDouble(jdof + did);
	            }
			}
		} else if (line.StartsWith("ADDEDMASS") || line.StartsWith("DAMPING")) {
			bool am = line.StartsWith("ADDEDMASS");
			String sarea = am ? "ADDEDMASS" : "DAMPING";
			for (int ib0 = 0; ib0 < hd().Nb; ++ib0) {
				for (int ib1 = 0; ib1 < hd().Nb; ++ib1) {
					for (int ifr = 0; ifr < hd().Nf; ++ifr) {
			            for (int idof = 0; idof < 6; ++idof) {
			                f.Load(in.GetLine());
			                int did = 0;
			                if (idof == 0) {
			                    did = 3;
			                    int itb0 = f.GetInt(0);
			                    if (itb0 - 1 != ib0)
			                        throw Exc(Format(t_("[%d] Body # does not match in '%s' %d<>%d"), in.GetLineNumber(), sarea, itb0, ib0+1));
			                    int itb1 = f.GetInt(1);
			                    if (itb1 - 1 != ib1)
			                        throw Exc(Format(t_("[%d] Body # does not match in '%s' %d<>%d"), in.GetLineNumber(), sarea, itb1, ib1+1));
			                    int itfr = f.GetInt(2);
			                    if (itfr - 1 != ifr)
			                        throw Exc(Format(t_("[%d] Frequency # does not match in '%s' %d<>%d"), in.GetLineNumber(), sarea, itfr, ifr+1));
			                } else
			                    did = 0;
			                for (int jdof = 0; jdof < 6; ++jdof) {
			                    if (am)
			                    	hd().A[ifr](6*ib0 + idof, 6*ib1 + jdof) = f.GetDouble(jdof + did);
			                    else
			                        hd().B[ifr](6*ib0 + idof, 6*ib1 + jdof) = f.GetDouble(jdof + did);
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
	                  	for (int idof = 0; idof < 6; ++idof) {
		                    int itb = f.GetInt(0);
		                    if (itb - 1 != ib)
		                        throw Exc(Format(t_("[%d] Body # does not match in 'FORCERAO' %d<>%d"), in.GetLineNumber(), itb, ib+1));
		                    int ith = f.GetInt(1);
		                    if (ith - 1 != ih)
		                        throw Exc(Format(t_("[%d] Heading # does not match in 'FORCERAO' %d<>%d"), in.GetLineNumber(), ith, ih+1));
		                    int itfr = f.GetInt(2);
							if (itfr - 1 != ifr)
								throw Exc(Format(t_("[%d] Frequency # does not match in 'FORCERAO' %d<>%d"), in.GetLineNumber(), itfr, ifr+1));
			                hd().ex.ma[ih](ifr, idof + 6*ib) = f.GetDouble(idof + 3);
	                  	}
	                  	f.Load(in.GetLine());
	                  	for (int idof = 0; idof < 6; ++idof) 
	                       	hd().ex.ph[ih](ifr, idof + 6*ib) = -f.GetDouble(idof)*M_PI/180;
	                    for (int idof = 0; idof < 6; ++idof) {   	
		                    hd().ex.re[ih](ifr, idof + 6*ib) = hd().ex.ma[ih](ifr, idof + 6*ib)*cos(hd().ex.ph[ih](ifr, idof*ib));
		       				hd().ex.im[ih](ifr, idof + 6*ib) = hd().ex.ma[ih](ifr, idof + 6*ib)*sin(hd().ex.ph[ih](ifr, idof*ib));
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
	
	String line; 
	FieldSplit f(in);
	int pos;
	
	FileInLine::Pos fpos = in.GetPos();
	
	hd().Nb = 0;
	while(!in.IsEof()) {
		line = TrimBoth(in.GetLine());
		
		if ((pos = line.FindAfter("F O R   S T R U C T U R E")) >= 0) {
			int ib = ScanInt(line.Mid(pos)); 
			hd().Nb = max(ib, hd().Nb);
		}
	}
	if (hd().Nb == 0)
		throw Exc(t_("Number of bodies not found"));
					
	in.Seek(fpos);			
				
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
	
	hd().Vo.SetCount(hd().Nb, Null);
	hd().cg.setConstant(3, hd().Nb, Null);
	hd().cb.setConstant(3, hd().Nb, Null);
	hd().C.SetCount(hd().Nb);
	for (int ib = 0; ib < hd().Nb; ++ib) 
		hd().C[ib].setConstant(6, 6, Null);
	
	int idb;	
	
	while(!in.IsEof()) {
		line = TrimBoth(in.GetLine());
		
		if ((pos = line.FindAfter("S T R U C T U R E")) >= 0) {
			idb = ScanInt(line.Mid(pos)) - 1; 
			if (idb >= hd().Nb)
				throw Exc(Format(t_("[%d] Wrong body %d"), in.GetLineNumber(), idb));
		} else if (line.StartsWith("CENTRE OF GRAVITY")) {
			f.Load(line);
			hd().cg(0, idb) = f.GetDouble(3);
			hd().cg(1, idb) = f.GetDouble(4);
			hd().cg(2, idb) = f.GetDouble(5);
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
				if (hd().w.GetCount() == 0 || newparagraph) {
					w = f.GetDouble(2);
					newparagraph = false;
				} else
					w = f.GetDouble(1);
				hd().w << w;
				hd().T << 2*M_PI/w;
			}
			hd().Nf = hd().w.GetCount();
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
				for (int i = idini; i < f.GetCount(); ++i)
					hd().head << f.GetDouble(i);
				idini = 0;
			}
			hd().Nh = hd().head.GetCount();
			break;
		}
	}
	if (IsNull(hd().Nf))
		throw Exc(t_("Number of frequencies not found"));
	if (IsNull(hd().Nh))
		throw Exc(t_("Number of headings not found"));
					
	hd().A.SetCount(hd().Nf);
	hd().B.SetCount(hd().Nf);
	for (int ifr = 0; ifr < hd().Nf; ++ifr) {
		hd().A[ifr].setConstant(hd().Nb*6, hd().Nb*6, Null);
		hd().B[ifr].setConstant(hd().Nb*6, hd().Nb*6, Null);
	}
	
	hd().Initialize_Forces(hd().fk);
	hd().Initialize_Forces(hd().sc);
	hd().Initialize_RAO();
	
	while(!in.IsEof()) {
		line = TrimBoth(in.GetLine());
		
		if ((pos = line.FindAfter("S T R U C T U R E")) >= 0) {
			idb = ScanInt(line.Mid(pos)) - 1; 
			if (idb >= hd().Nb)
				throw Exc(Format(t_("[%d] Wrong body %d"), in.GetLineNumber(), idb));
		} else if (line.Find("STIFFNESS MATRIX AT THE CENTRE OF GRAVITY") >= 0) {
			in.GetLine(3);
			line = in.GetLine();
			pos = line.FindAfter("=");
			if (pos < 0)
				throw Exc(Format(t_("[%d] Format error, '=' not found"), in.GetLineNumber()));
			f.Load(line.Mid(pos));
	       	hd().C[idb](2, 2) = f.GetDouble(0);
	       	hd().C[idb](2, 3) = f.GetDouble(1);
	       	hd().C[idb](2, 4) = f.GetDouble(2);
			line = in.GetLine();
			pos = line.FindAfter("=");
			if (pos < 0)
				throw Exc(Format(t_("[%d] Format error, '=' not found"), in.GetLineNumber()));
			f.Load(line.Mid(pos));
	       	hd().C[idb](3, 2) = f.GetDouble(0);
			hd().C[idb](3, 3) = f.GetDouble(1);
			hd().C[idb](3, 4) = f.GetDouble(2);
			line = in.GetLine();
			pos = line.FindAfter("=");
			if (pos < 0)
				throw Exc(Format(t_("[%d] Format error, '=' not found"), in.GetLineNumber()));
			f.Load(line.Mid(pos));
	       	hd().C[idb](4, 2) = f.GetDouble(0);
			hd().C[idb](4, 3) = f.GetDouble(1);
			hd().C[idb](4, 4) = f.GetDouble(2);
		} else if (line.StartsWith("MESH BASED DISPLACEMENT")) {
			pos = line.FindAfter("=");
			if (pos < 0)
				throw Exc(Format(t_("[%d] = not found"), in.GetLineNumber()));
			hd().Vo[idb] = ScanDouble(line.Mid(pos));
		} else if (line.StartsWith("POSITION OF THE CENTRE OF BUOYANCY")) {
			f.Load(line);
			hd().cb(0, idb) = f.GetDouble(8);	
			hd().cb(1, idb) = f.Load(in.GetLine()).GetDouble(2);
			hd().cb(2, idb) = f.Load(in.GetLine()).GetDouble(2);
		} else if (line.StartsWith("FROUDE KRYLOV + DIFFRACTION FORCES-VARIATION WITH WAVE PERIOD/FREQUENCY") ||
				   line.StartsWith("DIFFRACTION FORCES-VARIATION WITH WAVE PERIOD/FREQUENCY") ||
				   line.StartsWith("R.A.O.S-VARIATION WITH WAVE PERIOD/FREQUENCY")) {
			Hydro::Forces *pfrc;
			if (line.StartsWith("FROUDE KRYLOV + DIFFRACTION FORCES-VARIATION WITH WAVE PERIOD/FREQUENCY")) 
				pfrc = &hd().fk;
			else if (line.StartsWith("DIFFRACTION FORCES-VARIATION WITH WAVE PERIOD/FREQUENCY")) 
				pfrc = &hd().sc;
			else if (line.StartsWith("R.A.O.S-VARIATION WITH WAVE PERIOD/FREQUENCY")) 
				pfrc = &hd().rao;
			Hydro::Forces &frc = *pfrc; 
			
			in.GetLine(5);
			line = in.GetLine();
			while(!in.IsEof()) {
				if (line[0] == '1')
					break;
				f.Load(in.GetLine());
				double heading = f.GetDouble(2);
				int idh = FindIndex(hd().head, heading);
				if (idh < 0)
					throw Exc(Format(t_("[%d] Heading %f not found"), in.GetLineNumber(), heading));
				int dd = 1;
				for (int i = 0; i < hd().Nf; ++i) {
					double freq = f.GetDouble(1);
					int ifr = FindIndexDelta(hd().w, freq, 0.001);
					if (ifr < 0)
						throw Exc(Format(t_("[%d] Frequency %f not found"), in.GetLineNumber(), freq));
					for (int idof = 0; idof < 6; ++idof) {
						frc.ma[idh](ifr, idof + 6*idb) = f.GetDouble(2 + dd + idof*2);
						frc.ph[idh](ifr, idof + 6*idb) = f.GetDouble(2 + dd + idof*2 + 1);
						frc.re[idh](ifr, idof + 6*idb) = frc.ma[idh](ifr, idof + 6*idb)*cos(frc.ph[idh](ifr, idof + 6*idb));
			       		frc.im[idh](ifr, idof + 6*idb) = frc.ma[idh](ifr, idof + 6*idb)*sin(frc.ph[idh](ifr, idof + 6*idb));
					}
					dd = 0;
					line = in.GetLine();
					f.Load(line);
				}
			}
		} else if (line.Find("WAVE PERIOD") >= 0 && line.Find("WAVE FREQUENCY") >= 0) {
			f.Load(line);
			double freq = f.GetDouble(7);
			int ifr = FindIndexDelta(hd().w, freq, 0.001);
			if (ifr < 0)
				throw Exc(t_(Format(t_("[%d] Frequency %f not found"), in.GetLineNumber(), freq)));

			in.GetLine(2);
			if (TrimBoth(in.GetLine()) == "ADDED  MASS") {
				in.GetLine(5);
			
				for (int idof = 0; idof < 6; ++idof) {
					in.GetLine();
					f.Load(in.GetLine());
					for (int jdof = 0; jdof < 6; ++jdof) 
						hd().A[ifr](6*idb + idof, 6*idb + jdof) = f.GetDouble(1 + jdof);
				}
				in.GetLine(8);
				for (int idof = 0; idof < 6; ++idof) {
					in.GetLine();
					f.Load(in.GetLine());
					for (int jdof = 0; jdof < 6; ++jdof) 
						hd().B[ifr](6*idb + idof, 6*idb + jdof) = f.GetDouble(1 + jdof);
				}
			}
		}
	}
		
	return true;
}
	
void Aqwa::Save(String file) {
	throw Exc("Option not implemented");
}		
