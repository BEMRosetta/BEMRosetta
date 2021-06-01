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
		BEMData::Print("\n\n" + Format(t_("Loading '%s'"), file));

		BEMData::Print("\n- " + S(t_("LIS file")));
		if (!Load_LIS()) {
			BEMData::PrintWarning(S(": **") + t_("Not found") + "**");
			BEMData::Print("\n- " + S(t_("AH1 file")));
			if (!Load_AH1()) {
				BEMData::PrintWarning(S(": **") + t_("Not found") + "**");
				hd().Nh = hd().Nf = 0;
			//	throw Exc(t_("No .AH1 or .LIS file found"));
			}
		}
		//if (IsNull(hd().Nb))
		//	return false;
		
		BEMData::Print("\n- " + S(t_("QTF file")));
		if (!Load_QTF()) 
			BEMData::PrintWarning(S(": **") + t_("Not found") + "**");
		
		if (IsNull(hd().Nb))
			return false;
	
		hd().dof.Clear();	hd().dof.SetCount(hd().Nb, 0);
		for (int i = 0; i < hd().Nb; ++i)
			hd().dof[i] = 6;
	} catch (Exc e) {
		BEMData::PrintError(Format("\n%s: %s", t_("Error"), e));
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
				throw Exc(in.Str() + "\n"  + Format(t_("Wrong body %d"), idb));
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
			idb = ScanInt(line.Mid(pos));
			if (!IsNull(idb)) {
				idb -= 1; 
				if (idb >= hd().Nb)
					throw Exc(in.Str() + "\n"  + Format(t_("Wrong body %d"), idb));
			}
		} else if (line.Find("STIFFNESS MATRIX AT THE CENTRE OF GRAVITY") >= 0 ||
				   line.Find("TOTAL HYDROSTATIC STIFFNESS") >= 0) {
			hd().C[idb].setConstant(6, 6, 0);
			in.GetLine(3);
			line = in.GetLine();
			if (Trim(line).StartsWith("Z"))
				line = in.GetLine();
			pos = line.FindAfter("=");
			if (pos < 0)
				throw Exc(in.Str() + "\n"  + t_("Format error, '=' not found"));
			f.Load(line.Mid(pos));
			hd().C[idb](2, 2) = f.GetDouble(0);
			hd().C[idb](2, 3) = f.GetDouble(1);
			hd().C[idb](2, 4) = f.GetDouble(2);
			if (f.size() > 3)
				hd().C[idb](2, 5) = f.GetDouble(3);	
			line = in.GetLine();
			pos = line.FindAfter("=");
			if (pos < 0)
				throw Exc(in.Str() + "\n"  + t_("Format error, '=' not found"));
			f.Load(line.Mid(pos));
			hd().C[idb](3, 2) = f.GetDouble(0);
			hd().C[idb](3, 3) = f.GetDouble(1);
			hd().C[idb](3, 4) = f.GetDouble(2);
			if (f.size() > 3)
				hd().C[idb](3, 5) = f.GetDouble(3);	
			line = in.GetLine();
			pos = line.FindAfter("=");
			if (pos < 0)
				throw Exc(in.Str() + "\n"  + t_("Format error, '=' not found"));
			f.Load(line.Mid(pos));
			hd().C[idb](4, 2) = f.GetDouble(0);
			hd().C[idb](4, 3) = f.GetDouble(1);
			hd().C[idb](4, 4) = f.GetDouble(2);
			if (f.size() > 3)
	       		hd().C[idb](4, 5) = f.GetDouble(3);	
		} else if (line.Find("STIFFNESS MATRIX") >= 0 && IsNull(hd().C[idb](0,0))) {	// 2nd option to get stiffness matrix
			hd().C[idb].setConstant(6, 6, 0);
			in.GetLine(6);
			for (int r = 0; r < 6; ++r) {
				f.Load(in.GetLine());
				for (int c = 0; c < 6; ++c) 
					hd().C[idb](r, c) = f.GetDouble(c + 1);
				in.GetLine();
			}
		} else if (line.StartsWith("MESH BASED DISPLACEMENT")) {
			pos = line.FindAfter("=");
			if (pos < 0)
				throw Exc(in.Str() + "\n"  + t_("= not found"));
			hd().Vo[idb] = ScanDouble(line.Mid(pos));
		} else if (line.StartsWith("POSITION OF THE CENTRE OF BUOYANCY")) {
			f.Load(line);
			hd().cb(0, idb) = f.GetDouble(8);	
			hd().cb(1, idb) = f.Load(in.GetLine()).GetDouble(2);
			hd().cb(2, idb) = f.Load(in.GetLine()).GetDouble(2);
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
				f.Load(in.GetLine());
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
						frc.ma[idh](ifr, idf + 6*idb) = f.GetDouble(2 + dd + idf*2);
						frc.ph[idh](ifr, idf + 6*idb) = -f.GetDouble(2 + dd + idf*2 + 1)*M_PI/180; // Negative to follow Wamit
						frc.re[idh](ifr, idf + 6*idb) = frc.ma[idh](ifr, idf + 6*idb)*cos(frc.ph[idh](ifr, idf + 6*idb));
			       		frc.im[idh](ifr, idf + 6*idb) = frc.ma[idh](ifr, idf + 6*idb)*sin(frc.ph[idh](ifr, idf + 6*idb));
					}
					dd = 0;
					line = in.GetLine();
					f.Load(line);
				}
			}
		} else if (line.Find("WAVE PERIOD") >= 0 && line.Find("WAVE FREQUENCY") >= 0) {
			f.Load(line);
			double freq = f.GetDouble(7);
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
						hd().A[6*idb + idf][6*idb + jdf][ifr] = f.GetDouble(1 + jdf);
				}
				in.GetLine(8);
				for (int idf = 0; idf < 6; ++idf) {
					in.GetLine();
					f.Load(in.GetLine());
					if (f.GetText(0) != textDOF[idf])
						throw Exc(in.Str() + "\n"  + Format(t_("Expected %s data, found '%s'"), textDOF[idf], f.GetText())); 
					for (int jdf = 0; jdf < 6; ++jdf) 
						hd().B[6*idb + idf][6*idb + jdf][ifr] = f.GetDouble(1 + jdf);
				}
			}
		}
	}
	//hd().Initialize_Forces(hd().ex);
	//hd().GetFexFromFscFfk();
		
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
		//if (Nh != hd().Nh)
		//	throw Exc(in.Str() + "\n"  + Format(t_("%d headings found when num are %d"), Nh, hd().Nh));
		int Nf = f.GetInt(2);
		//if (Nf != hd().Nf)
		//	throw Exc(in.Str() + "\n"  + Format(t_("%d frequencies found when num are %d"), Nf, hd().Nf));
		
		int col = 3;
		int ih = 0;
		while (!in.IsEof()) {		// Check headings
			while (col < f.size() && ih < Nh) {
				double head = f.GetDouble(col++);
				FindAddRatio(hd().qtfhead, head, 0.001);
				//if (FindRatio(hd().head, head, 0.001) < 0)	
				//	throw Exc(in.Str() + "\n"  + Format(t_("%f heading not found"), head));
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
				//if (FindRatio(hd().w, w, 0.001) < 0)	
				//	throw Exc(in.Str() + "\n"  + Format(t_("%f frequency not found"), w));
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
				qtfsum.fim[idof] = f.GetDouble(idof);
		
			for (int idof = 0; idof < 6; ++idof) {
				qtfsum.fma[idof] = sqrt(sqr(qtfsum.fre[idof]) + sqr(qtfsum.fim[idof])); 
				qtfsum.fph[idof] = atan2(qtfsum.fim[idof], qtfsum.fre[idof]);
			}
		}
	}
	return true;
}

void Aqwa::Save(String ) {
	throw Exc("Option not implemented");
}		
