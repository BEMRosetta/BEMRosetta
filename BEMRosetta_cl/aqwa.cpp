// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"

#include "FastOut.h"

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
	LineParser f(in);
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
		FindAdd(hd().head, FixHeading_180(f.GetDouble(i)));
	Sort(hd().head);
	
	hd().names.SetCount(hd().Nb);
	hd().cg.setConstant(3, hd().Nb, NaNDouble);
	hd().c0.setConstant(3, hd().Nb, NaNDouble);
	hd().C.SetCount(hd().Nb);
	for (int ib = 0; ib < hd().Nb; ++ib) 
		hd().C[ib].setConstant(6, 6, NaNDouble); 
	hd().Initialize_AB(hd().A);
	hd().Initialize_AB(hd().B);
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
	                    VectorXd ma(6), ph(6);
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
			                ma(idf) = f.GetDouble(idf + 3);
	                  	}
	                  	f.Load(in.GetLine());
	                  	for (int idf = 0; idf < 6; ++idf) 
	                       	ph(idf) = -f.GetDouble(idf)*M_PI/180;
	                  	for (int idf = 0; idf < 6; ++idf) 
	                       	hd().ex.force[ih](ifr, idf + 6*ib) = std::complex<double>(ma(idf), ph(idf));
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
	LineParser f(in);
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
	

	factorMass = 1, factorLength = 1;
				
	in.SeekPos(fpos);			
				
	while(!in.IsEof()) {
		line = TrimBoth(in.GetLine());
		f.Load(line);
		
		if (line.StartsWith("WATER  DEPTH")) 
			hd().h = f.GetDouble(f.LAST)*factorLength;
		else if (line.StartsWith("DENSITY  OF  WATER")) 
			hd().rho = f.GetDouble(f.LAST)*factorMass;
		else if (line.StartsWith("ACCELERATION  DUE  TO GRAVITY")) {
			hd().g = f.GetDouble(f.LAST)*factorLength;	
			break;
		} else if ((pos = line.FindAfter("Unit System :")) >= 0) {
			String system = Trim(line.Mid(pos));
			if (system.Find("Metric") < 0)
				throw Exc(in.Str() + "\n" + t_("Only metric system is supported"));
			if (system.Find("kg") > 0)
				factorMass = 1;
			else if (system.Find("tonne") > 0)
				factorMass = 1000;
			else 
				throw Exc(in.Str() + "\n" + t_("Unknown mass unit"));
			if (system.Find("m ") > 0)
				factorLength = 1;
			else if (system.Find("km ") > 0)
				factorLength = 1000;
			else 
				throw Exc(in.Str() + "\n" + t_("Unknown length unit"));
		}
	}
	hd().names.SetCount(hd().Nb);
	hd().Vo.SetCount(hd().Nb, NaNDouble);
	hd().cg.setConstant(3, hd().Nb, NaNDouble);
	hd().c0.setConstant(3, hd().Nb, NaNDouble);
	hd().cb.setConstant(3, hd().Nb, NaNDouble);
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
			double mass = f.GetDouble(2)*factorMass;
			hd().M[ib](0, 0) = hd().M[ib](1, 1) = hd().M[ib](2, 2) = mass;
		} else if (line.StartsWith("CENTRE OF GRAVITY")) {
			f.Load(line);
			hd().cg(0, ib) = f.GetDouble(3)*factorLength;
			hd().cg(1, ib) = f.GetDouble(4)*factorLength;
			hd().cg(2, ib) = f.GetDouble(5)*factorLength;
			hd().c0 = clone(hd().cg);
		} else if (line.StartsWith("INERTIA MATRIX")) {
			Eigen::MatrixXd &M = hd().M[ib];
			f.Load(line);
			M(3, 3) = f.GetDouble(2)*factorMass;
			M(3, 4) = f.GetDouble(3)*factorMass;
			M(3, 5) = f.GetDouble(4)*factorMass;
			in.GetLine();
			f.Load(TrimBoth(in.GetLine()));
			M(4, 3) = f.GetDouble(0)*factorMass;
			M(4, 4) = f.GetDouble(1)*factorMass;
			M(4, 5) = f.GetDouble(2)*factorMass;
			in.GetLine();
			f.Load(TrimBoth(in.GetLine()));
			M(5, 3) = f.GetDouble(0)*factorMass;
			M(5, 4) = f.GetDouble(1)*factorMass;
			M(5, 5) = f.GetDouble(2)*factorMass;
			double mass = M(0, 0);
			double cx = mass*hd().cg(0, ib);
			double cy = mass*hd().cg(1, ib);
			double cz = mass*hd().cg(2, ib);
			M(1, 5) = M(5, 1) =  cx;
			M(2, 4) = M(4, 2) = -cx;
			M(0, 5) = M(5, 0) = -cy;
			M(2, 3) = M(3, 2) =  cy;
			M(0, 4) = M(4, 0) =  cz;
			M(1, 3) = M(3, 1) = -cz;
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
					FindAdd(hd().head, FixHeading_180(f.GetDouble(i)));
				idini = 0;
			}
			hd().Nh = hd().head.size();
			break;
		}
	}
	Sort(hd().head);
	if (IsNull(hd().Nf))
		throw Exc(t_("Number of frequencies not found"));
	if (IsNull(hd().Nh))
		throw Exc(t_("Number of headings not found"));
					
	hd().Initialize_AB(hd().A);
	hd().Initialize_AB(hd().B);
	
	hd().Initialize_Forces(hd().ex);
	hd().Initialize_Forces(hd().fk);
	hd().Initialize_Forces(hd().sc);
	hd().Initialize_Forces(hd().rao);
	
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
				   line.Find("TOTAL HYDROSTATIC STIFFNESS") >= 0) {		// One place to get the hydrostatic stiffness
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
			
			hd().C[ib] *= factorMass;
		} else if (Trim(line) == "STIFFNESS MATRIX") {		// Other place to get the hydrostatic stiffness, but less reliable
			if (hd().C[ib](2, 2) != 0) {					// ... but less reliable
				MatrixXd C;
				C.setConstant(6, 6, 0);
				in.GetLine(6);
				for (int r = 0; r < 6; ++r) {
					f.Load(in.GetLine());
					for (int c = 0; c < 6; ++c) 
						C(r, c) = f.GetDouble(c + 1);
					in.GetLine();
				}
				if (C(0, 0) == 0 && C(1, 1) == 0)			// Only use if ASTF additional hydrostatics has not been used: surge = sway = 0
					hd().C[ib] = C;							
			}
		} else if (line.StartsWith("MESH BASED DISPLACEMENT")) {
			pos = line.FindAfter("=");
			if (pos < 0)
				throw Exc(in.Str() + "\n"  + t_("= not found"));
			hd().Vo[ib] = ScanDouble(line.Mid(pos));
		} else if (line.StartsWith("POSITION OF THE CENTRE OF BUOYANCY")) {
			f.Load(line);
			hd().cb(0, ib) = f.GetDouble(8)*factorLength;	
			hd().cb(1, ib) = f.Load(in.GetLine()).GetDouble(2)*factorLength;
			hd().cb(2, ib) = f.Load(in.GetLine()).GetDouble(2)*factorLength;
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
			else //if (line.StartsWith("R.A.O.S-VARIATION WITH WAVE PERIOD/FREQUENCY")) // Only possible option
				pfrc = &hd().rao;
			Hydro::Forces &frc = *pfrc; 
			
			in.GetLine(5);
			line = in.GetLine();
			while(!in.IsEof()) {
				if (line[0] == '1')
					break;
				
				static const UVector<int> separatorsh = {8,16,26,36,44,54,62,72,80,90,98,108,116,126};
				f.Load(in.GetLine(), separatorsh);

				double heading = FixHeading_180(f.GetDouble(2));
				int idh = FindClosest(hd().head, heading);
				if (idh < 0)
					throw Exc(in.Str() + "\n"  + Format(t_("Heading %f is unknown"), heading));
				int dd = 1;
				double factorM = pfrc == &hd().rao ? 1 : factorMass;
				for (int ifr = 0; ifr < hd().Nf; ++ifr) {
					double freq = f.GetDouble(1);
					int ifrr = FindClosest(hd().w, freq);
					if (ifrr < 0)
						throw Exc(in.Str() + "\n"  + Format(t_("Frequency %f is unknown"), freq));
					for (int idf = 0; idf < 6; ++idf) {
						double factor;
						if (pfrc != &hd().rao || idf < 3)
							factor = 1;
						else
							factor = M_PI/180;	// Only for RAO rotations
						
						frc.force[idh](ifr, idf + 6*ib) = std::polar<double>(f.GetDouble(2 + dd + idf*2)*factor, 
																			-ToRad(f.GetDouble(2 + dd + idf*2 + 1))*factorM); // Negative to follow Wamit 
					}
					dd = 0;
					line = in.GetLine();
					static const UVector<int> separators = {8,16,36,44,54,62,72,80,90,98,108,116,126};
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
			line = in.GetLine();
			int ib2 = -1;
			if (TrimBoth(line) == "ADDED  MASS")
				ib2 = ib;
			else if (line.Find("ADDED MASS FOR FORCE") > 0) {
				int id = line.FindAfter("STR#");
				if (id > 0) {
					ib2 = ScanInt(line.Mid(id));
					if (!IsNull(ib2))
						ib2--;
				}
			}
			if (ib2 >= 0) {
				in.GetLine(5);
			
				for (int idf = 0; idf < 6; ++idf) {
					in.GetLine();
					f.Load(in.GetLine());
					if (f.GetText(0) != textDOF[idf])
						throw Exc(in.Str() + "\n"  + Format(t_("Expected %s data, found '%s'"), textDOF[idf], f.GetText())); 
					for (int jdf = 0; jdf < 6; ++jdf) 
						hd().A[6*ib + idf][6*ib2 + jdf][ifr] = f.GetDouble(1 + jdf)*factorMass;
				}
				in.GetLine(8);
				for (int idf = 0; idf < 6; ++idf) {
					in.GetLine();
					f.Load(in.GetLine());
					if (f.GetText(0) != textDOF[idf])
						throw Exc(in.Str() + "\n"  + Format(t_("Expected %s data, found '%s'"), textDOF[idf], f.GetText())); 
					for (int jdf = 0; jdf < 6; ++jdf) 
						hd().B[6*ib + idf][6*ib2 + jdf][ifr] = f.GetDouble(1 + jdf)*factorMass;
				}
			}
		} else if (line.Find("H Y D R O D Y N A M I C   P A R A M E T E R S   A T   L O W   &   H I G H") >= 0) {	// To capture A0 and Ainf
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
									hd().A0.setConstant(6*hd().Nb, 6*hd().Nb, NaNDouble);
								else if (freq > 90 && hd().Ainf.size() == 0)
									hd().Ainf.setConstant(6*hd().Nb, 6*hd().Nb, NaNDouble);
								break;
							}
						}
						while(!in.IsEof() && !(line[0] == '1')) {
							line = in.GetLine();
							int ib2 = -1;
							if (TrimBoth(line) == "ADDED MASS")
								ib2 = ib;
							else if (line.Find("ADDED MASS FOR FORCE") > 0) {
								int id = line.FindAfter("STR#");
								if (id > 0) {
									ib2 = ScanInt(line.Mid(id));
									if (!IsNull(ib2))
										ib2--;
								}
							}
							if (line.Find("ADDED MASS") >= 0) {
								while(!in.IsEof() && !(line[0] == '1')) {
									f.GetLine();
									if (f.size() == 7 && f.GetText(0) == "X") {
										for (int idf = 0; idf < 6; ++idf) {
											for (int jdf = 0; jdf < 6; ++jdf) {
												if (freq < 0.1)			// 0.001
													hd().A0(6*ib + idf, 6*ib2 + jdf) = f.GetDouble(jdf + 1)*factorMass;
												else if (freq > 90)		// 100.
													hd().Ainf(6*ib + idf, 6*ib2 + jdf) = f.GetDouble(jdf + 1)*factorMass;
											}
											f.GetLine(); 
											f.GetLine();
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
		} else if (Trim(line) == "FREQUENCY INDEPENDENT DAMPING") {
			if (hd().Dlin.size() == 0)
				hd().Dlin = Eigen::MatrixXd::Zero(6*hd().Nb, 6*hd().Nb);
			in.GetLine(6);
			for (int idof = 0; idof < 6; ++idof) {
				f.GetLine();
				for (int jdof = 0; jdof < 6; ++jdof) 
					hd().Dlin(6*ib + idof, 6*ib + jdof) = f.GetDouble(jdof + 1)*factorMass;
				f.GetLine(); 
			}
		} else if (line.Find("W A V E - D R I F T   L O A D S ") >= 0) {
			if (hd().md.size() == 0) {
				hd().mdhead.resize(hd().head.size());
				for (int ih = 0; ih < hd().head.size(); ++ih)
					hd().mdhead[ih] = std::complex<double>(hd().head[ih], hd().head[ih]);
				Hydro::Initialize_MD(hd().md, hd().Nb, int(hd().mdhead.size()), hd().Nf);
			}
			int id;
			while(!in.IsEof() && (id = line.FindAfter("S T R U C T U R E")) < 0) 
				line = in.GetLine();
			int ib = ScanInt(line.Mid(id));
			if (ib < 1 || ib > hd().Nb)
				throw Exc(in.Str() + "\n"  + Format(t_("Wrong body id %d found"), ib)); 
			ib -= 1;
			
			while(!in.IsEof() && (id = f.GetText().FindAfter("(RADIANS/SEC)")) < 0)
				f.GetLine();
			line = f.GetText().Mid(id);				// Because of joined fields like this: 'DUE TO (RADIANS/SEC)-120.0'
			f.Load(line);
			
			UVector<int> idhblock(f.size());
			for (int ih = 0; ih < f.size(); ++ih) {
				double heading = FixHeading_180(f.GetDouble(ih));	
				int id = FindClosest(hd().head, heading);
				if (id < 0)
					throw Exc(in.Str() + "\n"  + Format(t_("Heading %f is unknown"), heading));
				idhblock[ih] = id;
			}

			while(!in.IsEof() && Trim(line) != "DRIFT") 
				line = in.GetLine();
			in.GetLine();
			line = ToLower(Trim(in.GetLine()));	
			int idof;
			if (line.StartsWith("surge"))
				idof = 0;
			else if (line.StartsWith("sway"))
				idof = 1;
			else if (line.StartsWith("heave"))
				idof = 2;
			else if (line.StartsWith("roll"))
				idof = 3;
			else if (line.StartsWith("pitch"))
				idof = 4;
			else if (line.StartsWith("yaw"))
				idof = 5;
			else 
				throw Exc(in.Str() + "\n"  + Format(t_("Wrong DOF %s found"), line));

			for (int idf = 0; idf < hd().Nf; ++idf) {
				f.GetLine();
				for (int ih = 0; ih < idhblock.size(); ++ih) 
					hd().md[ib][idhblock[ih]][idof](idf) = f.GetDouble(1 + ih)*factorMass;
			}
		}
	}
	if (hd().IsLoadedMD()) {
		if (IsNum(hd().md[0][0][2][0]))			
			hd().mdtype = 9;				// Pressure integration/Near field
		else									
			hd().mdtype = 8;				// Momentum conservation/Far field
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
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
	
	f.GetLine();	// AQWA version	
	f.GetLine();
	
	int Nb = f.GetInt(0);
	int Nh = f.GetInt(1);
	int Nf = f.GetInt(2);	
	
	if (!IsNull(hd().Nb) && hd().Nb < Nb)
		throw Exc(in.Str() + "\n"  + Format(t_("Number of bodies loaded is lower than previous (%d != %d)"), hd().Nb, Nb));
	hd().Nb = Nb;
	if (hd().names.IsEmpty())
		hd().names.SetCount(hd().Nb);
			
    UVector<double> _qw;
    UArray<std::complex<double>> _qh;
    
	int col = 3;
	int ih = 0;
	while (!in.IsEof()) {		// Check headings
		while (col < f.size() && ih < Nh) {
			double head = f.GetDouble(col++);
			_qh << std::complex<double>(head, head);
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
			_qw << w;
			ifr++;
		}
		if (ifr >= Nf)
			break;
		f.Load(in.GetLine());
		col = 0;
	}

	if (_qw.size() != Nf)
		throw Exc(in.Str() + "\n"  + Format(t_("Wrong number of frequencies %d found in qtf header. They should have to be %d"), _qw.size(), Nf));
	if (_qh.size() != Nh)
		throw Exc(in.Str() + "\n"  + Format(t_("Wrong number of headings %d found in qtf header. They should have to be %d"), _qh.size(), Nh));
						
	Copy(_qw, hd().qw);
	Copy(_qh, hd().qh);
	
	Hydro::Initialize_QTF(hd().qtfsum, Nb, Nh, Nf);
	Hydro::Initialize_QTF(hd().qtfdif, Nb, Nh, Nf);
	
	int nrows = Nb*Nh*Nf*Nf;
	
	for (int i = 0; i < nrows; ++i) {
		f.Load(in.GetLine());
		int ib = f.GetInt(0)-1;
		if (ib >= Nb)
			throw Exc(in.Str() + "\n"  + Format(t_("Body id %d higher than number of bodies"), ib+1, Nb));
		int ih = f.GetInt(1)-1;
		if (ih >= Nh)
			throw Exc(in.Str() + "\n"  + Format(t_("Heading id %d higher than number of headings"), ih+1, Nh));
		int ifr1 = f.GetInt(2)-1;
		if (ifr1 >= Nf)
			throw Exc(in.Str() + "\n"  + Format(t_("Frequency id %d higher than number of frequencies"), ifr1+1, Nf));
		int ifr2 = f.GetInt(3)-1;
		if (ifr2 >= Nf)
			throw Exc(in.Str() + "\n"  + Format(t_("Frequency id %d higher than number of frequencies"), ifr2+1, Nf));

		for (int idf = 0; idf < 6; ++idf) 
			hd().qtfdif[ib][ih][idf](ifr1, ifr2).real(f.GetDouble(4 + idf)*factorMass);
			
        f.Load(in.GetLine());
        for (int idf = 0; idf < 6; ++idf)
            hd().qtfdif[ib][ih][idf](ifr1, ifr2).imag(-f.GetDouble(idf)*factorMass);	// Negative to follow Wamit
        
		f.Load(in.GetLine());
        for (int idf = 0; idf < 6; ++idf)
            hd().qtfsum[ib][ih][idf](ifr1, ifr2).real(f.GetDouble(idf)*factorMass);
        
        f.Load(in.GetLine());
        for (int idf = 0; idf < 6; ++idf)
            hd().qtfsum[ib][ih][idf](ifr1, ifr2).imag(-f.GetDouble(idf)*factorMass);
	}
	
	return true;
}

bool Aqwa::Save(String file, Function <bool(String, int)> Status) {
	try {
		BEM::Print("\n\n" + Format(t_("Saving '%s'"), file));

		if (hd().IsLoadedQTF(true) || hd().IsLoadedQTF(false)) {
			BEM::Print("\n- " + S(t_("QTF file")));
			Save_QTF(ForceExt(file, ".qtf"), Status);
		}
	} catch (Exc e) {
		BEM::PrintError(Format("\n%s: %s", t_("Error"), e));
		hd().lastError = e;
		return false;
	}
	return true;
}		

void Aqwa::Save_QTF(String file, Function <bool(String, int)> Status) {
	FileOut out(file);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'"), file));
	
	out << "AQTF-2.0 :                                                            \n";
	out << Format(" %2d %2d %3d    ", hd().Nb, int(hd().qh.size()), int(hd().qw.size()));
	int icol = 0;
	for (int ih = 0; ih < hd().qh.size(); ++ih) {
		if (hd().qh[ih].real() != hd().qh[ih].imag())
			continue;
		if (icol > 5) {
			out << "\n              ";
			icol = 0;
		} 
		out << Format("   % 9.5f", hd().qh[ih].real());
		icol++;
	}
	out << "\n               ";
	icol = 0;
	for (int ifr = 0; ifr < hd().qw.size(); ++ifr) {
		if (icol > 5) {
			out << "\n               ";
			icol = 0;
		}
		out << Format("   %9.7f", hd().qw[ifr]);
		icol++;
	}	

	int num = int(hd().Nb*hd().qh.size());
	int inum = 0;
	for (int ib = 0; ib < hd().Nb; ++ib)
        for (int ih = 0, realih = 0; ih < hd().qh.size(); ++ih) {
            inum++;
            if (Status && !Status(Format("Saving %s", file), (100*ih)/int(hd().qh.size())))
				throw Exc(t_("Stop by user"));
            
	        if (hd().qh[ih].real() != hd().qh[ih].imag())
				continue; 
	        for (int ifr1 = 0; ifr1 < hd().qw.size(); ++ifr1) 
				for (int ifr2 = 0; ifr2 < hd().qw.size(); ++ifr2) {
					bool nosum = hd().qtfsum.size() <= ib || hd().qtfsum[ib].size() <= ih || hd().qtfsum[ib][ih][0].rows() <= ifr1 || hd().qtfsum[ib][ih][0].cols() <= ifr2;
			        bool nodif = hd().qtfdif.size() <= ib || hd().qtfdif[ib].size() <= ih || hd().qtfdif[ib][ih][0].rows() <= ifr1 || hd().qtfdif[ib][ih][0].cols() <= ifr2;
			
					out << Format("\n %2d %2d %3d %3d ", ib+1, realih+1, ifr1+1, ifr2+1);
					for (int idf = 0; idf < 6; ++idf) 
						out << Format(" % 6.4E", nodif ? 0 : hd().F_dim(hd().qtfdif[ib][ih][idf](ifr1, ifr2), idf).real());
					out << "\n               ";
					for (int idf = 0; idf < 6; ++idf) 
						out << Format(" % 6.4E", nodif ? 0 : hd().F_dim(hd().qtfdif[ib][ih][idf](ifr1, ifr2), idf).imag());
					out << "\n               ";
					for (int idf = 0; idf < 6; ++idf) 
						out << Format(" % 6.4E", nosum ? 0 : hd().F_dim(hd().qtfsum[ib][ih][idf](ifr1, ifr2), idf).real());
					out << "\n               ";
					for (int idf = 0; idf < 6; ++idf) 
						out << Format(" % 6.4E", nosum ? 0 : hd().F_dim(hd().qtfsum[ib][ih][idf](ifr1, ifr2), idf).imag());
				}
			realih++;	
        }
}
	
	
bool AQWACase::Load(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	solver = AQWA;
	
	String line;
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
	
	bodies.SetCount(1);
	BEMBody &body = bodies[0];
	
	body.meshFile = fileName;
	body.ndof = 6;
	Resize(body.dof, 6, true);
	body.mass = MatrixXd::Zero(6, 6);
	
	UVector<double> hrtz, head;
	
	while (!in.IsEof()) {
		f.GetLine();
		
		if (f.size() == 1) {
			if (f.GetText(0) == "MATE") {
				f.GetLine();
				body.mass(0, 0) = body.mass(1, 1) = body.mass(2, 2) = f.GetDouble(2);
			} else if (f.GetText(0) == "GEOM") {
				while (true) {
					f.GetLine();
					if (f.IsEof())
						break;
					if (f.size() == 1 && f.GetText(0) == "END")
						break;
					if (f.size() > 1 && f.GetText(0) == "1PMAS") {
						body.mass(3, 3) = f.GetDouble(2);
						body.mass(3, 4) = body.mass(4, 3) = f.GetDouble(3);
						body.mass(3, 5) = body.mass(5, 3) = f.GetDouble(4);
						body.mass(4, 4) = f.GetDouble(5);
						body.mass(4, 5) = body.mass(5, 4) = f.GetDouble(6);
						body.mass(5, 5) = f.GetDouble(7);
					}
				}
			} else if (f.GetText(0) == "WFS1") {
				const UVector<int> pos = {0, 10, 15, 20, 30, 40, 50, 60, 70, 80};
				f.GetLineFields(pos);
				for (int r = 0; r < 6; ++r) {
					if (f.GetText(0) != "ASTF")		// Additional Hydrostatic Stiffness Matrix (sometimes due to mooring)
						break;
					if (f.GetInt(2) != r+1)
						throw Exc(Format("Wrong row id '%s' in ASTF", f.GetText(1)));
					for (int c = 0; c < 6; ++c) 
						body.Cadd(r, c) = f.GetDouble(3+c);
					f.GetLineFields(pos);
				}
				for (int r = 0; r < 6; ++r) {
					if (f.GetText(0) != "SSTF")		//  Additional Structural Stiffness Matrix (normally due to mooring)
						break;
					if (f.GetInt(2) != r+1)
						throw Exc(Format("Wrong row id '%s' in SSTF", f.GetText(1)));
					for (int c = 0; c < 6; ++c) 
						body.Cext(r, c) = f.GetDouble(3+c);
					f.GetLineFields(pos);
				}
				for (int r = 0; r < 6; ++r) {
					if (f.GetText(0) != "FIDP")		// Frequency independent Damping Matrix
						break;
					if (f.GetInt(2) != r+1)
						throw Exc(Format("Wrong row id '%s' in FIDP", f.GetText(1)));
					for (int c = 0; c < 6; ++c) 
						body.Dlin(r, c) = f.GetDouble(3+c);
					f.GetLineFields(pos);
				}
				for (int r = 0; r < 6; ++r) {
					if (f.GetText(0) != "FIAM")		// Frequency independent Added Mass Matrix
						break;
					if (f.GetInt(2) != r+1)
						throw Exc(Format("Wrong row id '%s' in FIAM", f.GetText(1)));
					for (int c = 0; c < 6; ++c) 
						body.Aadd(r, c) = f.GetDouble(3+c);
					f.GetLineFields(pos);
				}	
			}
		} else if (f.size() == 2) {
			if (f.GetText(0) == "DPTH")
				h = f.GetDouble(1);
			else if (f.GetText(0) == "DENS")
				rho = f.GetDouble(1);
			else if (f.GetText(0) == "ACCG")
				g = f.GetDouble(1);	
		} else if (f.size() == 4) {
			if (f.GetText(0) == "1HRTZ")
				hrtz << f.GetDouble(3);
			else if (f.GetText(0) == "1DIRN") 
				FindAddDelta(head, FixHeading_180(f.GetDouble(3)), 0.01);
			else if (f.GetInt_nothrow(0) == 198000) {
				body.cg[0] = f.GetDouble(1);
				body.cg[1] = f.GetDouble(2);
				body.cg[2] = f.GetDouble(3);
				body.c0 = clone(body.cg);
			}
		}
	}
	Sort(head);
	
	Nf = hrtz.size();
	if (Nf <= 0)
		throw Exc("Number of frequencies should have to be higher than zero");
	Nh = head.size();
	if (Nh <= 0)
		throw Exc("Number of headings should have to be higher than zero");
	
	minF = fround(hrtz[0]*2*M_PI, 5);				// 5 decimals is enough to filter conversion
	maxF = fround(hrtz[hrtz.size()-1]*2*M_PI, 5);
	
	minH = head[0];
	maxH = head.Top();
	
	return true;
}

			
String FastOut::LoadLis(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return t_("Impossible to open file");
	
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
	
	Clear();
	
	parameters << "time";
	units << "s";	

	auto Add3_dof =[&](String par) {
		parameters << (par + "_x");		parameters << (par + "_y");		parameters << (par + "_z");
	};
	auto Add_3dof =[&](String par) {
		parameters << (par + "_rx");	parameters << (par + "_ry");	parameters << (par + "_rz");
	};
	auto AddUnits6 =[&](String par) {
		UVector<UVector<String>> sunits = {{"POSITION", "m", "deg"}, {"VELOCITY", "m/s", "deg/s"}, {"ACCEL", "m/s2", "deg/s2"},
			{"GRAVITY", "N", "Nm"}, {"HYDROSTATIC", "N", "Nm"}, {"MORISON", "N", "Nm"}, {"KRYLOV", "N", "Nm"}, {"SLAM", "N", "Nm"},
			{"DIFFRACTION", "N", "Nm"}, {"MOORING", "N", "Nm"}, {"DAMPING", "N", "Nm"}, {"GYROSCOPIC", "N", "Nm"},
			{"FORCE", "N", "Nm"}, {"DRAG", "N", "Nm"}, {"WIND", "m/s", "m/s"}, {"ERROR", "%", "%"}};
		for (int i = 0; i < sunits.size(); ++i) {
			if (PatternMatch("*" + sunits[i][0] + "*", par)) {
				units << sunits[i][1];	units << sunits[i][1];	units << sunits[i][1];
				units << sunits[i][2];	units << sunits[i][2];	units << sunits[i][2];
				return;
			}
		}
		units << "";	units << "";	units << "";	units << "";	units << "";	units << "";
	};
	double factorMass = 1, factorLength = 1;
	
	try {
		String line;
		line = in.GetLine();
		
		if (!Trim(line).StartsWith("*********1*********2*********3"))
			return t_("Format error in AQWA file");	// To detect AQWA format
		
		FileInLine::Pos fpos;
	
		bool endwhile = false;
		while(!f.IsEof() && !endwhile) {
			fpos = in.GetPos();
			line = Trim(in.GetLine());
			int pos;
			if (line.StartsWith("*")) {
				if ((pos = line.FindAfter("Unit System :")) >= 0) {
					String system = Trim(line.Mid(pos));
					if (system.Find("Metric") < 0)
						return(in.Str() + "\n" + t_("Only metric system is supported"));
					if (system.Find("kg") > 0)
						factorMass = 1;
					else if (system.Find("tonne") > 0)
						factorMass = 1000;
					else 
						return(in.Str() + "\n" + t_("Unknown mass unit"));
					if (system.Find("m ") > 0)
						factorLength = 1;
					else if (system.Find("km ") > 0)
						factorLength = 1000;
					else 
						return(in.Str() + "\n" + t_("Unknown length unit"));
				}
			} else if (line.StartsWith("RECORD NO.")) {
				bool inparameters = true;
				while(!f.IsEof()) {		
					String str = f.GetLine();
					if (f.GetCount() > 3) {
						String par = Trim(str.Mid(20, 47-20));
						int id;
						if ((id = par.FindAfter("NODE")) > 0) {
							inparameters = false;
							String node = Trim(par.Mid(id));
							Add3_dof(S("Node_pos_") + node); units << "m";		units << "m";		units << "m";
							Add3_dof(S("Node_vel_") + node); units << "m/s";	units << "m/s";		units << "m/s";
							Add3_dof(S("Node_acc_") + node); units << "m/s2";	units << "m/s2";	units << "m/s2";
						} else if (f.GetText(0) == "FORCE") {
							inparameters = false;
							String line = f.GetText(2);
							Add3_dof(S("Force_") + line);	   units << "N";	units << "N";		units << "N";
							parameters << ("Tension_" << line);units << "N";
						} else if (f.GetText(0) == "TENSION") {
							inparameters = false;
							while(!f.IsEof() && f.size() >= 3) {		
								String line = f.GetText(2);
								do {
									String joint;
									if (f.size() == 6)
									 	joint = Trim(f.GetText(4));
									 else
									    joint = Trim(f.GetText(1)); 
									parameters << ("Tension_" + line + "_" + joint);	units << "N";
									f.GetLine();
								} while (f.size() == 3);
							}
							endwhile = true;
							break;
						} else if (inparameters) {
							par = Replace(par, " ", "_");
							Add3_dof(par);
							Add_3dof(par);
							AddUnits6(par);
						}
					}
				}
			} else if (line.StartsWith("WAVE AMPLITUDE")) {
				f.Load(line);
				Hs = 2*f.GetDouble(f.size()-1);
			} else if (line.StartsWith("WAVE PERIOD")) {
				f.Load(line);
				Tp = f.GetDouble(f.size()-1);
			} else if (line.StartsWith("WAVE DIRECTION")) {
				f.Load(line);
				heading = f.GetDouble(f.size()-1);
			}
		}
		if (parameters.size() != units.size()) 
			return t_("Number of parameters and units do not match");
		
		if (factorMass == 1000) {
			for (String &s : units) {
				s = Replace(s, "N", "kN");
				s = Replace(s, "kg", "t");
			}
		}
		if (factorLength == 1000) {
			for (String &s : units)
				s = Replace(s, "m", "km");
		}
		
		dataOut.SetCount(parameters.size());
		
		in.SeekPos(fpos);
		
		while(!f.IsEof()) {
			line = f.GetLine();
			if (line.StartsWith(" RECORD NO.")) {
				int c = 0;
				while(!f.IsEof()) {		// Look for time
					String str = Trim(f.GetLine());
					if (!str.IsEmpty()) {
						double time = ScanDouble(str);
						dataOut[c++] << time;
						break;
					}
				}
				while(!f.IsEof()) {
					line = f.GetLine();
					if (line.StartsWith("1"))
						break;
					f.LoadFields(line, {20, 49, 62, 75, 90, 103, 116});
					if (f.GetCount() > 3) {
						if (TrimLeft(line).StartsWith("***")) 
							;		// System warning
						else if (TrimLeft(f.GetText(0)).StartsWith("POSITION NODE")) {
							f.Load(line);
							dataOut[c++] << f.GetDouble(3);	dataOut[c++] << f.GetDouble(4);	dataOut[c++] << f.GetDouble(5);
							f.Load(f.GetLine());
							dataOut[c++] << f.GetDouble(1);	dataOut[c++] << f.GetDouble(2);	dataOut[c++] << f.GetDouble(3);
							f.Load(f.GetLine());
							dataOut[c++] << f.GetDouble(1);	dataOut[c++] << f.GetDouble(2);	dataOut[c++] << f.GetDouble(3);
						} else if (TrimLeft(f.GetText(0)).StartsWith("FORCE")) {
							f.LoadFields(line, {48, 61, 75, 88, 91, 103});
							dataOut[c++] << f.GetDouble(0);	dataOut[c++] << f.GetDouble(1);	dataOut[c++] << f.GetDouble(2);	dataOut[c++] << f.GetDouble(4);
						} else if (TrimLeft(f.GetText(0)).StartsWith("TENSION")) {
							f.Load(line);
							dataOut[c++] << f.GetDouble(5);
						} else if (Trim(line).StartsWith("S#")) {
							f.Load(line);
							dataOut[c++] << f.GetDouble(2);
						} else {
							dataOut[c++] << f.GetDouble(1);	dataOut[c++] << f.GetDouble(2);	dataOut[c++] << f.GetDouble(3);
							dataOut[c++] << f.GetDouble(4);	dataOut[c++] << f.GetDouble(5);	dataOut[c++] << f.GetDouble(6);
						}
					}
				}
			}
		}		
	} catch (Exc e) {
		return t_("Parsing error: ") + e;
	}			
	
	if (dataOut.IsEmpty()) 
		return Format(t_("Problem reading '%s'"), fileName); 
	
	if (!IsNull(Hs)) {
		parameters << "Hs";
		units << "m";
		UVector<double> &data = dataOut.Add();
		data.SetCount(First(dataOut).size());
		for (int i = 0; i < data.size(); ++i)
			data[i] = Hs;
	}
	if (!IsNull(Tp)) {
		parameters << "Tp";
		units << "s";
		UVector<double> &data = dataOut.Add();
		data.SetCount(First(dataOut).size());
		for (int i = 0; i < data.size(); ++i)
			data[i] = Tp;
	}
	if (!IsNull(Hs) && !IsNull(Tp)) {
		parameters << "Wave1Elev";
		units << "m";
		UVector<double> &data = dataOut.Add();
		UVector<double> &time = First(dataOut);
		data.SetCount(First(dataOut).size());
		double A = Hs/2;
		double w = 2*M_PI/Tp;
		for (int i = 0; i < data.size(); ++i)
			data[i] = A*cos(w*time[i]);
	}
	return "";	
}