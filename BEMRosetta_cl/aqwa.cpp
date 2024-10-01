// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"

#include "FastOut.h"

const char *textDOF[] = {"X", "Y", "Z", "RX", "RY", "RZ"};

String Aqwa::Load(String file, Function <bool(String, int)> Status, double) {
	dt.file = file;
	dt.name = GetFileTitle(file);
	if (ToLower(dt.name) == "analysis")
		dt.name = GetFileTitle(GetFileFolder(file));	// The folder names the model
	dt.dimen = true;
	dt.len = 1;
	dt.solver = Hydro::AQWA;
	dt.Nb = 0;
	
	try {
		BEM::Print("\n\n" + Format(t_("Loading '%s'"), file));
		
		double factorMass = 1;
		
		if (ToLower(GetFileExt(file)) == ".lis") {
			BEM::Print("\n- " + S(t_("LIS file")));
			if (!Load_LIS(factorMass, Status)) {
				BEM::PrintWarning(S(": ** LIS file ") + t_("Not found") + "**");
				BEM::Print("\n- " + S(t_("AH1 file")));
				if (!Load_AH1()) {
					BEM::PrintWarning(S(": ** AH1 file ") + t_("Not found") + "**");
					if (!AQWABody::LoadDat(dt.msh, *this, file).IsEmpty()) 
						dt.Nh = dt.Nf = 0;
					
				}
			}
		} else {
			BEM::Print("\n- " + S(t_("AH1 file")));
			if (!Load_AH1()) {
				BEM::PrintWarning(S(": ** AH1 file ") + t_("Not found") + "**");
				BEM::Print("\n- " + S(t_("LIS file")));
				if (!Load_LIS(factorMass, Status)) {
					BEM::PrintWarning(S(": ** LIS file ") + t_("Not found") + "**");
					if (!AQWABody::LoadDat(dt.msh, *this, file).IsEmpty()) 
						dt.Nh = dt.Nf = 0;
					
				}
			}
		}
		
		//if (IsNull(dt.Nb))
		//	return false;
		
		BEM::Print("\n- " + S(t_("QTF file")));
		if (!Load_QTF(factorMass)) 
			BEM::Print(S(": ** QTF file ") + t_("Not found") + "**");
		
		if (IsNull(dt.Nh) || dt.Nh <= 0) 
			return t_("No data found");
		
		for (int ib = 0; ib < dt.Nb; ++ib)					// Translates all bodies phase to 0,0
			AddWave(ib, -dt.msh[ib].dt.c0.x, -dt.msh[ib].dt.c0.y, dt.g);	// Phase is translated -c0_x,-c0_y
	
		dt.x_w = dt.y_w = 0;
			
	} catch (Exc e) {
		return e;
	}
	
	return String();
}

bool Aqwa::Load_AH1() {
	String fileName = ForceExtSafer(dt.file, ".AH1");
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;

	dt.Nb = dt.Nh = dt.Nf = Null;
	dt.head.Clear();
	dt.w.Clear();
	dt.rho = dt.g = dt.h = Null;
	//dt.dataFromW = true;
	
	String line;
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
	while(!in.IsEof()) {
		line = in.GetLine();
		
		if (TrimBoth(line) == "*")
			break;
	}
	f.Load(in.GetLine());
	dt.Nb = f.GetInt(0);
	dt.Nh = f.GetInt(1);
	dt.Nf = f.GetInt(2);
	
	int i0 = 3;
	do {
		for (int i = i0; i < f.size(); ++i)
			dt.head << f.GetDouble(i);

		if (dt.head.size() == dt.Nh)
			break;
		else if (dt.head.size() > dt.Nh)
			throw Exc(in.Str() + "\n"  + Format(t_("Number of headings do not match %d<>%d"), dt.Nh, dt.head.size()));
		f.GetLine();
		i0 = 0;
	} while(!in.IsEof());
	
	//dt.names.SetCount(dt.Nb);
	dt.msh.SetCount(dt.Nb);
	//dt.cg.setConstant(3, dt.Nb, NaNDouble);
	//dt.c0.setConstant(3, dt.Nb, NaNDouble);
	//dt.C.SetCount(dt.Nb);
	for (int ib = 0; ib < dt.Nb; ++ib) 
		dt.msh[ib].dt.C.setConstant(6, 6, NaNDouble); 
	Initialize_AB(dt.A);
	Initialize_AB(dt.B);
	Initialize_Forces(dt.ex);
	
	while(!in.IsEof() && dt.w.size() < dt.Nf) {
		f.GetLine();
		
		for (int i = 0; i < f.size(); ++i) 
			dt.w << f.GetDouble(i);
	}	
	if (dt.Nf != dt.w.size())
		throw Exc(in.Str() + "\n"  + Format(t_("Number of frequencies do not match %d<>%d"), dt.Nf, dt.w.size()));

	while(!in.IsEof()) {
		line = TrimBoth(in.GetLine());
		if (line.StartsWith("GENERAL")) {
			f.Load(in.GetLine());
			dt.h   = f.GetDouble(1);
			dt.rho = f.GetDouble(2);
			dt.g   = f.GetDouble(3);
		} else if (line.StartsWith("COG")) {
			for (int i = 0; i < dt.Nb; ++i) {
				f.Load(in.GetLine());
				int ib = f.GetInt(0) - 1;
				if (ib < 0 || ib >= dt.Nb)
					throw Exc(in.Str() + "\n"  + Format(t_("Unknown body found in COG %d (%d)"), ib+1, dt.Nb));
				
				dt.msh[ib].dt.cg.x = f.GetDouble(1);
				dt.msh[ib].dt.cg.y = f.GetDouble(2);
				dt.msh[ib].dt.cg.z = f.GetDouble(3);
				dt.msh[ib].dt.c0 = clone(dt.msh[ib].dt.cg);
			}
		} else if (line.StartsWith("COB")) {
			for (int i = 0; i < dt.Nb; ++i) {
				f.Load(in.GetLine());
				int ib = f.GetInt(0) - 1;
				if (ib < 0 || ib >= dt.Nb)
					throw Exc(in.Str() + "\n"  + Format(t_("Unknown body found in COB %d (%d)"), ib+1, dt.Nb));
				
				dt.msh[ib].dt.cb.x = f.GetDouble(1);
				dt.msh[ib].dt.cb.y = f.GetDouble(2);
				dt.msh[ib].dt.cb.z = f.GetDouble(3);
			}
		} else if (line.StartsWith("DISPLACEMENT")) {
			for (int i = 0; i < dt.Nb; ++i) {
				f.Load(in.GetLine());
				int ib = f.GetInt(0) - 1;
				if (ib < 0 || ib >= dt.Nb)
					throw Exc(in.Str() + "\n"  + Format(t_("Unknown body found in COB %d (%d)"), ib+1, dt.Nb));
				
				dt.msh[ib].dt.Vo = f.GetDouble(1);
			}
		} else if (line.StartsWith("HYDSTIFFNESS")) {
			for (int ib = 0; ib < dt.Nb; ++ib) {
				dt.msh[ib].dt.C.setConstant(6, 6, 0);
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
	                    dt.msh[ib].dt.C(idf, jdf) = f.GetDouble(jdf + did);
	            }
			}
		} else if (line.StartsWith("MASS")) {
			for (int ib = 0; ib < dt.Nb; ++ib) {
				dt.msh[ib].dt.M.setConstant(6, 6, 0);
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
	                    dt.msh[ib].dt.M(idf, jdf) = f.GetDouble(jdf + did);
	            }
			}
		} else if (line == "ADDEDMASS" || line == "DAMPING") {
			bool am = line == "ADDEDMASS";
			String sarea = am ? "ADDEDMASS" : "DAMPING";
			for (int ib0 = 0; ib0 < dt.Nb; ++ib0) {
				for (int ib1 = 0; ib1 < dt.Nb; ++ib1) {
					for (int ifr = 0; ifr < dt.Nf; ++ifr) {
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
			                }
			                for (int jdf = 0; jdf < 6; ++jdf) {
			            		if (am)
			                    	dt.A[6*ib0 + idf][6*ib1 + jdf][ifr] = f.GetDouble(jdf + did);
			                    else
			                        dt.B[6*ib0 + idf][6*ib1 + jdf][ifr] = f.GetDouble(jdf + did);
			                }
			            }
					}
				}
			}
		} else if (line.StartsWith("ADDEDMASS_LF")) {
			dt.A0.setConstant(6*dt.Nb, 6*dt.Nb, NaNDouble);
			for (int ib0 = 0; ib0 < dt.Nb; ++ib0) {
				for (int ib1 = 0; ib1 < dt.Nb; ++ib1) {
		            for (int idf = 0; idf < 6; ++idf) {
		                f.Load(in.GetLine());
		                int did = 0;
		        		if (idf == 0) {
		                    did = 2;
		                    int itb0 = f.GetInt(0);
		            		if (itb0 - 1 != ib0)
		                        throw Exc(in.Str() + "\n"  + Format(t_("Body # does not match in ADDEDMASS_LF %d<>%d"), itb0, ib0+1));
		                    int itb1 = f.GetInt(1);
		            		if (itb1 - 1 != ib1)
		                        throw Exc(in.Str() + "\n"  + Format(t_("Body # does not match in ADDEDMASS_LF %d<>%d"), itb1, ib1+1));
		                }
		                for (int jdf = 0; jdf < 6; ++jdf) 
		                    dt.A0(6*ib0 + idf, 6*ib1 + jdf) = f.GetDouble(jdf + did);
		            }
				}
			}
		} else if (line.StartsWith("ADDEDMASS_HF")) {
			dt.Ainf.setConstant(6*dt.Nb, 6*dt.Nb, NaNDouble);
			for (int ib0 = 0; ib0 < dt.Nb; ++ib0) {
				for (int ib1 = 0; ib1 < dt.Nb; ++ib1) {
		            for (int idf = 0; idf < 6; ++idf) {
		                f.Load(in.GetLine());
		                int did = 0;
		        		if (idf == 0) {
		                    did = 2;
		                    int itb0 = f.GetInt(0);
		            		if (itb0 - 1 != ib0)
		                        throw Exc(in.Str() + "\n"  + Format(t_("Body # does not match in ADDEDMASS_LF %d<>%d"), itb0, ib0+1));
		                    int itb1 = f.GetInt(1);
		            		if (itb1 - 1 != ib1)
		                        throw Exc(in.Str() + "\n"  + Format(t_("Body # does not match in ADDEDMASS_LF %d<>%d"), itb1, ib1+1));
		                }
		                for (int jdf = 0; jdf < 6; ++jdf) 
		                    dt.Ainf(6*ib0 + idf, 6*ib1 + jdf) = f.GetDouble(jdf + did);
		            }
				}
			}
		} else if (line.StartsWith("FORCERAO")) {
	        for (int ib = 0; ib < dt.Nb; ++ib) {
	            for (int ih = 0; ih < dt.Nh; ++ih) {
	                for (int ifr = 0; ifr < dt.Nf; ++ifr) {
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
	                       	dt.ex[ib][ih](ifr, idf) = std::complex<double>(ma(idf), ph(idf));
	                }
	            }
	        }
		}
	}
	return true;
}

bool Aqwa::Load_LIS(double &factorMass, Function <bool(String, int)> Status) {
	String fileName = ForceExtSafer(dt.file, ".LIS");
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;

	Status(t_("Loading"), 0);
	
	dt.Nb = dt.Nh = dt.Nf = Null;
	dt.head.Clear();
	dt.w.Clear();
	//dt.T.Clear();
	dt.rho = dt.g = dt.h = Null;
	//dt.dataFromW = true;
	
	String line; 
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
	int pos;
	
	FileInLine::Pos fpos = in.GetPos();
	
	Upp::Index<int> lidGroups;
	
	dt.Nb = 0;
	while(!in.IsEof()) {
		line = Trim(in.GetLine());
		
		if ((pos = line.FindAfter("F O R   S T R U C T U R E")) >= 0) {
			int ib = ScanInt(line.Mid(pos)); 
			dt.Nb = max(ib, dt.Nb);
		} else if ((pos = line.FindAfter("ELEMENT GROUP NUMBER")) >= 0) {		// Panels in the lid are not loaded
			pos = line.FindAfter("=", pos);
			int idGroup	 = ScanInt(line.Mid(pos));
			lidGroups.FindAdd(idGroup);
		}
	}
	if (dt.Nb == 0)
		throw Exc(t_("Number of bodies not found"));
	
	factorMass = 1;
	double factorLength = 1;
	
	dt.msh.SetCount(dt.Nb);
	for (int ib = 0; ib < dt.Nb; ++ib) 
		dt.msh[ib].dt.C.setConstant(6, 6, 0);
	for (int ib = 0; ib < dt.Nb; ++ib) 
		dt.msh[ib].dt.M.setConstant(6, 6, 0);
	
	double nlines = in.GetLineNumber();
	double step = 0.01;
	int lastIdPot = -1;
	
	in.SeekPos(fpos);			
	
	UArray<Upp::Index<int>> panelIDs(dt.Nb);
	UArray<Upp::Index<int>> nodeIDs(dt.Nb);
		
	while(!in.IsEof()) {
		line = TrimBoth(in.GetLine());
		f.Load(line);
		
		if (line.StartsWith("WATER  DEPTH")) 
			dt.h = f.GetDouble(f.LAST)*factorLength;
		else if (line.StartsWith("DENSITY  OF  WATER")) 
			dt.rho = f.GetDouble(f.LAST)*factorMass;
		else if (line.StartsWith("ACCELERATION  DUE  TO GRAVITY")) {
			dt.g = f.GetDouble(f.LAST)*factorLength;	
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
		} else if ((pos = line.FindAfter("C O O R D I N A T E   D A T A")) >= 0) {
			if (dt.msh.IsEmpty()) 
				dt.msh.SetCount(dt.Nb);
			
			f.GetLine(7);
			while (!in.IsEof()) {
				if (f.GetLine()[0] == '1')
					break;
				int ib = f.GetInt(1) - 1;
				nodeIDs[ib] << f.GetInt(2);		// Node number
				Point3D &p = dt.msh[ib].dt.mesh.nodes.Add();
				p.x = f.GetDouble(3)*factorLength;
				p.y = f.GetDouble(4)*factorLength;
				p.z = f.GetDouble(5)*factorLength;
			}
		} else if ((pos = line.FindAfter("E L E M E N T   T O P O L O G Y   F O R   S T R U C T U R E")) >= 0) {
			if (dt.msh.IsEmpty()) 
				dt.msh.SetCount(dt.Nb);
			
			int ib = ScanInt(line.Mid(pos, 6)); 
			if (IsNull(ib))
				throw Exc(in.Str() + "\n"  + t_("Bad body id"));
			ib--;
			f.GetLine(7);
			while (!in.IsEof()) {
				if (f.GetLine()[0] == '1')
					break;
				
				String type = f.GetText(1);
				if (type == "QPPL" || type == "TPPL") {
					if (!lidGroups.IsEmpty()) {
						int idGroup = f.GetInt(8);
						if (lidGroups.Find(idGroup) >= 0)		// Panels in the lid are not loaded
							continue;
					}
					Panel &p = dt.msh[ib].dt.mesh.panels.Add();
					panelIDs[ib] << f.GetInt(0);
					
					int id;
					id = nodeIDs[ib].Find(f.GetInt(2));
					if (id < 0)
						throw Exc(in.Str() + "\n"  + t_("Node number #1 not found"));
					p.id[0] = id;
					
					id = nodeIDs[ib].Find(f.GetInt(3));
					if (id < 0)
						throw Exc(in.Str() + "\n"  + t_("Node number #2 not found"));
					p.id[1] = id;
					
					id = nodeIDs[ib].Find(f.GetInt(4));
					if (id < 0)
						throw Exc(in.Str() + "\n"  + t_("Node number #3 not found"));
					p.id[2] = id;
					
					id = f.GetInt(5);
					if (id == 0)
						p.id[3] = p.id[0];
					else {	
						id = nodeIDs[ib].Find(id);
						if (id < 0)
							throw Exc(in.Str() + "\n"  + t_("Node number #4 not found"));
						p.id[3] = id;
					}
				}
			} 
		} else if ((pos = line.FindAfter("M O O R I N G   S T I F F N E S S   F O R   S T R U C T U R E")) >= 0) {			// LIBRIUM
			for (int ib = 0; ib < dt.Nb; ++ib) 
				dt.msh[ib].dt.Cmoor.setConstant(6, 6, 0);
			
			int iib = ScanInt(line.Mid(pos));
			if (iib < 1 || iib > dt.Nb)
				throw Exc(in.Str() + "\n"  + t_("Bad body id"));
			iib--;
			f.GetLine(10); 
			for (int idof = 0; idof < 6; ++idof) {
				f.GetLine();
				for (int jdof = 0; jdof < 6; ++jdof) 
					dt.msh[iib].dt.Cmoor(idof, jdof) = f.GetDouble(jdof + 1)*factorMass;
				f.GetLine(); 
			}
		} else if ((pos = line.FindAfter("H Y D R O S T A T I C   S T I F F N E S S   O F   S T R U C T U R E")) >= 0) {	// LIBRIUM	
			int iib = ScanInt(line.Mid(pos));
			if (iib < 1 || iib > dt.Nb)
				throw Exc(in.Str() + "\n"  + t_("Bad body id"));
			iib--;
			f.GetLine(7); 
			for (int idof = 0; idof < 6; ++idof) {
				f.GetLine(3);
				for (int jdof = 0; jdof < 6; ++jdof) 
					dt.msh[iib].dt.C(idof, jdof) = f.GetDouble(jdof + 1)*factorMass;
			}
		} 
	}
	
	if (!dt.msh.IsEmpty()) {
		for (int ib = 0; ib < dt.Nb; ++ib)
			dt.msh[ib].dt.SetCode(Body::AQWA_LIS);
	}
			
	int ib;	
	
	while(!in.IsEof()) {
		line = TrimBoth(in.GetLine());
		
		if ((pos = line.FindAfter("S T R U C T U R E")) >= 0) {
			ib = ScanInt(line.Mid(pos)) - 1; 
			if (ib >= dt.Nb)
				throw Exc(in.Str() + "\n"  + Format(t_("Wrong body %d"), ib));
		} else if (line.StartsWith("PMAS")) {
			f.Load(line);
			double mass = f.GetDouble(2)*factorMass;
			dt.msh[ib].dt.M(0, 0) = dt.msh[ib].dt.M(1, 1) = dt.msh[ib].dt.M(2, 2) = mass;
		} else if (line.StartsWith("CENTRE OF GRAVITY")) {
			f.Load(line);
			dt.msh[ib].dt.cg.x = f.GetDouble(3)*factorLength;
			dt.msh[ib].dt.cg.y = f.GetDouble(4)*factorLength;
			dt.msh[ib].dt.cg.z = f.GetDouble(5)*factorLength;
			dt.msh[ib].dt.c0 = clone(dt.msh[ib].dt.cg);
		} else if (line.StartsWith("INERTIA MATRIX")) {
			Eigen::MatrixXd &M = dt.msh[ib].dt.M;
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
		} else if (IsNull(dt.Nf) && line.Find("W A V E   F R E Q U E N C I E S / P E R I O D S   A N D   D I R E C T I O N S") >= 0) {
			in.GetLine(5);
			dt.Nf = 0;
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
				if (dt.w.size() == 0 || newparagraph) {
					w = f.GetDouble(2);
					newparagraph = false;
				} else
					w = f.GetDouble(1);
				dt.w << w;
			}
			dt.Nf = dt.w.size();
					
		} else if (IsNull(dt.Nh) && line.Find("DIRECTIONS") >= 0) {
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
					//FindAdd(dt.head, FixHeading_0_360(f.GetDouble(i)));
					dt.head << f.GetDouble(i);
				idini = 0;
			}
			dt.Nh = dt.head.size();
			break;
		} 
	}
	//Sort(dt.head);
	
	if (IsNull(dt.Nf))
		dt.Nf = 0;
	if (IsNull(dt.Nh))
		dt.Nh = 0;
	
	Initialize_AB(dt.A);
	Initialize_AB(dt.B);
	
	Initialize_Forces(dt.ex);
	Initialize_Forces(dt.fk);
	Initialize_Forces(dt.sc);
	Initialize_Forces(dt.rao);
	
	UVector<double> bouyancyForce(dt.Nb, Null);
	
	int ifrPot = -1, prevTrans = Null;
	
	while(!in.IsEof()) {
		line = TrimBoth(in.GetLine());
		
		if (Status) {
			double adv = in.GetLineNumber()/nlines;
			if (adv > step) {
				step += 0.01;
				if (Status && !Status("Loading", int(100*adv)))
					throw Exc(t_("Stop by user"));
			}
		}
		
		if (line.FindAfter("* P O T E N T I A L S *") >= 0) {	
			line = in.GetLine(2);
			pos = line.FindAfter("S T R U C T U R E");
			ib = ScanInt(line.Mid(pos)) - 1; 
			if (ib < 0 || ib >= dt.Nb)
				throw Exc(in.Str() + "\n"  + Format(t_("Wrong body %d"), ib));
			
			if (!IsLoadedPotsRad(ib))
				Initialize_PotsRad(); 				// Initialise potentials
			
			int trans;
			line = in.GetLine(4);
			if (line.Find("TRANSLATIONAL FREEDOMS") >= 0)
				trans = true;
			else if (line.Find("ROTATIONAL FREEDOMS") >= 0)
				trans = false;
			else
				trans = Null;
			
			if (!IsNull(trans)) {
				if (dt.Nb == 1)
					in.GetLine(7);
				else {
					ib = ScanInt(Trim(line).Right(3)) - 1;
					if (ib < 0 || ib >= dt.Nb)
						throw Exc(in.Str() + "\n"  + Format(t_("Wrong body %d"), ib));
				
					in.GetLine(6);
				}
				if (trans) {
				 	if (IsNull(prevTrans) || !prevTrans) {
						prevTrans = true;
						if (ib == 0)
							++ifrPot;
					}
				} else
					prevTrans = false;
				
				while (!in.IsEof()) {
					line = in.GetLine();
					if (line[0] == '1')
						break;
					
					if (line.GetCount() < 30)
						continue;
					
					f.Load(line.Mid(20));		// Removed the initial 1:( X, Y,Z)
					
					if (f.size() == 8) {
						int idPanel = f.GetInt(0);
						
						if (lastIdPot >= 0 && lastIdPot+1 < panelIDs[ib].size() && panelIDs[ib][lastIdPot+1] == idPanel)
							lastIdPot++;
						else
							lastIdPot = panelIDs[ib].Find(idPanel);
						
						if (lastIdPot >= 0) {
							UArray<UArray<std::complex<double>>> &pan = dt.pots_rad[ib][lastIdPot];
							if (trans) {
								pan[0][ifrPot] += std::polar<double>(f.GetDouble(2), ToRad(f.GetDouble(3)));
								pan[1][ifrPot] += std::polar<double>(f.GetDouble(4), ToRad(f.GetDouble(5)));
								pan[2][ifrPot] += std::polar<double>(f.GetDouble(6), ToRad(f.GetDouble(7)));
							} else {
								pan[3][ifrPot] += std::polar<double>(f.GetDouble(2), ToRad(f.GetDouble(3)));
								pan[4][ifrPot] += std::polar<double>(f.GetDouble(4), ToRad(f.GetDouble(5)));
								pan[5][ifrPot] += std::polar<double>(f.GetDouble(6), ToRad(f.GetDouble(7)));
							}
						}
					}
				}
			}
		} else if (line.Find("G L O B A L   A D D I T I O N A L   S T R U C T U R E   S T I F F N E S S   M A T R I X") >= 0) {
			f.GetLine(8);
			int iib = f.GetInt(1);
			if (iib < 1 || iib > dt.Nb)
				throw Exc(in.Str() + "\n"  + t_("Bad body id"));
			iib--;
			if (dt.msh[iib].dt.Cmoor.size() == 0)
				dt.msh[iib].dt.Cmoor = Eigen::MatrixXd::Zero(6, 6);
			in.GetLine(4);
			for (int idof = 0; idof < 6; ++idof) {
				f.GetLine();
				for (int jdof = 0; jdof < 6; ++jdof) 
					dt.msh[iib].dt.Cmoor(idof, jdof) = f.GetDouble(jdof + 1)*factorMass;
				f.GetLine(); 
			}
		} else if ((pos = line.FindAfter("S T R U C T U R E")) >= 0) {
			ib = ScanInt(line.Mid(pos));
			if (!IsNull(ib)) {
				ib -= 1; 
				if (ib >= dt.Nb)
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
			dt.msh[ib].dt.C(2, 2) = f.GetDouble(0);
			dt.msh[ib].dt.C(2, 3) = f.GetDouble(1);
			dt.msh[ib].dt.C(2, 4) = f.GetDouble(2);
			if (f.size() > 3)
				dt.msh[ib].dt.C(2, 5) = f.GetDouble(3);	
			line = in.GetLine();
			pos = line.FindAfter("=");
			if (pos < 0)
				throw Exc(in.Str() + "\n"  + t_("Format error, '=' not found"));
			f.Load(line.Mid(pos));
			dt.msh[ib].dt.C(3, 2) = f.GetDouble(0);
			dt.msh[ib].dt.C(3, 3) = f.GetDouble(1);
			dt.msh[ib].dt.C(3, 4) = f.GetDouble(2);
			if (f.size() > 3)
				dt.msh[ib].dt.C(3, 5) = f.GetDouble(3);	
			line = in.GetLine();
			pos = line.FindAfter("=");
			if (pos < 0)
				throw Exc(in.Str() + "\n"  + t_("Format error, '=' not found"));
			f.Load(line.Mid(pos));
			dt.msh[ib].dt.C(4, 2) = f.GetDouble(0);
			dt.msh[ib].dt.C(4, 3) = f.GetDouble(1);
			dt.msh[ib].dt.C(4, 4) = f.GetDouble(2);
			if (f.size() > 3)
	       		dt.msh[ib].dt.C(4, 5) = f.GetDouble(3);	
			
			dt.msh[ib].dt.C *= factorMass;
		} else if (Trim(line) == "STIFFNESS MATRIX") {		// Other place to get the hydrostatic stiffness, 
			if (dt.msh[ib].dt.C(2, 2) == 0) {					// ... but less reliable
				MatrixXd C;
				C.setConstant(6, 6, 0);
				in.GetLine(6);
				for (int r = 0; r < 6; ++r) {
					f.Load(in.GetLine());
					for (int c = 0; c < 6; ++c) 
						C(r, c) = f.GetDouble(c + 1)*factorMass;
					in.GetLine();
				}
				if (C(0, 0) == 0 && C(1, 1) == 0)			// Only use if ASTF additional hydrostatics has not been used: surge = sway = 0
					dt.msh[ib].dt.C = C;							
			}
		} else if (line.StartsWith("BUOYANCY FORCE")) {
			pos = line.FindAfter("=");
			if (pos < 0)
				throw Exc(in.Str() + "\n"  + t_("= not found"));
			bouyancyForce[ib] = ScanDouble(line.Mid(pos));
		} else if (line.StartsWith("MESH BASED DISPLACEMENT")) {
			pos = line.FindAfter("=");
			if (pos < 0)
				throw Exc(in.Str() + "\n"  + t_("= not found"));
			dt.msh[ib].dt.Vo = ScanDouble(line.Mid(pos));
		} else if (line.StartsWith("POSITION OF THE CENTRE OF BUOYANCY")) {
			f.Load(line);
			dt.msh[ib].dt.cb.x = f.GetDouble(8)*factorLength;	
			dt.msh[ib].dt.cb.y = f.Load(in.GetLine()).GetDouble(2)*factorLength;
			dt.msh[ib].dt.cb.z = f.Load(in.GetLine()).GetDouble(2)*factorLength;
		} else if (line.StartsWith("FROUDE KRYLOV + DIFFRACTION FORCES-VARIATION WITH WAVE PERIOD/FREQUENCY") ||
				   line.StartsWith(				 "FROUDE KRYLOV FORCES-VARIATION WITH WAVE PERIOD/FREQUENCY") ||
				   line.StartsWith(  			   "DIFFRACTION FORCES-VARIATION WITH WAVE PERIOD/FREQUENCY") ||
				   line.StartsWith(             			  "R.A.O.S-VARIATION WITH WAVE PERIOD/FREQUENCY")) {
			Hydro::Forces *pfrc;
			if (line.StartsWith("FROUDE KRYLOV + DIFFRACTION FORCES-VARIATION WITH WAVE PERIOD/FREQUENCY"))
				pfrc = &dt.ex;
			else if (line.StartsWith("FROUDE KRYLOV FORCES-VARIATION WITH WAVE PERIOD/FREQUENCY")) 
				pfrc = &dt.fk;
			else if (line.StartsWith("DIFFRACTION FORCES-VARIATION WITH WAVE PERIOD/FREQUENCY")) 
				pfrc = &dt.sc;
			else //if (line.StartsWith("R.A.O.S-VARIATION WITH WAVE PERIOD/FREQUENCY")) // Only possible option
				pfrc = &dt.rao;
			Hydro::Forces &frc = *pfrc; 
			
			in.GetLine(5);
			line = in.GetLine();
			while(!in.IsEof()) {
				if (line[0] == '1')
					break;
				
				static const UVector<int> separatorsh = {8,16,26,36,44,54,62,72,80,90,98,108,116,126};
				f.Load(in.GetLine(), separatorsh);

				int idh = FindClosest(dt.head, f.GetDouble(2));
				if (idh < 0)
					throw Exc(in.Str() + "\n"  + Format(t_("Heading %f is unknown"), f.GetDouble(2)));
				int dd = 1;
				for (int ifr = 0; ifr < dt.Nf; ++ifr) {
					double freq = f.GetDouble(1);
					int ifrr = FindClosest(dt.w, freq);
					if (ifrr < 0)
						throw Exc(in.Str() + "\n"  + Format(t_("Frequency %f is unknown"), freq));
					for (int idf = 0; idf < 6; ++idf) {
						double factorM;
						if (pfrc != &dt.rao || idf < 3)
							factorM = factorMass;
						else
							factorM = factorMass*M_PI/180;	// Only for RAO rotations, conversion from degrees to radians
								
						pos = 2 + dd + idf*2;
						frc[ib][idh](ifr, idf) = std::polar<double>(f.GetDouble(pos)*factorM, 
															 -ToRad(f.GetDouble(pos + 1))); // Negative to follow Wamit 
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
			int ifr = FindClosest(dt.w, freq);
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
						dt.A[6*ib + idf][6*ib2 + jdf][ifr] = f.GetDouble(1 + jdf)*factorMass;
				}
				in.GetLine(8);
				for (int idf = 0; idf < 6; ++idf) {
					in.GetLine();
					f.Load(in.GetLine());
					if (f.GetText(0) != textDOF[idf])
						throw Exc(in.Str() + "\n"  + Format(t_("Expected %s data, found '%s'"), textDOF[idf], f.GetText())); 
					for (int jdf = 0; jdf < 6; ++jdf) 
						dt.B[6*ib + idf][6*ib2 + jdf][ifr] = f.GetDouble(1 + jdf)*factorMass;
				}
			}
		} else if (line.Find("H Y D R O D Y N A M I C   P A R A M E T E R S   A T   L O W   &   H I G H") >= 0) {	// To capture A0 and Ainf
			int iib = -1;
			double freq = -1;
			while(!in.IsEof() && !(line[0] == '1')) {
				line = in.GetLine();
				int id = line.FindAfter("F O R   S T R U C T U R E");
				if (id >= 0) {
					iib = ScanInt(line.Mid(id)) - 1;
					for (int ii = 0; ii < 2; ++ii) {
						while(!in.IsEof() && !(line[0] == '1')) {
							line = in.GetLine();
							int idr = line.FindAfter("WAVE FREQUENCY");
							if (idr >= 0) {
								freq = ScanDouble(line.Mid(idr + 3));
								if (freq < 0.1 && dt.A0.size() == 0)
									dt.A0.setConstant(6*dt.Nb, 6*dt.Nb, NaNDouble);
								else if (freq > 90 && dt.Ainf.size() == 0)
									dt.Ainf.setConstant(6*dt.Nb, 6*dt.Nb, NaNDouble);
								break;
							}
						}
						while(!in.IsEof() && !(line[0] == '1')) {
							line = in.GetLine();
							int ib2 = -1;
							if (TrimBoth(line) == "ADDED MASS")
								ib2 = iib;
							else if (line.Find("ADDED MASS FOR FORCE") > 0) {
								id = line.FindAfter("STR#");
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
													dt.A0(6*iib + idf, 6*ib2 + jdf) = f.GetDouble(jdf + 1)*factorMass;
												else if (freq > 90)		// 100.
													dt.Ainf(6*iib + idf, 6*ib2 + jdf) = f.GetDouble(jdf + 1)*factorMass;
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
			if (dt.msh[ib].dt.Dlin.size() == 0)
				dt.msh[ib].dt.Dlin = Eigen::MatrixXd::Zero(6, 6);
			in.GetLine(6);
			for (int idof = 0; idof < 6; ++idof) {
				f.GetLine();
				for (int jdof = 0; jdof < 6; ++jdof) 
					dt.msh[ib].dt.Dlin(idof, jdof) = f.GetDouble(jdof + 1)*factorMass;
				f.GetLine(); 
			}
		} else if (line.Find("W A V E - D R I F T   L O A D S ") >= 0) {
			if (dt.md.size() == 0) {
				dt.mdhead.resize(dt.head.size());
				for (int ih = 0; ih < dt.head.size(); ++ih)
					dt.mdhead[ih] = std::complex<double>(dt.head[ih], dt.head[ih]);
				Hydro::Initialize_MD(dt.md, dt.Nb, int(dt.mdhead.size()), dt.Nf);
			}
			int id;
			while(!in.IsEof() && (id = line.FindAfter("S T R U C T U R E")) < 0) 
				line = in.GetLine();
			int iib = ScanInt(line.Mid(id));
			if (iib < 1 || iib > dt.Nb)
				throw Exc(in.Str() + "\n"  + Format(t_("Wrong body id %d found"), ib)); 
			iib -= 1;
			
			while(!in.IsEof() && (id = f.GetText().FindAfter("(RADIANS/SEC)")) < 0)
				f.GetLine();
			line = f.GetText().Mid(id);				// Because of joined fields like this: 'DUE TO (RADIANS/SEC)-120.0'
			f.Load(line);
			
			UVector<int> idhblock(f.size());
			for (int ih = 0; ih < f.size(); ++ih) {
				id = FindClosest(dt.head, f.GetDouble(ih));
				if (id < 0)
					throw Exc(in.Str() + "\n"  + Format(t_("Heading %f is unknown"), f.GetDouble(ih)));
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

			for (int idf = 0; idf < dt.Nf; ++idf) {
				f.GetLine();
				for (int ih = 0; ih < idhblock.size(); ++ih) 
					dt.md[iib][idhblock[ih]][idof](idf) = f.GetDouble(1 + ih)*factorMass;
			}
		}
	}
	if (IsLoadedMD()) {
		if (IsNum(dt.md[0][0][2][0]))			
			dt.mdtype = 9;				// Pressure integration/Near field
		else									
			dt.mdtype = 8;				// Momentum conservation/Far field
	}
	if (IsLoadedQTF(true)) {
		if (IsNum(dt.qtfsum[0][0][2](0, 0)))			
			dt.qtftype = 9;				// Pressure integration/Near field
		else									
			dt.qtftype = 8;				// Momentum conservation/Far field
	} else if (IsLoadedQTF(false)) {
		if (IsNum(dt.qtfdif[0][0][2](0, 0)))			
			dt.qtftype = 9;				// Pressure integration/Near field
		else									
			dt.qtftype = 8;				// Momentum conservation/Far field
	}
		
	if (!dt.msh.IsEmpty()) {
		// Removes mooring, Morison and other points unrelated with panels
		dt.symX = dt.symY = false;
		for (Body &m : dt.msh) {
			if (IsNull(m.dt.Vo) && !IsNull(bouyancyForce[ib]))
				m.dt.Vo = bouyancyForce[ib]/dt.rho/dt.g;
			
			m.dt.mesh.GetPanelParams();
			Surface::RemoveDuplicatedPointsAndRenumber(m.dt.mesh.panels, m.dt.mesh.nodes);
			
			m.dt.under.CutZ(m.dt.mesh, -1);
			m.dt.under.GetVolume();
			Point3D cb = m.dt.under.GetCentreOfBuoyancy();
			if (!IsNull(m.dt.cb)) {								// Guess the symmetries...
				if (abs(m.dt.cb.x) < 0.001 && abs(cb.x) > 0.001)// ...comparing real cb vs the obtained through the mesh
					dt.symX = true;
				if (abs(m.dt.cb.y) < 0.001 && abs(cb.y) > 0.001)
					dt.symY = true;
			} else if (!IsNull(m.dt.Vo)) {						// Guess the symmetries...
				double ratio = m.dt.Vo/m.dt.under.volume;		// ...depending on the real displacement vs the obtained through the mesh
				if (EqualRatio(ratio, 4., 0.01))					// Double symmetry
					dt.symX = dt.symY = true;					
				else if (EqualRatio(ratio, 2., 0.01)) {				// Simple symmetry
					if (abs(cb.x) > 0.001 && abs(cb.y) > 0.001)
						dt.symY = true;
					else if (abs(cb.y) > 0.001 && abs(cb.x) > 0.001)
						dt.symX = true;
				}
			}
		}
		
		// Checks the symmetry
		/*dt.symX = dt.symY = false;
		
		Surface full = clone(dt.msh[0].dt.mesh);	
		for (int i = 1; i < dt.msh.size(); ++i)
			full.Append(dt.msh[i].dt.mesh);
		
		full.GetPanelParams();
		VolumeEnvelope env = full.GetEnvelope();
	
		if (abs(env.minX) < EPS_LEN || abs(env.maxX) < EPS_LEN) {
			double xProjectionPos = full.GetAreaXProjection(true, false);
			double xProjectionNeg = full.GetAreaXProjection(false, true);
			if (abs(1-abs(xProjectionPos/xProjectionNeg)) > 0.2)
				dt.symX = true;
		}
		if (abs(env.minY) < EPS_LEN || abs(env.maxY) < EPS_LEN) {
			double yProjectionPos = full.GetAreaYProjection(true, false);
			double yProjectionNeg = full.GetAreaYProjection(false, true);
			if (abs(1-abs(yProjectionPos/yProjectionNeg)) > 0.2)
				dt.symY = true;
		}*/
		
		
	}
	if (!dt.pots_rad.IsEmpty()) {		// Transform potentials to Wamit
		for (int iib = 0; iib < dt.Nb; ++iib) 
			for (int ip = 0; ip < dt.pots_rad[iib].size(); ++ip) 
				for (int idf = 0; idf < 6; ++idf)
					for (int ifr = 0; ifr < dt.Nf; ++ifr) {	
						auto &d = dt.pots_rad[iib][ip][idf][ifr];
						d = std::complex<double>(d.imag()/dt.w[ifr], d.real()/dt.w[ifr]);
					}	
	}					
	
	return true;
}

bool Aqwa::Load_QTF(double factorMass) {
	String fileName = ForceExtSafer(dt.file, ".QTF");
	FileInLine in(fileName);
	if (!in.IsOpen()) {
		fileName = AFX(GetFileFolder(fileName), "analysis.qtf"); 
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
	
	if (!IsNull(dt.Nb) && dt.Nb < Nb)
		throw Exc(in.Str() + "\n"  + Format(t_("Number of bodies loaded is lower than previous (%d != %d)"), dt.Nb, Nb));
	dt.Nb = Nb;
	if (dt.msh.IsEmpty())
		dt.msh.SetCount(dt.Nb);
			
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
						
	::Copy(_qw, dt.qw);
	::Copy(_qh, dt.qhead);
	
	Hydro::Initialize_QTF(dt.qtfsum, Nb, Nh, Nf);
	Hydro::Initialize_QTF(dt.qtfdif, Nb, Nh, Nf);
	
	int nrows = Nb*Nh*Nf*Nf;
	
	for (int i = 0; i < nrows; ++i) {
		f.Load(in.GetLine());
		int ib = f.GetInt(0)-1;
		if (ib >= Nb)
			throw Exc(in.Str() + "\n"  + Format(t_("Body id %d higher than number of bodies"), ib+1, Nb));
		ih = f.GetInt(1)-1;
		if (ih >= Nh)
			throw Exc(in.Str() + "\n"  + Format(t_("Heading id %d higher than number of headings"), ih+1, Nh));
		int ifr1 = f.GetInt(2)-1;
		if (ifr1 >= Nf)
			throw Exc(in.Str() + "\n"  + Format(t_("Frequency id %d higher than number of frequencies"), ifr1+1, Nf));
		int ifr2 = f.GetInt(3)-1;
		if (ifr2 >= Nf)
			throw Exc(in.Str() + "\n"  + Format(t_("Frequency id %d higher than number of frequencies"), ifr2+1, Nf));

		for (int idf = 0; idf < 6; ++idf) 
			dt.qtfdif[ib][ih][idf](ifr1, ifr2).real(f.GetDouble(4 + idf)*factorMass);
			
        f.Load(in.GetLine());
        for (int idf = 0; idf < 6; ++idf)
            dt.qtfdif[ib][ih][idf](ifr1, ifr2).imag(-f.GetDouble(idf)*factorMass);	// Negative to follow Wamit
        
		f.Load(in.GetLine());
        for (int idf = 0; idf < 6; ++idf)
            dt.qtfsum[ib][ih][idf](ifr1, ifr2).real(f.GetDouble(idf)*factorMass);
        
        f.Load(in.GetLine());
        for (int idf = 0; idf < 6; ++idf)
            dt.qtfsum[ib][ih][idf](ifr1, ifr2).imag(-f.GetDouble(idf)*factorMass);
	}
	
	return true;
}

void Aqwa::Save(String file, Function <bool(String, int)> Status) const {
	BEM::Print("\n\n" + Format(t_("Saving '%s'"), file));

	if (IsLoadedQTF(true) || IsLoadedQTF(false)) {
		BEM::Print("\n- " + S(t_("QTF file")));
		Save_QTF(ForceExt(file, ".qtf"), Status);
	}
}		

void Aqwa::Save_QTF(String file, Function <bool(String, int)> Status) const {
	FileOut out(file);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'"), file));
	
	out << "AQTF-2.0 :                                                            \n";
	out << Format(" %2d %2d %3d    ", dt.Nb, int(dt.qhead.size()), int(dt.qw.size()));
	int icol = 0;
	for (int ih = 0; ih < dt.qhead.size(); ++ih) {
		if (dt.qhead[ih].real() != dt.qhead[ih].imag())
			continue;
		if (icol > 5) {
			out << "\n              ";
			icol = 0;
		} 
		out << Format("   % 9.5f", dt.qhead[ih].real());
		icol++;
	}
	out << "\n               ";
	icol = 0;
	for (int ifr = 0; ifr < dt.qw.size(); ++ifr) {
		if (icol > 5) {
			out << "\n               ";
			icol = 0;
		}
		out << Format("   %9.7f", dt.qw[ifr]);
		icol++;
	}	

	//int num = int(dt.Nb*dt.qh.size());
	//int inum = 0;
	for (int ib = 0; ib < dt.Nb; ++ib)
        for (int ih = 0, realih = 0; ih < dt.qhead.size(); ++ih) {
            //inum++;
    		if (Status && !Status(Format("Saving %s", file), (100*ih)/int(dt.qhead.size())))
				throw Exc(t_("Stop by user"));
            
			if (dt.qhead[ih].real() != dt.qhead[ih].imag())
				continue; 
	        for (int ifr1 = 0; ifr1 < dt.qw.size(); ++ifr1) 
				for (int ifr2 = 0; ifr2 < dt.qw.size(); ++ifr2) {
					bool nosum = dt.qtfsum.size() <= ib || dt.qtfsum[ib].size() <= ih || dt.qtfsum[ib][ih][0].rows() <= ifr1 || dt.qtfsum[ib][ih][0].cols() <= ifr2;
			        bool nodif = dt.qtfdif.size() <= ib || dt.qtfdif[ib].size() <= ih || dt.qtfdif[ib][ih][0].rows() <= ifr1 || dt.qtfdif[ib][ih][0].cols() <= ifr2;
			
					out << Format("\n %2d %2d %3d %3d ", ib+1, realih+1, ifr1+1, ifr2+1);
					for (int idf = 0; idf < 6; ++idf) 
						out << Format(" % 6.4E", nodif ? 0 : F_dim(dt.qtfdif[ib][ih][idf](ifr1, ifr2), idf).real());
					out << "\n               ";
					for (int idf = 0; idf < 6; ++idf) 
						out << Format(" % 6.4E", nodif ? 0 : F_dim(dt.qtfdif[ib][ih][idf](ifr1, ifr2), idf).imag());
					out << "\n               ";
					for (int idf = 0; idf < 6; ++idf) 
						out << Format(" % 6.4E", nosum ? 0 : F_dim(dt.qtfsum[ib][ih][idf](ifr1, ifr2), idf).real());
					out << "\n               ";
					for (int idf = 0; idf < 6; ++idf) 
						out << Format(" % 6.4E", nosum ? 0 : F_dim(dt.qtfsum[ib][ih][idf](ifr1, ifr2), idf).imag());
				}
			realih++;	
        }
}
			
String FastOut::Load_LIS(String file) {
	FileInLine in(file);
	if (!in.IsOpen())
		return t_("Impossible to open file");
	
	fileName = file;
	
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
							line = f.GetText(2);
							Add3_dof(S("Force_") + line);	   units << "N";	units << "N";		units << "N";
							parameters << ("Tension_" << line);units << "N";
						} else if (f.GetText(0) == "TENSION") {
							inparameters = false;
							while(!f.IsEof() && f.size() >= 3) {		
								line = f.GetText(2);
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

void Aqwa::SaveCaseDat(String folder, int numThreads, bool withPotentials, bool withQTF, bool x0z, bool y0z) const {
	String file = AFX(folder, "Analysis.dat");	
	
	int nNodes, nPanels;
	Body::SaveAs(dt.msh, file, Body::AQWA_DAT, Body::UNDERWATER, Bem().rho, Bem().g, y0z, x0z, nNodes, nPanels,
		dt.w, dt.head, withQTF, withPotentials, dt.h, numThreads);
}