// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2024, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"


String AQWABody::LoadLis(UArray<Body> &msh, String fileName, double g, bool &y0z, bool &x0z) {
	y0z = x0z = false;
	
	BEM bem;
	bem.g = g;	// This value is necessary, but discarded
	
	try {
		bem.LoadBEM(fileName, Null, false);
		
		Hydro &hy = bem.hydros[0];

		msh = pick(hy.dt.msh);
		y0z = hy.dt.symX;
		x0z = hy.dt.symY;
	} catch (Exc e) {
		return t_("Parsing error: ") + e;
	}
		
	return String();
}

String AQWABody::LoadDatANSYSTOAQWA(UArray<Body> &mesh, Hydro &hy, String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		return t_(Format("Impossible to open '%s'", fileName));
		
	hy.dt.symX = hy.dt.symY = false;
	
	hy.dt.w.Clear();
	hy.dt.head.Clear();
	mesh.Clear();

	UArray<Upp::Index<int>> ids;
	try {
		String line;
		line = in.GetLine();
		
		if (!line.StartsWith("******************************"))
			return t_("Format error in ANSYSTOAQWA .dat mesh file");	// To detect ANSYSTOAQWA format
		
		bool ansysFound = false;
		for (int i = 0; i < 3; ++i) {
			line = in.GetLine();
			if (line.Find("ANSYS") >= 0) {
				ansysFound = true;
				break;
			}
		}
		if (!ansysFound)
			return t_("Format error in ANSYSTOAQWA .dat mesh file");	// To detect ANSYSTOAQWA format
		
		LineParser f(in);
		
		int ib = -1;
		int deck = -1;
		while(!f.IsEof()) {
			line = f.GetLine();
			if (line.GetCount() >= 6)
				deck = ScanInt(line.Mid(4, 2));
			
			if (deck == 1) {
				f.LoadFields(line, {0, 6, 12, 20, 30, 40});
				if (line.Find("STRC") >= 0) {
					ib = f.GetInt(2)-1;
					mesh.SetCount(ib+1);
					ids.SetCount(ib+1);
					mesh[ib].dt.fileName = fileName;
					mesh[ib].dt.SetCode(Body::AQWA_DAT);
				} else if (f.GetCount() == 6) {
					if (ib == -1) {	// STRC not found, one body only
						ib = 0;
						mesh.SetCount(ib+1);
						ids.SetCount(ib+1);
						mesh[ib].dt.fileName = fileName;
						mesh[ib].dt.SetCode(Body::AQWA_DAT);
					}
					int id = f.GetInt(1);
					if (id == 99999) {
						mesh[ib].dt.cg.x = f.GetDouble(3);
						mesh[ib].dt.cg.y = f.GetDouble(4);
						mesh[ib].dt.cg.z = f.GetDouble(5);
						mesh[ib].dt.c0 = clone(mesh[ib].dt.cg);		// In AQWA, cg == c0
					} else {
						ids[ib] << id;
						Point3D &node = mesh[ib].dt.mesh.nodes.Add();
						node.x = f.GetDouble(3);
						node.y = f.GetDouble(4);
						node.z = f.GetDouble(5);
					}
				}
			} else if (deck == 2) {
				f.LoadFields(line, {4, 6, 11, 15, 25, 32, 39, 46, 53});
				
				if (line.Find("ELM") >= 0) 
					ib = ScanInt(line.Mid(13, 2))-1;
				else if (f.GetCount() > 6) {	
					String code = f.GetText(1);
					if (code == "QPPL" || code == "TPPL") {
						Panel &panel = mesh[ib].dt.mesh.panels.Add();
						int id;
						id = f.GetInt(4);
						if ((id = ids[ib].Find(id)) < 0)
							throw Exc(in.Str() + "\n"  + t_("id 1 not found"));
						panel.id[0] = id;
						id = f.GetInt(5);
						if ((id = ids[ib].Find(id)) < 0)
							throw Exc(in.Str() + "\n"  + t_("id 2 not found"));
						panel.id[1] = id;
						id = f.GetInt(6);
						if ((id = ids[ib].Find(id)) < 0)
							throw Exc(in.Str() + "\n"  + t_("id 3 not found"));
						panel.id[2] = id;
						
						if (code == "QPPL") {
							id =  f.GetInt(7);
							if ((id = ids[ib].Find(id)) < 0)
								throw Exc(in.Str() + "\n"  + t_("id 4 not found"));
							panel.id[3] = id;
						} else if (code == "TPPL") 
							panel.id[3] = panel.id[0];
					}
				}
			} else if (deck == 3) {
				if (f.GetCount() == 3) {
					int ib = f.GetInt(1)-1;
					if (ib >= mesh.size())
						throw Exc(in.Str() + "\n"  + Format(t_("Body %d not found"), ib+1));
					if (mesh[ib].dt.M.size() != 36)
						mesh[ib].dt.M = MatrixXd::Zero(6, 6);
					mesh[ib].dt.M(0, 0) = mesh[ib].dt.M(1, 1) = mesh[ib].dt.M(2, 2) = f.GetDouble(2);
				}	
			} else if (deck == 4) {
				String str = f.GetText(0);
				if (str.Find("PMAS") > 0) {
					int ib = f.GetInt(1)-1;
					if (ib >= mesh.size())
						throw Exc(in.Str() + "\n"  + Format(t_("Body %d not found"), ib+1));
					if (mesh[ib].dt.M.size() != 36)
						mesh[ib].dt.M = MatrixXd::Zero(6, 6);
					mesh[ib].dt.M(3, 3) = f.GetDouble(2);
					mesh[ib].dt.M(3, 4) = mesh[ib].dt.M(4, 3) = f.GetDouble(3);
					mesh[ib].dt.M(3, 5) = mesh[ib].dt.M(5, 3) = f.GetDouble(4);
					mesh[ib].dt.M(4, 4) = f.GetDouble(5);
					mesh[ib].dt.M(4, 5) = mesh[ib].dt.M(5, 4) = f.GetDouble(6);
					mesh[ib].dt.M(5, 5) = f.GetDouble(7);
				}
			} else if (deck == 5) {
				f.LoadFields(line, {6, 10, 21});
				if (f.GetCount() == 2) {
					if (f.GetText(0) == "DPTH")
						hy.dt.h = f.GetDouble(1);
					else if (f.GetText(0) == "DENS")
						hy.dt.rho = f.GetDouble(1);
					else if (f.GetText(0) == "ACCG")
						hy.dt.g = f.GetDouble(1);
				}
			} else if (deck == 6) {
				if (f.GetText(1) == "FDR1") 
					ib = 0;
				else if (f.GetText(1).StartsWith("FDR"))
					ib = -1;	
			    else if (ib == 0) {
				    if (f.GetText(0).Find("FREQ") >= 0) {
				        if (f.GetCount() == 4)
				        	LinSpaced(hy.dt.w, f.GetInt(1), f.GetDouble(2), f.GetDouble(3)); 
				        else {
				            for (int i = 3; i < f.size(); ++i)
				            	 hy.dt.w << f.GetDouble(i);
				        }
				    } else if (f.GetText(0).Find("PERD") >= 0) {
				        if (f.GetCount() == 4) {
				            UVector<double> T;
				        	LinSpaced(T, f.GetInt(1), f.GetDouble(2), f.GetDouble(3));
				        	hy.dt.w.SetCount(T.size());
				        	for (int i = 0; i < T.size(); ++i)
				        		hy.dt.w[i] = 2*M_PI/T[i];
				        }
				    } else if (f.GetText(0).Find("DIRN") >= 0) {
				        if (f.GetCount() == 4 && (f.GetDouble(1) != f.GetDouble(2)))
				        	LinSpaced(hy.dt.head, f.GetInt(1), f.GetDouble(2), f.GetDouble(3)); 
				        else {
				            for (int i = 3; i < f.size(); ++i)
				            	 hy.dt.head << f.GetDouble(i);
				        }
				    }
			    }
			} 
		}
		
		hy.dt.Nf = hy.dt.w.size();
		hy.dt.Nh = hy.dt.head.size();
		hy.dt.Nb = mesh.size();
		
		// Removes mooring, Morison and other points unrelated with panels
		for (Body &m : mesh) 
			Surface::RemoveDuplicatedPointsAndRenumber(m.dt.mesh.panels, m.dt.mesh.nodes);
			
	} catch (Exc e) {
		return t_("Parsing error: ") + e;
	}
	
	return String();	
}

String AQWABody::LoadDat(UArray<Body> &mesh, Hydro &hy, String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		return t_(Format("Impossible to open '%s'", fileName));
	
	//double factorMass = 1;
	double factorLength = 1;
	
	hy.dt.symX = hy.dt.symY = false;
	
	hy.dt.w.Clear();
	hy.dt.head.Clear();
	mesh.Clear();
	
	UArray<Upp::Index<int>> ids;
	try {
		String line;
		line = in.GetLine();	
		
		if (!line.StartsWith("*********1*********2*********3"))
			return t_("Format error in AQWA .dat mesh file");	// To detect AQWA format

		LineParser f(in);
		
		int deck = -1;
		while(!f.IsEof()) {
			line = f.GetLine();
			int pos;
			if (line.StartsWith("*")) {
				if ((pos = line.FindAfter("DECK")) > 0)
					deck = ScanInt(line.Mid(pos));
				else if ((pos = line.FindAfter("Unit System :")) >= 0) {
					String system = Trim(line.Mid(pos));
					if (system.Find("Metric") < 0)
						throw Exc(in.Str() + "\n" + t_("Only metric system is supported"));
					/*if (system.Find("kg") > 0)
						factorMass = 1;
					else if (system.Find("tonne") > 0)
						factorMass = 1000;
					else 
						throw Exc(in.Str() + "\n" + t_("Unknown mass unit"));*/
					if (system.Find("m ") > 0)
						factorLength = 1;
					else if (system.Find("km ") > 0)
						factorLength = 1000;
					else 
						throw Exc(in.Str() + "\n" + t_("Unknown length unit"));
				}
			}
			
			if (deck == 1) {
				f.LoadFields(line, {4, 6, 20, 30, 40});
				
				int ib = f.GetInt_nothrow(0);
				if (!IsNull(ib)) {
					ib--;
					if (ib+1 > mesh.size()) {
						mesh.SetCount(ib+1);
						ids.SetCount(ib+1);
						mesh[ib].dt.fileName = fileName;
						mesh[ib].dt.SetCode(Body::AQWA_DAT);
					}
					int id = f.GetInt(1);
					if (id >= 98000 && id < 99000) {
						mesh[ib].dt.cg.x = f.GetDouble(2);
						mesh[ib].dt.cg.y = f.GetDouble(3);
						mesh[ib].dt.cg.z = f.GetDouble(4);
						mesh[ib].dt.c0 = clone(mesh[ib].dt.cg);		// In AQWA, cg == c0
					} else {
						ids[ib] << id;
						Point3D &node = mesh[ib].dt.mesh.nodes.Add();
						node.x = f.GetDouble(2)*factorLength;
						node.y = f.GetDouble(3)*factorLength;
						node.z = f.GetDouble(4)*factorLength;
					}
				}
			} else if (deck == 2) {
				f.LoadFields(line, {4, 6, 11, 15, 24, 31, 38, 45, 52});
				
				int ib = f.GetInt_nothrow(0);
				if (IsNull(ib)) {
					int ps;
					if ((ps = line.FindAfter("ELM")) >= 0) {
						ib = ScanInt(line.Mid(ps, 3))-1;
						String name = Trim(line.Mid(ps+3));
						mesh[ib].dt.name = name;
					} else {
						if ((ps = line.FindAfter("SYMX")) >= 0) 
							hy.dt.symY = true;
						if ((ps = line.FindAfter("SYMY")) >= 0) 
							hy.dt.symX = true;
					}
				} else {
					ib--;
					String code = f.GetText(1);
					if (code == "QPPL" || code == "TPPL") {
						Panel &panel = mesh[ib].dt.mesh.panels.Add();
						int id;
						id = f.GetInt(4);
						if ((id = ids[ib].Find(id)) < 0)
							throw Exc(in.Str() + "\n"  + t_("id 1 not found"));
						panel.id[0] = id;
						id = f.GetInt(5);
						if ((id = ids[ib].Find(id)) < 0)
							throw Exc(in.Str() + "\n"  + t_("id 2 not found"));
						panel.id[1] = id;
						id = f.GetInt(6);
						if ((id = ids[ib].Find(id)) < 0)
							throw Exc(in.Str() + "\n"  + t_("id 3 not found"));
						panel.id[2] = id;
						
						if (code == "QPPL") {
							id =  f.GetInt(7);
							if ((id = ids[ib].Find(id)) < 0)
								throw Exc(in.Str() + "\n"  + t_("id 4 not found"));
							panel.id[3] = id;
						} else if (code == "TPPL") 
							panel.id[3] = panel.id[0];
					}
				}
			} else if (deck == 3) {
				if (f.IsInt(0)) {
					int ib = f.GetInt(0)-1;
					if (ib >= mesh.size())
						throw Exc(in.Str() + "\n"  + Format(t_("Body %d not found"), ib+1));
					if (f.GetInt(1) >= 98000) {
						if (mesh[ib].dt.M.size() != 36)
							mesh[ib].dt.M = MatrixXd::Zero(6, 6);
						mesh[ib].dt.M(0, 0) = mesh[ib].dt.M(1, 1) = mesh[ib].dt.M(2, 2) = f.GetDouble(2);
					}
				}	
			} else if (deck == 4) {
				f.LoadFields(line, {1, 15, 21, 30, 41, 51, 61, 71});
				
				String str = f.GetText(0);
				if (str.Find("PMAS") > 0) {
					str.Replace("PMAS", "");
					int ib = ScanInt(str);
					if (!IsNull(ib) && f.GetInt(1) >= 98000) {
						ib--;
						if (ib >= mesh.size())
							throw Exc(in.Str() + "\n"  + Format(t_("Body %d not found"), ib+1));
						if (mesh[ib].dt.M.size() != 36)
							mesh[ib].dt.M = MatrixXd::Zero(6, 6);
						mesh[ib].dt.M(3, 3) = f.GetDouble(2);
						mesh[ib].dt.M(3, 4) = mesh[ib].dt.M(4, 3) = f.GetDouble(3);
						mesh[ib].dt.M(3, 5) = mesh[ib].dt.M(5, 3) = f.GetDouble(4);
						mesh[ib].dt.M(4, 4) = f.GetDouble(5);
						mesh[ib].dt.M(4, 5) = mesh[ib].dt.M(5, 4) = f.GetDouble(6);
						mesh[ib].dt.M(5, 5) = f.GetDouble(7);
					}
				}
			} else if (deck == 5) {
				if (f.GetText(0) == "DPTH")
					hy.dt.h = f.GetDouble(1);
				else if (f.GetText(0) == "DENS")
					hy.dt.rho = f.GetDouble(1);
				else if (f.GetText(0) == "ACCG")
					hy.dt.g = f.GetDouble(1);
			} else if (deck == 6) {
			    if (f.GetText(0) == "1HRTZ") 
			     	hy.dt.w << 2*M_PI*f.GetDouble(3);
			    else if (f.GetText(0) == "1DIRN") 
					hy.dt.head << f.GetDouble(3);
			} 
		}
		
		hy.dt.Nf = hy.dt.w.size();
		hy.dt.Nh = hy.dt.head.size();
		hy.dt.Nb = mesh.size();
		
		// Removes mooring, Morison and other points unrelated with panels
		for (Body &m : mesh) 
			Surface::RemoveDuplicatedPointsAndRenumber(m.dt.mesh.panels, m.dt.mesh.nodes);
		
	} catch (Exc e) {
		return t_("Parsing error: ") + e;
	}
	
	return String();
}

void AQWABody::SaveDat(String fileName, const UArray<Body> &mesh, const UArray<Surface> &surfs, double rho, double g, bool y0z, bool x0z,
			const UVector<double> &w, const UVector<double> &head, bool getQTF, bool getPotentials, double h, int numThreads) {
	FileOut ret(fileName);
	if (!ret.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'\n"), fileName));
	
	if (IsNull(numThreads) || numThreads <= 0)
		numThreads = 8;
	
	ASSERT(mesh.size() == surfs.size());
	
	Time t = GetSysTime();
	String st = Format("%02d/%02d/%4d %02d:%02d:%02d", t.day, t.month, t.year, t.hour, t.minute, t.second);
	
	double factorMass = 1;
	double factorLength = 1;
	
	
	UVector<double> hrtz(w.size());
	for (int i = 0; i < w.size(); ++i)
		hrtz[i] = w[i]/2/M_PI;
	
	ret 
	<< "*********1*********2*********3*********4*********5*********6*********7*********8" << "\n"
	<< "*2345678901234567890123456789012345678901234567890123456789012345678901234567890" << "\n"
	<< "********************************************************************************" << "\n"
	<< "************************** File generated by BEMRosetta ************************" << "\n"
	<< "********************************************************************************" << "\n"
	<< "* Project Title                   : " << "\n"
	<< "* Project Reference               : " << "\n"
	<< "* Project Author                  : " << "\n"
	<< "* Project Description             : " << "\n"
	<< "* Date of Creation                : " << st << "\n"
	<< "* Last Modified                   : " << st << "\n"
	<< "* This analysis file written at   : " << st << "\n";
	String kg = factorMass == 1 ? "kg" : "tonne";
	String N = factorMass == 1 ? "N" : "kN";
	String m = factorLength == 1 ? "m" : "km";
	ret << Format("* Hydrodynamic Solver Unit System : Metric: %s, %s [%s]\n", kg, m, N);
	ret
	<< "********************************************************************************" << "\n"
	<< "*********************************** DECK  0 ************************************" << "\n"
	<< "********************************************************************************" << "\n"
	<< "JOB AQWA  LINE" << "\n"
	<< "TITLE               " << "\n"
	<< "NUM_CORES         " << numThreads << "\n"
	<< Format("OPTIONS %s%s%s", getQTF ? "AQTF " : "", "GOON ", getQTF ? "CQTF " : "") << "\n"
	<< "OPTIONS AHD1 " << (getPotentials ? "PRPT PRPR" : "") <<"\n"
	<< Format("OPTIONS %s%s LHFR REST END", getQTF ? "NQTF " : "", "LDRG") << "\n"
	<< "RESTART  1  5" << "\n"
	<< "********************************************************************************" << "\n"
	<< "*********************************** DECK  1 ************************************" << "\n"
	<< "********************************************************************************" << "\n"
	<< "          COOR" << "\n"
	<< "      NOD5" << "\n";

	UVector<int> firstidbody(surfs.size());
	firstidbody[0] = 0;
	for (int ib = 1; ib < surfs.size(); ++ib) 
		firstidbody[ib] = firstidbody[ib-1] + surfs[ib-1].nodes.size();	
		
	for (int ib = 0; ib < surfs.size(); ++ib) {
		ret << Format("      STRC       %2d\n", ib+1);
		
		const Surface &surf = surfs[ib];
		
		for (int in = 0; in < surf.nodes.size(); ++in) {
			const Point3D &p = surf.nodes[in];
			ret << Format("%6d%5d         %s%s%s\n", ib+1, firstidbody[ib] + in +1, 
						Replace(FDS(p.x/factorLength, 10, true), "E", "e"), 
						Replace(FDS(p.y/factorLength, 10, true), "E", "e"), 
						Replace(FDS(p.z/factorLength, 10, true), "E", "e"));
		}
		const Point3D &cg = mesh[ib].dt.cg;
		ret << Format("%6d%5d        %s%s%s\n", ib+1, 98000+ib, 
						FDS(cg.x/factorLength, 10, true), FDS(cg.y/factorLength, 10, true), FDS(cg.z/factorLength, 10, true));	
	}
	ret << " END" << "\n"
	<< "********************************************************************************" << "\n"
	<< "*********************************** DECK  2 ************************************" << "\n"
	<< "********************************************************************************" << "\n";
	
	auto IsDry = [](const Surface &surf, int ip) {
		const Panel &p = surf.panels[ip];
		for (int id = 0; id < 4; ++id) {
			int pid = p.id[id];
			if (surf.nodes[pid].z < -EPS_LEN)
				return false;
		}
		return true;
	};
	
	int ipall = 1;
	
	for (int ib = 0; ib < surfs.size(); ++ib) {
		ret << "********************************************************************************\n";
		if (ib == 0)
			ret << Format("          ELM%d      Body%d\n", ib+1, ib+1);
		if (y0z || x0z) {
			ret << "     ";
			if (y0z)
				ret << " SYMY";
			if (x0z)
				ret << " SYMX";
			ret << "\n";
		}
		ret
		<< "      SEAG          ( 81, 51,-270.24102, 339.52305,-230.62436, 230.62436)" << "\n"
		<< "      ZLWL          (        0.)" << "\n";
		
		const Surface &surf = surfs[ib];
		
		for (int ip = 0; ip < surf.panels.size(); ++ip) {
			const Panel &p = surf.panels[ip];
			
			String diff = IsDry(surf, ip) ? "    " : "DIFF";
			
			if (p.IsTriangle())
				ret << Format("%6d%s %s%5d(1)(%5d)(%5d)(%5d)  WB Elem Id:%5d Aqwa Elem No:%5d\n", ib+1, "TPPL", diff, 14+ib, 
							p.id[0]+1+firstidbody[ib], p.id[1]+1+firstidbody[ib], p.id[2]+1+firstidbody[ib], ipall, ip+1);
			else
				ret << Format("%6d%s %s%5d(1)(%5d)(%5d)(%5d)(%5d)  WB Elem Id:%5d Aqwa Elem No:%5d\n", ib+1, "QPPL", diff, 14+ib, 
							p.id[0]+1+firstidbody[ib], p.id[1]+1+firstidbody[ib], p.id[2]+1+firstidbody[ib], p.id[3]+1+firstidbody[ib], ipall, ip+1);
		     
			ipall++;
		}
		ret << Format("%6d%s     %5d(1)(%5d)(%5d)(%5d)\n", ib+1, "PMAS", 50+ib, 98000+ib, 98000+ib, 98000+ib);
		ret
		<< " END\n"
		<< "********************************************************************************\n";
	}
	ret << "          FINI\n";
	ret
	<< "********************************************************************************\n"
	<< "*********************************** DECK  3 ************************************\n"
	<< "********************************************************************************\n"
	<< "          MATE\n";
		
	for (int ib = 0; ib < mesh.size(); ++ib)
		ret << Format("    %2d         %5d  %s", ib+1, 98000+ib, FDS(mesh[ib].GetMass(), 8, true));

	ret	<< "\n END\n";	
		
	ret
	<< "********************************************************************************\n"
	<< "*********************************** DECK  4 ************************************\n"
	<< "********************************************************************************\n"
	<< "          GEOM\n";

	for (int ib = 0; ib < mesh.size(); ++ib)
		ret << Format("%6d%s     %5d %s %s %s %s %s %s", ib+1, "PMAS", 98000+ib, 
			FDS(mesh[ib].dt.M(3, 3), 9, true), FDS(mesh[ib].dt.M(3, 4), 9, true), FDS(mesh[ib].dt.M(3, 5), 9, true),
			FDS(mesh[ib].dt.M(4, 4), 9, true), FDS(mesh[ib].dt.M(4, 5), 9, true), FDS(mesh[ib].dt.M(5, 5), 9, true));

	ret	<< "\n END\n";
	
	ret
	<< "********************************************************************************\n"
	<< "*********************************** DECK  5 ************************************\n"
	<< "********************************************************************************\n"
	<< "          GLOB\n";
	
	ret 
	<<	Format("      DPTH%s\n", FDS(h/factorLength, 10, true))
	<<	Format("      DENS%s\n", FDS(rho/factorMass, 10, true))
	<<	Format("      ACCG%s\n", FDS(g/factorLength, 10, true));

	ret	<< " END\n";

	ret
	<< "********************************************************************************\n"
	<< "*********************************** DECK  6 ************************************\n"
	<< "********************************************************************************\n";
	
	UVector<double> head180;
	for (int i = 0; i < head.size(); ++i)
		FindAdd(head180, FixHeading_180(head[i]));
	Sort(head180);
	if (First(head180) != -180)
		head180.Insert(0, -180);
	if (Last(head180) != 180)
		head180.Insert(0, 180);
	
	for (int ib = 0; ib < mesh.size(); ++ib) {
		ret << "********************************************************************************\n";
		ret << Format("          FDR%d\n", ib+1);
		for (int iw = 0; iw < hrtz.size(); ++iw) 
			ret << Format("    %2d%s  %3d  %3d %s\n", ib+1, "HRTZ", iw+1, iw+1, Replace(FDS(hrtz[iw], 9, true), "E", "e"));
		for (int ih = 0; ih < head180.size(); ++ih) 
			ret << Format("    %2d%s   %2d   %2d %s\n", ib+1, "DIRN", ih+1, ih+1, FDS(head180[ih], 9, true));	
		ret	<< " END\n";
		ret << "********************************************************************************\n";
	}
	ret
	<< "********************************************************************************\n"
	<< "*********************************** DECK  7 ************************************\n"
	<< "********************************************************************************\n"
	<< "********************************************************************************\n"
	<< "          WFS1\n"
	<< " END\n"
	<< "********************************************************************************\n"
	<< "********************************************************************************\n"
	<< "*********************************** DECK  8 ************************************\n"
	<< "********************************************************************************\n"
	<< "********************************************************************************\n"
	<< "          DRC1\n"
	<< "* No data defined for this structure\n"
	<< " END\n"
	<< "********************************************************************************\n"
	<< "********************************************************************************\n"
	<< "          NONE\n"
	<< "          NONE\n"
	<< "          NONE\n"
	<< "          NONE\n"
	<< "*********************************** DECK 13 ************************************\n"
	<< "********************************************************************************\n"
	<< "          NONE\n"
	<< "          NONE\n"
	<< "          NONE\n"
	<< "          NONE\n"
	<< "          NONE\n"
	<< "          NONE";
}
