// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"


String AQWAMesh::Load_LIS(UArray<Mesh> &msh, String fileName, double g, bool &y0z, bool &x0z) {
	y0z = x0z = false;
	
	BEM bem;
	bem.g = 9.8;	// This value is necessary, but discarded
	
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

String AQWAMesh::LoadDat(UArray<Mesh> &mesh, String fileName, bool &y0z, bool &x0z) {
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		return t_("Impossible to open file");
	
	//double factorMass = 1;
	double factorLength = 1;
	
	y0z = x0z = false;
	
	UArray<Upp::Index<int>> ids;
	try {
		String line;
		line = in.GetLine();	
		
		if (!line.StartsWith("*********1*********2*********3"))
			return t_("Format error in AQWA .dat mesh file");	// To detect AQWA format

		LineParser f(in);
		
		int deck = -1;
		while(!f.IsEof()) {
			line = in.GetLine();
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
						mesh[ib].dt.SetCode(Mesh::AQWA_DAT);
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
							x0z = true;
						if ((ps = line.FindAfter("SYMY")) >= 0) 
							y0z = true;
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
			} else if (deck == 3) 
				break;
		}
		
		// Removes mooring, Morison and other points unrelated with panels
		for (Mesh &m : mesh) 
			Surface::RemoveDuplicatedPointsAndRenumber(m.dt.mesh.panels, m.dt.mesh.nodes);
		
	} catch (Exc e) {
		return t_("Parsing error: ") + e;
	}
	
	return String();
}

void AQWAMesh::SaveDat(String fileName, const UArray<Mesh> &meshes, const UArray<Surface> &surf, double rho, double g, bool y0z, bool x0z) {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'\n"), fileName));	
	
	SaveDat(out, meshes, surf, rho, g, y0z, x0z);
}

void AQWAMesh::SaveDat(Stream &ret, const UArray<Mesh> &meshes, const UArray<Surface> &surfs, double rho, double g, bool y0z, bool x0z) {
	ASSERT(meshes.size() == surfs.size());
	
	Time t = GetSysTime();
	String st = Format("%02d/%02d/%4d %02d:%02d:%02d", t.day, t.month, t.year, t.hour, t.minute, t.second);
	
	bool getQTF = false;
	double h = 300;
	double factorMass = 1;
	double factorLength = 1;
	int numCores = 4;
	
	VectorXd hrtz = VectorXd::LinSpaced(20, 0.1, 2.5)/2/M_PI;
	VectorXd head = VectorXd::LinSpaced(8, -135, 180);
	
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
	<< "NUM_CORES         " << numCores << "\n"
	<< Format("OPTIONS %s%s%s", getQTF ? "AQTF " : "", "GOON ", getQTF ? "CQTF " : "") << "\n"
	<< "OPTIONS AHD1 " << "\n"
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
						FDS(p.x/factorLength, 10, true), FDS(p.y/factorLength, 10, true), FDS(p.z/factorLength, 10, true));
		}
		const Point3D &cg = meshes[ib].dt.cg;
		ret << Format("%6d%5d         %s%s%s\n", ib+1, 98000+ib, 
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
		
	for (int ib = 0; ib < meshes.size(); ++ib)
		ret << Format("    %2d         %5d  %s", ib+1, 9800+ib, FDS(meshes[ib].GetMass(), 8, true));

	ret	<< "\n END\n";	
		
	ret
	<< "********************************************************************************\n"
	<< "*********************************** DECK  4 ************************************\n"
	<< "********************************************************************************\n"
	<< "          GEOM\n";

	for (int ib = 0; ib < meshes.size(); ++ib)
		ret << Format("%6d%s     %5d%s%s%s%s%s%s", ib+1, "PMAS", 9800+ib, 
			FDS(meshes[ib].dt.M(3, 3), 10, true), FDS(meshes[ib].dt.M(3, 4), 10, true), FDS(meshes[ib].dt.M(3, 5), 10, true),
			FDS(meshes[ib].dt.M(4, 4), 10, true), FDS(meshes[ib].dt.M(4, 5), 10, true), FDS(meshes[ib].dt.M(5, 5), 10, true));

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
	
	for (int ib = 0; ib < meshes.size(); ++ib) {
		ret << "********************************************************************************\n";
		ret << Format("          FDR%d\n", ib+1);
		for (int iw = 0; iw < hrtz.size(); ++iw) 
			ret << Format("    %2d%s   %2d   %2d %s\n", ib+1, "HRTZ", iw+1, iw+1, FDS(hrtz[iw], 9, true));
		for (int ih = 0; ih < head.size(); ++ih) 
			ret << Format("    %2d%s   %2d   %2d %s\n", ib+1, "DIRN", ih+1, ih+1, FDS(head[ih], 9, true));	
		ret	<< " END\n";
		ret << "********************************************************************************\n";
	}
}
