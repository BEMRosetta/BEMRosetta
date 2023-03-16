// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"


String AQWAMesh::LoadDat(UArray<Mesh> &mesh, String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		return Format(t_("Impossible to open file '%s'"), fileName);
	
	double factorMass = 1;
	double factorLength = 1;
	
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
			
			if (deck == 1) {
				f.LoadFields(line, {4, 6, 20, 30, 40});
				
				int ib = f.GetInt_nothrow(0);
				if (!IsNull(ib)) {
					ib--;
					if (ib+1 > mesh.size()) {
						mesh.SetCount(ib+1);
						ids.SetCount(ib+1);
						mesh[ib].fileName = fileName;
						mesh[ib].SetCode(Mesh::AQWA_DAT);
					}
					int id = f.GetInt(1);
					if (id >= 98000 && id < 99000) {
						mesh[ib].cg.x = f.GetDouble(2);
						mesh[ib].cg.y = f.GetDouble(3);
						mesh[ib].cg.z = f.GetDouble(4);
						mesh[ib].c0 = clone(mesh[ib].cg);		// In AQWA, cg == c0
					} else {
						ids[ib] << id;
						Point3D &node = mesh[ib].mesh.nodes.Add();
						node.x = f.GetDouble(2)*factorLength;
						node.y = f.GetDouble(3)*factorLength;
						node.z = f.GetDouble(4)*factorLength;
					}
				}
			} else if (deck == 2) {
				f.LoadFields(line, {4, 6, 11, 15, 24, 31, 38, 45, 52});
				
				int ib = f.GetInt_nothrow(0);
				if (IsNull(ib)) {
					int pos = line.FindAfter("ELM");
					if (pos > 0) {
						ib = ScanInt(line.Mid(pos, 3))-1;
						String name = Trim(line.Mid(pos+3));
						mesh[ib].name = name;
					}
				} else {
					ib--;
					String code = f.GetText(1);
					if (code == "QPPL" || code == "TPPL") {
						Panel &panel = mesh[ib].mesh.panels.Add();
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
	} catch (Exc e) {
		return t_("Parsing error: ") + e;
	}
	
	return String();
}

String AQWAMesh::SaveDat(const UArray<Surface> &surf, bool y0z, bool x0z) {
	String ret;
		
	return ret;
}

void AQWAMesh::SaveDat(String fileName, const UArray<Surface> &surf, bool y0z, bool x0z) {
	
	
	
}