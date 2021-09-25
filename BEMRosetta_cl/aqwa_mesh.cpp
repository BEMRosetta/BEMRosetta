// SPDX-License-Identifier: GPL-3.0-or-later
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"


String MeshData::LoadDatAQWA(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		return Format(t_("Impossible to open file '%s'"), fileName);
	
	this->fileName = fileName;
	SetCode(MeshData::AQWA_DAT);
	
	mesh.Clear();
	
	Index<int> ids;
	try {
		String line;
		line = in.GetLine();	
		
		if (!line.StartsWith("*********1*********2*********3"))
			return t_("Format error in AQWA .dat mesh file");	// To detect AQWA format

		bool done = false;
		while(!in.IsEof()) {
			line = in.GetLine();
			
			if (line.Find("DECK  1") >= 0) {
				done = true;
				break;
			}
		}
		if (!done)
			throw Exc(t_("Format error in AQWA .dat mesh file. DECK 1 not found"));
		
		done = false;
		while(!in.IsEof()) {
			line = in.GetLine();
			
			if (line.GetCount() > 5 && line[5] == '1') {		// First body?
				int id = ScanInt(line.Mid(6, 5));
				if (id == 98000) {
					cg.x = ScanDouble(line.Mid(20, 10));
					cg.y = ScanDouble(line.Mid(30, 10));
					cg.z = ScanDouble(line.Mid(40, 10));
				} else {
					ids << id;
					Point3D &node = mesh.nodes.Add();
					node.x = ScanDouble(line.Mid(20, 10));
					node.y = ScanDouble(line.Mid(30, 10));
					node.z = ScanDouble(line.Mid(40, 10));
				}
			} else if (line.StartsWith(" END")) {
				done = true;
				break;
			}
		}
		if (!done)
			throw Exc(t_("Format error in AQWA .dat mesh file. Points list END not found"));
		
		done = false;
		while(!in.IsEof()) {
			line = in.GetLine();
			
			if (line.GetCount() > 10 && line[5] == '1') {
				String code = line.Mid(6, 4);
				bool diff = true;//line.Mid(11, 4) == "DIFF";
				
				if (code == "QPPL" || code == "TPPL") {
					Panel &panel = mesh.panels.Add();
					int id;
					id = ScanInt(line.Mid(24, 5));
					if ((id = ids.Find(id)) < 0)
						throw Exc(in.Str() + "\n"  + t_("id 1 not found"));
					panel.id[0] = id;
					id = ScanInt(line.Mid(31, 5));
					if ((id = ids.Find(id)) < 0)
						throw Exc(in.Str() + "\n"  + t_("id 2 not found"));
					panel.id[1] = id;
					id = ScanInt(line.Mid(38, 5));
					if ((id = ids.Find(id)) < 0)
						throw Exc(in.Str() + "\n"  + t_("id 3 not found"));
					panel.id[2] = id;
					
					if (diff) {
						if (code == "QPPL") {
							id = ScanInt(line.Mid(45, 5));
							if ((id = ids.Find(id)) < 0)
								throw Exc(in.Str() + "\n"  + t_("id 4 not found"));
							panel.id[3] = id;
						} else if (code == "TPPL") 
							panel.id[3] = panel.id[0];
					}
				}
			} else if (line.StartsWith(" END")) {
				done = true;
				break;
			}
		}
		if (!done)
			throw Exc(t_("Format error in AQWA .dat mesh file. Panels list END not found"));
			
	} catch (Exc e) {
		return t_("Parsing error: ") + e;
	}
	
	return String();
}

