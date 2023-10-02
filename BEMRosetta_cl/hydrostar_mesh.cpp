// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2023, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"


String HydrostarMesh::LoadHst(UArray<Mesh> &mesh, String fileName, bool &y0z, bool &x0z) {
	Mesh &msh = mesh.Add();
	String ret = static_cast<HydrostarMesh&>(msh).LoadHst0(fileName, y0z, x0z);
	
	if (!ret.IsEmpty() && !ret.StartsWith(t_("Parsing error: "))) {
		mesh.Clear();
		return ret;
	}
	
	return ret;
}
	
String HydrostarMesh::LoadHst0(String fileName, bool &y0z, bool &x0z) {
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		return Format(t_("Impossible to open file '%s'"), fileName);
	
	this->fileName = fileName;
	SetCode(Mesh::HYDROSTAR_HST);
	x0z = y0z = false;
	
	try {
		String line;
		LineParser f(in);	
		f.IsSeparator = IsTabSpace;
				
		mesh.Clear();
		UIndex<int> idnodes;
		
		while(true) {
			f.GetLine_discard_empty();
			if (f.IsEof())
				break;	
			
			if (f.GetText(0) == "SYMMETRY_BODY") {
				int id1 = f.GetInt(1), 
					id2 = f.GetInt(2); 
				if (id1 == 1 && id2 == 1) {
					x0z = true;
					y0z = false;
				} else if (id1 == 1 && id2 == 2) {
					x0z = true;
					y0z = true;
				} 
			} else if (f.GetText(0) == "COORDINATES") {
				while(!f.IsEof()) {
					f.GetLine_discard_empty();
					if (f.GetText(0) == "ENDCOORDINATES")
						break;
			
					idnodes << f.GetInt(0);	
					Point3D &node = mesh.nodes.Add();
					node.x = f.GetDouble(1);
					node.y = f.GetDouble(2);
					node.z = f.GetDouble(3); 
				}
			} else if (f.GetText(0) == "PANEL") {
				while(!f.IsEof()) {
					f.GetLine_discard_empty();
					if (f.GetText(0) == "ENDPANEL")
						break;
					
					Panel &panel = mesh.panels.Add();
					panel.id[0] = idnodes.Find(f.GetInt(1));	
					panel.id[1] = idnodes.Find(f.GetInt(2));	
					panel.id[2] = idnodes.Find(f.GetInt(3));	
					if (f.size() >= 5)
						panel.id[3] = idnodes.Find(f.GetInt(4));	
					else
						panel.id[3] = panel.id[2];	
				}
			}
		}	
	} catch (Exc e) {
		return t_("Parsing error: ") + e;
	}
	
	return String();
}
