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
	
	try {
		String line;
		LineParser f(in);	
		f.IsSeparator = IsTabSpace;
			
		line = in.GetLine();	
		f.Load(line);
	
		int nnodes;
	
		if (f.size() != 2 || !f.IsInt(0) || (nnodes = f.GetInt(0)) < 10 || !f.IsInt(1) || f.GetInt(1) < 10)
			return t_("Format error in Nemoh .dat mesh file");	// To detect Salome format
		
		mesh.Clear();
		UIndex<int> idnodes;
		
		while(true) {
			f.GetLine();
			if (f.IsEof())
				break;	

			int id = f.GetInt(0);	
			if (id > nnodes)
				throw Exc(t_("Found more nodes that the indicated inthe first row"));
			
			idnodes << (id-1);	
			Point3D &node = mesh.nodes.Add();
			node.x = f.GetDouble(1);
			node.y = f.GetDouble(2);
			node.z = f.GetDouble(3); 
			if (id == nnodes)
				break;
		}
		while(!f.IsEof()) {
			f.GetLine();	
			
			int type = f.GetInt(1);	
			if (type == 203 || type == 204) {
				Panel &panel = mesh.panels.Add();
				panel.id[0] = idnodes.Find(f.GetInt(2)-1);	
				panel.id[1] = idnodes.Find(f.GetInt(3)-1);	
				panel.id[2] = idnodes.Find(f.GetInt(4)-1);	
				if (type == 204) 
					panel.id[3] = f.GetInt(5)-1;		
				else
					panel.id[3] = panel.id[0];
			}
		}	
	} catch (Exc e) {
		return t_("Parsing error: ") + e;
	}
	
	return String();
}
