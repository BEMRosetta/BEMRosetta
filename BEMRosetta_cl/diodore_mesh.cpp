// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"


String DiodoreMesh::LoadDat(UArray<Mesh> &_mesh, String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		return t_("Impossible to open file");
	
	String line;
	LineParser f(in);	
	f.IsSeparator = IsTabSpace;

	try {		
		FileInLine::Pos fpos = in.GetPos();
		
		f.GetLine();		// To check if format is right
		
		if (f.size() == 1) {
			if (f.GetText(0) != "*NODE")
				return t_("Format error in Diodore .txt mesh file");	
		} else {
			f.GetInt(0);
			f.GetDouble(1);
			f.GetDouble(2);
			f.GetDouble(3);
		}
		
		in.SeekPos(fpos);
	} catch (...) {
		return t_("Format error");
	}

	UIndex<int> idnodes;
	try {		
		Mesh &msh = _mesh.Add();
		msh.dt.fileName = fileName;
		msh.dt.SetCode(Mesh::DIODORE_DAT);

		while(true) {
			line = in.GetLine();
			if (in.IsEof())
				break;	
			f.Load(line);
			if (f.GetText(0) == "*RETURN" || f.GetText(0) == "*TRIANGLE" || f.GetText(0) == "*QUADRANGLE")
				break;
			else if (f.GetText(0) == "*NODE")
				;
			else {
				idnodes << (f.GetInt(0)-1);		// Node ids may jump
				Point3D &node = msh.dt.mesh.nodes.Add();
				node.x = f.GetDouble(1);
				node.y = f.GetDouble(2);
				node.z = f.GetDouble(3); 
			}
		}
		for (int i = 0; i < 2; ++i) {	// Twice: Quads and triangles
			while(true) {
				line = in.GetLine();	
				if (in.IsEof())
					break;
				f.Load(line);
				if (f.GetText(0) == "*RETURN")
					break;
				else if (f.GetText(0) == "*TRIANGLE" || f.GetText(0) == "*QUADRANGLE") 
					;
				else {
					Panel &panel = msh.dt.mesh.panels.Add();
					panel.id[0] = idnodes.Find(f.GetInt(1)-1);		// Reassign the ids to be consecutive
					panel.id[1] = idnodes.Find(f.GetInt(2)-1);	
					panel.id[2] = idnodes.Find(f.GetInt(3)-1);
					if (f.GetCount() >= 5)	
						panel.id[3] = idnodes.Find(f.GetInt(4)-1);	
					else
						panel.id[3] = panel.id[0];
				}
			}	
		}
	} catch (Exc e) {
		return t_("Parsing error: ") + e;
	}
	
	return String();
}

void DiodoreMesh::SaveDat(String fileName, const Surface &surf) {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(t_("Impossible to open file"));	
	
	const UVector<Panel> &panels = surf.panels;
	const UVector<Point3D> &nodes = surf.nodes;
	
	for (int i = 0; i < nodes.size(); ++i) {
		const Point3D &node = nodes[i];
		out << Format("  %8d    % 014.7E   %014.7E   % 014.7E\n", i+1, node.x, node.y, node.z);
	}
	out << "*RETURN\n";
	
	int panelId = 1;
	for (const Panel &panel : panels) {		// First Quads
		if (panel.IsTriangle())
			continue;
		out << Format("  %8d", panelId++);
		for (int i = 0; i < 4; ++i)
			out << Format("   %8d", panel.id[i]+1);
		out << "\n";
	}
	out << "*RETURN\n";
	for (const Panel &panel : panels) {		// Second triangles
		if (!panel.IsTriangle())
			continue;
		out << Format("  %8d", panelId++);
		for (int i = 0; i < 3; ++i)
			out << Format("   %8d", panel.id[i]+1);
		out << "\n";
	}
	out << "*RETURN";	
}