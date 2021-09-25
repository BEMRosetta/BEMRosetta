// SPDX-License-Identifier: GPL-3.0-or-later
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"


String MeshData::LoadDatNemoh(String fileName, bool &x0z) {
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		return Format(t_("Impossible to open file '%s'"), fileName);
	
	this->fileName = fileName;
	SetCode(MeshData::NEMOH_DAT);
	
	try {
		String line;
		FieldSplit f(in);	
		f.IsSeparator = IsTabSpace;
			
		line = in.GetLine();	
		f.Load(line);
	
		if (!f.IsInt(0) || f.GetInt(0) != 2)
			return t_("Format error in Nemoh .dat mesh file");	// To detect Nemoh format
		
		if (f.GetInt(1) == 1)
			x0z = true;
		
		mesh.Clear();
		
		while(true) {
			line = in.GetLine();
			if (in.IsEof())
				break;	
			f.Load(line);
			int id = f.GetInt(0);	
			if (id == 0)
				break;
			Point3D &node = mesh.nodes.Add();
			node.x = f.GetDouble(1);
			node.y = f.GetDouble(2);
			node.z = f.GetDouble(3); 
		}
		while(true) {
			line = in.GetLine();	
			if (in.IsEof())
				break;
			f.Load(line);
			int id0 = f.GetInt(0);	
			if (id0 == 0)
				break;
			Panel &panel = mesh.panels.Add();
			panel.id[0] = id0-1;
			panel.id[1] = f.GetInt(1)-1;	
			panel.id[2] = f.GetInt(2)-1;	
			panel.id[3] = f.GetInt(3)-1;	
		}	
	} catch (Exc e) {
		return t_("Parsing error: ") + e;
	}
	
	return String();
}

void MeshData::SaveDatNemoh(String fileName, const Surface &surf, bool x0z) {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'\n"), fileName));	
	
	const Vector<Panel> &panels = surf.panels;
	const Vector<Point3D> &nodes = surf.nodes;
	
	out << "    2   " << (x0z ? "1" : "0") << "\n";
	for (int i = 0; i < nodes.size(); ++i) {
		const Point3D &node = nodes[i];
		out << Format("  %8d    % 014.7E   %0 14.7E   % 014.7E\n", i+1, node.x, node.y, node.z);
	}
	out << Format("  %8d    % 014.7E   %0 14.7E   % 014.7E\n", 0, 0, 0, 0);
	
	for (int i = 0; i < panels.size(); ++i) {
		const Panel &panel = panels[i];
		out << Format("  %8d   %8d   %8d   %8d\n", panel.id[0]+1, panel.id[1]+1, panel.id[2]+1, panel.id[3]+1);
	}
	out << Format("  %8d   %8d   %8d   %8d\n", 0, 0, 0, 0);
}

void MeshData::SavePreMeshNemoh(String fileName, const Surface &surf) {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'\n"), fileName));	
	
	const Vector<Panel> &panels = surf.panels;
	const Vector<Point3D> &nodes = surf.nodes;
	
	out << nodes.size() << "\n";
	out << panels.size() << "\n";
	
	for (int i = 0; i < nodes.size(); ++i) {
		const Point3D &node = nodes[i];
		out << Format("%014.7E   %014.7E   % 014.7E\n", node.x, node.y, node.z);
	}
	
	for (int i = 0; i < panels.size(); ++i) {
		const Panel &panel = panels[i];
		out << Format("%8d   %8d   %8d   %8d\n", panel.id[0]+1, panel.id[1]+1, panel.id[2]+1, panel.id[3]+1);
	}
}

