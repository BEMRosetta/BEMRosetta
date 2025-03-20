// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"


String OffBody::LoadOff(UArray<Body> &_mesh, String fileName) {
	int a, b, c, d;
	float x, y, z;
	
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		return t_(Format("Impossible to open '%s'", fileName));
	
	String line;
	LineParser f(in);	
	f.IsSeparator = IsTabSpace;
	
	while(!in.IsEof()) {
		line = in.GetLine();
		if (line.StartsWith("OFF"))
			break;
	}
	f.GetLine();
	int np = f.GetInt(0);
	int nf = f.GetInt(1);
	int num = f.GetInt(2);
	
	try {		
		Body &msh = _mesh.Add();
		msh.dt.fileName = fileName;
		msh.dt.SetCode(Body::GEOMVIEW_OFF);
		
		msh.dt.mesh.nodes.SetCount(np);
		msh.dt.mesh.panels.SetCount(nf);
		
		for (int n = 0; n < np; n++) {
			f.GetLine();
			
			Point3D &p = msh.dt.mesh.nodes[n];
			p.x = f.GetDouble(0);
			p.y = f.GetDouble(1);
			p.z = f.GetDouble(2);
		}
		for (int n = 0; n < nf; n++) {
			f.GetLine();
			int nm = f.GetInt(0);
			
			Panel &pan = msh.dt.mesh.panels[n];
			pan.id[0] = f.GetInt(1);
			pan.id[1] = f.GetInt(2);
			pan.id[2] = f.GetInt(3);
			if (nm == 4)
				pan.id[3] = f.GetInt(4);
			else
				pan.id[3] = pan.id[0];
			for (int i = 0; i < 4; ++i) {
				int id = pan.id[i];
				if (id < 0 || id >= np)
					throw Exc(in.Str() + "\n"  + t_("Wrong vertex identifier"));
			}
		}		
	} catch (Exc e) {
		return t_("Parsing error: ") + e;
	}
	
	return String();
}

void OffBody::SaveOff(String fileName, const Surface &surf) {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(t_(Format("Impossible to open '%s'", fileName)));

	const UVector<Panel> &panels = surf.panels;
	const UVector<Point3D> &nodes = surf.nodes;
		
	out << "OFF\n";
	
	out << Format("%d %d %d\n", nodes.size(), panels.size(), 0);

	for (int n = 0; n < nodes.size(); n++) 
		out << Format("%f %f %f\n", nodes[n].x, nodes[n].y, nodes[n].z);

	for (int n = 0; n < panels.size(); n++) {
		if (panels[n].IsTriangle())
			out << Format("3 %d %d %d\n", panels[n].id[0], panels[n].id[1], panels[n].id[2]);
		else
			out << Format("4 %d %d %d %d\n", panels[n].id[0], panels[n].id[1], panels[n].id[2], panels[n].id[3]);
	}
}