// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2023, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"


String HydrostarBody::LoadHst(UArray<Body> &mesh, String fileName, bool &y0z, bool &x0z) {
	Body &msh = mesh.Add(),
		 &damp = mesh.Add();
	String ret = LoadHst0(fileName, msh, damp, y0z, x0z);
	
	if (!ret.IsEmpty() && !ret.StartsWith(t_("Parsing error: "))) {
		mesh.Clear();
		return ret;
	}
	if (damp.dt.mesh.GetNumPanels() == 0) {
		mesh.SetCount(1);
	}
	return ret;
}
	
String HydrostarBody::LoadHst0(String fileName, Body &mesh, Body &damp, bool &y0z, bool &x0z) {
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		return t_(Format("Impossible to open '%s'", fileName));
	
	mesh.dt.fileName = damp.dt.fileName = fileName;
	mesh.dt.SetCode(Body::HYDROSTAR_HST);
	damp.dt.SetCode(Body::HYDROSTAR_HST);
	
	x0z = y0z = false;
	
	try {
		String line;
		LineParser f(in);	
		f.IsSeparator = IsTabSpace;
				
		mesh.dt.mesh.Clear();
		UIndex<int> idnodes;
		
		while(true) {
			f.GetLine_discard_empty();
			if (f.IsEof())
				break;	
			
			if (f.GetText(0) == "ZONEDAMPING" && f.size() == 8) {	// Only rectangular damping panels
				double xmin = f.GetDouble(1),
					   xmax = f.GetDouble(2),
					   dltx = f.GetDouble(3), 
					   ymin = f.GetDouble(4),
					   ymax = f.GetDouble(5),
					   dlty = f.GetDouble(6);
				
				Body surf;
				surf.dt.mesh.AddFlatRectangle(xmax-xmin, ymax-ymin, dltx, dlty); 
				surf.dt.mesh.Translate(xmin, ymin, 0);
				
				damp.Append(surf.dt.mesh, Null, Null);
				
			} else if (f.GetText(0) == "SYMMETRY_BODY") {
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
					Point3D &node = mesh.dt.mesh.nodes.Add();
					node.x = f.GetDouble(1);
					node.y = f.GetDouble(2);
					node.z = f.GetDouble(3); 
				}
			} else if (f.GetText(0) == "PANEL") {
				int panelType = f.GetInt(2);
				int ibegin = panelType == 1 ? 1 : 0;
				
				while(!f.IsEof()) {
					f.GetLine_discard_empty();
					if (f.GetText(0) == "ENDPANEL")
						break;
					
					Panel &panel = mesh.dt.mesh.panels.Add();
					panel.id[0] = idnodes.Find(f.GetInt(ibegin));	
					panel.id[1] = idnodes.Find(f.GetInt(ibegin+1));	
					panel.id[2] = idnodes.Find(f.GetInt(ibegin+2));	
					if (f.size() >= 4+ibegin)
						panel.id[3] = idnodes.Find(f.GetInt(ibegin+3));	
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

void HydrostarBody::SaveHST(String fileName, const Surface &surf, bool y0z, bool x0z) {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'"), fileName));	
	
	out << "NBBODY 1\n";
	out << "SYMMETRY_BODY " << (x0z ? 1 : 0) << " " << (y0z ? 2 : 1) << "\n";
	out << " COORDINATES TYPE 0\n";
	for (int i = 0; i < surf.nodes.size(); ++i)
		out << Format("   %d   %f   %f   %f\n", i+1, surf.nodes[i].x, surf.nodes[i].y, surf.nodes[i].z);
	out << " ENDCOORDINATES\n";
	out << " PANEL TYPE 1\n";
	for (int i = 0; i < surf.panels.size(); ++i)
		out << Format("   %d   %d   %d   %d   %d\n", i+1, surf.panels[i].id[0]+1, surf.panels[i].id[1]+1, 
														  surf.panels[i].id[2]+1, surf.panels[i].id[3]+1);
	out << " ENDPANEL\n";	
	out << " ENDFILE";	
}