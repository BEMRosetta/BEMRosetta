// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"


String HAMSBody::LoadPnl(UArray<Body> &_mesh, String fileName, bool &y0z, bool &x0z) {
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		return t_(Format("Impossible to open '%s'", fileName));
	
	Body &msh = _mesh.Add();
	msh.dt.fileName = fileName;
	msh.dt.SetCode(Body::HAMS_PNL);
	
	LineParser f(in);	
	f.IsSeparator = IsTabSpace;
		
	Upp::Index<int> ids;
	try {
		String line;
		line = in.GetLine();	
		
		if (line.Find("Mesh File") < 0)
			return t_("File is not in HAMS .pnl format");	// To detect HAMS format

		bool done = false;
		while(!in.IsEof()) {
			line = in.GetLine();
			
			if (line.Find("Number of Panels") >= 0) {
				done = true;
				break;
			}
		}
		if (!done)
			throw Exc(t_("Format error in HAMS .dat mesh file. Number of Panels, Nodes, X-Symmetry and Y-Symmetry not found"));
		
		f.Load(in.GetLine());
		//int numPanels = f.GetInt(0);
		//int numNodes = f.GetInt(1);
		y0z = f.GetInt(2) == 1;
		x0z = f.GetInt(3) == 1;
		
		done = false;
		while(!in.IsEof()) {
			line = in.GetLine();
			
			if (line.Find("Start Definition of Node Coordinates") >= 0) {
				done = true;
				break;
			}
		}
		if (!done)
			throw Exc(t_("Format error in HAMS .dat mesh file. Start Definition of Node Coordinates not found"));
		
		done = false;
		while(!in.IsEof()) {
			line = in.GetLine();
			
			if (line.Find("End") >= 0) {
				done = true;
				break;
			} else {
				f.Load(line);
				ids << f.GetInt(0);
				Point3D &node = msh.dt.mesh.nodes.Add();
				node.x = f.GetDouble(1);
				node.y = f.GetDouble(2);
				node.z = f.GetDouble(3);
			} 
		}
		if (!done)
			throw Exc(t_("Format error in HAMS .dat mesh file. Points list End not found"));
		
		done = false;
		while(!in.IsEof()) {
			line = in.GetLine();
			
			if (line.Find("Start Definition of Node Relations") >= 0) {
				done = true;
				break;
			}
		}
		if (!done)
			throw Exc(t_("Format error in HAMS .dat mesh file. Start Definition of Node Relations not found"));		
		
		done = false;
		while(!in.IsEof()) {
			line = in.GetLine();
			
			 if (line.Find("End") > 0) {
				done = true;
				break;
			} else {
				f.Load(line);
				
				int numVert = f.GetInt(1);

				Panel &panel = msh.dt.mesh.panels.Add();
				int id;
				id = f.GetInt(2);
				if ((id = ids.Find(id)) < 0)
					throw Exc(in.Str() + "\n"  + t_("id 1 not found"));
				panel.id[0] = id;
				id = f.GetInt(3);
				if ((id = ids.Find(id)) < 0)
					throw Exc(in.Str() + "\n"  + t_("id 2 not found"));
				panel.id[1] = id;
				id = f.GetInt(4);
				if ((id = ids.Find(id)) < 0)
					throw Exc(in.Str() + "\n"  + t_("id 3 not found"));
				panel.id[2] = id;
					
				if (numVert == 4) {
					id = f.GetInt(5);
					if ((id = ids.Find(id)) < 0)
						throw Exc(in.Str() + "\n"  + t_("id 4 not found"));
					panel.id[3] = id;
				} else if (numVert == 3) 
					panel.id[3] = panel.id[0];
			}
		}
		if (!done)
			throw Exc(t_("Format error in HAMS .pnl mesh file. Panels list End not found"));
			
	} catch (Exc e) {
		return t_("Parsing error: ") + e;
	}
	
	return String();
}

void HAMSBody::SavePnl(String fileName, const Surface &surf, bool y0z, bool x0z) {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'"), fileName));	
	
	char type = surf.IsWaterPlaneMesh();
	if (type == 'x')
		throw Exc(t_("Mesh contains both waterplane and hull panels. Please separate them before saving."));
	
	if (y0z && x0z)
		throw Exc(t_("HAMS only allows one symmetry at a time, either X or Y."));
	
	
	const UVector<Panel> &panels = surf.panels;
	const UVector<Point3D> &nodes = surf.nodes;
	
	if (type == 'y')
		out << "    --------------Waterplane Mesh File---------------";
	else
		out << "    --------------Hull Mesh File---------------";
	
	out	<< "\n \n    # Number of Panels, Nodes, X-Symmetry and Y-Symmetry\n";
	out << Format("       %5d       %5d       %5d       %5d\n \n", panels.size(), nodes.size(), 
																	y0z ? 1 : 0, x0z ? 1 : 0);
	out << "    # Start Definition of Node Coordinates     ! node_number   x   y   z\n";
	for (int in = 0; in < nodes.size(); ++in) {
		const Point3D &p = nodes[in]; 
		out << Format("%5d         % .6f         % .6f          % .6f\n", in+1, p.x, p.y, p.z);
	}
	out << "    # End Definition of Node Coordinates\n \n";
	out << "  # Start Definition of Node Relations   ! panel_number  number_of_vertices   Vertex1_ID   Vertex2_ID   Vertex3_ID   (Vertex4_ID\n";
	
	for (int ip = 0; ip < panels.size(); ++ip) {
		const Panel &pan = panels[ip];
		int num = pan.GetNumNodes();

		out << Format("%5d    %5d         %5d         %5d         %5d", ip+1, num, 
						pan.id[0]+1, pan.id[1]+1, pan.id[2]+1);
		if (num == 4)
			out << Format("         %5d", pan.id[3]+1);
		out << "\n";
	}	 
	out << "    # End Definition of Node Relations\n \n";
	if (type == 'y')
		out << "    --------------End Waterplane Mesh File---------------\n";
	else
		out << "    --------------End Hull Mesh File---------------\n";
}
