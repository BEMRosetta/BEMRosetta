// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"


String NemohMesh::LoadDatFS(UArray<Mesh> &_mesh, String fileName, bool &x0z) {
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		return Format(t_("Impossible to open file '%s'"), fileName);
	
	try {
		String line;
		LineParser f(in);	
		f.IsSeparator = IsTabSpace;
			
		line = in.GetLine();	
		f.Load(line);
	
		if (f.size() != 4 || !f.IsInt(0) || (f.GetInt(0) != 1 && f.GetInt(0) != 0) || !f.IsInt(1) || !f.IsInt(2) || !f.IsInt(3))
			return t_("Format error in Nemoh .dat mesh file");	// To detect Nemoh format
		
		if (f.GetInt(0) == 1)
			x0z = true;
		
		Mesh &msh = _mesh.Add();
		msh.fileName = fileName;
		msh.SetCode(Mesh::NEMOHFS_DAT);

		bool isfirst = true;		
		while(true) {
			line = in.GetLine();
			if (in.IsEof())
				break;	
			f.Load(line);
			int id = f.GetInt(0);	
			if (id == 0)
				break;
			if (isfirst && id != 1)
				return t_("Format error in Nemoh .dat mesh file");	// To detect Nemoh format
			
			Point3D &node = msh.mesh.nodes.Add();
			node.x = f.GetDouble(1);
			node.y = f.GetDouble(2);
			node.z = f.GetDouble(3); 
			
			if (isfirst)
				isfirst = false;
		}
		while(true) {
			line = in.GetLine();	
			if (in.IsEof())
				break;
			f.Load(line);
			if (f.size() != 4)
				break;
			Panel &panel = msh.mesh.panels.Add();
			panel.id[0] = f.GetInt(0)-1;
			panel.id[1] = f.GetInt(1)-1;	
			panel.id[2] = f.GetInt(2)-1;	
			panel.id[3] = f.GetInt(3)-1;	
		}	
		while(true) {
			line = in.GetLine();	
			if (in.IsEof())
				break;
			f.Load(line);
			if (f.size() != 2)
				break;
			int id0 = f.GetInt(0);	
			if (id0 == 0)
				break;
			LineSegment &seg = msh.mesh.segments.Add();
			seg.inode0 = id0-1;
			seg.inode1 = f.GetInt(1)-1;	
		}	
	} catch (Exc e) {
		return t_("Parsing error: ") + e;
	}
	
	return String();
}	
		
String NemohMesh::LoadDat(UArray<Mesh> &mesh, String fileName, bool &x0z) {
	Mesh &msh = mesh.Add();
	String ret = static_cast<NemohMesh&>(msh).LoadDat0(fileName, x0z);
	
	if (!ret.IsEmpty() && !ret.StartsWith(t_("Parsing error: "))) {
		mesh.Clear();
		Mesh::DecrementIdCount();
		return ret;
	}
	MatrixXd cg_(3, 1), cb_(3, 1);
	UVector<double> Vo;
	if (!Nemoh::Load_Hydrostatics_static(GetFileFolder(fileName), 1, cg_, cb_, Vo))
		return ret;
	
	msh.cg.x = cg_(0, 0);		//Supposed 1 body
	msh.cg.y = cg_(1, 0);		
	msh.cg.z = cg_(2, 0);		
	msh.cb.x = cb_(0, 0);		//Supposed 1 body
	msh.cb.y = cb_(1, 0);		
	msh.cb.z = cb_(2, 0);
	
	return ret;
}
	
String NemohMesh::LoadDat0(String fileName, bool &x0z) {
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		return Format(t_("Impossible to open file '%s'"), fileName);
	
	this->fileName = fileName;
	SetCode(Mesh::NEMOH_DAT);
	
	try {
		String line;
		LineParser f(in);	
		f.IsSeparator = IsTabSpace;
			
		line = in.GetLine();	
		f.Load(line);
	
		if (f.size() != 2 || !f.IsInt(0) || f.GetInt(0) != 2)
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
		while(!in.IsEof()) {
			line = in.GetLine();	

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

void NemohMesh::SaveDat(const Mesh &mesh, String fileName, const Surface &surf, bool x0z, int &npanels) {
	SaveDat0(fileName, surf, x0z, npanels);
	
	MatrixXd cg_(3, 1), cb_(3, 1);
	UVector<double> Vo;
	
	if (!IsNull(mesh.cg) && !IsNull(mesh.cb)) {
		cg_(0, 0) = mesh.cg.x;		//Supposed 1 body
		cg_(1, 0) = mesh.cg.y;		
		cg_(2, 0) = mesh.cg.z;		
		cb_(0, 0) = mesh.cb.x;		//Supposed 1 body
		cb_(1, 0) = mesh.cb.y;		
		cb_(2, 0) = mesh.cb.z;	
		Nemoh::Save_Hydrostatics_static(GetFileFolder(fileName), 1, cg_, cb_, Vo);
	}
}

void NemohMesh::SaveDat0(String fileName, const Surface &surf, bool x0z, int &npanels) {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'\n"), fileName));	
	
	npanels = 0;
	
	const UVector<Panel> &panels = surf.panels;
	const UVector<Point3D> &nodes = surf.nodes;
	
	out << "    2   " << (x0z ? "1" : "0") << "\n";
	for (int i = 0; i < nodes.size(); ++i) {
		const Point3D &node = nodes[i];
		out << Format("  %8d    % 014.7E   %0 14.7E   % 014.7E\n", i+1, node.x, node.y, node.z);
	}
	out << Format("  %8d    % 014.7E   %0 14.7E   % 014.7E\n", 0, 0, 0, 0);
	
	for (int i = 0; i < panels.size(); ++i) {
		const Panel &panel = panels[i];
		if (panel.surface0+panel.surface1 >= 1.0E-07) {
			out << Format("  %8d   %8d   %8d   %8d\n", panel.id[0]+1, panel.id[1]+1, panel.id[2]+1, panel.id[3]+1);
			npanels++;
		}
	}
	out << Format("  %8d   %8d   %8d   %8d\n", 0, 0, 0, 0);
}

void NemohMesh::SavePreMesh(String fileName, const Surface &surf) {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'\n"), fileName));	
	
	const UVector<Panel> &panels = surf.panels;
	const UVector<Point3D> &nodes = surf.nodes;
	
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

void NemohMesh::SaveKH(String fileName) const {
	Nemoh::Save_6x6(C, fileName);
}
/*
void NemohMesh::SaveInertia(String fileName) const {
	Nemoh::Save_Inertia_static(M, fileName);
}*/

