// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2021, the BEMRosetta author and contributors
#include "BEMRosetta.h"

int Mesh::idCount = 0;


String Mesh::Load(String file, double rho, double g, bool cleanPanels) {
	bool y0z, x0z;
	return Load(file, rho, g, cleanPanels, y0z, x0z);
}
	
String Mesh::Load(String file, double rho, double g, bool cleanPanels, bool &y0z, bool &x0z) {
	String ext = ToLower(GetFileExt(file));
	String ret;
	y0z = x0z = false;
	if (ext == ".dat") {
		ret = static_cast<NemohMesh&>(*this).LoadDat(file, x0z);
		if (!ret.IsEmpty() && !ret.StartsWith(t_("Parsing error: "))) {
			ret = static_cast<WamitMesh &>(*this).LoadDat(file);
			if (!ret.IsEmpty() && !ret.StartsWith(t_("Parsing error: "))) 
				ret = static_cast<AQWAMesh &>(*this).LoadDat(file);
		}
	} else if (ext == ".gdf") 
		ret = static_cast<WamitMesh &>(*this).LoadGdf(file, y0z, x0z); 
	else if (ext == ".pnl") 
		ret = static_cast<HAMSMesh&>(*this).LoadPnl(file, y0z, x0z); 
	else if (ext == ".stl") {
		bool isText;
		try {
			LoadStl(file, mesh, isText, header);
		} catch(Exc e) {
			return std::move(e);
		}
		SetCode(isText ? Mesh::STL_TXT : Mesh::STL_BIN);
	} else
		ret = Format(t_("Unknown mesh file format '%s'"), file);	
	
	if (!ret.IsEmpty())
		return ret;
	
	ret = mesh.CheckErrors();
	if (!ret.IsEmpty())
		return ret;
	
	fileName = file;
	name = InitCaps(GetFileTitle(file));
	
	if (y0z)
		mesh.DeployXSymmetry();
	if (x0z)
		mesh.DeployYSymmetry();	
	
	if (cleanPanels) {
		Surface::RemoveDuplicatedPanels(mesh.panels);
		Surface::RemoveDuplicatedPointsAndRenumber(mesh.panels, mesh.nodes);
		Surface::RemoveDuplicatedPanels(mesh.panels);
	}
	
	if (!IsNull(rho))
		AfterLoad(rho, g, false, true);
	
	return String();
}

void Mesh::SaveAs(String file, MESH_FMT type, double g, MESH_TYPE meshType, bool symX, bool symY, 
						int &nNodes, int &nPanels) {
	Surface surf;
	if (meshType == UNDERWATER) 
		surf = clone(under);
	else
		surf = clone(mesh);
	
	if (type == UNKNOWN) {
		String ext = ToLower(GetFileExt(file));
		
		if (ext == ".gdf")
			type = WAMIT_GDF;
		else if (ext == ".dat")
			type = NEMOH_DAT;
		else if (ext == ".")
			type = NEMOH_PRE;
		else if (ext == ".pnl")
			type = HAMS_PNL;
		else if (ext == ".stl")
			type = STL_TXT;
		else
			throw Exc(Format(t_("Conversion to type of file '%s' not supported"), file));
	}
	
	if (symX && (type == WAMIT_GDF || type == HAMS_PNL)) {
		Surface nsurf;
		nsurf.CutX(surf);
		surf = pick(nsurf);
	}
	if (symY && (type == WAMIT_GDF || type == NEMOH_DAT || type == NEMOH_PRE || type == HAMS_PNL)) {
		Surface nsurf;
		nsurf.CutY(surf);
		surf = pick(nsurf);
	}
	if (meshType == UNDERWATER || symX || symY) {// Some healing before saving
		Surface::RemoveDuplicatedPanels(surf.panels);
		Surface::RemoveDuplicatedPointsAndRenumber(surf.panels, surf.nodes);
		Surface::RemoveDuplicatedPanels(surf.panels);
		Surface::DetectTriBiP(surf.panels);
	}
	
	if (surf.panels.IsEmpty())
		throw Exc(t_("Model is empty. No panels found"));
	
	nNodes = surf.nodes.size();
	nPanels = surf.panels.size();
	
	if (type == WAMIT_GDF) 
		static_cast<const WamitMesh &>(*this).SaveGdf(file, surf, g, symX, symY);	
	else if (type == NEMOH_DAT) 
		static_cast<NemohMesh&>(*this).SaveDat(file, surf, symY);
	else if (type == NEMOH_PRE) 
		static_cast<NemohMesh&>(*this).SavePreMesh(file, surf);
	else if (type == HAMS_PNL)		
		static_cast<HAMSMesh&>(*this).SavePnl(file, surf, symX, symY);	// Only one symmetry really available
	else if (type == STL_BIN)		
		SaveStlBin(file, surf);
	else if (type == STL_TXT)		
		SaveStlTxt(file, surf);
	else
		throw Exc(t_("Unknown mesh file type"));
}

String Mesh::Heal(bool basic, double rho, double g, Function <bool(String, int pos)> Status) {
	String ret = mesh.Heal(basic, Status);
	if (!ret.IsEmpty())
		return ret;
	
	AfterLoad(rho, g, false, false);
	
	return String();
}

void Mesh::Orient() {
	mesh.Orient();
}

void Mesh::Join(const Surface &orig, double rho, double g) {
	mesh.Join(orig);
	
	AfterLoad(rho, g, false, true);
}

void Mesh::Reset(double rho, double g) {
	mesh = clone(mesh0);
	cg = clone(cg0);
	AfterLoad(rho, g, false, false);
}
	
void Mesh::AfterLoad(double rho, double g, bool onlyCG, bool isFirstTime) {
	if (!onlyCG) {
		mesh.GetPanelParams();
		mesh.GetLimits();
		mesh.GetSurface();
		mesh.GetVolume();
		
		under.CutZ(mesh, -1);
		under.GetPanelParams();
		xProjectionPos = under.GetSurfaceXProjection(true, false);
		xProjectionNeg = under.GetSurfaceXProjection(false, true);
		yProjectionPos = under.GetSurfaceYProjection(true, false);
		yProjectionNeg = under.GetSurfaceYProjection(false, true);
		zProjectionPos = under.GetSurfaceZProjection(true, false);
		zProjectionNeg = under.GetSurfaceZProjection(false, true);
		under.GetSurface();
		under.GetVolume();
		
		if (IsNull(mass))
			mass = under.volume*rho;
		cb = under.GetCenterOfBuoyancy();
	}
	if (isFirstTime) {
		mesh0 = clone(mesh);
		cg0 = clone(cg);
	}
	under.GetHydrostaticStiffness(C, cb, rho, cg, mass, g);
}

void Mesh::Report(double rho) const {
	BEMData::Print("\n\n" + Format(t_("Mesh file '%s'"), fileName));
	
	BEMData::Print(S("\n") + Format(t_("Limits [m] (%f - %f, %f - %f, %f - %f)"), 
			mesh.env.minX, mesh.env.maxX, mesh.env.minY, mesh.env.maxY, mesh.env.minZ, mesh.env.maxZ));
	BEMData::Print(S("\n") + Format(t_("Water-plane area. Surface projection Z-axis [m2] %f - %f = %f"), -zProjectionPos, zProjectionNeg, zProjectionPos + zProjectionNeg));
	BEMData::Print(S("\n") + Format(t_("Surface projection X-axis [m2] %f - %f = %f"), -xProjectionPos, xProjectionNeg, xProjectionPos + xProjectionNeg));
	BEMData::Print(S("\n") + Format(t_("Surface projection Y-axis [m2] %f - %f = %f"), -yProjectionPos, yProjectionNeg, yProjectionPos + yProjectionNeg));
	BEMData::Print(S("\n") + Format(t_("Surface [m2] %f"), mesh.surface));
	BEMData::Print(S("\n") + Format(t_("Volume [m3] %f"), mesh.volume));
	BEMData::Print(S("\n") + Format(t_("Underwater surface [m2] %f"), under.surface));
	BEMData::Print(S("\n") + Format(t_("Underwater volume [m3] %f"), under.volume));
	BEMData::Print(S("\n") + Format(t_("Displacement [tm] %f"), under.volume*rho/1000));
	BEMData::Print(S("\n") + Format(t_("Center of buoyancy [m] (%f, %f, %f)"), cb.x, cb.y, cb.z));
	
	BEMData::Print(S("\n") + Format(t_("Loaded %d panels and %d nodes"), mesh.panels.size(), mesh.nodes.size()));
}

bool Mesh::IsSymmetricX() {
	return abs(cb.x)/abs(mesh.env.maxX - mesh.env.minX) < 0.001 && abs(mesh.env.maxX + mesh.env.minX) < 0.001;
}

bool Mesh::IsSymmetricY() {
	return abs(cb.y)/abs(mesh.env.maxY - mesh.env.minY) < 0.001 && abs(mesh.env.maxY + mesh.env.minY) < 0.001;
}
