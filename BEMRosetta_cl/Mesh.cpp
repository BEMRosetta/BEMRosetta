// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
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
	} else if (ext == ".msh") {
		try {
			LoadTDynMsh(file, mesh);
		} catch(Exc e) {
			return std::move(e);
		}
		SetCode(Mesh::MSH_TDYN);
	} else if (ext == ".mesh") {
		try {
			LoadMesh(file, mesh, mass, cg);
		} catch(Exc e) {
			return std::move(e);
		}
		SetCode(Mesh::BEM_MESH);
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
		else if (ext == ".mesh")
			type = BEM_MESH;
		else
			throw Exc(Format(t_("Conversion to file type '%s' not supported"), file));
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
		SaveStlBin(file, surf, 1000);
	else if (type == STL_TXT)		
		SaveStlTxt(file, surf, 1000);
	else if (type == BEM_MESH)		
		SaveMesh(file, surf, mass, cg);
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
	if (mesh.IsEmpty())
		return;
	if (!onlyCG) {
		mesh.GetPanelParams();
		mesh.GetEnvelope();
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
		cgZ0surface = under.GetSurfaceZProjectionCG();
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
	under.GetHydrostaticStiffness(C, c0, cg, cb, rho, g, mass);
}

void Mesh::Report(double rho) const {
	BEM::Print("\n\n" + Format(t_("Mesh file '%s'"), fileName));
	
	BEM::Print(S("\n") + Format(t_("Limits [m] (%f - %f, %f - %f, %f - %f)"), 
			mesh.env.minX, mesh.env.maxX, mesh.env.minY, mesh.env.maxY, mesh.env.minZ, mesh.env.maxZ));
	BEM::Print(S("\n") + Format(t_("Water-plane area. Surface projection Z-axis [m2] %f - %f = %f"), -zProjectionPos, zProjectionNeg, zProjectionPos + zProjectionNeg));
	BEM::Print(S("\n") + Format(t_("Surface projection X-axis [m2] %f - %f = %f"), -xProjectionPos, xProjectionNeg, xProjectionPos + xProjectionNeg));
	BEM::Print(S("\n") + Format(t_("Surface projection Y-axis [m2] %f - %f = %f"), -yProjectionPos, yProjectionNeg, yProjectionPos + yProjectionNeg));
	BEM::Print(S("\n") + Format(t_("Surface [m2] %f"), mesh.surface));
	BEM::Print(S("\n") + Format(t_("Volume [m3] %f"), mesh.volume));
	BEM::Print(S("\n") + Format(t_("Underwater surface [m2] %f"), under.surface));
	BEM::Print(S("\n") + Format(t_("Underwater volume [m3] %f"), under.volume));
	BEM::Print(S("\n") + Format(t_("Displacement [tm] %f"), under.volume*rho/1000));
	BEM::Print(S("\n") + Format(t_("Centre of buoyancy [m] (%f, %f, %f)"), cb.x, cb.y, cb.z));
	
	BEM::Print(S("\n") + Format(t_("Loaded %d panels and %d nodes"), mesh.panels.size(), mesh.nodes.size()));
}

bool Mesh::IsSymmetricX() {
	return abs(cb.x)/abs(mesh.env.maxX - mesh.env.minX) < 0.001 && abs(mesh.env.maxX + mesh.env.minX) < 0.001;
}

bool Mesh::IsSymmetricY() {
	return abs(cb.y)/abs(mesh.env.maxY - mesh.env.minY) < 0.001 && abs(mesh.env.maxY + mesh.env.minY) < 0.001;
}

double Mesh::GMroll(double rho, double g) const {
	return C(3, 3)/(rho*g*under.volume);
}

double Mesh::GMpitch(double rho, double g) const {
	return C(4, 4)/(rho*g*under.volume);
}
		
void Mesh::GZ(double from, double to, double delta, double angleCalc, double rho, double g,
				Vector<double> &dataangle, Vector<double> &datagz) {
	Vector<double> dataMoment, vol, disp, wett, wplane, draft;
	Vector<Point3D> dcb, dcg;
	GZ(from, to, delta, angleCalc, rho, g, Null, dataangle, datagz, dataMoment, 
		vol, disp, wett, wplane, draft, dcb, dcg);
}

void Mesh::GZ(double from, double to, double delta, double angleCalc, double rho, double g,
	Function <bool(String, int pos)> Status, 
	Vector<double> &dataangle, Vector<double> &datagz, Vector<double> &dataMoment,
	Vector<double> &vol, Vector<double> &disp, Vector<double> &wett, Vector<double> &wplane,
	Vector<double> &draft, Vector<Point3D> &dcb, Vector<Point3D> &dcg) {
	
	dataangle.Clear();
	datagz.Clear();
	dataMoment.Clear();
	vol.Clear();
	disp.Clear();
	wett.Clear();
	wplane.Clear();
	draft.Clear();
	dcb.Clear();
	dcg.Clear();
	
	Surface base0 = clone(mesh);
	Point3D cg0 = clone(cg);
	Point3D c00 = clone(c0);

	base0.Rotate(0, 0, angleCalc, c00.x, c00.y, c00.z);
	cg0.Rotate(0, 0, angleCalc, c00.x, c00.y, c00.z);
	
	double dz = 0.1;
	for (double angle = from; angle <= to; angle += delta) {
		if (Status && !Status("", int(100*(angle - from)/(to - from))))
			throw Exc(t_("Cancelled by the user"));
				
		Surface base = clone(base0);
		Point3D cg = clone(cg0);
		
		base.Rotate(0, angle, 0, c00.x, c00.y, c00.z);
		cg.Rotate(0, angle, 0, c00.x, c00.y, c00.z);
		
		Surface under;
		if (!base.TranslateArchimede(mass, rho, dz, under))
			throw Exc(t_("Problem obtaining GZ"));
		
		cg.Translate(0, 0, dz);
		
		Point3D cb = under.GetCenterOfBuoyancy();
		
		Force6 fcb, fcg;
		under.GetHydrostaticForceCB(fcb, c00, cb, rho, g);	
		Surface::GetMassForce(fcg, c00, cg, mass, g);
	
		double moment = -(fcg.ry + fcb.ry);
		double gz = moment/mass/g;
		
		dataangle << angle;
		datagz << gz;
		dataMoment << moment;
		vol << under.volume;
		disp << under.volume*rho;
		wett << under.GetSurface();
		wplane << under.GetWaterPlaneArea();
		draft << under.GetEnvelope().minZ;
		dcb << cb;
		dcg << cg;
	}	
}
	