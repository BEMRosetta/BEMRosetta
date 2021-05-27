#include "BEMRosetta.h"

int MeshData::idCount = 0;

String MeshData::Load(String file) {
	bool y0z, x0z;
	return Load(file, Null, Null, false, y0z, x0z);	
}

String MeshData::Load(String file, bool &y0z, bool &x0z) {
	return Load(file, Null, Null, false, y0z, x0z);	
}

String MeshData::Load(String file, double rho, double g, bool cleanPanels) {
	bool y0z, x0z;
	return Load(file, rho, g, cleanPanels, y0z, x0z);
}
	
String MeshData::Load(String file, double rho, double g, bool cleanPanels, bool &y0z, bool &x0z) {
	String ext = ToLower(GetFileExt(file));
	String ret;
	y0z = x0z = false;
	if (ext == ".dat") {
		ret = LoadDatNemoh(file, x0z);
		if (!ret.IsEmpty() && !ret.StartsWith(t_("Parsing error: "))) {
			ret = LoadDatWamit(file);
			if (!ret.IsEmpty() && !ret.StartsWith(t_("Parsing error: "))) 
				ret = LoadDatAQWA(file);
		}
	} else if (ext == ".gdf") 
		ret = LoadGdfWamit(file, y0z, x0z); 
	else if (ext == ".pnl") 
		ret = LoadPnlHAMS(file, y0z, x0z); 
	else if (ext == ".stl") {
		bool isText;
		try {
			LoadStl(file, mesh, isText, header);
		} catch(Exc e) {
			return std::move(e);
		}
		SetCode(isText ? MeshData::STL_TXT : MeshData::STL_BIN);
	} else
		throw Exc(Format(t_("Unknown MESH file extension in '%s'"), file));	
	
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
		AfterLoad(rho, g, false);
	
	return String();
}

void MeshData::SaveAs(String file, MESH_FMT type, double g, MESH_TYPE meshType, bool symX, bool symY) {
	Surface surf;
	if (meshType == UNDERWATER) 
		surf = clone(under);
	else
		surf = clone(mesh);
	
	if (symX && (type == WAMIT_GDF || type == HAMS_PNL)) {
		Surface nsurf;
		nsurf.CutX(surf);
		surf = pick(nsurf);
	}
	if (symY && (type == WAMIT_GDF || type == NEMOH_DAT || type == HAMS_PNL)) {
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
	
	if (type == WAMIT_GDF) 
		SaveGdfWamit(file, surf, g, symX, symY);	
	else if (type == NEMOH_DAT) 
		SaveDatNemoh(file, surf, symY);
	else if (type == NEMOH_PRE) 
		SavePreMeshNemoh(file, surf);
	else if (type == HAMS_PNL)		
		SavePnlHAMS(file, surf, false, symY);	// Only one symmetry really available
	else if (type == STL_BIN)		
		SaveStlBin(file, surf);
	else if (type == STL_TXT)		
		SaveStlTxt(file, surf);
	else
		throw Exc(t_("Unknown mesh file type"));
}

String MeshData::Heal(bool basic, Function <void(String, int pos)> Status) {
	String ret = mesh.Heal(basic, Status);
	if (!ret.IsEmpty())
		return ret;
	
	return String();
}

void MeshData::Orient() {
	mesh.Orient();
}

void MeshData::Join(const Surface &orig, double rho, double g) {
	mesh.Join(orig);
	
	AfterLoad(rho, g, false);
}
	
void MeshData::AfterLoad(double rho, double g, bool onlyCG) {
	if (!onlyCG) {
		mesh.GetPanelParams();
		mesh.GetLimits();
		mesh.GetSurface();
		mesh.GetVolume();
		
		under.CutZ(mesh, -1);
		under.GetPanelParams();
		waterPlaneArea = under.GetWaterPlaneArea();
		under.GetSurface();
		under.GetVolume();
		
		if (IsNull(mass))
			mass = under.volume*rho;
		cb = under.GetCenterOfBuoyancy();
	}
	under.GetHydrostaticStiffness(C, cb, rho, cg, mass, g);
}

void MeshData::Report(double rho) {
	BEMData::Print("\n\n" + Format(t_("Loaded mesh '%s'"), fileName));
	
	BEMData::Print(S("\n") + Format(t_("Limits [m] (%f - %f, %f - %f, %f - %f)"), 
			mesh.env.minX, mesh.env.maxX, mesh.env.minY, mesh.env.maxY, mesh.env.minZ, mesh.env.maxZ));
	BEMData::Print(S("\n") + Format(t_("Water-plane area [m2] %f"), waterPlaneArea));
	BEMData::Print(S("\n") + Format(t_("Surface [m2] %f"), mesh.surface));
	BEMData::Print(S("\n") + Format(t_("Volume [m3] %f"), mesh.volume));
	BEMData::Print(S("\n") + Format(t_("Underwater surface [m2] %f"), under.surface));
	BEMData::Print(S("\n") + Format(t_("Underwater volume [m3] %f"), under.volume));
	BEMData::Print(S("\n") + Format(t_("Displacement [tm] %f"), under.volume*rho/1000));
	BEMData::Print(S("\n") + Format(t_("Center of buoyancy [m] (%f, %f, %f)"), cb.x, cb.y, cb.z));
	
	BEMData::Print(S("\n") + Format(t_("Loaded %d panels and %d nodes"), mesh.panels.size(), mesh.nodes.size()));
}

bool MeshData::IsSymmetricX() {
	return abs(cb.x)/abs(mesh.env.maxX - mesh.env.minX) < 0.001 && abs(mesh.env.maxX + mesh.env.minX) < 0.001;
}

bool MeshData::IsSymmetricY() {
	return abs(cb.y)/abs(mesh.env.maxY - mesh.env.minY) < 0.001 && abs(mesh.env.maxY + mesh.env.minY) < 0.001;
}

void MeshData::SaveHST(String fileName, double rho, double g) const {
	Wamit::Save_hst_static(C, fileName, rho, g);
}