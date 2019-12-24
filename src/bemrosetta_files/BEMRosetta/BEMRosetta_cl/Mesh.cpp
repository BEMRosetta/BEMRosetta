#include "BEMRosetta.h"

int MeshData::idCount = 0;
	
String MeshData::Load(String file, double rho, double g) {
	String ext = ToLower(GetFileExt(file));
	String ret;
	bool y0z = false, x0z = false;
	if (ext == ".dat") {
		ret = LoadDatNemoh(file, x0z);
		if (!ret.IsEmpty() && !ret.StartsWith(t_("Parsing error: "))) 
			ret = LoadDatWamit(file);
	} else if (ext == ".gdf") 
		ret = LoadGdfWamit(file, y0z, x0z); 
	else if (ext == ".stl") {
		bool isText;
		ret = LoadStlTxt(file, isText);
		if (!ret.IsEmpty() && !isText) 
			ret = LoadStlBin(file);
		else
			ret = String(); 
	} else
		throw Exc(Format(t_("Unknown MESH file extension in '%s'"), file));	
	
	if (!ret.IsEmpty())
		return ret;
	
	ret = mesh.CheckErrors();
	if (!ret.IsEmpty())
		return ret;
		
	if (y0z)
		mesh.DeployXSymmetry();
	if (x0z)
		mesh.DeployYSymmetry();	
	
	cg = cg0 = Point3D(0, 0, 0);
	mass = Null;
	
	mesh.nodes = clone(mesh.nodes0);
	
	AfterLoad(rho, g, false);
	
	return String();
}

void MeshData::SaveAs(String file, MESH_FMT type, double g, MESH_TYPE meshType) {
	Vector<Panel> panelsRw;
	Vector<Point3D> nodesRw;
	if (meshType == UNDERWATER) {		// Some healing before saving
		panelsRw = clone(under.panels);
		nodesRw = clone(under.nodes);
		Surface::RemoveDuplicatedPanels(panelsRw);
		Surface::RemoveDuplicatedPointsAndRenumber(panelsRw, nodesRw);
		Surface::RemoveDuplicatedPanels(panelsRw);
	}
	const Vector<Panel> &panels = meshType != UNDERWATER ? mesh.panels : panelsRw;
	const Vector<Point3D> &nodes = [&]()->const Vector<Point3D> & {
		switch(meshType) {
		case 0:		return mesh.nodes0;
		case 1:		return mesh.nodes;
		default:	return nodesRw;
		}
	}();
	if (panels.IsEmpty())
		throw Exc(t_("Model is empty. No panels found"));
		
	if (type == UNKNOWN) {
		String ext = ToLower(GetFileExt(file));
		
		if (ext == ".gdf")
			type = MeshData::WAMIT_GDF;
		else if (ext == ".dat")
			type = MeshData::NEMOH_DAT;
		else if (ext == ".")
			type = MeshData::NEMOH_PRE;
		else if (ext == ".stl")
			type = MeshData::STL_TXT;
		else
			throw Exc(Format(t_("Conversion to type of file '%s' not supported"), file));
	}
	
	bool y0z = false, x0z = false;	// Pending
	
	if (type == WAMIT_GDF) 
		SaveGdfWamit(file, panels, nodes, g, y0z, x0z);	
	else if (type == NEMOH_DAT) 
		SaveDatNemoh(file, panels, nodes, x0z);
	else if (type == NEMOH_PRE) 
		SavePreMeshNemoh(file, panels, nodes);
	else if (type == STL_BIN)		
		SaveStlBin(file, panels, nodes);
	else if (type == STL_TXT)		
		SaveStlTxt(file, panels, nodes);
	else
		throw Exc(t_("Unknown mesh file type"));
}

String MeshData::Heal(Function <void(String, int pos)> Status) {
	String ret = mesh.Heal(Status);
	if (!ret.IsEmpty())
		return ret;
	
	return String();
}

void MeshData::AfterLoad(double rho, double g, bool onlyCG) {
	if (!onlyCG) {
		mesh.GetPanelParams();
		mesh.GetLimits();
		mesh.GetSurface();
		mesh.GetVolume();
		
		under.Underwater(mesh);
		under.GetPanelParams();
		waterPlaneArea = under.GetWaterPlaneArea();
		under.GetSurface();
		under.GetVolume();
		
		if (IsNull(mass))
			mass = under.volume*rho;
		cb = under.GetCenterOfBuoyancy();
	}
	under.GetHydrostaticStiffness(C, cb, rho, cg, mass, g, 0);
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
	
	BEMData::Print(S("\n") + Format(t_("Loaded %d panels and %d nodes"), mesh.panels.GetCount(), mesh.nodes.GetCount()));
}
