#include "BEMRosetta.h"

	
String MeshData::Load(String fileName, double rho, double g) {
	String ext = ToLower(GetFileExt(fileName));
	String ret;
	bool y0z = false, x0z = false;
	if (ext == ".dat") {
		ret = LoadDatNemoh(fileName, x0z);
		if (!ret.IsEmpty()) 
			ret = LoadDatWamit(fileName);
		else
			ret = String();
	} else if (ext == ".gdf") 
		ret = LoadGdfWamit(fileName, y0z, x0z); 
	else if (ext == ".stl") {
		bool isText;
		ret = LoadStlTxt(fileName, isText);
		if (!ret.IsEmpty() && !isText) 
			ret = LoadStlBin(fileName);
		else
			ret = String(); 
	} else
		return t_("Unknown file extension");	
	
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

void MeshData::SaveAs(String file, MESH_FMT type, double g, int meshType) {
	const Vector<Panel> &panels = meshType < 2 ? mesh.panels : under.panels;
	const Vector<Point3D> &nodes = [&]()->const Vector<Point3D> & {
		switch(meshType) {
		case 0:		return mesh.nodes0;
		case 1:		return mesh.nodes;
		default:	return under.nodes;
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
	under.GetHydrostaticStiffness(c, cb, rho, cg, mass, g, 0);
}

void MeshData::Report(double rho) {
	BEMData::Print("\n\n" + Format(t_("Loaded mesh '%s'"), file));
	
	BEMData::Print(x_("\n") + Format(t_("Limits [m] (%f - %f, %f - %f, %f - %f)"), 
			mesh.env.minX, mesh.env.maxX, mesh.env.minY, mesh.env.maxY, mesh.env.minZ, mesh.env.maxZ));
	BEMData::Print(x_("\n") + Format(t_("Water-plane area [m2] %f"), waterPlaneArea));
	BEMData::Print(x_("\n") + Format(t_("Surface [m2] %f"), mesh.surface));
	BEMData::Print(x_("\n") + Format(t_("Volume [m3] %f"), mesh.volume));
	BEMData::Print(x_("\n") + Format(t_("Underwater surface [m2] %f"), under.surface));
	BEMData::Print(x_("\n") + Format(t_("Underwater volume [m3] %f"), under.volume));
	BEMData::Print(x_("\n") + Format(t_("Displacement [tm] %f"), under.volume*rho/1000));
	BEMData::Print(x_("\n") + Format(t_("Center of buoyancy [m] (%f, %f, %f)"), cb.x, cb.y, cb.z));
	
	BEMData::Print(x_("\n") + Format(t_("Loaded %d panels and %d nodes"), mesh.panels.GetCount(), mesh.nodes.GetCount()));
}