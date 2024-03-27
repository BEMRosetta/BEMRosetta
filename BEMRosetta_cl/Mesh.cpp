// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"

int Mesh::idCount = 0;


void Mesh::Copy(const Mesh &msh) {
	xProjectionPos = msh.xProjectionPos;
	xProjectionNeg = msh.xProjectionNeg;
	yProjectionPos = msh.yProjectionPos;
	yProjectionNeg = msh.yProjectionNeg;
	zProjectionPos = msh.zProjectionPos;
	zProjectionNeg = msh.zProjectionNeg;
	
	cgZ0surface = clone(msh.cgZ0surface);
	cb = clone(msh.cb);
	cg = clone(msh.cg);
	cg0 = clone(msh.cg0);
	c0 = clone(msh.c0);
	
	M = clone(msh.M);
	C = clone(msh.C);
	
	name = msh.name;
	fileName = msh.fileName;
	header = msh.header;
	
	mesh = clone(msh.mesh);
	under = clone(msh.under);
	mesh0 = clone(msh.mesh0);
	
	code = msh.code;
	id = msh.id;
}

String Mesh::Load(Mesh &mesh, String file, double rho, double g, bool cleanPanels, double grid, double eps) {
	bool y0z, x0z;
	return Load(mesh, file, rho, g, cleanPanels, grid, eps, y0z, x0z);
}
	
String Mesh::Load(Mesh &mesh, String file, double rho, double g, bool cleanPanels, double grid, double eps, bool &y0z, bool &x0z) {
	UArray<Mesh> msh;
	String ret = Load(msh, file, rho, g, cleanPanels, grid, eps, y0z, x0z);
	if (!ret.IsEmpty())
		return ret;
	mesh = pick(First(msh));
	return ret;
}
	
String Mesh::Load(UArray<Mesh> &mesh, String file, double rho, double g, bool cleanPanels, double grid, double eps) {
	bool y0z, x0z;
	return Load(mesh, file, rho, g, cleanPanels, grid, eps, y0z, x0z);
}
	
String Mesh::Load(UArray<Mesh> &mesh, String file, double rho, double g, bool cleanPanels, double grid, double eps, bool &y0z, bool &x0z) {
	String ext = ToLower(GetFileExt(file));
	String ret;
	y0z = x0z = false;
	
	mesh.Clear();
	
	if (ext == ".dat") {
		ret = NemohMesh::LoadDat(mesh, file, x0z);
		if (!ret.IsEmpty() && !ret.StartsWith(t_("Parsing error: "))) {
			ret = NemohMesh::LoadDatFS(mesh, file, x0z);
			if (!ret.IsEmpty() && !ret.StartsWith(t_("Parsing error: "))) {
				ret = SalomeMesh::LoadDat(mesh, file);
				if (!ret.IsEmpty() && !ret.StartsWith(t_("Parsing error: "))) {
					ret = WamitMesh::LoadDat(mesh, file);
					if (!ret.IsEmpty() && !ret.StartsWith(t_("Parsing error: "))) {
						ret = DiodoreMesh::LoadDat(mesh, file);
						if (!ret.IsEmpty() && !ret.StartsWith(t_("Parsing error: "))) 	
							ret = AQWAMesh::LoadDat(mesh, file, y0z, x0z);
					}
				}
			}
		}
	} else if (ext == ".lis")
		ret = AQWAMesh::Load_LIS(mesh, file, g, y0z, x0z);
	else if (ext == ".txt") 
		ret = DiodoreMesh::LoadDat(mesh, file); 
	else if (ext == ".gdf") 
		ret = WamitMesh::LoadGdf(mesh, file, y0z, x0z); 
	else if (ext == ".pnl") 
		ret = HAMSMesh::LoadPnl(mesh, file, y0z, x0z); 
	else if (ext == ".hst") 
		ret = HydrostarMesh::LoadHst(mesh, file, y0z, x0z); 
	else if (ext == ".stl") {
		bool isText;
		Mesh &m = mesh.Add();
		try {
			LoadStl(file, m.mesh, isText, m.header);
		} catch(Exc e) {
			return std::move(e);
		}
		m.SetCode(isText ? Mesh::STL_TXT : Mesh::STL_BIN);
	} else if (ext == ".msh") {
		Mesh &m = mesh.Add();
		try {
			LoadTDynMsh(file, m.mesh);
		} catch(Exc e) {
			return std::move(e);
		}
		m.SetCode(Mesh::MSH_TDYN);
	} else if (ext == ".mesh" || ext == ".bem") {
		Mesh &m = mesh.Add();
		try {
			m.mesh.Load(file);
		} catch(Exc e) {
			return std::move(e);
		}
		m.SetCode(Mesh::BEM_MESH);
	} else
		ret = Format(t_("Unknown mesh file format '%s'"), file);	
	
	if (!ret.IsEmpty())
		return ret;
	
	for (Mesh &m : mesh) {
		ret = m.mesh.CheckErrors();
		if (!ret.IsEmpty())
			return ret;
		
		m.fileName = file;
		if (m.name.IsEmpty())
			m.name = InitCaps(GetFileTitle(file));
		
		if (y0z)
			m.mesh.DeployXSymmetry();
		if (x0z)
			m.mesh.DeployYSymmetry();	
		
		if (cleanPanels) 
			m.mesh.Heal(true, grid, eps);
		
		if (!IsNull(rho))
			m.AfterLoad(rho, g, false, true);
	}
	return String();
}

void Mesh::SaveAs(const UArray<Mesh*> &meshes, String file, MESH_FMT type, MESH_TYPE meshType, double rho, double g, bool symX, bool symY, 
						int &nNodes, int &nPanels) {
	UArray<Surface> surfs(meshes.size());
	nNodes = nPanels = 0;
	
	for (int i = 0; i < meshes.size(); ++i) {
		Surface &surf = surfs[i];
		if (meshType == UNDERWATER) 
			surf = clone(meshes[i]->under);
		else {
			if (type == AQWA_DAT) {		// Appends dry and wet sides. This way there are no panels between dry and wet side
				surf = clone(meshes[i]->under);		// First the wet
				Surface dry;	
				dry.CutZ(meshes[i]->mesh, 1);
				surf.Append(dry);					// Next the dry
			} else
				surf = clone(meshes[i]->mesh);
		}
		
		if (symX && (type == WAMIT_GDF || type == HAMS_PNL || type == DIODORE_DAT || type == AQWA_DAT)) {
			Surface nsurf;
			nsurf.CutX(surf);
			surf = pick(nsurf);
		}
		if (symY && (type == WAMIT_GDF || type == NEMOH_DAT || type == NEMOH_PRE || 
					 type == HAMS_PNL || type == DIODORE_DAT || type == AQWA_DAT)) {
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
		if (surf.panels.IsEmpty() && surf.lines.IsEmpty())
			throw Exc(t_("Model is empty. No data found"));

		nNodes += surf.nodes.size();
		nPanels += surf.panels.size();
	}
	
	if (type == WAMIT_GDF) 
		WamitMesh::SaveGdf(file, First(surfs), g, symX, symY);	
	else if (type == NEMOH_DAT) 
		NemohMesh::SaveDat(*(First(meshes)), file, First(surfs), symY, nPanels);
	else if (type == NEMOH_PRE) 
		NemohMesh::SavePreMesh(file, First(surfs));
	else if (type == HAMS_PNL)		
		HAMSMesh::SavePnl(file, First(surfs), symX, symY);	// Only one symmetry is really available
	else if (type == AQWA_DAT)		
		AQWAMesh::SaveDat(file, meshes, surfs, rho, g, symX, symY);	
	else if (type == DIODORE_DAT) 
		DiodoreMesh::SaveDat(file, First(surfs));
	else if (type == STL_BIN)		
		SaveStlBin(file, First(surfs));
	else if (type == STL_TXT)		
		SaveStlTxt(file, First(surfs));
	else if (type == BEM_MESH)		
		First(surfs).Save(file);
	else
		throw Exc(t_("Unknown mesh file type"));
}

String Mesh::Heal(bool basic, double rho, double g, double grid, double eps, Function <bool(String, int pos)> Status) {
	String ret = mesh.Heal(basic, grid, eps, Status);
	if (!ret.IsEmpty())
		return ret;
	
	AfterLoad(rho, g, false, false);
	
	return String();
}

void Mesh::Orient() {
	mesh.Orient();
}

void Mesh::Append(const Surface &orig, double rho, double g) {
	mesh.Append(orig);
	
	AfterLoad(rho, g, false, true);
}

void Mesh::Reset(double rho, double g) {
	mesh = clone(mesh0);
	cg = clone(cg0);
	AfterLoad(rho, g, false, false);
}
	
void Mesh::AfterLoad(double rho, double g, bool onlyCG, bool isFirstTime, bool massBuoy) {
	if (mesh.IsEmpty()) {
		if (!mesh.lines.IsEmpty())
			mesh.GetEnvelope();
		return;
	}
	if (!onlyCG) {
		mesh.GetPanelParams();
		mesh.GetEnvelope();
		mesh.GetArea();
		mesh.GetVolume();
		
		under.CutZ(mesh, -1);
		under.GetPanelParams();
		xProjectionPos = under.GetAreaXProjection(true, false);
		xProjectionNeg = under.GetAreaXProjection(false, true);
		yProjectionPos = under.GetAreaYProjection(true, false);
		yProjectionNeg = under.GetAreaYProjection(false, true);
		zProjectionPos = under.GetAreaZProjection(true, false);
		zProjectionNeg = under.GetAreaZProjection(false, true);
		cgZ0surface = under.GetAreaZProjectionCG();
		under.GetArea();
		under.GetVolume();
		
		if (M.size() != 36)
			M = MatrixXd::Zero(6,6);
		if (GetMass() == 0 && !IsNull(rho))
			SetMass(under.volume*rho);
		cb = under.GetCentreOfBuoyancy();
	}
	if (isFirstTime) {
		mesh0 = clone(mesh);
		cg0 = clone(cg);
	}
	if (!IsNull(rho) && !IsNull(g))
		under.GetHydrostaticStiffness(C, c0, cg, cb, rho, g, GetMass(), massBuoy);
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
	if (!IsNum(rho) || !IsNum(g) || under.volume == 0 || C.size() == 0)
		return Null;
	return C(3, 3)/(rho*g*under.volume);
}

double Mesh::GMpitch(double rho, double g) const {
	if (!IsNum(rho) || !IsNum(g) || under.volume == 0 || C.size() == 0)
		return Null;
	return C(4, 4)/(rho*g*under.volume);
}
		
void Mesh::GZ(double from, double to, double delta, double angleCalc, double rho, double g,
  double tolerance, UVector<double> &dataangle, UVector<double> &datagz, String &error) {
	UVector<double> dataMoment, vol, disp, wett, wplane, draft;
	UVector<Point3D> dcb, dcg;
	GZ(from, to, delta, angleCalc, rho, g, tolerance, Null, dataangle, datagz, dataMoment, 
		vol, disp, wett, wplane, draft, dcb, dcg, error);
}

void Mesh::GZ(double from, double to, double delta, double angleCalc, double rho, double g, double tolerance,
	Function <bool(String, int pos)> Status, 
	UVector<double> &dataangle, UVector<double> &datagz, UVector<double> &dataMoment,
	UVector<double> &vol, UVector<double> &disp, UVector<double> &wett, UVector<double> &wplane,
	UVector<double> &draft, UVector<Point3D> &dcb, UVector<Point3D> &dcg, String &error) {
	
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

	base0.Rotate(0, 0, ToRad(angleCalc), c0.x, c0.y, c0.z);
	cg0.Rotate(0, 0, ToRad(angleCalc), c0.x, c0.y, c0.z);
	
	double dz = 0.1;
	for (double angle = from; angle <= to; angle += delta) {
		double progress =  to == from ? 1 : (angle - from)/(to - from);
		if (Status && !Status("", int(100*progress)))
			throw Exc(t_("Cancelled by the user"));
				
		Surface base = clone(base0);
		Point3D cg = clone(cg0);
		
		base.Rotate(0, ToRad(angle), 0, c0.x, c0.y, c0.z);
		cg.Rotate(0, ToRad(angle), 0, c0.x, c0.y, c0.z);
		
		if (GetMass() == 0)
			throw Exc(t_("Problem obtaining GZ. Mass is zero"));
		
		Surface under;
		if (!base.TranslateArchimede(GetMass(), rho, dz, under))
			throw Exc(t_("Problem obtaining GZ"));
		
		cg.Translate(0, 0, dz);
		
		if (under.VolumeMatch(tolerance, tolerance) < 0) {
			if (!error.IsEmpty())
				error << "\n";
			error << Format("Around %.2f, angle %.2f", angleCalc, angle);
			
			dataangle << angle;
			datagz << Null;
			dataMoment << Null;
			vol << Null;
			disp << Null;
			wett << Null;
			wplane << Null;
			draft << Null;
			dcb << Null;
			dcg << cg;
			
		} else {
			Point3D cb = under.GetCentreOfBuoyancy();
			Force6D fcb = under.GetHydrostaticForceCB(c0, cb, rho, g);
			//base.GetPanelParams();	
			//Force6D fcb = base.GetHydrodynamicForce(c0, true, 
			//				[&](double x, double y)->double {return 0;}, 
			//				[&](double x, double y, double z, double et)->double {return z > 0 ? 0 : rho*g*z;});
			Force6D fcg = Surface::GetMassForce(c0, cg, GetMass(), g);
		
			double moment = -(fcg.r.y + fcb.r.y);
			double gz = moment/GetMass()/g;
			
			dataangle << angle;
			datagz << gz;
			dataMoment << moment;
			vol << under.volume;
			disp << under.volume*rho;
			wett << under.GetArea();
			wplane << under.GetWaterPlaneArea();
			draft << under.GetEnvelope().minZ;
			dcb << cb;
			dcg << cg;
		}
	}	
}

void Mesh::Move(double dx, double dy, double dz, double ax, double ay, double az, double rho, double g, bool setnewzero) {
	mesh = clone(mesh0);
	cg = clone(cg0);					
	mesh.TransRot(dx, dy, dz, ax, ay, az, c0.x, c0.y, c0.z);
	cg.TransRot(dx, dy, dz, ax, ay, az, c0.x, c0.y, c0.z);
	AfterLoad(rho, g, false, setnewzero);	
}

void Mesh::Move(const double *pos, double rho, double g, bool setnewzero) {
	Move(pos[0], pos[1], pos[2], pos[3], pos[4], pos[5], rho, g, setnewzero);	
}

void Mesh::Move(const float *pos, double rho, double g, bool setnewzero) {
	Move(pos[0], pos[1], pos[2], pos[3], pos[4], pos[5], rho, g, setnewzero);	
}

void Mesh::SetMass(double m) {
	if (M.size() != 36)
		M = MatrixXd::Zero(6,6);	
	double oldm = M(0,0);
	if (oldm != 0) {
		double factor = m/oldm;
		M.array() *= factor;
	} else
		M(0,0) = M(1,1) = M(2,2) = m;
}