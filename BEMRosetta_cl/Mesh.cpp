// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"

	
// enum MESH_FMT 			    	  {WAMIT_GDF,  WAMIT_DAT,   NEMOH_DAT,   NEMOHFS_DAT,   NEMOH_PRE,      AQWA_DAT,   AQWA LIS, HAMS_PNL,  STL_BIN,     STL_TXT,   EDIT,  MSH_TDYN,   BEM_MESH, DIODORE_DAT,   HYDROSTAR_HST,    ORCA_OWR, 	   MIKE21_GRD	 CAPY_NC, 		OBJ,    ORCAFLEX_YML,   UNKNOWN, NUMMESH};	
const char *Body::meshStr[]         = {"Wamit.gdf","Wamit.dat",	"Nemoh.dat", "NemohFS.dat", "Nemoh premesh","AQWA.dat", "AQWA.lis","HAMS.pnl","STL.Binary","STL.Text","Edit","TDyn.msh", "BEMR",   "Diodore.dat", "HydroStar.hst", "OrcaWave.owr", "MIKE21.grd", "Capytaine.nc","Obj",  "OrcaFlex.yml", "Unknown"};	
const bool Body::meshCanSave[] 		= {true, 	   false,	    true,		 false,			false, 		    true,		false,	   true,	   true,		true,	   false, false, 	  true, 	true,		   false,   	   false, 		   true, 		 false, 	    false,  false, 	        false};       
const char *Body::meshExt[]	  		= {"*.gdf",    "*.dat",	 	"*.dat",	 "*.dat", 		"",		        "*.dat",	"*.lis",   "*.pnl",   "*.stl",     "*.stl",    "",	  "*.msh",   "*.bemr",  "*.dat", 	  "*.hst", 	   	   "*.owr",		   "*.grd", 	 "*.nc", 	    "*.obj","*.yml",        "*.*"};       

int Body::idCount = 0;


String Body::GetMeshExt() {
	String ret;
	
	for (int i = 0; i < MESH_FMT::NUMMESH; ++i) {
		if (ret.FindAfter(meshExt[i]) < 0) {
			if (!ret.IsEmpty())
				ret << " ";
			ret << meshExt[i];		
		}
	}
	return ret;
}
	
void Body::Copy(const Body &msh) {
	dt.Copy(msh.dt);
	cdt.Copy(msh.cdt);
}

void Body::Data::Copy(const Body::Data &msh) {
	projectionPos = clone(msh.projectionPos);
	projectionNeg = clone(msh.projectionNeg);
	
	cgZ0surface = clone(msh.cgZ0surface);
	cb = clone(msh.cb);
	cg = clone(msh.cg);
	cg0 = clone(msh.cg0);
	c0 = clone(msh.c0);
	Vo = msh.Vo;
	
	M = clone(msh.M);
	C = clone(msh.C);
	Cmoor = clone(msh.Cmoor);
	Cadd = clone(msh.Cadd);
	Dlin = clone(msh.Dlin);
	Dquad = clone(msh.Dquad);
	Aadd = clone(msh.Aadd);
	
	name = msh.name;
	fileName = msh.fileName;
	fileHeader = msh.fileHeader;
	
	mesh = clone(msh.mesh);
	under = clone(msh.under);
	mesh0 = clone(msh.mesh0);
		
	SetCode(msh.GetCode());
	SetId(msh.GetId());
}

void Body::ControlData::Copy(const Body::ControlData &msh) {
	controlPointsA = clone(msh.controlPointsA);
	controlPointsB = clone(msh.controlPointsB);
	controlPointsC = clone(msh.controlPointsC);
	controlLoads = clone(msh.controlLoads);
	
	controlPointsA0 = clone(msh.controlPointsA0);
	controlPointsB0 = clone(msh.controlPointsB0);
	controlPointsC0 = clone(msh.controlPointsC0);
	controlLoads0 = clone(msh.controlLoads0);
	damagedBodies = clone(msh.damagedBodies);
}

String Body::Load(Body &mesh, String file, double rho, double g, bool cleanPanels, double grid, double eps) {
	bool y0z, x0z;
	return Load(mesh, file, rho, g, cleanPanels, grid, eps, y0z, x0z);
}
	
String Body::Load(Body &mesh, String file, double rho, double g, bool cleanPanels, double grid, double eps, bool &y0z, bool &x0z) {
	UArray<Body> msh;
	String ret = Load(msh, file, rho, g, cleanPanels, grid, eps, y0z, x0z);
	if (!ret.IsEmpty())
		return ret;
	mesh = pick(First(msh));
	return ret;
}
	
String Body::Load(UArray<Body> &mesh, String file, double rho, double g, bool cleanPanels, double grid, double eps) {
	bool y0z, x0z;
	return Load(mesh, file, rho, g, cleanPanels, grid, eps, y0z, x0z);
}
	
String Body::Load(UArray<Body> &mesh, String file, double rho, double g, bool cleanPanels, double grid, double eps, bool &y0z, bool &x0z) {
	String ext = ToLower(GetFileExt(file));
	String ret;
	y0z = x0z = false;
	
	mesh.Clear();
	
	if (ext == ".dat") {
		ret = NemohBody::LoadDat(mesh, file, x0z);
		if (!ret.IsEmpty() && !ret.StartsWith(t_("Parsing error: "))) {
			ret = NemohBody::LoadDatFS(mesh, file, x0z);
			if (!ret.IsEmpty() && !ret.StartsWith(t_("Parsing error: "))) {
				ret = SalomeBody::LoadDat(mesh, file);
				if (!ret.IsEmpty() && !ret.StartsWith(t_("Parsing error: "))) {
					ret = WamitBody::LoadDat(mesh, file);
					if (!ret.IsEmpty() && !ret.StartsWith(t_("Parsing error: "))) {
						ret = DiodoreBody::LoadDat(mesh, file);
						if (!ret.IsEmpty() && !ret.StartsWith(t_("Parsing error: "))) { 	
							Hydro hy;
							ret = AQWABody::LoadDat(mesh, hy, file);
							if (!ret.IsEmpty() && !ret.StartsWith(t_("Parsing error: "))) 
								ret = AQWABody::LoadDatANSYSTOAQWA(mesh, hy, file);	
							y0z = hy.dt.symX;
							x0z = hy.dt.symY;
						}
					}
				}
			}
		}
	} else if (ext == ".lis")
		ret = AQWABody::LoadLis(mesh, file, g, y0z, x0z);
	else if (ext == ".yml") {
		OrcaWave orca;
		ret = orca.Load(file);
		if (ret.IsEmpty()) {
			mesh.Append(orca.dt.msh);
			y0z = orca.dt.symX;
			x0z = orca.dt.symY;
		}
	}
#ifdef PLATFORM_WIN32	
	else if (ext == ".owr")
		ret = ORCABody::Load_OWR(mesh, file, g, y0z, x0z);
#endif	
	else if (ext == ".txt") 
		ret = DiodoreBody::LoadDat(mesh, file); 
	else if (ext == ".gdf") 
		ret = WamitBody::LoadGdf(mesh, file, y0z, x0z); 
	else if (ext == ".pnl") 
		ret = HAMSBody::LoadPnl(mesh, file, y0z, x0z); 
	else if (ext == ".hst") 
		ret = HydrostarBody::LoadHst(mesh, file, y0z, x0z); 
	else if (ext == ".nc") 
		ret = CapyBody::Load_NC(mesh, file, g);
	else if (ext == ".stl") {
		bool isText;
		Body &m = mesh.Add();
		try {
			LoadStl(file, m.dt.mesh, isText, m.dt.fileHeader);
		} catch(Exc e) {
			return std::move(e);
		}
		m.dt.SetCode(isText ? Body::STL_TXT : Body::STL_BIN);
	} else if (ext == ".msh") {
		Body &m = mesh.Add();
		try {
			LoadTDynMsh(file, m.dt.mesh);
		} catch(Exc e) {
			if (e.StartsWith(t_("Parsing error: "))) {
				try {
					LoadGMSH(file, m.dt.mesh);
				} catch(Exc e) {
					return std::move(e);
				}	
			} else
				return std::move(e);
		}
		m.dt.SetCode(Body::MSH_TDYN);
	} else if (ext == ".grd") {
		Body &m = mesh.Add();
		try {
			LoadGRD(file, m.dt.mesh, y0z, x0z);
		} catch(Exc e) {
			return std::move(e);
		}
		m.dt.SetCode(Body::MIKE21_GRD);
	} else if (ext == ".obj") {
		Body &m = mesh.Add();
		try {
			LoadOBJ(file, m.dt.mesh);
		} catch(Exc e) {
			return std::move(e);
		}
		m.dt.SetCode(Body::OBJ);
	} else if (ext == ".mesh" || ext == ".bem") {
		Body &m = mesh.Add();
		try {
			m.dt.mesh.Load(file);
		} catch(Exc e) {
			return std::move(e);
		}
		m.dt.SetCode(Body::BEM_MESH);
	} else
		ret = Format(t_("Unknown mesh file format '%s'"), file);	
	
	if (!ret.IsEmpty())
		return ret;
	
	for (Body &m : mesh) {
		if (IsNull(m.dt.c0))
			m.dt.c0 = Point3D(0, 0, 0);
		
		ret = m.dt.mesh.CheckErrors();
		if (!ret.IsEmpty())
			return ret;
		
		m.dt.fileName = file;
		if (m.dt.name.IsEmpty())
			m.dt.name = InitCaps(GetFileTitle(file));
		
		if (y0z)
			m.dt.mesh.DeployXSymmetry();
		if (x0z)
			m.dt.mesh.DeployYSymmetry();	
		
		if (cleanPanels) 
			m.dt.mesh.Heal(true, grid, eps);
		
		m.AfterLoad(rho, g, false, true);
	
	}
	return String();
}

void Body::SaveAs(const UArray<Body> &meshes, String file, MESH_FMT type, MESH_TYPE meshType, double rho, double g, bool symX, bool symY, 
				int &nNodes, int &nPanels, const UVector<double> &w, const UVector<double> &head, bool getQTF, bool getPotentials, double h, int numCores) {
	UArray<Surface> surfs(meshes.size());
	nNodes = nPanels = 0;
	
	for (int i = 0; i < meshes.size(); ++i) {
		Surface &surf = surfs[i];
		if (meshType == UNDERWATER) 
			surf = clone(meshes[i].dt.under);
		else {
			if (type == AQWA_DAT) {		// Appends dry and wet sides. This way there are no panels between dry and wet side
				surf = clone(meshes[i].dt.under);	// First the wet
				Surface dry;	
				dry.CutZ(meshes[i].dt.mesh, 1);
				surf.Append(dry);					// Next the dry
			} else
				surf = clone(meshes[i].dt.mesh);
		}
		
		if (symX && (type == WAMIT_GDF || type == HAMS_PNL || type == DIODORE_DAT || type == AQWA_DAT || type == MIKE21_GRD)) {
			Surface nsurf;
			nsurf.CutX(surf);
			surf = pick(nsurf);
		}
		if (symY && (type == WAMIT_GDF || type == NEMOH_DAT || type == NEMOH_PRE || 
					 type == HAMS_PNL || type == DIODORE_DAT || type == AQWA_DAT || type == MIKE21_GRD)) {
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
			throw Exc(t_("Impossible to save mesh. No mesh found"));

		nNodes += surf.nodes.size();
		nPanels += surf.panels.size();
	}
	
	if (type == WAMIT_GDF) 
		WamitBody::SaveGdf(file, First(surfs), g, symX, symY);	
	else if (type == NEMOH_DAT) 
		NemohBody::SaveDat(meshes, file, First(surfs), symY, nPanels);
	else if (type == NEMOH_PRE) 
		NemohBody::SavePreBody(file, First(surfs));
	else if (type == HAMS_PNL)		
		HAMSBody::SavePnl(file, First(surfs), symX, symY);	// Only one symmetry is really available
	else if (type == AQWA_DAT)		
		AQWABody::SaveDat(file, meshes, surfs, rho, g, symX, symY, w, head, getQTF, getPotentials, h, numCores);	
	else if (type == DIODORE_DAT) 
		DiodoreBody::SaveDat(file, First(surfs));
	else if (type == STL_BIN)		
		SaveStlBin(file, First(surfs));
	else if (type == STL_TXT)		
		SaveStlTxt(file, First(surfs));
	else if (type == BEM_MESH)		
		First(surfs).Save(file);
	else if (type == MIKE21_GRD)		
		SaveGRD(file, First(surfs), g, symX, symY);
	else
		throw Exc(t_("Unknown mesh file type"));
}

String Body::Heal(bool basic, double rho, double g, double grid, double eps, Function <bool(String, int pos)> Status) {
	String ret = dt.mesh.Heal(basic, grid, eps, Status);
	if (!ret.IsEmpty())
		return ret;
	
	AfterLoad(rho, g, false, false);
	
	return String();
}

void Body::Orient() {
	dt.mesh.Orient();
}

void Body::Append(const Surface &orig, double rho, double g) {
	dt.mesh.Append(orig);
	
	AfterLoad(rho, g, false, true);
}

void Body::Reset(double rho, double g) {
	dt.mesh = clone(dt.mesh0);
	if (!IsNull(dt.cg0))
		dt.cg = clone(dt.cg0);
	cdt.controlPointsA = clone(cdt.controlPointsA0);
	cdt.controlPointsB = clone(cdt.controlPointsB0);
	cdt.controlPointsC = clone(cdt.controlPointsC0);
	cdt.controlLoads = clone(cdt.controlLoads0);
	
	for (Body *b : cdt.damagedBodies)
		if (b->IsValid())
			b->Reset(rho, g);
		
	AfterLoad(rho, g, false, false);
}
	
void Body::AfterLoad(double rho, double g, bool onlyCG, bool isFirstTime, bool massBuoy, bool reZero) {
	if (dt.mesh.IsEmpty()) {
		if (!dt.mesh.lines.IsEmpty())
			dt.mesh.GetEnvelope();
		return;
	}
	if (!onlyCG) {
		dt.mesh.GetPanelParams();
		dt.mesh.GetEnvelope();
		dt.mesh.GetArea();
		dt.mesh.GetVolume();
		
		dt.under.CutZ(dt.mesh, -1);
		dt.under.GetPanelParams();
		dt.projectionPos.x = dt.under.GetAreaXProjection(true, false);
		dt.projectionNeg.x = dt.under.GetAreaXProjection(false, true);
		dt.projectionPos.y = dt.under.GetAreaYProjection(true, false);
		dt.projectionNeg.y = dt.under.GetAreaYProjection(false, true);
		dt.projectionPos.z = dt.under.GetAreaZProjection(true, false);
		dt.projectionNeg.z = dt.under.GetAreaZProjection(false, true);
		dt.cgZ0surface = dt.under.GetAreaZProjectionCG();
		dt.under.GetArea();
		dt.under.GetVolume();
		
		if (dt.M.size() != 36)
			dt.M = MatrixXd::Zero(6,6);
		if (GetMass() == 0 && !IsNull(rho))
			SetMass(dt.under.volume*rho);
		dt.cb = dt.under.GetCentreOfBuoyancy();
		if (IsNull(dt.Vo))
			dt.Vo = dt.under.volume;
	}
	if (isFirstTime) 
		IncrementIdCount();
	
	if (isFirstTime || reZero) {
		dt.mesh0 = clone(dt.mesh);
		dt.cg0 = clone(dt.cg);
		cdt.controlPointsA0 = clone(cdt.controlPointsA);
		cdt.controlPointsB0 = clone(cdt.controlPointsB);
		cdt.controlPointsC0 = clone(cdt.controlPointsC);
		cdt.controlLoads0 = clone(cdt.controlLoads);
	}
	if (!onlyCG && !IsNull(rho) && !IsNull(g) && !IsNull(dt.cg) && !IsNull(dt.cb))
		dt.under.GetHydrostaticStiffness(dt.C, dt.c0, dt.cg, dt.cb, rho, g, GetMass(), massBuoy);
	
	for (Body *b : cdt.damagedBodies)
		if (b->IsValid())
			b->AfterLoad(rho, g, onlyCG,isFirstTime, massBuoy, reZero);
}

void Body::Report(double rho) const {
	BEM::Print("\n\n" + Format(t_("Body file '%s'"), dt.fileName));
	
	BEM::Print(S("\n") + Format(t_("Limits [m] (%f - %f, %f - %f, %f - %f)"), 
			dt.mesh.env.minX, dt.mesh.env.maxX, dt.mesh.env.minY, dt.mesh.env.maxY, dt.mesh.env.minZ, dt.mesh.env.maxZ));
	BEM::Print(S("\n") + Format(t_("Water-plane area. Surface projection Z-axis [m2] %f - %f = %f"), -dt.projectionPos.z, dt.projectionNeg.z, dt.projectionPos.z + dt.projectionNeg.z));
	BEM::Print(S("\n") + Format(t_("Surface projection X-axis [m2] %f - %f = %f"), -dt.projectionPos.x, dt.projectionNeg.x, dt.projectionPos.x + dt.projectionNeg.x));
	BEM::Print(S("\n") + Format(t_("Surface projection Y-axis [m2] %f - %f = %f"), -dt.projectionPos.y, dt.projectionNeg.y, dt.projectionPos.y + dt.projectionNeg.y));
	BEM::Print(S("\n") + Format(t_("Surface [m2] %f"), dt.mesh.surface));
	BEM::Print(S("\n") + Format(t_("Volume [m³] %f"), dt.mesh.volume));
	BEM::Print(S("\n") + Format(t_("Underwater surface [m2] %f"), dt.under.surface));
	BEM::Print(S("\n") + Format(t_("Underwater volume [m³] %f"), dt.under.volume));
	BEM::Print(S("\n") + Format(t_("Displacement [tm] %f"), dt.under.volume*rho/1000));
	BEM::Print(S("\n") + Format(t_("Centre of buoyancy [m] (%f, %f, %f)"), dt.cb.x, dt.cb.y, dt.cb.z));
	
	BEM::Print(S("\n") + Format(t_("Loaded %d panels and %d nodes"), dt.mesh.panels.size(), dt.mesh.nodes.size()));
}

bool Body::IsSymmetricX() {
	return abs(dt.cb.x)/abs(dt.mesh.env.maxX - dt.mesh.env.minX) < 0.001 && abs(dt.mesh.env.maxX + dt.mesh.env.minX) < 0.001;
}

bool Body::IsSymmetricY() {
	return abs(dt.cb.y)/abs(dt.mesh.env.maxY - dt.mesh.env.minY) < 0.001 && abs(dt.mesh.env.maxY + dt.mesh.env.minY) < 0.001;
}

double Body::GMroll(double rho, double g) const {
	if (!IsNum(rho) || !IsNum(g) || dt.under.volume == 0 || dt.C.size() == 0)
		return Null;
	return dt.C(3, 3)/(rho*g*dt.under.volume);
}

double Body::GMpitch(double rho, double g) const {
	if (!IsNum(rho) || !IsNum(g) || dt.under.volume == 0 || dt.C.size() == 0)
		return Null;
	return dt.C(4, 4)/(rho*g*dt.under.volume);
}
		
void Body::GZ(double from, double to, double delta, double angleCalc, double rho, double g,
  double tolerance, UVector<double> &dataangle, UVector<double> &datagz, String &error) {
	UVector<double> dataMoment, vol, disp, wett, wplane, draft;
	UVector<Point3D> dcb, dcg;
	UArray<UVector<double>> zA, zB, zC;
	GZ(from, to, delta, angleCalc, rho, g, tolerance, Null, dataangle, datagz, dataMoment, 
		vol, disp, wett, wplane, draft, dcb, dcg, error, zA, zB, zC);
}

void Body::GZ(double from, double to, double delta, double angleCalc, double rho, double g, double tolerance,
	Function <bool(String, int pos)> Status, 
	UVector<double> &dataangle, UVector<double> &datagz, UVector<double> &dataMoment,
	UVector<double> &vol, UVector<double> &disp, UVector<double> &wett, UVector<double> &wplane,
	UVector<double> &draft, UVector<Point3D> &dcb, UVector<Point3D> &dcg, String &error, 
	UArray<UVector<double>> &zA, UArray<UVector<double>> &zB, UArray<UVector<double>> &zC) {
	
	if (GetMass() == 0)
		throw Exc(t_("Problem obtaining GZ. Mass is zero"));
	
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
	zA.Clear();
	zB.Clear();
	zC.Clear();
	
	Body base0 = clone(*this);

	base0.Rotate(0, 0, ToRad(angleCalc), dt.c0.x, dt.c0.y, dt.c0.z);
	
	zA.SetCount(cdt.controlPointsA.size());
	zB.SetCount(cdt.controlPointsB.size());
	zC.SetCount(cdt.controlPointsC.size());
	
	double dz = 0.1;
	for (double angle = from; angle <= to; angle += delta) {
		double progress =  to == from ? 1 : (angle - from)/(to - from);
		if (Status && !Status("", int(100*progress)))
			throw Exc(t_("Cancelled by the user"));
				
		Body base = clone(base0);
		UVector<Body> damaged;	// Handles a copy of the damaged bodies
		base.cloneDamaged(damaged);
			
		base.Rotate(0, ToRad(angle), 0, dt.c0.x, dt.c0.y, dt.c0.z);
		
		Point3D ccb;
		double allvol;
		if (!base.TranslateArchimede(rho, tolerance, dz, ccb, allvol)) {
			//throw Exc(t_("Problem obtaining GZ"));
		
		//if (uunder.VolumeMatch(tolerance, tolerance) < 0) {
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
			dcg << base.dt.cg;
			for (auto &z : zA)
				z << Null;
			for (auto &z : zB)
				z << Null;
			for (auto &z : zC)
				z << Null;
		} else {
			//Point3D ccb = uunder.GetCentreOfBuoyancy();
			Force6D fcb = Force6D::Zero();//Surface::GetHydrostaticForceCB(dt.c0, ccb, uunder.volume, rho, g);
			fcb.Add(Vector3D(0, 0, allvol*rho*g), ccb, dt.c0);
			
			Force6D fcg = Surface::GetMassForce(dt.c0, base.dt.cg, GetMass(), g);
			for (const auto &d : base.cdt.controlLoads)
				fcg += Surface::GetMassForce(dt.c0, d.p, d.mass, g);
			
			double moment = -(fcg.r.y + fcb.r.y);
			double gz = moment/GetMass_all()/g;
			
			dataangle << angle;
			datagz << gz;
			dataMoment << moment;
			vol << allvol;
			disp << allvol*rho;
			wett << base.dt.under.GetArea();
			wplane << base.dt.under.GetWaterPlaneArea();
			draft << base.dt.under.GetEnvelope().minZ;
			dcb << ccb;
			dcg << base.GetCG_all();
			for (int i = 0; i < zA.size(); ++i)
				zA[i] << base.cdt.controlPointsA[i].p.z;
			for (int i = 0; i < zB.size(); ++i)
				zB[i] << base.cdt.controlPointsB[i].p.z;
			for (int i = 0; i < zC.size(); ++i)
				zC[i] << base.cdt.controlPointsC[i].p.z;
		}
	}	
}

void Body::Move(double dx, double dy, double dz, double ax, double ay, double az, double rho, double g, bool setnewzero) {
	dt.mesh = clone(dt.mesh0);
	if (!IsNull(dt.cg0))
		dt.cg = clone(dt.cg0);					
	dt.mesh.TransRot(dx, dy, dz, ax, ay, az, dt.c0.x, dt.c0.y, dt.c0.z);
	if (!IsNull(dt.cg))
		dt.cg.TransRot(dx, dy, dz, ax, ay, az, dt.c0.x, dt.c0.y, dt.c0.z);
	
	for (auto &d : cdt.controlPointsA)
		d.p.TransRot(dx, dy, dz, ax, ay, az, dt.c0.x, dt.c0.y, dt.c0.z);
	for (auto &d : cdt.controlPointsB)
		d.p.TransRot(dx, dy, dz, ax, ay, az, dt.c0.x, dt.c0.y, dt.c0.z);
	for (auto &d : cdt.controlPointsC)
		d.p.TransRot(dx, dy, dz, ax, ay, az, dt.c0.x, dt.c0.y, dt.c0.z);
	for (auto &d : cdt.controlLoads)
		d.p.TransRot(dx, dy, dz, ax, ay, az, dt.c0.x, dt.c0.y, dt.c0.z);
	
	AfterLoad(rho, g, false, setnewzero);	
}

void Body::Move(const double *pos, double rho, double g, bool setnewzero) {
	Move(pos[0], pos[1], pos[2], pos[3], pos[4], pos[5], rho, g, setnewzero);	
}

void Body::Move(const float *pos, double rho, double g, bool setnewzero) {
	Move(pos[0], pos[1], pos[2], pos[3], pos[4], pos[5], rho, g, setnewzero);	
}

void Body::SetMass(double m) {
	if (dt.M.size() != 36)
		dt.M = MatrixXd::Zero(6,6);	
	double oldm = dt.M(0,0);
	if (oldm != 0) {
		double factor = m/oldm;
		dt.M.array() *= factor;
	} else
		dt.M(0,0) = dt.M(1,1) = dt.M(2,2) = m;
}

void Body::Translate(double dx, double dy, double dz) {
	dt.mesh.Translate(dx, dy, dz);
	if (!IsNull(dt.cg))
		dt.cg.Translate(dx, dy, dz);
	for (auto &d : cdt.controlPointsA)
		d.p.Translate(dx, dy, dz);
	for (auto &d : cdt.controlPointsB)
		d.p.Translate(dx, dy, dz);
	for (auto &d : cdt.controlPointsC)
		d.p.Translate(dx, dy, dz);
	for (auto &d : cdt.controlLoads)
		d.p.Translate(dx, dy, dz);
	for (Body *b : cdt.damagedBodies)
		if (b->IsValid())
			b->Translate(dx, dy, dz);
}

void Body::Rotate(double a_x, double a_y, double a_z, double c_x, double c_y, double c_z) {
	dt.mesh.Rotate(a_x, a_y, a_z, c_x, c_y, c_z);
	if (!IsNull(dt.cg))
		dt.cg.Rotate(a_x, a_y, a_z, c_x, c_y, c_z);
	for (auto &d : cdt.controlPointsA)
		d.p.Rotate(a_x, a_y, a_z, c_x, c_y, c_z);
	for (auto &d : cdt.controlPointsB)
		d.p.Rotate(a_x, a_y, a_z, c_x, c_y, c_z);
	for (auto &d : cdt.controlPointsC)
		d.p.Rotate(a_x, a_y, a_z, c_x, c_y, c_z);
	for (auto &d : cdt.controlLoads)
		d.p.Rotate(a_x, a_y, a_z, c_x, c_y, c_z);
	for (Body *b : cdt.damagedBodies)
		if (b->IsValid())
			b->Rotate(a_x, a_y, a_z, c_x, c_y, c_z);
}	

bool Body::TranslateArchimede(double rho, double tolerance, double &dz) {
	Point3D cb;
	double allvol;
	return TranslateArchimede(rho, tolerance, dz, cb, allvol);
}

bool Body::TranslateArchimede(double rho, double tolerance, double &dz, Point3D &cb, double &allvol) {
	UVector<Surface *> damaged;
	for (Body *b : cdt.damagedBodies)
		if (b->IsValid())
			damaged << &(b->dt.mesh);
	
	if (!dt.mesh.TranslateArchimede(GetMass_all(), rho, Bem().volError/100., damaged, tolerance, dz, cb, allvol))
		return false;
	if (!IsNull(dt.cg))
		dt.cg.Translate(0, 0, dz);
	for (auto &d : cdt.controlPointsA)
		d.p.Translate(0, 0, dz);
	for (auto &d : cdt.controlPointsB)
		d.p.Translate(0, 0, dz);
	for (auto &d : cdt.controlPointsC)
		d.p.Translate(0, 0, dz);
	for (auto &d : cdt.controlLoads)
		d.p.Translate(0, 0, dz);
	for (auto &b : cdt.damagedBodies)
		if (b->IsValid())
			b->Translate(0, 0, dz);
	return true;
}

void Body::PCA(double &yaw) {
	//Value3D ax1, ax2, ax3;
	//dt.mesh.PrincipalComponents(ax1, ax2, ax3);
	//yaw = -atan(abs(ax1.y/ax1.x));
	
	yaw = dt.mesh.YawMainAxis();
	
	Rotate(0, 0, yaw, dt.c0.x, dt.c0.y, dt.c0.z);		// Adjusts the yaw error
}
	
bool Body::Archimede(double rho, double g, double tolerance, double &roll, double &pitch, double &dz) {
	Point3D cb_dmg;
	double resroll, respitch;
	Force6D fcg, fcb;
	double vol_dmg;
	
	auto Residual = [&](double roll, double pitch, double &resroll, double &respitch) {
		Body base = clone(*this);		// Handles a copy of the body
		UVector<Body> damaged;	// Handles a copy of the damaged bodies
		base.cloneDamaged(damaged);
		
		base.Rotate(roll, pitch, 0, dt.c0.x, dt.c0.y, dt.c0.z);
		base.TranslateArchimede(rho, tolerance, dz, cb_dmg, vol_dmg);	// All volumes are included
		
		UVector<Point3D> cgs;											// All loads are included
		UVector<double> masses;
		cgs << base.dt.cg;
		masses << base.GetMass();
		for (ControlData::ControlLoad &c : cdt.controlLoads) {
			if (c.loaded) {
				cgs << c.p;
				masses << c.mass;
			}
		}
		fcg = Surface::GetMassForce(base.dt.c0, cgs, masses, g);
		fcb = Surface::GetHydrostaticForceCB(base.dt.c0, cb_dmg, vol_dmg, rho, g);
	
		resroll  = fcb.r.x + fcg.r.x;		// ∑ Froll = 0	
		respitch = fcb.r.y + fcg.r.y;		// ∑ Fpitch = 0	
	};
	
	// First, vessel is put into the water
	roll = pitch = 0;
	Residual(roll, pitch, resroll, respitch);
	
	// First delta is a fraction of the angle between cg and cb
	Value3D delta_cbcg = GetCG_all() - cb_dmg;
	
	double droll  = 0.25*abs(M_PI/2 - abs(atan2(delta_cbcg.z, delta_cbcg.y)));
	if (Sign(droll) != Sign(resroll))
		droll = -droll;
	double dpitch = 0.25*abs(M_PI/2 - abs(atan2(delta_cbcg.z, delta_cbcg.x)));
	if (Sign(dpitch) != Sign(respitch))
		dpitch = -dpitch;
	
	// Iterate rotating around cg	
	int nIter = 0;	
	int maxIter = 500;
	for (; nIter < maxIter; ++nIter) {
		if (abs(droll) < M_PI/1000000 && abs(dpitch) < M_PI/1000000) 		// > 0.00018 deg
			break;
		
		roll += droll;
		pitch += dpitch;
		
		LOG(Format("Archimede roll: %f %f %.0f pitch: %f %f %.0f", roll, droll, resroll, pitch, dpitch, respitch));
		
		double nresroll, nrespitch;
		Residual(roll, pitch, nresroll, nrespitch);
		
		if (Sign(resroll) != Sign(nresroll)) 	// We have crossed a zero
			droll /= 2;								// Change the sign reducing the delta
		resroll = nresroll;
		if (Sign(droll) != Sign(resroll))
			droll = -droll;
		
		if (Sign(respitch) != Sign(nrespitch)) 
			dpitch /= 2;				
		respitch = nrespitch;
		if (Sign(dpitch) != Sign(respitch))
			dpitch = -dpitch;
	}
		
	LOG(Format("Archimede NumIter: %d", nIter));
	
	if (nIter >= maxIter)
		return false;
	
	Rotate(roll, pitch, 0, dt.c0.x, dt.c0.y, dt.c0.z);	// Moves just the necessary
	TranslateArchimede(rho, tolerance, dz);
	Value3D delta = dt.cg0 - dt.cg;
	Translate(delta.x, delta.y, 0);
	
	return true;
}

double Body::GetMass_all() const	{
	double mass = GetMass();
	for (const auto &d : cdt.controlLoads)
		if (d.loaded && !IsNull(d.mass)) 
			mass += d.mass;
			
	return mass;
}

Point3D Body::GetCG_all() const {
	double mass = GetMass();
	Point3D cg = dt.cg*mass;
	for (const auto &d : cdt.controlLoads) {
		if (d.loaded && !IsNull(d.p) && !IsNull(d.mass)) {
			cg.x += d.mass*d.p.x;
			cg.y += d.mass*d.p.y;
			cg.z += d.mass*d.p.z;
			mass += d.mass;
		}
	}
	cg /= mass;	

	return cg;		
}

Point3D Body::GetCB_all() const {
	double allvol = dt.under.volume;
	Point3D cb = dt.under.GetCentreOfBuoyancy()*dt.under.volume;
	for (Body *pb : cdt.damagedBodies) {
		const Body &b = *pb;
		if (b.dt.under.VolumeMatch(Bem().volError, Bem().volError) < 0)
			return Null;
		double vv = b.dt.under.volume;
		if (vv > 0) {
			Point3D ccb = b.dt.under.GetCentreOfBuoyancy();
			cb -= ccb*vv;
			allvol -= vv;
		}
	}
	if (allvol <= 0)
		return Null;
	
	cb /= allvol;
	return cb;
}
	
void Body::Jsonize(JsonIO &json) {
	json
		("projectionPos", dt.projectionPos)
		("projectionNeg", dt.projectionNeg)
		("cgZ0surface", dt.cgZ0surface)
		("cb", dt.cb)
		("cg", dt.cg)
		("cg0", dt.cg0)
		("c0", dt.c0)
		("cb", dt.cb)
		("Vo", dt.Vo)
		("M", dt.M)
		("C", dt.C)
		("Cmoor", dt.Cmoor)
		("Cadd", dt.Cadd)
		("Dlin", dt.Dlin)
		("Dquad", dt.Dquad)
		("Aadd", dt.Aadd)
		("name", dt.name)
		("fileName", dt.fileName)
		("fileHeader", dt.fileHeader)
		("mesh", dt.mesh)
		("under", dt.under)
		("mesh0", dt.mesh0)
		("ControlData", cdt)
	;
}

	