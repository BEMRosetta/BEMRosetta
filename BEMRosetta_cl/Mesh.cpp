// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"

	
// enum MESH_FMT 			    	  {WAMIT_GDF,  WAMIT_DAT,   NEMOH_DAT,   NEMOHFS_DAT,   NEMOH_PRE,      AQWA_DAT,   AQWA LIS, HAMS_PNL,  STL_BIN,     STL_TXT,   EDIT,  MSH_TDYN,   BEM_MESH, DIODORE_DAT,   HYDROSTAR_HST,    ORCA_OWR, 	   MIKE21_GRD	UNKNOWN, NUMMESH};	
const char *Body::meshStr[]         = {"Wamit.gdf","Wamit.dat",	"Nemoh.dat", "NemohFS.dat", "Nemoh premesh","AQWA.dat", "AQWA.lis","HAMS.pnl","STL.Binary","STL.Text","Edit","TDyn.msh", "BEMR",   "Diodore.dat", "HydroStar.hst", "OrcaWave.owr", "MIKE21.grd", "Unknown"};	
const bool Body::meshCanSave[] 		= {true, 	   false,	    true,		 false,			false, 		    true,		false,	   true,	   true,		true,	   false, false, 	  true, 	true,		   false,   	   false, 		   true, 		 false};       
const char *Body::meshExt[]	  		= {"*.gdf",    "*.dat",	 	"*.dat",	 "*.dat", 		"",		        "*.dat",	"*.lis",   "*.pnl",   "*.stl",     "*.stl",    "",	  "*.msh",   "*.bemr",  "*.dat", 	  "*.hst", 	   	   "*.owr",		   "*.grd", 	 "*.*"};       

int Body::idCount = 0;


void Body::Copy(const Body &msh) {
	dt.projectionPos = clone(msh.dt.projectionPos);
	dt.projectionNeg = clone(msh.dt.projectionNeg);
	
	dt.cgZ0surface = clone(msh.dt.cgZ0surface);
	dt.cb = clone(msh.dt.cb);
	dt.cg = clone(msh.dt.cg);
	dt.cg0 = clone(msh.dt.cg0);
	dt.c0 = clone(msh.dt.c0);
	dt.Vo = msh.dt.Vo;
	
	dt.M = clone(msh.dt.M);
	dt.C = clone(msh.dt.C);
	dt.Cmoor = clone(msh.dt.Cmoor);
	dt.Cadd = clone(msh.dt.Cadd);
	//Dlin = clone(msh.Dlin);
	//Dquad = clone(msh.Dquad);
	dt.Aadd = clone(msh.dt.Aadd);
	
	dt.name = msh.dt.name;
	dt.fileName = msh.dt.fileName;
	dt.fileHeader = msh.dt.fileHeader;
	dt.lidFile = msh.dt.lidFile;
	
	dt.mesh = clone(msh.dt.mesh);
	dt.under = clone(msh.dt.under);
	dt.mesh0 = clone(msh.dt.mesh0);
	
	dt.SetCode(msh.dt.GetCode());
	dt.SetId(msh.dt.GetId());
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
						if (!ret.IsEmpty() && !ret.StartsWith(t_("Parsing error: "))) 	
							ret = AQWABody::LoadDat(mesh, file, y0z, x0z);
					}
				}
			}
		}
	} else if (ext == ".lis")
		ret = AQWABody::Load_LIS(mesh, file, g, y0z, x0z);
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
		
		if (!IsNull(rho))
			m.AfterLoad(rho, g, false, true);
		
		m.IncrementIdCount();
	}
	return String();
}

void Body::SaveAs(const UArray<Body> &meshes, String file, MESH_FMT type, MESH_TYPE meshType, double rho, double g, bool symX, bool symY, 
						int &nNodes, int &nPanels) {
	UArray<Surface> surfs(meshes.size());
	nNodes = nPanels = 0;
	
	for (int i = 0; i < meshes.size(); ++i) {
		Surface &surf = surfs[i];
		if (meshType == UNDERWATER) 
			surf = clone(meshes[i].dt.under);
		else {
			if (type == AQWA_DAT) {		// Appends dry and wet sides. This way there are no panels between dry and wet side
				surf = clone(meshes[i].dt.under);		// First the wet
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
			throw Exc(t_("Model is empty. No data found"));

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
		AQWABody::SaveDat(file, meshes, surfs, rho, g, symX, symY);	
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
	dt.cg = clone(dt.cg0);
	AfterLoad(rho, g, false, false);
}
	
void Body::AfterLoad(double rho, double g, bool onlyCG, bool isFirstTime, bool massBuoy) {
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
	}
	if (isFirstTime) {
		dt.mesh0 = clone(dt.mesh);
		dt.cg0 = clone(dt.cg);
	}
	if (!IsNull(rho) && !IsNull(g) && !IsNull(dt.cg) && !IsNull(dt.cb))
		dt.under.GetHydrostaticStiffness(dt.C, dt.c0, dt.cg, dt.cb, rho, g, GetMass(), massBuoy);
}

void Body::Report(double rho) const {
	BEM::Print("\n\n" + Format(t_("Body file '%s'"), dt.fileName));
	
	BEM::Print(S("\n") + Format(t_("Limits [m] (%f - %f, %f - %f, %f - %f)"), 
			dt.mesh.env.minX, dt.mesh.env.maxX, dt.mesh.env.minY, dt.mesh.env.maxY, dt.mesh.env.minZ, dt.mesh.env.maxZ));
	BEM::Print(S("\n") + Format(t_("Water-plane area. Surface projection Z-axis [m2] %f - %f = %f"), -dt.projectionPos.z, dt.projectionNeg.z, dt.projectionPos.z + dt.projectionNeg.z));
	BEM::Print(S("\n") + Format(t_("Surface projection X-axis [m2] %f - %f = %f"), -dt.projectionPos.x, dt.projectionNeg.x, dt.projectionPos.x + dt.projectionNeg.x));
	BEM::Print(S("\n") + Format(t_("Surface projection Y-axis [m2] %f - %f = %f"), -dt.projectionPos.y, dt.projectionNeg.y, dt.projectionPos.y + dt.projectionNeg.y));
	BEM::Print(S("\n") + Format(t_("Surface [m2] %f"), dt.mesh.surface));
	BEM::Print(S("\n") + Format(t_("Volume [m3] %f"), dt.mesh.volume));
	BEM::Print(S("\n") + Format(t_("Underwater surface [m2] %f"), dt.under.surface));
	BEM::Print(S("\n") + Format(t_("Underwater volume [m3] %f"), dt.under.volume));
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
	GZ(from, to, delta, angleCalc, rho, g, tolerance, Null, dataangle, datagz, dataMoment, 
		vol, disp, wett, wplane, draft, dcb, dcg, error);
}

void Body::GZ(double from, double to, double delta, double angleCalc, double rho, double g, double tolerance,
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
	
	Surface base0 = clone(dt.mesh);
	Point3D ccg0 = clone(dt.cg);

	base0.Rotate(0, 0, ToRad(angleCalc), dt.c0.x, dt.c0.y, dt.c0.z);
	ccg0.Rotate(0, 0, ToRad(angleCalc), dt.c0.x, dt.c0.y, dt.c0.z);
	
	double dz = 0.1;
	for (double angle = from; angle <= to; angle += delta) {
		double progress =  to == from ? 1 : (angle - from)/(to - from);
		if (Status && !Status("", int(100*progress)))
			throw Exc(t_("Cancelled by the user"));
				
		Surface base = clone(base0);
		Point3D ccg = clone(ccg0);
		
		base.Rotate(0, ToRad(angle), 0, dt.c0.x, dt.c0.y, dt.c0.z);
		ccg.Rotate(0, ToRad(angle), 0, dt.c0.x, dt.c0.y, dt.c0.z);
		
		if (GetMass() == 0)
			throw Exc(t_("Problem obtaining GZ. Mass is zero"));
		
		Surface uunder;
		if (!base.TranslateArchimede(GetMass(), rho, dz, uunder))
			throw Exc(t_("Problem obtaining GZ"));
		
		ccg.Translate(0, 0, dz);
		
		if (uunder.VolumeMatch(tolerance, tolerance) < 0) {
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
			dcg << ccg;
			
		} else {
			Point3D ccb = uunder.GetCentreOfBuoyancy();
			Force6D fcb = uunder.GetHydrostaticForceCB(dt.c0, ccb, rho, g);
			//base.GetPanelParams();	
			//Force6D fcb = base.GetHydrodynamicForce(c0, true, 
			//				[&](double x, double y)->double {return 0;}, 
			//				[&](double x, double y, double z, double et)->double {return z > 0 ? 0 : rho*g*z;});
			Force6D fcg = Surface::GetMassForce(dt.c0, ccg, GetMass(), g);
		
			double moment = -(fcg.r.y + fcb.r.y);
			double gz = moment/GetMass()/g;
			
			dataangle << angle;
			datagz << gz;
			dataMoment << moment;
			vol << uunder.volume;
			disp << uunder.volume*rho;
			wett << uunder.GetArea();
			wplane << uunder.GetWaterPlaneArea();
			draft << uunder.GetEnvelope().minZ;
			dcb << ccb;
			dcg << ccg;
		}
	}	
}

void Body::Move(double dx, double dy, double dz, double ax, double ay, double az, double rho, double g, bool setnewzero) {
	dt.mesh = clone(dt.mesh0);
	dt.cg = clone(dt.cg0);					
	dt.mesh.TransRot(dx, dy, dz, ax, ay, az, dt.c0.x, dt.c0.y, dt.c0.z);
	dt.cg.TransRot(dx, dy, dz, ax, ay, az, dt.c0.x, dt.c0.y, dt.c0.z);
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
		("lidFile", dt.lidFile)
		("mesh", dt.mesh)
		("under", dt.under)
		("mesh0", dt.mesh0)
	;
}

	