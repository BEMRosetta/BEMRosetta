// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2026, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include <SysInfo/Crash.h>

BEM &Bem() 		{static BEM bem;		return bem;}

#ifdef PLATFORM_WIN32
#include "orca.h"


Function<bool(String)> Orca::WhenPrint = [](String str)->bool {
	BEM::Print(F("\n%s", str)); 
	return 0;
};

Time Orca::startCalc = Null, Orca::beginNoLicense = Null, Orca::lastLog;
int64 Orca::noLicenseTime = 0;

#endif


#include "FastOut.h"
#include "libbemrosetta.h"

BMR_Data &BMR() {
	static BMR_Data dll;
	return dll;
}

const char *_BMR_GetLastError() noexcept {
	if (BMR().errorStr.IsEmpty())
		return nullptr;
	return BMR().errorStr;
}

void _BMR_ClearLastError() noexcept {
	BMR().errorStr.Clear();
}

#ifndef flagBEMR_TEST_DLL


void _BMR_SetErrorHandler(void (*error_callback)(const char *, void *), void *user_data) noexcept {
	BMR().error_callback = error_callback;
	BMR().user_data = user_data;
}

error_callback_t _BMR_GetErrorHandlerCallback() noexcept {
	return BMR().error_callback;
}

void *_BMR_GetErrorHandlerData() noexcept {
	return BMR().user_data;
}

void _BMR_Init() noexcept {
	BMR();
}			

const char *_BMR_Version() noexcept {
	static String version;
	version << __DATE__ << ", " << __TIME__;
	return version;	
}

void _BMR_EnablePrint(bool print) noexcept {
	Bem().print = print;
}

void _BMR_Wamit_V6s_Set(int force) noexcept {
	Bem().opForceV6 = force;
}

void _BMR_AQWA_ShowCalculationDialog_Set(int show) noexcept {
	Bem().opNoWind = show;
}

int _BMR_Mesh_Load(const char *file) noexcept {
	try {
		if (!FileExists(file))
			throw Exc(F(t_("File '%s' not found"), file)); 
		
		_BMR_Mesh_Id_Set(Bem().surfs.size());						
		Bem().LoadBody(file, BMR().echo ? BMR().Status : BMR().NoPrint, false, false, BMR().meshid);
	} catch(Exc err) {
		BMR().errorStr = err;
		return NullInt;
	}
	BMR().errorStr.Clear();
	return BMR().meshid;
}

bool _BMR_Mesh_Report() noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		Bem().surfs[BMR().meshid].Report(Bem().rho);
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

void _BMR_Mesh_Clear() noexcept {
	Bem().surfs.Clear();
	BMR().meshid = -1;
}

bool _BMR_Mesh_Id_Set(int id) noexcept {
	try {
		if (IsNull(id) || id < 0)
			throw Exc(F(t_("Invalid id %d"), id));
		
		if (id >= Bem().surfs.size())
			Bem().surfs.SetCount(id+1);
		
		BMR().meshid = id;
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

int _BMR_Mesh_Id_Get() noexcept {
	return BMR().meshid;
}

bool _BMR_Mesh_Name_Set(const char *name) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		
		Bem().surfs[BMR().meshid].dt.name = name;
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

const char *_BMR_Mesh_Name_Get() noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		
	} catch(Exc err) {
		BMR().errorStr = err;
		return nullptr;
	}
	BMR().errorStr.Clear();
	return Bem().surfs[BMR().meshid].dt.name;;
}
		
bool _BMR_Mesh_Save(const char *file, const char *format, int symX, int symY) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		
		Body::MESH_FMT meshFmt = Body::GetCodeBodyStr(format);
		if (!Body::meshInfo[meshFmt].canSave)
			throw Exc(F(t_("Saving format '%s' is not implemented"), format));						
		
		Body::SaveAs(Bem().surfs[BMR().meshid], file, meshFmt, Body::ALL, Bem().rho, Bem().g, symX, symY);					
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

bool _BMR_Mesh_Translate(double x, double y, double z) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		msh.dt.mesh.Translate(x, y, z);
		msh.dt.spline.Translate(Point3D(x, y, z));
		msh.dt.cg.Translate(x, y, z);
		msh.AfterLoad(Bem().rho, Bem().g, false, false);	
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

bool _BMR_Mesh_Rotate(double ax, double ay, double az, double cx, double cy, double cz) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		msh.dt.mesh.Rotate(ToRad(ax), ToRad(ay), ToRad(az), cx, cy, cz);	
		msh.dt.spline.Rotate(Point3D(ToRad(ax), ToRad(ay), ToRad(az)), Point3D(cx, cy, cz));
		msh.dt.cg.Rotate(ToRad(ax), ToRad(ay), ToRad(az), cx, cy, cz);
		msh.AfterLoad(Bem().rho, Bem().g, false, false);	
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

bool _BMR_Mesh_Cg_Set(double x, double y, double z) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		msh.dt.cg = Point3D(x, y, z);
		msh.AfterLoad(Bem().rho, Bem().g, true, false);
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

bool _BMR_Mesh_Cg_Get(double *x, double *y, double *z) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		*x = msh.dt.cg.x;
		*y = msh.dt.cg.y;
		*z = msh.dt.cg.z;
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

bool _BMR_Mesh_C0_Set(double x, double y, double z) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		msh.dt.c0 = Point3D(x, y, z);
		msh.AfterLoad(Bem().rho, Bem().g, true, false);
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

bool _BMR_Mesh_C0_Get(double *x, double *y, double *z) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		*x = msh.dt.c0.x;
		*y = msh.dt.c0.y;
		*z = msh.dt.c0.z;
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

bool _BMR_Mesh_Mass_Set(double mass) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		msh.SetMass(mass);
		msh.AfterLoad(Bem().rho, Bem().g, true, false);
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

bool _BMR_Mesh_Inertia_Set(const double *data, const int dim[2]) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		if (dim[0] != 6 || dim[1] != 6)
			throw Exc(F(t_("Matrix dimensions (%d,%d) are not correct"), dim[0], dim[1]));
		msh.dt.M.resize(6, 6);
		for (int r = 0; r < 6; ++r)
			for (int c = 0; c < 6; ++c)
				msh.dt.M(r, c) = data[r*6 + c];
			
		msh.AfterLoad(Bem().rho, Bem().g, true, false);
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

bool _BMR_Mesh_LinearDamping_Set(const double *data, const int dim[2]) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		if (dim[0] != 6 || dim[1] != 6)
			throw Exc(F(t_("Matrix dimensions (%d,%d) are not correct"), dim[0], dim[1]));
		msh.dt.Dlin.resize(6, 6);
		for (int r = 0; r < 6; ++r)
			for (int c = 0; c < 6; ++c)
				msh.dt.Dlin(r, c) = data[r*6 + c];
			
		msh.AfterLoad(Bem().rho, Bem().g, true, false);
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

bool _BMR_Mesh_MooringStiffness_Set(const double *data, const int dim[2]) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		if (dim[0] != 6 || dim[1] != 6)
			throw Exc(F(t_("Matrix dimensions (%d,%d) are not correct"), dim[0], dim[1]));
		msh.dt.Cmoor.resize(6, 6);
		for (int r = 0; r < 6; ++r)
			for (int c = 0; c < 6; ++c)
				msh.dt.Cmoor(r, c) = data[r*6 + c];
			
		msh.AfterLoad(Bem().rho, Bem().g, true, false);
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

bool _BMR_Mesh_Reset() noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		Bem().surfs[BMR().meshid].Reset(Bem().rho, Bem().g);
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

int _BMR_Mesh_Duplicate() noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		
		Body &msh = Bem().surfs[BMR().meshid];
		_BMR_Mesh_Id_Set(Bem().surfs.size());
		Bem().surfs[BMR().meshid] = clone(msh);
	} catch(Exc err) {
		BMR().errorStr = err;
		return NullInt;
	}
	BMR().errorStr.Clear();
	return BMR().meshid;
}

int _BMR_Mesh_GetWaterPlane() noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		Bem().AddWaterSurface(BMR().meshid, 'e', 1, false);
		BMR().meshid = Bem().surfs.size() - 1;
	} catch(Exc err) {
		BMR().errorStr = err;
		return NullInt;
	}
	BMR().errorStr.Clear();
	return BMR().meshid;
}

int _BMR_Mesh_GetHull() noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		Bem().AddWaterSurface(BMR().meshid, 'r', 1, false);
		BMR().meshid = Bem().surfs.size() - 1;
	} catch(Exc err) {
		BMR().errorStr = err;
		return NullInt;
	}
	BMR().errorStr.Clear();
	return BMR().meshid;
}

int _BMR_Mesh_FillWaterplane(double ratio, int quads) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		if (ratio < 0 || ratio > 100)
			throw Exc(F(t_("Wrong mesh ratio %s"), ratio));
		Bem().AddWaterSurface(BMR().meshid, 'f', ratio, quads);
		BMR().meshid = Bem().surfs.size() - 1;
	} catch(Exc err) {
		BMR().errorStr = err;
		return NullInt;
	}
	BMR().errorStr.Clear();
	return BMR().meshid;
}
		
int _BMR_Mesh_GetControlSurface(double distance, double ratio, int quads, int bottom, int top) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		if (ratio < 0 || ratio > 100)
			throw Exc(F(t_("Wrong mesh ratio %s"), ratio));
		if (!bottom && !top)
			throw Exc(t_("Some surface has to be indicated, bottom/side or top"));
		UVector<int> ids;
		ids << BMR().meshid;
		Bem().GetCS(ids, distance, ratio, quads, bottom, top);
		BMR().meshid = Bem().surfs.size() - 1;
	} catch(Exc err) {
		BMR().errorStr = err;
		return NullInt;
	}
	BMR().errorStr.Clear();
	return BMR().meshid;
}

				
bool _BMR_Mesh_Volume_Get(double *vx, double *vy, double *vz) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		*vx = msh.dt.mesh.volumex;
		*vy = msh.dt.mesh.volumey;
		*vz = msh.dt.mesh.volumez;	
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

bool _BMR_Mesh_UnderwaterVolume_Get(double *vx, double *vy, double *vz) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		*vx = msh.dt.under.volumex;
		*vy = msh.dt.under.volumey;
		*vz = msh.dt.under.volumez;
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

double _BMR_Mesh_Surface_Get() noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
	} catch(Exc err) {
		BMR().errorStr = err;
		return NullDouble;
	}
	BMR().errorStr.Clear();
	return Bem().surfs[BMR().meshid].dt.mesh.surface;
}

double _BMR_Mesh_UnderwaterSurface_Get() noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
	} catch(Exc err) {
		BMR().errorStr = err;
		return NullDouble;
	}
	BMR().errorStr.Clear();
	return Bem().surfs[BMR().meshid].dt.under.surface;
}

bool _BMR_Mesh_Centre_Volume_Get(double *x, double *y, double *z) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		Point3D cg = msh.dt.mesh.GetCentreOfBuoyancy();
		*x = cg.x;
		*y = cg.y;
		*z = cg.z;
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

bool _BMR_Mesh_Centre_Surface_Get(double *x, double *y, double *z) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		Point3D cg = msh.dt.mesh.GetCentreOfGravity_Surface();
		*x = cg.x;
		*y = cg.y;
		*z = cg.z;
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

bool _BMR_Mesh_HydrostaticStiffness_Get(double **data, int dim[2]) noexcept {
	static UVector<double> d;
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		dim[0] = (int)msh.dt.C.rows();
		dim[1] = (int)msh.dt.C.cols();
		CopyRowMajor(msh.dt.C, d);
		*data = d.begin();
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

int _BMR_Mesh_NumPanels_Get() noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		
		BMR().errorStr.Clear();
		return msh.dt.mesh.panels.size();
	} catch(Exc err) {
		BMR().errorStr = err;
		return NullInt;
	}
}

bool _BMR_Mesh_VolumeEnvelope_Get(double *minx, double *maxx, double *miny, double *maxy, double *minz, double *maxz) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No mesh is loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		
		*minx = msh.dt.mesh.env.minX;
		*maxx = msh.dt.mesh.env.maxX;
		*miny = msh.dt.mesh.env.minY;
		*maxy = msh.dt.mesh.env.maxY;
		*minz = msh.dt.mesh.env.minZ;
		*maxz = msh.dt.mesh.env.maxZ;
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

void _BMR_Bem_Clear() noexcept {
	Bem().hydros.Clear();
	BMR().bemid = -1;
}

int _BMR_Bem_New() noexcept {
	try {
		_BMR_Bem_Id_Set(Bem().hydros.size());
	} catch(Exc err) {
		BMR().errorStr = err;
		return NullInt;
	}
	BMR().errorStr.Clear();
	return BMR().bemid;
}

bool _BMR_Bem_Name_Set(const char *name) noexcept {
	try {
		if (Bem().hydros.IsEmpty()) 
			throw Exc(t_("No bem case is loaded"));
		
		Bem().hydros[BMR().bemid].dt.name = name;
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

const char *_BMR_Bem_Name_Get() noexcept {
	try {
		if (Bem().hydros.IsEmpty()) 
			throw Exc(t_("No bem case is loaded"));
		
	} catch(Exc err) {
		BMR().errorStr = err;
		return nullptr;
	}
	BMR().errorStr.Clear();
	return Bem().hydros[BMR().bemid].dt.name;
}

int _BMR_Bem_Load(const char *file) noexcept {
	try {
		Bem().LoadBEM(file);
		_BMR_Bem_Id_Set(Bem().hydros.size()-1);
	} catch(Exc err) {
		BMR().errorStr = err;
		return NullInt;
	}
	BMR().errorStr.Clear();
	return BMR().bemid;
}

bool _BMR_Bem_Save(const char *file) noexcept {
	try {
		Bem().hydros[BMR().bemid].SaveAs(file, Null, Hydro::UNKNOWN, Null);
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}
							
bool _BMR_Bem_Id_Set(int id) noexcept {
	try {
		if (IsNull(id) || id < 0)
			throw Exc(F(t_("Invalid id %d"), id));
		
		if (id >= Bem().hydros.size())
			Bem().hydros.SetCount(id+1);
		
		BMR().bemid = id;
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

int _BMR_Bem_Id_Get() noexcept {
	return BMR().bemid;
}

int _BMR_Bem_size() noexcept {
	int ret = -1;
	try {
		if (BMR().bemid < 0 || Bem().hydros.size() <= BMR().bemid) 
			throw Exc(F(t_("Wrong %d BEM case"), BMR().bemid));
		
		ret = Bem().hydros.size();
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return ret;
}

bool _BMR_Bem_Support(const char *solver, int *irregular, int *autoIrregular, int *middle7, int *far8, int *near9, int *autoCS, int *multibody) noexcept {
	for (int i = 0; i < Hydro::NUMBEM; ++i) {
		const Hydro::BEMInfo &info = Hydro::bemInfo[i];
		if (ToLower(info.str) == ToLower(solver)) {
			*multibody = info.multibody;
			String sqtf = F(info.qtf);
			*middle7 = sqtf.Find('7') >= 0;
			*far8 = sqtf.Find('8') >= 0;
			*near9 = sqtf.Find('9') >= 0;
			*irregular = info.irregular;
			*autoIrregular = info.autoIrregular;
			*autoCS = info.autoCS;
			return true;
		}
	}
	BMR().errorStr = F("Solver '%s' not supported", solver);
	return false;
}

bool _BMR_Bem_WaveOrigin_Set(double x, double y) noexcept {
	try {
		if (Bem().hydros.IsEmpty()) 
			throw Exc(t_("No bem case is loaded"));
		
		Bem().hydros[BMR().bemid].dt.x_w = x;
		Bem().hydros[BMR().bemid].dt.y_w = y;
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

bool _BMR_Bem_depth_Set(double h) noexcept {
	try {
		if (Bem().hydros.IsEmpty()) 
			throw Exc(t_("No bem case is loaded"));
		
		if (IsNull(h))
			throw Exc(F(t_("Wrong depth '%f'"), h));
		
		Bem().hydros[BMR().bemid].dt.h = h;
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

bool _BMR_Bem_g_Set(double g) noexcept {
	try {
		if (Bem().hydros.IsEmpty()) 
			throw Exc(t_("No bem case is loaded"));
		
		if (IsNull(g) || g < 0)
			throw Exc(F(t_("Wrong gravity '%f'"), g));
		
		Bem().hydros[BMR().bemid].dt.g = g;
	} catch(Exc err) {
		BMR().errorStr = err;
		return true;
	}
	BMR().errorStr.Clear();
	return false;
}

bool _BMR_Bem_rho_Set(double rho) noexcept {
	try {
		if (Bem().hydros.IsEmpty()) 
			throw Exc(t_("No bem case is loaded"));
		
		if (IsNull(rho) || rho < 0)
			throw Exc(F(t_("Wrong density '%f'"), rho));
		
		Bem().hydros[BMR().bemid].dt.rho = rho;
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

bool _BMR_Bem_w_Set(const double *w, int dim) noexcept {
	try {
		if (Bem().hydros.IsEmpty()) 
			throw Exc(t_("No bem case is loaded"));
		
		Hydro &hy = Bem().hydros[BMR().bemid];
		
		hy.dt.w.SetCount(dim);
		Copy(w, (size_t)dim, hy.dt.w);
		hy.dt.Nf = dim;
		
		hy.SortFrequencies();
		
		hy.dt.A.Clear();
		hy.dt.Ainf_w.Clear();
		hy.dt.A_P.Clear();
		hy.dt.B.Clear();
		hy.dt.B_H.Clear();
		hy.dt.B_P.Clear();
		
		hy.dt.ex.Clear();
		hy.dt.sc.Clear();
		hy.dt.fk.Clear();
		hy.dt.sc_pot.Clear();
		hy.dt.fk_pot.Clear();
		hy.dt.fk_pot_bmr.Clear();
		hy.dt.rao.Clear();
		
		hy.dt.qw = hy.Get_w();
		hy.dt.qtfsum.Clear();
		hy.dt.qtfdif.Clear();
		hy.dt.md.Clear();
		hy.dt.pots_rad.Clear();
		hy.dt.pots_dif.Clear();
		hy.dt.pots_inc.Clear();
		hy.dt.pots_inc_bmr.Clear();
		hy.dt.Apan.Clear();
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

bool _BMR_Bem_w_Get(double **data, int dim[1]) noexcept {
	try {
		if (Bem().hydros.IsEmpty()) 
			throw Exc(t_("No bem case is loaded"));
		
		Hydro &hy = Bem().hydros[BMR().bemid];
		
		dim[0] = hy.dt.Nf;
		*data = hy.dt.w.begin();
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

int _BMR_Bem_w_size() noexcept {
	double ret;
	try {
		if (Bem().hydros.IsEmpty()) 
			throw Exc(t_("No bem case is loaded"));
		
		Hydro &hy = Bem().hydros[BMR().bemid];
		
		ret = hy.dt.Nf;
	} catch(Exc err) {
		BMR().errorStr = err;
		return -1;
	}
	BMR().errorStr.Clear();
	return ret;
}

bool _BMR_Bem_headings_Set(const double *head, int dim) noexcept {
	try {
		if (Bem().hydros.IsEmpty()) 
			throw Exc(t_("No bem case is loaded"));
		
		Hydro &hy = Bem().hydros[BMR().bemid];
		
		hy.dt.head.SetCount(dim);
		Copy(head, (size_t)dim, hy.dt.head);
		hy.dt.Nh = hy.dt.head.size();
		
		hy.SortHeadings(BasicBEM::HEAD_0_360, BasicBEM::HEAD_0_360, BasicBEM::HEAD_0_360);
		
		hy.dt.ex.Clear();
		hy.dt.sc.Clear();
		hy.dt.fk.Clear();
		hy.dt.sc_pot.Clear();
		hy.dt.fk_pot.Clear();
		hy.dt.fk_pot_bmr.Clear();
		hy.dt.rao.Clear();									
		
		hy.dt.qtfsum.Clear();
		hy.dt.qtfdif.Clear();
		hy.dt.qhead = VectorXcd();

		hy.dt.mdhead = VectorXcd();
		hy.dt.md.Clear();
											
		hy.dt.pots_dif.Clear();
		hy.dt.pots_inc.Clear();
		hy.dt.pots_inc_bmr.Clear();						
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

bool _BMR_Bem_headings_Get(double **data, int dim[1]) noexcept {
	try {
		if (Bem().hydros.IsEmpty()) 
			throw Exc(t_("No bem case is loaded"));
		
		Hydro &hy = Bem().hydros[BMR().bemid];
		
		dim[0] = hy.dt.Nh;
		*data = hy.dt.head.begin();
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

int _BMR_Bem_headings_size() noexcept {
	double ret;
	try {
		if (Bem().hydros.IsEmpty()) 
			throw Exc(t_("No bem case is loaded"));
		
		Hydro &hy = Bem().hydros[BMR().bemid];
		
		ret = hy.dt.Nh;
	} catch(Exc err) {
		BMR().errorStr = err;
		return -1;
	}
	BMR().errorStr.Clear();
	return ret;
}

int _BMR_Bem_Duplicate() noexcept {
	try {
		if (Bem().hydros.IsEmpty()) 
			throw Exc(t_("No bem case is loaded"));
		
		Hydro &hy = Bem().hydros[BMR().bemid];
		_BMR_Bem_Id_Set(Bem().hydros.size());
		Bem().hydros[BMR().bemid] = clone(hy);									
	} catch(Exc err) {
		BMR().errorStr = err;
		return NullInt;
	}
	BMR().errorStr.Clear();
	return BMR().bemid;
}

static bool BMR_Bem_Body_Id_Allocate(int id) noexcept {
	try {
		if (Bem().hydros.IsEmpty()) 
			throw Exc(t_("No bem case is loaded"));
							
		if (IsNull(id) || id < 0)
			throw Exc(F(t_("Invalid id %d"), id));
		
		Hydro &hy = Bem().hydros[BMR().bemid];
		if (id >= hy.dt.msh.size()) {
			hy.dt.msh.SetCount(id+1);
			hy.dt.lids.SetCount(id+1);
			hy.dt.css.SetCount(id+1);
			hy.dt.Nb = hy.dt.msh.size();
		}
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

bool _BMR_Bem_Mesh_Load(int idBody, int idMesh) noexcept {
	try {
		if (BMR().bemid < 0 || Bem().hydros.size() <= BMR().bemid) 
			throw Exc(F(t_("Wrong %d BEM case"), BMR().bemid));
		
		BMR_Bem_Body_Id_Allocate(idBody);
		
		if (Bem().surfs.size() < idMesh) 
			throw Exc(F(t_("Id %d is not loaded"), idMesh));
		
		Hydro &hy = Bem().hydros[BMR().bemid];
		hy.dt.msh[idBody] = clone(Bem().surfs[idMesh]);
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

int _BMR_Bem_Mesh_size() noexcept {
	int ret = -1;
	try {
		if (BMR().bemid < 0 || Bem().hydros.size() <= BMR().bemid) 
			throw Exc(F(t_("Wrong %d BEM case"), BMR().bemid));
		
		Hydro &hy = Bem().hydros[BMR().bemid];
		ret = hy.dt.msh.size();
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return ret;
}

bool _BMR_Bem_Lid_Load(int idBody, int idMesh) noexcept {
	try {
		if (BMR().bemid < 0 || Bem().hydros.size() <= BMR().bemid) 
			throw Exc(F(t_("Wrong %d BEM case"), BMR().bemid));
		
		BMR_Bem_Body_Id_Allocate(idBody);
		
		if (Bem().surfs.size() < idMesh) 
			throw Exc(F(t_("Id %d is not loaded"), idMesh));
		
		Hydro &hy = Bem().hydros[BMR().bemid];
		hy.dt.lids[idBody] = clone(Bem().surfs[idMesh]);
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

bool _BMR_Bem_ControlSurface_Load(int idBody, int idMesh) noexcept {
	try {
		if (BMR().bemid < 0 || Bem().hydros.size() <= BMR().bemid) 
			throw Exc(F(t_("Wrong %d BEM case"), BMR().bemid));
		
		BMR_Bem_Body_Id_Allocate(idBody);
		
		if (Bem().surfs.size() < idMesh) 
			throw Exc(F(t_("Id %d is not loaded"), idMesh));
		
		Hydro &hy = Bem().hydros[BMR().bemid];
		hy.dt.css[idBody] = clone(Bem().surfs[idMesh]);
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();
	return true;
}

bool _BMR_Bem_SaveCase(const char *folder, const char *solver, bool x0z, bool y0z, 
		bool irregular, bool autoIrregular, const char *qtfType, bool autoQTF, 
		bool bin, int numCases, int numThreads, bool withPotentials, bool withMesh) noexcept {
	try {
		if (BMR().bemid < 0 || Bem().hydros.size() <= BMR().bemid) 
			throw Exc(F(t_("Wrong %d BEM case"), BMR().bemid));
		Hydro &hy = Bem().hydros[BMR().bemid];
		
		UVector<String> candidates;
		int type, icase;
		String lsolver = ToLower(Replace(solver, " ", ""));
		for (type = 0; type < Hydro::NUMBEM; ++type) {
			if (Hydro::bemInfo[type].caseCanSave) {
				String fmt = Hydro::GetBemStrCase(static_cast<Hydro::BEM_FMT>(type));
				if (ToLower(Replace(fmt, " ", "")) == lsolver) {				// Found the same
					candidates.Clear();
					candidates << fmt;
					icase = type;
					break;
				} else if (ToLower(Replace(fmt, " ", "")).Find(lsolver) >= 0) {	// Found similar
					candidates << fmt;
					icase = type;
				}
			}
		}

		if (candidates.IsEmpty())
			throw Exc(F(t_("Unsupported format %s for saving case"), solver));
		if (candidates.size() > 1) {
			String ret;
			for (int i = 0; i < candidates.size(); ++i) {
				if (i > 0)
					ret << ", ";
				ret << candidates[i];
			}
			throw Exc(F(t_("Format %s has more than one option: %s"), solver, ret));
		}

		int iqtfType = 0;
		String sqtfType = ToLower(qtfType);
		if (sqtfType.Find("control") >= 0 || sqtfType.Find("middle") >= 0)
			iqtfType = 7;
		else if (sqtfType.Find("momentum") >= 0 || sqtfType.Find("far") >= 0)
			iqtfType = 8;
		else if (sqtfType.Find("pressure") >= 0 || sqtfType.Find("near") >= 0)
			iqtfType = 9;
	
		if (sqtfType.Find("drift") >= 0 || sqtfType.Find("md") >= 0)
			iqtfType += 10;
		
		UVector<String> errors = hy.Check(static_cast<Hydro::BEM_FMT>(icase), irregular, autoIrregular, iqtfType, autoQTF);
		if (!errors.IsEmpty()) {
			String str;
			if (errors.size() == 1)
				str << "\n " << errors[0];
			else {
				for (int i = 0; i < errors.size(); ++i)
				 	str << "\n- " << errors[i];
			}
			BEM::Print(F(t_("\nProblems found in data: %s"), str));
		}		
		UVector<bool> listDOF(6, true);
		UVector<Point3D> dummy;
		hy.SaveCase(folder, static_cast<Hydro::BEM_FMT>(icase), x0z, y0z, 
				irregular, autoIrregular, iqtfType, autoQTF, 
				bin, numCases, numThreads, withPotentials, withMesh, listDOF, dummy);
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();	
	return true;
}
	
bool _BMR_Bem_FroudeKrylov_Calc() noexcept {
	try {
		if (BMR().bemid < 0 || Bem().hydros.size() <= BMR().bemid) 
			throw Exc(F(t_("Wrong %d BEM case"), BMR().bemid));
		Hydro &hy = Bem().hydros[BMR().bemid];
		
		hy.GetPotentialsIncident();
		hy.GetForcesFromPotentials(hy.dt.pots_inc_bmr, hy.dt.fk);	
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();	
	return true;
}

bool _BMR_Bem_Diffraction_Reset() noexcept {
	try {
		if (BMR().bemid < 0 || Bem().hydros.size() <= BMR().bemid) 
			throw Exc(F(t_("Wrong %d BEM case"), BMR().bemid));
		Hydro &hy = Bem().hydros[BMR().bemid];
		
		hy.Initialize_Forces(hy.dt.sc, -1, std::complex<double>(0, 0));
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();	
	return true;	
}

bool _BMR_Bem_Excitation_Reset() noexcept {
	try {
		if (BMR().bemid < 0 || Bem().hydros.size() <= BMR().bemid) 
			throw Exc(F(t_("Wrong %d BEM case"), BMR().bemid));
		Hydro &hy = Bem().hydros[BMR().bemid];
		
		hy.Initialize_Forces(hy.dt.ex, -1, std::complex<double>(0, 0));
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();	
	return true;	
}

bool _BMR_Bem_Excitation_GetFromDiffFK() noexcept {	
	try {
		if (BMR().bemid < 0 || Bem().hydros.size() <= BMR().bemid) 
			throw Exc(F(t_("Wrong %d BEM case"), BMR().bemid));
		Hydro &hy = Bem().hydros[BMR().bemid];
		
		hy.GetFexFromFscFfk();
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();	
	return true;	
}

bool _BMR_Bem_AddedMass_Reset() noexcept {
	try {
		if (BMR().bemid < 0 || Bem().hydros.size() <= BMR().bemid) 
			throw Exc(F(t_("Wrong %d BEM case"), BMR().bemid));
		Hydro &hy = Bem().hydros[BMR().bemid];
		
		hy.Initialize_AB(hy.dt.A, 0);
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();	
	return true;	
}

bool _BMR_Bem_RadiationDamping_Reset() noexcept {
	try {
		if (BMR().bemid < 0 || Bem().hydros.size() <= BMR().bemid) 
			throw Exc(F(t_("Wrong %d BEM case"), BMR().bemid));
		Hydro &hy = Bem().hydros[BMR().bemid];
		
		hy.Initialize_AB(hy.dt.B, 0);
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();	
	return true;	
}

bool _BMR_Bem_AB_Set(bool a, int idBodyRow, int idBodyCol, int idBemFrom, int idBodyFromRow, int idBodyFromCol) noexcept {
	try {
		if (BMR().bemid < 0 || Bem().hydros.size() <= BMR().bemid) 
			throw Exc(F(t_("Wrong %d BEM case"), BMR().bemid));
		Hydro &hy = Bem().hydros[BMR().bemid];

		if (idBodyRow < 0 || hy.dt.Nb <= idBodyRow)
			throw Exc(F(t_("Wrong id for body %d of BEM case %d"), idBodyRow, BMR().bemid));
		if (idBodyCol < 0 || hy.dt.Nb <= idBodyCol)
			throw Exc(F(t_("Wrong id for body %d of BEM case %d"), idBodyCol, BMR().bemid));
						
		if (idBemFrom < 0 || Bem().hydros.size() < idBemFrom) 
			throw Exc(F(t_("Wrong %d BEM case"), idBemFrom));
		
		Hydro &hyFrom = Bem().hydros[idBemFrom];
		
		if (idBodyFromRow < 0 || hyFrom.dt.Nb <= idBodyFromRow)
			throw Exc(F(t_("Wrong id for body %d of BEM case %d"), idBodyFromRow, idBemFrom));
		if (idBodyFromCol < 0 || hyFrom.dt.Nb <= idBodyFromCol)
			throw Exc(F(t_("Wrong id for body %d of BEM case %d"), idBodyFromCol, idBemFrom));
		
		if (hy.dt.Nf != hyFrom.dt.Nf)
			throw Exc(F(t_("Nomber of frequencies of BEM cases %d and %d do not match (%d != %d)"), BMR().bemid, idBemFrom, hy.dt.Nf, hyFrom.dt.Nf));

		if (a)
			hy.Set_AB(hy.dt.A, idBodyRow, idBodyCol, hyFrom.dt.A, idBodyFromRow, idBodyFromCol);
		else
			hy.Set_AB(hy.dt.B, idBodyRow, idBodyCol, hyFrom.dt.B, idBodyFromRow, idBodyFromCol);
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();	
	return true;
}
	
bool _BMR_Bem_AddedMass_Set(int idBodyRow, int idBodyCol, int idBemFrom, int idBodyFromRow, int idBodyFromCol) noexcept {
	return _BMR_Bem_AB_Set(true, idBodyRow, idBodyCol, idBemFrom, idBodyFromRow, idBodyFromCol);
}

bool _BMR_Bem_RadiationDamping_Set(int idBodyRow, int idBodyCol, int idBemFrom, int idBodyFromRow, int idBodyFromCol) noexcept {
	return _BMR_Bem_AB_Set(false, idBodyRow, idBodyCol, idBemFrom, idBodyFromRow, idBodyFromCol);
}

int _BMR_Bem_MapToMesh(int idBody, int idMesh, double tolerance, bool rad, bool diff, bool inc, bool relatedToBody) noexcept {
	int newId;
	try {
		if (BMR().bemid < 0 || Bem().hydros.size() <= BMR().bemid) 
			throw Exc(F(t_("Wrong %d BEM case"), BMR().bemid));
		Hydro &hy = Bem().hydros[BMR().bemid];

		if (idBody < 0 || hy.dt.Nb <= idBody)
			throw Exc(F(t_("Wrong id for body %d of BEM case %d"), idBody, BMR().bemid));
		
		if (idMesh < 0 || idMesh >= Bem().surfs.size())
			throw Exc(F(t_("Wrong id for mesh %d"), idMesh));
		
		UVector<int> idms = {idMesh};

		hy.MapMeshes(Bem().hydros, idBody, idms, true, relatedToBody, tolerance, rad, diff, inc);
		
		newId = Bem().hydros.size()-1;
		
		Bem().hydros[newId].FillWithPotentials();
		
	} catch(Exc err) {
		BMR().errorStr = err;
		return -1;
	}
	BMR().errorStr.Clear();	
	return newId;
}		

bool _BMR_Bem_FroudeKrylov_Set(int idBody, int idBemFrom, int idBodyFrom) noexcept {
	try {
		if (BMR().bemid < 0 || Bem().hydros.size() <= BMR().bemid) 
			throw Exc(F(t_("Wrong %d BEM case"), BMR().bemid));
		Hydro &hy = Bem().hydros[BMR().bemid];

		if (idBody < 0 || hy.dt.Nb <= idBody)
			throw Exc(F(t_("Wrong id for body %d of BEM case %d"), idBody, BMR().bemid));
					
		if (idBemFrom < 0 || Bem().hydros.size() < idBemFrom) 
			throw Exc(F(t_("Wrong %d BEM case"), idBemFrom));
		Hydro &hyFrom = Bem().hydros[idBemFrom];
		
		if (idBodyFrom < 0 || hyFrom.dt.Nb <= idBodyFrom)
			throw Exc(F(t_("Wrong id for body %d of BEM case %d"), idBodyFrom, idBemFrom));
		
		if (hy.dt.Nf != hyFrom.dt.Nf)
			throw Exc(F(t_("Nomber of frequencies of BEM cases %d and %d do not match (%d != %d)"), BMR().bemid, idBemFrom, hy.dt.Nf, hyFrom.dt.Nf));

		hy.Set_Force(hy.dt.fk, idBody, hyFrom.dt.fk, idBodyFrom,  hyFrom.dt.head);
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();	
	return true;
}

bool _BMR_Bem_Diffraction_Set(int idBody, int idBemFrom, int idBodyFrom) noexcept {
	try {
		if (BMR().bemid < 0 || Bem().hydros.size() <= BMR().bemid) 
			throw Exc(F(t_("Wrong %d BEM case"), BMR().bemid));
		Hydro &hy = Bem().hydros[BMR().bemid];

		if (idBody < 0 || hy.dt.Nb <= idBody)
			throw Exc(F(t_("Wrong id for body %d of BEM case %d"), idBody, BMR().bemid));
					
		if (idBemFrom < 0 || Bem().hydros.size() < idBemFrom) 
			throw Exc(F(t_("Wrong %d BEM case"), idBemFrom));
		Hydro &hyFrom = Bem().hydros[idBemFrom];
		
		if (idBodyFrom < 0 || hyFrom.dt.Nb <= idBodyFrom)
			throw Exc(F(t_("Wrong id for body %d of BEM case %d"), idBodyFrom, idBemFrom));
		
		if (hy.dt.Nf != hyFrom.dt.Nf)
			throw Exc(F(t_("Nomber of frequencies of BEM cases %d and %d do not match (%d != %d)"), BMR().bemid, idBemFrom, hy.dt.Nf, hyFrom.dt.Nf));

		hy.Set_Force(hy.dt.sc, idBody, hyFrom.dt.sc, idBodyFrom, hyFrom.dt.head);
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();	
	return true;
}

bool _BMR_Bem_Excitation_Set(int idBody, int idBemFrom, int idBodyFrom) noexcept {
	try {
		if (BMR().bemid < 0 || Bem().hydros.size() <= BMR().bemid) 
			throw Exc(F(t_("Wrong %d BEM case"), BMR().bemid));
		Hydro &hy = Bem().hydros[BMR().bemid];

		if (idBody < 0 || hy.dt.Nb <= idBody)
			throw Exc(F(t_("Wrong id for body %d of BEM case %d"), idBody, BMR().bemid));
					
		if (idBemFrom < 0 || Bem().hydros.size() < idBemFrom) 
			throw Exc(F(t_("Wrong %d BEM case"), idBemFrom));
		Hydro &hyFrom = Bem().hydros[idBemFrom];
		
		if (idBodyFrom < 0 || hyFrom.dt.Nb <= idBodyFrom)
			throw Exc(F(t_("Wrong id for body %d of BEM case %d"), idBodyFrom, idBemFrom));
		
		if (hy.dt.Nf != hyFrom.dt.Nf)
			throw Exc(F(t_("Nomber of frequencies of BEM cases %d and %d do not match (%d != %d)"), BMR().bemid, idBemFrom, hy.dt.Nf, hyFrom.dt.Nf));

		hy.Set_Force(hy.dt.ex, idBody, hyFrom.dt.ex, idBodyFrom,  hyFrom.dt.head);
	} catch(Exc err) {
		BMR().errorStr = err;
		return false;
	}
	BMR().errorStr.Clear();	
	return true;
}
			
bool _BMR_FAST_Load(const char *filename) noexcept {
	try {
		String ret = BMR().fast.Load(filename, Null);
		if (ret.IsEmpty()) {
			BMR().errorStr.Clear();	
			return true;
		} else {
			BMR().errorStr = ret;
			return false;
		}
	} catch (...) {
		BMR().errorStr = "Unknown error in BMR_FAST_Load()";
		return false;
	}
}

const char *_BMR_FAST_GetParameterName(int id) noexcept {
	static String ret;
	try {
		BMR().errorStr.Clear();	
		return ret = BMR().fast.GetParameter(id);
	} catch (...) {
		BMR().errorStr = "Unknown error in BMR_FAST_GetParameterName()";
		return ret = "";
	}
}

const char *_BMR_FAST_GetUnitName(int id) noexcept {
	static String ret;
	try {
		BMR().errorStr.Clear();	
		return ret = BMR().fast.GetUnit(id);
	} catch (...) {
		BMR().errorStr = "Unknown error in BMR_FAST_GetUnitName()";
		return ret = "";
	}
}

int _BMR_FAST_GetParameterId(const char *name) noexcept {
	try {
		UVector<int> p = BMR().fast.FindParameterMatch(name);
		if (p.IsEmpty())
			return NullInt;
		else {
			BMR().errorStr.Clear();	
			return p[0];
		}
	} catch (...) {
		Cout() << "Unknown error in BMR_FAST_GetParameterCount()";
		return NullInt;
	}
}

int _BMR_FAST_GetParameterCount() noexcept {
	try {
		BMR().errorStr.Clear();	
		return BMR().fast.GetParameterCount();
	} catch (...) {
		BMR().errorStr = "Unknown error in BMR_FAST_GetParameterCount()";
		return NullInt;
	}
}

int _BMR_FAST_GetLen() noexcept {
	try {
		BMR().errorStr.Clear();	
		return BMR().fast.GetNumData();
	} catch (...) {
		BMR().errorStr = "Unknown error in BMR_FAST_GetLen()";
		return NullInt;
	}
}

double _BMR_FAST_GetTimeStart() noexcept {
	try {
		BMR().errorStr.Clear();	
		return BMR().fast.GetTimeStart();
	} catch (...) {
		BMR().errorStr = "Unknown error in BMR_FAST_GetTimeStart()";
		return NullInt;
	}
}

double _BMR_FAST_GetTimeEnd() noexcept {
	try {
		BMR().errorStr.Clear();	
		return BMR().fast.GetTimeEnd();
	} catch (...) {
		BMR().errorStr = "Unknown error in BMR_FAST_GetTimeEnd()";
		return NullDouble;
	}
}

double _BMR_FAST_GetTime(int idtime) noexcept {
	return _BMR_FAST_GetData(idtime, 0);
}

double _BMR_FAST_GetData(int idtime, int idparam) noexcept {
	if (idtime < 0) {
		BMR().errorStr = "Error in BMR_FAST_GetData() idtime < 0";
		return NullDouble;
	}
	if (idtime >= BMR().fast.GetNumData()) {
		BMR().errorStr = "Error in BMR_FAST_GetData() idtime >= time";
		return NullDouble;
	}
	if (idparam < 0) {
		BMR().errorStr = "Error in BMR_FAST_GetData() idparam < 0";
		return NullDouble;
	}
	if (idparam >= BMR().fast.GetParameterCount()) {
		BMR().errorStr = "Error in BMR_FAST_GetData() idparam >= num_params";
		return NullDouble;
	}
	BMR().errorStr.Clear();		
	return BMR().fast.GetVal(idtime, idparam);
}

static void BMR_FAST_GetData(int idparam, int idbegin, int idend, VectorXd &data) {
	if (idparam < 0) 
		throw Exc("idparam < 0");
	if (idparam >= BMR().fast.GetParameterCount()) 
		throw Exc("idparam >= num_params");

	if (idbegin < 0)
		idbegin = 0;
	if (idend < 0)
		idend = BMR().fast.GetNumData()-1;
	
	if (idbegin > idend) 
		throw Exc("idbegin > idend");
		
	data = BMR().fast.GetVector(idparam).segment(idbegin, idend - idbegin + 1);
}

bool _BMR_FAST_GetArray(int idparam, int idbegin, int idend, double **data, int *num) noexcept {
	static VectorXd v;
	
	try {
		BMR_FAST_GetData(idparam, idbegin, idend, v);
		
		*num = int(v.size());
		*data = v.data();
		
		BMR().errorStr.Clear();	
		return true;
	} catch (Exc e) {
		BMR().errorStr = F("Error in BMR_FAST_GetArray(): %s", e);
	} catch (...) {
		BMR().errorStr = "Unknown error in BMR_FAST_GetArray()";
	}
	return false;	
}

double _BMR_FAST_GetAvg(int idparam, int idbegin, int idend) noexcept {
	try {
		VectorXd data;
		
		BMR_FAST_GetData(idparam, idbegin, idend, data);
		
		BMR().errorStr.Clear();	
		return data.mean();
	} catch (Exc e) {
		BMR().errorStr = F("Error in BMR_FAST_GetAvg(): %s", e);
	} catch (...) {
		BMR().errorStr = "Unknown error in BMR_FAST_GetAvg()";
	}
	return NullDouble;
}

double _BMR_FAST_GetMax(int idparam, int idbegin, int idend) noexcept {
	try {
		VectorXd data;
		
		BMR_FAST_GetData(idparam, idbegin, idend, data);
		
		BMR().errorStr.Clear();	
		return data.maxCoeff();
	} catch (Exc e) {
		BMR().errorStr = F("Error in BMR_FAST_GetMax(): %s", e);
	} catch (...) {
		BMR().errorStr = "Unknown error in BMR_FAST_GetAvg()";
	}
	return NullDouble;
}

double _BMR_FAST_GetMin(int idparam, int idbegin, int idend) noexcept {
	try {
		VectorXd data;
		
		BMR_FAST_GetData(idparam, idbegin, idend, data);
		
		BMR().errorStr.Clear();	
		return data.minCoeff();
	} catch (Exc e) {
		BMR().errorStr = F("Error in BMR_FAST_GetMin(): %s", e);
	} catch (...) {
		BMR().errorStr = "Unknown error in BMR_FAST_GetAvg()";
	}
	return NullDouble;
}

bool _BMR_FAST_LoadFile(const char *file) noexcept {
	BMR().fastFileStr = LoadFile(file);
	BMR().fastFileName = file;
	return !BMR().fastFileStr.IsEmpty();
}

bool _BMR_FAST_SaveFile(const char *file) noexcept {
	bool ret;
	try {
		String sfile(file);
		if (!sfile.IsEmpty()) 
			BMR().fastFileName = sfile;
		ret = SaveFile(BMR().fastFileName, BMR().fastFileStr);
		
	} catch (Exc e) {
		BMR().errorStr = e;
		return false;
	}
	BMR().errorStr.Clear();	
	return ret;	
}

bool _BMR_FAST_SetVar(const char *name, const char *paragraph, const char *value) noexcept {
	try {
		SetFASTVar(BMR().fastFileStr, name, value, paragraph);
	} catch (Exc e) {
		BMR().errorStr = e;
		return false;
	}
	BMR().errorStr.Clear();	
	return true;
}

const char *_BMR_FAST_GetVar(const char *name, const char *paragraph) noexcept {
	static String ret;

	try {
		ret = GetFASTVar(BMR().fastFileStr, name, paragraph);
	} catch (Exc e) {
		BMR().errorStr = e;
		return ret = "";
	}
	if (IsVoid(ret))
		return ret = "";

	BMR().errorStr.Clear();	
	return ret;
}

double _BMR_DemoVectorPyC(const double *v, int num) noexcept {
    double res = 0;
    for (int i = 0; i < num; ++i) 
        res += v[i];
    return res;
}


#endif

#if defined(flagBEMR_TEST_DLL_INTERNAL) || defined(flagBEMR_TEST_DLL)
void BEM_Throw() {
	if (_BMR_GetLastError())
		throw Exc(_BMR_GetLastError());
}
#endif
//#endif

//#endif