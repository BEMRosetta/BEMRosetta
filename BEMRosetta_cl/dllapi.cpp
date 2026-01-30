// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2026, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include <SysInfo/Crash.h>

BEM &Bem() 		{static BEM bem;		return bem;}

#ifdef PLATFORM_WIN32
#include "orca.h"

Function<bool(String, int, const Time &)> Orca::WhenWave = [](String str, int perc, const Time &et)->bool {
	if (IsNull(et))
		BEM::Print("\nCompleted 0%"); 
	else
		BEM::Print(Format("\nCompleted %d%%. Et: %", perc, et)); 
	return 0;
};

Function<bool(String)> Orca::WhenPrint = [](String str)->bool {
	BEM::Print(Format("\n%s", str)); 
	return 0;
};

Time Orca::startCalc = Null, Orca::beginNoLicense = Null, Orca::lastLog;
int64 Orca::noLicenseTime = 0;

#endif

//#if defined(flagBEMR_DLL) || defined(flagBEMR_TEST_DLL_INTERNAL) || defined(flagBEMR_TEST_DLL) || defined(flagBEMR_CL)

#include "FastOut.h"
#include "export.h"


DLL_Data &DLL() {
	static DLL_Data dll;
	return dll;
}

void DLL_Init() noexcept {
	DLL();
}

const char *DLL_GetLastError() noexcept {
	if (DLL().errorStr.IsEmpty())
		return nullptr;
	return DLL().errorStr;
}
														
																												
void DLL_NoPrint() noexcept {
	CoutStreamX::NoPrint();
}
	
const char *DLL_Version() noexcept {
	static String version;
	version << __DATE__ << ", " << __TIME__;
	return version;	
}

void DLL_Echo(const char *str) noexcept {
	try {
		CoutX() << str;
	} catch(Exc err) {
		DLL().errorStr = err;
		return;
	}
	DLL().errorStr.Clear();
}


void DLL_Mesh_Input(const char *file) noexcept {
	try {
		if (!FileExists(file))
			throw Exc(Format(t_("File '%s' not found"), file)); 
								
		BEM::Print("\n");
		Bem().LoadBody(file, DLL().echo ? DLL().Status : DLL().NoPrint, false, false);		// Doesn't work for multibody .dat
		DLL().meshid = Bem().surfs.size() - 1;
	} catch(Exc err) {
		DLL().errorStr = err;
		return;
	}
	DLL().errorStr.Clear();
}

void DLL_Mesh_Report() noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Bem().surfs[DLL().meshid].Report(Bem().rho);
	} catch(Exc err) {
		DLL().errorStr = err;
		return;
	}
	DLL().errorStr.Clear();
}

void DLL_Mesh_Clear() noexcept {
	try {
		Bem().surfs.Clear();
		DLL().meshid = -1;
		BEM::Print("\n" + S(t_("Body data cleared")));
	} catch(Exc err) {
		DLL().errorStr = err;
		return;
	}
	DLL().errorStr.Clear();
}

void DLL_Mesh_SetId(int id) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
							
		if (IsNull(id) || id < 0 || id > Bem().surfs.size()-1)
			throw Exc(Format(t_("Invalid id %d"), id));
		
		DLL().meshid = id;
		BEM::Print("\n" + Format(t_("Body active model id is %d"), DLL().meshid));	
	} catch(Exc err) {
		DLL().errorStr = err;
		return;
	}
	DLL().errorStr.Clear();
}

void DLL_Mesh_Convert(const char *file, const char *format, int symX, int symY) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		
		Body::MESH_FMT meshFmt = Body::GetCodeBodyStr(format);
		if (!Body::meshInfo[meshFmt].canSave)
			throw Exc(Format(t_("Saving format '%s' is not implemented"), format));						
		
		Body::SaveAs(Bem().surfs[DLL().meshid], file, meshFmt, Body::ALL, Bem().rho, Bem().g, symX, symY);
		BEM::Print("\n" + Format(t_("Model id %d saved as '%s'"), DLL().meshid, file));						
	} catch(Exc err) {
		DLL().errorStr = err;
		return;
	}
	DLL().errorStr.Clear();
}

void DLL_Mesh_Translate(double x, double y, double z) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Body &msh = Bem().surfs[DLL().meshid];
		msh.dt.mesh.Translate(x, y, z);
		msh.dt.cg.Translate(x, y, z);
		msh.AfterLoad(Bem().rho, Bem().g, false, false);	
		BEM::Print("\n" + Format(t_("Body id %d translated %f, %f, %f"), DLL().meshid, x, y, z)); 
	} catch(Exc err) {
		DLL().errorStr = err;
		return;
	}
	DLL().errorStr.Clear();
}

void DLL_Mesh_Rotate(double ax, double ay, double az, double cx, double cy, double cz) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Body &msh = Bem().surfs[DLL().meshid];
		msh.dt.mesh.Rotate(ToRad(ax), ToRad(ay), ToRad(az), cx, cy, cz);	
		msh.dt.cg.Rotate(ToRad(ax), ToRad(ay), ToRad(az), cx, cy, cz);
		msh.AfterLoad(Bem().rho, Bem().g, false, false);	
		BEM::Print("\n" + Format(t_("Body id %d rotated angles %f, %f, %f around centre %f, %f, %f"), DLL().meshid, ax, ay, az, cx, cy, cz));
	} catch(Exc err) {
		DLL().errorStr = err;
		return;
	}
	DLL().errorStr.Clear();
}

void DLL_Mesh_Cg(double x, double y, double z) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Body &msh = Bem().surfs[DLL().meshid];
		msh.dt.cg = Point3D(x, y, z);
		msh.AfterLoad(Bem().rho, Bem().g, true, false);
		BEM::Print("\n" + Format(t_("CG is %f, %f, %f"), x, y, z));	
	} catch(Exc err) {
		DLL().errorStr = err;
		return;
	}
	DLL().errorStr.Clear();
}

void DLL_Mesh_C0(double x, double y, double z) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Body &msh = Bem().surfs[DLL().meshid];
		msh.dt.c0 = Point3D(x, y, z);
		msh.AfterLoad(Bem().rho, Bem().g, true, false);
		BEM::Print("\n" + Format(t_("C0 is %f, %f, %f"), x, y, z));	
	} catch(Exc err) {
		DLL().errorStr = err;
		return;
	}
	DLL().errorStr.Clear();
}

void DLL_Mesh_Mass(double mass) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Body &msh = Bem().surfs[DLL().meshid];
		msh.SetMass(mass);
		msh.AfterLoad(Bem().rho, Bem().g, true, false);
		BEM::Print("\n" + Format(t_("Mass is %f"), mass));
	} catch(Exc err) {
		DLL().errorStr = err;
		return;
	}
	DLL().errorStr.Clear();
}

void DLL_Mesh_Reset() noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Bem().surfs[DLL().meshid].Reset(Bem().rho, Bem().g);
		BEM::Print("\n" + Format(t_("Body id %d position is reset"), DLL().meshid));
	} catch(Exc err) {
		DLL().errorStr = err;
		return;
	}
	DLL().errorStr.Clear();
}

void DLL_Mesh_GetWaterPlane() noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Bem().AddWaterSurface(DLL().meshid, 'e', 1, false);
		DLL().meshid = Bem().surfs.size() - 1;
		BEM::Print("\n" + Format(t_("Body id %d waterplane is got"), DLL().meshid));
	} catch(Exc err) {
		DLL().errorStr = err;
		return;
	}
	DLL().errorStr.Clear();
}

void DLL_Mesh_GetHull() noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Bem().AddWaterSurface(DLL().meshid, 'r', 1, false);
		DLL().meshid = Bem().surfs.size() - 1;
		BEM::Print("\n" + Format(t_("Body id %d hull is got"), DLL().meshid));
	} catch(Exc err) {
		DLL().errorStr = err;
		return;
	}
	DLL().errorStr.Clear();
}

void DLL_Mesh_FillWaterplane(double ratio) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		if (ratio < 0 || ratio > 100)
			throw Exc(Format(t_("Wrong mesh ratio %s"), ratio));
		Bem().AddWaterSurface(DLL().meshid, 'f', ratio, false);
		DLL().meshid = Bem().surfs.size() - 1;
		BEM::Print("\n" + Format(t_("Body id %d lid is got"), DLL().meshid));
	} catch(Exc err) {
		DLL().errorStr = err;
		return;
	}
	DLL().errorStr.Clear();
}
						
void DLL_Mesh_GetVolume(double *vx, double *vy, double *vz) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Body &msh = Bem().surfs[DLL().meshid];
		*vx = msh.dt.mesh.volumex;
		*vy = msh.dt.mesh.volumey;
		*vz = msh.dt.mesh.volumez;	
	} catch(Exc err) {
		DLL().errorStr = err;
		return;
	}
	DLL().errorStr.Clear();
}

void DLL_Mesh_GetUnderwaterVolume(double *vx, double *vy, double *vz) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Body &msh = Bem().surfs[DLL().meshid];
		*vx = msh.dt.under.volumex;
		*vy = msh.dt.under.volumey;
		*vz = msh.dt.under.volumez;
	} catch(Exc err) {
		DLL().errorStr = err;
		return;
	}
	DLL().errorStr.Clear();
}

double DLL_Mesh_GetSurface() noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
	} catch(Exc err) {
		DLL().errorStr = err;
		return Null;
	}
	DLL().errorStr.Clear();
	return Bem().surfs[DLL().meshid].dt.mesh.surface;
}

double DLL_Mesh_GetUnderwaterSurface() noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
	} catch(Exc err) {
		DLL().errorStr = err;
		return Null;
	}
	DLL().errorStr.Clear();
	return Bem().surfs[DLL().meshid].dt.under.surface;
}

void DLL_Mesh_GetCentreOfGravity(double *x, double *y, double *z) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Body &msh = Bem().surfs[DLL().meshid];
		*x = msh.dt.cb.x;
		*y = msh.dt.cb.y;
		*z = msh.dt.cb.z;
	} catch(Exc err) {
		DLL().errorStr = err;
		return;
	}
	DLL().errorStr.Clear();
}

void DLL_Mesh_GetCentreOfGravity_Volume(double *x, double *y, double *z) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Body &msh = Bem().surfs[DLL().meshid];
		Point3D cg = msh.dt.mesh.GetCentreOfBuoyancy();
		*x = cg.x;
		*y = cg.y;
		*z = cg.z;
	} catch(Exc err) {
		DLL().errorStr = err;
		return;
	}
	DLL().errorStr.Clear();
}

void DLL_Mesh_GetCentreOfGravity_Surface(double *x, double *y, double *z) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Body &msh = Bem().surfs[DLL().meshid];
		Point3D cg = msh.dt.mesh.GetCentreOfGravity_Surface();
		*x = cg.x;
		*y = cg.y;
		*z = cg.z;
	} catch(Exc err) {
		DLL().errorStr = err;
		return;
	}
	DLL().errorStr.Clear();
}

void DLL_Mesh_GetHydrostaticStiffness(double **data, int dim[2]) noexcept {
	static UVector<double> d;
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Body &msh = Bem().surfs[DLL().meshid];
		dim[0] = msh.dt.C.rows();
		dim[1] = msh.dt.C.cols();
		CopyRowMajor(msh.dt.C, d);
		*data = d.begin();
	} catch(Exc err) {
		DLL().errorStr = err;
		return;
	}
	DLL().errorStr.Clear();
}
								
									



int DLL_FAST_Load(const char *filename) noexcept {
	try {
		String ret = DLL().fast.Load(filename, Null);
		if (ret.IsEmpty())
			return 1;
		else {
			CoutX() << Format("Error: %s", ret);
			return 0;
		}
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_Load()";
		return 0;
	}
}

const char *DLL_FAST_GetParameterName(int id) noexcept {
	static String ret;
	try {
		return ret = DLL().fast.GetParameter(id);
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetParameterName()";
		return ret = "Error";
	}
}

const char *DLL_FAST_GetUnitName(int id) noexcept {
	static String ret;
	try {
		return ret = DLL().fast.GetUnit(id);
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetUnitName()";
		return ret = "Error";
	}
}

int DLL_FAST_GetParameterId(const char *name) noexcept {
	try {
		UVector<int> p = DLL().fast.FindParameterMatch(name);
		if (p.IsEmpty())
			return -1;
		else
			return p[0];
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetParameterCount()";
		return Null;
	}
}

int DLL_FAST_GetParameterCount() noexcept {
	try {
		return DLL().fast.GetParameterCount();
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetParameterCount()";
		return Null;
	}
}

int DLL_FAST_GetLen() noexcept {
	try {
		return DLL().fast.GetNumData();
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetLen()";
		return Null;
	}
}

double DLL_FAST_GetTimeStart() noexcept {
	try {
		return DLL().fast.GetTimeStart();
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetTimeStart()";
		return Null;
	}
}

double DLL_FAST_GetTimeEnd() noexcept {
	try {
		return DLL().fast.GetTimeEnd();
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetTimeEnd()";
		return Null;
	}
}

double DLL_FAST_GetTime(int idtime) noexcept {
	return DLL_FAST_GetData(idtime, 0);
}

double DLL_FAST_GetData(int idtime, int idparam) noexcept {
	if (idtime < 0) {
		CoutX() << "Error in DLL_FAST_GetData() idtime < 0";
		return Null;
	}
	if (idtime >= DLL().fast.GetNumData()) {
		CoutX() << "Error in DLL_FAST_GetData() idtime >= time";
		return Null;
	}
	if (idparam < 0) {
		CoutX() << "Error in DLL_FAST_GetData() idparam < 0";
		return Null;
	}
	if (idparam >= DLL().fast.GetParameterCount()) {
		CoutX() << "Error in DLL_FAST_GetData() idparam >= num_params";
		return Null;
	}
		
	return DLL().fast.GetVal(idtime, idparam);
}

static void DLL_FAST_GetData(int idparam, int idbegin, int idend, VectorXd &data) {
	if (idparam < 0) 
		throw Exc("idparam < 0");
	if (idparam >= DLL().fast.GetParameterCount()) 
		throw Exc("idparam >= num_params");

	if (idbegin < 0)
		idbegin = 0;
	if (idend < 0)
		idend = DLL().fast.GetNumData()-1;
	
	if (idbegin > idend) 
		throw Exc("idbegin > idend");
		
	data = DLL().fast.GetVector(idparam).segment(idbegin, idend - idbegin + 1);
}

int DLL_FAST_GetArray(int idparam, int idbegin, int idend, double **data, int *num) noexcept {
	static VectorXd v;
	
	try {
		DLL_FAST_GetData(idparam, idbegin, idend, v);
		
		*num = int(v.size());
		*data = v.data();
		return 1;
	} catch (Exc e) {
		CoutX() << Format("Error in DLL_FAST_GetArray(): %s", e);
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetArray()";
	}
	return Null;	
}

double DLL_FAST_GetAvg(int idparam, int idbegin, int idend) noexcept {
	try {
		VectorXd data;
		
		DLL_FAST_GetData(idparam, idbegin, idend, data);
		
		return data.mean();
	} catch (Exc e) {
		CoutX() << Format("Error in DLL_FAST_GetAvg(): %s", e);
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetAvg()";
	}
	return Null;
}

double DLL_FAST_GetMax(int idparam, int idbegin, int idend) noexcept {
	try {
		VectorXd data;
		
		DLL_FAST_GetData(idparam, idbegin, idend, data);
		
		return data.maxCoeff();
	} catch (Exc e) {
		CoutX() << Format("Error in DLL_FAST_GetMax(): %s", e);
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetAvg()";
	}
	return Null;
}

double DLL_FAST_GetMin(int idparam, int idbegin, int idend) noexcept {
	try {
		VectorXd data;
		
		DLL_FAST_GetData(idparam, idbegin, idend, data);
		
		return data.minCoeff();
	} catch (Exc e) {
		CoutX() << Format("Error in DLL_FAST_GetMin(): %s", e);
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetAvg()";
	}
	return Null;
}

int DLL_IsNull(double val) noexcept {return IsNull(val);}


int DLL_FAST_LoadFile(const char *file) noexcept {
	DLL().fastFileStr = LoadFile(file);
	DLL().fastFileName = file;
	return !DLL().fastFileStr.IsEmpty();
}

int DLL_FAST_SaveFile(const char *file) noexcept {
	bool ret;
	try {
		String sfile(file);
		if (!sfile.IsEmpty()) 
			DLL().fastFileName = sfile;
		ret = SaveFile(DLL().fastFileName, DLL().fastFileStr);
		
	} catch (Exc e) {
		SetConsoleColor(CONSOLE_COLOR::LTYELLOW);
		CoutX() << "\n" << "Error: " << e;
		SetConsoleColor(CONSOLE_COLOR::PREVIOUS);
		return false;
	}
	return ret;	
}

int DLL_FAST_SetVar(const char *name, const char *paragraph, const char *value) noexcept {
	try {
		SetFASTVar(DLL().fastFileStr, name, value, paragraph);
	} catch (Exc e) {
		SetConsoleColor(CONSOLE_COLOR::LTYELLOW);
		CoutX() << "\n" << "Error: " << e;
		SetConsoleColor(CONSOLE_COLOR::PREVIOUS);
		return false;
	}
	return true;
}

const char *DLL_FAST_GetVar(const char *name, const char *paragraph) noexcept {
	static String ret;

	try {
		ret = GetFASTVar(DLL().fastFileStr, name, paragraph);
	} catch (Exc e) {
		SetConsoleColor(CONSOLE_COLOR::LTYELLOW);
		CoutX() << "\n" << "Error: " << e;
		SetConsoleColor(CONSOLE_COLOR::PREVIOUS);
		return ret = "";
	}
	if (IsVoid(ret))
		return "";
	return ret;
}

double DLL_DemoVectorPyC(const double *v, int num) noexcept {
    double res = 0;
    for (int i = 0; i < num; ++i) 
        res += v[i];
    return res;
}


//#endif
//#endif

//#endif