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
	BEM::Print(F("\n%s", str)); 
	return 0;
};

Time Orca::startCalc = Null, Orca::beginNoLicense = Null, Orca::lastLog;
int64 Orca::noLicenseTime = 0;

#endif


//#if defined(flagBEMR_DLL) || defined(flagBEMR_TEST_BMR_INTERNAL) || defined(flagBEMR_TEST_DLL) || defined(flagBEMR_CL)

#include "FastOut.h"
#include "export.h"

BMR_Data &BMR() {
	static BMR_Data dll;
	return dll;
}

const char *BMR_GetLastError() noexcept {
	if (BMR().errorStr.IsEmpty())
		return nullptr;
	return BMR().errorStr;
}

#ifndef flagBEMR_TEST_DLL

void BMR_Init() noexcept {
	BMR();
	BMR_Bem_Id_Set(0);
	BMR_Mesh_Id_Set(0);
}			
																												
void BMR_NoPrint() noexcept {
	CoutStreamX::NoPrint();
}
	
const char *BMR_Version() noexcept {
	static String version;
	version << __DATE__ << ", " << __TIME__;
	return version;	
}

void BMR_Echo(const char *str) noexcept {
	try {
		CoutX() << str;
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}


void BMR_Mesh_Load(const char *file) noexcept {
	try {
		if (!FileExists(file))
			throw Exc(F(t_("File '%s' not found"), file)); 
								
		Bem().LoadBody(file, BMR().echo ? BMR().Status : BMR().NoPrint, false, false, BMR().meshid);
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Mesh_Report() noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Bem().surfs[BMR().meshid].Report(Bem().rho);
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Mesh_Clear() noexcept {
	try {
		Bem().surfs.Clear();
		BMR_Mesh_Id_Set(0);
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Mesh_Id_Set(int id) noexcept {
	try {
		if (IsNull(id) || id < 0)
			throw Exc(F(t_("Invalid id %d"), id));
		
		if (id >= Bem().surfs.size())
			Bem().surfs.SetCount(id+1);
		
		BMR().meshid = id;
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

int BMR_Mesh_Id_Get() noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		
		return BMR().meshid;
	} catch(Exc err) {
		BMR().errorStr = err;
		return -1;
	}
	BMR().errorStr.Clear();
}


void BMR_Mesh_Save(const char *file, const char *format, int symX, int symY) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		
		Body::MESH_FMT meshFmt = Body::GetCodeBodyStr(format);
		if (!Body::meshInfo[meshFmt].canSave)
			throw Exc(F(t_("Saving format '%s' is not implemented"), format));						
		
		Body::SaveAs(Bem().surfs[BMR().meshid], file, meshFmt, Body::ALL, Bem().rho, Bem().g, symX, symY);					
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Mesh_Translate(double x, double y, double z) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		msh.dt.mesh.Translate(x, y, z);
		msh.dt.cg.Translate(x, y, z);
		msh.AfterLoad(Bem().rho, Bem().g, false, false);	
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Mesh_Rotate(double ax, double ay, double az, double cx, double cy, double cz) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		msh.dt.mesh.Rotate(ToRad(ax), ToRad(ay), ToRad(az), cx, cy, cz);	
		msh.dt.cg.Rotate(ToRad(ax), ToRad(ay), ToRad(az), cx, cy, cz);
		msh.AfterLoad(Bem().rho, Bem().g, false, false);	
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Mesh_Cg_Set(double x, double y, double z) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		msh.dt.cg = Point3D(x, y, z);
		msh.AfterLoad(Bem().rho, Bem().g, true, false);
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Mesh_C0_Set(double x, double y, double z) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		msh.dt.c0 = Point3D(x, y, z);
		msh.AfterLoad(Bem().rho, Bem().g, true, false);
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Mesh_Mass_Set(double mass) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		msh.SetMass(mass);
		msh.AfterLoad(Bem().rho, Bem().g, true, false);
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Mesh_Inertia_Set(const double *data, const int dim[2]) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
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
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Mesh_LinearDamping_Set(const double *data, const int dim[2]) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
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
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Mesh_MooringStiffness_Set(const double *data, const int dim[2]) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
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
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Mesh_Reset() noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Bem().surfs[BMR().meshid].Reset(Bem().rho, Bem().g);
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

int BMR_Mesh_Duplicate() noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		
		Body &msh = Bem().surfs[BMR().meshid];
		BMR_Mesh_Id_Set(Bem().surfs.size());
		Bem().surfs[BMR().meshid] = clone(msh);
	} catch(Exc err) {
		BMR().errorStr = err;
		return -1;
	}
	BMR().errorStr.Clear();
	return BMR().meshid;
}

int BMR_Mesh_GetWaterPlane() noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Bem().AddWaterSurface(BMR().meshid, 'e', 1, false);
		BMR().meshid = Bem().surfs.size() - 1;
	} catch(Exc err) {
		BMR().errorStr = err;
		return -1;
	}
	BMR().errorStr.Clear();
	return BMR().meshid;
}

int BMR_Mesh_GetHull() noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Bem().AddWaterSurface(BMR().meshid, 'r', 1, false);
		BMR().meshid = Bem().surfs.size() - 1;
	} catch(Exc err) {
		BMR().errorStr = err;
		return -1;
	}
	BMR().errorStr.Clear();
	return BMR().meshid;
}

int BMR_Mesh_FillWaterplane(double ratio, int quads) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		if (ratio < 0 || ratio > 100)
			throw Exc(F(t_("Wrong mesh ratio %s"), ratio));
		Bem().AddWaterSurface(BMR().meshid, 'f', ratio, quads);
		BMR().meshid = Bem().surfs.size() - 1;
	} catch(Exc err) {
		BMR().errorStr = err;
		return -1;
	}
	BMR().errorStr.Clear();
	return BMR().meshid;
}
		
int BMR_Mesh_GetControlSurface(double distance, double ratio, int quads) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		if (ratio < 0 || ratio > 100)
			throw Exc(F(t_("Wrong mesh ratio %s"), ratio));
		UVector<int> ids;
		ids << BMR().meshid;
		Bem().GetCS(ids, distance, ratio, quads);
		BMR().meshid = Bem().surfs.size() - 1;
	} catch(Exc err) {
		BMR().errorStr = err;
		return -1;
	}
	BMR().errorStr.Clear();
	return BMR().meshid;
}

				
void BMR_Mesh_Volume_Get(double *vx, double *vy, double *vz) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		*vx = msh.dt.mesh.volumex;
		*vy = msh.dt.mesh.volumey;
		*vz = msh.dt.mesh.volumez;	
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Mesh_UnderwaterVolume_Get(double *vx, double *vy, double *vz) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		*vx = msh.dt.under.volumex;
		*vy = msh.dt.under.volumey;
		*vz = msh.dt.under.volumez;
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

double BMR_Mesh_Surface_Get() noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
	} catch(Exc err) {
		BMR().errorStr = err;
		return Null;
	}
	BMR().errorStr.Clear();
	return Bem().surfs[BMR().meshid].dt.mesh.surface;
}

double BMR_Mesh_UnderwaterSurface_Get() noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
	} catch(Exc err) {
		BMR().errorStr = err;
		return Null;
	}
	BMR().errorStr.Clear();
	return Bem().surfs[BMR().meshid].dt.under.surface;
}

void BMR_Mesh_Centre_Volume_Get(double *x, double *y, double *z) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		Point3D cg = msh.dt.mesh.GetCentreOfBuoyancy();
		*x = cg.x;
		*y = cg.y;
		*z = cg.z;
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Mesh_Centre_Surface_Get(double *x, double *y, double *z) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		Point3D cg = msh.dt.mesh.GetCentreOfGravity_Surface();
		*x = cg.x;
		*y = cg.y;
		*z = cg.z;
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Mesh_HydrostaticStiffness_Get(double **data, int dim[2]) noexcept {
	static UVector<double> d;
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		dim[0] = (int)msh.dt.C.rows();
		dim[1] = (int)msh.dt.C.cols();
		CopyRowMajor(msh.dt.C, d);
		*data = d.begin();
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Mesh_NumPanels_Get(int *num) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		
		*num = msh.dt.mesh.panels.size();
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Mesh_VolumeEnvelope_Get(double *minx, double *maxx, double *miny, double *maxy, double *minz, double *maxz) noexcept {
	try {
		if (Bem().surfs.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		Body &msh = Bem().surfs[BMR().meshid];
		
		*minx = msh.dt.mesh.env.minX;
		*maxx = msh.dt.mesh.env.maxX;
		*miny = msh.dt.mesh.env.minY;
		*maxy = msh.dt.mesh.env.maxY;
		*minz = msh.dt.mesh.env.minZ;
		*maxz = msh.dt.mesh.env.maxZ;
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Bem_Clear() noexcept {
	Bem().hydros.Clear();
	BMR_Bem_Id_Set(0);
}

int BMR_Bem_New() noexcept {
	BMR_Bem_Id_Set(Bem().hydros.size());
	return BMR().bemid;
}

void BMR_Bem_Id_Set(int id) noexcept {
	try {
		if (IsNull(id) || id < 0)
			throw Exc(F(t_("Invalid id %d"), id));
		
		if (id >= Bem().hydros.size())
			Bem().hydros.SetCount(id+1);
		
		BMR().bemid = id;
		BMR_Bem_Body_Id_Set(0);
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

int BMR_Bem_Id_Get() noexcept {
	try {
		if (Bem().hydros.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		
		return BMR().bemid;
	} catch(Exc err) {
		BMR().errorStr = err;
		return -1;
	}
	BMR().errorStr.Clear();
}

void BMR_Bem_Support(const char *solver, int *irregular, int *autoIrregular, int *middle7, int *far8, int *near9, int *autoCS, int *multibody) noexcept {
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
		}
	}
}

void BMR_Bem_depth_Set(double h) noexcept {
	try {
		if (Bem().hydros.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		
		if (IsNull(h))
			throw Exc(F(t_("Wrong depth '%f'"), h));
		
		Bem().hydros[BMR().bemid].dt.h = h;
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Bem_g_Set(double g) noexcept {
	try {
		if (Bem().hydros.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		
		if (IsNull(g) || g < 0)
			throw Exc(F(t_("Wrong gravity '%f'"), g));
		
		Bem().hydros[BMR().bemid].dt.g = g;
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Bem_rho_Set(double rho) noexcept {
	try {
		if (Bem().hydros.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		
		if (IsNull(rho) || rho < 0)
			throw Exc(F(t_("Wrong density '%f'"), rho));
		
		Bem().hydros[BMR().bemid].dt.rho = rho;
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Bem_w_Set(const double *w, int dim) noexcept {
	try {
		if (Bem().hydros.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		
		Hydro &hy = Bem().hydros[BMR().bemid];
		
		hy.dt.w.SetCount(dim);
		Copy(w, dim, hy.dt.w);
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
		hy.dt.Apan = Eigen::Tensor<double, 5>();
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Bem_headings_Set(const double *head, int dim) noexcept {
	try {
		if (Bem().hydros.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		
		Hydro &hy = Bem().hydros[BMR().bemid];
		
		hy.dt.head.SetCount(dim);
		Copy(head, dim, hy.dt.head);
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
		return;
	}
	BMR().errorStr.Clear();
}

int BMR_Bem_Duplicate() noexcept {
	try {
		if (Bem().hydros.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		
		Hydro &hy = Bem().hydros[BMR().bemid];
		BMR_Bem_Id_Set(Bem().hydros.size());
		Bem().hydros[BMR().bemid] = clone(hy);									
	} catch(Exc err) {
		BMR().errorStr = err;
		return -1;
	}
	BMR().errorStr.Clear();
	return BMR().bemid;
}

void BMR_Bem_Body_Id_Set(int id) noexcept {
	try {
		if (Bem().hydros.IsEmpty()) 
			throw Exc(t_("No file loaded"));
							
		if (IsNull(id) || id < 0)
			throw Exc(F(t_("Invalid id %d"), id));
		
		Hydro &hy = Bem().hydros[BMR().bemid];
		if (id >= hy.dt.msh.size()) {
			hy.dt.msh.SetCount(id+1);
			hy.dt.Nb = hy.dt.msh.size();
		}
		BMR().bembodyid = id;
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

int BMR_Bem_Body_Id_Get() noexcept {
	try {
		if (Bem().hydros.IsEmpty()) 
			throw Exc(t_("No file loaded"));
		
		return BMR().bembodyid;
	} catch(Exc err) {
		BMR().errorStr = err;
		return -1;
	}
	BMR().errorStr.Clear();
}

void BMR_Bem_Body_LoadMesh(const char *file) noexcept {
	try {
		if (Bem().hydros.size() < BMR().bemid) 
			throw Exc(F(t_("Model %d is not set"), BMR().bemid));
		Hydro &hy = Bem().hydros[BMR().bemid];
		hy.dt.msh.SetCount(max(BMR().bembodyid+1, hy.dt.msh.size()));
		Body::Load(hy.dt.msh[BMR().bembodyid], file, Bem().rho, Bem().g, Null, Null, false);
		hy.dt.Nb = hy.dt.msh.size();
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Bem_Body_LoadLid(const char *file) noexcept {
	try {
		if (Bem().hydros.size() < BMR().bemid) 
			throw Exc(F(t_("Model %d is not set"), BMR().bemid));
		Hydro &hy = Bem().hydros[BMR().bemid];
		hy.dt.lids.SetCount(max(BMR().bembodyid+1, hy.dt.lids.size()));
		Body::Load(hy.dt.lids[BMR().bembodyid], file, Bem().rho, Bem().g, Null, Null, false);
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Bem_Body_LoadControlSurface(const char *file) noexcept {
	try {
		if (Bem().hydros.size() < BMR().bemid) 
			throw Exc(F(t_("Model %d is not set"), BMR().bemid));
		Hydro &hy = Bem().hydros[BMR().bemid];
		hy.dt.css.SetCount(max(BMR().bembodyid+1, hy.dt.css.size()));	
		Body::Load(hy.dt.css[BMR().bembodyid], file, Bem().rho, Bem().g, Null, Null, false);
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Bem_Body_LoadMeshFromMesh(int idsurf) noexcept {
	try {
		if (Bem().hydros.size() < BMR().bemid) 
			throw Exc(F(t_("Model %d is not set"), BMR().bemid));
		Hydro &hy = Bem().hydros[BMR().bemid];
		hy.dt.msh.SetCount(max(BMR().bembodyid+1, hy.dt.msh.size()));
		if (Bem().surfs.size() < idsurf) 
			throw Exc(F(t_("Id %d is not loaded"), idsurf));
		hy.dt.msh[BMR().bembodyid] = clone(Bem().surfs[idsurf]);
		hy.dt.Nb = hy.dt.msh.size();
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Bem_Body_LoadLidFromMesh(int idsurf) noexcept {
	try {
		if (Bem().hydros.size() < BMR().bemid) 
			throw Exc(F(t_("Model %d is not set"), BMR().bemid));
		Hydro &hy = Bem().hydros[BMR().bemid];
		hy.dt.lids.SetCount(max(BMR().bembodyid+1, hy.dt.lids.size()));
		if (Bem().surfs.size() < idsurf) 
			throw Exc(F(t_("Id %d is not loaded"), idsurf));
		hy.dt.lids[BMR().bembodyid] = clone(Bem().surfs[idsurf]);
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Bem_Body_LoadControlSurfaceFromMesh(int idsurf) noexcept {
	try {
		if (Bem().hydros.size() < BMR().bemid) 
			throw Exc(F(t_("Model %d is not set"), BMR().bemid));
		Hydro &hy = Bem().hydros[BMR().bemid];
		hy.dt.css.SetCount(max(BMR().bembodyid+1, hy.dt.css.size()));
		if (Bem().surfs.size() < idsurf) 
			throw Exc(F(t_("Id %d is not loaded"), idsurf));
		hy.dt.css[BMR().bembodyid] = clone(Bem().surfs[idsurf]);
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Bem_Body_C0_Set(double x, double y, double z) noexcept {
	try {
		if (Bem().hydros.size() < BMR().bemid) 
			throw Exc(F(t_("Model %d is not set"), BMR().bemid));
		Hydro &hy = Bem().hydros[BMR().bemid];
		hy.dt.msh[BMR().bembodyid].dt.c0 = Point3D(x, y, z);
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Bem_Body_Cg_Set(double x, double y, double z) noexcept {
	try {
		if (Bem().hydros.size() < BMR().bemid) 
			throw Exc(F(t_("Model %d is not set"), BMR().bemid));
		Hydro &hy = Bem().hydros[BMR().bemid];
		hy.dt.msh[BMR().bembodyid].dt.cg = Point3D(x, y, z);
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}

void BMR_Bem_Body_Inertia_Set(const double *data, const int dim[2]) noexcept {
	try {
		if (Bem().hydros.size() < BMR().bemid) 
			throw Exc(F(t_("Model %d is not set"), BMR().bemid));
		Hydro &hy = Bem().hydros[BMR().bemid];
		if (dim[0] != 6 || dim[1] != 6)
			throw Exc(F(t_("Matrix dimensions (%d,%d) are not correct"), dim[0], dim[1]));
		hy.dt.msh[BMR().bembodyid].dt.M.resize(6, 6);
		for (int r = 0; r < 6; ++r)
			for (int c = 0; c < 6; ++c)
				hy.dt.msh[BMR().bembodyid].dt.M(r, c) = data[r*6 + c];
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();	
}

void BMR_Bem_Body_Name_Set(const char *name) noexcept {
	try {
		if (Bem().hydros.size() < BMR().bemid) 
			throw Exc(F(t_("Model %d is not set"), BMR().bemid));
		Hydro &hy = Bem().hydros[BMR().bemid];
		hy.dt.msh[BMR().bembodyid].dt.name = name;
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();
}



void BMR_Bem_SaveCase(const char *folder, const char *solver, bool x0z, bool y0z, 
		bool irregular, bool autoIrregular, const char *qtfType, bool autoQTF, 
		bool bin, int numCases, int numThreads, bool withPotentials, bool withMesh) noexcept {
	try {
		if (Bem().hydros.size() < BMR().bemid) 
			throw Exc(F(t_("Model %d is not set"), BMR().bemid));
		Hydro &hy = Bem().hydros[BMR().bemid];
		
		UVector<String> candidates;
		int type, icase;
		String lsolver = ToLower(Replace(solver, " ", ""));
		for (type = 0; type < Hydro::NUMBEM; ++type) {
			if (Hydro::bemInfo[type].caseCanSave) {
				String fmt = Hydro::GetBemStrCase(static_cast<Hydro::BEM_FMT>(type));
				if (ToLower(Replace(fmt, " ", "")).Find(ToLower(lsolver)) >= 0) {
					candidates << fmt;
					icase = type;
				}
			}
		}
		if (candidates.IsEmpty())
			throw Exc(F(t_("Unknown format %s"), solver));
		if (candidates.size() > 1) {
			String ret;
			for (int i = 0; i < candidates.size(); ++i) {
				if (i > 0)
					ret << ", ";
				ret << candidates[i];
			}
			throw Exc(F(t_("Format %s can be confused with %s"), solver, ret));
		}
		int iqtfType = -1;
		String sqtfType = ToLower(qtfType);
		if (sqtfType.Find("control") >= 0 || sqtfType.Find("middle") >= 0)
			iqtfType = 7;
		else if (sqtfType.Find("momentum") >= 0 || sqtfType.Find("far") >= 0)
			iqtfType = 8;
		else if (sqtfType.Find("pressure") >= 0 || sqtfType.Find("near") >= 0)
			iqtfType = 9;
		
		UVector<String> errors = hy.Check(static_cast<Hydro::BEM_FMT>(icase));
		if (!errors.IsEmpty()) {
			String str;
			if (errors.size() == 1)
				str << "\n " << errors[0];
			else {
				for (int i = 0; i < errors.size(); ++i)
				 	str << "\n- " << errors[i];
			}
			throw Exc(F(t_("Problems found in data: %s"), str));
		}
		UVector<bool> listDOF(6, true);
		UVector<Point3D> dummy;
		hy.SaveCase(folder, static_cast<Hydro::BEM_FMT>(icase), x0z, y0z, 
				irregular, autoIrregular, iqtfType, autoQTF, 
				bin, numCases, numThreads, withPotentials, withMesh, listDOF, dummy);
	} catch(Exc err) {
		BMR().errorStr = err;
		return;
	}
	BMR().errorStr.Clear();	
}
	
	
int BMR_FAST_Load(const char *filename) noexcept {
	try {
		String ret = BMR().fast.Load(filename, Null);
		if (ret.IsEmpty())
			return 1;
		else {
			CoutX() << F("Error: %s", ret);
			return 0;
		}
	} catch (...) {
		CoutX() << "Unknown error in BMR_FAST_Load()";
		return 0;
	}
}

const char *BMR_FAST_GetParameterName(int id) noexcept {
	static String ret;
	try {
		return ret = BMR().fast.GetParameter(id);
	} catch (...) {
		CoutX() << "Unknown error in BMR_FAST_GetParameterName()";
		return ret = "Error";
	}
}

const char *BMR_FAST_GetUnitName(int id) noexcept {
	static String ret;
	try {
		return ret = BMR().fast.GetUnit(id);
	} catch (...) {
		CoutX() << "Unknown error in BMR_FAST_GetUnitName()";
		return ret = "Error";
	}
}

int BMR_FAST_GetParameterId(const char *name) noexcept {
	try {
		UVector<int> p = BMR().fast.FindParameterMatch(name);
		if (p.IsEmpty())
			return -1;
		else
			return p[0];
	} catch (...) {
		CoutX() << "Unknown error in BMR_FAST_GetParameterCount()";
		return Null;
	}
}

int BMR_FAST_GetParameterCount() noexcept {
	try {
		return BMR().fast.GetParameterCount();
	} catch (...) {
		CoutX() << "Unknown error in BMR_FAST_GetParameterCount()";
		return Null;
	}
}

int BMR_FAST_GetLen() noexcept {
	try {
		return BMR().fast.GetNumData();
	} catch (...) {
		CoutX() << "Unknown error in BMR_FAST_GetLen()";
		return Null;
	}
}

double BMR_FAST_GetTimeStart() noexcept {
	try {
		return BMR().fast.GetTimeStart();
	} catch (...) {
		CoutX() << "Unknown error in BMR_FAST_GetTimeStart()";
		return Null;
	}
}

double BMR_FAST_GetTimeEnd() noexcept {
	try {
		return BMR().fast.GetTimeEnd();
	} catch (...) {
		CoutX() << "Unknown error in BMR_FAST_GetTimeEnd()";
		return Null;
	}
}

double BMR_FAST_GetTime(int idtime) noexcept {
	return BMR_FAST_GetData(idtime, 0);
}

double BMR_FAST_GetData(int idtime, int idparam) noexcept {
	if (idtime < 0) {
		CoutX() << "Error in BMR_FAST_GetData() idtime < 0";
		return Null;
	}
	if (idtime >= BMR().fast.GetNumData()) {
		CoutX() << "Error in BMR_FAST_GetData() idtime >= time";
		return Null;
	}
	if (idparam < 0) {
		CoutX() << "Error in BMR_FAST_GetData() idparam < 0";
		return Null;
	}
	if (idparam >= BMR().fast.GetParameterCount()) {
		CoutX() << "Error in BMR_FAST_GetData() idparam >= num_params";
		return Null;
	}
		
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

int BMR_FAST_GetArray(int idparam, int idbegin, int idend, double **data, int *num) noexcept {
	static VectorXd v;
	
	try {
		BMR_FAST_GetData(idparam, idbegin, idend, v);
		
		*num = int(v.size());
		*data = v.data();
		return 1;
	} catch (Exc e) {
		CoutX() << F("Error in BMR_FAST_GetArray(): %s", e);
	} catch (...) {
		CoutX() << "Unknown error in BMR_FAST_GetArray()";
	}
	return Null;	
}

double BMR_FAST_GetAvg(int idparam, int idbegin, int idend) noexcept {
	try {
		VectorXd data;
		
		BMR_FAST_GetData(idparam, idbegin, idend, data);
		
		return data.mean();
	} catch (Exc e) {
		CoutX() << F("Error in BMR_FAST_GetAvg(): %s", e);
	} catch (...) {
		CoutX() << "Unknown error in BMR_FAST_GetAvg()";
	}
	return Null;
}

double BMR_FAST_GetMax(int idparam, int idbegin, int idend) noexcept {
	try {
		VectorXd data;
		
		BMR_FAST_GetData(idparam, idbegin, idend, data);
		
		return data.maxCoeff();
	} catch (Exc e) {
		CoutX() << F("Error in BMR_FAST_GetMax(): %s", e);
	} catch (...) {
		CoutX() << "Unknown error in BMR_FAST_GetAvg()";
	}
	return Null;
}

double BMR_FAST_GetMin(int idparam, int idbegin, int idend) noexcept {
	try {
		VectorXd data;
		
		BMR_FAST_GetData(idparam, idbegin, idend, data);
		
		return data.minCoeff();
	} catch (Exc e) {
		CoutX() << F("Error in BMR_FAST_GetMin(): %s", e);
	} catch (...) {
		CoutX() << "Unknown error in BMR_FAST_GetAvg()";
	}
	return Null;
}

int BMR_IsNull(double val) noexcept {return IsNull(val);}


int BMR_FAST_LoadFile(const char *file) noexcept {
	BMR().fastFileStr = LoadFile(file);
	BMR().fastFileName = file;
	return !BMR().fastFileStr.IsEmpty();
}

int BMR_FAST_SaveFile(const char *file) noexcept {
	bool ret;
	try {
		String sfile(file);
		if (!sfile.IsEmpty()) 
			BMR().fastFileName = sfile;
		ret = SaveFile(BMR().fastFileName, BMR().fastFileStr);
		
	} catch (Exc e) {
		SetConsoleColor(CONSOLE_COLOR::LTYELLOW);
		CoutX() << "\n" << "Error: " << e;
		SetConsoleColor(CONSOLE_COLOR::PREVIOUS);
		return false;
	}
	return ret;	
}

int BMR_FAST_SetVar(const char *name, const char *paragraph, const char *value) noexcept {
	try {
		SetFASTVar(BMR().fastFileStr, name, value, paragraph);
	} catch (Exc e) {
		SetConsoleColor(CONSOLE_COLOR::LTYELLOW);
		CoutX() << "\n" << "Error: " << e;
		SetConsoleColor(CONSOLE_COLOR::PREVIOUS);
		return false;
	}
	return true;
}

const char *BMR_FAST_GetVar(const char *name, const char *paragraph) noexcept {
	static String ret;

	try {
		ret = GetFASTVar(BMR().fastFileStr, name, paragraph);
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

double BMR_DemoVectorPyC(const double *v, int num) noexcept {
    double res = 0;
    for (int i = 0; i < num; ++i) 
        res += v[i];
    return res;
}


#endif

#if defined(flagBEMR_TEST_DLL_INTERNAL) || defined(flagBEMR_TEST_DLL)
void BEM_Throw() {
	if (BMR_GetLastError())
		throw Exc(BMR_GetLastError());
}
#endif
//#endif

//#endif