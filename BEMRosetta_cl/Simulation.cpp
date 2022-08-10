// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "Simulation.h"


const char *Simulation::strCalculation[] = {t_("None"), t_("Static"), t_("Dynamic-Static"), t_("Dynamic"), t_("OpenFAST"), t_("Unknown")};
	
Simulation::~Simulation() {
	String outFile = AppendFileNameX(folder, "DLLData.out");
	if (out.GetParamCount() > 1)
		out.Save(outFile, ".out");
	else
		DeleteFile(outFile);
}

void Simulation::Load(const String &datfile, int stiffMod, int dllForce, 
			double _rho, double _g, double c0x, double c0y, double c0z, const BEM &bem) {
	rho = _rho;
	g = _g;
	
	String strFile = LoadFile(datfile);
	if (IsNull(strFile)) 
		throw Exc(Format("File '%s' not found", datfile));
	
	folder = GetFileFolder(datfile);
	
	if (stiffMod < 3)
		calculation = OPENFAST;
	else {
		String strcalc = GetFASTVar(strFile, "Hydrostatics", "");
		if (IsNull(strcalc) || strcalc == "Static") 
			calculation = STATIC;
		else if (strcalc == "None")
			calculation = NONE;
		else if (strcalc == "DynamicStatic")
			calculation = DYN_STATIC;
		else if (strcalc == "Dynamic")
			calculation = DYNAMIC;
		else 
			throw Exc(Format("Unknown Hydrostatics '%s'", strcalc));
	
		if (calculation != NONE) {
			String meshfile = GetFASTVar(strFile, "MeshFile", "");
			if (IsNull(meshfile)) 
				throw Exc("Variable 'meshfile' not found");
			meshfile.Replace("\"", "");
			
			String ret = mesh.Load(meshfile, rho, g, false);
			if (!ret.IsEmpty()) 
				throw Exc(Format("Error loading mesh: %s", ret));
		
			mesh.c0.Set(c0x, c0y, c0z);		// Reference system
			mesh.mass = 0;
			mesh.AfterLoad(rho, g, false, false);
			
			stiff = clone(mesh.C);
			cb = mesh.under.GetCenterOfBuoyancy();
			Force6D forceb0 = mesh.under.GetHydrostaticForceCB(mesh.c0, cb, rho, g);
			forceb = forceb0.ToVector();		
		}
	}
	UVector<UVector<String>> sforces = GetFASTArray(strFile, "NCoupledForces");
	for (int i = 0; i < sforces.size(); ++i) {
		if (sforces[i].size() < 9) 
			throw Exc("Coupled forces have to include six components");
		auto &f = coupledForces.Add();
		f.Set(ScanDouble(sforces[i][0]), ScanDouble(sforces[i][1]), ScanDouble(sforces[i][2]),
			  sforces[i][3], sforces[i][4], sforces[i][5], 
			  sforces[i][6], sforces[i][7], sforces[i][8]);
	}

	UVector<UVector<String>> sfixed = GetFASTArray(strFile, "NFixedForces");
	for (int i = 0; i < sfixed.size(); ++i) {
		if (sfixed[i].size() < 9) 
			throw Exc("\nFixed forces have to include six components");
		auto &f = fixedForces.Add();
		f.Set(ScanDouble(sfixed[i][0]), ScanDouble(sfixed[i][1]), ScanDouble(sfixed[i][2]),
			  sfixed[i][3], sfixed[i][4], sfixed[i][5], 
			  sfixed[i][6], sfixed[i][7], sfixed[i][8]);
	}
	
	if (dllForce > 0) {		
		dampingCentre = GetFASTMatrix(strFile, "DampingCentre", 1, 3).row(0);
		linearDamping = GetFASTMatrix(strFile, "AddBLin", 6, 6);
		if (linearDamping.isZero(0))
			linearDamping.resize(0,0);
		quadDamping   = GetFASTMatrix(strFile, "AddBQuad", 6, 6);
		if (quadDamping.isZero(0))
			quadDamping.resize(0,0);
	}
	
	auto AddVar = [&] (String name, String units, bool isTime = false) -> int {
		if (isTime)
			return out.AddParam(name, units);
		else {
			String res = ToLower(GetFASTVar(strFile, name, "OUTPUT CHANNELS"));
			if (!IsNull(res) && res == "true") 
				return out.AddParam(name, units);
		}
		return -1;
	};
		
	output.time = AddVar("time", "s", true);

	output.ptfmSurge= AddVar("ptfmSurge", "m");
	output.ptfmSway = AddVar("ptfmSway",  "m");
	output.ptfmHeave= AddVar("ptfmHeave", "m");
	output.ptfmRoll = AddVar("ptfmRoll",  "deg");
	output.ptfmPitch= AddVar("ptfmPitch", "deg");
	output.ptfmYaw  = AddVar("ptfmYaw",   "deg");
	
	output.wave1Elev= AddVar("wave1Elev", "m");
	
	output.ptfmSurge_cd= AddVar("ptfmSurge_cd", "m");
	output.ptfmSway_cd = AddVar("ptfmSway_cd",  "m");
	output.ptfmHeave_cd= AddVar("ptfmHeave_cd", "m");
	output.ptfmRoll_cd = AddVar("ptfmRoll_cd",  "deg");
	output.ptfmPitch_cd= AddVar("ptfmPitch_cd", "deg");
	output.ptfmYaw_cd  = AddVar("ptfmYaw_cd",   "deg");
	
	output.ptfmTVxi_cd = AddVar("ptfmTVxi_cd", "m/s");
	output.ptfmTVyi_cd = AddVar("ptfmTVyi_cd", "m/s");
	output.ptfmTVzi_cd = AddVar("ptfmTVzi_cd", "m/s");
	output.ptfmRVxi_cd = AddVar("ptfmRVxi_cd", "deg/s");
	output.ptfmRVyi_cd = AddVar("ptfmRVyi_cd", "deg/s");
	output.ptfmRVzi_cd = AddVar("ptfmRVzi_cd", "deg/s");
	
	output.ptfmTAxi_cd = AddVar("ptfmTAxi_cd", "m/s^2");
	output.ptfmTAyi_cd = AddVar("ptfmTAyi_cd", "m/s^2");
	output.ptfmTAzi_cd = AddVar("ptfmTAzi_cd", "m/s^2");
	output.ptfmRAxi_cd = AddVar("ptfmRAxi_cd", "deg/s^2");
	output.ptfmRAyi_cd = AddVar("ptfmRAyi_cd", "deg/s^2");
	output.ptfmRAzi_cd = AddVar("ptfmRAzi_cd", "deg/s^2");
	
	output.ptfmCdFx	= AddVar("ptfmCdFx", "N");
	output.ptfmCdFy = AddVar("ptfmCdFy", "N");
	output.ptfmCdFz = AddVar("ptfmCdFz", "N");
	output.ptfmCdMx = AddVar("ptfmCdMx", "N-m");
	output.ptfmCdMy = AddVar("ptfmCdMy", "N-m");
	output.ptfmCdMz = AddVar("ptfmCdMz", "N-m");
	
	output.ptfmCd2Fx = AddVar("ptfmCd2Fx", "N");
	output.ptfmCd2Fy = AddVar("ptfmCd2Fy", "N");
	output.ptfmCd2Fz = AddVar("ptfmCd2Fz", "N");
	output.ptfmCd2Mx = AddVar("ptfmCd2Mx", "N-m");
	output.ptfmCd2My = AddVar("ptfmCd2My", "N-m");
	output.ptfmCd2Mz = AddVar("ptfmCd2Mz", "N-m");
	
	output.ptfmStiffFx = AddVar("ptfmStiffFx", "N");
	output.ptfmStiffFy = AddVar("ptfmStiffFy", "N");
	output.ptfmStiffFz = AddVar("ptfmStiffFz", "N");
	output.ptfmStiffMx = AddVar("ptfmStiffMx", "N-m");
	output.ptfmStiffMy = AddVar("ptfmStiffMy", "N-m");
	output.ptfmStiffMz = AddVar("ptfmStiffMz", "N-m");
	
	output.ptfmCBx = AddVar("ptfmCBx", "m");
	output.ptfmCBy = AddVar("ptfmCBy", "m");
	output.ptfmCBz = AddVar("ptfmCBz", "m");
	
	output.ptfmVol = AddVar("ptfmVol", "m^3");
}

Force6D Simulation::CalcStiff(double time, const float *pos, double volTolerance, SeaWaves &waves) {
	Force6D f6;
	
	if (calculation == Simulation::NONE)
		f6.Reset();
	else if (calculation == Simulation::STATIC)
		f6 = CalcStiff_Static(pos);
	else if (calculation == Simulation::DYN_STATIC)
		f6 = CalcStiff_DynamicStatic(time, pos, volTolerance);
	else if (calculation == Simulation::DYNAMIC)
		f6 = CalcStiff_Dynamic(time, pos, volTolerance, waves);

	out.SetNextTime(time);

	if (output.ptfmCBx >= 0 || output.ptfmCBy >= 0 || output.ptfmCBz >= 0) {
		if (calculation != Simulation::STATIC)
			cb = mesh.under.GetCenterOfBuoyancy();
	}
	if (output.ptfmVol >= 0) {
		if (calculation != Simulation::STATIC)
			mesh.under.GetVolume();
	}
	if (output.wave1Elev >= 0) {
		if (calculation != Simulation::STATIC)
			out.SetVal(output.wave1Elev, waves.ZSurf(0, 0, time));
	}
	
	out.SetVal(output.ptfmSurge, pos[0]);
	out.SetVal(output.ptfmSway,  pos[1]);
	out.SetVal(output.ptfmHeave, pos[2]);
	out.SetVal(output.ptfmRoll,  ToDeg(pos[3]));
	out.SetVal(output.ptfmPitch, ToDeg(pos[4]));
	out.SetVal(output.ptfmYaw,   ToDeg(pos[5]));
	
	out.SetVal(output.ptfmStiffFx, f6.t.x);
	out.SetVal(output.ptfmStiffFy, f6.t.y);
	out.SetVal(output.ptfmStiffFz, f6.t.z);
	out.SetVal(output.ptfmStiffMx, f6.r.x);
	out.SetVal(output.ptfmStiffMy, f6.r.y);
	out.SetVal(output.ptfmStiffMz, f6.r.z);

	out.SetVal(output.ptfmCBx, 	cb.x);
	out.SetVal(output.ptfmCBy, 	cb.y);
	out.SetVal(output.ptfmCBz,  cb.z);

	out.SetVal(output.ptfmVol,  mesh.under.volume);
			
	return f6;
}

Force6D Simulation::CalcStiff_Static(const float *pos) {
	VectorXd p = C6ToVector(pos);

	VectorXd res = forceb - stiff*p;
	return Force6D(res);
}

Force6D Simulation::CalcStiff_DynamicStatic(double time, const float *pos, double volTolerance) {
	mesh.Move(pos, rho, g, false);

	Point3D p(pos[0], pos[1], pos[2]);
	Force6D f6 = mesh.under.GetHydrostaticForce(p - mesh.c0, rho, g);
	//Point3D cb = mesh.under.GetCenterOfBuoyancy();
	//Force6D f6 = mesh.under.GetHydrostaticForceCB(p - mesh.c0, cb, rho, g);
	
	if (mesh.under.VolumeMatch(volTolerance, volTolerance) == -2) 
		throw Exc("Error: Mesh opened in the waterline");
	
	return f6;
}

Force6D Simulation::CalcStiff_Dynamic(double time, const float *pos, double volTolerance, 
			SeaWaves &waves) {
	mesh.Move(pos, rho, g, false);
	
	bool clip = true;
	
	Point3D p(pos[0], pos[1], pos[2]);	
	Force6D f6 = mesh.mesh.GetHydrodynamicForce(p - mesh.c0, clip,  
		[&](double x, double y)->double {return waves.ZSurf(x, y, time);},
		[&](double x, double y, double z, double et)->double {
			return -waves.Pressure(x, y, waves.ZWheelerStretching(z, et), time);}
	);
	//if (mesh.under.VolumeMatch(volTolerance, volTolerance) == -2) 
	//	throw Exc("Error: Mesh opened in the waterline");
	
	return f6;
}
	
Force6D Simulation::CalcForces(double time, const float *pos, const float *vel, const float *acc) {
	Force6D f6;
	f6.Reset();
	
	Affine3d aff;
	GetTransform(aff, pos[0], pos[1], pos[2], pos[3], pos[4], pos[5], mesh.c0.x, mesh.c0.y, mesh.c0.z);
	
	Point3D meshc0 = clone(mesh.c0);
	Point3D dcentre = clone(dampingCentre);
	
	meshc0.TransRot(aff);
	dcentre.TransRot(aff);
	
	Point3D rpq = dcentre - meshc0; 

	Value6D p(pos);
	p.t += rpq;
	
	Velocity6D v(vel);
	v.Translate(rpq);

	Acceleration6D a(acc);
	a.Translate(rpq, v);			
		
	if (fixedForces.size() > 0) {
		for (int i = 0; i < fixedForces.size(); ++i) {
			ForceVector f(fixedForces[i].p.x, fixedForces[i].p.y, fixedForces[i].p.z,
				fixedForces[i].tx.LinearInterpolate(time),
				fixedForces[i].ty.LinearInterpolate(time),
				fixedForces[i].tz.LinearInterpolate(time),
				fixedForces[i].rx.LinearInterpolate(time),
				fixedForces[i].ry.LinearInterpolate(time),
				fixedForces[i].rz.LinearInterpolate(time));
			f6.Add(f, meshc0);
		}
	}
	if (coupledForces.size() > 0) {
		for (int i = 0; i < coupledForces.size(); ++i) {
			ForceVector f(coupledForces[i].p.x, coupledForces[i].p.y, coupledForces[i].p.z,
				coupledForces[i].tx.LinearInterpolate(time),
				coupledForces[i].ty.LinearInterpolate(time),
				coupledForces[i].tz.LinearInterpolate(time),
				coupledForces[i].rx.LinearInterpolate(time),
				coupledForces[i].ry.LinearInterpolate(time),
				coupledForces[i].rz.LinearInterpolate(time));
			f6.Add(f.TransRot(aff), meshc0);
		}
	}
	if (linearDamping.size() == 36 || quadDamping.size() == 36) {
		VectorXd vvel = v.ToVector();
		VectorXd vfdamp = VectorXd::Zero(6);
		if (linearDamping.size() == 36) 
			vfdamp = -linearDamping*vvel;			// Damping is negative
		if (quadDamping.size() == 36) {
			VectorXd vvel2(6);
			for (int i = 0; i < 6; ++i)				// vvel^2
				vvel2[i] = vvel[i]*abs(vvel[i]);	// Always with the direction of speed
			vfdamp -= quadDamping*vvel2;			// Damping is negative
		}
		for (int i = 0; i < 6; ++i) {			
			if (Sign(vfdamp[i]) == Sign(vvel[i]))	// Negative damping is forgiven
				vfdamp[i] = 0;
		}
		f6.Add(Force6D(vfdamp), dcentre, meshc0);
	}
	
	out.SetNextTime(time);

	out.SetVal(output.ptfmSurge, pos[0]);
	out.SetVal(output.ptfmSway,  pos[1]);
	out.SetVal(output.ptfmHeave, pos[2]);
	out.SetVal(output.ptfmRoll,  ToDeg(pos[3]));
	out.SetVal(output.ptfmPitch, ToDeg(pos[4]));
	out.SetVal(output.ptfmYaw,   ToDeg(pos[5]));
	
	out.SetVal(output.ptfmSurge_cd, p.t.x);
	out.SetVal(output.ptfmSway_cd,  p.t.y);
	out.SetVal(output.ptfmHeave_cd, p.t.z);
	out.SetVal(output.ptfmRoll_cd,  ToDeg(p.r.x));
	out.SetVal(output.ptfmPitch_cd, ToDeg(p.r.y));
	out.SetVal(output.ptfmYaw_cd,   ToDeg(p.r.z));

	out.SetVal(output.ptfmTVxi_cd,  v.t.x);
	out.SetVal(output.ptfmTVyi_cd,  v.t.y);
	out.SetVal(output.ptfmTVzi_cd,  v.t.z);
	out.SetVal(output.ptfmRVxi_cd,  ToDeg(v.r.x));
	out.SetVal(output.ptfmRVyi_cd,  ToDeg(v.r.y));
	out.SetVal(output.ptfmRVzi_cd,  ToDeg(v.r.z));

	out.SetVal(output.ptfmTAxi_cd,  a.t.x);
	out.SetVal(output.ptfmTAyi_cd,  a.t.y);
	out.SetVal(output.ptfmTAzi_cd,  a.t.z);
	out.SetVal(output.ptfmRAxi_cd,  ToDeg(a.r.x));
	out.SetVal(output.ptfmRAyi_cd,  ToDeg(a.r.y));
	out.SetVal(output.ptfmRAzi_cd,  ToDeg(a.r.z));

	out.SetVal(output.ptfmCdFx, 0);
	out.SetVal(output.ptfmCdFy, 0);
	out.SetVal(output.ptfmCdFz, 0);
	out.SetVal(output.ptfmCdMx, 0);
	out.SetVal(output.ptfmCdMy, 0);
	out.SetVal(output.ptfmCdMz, 0);

	out.SetVal(output.ptfmCd2Fx, 0);
	out.SetVal(output.ptfmCd2Fy, 0);
	out.SetVal(output.ptfmCd2Fz, 0);
	out.SetVal(output.ptfmCd2Mx, 0);
	out.SetVal(output.ptfmCd2My, 0);
	out.SetVal(output.ptfmCd2Mz, 0);

	return f6;
}

