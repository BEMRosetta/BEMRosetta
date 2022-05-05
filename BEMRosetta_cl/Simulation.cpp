// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "Simulation.h"

Simulation::~Simulation() {
	if (out.GetParamCount() > 0)
		out.Save(AppendFileNameX(folder, "DLLData.out"));
}

void Simulation::Load(const String &datfile, double _rho, double _g, double c0x, double c0y, double c0z) {
	rho = _rho;
	g = _g;
	
	String strFile = LoadFile(datfile);
	if (IsNull(strFile)) 
		throw Exc(Format("File '%s' not found", datfile));
		
	String strcalc = GetFASTVar(strFile, "Hydrostatics", "");
	if (IsNull(strcalc) || strcalc == "Static") 
		calculation = 0;
	else if (strcalc == "DynamicStatic")
		calculation = 1;
	else 
		throw Exc(Format("Unknown Hydrostatics '%s'", strcalc));
	
	folder = GetFileFolder(datfile);
	
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
	Point3D cb = mesh.under.GetCenterOfBuoyancy();
	Force6D forceb0;
	mesh.under.GetHydrostaticForceCB(forceb0, mesh.c0, cb, rho, g);
	forceb = forceb0.ToVector();		
	
	UVector<UVector<String>> sforces = GetFASTArray(strFile, "NCoupledForces");
	for (int i = 0; i < sforces.size(); ++i) {
		if (sforces[i].size() < 9) 
			throw Exc("Coupled forces have to include six components");
		auto &f = coupledForces.Add();
		f.Set(ScanDouble(sforces[i][0]), ScanDouble(sforces[i][1]), ScanDouble(sforces[i][2]),
			  ScanDouble(sforces[i][3]), ScanDouble(sforces[i][4]), ScanDouble(sforces[i][5]), 
			  ScanDouble(sforces[i][6]), ScanDouble(sforces[i][7]), ScanDouble(sforces[i][8]));
	}

	UVector<UVector<String>> sfixed = GetFASTArray(strFile, "NFixedForces");
	for (int i = 0; i < sfixed.size(); ++i) {
		if (sfixed[i].size() < 9) 
			throw Exc("\nFixed forces have to include six components");
		auto &f = fixedForces.Add();
		f.Set(ScanDouble(sfixed[i][0]), ScanDouble(sfixed[i][1]), ScanDouble(sfixed[i][2]),
			  ScanDouble(sfixed[i][3]), ScanDouble(sfixed[i][4]), ScanDouble(sfixed[i][5]), 
			  ScanDouble(sfixed[i][6]), ScanDouble(sfixed[i][7]), ScanDouble(sfixed[i][8]));
	}
			
	dampingCentre = GetFASTMatrix(strFile, "DampingCentre", 1, 3).row(0);
	linearDamping = GetFASTMatrix(strFile, "AddBLin", 6, 6);
	if (linearDamping.isZero(0))
		linearDamping.resize(0,0);
	quadDamping   = GetFASTMatrix(strFile, "AddBQuad", 6, 6);
	if (quadDamping.isZero(0))
		quadDamping.resize(0,0);
	addedMass     = GetFASTMatrix(strFile, "AddAmass", 6, 6);
	if (addedMass.isZero(0))
		addedMass.resize(0,0);
	
	auto AddVar = [&] (String name, String units, bool force = false) -> int {
		if (force)
			return out.AddParam(name, units);
		else {
			String res = ToLower(GetFASTVar(strFile, name, "OUTPUT CHANNELS"));
			if (!IsNull(res) && res == "true") {
				return out.AddParam(name, units);
			}
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
	
	output.PtfmTVxi = AddVar("PtfmTVxi", "m/s");
	output.PtfmTVyi = AddVar("PtfmTVyi", "m/s");
	output.PtfmTVzi = AddVar("PtfmTVzi", "m/s");
	output.PtfmRVxi = AddVar("PtfmRVxi", "deg/s");
	output.PtfmRVyi = AddVar("PtfmRVyi", "deg/s");
	output.PtfmRVzi = AddVar("PtfmRVzi", "deg/s");
	
	output.PtfmTAxi = AddVar("PtfmTAxi", "m/s^2");
	output.PtfmTAyi = AddVar("PtfmTAyi", "m/s^2");
	output.PtfmTAzi = AddVar("PtfmTAzi", "m/s^2");
	output.PtfmRAxi = AddVar("PtfmRAxi", "deg/s^2");
	output.PtfmRAyi = AddVar("PtfmRAyi", "deg/s^2");
	output.PtfmRAzi = AddVar("PtfmRAzi", "deg/s^2");
	
}

Force6D Simulation::CalcForces_Static(const float *pos) {
	VectorXd p = C6ToVector(pos);

	VectorXd res = forceb - stiff*p;
	return Force6D(res);
}

Force6D Simulation::CalcForces_Dynamic(const float *pos, double volTolerance) {
	mesh.Move(pos, rho, g, false);

	Force6D b;
	Point3D cb = mesh.under.GetCenterOfBuoyancy();
	mesh.under.GetHydrostaticForceCB(b, mesh.c0, cb, rho, g);
	
	if (mesh.under.VolumeMatch(volTolerance, volTolerance) == -2) 
		throw Exc("Error: Mesh opened in the waterline");

	return b;
}

Force6D Simulation::CalcForces(double time, const float *pos, const float *vel, const float *acc) {
	static double lasttime = -1;
	
	Force6D f6;
	f6.Reset();
	
	if (lasttime == -1 || lasttime != time) {
		out.AddVal(output.time,      time);
		
		out.AddVal(output.ptfmSurge, pos[0]);
		out.AddVal(output.ptfmSway,  pos[1]);
		out.AddVal(output.ptfmHeave, pos[2]);
		out.AddVal(output.ptfmRoll,  ToDeg(pos[3]));
		out.AddVal(output.ptfmPitch, ToDeg(pos[4]));
		out.AddVal(output.ptfmYaw,   ToDeg(pos[5]));
		
		out.AddVal(output.PtfmTVxi,  vel[0]);
		out.AddVal(output.PtfmTVyi,  vel[1]);
		out.AddVal(output.PtfmTVzi,  vel[2]);
		out.AddVal(output.PtfmRVxi,  ToDeg(vel[3]));
		out.AddVal(output.PtfmRVyi,  ToDeg(vel[4]));
		out.AddVal(output.PtfmRVzi,  ToDeg(vel[5]));
		
		out.AddVal(output.PtfmTAxi,  acc[0]);
		out.AddVal(output.PtfmTAyi,  acc[1]);
		out.AddVal(output.PtfmTAzi,  acc[2]);
		out.AddVal(output.PtfmRAxi,  ToDeg(acc[3]));
		out.AddVal(output.PtfmRAyi,  ToDeg(acc[4]));
		out.AddVal(output.PtfmRAzi,  ToDeg(acc[5]));
		
		lasttime = time;
	}
	
	if (fixedForces.size() > 0) {
		for (int i = 0; i < fixedForces.size(); ++i) 
			f6.Add(fixedForces[i], mesh.c0);
	}
	if (coupledForces.size() > 0) {
		for (int i = 0; i < coupledForces.size(); ++i) {
			ForceVector force = clone(coupledForces[i]);	
			force.TransRot(pos[0], pos[1], pos[2], 
						   pos[3], pos[4], pos[5], 
						   mesh.c0.x, mesh.c0.y, mesh.c0.z);
			f6.Add(force, mesh.c0);
		}
	}
	if (linearDamping.size() == 36 || quadDamping.size() == 36) {
		Velocity6D v;
		v.Set(vel);
		v.Translate(mesh.c0, dampingCentre);	
		VectorXd vvel = v.ToVector();
		VectorXd vfdamp = VectorXd::Zero(6);
		if (linearDamping.size() == 36) 
			vfdamp = -linearDamping*vvel;		// Damping is negative
		if (quadDamping.size() == 36) {
			VectorXd vvel2(6);
			for (int i = 0; i < 6; ++i)			// vvel^2
				vvel2[i] = vvel[i]*abs(vvel[i]);// Always with the direction of speed
			vfdamp -= quadDamping*vvel2;		// Damping is negative
		}
		for (int i = 0; i < 6; ++i) {			
			if (Sign(vfdamp[i]) == Sign(vvel[i]))// Negative damping is forgiven
				vfdamp[i] = 0;
		}
		f6.Add(Force6D(vfdamp), dampingCentre, mesh.c0);
	}
	if (addedMass.size() == 36) {
		Velocity6D v;
		v.Set(vel);
		Acceleration6D a;
		a.Set(acc);
		a.Translate(mesh.c0, dampingCentre, v);	
		VectorXd vacc = a.ToVector();
		f6.Add(Force6D(addedMass*vacc), dampingCentre, mesh.c0);
	}
	return f6;
}
