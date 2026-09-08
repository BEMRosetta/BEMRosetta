// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2025, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#include <STEM4U/SeaWaves.h>
#include <Surface/Surface.h>
#include <Surface/Translation.h>
#include "orca.h"

String OrcaWave::Load(String _file, double) {
	dt.file = _file;
	dt.name = GetFileTitle(dt.file);
	dt.dimen = true;
	dt.len = 1;
	dt.Nb = Null;
	
	try {
		BEM::Print("\n\n" + F(t_("Loading '%s'"), dt.file));

#ifdef PLATFORM_WIN32
		if (GetFileExt(dt.file) == ".owr")
			Load_OWR();
		else
#endif
			if (!Load_OF_YML())
				if (!Load_OW_YML())
					return t_("Unknown yml format");

		if (IsNull(dt.Nb))
			return t_("No data found");
	
	} catch (Exc e) {
		return e;
	}
	return String();
}

#ifdef PLATFORM_WIN32
void OrcaWave::Load_OWR() {
	String fileName = ForceExtSafer(dt.file, ".owr");
	
	Orca orca;
	
	orca.LoadWaveResults(fileName);
	
	dt.solver = Hydro::ORCAWAVE_OWR;
	dt.dimen = true;
	dt.x_w = dt.y_w = 0;
	dt.Nb = dt.Nf = dt.Nh = Null;
	
	orca.LoadParameters(*this, Null);
}
#endif

bool OrcaWave::Load_OW_YML() {
	String fileName = ForceExtSafer(dt.file, ".yml");
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		throw Exc(in.Str() + "\n" + t_("File not found or blocked"));
	
	dt.solver = Hydro::ORCAWAVE_YML;

	dt.dimen = true;
	
	dt.x_w = dt.y_w = 0;
	
	dt.Nf = dt.Nh = Null;
	dt.Nb = 0; 
	
	dt.g = 9.80665;		// Default value used when SI units

	dt.symX = dt.symY = false;
	
	char dataFrom = 'r';
	
	int ib = -1;
	Point3D bodyMeshPosition;
	Value3D bodyMeshAttitude;
	bool originCM = false;
	double mass;
	Matrix3d inertia;
	String bodyname;
	
	OrcaFactors factor;	

	YmlParser fy(in);

	FileInLine::Pos fpos = in.GetPos();
	
	auto GetMsh = [&]()->Body& {
		if (ib < 0)
			throw Exc(in.Str() + "\n" + t_("BodyName not found"));
		return dt.msh[ib];
	};
	auto GetCss = [&]()->Body& {
		if (ib < 0)
			throw Exc(in.Str() + "\n" + t_("BodyName not found"));
		return dt.css[ib];
	};
		
	while(fy.GetLine()) {
		if (fy.FirstIs("UnitsSystem")) {
			if (fy.GetVal() == "SI") {
				factor.mass = factor.force = 1000;
				factor.len = 1;	
				factor.Update();
			} else if (fy.GetVal() == "User") 
				;
			else
				throw Exc(in.Str() + "\n" + F(t_("Only SI and User units are supported. Read '%s'"), fy.GetVal()));
		} else if (fy.FirstIs("LengthUnits")) {
			factor.len = FactorLen(fy.GetVal());
			factor.Update();
		} else if (fy.FirstIs("MassUnits")) {
			factor.mass = FactorMass(fy.GetVal());
			factor.Update();
		} else if (fy.FirstIs("ForceUnits")) {
			factor.force = FactorForce(fy.GetVal());
			factor.Update();
		} else if (fy.FirstIs("g")) 
			dt.g = fy.GetDouble()*factor.len;
		else if (fy.FirstIs("WaterDepth")) 
			dt.h = fy.GetDouble()*factor.len;
		else if (fy.FirstIs("WaterDensity")) 
			dt.rho = fy.GetDouble()*factor.rho;
		else if (fy.FirstIs("WavesReferredToBy")) {	
			String val = fy.GetVal();
			if (val.Find("(rad/s)") >= 0)
				dataFrom = 'r';
			else if (val.Find("(s)") >= 0)
				dataFrom = 's';
			else if (val.Find("(Hz)") >= 0)
				dataFrom = 'h';
			else
				throw Exc(in.Str() + "\n"  + F(t_("Unknown data in WavesReferredToBy: %s"), val));
		} else if (fy.FirstIs("PeriodOrFrequency")) {
			const UVector<UVector<String>> &mat = fy.GetMatrix();
			dt.w.SetCount(mat.size());
			for (int i = 0; i < mat.size(); ++i)
				dt.w[i] = ScanDouble(mat[i][0]);
			dt.Nf = mat.size();		
		} else if (fy.FirstIs("WaveHeading")) {
			const UVector<UVector<String>> &mat = fy.GetMatrix();
			dt.head.SetCount(mat.size());
			for (int i = 0; i < mat.size(); ++i)
				dt.head[i] = ScanDouble(mat[i][0]);
			dt.Nh = mat.size();		
		} else if (fy.FirstIs("Bodies")) {
			if (fy.FirstIs("BodyName")) {
				if (fy.Index() != dt.Nb)
					throw Exc(in.Str() + "\n" + t_("Failed body count"));
				ib = dt.Nb;
				dt.msh.SetCount(dt.Nb+1);
				bodyname = fy.GetVal();
				dt.Nb++;
			} else if (fy.FirstIs("BodyCentreOfMass")) {
				UVector<double> line = fy.GetVectorDouble();
				if (line.size() != 3)
					throw Exc(in.Str() + "\n" + t_("Incorrect BodyCentreOfMass"));
						
				GetMsh().dt.cg.x = line[0]*factor.len;
				GetMsh().dt.cg.y = line[1]*factor.len;
				GetMsh().dt.cg.z = line[2]*factor.len;
				
				if (IsNull(GetMsh().dt.c0))
					GetMsh().dt.c0 = clone(GetMsh().dt.cg);
			} else if (fy.FirstIs("BodyMeshPosition")) {
				UVector<double> line = fy.GetVectorDouble();
				if (line.size() != 3)
					throw Exc(in.Str() + "\n" + t_("Incorrect BodyMeshPosition"));
						
				bodyMeshPosition.x = line[0]*factor.len;
				bodyMeshPosition.y = line[1]*factor.len;
				bodyMeshPosition.z = line[2]*factor.len;
				GetMsh().dt.c0 = bodyMeshPosition;
			} else if (fy.FirstIs("BodyMeshAttitude")) {
				UVector<double> line = fy.GetVectorDouble();
				if (line.size() != 3)
					throw Exc(in.Str() + "\n" + t_("Incorrect BodyMeshAttitude"));
						
				bodyMeshAttitude.x = ToRad(line[0]);
				bodyMeshAttitude.y = ToRad(line[1]);
				bodyMeshAttitude.z = ToRad(line[2]);
			} else if (fy.FirstIs("BodyMeshFileName")) {
				String meshname = fy.GetVal();
				if (!FileExists(meshname))
					meshname = AFX(GetFileFolder(fileName), meshname);
			
				Body::Load(GetMsh(), meshname, rho_ndim(), g_ndim(), Null, Null, false);
				GetMsh().dt.fileName = meshname;
				GetMsh().dt.name = bodyname;
				GetMsh().Translate(bodyMeshPosition);
				GetMsh().Rotate(bodyMeshAttitude, bodyMeshPosition);
				GetMsh().AfterLoad(rho_ndim(), g_ndim(), false, true);
			} else if (fy.FirstIs("BodyControlSurfaceMeshFileName")) {
				dt.css.SetCount(dt.Nb);
				String meshname = fy.GetVal();
				if (!FileExists(meshname))
					meshname = AFX(GetFileFolder(fileName), meshname);
			
				Body::Load(GetCss(), meshname, rho_ndim(), g_ndim(), Null, Null, false);
				GetCss().dt.fileName = meshname;
				GetCss().dt.name = bodyname;
				GetCss().Translate(bodyMeshPosition);
				GetCss().Rotate(bodyMeshAttitude, bodyMeshPosition);
				GetCss().AfterLoad(rho_ndim(), g_ndim(), false, true);	
			} else if (fy.FirstIs("BodyMeshSymmetry")) {	
				String sym = fy.GetVal(); 
				//dt.symY = sym.Find("xz") >= 0;		// mesh file already includes symmetry. If not it would be deployed twice
				//dt.symX = sym.Find("yz") >= 0;
			} else if (fy.FirstIs("BodyCentreOfMassZRelativeToFreeSurface")) {
				double z = fy.GetDouble()*factor.len;
				GetMsh().dt.cg.x = bodyMeshPosition.x;
				GetMsh().dt.cg.y = bodyMeshPosition.y;
				GetMsh().dt.cg.z = z;
			} else if (fy.FirstIs("BodyUserOrigin")) {
				UVector<double> line = fy.GetVectorDouble();
				if (line.size() != 3)
					throw Exc(in.Str() + "\n" + t_("Incorrect BodyUserOrigin"));
				
				GetMsh().dt.c0.x = line[0]*factor.len;
				GetMsh().dt.c0.y = line[1]*factor.len;
				GetMsh().dt.c0.z = line[2]*factor.len;	
			} else if (fy.FirstMatch("BodyRadiiOfGyrationRx*")) {	
				UVector<UVector<double>> rad = fy.GetMatrixDouble(true);
				if (rad.size() != 3)
					throw Exc(in.Str() + "\n" + t_("Incorrect BodyRadiiOfGyration matrix"));
				
				GetMsh().SetMass(GetMsh().dt.under.volume*dt.rho);
				
				MatrixXd &M = GetMsh().dt.M;
				
				for (int r = 0; r < 3; ++r) {
					if (rad[r].size() != 3)
						throw Exc(in.Str() + "\n" + t_("Incorrect BodyRadiiOfGyration matrix"));
					for (int c = 0; c < 3; ++c) 
						M(r+3, c+3) = M(0, 0)*sqr(rad[r][c]*factor.len);
				}				
				
			} else if (fy.FirstMatch("BodyExternalStiffnessMatrixx*")) {
				UVector<UVector<double>> mat = fy.GetMatrixDouble(true);
				if (mat.size() != 6)
					throw Exc(in.Str() + "\n" + t_("Incorrect BodyExternalStiffnessMatrix"));
				GetMsh().dt.Cadd.resize(6, 6);
				
				for (int r = 0; r < 6; ++r) {
					if (mat[r].size() != 6)
						throw Exc(in.Str() + "\n" + t_("Incorrect BodyExternalStiffnessMatrix"));
					for (int c = 0; c < 6; ++c) 
						GetMsh().dt.Cadd(r, c) = mat[r][c]*factor.C(r, c);
				}
			} else if (fy.FirstMatch("BodyExternalDampingMatrixx*")) {
				UVector<UVector<double>> mat = fy.GetMatrixDouble(true);
				if (mat.size() != 6)
					throw Exc(in.Str() + "\n" + t_("Incorrect BodyExternalDampingMatrix"));
				GetMsh().dt.Dlin.resize(6, 6);
				
				for (int r = 0; r < 6; ++r) {
					if (mat[r].size() != 6)
						throw Exc(in.Str() + "\n" + t_("Incorrect BodyExternalDampingMatrix"));
					for (int c = 0; c < 6; ++c)
						GetMsh().dt.Dlin(r, c) = mat[r][c]*factor.Dlin(r, c);
				}
			} else if (fy.FirstIs("BodyMass"))
				mass = fy.GetDouble();
			else if (fy.FirstMatch("BodyInertiaTensorRx*")) {
				UVector<UVector<double>> mat = fy.GetMatrixDouble(true);
				if (mat.size() != 3)
					throw Exc(in.Str() + "\n" + t_("Incorrect BodyInertiaTensorR"));
				for (int r = 0; r < 3; ++r) {
					if (mat[r].size() != 3)
						throw Exc(in.Str() + "\n" + t_("Incorrect BodyInertiaTensorR"));
					for (int c = 0; c < 3; ++c) 
						inertia(r, c) = mat[r][c]*factor.M(r, c);
				}
			} else if (fy.FirstIs("BodyInertiaTensorOriginType")) {	
				String origin = fy.GetVal();
				if (origin == "User specified")
					;	// Wait for BodyInertiaTensorUserOrigin
				else if (origin == "Centre of mass") {
					GetMsh().dt.M = Eigen::MatrixXd::Zero(6, 6);
					GetMsh().dt.M(0, 0) = GetMsh().dt.M(1, 1) = GetMsh().dt.M(2, 2) = fy.GetDouble();
					GetMsh().dt.M.block(3, 3, 3, 3) = inertia;
				} else if (origin == "Body origin") {
					Surface::TranslateInertia33(inertia, mass, GetMsh().dt.cg, GetMsh().dt.c0, GetMsh().dt.cg);	// Assemble inertia at CG
					GetMsh().dt.M = Eigen::MatrixXd::Zero(6, 6);
					GetMsh().dt.M(0, 0) = GetMsh().dt.M(1, 1) = GetMsh().dt.M(2, 2) = mass;
					GetMsh().dt.M.block(3, 3, 3, 3) = inertia;	
					Surface::TranslateInertia66(GetMsh().dt.M, GetMsh().dt.cg, GetMsh().dt.cg, GetMsh().dt.c0);	// Translate it to C0
				} else
					throw Exc(in.Str() + "\n" + t_("Incorrect BodyInertiaTensorOriginType"));
			} else if (fy.FirstIs("BodyInertiaTensorUserOrigin")) {
				UVector<double> line = fy.GetVectorDouble();
				if (line.size() != 3)
					throw Exc(in.Str() + "\n" + t_("Incorrect BodyInertiaTensorUserOrigin"));
				Vector3d ref(line[0], line[1], line[2]);
				Surface::TranslateInertia33(inertia, mass, GetMsh().dt.cg, ref, GetMsh().dt.cg);				// Assemble inertia at CG
				GetMsh().dt.M = Eigen::MatrixXd::Zero(6, 6);
				GetMsh().dt.M(0, 0) = GetMsh().dt.M(1, 1) = GetMsh().dt.M(2, 2) = mass;
				GetMsh().dt.M.block(3, 3, 3, 3) = inertia;	
				Surface::TranslateInertia66(GetMsh().dt.M, GetMsh().dt.cg, GetMsh().dt.cg, GetMsh().dt.c0);	// Translate it to C0
			} else if (fy.FirstIs("BodyExternalDampingMatrixOriginType")) {	
				String origin = fy.GetVal();
				if (origin == "User specified")
					;	// Wait for BodyExternalDampingMatrixUserOrigin
				else if (origin == "Centre of mass")
					TranslateMatrix6(GetMsh().dt.Dlin, Eigen::Vector3d(GetMsh().dt.cg), Eigen::Vector3d(GetMsh().dt.c0));	
				else if (origin == "Body origin")
					;
				else
					throw Exc(in.Str() + "\n" + t_("Incorrect BodyExternalDampingMatrixOriginType"));
			} else if (fy.FirstIs("BodyExternalDampingMatrixUserOrigin")) {
				UVector<double> line = fy.GetVectorDouble();
				if (line.size() != 3)
					throw Exc(in.Str() + "\n" + t_("Incorrect BodyExternalDampingMatrixUserOrigin"));
				Vector3d ref(line[0], line[1], line[2]);
				TranslateMatrix6(GetMsh().dt.Dlin, ref, Eigen::Vector3d(GetMsh().dt.c0));
			} else if (fy.FirstIs("BodyExternalStiffnessMatrixOriginType")) {	
				String origin = fy.GetVal();
				if (origin == "User specified")
					;	// Wait for BodyExternalStiffnessMatrixUserOrigin
				else if (origin == "Centre of mass")
					TranslateMatrix6(GetMsh().dt.Cadd, Eigen::Vector3d(GetMsh().dt.cg), Eigen::Vector3d(GetMsh().dt.c0));	
				else if (origin == "Body origin")
					;
				else
					throw Exc(in.Str() + "\n" + t_("Incorrect BodyExternalDampingMatrixOriginType"));
			} else if (fy.FirstIs("BodyExternalStiffnessMatrixUserOrigin")) {
				UVector<double> line = fy.GetVectorDouble();
				if (line.size() != 3)
					throw Exc(in.Str() + "\n" + t_("Incorrect BodyExternalStiffnessMatrixUserOrigin"));
				Vector3d ref(line[0], line[1], line[2]);
				TranslateMatrix6(GetMsh().dt.Cadd, ref, Eigen::Vector3d(GetMsh().dt.c0));
			}
		} else if (fy.FirstMatch("FieldPointX*")) {
			UVector<UVector<double>> dat = fy.GetMatrixDouble(false);
			listPointsTemp.SetCount(dat.size());
			for (int r = 0; r < dat.size(); ++r) {
				if (dat[r].size() != 3)
					throw Exc(in.Str() + "\n" + t_("FieldPoint have to contain x, y, z"));
				for (int c = 0; c < 3; ++c) 
					listPointsTemp[r][c] = dat[r][c];
			}
		}
	}

	if (dataFrom == 'h') {
		for (double &w : dt.w)
			w *= 2*M_PI;	
	} else if (dataFrom == 's') {
		for (double &w : dt.w)
			w = 2*M_PI/w;	
	}

	if (dt.Nb == 0)
		throw Exc(t_("Incorrect .yml format"));
	
	for (ib = 0; ib < dt.css.size(); ++ib)
		dt.css[ib].dt.c0 = clone(dt.msh[ib].dt.c0);
	
	return true;
}

bool OrcaWave::Load_OF_YML() {
	String fileName = ForceExtSafer(dt.file, ".yml");
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		throw Exc(in.Str() + "\n" + t_("File not found or blocked"));
	
	dt.solver = Hydro::ORCAFLEX_YML;
	
	char dataFrom = 'r';
	dt.dimen = true;
	
	dt.x_w = dt.y_w = 0;
	
	dt.Nf = dt.Nh = Null;
	dt.Nb = 0; 
	
	dt.g = 9.80665;		// Default value used when SI units
	
	int ib = -1;

	UArray<Point3D> c0s;	

	auto Origin = [&] (int ib, const UVector<String> &norig) {
		UVector<String> snorig(norig.size());
		for (int i = 0; i < norig.size(); ++i)
			snorig[i] = Trim(norig[i]);
		if (IsNull(c0s[ib])) 
			c0s[ib].Set(ScanDouble(snorig[0]), ScanDouble(snorig[1]), ScanDouble(snorig[2]));
		else if (snorig[0] == "~" || snorig[1] == "~" || snorig[2] == "~")
			;
		else if (c0s[ib].x != ScanDouble(snorig[0]) || c0s[ib].y != ScanDouble(snorig[1]) || c0s[ib].z != ScanDouble(snorig[2]))
			throw Exc(in.Str() + "\n"  + F(t_("RAOOrigin for body %d (%s, %s, %s) diferent than the previously set (%f, %f, %f)"), 
							ib, snorig[0], snorig[1], snorig[2], c0s[ib].x, c0s[ib].y, c0s[ib].z));
	};
	auto Phase = [&] (const UVector<String> &norig) {
		UVector<String> snorig(norig.size());
		for (int i = 0; i < norig.size(); ++i)
			snorig[i] = Trim(norig[i]);
		if (snorig[0] == "~" || snorig[1] == "~" || snorig[2] == "~")
			;
		else if (!IsEqualRange<UVector<String>>({"0","0","0"}, snorig))
			throw Exc(in.Str() + "\n"  + F(t_("Only PhaseOrigin:[0,0,0] is supported. Read (%s, %s, %s)"), snorig[0], snorig[1], snorig[2]));
	};
		
	OrcaFactors factor;	

	YmlParser fy(in);

	FileInLine::Pos fpos = in.GetPos();

	while(fy.GetLine()) {
		if (fy.FirstIs("LoadRAOCalculationMethod")) 
			return false;
		else if (fy.FirstIs("General")) {
			if (fy.FirstIs("UnitsSystem")) {
				if (fy.GetVal() == "SI") {
					factor.mass = factor.force = 1000;
					factor.len = 1;	
					factor.Update();
				} else if (fy.GetVal() == "User") 
					;
				else
					throw Exc(in.Str() + "\n" + F(t_("Only SI and User units are supported. Read '%s'"), fy.GetVal()));
			} else if (fy.FirstIs("LengthUnits")) {
				factor.len = FactorLen(fy.GetVal());
				factor.Update();
			} else if (fy.FirstIs("MassUnits")) {
				factor.mass = FactorMass(fy.GetVal());
				factor.Update();
			} else if (fy.FirstIs("ForceUnits")) {
				factor.force = FactorForce(fy.GetVal());
				factor.Update();
			} else if (fy.FirstIs("g")) 
				dt.g = ScanDouble(fy.GetVal())*factor.len;
		} else if (fy.FirstIs("Environment")) {
			if (fy.FirstIs("WaterSurfaceZ") && fy.GetVal() != "0") 
				throw Exc(in.Str() + "\n" + F(t_("Only WaterSurfaceZ 0 is supported. Read '%s'"), fy.GetVal()));
			else if (fy.FirstIs("SeabedOriginDepth"))
				dt.h = ScanDouble(fy.GetVal())*factor.len;
		} else if (fy.FirstIs("VesselTypes")) {
			if (fy.FirstIs("Name")) {
				if (fy.Index() != dt.Nb)
					throw Exc(in.Str() + "\n" + t_("Failed body count"));
				dt.msh.SetCount(dt.Nb+1);
				dt.msh[dt.Nb].dt.name = fy.GetVal();
				dt.msh[dt.Nb].dt.SetCode(Body::ORCAFLEX_YML);
				dt.Nb++;
				c0s << Null;
			} else if (fy.FirstIs("WavesReferredToBy") && fy.Index() == 0) {		// Only for the first body
				String val = fy.GetVal();
				if (val.Find("(rad/s)") >= 0)
					dataFrom = 'r';
				else if (val.Find("(s)") >= 0)
					dataFrom = 's';
				else if (val.Find("(Hz)") >= 0)
					dataFrom = 'h';
				else
					throw Exc(in.Str() + "\n"  + F(t_("Unknown data in WavesReferredToBy: %s"), val));
			
			} else if (fy.FirstIs("SurgePositive")) {	
				if (fy.GetVal() != "forward")
					throw Exc(in.Str() + "\n"  + F(t_("Only SurgePositive: 'forward' is supported. Read '%s'"), fy.GetVal()));
			} else if (fy.FirstIs("SwayPositive")) {	
				if (fy.GetVal() != "port")
					throw Exc(in.Str() + "\n"  + F(t_("Only SwayPositive: 'port' is supported. Read '%s'"), fy.GetVal()));
			} else if (fy.FirstIs("HeavePositive")) {	
				if (fy.GetVal() != "up")
					throw Exc(in.Str() + "\n"  + F(t_("Only HeavePositive: 'up' is supported. Read '%s'"), fy.GetVal()));
			} else if (fy.FirstIs("RollPositiveStarboard")) {	
				if (fy.GetVal() != "down")
					throw Exc(in.Str() + "\n"  + F(t_("Only RollPositiveStarboard: 'down' is supported. Read '%s'"), fy.GetVal()));
			} else if (fy.FirstIs("PitchPositiveBowe")) {	
				if (fy.GetVal() != "down")
					throw Exc(in.Str() + "\n"  + F(t_("Only PitchPositiveBowe: 'down' is supported. Read '%s'"), fy.GetVal()));
			} else if (fy.FirstIs("YawPositiveBow")) {	
				if (fy.GetVal() != "port")
					throw Exc(in.Str() + "\n"  + F(t_("Only YawPositiveBow: 'port' is supported. Read '%s'"), fy.GetVal()));
			} else if (fy.FirstIs("QTFConventionsRotationOrder")) {	
				if (fy.GetVal() != "RzRyRx")
					throw Exc(in.Str() + "\n"  + F(t_("Only QTFConventionsRotationOrder: 'RzRyRx is supported. Read '%s'"), fy.GetVal()));
			
			} else if (fy.FirstIs("Draughts")) {
				ib = fy.GetIndex()[1];
				if (fy.FirstIs("DisplacementRAOs")) {
					if (fy.FirstIs("RAOOrigin")) 
						Origin(ib, fy.GetVector());
					else if (fy.FirstIs("PhaseOrigin")) 
						Phase(fy.GetVector());
					else if (fy.FirstIs("RAOs")) {
						if (fy.FirstIs("RAODirection") && fy.GetIndex()[1] == 0) {		// Only for the first body
							FindAdd(dt.head, ScanDouble(fy.GetVal()));
						}
					}
				} else if (fy.FirstIs("LoadRAOs")) {
					if (fy.FirstIs("RAOOrigin")) 
						Origin(ib, fy.GetVector());	
					else if (fy.FirstIs("PhaseOrigin")) 
						Phase(fy.GetVector());
					else if (fy.FirstIs("RAOs")) {
						if (fy.FirstIs("RAODirection") && fy.GetIndex()[1] == 0) {		// Only for the first body
							FindAdd(dt.head, ScanDouble(fy.GetVal()));
						}
					}
				} else if (fy.FirstIs("WaveDrift")) {
					if (fy.FirstIs("RAOOrigin")) {
						;//Origin(ib, fy.GetVector());	// If no QTF is set to 0,0,0
					} else if (fy.FirstIs("RAOs")) {
						if (fy.FirstIs("RAODirection") && fy.GetIndex()[1] == 0) {		// Only for the first body
							FindAdd(dt.head, ScanDouble(fy.GetVal()));
						}
					}
				} else if (fy.FirstIs("SumFrequencyQTFs")) {
					if (fy.FirstIs("RAOOrigin")) 
						;//Origin(ib, fy.GetVector());	
					else if (fy.FirstIs("PhaseOrigin")) 
						Phase(fy.GetVector());
				} else if (fy.FirstIs("OtherDampingOrigin")) {
					Origin(ib, fy.GetVector());	
					//if (!IsEqualRange<UVector<double>>({0,0,0}, fy.GetVectorDouble()))
					//	throw Exc(in.Str() + "\n"  + F(t_("Only OtherDampingOrigin:[0,0,0] is supported. Read '%s'"), fy.StrVar()));
				} else if (fy.FirstIs("ReferenceOrigin")) {
					Origin(ib, fy.GetVector());	
					//if (!IsEqualRange<UVector<double>>({0,0,0}, fy.GetVectorDouble()))
					//	throw Exc(in.Str() + "\n"  + F(t_("Only ReferenceOrigin:[0,0,0] is supported. Read '%s'"), fy.StrVar()));
				//} else if (fy.FirstIs("ReferenceOriginDatumPosition")) {
				//	if (!IsEqualRange<UVector<double>>({0,0,0}, fy.GetVectorDouble()))
				//		throw Exc(in.Str() + "\n"  + F(t_("Only ReferenceOriginDatumPosition:[0,0,0] is supported. Read '%s'"), fy.StrVar()));
				} else if (fy.FirstIs("FrequencyDependentAddedMassAndDamping")) {
					if (fy.FirstIs("AMDPeriodOrFrequency")) {
						if (((dataFrom == 'r' || dataFrom == 'h') && fy.GetVal() != "Infinity") ||
							 (dataFrom == 's' && fy.GetVal() != "0")) {
							if (dt.w.size() != fy.Index() - 1)		// -1 because Infinity is the first
								throw Exc(in.Str() + "\n" + t_("Failed frequencies count"));			
							dt.w << ScanDouble(fy.GetVal());
						}
					}
				}
			}
		} else if (fy.FirstIs("MultibodyGroups")) {
			if (fy.FirstIs("Bodies")) {
				if (fy.FirstIs("Name")) 
					ib = fy.Index();
			} else if (fy.FirstIs("MultibodyAddedMassAndDamping")) {
				if (fy.FirstIs("AMDPeriodOrFrequency")) {
					if (fy.GetVal() != "Infinity") {
						if (dt.w.size() != fy.Index() - 1)
							throw Exc(in.Str() + "\n" + t_("Failed frequencies count"));			
						dt.w << ScanDouble(fy.GetVal());
					}
				}
			}
		}
	}
	
	factor.Update();
	
	if (dt.Nb == 0)
		throw Exc(F("\n") + t_("No body found"));
	
	dt.msh.SetCount(dt.Nb);
	
	for (int iib = 0; iib < dt.Nb; ++iib)				
		for (int idf = 0; idf < 3; ++idf)
			dt.msh[iib].dt.c0[idf] = c0s[iib][idf];
	
	for (int iib = 0; iib < dt.Nb; ++iib) 
		dt.msh[iib].dt.C.setConstant(6, 6, 0);
	for (int iib = 0; iib < dt.Nb; ++iib) 
		dt.msh[iib].dt.M.setConstant(6, 6, 0);

	dt.Nf = dt.w.size();
	dt.Nh = dt.head.size();
	
	dt.Ainf.setConstant(dt.Nb*6, dt.Nb*6, NaNDouble);
	
	Initialize_AB(dt.A);
	Initialize_AB(dt.B);
		
	if (dataFrom == 'h') {
		for (int i = 0; i < dt.w.size(); ++i)
			dt.w[i] *= 2*M_PI;	
	} else if (dataFrom == 's'){
		for (int i = 0; i < dt.w.size(); ++i)
			dt.w[i] = 2*M_PI/dt.w[i];	
	}	

	Initialize_Forces(dt.ex);
	Initialize_Forces(dt.rao);

	int rrow = -1, ccol = -1;
	int idf = -1;
	bool diffFullQTF = false;
	
	in.SeekPos(fpos);
	
	while(fy.GetLine()) {
		if (fy.FirstIs("Environment")) {
			if (fy.FirstIs("Density")) 
				dt.rho = ScanDouble(fy.GetVal())*factor.rho;		// In kg/m3
			else if (fy.FirstIs("WaterDepth")) { 
				String sh = fy.GetVal(); 
				if (ToLower(sh) == "infinite")
					dt.h = -1;
				else	
					dt.h = ScanDouble(sh)*factor.len;
			}
		} else if (fy.FirstIs("VesselTypes")) {
			if (fy.FirstIs("Name")) 
				ib = fy.GetIndex()[1];
			else if (fy.FirstIs("WireFrameSymmetry")) {
				String sym = fy.GetVal();
				dt.symY = sym.Find("xz") >= 0;
				dt.symX = sym.Find("yz") >= 0;
				//foundWireFrameSymmetry = true;
			} else if (fy.FirstMatch("Vertex*")) {
				UVector<UVector<double>> mat = fy.GetMatrixDouble();
				
				dt.msh[ib].dt.mesh.nodes.SetCount(mat.size());
				for (int inn = 0; inn < mat.size(); ++inn) 
					dt.msh[ib].dt.mesh.nodes[inn] = Point3D(mat[inn]);
			} else if (fy.FirstMatch("PanelVertexIndex*")) {
				UVector<UVector<double>> mat = fy.GetMatrixDouble();
				
				int numNodes = dt.msh[ib].dt.mesh.nodes.size();
				dt.msh[ib].dt.mesh.panels.SetCount(mat.size());
				for (int ip = 0; ip < mat.size(); ++ip) {
					const UVector<double> &pan = mat[ip];
					if (Min(pan) < 1 || Max(pan) > numNodes)
						throw Exc(in.Str() + "\n" + F(t_("Wrong panel %d (%d,%d,%d,%d)"), ip+1, pan[0], pan[1], pan[2], pan[3]));
					for (int i4 = 0; i4 < 4; ++i4)
						dt.msh[ib].dt.mesh.panels[ip].id[i4] = int(pan[i4])-1;
				}		
			} else if (fy.FirstIs("Draughts")) {
				Eigen::MatrixXd &inertia = dt.msh[ib].dt.M;
				if (fy.FirstIs("Mass")) 
					inertia(0, 0) = inertia(1, 1) = inertia(2, 2) = ScanDouble(fy.GetVal())*factor.mass;
				else if (fy.FirstMatch("MomentOfInertiaTensor*")) {			// Referred to cg
					UVector<UVector<double>> mat = fy.GetMatrixDouble();
					
					inertia(3, 3) = mat[0][0]*factor.M(3, 3);		// In kg
					inertia(3, 4) = mat[0][1]*factor.M(3, 4);
					inertia(3, 5) = mat[0][2]*factor.M(3, 5);
					inertia(4, 3) = mat[1][0]*factor.M(4, 3);
					inertia(4, 4) = mat[1][1]*factor.M(4, 4);
					inertia(4, 5) = mat[1][2]*factor.M(4, 5);
					inertia(5, 3) = mat[2][0]*factor.M(5, 3);
					inertia(5, 4) = mat[2][1]*factor.M(5, 4);
					inertia(5, 5) = mat[2][2]*factor.M(5, 5);
				} else if (fy.FirstMatch("HydrostaticStiffnessz*")) {
					UVector<UVector<double>> mat = fy.GetMatrixDouble();
					
					for (int r = 0; r < 3; ++r)				// Only heave, roll, pitch
						for (int c = 0; c < 3; ++c)
							dt.msh[ib].dt.C(r+2, c+2) = mat[r][c]*factor.K(r+2, c+2);
				} else if (fy.FirstIs("CentreOfMass")) {
					UVector<double> line = fy.GetVectorDouble();
					
					dt.msh[ib].dt.cg.x = line[0]*factor.len;
					dt.msh[ib].dt.cg.y = line[1]*factor.len;
					dt.msh[ib].dt.cg.z = line[2]*factor.len;
				} else if (fy.FirstIs("CentreOfBuoyancy")) {
					UVector<double> line = fy.GetVectorDouble();
					
					dt.msh[ib].dt.cb.x = line[0]*factor.len;
					dt.msh[ib].dt.cb.y = line[1]*factor.len;
					dt.msh[ib].dt.cb.z = line[2]*factor.len;
				} else if (fy.FirstIs("DisplacedVolume")) 
					dt.msh[ib].dt.Vo = ScanDouble(fy.GetVal());
				else if (fy.FirstIs("DisplacementRAOs")) {
					if (fy.FirstIs("RAOs")) {
						if (fy.FirstMatch("RAOPeriodOrFrequency*")) {	
							int idh = fy.Index();
							if (idh < 0 || idh >= dt.head.size())
								throw Exc(in.Str() + "\n" + t_("Wrong heading"));
								
							UVector<UVector<double>> mat = fy.GetMatrixDouble();
							
							if (mat.size() != dt.Nf || mat[0].size() != 13)
								throw Exc(in.Str() + "\n"  + F(t_("Wrong number of numbers in DisplacementRAOs matrix"), fy.GetText()));

							for (int ifr = 0; ifr < dt.Nf; ++ifr) 
								for (int idof = 0; idof < 6; ++idof) 
									dt.rao[ib][idh](ifr, idof) = std::polar<double>(mat[ifr][1 + 2*idof]*factor.RAO(idof), ToRad(mat[ifr][1 + 2*idof + 1]));
						}
					}
				}  else if (fy.FirstIs("LoadRAOs")) {
					if (fy.FirstIs("RAOs")) {
						if (fy.FirstMatch("RAOPeriodOrFrequency*")) {	
							int idh = fy.Index();
							if (idh < 0 || idh >= dt.head.size())
								throw Exc(in.Str() + "\n" + t_("Wrong heading"));
								
							UVector<UVector<double>> mat = fy.GetMatrixDouble();
							
							if (mat.size() != dt.Nf || mat[0].size() != 13)
								throw Exc(in.Str() + "\n"  + F(t_("Wrong number of numbers in LoadRAOs matrix"), fy.GetText()));

							for (int ifr = 0; ifr < dt.Nf; ++ifr) 
								for (int idof = 0; idof < 6; ++idof) 
									dt.ex[ib][idh](ifr, idof) = std::polar<double>(mat[ifr][1 + 2*idof]*factor.F(idof), ToRad(mat[ifr][1 + 2*idof + 1]));
						}
					}
				} else if (fy.FirstIs("WaveDriftQTFMethod"))
					diffFullQTF = fy.GetVal() == "Full QTFs";
				else if (!diffFullQTF && fy.FirstIs("WaveDrift")) {
					if (fy.FirstIs("RAOs")) {
						if (fy.FirstMatch("RAOPeriodOrFrequency*")) {	
							if (!IsLoadedMD()) {
								dt.mdhead.resize(dt.head.size());
								for (int ih = 0; ih < dt.head.size(); ++ih)
									dt.mdhead[ih] = std::complex<double>(dt.head[ih], dt.head[ih]);
								Hydro::Initialize_MD(dt.md, dt.Nb, int(dt.mdhead.size()), dt.Nf);
							}
							
							int idh = fy.Index();
							if (idh < 0 || idh >= dt.head.size())
								throw Exc(in.Str() + "\n" + t_("Wrong heading"));
								
							UVector<UVector<double>> mat = fy.GetMatrixDouble();
							
							if (mat.size() != dt.Nf || mat[0].size() != 7)
								throw Exc(in.Str() + "\n"  + F(t_("Wrong number of numbers in WaveDrift matrix"), fy.GetText()));

							for (int ifr = 0; ifr < dt.Nf; ++ifr) 
								for (int idof = 0; idof < 6; ++idof) 
									dt.md[ib][idh][idof](ifr) = mat[ifr][1 + idof]*factor.MD(idof);
						}
					}
				} else if (fy.FirstIs("SumFrequencyQTFs") || (diffFullQTF && fy.FirstIs("WaveDrift"))) {
					if (fy.FirstMatch("RAOPeriodOrFrequency*")) {
						UVector<UVector<double>> mat = fy.GetMatrixDouble();
						
						UArray<UArray<UArray<MatrixXcd>>> &q = diffFullQTF ? dt.qtfdif : dt.qtfsum;
						double phmult = !diffFullQTF ? 1 : -1;		// Difference is conjugate-symmetric
						
						if (!IsLoadedQTF(!diffFullQTF)) {		// Gets frequencies and headings
							::Copy(dt.w, dt.qw);
							
							UArray<std::complex<double>> qh;
							for (int row = 0; row < mat.size(); ++row) {
								if (mat[row].size() != 16)
									throw Exc(in.Str() + "\n"  + t_("Wrong data in list"));
								double h1 = mat[row][2],
									   h2 = mat[row][3];
									   
								FindAddDelta(qh, std::complex<double>(h1, h2), 0.0001);	
							}
							::Copy(qh, dt.qhead);
	
	
							Hydro::Initialize_QTF(q, dt.Nb, int(qh.size()), dt.Nf);
							dt.mdtype = dt.qtftype = 9;
						}
						diffFullQTF = false;
						
						for (int row = 0; row < mat.size(); ++row) {
							if (mat[row].size() != 16)
								throw Exc(in.Str() + "\n"  + t_("Wrong data in list"));
							double w1 = mat[row][0],
								   w2 = mat[row][1],
								   h1 = mat[row][2],
								   h2 = mat[row][3];
								   
							if (dataFrom == 'h') {
								w1 *= 2*M_PI;
								w2 *= 2*M_PI;
							} else if (dataFrom == 's') {
								w1 = 2*M_PI/w1;
								w2 = 2*M_PI/w2;
							}
							int ifr1 = FindDelta(dt.qw, w1, 0.0001),
								ifr2 = FindDelta(dt.qw, w2, 0.0001),
								ih = FindDelta(dt.qhead, std::complex<double>(h1, h2), 0.0001);	
							if (ifr1 < 0)
								throw Exc(in.Str() + "\n"  + F(t_("Wrong frequency '%d' in QTF"), w1));
							if (ifr2 < 0)
								throw Exc(in.Str() + "\n"  + F(t_("Wrong frequency '%d' in QTF"), w2));
							if (ih < 0)
								throw Exc(in.Str() + "\n"  + F(t_("Wrong head (%d,%d) in QTF"), h1, h2));
							for (int iidf = 0; iidf < 6; ++iidf) {
								double mag = mat[row][4+iidf*2]*factor.F(iidf);
								double ph  = ToRad(mat[row][4+iidf*2+1]);
								q[ib][ih][iidf](ifr1, ifr2) = std::polar<double>(mag, ph);
								q[ib][ih][iidf](ifr2, ifr1) = std::polar<double>(mag, ph*phmult);	// Diagonal
							}
						}
					}
				} else if (fy.FirstIs("FrequencyDependentAddedMassAndDamping")) {
					if (fy.FirstMatch("AddedMassMatrixX*")) {
						idf = fy.Index()-1;
						if (idf < -1 || idf >= dt.w.size())		// Infinity is the first
							throw Exc(in.Str() + "\n" + t_("Wrong frequency"));			
						
						UVector<UVector<double>> mat = fy.GetMatrixDouble();
						
						if (idf == -1) {
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) 
									dt.Ainf(r, c) = mat[r][c]*factor.A(r, c); 
							}
						} else {
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) 
									dt.A[r][c](idf) = mat[r][c]*factor.A(r, c); 
							}
						}
					} else if (fy.FirstMatch("DampingX*")) {
						idf = fy.Index()-1;
						if (idf < -1 || idf >= dt.w.size())		// Infinity is the first
							throw Exc(in.Str() + "\n" + t_("Wrong frequency"));			
						
						if (idf >= 0) {
							UVector<UVector<double>> mat = fy.GetMatrixDouble();
							
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) 
									dt.B[r][c](idf) = mat[r][c]*factor.B(r, c); 
							}
						}
					}
				}
			}
		} else if (fy.FirstIs("MultibodyGroups")) {
			if (fy.FirstIs("Bodies")) {
				if (fy.FirstIs("Name")) {
					ib = fy.Index();
					dt.msh[ib].dt.name = fy.GetVal();
				} else if (fy.FirstIs("DisplacedVolume")) 
					dt.msh[ib].dt.Vo = ScanDouble(fy.GetVal())*factor.len*factor.len*factor.len;
				else if (fy.FirstIs("CentreOfBuoyancy")) {
					UVector<double> line = fy.GetVectorDouble();
					
					dt.msh[ib].dt.cb.x = line[0]*factor.len;
					dt.msh[ib].dt.cb.y = line[1]*factor.len;
					dt.msh[ib].dt.cb.z = line[2]*factor.len;
				} else if (fy.FirstMatch("HydrostaticStiffnessz*")) {
					UVector<UVector<double>> mat = fy.GetMatrixDouble();
					
					for (int r = 0; r < 3; ++r)				// Only heave, roll, pitch
						for (int c = 0; c < 3; ++c)
							dt.msh[ib].dt.C(r+2, c+2) = mat[r][c]*factor.K(r, c);
				}
			} else if (fy.FirstIs("MultibodyAddedMassAndDamping")) {
				if (fy.FirstIs("AMDPeriodOrFrequency")) {
					idf = fy.Index()-1;
					if (idf < -1 || idf >= dt.w.size())		// Infinity is the first
						throw Exc(in.Str() + "\n" + t_("Wrong frequency"));	
						
				} else if (fy.FirstIs("Matrices")) {
					if (fy.FirstIs("Row")) 
						rrow = ScanInt(fy.GetVal())-1;
					else if (fy.FirstIs("Column")) 
						ccol = ScanInt(fy.GetVal())-1;
					else if (fy.FirstMatch("AddedMassX*")) {
						UVector<UVector<double>> mat = fy.GetMatrixDouble();
						
						if (idf == -1) {
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) {
									dt.Ainf(r+rrow*6, c+ccol*6) = mat[r][c]*factor.A(r, c); 
									if (rrow != ccol)		// Fill the symmetric
										dt.Ainf(r+ccol*6, c+rrow*6) = dt.Ainf(r+rrow*6, c+ccol*6); 
								}
							}
						} else {
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) {
									dt.A[r+rrow*6][c+ccol*6](idf) = mat[r][c]*factor.A(r, c); 
									if (rrow != ccol)		// Fill the symmetric
										dt.A[r+ccol*6][c+rrow*6](idf) = dt.A[r+rrow*6][c+ccol*6](idf); 
								}
							}
						}
					} else if (fy.FirstMatch("DampingX*")) {
						if (idf >= 0) {
							UVector<UVector<double>> mat = fy.GetMatrixDouble();
							
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) {
									dt.B[r+rrow*6][c+ccol*6](idf) = mat[r][c]*factor.B(r, c); 
									if (rrow != ccol)		// Fill the symmetric
										dt.B[r+ccol*6][c+rrow*6](idf) = dt.B[r+rrow*6][c+ccol*6](idf);
								}
							}
						}
					}
				}
			}
		}
	}
	
	if (dt.Nb == 0)
		throw Exc(t_("Incorrect .yml format"));
	
	// Inertia matrices have to be translated from cg to c0
	for (int iib = 0; iib < dt.Nb; ++iib)
		Surface::TranslateInertia66(dt.msh[iib].dt.M, dt.msh[iib].dt.cg, dt.msh[iib].dt.cg, dt.msh[iib].dt.c0);
	
	return true;
}

void OrcaWave::SaveCase_OW_YML(String folder, bool bin, int numThreads, bool withPotentials, bool withMesh, bool x0z, bool y0z, UVector<Point3D> &listPoints,
								bool irregular, bool autoIrregular, int qtfType, bool autoQTF) const {
	bool createQTFFreeSurface = true;
	bool onlyMeanDrift = qtfType > 10;
	
	qtfType %= 10;
	
	if (irregular && qtfType == 8)
		throw Exc(t_("The irregular frequencies removal cannot be set with the far field/momentum conservation method"));

	if (!onlyMeanDrift && qtfType == 8)
		throw Exc(t_("The Momentum Conservation method cannnot perform Full QTF calculation"));
				
	bool isv15 = false;
	
	String exeName = "bemrosetta_cl";
	if (bin) {
		String exe = GetExeFilePath();		
		String ext = GetFileExt(exe);
		String exe_cl = AFX(GetFileFolder(GetExeFilePath()), exeName + ext);
		if (FileExists(exe_cl))
			FileCopy(exe_cl, AFX(folder, exeName + ext));
		else {
			exeName = "bemrosetta";
			FileCopy(exe, AFX(folder, exeName + ext));
		}
	}
	if (IsNull(numThreads) || numThreads <= 0)
		numThreads = 8;
		
	String name = GetFileTitle(folder); 
	String fileYaml = AFX(folder, "OrcaWave.wave.yml");
	String fileBat  = AFX(folder, "OrcaWave_bat.bat");
	
	DirectoryCreate(folder);
	
	FileOut bat;
	if (!bat.Open(fileBat))
		throw Exc(F(t_("Impossible to open file '%s'"), fileBat));
	bat << BatchStart();
	const Point3D &c0 = dt.msh[0].dt.c0;
	bat << F("%s -orca -numtries 10 -numthread %d -rw \"%s\" \"%s\"", exeName, numThreads, "OrcaWave.wave.yml", "OrcaWave.flex.yml");
	bat << BatchEnd();
	
	FileOut	out;
	if (!out.Open(fileYaml))
		throw Exc(F(t_("Impossible to open file '%s'"), fileYaml));
	
	for (int ib = 0; ib < dt.Nb; ++ib) {
		const Body &b = dt.msh[ib];
		Body::SaveAs(b, AFX(folder, F("Body_%d.gdf", ib+1)), Body::WAMIT_GDF, Body::UNDERWATER, rho_ndim(), g_ndim(), y0z, x0z);
		if (irregular && !autoIrregular) {
			const Body &l = dt.lids[ib];
			Body::SaveAs(l, AFX(folder, F("Body_%d_lid.gdf", ib+1)), Body::WAMIT_GDF, Body::ALL, rho_ndim(), g_ndim(), y0z, x0z);			
		}
		if (qtfType == 7 && !autoQTF) {
			const Body &c = dt.css[ib];
			Body::SaveAs(c, AFX(folder, F("Body_%d_cs.gdf", ib+1)), Body::WAMIT_GDF, Body::ALL, rho_ndim(), g_ndim(), y0z, x0z);			
		}
	}
	
	String machineName, domain, ip4, ip6;
	GetNetworkInfo(machineName, domain, ip4, ip6);
	Time t = GetSysTime(); 
	
	out << 	"\%YAML 1.1\n"
			"# Type: Diffraction\n"
			"# Program: OrcaWave " << Orca::BEMRVersion() << "\n"
			"# File: " << fileYaml << "\n"
			"# Created: " << F("%02d:%02d", t.hour, t.minute) << " on " << F("%02d/%02d/%04d", t.day, t.month, t.year) << "\n"
			"# User: " << GetUserName() << "\n"
			"# Machine: " << machineName << "\n"
			"---\n";
			
	out << 	"# Model\n"
			"UnitsSystem: User\n"
			"LengthUnits: m\n"
			"MassUnits: kg\n"
			"ForceUnits: N\n"
			"g: " << F("%.5f", dt.g) << "\n";
	
	out << 	"# Calculation & output\n";
	out << 	"SolveType: ";
	if (!listPoints.IsEmpty())
		out << "Potential and source formulations";
	else if (qtfType > 0 ) {
		if (onlyMeanDrift)
			out << "Potential and source formulations";
		else
			out << "Full QTF calculation";	
	} else
		out << 	   "Potential formulation only";
	out <<  "\n";
	
	out << 	"LoadRAOCalculationMethod: Diffraction\n";
	if (qtfType == 7 || qtfType == 9)
		out << "QuadraticLoadPressureIntegration: " << (qtfType == 9 ? "Yes" : "No") << "\n";
	out <<	"QuadraticLoadControlSurface: " << (qtfType == 7 ? "Yes" : "No") << "\n"
			"QuadraticLoadMomentumConservation: " << (qtfType == 8 ? "Yes" : "No") << "\n";
	if (qtfType > 0) {		
		if (qtfType == 7)
			out << "PreferredQuadraticLoadCalculationMethod: Control surface\n";
		else if (qtfType == 8)
			out << "PreferredQuadraticLoadCalculationMethod: Momentum conservation\n"
				   "MomentumConservationNodeCount: 64\n";
		else if (qtfType == 9)
			out << "PreferredQuadraticLoadCalculationMethod: Pressure integration\n";
	}
	out << 	"HasResonanceDampingLid: No\n"
			"LengthTolerance: 100e-9\n"
			"WaterlineZTolerance: 1e-6\n"
			"WaterlineGapTolerance: 1e-6\n"
			"DivideNonPlanarPanels: No\n"
			"LinearSolverMethod: Direct LU\n"
			"OutputPanelPressures: " << (withPotentials ? "Yes" : "No") << "\n"
			//"OutputPanelVelocities: No\n"
			"OutputBodyWireFrames: " << (withMesh ? "Yes" : "No") << "\n"
			"OutputIntermediateResults: " << (withPotentials ? "Yes" : "No") << "\n"
			"ValidatePanelArrangement: No\n"
			"BodyVolumeWarningLevel: 1e-12\n"
			"PanelAspectRatioWarningLevel: 25\n"
			"PanelsPerWavelengthWarningLevel: 5\n"
			"ComputationStrategy: Optimised for memory\n";			// Also "Optimised for run time"
	
	out << 	"# Environment\n";
	//if (isv15)
		out << 	"WaterDepth: " << (dt.h > 0 ? F("%.2f", dt.h) : "Infinity") << "\n";
	//else {
	//	out << 	"SeabedOriginDepth: " << (dt.h > 0 ? F("%.2f", dt.h) : "Infinity") << "\n"
	//			"NominalDepth: ~\n"; 
	//}
	out <<  "WaterDensity: " << F("%.2f", dt.rho) << "\n"
			"WavesReferredToBy: frequency (rad/s)\n"
			"HasWaveSpectrumForDragLinearisation: No\n"
			"MorisonFluidVelocity: Undisturbed incident wave\n"
			"PeriodOrFrequency:\n";
	for (int ifr = 0; ifr < dt.Nf; ++ifr)
		out << "  - " << F("%.5f", dt.w[ifr]) << "\n";
	out << 	"WaveHeading:\n";
	for (int ih = 0; ih < dt.Nh; ++ih)
		out << "  - " << F("%.1f", dt.head[ih]) << "\n";
	
	if (qtfType > 0) 
		out << 	"QTFMinCrossingAngle: 0\n"
				"QTFMaxCrossingAngle: 180\n";
	if (!onlyMeanDrift && (qtfType == 7 || qtfType == 9)) 
		out <<	"QTFMinPeriodOrFrequency: 0\n"
				"QTFMaxPeriodOrFrequency: Infinity\n"
				"QTFFrequencyTypes: Both\n";
	
	out << 	"# Bodies\n"
			"Bodies:\n";
	for (int ib = 0; ib < dt.Nb; ++ib) {
		const Body &b = dt.msh[ib];
		const Body::Data &d = b.dt;
		
		double panelSize = d.under.avgFacetSideLen;
		double separationFromBody = panelSize*8;		// Orcina advises 5-10
		
		out <<	"  - BodyName: " << d.name << "\n"
				"    BodyMeshPosition: [0, 0, 0]\n"
				"    BodyMeshAttitude: [0, 0, 0]\n"
				"    BodyIncludedInAnalysis: Yes\n"
				"    BodyMeshFileName: " << F(".\\Body_%d.gdf", ib+1) << "\n"
				"    BodyMeshFormat: Wamit gdf\n"
				"    BodyMeshLengthUnits: m\n"
				"    BodyMeshSymmetry: ";
		if (y0z) {	
			if (x0z)
				out << 	"xz and yz planes";
			else
				out << 	"yz plane";
		} else if (x0z)
			out << 	"xz plane";
		else
			out << 	"None";
		out << 	"\n"
				"    BodyMeshDipolePanels: \n";
		if (!isv15) {
			if (!IsNull(d.c0))
				out <<	"    BodyOriginType: User specified\n"
						"    BodyUserOrigin: [" << F("%f, %f, %f", d.c0.x, d.c0.y, d.c0.z) << "]\n";
			else
				out <<	"    BodyOriginType: Free surface\n";
		}
		if (!irregular) 
			out <<	"    BodyAddInteriorSurfacePanels: No\n";
		else {
			if (autoIrregular)
				out <<	"    BodyAddInteriorSurfacePanels: Yes\n"
						"    BodyInteriorSurfacePanelMethod: Triangulation method\n";
			else
				throw Exc(t_("OrcaWave does not support irregular frequencies removal through user supplied lid"));
		}
		if (qtfType == 7) {
			if (!autoQTF) {
				out <<	"    BodyControlSurfaceType: Defined by mesh file\n"
						"    BodyControlSurfaceMeshFileName: " << F("Body_%d_cs.gdf", ib+1) << "\n"
						"    BodyControlSurfaceMeshFormat: Wamit gdf\n"
						"    BodyControlSurfaceMeshLengthUnits: m\n";				
			} else {
				out <<	"    BodyControlSurfaceType: Automatically generated\n"
						"    BodyControlSurfacePanelSize: " << panelSize << "\n"
						"    BodyControlSurfaceSeparationFromBody: " << separationFromBody << "\n"
						"    BodyControlSurfaceIncludeFreeSurface: Yes\n";	
			}
		}
    	out <<	"    BodyOrcaFlexImportSymmetry: Use global mesh symmetry\n"
    			"    BodyOrcaFlexImportLength: 1\n"
    			"    BodyHydrostaticIntegralMethod: Standard\n"
    			"    BodyHydrostaticStiffnessMethod: Displacement\n"
    			"    BodyInertiaSpecifiedBy: Matrix (for a general body)\n"
    			"    BodyCentreOfMass: [" << F("%f, %f, %f", d.cg.x, d.cg.y, d.cg.z) << "]\n"
    			"    BodyMass: " << b.GetMass() << "\n"
    			"    BodyInertiaTensorRx, BodyInertiaTensorRy, BodyInertiaTensorRz:\n";
		if (d.M.size() > 0)
			out << 	"      - [" << F("%f, %f, %f", Nvl(d.M(3, 3), 0.), Nvl(d.M(3, 4), 0.), Nvl(d.M(3, 5), 0.)) << "]\n"
					"      - [" << F("%f, %f, %f", Nvl(d.M(4, 3), 0.), Nvl(d.M(4, 4), 0.), Nvl(d.M(4, 5), 0.)) << "]\n"
					"      - [" << F("%f, %f, %f", Nvl(d.M(5, 3), 0.), Nvl(d.M(5, 4), 0.), Nvl(d.M(5, 5), 0.)) << "]\n";
		else
			out <<	"      - [0.001, 0, 	0]\n"
					"      - [0, 	 0.001,	0]\n"
					"      - [0, 	 0, 	0.001]\n";
    	out << 	"    BodyInertiaTensorOriginType: Body origin\n" //User specified\n"
    			//"    BodyInertiaTensorUserOrigin: [" << F("%f, %f, %f", 0, 0, 0) << "]\n" // d.c0.x, d.c0.y, d.c0.z) << "]\n"
    			"    BodyExternalStiffnessMatrixx, BodyExternalStiffnessMatrixy, BodyExternalStiffnessMatrixz, BodyExternalStiffnessMatrixRx, BodyExternalStiffnessMatrixRy, BodyExternalStiffnessMatrixRz:\n";
    	if (d.Cadd.size() == 36) {
    		for (int r = 0; r < 6; ++r) {
    			out << "      - [";
    			for (int c = 0; c < 6; ++c) {
    				if (c > 0)
    					out << ", ";
    				out << Nvl(d.Cadd(r, c), 0.);
    			}
    			out << "]\n";
    		}
		} else if (d.Cmoor.size() == 36) {
    		for (int r = 0; r < 6; ++r) {
    			out << "      - [";
    			for (int c = 0; c < 6; ++c) {
    				if (c > 0)
    					out << ", ";
    				out << Nvl(d.Cmoor(r, c), 0.);
    			}
    			out << "]\n";
    		}
		} else {
			for (int r = 0; r < 6; ++r) 
				out << "      - [0, 0, 0, 0, 0, 0]\n";
		}
    	out <<	"    BodyExternalStiffnessMatrixOriginType: Body origin\n"
    			"    BodyExternalDampingMatrixx, BodyExternalDampingMatrixy, BodyExternalDampingMatrixz, BodyExternalDampingMatrixRx, BodyExternalDampingMatrixRy, BodyExternalDampingMatrixRz:\n";
    	if (d.Dlin.size() > 0) {
    		for (int r = 0; r < 6; ++r) {
    			out << "      - [";
    			for (int c = 0; c < 6; ++c) {
    				if (c > 0)
    					out << ", ";
    				out << Nvl(d.Dlin(r, c), 0.);
    			}
    			out << "]\n";
    		}
		} else {
			for (int r = 0; r < 6; ++r) 
				out << "      - [0, 0, 0, 0, 0, 0]\n";
		}		
    	out <<	"    BodyExternalDampingMatrixOriginType: Body origin\n"
    			"    BodyConnectionParent: Free\n"
    			"    BodyIncreaseRollDampingToTarget: No\n"
    			"    BodyFixedDOFx: No\n"
    			"    BodyFixedDOFy: No\n"
    			"    BodyFixedDOFz: No\n"
    			"    BodyFixedDOFRx: No\n"
    			"    BodyFixedDOFRy: No\n"
    			"    BodyFixedDOFRz: No\n";
	}
	out <<	"# Field points\n";
	if (!listPoints.IsEmpty()) {
		out << "FieldPointX, FieldPointY, FieldPointZ:\n";
		for (int i = 0; i < listPoints.size(); ++i) 
			out << F("- [%8.6g, %8.6g, %8.6g]\n", listPoints[i].x, listPoints[i].y, listPoints[i].z);
	}
	out <<	"DetectAndSkipFieldPointsInsideBodies: Yes\n";
						
	if (!onlyMeanDrift && (qtfType == 7 || qtfType == 9)) {
		out << 	"# QTFs\n"
					"QTFCalculationMethod: Direct method\n";
		if(!createQTFFreeSurface) {
			out << 	"FreeSurfacePanelledZoneType: Defined by mesh file\n"
					"FreeSurfacePanelledZoneMeshFileName: \n"
					"FreeSurfacePanelledZoneMeshFormat: Wamit gdf\n"
					"FreeSurfacePanelledZoneMeshLengthUnits: m\n";
		} else {
			Surface first = clone(First(dt.msh).dt.mesh);
			first.GetSegments();
			double panelSize = first.GetAvgLenSegment();
			double rsb = 0; // Smallest bounding sphere radius of the bodies 
			for (int ib = 0; ib < dt.Nb; ++ib) {
				const VolumeEnvelope &env = dt.msh[ib].dt.mesh.env;
				rsb = max(rsb, env.Max());
			}
			double innerRadius = rsb + 10*panelSize;
			double Tmin = 4;
			double lambda_min = SeaWaves::WaveLength(max(Tmin, 2*M_PI/(2*Last(dt.w))), dt.h, Bem().g);
			double lambda_max = SeaWaves::WaveLength(2*M_PI/First(dt.w), dt.h, Bem().g);
			double roc = rsb;
			if (dt.h > 0)
				roc += min(dt.h, lambda_max);		// Outer Circle Radius
			else
				roc += lambda_max;
			//int numberOfAnnuli = max(1, int((roc - innerRadius)/2./lambda_min));
			//double radiusStep = (roc - innerRadius)/numberOfAnnuli;
			//int numberOfAzimuthalNodes = int(6*M_PI*innerRadius/lambda_min);
			int numberOfSegments = min(1000, int(6*M_PI*roc/lambda_min));

			out << 	"FreeSurfacePanelledZoneType: Automatically generated\n";
			out << 	F("FreeSurfacePanelledZonePanelSize: %.3f\n", panelSize*1.2);
			out << 	F("FreeSurfacePanelledZoneInnerRadius: %.3f\n", innerRadius);
			//out << 	F("FreeSurfaceQuadratureZoneNumberOfAnnuli: %d\n", numberOfAnnuli);
			//out << 	F("FreeSurfaceQuadratureZoneRadiusStep: %.3f\n", radiusStep);
			//out << 	"FreeSurfaceQuadratureZoneNumberOfRadialNodes: 6\n";
			//out << 	F("FreeSurfaceQuadratureZoneNumberOfAzimuthalNodes: %d\n", numberOfAzimuthalNodes);
			out << 	F("FreeSurfaceOuterCircleNumberOfSegments: %d\n", numberOfSegments);
			out << 	"FreeSurfaceAsymptoticZoneExpansionOrder: 24\n";
		}
	}
	out << 	"...\n";
}
