// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors

#ifndef _BEMRosetta_cl_Simulation_h_
#define _BEMRosetta_cl_Simulation_h_

class Simulation {
public:
	void Load(const String &strFile, double _rho, double _g, double c0x, double c0y, double c0z);

	Force6D CalcForces_Static(const float *pos);
	Force6D CalcForces_Dynamic(const float *pos, double volTolerance);
	Force6D CalcForces(const float *pos, const float *vel, const float *acc);

private:	
	MatrixXd stiff;
	VectorXd forceb;
	
	Point3D dampingCentre;
	MatrixXd linearDamping;
	MatrixXd quadDamping;
	MatrixXd addedMass;
	
	UArray<ForceVector> coupledForces, fixedForces;
	
	double rho, g;
	Mesh mesh;
};

#endif
