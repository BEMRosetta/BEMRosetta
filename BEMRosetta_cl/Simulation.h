// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors

#ifndef _BEMRosetta_cl_Simulation_h_
#define _BEMRosetta_cl_Simulation_h_

#include <STEM4U/SeaWaves.h>
#include "FastOut.h"

class Simulation {
public:
	~Simulation();
	
	void Load(const String &strFile, int stiffMod, int dllForce, double _rho, double _g, 
			double c0x, double c0y, double c0z, const BEM &bem);

	Force6D CalcStiff_Static(const float *pos);
	Force6D CalcStiff_DynamicStatic(double time, const float *pos, double volTolerance);
	Force6D CalcStiff_Dynamic(double time, const float *pos, double volTolerance, SeaWaves &waves);
	Force6D CalcStiff(double time, const float *pos, double volTolerance, SeaWaves &waves);
	Force6D CalcForces(double time, const float *pos, const float *vel, const float *acc);
	
	enum CALC_TYPE {NONE, STATIC, DYN_STATIC, DYNAMIC, OPENFAST, UNKNOWN};
	static const char *strCalculation[];
	CALC_TYPE calculation = UNKNOWN;

	MatrixXd stiff;
	Point3D cb;
	
	Mesh mesh;
	
	double rho, g;
	
private:	
	VectorXd forceb;
	
	Point3D dampingCentre;
	MatrixXd linearDamping;
	MatrixXd quadDamping;
	
	struct InterpolateVal {
		UVector<double> vx, vy;
		void Add(String str) {
			if (str.Find(':') < 0) 
				vy << ScanDouble(str);
			else {
				UVector<String> svals = Split(str, ':');
				for (int i = 0; i < svals.size(); ++i) {
					UVector<String> sxy = Split(svals[i], ',');
					if (sxy.size() != 2)
						throw Exc(Format("Format problem in '%s'", str));
					vx << ScanDouble(sxy[0]);
					vy << ScanDouble(sxy[1]);
				}
			}
		}
		double LinearInterpolate(double x) {
			if (vx.size() == 0)
				return vy[0];
			else
				return ::LinearInterpolate(x, vx, vy);
		}
	};
	struct ForceInterpolate {
		Point3D p;
		InterpolateVal tx, ty, tz, rx, ry, rz;
		void Set(double x, double y, double z, String stx, String sty, String stz, 
				String srx, String sry, String srz) {
			p.Set(x, y, z);	
			tx.Add(stx);
			ty.Add(sty);
			tz.Add(stz);
			rx.Add(srx);
			ry.Add(sry);
			rz.Add(srz);
		}
	};
	
	UArray<ForceInterpolate> coupledForces, fixedForces;
	
	FastOut out;
	struct Outputs {
		int time, 
			ptfmSurge, ptfmSway, ptfmHeave, ptfmRoll, ptfmPitch, ptfmYaw,
			wave1Elev,
			ptfmSurge_cd, ptfmSway_cd, ptfmHeave_cd, ptfmRoll_cd, ptfmPitch_cd, ptfmYaw_cd,
			ptfmTVxi_cd, ptfmTVyi_cd, ptfmTVzi_cd, ptfmRVxi_cd, ptfmRVyi_cd, ptfmRVzi_cd,
			ptfmTAxi_cd, ptfmTAyi_cd, ptfmTAzi_cd, ptfmRAxi_cd, ptfmRAyi_cd, ptfmRAzi_cd;
		int ptfmStiffFx, ptfmStiffFy, ptfmStiffFz, ptfmStiffMx, ptfmStiffMy, ptfmStiffMz;
		int ptfmCdFx, ptfmCdFy, ptfmCdFz, ptfmCdMx, ptfmCdMy, ptfmCdMz;
		int ptfmCd2Fx, ptfmCd2Fy, ptfmCd2Fz, ptfmCd2Mx, ptfmCd2My, ptfmCd2Mz;
		int ptfmCBx, ptfmCBy, ptfmCBz;
		int ptfmVol;
	} output;
	
	String folder;
};

#endif
