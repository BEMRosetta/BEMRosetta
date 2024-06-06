// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#ifndef _BEMRosetta_BEMRosetta_cl_functions_h_
#define _BEMRosetta_BEMRosetta_cl_functions_h_

#include <ScatterDraw/ScatterDraw.h>

String GetBEMRosettaDataFolder();

double GetKirfMaxT(const Eigen::VectorXd &w);
void GetTirf(Eigen::VectorXd &Tirf, int numT, double maxT);
	
void GetKirf(Eigen::VectorXd &Kirf, const Eigen::VectorXd &Tirf, const Eigen::VectorXd &w, const Eigen::VectorXd &B);
void GetKirf(Eigen::VectorXd &Kirf, const Eigen::VectorXd &Tirf, double w0, double dw, 	 const Eigen::VectorXd &B);	
double GetAinf_Kirf(Eigen::VectorXd &Kirf, const Eigen::VectorXd &w, const Eigen::VectorXd &A, const Eigen::VectorXd &B, int numT, double maxT);
double GetAinf_Kirf(Eigen::VectorXd &Kirf, double w0, double dw, const Eigen::VectorXd &A, const Eigen::VectorXd &B, int numT, double maxT);

void GetAinf_w(Eigen::VectorXd &Ainf_w, const Eigen::VectorXd &Kirf, const Eigen::VectorXd &Tirf, const Eigen::VectorXd &w, 
			const Eigen::VectorXd &A);
double GetAinf(const Eigen::VectorXd &Kirf, const Eigen::VectorXd &Tirf, const Eigen::VectorXd &w, 
			const Eigen::VectorXd &A);
double GetA0(const VectorXd &Kirf, const VectorXd &Tirf, double ainf);
void GetA(Eigen::VectorXd &A, const Eigen::VectorXd &Kirf, const Eigen::VectorXd &w, double ainf, double dt);

//double Fradiation2(double t, const Eigen::VectorXd &vel, const Eigen::VectorXd &irf, double dt);
double Fradiation(const Eigen::VectorXd &vel, const Eigen::VectorXd &irf, Eigen::Index iiter, double dt, Eigen::Index velSize = -1);

struct ParamDampedSin {
	double z0,
		   zDecay, 
		   ainf,			// Added mass at infinite frequency 
		   b, 
		   gamma, 
		   w_d, 			// Damped natural frequency
		   w02,				// Undamped natural frequency 
		   eta, 			// Damping ratio (fraction of the critical damping)
		   phi, 
		   t0;
	bool gett0, getz0;
	
	ParamDampedSin() {Clear();}
	
	void Clear() {
		z0 = 0;
		zDecay = ainf = gamma = w_d = 0.1;
		b = w02 = eta = phi = Upp::Null;
		gett0 = getz0 = false;
		t0 = 0;
	}
};

struct ParamDecay {
	double ainf, b, b2, av;
	Eigen::VectorXd Bdec, wdec;	// Real unknowns
	int numB;					// Num. of unknown B values
	double maxFreq;				// rad/s	max w for B
	Eigen::VectorXd Bspl, wspl;	// Processed ready to get Kirf
	Eigen::VectorXd B, wB;		// Real data, if available
	Eigen::VectorXd Kirf;
	SplineEquation splineB;
	
	bool getb, getb2, getav;
	
	ParamDecay() {Clear();}
	
	void Clear() {
		ainf = b = b2 = av = 0;
		numB = 0;			// Num. of unknown B values
		getb = getb2 = getav = false;
		maxFreq = 15;
		Bdec.resize(0);
		wdec.resize(0);
		Bspl.resize(0);
		wspl.resize(0);
		B.resize(0);
		wB.resize(0);
		Kirf.resize(0);
	}
};

double FitToDampedSin(const Eigen::VectorXd &y, double dt, double mass, double Kh, double g, ParamDampedSin &par);
double DampedSin(double x, double z0, double zDecay, double mass, double ainf, double b, double w_d, double t0, double phi);
void DampedSin(double x, double z0, double zDecay, double mass, double ainf, double b, double w_d, double t0, double phi,
		double &y, double &dy, double &d2y);
		
void FitToDecay(const Eigen::VectorXd &z, const Eigen::VectorXd &dz, const Eigen::VectorXd &d2z, 
			double dt, double mass, double Kh, double g, ParamDecay &param);
void Decay(double mass, double ainf, double av, double Kh, double b, double b2, double dt, double zDecay, double maxT, Eigen::VectorXd &z);

double FixHeading_0_360(double head);
double FixHeading_180(double head);
double FixHeading(double head, BEM::HeadingType range);
std::complex<double> FixHeading_180(const std::complex<double> &head);
std::complex<double> FixHeading_0_360(const std::complex<double> &head);
std::complex<double> FixHeading(const std::complex<double> &head, BEM::HeadingType range);

bool SortComplex(const std::complex<double>& a, const std::complex<double>& b);
	
void SetPhaseToMag(std::complex<double> &val, double arg);
void AddPhase(std::complex<double> &val, double arg);

double FactorLen(const char *units);
double FactorMass(const char *units);
double FactorForce(const char *units);
	
#endif
