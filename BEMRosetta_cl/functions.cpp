#include "BEMRosetta.h"
#include <STEM4U/Sundials.h>
#include <STEM4U/Integral.h>
#include <STEM4U/Utility.h>

#include "functions.h"

using namespace Eigen;


void GetKirfTirf(VectorXd &Kirf, VectorXd &Tirf, const VectorXd &w, const VectorXd &B, double dt, double &maxT) {
	ASSERT(B.size() >= 2);
	
    maxT = min(maxT, M_PI/(w[1] - w[0]));
    int numT = int(maxT/dt);
    
    Tirf = VectorXd::LinSpaced(numT, 0, maxT);
    Kirf.resize(numT);
	
	size_t Nf = B.size();
	
	VectorXd y(Nf);
    for (int it = 0; it < numT; ++it) {
		for (int iw = 0; iw < Nf; ++iw)
			y(iw) = B(iw)*cos(w(iw)*Tirf(it));
		Kirf(it) = Integral(w, y, SIMPSON_1_3)*2/M_PI;
	}
}

void GetKirfTirf(VectorXd &Kirf, VectorXd &Tirf, double w0, double dw, const VectorXd &B, double dt, double &maxT) {
	ASSERT(B.size() >= 2);
	
	maxT = min(maxT, M_PI/dw);
    int numT = int(maxT/dt);
	
    Tirf = VectorXd::LinSpaced(numT, 0, maxT);
    Kirf.resize(numT);
	
	size_t Nf = B.size();
	
	VectorXd y(Nf);
    for (int it = 0; it < numT; ++it) {
		for (int iw = 0; iw < Nf; ++iw)
			y(iw) = B(iw)*cos((w0 + iw*dw)*Tirf(it));
		Kirf(it) = Integral(y, dw, IntegralType::SIMPSON_1_3)*2/M_PI;
	}
}

void GetKirf(VectorXd &Kirf, const VectorXd &w, const VectorXd &B, double dt, double maxT) {
	VectorXd Tirf;
	GetKirfTirf(Kirf, Tirf, w, B, dt, maxT);
}

void GetKirf(VectorXd &Kirf, double w0, double dw, const VectorXd &B, double dt, double maxT) {
	VectorXd Tirf;
	GetKirfTirf(Kirf, Tirf, w0, dw, B, dt, maxT);
}

double GetAinf_Kirf(VectorXd &Kirf, const VectorXd &w, const VectorXd &A, const VectorXd &B, double dt, double maxT) {
	VectorXd Tirf;
	GetKirfTirf(Kirf, Tirf, w, B, dt, maxT);
	return GetAinf(Kirf, Tirf, w, A, dt, maxT);
}

void GetAinf_Kirf(double &Ainf, VectorXd &Kirf, double w0, double dw, const VectorXd &A, const VectorXd &B, double dt, double maxT) {
    VectorXd Tirf;
    GetKirfTirf(Kirf, Tirf, w0, dw, B, dt, maxT);
    
    int numT = int(maxT/dt);
	Eigen::Index Nf = B.size();
	
	VectorXd y(numT);
    Ainf = 0;
    for (int iw = 0; iw < Nf; ++iw) {
        double w = w0 + iw*dw;
        for (int it = 0; it < numT; ++it) 
        	y(it) = Kirf(it)*sin(w*Tirf(it));
        Ainf += A(iw) + Integral(y, dt, IntegralType::SIMPSON_1_3)/w;	// Ogilvie's formula
	}
	Ainf = Ainf/Nf;
}

double GetAinf(const VectorXd &Kirf, const VectorXd &Tirf, const VectorXd &w, const VectorXd &A, double dt, double maxT) {
    int numT = int(maxT/dt);
	Eigen::Index Nf = A.size();
	
	VectorXd y(numT);
    double Ainf = 0;
    for (int iw = 0; iw < Nf; ++iw) {
        for (int it = 0; it < numT; ++it) 
        	y(it) = Kirf(it)*sin(w(iw)*Tirf(it));
        Ainf += A(iw) + Integral(y, dt, IntegralType::SIMPSON_1_3)/w(iw);	// Ogilvie's formula
	}
	return Ainf/Nf;
}


/*double Fradiation2(double t, const VectorXd &vel, const VectorXd &irf, double dt) {
	Eigen::Index numV = int(t/dt);
	if (numV < 2)
		return 0;
	double ret = 0;
	for (Eigen::Index i = numV-2, idtau = 0; i >= 0 && idtau < irf.size()-1; --i, ++idtau) 
		ret += Avg(irf(idtau), irf(idtau+1))*Avg(vel(i), vel(i+1))*dt;
	return ret;
}*/	

double Fradiation(const VectorXd &vel, const VectorXd &irf, Eigen::Index iiter, double dt, Eigen::Index velSize) {
	if (irf.size() == 0)
		return 0;
	if (velSize < 0)
		velSize = vel.size();
	if (velSize < 2)
		return 0;
	
	Eigen::Index num = min(iiter, velSize, irf.size());
	VectorXd vvel = vel.segment(iiter-num, num);
	VectorXd iirf = irf.head(num).reverse();
	VectorXd cont = vvel.array()*iirf.array();
	
	return Integral(cont, dt, SIMPSON_1_3);
}	

double DampedSin(double x, double z0, double zDecay, double mass, double ainf, double b, double w_d, double t0, double phi) {
	double gamma = b/2/(mass + ainf);
	return z0 + zDecay*exp(-gamma*(x-t0))*cos(w_d*(x-t0) + phi);
}

void DampedSin(double x, double z0, double zDecay, double mass, double ainf, double b, double w_d, double t0, double phi,
		double &y, double &dy, double &d2y) {
	double gamma = b/2/(mass + ainf);
	y = z0 + zDecay*exp(-gamma*(x-t0))*cos(w_d*(x-t0) + phi);
	dy = -(zDecay/(2*mass))*exp(-gamma*(x-t0))*(2*mass*w_d*sin(w_d*(x-t0) + phi) + b*cos(w_d*(x-t0) + phi));
	d2y = (zDecay/sqr(2*mass))*exp(-gamma*(x-t0))*(4*b*mass*w_d*sin(w_d*(x-t0) + phi) + (b*b - sqr(2*mass*w_d))*cos(w_d*(x-t0) + phi));
}

double FitToDampedSin(const VectorXd &y, double dt, double mass, double Kh, double g, ParamDampedSin &par) {
	int num = int(y.size());
	
	double w0 = sqrt(Kh/mass);
	double z0 = y.mean();
	double zDecay = y[0] - z0;
	
	int numCoeff = 3 + par.gett0 + par.getz0;
	VectorXd coeff(numCoeff); 
	int idc = 0;
	coeff[idc++] = zDecay;
	coeff[idc++] = par.gamma;
	coeff[idc++] = w0;
	if (par.getz0)
		coeff[idc++] = par.z0;
	if (par.gett0)
		coeff[idc++] = par.t0;
	
	if (!NonLinearOptimization(coeff, num, [&](const VectorXd &coeff, VectorXd &residual)->int {
		for(int i = 0; i < num; i++) {
			double x = i*dt;
			int idc = 0;
			double zDecay = coeff[idc++];
			double gamma  = coeff[idc++];
			double w_d 	  = abs(coeff[idc++]);	
			double z0 	  = par.getz0 ? coeff[idc++] : par.z0;
			double t0 	  = par.gett0 ? coeff[idc++] : par.t0;
			double w02 	  = sqrt(sqr(w_d) + sqr(gamma));
			double ainf   = Kh/sqr(w02) - mass; 
			double b 	  = 2*(mass + ainf)*gamma;
			double phi    = -atan2(b, 2*(mass + ainf)*w02);
			residual[i]   = DampedSin(x, z0, zDecay, mass, ainf, b, w_d, t0, phi) - y[i];
		}
		return 0;	
	}))
		throw Exc("Problem fitting damped sinusoidal");					
	
	idc = 0;
	par.zDecay	= coeff[idc++];
	par.gamma   = coeff[idc++];
	par.w_d     = abs(coeff[idc++]);				// Damped natural frequency
	if (par.getz0)
		par.z0	= coeff[idc++];
	if (par.gett0)
		par.t0	= coeff[idc++];
	par.w02 = sqrt(sqr(par.w_d) + sqr(par.gamma));	// Undamped natural frequency
	par.eta = par.gamma/w0;							// Damping ratio (fraction of the critical damping)
	par.ainf = Kh/sqr(par.w02) - mass;            	// Added mass at infinite frequency
	par.b = 2*(mass + par.ainf)*par.gamma;	
	par.phi = -atan2(par.b, 2*(mass + par.ainf)*par.w02);
	
	VectorXd fitted(num);
	for(int i = 0; i < num; i++)
		fitted[i] = DampedSin(i*dt, par.z0, par.zDecay, mass, par.ainf, par.b, par.w_d, par.t0, par.phi);
	double r2 = R2(y, fitted, Null);
	if (r2 < 0.8)
		throw Exc("Low quality damped sinusoidal fitting");
	
	return r2;
}

void Decay(double mass, double ainf, double av, double Kh, double b, double b2, double dt, double zDecay, double maxT, VectorXd &z ) {
	Upp::Array<VectorXd> res, dres;
	VectorXd y(2), dy(2);
	y[0] = zDecay;
	y[1] = dy[0] = 0;
	dy[1] = -Kh*zDecay/(mass + ainf);
	SolveDAE(y, dy, dt, maxT, res, dres, [&](double t, Eigen::Index iiter, const double y[], const double dy[], double residual[])->bool {
		t = min(t, maxT);
		int iter = int(maxT/dt);
		res[1](iter) = y[1];
		residual[0] = (mass + ainf + av*abs(y[1]))*dy[1] + Kh*y[0] + b*y[1] + b2*y[1]*abs(y[1]);
		residual[1] = y[1] - dy[0];
		return true;
	});
	z = pick(res[0]);	
}

void FitToDecay(const VectorXd &z, const VectorXd &dz, const VectorXd &d2z, 
			double dt, double mass, double Kh, double g, ParamDecay &par) {
				
	int numCoeff = 1 + par.getav + par.getb + par.getb2;
	if (par.numB > 0) {
		ASSERT(par.numB > 2);
		numCoeff += par.numB - 2;
	}
	VectorXd coeff = VectorXd::Constant(numCoeff, 1.);
	int idc = 0;
	coeff[idc++] = par.ainf; 	
	if (par.getav)	
		coeff[idc++] = par.av; 
	if (par.getb) 	
		coeff[idc++] = par.b; 		
	if (par.getb2) 
		coeff[idc] = par.b2; 	
	if (par.numB > 0) {
		par.wdec = VectorXd::LinSpaced(par.numB, 0, par.maxFreq);
		par.Bdec.resize(par.numB);
		par.Bdec(0) = par.Bdec(par.numB-1) = 0;
		par.wspl = VectorXd::LinSpaced(par.numB, 0, par.maxFreq);
		par.Bspl.resize(par.wspl.size());
	}

	if (!NonLinearOptimization(coeff, z.size(), [&](const VectorXd &x, VectorXd &err)->int {
		int idc = 0;
		double ainf = abs(x[idc++]);
		double av   = par.getav ? abs(x[idc++]) : 0;
		double b    = par.getb  ? x[idc++] : 0;
		double b2   = par.getb2 ? x[idc++] : 0; 
		if (par.numB > 0) {
			for (int i = 0; i < par.numB-2; ++i)
				par.Bdec(i+1) = abs(x[i+idc]);

			par.splineB.Init(par.wdec, par.Bdec);
			for (int i = 0; i < par.wspl.size(); ++i) 
				par.Bspl(i) = par.splineB.f(par.wspl[i]);

			GetKirf(par.Kirf, par.wspl, par.Bspl, dt);
		}
		for(Eigen::Index i = 0; i < z.size(); i++) {
			double Frad = Fradiation(dz, par.Kirf, i, dt);
			err[i] = (mass + ainf + av*abs(dz[i]))*d2z[i] + Frad + Kh*z[i] + b*dz[i] + b2*dz[i]*abs(dz[i]);
		}
		return true;	
	}, Null, Null, 1000))
		Cout() << "\nNLO problem";

	idc = 0;
	par.ainf = abs(coeff[idc++]);	
	par.av   = par.getav ? abs(coeff[idc++]) : 0;
	par.b    = par.getb  ? coeff[idc++] : 0;
	par.b2   = par.getb2 ? coeff[idc]   : 0;
}
