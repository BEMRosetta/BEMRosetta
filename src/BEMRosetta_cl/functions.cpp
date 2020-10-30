#include "BEMRosetta.h"
#include <STEM4U/Integral.h>

#include "functions.h"

using namespace Eigen;


static void GetKirfTirf(VectorXd &Kirf, VectorXd &Tirf, const VectorXd &w, const VectorXd &B, double dt, double &maxT) {
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
		Kirf(it) = Integral<VectorXd, double>(y, w)*2/M_PI;
	}
}

static void GetKirfTirf(VectorXd &Kirf, VectorXd &Tirf, double w0, double dw, const VectorXd &B, double dt, double &maxT) {
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

void GetAinf_Kirf(double &Ainf, VectorXd &Kirf, const VectorXd &w, const VectorXd &A, const VectorXd &B, double dt, double maxT) {
	VectorXd Tirf;
	GetKirfTirf(Kirf, Tirf, w, B, dt, maxT);

    int numT = int(maxT/dt);
	Eigen::Index Nf = B.size();
	
	VectorXd y(numT);
    Ainf = 0;
    for (int iw = 0; iw < Nf; ++iw) {
        for (int it = 0; it < numT; ++it) 
        	y(it) = Kirf(it)*sin(w(iw)*Tirf(it));
        Ainf += A(iw) + Integral(y, dt, IntegralType::SIMPSON_1_3)/w(iw);	// Ogilvie's formula
	}
	Ainf = Ainf/Nf;
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


double Fradiation2(double t, const VectorXd &vel, const VectorXd &irf, double dt) {
	Eigen::Index numV = int(t/dt);
	if (numV < 2)
		return 0;
	double ret = 0;
	for (Eigen::Index i = numV-2, idtau = 0; i >= 0 && idtau < irf.size()-1; --i, ++idtau) 
		ret += Avg(irf(idtau), irf(idtau+1))*Avg(vel(i), vel(i+1))*dt;
	return ret;
}	

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

