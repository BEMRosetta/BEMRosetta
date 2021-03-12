#ifndef _BEMRosetta_BEMRosetta_cl_functions_h_
#define _BEMRosetta_BEMRosetta_cl_functions_h_

void GetKirfTirf(Eigen::VectorXd &Kirf, Eigen::VectorXd &Tirf, const Eigen::VectorXd &w, const Eigen::VectorXd &B, double dt, double &maxT);
void GetKirfTirf(Eigen::VectorXd &Kirf, Eigen::VectorXd &Tirf, double w0, double dw, 	 const Eigen::VectorXd &B, double dt, double &maxT);	
void GetKirf(Eigen::VectorXd &Kirf, const Eigen::VectorXd &w, const Eigen::VectorXd &B, double dt, double maxT = 30);
void GetAinf_Kirf(double &Ainf, Eigen::VectorXd &Kirf, const Eigen::VectorXd &w, const Eigen::VectorXd &A, const Eigen::VectorXd &B, double dt, double maxT = 30);
void GetKirf(Eigen::VectorXd &Kirf, double w0, double dw, const Eigen::VectorXd &B, double dt, double maxT = 30);
void GetAinf_Kirf(double &Ainf, Eigen::VectorXd &Kirf, double w0, double dw, const Eigen::VectorXd &A, const Eigen::VectorXd &B, double dt, double maxT = 30);
//double Fradiation2(double t, const Eigen::VectorXd &vel, const Eigen::VectorXd &irf, double dt);
double Fradiation(const Eigen::VectorXd &vel, const Eigen::VectorXd &irf, Eigen::Index iiter, double dt, Eigen::Index velSize = -1);
	
#endif
