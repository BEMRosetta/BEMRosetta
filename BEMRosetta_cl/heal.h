// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#ifndef _BEMRosetta_cl_heal_h_
#define _BEMRosetta_cl_heal_h_

class HealBEM {
public:
	void Heal(bool zremoval, bool thinremoval, bool decayingTail, bool haskind);
	bool Load(const Eigen::VectorXd &w, const Eigen::VectorXd &A, const Eigen::VectorXd &B, int numT, double maxT, const Eigen::MatrixXd &ex_hf);
	void Save(const Eigen::VectorXd &w, Eigen::VectorXd &A, Eigen::VectorXd &Ainfw, double &ainf, Eigen::VectorXd &B, 
				Eigen::VectorXd &Tirf, Eigen::VectorXd &Kinf);
	void Reset(const Eigen::VectorXd &w, Eigen::VectorXd &A, Eigen::VectorXd &Ainfw, double &ainf, Eigen::VectorXd &B, 
				Eigen::VectorXd &Tirf, Eigen::VectorXd &Kinf);
			   				
	Eigen::VectorXd w, A, B;
	Eigen::VectorXd fB, fA, fAinf;
	double fainf;
	static String filterType;
	Eigen::VectorXd Tirf, fKirf;
	Eigen::MatrixXd ex_hf;
	
	String title;
	String sdof;

	UVector<Pointf> aoi, maxB;
	
	// Parameters for IRF calculation	
	double maxT;
	int numT;
	
	// 1. AreaOfInterest
	double  aoix0, 			// Value x from which the AOI begins
			aoidx,			// AOI width 
			aoidy;			// AOI height (aoiy0 = 0)
	int idaoix0,			// index of aoix0
		idaoixMx,			// index of aoix0 + aoidx
		idaoiyMx;			// Index of the max value in y
		
private:
	bool AreaOfInterest(double percMin, 
						double percMax, 
						double &aoix0, 		// Value x from which the AOI begins
						double &aoidx, 		// AOI width 
						double &aoidy, 		// AOI height (aoiy0 = 0)
						int &idaoix0, 		// index of aoix0
						int &idaoixMx, 		// index of aoix0 + aoidx
						int &idaoiyMx);
					 
	Upp::Index<int> SpineRemovalRight(int idpk, double maxDer);
	Upp::Index<int> SpineRemovalLeft(int idpk, double maxDer);
	void SpineRemoval(int idpk, double maxDer, UVector<bool> &idToVoid);
	
	static void CubicFromEnds(double x0, double y0, double p0, double x1, double y1, double p1,
				   			double &a, double &b, double &c, double &d);
	static void ScrimTape(Eigen::VectorXd &x, Eigen::VectorXd &y, int from, int to);
};


#endif
