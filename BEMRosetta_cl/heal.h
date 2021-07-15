#ifndef _BEMRosetta_cl_heal_h_
#define _BEMRosetta_cl_heal_h_

class HealBEM {
public:
	void Heal();
	void Load(const VectorXd &w, const VectorXd &A, const VectorXd &B, double maxT, int num);
	void Save(const VectorXd &w, VectorXd &A, VectorXd &Ainfw, double &ainf, VectorXd &B, 
				VectorXd &Tirf, VectorXd &Kinf);
			   				
	VectorXd w, A, B;
	VectorXd fB, fA, fAinf, fAinfw;
	double fainf;
	static String filterType;
	VectorXd Tirf, fKirf;
	String title;
	String sdof;

	Vector<Pointf> aoi, maxB;
	
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
	void AreaOfInterest(double percMin, 
						double percMax, 
						double &aoix0, 		// Value x from which the AOI begins
						double &aoidx, 		// AOI width 
						double &aoidy, 		// AOI height (aoiy0 = 0)
						int &idaoix0, 		// index of aoix0
						int &idaoixMx, 		// index of aoix0 + aoidx
						int &idaoiyMx);
					 
	Upp::Index<int> SpineRemovalRight(int idpk, double maxDer);
	Upp::Index<int> SpineRemovalLeft(int idpk, double maxDer);
	void SpineRemoval(int idpk, double maxDer, Vector<bool> &idToVoid);
	
	static void CubicFromEnds(double x0, double y0, double p0, double x1, double y1, double p1,
				   			double &a, double &b, double &c, double &d);
	static void ScrimTape(VectorXd &x, VectorXd &y, int from, int to);
};


#endif
