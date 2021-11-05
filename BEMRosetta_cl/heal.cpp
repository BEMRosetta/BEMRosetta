// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2021, the BEMRosetta author and contributors
#include <Core/Core.h>
#include <Functions4U/Functions4U.h>
#include <STEM4U/Integral.h>
#include <STEM4U/Sundials.h>
#include <STEM4U/Butterworth.h>
#include <STEM4U/LocalFitting.h>
#include <STEM4U/Utility.h>

using namespace Upp;
using namespace Eigen;

#include <BEMRosetta_cl/BEMRosetta.h>
#include <BEMRosetta_cl/functions.h>
#include <BEMRosetta_cl/heal.h>


String HealBEM::filterType = "full";

// Gets the x and y limits of the area of interest of B(w). The obtained rectangle will serve
// as a reference for all the algorithms to be implemented later
bool HealBEM::AreaOfInterest(
	double percMin, 
	double percMax, 
	double &aoix0, 		// Value x from which the AOI begins
	double &aoidx, 		// AOI width 
	double &aoidy, 		// AOI height (aoiy0 = 0)
	int &idaoix0, 		// index of aoix0
	int &idaoixMx, 		// index of aoix0 + aoidx
	int &idaoiyMx) {	// Index of the max value in y
				
	VectorXd accum(B.size());
	double srate = w[1] - w[0];
	
	accum[0] = B[0];								// Accumulated values per index
	for (int i = 1; i < B.size(); ++i) 
		accum[i] = accum[i-1] + B[i];
	
	double tot = accum[accum.size()-1];
	
	idaoix0 = Null;
	for (int i = 0; i < B.size(); ++i) {			// The area of interest goes
		if (IsNull(idaoix0) && accum[i] > tot*percMin)	// from minimum percentile	
			idaoix0 = i;
		if (accum[i] > tot*percMax) {					// to maximum percentile
			idaoixMx = i;								// idaoixMx initial guess
			break;
		}
	}
	if (IsNull(idaoix0))
		return false;
	
	///////////////////////////////////////////////////////
	// REPLACE THIS ALGORITHM WITH STRONG FILTER instead of these averages
	///////////////////////////////////////////////////////
	
	// Second pass to calculate idaoixMx by somewhat reducing the effect of the mini-peaks on the right-hand side.
	double maxB = B.maxCoeff();						// Preliminary (unfiltered) max
	int iddeltax = int(w.size()/20);				// Divide x (w) in 20 sections
	for (int i = idaoixMx; i >= idaoix0; i -= iddeltax) {				// From initial idaoixMx guess
		int idfrom = max(0, i-iddeltax/2);
		int iddeltaxx = min(iddeltax, int(B.size()) - idfrom - 1);
		double avg = B.segment(idfrom, iddeltaxx).mean();	// Section avg
		if (avg > maxB*0.05) {											// If section avg is relevant (> 5% Max)
			idaoixMx = i;												// it may be the max relevant x
			break;														
		}
	}
	
	aoix0 = w[idaoix0];
	aoidx = w[idaoixMx] - w[idaoix0];
	
	VectorXd fB;
	double Tcut = aoidx*0.2;		// The maximum value is taken from the strongly filtered series
	VectorXd num, den;
	ButterLowPass(3, 2*srate/Tcut, num, den);
	Filtfilt(B, num, den, fB);
	
	aoidy = fB.maxCoeff(&idaoiyMx);	// Index of the maximum
	
	return true;
}

// From the peak to the right, it detects the window were the slope is higher than maxDer 
Upp::Index<int> HealBEM::SpineRemovalRight(int idpk, double maxDer) {
	Upp::Index<int> ret;	
	while (true) {
		int k;
		for (k = 1; k < 8 && idpk+k < w.size(); ++k) {
			double der = abs((B[idpk]-B[idpk+k])/(w[idpk]-w[idpk+k]));
			if (der > maxDer) {
				for (int l = 0; l <= k; ++l)
					ret.FindAdd(idpk+l);
				idpk += k;
			} else
				return ret;
		}
		if (idpk+k >= w.size())
			return ret;
	}
}

// From the peak to the left, it detects the window were the slope is higher than maxDer 
Upp::Index<int> HealBEM::SpineRemovalLeft(int idpk, double maxDer) {
	Upp::Index<int> ret;	
	while (true) {
		int k;
		for (k = 1; k < 8 && idpk-k >= 0; ++k) {
			double der = abs((B[idpk]-B[idpk-k])/(w[idpk]-w[idpk-k]));
			if (der > maxDer) {
				for (int l = 0; l <= k; ++l)
					ret.FindAdd(idpk-l);
				idpk -= k;
			} else
				return ret;
		}
		if (idpk-k < 0)
			return ret;
	}
}

// Removes both sides of a spine
void HealBEM::SpineRemoval(int idpk, double maxDer, Vector<bool> &idToVoid) {
	int idpk0 = max(0, idpk - 3);
	int idpkE = min(int(w.size())-1, idpk + 3);	
	Upp::Index<int> maxLeft, maxRight;
	double maxheightLeft = 0, maxheightRight = 0;
	for (int id = idpk0; id <= idpkE; ++id) {
		double height;
		Upp::Index<int> left  = SpineRemovalLeft (id, maxDer);
		if (left.size() > 0) {
			height = abs(B[left[0]] - B[left.Top()]);
			if (height > maxheightLeft) {
				maxLeft  = pick(left);
				maxheightLeft = height;
			}
		}
		Upp::Index<int> right = SpineRemovalRight(id, maxDer);
		if (right.size() > 0) {
			height = abs(B[right[0]] - B[right.Top()]);
			if (height > maxheightRight) {
				maxRight = pick(right);
				maxheightRight = height;
			}
		}
	}
	double maxheight = min(maxheightLeft, maxheightRight);		// Both sides of spine have to have the same length
	if (maxheight == 0)
		return;
	
	idToVoid[idpk] = false;
		
	int numleft = 0;
	for (int i = idpk; i >= 0 && abs(B[idpk] - B[i]) <= maxheight; --i) // Effectively signs to remove both sides of the spine
		numleft++;
	numleft = int(1.5*numleft);
	for (; numleft > 0 && idpk - numleft >= 0; numleft--)
		idToVoid[idpk - numleft] = false;
	
	int numright = 0;
	for (int i = idpk; i < B.size() && abs(B[idpk] - B[i]) <= maxheight; ++i) 
		numright++;
	numright = int(2*numright);
	for (; numright > 0 && idpk + numright < B.size(); numright--)
		idToVoid[idpk + numright] = false;
}

// Based on two points coordinates and the slope on them, ot returs the coefficients of a cubic
// complying x, y and slope in both ends
void HealBEM::CubicFromEnds(double x0, double y0, double p0, double x1, double y1, double p1,
				   			double &a, double &b, double &c, double &d) {
	double x0_2 = x0*x0; 
	double x0_3 = x0_2*x0;
	double x1_2 = x1*x1; 
	double x1_3 = x1_2*x1;
	double deltax = x0-x1;
	double den = pow(deltax, 3);
	a = (p0*deltax + p1*deltax + 2*(y1-y0))/den;
	b = (-p0*(x0_2+x0*x1-2*x1_2) + p1*(-2*x0_2+x0*x1+x1_2) + 3*(x0+x1)*(y0-y1))/den;
	c = (p1*x0*(x0_2+x0*x1-2*x1_2) - x1*(p0*(-2*x0_2+x0*x1+x1_2) + 6*x0*(y0-y1)))/den;
	d = (p0*x0*x1_2*(x1-x0) + p1*x0_2*x1*(x1-x0) + x0_3*y1 - 3*x0_2*x1*y1 + 3*x0*x1_2*y0 - x1_3*y0)/den;
}


// Some scrimtape maintaining coordinates and slope of patch ends
void HealBEM::ScrimTape(VectorXd &x, VectorXd &y, int from, int to) {
	if (from <= 0)
		from = 1;
	if (to >= x.size()-1)
		to = int(x.size())-2;
	if (from >= to)
		return;
	double p0 = (y[from]-y[from-1])/(x[from]-x[from-1]);
	double p1 = (y[to+1]-y[to])/(x[to+1]-x[to]);
	double a, b, c, d;
	CubicFromEnds(x[from], y[from], p0, x[to], y[to], p1, a, b, c, d);
	for (int j = from+1; j < to; ++j) {
        double s = x[j];
        double s_2 = s*s;
		y[j] = a*s*s_2 + b*s_2 + c*s + d;
    }
}
				   

bool IsNull(const VectorXd &data) {
	if (data.size() == 0)
		return true;
	
	for (int i = 0; i < data.size(); ++i)
		if (!IsNum(data[i]))
			return true;
	return false;
}

bool HealBEM::Load(const VectorXd &w, const VectorXd &A, const VectorXd &B, int numT, double maxT) {
	if (A.size() == 0 || B.size() == 0)
		return false;
	this->w = clone(w);
	this->B = clone(B);
	this->A = clone(A);
	this->numT = numT;
	this->maxT = maxT;
	return true;
}

void HealBEM::Reset(const VectorXd &w, VectorXd &A, VectorXd &Ainfw, double &ainf, VectorXd &B, 
		VectorXd &Tirf, VectorXd &Kirf) {
	A.resize(0);
	Ainfw.resize(0);
	ainf = 0;
	B.resize(0);
	Kirf.resize(0);
}
	
void HealBEM::Save(const VectorXd &w, VectorXd &A, VectorXd &Ainfw, double &ainf, VectorXd &B, 
		VectorXd &Tirf, VectorXd &Kirf) {
	Tirf = this->Tirf;
	Kirf = this->fKirf;
	ainf = this->fainf;
	
	for (int iw = 0; iw < w.size(); ++iw) {
		A(iw)  = LinearInterpolate(w[iw], this->w, this->fA);
		Ainfw(iw) = LinearInterpolate(w[iw], this->w, this->fAinf);
		B(iw)  = LinearInterpolate(w[iw], this->w, this->fB);
	}	
}
	
void HealBEM::Heal(bool zremoval, bool thinremoval, bool decayingTail) {
	// Removes NaN, Inf, duplicated (or nearly) w, sorts by w 
	CleanNANDupXSort(w, A, B, w, A, B);
	double srate = GetSampleRate(w, 4, .8);	// Gets the most probable sample rate, or the average if the most probable probability is lower than 0.8
	Resample(w, A, B, w, A, B, srate);	
	
	
	// Filtered process, with irregular frequencies removal. 
	
	// 1. AreaOfInterest
	if (!AreaOfInterest(0.01, 0.98, aoix0, aoidx, aoidy, idaoix0, idaoixMx, idaoiyMx))
		return;


	// Removes spikes	
	// Different ways
	if (filterType == "butterworth") {		  // Just a filter
		double Tcut = aoidx*0.2;				// Strong filter
		
		VectorXd num, den;
		ButterLowPass(3, 2*srate/Tcut, num, den);
		Filtfilt(B, num, den, fB);
	} else if (filterType == "localfitting") {// Just a filter
		double Tcut = aoidx*0.08;				// Light filter
		
		VectorXd num, den;
		ButterLowPass(3, 2*srate/Tcut, num, den);
		Filtfilt(B, num, den, fB);
		
		VectorXd nw, nB;
		CleanOutliers(w, B, fB, nw, nB, 0.4);	// Removed values +- 40 out of filtered B
		CleanCondition(nw, nB, nw, nB, [&](int id) {return nB[id] > 0;});	// Removed negatives
		
		LocalFitting(nw, nB, w, fB, 2, aoidx*0.25, false);		// Local fitting, quadratic
	} else if (filterType == "full") {// The complex way
		
		Vector<bool>idToRemove;						// Signals with false if a value of B should be removed
				
		Resize(idToRemove, int(B.size()), true);	// Initialisation, resets points to be cleaned
	
		if (zremoval) {		// Removes sign of Zorro. Basic algorithm, but functional
			Vector<int64> allupperPk, alllowerPk;
			FindPeaks(B, allupperPk, alllowerPk);		// Finds all lower and upper peaks, even the tiniest
			
			Vector<int> idZ;
			Vector<double> heightZ;
			
			for (int ip = 0; ip < allupperPk.size(); ++ip) {	// Max peaks
				int i = int(allupperPk[ip]);	// Detects too much step between a peak and adjacent values
				if (abs(B[i] - B[i+1]) > 0.2*aoidy || abs(B[i] - B[i-1]) > 0.2*aoidy) {	// 20% of aoidy
					idZ << i;
					heightZ << max(abs(B[i] - B[i+1]), abs(B[i] - B[i-1]));
				}
			}
			for (int ip = 0; ip < alllowerPk.size(); ++ip) {	// Min peaks
				int i = int(alllowerPk[ip]);	// Detects too much step between a peak and adjacent values
				if (abs(B[i] - B[i+1]) > 0.2*aoidy || abs(B[i] - B[i-1]) > 0.2*aoidy) {	// 20% of aoidy
					idZ << i;
					heightZ << max(abs(B[i] - B[i+1]), abs(B[i] - B[i-1]));
				}
			}	
			// Bad peak detected. Now it's chosen window to clean around
			for (int ip = 0; ip < idZ.size(); ++ip) {
				int i = idZ[ip];
				idToRemove[i] = false;	// Removed, but total width to remove depends on peak height until a max
				double deltaw =  min(0.2*aoidx, 0.1*(heightZ[ip]/aoidy)*aoidx); // Window to clean
				for (int j = i+1; j < w.size() && w[j] - w[i] < deltaw; ++j)
					 idToRemove[j] = false;
				for (int j = i-1; j >= 0 && w[i] - w[j] < deltaw; --j)
					 idToRemove[j] = false;
			}
			// Recalculates aoidy and idaoiyMx after sign of Zorro removal
			aoidy = 0;
			for (int i = 0; i < w.size(); ++i) {
				if (w[i] > aoix0 && w[i] < aoix0 + aoidx) {
					if (idToRemove[i] && B[i] > aoidy) {
						idaoiyMx = i;
						aoidy = B[i];
					}
				}
			}
		}
		// Removes spines

		if (thinremoval) {	// This option can be avoided
			double maxDer = 6*abs(aoidy/aoidx); // Max slope. Less sharp spine detection to get and remove bear claws
			
			Vector<int64> upperPk, lowerPk;
			FindPeaks(w, B, 0.1*aoidx, upperPk, lowerPk);	// Considers 0.1*aoidx as the window sensitivity to detect
			
			// Bear claws are removed after the max value idaoiyMx, not before
			for (int i = 0; i < upperPk.size(); ++i)
				if (upperPk[i] > idaoiyMx+3 && idToRemove[int(upperPk[i])]) {	// Only removed after max
					SpineRemoval(int(upperPk[i]), maxDer, idToRemove);
					//upper << Pointf(w[upperPk[i]], B[upperPk[i]]);
				}
			
			for (int i = 0; i < lowerPk.size(); ++i) 
				if (lowerPk[i] > idaoiyMx+3 && idToRemove[int(lowerPk[i])]) {	// Only removed after max
					SpineRemoval(int(lowerPk[i]), maxDer, idToRemove);	
					//lower << Pointf(w[lowerPk[i]], B[lowerPk[i]]);
				}
		}
		
		// Cleans negatives
		
		for (int i = 0; i < B.size(); ++i)
			if (B[i] < 0)
				idToRemove[i] = false;		
		
		
		// The points to be deleted have been detected, but a new B has to be created
	
		if (false) {		// Just shows the results after spine removal. No new B creation			
			fB = B;
			for (int i = 0; i < B.size(); ++i) 
				if (!idToRemove[i])
					fB[i] = Null;	// Null
		} else {
			// Vectors to signal the begin and end of each area to be cleaned. size() is the
			// number of scratches
			Vector<int> fromScratch,		// Id of the first value of each scratch 
						toScratch;			// Id of the list value of each scratch
			for (int i = 0; i < idToRemove.size(); ++i) {
				bool equal = fromScratch.size() == toScratch.size();
				if (!idToRemove[i]) {
					if (equal)
						fromScratch << i;
				} else {
					if (!equal)
						toScratch << i-1;
				}
			}
			if (fromScratch.size() != toScratch.size())
				toScratch << idToRemove.size() - 1;
			
			
			// The base filtered B fB is a soft version of B
			VectorXd num, den;
			double Tcutsoft = aoidx*0.1;
			ButterLowPass(3, 2*srate/Tcutsoft, num, den);
			Filtfilt(B, num, den, fB);
			
			// f2B is a more strongly filtered version
			VectorXd f2B;
			double Tcuthard = aoidx*0.4;
			ButterLowPass(3, 2*srate/Tcuthard, num, den);
			Filtfilt(B, num, den, f2B);			
			
			for (int i = 0; i < fromScratch.size(); ++i) {	// For each scratch
				int ifrom = max(0, fromScratch[i]-1);
				int ito = min(int(B.size())-1, toScratch[i]+1);
				
				double x0 = w[ifrom];		// w begin of scratch
				double y0 = fB[ifrom];
				double x1 = w[ito];			// w end of scratch
				double y1 = fB[ito];

				if (x0 >= aoix0 + aoidx) 	// Zone out of aoi is not fixed
					break;
									
				if ((x1 - x0) < aoidx*0.1) {// Small injuries are solved with a plaster
					for (int j = fromScratch[i]; j <= toScratch[i]; ++j) {
						double factorx = (w[j] - x0)/(x1 - x0);
						fB[j] = y0 + (y1 - y0)*factorx;	// made with a linear connection between scratch ends
					}
				/*} else {
					double f2y0 = f2B[ifrom];
					double f2yE = f2B[ito];
					
					for (int j = fromScratch[i]; j <= toScratch[i]; ++j) {
						double factorx = (w[j] - x0)/(x1 - x0);
						
						// Linear morphing of f2B
						fB[j] = f2B[j] - (f2y0 - y0) - factorx*((y0-y1) - (f2y0-f2yE));
					}
			*/  } else {						// Big injuries require a graft
					double a, b, c, d;
					
					double p0 = 0;
					int num = 0;
					for (int j = fromScratch[i] - 1; j >= 0 && num < 4; --j) {
						p0 += (fB[j+1]-fB[j])/(w[j+1]-w[j]);
						p0 += (f2B[j+1]-f2B[j])/(w[j+1]-w[j]);
						num++;
					}
					p0 /= (2*num);	// Begin slope is the avg of the slope in 4 previous points for fB and f2B
					double p1 = 0;	// End slope is zero (Could this be improved?)
					
					CubicFromEnds(x0, y0, p0, x1, y1, p1, a, b, c, d);	// Cubic from end points and slope
					for (int j = fromScratch[i]; j <= toScratch[i]; ++j) {
				        double x = w[j];
				        double x_2 = x*x;
						fB[j] = a*x*x_2 + b*x_2 + c*x + d;
				    }
				}
				// Some scrimtape +- 3 points around patch
				ScrimTape(w, fB, fromScratch[i]-3, fromScratch[i]+3);	// Begin of the repaired scratch
				ScrimTape(w, fB, toScratch[i]-3, toScratch[i]+3);		// End of the repaired scratch
			}
	 
	 		if (decayingTail && aoidx > 0) {		// Decaying right tail
	 			int fromid, toid;
	 			double fromw = aoix0 + 1*aoidx;		// Begin and 
	 			double tow   = aoix0 + 1.5*aoidx;	// end of decaying. After end, B is zero

				RealExponentEquation real;		// Decaying function B(w) = a*w^b
	 			
	 			while (true) {
		 			for (int i = 0; i < w.size() && w[fromid = i] < fromw; ++i) 		// id of fromw
		 				;
		 			for (int i = fromid; i < w.size() && w[toid = i] < tow; ++i) 		// id of tow
		 				;	
		 			
		 			if (fromid == toid)
		 				break;
		 			
		 			VectorXd fitw =  w.segment(fromid, toid-fromid+1);		// Gets the fragment of original B to be fit
		 			VectorXd fitB = fB.segment(fromid, toid-fromid+1);
		 	
		 			EigenVector fit(fitw, fitB);	// Eigen is the name of the library to do the fitting		
					
					double r2;
					real.Fit(fit, r2);				// Fitting
					if (!IsNull(r2)) {				// Fitting OK
						title += Format("R2:%.3f", r2);
						break;
					}
					
					fromw += 0.1*aoidx;				// No fitting..., B is not decaying. Repeat moving fromw to the right
	 				tow   += 0.1*aoidx;
	 			}
								
				for (int i = idaoixMx; i < w.size(); ++i)
					fB[i] = real.f(w[i]);			// Replaces old with decaying values from B(w) = a*w^b 
				
				ScrimTape(w, fB, idaoixMx-3, idaoixMx+3);	// Scrimtape to fix the patch +- 3 points around
	 		}
		}

		// Interval to get Ainf
		int fromA;
		for (int i = 0; i < w.size() && w[fromA = i] < 0.1*aoidx; ++i) 	// From 0.1*aoidx
		 	;
		int toA = idaoiyMx;												// To the max value idaoiyMx
		
		
		// IRF and Ainf obtained from filtered B fB. A used directly (should it be softly
		// filtered?
		GetTirf(Tirf, numT, maxT);
		GetKirf(fKirf, Tirf, w, fB);
		GetAinf_w(fAinf, fKirf, Tirf, w, A);
		
		if (toA > fromA) {
			VectorXd tmp = fAinf.segment(fromA, toA-fromA);
			fainf = tmp.mean();		// New clean Ainf
			
			double dt = maxT/(numT-1);
			// New fA obtained from ainf and fB
			GetA(fA, fKirf, w, fainf, dt);
			
			fAinf = VectorXd::Constant(w.size(), fainf);
		} else 
			fA = A;
	}
}

 
