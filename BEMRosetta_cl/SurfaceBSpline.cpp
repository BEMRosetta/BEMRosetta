#include <Core/Core.h>
#include <Surface/Surface.h>
#include "SurfaceBSpline.h"

using namespace Upp;



int FindSpan(double u, const UVector<double>& knots, int n, int p) {
	int m = knots.size() - 1;
	if (u >= knots[m]) 
		return n;
	if (u <= knots[0]) 
		return p;
	
	int low = 0, high = m, mid = (low + high)/2;
	while (u < knots[mid] || u >= knots[mid + 1]) {
		if (u < knots[mid]) 
			high = mid;
		else 
			low = mid;
		mid = (low + high)/2;
	}
	return mid;
}

void BasisFunctions(double u, int span, const UVector<double>& knots, int p, double* N) {
	double left[20], right[20];
	N[0] = 1;
	
	for (int j = 1; j <= p; j++) {
		left[j] = u - knots[span + 1 - j];
		right[j] = knots[span + j] - u;
		double saved = 0.0;
		
		for (int r = 0; r < j; r++) {
			double temp = N[r];
			double denom = right[r+1] + left[j-r];
			if (denom > 1e-10) {
				N[r] = saved + right[r+1] * temp / denom;
				saved = left[j-r] * temp / denom;
			} else {
				N[r] = saved;
				saved = 0.0;
			}
		}
		N[j] = saved;
	}
}

void BasisDerivatives(double u, int span, const UVector<double>& knots, int p, double* dN) {
	// dN[i] = derivative of N_{i,p} at u
	// Uses formula: dN_{i,p}/du = p/(knots[i+p]-knots[i]) * N_{i,p-1} - p/(knots[i+p+1]-knots[i+1]) * N_{i+1,p-1}
	
	double left[20], right[20];
	double N[2][20];  // N[0] for degree p-1, N[1] for degree p
	
	N[0][0] = 1;
	
	for (int j = 1; j <= p; j++) {			// Compute basis functions up to degree p
		left[j] = u - knots[span + 1 - j];
		right[j] = knots[span + j] - u;
		double saved = 0;
		
		for (int r = 0; r < j; r++) {
			double temp = N[(j-1)%2][r];
			double denom = right[r+1] + left[j-r];
			if (denom > 1e-10) {
				N[j%2][r] = saved + right[r+1]*temp/denom;
				saved = left[j-r]*temp/denom;
			} else {
				N[j%2][r] = saved;
				saved = 0;
			}
		}
		N[j%2][j] = saved;
	}
	
	for (int i = 0; i <= p; i++) {		// Now compute derivatives
		double dval = 0.0;
		
		if (i > 0) {					// First term: p/(knots[span+i] - knots[span+i-p]) * N_{span+i-p, p-1}
			int idx1 = span + i;
			int idx2 = span + i - p;
			if (idx2 < 0) idx2 = 0;
			if (idx1 < knots.size()) {
				double denom = knots[idx1] - knots[idx2];
				if (denom > 1e-10)
					dval += p * N[(p-1)%2][i-1] / denom;
			}
		}
		if (i < p) {					// Second term: -p/(knots[span+i+p+1] - knots[span+i+1]) * N_{span+i+1-p, p-1}
			int idx1 = span + i + p + 1;
			int idx2 = span + i + 1;
			if (idx1 >= knots.size()) 
				idx1 = knots.size() - 1;
			double denom = knots[idx1] - knots[idx2];
			if (denom > 1e-10)
				dval -= p * N[(p-1)%2][i]/denom;
		}
		
		dN[i] = dval;
	}
}

Point3D BSplinePatch::EvaluateSurface(const UVector<Point3D>& cpts, int nU, int nV,
                     const UVector<double>& uKnots, const UVector<double>& vKnots,
                     int p, int q, double u, double v) {
	double umin = uKnots[p];
	double umax = uKnots[uKnots.size() - p - 1];
	double vmin = vKnots[q];
	double vmax = vKnots[vKnots.size() - q - 1];
	
	u = clamp(u, umin, umax);
	v = clamp(v, vmin, vmax);
	
	int uspan = FindSpan(u, uKnots, nU - 1, p);
	int vspan = FindSpan(v, vKnots, nV - 1, q);
	
	double Nu[20], Nv[20];
	BasisFunctions(u, uspan, uKnots, p, Nu);
	BasisFunctions(v, vspan, vKnots, q, Nv);
	
	Point3D result(0, 0, 0);
	for (int l = 0; l <= q; l++) {
		int vidx = clamp(vspan - q + l, 0, nV - 1);
		for (int k = 0; k <= p; k++) {
			int uidx = clamp(uspan - p + k, 0, nU - 1);
			result += cpts[vidx * nU + uidx]*(Nu[k]*Nv[l]);
		}
	}
	return result;
}

void BSplinePatch::EvaluateDerivatives(double u, double v, Point3D& Su, Point3D& Sv) const {
	double du = 0.001*(GetUMax() - GetUMin());
	double dv = 0.001*(GetVMax() - GetVMin());
	
	Point3D pu  = Evaluate(u - du, v);
	Point3D puu = Evaluate(u + du, v);
	Su = (puu - pu)/(2*du);
	
	Point3D pv  = Evaluate(u, v - dv);
	Point3D pvv = Evaluate(u, v + dv);
	Sv = (pvv - pv)/(2*dv);
}

double BSplinePatch::EvaluateAreaElement(double u, double v) const {
	Point3D Su, Sv;
	EvaluateDerivatives(u, v, Su, Sv);
	return cross(Su, Sv).Length();
}

double BSplinePatch::ComputeArea(int order) const {
	static const double gp[5][5] = {		// Gauss-Legendre points
		{0},
		{-0.577350269189626, 0.577350269189626},
		{-0.774596669241483, 0.0, 0.774596669241483},
		{-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053},
		{-0.906179845938664, -0.538469310105683, 0.0, 0.538469310105683, 0.906179845938664}
	};
	static const double gw[5][5] = {		// Gauss-Legendre weights
		{2},
		{1.0, 1.0},
		{0.555555555555556, 0.888888888888889, 0.555555555555556},
		{0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454},
		{0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189}
	};
	
	order = clamp(order, 1, 5);
	int n = order;
	
	double umin = GetUMin(), umax = GetUMax();
	double vmin = GetVMin(), vmax = GetVMax();
	double area = 0;
	
	for (int j = 0; j < n; j++) {
		double v = 0.5*((vmax - vmin)*gp[n-1][j] + (vmax + vmin));
		double wv = gw[n-1][j];
		
		for (int i = 0; i < n; i++) {
			double u = 0.5*((umax - umin)*gp[n-1][i] + (umax + umin));
			double wu = gw[n-1][i];
			
			area += wu*wv*EvaluateAreaElement(u, v);
		}
	}
	return 0.25*(umax - umin) * (vmax - vmin)*area;	// Jacobian of transformation
}

void BSplinePatch::GetBoundingBox(Point3D& mn, Point3D& mx) const {
	if (controlPoints.IsEmpty()) {
		mn = mx = Point3D(); 
		return; 
	}
	mn = mx = controlPoints[0];
	for (const auto& cp : controlPoints) {
		mn.x = min(mn.x, cp.x); mx.x = max(mx.x, cp.x);
		mn.y = min(mn.y, cp.y); mx.y = max(mx.y, cp.y);
		mn.z = min(mn.z, cp.z); mx.z = max(mx.z, cp.z);
	}
}

UVector<double> MakeUniformKnots(int nSpans, int order) {
	int nKnots = nSpans + 2*order - 1;
	UVector<double> knots;
	knots.SetCount(nKnots);
	for (int i = 0; i < nKnots; i++) {
		if (i < order) 
			knots[i] = -1;
		else if (i >= nSpans + order - 1) 
			knots[i] = 1;
		else 
			knots[i] = -1. + 2.*(i - order + 1.)/nSpans;
	}
	return knots;
}

/* 	R	Radius of cylinder
	H	Height of cylinder
	nU	Number of panel subdivisions (spans/segments) in circumferential direction
	nV	Number of panel subdivisions in vertical direction
	kU	B-spline order in circumferential direction (4 = cubic)
	kV	B-spline order in vertical direction (4 = cubic)
*/
SurfaceBSpline MakeCylinder(double R, double H, int nU = 4, int nV = 4, int kU = 4, int kV = 4) {
	SurfaceBSpline mgr;
	BSplinePatch patch;

	patch.nug = nU;
	patch.nvg = nV;
	patch.kug = kU;
	patch.kvg = kV;
	patch.uKnots = MakeUniformKnots(nU, kU);
	patch.vKnots = MakeUniformKnots(nV, kV);
	
	int nBasisU = patch.GetNUBasis();
	int nBasisV = patch.GetNVBasis();
	patch.controlPoints.SetCount(nBasisU*nBasisV);
	
	for (int j = 0; j < nBasisV; j++) {
		double v = -1. + 2.*j/(nBasisV - 1.);
		double z = v*H/2.;
		for (int i = 0; i < nBasisU; i++) {
			double theta = 2.*M_PI*i/nBasisU;
			patch.controlPoints[j*nBasisU + i] = Point3D(R*cos(theta), R*sin(theta), z);
		}
	}
	mgr.Add(pick(patch));
	return mgr;
}

/* 	R	Radius of sphere
	nU	Number of panel subdivisions in circumferential direction
	nV	Number of panel subdivisions in vertical direction
	kU	B-spline order in circumferential direction (4 = cubic)
	kV	B-spline order in vertical direction (4 = cubic)
*/
SurfaceBSpline MakeSphere(double R, int nU = 6, int nV = 6, int kU = 4, int kV = 4) {
	SurfaceBSpline mgr;
	BSplinePatch patch;
	
	patch.nug = nU;
	patch.nvg = nV;
	patch.kug = kU;
	patch.kvg = kV;
	patch.uKnots = MakeUniformKnots(nU, kU);
	patch.vKnots = MakeUniformKnots(nV, kV);
	
	int nBasisU = patch.GetNUBasis();
	int nBasisV = patch.GetNVBasis();
	patch.controlPoints.SetCount(nBasisU*nBasisV);
	
	for (int j = 0; j < nBasisV; j++) {
		double v = -1. + 2.*j/(nBasisV - 1.);
		double phi = M_PI*(v + 1.)/2.;
		for (int i = 0; i < nBasisU; i++) {
			double u = -1. + 2.*i/(nBasisU - 1.);
			double theta = M_PI*(u + 1.0);
			patch.controlPoints[j * nBasisU + i] = Point3D(
				R*sin(phi)*cos(theta),
				R*sin(phi)*sin(theta),
				R*cos(phi)
			);
		}
	}
	
	mgr.Add(pick(patch));
	return mgr;
}

/* a	Semi-axis in X direction (length/2)
	b	Semi-axis in Y direction (width/2)
	c	Semi-axis in Z direction (height/2)
*/
SurfaceBSpline MakeEllipsoid(double a, double b, double c, int nU = 6, int nV = 6, int kU = 4, int kV = 4) {
	SurfaceBSpline mgr;
	BSplinePatch patch;
	
	patch.nug = nU;
	patch.nvg = nV;
	patch.kug = kU;
	patch.kvg = kV;
	patch.uKnots = MakeUniformKnots(nU, kU);
	patch.vKnots = MakeUniformKnots(nV, kV);
	
	int nBasisU = patch.GetNUBasis();
	int nBasisV = patch.GetNVBasis();
	patch.controlPoints.SetCount(nBasisU*nBasisV);
	
	for (int j = 0; j < nBasisV; j++) {
		double v = -1. + 2.*j/(nBasisV - 1.);
		double phi = M_PI*(v + 1.)/2.;
		double sinPhi = sin(phi);
		double cosPhi = cos(phi);
		for (int i = 0; i < nBasisU; i++) {
			double u = -1. + 2.*i/(nBasisU - 1.);
			double theta = M_PI*(u + 1.);
			patch.controlPoints[j * nBasisU + i] = Point3D(
				a*sinPhi*cos(theta),
				b*sinPhi*sin(theta),
				c*cosPhi
			);
		}
	}
	
	mgr.Add(pick(patch));
	return mgr;
}

/* 	L	Length of hull (longitudinal, X direction)
	B	Beam (maximum width, Y direction)
	D	Draft (depth below waterline, Z direction)
*/
SurfaceBSpline MakeWigley(double L, double B, double D, int nU = 8, int nV = 4, int kU = 4, int kV = 4) {
	SurfaceBSpline mgr;
	BSplinePatch patch;
	
	patch.nug = nU;
	patch.nvg = nV;
	patch.kug = kU;
	patch.kvg = kV;
	patch.uKnots = MakeUniformKnots(nU, kU);
	patch.vKnots = MakeUniformKnots(nV, kV);
	
	int nBasisU = patch.GetNUBasis();
	int nBasisV = patch.GetNVBasis();
	patch.controlPoints.SetCount(nBasisU*nBasisV);
	
	for (int j = 0; j < nBasisV; j++) {
		double v = -1. + 2.*j/(nBasisV - 1);
		double z = v*D/2. - D/2.;
		double zFactor = 1. - (z/D)*(z/D);
		for (int i = 0; i < nBasisU; i++) {
			double u = -1. + 2.*i/(nBasisU - 1);
			double x = u*L/2.;
			double xFactor = 1. - (2.*x/L)*(2.*x/L);
			double y = (B/2.)*xFactor*zFactor;
			patch.controlPoints[j*nBasisU + i] = Point3D(x, y, z);
		}
	}
	
	mgr.Add(pick(patch));
	return mgr;
}

/* Barge with rounded bilge (3 patches)
	L	Length of barge (longitudinal, X direction)	100.0
	B	Beam (total width, Y direction)	20.0
	D	Draft (depth below waterline, Z direction)	10.0
	R	Bilge radius (rounded corner at bottom)
*/
SurfaceBSpline MakeBarge(double L, double B, double D, double R, int nU = 6, int nV = 4, int kU = 4, int kV = 4) {
	SurfaceBSpline mgr;
	
	{	// Patch 0: Flat bottom
		BSplinePatch patch;
		
		patch.nug = nU;
		patch.nvg = nV;
		patch.kug = kU;
		patch.kvg = kV;
		patch.uKnots = MakeUniformKnots(nU, kU);
		patch.vKnots = MakeUniformKnots(nV, kV);
		
		int nBasisU = patch.GetNUBasis();
		int nBasisV = patch.GetNVBasis();
		patch.controlPoints.SetCount(nBasisU*nBasisV);
		
		for (int j = 0; j < nBasisV; j++) {
			for (int i = 0; i < nBasisU; i++) {
				double u = -1. + 2.*i/(nBasisU - 1);
				double x = u*L/2.;
				patch.controlPoints[j*nBasisU + i] = Point3D(x, 0., -D + R);
			}
		}
		mgr.Add(pick(patch));
	}
	{	// Patch 1: Bilge (rounded)
		BSplinePatch patch;
		
		patch.nug = nU;
		patch.nvg = nV;
		patch.kug = kU;
		patch.kvg = kV;
		patch.uKnots = MakeUniformKnots(nU, kU);
		patch.vKnots = MakeUniformKnots(nV, kV);
		
		int nBasisU = patch.GetNUBasis();
		int nBasisV = patch.GetNVBasis();
		patch.controlPoints.SetCount(nBasisU*nBasisV);
		
		for (int j = 0; j < nBasisV; j++) {
			double v = -1. + 2.*j/(nBasisV - 1);
			double theta = M_PI/2.*(v + 1.)/2.;
			double y = R*(1. - cos(theta));
			double z = -D + R + R*sin(theta);
			for (int i = 0; i < nBasisU; i++) {
				double u = -1. + 2.*i/(nBasisU - 1);
				double x = u*L/2.;
				patch.controlPoints[j*nBasisU + i] = Point3D(x, y, z);
			}
		}
		mgr.Add(pick(patch));
	}
	{	// Patch 2: Side
		BSplinePatch patch;
		
		patch.nug = nU;
		patch.nvg = nV;
		patch.kug = kU;
		patch.kvg = kV;
		patch.uKnots = MakeUniformKnots(nU, kU);
		patch.vKnots = MakeUniformKnots(nV, kV);
		
		int nBasisU = patch.GetNUBasis();
		int nBasisV = patch.GetNVBasis();
		patch.controlPoints.SetCount(nBasisU*nBasisV);
		
		for (int j = 0; j < nBasisV; j++) {
			double v = -1. + 2.*j/(nBasisV - 1);
			double y = R + v*(B/2. - R);
			for (int i = 0; i < nBasisU; i++) {
				double u = -1. + 2.*i/(nBasisU - 1);
				double x = u*L/2.;
				patch.controlPoints[j*nBasisU + i] = Point3D(x, y, -D + R);
			}
		}
		mgr.Add(pick(patch));
	}
	
	return mgr;
}

void SurfaceBSpline::DeployXSymmetry() {
	int sz = patches.size();
	patches.SetCount(2*sz);
	for (int i = 0; i < sz; ++i) {
		BSplinePatch &patch = patches[i + sz];
		patch = clone(patches[i]);
		for (Point3D &p : patch.controlPoints)
			p.x = -p.x;
		
    	int nU = patch.GetNUBasis();
    	int nV = patch.GetNVBasis();
    
   		for (int iV = 0; iV < nV; iV++) {
        	for (int iU = 0; iU < nU/2; iU++) {
           		int leftIdx  = iV*nU + iU;      
            	int rightIdx = iV*nU + (nU - 1 - iU);
            	Swap(patch.controlPoints[leftIdx], patch.controlPoints[rightIdx]);	
        	}
        }
	}
}

void SurfaceBSpline::DeployYSymmetry() {
	int sz = patches.size();
	patches.SetCount(2*sz);
	for (int i = 0; i < sz; ++i) {
		BSplinePatch &patch = patches[i + sz];
		patch = clone(patches[i]);
		for (Point3D &p : patch.controlPoints)
			p.y = -p.y;
		
		int nU = patch.GetNUBasis();
    	int nV = patch.GetNVBasis();
    
   		for (int iV = 0; iV < nV; iV++) {
        	for (int iU = 0; iU < nU/2; iU++) {
           		int leftIdx  = iV*nU + iU;      
            	int rightIdx = iV*nU + (nU - 1 - iU);
            	Swap(patch.controlPoints[leftIdx], patch.controlPoints[rightIdx]);	
        	}
        }
	}	
}

void SurfaceBSpline::CutX(bool leavePositive) {
	for (int i = patches.size()-1; i >= 0; --i) {
		BSplinePatch &patch = patches[i];
		bool arenegative = true, arepositive = true;
		for (Point3D &p : patch.controlPoints) {
			if (p.x > -EPS_LEN)
				arenegative = false;
			if (p.x <  EPS_LEN)
				arepositive = false;
		}
		if ((arenegative && leavePositive) || (arepositive && !leavePositive))	// Remove the patch if ALL control points comply
			patches.Remove(i);
	}
}

void SurfaceBSpline::CutY(bool leavePositive) {
	for (int i = patches.size()-1; i >= 0; --i) {
		BSplinePatch &patch = patches[i];
		bool arenegative = true, arepositive = true;
		for (Point3D &p : patch.controlPoints) {
			if (p.y > -EPS_LEN)
				arenegative = false;
			if (p.y <  EPS_LEN)
				arepositive = false;
		}
		if ((arenegative && leavePositive) || (arepositive && !leavePositive))	// Remove the patch if ALL control points comply
			patches.Remove(i);
	}
}

void SurfaceBSpline::CutZ(bool leavePositive) {
	for (int i = patches.size()-1; i >= 0; --i) {
		BSplinePatch &patch = patches[i];
		bool arenegative = true, arepositive = true;
		for (Point3D &p : patch.controlPoints) {
			if (p.z > -EPS_LEN)
				arenegative = false;
			if (p.z <  EPS_LEN)
				arepositive = false;
		}
		if ((arenegative && leavePositive) || (arepositive && !leavePositive))	// Remove the patch if ALL control points comply
			patches.Remove(i);
	}
}

void SurfaceBSpline::GetHull() {
	for (int i = patches.size()-1; i >= 0; --i) {
		BSplinePatch &patch = patches[i];
		bool ishull = false;
		for (Point3D &p : patch.controlPoints) {
			if (p.z < -EPS_LEN) {		// If just a point is underwater, it is hull
				ishull = true;
				break;
			}
		}
		if (!ishull)	
			patches.Remove(i);
	}	
}

void SurfaceBSpline::GetLid() {
	for (int i = patches.size()-1; i >= 0; --i) {
		BSplinePatch &patch = patches[i];
		bool islid = true;
		for (Point3D &p : patch.controlPoints) {
			if (p.z > EPS_LEN || p.z < -EPS_LEN) {
				islid = false;
				break;
			}
		}
		if (!islid)	
			patches.Remove(i);
	}	
}

	
void SurfaceBSpline::RoundClosest(double grid, double eps) {
	for (BSplinePatch &patch : patches) {
		for (Point3D &p : patch.controlPoints) {
			p.x = Upp::RoundClosest(p.x, grid, eps);
			p.y = Upp::RoundClosest(p.y, grid, eps);
			p.z = Upp::RoundClosest(p.z, grid, eps);
		}
	}
}

int SurfaceBSpline::FitToZ0(double zTolerance) {
	zTolerance = abs(zTolerance);
	for (BSplinePatch &patch : patches) {
		for (const Point3D &p : patch.controlPoints) {
			if (p.z > zTolerance)	// Just works if mesh is underwater
				return 0;
		}
	}
	int ret = 0;
	for (BSplinePatch &patch : patches) {
		for (Point3D &p : patch.controlPoints) {
			if (Between(p.z, -zTolerance, zTolerance)) {
				p.z = 0;
				ret++;
			}
		}
	}
	return ret;
}

String SurfaceBSpline::Heal(double grid, double eps) {
	String ret;
	
	RoundClosest(grid, eps);
	double zTolerance = -0.1;
	int num0 = FitToZ0(zTolerance);
	if (num0 > 0) 
		ret << "\n" << F(t_("Fitted to Z=0 %d points"), num0);
	
	return ret;
}

bool SurfaceBSpline::SaveGdf(const String& fileName, double g, bool symX, bool symY, bool iscsf) const {
	if (iscsf)
		ForceExt(fileName, ".csf");
	
	FileOut out(fileName);
	if (!out.IsOpen()) 
		return false;
	
	out << "BEMRosetta GDF mesh file export\n";
	if (!iscsf)
		out << F("%16<s ULEN GRAV\n", F("%d %12f", 1, g));
	else
		out << "ILOWHICSF=1\n";
	out << F("%16<s ISX ISY\n", F("%d %d", (symX ? 1 : 0), (symY ? 1 : 0)));
	out << F("%16<s NPATCH IGDEF\n", F("%d 1", size()));
	
	for (const BSplinePatch& p : patches) {
		out << p.nug << " " << p.nvg << "\n";
		out << p.kug << " " << p.kvg << "\n";
		
		for (int i = 0; i < p.uKnots.size(); i++)
			out << p.uKnots[i] << ((i % 4 == 3) ? "\n" : " ");
		if (p.uKnots.size() % 4 != 0) 
			out << "\n";
		
		for (int i = 0; i < p.vKnots.size(); i++)
			out << p.vKnots[i] << ((i % 4 == 3) ? "\n" : " ");
		if (p.vKnots.size() % 4 != 0) 
			out << "\n";
		
		for (const Point3D& cp : p.controlPoints)
			out << F("%10.5f %10.5f %10.5f", cp.x, cp.y, cp.z) << "\n";
	}	
	return true;
}

void SurfaceBSpline::Tessellate(int nu, int nv, Surface& surf) {
	for (const BSplinePatch& patch : patches) {
		Surface s;
		
		double umin = patch.GetUMin();
		double umax = patch.GetUMax();
		double vmin = patch.GetVMin();
		double vmax = patch.GetVMax();
		
		double du = (umax - umin)/nu;
		double dv = (vmax - vmin)/nv;
		
		int nU = patch.GetNUBasis();
		int nV = patch.GetNVBasis();
		int p  = patch.GetUDegree();
		int q  = patch.GetVDegree();
		
		s.nodes.SetCount((nu + 1)*(nv + 1));
		for (int j = 0; j <= nv; j++) {
			double v = (j == nv) ? vmax : vmin + j*dv;
			for (int i = 0; i <= nu; i++) {
				double u = (i == nu) ? umax : umin + i*du;
				s.nodes[j*(nu + 1) + i] = patch.Evaluate(u, v);
			}
		}
		
		s.panels.SetCount(nu*nv);
		int ipanel = 0;
		for (int j = 0; j < nv; j++) {
			for (int i = 0; i < nu; i++) {
				s.panels[ipanel].id[0] = j*(nu + 1) + i;
				s.panels[ipanel].id[1] = (j + 1)*(nu + 1) + i;
				s.panels[ipanel].id[2] = (j + 1)*(nu + 1) + i + 1;
				s.panels[ipanel].id[3] = j*(nu + 1) + i + 1;
				ipanel++;
			}
		}
		surf.Append(s);
	}
	surf.GetPanelParams();
	Surface::RemoveDuplicatedPointsAndRenumber(surf.panels, surf.nodes, surf.segments);
}