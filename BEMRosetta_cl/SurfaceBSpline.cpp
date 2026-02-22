#include <Core/Core.h>
#include <Surface/Surface.h>
#include "SurfaceBSpline.h"

using namespace Upp;



int FindSpan(double u, const UVector<double>& knots, int n, int p) {
	int m = knots.GetCount() - 1;
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
			if (idx1 < knots.GetCount()) {
				double denom = knots[idx1] - knots[idx2];
				if (denom > 1e-10)
					dval += p * N[(p-1)%2][i-1] / denom;
			}
		}
		if (i < p) {					// Second term: -p/(knots[span+i+p+1] - knots[span+i+1]) * N_{span+i+1-p, p-1}
			int idx1 = span + i + p + 1;
			int idx2 = span + i + 1;
			if (idx1 >= knots.GetCount()) 
				idx1 = knots.GetCount() - 1;
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
	double umax = uKnots[uKnots.GetCount() - p - 1];
	double vmin = vKnots[q];
	double vmax = vKnots[vKnots.GetCount() - q - 1];
	
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
			result = result + cpts[vidx * nU + uidx]*(Nu[k]*Nv[l]);
		}
	}
	return result;
}


bool SurfaceBSpline::ParseGDF(const String& filename, SurfaceBSpline& manager) {
	FileIn file(filename);
	if (!file.IsOpen()) {
		Cout() << "Cannot open: " << filename << "\n";
		return false;
	}
	
	auto readLine = [&]() -> String {
		String line = file.GetLine();
		int pos = min(line.Find('#'), line.Find('!'));
		if (pos >= 0) line = line.Left(pos);
		return TrimBoth(line);
	};
	
	auto read2int = [&](int& a, int& b) -> bool {
		String s = readLine();
		CParser p(s);
		try { a = p.ReadInt(); b = p.ReadInt(); return true; }
		catch (...) { return false; }
	};
	
	auto read2double = [&](double& a, double& b) -> bool {
		String s = readLine();
		CParser p(s);
		try { a = p.ReadDouble(); b = p.ReadDouble(); return true; }
		catch (...) { return false; }
	};
	
	String header = readLine();
	Cout() << "Header: " << header << "\n";
	
	double ulen, grav;
	int isx, isy, npatch, igdef;
	
	if (!read2double(ulen, grav)) return false;
	if (!read2int(isx, isy)) return false;
	if (!read2int(npatch, igdef)) return false;
	
	Cout() << "NPATCH=" << npatch << ", IGDEF=" << igdef << "\n";
	
	if (igdef != 1) {
		Cout() << "Only IGDEF=1 supported\n";
		return false;
	}
	
	for (int i = 0; i < npatch; i++) {
		BSplinePatch patch;
		patch.id = i;
		
		if (!read2int(patch.nug, patch.nvg)) return false;
		if (!read2int(patch.kug, patch.kvg)) return false;
		
		int nua = patch.nug + 2*patch.kug - 1;
		int nva = patch.nvg + 2*patch.kvg - 1;
		int nb = patch.GetNUBasis() * patch.GetNVBasis();
		
		patch.uKnots.SetCount(nua);
		int read = 0;
		while (read < nua) {
			String s = readLine();
			CParser p(s);
			while (!p.IsEof() && read < nua)
				patch.uKnots[read++] = p.ReadDouble();
		}
		
		patch.vKnots.SetCount(nva);
		read = 0;
		while (read < nva) {
			String s = readLine();
			CParser p(s);
			while (!p.IsEof() && read < nva)
				patch.vKnots[read++] = p.ReadDouble();
		}
		
		patch.controlPoints.SetCount(nb);
		for (int j = 0; j < nb; j++) {
			String s = readLine();
			CParser p(s);
			double x = p.ReadDouble();
			double y = p.ReadDouble();
			double z = p.ReadDouble();
			patch.controlPoints[j] = Point3D(x, y, z);
		}
		
		manager.Add(pick(patch));
	}
	
	return true;
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

	patch.id = 0;
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
	patch.id = 0;
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
	patch.id = 0;
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
	patch.id = 0;
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
		patch.id = 0;
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
		patch.id = 1;
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
		patch.id = 2;
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


bool SaveSurfaceBSpline(const SurfaceBSpline& mgr, const String& filename,
                         const String& header = "Generated GDF",
                         double ulen = 1.0, double grav = 9.80665,
                         int isx = 0, int isy = 0) {
	FileOut file(filename);
	if (!file.IsOpen()) return false;
	
	file << header << "\n";
	file << ulen << " " << grav << "\n";
	file << isx << " " << isy << "\n";
	file << mgr.GetCount() << " 1\n";  // IGDEF=1
	
	for (const BSplinePatch& p : mgr.patches) {
		file << p.nug << " " << p.nvg << "\n";
		file << p.kug << " " << p.kvg << "\n";
		
		// U knots
		for (int i = 0; i < p.uKnots.GetCount(); i++)
			file << p.uKnots[i] << ((i % 4 == 3) ? "\n" : " ");
		if (p.uKnots.GetCount() % 4 != 0) 
			file << "\n";
		
		// V knots
		for (int i = 0; i < p.vKnots.GetCount(); i++)
			file << p.vKnots[i] << ((i % 4 == 3) ? "\n" : " ");
		if (p.vKnots.GetCount() % 4 != 0) 
			file << "\n";
		
		// Control points
		for (const Point3D& cp : p.controlPoints)
			file << Format("%.6f %.6f %.6f", cp.x, cp.y, cp.z) << "\n";
	}
	
	file.Close();
	return true;
}


void SurfaceBSpline::Tessellate(const SurfaceBSpline& manager, int nu, int nv, UArray<FlatPanel>& panels) {
	for (const BSplinePatch& patch : manager.patches) {
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
		
		UVector<Point3D> grid;
		grid.SetCount((nu + 1)*(nv + 1));
		
		for (int j = 0; j <= nv; j++) {
			double v = (j == nv) ? vmax : vmin + j*dv;
			for (int i = 0; i <= nu; i++) {
				double u = (i == nu) ? umax : umin + i*du;
				grid[j*(nu + 1) + i] = patch.Evaluate(u, v);
			}
		}
		
		// Create panels
		for (int j = 0; j < nv; j++) {
			for (int i = 0; i < nu; i++) {
				FlatPanel panel;
				panel.v[0] = grid[j*(nu + 1) + i];
				panel.v[1] = grid[(j + 1)*(nu + 1) + i];
				panel.v[2] = grid[(j + 1)*(nu + 1) + i + 1];
				panel.v[3] = grid[j*(nu + 1) + i + 1];
				
				panels.Add(panel);
			}
		}
	}
}