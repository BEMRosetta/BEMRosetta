#ifndef _KK_SurfaceBSpline_h_
#define _KK_SurfaceBSpline_h_

namespace Upp {
	
class BSplinePatch : Moveable<BSplinePatch> {
public:
	int id;
	int nug, nvg;
	int kug, kvg;
	UVector<double> uKnots, vKnots;
	UVector<Point3D> controlPoints;

	BSplinePatch() {}
	BSplinePatch(const BSplinePatch &patch, int)		{Copy(patch);}
	BSplinePatch(const BSplinePatch &patch)			{Copy(patch);}
	BSplinePatch& operator=(const BSplinePatch &patch) {Copy(patch); return *this;};
	BSplinePatch& operator=(BSplinePatch&&) noexcept = default;
	void Copy(const BSplinePatch &patch) {
		id = patch.id;
		nug = patch.nug;
		nvg = patch.nvg;
		kug = patch.kug;
		kvg = patch.kvg;
		uKnots = clone(patch.uKnots);
		vKnots = clone(patch.vKnots);
		controlPoints = clone(patch.controlPoints);
	}
		
	int GetNUBasis() const { return nug + kug - 1; }
	int GetNVBasis() const { return nvg + kvg - 1; }
	double GetUMin() const { return uKnots[kug - 1]; }
	double GetUMax() const { return uKnots[uKnots.GetCount() - kug]; }
	double GetVMin() const { return vKnots[kvg - 1]; }
	double GetVMax() const { return vKnots[vKnots.GetCount() - kvg]; }
	int GetUDegree() const { return kug - 1; }
	int GetVDegree() const { return kvg - 1; }
	
	Point3D GetControlPoint(int i, int j) const {return controlPoints[j * GetNUBasis() + i];}
	void SetControlPoint(int i, int j, const Point3D& pt) {controlPoints[j * GetNUBasis() + i] = pt;}
	
	static Point3D EvaluateSurface(const UVector<Point3D>& cpts, int nU, int nV,
                     const UVector<double>& uKnots, const UVector<double>& vKnots,
                     int p, int q, double u, double v);
                         
	Point3D Evaluate(double u, double v) const {
		return EvaluateSurface(controlPoints, GetNUBasis(), GetNVBasis(),
			uKnots, vKnots, GetUDegree(), GetVDegree(), u, v);
	}
	
	void EvaluateDerivatives(double u, double v, Point3D& Su, Point3D& Sv) const {
		double du = 0.001*(GetUMax() - GetUMin());
		double dv = 0.001*(GetVMax() - GetVMin());
		
		Point3D pu  = Evaluate(u - du, v);
		Point3D puu = Evaluate(u + du, v);
		Su = (puu - pu)/(2*du);
		
		Point3D pv  = Evaluate(u, v - dv);
		Point3D pvv = Evaluate(u, v + dv);
		Sv = (pvv - pv)/(2*dv);
	}
	
	double EvaluateAreaElement(double u, double v) const {
		Point3D Su, Sv;
		EvaluateDerivatives(u, v, Su, Sv);
		return cross(Su, Sv).Length();
	}

	double ComputeArea(int order = 3) const {
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
	
	void Scale(double sx, double sy, double sz) {
		for (auto& cp : controlPoints) {
			cp.x *= sx;
			cp.y *= sy;
			cp.z *= sz;
		}
	}
	void Scale(double s) { Scale(s, s, s); }
	
	void Translate(const Point3D& offset) {
		for (auto& cp : controlPoints) {
			cp = cp + offset;
		}
	}
	
	void RotateX(double angle) {
		double c = cos(angle), s = sin(angle);
		for (auto& cp : controlPoints) {
			double y = cp.y*c - cp.z*s;
			double z = cp.y*s + cp.z*c;
			cp.y = y; 
			cp.z = z;
		}
	}
	void RotateY(double angle) {
		double c = cos(angle), s = sin(angle);
		for (auto& cp : controlPoints) {
			double x =  cp.x*c + cp.z*s;
			double z = -cp.x*s + cp.z*c;
			cp.x = x; 
			cp.z = z;
		}
	}
	void RotateZ(double angle) {
		double c = cos(angle), s = sin(angle);
		for (auto& cp : controlPoints) {
			double x = cp.x*c - cp.y*s;
			double y = cp.x*s + cp.y*c;
			cp.x = x; 
			cp.y = y;
		}
	}
	
	void GetBoundingBox(Point3D& mn, Point3D& mx) const {
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
};

class SurfaceBSpline {
public:
	UVector<BSplinePatch> patches;
	
	void Add(BSplinePatch&& p) 					{patches.Add(pick(p));}
	int GetCount() const 						{return patches.GetCount();}
	BSplinePatch& operator[](int i) 			{return patches[i];}
	const BSplinePatch& operator[](int i) const {return patches[i];}
	
	// Batch transformations
	void Scale(double sx, double sy, double sz) {
		for (auto& p : patches) 
			p.Scale(sx, sy, sz);
	}
	void Scale(double s) { Scale(s, s, s); }
	void Translate(const Point3D& offset) {
		for (auto& p : patches) 
			p.Translate(offset);
	}
	void RotateX(double angle) { 
		for (auto& p : patches) 
			p.RotateX(angle); 
	}
	void RotateY(double angle) { 
		for (auto& p : patches) 
			p.RotateY(angle); 
	}
	void RotateZ(double angle) { 
		for (auto& p : patches) 
			p.RotateZ(angle); 
	}
	
	double ComputeTotalArea(int order = 3) const {
		double total = 0;
		for (const auto& p : patches) 
			total += p.ComputeArea(order);
		return total;
	}

	void GetBoundingBox(Point3D& mn, Point3D& mx) const {
		if (patches.IsEmpty()) { 
			mn = mx = Point3D(); 
			return; 
		}
		patches[0].GetBoundingBox(mn, mx);
		for (int i = 1; i < patches.GetCount(); i++) {
			Point3D pmin, pmax;
			patches[i].GetBoundingBox(pmin, pmax);
			mn.x = min(mn.x, pmin.x); mx.x = max(mx.x, pmax.x);
			mn.y = min(mn.y, pmin.y); mx.y = max(mx.y, pmax.y);
			mn.z = min(mn.z, pmin.z); mx.z = max(mx.z, pmax.z);
		}
	}

	static bool ParseGDF(const String& filename, SurfaceBSpline& manager);
	struct FlatPanel {
		Point3D v[4];
	};
	static void Tessellate(const SurfaceBSpline& manager, int nu, int nv, UArray<FlatPanel>& panels);
};





}

#endif
