#ifndef _KK_SurfaceBSpline_h_
#define _KK_SurfaceBSpline_h_

namespace Upp {
	
class BSplinePatch : Moveable<BSplinePatch> {
public:
	int nug, nvg;
	int kug, kvg;
	UVector<double> uKnots, vKnots;
	UVector<Point3D> controlPoints;

	BSplinePatch() {}
	BSplinePatch(const BSplinePatch &patch, int)		{Copy(patch);}
	BSplinePatch(const BSplinePatch &patch)				{Copy(patch);}
	BSplinePatch& operator=(const BSplinePatch &patch) 	{Copy(patch); return *this;};
	BSplinePatch& operator=(BSplinePatch&&) noexcept = default;
	void Copy(const BSplinePatch &patch) {
		nug = patch.nug;
		nvg = patch.nvg;
		kug = patch.kug;
		kvg = patch.kvg;
		uKnots = clone(patch.uKnots);
		vKnots = clone(patch.vKnots);
		controlPoints = clone(patch.controlPoints);
	}
	BSplinePatch(BSplinePatch &&patch) noexcept {
		nug = patch.nug;
		nvg = patch.nvg;
		kug = patch.kug;
		kvg = patch.kvg;
		uKnots = pick(patch.uKnots);
		vKnots = pick(patch.vKnots);
		controlPoints = pick(patch.controlPoints);		
	}
		
	int GetNUBasis() const {return nug + kug - 1;}
	int GetNVBasis() const {return nvg + kvg - 1;}
	double GetUMin() const {return uKnots[kug - 1];}
	double GetUMax() const {return uKnots[uKnots.GetCount() - kug];}
	double GetVMin() const {return vKnots[kvg - 1];}
	double GetVMax() const {return vKnots[vKnots.GetCount() - kvg];}
	int GetUDegree() const {return kug - 1;}
	int GetVDegree() const {return kvg - 1;}
	
	Point3D GetControlPoint(int i, int j) const {return controlPoints[j * GetNUBasis() + i];}
	void SetControlPoint(int i, int j, const Point3D& pt) {controlPoints[j * GetNUBasis() + i] = pt;}
	
	static Point3D EvaluateSurface(const UVector<Point3D>& cpts, int nU, int nV,
                     const UVector<double>& uKnots, const UVector<double>& vKnots,
                     int p, int q, double u, double v);
                         
	Point3D Evaluate(double u, double v) const {
		return EvaluateSurface(controlPoints, GetNUBasis(), GetNVBasis(),
			uKnots, vKnots, GetUDegree(), GetVDegree(), u, v);
	}
	
	void EvaluateDerivatives(double u, double v, Point3D& Su, Point3D& Sv) const;
	double EvaluateAreaElement(double u, double v) const;
	double ComputeArea(int order = 3) const;
	
	void Scale(const Value3D& scale) {
		for (auto& cp : controlPoints)
			cp = cp.dot(scale);
	}
	void Scale(const Value3D& scale, const Point3D &c0) {
		for (auto& cp : controlPoints)
			cp.Translate((cp - c0)*scale); 
	}
	
	void Translate(const Point3D& offset) {
		for (auto& cp : controlPoints)
			cp.Translate(offset);
	}
	void Rotate(const Value3D& angle, const Point3D& c0, RotationOrder order = RotationOrder::XYZ) {
		for (auto& cp : controlPoints)
			cp.Rotate(angle, c0, order);
	}
	void TransRot(const Point3D& offset, const Value3D& angle, const Point3D& c0, RotationOrder order = RotationOrder::XYZ) {
		for (auto& cp : controlPoints)
			cp.TransRot(offset, angle, c0, order);
	}
		
	void GetBoundingBox(Point3D& mn, Point3D& mx) const;
	
	void Jsonize(JsonIO &json) {
		json
			("nug", nug)
			("nvg", nvg)
			("kug", kug)
			("kvg", kvg)
			("uKnots", uKnots)
			("vKnots", vKnots)
			("controlPoints", controlPoints)
		;
	}
};

class SurfaceBSpline : Moveable<SurfaceBSpline> {
public:
	UVector<BSplinePatch> patches;

	SurfaceBSpline() {}
	SurfaceBSpline(const SurfaceBSpline &patch, int)		{Copy(patch);}
	SurfaceBSpline(const SurfaceBSpline &patch)				{Copy(patch);}
	SurfaceBSpline& operator=(const SurfaceBSpline &patch) 	{Copy(patch); return *this;};
	SurfaceBSpline& operator=(SurfaceBSpline&&) noexcept = default;
	void Copy(const SurfaceBSpline &surf) 					{patches = clone(surf.patches);}
	SurfaceBSpline(SurfaceBSpline &&surf) noexcept 			{patches = pick(surf.patches);}
	
	void Append(const SurfaceBSpline &surf) 	{patches.Append(surf.patches);}
		
	void Add(BSplinePatch&& p) 					{patches.Add(pick(p));}
	int size() const 							{return patches.size();}
	bool IsEmpty() const						{return patches.IsEmpty();}
	BSplinePatch& operator[](int i) 			{return patches[i];}
	const BSplinePatch& operator[](int i) const {return patches[i];}
	
	void DeployXSymmetry();
	void DeployYSymmetry();
	void CutX(bool leavePositive = true);
	void CutY(bool leavePositive = true);
	void CutZ(bool leavePositive = true);

	String Heal(double grid, double eps);
	
	void Scale(const Value3D& scale) {
		for (auto& p : patches) 
			p.Scale(scale);
	}
	void Scale(const Value3D& scale, const Point3D &c0) {
		for (auto& p : patches) 
			p.Scale(scale, c0);
	}
		
	void Translate(const Value3D& offset) {
		for (auto& p : patches) 
			p.Translate(offset);
	}
	void Rotate(const Value3D& angle, const Point3D& centre) {
		for (auto& p : patches) 
			p.Rotate(angle, centre);
	}
	void TransRot(const Value3D& offset, const Value3D& angle, const Point3D& centre) {
		for (auto& p : patches) 
			p.TransRot(offset, angle, centre);
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
	void Tessellate(int nu, int nv, Surface& surf);
	bool SaveGdf(const String& filename, double g, bool symx, bool symy, bool iscsf) const;
	
	void Jsonize(JsonIO &json) {
		json
			("patches", patches)
		;
	}
	
private:
	void RoundClosest(double grid, double eps);
	int FitToZ0(double zTolerance);
};


}

#endif
