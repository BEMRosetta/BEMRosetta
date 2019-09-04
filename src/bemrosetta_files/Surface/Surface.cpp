#include <Core/Core.h>

using namespace Upp;

#include <Surface/Surface.h>
#include <plugin/Eigen/Eigen.h>

using namespace Eigen;

Point3D GetCentroid(const Point3D &a, const Point3D &b) {
	return Point3D(avg(a.x, b.x), avg(a.y, b.y), avg(a.z, b.z));	
}

Point3D GetCentroid(const Point3D &a, const Point3D &b, const Point3D &c) {
	return Point3D(avg(a.x, b.x,c.x), avg(a.y, b.y, c.y), avg(a.z, b.z, c.z));	
}

Vector3D GetNormal(const Point3D &a, const Point3D &b, const Point3D &c) {
	return Vector3D((a - b) % (b - c)).Normalize();
}

Point3D Intersection(const Vector3D &lineVector, const Point3D &linePoint, const Point3D &planePoint, const Vector3D &planeNormal) {
	Vector3D diff = planePoint - linePoint;
	double prod1 = diff.dot(planeNormal);
	double prod2 = lineVector.dot(planeNormal);
	if (abs(prod2) < EPS)
		return Null;
	double factor = prod1/prod2;
	return linePoint + lineVector*factor;	
}

void Point3D::MoveTo(const Point3D &pos, double _x, double _y, double _z, const Affine3d &quat) {
	x = pos.x + _x;
	y = pos.y + _y;
	z = pos.z + _z;

	Vector3d pnt0(x, y, z);	
	Vector3d pnt = quat * pnt0;

	x = pnt[0];
	y = pnt[1];
	z = pnt[2];
}

void Point3D::MoveTo(const Point3D &pos, double _x, double _y, double _z, double a_x, double a_y, double a_z, double c_x, double c_y, double c_z) {
	Affine3d aff;
	GetTransform(aff, a_x, a_y, a_z, c_x, c_y, c_z);
	MoveTo(pos, _x, _y, _z, aff);
}

void GetTransform(Affine3d &aff, double a_x, double a_y, double a_z, double c_x, double c_y, double c_z) {
	Vector3d c(c_x, c_y, c_z);	
	aff =	Translation3d(c) *
			AngleAxisd(a_x*M_PI/180, Vector3d::UnitX()) *
		 	AngleAxisd(a_y*M_PI/180, Vector3d::UnitY()) *
		 	AngleAxisd(a_z*M_PI/180, Vector3d::UnitZ()) *
		 	Translation3d(-c);
}

Point3D Segment3D::IntersectionPlaneX(double x) {
	if (from.x >= x && to.x >= x)
		return Point3D(true);
	if (from.x <= x && to.x <= x)
		return Point3D(false);
	
	double factor = (x - from.x)/(to.x - from.x);
	return Point3D(x, from.y + (to.y - from.y)*factor, from.z + (to.z - from.z)*factor);
}

Point3D Segment3D::IntersectionPlaneY(double y) {
	if (from.y >= y && to.y >= y)
		return Point3D(true);
	if (from.y <= y && to.y <= y)
		return Point3D(false);
	
	double factor = (y - from.y)/(to.y - from.y);
	return Point3D(from.x + (to.x - from.x)*factor, y, from.z + (to.z - from.z)*factor);
}

Point3D Segment3D::IntersectionPlaneZ(double z) {
	if (from.z >= z && to.z >= z)
		return Point3D(true);
	if (from.z <= z && to.z <= z)
		return Point3D(false);
	
	double factor = (z - from.z)/(to.z - from.z);
	return Point3D(from.x + (to.x - from.x)*factor, from.y + (to.y - from.y)*factor, z);
}

Point3D Segment3D::Intersection(const Point3D &planePoint, const Vector3D &planeNormal) {
	Vector3D vector = Vector();
	Vector3D diff = planePoint - from;
	double prod1 = diff.dot(planeNormal);
	double prod2 = vector.dot(planeNormal);
	if (abs(prod2) < EPS)
		return Null;
	double factor = prod1/prod2;
	if (factor >= 1)
		return Point3D(true);
	if (factor <= 0)
		return Point3D(false);
	return from + vector*factor;	
}

void Surface::Clear() {
	nodes0.Clear();
	nodes.Clear();
	panels.Clear();
	skewed.Clear();
	segWaterlevel.Clear();
	segTo1panel.Clear();
	segTo3panel.Clear();
	segments.Clear();
	x = y = z = a_x = a_y = a_z = 0;
}

Surface::Surface(const Surface &orig, int) {
	healing = orig.healing;
	numTriangles = orig.numTriangles;
	numBiQuads = orig.numBiQuads;
	numMonoQuads = orig.numMonoQuads;
	
	panels = clone(orig.panels);
	nodes0 = clone(orig.nodes0);
	skewed = clone(orig.skewed);
	segWaterlevel = clone(orig.segWaterlevel);
	segTo1panel = clone(orig.segTo1panel);
	segTo3panel = clone(orig.segTo3panel);
	segments = clone(orig.segments);
	
	env = clone(orig.env);
	
	surface = orig.surface;
	volume = orig.volume;
}

bool Surface::IsEmpty() {
	return nodes0.IsEmpty();
}

int Surface::FixSkewed() {	
	int num = 0;
	for (int i = 0; i < panels.GetCount(); ++i) {
		int &id0 = panels[i].id[0];
		int &id1 = panels[i].id[1];
		int &id2 = panels[i].id[2];
		int &id3 = panels[i].id[3];
		Point3D &p0 = nodes0[id0];
		Point3D &p1 = nodes0[id1];
		Point3D &p2 = nodes0[id2];
		Point3D &p3 = nodes0[id3];
		if (id0 != id3) {		// Is not triangular 
			Vector3D normal301 = GetNormal(p3, p0, p1);
			Vector3D normal012 = GetNormal(p0, p1, p2);
			Vector3D normal123 = GetNormal(p1, p2, p3);
			Vector3D normal230 = GetNormal(p2, p3, p0);
			double d0  = normal301.Manhattan();
			double d01 = normal301.Manhattan(normal012);
			double d02 = normal301.Manhattan(normal123);
			double d03 = normal301.Manhattan(normal230);
			
			int numg = 0;
			if (d0 < d01)
				numg++;
			if (d0 < d02)
				numg++;
			if (d0 < d03)
				numg++;	 
			if (numg > 1) {
				skewed << Segment3D(p0, p1) << Segment3D(p1, p2) << Segment3D(p2, p3) << Segment3D(p3, p0);
				
				if (d0 < d01)
					Swap(panels[i].id[1], panels[i].id[2]);
				else if (d0 < d02)
					Swap(panels[i].id[2], panels[i].id[3]);
				num++;
			}
		}
	}	
	return num;
}

void Surface::DetectTriBiP(int &numTri, int &numBi, int &numP) {
	numTri = numBi = numP = 0;
	for (int i = panels.GetCount()-1; i >= 0; --i) {
		Panel &panel = panels[i];
		Upp::Index<int> ids;
		ids.FindAdd(panel.id[0]);
		ids.FindAdd(panel.id[1]);
		ids.FindAdd(panel.id[2]);
		ids.FindAdd(panel.id[3]);
		if (ids.GetCount() == 4)
			;
		else if (ids.GetCount() == 3) {
			numTri++;
			panel.id[0] = ids[0];
			panel.id[1] = ids[1];
			panel.id[2] = ids[2];
			panel.id[3] = ids[0];
		} else if (ids.GetCount() == 2) {
			numBi++;
			panels.Remove(i, 1);
		} else {
			numP++;
			panels.Remove(i, 1);
		}
	}
}

int Surface::RemoveDuplicatedPanels() {		
	int num = 0;
	for (int i = 0; i < panels.GetCount()-1; ++i) {
		Panel &panel = panels[i];
		for (int j = panels.GetCount()-1; j >= i+1; --j) {
			if (panel == panels[j]) {
				num++;
				panels.Remove(j, 1);
			}
		}
	}
	return num;
}

int Surface::RemoveDuplicatedPointsAndRenumber() {
	int num = 0;
	
	// Detect duplicate points in nodes
	double similThres = 0.00001;
	Upp::Index<int> similars, goods;
	for (int i = 0; i < nodes0.GetCount()-1; ++i) {
		if (similars.Find(i) >= 0)
			continue;
		for (int j = i+1; j < nodes0.GetCount(); ++j) {
			if (nodes0[i].IsSimilar(nodes0[j], similThres)) {
				similars.Add(j);
				goods.Add(i);
				num++;
			}
		}
	}
	// Replace duplicated points with good ones in panels
	for (int i = 0; i < panels.GetCount(); ++i) {
		int numP = PanelGetNumNodes(i);
		for (int j = 0; j < numP; ++j) {
			int &id = panels[i].id[j];
			int pos = similars.Find(id);
			if (pos >= 0)
				id = goods[pos];
		}
	}
	
	// Find duplicated and unused nodes
	Vector<int> newId;
	newId.SetCount(nodes0.GetCount());
	int avId = 0;
	for (int i = 0; i < nodes0.GetCount(); ++i) {
		bool found = false;
		for (int ip = 0; ip < panels.GetCount() && !found; ++ip) {
			int numP = PanelGetNumNodes(ip);
			for (int j = 0; j < numP; ++j) {
				if (panels[ip].id[j] == i) {
					found = true;
					break;
				}
			}
		}
		if (!found)
			newId[i] = Null;	// Remove unused nodes
		else if (similars.Find(i) < 0) {
			newId[i] = avId;
			avId++;
		} else
			newId[i] = Null;
	}
	// Remove duplicated nodes
	for (int i = nodes0.GetCount()-1; i >= 0; --i) {
		if (IsNull(newId[i]))
			nodes0.Remove(i, 1);
	}
	// Renumber panels
	for (int i = 0; i < panels.GetCount(); ++i) {
		int numP = PanelGetNumNodes(i);
		for (int j = 0; j < numP; ++j) {
			int& id = panels[i].id[j];
			id = newId[id];
		}
	}
	return num;
}
	
void Surface::AddSegment(int inode0, int inode1, int ipanel) {
	for (int i = 0; i < segments.GetCount(); ++i) {
		if ((segments[i].inode0 == inode0 && segments[i].inode1 == inode1) ||
			(segments[i].inode1 == inode0 && segments[i].inode0 == inode1)) {
			segments[i].panels << ipanel;
			return;
		}
	}
	Segment &sg = segments.Add();
	sg.inode0 = inode0;
	sg.inode1 = inode1;
	sg.panels << ipanel;
}

void Surface::AnalyseSegments(double zTolerance) {
	segments.Clear();
	
	for (int i = 0; i < panels.GetCount(); ++i) {
		AddSegment(panels[i].id[0], panels[i].id[1], i);
		AddSegment(panels[i].id[1], panels[i].id[2], i);
		if (IsPanelTriangle(i)) 
			AddSegment(panels[i].id[2], panels[i].id[0], i);
		else {
			AddSegment(panels[i].id[2], panels[i].id[3], i);
			AddSegment(panels[i].id[3], panels[i].id[0], i);
		}
	}
	
	for (int i = 0; i < segments.GetCount(); ++i) {
		int inode0 = segments[i].inode0;
		int inode1 = segments[i].inode1;
		
		if (inode0 >= nodes0.GetCount())
			throw Exc(Format(t_("Node %d is pointing out of scope"), inode0+1));	
		if (inode1 >= nodes0.GetCount())
			throw Exc(Format(t_("Node %d is pointing out of scope"), inode1+1));
		
		int num = segments[i].panels.GetCount();
				
		if (num == 1) {
			if (nodes0[segments[i].inode0].z >= zTolerance && nodes0[segments[i].inode1].z >= zTolerance)
				segWaterlevel << Segment3D(nodes0[segments[i].inode0], nodes0[segments[i].inode1]);
			else
				segTo1panel << Segment3D(nodes0[segments[i].inode0], nodes0[segments[i].inode1]);
		}
		if (num > 2)
			segTo3panel << Segment3D(nodes0[segments[i].inode0], nodes0[segments[i].inode1]);
	}
}

bool Surface::ReorientPanels() {
	numUnprocessed = -1;
	
	// Get the lowest panel with normal non horizontal
	int iLowSeg = Null;
	int iLowPanel;
	double zLowSeg = DBL_MAX;
	for (int i = 0; i < segments.GetCount(); ++i) {
		const Segment &seg = segments[i];
		if (seg.panels.GetCount() == 2) {
			for (int ip = 0; ip < seg.panels.GetCount(); ++ip) {
				if (panels[seg.panels[ip]].normal0.z != 0) {
					double zz = max(nodes0[seg.inode0].z, nodes0[seg.inode1].z);
					if (zz < zLowSeg) {
						zLowSeg = zz;
						iLowSeg = i;
						iLowPanel = ip;
					}
				}
			}
		}
	}
	if (IsNull(iLowSeg))
		return false;
	
	// Reorient lowest panel downwards to be the seed
	int ip = segments[iLowSeg].panels[iLowPanel];
	if (panels[ip].normal0.z > 0)
		ReorientPanel(ip);
	
	Vector<int> panelStack;
	Upp::Index<int> panelProcessed;
	
	panelStack << ip;
	while (!panelStack.IsEmpty()) {
		int id = panelStack.GetCount() - 1;
		int ipp = panelStack[id];
		panelStack.Remove(id, 1);
		panelProcessed << ipp;
		
		for (int is = 0; is < segments.GetCount(); ++is) {
			const Upp::Index<int> &segPanels = segments[is].panels;
			if (segPanels.Find(ipp) >= 0) {
				for (int i = 0; i < segPanels.GetCount(); ++i) {
					int ipadyac = segPanels[i];
					if (ipadyac != ipp && panelProcessed.Find(ipadyac) < 0) {
						panelStack << ipadyac;
						if (!SameOrderPanel(ipp, ipadyac, segments[is].inode0, segments[is].inode1))
							ReorientPanel(ipadyac);
					}
				}
			}
		}
	}
	
	numUnprocessed = panels.GetCount() - panelProcessed.GetCount();

	return true;
}

void Surface::ReorientPanel(int ip) {
	panels[ip].Swap();
	panels[ip].normal0.Mirror();
	if (panels[ip].IsTriangle()) 
		panels[ip].normal1.Mirror();
}

bool Panel::FirstNodeIs0(int in0, int in1) const {
	if (IsTriangle()) {
		if ((id[0] == in0 && id[1] == in1) ||
			(id[1] == in0 && id[2] == in1) ||
			(id[2] == in0 && id[0] == in1))
			return true;
		else
			return false;
	} else {
		if ((id[0] == in0 && id[1] == in1) ||
			(id[1] == in0 && id[2] == in1) ||
			(id[2] == in0 && id[3] == in1) ||
			(id[3] == in0 && id[0] == in1))
			return true;
		else
			return false;
	}
}

void Panel::RedirectTriangles() {
	int shift = 0;
	if (id[0] == id[1])
		shift = -1;
	else if (id[1] == id[2])
		shift = -2;
	else if (id[2] == id[3])
		shift = 1;
	else
		return;
	ShiftNodes(shift);
}

void Panel::ShiftNodes(int shift) {
	if (shift == 1) {
		id[1] = id[0];
		id[2] = id[1];
		id[3] = id[2];
		id[0] = id[3];
	} else if (shift == -1) { 
		id[0] = id[1];
		id[1] = id[2];
		id[2] = id[3];
		id[3] = id[0];
	} else if (shift == -2) { 
		id[0] = id[2];
		id[1] = id[3];
		id[2] = id[0];
		id[3] = id[1];
	} else
		throw t_("ShiftNodes value not implemented");
}

double Panel::GetSurface(const Point3D &p0, const Point3D &p1, const Point3D &p2) {
	double l01 = p0.Distance(p1);
	double l12 = p1.Distance(p2);
	double l02 = p0.Distance(p2);

	double s = (l01 + l12 + l02)/2;
	return sqrt(max(s*(s - l01)*(s - l12)*(s - l02), 0.)); 
}

bool Surface::SameOrderPanel(int ip0, int ip1, int in0, int in1) {
	bool first0in0 = panels[ip0].FirstNodeIs0(in0, in1);
	bool first1in0 = panels[ip1].FirstNodeIs0(in0, in1);
	
	return first0in0 != first1in0;
}

String Surface::Heal(Function <void(String, int pos)> Status) {
	String ret;
	
	Status(t_("Detecting triangles and wrong panels"), 50);
	DetectTriBiP(numTriangles, numBiQuads, numMonoQuads);
	if (numTriangles > 0)
		ret << "\n" << Format(t_("Found %d triangles"), numTriangles);
	if (numBiQuads > 0)
		ret << "\n" << Format(t_("Removed %d 2 points quads"), numBiQuads);
	if (numMonoQuads > 0)
		ret << "\n" << Format(t_("Removed %d 1 points quads"), numMonoQuads);
	
	Status(t_("Removing duplicated panels (pass 1)"), 55);
	numDupPan = RemoveDuplicatedPanels();
	
	Status(t_("Fixing skewed panels"), 60);
	numSkewed = FixSkewed();
	if (numSkewed > 0) 
		ret << "\n" << Format(t_("Fixed %d skewed panels"), numSkewed);

	Status(t_("Removing duplicated points"), 65);
	numDupP = RemoveDuplicatedPointsAndRenumber();
	if (numDupP > 0) 
		ret << "\n" << Format(t_("Removed %d duplicated points"), numDupP);	

	Status(t_("Removing duplicated panels (pass 2)"), 70);
	numDupPan += RemoveDuplicatedPanels();	// Second time after duplicated points
	if (numDupPan > 0) 
		ret << "\n" << Format(t_("Removed %d duplicated panels"), numDupPan);

	Status(t_("Analysing water tightness"), 75);
	double zTolerance = -0.1;
	AnalyseSegments(zTolerance);
	ret << "\n" << Format(t_("%d segments, %d water level, %d water leak and %d multipanel"), 
								segments.GetCount(), segWaterlevel.GetCount(), 
								segTo1panel.GetCount(), segTo3panel.GetCount());
	
	Status(t_("Reorienting panels water side"), 80);
	if (!ReorientPanels())
		ret << "\n" << t_("Failed to reorient panels to water side");
	else if (numUnprocessed > 0)
		ret << "\n" << Format(t_("%d panels not reoriented. Body contains separated surfaces"), numUnprocessed);
	
	healing = true;
	
	return ret;
}

void Surface::GetLimits() {
	env.maxX = env.maxY = env.maxZ = -DBL_MAX; 
	env.minX = env.minY = env.minZ = DBL_MAX;
	for (int i = 0; i < nodes.GetCount(); ++i) {
		env.maxX = max(env.maxX, nodes[i].x);
		env.minX = min(env.minX, nodes[i].x);
		env.maxY = max(env.maxY, nodes[i].y);
		env.minY = min(env.minY, nodes[i].y);
		env.maxZ = max(env.maxZ, nodes[i].z);
		env.minZ = min(env.minZ, nodes[i].z);
	}
}

void Surface::GetPanelParams(Panel &panel) {
	panel.RedirectTriangles();
	
	const Point3D &p0 = nodes[panel.id[0]];
	const Point3D &p1 = nodes[panel.id[1]];
	const Point3D &p2 = nodes[panel.id[2]];
	const Point3D &p3 = nodes[panel.id[3]];
	
	panel.surface0 = panel.GetSurface(p0, p1, p2);
	panel.centroid0 = GetCentroid(p0, p1, p2);
	panel.normal0 = GetNormal(p0, p1, p2);
	if (!panel.IsTriangle()) {
		panel.surface1 = panel.GetSurface(p2, p3, p0);
		panel.centroid1 = GetCentroid(p2, p3, p0);
		panel.normal1 = GetNormal(p2, p3, p0);
		double surf = panel.surface0 + panel.surface1;
		panel.centroidPaint.x = (panel.centroid0.x*panel.surface0 + panel.centroid1.x*panel.surface1)/surf;
		panel.centroidPaint.y = (panel.centroid0.y*panel.surface0 + panel.centroid1.y*panel.surface1)/surf;
		panel.centroidPaint.z = (panel.centroid0.z*panel.surface0 + panel.centroid1.z*panel.surface1)/surf;
		panel.normalPaint.x = (panel.normal0.x*panel.surface0 + panel.normal1.x*panel.surface1)/surf;
		panel.normalPaint.y = (panel.normal0.y*panel.surface0 + panel.normal1.y*panel.surface1)/surf;
		panel.normalPaint.z = (panel.normal0.z*panel.surface0 + panel.normal1.z*panel.surface1)/surf;
		panel.normalPaint.Normalize();
	} else {
		panel.surface1 = 0;
		panel.centroidPaint = panel.centroid1 = panel.centroid0;
		panel.normalPaint = panel.normal1 = panel.normal0;
	}
}

void Surface::GetPanelParams() {
	for (int ip = 0; ip < panels.GetCount(); ++ip) {
		Panel &panel = panels[ip];
		GetPanelParams(panel);
	}	
}

void Surface::GetSurface() {
	surface = 0;
	for (int ip = 0; ip < panels.GetCount(); ++ip) 
		surface += panels[ip].surface0 + panels[ip].surface1;
	avgFacetSideLen  = sqrt(surface/panels.GetCount());
}

double Surface::GetWaterPlaneArea() {
	double area = 0;
	
	for (int ip = 0; ip < panels.GetCount(); ++ip) {
		Panel &panel = panels[ip];
		
		area += -panel.surface0*panel.normal0.z;
		
		if (!panel.IsTriangle()) 
			area += -panel.surface1*panel.normal1.z;
	}
	return area;
}

void Surface::GetVolume() {
	double volx = 0, voly = 0, volz = 0;
	
	for (int ip = 0; ip < panels.GetCount(); ++ip) {
		Panel &panel = panels[ip];
		
		volx += panel.surface0*panel.normal0.x*panel.centroid0.x;
		voly += panel.surface0*panel.normal0.y*panel.centroid0.y;
		volz += panel.surface0*panel.normal0.z*panel.centroid0.z;
		
		if (!panel.IsTriangle()) {
			volx += panel.surface1*panel.normal1.x*panel.centroid1.x;
			voly += panel.surface1*panel.normal1.y*panel.centroid1.y;
			volz += panel.surface1*panel.normal1.z*panel.centroid1.z;
		}
	}
	volume = avg(volx, voly, volz);
}
	
Point3D Surface::GetCenterOfBuoyancy() {
	double xb = 0, yb = 0, zb = 0;
	
	for (int ip = 0; ip < panels.GetCount(); ++ip) {
		Panel &panel = panels[ip];
		
		xb += panel.surface0*panel.normal0.x*sqr(panel.centroid0.x);
		yb += panel.surface0*panel.normal0.y*sqr(panel.centroid0.y);
		zb += panel.surface0*panel.normal0.z*sqr(panel.centroid0.z);
		
		if (!panel.IsTriangle()) {
			xb += panel.surface1*panel.normal1.x*sqr(panel.centroid1.x);
			yb += panel.surface1*panel.normal1.y*sqr(panel.centroid1.y);
			zb += panel.surface1*panel.normal1.z*sqr(panel.centroid1.z);
		}
	}
	
	xb /= 2*volume;
	yb /= 2*volume;
	zb /= 2*volume;
	
	return Point3D(xb, yb, zb);
}

void Surface::GetHydrostaticStiffness(MatrixXd &c, const Point3D &cb, double rho, 
					const Point3D &cg, double mass, double g, double zTolerance) {
	c.setConstant(6, 6, 0);
	
	if (volume < 0.00001)
		return;
	
	if (IsNull(mass))
		mass = rho*volume;
		
	for (int ip = 0; ip < panels.GetCount(); ++ip) {
		Panel &panel = panels[ip];
		const Point3D &p0 = nodes[panel.id[0]];
		const Point3D &p1 = nodes[panel.id[1]];
		const Point3D &p2 = nodes[panel.id[2]];
		const Point3D &p3 = nodes[panel.id[3]];
		
		if (p0.z <= zTolerance && p1.z <= zTolerance && p2.z <= zTolerance && 
			(panel.IsTriangle() || p3.z <= zTolerance)) {
			double momentz0 = panel.normal0.z*panel.surface0;
			double momentz1 = panel.normal1.z*panel.surface1;
			double x0 = panel.centroid0.x - cg.x;
			double y0 = panel.centroid0.y - cg.y;
			double x1 = panel.centroid1.x - cg.x;
			double y1 = panel.centroid1.y - cg.y;
			c(2, 2) -= (momentz0 + momentz1);
            c(2, 3) -= (y0*momentz0 + y1*momentz1);
            c(2, 4) += (x0*momentz0 + x1*momentz1);
            c(3, 3) -= (y0*y0*momentz0 + y1*y1*momentz1);
            c(3, 4) += (x0*y0*momentz0 + x1*y1*momentz1);
            c(4, 4) -= (x0*x0*momentz0 + x1*x1*momentz1);
		}
	}

	c(2, 2) = c(2, 2)*rho*g;
	c(2, 3) = c(2, 3)*rho*g;
	c(2, 4) = c(2, 4)*rho*g;
	c(3, 4) = c(3, 4)*rho*g;
	c(3, 3) = c(3, 3)*rho*g - rho*g*volume*cb.z + mass*g*cg.z;
	c(4, 4) = c(4, 4)*rho*g - rho*g*volume*cb.z + mass*g*cg.z;
	c(3, 5) = rho*g*volume*cb.x - mass*g*cg.x;
	c(4, 5) = rho*g*volume*cb.y - mass*g*cg.y;
	
	c(3, 2) = c(2, 3);
	c(4, 2) = c(2, 4);
	c(4, 3) = c(3, 4); 
}


inline static void AddSegZero(Vector<Segment3D> &seg, const Point3D &p0, const Point3D &p1, 
			const Point3D &p2, const Point3D &p3) {
	if (abs(p0.z - 0) <= 0.01 && abs(p1.z - 0) <= 0.01)
		seg << Segment3D(p0, p1);
	if (abs(p1.z - 0) <= 0.01 && abs(p2.z - 0) <= 0.01)
		seg << Segment3D(p1, p2);
	if (abs(p2.z - 0) <= 0.01 && abs(p3.z - 0) <= 0.01)
		seg << Segment3D(p2, p3);
	if (abs(p3.z - 0) <= 0.01 && abs(p0.z - 0) <= 0.01)
		seg << Segment3D(p3, p0);
}

void Surface::Underwater(const Surface &orig) {
	nodes = clone(orig.nodes);
	panels.Clear();
	Fixed();
	
	segWaterlevel.Clear();
	
	for (int ip = 0; ip < orig.panels.GetCount(); ++ip) {
		const int &id0 = orig.panels[ip].id[0];
		const int &id1 = orig.panels[ip].id[1];
		const int &id2 = orig.panels[ip].id[2];
		const int &id3 = orig.panels[ip].id[3];
		const Point3D &p0 = nodes[id0];
		const Point3D &p1 = nodes[id1];
		const Point3D &p2 = nodes[id2];
		const Point3D &p3 = nodes[id3];	
		
		if (p0.z >= 0 && p1.z >= 0 && p2.z >= 0 && p3.z >= 0) {
			AddSegZero(segWaterlevel, p0, p1, p2, p3);
		} else if (p0.z <= 0 && p1.z <= 0 && p2.z <= 0 && p3.z <= 0) {
			AddSegZero(segWaterlevel, p0, p1, p2, p3);
			panels << Panel(orig.panels[ip]);
		} else {
			const int *origPanelid = orig.panels[ip].id;
			Vector<int> nodeFrom, nodeTo;
			Segment3D segWL;
			segWL.from = segWL.to = Null;
			const int ids[] = {0, 1, 2, 3, 0};
			for (int i = 0; i < 4; ++i) {
				Point3D from(nodes[origPanelid[ids[i]]]);
				Point3D to(nodes[origPanelid[ids[i+1]]]);
				Segment3D seg(from, to);
				if (abs(from.z - 0) <= 0.01 && abs(to.z - 0) <= 0.01) {
					nodeFrom << origPanelid[ids[i]];
					nodeTo << origPanelid[ids[i+1]];
					segWL.from = from;
					segWL.to = to;	
				} else if (from.z <= 0 && to.z <= 0) {
					nodeFrom << origPanelid[ids[i]];
					nodeTo << origPanelid[ids[i+1]];
				} else if (from.z >= 0 && to.z >= 0) 
					;
				else {
					Point3D inter = seg.IntersectionPlaneZ(0);
					if (from.z < 0) {
						nodeFrom << origPanelid[ids[i]];
						nodes << inter;
						nodeTo << nodes.GetCount() - 1;
					} else {
						nodeTo << origPanelid[ids[i+1]];
						nodes << inter;
						nodeFrom << nodes.GetCount() - 1;
					}
					if (IsNull(segWL.from))
						segWL.from = inter;
					else if (IsNull(segWL.to))
						segWL.to = inter;
				}
			}
			
			if (!IsNull(segWL))
				segWaterlevel << segWL;
			
			int pos = -1, nFrom, nTo;
			for (int i = 0; i < nodeFrom.GetCount(); ++i) {
				int i_1 = i + 1;
				if (i_1 >= nodeFrom.GetCount())
					i_1 = 0;
				if (nodeTo[i] != nodeFrom[i_1]) {
					pos = i+1;
					nFrom = nodeTo[i];
					nTo = nodeFrom[i_1];
					break;
				}
			}
			if (pos == nodeTo.GetCount()) {
				nodeFrom << nFrom;
				nodeTo << nTo;
			} else if (pos >= 0) {
				nodeFrom.Insert(pos, nFrom);		
				nodeTo.Insert(pos, nTo);
			}
			Panel panel;
			if (nodeFrom.GetCount() == 3) {
				panel.id[0] = nodeFrom[0];
				panel.id[1] = nodeFrom[1];
				panel.id[2] = nodeFrom[2];
				panel.id[3] = nodeFrom[2];
			} else if (nodeFrom.GetCount() == 4) {
				panel.id[0] = nodeFrom[0];
				panel.id[1] = nodeFrom[1];
				panel.id[2] = nodeFrom[2];
				panel.id[3] = nodeFrom[3];
			} else if (nodeFrom.GetCount() == 5) {
				panel.id[0] = nodeFrom[0];
				panel.id[1] = nodeFrom[1];
				panel.id[2] = nodeFrom[2];
				panel.id[3] = nodeFrom[3];
				Panel panel2;
				panel2.id[0] = nodeFrom[0];
				panel2.id[1] = nodeFrom[3];
				panel2.id[2] = nodeFrom[4];
				panel2.id[3] = nodeFrom[4];
				panels << panel2;
			}
			panels << panel;
		}
	}
}

bool Surface::IsMoved(double _x, double _y, double _z, double _a_x, double _a_y, double _a_z) const {
	return x != _x || y != _y || z != _z || a_x != _a_x || a_y != _a_y || a_z != _a_z;
}

void Surface::MoveTo() {
	if (IsFixed())
		throw(t_("Error moving a fixed surface"));
	
	nodes.SetCount(nodes0.GetCount());
	
	Affine3d quat;
	GetTransform(quat, a_x, a_y, a_z, c_x, c_y, c_z);
	
	for (int i = 0; i < nodes0.GetCount(); ++i) 
		nodes[i].MoveTo(nodes0[i], x, y, z, quat); 
}

void Surface::MoveTo(double _x, double _y, double _z, double _a_x, double _a_y, double _a_z, double _c_x, double _c_y, double _c_z) {
	x = _x;
	y = _y;
	z = _z;
	a_x = _a_x;
	a_y = _a_y;
	a_z = _a_z;
	c_x = _c_x;
	c_y = _c_y;
	c_z = _c_z;
	
	MoveTo();
}

void Surface::DeployXSymmetry() {
	int nnodes = nodes0.GetCount();
	for (int i = 0; i < nnodes; ++i) {
		Point3D 	  &dest = nodes0.Add();
		const Point3D &orig = nodes0[i];
		dest.x = -orig.x;
		dest.y =  orig.y;
		dest.z =  orig.z;
	}
	int npanels = panels.GetCount();
	for (int i = 0; i < npanels; ++i) {
		Panel 		&dest = panels.Add();
		const Panel &orig = panels[i];
		dest.id[0] = orig.id[3] + nnodes;
		dest.id[1] = orig.id[2] + nnodes;
		dest.id[2] = orig.id[1] + nnodes;
		dest.id[3] = orig.id[0] + nnodes;
	}
}

void Surface::DeployYSymmetry() {
	int nnodes = nodes0.GetCount();
	for (int i = 0; i < nnodes; ++i) {
		Point3D 	  &dest = nodes0.Add();
		const Point3D &orig = nodes0[i];
		dest.x =  orig.x;
		dest.y = -orig.y;
		dest.z =  orig.z;
	}
	int npanels = panels.GetCount();
	for (int i = 0; i < npanels; ++i) {
		Panel 		&dest = panels.Add();
		const Panel &orig = panels[i];
		dest.id[0] = orig.id[3] + nnodes;
		dest.id[1] = orig.id[2] + nnodes;
		dest.id[2] = orig.id[1] + nnodes;
		dest.id[3] = orig.id[0] + nnodes;
	}
}

void VolumeEnvelope::MixEnvelope(VolumeEnvelope &env) {
	maxX = maxNotNull(env.maxX, maxX);
	minX = minNotNull(env.minX, minX);
	maxY = maxNotNull(env.maxY, maxY);
	minY = minNotNull(env.minY, minY);
	maxZ = maxNotNull(env.maxZ, maxZ);
	minZ = minNotNull(env.minZ, minZ);
}
