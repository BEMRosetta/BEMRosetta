#include <Core/Core.h>

using namespace Upp;

#include <Surface/Surface.h>

Point3D GetCentroid(const Point3D &a, const Point3D &b) {
	return Point3D(avg(a.x, b.x), avg(a.y, b.y), avg(a.z, b.z));	
}

Point3D GetCentroid(const Point3D &a, const Point3D &b, const Point3D &c) {
	return Point3D(avg(a.x, b.x,c.x), avg(a.y, b.y, c.y), avg(a.z, b.z, c.z));	
}

Vector3D GetNormal(const Point3D &a, const Point3D &b, const Point3D &c) {
	return Vector3D((a - b) % (b - c)).Normalize();
}

int Surface::FixSkewed() {	
	double similThres = 0.00001;
	int num = 0;
	for (int i = 0; i < panels.GetCount(); ++i) {
		int &id0 = panels[i].id[0];
		int &id1 = panels[i].id[1];
		int &id2 = panels[i].id[2];
		int &id3 = panels[i].id[3];
		Point3D &p0 = nodes[id0];
		Point3D &p1 = nodes[id1];
		Point3D &p2 = nodes[id2];
		Point3D &p3 = nodes[id3];
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
		Index<int> ids;
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
	Index<int> similars, goods;
	for (int i = 0; i < nodes.GetCount()-1; ++i) {
		if (similars.Find(i) >= 0)
			continue;
		for (int j = i+1; j < nodes.GetCount(); ++j) {
			if (nodes[i].IsSimilar(nodes[j], similThres)) {
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
	newId.SetCount(nodes.GetCount());
	int avId = 0;
	for (int i = 0; i < nodes.GetCount(); ++i) {
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
	for (int i = nodes.GetCount()-1; i >= 0; --i) {
		if (IsNull(newId[i]))
			nodes.Remove(i, 1);
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
	
void Surface::AddSegment(Vector<Segment> &segments, int inode0, int inode1, int ipanel) {
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

void Surface::AnalyseSegments(Vector<Segment> &segments, double zTolerance) {
	segments.Clear();
	
	for (int i = 0; i < panels.GetCount(); ++i) {
		AddSegment(segments, panels[i].id[0], panels[i].id[1], i);
		AddSegment(segments, panels[i].id[1], panels[i].id[2], i);
		if (IsPanelTriangle(i)) 
			AddSegment(segments, panels[i].id[2], panels[i].id[0], i);
		else {
			AddSegment(segments, panels[i].id[2], panels[i].id[3], i);
			AddSegment(segments, panels[i].id[3], panels[i].id[0], i);
		}
	}
	
	for (int i = 0; i < segments.GetCount(); ++i) {
		int num = segments[i].panels.GetCount();
		if (y0z && nodes[segments[i].inode0].x == 0 && nodes[segments[i].inode1].x == 0)
			num *= 2;
		if (x0z && nodes[segments[i].inode0].y == 0 && nodes[segments[i].inode1].y == 0)
			num *= 2;
				
		if (num == 1) {
			if (nodes[segments[i].inode0].z >= zTolerance && nodes[segments[i].inode1].z >= zTolerance)
				segWaterlevel << Segment3D(nodes[segments[i].inode0], nodes[segments[i].inode1]);
			else
				segTo1panel << Segment3D(nodes[segments[i].inode0], nodes[segments[i].inode1]);
		}
		if (num > 2)
			segTo3panel << Segment3D(nodes[segments[i].inode0], nodes[segments[i].inode1]);
	}
}

bool Surface::ReorientPanels(const Vector<Segment> &segments, int &numUnprocessed) {
	// Get lowest panel with normal not horizontal
	int iLowSeg = Null;
	int iLowPanel;
	double zLowSeg = DBL_MAX;
	for (int i = 0; i < segments.GetCount(); ++i) {
		const Segment &seg = segments[i];
		if (seg.panels.GetCount() == 2) {
			for (int ip = 0; ip < seg.panels.GetCount(); ++ip) {
				if (panels[seg.panels[ip]].normal.Dz() != 0) {
					double z = max(nodes[seg.inode0].z, nodes[seg.inode1].z);
					if (z < zLowSeg) {
						zLowSeg = z;
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
	if (panels[ip].normal.Dz() > 0)
		ReorientPanel(ip);
	
	Vector<int> panelStack;
	Index<int> panelProcessed;
	
	panelStack << ip;
	while (!panelStack.IsEmpty()) {
		int id = panelStack.GetCount() - 1;
		int ip = panelStack[id];
		panelStack.Remove(id, 1);
		panelProcessed << ip;
		
		for (int is = 0; is < segments.GetCount(); ++is) {
			const Index<int> &segPanels = segments[is].panels;
			if (segPanels.Find(ip) >= 0) {
				for (int i = 0; i < segPanels.GetCount(); ++i) {
					int ipadyac = segPanels[i];
					if (ipadyac != ip && panelProcessed.Find(ipadyac) < 0) {
						panelStack << ipadyac;
						if (!SameOrderPanel(ip, ipadyac, segments[is].inode0, segments[is].inode1))
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
	panels[ip].normal.Mirror(panels[ip].normal.from);
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


bool Surface::SameOrderPanel(int ip0, int ip1, int in0, int in1) {
	bool first0in0 = panels[ip0].FirstNodeIs0(in0, in1);
	bool first1in0 = panels[ip1].FirstNodeIs0(in0, in1);
	
	return first0in0 != first1in0;
}

String Surface::Heal(Function <void(String, int pos)> Status) {
	String ret;
	
	Status(t_("Detecting triangles and wrong panels"), 50);
	int numTri, numBi, numP;
	DetectTriBiP(numTri, numBi, numP);
	if (numTri > 0)
		ret << "\n" << Format(t_("Found %d triangles"), numTri);
	if (numBi > 0)
		ret << "\n" << Format(t_("Removed %d 2 points quads"), numBi);
	if (numP > 0)
		ret << "\n" << Format(t_("Removed %d 1 points quads"), numP);
	
	Status(t_("Removing duplicated panels (pass 1)"), 55);
	int numDup = RemoveDuplicatedPanels();
	
	Status(t_("Fixing skewed panels"), 60);
	int numSkewed = FixSkewed();
	if (numSkewed > 0) 
		ret << "\n" << Format(t_("Fixed %d skewed panels"), numSkewed);

	Status(t_("Removing duplicated points"), 65);
	int numDupP = RemoveDuplicatedPointsAndRenumber();
	if (numDupP > 0) 
		ret << "\n" << Format(t_("Removed %d duplicated points"), numDupP);	

	Status(t_("Removing duplicated panels (pass 2)"), 70);
	numDup += RemoveDuplicatedPanels();	// Second time after duplicated points
	if (numDup > 0) 
		ret << "\n" << Format(t_("Removed %d duplicated panels"), numDup);

	Status(t_("Analysing water tightness"), 75);
	Vector<Segment> segments;
	double zTolerance = -0.1;
	AnalyseSegments(segments, zTolerance);
	ret << "\n" << Format(t_("%d segments, %d water level, %d water leak and %d multipanel"), 
								segments.GetCount(), segWaterlevel.GetCount(), 
								segTo1panel.GetCount(), segTo3panel.GetCount());
	
	GetNormals();
	
	Status(t_("Reorienting panels water side"), 80);
	int numUnprocessed;
	if (!ReorientPanels(segments, numUnprocessed))
		ret << "\n" << t_("Failed to reorient panels to water side");
	else if (numUnprocessed > 0)
		ret << "\n" << Format(t_("%d panels not reoriented. Body contains separated surfaces"), numUnprocessed);
	
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

void Surface::GetNormals() {
	for (int ip = 0; ip < panels.GetCount(); ++ip) {
		Panel &panel = panels[ip];
		Point3D p0 = nodes[panel.id[0]];
		Point3D p1 = nodes[panel.id[1]];
		Point3D p2 = nodes[panel.id[2]];
		Point3D p3 = nodes[panel.id[3]];
		
		Point3D from = GetCentroid(p0, p2);		
		Point3D pnormal = GetNormal(p0, p1, p2);
		panel.normal.Set(from, pnormal, 1);
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
