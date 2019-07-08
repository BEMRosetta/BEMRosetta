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
	return num;
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
		for (int j = 0; j < 4; ++j) {
			int &id = panels[i].id[j];
			int pos = similars.Find(id);
			if (pos >= 0)
				id = goods[pos];
		}
	}
	
	// Find duplicates
	Vector<int> newId;
	newId.SetCount(nodes.GetCount());
	int avId = 0;
	for (int i = 0; i < nodes.GetCount(); ++i) {
		if (similars.Find(i) < 0) {
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
		for (int j = 0; j < 4; ++j) {
			int& id = panels[i].id[j];
			id = newId[id];
		}
	}
	
	return num;
}

String Surface::Heal() {
	String ret;
	
	int numDup = RemoveDuplicatedPanels();
	
	int numSkewed = FixSkewed();
	if (numSkewed > 0) 
		ret << "\n" << Format(t_("Fixed %d skewed panels"), numSkewed);

	int numDupP = RemoveDuplicatedPointsAndRenumber();
	if (numDupP > 0) 
		ret << "\n" << Format(t_("Removed %d duplicated points"), numDupP);	

	numDup += RemoveDuplicatedPanels();	// Second time after duplicated points
	if (numDup > 0) 
		ret << "\n" << Format(t_("Removed %d duplicated panels"), numDup);

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
	normals.SetCount(panels.GetCount());
	for (int ip = 0; ip < panels.GetCount(); ++ip) {
		Panel &panel = panels[ip];
		Point3D p0 = nodes[panel.id[0]];
		Point3D p1 = nodes[panel.id[1]];
		Point3D p2 = nodes[panel.id[2]];
		Point3D p3 = nodes[panel.id[3]];
		
		Point3D from = GetCentroid(p0, p2);		
		Point3D pnormal = GetNormal(p0, p1, p2);
		normals[ip].Set(from, pnormal, 1);
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
