#include <Core/Core.h>

using namespace Upp;

#include <Surface/Surface.h>

Point3D GetCentroid(Point3D &a, Point3D &b) {
	return Point3D(avg(a.x, b.x), avg(a.y, b.y), avg(a.z, b.z));	
}

Point3D GetCentroid(Point3D &a, Point3D &b, Point3D &c) {
	return Point3D(avg(a.x, b.x,c.x), avg(a.y, b.y, c.y), avg(a.z, b.z, c.z));	
}

Vector3D GetNormal(Point3D &a, Point3D &b, Point3D &c) {
	return Vector3D((a - b) % (b - c)).Normalize();
}

int Surface::FixSkewed() {	
	int num = 0;
	for (int i = 0; i < panels.GetCount(); ++i) {
		Vector3D normal012 = GetNormal(nodes[panels[i].id[0]], nodes[panels[i].id[1]], nodes[panels[i].id[2]]);
		Vector3D normal123 = GetNormal(nodes[panels[i].id[1]], nodes[panels[i].id[2]], nodes[panels[i].id[3]]);
		double dman = normal012.Manhattan(normal123);
		double d0   = normal012.Manhattan();
		if (dman > d0) {
			num++;
			Swap(panels[i].id[0], panels[i].id[2]);
		}
	}	
	return num;
}

int Surface::RemoveDuplicatedPanels() {		
	int num = 0;
	for (int i = 0; i < panels.GetCount()-1; ++i) {
		Panel &panel = panels[i];
		for (int j = i+1; j < panels.GetCount(); ++j) {
			if (panel == panels[j]) {
				num++;
				panels.Remove(j, 1);
			}
		}
	}
	return num;
}

int Surface::RemoveDuplicatedPoints() {
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
	// Remove duplicate points in panels
	for (int i = 0; i < panels.GetCount(); ++i) {
		for (int j = 0; j < 4; ++j) {
			int id = panels[i].id[j];
			int pos = similars.Find(id);
			if (pos >= 0)
				panels[i].id[j] = goods[id];
		}
	}
	// Renumber nodes
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
	// Remove nodes
	for (int i = nodes.GetCount()-1; i >= 0; --i) {
		if (IsNull(newId[i]))
			nodes.Remove(i, 1);
	}
	// Renumber panels
	for (int i = 0; i < panels.GetCount(); ++i) {
		for (int j = 0; j < 4; ++j) {
			int id = panels[i].id[j];
			panels[i].id[j] = newId[id];
		}
	}
	return num;
}
		
String Surface::Heal() {						// Clean surface
	return "";
	
	String ret;
	
	int numDup = RemoveDuplicatedPanels();
	if (numDup > 0) 
		ret << Format(t_("Removed %d duplicated panels"), numDup) << "\n";
	
	int numSkewed = FixSkewed();
	if (numSkewed > 0) 
		ret << Format(t_("Fixed %d skewed panels"), numSkewed) << "\n";
	 
	int numDupP = RemoveDuplicatedPoints();
	if (numDupP > 0) 
		ret << Format(t_("Removed %d duplicated points"), numDupP) << "\n";	

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

void VolumeEnvelope::MixEnvelope(VolumeEnvelope &env) {
	maxX = maxNotNull(env.maxX, maxX);
	minX = minNotNull(env.minX, minX);
	maxY = maxNotNull(env.maxY, maxY);
	minY = minNotNull(env.minY, minY);
	maxZ = maxNotNull(env.maxZ, maxZ);
	minZ = minNotNull(env.minZ, minZ);
}