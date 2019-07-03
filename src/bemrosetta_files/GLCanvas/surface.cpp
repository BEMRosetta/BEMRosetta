#include <Core/Core.h>

using namespace Upp;

#include "surface.h"

void Surface::Compact(const Index<int> &nIds) {	// Remove nIds nodes and renumber in panels
	
	
	
}

void Surface::Heal() {						// Clean surface
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
	Compact(similars);



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