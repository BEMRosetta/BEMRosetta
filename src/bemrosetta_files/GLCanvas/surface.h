#ifndef _GLCanvas_surface_h_
#define _GLCanvas_surface_h_

class Panel : public Moveable<Panel> {
public:
	int id[4];

	String ToString() const { return FormatInt(id[0]) + "," + FormatInt(id[1]) + "," + FormatInt(id[2]) + "," + FormatInt(id[3]); }
};

class Point3D : public Moveable<Point3D> {
public:
	double x, y, z;

	String ToString() const { return FormatDouble(x) + "," + FormatDouble(y) + "," + FormatDouble(z); }
};


template <typename T>
inline T const& maxNotNull(T const& a, T const& b) {
	if (IsNull(a))
		return b;
	else if (IsNull(b))
		return a;
	else
    	return a > b ? a : b;
}

template <typename T>
inline T const& minNotNull(T const& a, T const& b) {
	if (IsNull(a))
		return b;
	else if (IsNull(b))
		return a;
	else
    	return a < b ? a : b;
}


class VolumeEnvelope {
public:
	VolumeEnvelope() {Reset();}
	void Reset() {
		maxX = minX = maxY = minY = maxZ = minZ = Null;
	}
	
	void MixEnvelope(VolumeEnvelope &env) {
		maxX = maxNotNull(env.maxX, maxX);
		minX = minNotNull(env.minX, minX);
		maxY = maxNotNull(env.maxY, maxY);
		minY = minNotNull(env.minY, minY);
		maxZ = maxNotNull(env.maxZ, maxZ);
		minZ = minNotNull(env.minZ, minZ);
	}
	
	double maxX, minX, maxY, minY, maxZ, minZ;
};

class Surface {
public:
	Surface() : x0z(false), y0z(false) {}
	Vector<Point3D> nodes;
	Vector<Panel> panels;
	bool x0z, y0z;
	String file;
	VolumeEnvelope env;
	/*
	bool Check() {
		for (int i = 0; i < panels.GetCount(); ++i)
			if (!panels[i].id[0] || !panels[i].id[1] || !panels[i].id[2] || !panels[i].id[3])
				return false;
		return true;
	}*/
	void GetLimits() {
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
	String GetLastError()	{return lastError;}
	String lastError;
	
private:
	inline bool CheckId(int id) {return id >= 0 && id < nodes.GetCount()-1;}
};

#endif
