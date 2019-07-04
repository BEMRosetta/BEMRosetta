#ifndef _GLCanvas_surface_h_
#define _GLCanvas_surface_h_


template<class T>
inline T avg(T a, T b) 			{return T(a+b)/2;}
template<class T>
inline T avg(T a, T b, T c)		{return T(a+b+c)/3;}

template<class T>
void Sort(T a, T b, T c, T d) {
	if (a > b) 
		Swap(a, b);
	if (c > d) 
		Swap(c, d);
	if (a > c) 
		Swap(a, c);
	if (b > d) 
		Swap(b, d);
	if (b > c) 
		Swap(b, c);
}

template<class T>
void Sort(T a, T b, T c) {
	if (a > b) 
		Swap(a, b);
	if (a > c) 
		Swap(a, c);
	if (b > c) 
		Swap(b, c);
}
	
class Panel : public Moveable<Panel> {
public:
	int id[4];

	bool operator==(const Panel &p) {
		int id0 = id[0], id1 = id[1], id2 = id[2], id3 = id[3];
		int pid0 = p.id[0], pid1 = p.id[1], pid2 = p.id[2], pid3 = p.id[3];
		if (id0 + id1 + id2 + id3 != pid0 + pid1 + pid2 + pid3)
			return false;
		Sort(id0, id1, id2, id3);
		Sort(pid0, pid1, pid2, pid3);
		if (id0 == pid0 && id1 == pid1 && id2 == pid2 && id3 == pid3)
			return true;
		return false;
	}
	
	String ToString() const { return FormatInt(id[0]) + "," + FormatInt(id[1]) + "," + FormatInt(id[2]) + "," + FormatInt(id[3]); }
};

class Point3D : public Moveable<Point3D> {
public:
	double x, y, z;

	Point3D() {}
	Point3D(const Point3D &p) : x(p.x), y(p.y), z(p.z) {}
	Point3D(double x, double y, double z) : x(x), y(y), z(z) {}
	
	String ToString() const { return FormatDouble(x) + "," + FormatDouble(y) + "," + FormatDouble(z); }
	
	bool IsSimilar(const Point3D &p, double similThres) {
		if (abs(p.x - x) < similThres && abs(p.y - y) < similThres && abs(p.z - z) < similThres)
			return true;
		return false;
	}
	Point3D operator%(const Point3D &p) {return Point3D(y*p.z - z*p.y, z*p.x - x*p.z, x*p.y - y*p.x);}
	Point3D operator+(const Point3D &p) {return Point3D(x+p.x, y+p.y, z+p.z);}
	Point3D operator-(const Point3D &p) {return Point3D(x-p.x, y-p.y, z-p.z);}
	double GetLength()	{return sqrt(x*x + y*y + z*z);}
	Point3D &Normalize() {
		double length = GetLength();
		
		if (length == 0) 
			x = y = z = 0;
		else {
		    x = x/length;
		    y = y/length;
		    z = z/length;
		}
		return *this;
	}
	double Distance(const Point3D &p) 	{return sqrt(sqr(x-p.x) + sqr(y-p.y) + sqr(z-p.z));}
	double Manhattan(const Point3D &p) 	{return (x-p.x) + (y-p.y) + (z-p.z);}
	double Manhattan() 					{return x + y + z;}
};

typedef Point3D Vector3D;

class Line3D : public Moveable<Line3D> {
public:
	Point3D from, to;
	
	Line3D() {}
	Line3D(const Point3D &from, const Point3D &to) : from(from), to(to) {}
	Line3D(const Point3D &from, const Vector3D &normal, double length) : from(from) {
		to = Point3D(from.x + length*normal.x, from.y + length*normal.y, from.z + length*normal.z);
	}
};

Point3D GetCentroid(Point3D &A, Point3D &B);
Point3D GetCentroid(Point3D &A, Point3D &B, Point3D &C);
Vector3D GetNormal(Point3D &A, Point3D &B, Point3D &C);

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
	void Reset() {maxX = minX = maxY = minY = maxZ = minZ = Null;}
	
	void MixEnvelope(VolumeEnvelope &env);
	
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
	
	String Heal();
	void GetLimits(); 
	String GetLastError()	{return lastError;}
	String lastError;
	
private:
	inline bool CheckId(int id) {return id >= 0 && id < nodes.GetCount()-1;}
	
	int FixSkewed();
	int RemoveDuplicatedPanels();
	int RemoveDuplicatedPoints();
};

#endif
