#ifndef _GLCanvas_surface_h_
#define _GLCanvas_surface_h_


template<class T>
inline T avg(T a, T b) 			{return T(a+b)/2;}
template<class T>
inline T avg(T a, T b, T c)		{return T(a+b+c)/3;}

const double EPS = 0.00001;
 
template<class T>
void Sort(T& a, T& b, T& c, T& d) {
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
void Sort(T& a, T& b, T& c) {
	if (a > b) 
		Swap(a, b);
	if (a > c) 
		Swap(a, c);
	if (b > c) 
		Swap(b, c);
}
	
class Point3D : public Moveable<Point3D> {
public:
	double x, y, z;

	Point3D() {}
	Point3D(const Nuller&) {x = Null;}
	Point3D(const Point3D &p) : x(p.x), y(p.y), z(p.z) {}
	Point3D(double x, double y, double z) : x(x), y(y), z(z) {}
	
	void Set(const Point3D &p) 				{this->x = p.x; this->y = p.y; this->z = p.z;}
	void Set(double x, double y, double z) 	{this->x = x; this->y = y; this->z = z;}
	
	String ToString() const { return FormatDouble(x) + "," + FormatDouble(y) + "," + FormatDouble(z); }
	
	bool IsSimilar(const Point3D &p, double similThres) {
		if (abs(p.x - x) < similThres && abs(p.y - y) < similThres && abs(p.z - z) < similThres)
			return true;
		return false;
	}
	
	friend Point3D operator%(const Point3D& a, const Point3D& b)  {return Point3D(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);}
	friend Point3D operator+(const Point3D& a, const Point3D& b)  {return Point3D(a.x+b.x, a.y+b.y, a.z+b.z);}
	friend Point3D operator-(const Point3D& a, const Point3D& b)  {return Point3D(a.x-b.x, a.y-b.y, a.z-b.z);}

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
	double Manhattan(const Point3D &p) 	{return abs(x-p.x) + abs(y-p.y) + abs(z-p.z);}
	double Manhattan() 					{return abs(x) + abs(y) + abs(z);}
	
	void SimX() {x = -x;}
	void SimY() {y = -y;}
	void SimZ() {z = -z;}
	
	void Mirror(const Point3D &p0) {
		x = 2*p0.x - x;
		y = 2*p0.y - y;
		z = 2*p0.z - z;
	}
};

typedef Point3D Vector3D;

class Segment3D : public Moveable<Segment3D> {
public:
	Point3D from, to;
	
	Segment3D() {}
	Segment3D(const Nuller&) {SetNull();}
	Segment3D(const Point3D &from, const Point3D &to) : from(from), to(to) {}
	Segment3D(const Point3D &from, const Vector3D &normal, double length) : from(from) {
		to = Point3D(from.x + length*normal.x, from.y + length*normal.y, from.z + length*normal.z);
	}
	void SetNull() {from = Null;}
		
	void Set(const Point3D &from, const Point3D &to) {
		this->from.Set(from);
		this->to.Set(to);
	}
	void Set(const Point3D &from, const Vector3D &normal, double length) {
		this->from.Set(from);
		to.Set(from.x + length*normal.x, from.y + length*normal.y, from.z + length*normal.z);
	}
	void SimX() {
		from.SimX();
		to.SimX();
	}
	void SimY() {
		from.SimY();
		to.SimY();
	}
	void SimZ() {
		from.SimZ();
		to.SimZ();
	}
	double Length() {return from.Distance(to);}
	double Dx()		{return to.x - from.x;}
	double Dy()		{return to.y - from.y;}	
	double Dz()		{return to.z - from.z;}
	
	void Mirror(const Point3D &p0) {
		from.Mirror(p0);
		to.Mirror(p0);
	}
};

bool MinSegment(const Point3D& p1, const Point3D& p2, const Point3D& p3, const Point3D& p4, Segment3D &ret, double &mua, double &mub);
bool Cross(const Point3D& p1, const Point3D& p2, const Point3D& p3, const Point3D& p4);
Point3D GetCentroid(const Point3D &a, const Point3D &b);
Point3D GetCentroid(const Point3D &a, const Point3D &b, const Point3D &c);
Vector3D GetNormal(const Point3D &a, const Point3D &b, const Point3D &c);

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

class Panel : public Moveable<Panel> {
public:
	int id[4];
	Segment3D normal;

	bool operator==(const Panel &p) const {
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
	inline void Swap() {
		if (IsTriangle())
			::Swap(id[1], id[2]);
		else
			::Swap(id[1], id[3]);
	}
	inline bool IsTriangle() const	{return id[0] == id[3];}
	inline int GetNumNodes() const	{return IsTriangle() ? 3 : 4;}
	bool FirstNodeIs0(int in0, int in1) const;
	
	String ToString() const { return FormatInt(id[0]) + "," + FormatInt(id[1]) + "," + FormatInt(id[2]) + "," + FormatInt(id[3]); }
};

class Segment : public Moveable<Segment> {
public:
	int inode0, inode1;
	Index<int> panels;
};

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
	
	Vector<Segment3D> skewed;
	Vector<Segment3D> segWaterlevel, segTo1panel, segTo3panel;
	
	bool x0z, y0z;
	String file;
	VolumeEnvelope env;
	
	String Heal(Function <void(String, int pos)> Status);
	void GetLimits(); 
	void GetNormals();
	String GetLastError()	{return lastError;}
	String lastError;

private:
	inline bool CheckId(int id) {return id >= 0 && id < nodes.GetCount()-1;}
	
	void DetectTriBiP(int &numTri, int &numBi, int &numP);
	int FixSkewed();
	int RemoveDuplicatedPanels();
	int RemoveDuplicatedPointsAndRenumber();
	void AnalyseSegments(Vector<Segment> &segments, double zTolerance);
	void AddSegment(Vector<Segment> &segments, int ip0, int ip1, int ipanel);
	bool ReorientPanels(const Vector<Segment> &segments, int &numUnprocessed);
	void ReorientPanel(int ip);
	bool SameOrderPanel(int ip0, int ip1, int in0, int in1);
	int PanelGetNumNodes(int ip) 	{return panels[ip].GetNumNodes();}
	bool IsPanelTriangle(int ip) 	{return panels[ip].IsTriangle();}
};

#endif
