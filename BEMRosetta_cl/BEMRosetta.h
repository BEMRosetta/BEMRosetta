#ifndef _BEM_Rosetta_BEM_Rosetta_h_
#define _BEM_Rosetta_BEM_Rosetta_h_

#include <Core/Core.h>
#include <SysInfo/SysInfo.h>

using namespace Upp;

#include <plugin/Eigen/Eigen.h>

using namespace Eigen;


class Hydro {
public:
	enum BEM_SOFT {WAMIT, FAST_WAMIT, WAMIT_1_3, NEMOH, SEAFEM_NEMOH};
	
	void Save(String file, BEM_SOFT type);
	void Report();
	Hydro() : g(Null), h(Null), rho(Null), Nb(Null), Nf(Null), Nh(Null) {}
	virtual ~Hydro() {}	

	static Function <void(String)> Print, PrintError;	

	static bool MatchCoeffStructure(Upp::Array<Hydro> &hydro, String &strError);
	
	String GetCodeStr()	{
		switch (code) {
		case WAMIT: 		return "Wamit";
		case WAMIT_1_3: 	return "Wamit.1.3";
		case FAST_WAMIT: 	return "FAST-Wamit";
		case NEMOH:			return "Nemoh";
		case SEAFEM_NEMOH:	return "SeaFEM-Nemoh";
		}
		return "Unknown";
	}
	
	inline bool IsAvailableDOF(int ib, int idof) {
		return (Awinf.size() > 0 && !IsNaN(Awinf((ib+1)*idof, (ib+1)*idof))) || 
			   (!A.IsEmpty() && A[0].size() > 0 && !IsNaN(A[0]((ib+1)*idof, (ib+1)*idof)));
	}

	String file;        	// BEM output filename
	String name;
    double g;           	// gravity
    double h;           	// water depth
   	double rho;        		// density
    int Nb;          		// number of bodies
    int Nf;          		// number of wave frequencies
    int Nh;          		// number of wave headings
 	double len;
 	
	Upp::Array<MatrixXd> A;	// [Nf](6*Nb, 6*Nb)		Added mass
    MatrixXd Awinf;        	// (6*Nb, 6*Nb)        	Infinite frequency added mass
    MatrixXd Aw0;        	// (6*Nb, 6*Nb)       	Infinite period added mass
    Upp::Array<MatrixXd> B; // [Nf](6*Nb, 6*Nb)    	Radiation damping
    Vector<double> head;	// [Nh]                	Wave headings (deg)
    Vector<String> names;  	// {Nb}                	Body names
    Upp::Array<MatrixXd> C; // [Nb](6, 6)          	Hydrostatic restoring coefficients:
    MatrixXd cb;          	// (3,Nb)            	Center of buoyancy
    MatrixXd cg;          	// (3,Nb)     			Center of gravity
    BEM_SOFT code;        	// BEM_SOFT				BEM code (WAMIT, AQWA, or NEMOH)
    Vector<int> dof;      	// [Nb]             	Degrees of freedom for each body 
    Vector<int> dofOrder;	//						Order of DOF
    
    struct Forces {
    	Upp::Array<MatrixXd> ma, ph;   	// [Nh](Nf, 6*Nb) 	Magnitude and phase
    	Upp::Array<MatrixXd> re, im;	// [Nh](Nf, 6*Nb)	Real and imaginary components
    };
    
    Forces ex; 				// Excitation
    Forces sc; 				// Diffraction scattering
    Forces fk; 				// Froude-Krylov
    
    Vector<double> T; 		// [Nf]    				Wave periods
    Vector<double> w;      	// [Nf]               	Wave frequencies
    Vector<double> Vo;    	// [Nb]             	Displaced volume
    
    void Dimensionalize();
    void Normalize();
    
    static String C_units(int i, int j);
    
protected:
	void AfterLoad();
	void Normalize_Forces(Forces &f);
	void Dimensionalize_Forces(Forces &f);
	void Initialize_Forces();
	void Initialize_Forces(Forces &f);
	inline int GetK_AB(int i, int j);
	inline int GetK_F(int i);
	inline int GetK_C(int i, int j);
	
	void GetBodyDOF();
	
	String lastError;

private:
	static const char *strDOF[6];
	static const char *strDOFAbrev[6];
	static String C_units_base(int i, int j);
		
public:
	static String StrDOF(int i) {
		int ib = i/6 + 1;
		int idof = i - (ib - 1)*6;
		return Format("%d.%s", ib, strDOF[idof]);
	}
	
	static int DOFStr(String &str) {
		for (int i = 0; i < 6; ++i)
			if (strDOF[i] == ToLower(str))
				return i;
		return -1;
	}
	
	static void DOFFromStr(String str, int &ib, int &idof) {
		int pos = str.Find(".");
		ib = ScanInt(str.Left(pos-1));
		String sdof = str.Mid(pos+1);
		idof = DOFStr(sdof);	
	}
	
	static String StrDOFAbrev(int i) {
		int nb = i/6 + 1;
		int ni = i - (nb - 1)*6;
		return Format("%d%s", nb, strDOFAbrev[ni]);
	}
	
	bool IsLoadedA() 	{return A.GetCount() > 0;}
	bool IsLoadedAwinf(){return Awinf.size() > 0;}
	bool IsLoadedAw0()	{return Aw0.size() > 0;}
	bool IsLoadedB() 	{return B.GetCount() > 0;}
	bool IsLoadedC()	{return C.size() > 0;}
	bool IsLoadedFex() 	{return ex.ma.GetCount() > 0;}
	bool IsLoadedFsc() 	{return sc.ma.GetCount() > 0;}
	bool IsLoadedFfk() 	{return fk.ma.GetCount() > 0;}
	
	const Vector<int> &GetOrder()		{return dofOrder;}
	void SetOrder(Vector<int> &order)	{dofOrder = pick(order);}
	
	String GetLastError()	{return lastError;}
};


class Wamit : public Hydro {
public:
	bool Load(String file, double rho = 1000);
	void Save(String file);
	
private:	
	bool Load_out();
	void Load_A(FileIn &in, MatrixXd &A);
	bool Load_Scattering(String fileName);
	bool Load_FK(String fileName);
};

class Fast : public Hydro {
public:
	Fast() : WaveNDir(Null), WaveDirRange(Null) {}
	bool Load(String file, double g = 9.81);
	void Save(String file, bool isFast);
	
private:
	bool Load_dat();	
	bool Load_1(String fileName);
	bool Load_3(String fileName);
	bool Load_hst(String fileName);
	
	void Save_dat(String fileName, bool force);	
	void Save_1(String fileName);
	void Save_3(String fileName);
	void Save_hst(String fileName);
	
	String hydroFolder;
	bool readW;
	int WaveNDir;
	double WaveDirRange;
};

class Nemoh : public Hydro {
public:
	bool Load(String file, double rho = Null);
	void Save(String file);
	
private:
	String folder;
	
	bool Load_Cal(String fileName);
	bool Load_Inf(String fileName);
	bool Load_Hydrostatics();
	bool Load_KH();
	bool Load_Radiation(String fileName);
	bool Load_Excitation(String folder);
	bool Load_Diffraction(String folder);
	bool Load_FroudeKrylov(String folder);
	bool Load_Forces(Forces &f, String nfolder, String fileName, String textDelim);
	bool Load_IRF(String fileName);
};

template <class T>
void LinSpaced(Vector<T> &v, int n, T min, T max) {
	ASSERT(n > 0);
	v.SetCount(n);
	if (n == 1)
		v[0] = min;
	else {
		for (int i = 0; i < n; ++i)
			v[i] = min + (max - min)*i/(n - 1);
	}
}

int IsTabSpace(int c);

class FieldSplit {
public:
	void Load(String line) {
		this->line = line;
		fields = Split(line, IsTabSpace, true);
	}
	void LoadWamitJoinedFields(String line) {		// Trick for "glued" fields in Wamit
		this->line = line;
		fields.Clear();
		Vector<String> prefields = Split(line, IsTabSpace, true);
		for (int id = 0; id < prefields.GetCount(); ++id) {
			String s = prefields[id];
			String ns;
			for (int i = 0; i < s.GetCount(); ++i) {	
				int c = s[i];
				if (c == '-') {
					if (i == 0)
						ns.Cat(c);
					else if (s[i-1] == 'E')
						ns.Cat(c);
					else {
						fields << ns;
						ns.Clear();
						ns.Cat(c);
					}
				} else
					ns.Cat(c);
			}
			fields << ns;
		}
	}
	String GetText(int i) {
		CheckId(i);
		return fields[i];
	}
	int GetInt(int i) {
		CheckId(i);
		int res = ScanInt(fields[i]);
		if (IsNull(res))
			throw Exc(Format("Bad integer '%s' in field #%d, line '%s'", fields[i], i+1, line));
		return res;
	}
	double GetDouble(int i) {
		CheckId(i);
		double res = ScanDouble(fields[i]);
		if (IsNull(res))
			throw Exc(Format("Bad double '%s' in field #%d, line '%s'", fields[i], i+1, line));
		return res;
	}
private:
	String line;
	Vector<String> fields;
	
	void CheckId(int i) {
		if (i >= fields.GetCount())
			throw Exc(Format("Field #%d not found in line '%s'", i+1, line));
	}
};

template <class Range, class V>
void FindAdd(Range& r, const V& value, int from = 0)
{
	for(int i = from; i < r.GetCount(); i++)
		if(r[i] == value) 
			return;
	r.Add(value);
}

String FormatWam(double d);

#define TV(v, i)		TV_(v, i, #v, #i, __FILE__, __LINE__)

template <class T>
inline T& TV_(T &v, int i, const char *strvar, const char *strindex, const char *strfile, int strline) {
	if (i < 0 || i >= v.GetCount()) 
		throw Exc(Format("Array %s[%s] is out of limits in file '%s', line '%d'", strvar, strindex, strfile, strline));
	return v;
}
	
#endif
