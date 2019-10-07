#ifndef _BEM_Rosetta_BEM_Rosetta_h_
#define _BEM_Rosetta_BEM_Rosetta_h_

#include <Core/Core.h>
#include <SysInfo/SysInfo.h>
#include <Surface/Surface.h>

using namespace Upp;


class BEMData;
class FieldSplit;

void ConsoleMain(const Vector<String>& command, bool gui);
void SetBuildInfo(String &str);
	
class Hydro {
public:
	enum BEM_SOFT {WAMIT, FAST_WAMIT, WAMIT_1_3, NEMOH, SEAFEM_NEMOH, AQWA, FOAMM, UNKNOWN};
	
	void SaveAs(String file, BEM_SOFT type = UNKNOWN);
	void Report();
	Hydro(BEMData &_bem) : g(Null), h(Null), rho(Null), len(Null), Nb(Null), Nf(Null), Nh(Null), dataFromW(true), bem(&_bem) {}
	virtual ~Hydro() {}	
	
	String GetCodeStr()	const {
		switch (code) {
		case WAMIT: 		return t_("Wamit");
		case WAMIT_1_3: 	return t_("Wamit.1.3");
		case FAST_WAMIT: 	return t_("FAST-Wamit");
		case NEMOH:			return t_("Nemoh");
		case SEAFEM_NEMOH:	return t_("SeaFEM-Nemoh");
		case AQWA:			return t_("AQWA");
		case FOAMM:			return t_("FOAMM");
		case UNKNOWN:		return t_("Unknown");
		}
		return t_("Unknown");
	}
	
	String GetCodeStrAbr() const {
		switch (code) {
		case WAMIT: 		return t_("W.o");
		case WAMIT_1_3: 	return t_("W.1");
		case FAST_WAMIT: 	return t_("FST");
		case NEMOH:			return t_("Nmh");
		case SEAFEM_NEMOH:	return t_("SFM");
		case AQWA:			return t_("AQW");
		case FOAMM:			return t_("FMM");
		case UNKNOWN:		return t_("¿?");
		}
		return t_("Unknown");
	}
		
	inline bool IsAvailableDOF(int ib, int idf) {
		if (dof.IsEmpty())
			return false;
		if (dof[ib] <= idf)
			return false;
		return (Awinf.size() > 0 && !IsNull(Awinf((ib+1)*idf, (ib+1)*idf))) || 
			   (!A.IsEmpty() && A[0].size() > 0 && !IsNull(A[0]((ib+1)*idf, (ib+1)*idf)));
	}

	String file;        	// BEM output file name
	String name;
    double g;           	// gravity
    double h;           	// water depth
   	double rho;        		// density
   	double len;				// Length scale
   	int dimen;				// false if data is dimensionless
    int Nb;          		// number of bodies
    int Nf;          		// number of wave frequencies
    int Nh;          		// number of wave headings
 	
	Upp::Array<MatrixXd> A;					// [Nf](6*Nb, 6*Nb)	Added mass
    MatrixXd Awinf;        					// (6*Nb, 6*Nb) 	Infinite frequency added mass
    MatrixXd Aw0;        					// (6*Nb, 6*Nb)  	Infinite period added mass
    Upp::Array<MatrixXd> B; 				// [Nf](6*Nb, 6*Nb)	Radiation damping
    Vector<double> head;					// [Nh]             Wave headings (deg)
    Vector<String> names;  					// {Nb}             Body names
    Upp::Array<MatrixXd> C;				 	// [Nb](6, 6)		Hydrostatic restoring coefficients:
    MatrixXd cb;          					// (3,Nb)           Centre of buoyancy
    MatrixXd cg;          					// (3,Nb)     		Centre of gravity
    BEM_SOFT code;        					// BEM_SOFT			BEM code 
    Vector<int> dof;      					// [Nb]            	Degrees of freedom for each body 
    Vector<int> dofOrder;					//					Order of DOF
    
    Upp::Array<MatrixXd> Kirf;				// [Nt](6*Nb, 6*Nb)	Radiation impulse response function IRF
    Vector<double> Tirf;	  				// [Nt]				Time-window for the calculation of the IRF
    	
    struct Forces {
    	Upp::Array<MatrixXd> ma, ph;   		// [Nh](Nf, 6*Nb) 	Magnitude and phase
    	Upp::Array<MatrixXd> re, im;		// [Nh](Nf, 6*Nb)	Real and imaginary components
    };
    
    Forces ex; 								// Excitation
    Forces sc;			 					// Diffraction scattering
    Forces fk; 								// Froude-Krylov
    
  	typedef struct Forces RAO;
   
   	RAO rao;
    
    struct StateSpace {
	    Upp::Array<std::complex<double>> TFSResponse;
		Upp::Array<std::complex<double>> Z;
		MatrixXd A_ss;
		VectorXd B_ss;
		VectorXd C_ss;
		VectorXd ssFrequencies, ssFreqRange, ssFrequencies_index;
		double ssMAE;
    };
    Upp::Array<Upp::Array<StateSpace>> sts;	// (6*Nb, 6*Nb)		State space data
    
    Vector<double> T; 						// [Nf]    			Wave periods
    Vector<double> w;     		 			// [Nf]             Wave frequencies
    bool dataFromW;
    Vector<double> Vo;   		 			// [Nb]             Displaced volume
    		
    void Dimensionalize();
    void Normalize();
    
    static String C_units(int i, int j);
    
	int GetHeadId(double hd) const {
		for (int i = 0; i < head.GetCount(); ++i) {
			if (EqualRatio(head[i], hd, 0.01))
				return i;
		}
		return -1;
	}
	
	bool AfterLoad(Function <bool(String, int)> Status);
	void Normalize_Forces(Forces &f);
	void Dimensionalize_Forces(Forces &f);
	void Initialize_Forces();
	void Initialize_Forces(Forces &f);
	void Add_Forces(Forces &to, const Hydro &hydro, const Forces &from);
	void Initialize_RAO();
	void GetFexFromFscFfk();
	void InitializeSts();
		
	static int GetK_AB(int i, int j) {
		while (i > 5)
			i -= 6;
		while (j > 5)
			j -= 6;
		if      ((i == 0 || i == 1 || i == 2) && (j == 0 || j == 1 || j == 2))
			return 3;
		else if ((i == 3 || i == 4 || i == 5) && (j == 3 || j == 4 || j == 5))
			return 5;
		else
			return 4;
	}
	
	static int GetK_F(int i) {
		while (i > 5)
			i -= 6;
		if (i < 3)
			return 2;
		else
			return 3;
	}
	
	static int GetK_C(int i, int j) {
		if (i == 2 && j == 2)
			return 2;	
		else if (i < 3) 
			return 3;
		else
			return 4;
	}
	
	static int GetK_RAO(int i) {
		if (i < 3)
			return 0;	
		else
			return 1;
	}
	
	void GetBodyDOF();
	
	int GetIrregularHead() const;	
	int GetIrregularFreq() const;	
	
	String lastError;
	
	double g_dim()		const;
	double g_ndim()		const;
	double rho_dim()	const;
	double rho_ndim()	const;
	double g_rho_dim()  const;
	double g_rho_ndim() const;
	
	double A_dim(int ifr, int idf, int jdf) 	const {return dimen  ? A[ifr](idf, jdf)*g_rho_dim()/g_rho_ndim() : A[ifr](idf, jdf)*(rho_dim()*pow(len, GetK_AB(idf, jdf)));}
	double A_ndim(int ifr, int idf, int jdf) 	const {return !dimen ? A[ifr](idf, jdf) : A[ifr](idf, jdf)/(rho_ndim()*pow(len, GetK_AB(idf, jdf)));}
	double A_(bool ndim, int ifr, int idf, int jdf) const {return ndim ? A_ndim(ifr, idf, jdf) : A_dim(ifr, idf, jdf);}
	double Aw0_dim(int idf, int jdf)   		 	const {return dimen  ? Aw0(idf, jdf)*g_rho_dim()/g_rho_ndim()    : Aw0(idf, jdf)  *(rho_dim()*pow(len, GetK_AB(idf, jdf)));}
	double Aw0_ndim(int idf, int jdf)  		 	const {return !dimen ? Aw0(idf, jdf)    : Aw0(idf, jdf)  /(rho_ndim()*pow(len, GetK_AB(idf, jdf)));}
	double Aw0_(bool ndim, int idf, int jdf) 	const {return ndim ? Aw0_ndim(idf, jdf) : Aw0_dim(idf, jdf);}
	double Awinf_dim(int idf, int jdf) 		 	const {return dimen  ? Awinf(idf, jdf)*g_rho_dim()/g_rho_ndim()  : Awinf(idf, jdf)*(rho_dim()*pow(len, GetK_AB(idf, jdf)));}
	double Awinf_ndim(int idf, int jdf)		 	const {return !dimen ? Awinf(idf, jdf)  : Awinf(idf, jdf)/(rho_ndim()*pow(len, GetK_AB(idf, jdf)));}
	double Awinf_(bool ndim, int idf, int jdf) 	const {return ndim ? Awinf_ndim(idf, jdf) : Awinf_dim(idf, jdf);}
	
	double B_dim(int ifr, int idf, int jdf)  	   const {return dimen  ? B[ifr](idf, jdf)*g_rho_dim()/g_rho_ndim() : B[ifr](idf, jdf)*(rho_dim()*pow(len, GetK_AB(idf, jdf))*w[ifr]);}
	double B_ndim(int ifr, int idf, int jdf) 	   const {return !dimen ? B[ifr](idf, jdf)*g_rho_ndim()/g_rho_dim() : B[ifr](idf, jdf)/(rho_ndim()*pow(len, GetK_AB(idf, jdf))*w[ifr]);}
	double B_(bool ndim, int ifr, int idf, int jdf)const {return ndim ? B_ndim(ifr, idf, jdf) : B_dim(ifr, idf, jdf);}	
	double C_dim(int ib, int idf, int jdf)   	   const {return dimen  ? C[ib](idf, jdf)*g_rho_dim()/g_rho_ndim()  : C[ib](idf, jdf)*(g_rho_dim()*pow(len, GetK_C(idf, jdf)));}
	double C_ndim(int ib, int idf, int jdf)  	   const {return !dimen ? C[ib](idf, jdf)  : C[ib](idf, jdf)/(g_rho_ndim()*pow(len, GetK_C(idf, jdf)));}
	double C_(bool ndim, int ib, int idf, int jdf) const {return ndim ? C_ndim(ib, idf, jdf) : C_dim(ib, idf, jdf);}

	double F_ma_dim(const Forces &f, int _h, int ifr, int idf)  	   const {return dimen ? f.ma[_h](ifr, idf)*g_rho_dim()/g_rho_ndim()  : f.ma[_h](ifr, idf)*(g_rho_dim()*pow(len, GetK_F(idf)));}
	double F_re_dim(const Forces &f, int _h, int ifr, int idf)  	   const {return dimen ? f.re[_h](ifr, idf)*g_rho_dim()/g_rho_ndim()  : f.re[_h](ifr, idf)*(g_rho_dim()*pow(len, GetK_F(idf)));}
	double F_im_dim(const Forces &f, int _h, int ifr, int idf)  	   const {return dimen ? f.im[_h](ifr, idf)*g_rho_dim()/g_rho_ndim()  : f.im[_h](ifr, idf)*(g_rho_dim()*pow(len, GetK_F(idf)));}
	double F_ma_ndim(const Forces &f, int _h, int ifr, int idf) 	   const {return !dimen ? f.ma[_h](ifr, idf) : f.ma[_h](ifr, idf)/(g_rho_ndim()*pow(len, GetK_F(idf)));}
	double F_re_ndim(const Forces &f, int _h, int ifr, int idf) 	   const {return !dimen ? f.re[_h](ifr, idf) : f.re[_h](ifr, idf)/(g_rho_ndim()*pow(len, GetK_F(idf)));}
	double F_im_ndim(const Forces &f, int _h, int ifr, int idf) 	   const {return !dimen ? f.im[_h](ifr, idf) : f.im[_h](ifr, idf)/(g_rho_ndim()*pow(len, GetK_F(idf)));}
	double F_ma_(bool ndim, const Forces &f, int _h, int ifr, int idf) const {return ndim ? F_ma_ndim(f, _h, ifr, idf) : F_ma_dim(f, _h, ifr, idf);}
	double F_re_(bool ndim, const Forces &f, int _h, int ifr, int idf) const {return ndim ? F_re_ndim(f, _h, ifr, idf) : F_re_dim(f, _h, ifr, idf);}
	double F_im_(bool ndim, const Forces &f, int _h, int ifr, int idf) const {return ndim ? F_im_ndim(f, _h, ifr, idf) : F_im_dim(f, _h, ifr, idf);}
	
	double R_ma_dim(const Forces &f, int _h, int ifr, int idf)  	   const {return dimen ? f.ma[_h](ifr, idf)*g_rho_dim()/g_rho_ndim()  : f.ma[_h](ifr, idf)*(g_rho_dim()*pow(len, GetK_RAO(idf)));}
	double R_re_dim(const Forces &f, int _h, int ifr, int idf)  	   const {return dimen ? f.re[_h](ifr, idf)*g_rho_dim()/g_rho_ndim()  : f.re[_h](ifr, idf)*(g_rho_dim()*pow(len, GetK_RAO(idf)));}
	double R_im_dim(const Forces &f, int _h, int ifr, int idf)  	   const {return dimen ? f.im[_h](ifr, idf)*g_rho_dim()/g_rho_ndim()  : f.im[_h](ifr, idf)*(g_rho_dim()*pow(len, GetK_RAO(idf)));}
	double R_ma_ndim(const Forces &f, int _h, int ifr, int idf) 	   const {return !dimen ? f.ma[_h](ifr, idf) : f.ma[_h](ifr, idf)/(g_rho_ndim()*pow(len, GetK_RAO(idf)));}
	double R_re_ndim(const Forces &f, int _h, int ifr, int idf) 	   const {return !dimen ? f.re[_h](ifr, idf) : f.re[_h](ifr, idf)/(g_rho_ndim()*pow(len, GetK_RAO(idf)));}
	double R_im_ndim(const Forces &f, int _h, int ifr, int idf) 	   const {return !dimen ? f.im[_h](ifr, idf) : f.im[_h](ifr, idf)/(g_rho_ndim()*pow(len, GetK_RAO(idf)));}
	double R_ma_(bool ndim, const Forces &f, int _h, int ifr, int idf) const {return ndim ? R_ma_ndim(f, _h, ifr, idf) : R_ma_dim(f, _h, ifr, idf);}
	double R_re_(bool ndim, const Forces &f, int _h, int ifr, int idf) const {return ndim ? R_re_ndim(f, _h, ifr, idf) : R_re_dim(f, _h, ifr, idf);}
	double R_im_(bool ndim, const Forces &f, int _h, int ifr, int idf) const {return ndim ? R_im_ndim(f, _h, ifr, idf) : R_im_dim(f, _h, ifr, idf);}
	
private:
	static const char *strDOF[6];
	static const char *strDOFAbrev[6];
	static String C_units_base(int i, int j);
	BEMData *bem;
	
public:
	static String StrBDOF(int i) {
		int ib = i/6 + 1;
		int idf = i - (ib - 1)*6;
		return Format("%d.%s", ib, strDOF[idf]);
	}
	
	static String StrBDOFFull(int i) {
		int ib = i/6 + 1;
		int idf = i - (ib - 1)*6;
		return Format("Body #%d. DoF: %s", ib, strDOF[idf]);
	}
	
	static String StrBDOF(int i, int j) {
		if (i != j) {
			int ib = i/6 + 1;
			int idf = i - (ib - 1)*6;
			int jb = j/6 + 1;
			int jdf = j - (jb - 1)*6;
			if (ib != jb)
				return Format("%d.%s_%d.%s", ib, strDOF[idf], jb, strDOF[jdf]);
			else
				return Format("%d.%s_%s", ib, strDOFAbrev[idf], strDOFAbrev[jdf]);
		} else
			return StrBDOF(i);
	}

	static String StrBDOFFull(int i, int j) {
		if (i != j) {
			int ib = i/6 + 1;
			int idf = i - (ib - 1)*6;
			int jb = j/6 + 1;
			int jdf = j - (jb - 1)*6;
			if (ib != jb)
				return Format("Body #%d, DoF: %s. Body #%d, DoF: %s", ib, strDOF[idf], jb, strDOF[jdf]);
			else
				return Format("Body #%d. DoF: %s, DoF: %s", ib, strDOF[idf], strDOF[jdf]);
		} else
			return StrBDOF(i);
	}
		
	static const char *StrDOF_base(int i) {
		return strDOF[i];
	}
	
	static const char *StrDOFAbrev_base(int i) {
		return strDOFAbrev[i];
	}
		
	static String StrBDOFAbrev(int i) {
		int nb = i/6 + 1;
		int ni = i - (nb - 1)*6;
		return Format("%d%s", nb, strDOFAbrev[ni]);
	}
	
	static int DOFStr(const String &str) {
		for (int i = 0; i < 6; ++i)
			if (strDOF[i] == ToLower(str))
				return i;
		return -1;
	}
	
	static void DOFFromStr(const String str, int &ib, int &idf) {
		int pos = str.Find(".");
		ib = ScanInt(str.Left(pos))-1;
		String sdof = str.Mid(pos+1);
		idf = DOFStr(sdof);	
	}
	
	static int DOFStrAbrev(const String &str) {
		for (int i = 0; i < 6; ++i)
			if (strDOFAbrev[i] == ToLower(str))
				return i;
		return -1;
	}
	
	const BEMData &GetBEMData() const {return *bem;}
	
	bool IsLoadedA() 	 const {return A.GetCount() > 0;}
	bool IsLoadedAwinf() const {return Awinf.size() > 0;}
	bool IsLoadedAw0()	 const {return Aw0.size() > 0;}
	bool IsLoadedB() 	 const {return B.GetCount() > 0;}
	bool IsLoadedC()	 const {return C.size() > 0;}
	bool IsLoadedFex() 	 const {return ex.ma.GetCount() > 0;}
	bool IsLoadedFsc() 	 const {return sc.ma.GetCount() > 0;}
	bool IsLoadedFfk() 	 const {return fk.ma.GetCount() > 0;}
	bool IsLoadedRAO() 	 const {return rao.ma.GetCount() > 0;}
	bool IsLoadedForce(Forces &f) const {return f.ma.GetCount() > 0;}
	bool IsLoadedStateSpace()	  const {return sts.GetCount() > 0;}
	
	void RemoveThresDOF_A(double thres);
	void RemoveThresDOF_B(double thres);
	void RemoveThresDOF_Force(Forces &f, double thres);
	
	void Compare_rho(Hydro &a);
	void Compare_g(Hydro &a);
	void Compare_h(Hydro &a);
	void Compare_w(Hydro &a);
	void Compare_head(Hydro &a);
	void Compare_Nb(Hydro &a);
	void Compare_A(Hydro &a);
	void Compare_B(Hydro &a);
	void Compare_C(Hydro &a);
	void Compare_cg(Hydro &a);
	
	//const Vector<int> &GetOrder()		{return dofOrder;}
	void SetOrder(Vector<int> &order)	{dofOrder = pick(order);}
	
	int GetW0();
	void A0();
		
	void K_IRF(double maxT = 120, int numT = 1000);
	void Ainf();
	
	void Join(const Vector<Hydro *> &hydrosp);
	
	String S_g()	const {return IsNull(g)   ? S("unknown") : Format("%.3f", g);}
	String S_h()	const {return IsNull(h)   ? S("unknown") : (h < 0 ? S(t_("INFINITY")) : Format("%.1f", h));}
	String S_rho() 	const {return IsNull(rho) ? S("unknown") : Format("%.3f", rho);}
	String S_len() 	const {return IsNull(len) ? S("unknown") : Format("%.1f", len);}

	String GetLastError()	{return lastError;}
};

		
class HydroData {
public:
	HydroData() : data(0) {}
	HydroData(BEMData &bem, Hydro *_data = 0) {
		if (!_data) {
			manages = true;
			data = new Hydro(bem);
		} else {
			manages = false;
			data = _data;
		}
	}
	Hydro &operator()()				{return *data;}
	const Hydro &operator()() const	{return *data;}
	virtual ~HydroData() {
		if (manages)
			delete data;
	}
	
private:
	Hydro *data;	
	bool manages;
};

class HydroClass {
public:
	HydroClass()							{}
	HydroClass(BEMData &bem, Hydro *hydro = 0) : hd(bem, hydro)	{}
	virtual ~HydroClass()					{}
	
	HydroData hd;	
};

class MeshData {
public:
	enum MESH_FMT {WAMIT_GDF, WAMIT_DAT, NEMOH_DAT, STL_BIN, STL_TXT, UNKNOWN};
	enum MESH_TYPE {ORIGINAL, MOVED, UNDERWATER};
	
	String GetCodeStr()	const {
		switch (code) {
		case WAMIT_GDF: 	return t_("Wamit.gdf");
		case WAMIT_DAT: 	return t_("Wamit.dat");
		case NEMOH_DAT: 	return t_("Nemoh.dat");
		case STL_BIN: 		return t_("Binary.stl");
		case STL_TXT: 		return t_("Text.stl");
		case UNKNOWN:		return t_("Unknown");
		}
		return t_("Unknown");
	}
	void SetCode(MESH_FMT _code){code = _code;}
	MESH_FMT GetCode()			{return code;}

	String Load(String fileName, double rho, double g);
	String LoadDatNemoh(String fileName, bool &x0z);
	String LoadDatWamit(String fileName);
	String LoadGdfWamit(String fileName, bool &y0z, bool &x0z);
	String LoadStlTxt(String fileName, bool &isText);
	String LoadStlBin(String fileName);
	
	String Heal(Function <void(String, int pos)> Status);
		
	void AfterLoad(double rho, double g, bool onlyCG);

	void SaveAs(String fileName, MESH_FMT type, double g, MESH_TYPE meshType);
	static void SaveDatNemoh(String fileName, const Vector<Panel> &panels, const Vector<Point3D> &nodes, bool x0z);
	static void SaveGdfWamit(String fileName, const Vector<Panel> &panels, const Vector<Point3D> &nodes, double g, bool y0z, bool x0z);
	static void SaveStlTxt(String fileName, const Vector<Panel> &panels, const Vector<Point3D> &nodes);
	static void SaveStlBin(String fileName, const Vector<Panel> &panels, const Vector<Point3D> &nodes);
	
	void Report(double rho);
	
	double waterPlaneArea; 
	Point3D cb;
	Point3D cg, cg0;
	double mass;
	MatrixXd c;
	
	String file;
	String header;
	
	Surface mesh, under;
	
private:
	MESH_FMT code;
};

class FileInLine : public FileIn {
public:
	FileInLine(String data) : FileIn(data), line(0) {};
	String GetLine() {
		line++;	
		return FileIn::GetLine();
	}
	void GetLine(int num) {
		for (int i = 0; i < num; ++i)
			GetLine();
	}
	int GetLineNumber()	const {return line;}
	
	struct Pos {
		Pos() : byt(0), line(0) {}
		int64 byt;
		int line;
	};
	
	Pos GetPos() {
		Pos ret;
		ret.byt = FileIn::GetPos();
		ret.line = line;
		return ret;
	}
	
	void SeekPos(Pos &_pos) {
		FileIn::Seek(_pos.byt);
		line = _pos.line;
	}
	
private:
	int line;
};

class FileInData : public FileIn {
public:
	FileInData(const char *fn) : FileIn(fn)	{}
	FileInData()                          	{}
	
	void Read(void *data, int sz) {
		size_t len = Get64(data, sz);
		if (len != static_cast<size_t>(sz))
			throw Exc(Format("Data not loaded in Read(%d)", sz));
	}
	
	template <class T>
	T Read() {
		T data;
		Read(&data, sizeof(T));
		return data;
	}
};

class FileOutData : public FileOut {
public:
	FileOutData(const char *fn) : FileOut(fn)	{}
	FileOutData()                          		{}
	
	template <class T>
	void Write(T data) {
		Put64(&data, sizeof(T));
	}
};

class Wamit : public HydroClass {
public:
	Wamit(BEMData &bem, Hydro *hydro = 0) : HydroClass(bem, hydro) {}
	bool Load(String file);
	void Save(String file);
	virtual ~Wamit()	{}
	
	bool LoadGdfMesh(String file);
	bool LoadDatMesh(String file);
	void SaveGdfMesh(String fileName);
	
protected:
	bool Load_out();							
	void Load_A(FileInLine &in, MatrixXd &A);
	bool Load_Scattering(String fileName);
	bool Load_FK(String fileName);

	bool Load_1(String fileName);				
	bool Load_3(String fileName);
	bool Load_hst(String fileName);
	bool Load_4(String fileName);
	
	void Save_1(String fileName);
	void Save_3(String fileName);
	void Save_hst(String fileName);
	void Save_4(String fileName);
};

class Foamm : public HydroClass {
public:
	Foamm(BEMData &bem, Hydro *hydro = 0) : HydroClass(bem, hydro) {}
	bool Load(String file);
	void Get_Each(int ibody, int idof, int jdof, double from, double to, const Vector<double> &freqs, Function <bool(String, int)> Status, Function <void(String)> FOAMMMessage);
	void Get(const Vector<int> &ibs, const Vector<int> &idofs, const Vector<int> &jdofs,
		const Vector<double> &froms, const Vector<double> &tos, const Vector<Vector<double>> &freqs, 
		Function <bool(String, int)> Status, Function <void(String)> FOAMMMessage);
	virtual ~Foamm()	{}
	
protected:
	bool Load_mat(String fileName, int ib, int jb, bool loadCoeff);
};

class Fast : public Wamit {
public:
	Fast(BEMData &bem, Hydro *hydro = 0) : Wamit(bem, hydro), WaveNDir(Null), WaveDirRange(Null) {}
	bool Load(String file, double g = 9.81);
	void Save(String file);
	virtual ~Fast()	{}
	
private:
	bool Load_HydroDyn();	
	void Save_HydroDyn(String fileName, bool force);	
	
	String hydroFolder;
	int WaveNDir;
	double WaveDirRange;
};

class NemohBody : public Moveable<NemohBody> {
public:
	NemohBody() {surge = sway = heave = roll = pitch = yaw = false;}
	NemohBody(const NemohBody &d) : meshFile(d.meshFile), 
		surge(d.surge), sway(d.sway), heave(d.heave), roll(d.roll), pitch(d.pitch), yaw(d.yaw),
		cx(d.cx), cy(d.cy), cz(d.cz), npoints(d.npoints), npanels(d.npanels) {}
	
	String meshFile;
	bool surge, sway, heave, roll, pitch, yaw;
	double cx, cy, cz;	
	int npoints, npanels;
	
	int GetNDOF() const;
};

class NemohCal {
public:	
	double rho = Null, g = Null, h = Null;
	double xeff, yeff;
	Vector<NemohBody> bodies;
	int Nf = Null;
	double minF, maxF;
	int Nh = Null;
	double minH, maxH;
	bool irf;
	double irfStep, irfDuration;	
	bool showPressure;
	int nFreeX, nFreeY;
	double domainX, domainY;
	int nKochin;
	double minK, maxK;
	
	bool Load(String fileName);
	void Save_Cal(String folder, int _nf, double _minf, double _maxf) const;
	void SaveFolder(String folder, bool bin, int numCases, const BEMData &bem) const;
	Vector<String> Check();
	
private:
	static int GetNumArgs(const FieldSplit &f);
	void LoadFreeSurface(const FileInLine &in, const FieldSplit &f);
	void LoadKochin(const FileInLine &in, const FieldSplit &f);

	void CreateId(String folder) const;
	void CreateBat(String folder, String batname, String caseFolder, bool bin, String preName, String solvName, String postName) const;
	void CreateInput(String folder) const;
};

class Nemoh : public HydroClass {
public:
	Nemoh(BEMData &bem, Hydro *hydro = 0) : HydroClass(bem, hydro) {}
	bool Load(String file, double rho = Null);
	void Save(String file);
	virtual ~Nemoh()	{}
	
	bool LoadDatMesh(String file);
	void SaveDatMesh(String file); 
	
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
	bool Load_Forces(Hydro::Forces &f, String nfolder, String fileName, String textDelim);
	bool Load_IRF(String fileName);
};

class Aqwa : public HydroClass {
public:
	Aqwa(BEMData &bem, Hydro *hydro = 0) : HydroClass(bem, hydro) {}
	bool Load(String file, double rho = Null);
	void Save(String file);
	virtual ~Aqwa()	{}
	
private:
	bool Load_AH1();
	bool Load_LIS();
};

template <class Range, class T>
void LinSpaced(Range &v, int n, T min, T max) {
	ASSERT(n > 0);
	v.SetCount(n);
	if (n == 1)
		v[0] = min;
	else {
		for (int i = 0; i < n; ++i)
			v[i] = min + ((max - min)*i)/(n - 1);
	}
}

int IsTabSpace(int c);

class FieldSplit {
public:
	const int FIRST = 0;
	const int LAST = Null;
	
	FieldSplit(FileInLine &_in) {in = &_in;}
	
	FieldSplit& Load(String _line) {
		line = _line;
		fields = Split(line, IsTabSpace, true);
		return *this;
	}
	void LoadWamitJoinedFields(String _line) {		// Trick for "glued" fields in Wamit
		line = _line;
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
	String GetText() const {
		return line;
	}
	String GetText(int i) const {
		if (fields.IsEmpty())
			throw Exc(Format(t_("[line %d] No data available"), in->GetLineNumber()));
		if (IsNull(i))
			i = fields.GetCount()-1;
		CheckId(i);
		return fields[i];
	}
	int GetInt(int i) const {
		if (fields.IsEmpty())
			throw Exc(Format(t_("[line %d] No data available"), in->GetLineNumber()));
		if (IsNull(i))
			i = fields.GetCount()-1;
		CheckId(i);
		int res = ScanInt(fields[i]);
		if (IsNull(res))
			throw Exc(Format(t_("[line %d] Bad %s '%s' in field #%d, line\n'%s'"), in->GetLineNumber(), "integer", fields[i], i+1, line));
		return res;
	}
	double GetDouble(int i) const {
		if (fields.IsEmpty())
			throw Exc(Format(t_("[line %d] No data available"), in->GetLineNumber()));
		if (IsNull(i))
			i = fields.GetCount()-1;
		CheckId(i);
		double res = ScanDouble(fields[i]);
		if (IsNull(res))
			throw Exc(Format(t_("[line %d] Bad %s '%s' in field #%d, line\n'%s'"), in->GetLineNumber(), "double", fields[i], i+1, line));
		return res;
	}
	int GetCount() const {
		return fields.GetCount();
	}
	
private:
	String line;
	Vector<String> fields;
	FileInLine *in;
	
	void CheckId(int i) const {
		if (i >= fields.GetCount() || i < 0)
			throw Exc(Format(t_("[line %d] Field #%d not found in line\n'%s'"), in->GetLineNumber(), i+1, line));
	}
};

String FormatWam(double d);

class BEMData {
public:
	BEMData();
	
	Upp::Array<HydroClass> hydros;
	Upp::Array<MeshData> surfs;
	
	static Function <void(String)> Print, PrintWarning, PrintError;	
	
	Vector<double> head;	// Common models data
	int Nb = 0;				
	
	double depth, rho, g, length;
	int discardNegDOF;
	double thres;
	int calcAwinf;
	double maxTimeA;
	int numValsA;
	int onlyDiagonal;
	
	String nemohPathPreprocessor, nemohPathSolver, nemohPathPostprocessor, nemohPathNew, nemohPathGREN;
	bool experimental, experimentalFOAMM;
	String foammPath;
	
	void Load(String file, Function <bool(String, int pos)> Status);
	void Join(Vector<int> &ids, Function <bool(String, int)> Status);
	void LoadMesh(String file, Function <void(String, int pos)> Status);
	void HealingMesh(int id, Function <void(String, int pos)> Status);
	void UnderwaterMesh(int id, Function <void(String, int pos)> Status);
			
	bool LoadSerializeJson();
	bool StoreSerializeJson();
	bool ClearTempFiles();
	static String GetTempFilesFolder() {return AppendFileName(GetAppDataFolder(), "BEMRosetta", "Temp");}
	
	const String bemFilesExt = ".1 .3 .hst .4 .out .cal .inf .ah1 .lis .mat .dat";
	String bemFilesAst;
	
	void Jsonize(JsonIO &json) {
		json
			("depth", depth)
			("rho", rho)
			("g", g)
			("length", length)
			("discardNegDOF", discardNegDOF)
			("thres", thres)
			("calcAwinf", calcAwinf)
			("maxTimeA", maxTimeA)
			("numValsA", numValsA)
			("onlyDiagonal", onlyDiagonal)
			("nemohPathPreprocessor", nemohPathPreprocessor)
			("nemohPathSolver", nemohPathSolver)
			("nemohPathPostprocessor", nemohPathPostprocessor)
			("nemohPathGREN", nemohPathGREN)
			("nemohPathNew", nemohPathNew)
			("experimental", experimental)
			("experimentalFOAMM", experimentalFOAMM)
			("foammPath", foammPath)
		;
	}
};

template <class T>
bool OUTB(int id, T total) {
	if (id < 0	|| id >= static_cast<int>(total))
		return true;
	return false;
}

#endif
