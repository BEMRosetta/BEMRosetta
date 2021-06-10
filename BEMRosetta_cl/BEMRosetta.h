#ifndef _BEM_Rosetta_BEM_Rosetta_h_
#define _BEM_Rosetta_BEM_Rosetta_h_

#include <Core/Core.h>
#include <SysInfo/SysInfo.h>
#include <Surface/Surface.h>

using namespace Upp;


class BEMData;

void ConsoleMain(const Upp::Vector<String>& command, bool gui);
void SetBuildInfo(String &str);

enum {SURGE = 0, SWAY, HEAVE, ROLL, PITCH, YAW};

class Hydro {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	enum BEM_SOFT {WAMIT, FAST_WAMIT, WAMIT_1_3, NEMOH, SEAFEM_NEMOH, AQWA, FOAMM, BEMROSETTA, UNKNOWN};
	
	void SaveAs(String file, BEM_SOFT type = UNKNOWN, int qtfHeading = Null);
	void Report();
	Hydro(BEMData &_bem) : g(Null), h(Null), rho(Null), len(Null), Nb(Null), Nf(Null), Nh(Null), 
							dataFromW(Null), bem(&_bem) {id = idCount++;}
	virtual ~Hydro() noexcept {}	
	
	String GetCodeStr()	const {
		switch (code) {
		case WAMIT: 		return t_("Wamit");
		case WAMIT_1_3: 	return t_("Wamit.1.3");
		case FAST_WAMIT: 	return t_("FAST-Wamit");
		case NEMOH:			return t_("Nemoh");
		case SEAFEM_NEMOH:	return t_("SeaFEM-Nemoh");
		case AQWA:			return t_("AQWA");
		case FOAMM:			return t_("FOAMM");
		case BEMROSETTA:	return t_("BEMRosetta");
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
		case BEMROSETTA:	return t_("BMR");
		case UNKNOWN:		return t_("¿?");
		}
		return t_("Unknown");
	}
		
	inline bool IsAvailableDOF(int ib, int idf) {
		if (dof.IsEmpty())
			return false;
		//if (dof[ib] <= idf)
		//	return false;
		
		int i = ib*6+idf;
		
		if (Awinf.size() > 0)
			if (!IsNull(Awinf(i, i)))
				return true;
		
		if (!A.IsEmpty())
			if (A[i][i].size() > 0)
				if (!IsNull(A[i][i][0]))
					return true;
		
		return false;			   
		//return (Awinf.size() > 0 && !IsNull(Awinf((ib+1)*idf, (ib+1)*idf))) || 
		//	   (!A.IsEmpty() && A[0].size() > 0 && !IsNull(A[(ib+1)*idf][(ib+1)*idf][0]));
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
 	
	//Upp::Array<Eigen::MatrixXd> A;			// [Nf](6*Nb, 6*Nb)	Added mass
	Upp::Array<Upp::Array<Eigen::VectorXd>> A;	// [6*Nb][6*Nb][Nf]	Added mass
    Eigen::MatrixXd Awinf;        			// (6*Nb, 6*Nb) 	Infinite frequency added mass
    Eigen::MatrixXd Aw0;        			// (6*Nb, 6*Nb)  	Infinite period added mass
    //Upp::Array<Eigen::MatrixXd> B; 			// [Nf](6*Nb, 6*Nb)	Radiation damping
    Upp::Array<Upp::Array<Eigen::VectorXd>> B; 	// [6*Nb][6*Nb][Nf]	Radiation damping
    Upp::Vector<double> head;				// [Nh]             Wave headings (deg)
    Upp::Vector<String> names;  			// {Nb}             Body names
    Upp::Array<Eigen::MatrixXd> C;			// [Nb](6, 6)		Hydrostatic restoring coefficients:
    Upp::Array<Eigen::MatrixXd> M;			// [Nb](6, 6)		Mass and inertia matrix
    Eigen::MatrixXd cb;          			// (3,Nb)           Centre of buoyancy
    Eigen::MatrixXd cg;          			// (3,Nb)     		Centre of gravity
    BEM_SOFT code;        					// BEM_SOFT			BEM code 
    Upp::Vector<int> dof;      				// [Nb]            	Degrees of freedom for each body 
    Upp::Vector<int> dofOrder;				// [6*Nb]			DOF order
    
    //Upp::Array<Eigen::MatrixXd> Kirf;		// [Nt](6*Nb, 6*Nb)	Radiation impulse response function IRF
    Upp::Array<Upp::Array<Eigen::VectorXd>> Kirf;	// [6*Nb][6*Nb][Nt]	Radiation impulse response function IRF
    Eigen::VectorXd Tirf;	  				// [Nt]				Time-window for the calculation of the IRF
    
    int GetHeadId(double hd) const;
	
    struct Forces {
    	Upp::Array<Eigen::MatrixXd> ma, ph;	// [Nh](Nf, 6*Nb) 	Magnitude and phase
    	Upp::Array<Eigen::MatrixXd> re, im;	// [Nh](Nf, 6*Nb)	Real and imaginary components
    	void Jsonize(JsonIO &json) {
			json
				("ma", ma)
				("ph", ph)
				("re", re)
				("im", im)
			;
    	}
    };
    
    Forces ex; 								// Excitation
    Forces sc;			 					// Diffraction scattering
    Forces fk; 								// Froude-Krylov
    
  	typedef struct Forces RAO;
   
   	RAO rao;
    
    String description;

    struct StateSpace {
	    Upp::Array<std::complex<double>> TFS;
		Eigen::MatrixXd A_ss;
		Eigen::VectorXd B_ss;
		Eigen::VectorXd C_ss;
		Eigen::VectorXd ssFrequencies, ssFreqRange, ssFrequencies_index;
		double ssMAE = Null;
		
		void GetTFS(const Upp::Vector<double> &w);
		
		void Jsonize(JsonIO &json) {
			json
				("TFSResponse", TFS)
				("A_ss", A_ss)
				("B_ss", B_ss)
				("C_ss", C_ss)
				("ssFrequencies", ssFrequencies)
				("ssFreqRange", ssFreqRange)
				("ssFrequencies_index", ssFrequencies_index)
				("ssMAE", ssMAE)
			;
    	}
    };
    Upp::Array<Upp::Array<StateSpace>> sts;	// (6*Nb, 6*Nb)		State space data
    int dimenSTS;							// false if data is dimensionless
    String stsProcessor;
    
    struct QTF {
        QTF() {
            fre.SetCount(6, Null);
            fim.SetCount(6, Null);
            fma.SetCount(6, Null);
            fph.SetCount(6, Null);
        }
        void Set(int _ib, int _ih1, int _ih2, int _ifr1, int _ifr2) {
            ib = _ib;
            ih1 = _ih1;
            ih2 = _ih2;
            ifr1 = _ifr1;
            ifr2 = _ifr2;
        }
        
        int ib = -1;
        int ih1, ih2;
        int ifr1, ifr2;
        Upp::Vector<double> fre, fim, fma, fph;
 
		void Jsonize(JsonIO &json) {
			json
				("ib", ib)
				("ih1", ih1)
				("ih2", ih2)
				("ifr1", ifr1)
				("ifr2", ifr2)
				("fre", fre)
				("fim", fim)
				("fma", fma)
				("fph", fph)
			;
    	}
    };
    Upp::Array<QTF> qtfsum, qtfdif;
    Upp::Vector<double> qtfw, qtfT, qtfhead;
    bool qtfdataFromW;
    
    int GetQTFHeadId(double hd) const;
    static int GetQTFId(const Upp::Array<Hydro::QTF> &qtfList, int _ib, int _ih1, int _ih2, int _ifr1, int _ifr2);
	static void GetQTFList(const Upp::Array<Hydro::QTF> &qtfList, Upp::Vector<int> &ibL, Upp::Vector<int> &ih1L, Upp::Vector<int> &ih2L);
							
    Upp::Vector<double> T;					// [Nf]    			Wave periods
    Upp::Vector<double> w;		 			// [Nf]             Wave frequencies
    bool dataFromW;
    Upp::Vector<double> Vo;   				// [Nb]             Displaced volume
    		
    void Dimensionalize();
    void Normalize();
    
    static String C_units(int i, int j);
    
    void SetC(int ib, const Eigen::MatrixXd &K);
	
	bool AfterLoad(Function <bool(String, int)> Status);
	
	void Initialize_Forces();
	void Initialize_Forces(Forces &f, int _Nh = -1);
	void Normalize_Forces(Forces &f);
	void Dimensionalize_Forces(Forces &f);
	void Add_Forces(Forces &to, const Hydro &hydro, const Forces &from);
	void Symmetrize_Forces();
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
	
	Eigen::VectorXd Get_w() 					const {return Eigen::Map<const Eigen::VectorXd>(w, w.size());}
	Eigen::VectorXd Get_T() 					const {return Eigen::Map<const Eigen::VectorXd>(T, T.size());}
	
	double A_dim(int ifr, int idf, int jdf) 	const {return dimen  ? A[idf][jdf][ifr]*rho_dim()/rho_ndim() : A[idf][jdf][ifr]*(rho_dim()*pow(len, GetK_AB(idf, jdf)));}
	Eigen::VectorXd A_dim(int idf, int jdf) 	const {return dimen  ? A[idf][jdf]*(rho_dim()/rho_ndim()) : A[idf][jdf]*(rho_dim()*pow(len, GetK_AB(idf, jdf)));}
	double A_ndim(int ifr, int idf, int jdf) 	const {return !dimen ? A[idf][jdf][ifr]*(rho_ndim()/rho_dim()) : A[idf][jdf][ifr]/(rho_ndim()*pow(len, GetK_AB(idf, jdf)));}
	Eigen::VectorXd A_ndim(int idf, int jdf)	const {return !dimen ? A[idf][jdf]*(rho_ndim()/rho_dim()) : A[idf][jdf]*(1/(rho_ndim()*pow(len, GetK_AB(idf, jdf))));}
	double A_(bool ndim, int ifr, int idf, int jdf) const {return ndim ? A_ndim(ifr, idf, jdf) : A_dim(ifr, idf, jdf);}
	double Aw0_dim(int idf, int jdf)   		 	const {return dimen  ? Aw0(idf, jdf)*rho_dim()/rho_ndim() : Aw0(idf, jdf)  *(rho_dim()*pow(len, GetK_AB(idf, jdf)));}
	double Aw0_ndim(int idf, int jdf)  		 	const {return !dimen ? Aw0(idf, jdf)    : Aw0(idf, jdf)  /(rho_ndim()*pow(len, GetK_AB(idf, jdf)));}
	double Aw0_(bool ndim, int idf, int jdf) 	const {return ndim ? Aw0_ndim(idf, jdf) : Aw0_dim(idf, jdf);}
	double Awinf_dim(int idf, int jdf) 		 	const {return dimen  ? Awinf(idf, jdf)*rho_dim()/rho_ndim() : Awinf(idf, jdf)*(rho_dim()*pow(len, GetK_AB(idf, jdf)));}
	double Awinf_ndim(int idf, int jdf)		 	const {return !dimen ? Awinf(idf, jdf) : Awinf(idf, jdf)/(rho_ndim()*pow(len, GetK_AB(idf, jdf)));}
	double Awinf_(bool ndim, int idf, int jdf) 	const {return ndim ? Awinf_ndim(idf, jdf) : Awinf_dim(idf, jdf);}
	
	double B_dim(int ifr, int idf, int jdf)  	   const {return dimen  ? B[idf][jdf][ifr]*rho_dim()/rho_ndim() : B[idf][jdf][ifr]*(rho_dim()*pow(len, GetK_AB(idf, jdf))*w[ifr]);}
	Eigen::VectorXd B_dim(int idf, int jdf)  	   const;
	double B_ndim(int ifr, int idf, int jdf) 	   const {return !dimen ? B[idf][jdf][ifr]*(rho_ndim()/rho_dim()) : B[idf][jdf][ifr]/(rho_ndim()*pow(len, GetK_AB(idf, jdf))*w[ifr]);}
	Eigen::VectorXd B_ndim(int idf, int jdf) 	   const;
	double B_(bool ndim, int ifr, int idf, int jdf)const {return ndim ? B_ndim(ifr, idf, jdf) : B_dim(ifr, idf, jdf);}	
	
	double Kirf_dim(int ifr, int idf, int jdf)  	   const {return dimen ? Kirf[idf][jdf][ifr]*g_rho_dim()/g_rho_ndim()  : Kirf[idf][jdf][ifr]*(g_rho_dim()*pow(len, GetK_F(idf)));}
	double Kirf_ndim(int ifr, int idf, int jdf) 	   const {return !dimen ? Kirf[idf][jdf][ifr] : Kirf[idf][jdf][ifr]/(g_rho_ndim()*pow(len, GetK_F(idf)));}
	double Kirf_(bool ndim, int ifr, int idf, int jdf) const {return ndim ? Kirf_ndim(ifr, idf, jdf) : Kirf_dim(ifr, idf, jdf);}
	
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

	double F_dim(double f, int idf)  	   const {return dimen ? f*g_rho_dim()/g_rho_ndim() : f*(g_rho_dim()*pow(len, GetK_F(idf)));}
	double F_ndim(double f, int idf) 	   const {return !dimen ? f : f/(g_rho_ndim()*pow(len, GetK_F(idf)));}
	double F_(bool ndim, double f, int idf) const {return ndim ? F_ndim(f, idf) : F_dim(f, idf);}

	inline std::complex<double>Z(bool ndim, int ifr, int idf, int jdf) const {
		return std::complex<double>(B_(ndim, ifr, idf, jdf), w[ifr]*(A_(ndim, ifr, idf, jdf) - Awinf_(ndim, idf, jdf))/(!ndim ? 1. : w[ifr]));
	}
	
	std::complex<double> TFS_dim(int ifr, int idf, int jdf) 		const {return dimenSTS  ? sts[idf][jdf].TFS[ifr]*g_rho_dim()/g_rho_ndim() : sts[idf][jdf].TFS[ifr]*(rho_dim()*pow(len, GetK_AB(idf, jdf))*w[ifr]);}
	std::complex<double> TFS_ndim(int ifr, int idf, int jdf) 		const {return !dimenSTS ? sts[idf][jdf].TFS[ifr]*g_rho_ndim()/g_rho_dim() : sts[idf][jdf].TFS[ifr]/(rho_ndim()*pow(len, GetK_AB(idf, jdf))*w[ifr]);}
	std::complex<double> TFS_(bool ndim, int ifr, int idf, int jdf) const {return ndim ? TFS_ndim(ifr, idf, jdf) : TFS_dim(ifr, idf, jdf);}
	
	void SetId(int _id)			{id = _id;}
	int GetId()	const			{return id;}
	
	void Jsonize(JsonIO &json);
	
private:
	static const char *strDOF[6];
	static const char *strDOFAbrev[6];
	static String C_units_base(int i, int j);
	BEMData *bem;
		
	void Symmetrize_Forces_Each0(const Forces &f, Forces &newf, const Upp::Vector<double> &newHead, double h, int ih, int idb);
	void Symmetrize_ForcesEach(const Forces &f, Forces &newf, const Upp::Vector<double> &newHead, int newNh);
	int id;
	static int idCount;
	 
	static void GetOldAB(const Upp::Array<Eigen::MatrixXd> &oldAB, Upp::Array<Upp::Array<Eigen::VectorXd>> &AB);
	static void SetOldAB(Upp::Array<Eigen::MatrixXd> &oldAB, const Upp::Array<Upp::Array<Eigen::VectorXd>> &AB);
	
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
	
	bool IsLoadedA() 	 const {return !A.IsEmpty();}
	bool IsLoadedAwinf() const {return Awinf.size() > 0;}
	bool IsLoadedAw0()	 const {return Aw0.size() > 0;}
	bool IsLoadedB() 	 const {return !B.IsEmpty();}
	bool IsLoadedC()	 const {return !C.IsEmpty() && C[0].size() > 0 && !IsNull(C[0](0, 0));}
	bool IsLoadedFex() 	 const {return !ex.ma.IsEmpty();}
	bool IsLoadedFsc() 	 const {return !sc.ma.IsEmpty();}
	bool IsLoadedFfk() 	 const {return !fk.ma.IsEmpty();}
	bool IsLoadedRAO() 	 const {return !rao.ma.IsEmpty();}
	bool IsLoadedForce(const Forces &f) const {return !f.ma.IsEmpty();}
	bool IsLoadedStateSpace()	  const {return !sts.IsEmpty();}
	bool IsLoadedQTF() 	 const {return !qtfsum.IsEmpty();}
	
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
	
	const Vector<int> &GetOrder() const		{return dofOrder;}
	void SetOrder(Upp::Vector<int> &order)	{dofOrder = pick(order);}
	
	int GetW0();
	void Get3W0(int &id1, int &id2, int &id3);
	void A0();
		
	void K_IRF(double maxT = 120, int numT = 1000);
	double GetK_IRF_MaxT();
	void Ainf();
	
	void Join(const Upp::Vector<Hydro *> &hydrosp);
	
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
	virtual ~HydroData() noexcept {
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
	virtual ~HydroClass() noexcept			{}
	bool Load(String file);
	bool Save(String file);
	
	HydroData hd;	
};

class MeshData {
public:
	enum MESH_FMT {WAMIT_GDF, WAMIT_DAT, NEMOH_DAT, NEMOH_PRE, AQWA_DAT, HAMS_PNL, STL_BIN, STL_TXT, EDIT, UNKNOWN};
	enum MESH_TYPE {MOVED, UNDERWATER};
	
	MeshData() {
		id = idCount++;
		cg = Point3D(0, 0, 0);
	}
	
	String GetCodeStr()	const {
		switch (code) {
		case WAMIT_GDF: 	return t_("Wamit.gdf");
		case WAMIT_DAT: 	return t_("Wamit.dat");
		case NEMOH_DAT: 	return t_("Nemoh.dat");
		case NEMOH_PRE:		return t_("Nemoh premesh.");
		case AQWA_DAT: 		return t_("AQWA.dat");
		case HAMS_PNL: 		return t_("HAMS.pnl");
		case STL_BIN: 		return t_("Binary.stl");
		case STL_TXT: 		return t_("Text.stl");
		case EDIT: 			return t_("Edit");
		case UNKNOWN:		return t_("Unknown");
		}
		return t_("Unknown");
	}
	void SetCode(MESH_FMT _code){code = _code;}
	MESH_FMT GetCode()			{return code;}
	int GetId()	const			{return id;}

	String Load(String fileName, double rho, double g, bool cleanPanels);
	String Load(String fileName, double rho, double g, bool cleanPanels, bool &y0z, bool &x0z);
	String Load(String fileName, bool &y0z, bool &x0z);
	String Load(String fileName);
	String LoadDatNemoh(String fileName, bool &x0z);
	String LoadDatWamit(String fileName);
	String LoadGdfWamit(String fileName, bool &y0z, bool &x0z);
	String LoadDatAQWA(String fileName);
	String LoadPnlHAMS(String fileName, bool &y0z, bool &x0z);
	
	String Heal(bool basic, Function <void(String, int pos)> Status);
	void Orient();
	void Join(const Surface &orig, double rho, double g);
	void Image(int axis);
		
	void AfterLoad(double rho, double g, bool onlyCG);

	void SaveAs(String fileName, MESH_FMT type, double g, MESH_TYPE meshType, bool symX, bool symY);
	static void SaveDatNemoh(String fileName, const Surface &surf, bool x0z);
	static void SavePreMeshNemoh(String fileName, const Surface &surf);
	static void SaveGdfWamit(String fileName, const Surface &surf, double g, bool y0z, bool x0z);
	static void SavePnlHAMS(String fileName, const Surface &surf, bool y0z, bool x0z);
	
	void SaveHST(String fileName, double rho, double g) const; 
	
	void Report(double rho);
	
	bool IsSymmetricX();
	bool IsSymmetricY();
	
	double waterPlaneArea; 
	Point3D cb;
	Point3D cg;
	double mass = Null;
	Eigen::MatrixXd C;
	
	String name;
	String fileName;
	String header;
	
	Surface mesh, under;
	
private:
	MESH_FMT code;
	int id;
	static int idCount;
};

class Wamit : public HydroClass {
public:
	Wamit(BEMData &bem, Hydro *hydro = 0) : HydroClass(bem, hydro) {}
	bool Load(String file);
	void Save(String file, bool force_T = false, int qtfHeading = Null);
	void Save_out(String file, double g, double rho);
	virtual ~Wamit() noexcept {}
	
	bool LoadGdfMesh(String file);
	bool LoadDatMesh(String file);
	void SaveGdfMesh(String fileName);
	
	static void Save_hst_static(Eigen::MatrixXd C, String fileName, double rho, double g);
	
protected:
	void ProcessFirstColumn(Vector<double> &w, Vector<double> &T);
	
	bool Load_cfg(String fileName);
	int iperout = Null;
	bool Load_pot(String fileName);
	bool Load_gdf(String fileName);
	
	bool Load_out();							
	void Load_A(FileInLine &in, Eigen::MatrixXd &A);
	bool Load_Scattering(String fileName);
	bool Load_FK(String fileName);

	bool Load_1(String fileName);				
	bool Load_3(String fileName);
	bool Load_hst(String fileName);
	bool Load_4(String fileName);
	bool Load_12(String fileName, bool isSum);
	
	void Save_1(String fileName, bool force_T = false);
	void Save_3(String fileName, bool force_T = false);
	void Save_hst(String fileName);
	void Save_4(String fileName, bool force_T = false);
	void Save_12(String fileName, bool isSum, bool force_T = false, bool force_Deg = true, int qtfHeading = Null);

	void Save_A(FileOut &out, Function <double(int, int)> fun, const Eigen::MatrixXd &base, String wavePeriod);
	void Save_AB(FileOut &out, int ifr);
	void Save_Forces(FileOut &out, int ifr);
	void Save_RAO(FileOut &out, int ifr);
};

class Foamm : public HydroClass {
public:
	Foamm(BEMData &bem, Hydro *hydro = 0) : HydroClass(bem, hydro) {}
	bool Load(String file);
	void Get_Each(int ibody, int idf, int jdf, double from, double to, const Upp::Vector<double> &freqs, Function <bool(String, int)> Status, Function <void(String)> FOAMMMessage);
	void Get(const Upp::Vector<int> &ibs, const Upp::Vector<int> &idfs, const Upp::Vector<int> &jdfs,
		const Upp::Vector<double> &froms, const Upp::Vector<double> &tos, const Upp::Vector<Upp::Vector<double>> &freqs, 
		Function <bool(String, int)> Status, Function <void(String)> FOAMMMessage);
	virtual ~Foamm() noexcept {}
	
protected:
	bool Load_mat(String fileName, int ib, int jb, bool loadCoeff);
};

class Fast : public Wamit {
public:
	Fast(BEMData &bem, Hydro *hydro = 0) : Wamit(bem, hydro), WaveNDir(Null), WaveDirRange(Null) {}
	bool Load(String file, double g = 9.81);
	void Save(String file, int qtfHeading = Null);
	virtual ~Fast() noexcept {}
	
private:
	bool Load_HydroDyn();	
	void Save_HydroDyn(String fileName, bool force);
	bool Load_SS(String fileName);	
	void Save_SS(String fileName);
	
	String hydroFolder;
	int WaveNDir;
	double WaveDirRange;
};

class BemBody : public Moveable<BemBody> {
public:
	BemBody();
	BemBody(const BemBody &d) : meshFile(d.meshFile), dof(clone(d.dof)), c0(clone(d.c0)), 
			cg(clone(d.cg)), mass(clone(d.mass)), 
			linearDamping(clone(d.linearDamping)), quadraticDamping(clone(d.quadraticDamping)), 
			hydrostaticRestoring(clone(d.hydrostaticRestoring)), externalRestoring(clone(d.externalRestoring)), 
			npoints(d.npoints), npanels(d.npanels) {}
	
	String meshFile, lidFile;
	int ndof;
	Vector<bool> dof;
	Eigen::Vector3d c0;	
	
	int npoints, npanels;
	
	Eigen::Vector3d cg;
	Eigen::MatrixXd mass, linearDamping, quadraticDamping, hydrostaticRestoring, externalRestoring;
	
	int GetNDOF() const;
};


class BemCal {
public:
	BemCal() {bodies.SetCount(1);}
		
	void BeforeSave(String folderBase, int numCases, bool deleteFolder) const;
	bool IsDof(int ib, int idf) {return bodies[ib].dof[idf];}
	Upp::Vector<String> Check(bool oneBody);
	
	double h = Null;
	
	Upp::Vector<BemBody> bodies;

	int Nf = Null;
	double minF, maxF;
	int Nh = Null;
	double minH, maxH;
	
	enum Solver {NEMOH, NEMOHv115, CAPYTAINE, HAMS};
};

class HamsCal : public BemCal {
public:
	void SaveFolder(String folderBase, bool bin, int numCases, int numThreads, const BEMData &bem) const;
	bool Load(String fileName);
	bool LoadHydrostatic(String fileName);
		
   	Upp::Vector<String> Check();
   	
private:
	void SaveFolder0(String folderBase, bool bin, int numCases, const BEMData &bem, bool deleteFolder, int numThreads) const;
	static void OutMatrix(FileOut &out, String header, const Eigen::MatrixXd &mat);
	static void InMatrix(FieldSplit &f, Eigen::MatrixXd &mat);
		
	void Save_Hydrostatic(String folderInput) const;
	void Save_ControlFile(String folderInput, int _nf, double _minf, double _maxf,
							int numThreads) const;
	void Save_Bat(String folder, String batname, String caseFolder, bool bin, String solvName) const;
};

class NemohCal : public BemCal {
public:	
	double rho = Null, g = Null;
	double xeff, yeff;
	
	bool irf;
	double irfStep, irfDuration;	
	bool showPressure;
	int nFreeX, nFreeY;
	double domainX, domainY;
	int nKochin;
	double minK, maxK;
	
	bool Load(String fileName);
	void Save_Cal(String folder, int _nf, double _minf, double _maxf) const;
	void SaveFolder(String folder, bool bin, int numCases, const BEMData &bem, 
					int solver) const;
	Upp::Vector<String> Check();
	
private:
	static int GetNumArgs(const FieldSplit &f);
	void LoadFreeSurface(const FileInLine &in, const FieldSplit &f);
	void LoadKochin(const FileInLine &in, const FieldSplit &f);

	void Save_Id(String folder) const;
	void Save_Bat(String folder, String batname, String caseFolder, bool bin, String preName, String solvName, String postName) const;
	void Save_Input(String folder) const;
	
	void SaveFolder0(String folder, bool bin, int numCases, const BEMData &bem, 
					bool deleteFolder, int solver) const;
};

class Nemoh : public HydroClass {
public:
	Nemoh(BEMData &bem, Hydro *hydro = 0) : HydroClass(bem, hydro) {}
	bool Load(String file, double rho = Null);
	void Save(String file);
	virtual ~Nemoh() noexcept {}
	
	bool LoadDatMesh(String file);
	void SaveDatMesh(String file); 
	
private:
	String folder;
	NemohCal datacal;
	
	bool Load_Cal(String fileName);
	bool Load_Inf(String fileName);
	bool Load_Hydrostatics();
	bool Load_KH();
	bool Load_Radiation(String fileName);
	bool Load_Excitation(String folder);
	bool Load_Diffraction(String folder);
	bool Load_FroudeKrylov(String folder);
	bool Load_Forces(Hydro::Forces &f, String nfolder, String fileName);
	bool Load_IRF(String fileName);
};

class Aqwa : public HydroClass {
public:
	Aqwa(BEMData &bem, Hydro *hydro = 0) : HydroClass(bem, hydro) {}
	bool Load(String file, double rho = Null);
	void Save(String file);
	virtual ~Aqwa() noexcept {}
	
private:
	bool Load_AH1();
	bool Load_LIS();
	bool Load_QTF();
};


Upp::Vector<int> NumSets(int num, int numsets);	


class FieldSplitWamit: public FieldSplit {
public:
	FieldSplitWamit(FileInLine &_in) : FieldSplit(_in) {}
	
	void LoadWamitJoinedFields(String _line);// Trick for "glued" fields in Wamit
};

String FormatWam(double d);

class BEMData {
public:
	BEMData();
	
	Upp::Array<HydroClass> hydros;
	Upp::Array<MeshData> surfs;
	
	int GetHydroId(int id) {
		for (int i = 0; i < hydros.size(); ++i) {
			if (hydros[i].hd().GetId() == id)
				return i;
		}
		return -1;
	}
	int GetMeshId(int id) {	
		for (int i = 0; i < surfs.size(); ++i) {
			if (surfs[i].GetId() == id)
				return i;
		}
		return -1;
	}
		
	static Function <void(String)> Print, PrintWarning, PrintError;	
	
	Upp::Vector<double> headAll;	// Common models data
	int Nb = 0;				
	
	double depth, rho, g, len;
	int discardNegDOF;
	double thres;
	int calcAwinf;
	double maxTimeA;
	int numValsA;
	int onlyDiagonal;
	
	String nemohPathPreprocessor, nemohPathSolver, nemohPathPostprocessor, nemohPathNew, nemohPathGREN;
	bool experimental;
	String foammPath;
	String hamsPath;
	
	void Load(String file, Function <bool(String, int pos)> Status, bool checkDuplicated);
	HydroClass &Join(Upp::Vector<int> &ids, Function <bool(String, int)> Status);
	void Symmetrize(int ids);
	void A0(int ids);
	void Ainf(int ids, double maxT);
	
	void LoadMesh(String file, Function <void(String, int pos)> Status, bool cleanPanels, bool checkDuplicated);
	void HealingMesh(int id, bool basic, Function <void(String, int pos)> Status);
	void OrientSurface(int id, Function <void(String, int)> Status);
	void ImageMesh(int id, int axis);
	void UnderwaterMesh(int id, Function <void(String, int pos)> Status);
	void RemoveMesh(int id);
	void JoinMesh(int idDest, int idOrig);
	Upp::Vector<int> SplitMesh(int id, Function <void(String, int pos)> Status);
	
	void AddFlatPanel(double x, double y, double z, double size, double panWidthX, double panWidthY);
	void AddRevolution(double x, double y, double z, double size, Upp::Vector<Pointf> &vals);
	void AddPolygonalPanel(double x, double y, double z, double size, Upp::Vector<Pointf> &vals);
	void AddWaterSurface(int id, char c);
	
	bool LoadSerializeJson(bool &firstTime);
	bool StoreSerializeJson();
	bool ClearTempFiles();
	static String GetTempFilesFolder() {return AppendFileNameX(GetAppDataFolder(), "BEMRosetta", "Temp");}
	
	const String bemFilesExt = ".1 .3 .hst .4 .12s .12d .out .cal .tec .inf .ah1 .lis .qtf .mat .dat .bem";
	String bemFilesAst;
	
	void Jsonize(JsonIO &json) {
		json
			("depth", depth)
			("rho", rho)
			("g", g)
			("length", len)
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
			("foammPath", foammPath)
			("hamsPath", hamsPath)
		;
	}
};

template <class T>
bool OUTB(int id, T total) {
	if (id < 0	|| id >= int(total))
		return true;
	return false;
}

String GetFASTVar(const String &strFile, String varName, String paragraph = "");
void SetFASTVar(String &strFile, String varName, String value, String paragraph = "");
	
	
class Mooring : public DeepCopyOption<Mooring> {
public:
	Mooring() {}
	Mooring(const Mooring &mooring, int) {
		lineTypes = clone(mooring.lineTypes);
		lineProperties = clone(mooring.lineProperties);
		connections = clone(mooring.connections);
	}
	
	bool Load(String file);
	bool Save(String file);
	String Test();
	void Jsonize(JsonIO &json);
	
	struct LineType : Moveable<LineType> {
		String name;
		double mass, diameter;
		void Jsonize(JsonIO &json);
	};
	Vector<LineType> lineTypes;
	
	struct LineProperty : Moveable<LineProperty> {
		String name, nameType;
		double length;
		String from, to;
		void Jsonize(JsonIO &json);
	};
	Vector<LineProperty> lineProperties;	

	struct Connection : Moveable<Connection> {
		String name;
		bool vessel;
		double x, y, z;
		void Jsonize(JsonIO &json);
	};
	Vector<Connection> connections;	
};
	
#endif
