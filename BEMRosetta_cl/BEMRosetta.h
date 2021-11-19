// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2021, the BEMRosetta author and contributors
#ifndef _BEM_Rosetta_BEM_Rosetta_h_
#define _BEM_Rosetta_BEM_Rosetta_h_

#include <Core/Core.h>
#include <SysInfo/SysInfo.h>
#include <Surface/Surface.h>

using namespace Upp;


class BEMData;

bool ConsoleMain(const Upp::Vector<String>& command, bool gui, Function <bool(String, int pos)> Status);
void SetBuildInfo(String &str);

bool PrintStatus(String s, int d);

class Hydro : public DeepCopyOption<Hydro> {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	enum BEM_SOFT {WAMIT, FAST_WAMIT, WAMIT_1_3, HAMS_WAMIT, WADAM_WAMIT, NEMOH, SEAFEM_NEMOH, AQWA, FOAMM, BEMROSETTA, UNKNOWN};
	
	enum DOF {SURGE = 0, SWAY, HEAVE, ROLL, PITCH, YAW};

	void Copy(const Hydro &hyd);
	Hydro(const Hydro &hyd, int) {Copy(hyd);}
	bool SaveAs(String file, Function <bool(String, int)> Status = Null, BEM_SOFT type = UNKNOWN, int qtfHeading = Null);
	void Report() const;
	Hydro(BEMData &_bem) : g(Null), h(Null), rho(Null), len(Null), Nb(Null), Nf(Null), Nh(Null), 
							dataFromW(Null), bem(&_bem) {id = idCount++;}
	virtual ~Hydro() noexcept {}	
	
	String GetCodeStr()	const {
		switch (code) {
		case WAMIT: 		return t_("Wamit");
		case WAMIT_1_3: 	return t_("Wamit.1.2.3");
		case HAMS_WAMIT: 	return t_("HAMS.Wamit");
		case WADAM_WAMIT: 	return t_("Wadam.Wamit");
		case FAST_WAMIT: 	return t_("Wamit.FAST");
		case NEMOH:			return t_("Nemoh");
		case SEAFEM_NEMOH:	return t_("SeaFEM.Nemoh");
		case AQWA:			return t_("AQWA");
		case FOAMM:			return t_("FOAMM");
		case BEMROSETTA:	return t_("BEMRosetta");
		case UNKNOWN:		return t_("Unknown");
		}
		return t_("Unknown");
	}
	
	String GetCodeStrAbr() const {
		switch (code) {
		case WAMIT: 		return t_("Wm.o");
		case WAMIT_1_3: 	return t_("Wm.1");
		case HAMS_WAMIT: 	return t_("HAM");
		case WADAM_WAMIT: 	return t_("WDM");
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
		
		int i = ib*6+idf;
		
		if (Ainf.size() > 0)
			if (!IsNull(Ainf(i, i)))
				return true;
		
		if (!A.IsEmpty())
			if (A[i][i].size() > 0)
				if (!IsNull(A[i][i][0]))
					return true;
		
		return false;			   
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
 	
	Upp::Array<Upp::Array<Eigen::VectorXd>> A;		// [6*Nb][6*Nb][Nf]	Added mass
	Upp::Array<Upp::Array<Eigen::VectorXd>> Ainf_w;	// [6*Nb][6*Nb][Nf]	Infinite frequency added mass (w)
    Eigen::MatrixXd Ainf;        			// (6*Nb, 6*Nb) 	Infinite frequency added mass
    Eigen::MatrixXd A0;        				// (6*Nb, 6*Nb)  	Infinite period added mass

	Eigen::MatrixXd Dlin;      				// (6*Nb, 6*Nb) 	Additional linear damping

    Upp::Array<Upp::Array<Eigen::VectorXd>> B; 		// [6*Nb][6*Nb][Nf]	Radiation damping
    Upp::Vector<double> head;				// [Nh]             Wave headings (deg)
    Upp::Vector<String> names;  			// {Nb}             Body names
    Upp::Array<Eigen::MatrixXd> C;			// [Nb](6, 6)		Hydrostatic restoring coefficients:
    Upp::Array<Eigen::MatrixXd> M;			// [Nb](6, 6)		Mass and inertia matrix
    Eigen::MatrixXd cb;          			// (3,Nb)           Centre of buoyancy
    Eigen::MatrixXd cg;          			// (3,Nb)     		Centre of gravity
    Eigen::MatrixXd c0;          			// (3,Nb)     		Centre of rotation
    BEM_SOFT code;        					// BEM_SOFT			BEM code 
    Upp::Vector<int> dof;      				// [Nb]            	Degrees of freedom for each body 
    Upp::Vector<int> dofOrder;				// [6*Nb]			DOF order
    
    Upp::Array<Upp::Array<Eigen::VectorXd>> Kirf;	// [6*Nb][6*Nb][Nt]	Radiation impulse response function IRF
    Eigen::VectorXd Tirf;	  				// [Nt]				Time-window for the calculation of the IRF
    
    double GetMass(int ib) {return M[ib](0, 0);}
    
    int GetHeadId(double hd) const;
	
    struct Forces : public DeepCopyOption<Forces> {
        Forces() {}
        Forces(const Forces &f, int) {
            ma = clone(f.ma);		ph = clone(f.ph);
            re = clone(f.re);		im = clone(f.im);
        }
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
    	void Clear() {ma.Clear(); ph.Clear(); re.Clear(); im.Clear();}
    };
    
    Forces ex; 								// Excitation
    Forces sc;			 					// Diffraction scattering
    Forces fk; 								// Froude-Krylov
    
  	typedef struct Forces RAO;
   
   	RAO rao;
    
    String description;

    struct StateSpace : public DeepCopyOption<StateSpace> {
        StateSpace() {}
        StateSpace(const StateSpace &s, int) {
        	TFS = clone(s.TFS);
			A_ss = clone(s.A_ss);
			B_ss = clone(s.B_ss);
			C_ss = clone(s.C_ss);
			ssFrequencies = clone(s.ssFrequencies);
			ssFreqRange = clone(s.ssFreqRange);
			ssFrequencies_index = clone(s.ssFrequencies_index);
			ssMAE = s.ssMAE;
        }
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
    
    struct QTF : public DeepCopyOption<QTF> {
        QTF() {
            fre.SetCount(6, Null);
            fim.SetCount(6, Null);
            fma.SetCount(6, Null);
            fph.SetCount(6, Null);
        }
        QTF(const QTF &q, int) {
        	ib = q.ib;
        	ih1 = q.ih1;	ih2 = q.ih2;
        	ifr1 = q.ifr1;	ifr2 = q.ifr2;
        	fre = clone(q.fre);		fim = clone(q.fim);
        	fma = clone(q.fma); 	fph = clone(q.fph);
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
				("ib",  ib)
				("ih1", ih1)
				("ih2", ih2)
				("ifr1",ifr1)
				("ifr2",ifr2)
				("fre", fre)
				("fim", fim)
				("fma", fma)
				("fph", fph)
			;
    	}
    	void Clear() {
    		fre.Clear();
    		fim.Clear();
    		fma.Clear();
    		fph.Clear();
    	}
    };
    Upp::Array<QTF> qtfsum, qtfdif;
    Upp::Vector<double> qtfw, qtfT, qtfhead;
    bool qtfdataFromW;
    
    struct QTFCases : public DeepCopyOption<QTFCases>  {
    	Upp::Vector<int> ib, ih1, ih2;
    	
    	QTFCases() {}
    	QTFCases(const QTFCases &q, int) {
    		ib = clone(q.ib);
    		ih1 = clone(q.ih1);
    		ih2 = clone(q.ih2);
    	}
    	void Clear() {
    		ib.Clear();
    		ih1.Clear();
    		ih2.Clear();
    	}
    } qtfCases;
    
    int GetQTFHeadId(double hd) const;
    static int GetQTFId(int lastid, const Upp::Array<Hydro::QTF> &qtfList, 
    			const QTFCases &qtfCases, int _ib, int _ih1, int _ih2, int _ifr1, int _ifr2);
	static void GetQTFList(const Upp::Array<Hydro::QTF> &qtfList, QTFCases &qtfCases);
							
    Upp::Vector<double> T;					// [Nf]    			Wave periods
    Upp::Vector<double> w;		 			// [Nf]             Wave frequencies
    bool dataFromW;
    Upp::Vector<double> Vo;   				// [Nb]             Displaced volume
    		
    void Dimensionalize();
    void Normalize();
    
    static String C_units(int i, int j);
    
    void SetC(int ib, const Eigen::MatrixXd &K);
	
	bool AfterLoad(Function <bool(String, int)> Status = Null);
	
	void Initialize_Forces();
	void Initialize_Forces(Forces &f, int _Nh = -1);
	void GetMaPh(Forces &f);
	void Normalize_Forces(Forces &f);
	void Dimensionalize_Forces(Forces &f);
	void Add_Forces(Forces &to, const Hydro &hydro, const Forces &from);
	void Symmetrize_Forces(bool xAxis);
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
	
	double A_dim(int ifr, int idf, int jdf) 	const {return dimen  ? A[idf][jdf][ifr]*rho_dim()/rho_ndim()  : A[idf][jdf][ifr]*(rho_dim()*pow(len, GetK_AB(idf, jdf)));}
	Eigen::VectorXd A_dim(int idf, int jdf) 	const {return dimen  ? A[idf][jdf]    *(rho_dim()/rho_ndim()) : A[idf][jdf]*     (rho_dim()*pow(len, GetK_AB(idf, jdf)));}
	double A_ndim(int ifr, int idf, int jdf) 	const {return !dimen ? A[idf][jdf][ifr]*(rho_ndim()/rho_dim()) : A[idf][jdf][ifr]/(rho_ndim()*pow(len, GetK_AB(idf, jdf)));}
	Eigen::VectorXd A_ndim(int idf, int jdf)	const {return !dimen ? A[idf][jdf]*(rho_ndim()/rho_dim()) : A[idf][jdf]*(1/(rho_ndim()*pow(len, GetK_AB(idf, jdf))));}
	double A_(bool ndim, int ifr, int idf, int jdf) const {return ndim ? A_ndim(ifr, idf, jdf) : A_dim(ifr, idf, jdf);}
	
	double A0_dim(int idf, int jdf)   		 	const {return dimen  ? A0(idf, jdf)*rho_dim()/rho_ndim() : A0(idf, jdf)  *(rho_dim()*pow(len, GetK_AB(idf, jdf)));}
	double A0_ndim(int idf, int jdf)  		 	const {return !dimen ? A0(idf, jdf)      : A0(idf, jdf)  /(rho_ndim()*pow(len, GetK_AB(idf, jdf)));}
	double A0_(bool ndim, int idf, int jdf) 	const {return ndim   ? A0_ndim(idf, jdf) : A0(idf, jdf);}
	double Ainf_dim(int idf, int jdf) 		 	const {return dimen  ? Ainf(idf, jdf)*rho_dim()/rho_ndim() : Ainf(idf, jdf)*(rho_dim()*pow(len, GetK_AB(idf, jdf)));}
	Eigen::MatrixXd Ainf_dim(int ib) const;
	double Ainf_ndim(int idf, int jdf)		 	const {return !dimen ? Ainf(idf, jdf) : Ainf(idf, jdf)/(rho_ndim()*pow(len, GetK_AB(idf, jdf)));}
	double Ainf_(bool ndim, int idf, int jdf) 	const {return ndim   ? Ainf_ndim(idf, jdf) : Ainf_dim(idf, jdf);}
	
	double B_dim(int ifr, int idf, int jdf)  	   const {return dimen  ? B[idf][jdf][ifr]*rho_dim()/rho_ndim() : B[idf][jdf][ifr]*(rho_dim()*pow(len, GetK_AB(idf, jdf))*w[ifr]);}
	Eigen::VectorXd B_dim(int idf, int jdf)  	   const;
	double B_ndim(int ifr, int idf, int jdf) 	   const {return !dimen ? B[idf][jdf][ifr]*(rho_ndim()/rho_dim()) : B[idf][jdf][ifr]/(rho_ndim()*pow(len, GetK_AB(idf, jdf))*w[ifr]);}
	Eigen::VectorXd B_ndim(int idf, int jdf) 	   const;
	double B_(bool ndim, int ifr, int idf, int jdf)const {return ndim ? B_ndim(ifr, idf, jdf) : B_dim(ifr, idf, jdf);}	
	
	double Kirf_dim(int it, int idf, int jdf)  	   	  const {return dimen ? Kirf[idf][jdf][it]*g_rho_dim()/g_rho_ndim()  : Kirf[idf][jdf][it]*(g_rho_dim()*pow(len, GetK_F(idf)));}
	double Kirf_ndim(int it, int idf, int jdf) 	   	  const {return !dimen ? Kirf[idf][jdf][it] : Kirf[idf][jdf][it]/(g_rho_ndim()*pow(len, GetK_F(idf)));}
	Eigen::VectorXd Kirf_ndim(int idf, int jdf) 	  const {return !dimen ? Kirf[idf][jdf]     : Kirf[idf][jdf]/(g_rho_ndim()*pow(len, GetK_F(idf)));}
	double Kirf_(bool ndim, int it, int idf, int jdf) const {return ndim ? Kirf_ndim(it, idf, jdf) : Kirf_dim(it, idf, jdf);}
	
	double Ainf_w_dim(int ifr, int idf, int jdf) 		const {return dimen  ? Ainf_w[idf][jdf][ifr]*rho_dim()/rho_ndim() : Ainf_w[idf][jdf][ifr]*(rho_dim()*pow(len, GetK_AB(idf, jdf)));}
	double Ainf_w_ndim(int ifr, int idf, int jdf) 		const {return !dimen ? Ainf_w[idf][jdf][ifr]*(rho_ndim()/rho_dim()) : Ainf_w[idf][jdf][ifr]/(rho_ndim()*pow(len, GetK_AB(idf, jdf)));}
	double Ainf_w_(bool ndim, int ifr, int idf, int jdf)const {return ndim   ? Ainf_w_ndim(ifr, idf, jdf) : Ainf_w_dim(ifr, idf, jdf);}
	
	double C_dim(int ib, int idf, int jdf)   	   const {return dimen  ? C[ib](idf, jdf)*g_rho_dim()/g_rho_ndim()  : C[ib](idf, jdf)*(g_rho_dim()*pow(len, GetK_C(idf, jdf)));}
	Eigen::MatrixXd C_dim(int ib) const;
	void C_dim();	
	double C_ndim(int ib, int idf, int jdf)  	   const {return !dimen ? C[ib](idf, jdf)  : C[ib](idf, jdf)/(g_rho_ndim()*pow(len, GetK_C(idf, jdf)));}
	double C_(bool ndim, int ib, int idf, int jdf) const {return ndim ? C_ndim(ib, idf, jdf) : C_dim(ib, idf, jdf);}

	Eigen::MatrixXd Dlin_dim(int ib) const;
	
	double F_ma_dim(const Forces &f, int _h, int ifr, int idf)  	   const {return dimen ? f.ma[_h](ifr, idf)*g_rho_dim()/g_rho_ndim()  : f.ma[_h](ifr, idf)*(g_rho_dim()*pow(len, GetK_F(idf)));}
	double F_re_dim(const Forces &f, int _h, int ifr, int idf)  	   const {return dimen ? f.re[_h](ifr, idf)*g_rho_dim()/g_rho_ndim()  : f.re[_h](ifr, idf)*(g_rho_dim()*pow(len, GetK_F(idf)));}
	double F_im_dim(const Forces &f, int _h, int ifr, int idf)  	   const {return dimen ? f.im[_h](ifr, idf)*g_rho_dim()/g_rho_ndim()  : f.im[_h](ifr, idf)*(g_rho_dim()*pow(len, GetK_F(idf)));}
	void F_dim(Forces &f);
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
		return std::complex<double>(B_(ndim, ifr, idf, jdf), w[ifr]*(A_(ndim, ifr, idf, jdf) - Ainf_(ndim, idf, jdf))/(!ndim ? 1. : w[ifr]));
	}
	
	std::complex<double> TFS_dim(int ifr, int idf, int jdf) 		const {return dimenSTS  ? sts[idf][jdf].TFS[ifr]*g_rho_dim()/g_rho_ndim() : sts[idf][jdf].TFS[ifr]*(rho_dim()*pow(len, GetK_AB(idf, jdf))*w[ifr]);}
	std::complex<double> TFS_ndim(int ifr, int idf, int jdf) 		const {return !dimenSTS ? sts[idf][jdf].TFS[ifr]*g_rho_ndim()/g_rho_dim() : sts[idf][jdf].TFS[ifr]/(rho_ndim()*pow(len, GetK_AB(idf, jdf))*w[ifr]);}
	std::complex<double> TFS_(bool ndim, int ifr, int idf, int jdf) const {return ndim ? TFS_ndim(ifr, idf, jdf) : TFS_dim(ifr, idf, jdf);}
	
	void SetId(int _id)			{id = _id;}
	int GetId()	const			{return id;}
	
	void CheckNaN();
		
	void Jsonize(JsonIO &json);
	
private:
	static const char *strDOF[];
	static const char *strDOFAbrev[];
	static const char *strDataToPlot[];
	static String C_units_base(int i, int j);
	BEMData *bem;
		
	void Symmetrize_Forces_Each0(const Forces &f, Forces &newf, const Upp::Vector<double> &newHead, double h, int ih, int idb);
	void Symmetrize_ForcesEach(const Forces &f, Forces &newf, const Upp::Vector<double> &newHead, int newNh, bool xAxis);
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
	
	enum DataToShow {DATA_A, DATA_B, DATA_AINFW, DATA_K, DATA_FORCE_SC, DATA_FORCE_FK, DATA_FORCE_EX, DATA_RAO, DATA_STS, DATA_STS2};
	enum DataToPlot {PLOT_A, PLOT_AINF, PLOT_A0, PLOT_B, PLOT_AINFW, PLOT_K, PLOT_FORCE_SC_MA, PLOT_FORCE_SC_PH,
				 PLOT_FORCE_FK_MA, PLOT_FORCE_FK_PH, PLOT_FORCE_EX_MA, PLOT_FORCE_EX_PH, 
				 PLOT_RAO_MA, PLOT_RAO_PH, PLOT_Z_MA, PLOT_Z_PH, PLOT_KR_MA, PLOT_KR_PH, 
				 PLOT_TFS_MA, PLOT_TFS_PH};
	enum DataMatrix {MAT_K, MAT_A, MAT_DAMP_LIN};
				 
	static const char *StrDataToPlot(DataToPlot dataToPlot) {
		return strDataToPlot[dataToPlot];
	}
	
	const BEMData &GetBEMData() const {return *bem;}
	
	bool IsLoadedA() 	 const {return !A.IsEmpty();}
	bool IsLoadedAinf_w()const {return !Ainf_w.IsEmpty();}
	bool IsLoadedAinf()  const {return Ainf.size() > 0;}
	bool IsLoadedDlin()  const {return Dlin.size() > 0;}
	bool IsLoadedA0()	 const {return A0.size() > 0;}
	bool IsLoadedB() 	 const {return !B.IsEmpty();}
	bool IsLoadedC()	 const {return !C.IsEmpty() && C[0].size() > 0 && !IsNull(C[0](0, 0));}
	bool IsLoadedM()	 const {return !M.IsEmpty() && M[0].size() > 0 && !IsNull(M[0](0, 0));}
	bool IsLoadedFex() 	 const {return !ex.ma.IsEmpty()  && ex.ma[0].size() > 0;}
	bool IsLoadedFsc() 	 const {return !sc.ma.IsEmpty()  && sc.ma[0].size() > 0;}
	bool IsLoadedFfk() 	 const {return !fk.ma.IsEmpty()  && fk.ma[0].size() > 0;}
	bool IsLoadedRAO() 	 const {return !rao.ma.IsEmpty() && rao.ma[0].size() > 0;}
	bool IsLoadedForce(const Forces &f) const {return !f.ma.IsEmpty();}
	bool IsLoadedStateSpace()	  const {return !sts.IsEmpty();}
	bool IsLoadedQTF() 	 const {return !qtfsum.IsEmpty();}
	bool IsLoadedKirf()	 const {return !Kirf.IsEmpty();}
	
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
	void GetA0();
		
	void GetK_IRF(double maxT = 120, int numT = 1000);
	double GetK_IRF_MaxT() const;
	static double GetK_IRF_MaxT(const Vector<double> &w);
	void GetAinf();
	void GetAinf_w();
	void InitAinf_w();
	void GetOgilvieCompliance(bool zremoval, bool thinremoval, bool decayingTail);
	void GetTranslationTo(double xto, double yto, double zto);
		
	void Join(const Upp::Vector<Hydro *> &hydrosp);
	
	String S_g()	const {return IsNull(g)   ? S("-") : Format("%.3f", g);}
	String S_h()	const {return IsNull(h)   ? S("-") : (h < 0 ? S(t_("INFINITY")) : Format("%.1f", h));}
	String S_rho() 	const {return IsNull(rho) ? S("-") : Format("%.3f", rho);}
	String S_len() 	const {return IsNull(len) ? S("-") : Format("%.1f", len);}

	String GetLastError()	{return lastError;}
};

bool IsNum(const Hydro::Forces &f);
			
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

class Mesh {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	enum MESH_FMT 			    		  		{WAMIT_GDF,  WAMIT_DAT,  NEMOH_DAT,  NEMOH_PRE,      AQWA_DAT,  HAMS_PNL,  STL_BIN,     STL_TXT,   EDIT,  UNKNOWN};	
	static constexpr const char *meshStr[]    = {"Wamit.gdf","Wamit.dat","Nemoh.dat","Nemoh premesh","AQWA.dat","HAMS.pnl","STL.Binary","STL.Text","Edit","Unknown"};	
	static constexpr const bool meshCanSave[] = {true, 	     false,	     true,		 false,			 false,		true,	   false,		true,	   false, false};       
	
	enum MESH_TYPE {MOVED, UNDERWATER, ALL};
	
	Mesh() {
		id = idCount++;
		cg  = Point3D(0, 0, 0);
		cg0 = Point3D(0, 0, 0);
		c0  = Point3D(0, 0, 0);
	}
	const char *GetCodeMeshStr() const {
		return meshStr[code];
	}
	static const char *GetCodeMeshStr(MESH_FMT c) {
		if (c < 0 || c >= UNKNOWN)
			return "Unknown";
		return meshStr[c];
	}
	static MESH_FMT GetCodeMeshStr(String fmt) {
		fmt = ToLower(Trim(fmt));
		for (int i = 0; i < sizeof(meshStr)/sizeof(char *); ++i)
			if (fmt == ToLower(meshStr[i]))
				return static_cast<MESH_FMT>(i);
		return UNKNOWN;
	}
	static bool MeshCanSave(MESH_FMT c) {
		if (c < 0 || c >= UNKNOWN)
			return false;
		return meshCanSave[c];
	}
	void SetCode(MESH_FMT _code){code = _code;}
	MESH_FMT GetCode()			{return code;}
	int GetId()	const			{return id;}

	String Load(String fileName, double rho, double g, bool cleanPanels);
	String Load(String fileName, double rho, double g, bool cleanPanels, bool &y0z, bool &x0z);
	
	String Heal(bool basic, double rho, double g, Function <bool(String, int pos)> Status);
	void Orient();
	void Join(const Surface &orig, double rho, double g);
	void Image(int axis);
		
	void AfterLoad(double rho, double g, bool onlyCG, bool isFirstTime);
	void Reset(double rho, double g);

	void SaveAs(String fileName, MESH_FMT type, double g, MESH_TYPE meshType, bool symX, bool symY, int &nNodes, int &nPanels);
	void SaveAs(String fileName, MESH_FMT type, double g, MESH_TYPE meshType, bool symX, bool symY) {
		int nNodes, nPanels;
		SaveAs(fileName, type, g, meshType, symX, symY, nNodes, nPanels);
	}
	
	void Report(double rho) const;
	
	bool IsSymmetricX();
	bool IsSymmetricY();
	
	double xProjectionPos, xProjectionNeg, yProjectionPos, yProjectionNeg, zProjectionPos, zProjectionNeg; 
	Point3D cb = Null;
	Point3D cg, cg0, c0;
	double mass = Null;
	Eigen::MatrixXd C;
	
	String name;
	String fileName;
	String header;
	
	Surface mesh, under, mesh0;
	
private:
	MESH_FMT code;
	int id;
	static int idCount;
};


class NemohMesh : public Mesh {
public:
	String LoadDat(String fileName, bool &x0z);
	static void SaveDat(String fileName, const Surface &surf, bool x0z);
	static void SavePreMesh(String fileName, const Surface &surf);
	void SaveKH(String fileName) const; 
	
	virtual ~NemohMesh() noexcept {}
};

class HAMSMesh : public Mesh {
public:
	String LoadPnl(String filefCaseName, bool &y0z, bool &x0z);
	static void SavePnl(String fileName, const Surface &surf, bool y0z, bool x0z);
	
	virtual ~HAMSMesh() noexcept {}
};

class WamitMesh : public Mesh {
public:
	String LoadDat(String fileName);
	String LoadGdf(String fileName, bool &y0z, bool &x0z);
	static void SaveGdf(String fileName, const Surface &surf, double g, bool y0z, bool x0z);
	void SaveHST(String fileName, double rho, double g) const; 

	virtual ~WamitMesh() noexcept {}
};

class AQWAMesh : public Mesh {
public:
	String LoadDat(String fileName);
	
	virtual ~AQWAMesh() noexcept {}
};


class Wamit : public HydroClass {
public:
	Wamit(BEMData &bem, Hydro *hydro = 0) : HydroClass(bem, hydro) {}
	bool Load(String file, Function <bool(String, int)> Status);
	bool Save(String file, Function <bool(String, int)> Status, bool force_T = false, int qtfHeading = Null);
	bool Save_out(String file, double g, double rho);
	virtual ~Wamit() noexcept {}
	
	bool LoadGdfMesh(String file);
	bool LoadDatMesh(String file);
	void SaveGdfMesh(String fileName);
	
	static void Save_hst_static(const Eigen::MatrixXd &C, String fileName, double rho, double g);
	
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
	bool Load_12(String fileName, bool isSum, Function <bool(String, int)> Status);
	
	void Save_1(String fileName, bool force_T = false);
	void Save_3(String fileName, bool force_T = false);
	void Save_hst(String fileName);
	void Save_4(String fileName, bool force_T = false);
	void Save_12(String fileName, bool isSum, Function <bool(String, int)> Status,
				bool force_T = false, bool force_Deg = true, int qtfHeading = Null);

	void Save_A(FileOut &out, Function <double(int, int)> fun, const Eigen::MatrixXd &base, String wavePeriod);
	void Save_AB(FileOut &out, int ifr);
	void Save_Forces(FileOut &out, int ifr);
	void Save_RAO(FileOut &out, int ifr);
};

class HAMS : public Wamit {
public:
	HAMS(BEMData &bem, Hydro *hydro = 0) : Wamit(bem, hydro) {}
	bool Load(String file, Function <bool(String, int)> Status, double g = 9.81);
	virtual ~HAMS() noexcept {}
	
	bool Load_Settings(String settingsFile);
	bool Load_HydrostaticMesh(String fileName, double rhog);
};

class Fast : public Wamit {
public:
	Fast(BEMData &bem, Hydro *hydro = 0) : Wamit(bem, hydro), WaveNDir(Null), WaveDirRange(Null) {}
	bool Load(String file, Function <bool(String, int)> Status, double g = 9.81);
	bool Save(String file, Function <bool(String, int)> Status, int qtfHeading = Null);
	virtual ~Fast() noexcept {}
	
private:
	bool Load_HydroDyn(String fileName);	
	void Save_HydroDyn(String fileName, bool force);
	bool Load_SS(String fileName);	
	void Save_SS(String fileName);
	
	String hydroFolder;
	int WaveNDir;
	double WaveDirRange;
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

class BEMBody : public MoveableAndDeepCopyOption<BEMBody> {
public:
	BEMBody();
	BEMBody(const BEMBody &d, int) : meshFile(d.meshFile), dof(clone(d.dof)), c0(clone(d.c0)), 
			cg(clone(d.cg)), mass(clone(d.mass)), 
			linearDamping(clone(d.linearDamping)), quadraticDamping(clone(d.quadraticDamping)), 
			hydrostaticRestoring(clone(d.hydrostaticRestoring)), externalRestoring(clone(d.externalRestoring)), 
			ndof(d.ndof) {}
	
	String meshFile, lidFile;
	int ndof;
	Vector<bool> dof;
	Eigen::Vector3d c0;	
	Eigen::Vector3d cg;
	Eigen::MatrixXd mass, linearDamping, quadraticDamping, hydrostaticRestoring, externalRestoring;
	
	int GetNDOF() const;
};


class BEMCase {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	BEMCase() {bodies.SetCount(1);}
	BEMCase(const BEMCase &bcase) : bodies(clone(bcase.bodies)), h(bcase.h),
		Nf(bcase.Nf), minF(bcase.minF), maxF(bcase.maxF), Nh(bcase.Nh),
		minH(bcase.minH), maxH(bcase.maxH), rho(bcase.rho), g(bcase.g),
		xeff(bcase.xeff), yeff(bcase.yeff), irf(bcase.irf), 
		irfStep(bcase.irfStep), irfDuration(bcase.irfDuration),
		showPressure(bcase.showPressure), 
		nFreeX(bcase.nFreeX), nFreeY(bcase.nFreeY), 
		domainX(bcase.domainX), domainY(bcase.domainY),
		nKochin(bcase.nKochin), minK(bcase.minK), maxK(bcase.maxK), solver(bcase.solver) {} 
		
	void Load(String file, const BEMData &bem);
	void SaveFolder(String folder, bool bin, int numCases, int numThreads, const BEMData &bem, int solver) const;
	Vector<String> Check(int solver) const;
	
	void BeforeSave(String folderBase, int numCases, bool deleteFolder) const;
	
	bool IsDof(int ib, int idf) {return bodies[ib].dof[idf];}
	bool IsNemoh() {return solver <= CAPYTAINE;}
	
	double h = Null;
	
	Upp::Vector<BEMBody> bodies;

	int Nf = Null;
	double minF = Null, maxF = Null;
	int Nh = Null;
	double minH = Null, maxH = Null;
	
	double rho = Null, g = Null;
	double xeff = 0, yeff = 0;
	
	bool irf = false;
	double irfStep = 0.1, irfDuration = 100;	
	bool showPressure = false;
	int nFreeX = 0, nFreeY = 0;
	double domainX = 0, domainY = 0;
	int nKochin = 0;
	double minK = 0, maxK = 0;
	
	enum Solver 			   		 			  {NEMOH, NEMOHv115, CAPYTAINE, HAMS, AQWA, NUMSOLVERS} solver;
	static const char *solverStr[];
	static constexpr const bool solverCanSave[] = {true,  true, 	 true, 		true, false};
	
	virtual ~BEMCase() noexcept {}
};

class AQWACase : public BEMCase {
public:
	bool Load(String fileName);
	Upp::Vector<String> Check() const;	
	
	virtual ~AQWACase() noexcept {}
};

class HamsCase : public BEMCase {
public:
	bool Load(String fileName);
	void SaveFolder(String folder, bool bin, int numCases, int numThreads, const BEMData &bem, int solver) const;
	Upp::Vector<String> Check() const;
	
	bool LoadHydrostatic(String fileName);
	
	virtual ~HamsCase() noexcept {}
	
private:
	void SaveFolder0(String folderBase, bool bin, int numCases, const BEMData &bem, bool deleteFolder, int numThreads) const;
	static void OutMatrix(FileOut &out, String header, const Eigen::MatrixXd &mat);
	static void InMatrix(FieldSplit &f, Eigen::MatrixXd &mat);
		
	void Save_Hydrostatic(String folderInput) const;
	void Save_ControlFile(String folderInput, int _nf, double _minf, double _maxf,
							int numThreads) const;
	void Save_Settings(String folderInput, bool thereIsLid) const;
	void Save_Bat(String folder, String batname, String caseFolder, bool bin, String solvName, String meshName) const;
};

class NemohCase : public BEMCase {
public:	
	bool Load(String fileName);
	void SaveFolder(String folder, bool bin, int numCases, int numThreads, const BEMData &bem, int solver) const;
	Upp::Vector<String> Check() const;
	
	void Save_Cal(String folder, int _nf, double _minf, double _maxf, const Vector<int> &nodes, const Vector<int> &panels) const;
	
	virtual ~NemohCase() noexcept {}
	
private:
	static int GetNumArgs(const FieldSplit &f);
	void LoadFreeSurface(const FileInLine &in, const FieldSplit &f);
	void LoadKochin(const FileInLine &in, const FieldSplit &f);

	void Save_Id(String folder) const;
	void Save_Bat(String folder, String batname, String caseFolder, bool bin, String preName, String solvName, String postName) const;
	void Save_Mesh_cal(String folder, int ib, String meshFile, Mesh &mesh, int npanels, bool x0z, Eigen::Vector3d cg, double rho, double g) const;
	void Save_Mesh_bat(String folder, String caseFolder, const Vector<String> &meshes, String meshName, bool bin) const;
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
	
	bool Save_KH(String folder) const;
	static bool Save_KH_static(const Eigen::MatrixXd &C, String fileKH);
	
private:
	String folder;
	BEMCase dcase;
	
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
	
	HydroClass &Duplicate(int id);
		
	Upp::Array<HydroClass> hydros;
	Upp::Array<Mesh> surfs;
	
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
	
	int calcAinf, calcAinf_w;
	double maxTimeA;
	int numValsA;
	int onlyDiagonal;
	
	String nemohPathPreprocessor, nemohPathSolver, nemohPathPostprocessor, nemohPathNew, nemohPathGREN;
	bool experimental;
	String foammPath;
	String hamsPath, hamsMeshPath;
	int volWarning, volError;
	
	void LoadBEM(String file, Function <bool(String, int pos)> Status = Null, bool checkDuplicated = false);
	HydroClass &Join(Upp::Vector<int> &ids, Function <bool(String, int)> Status = Null);
	void Symmetrize(int id, bool xAxis);
	void A0(int id);
	void Kirf(int id, double maxT);
	void Ainf(int id);
	void Ainf_w(int id);
	void OgilvieCompliance(int id, bool zremoval, bool thinremoval, bool decayingTail);
	void TranslationTo(int id, double xto, double yto, double zto);
	
	void LoadMesh(String file, Function <bool(String, int pos)> Status, bool cleanPanels, bool checkDuplicated);
	void HealingMesh(int id, bool basic, Function <bool(String, int pos)> Status);
	void OrientSurface(int id, Function <bool(String, int)> Status);
	void ImageMesh(int id, int axis);
	void UnderwaterMesh(int id, Function <bool(String, int pos)> Status);
	void RemoveMesh(int id);
	void JoinMesh(int idDest, int idOrig);
	Upp::Vector<int> SplitMesh(int id, Function <bool(String, int pos)> Status);
	
	void AddFlatPanel(double x, double y, double z, double size, double panWidthX, double panWidthY);
	void AddRevolution(double x, double y, double z, double size, Upp::Vector<Pointf> &vals);
	void AddPolygonalPanel(double x, double y, double z, double size, Upp::Vector<Pointf> &vals);
	void AddWaterSurface(int id, char c);
	
	bool LoadSerializeJson();
	bool StoreSerializeJson();
	bool ClearTempFiles();
	static String GetTempFilesFolder() {return AppendFileNameX(GetAppDataFolder(), "BEMRosetta", "Temp");}
	
	const String bemFilesExt = ".1 .2 .3 .hst .4 .12s .12d .out .in .cal .tec .inf .ah1 .lis .qtf .mat .dat .bem";
	String bemFilesAst;
	
	void Jsonize(JsonIO &json) {
		if (json.IsLoading()) {
			volWarning = 1;
			volError = 10;
		}
		json
			("depth", depth)
			("rho", rho)
			("g", g)
			("length", len)
			//("discardNegDOF", discardNegDOF)
			//("thres", thres)
			("calcAwinf", calcAinf)
			("calcAwinfw", calcAinf_w)
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
			("hamsMeshPath", hamsMeshPath)
			("volWarning", volWarning)
			("volError", volError)
		;
	}
};

template <class T>
bool OUTB(int id, T total) {
	if (id < 0	|| id >= int(total))
		return true;
	return false;
}


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
	
String FormatDoubleEmpty(double val);
String FormatIntEmpty(int val);


String GetFASTVar(const String &strFile, String varName, String paragraph = "");
void SetFASTVar(String &strFile, String varName, String value, String paragraph = "");


class FASTFiles {
public:
	void Load(String file) {
		String path = GetFileFolder(file);
		
		fast.fileName = file;
		elastodyn.fileName = AppendFileNameX(path, fast.GetString("EDFile"));
		hydrodyn.fileName = AppendFileNameX(path, fast.GetString("HydroFile"));
	}
	void Save() {
		fast.Save();
		elastodyn.Save();
		hydrodyn.Save();
	}
	
private:
	class File {
	public:
		String fileName;
		String fileText;
		bool isChanged;
		
		void Save() const {
			if (!isChanged)
				return;
			if (!SaveFile(fileName, fileText))
				throw Exc(Format(t_("Impossible to save file '%s'"), fileName));
		}
		
		String GetString(String var) {
			if (fileText.IsEmpty()) {
				fileText = LoadFile(fileName);
				if (fileText.IsEmpty())
					throw Exc(Format(t_("Impossible to read file '%s'"), fileName));
			}
			String res;
			Vector<String> vars = Split(var, "/");
			if (vars.size() == 2)
				res = GetFASTVar(fileText, vars[1], vars[0]);
			else if (vars.size() == 1)		
				res = GetFASTVar(fileText, vars[0]);
			else
				throw Exc(Format(t_("Wrong variable '%s' in GetString"), var));
			
			if (res[0] == '\"')		// Remove quotes
				res = res.Mid(1);
			if (res[res.GetCount()-1] == '\"')
				res = res.Left(res.GetCount()-1);
			return res;
		}
		double GetDouble(String var) {
			double ddata = ScanDouble(GetString(var));
			if (IsNull(ddata))
				throw Exc(Format(t_("Wrong variable '%s' in GetDouble"), var));
			return ddata;
		}
		double GetInt(String var) {
			int ddata = ScanInt(GetString(var));
			if (IsNull(ddata))
				throw Exc(Format(t_("Wrong variable '%s' in GetInt"), var));
			return ddata;
		}
		bool GetBool(String var) {
			String data = ToLower(GetString(var));
			if (data == "true")
				return true;
			if (data == "false")
				return true;
			int idata = ScanInt(data);
			if (idata == 1)
				return true;
			if (idata == 0)
				return false;
			throw Exc(Format(t_("Wrong variable '%s' in GetBool"), var));
		}
		double GetMatrixVal(String var, int row, int col) {
			if (fileText.IsEmpty()) {
				fileText = LoadFile(fileName);
				if (fileText.IsEmpty())
					throw Exc(Format(t_("Impossible to read file '%s'"), fileName));
			}
			int posIni, posEnd;
			GetMatrixIds(var, row, col, posIni, posEnd);
			
			String data = fileText.Mid(posIni, posEnd-posIni);
			double ddata = ScanDouble(data);
			if (IsNull(ddata))
				throw Exc(Format(t_("Problem reading variable '%s' in GetMatrix %d, %d"), var, row, col));
			return ddata;
		}
		Eigen::MatrixXd GetMatrix(String var, int rows, int cols) {
			Eigen::MatrixXd ret(rows, cols);
			
			for (int r = 0; r < rows; ++r)
				for (int c = 0; c < rows; ++c)
					ret(r, c) = GetMatrixVal(var, r, c);
					
			return ret;
		}
		void SetMatrixVal(String var, int row, int col, double val) {
			int posIni, posEnd;
			GetMatrixIds(var, row, col, posIni, posEnd);
			
			int delta = posEnd-posIni-1;
			
			fileText = fileText.Left(posIni) + S(" ") + FormatDoubleSize(val, delta, true) + fileText.Mid(posEnd);
		}
		
		void SetString(String var, String val) {
			val = S("\"") + val + S("\"");
			SetString0(var, val);
		}
		void SetInt(String var, int val) {
			SetString0(var, FormatInt(val));
		}
		void SetDouble(String var, double val) {
			SetString0(var, FormatDoubleSize(val, 10));
		}
		void SetBool(String var, bool val) {
			SetString0(var, val ? "True" : "False");
		}
		
	private:
		void SetString0(String var, String val) {
			if (fileText.IsEmpty()) {
				fileText = LoadFile(fileName);
				if (fileText.IsEmpty())
					throw Exc(Format(t_("Impossible to read file '%s'"), fileName));
			}
			Vector<String> vars = Split(var, "/");
			if (vars.size() == 2) {
				isChanged = true;
				return SetFASTVar(fileText, vars[1], val, vars[0]);
			} else if (vars.size() == 1) {
				isChanged = true;	
				return SetFASTVar(fileText, vars[0], val);
			} else
				throw Exc(Format(t_("Wrong variable '%s' in SetString"), var));
		}
		void GetMatrixIds(String var, int row, int col, int &posIni, int &posEnd) {
			if (fileText.IsEmpty()) {
				fileText = LoadFile(fileName);
				if (fileText.IsEmpty())
					throw Exc(Format(t_("Impossible to read file '%s'"), fileName));
			}
			int id;
			if ((id = fileText.Find(var)) < 0)
				throw Exc(Format(t_("Wrong variable '%s' in GetMatrixIds"), var));
			if ((id = fileText.ReverseFindAfter("\n", id)) < 0)
				throw Exc(Format(t_("Problem reading variable '%s' in GetMatrixIds"), var));
			
			for (int i = 0; i < row; ++i) {
				if ((id = fileText.FindAfter("\n", id)) < 0)
					throw Exc(Format(t_("Problem reading variable '%s' row %d in GetMatrixIds"), var, row));
			}
			for (int ic = 0; ic <= col; ++ic) {
				posIni = id;
				while (id < fileText.GetCount() && IsSpace(fileText[id]))
					id++;
				while (id < fileText.GetCount() && !IsSpace(fileText[id]))
					id++;
				posEnd = id;
			}
			if (id == fileText.GetCount())
				throw Exc(Format(t_("Problem reading variable '%s' col %d in GetMatrixIds"), var, col));				
		}
	};

public:
	File fast, elastodyn, hydrodyn;
};

	
#endif
