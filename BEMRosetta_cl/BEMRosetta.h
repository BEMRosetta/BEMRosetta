// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#ifndef _BEM_Rosetta_BEM_Rosetta_h_
#define _BEM_Rosetta_BEM_Rosetta_h_

#include <Core/Core.h>
#include <SysInfo/SysInfo.h>
#include <ScatterDraw/DataSource.h>
#include <Surface/Surface.h>
#include <STEM4U/Mooring.h>

using namespace Upp;


class BEM;

bool ConsoleMain(const UVector<String>& command, bool gui, Function <bool(String, int pos)> Status);
void SetBuildInfo(String &str);
String GetSystemInfo();

bool PrintStatus(String s, int d);

class Hydro : public DeepCopyOption<Hydro> {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	enum BEM_FMT 							   {WAMIT, 		  WAMIT_1_3, 					FAST_WAMIT, 				 	HAMS_WAMIT,   WADAM_WAMIT,   NEMOH,   SEAFEM_NEMOH,   AQWA,   					  FOAMM,   BEMROSETTA, 	   	   UNKNOWN};
	static constexpr const char *bemStr[]    = {"Wamit .out", "Wamit .1.2.3.hst.789.ss.12", "FAST .dat.1.2.3.hst.789.ss.12","HAMS Wamit", "Wadam Wamit", "Nemoh", "SeaFEM Nemoh", "AQWA .lis .ah1 .qtf [W]", "FOAMM", "BEMRosetta .bem", "By extension"};
	static constexpr const bool bemCanSave[] = {true, 	      true,	     					true,		 			 	 	false,		  false,		  false,   false, 		   true,  					  false,   true,			   true};       
	static constexpr const char *bemExt[]	 = {"*.out", 	  "*.1",	     				"*.1",		 			 	 	"",		   	  "",		      "",      "", 		   	   "*.qtf", 				  "",      "*.bem",		   	   "*.*"};       
	
	static void ResetIdCount()		{idCount = 0;}
	
	static const char *GetBemStr(BEM_FMT c) {
		if (c < 0 || c > UNKNOWN)
			return "Unknown";
		return bemStr[c];
	}
	static int GetBemStrCount() {return sizeof(bemStr)/sizeof(char *);}
	
	static BEM_FMT GetCodeBemStr(String fmt) {
		fmt = ToLower(Trim(fmt));
		for (int i = 0; i < GetBemStrCount(); ++i)
			if (fmt == ToLower(bemStr[i]))
				return static_cast<BEM_FMT>(i);
		return UNKNOWN;
	}
		
	void Copy(const Hydro &hyd);
	Hydro(const Hydro &hyd, int) {Copy(hyd);}
	bool SaveAs(String file, Function <bool(String, int)> Status = Null, BEM_FMT type = UNKNOWN, int qtfHeading = Null);
	void Report() const;
	Hydro(const BEM &_bem) : g(Null), h(Null), rho(Null), len(Null), Nb(Null), Nf(Null), Nh(Null), 
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
			if (IsNum(Ainf(i, i)))
				return true;
		
		if (!A.IsEmpty())
			if (A[i][i].size() > 0)
				if (IsNum(A[i][i][0]))
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
 	
	UArray<UArray<VectorXd>> A;		// [6*Nb][6*Nb][Nf]	Added mass
	UArray<UArray<VectorXd>> Ainf_w;// [6*Nb][6*Nb][Nf]	Infinite frequency added mass (w)
    MatrixXd Ainf;        			// (6*Nb, 6*Nb) 	Infinite frequency added mass
    MatrixXd A0;        			// (6*Nb, 6*Nb)  	Infinite period added mass

	MatrixXd linearDamping;			// (6*Nb, 6*Nb) 	Additional linear damping	ALWAYS DIMENSIONAL

    UArray<UArray<VectorXd>> B; 	// [6*Nb][6*Nb][Nf]	Radiation damping
    UVector<double> head;			// [Nh]             Wave headings (deg)
    UVector<String> names;  		// {Nb}             Body names
    UArray<MatrixXd> C;				// [Nb](6, 6)		Hydrostatic restoring coefficients:
    UArray<MatrixXd> M;				// [Nb](6, 6)		Mass and inertia matrix		ALWAYS DIMENSIONAL
    MatrixXd cb;          			// (3,Nb)           Centre of buoyancy
    MatrixXd cg;          			// (3,Nb)     		Centre of gravity
    MatrixXd c0;          			// (3,Nb)     		Centre of motion
    BEM_FMT code;        			// BEM_FMT			BEM code 
    UVector<int> dof;      			// [Nb]            	Degrees of freedom for each body 
    
    UArray<UArray<VectorXd>> Kirf;	// [6*Nb][6*Nb][Nt]	Radiation impulse response function IRF
    VectorXd Tirf;	  				// [Nt]				Time-window for the calculation of the IRF
    
    double GetMass(int ib) const	{return M[ib](0, 0);}
    
    int GetHeadId(double hd) const;
    int GetHeadIdMD(const std::complex<double> &h) const;
	
    struct Forces : public DeepCopyOption<Forces> {
        Forces() {}
        Forces(const Forces &f, int) {
            force = clone(f.force);
        }
    	UArray<MatrixXcd> force;	// [Nh](Nf, 6*Nb) 	
    
    	void Jsonize(JsonIO &json) {
			json
				("force", force)
			;
    	}
    	void Clear() {force.Clear();}
    };
    
    Forces ex; 								// Excitation
    Forces sc;			 					// Diffraction scattering
    Forces fk; 								// Froude-Krylov
    
    enum FORCE {NONE, SCATTERING, FK, ALL, QTFSUM, QTFDIF};
        
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
	    UArray<std::complex<double>> TFS;
		MatrixXd A_ss;
		VectorXd B_ss;
		VectorXd C_ss;
		VectorXd ssFrequencies, ssFreqRange, ssFrequencies_index;
		double ssMAE = Null;
		
		void GetTFS(const UVector<double> &w);
		
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
    UArray<UArray<StateSpace>> sts;			// (6*Nb, 6*Nb)		State space data
    int dimenSTS;							// false if data is dimensionless
    String stsProcessor;
    
    
    VectorXd  qw;		 							 // [Nf]             Wave frequencies
    VectorXcd qh;									 // [Nh]             Wave headings
    UArray<UArray<UArray<MatrixXcd>>> qtfsum, qtfdif;// [Nb][Nh][6](Nf, Nf)	
    bool qtfdataFromW = true;
     
    void InitQTF(UArray<UArray<UArray<MatrixXcd>>> &qtf, int nb, int nh, int nf) {
        qtf.SetCount(nb);
        for (int ib = 0; ib < nb; ++ib) {
            qtf[ib].SetCount(nh);
        	for (int ih = 0; ih < nh; ++ih) {    
        		qtf[ib][ih].SetCount(6);
        		for (int idf = 0; idf < 6; ++idf) 
        			qtf[ib][ih][idf].resize(nf, nf);
        	}
        }
    }
    
    VectorXcd mdhead;							// [Nh]             Wave headings
	UArray<UArray<UArray<VectorXd>>> md;		// [Nb][Nh][6](Nf)	
	int mdtype = 0;
	
    static void InitMD(UArray<UArray<UArray<VectorXd>>> &md, int nb, int nh, int nf) {
        md.SetCount(nb);
        for (int ib = 0; ib < nb; ++ib) {
            md[ib].SetCount(nh);
        	for (int ih = 0; ih < nh; ++ih) {    
        		md[ib][ih].SetCount(6);
        		for (int idf = 0; idf < 6; ++idf) 
        			md[ib][ih][idf].setConstant(nf, NaNDouble);
        	}
        }
    }
    					
    UVector<double> T;					// [Nf]    			Wave periods
    UVector<double> w;		 			// [Nf]             Wave frequencies
    int dataFromW = Null;
    UVector<double> Vo;   				// [Nb]             Displaced volume
    		
    void Dimensionalize();
    void Normalize();
    
    static String C_units(int i, int j);
    
    void SetC(int ib, const MatrixXd &K);
	
	bool AfterLoad(Function <bool(String, int)> Status = Null);
	
	void Initialize_Forces();
	void Initialize_Forces(Forces &f, int _Nh = -1);
	void Normalize_Forces(Forces &f);
	void Dimensionalize_Forces(Forces &f);
	void Add_Forces(Forces &to, const Hydro &hydro, const Forces &from);
	void Symmetrize_Forces(bool xAxis);
	void Symmetrize();
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
	
	VectorXd Get_w() 							const {return Map<const VectorXd>(w, w.size());}
	VectorXd Get_T()		 					const {return Map<const VectorXd>(T, T.size());}
	
	double A_dim(int ifr, int idf, int jdf) 	const {return dimen  ? A[idf][jdf][ifr]*rho_dim()/rho_ndim()  : A[idf][jdf][ifr]*(rho_dim()*pow(len, GetK_AB(idf, jdf)));}
	VectorXd A_dim(int idf, int jdf) 			const {return dimen  ? A[idf][jdf]    *(rho_dim()/rho_ndim()) : A[idf][jdf]*     (rho_dim()*pow(len, GetK_AB(idf, jdf)));}
	double A_ndim(int ifr, int idf, int jdf) 	const {return !dimen ? A[idf][jdf][ifr]/**(rho_ndim()/rho_dim())*/ : A[idf][jdf][ifr]/(rho_ndim()*pow(len, GetK_AB(idf, jdf)));}
	VectorXd A_ndim(int idf, int jdf)			const {return !dimen ? A[idf][jdf]/**(rho_ndim()/rho_dim())*/ : A[idf][jdf]*(1/(rho_ndim()*pow(len, GetK_AB(idf, jdf))));}
	double A_(bool ndim, int ifr, int idf, int jdf) const {return ndim ? A_ndim(ifr, idf, jdf) : A_dim(ifr, idf, jdf);}
	MatrixXd A_(bool ndim, int ifr, int ib) 	const;
	
	double A0_dim(int idf, int jdf)   		 	const {return dimen  ? A0(idf, jdf)*rho_dim()/rho_ndim() : A0(idf, jdf)  *(rho_dim()*pow(len, GetK_AB(idf, jdf)));}
	double A0_ndim(int idf, int jdf)  		 	const {return !dimen ? A0(idf, jdf)      : A0(idf, jdf)  /(rho_ndim()*pow(len, GetK_AB(idf, jdf)));}
	double A0_(bool ndim, int idf, int jdf) 	const {return ndim   ? A0_ndim(idf, jdf) : A0_dim(idf, jdf);}
	double Ainf_dim(int idf, int jdf) 		 	const {return dimen  ? Ainf(idf, jdf)*rho_dim()/rho_ndim() : Ainf(idf, jdf)*(rho_dim()*pow(len, GetK_AB(idf, jdf)));}
	MatrixXd Ainf_(bool ndim, int ib) const;
	double Ainf_ndim(int idf, int jdf)		 	const {return !dimen ? Ainf(idf, jdf) : Ainf(idf, jdf)/(rho_ndim()*pow(len, GetK_AB(idf, jdf)));}
	double Ainf_(bool ndim, int idf, int jdf) 	const {return ndim   ? Ainf_ndim(idf, jdf) : Ainf_dim(idf, jdf);}
	
	double B_dim(int ifr, int idf, int jdf)  	const {return dimen  ? B[idf][jdf][ifr]*rho_dim()/rho_ndim() : B[idf][jdf][ifr]*(rho_dim()*pow(len, GetK_AB(idf, jdf))*w[ifr]);}
	VectorXd B_dim(int idf, int jdf)  	   		const;
	double B_ndim(int ifr, int idf, int jdf) 	const {return !dimen ? B[idf][jdf][ifr]/**(rho_ndim()/rho_dim())*/ : B[idf][jdf][ifr]/(rho_ndim()*pow(len, GetK_AB(idf, jdf))*w[ifr]);}
	VectorXd B_ndim(int idf, int jdf) 	   		const;
	double B_(bool ndim, int ifr, int idf, int jdf)const {return ndim ? B_ndim(ifr, idf, jdf) : B_dim(ifr, idf, jdf);}	
	MatrixXd B_(bool ndim, int ifr, int ib) 	const;
	
	double Kirf_dim(int it, int idf, int jdf)  	   	  const {return dimen ? Kirf[idf][jdf][it]*g_rho_dim()/g_rho_ndim()  : Kirf[idf][jdf][it]*(g_rho_dim()*pow(len, GetK_F(idf)));}
	double Kirf_ndim(int it, int idf, int jdf) 	   	  const {return !dimen ? Kirf[idf][jdf][it] : Kirf[idf][jdf][it]/(g_rho_ndim()*pow(len, GetK_F(idf)));}
	VectorXd Kirf_ndim(int idf, int jdf) 	 		  const {return !dimen ? Kirf[idf][jdf]     : Kirf[idf][jdf]/(g_rho_ndim()*pow(len, GetK_F(idf)));}
	double Kirf_(bool ndim, int it, int idf, int jdf) const {return ndim ? Kirf_ndim(it, idf, jdf) : Kirf_dim(it, idf, jdf);}
	
	double Ainf_w_dim(int ifr, int idf, int jdf) 		const {return dimen  ? Ainf_w[idf][jdf][ifr]*rho_dim()/rho_ndim() : Ainf_w[idf][jdf][ifr]*(rho_dim()*pow(len, GetK_AB(idf, jdf)));}
	double Ainf_w_ndim(int ifr, int idf, int jdf) 		const {return !dimen ? Ainf_w[idf][jdf][ifr]/**(rho_ndim()/rho_dim())*/ : Ainf_w[idf][jdf][ifr]/(rho_ndim()*pow(len, GetK_AB(idf, jdf)));}
	double Ainf_w_(bool ndim, int ifr, int idf, int jdf)const {return ndim   ? Ainf_w_ndim(ifr, idf, jdf) : Ainf_w_dim(ifr, idf, jdf);}
	
	double C_dim(int ib, int idf, int jdf)   	   const {return dimen  ? C[ib](idf, jdf)/**g_rho_dim()/g_rho_ndim()*/  : C[ib](idf, jdf)*(g_rho_dim()*pow(len, GetK_C(idf, jdf)));}
	MatrixXd C_(bool ndim, int ib) 				   const;
	void C_dim();	
	double C_ndim(int ib, int idf, int jdf)  	   const {return !dimen ? C[ib](idf, jdf)  : C[ib](idf, jdf)/(g_rho_ndim()*pow(len, GetK_C(idf, jdf)));}
	double C_(bool ndim, int ib, int idf, int jdf) const {return ndim ? C_ndim(ib, idf, jdf) : C_dim(ib, idf, jdf);}

	double Md_dim(int idf, int ih, int ifr)  const {
		int ib = int(idf/6);
		idf -= 6*ib;
		return dimen ? md[ib][ih][idf](ifr)*g_rho_ndim()/g_rho_dim()  : md[ib][ih][idf](ifr)*(g_rho_dim()*pow(len, GetK_F(idf)));}
	double Md_ndim(int idf, int ih, int ifr) const {
		int ib = int(idf/6);
		idf -= 6*ib;
		return !dimen ? md[ib][ih][idf](ifr) : md[ib][ih][idf](ifr)/(g_rho_ndim()*pow(len, GetK_F(idf)));}
	double Md_(bool ndim, int idf, int ih, int ifr) const {return ndim ? Md_ndim(idf, ih, ifr) : Md_dim(idf, ih, ifr);}	// idf: body, jdf: heading, [Nb][Nh][6](Nf)
	
	MatrixXd Dlin_dim(int ib) const;
	
	std::complex<double> F_dim(const Forces &f, int _h, int ifr, int idf)  const {
		return dimen ? f.force[_h](ifr, idf)*g_rho_ndim()/g_rho_dim()  : f.force[_h](ifr, idf)*(g_rho_dim()*pow(len, GetK_F(idf)));}
	void F_dim(Forces &f);
	std::complex<double> F_ndim(const Forces &f, int _h, int ifr, int idf) const {
		if (!dimen) 
			return f.force[_h](ifr, idf); 
		else 
			return f.force[_h](ifr, idf)/(g_rho_ndim()*pow(len, GetK_F(idf)));
	}
	std::complex<double> F_(bool ndim, const Forces &f, int _h, int ifr, int idf) const {
		return ndim ? F_ndim(f, _h, ifr, idf) : F_dim(f, _h, ifr, idf);
	}
	std::complex<double> R_dim(const Forces &f,  int _h, int ifr, int idf) const {
		return dimen ? f.force[_h](ifr, idf)*g_rho_ndim()/g_rho_dim()  : f.force[_h](ifr, idf)*(g_rho_dim()*pow(len, GetK_RAO(idf)));
	}
	std::complex<double> R_ndim(const Forces &f, int _h, int ifr, int idf) const {
		if (!dimen) 
			return f.force[_h](ifr, idf); 
		else 
			return f.force[_h](ifr, idf)/(g_rho_ndim()*pow(len, GetK_RAO(idf)));
	}
	std::complex<double> R_(bool ndim, const Forces &f, int _h, int ifr, int idf) const {
		return ndim ? R_ndim(f, _h, ifr, idf) : R_dim(f, _h, ifr, idf);
	}
	
	template <class T>
	T F_dim(T f, int idf)  	      const {return  dimen ? f*g_rho_dim()/g_rho_ndim() : f*(g_rho_dim()*pow(len, GetK_F(idf)));}
	template <class T>
	T F_ndim(T f, int idf) 	      const {return !dimen ? f : f/(g_rho_ndim()*pow(len, GetK_F(idf)));}
	template <class T>
	T F_(bool ndim, T f, int idf) const {return   ndim ? F_ndim(f, idf) : F_dim(f, idf);}

	VectorXcd F_(bool ndim, const Forces &f, int _h, int ifr) const;

	inline std::complex<double> Z(bool ndim, int ifr, int idf, int jdf) const {
		return std::complex<double>(B_(ndim, ifr, idf, jdf), w[ifr]*(A_(ndim, ifr, idf, jdf) - Ainf_(ndim, idf, jdf))/(!ndim ? 1. : w[ifr]));
	}
	
	std::complex<double> TFS_dim(int ifr, int idf, int jdf) 		const {return dimenSTS  ? sts[idf][jdf].TFS[ifr]*g_rho_dim()/g_rho_ndim() : sts[idf][jdf].TFS[ifr]*(rho_dim()*pow(len, GetK_AB(idf, jdf))*w[ifr]);}
	std::complex<double> TFS_ndim(int ifr, int idf, int jdf) 		const {return !dimenSTS ? sts[idf][jdf].TFS[ifr]/**g_rho_ndim()/g_rho_dim()*/ : sts[idf][jdf].TFS[ifr]/(rho_ndim()*pow(len, GetK_AB(idf, jdf))*w[ifr]);}
	std::complex<double> TFS_(bool ndim, int ifr, int idf, int jdf) const {return ndim ? TFS_ndim(ifr, idf, jdf) : TFS_dim(ifr, idf, jdf);}
	
	double Tdof(int ib, int idf) const;
	double Theave(int ib) const;
	double Troll(int ib) const;
	double Tpitch(int ib) const;
	double GM(int ib, int idf) const;
	double GMroll(int ib) const;
	double GMpitch(int ib) const;
	
	void SetId(int _id)			{id = _id;}
	int GetId()	const			{return id;}
	
	void CheckNaN();
		
	void Jsonize(JsonIO &json);
	
private:
	static const char *strDataToPlot[];
	static String C_units_base(int i, int j);
	const BEM *bem;
		
	void Symmetrize_Forces_Each0(const Forces &f, Forces &newf, const UVector<double> &newHead, double h, int ih, int idb);
	void Symmetrize_ForcesEach(const Forces &f, Forces &newf, const UVector<double> &newHead, int newNh, bool xAxis);
	int id;
	static int idCount;
	 
	static void GetOldAB(const UArray<MatrixXd> &oldAB, UArray<UArray<VectorXd>> &AB);
	static void SetOldAB(UArray<MatrixXd> &oldAB, const UArray<UArray<VectorXd>> &AB);
	
	void ResetForces1st(Hydro::FORCE force);
	
public:
	enum DataToShow {DATA_A, DATA_B, DATA_AINFW, DATA_K, DATA_FORCE_SC, DATA_FORCE_FK, DATA_FORCE_EX, DATA_RAO, DATA_STS, DATA_STS2, DATA_MD};
	enum DataToPlot {PLOT_A, PLOT_AINF, PLOT_A0, PLOT_B, PLOT_AINFW, PLOT_KIRF, PLOT_FORCE_SC_1, PLOT_FORCE_SC_2,
				 PLOT_FORCE_FK_1, PLOT_FORCE_FK_2, PLOT_FORCE_EX_1, PLOT_FORCE_EX_2, 
				 PLOT_RAO_1, PLOT_RAO_2, PLOT_Z_1, PLOT_Z_2, PLOT_KR_1, PLOT_KR_2, 
				 PLOT_TFS_1, PLOT_TFS_2, PLOT_MD};
	enum DataMatrix {MAT_K, MAT_A, MAT_DAMP_LIN, MAT_M};
				 
	static const char *StrDataToPlot(DataToPlot dataToPlot) {
		return strDataToPlot[dataToPlot];
	}
	
	const BEM &GetBEM() const {return *bem;}
	
	bool IsLoadedA	   (int i = 0, int j = 0)const {return A.size() > i && A[i].size() > j && A[i][j].size() > 0 && IsNum(A[i][j][0]);}
	bool IsLoadedAinf_w(int i = 0, int j = 0)const {return Ainf_w.size() > i && Ainf_w[i].size() > j && Ainf_w[i][j].size() > 0 && IsNum(Ainf_w[i][j][0]);}
	bool IsLoadedAinf  (int i = 0, int j = 0)const {return Ainf.rows() > i && Ainf.cols() > j && IsNum(Ainf(i, j));}
	bool IsLoadedA0	   (int i = 0, int j = 0)const {return A0.rows() > i && A0.cols() > j && IsNum(A0(i, j));}
	bool IsLoadedLinearDamping()  			 const {return linearDamping.size() > 0;}
	bool IsLoadedB	   (int i = 0, int j = 0)const {return B.size() > i && B[i].size() > j && B[i][j].size() > 0 && IsNum(B[i][j][0]);}
	bool IsLoadedC(int ib = 0, int idf = 0, int jdf = 0)	const {return C.size() > ib && C[ib].rows() > idf && C[ib].cols() > jdf && IsNum(C[ib](idf, jdf));}
	bool IsLoadedM(int ib = 0, int idf = 0, int jdf = 0)	const {return M.size() > ib && M[ib].rows() > idf && M[ib].cols() > jdf && IsNum(M[ib](idf, jdf));}
	
	bool IsLoadedFex(int idf = 0, int ih = 0)const {return IsLoadedForce(ex, idf, ih);}
	bool IsLoadedFsc(int idf = 0, int ih = 0)const {return IsLoadedForce(sc, idf, ih);}
	bool IsLoadedFfk(int idf = 0, int ih = 0)const {return IsLoadedForce(fk, idf, ih);}
	bool IsLoadedRAO(int idf = 0, int ih = 0)const {return IsLoadedForce(rao,idf, ih);}
	bool IsLoadedForce(const Forces &f, int idf = 0, int ih = 0)
											 const {return f.force.size() > ih && f.force[ih].cols() > idf && IsNum(f.force[ih](0, idf));}
	
	bool IsLoadedStateSpace()	  			 const {return !sts.IsEmpty();}
	bool IsLoadedQTF(bool isSum) 			 const {return isSum ? !qtfsum.IsEmpty() : !qtfdif.IsEmpty();}
	bool IsLoadedMD(int ib = 0, int ih = 0)	 const {return md.size() > ib && md[ib].size() > ih && md[ib][ih].size() == 6 && md[ib][ih][0].size() > 0 && IsNum(md[ib][ih][0](0));}
	static bool IsLoadedMD(const UArray<UArray<UArray<VectorXd>>> &mD, int ib = 0, int ih = 0) {return mD.size() > ib && mD[ib].size() > ih && mD[ib][ih].size() == 6 && mD[ib][ih][0].size() > 0 && IsNum(mD[ib][ih][0](0));}
	bool IsLoadedKirf(int idf=0,int jdf = 0) const {return Kirf.size() > idf && Kirf[idf].size() > jdf && Kirf[idf][jdf].size() > 0 && IsNum(Kirf[idf][jdf][0]);}
	bool IsLoadedVo()	 	 				 const {return Vo.size() == 3 && IsNum(Vo[0]); }
	
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
	
//	const UVector<int> &GetOrder() const	{return dofOrder;}
//	void SetOrder(UVector<int> &order)		{dofOrder = pick(order);}
	
	int GetW0();
	void Get3W0(int &id1, int &id2, int &id3);
	void GetA0();
		
	void GetK_IRF(double maxT = 120, int numT = 1000);
	double GetK_IRF_MaxT() const;
	static double GetK_IRF_MaxT(const UVector<double> &w);
	void GetAinf();
	void GetAinf_w();
	void GetRAO();
	static VectorXcd GetRAO(double w, const MatrixXd &Aw, const MatrixXd &Bw, const VectorXcd &Fwh, 
				const MatrixXd &C, const MatrixXd &M, const MatrixXd &D, const MatrixXd &D2);
	void InitAinf_w();
	void GetOgilvieCompliance(bool zremoval, bool thinremoval, bool decayingTail, bool haskind, UVector<int> &vidof, UVector<int> &vjdof);
	void GetTranslationTo(double xto, double yto, double zto);
	
	void DeleteFrequencies(const UVector<int> &idFreq);
	void DeleteFrequenciesQTF(const UVector<int> &idFreqQTF);
	void DeleteHeadings(const UVector<int> &idHead);
	void DeleteHeadingsMD(const UVector<int> &idHead);
	void DeleteHeadingsQTF(const UVector<int> &idHeadQTF);
	void ResetForces(Hydro::FORCE force, bool forceMD, Hydro::FORCE forceQtf);
	void MultiplyDOF(double factor, const UVector<int> &idDOF, bool a, bool b, bool diag, bool f, bool md, bool qtf);
	void SwapDOF(int ib, int idof1, int idof2);
	
	void SymmetrizeDOF();
	
	void FillFrequencyGapsABForces(bool zero, int maxFreq);
	void FillFrequencyGapsQTF(bool zero, int maxFreq);
	
	void CopyQTF_MD();
	
	void Join(const UVector<Hydro *> &hydrosp);
	
	String S_g()	const {return !IsNum(g)   ? S("-") : Format("%.3f", g);}
	String S_h()	const {return !IsNum(h)   ? S("-") : (h < 0 ? S(t_("INFINITY")) : Format("%.1f", h));}
	String S_rho() 	const {return !IsNum(rho) ? S("-") : Format("%.3f", rho);}
	String S_len() 	const {return !IsNum(len) ? S("-") : Format("%.1f", len);}

	String GetLastError()	{return lastError;}
	
	static String K_units(bool ndim, int r, int c) {
		if (ndim) {
			if (r < 3 && c < 3)
				return "m^2";
			if (r < 3)
				return "m^3/rad";
			if (c < 3)
				return "m^3";
			return "m^4/rad";			
		} else {
			if (r < 3 && c < 3)
				return "N/m";
			if (r < 3)
				return "N/rad";
			if (c < 3)
				return "N";
			return "N-m/rad";
		}
	}
	static String Kirf_units(bool ndim, int r, int c) {
		return K_units(ndim, r, c);
	}
	static String B_units(bool ndim, int r, int c) {
		if (ndim) {
			if (r < 3 && c < 3)
				return "m^3-rad/s";
			if (r < 3)
				return "m^4/s";
			if (c < 3)
				return "m^4/s/rad";
			return "m^5/s";
		} else {
			if (r < 3 && c < 3)
				return "N/(m/s)";
			if (r < 3)
				return "N/(rad/s)";
			if (c < 3)
				return "N-s";
			return "N-m/(rad/s)";
		}
	}
	static String A_units(bool ndim, int r, int c) {
		if (ndim) {
			if (r < 3 && c < 3)
				return "m^3";
			if (r < 3)
				return "m^4/rad";
			if (c < 3)
				return "m^4";
			return "m^5/rad";			
		} else {
			if (r < 3 && c < 3)
				return "kg";
			if (r < 3)
				return "kg-m/rad";
			if (c < 3)
				return "kg-m";
			return "kg-m^2/rad";
		}
	}
	static String M_units(int r, int c) {
		if (r < 3 && c < 3)
			return "kg";
		else
			return "kg-m";
	}
	static String F_units(bool ndim, int r) {
		if (ndim) {
			if (r < 3)
				return "m^2";
			return "m^3";
		} else {
			if (r < 3)
				return "N/m";
			return "N-m/m";
		}
	}
	static String MD_units(bool ndim, int r) {
		if (ndim) {
			if (r < 3)
				return t_("m");
			return t_("m2");
		} else {
			if (r < 3)
				return t_("N/m2");
			return t_("N-m/m2");
		}
	} 
};

bool IsNum(const Hydro::Forces &f);
			
class HydroData {
public:
	HydroData() : data(0) {}
	HydroData(const BEM &bem, Hydro *_data = 0) {
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
	HydroClass(const BEM &bem, Hydro *hydro = 0) : hd(bem, hydro)	{}
	virtual ~HydroClass() noexcept			{}
	bool Load(String file);
	bool Save(String file);
	
	HydroData hd;	
};

class Mesh {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	enum MESH_FMT 			    		  		{WAMIT_GDF,  WAMIT_DAT,  NEMOH_DAT,  NEMOHFS_DAT,   NEMOH_PRE,      AQWA_DAT,  HAMS_PNL,  STL_BIN,     STL_TXT,   EDIT,  MSH_TDYN,   BEM_MESH, DIODORE_DAT,   UNKNOWN};	
	static constexpr const char *meshStr[]    = {"Wamit.gdf","Wamit.dat","Nemoh.dat","NemohFS.dat", "Nemoh premesh","AQWA.dat","HAMS.pnl","STL.Binary","STL.Text","Edit","TDyn.msh", "BEMR",   "DIODORE.dat", "Unknown"};	
	static constexpr const bool meshCanSave[] = {true, 	     false,	     true,		 false,			false, 		    false,		true,	   true,		true,	   false, false, 	  true, 	true,		   false};       
	static constexpr const char *meshExt[]	  = {"*.gdf", 	 "*.dat",	 "*.dat",	 "*.dat", 		"",		        "*.dat",	"*.pnl",   "*.stl", 	"*.stl",   "",    "*.msh",	  "*.bemr", "*.dat", 	   "*.*"};       
	
	enum MESH_TYPE {MOVED, UNDERWATER, ALL};
	
	Mesh() {
		id = idCount++;
		cg  = Point3D(0, 0, 0);
		cg0 = Point3D(0, 0, 0);
		c0  = Point3D(0, 0, 0);
	}
	bool IsEmpty() {return mesh.IsEmpty();}
	static void ResetIdCount()		{idCount = 0;}
	
	const char *GetMeshStr() const {
		return meshStr[code];
	}
	static const char *GetMeshStr(MESH_FMT c) {
		if (c < 0 || c > UNKNOWN)
			return "Unknown";
		return meshStr[c];
	}
	static int GetMeshStrCount() {return sizeof(meshStr)/sizeof(char *);}
	
	static MESH_FMT GetCodeMeshStr(String fmt) {
		fmt = ToLower(Trim(fmt));
		for (int i = 0; i < GetMeshStrCount(); ++i)
			if (fmt == ToLower(meshStr[i]))
				return static_cast<MESH_FMT>(i);
		return UNKNOWN;
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
	void Move(double dx, double dy, double dz, double ax, double ay, double az, 
			  double rho, double g, bool setnewzero);
	void Move(const double *pos, double rho, double g, bool setnewzero);
	void Move(const float *pos, double rho, double g, bool setnewzero);
		
	void AfterLoad(double rho, double g, bool onlyCG, bool isFirstTime);
	void Reset(double rho, double g);

	void GZ(double from, double to, double delta, double angle, double rho, double g, 
		double tolerance,
		Function <bool(String, int pos)> Status, 
		UVector<double> &dataangle, UVector<double> &dataGZ, UVector<double> &dataMoment,
		UVector<double> &vol, UVector<double> &disp, UVector<double> &wett, UVector<double> &wplane,
		UVector<double> &draft, UVector<Point3D> &dcb, UVector<Point3D> &dcg, String &error);
	void GZ(double from, double to, double delta, double angleCalc, double rho, double g,
		double tolerance, UVector<double> &dataangle, UVector<double> &datagz, String &error);

	double GMroll(double rho, double g) const;
	double GMpitch(double rho, double g) const;
	
	void SaveAs(String fileName, MESH_FMT type, double g, MESH_TYPE meshType, bool symX, bool symY, int &nNodes, int &nPanels);
	void SaveAs(String fileName, MESH_FMT type, double g, MESH_TYPE meshType, bool symX, bool symY) {
		int nNodes, nPanels;
		SaveAs(fileName, type, g, meshType, symX, symY, nNodes, nPanels);
	}
	
	void Report(double rho) const;
	
	bool IsSymmetricX();
	bool IsSymmetricY();
	
	double xProjectionPos, xProjectionNeg, yProjectionPos, yProjectionNeg, zProjectionPos, zProjectionNeg; 
	Pointf cgZ0surface = Null;
	Point3D cb = Null;
	Point3D cg, cg0, c0;
	double mass = Null;
	MatrixXd C;
	
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
	String LoadDatFS(String fileName, bool &x0z);
	void SaveDat(String fileName, const Surface &surf, bool x0z) const;
	static void SavePreMesh(String fileName, const Surface &surf);
	void SaveKH(String fileName) const; 
	
	virtual ~NemohMesh() noexcept {}

private:
	String LoadDat0(String fileName, bool &x0z);
	void SaveDat0(String fileName, const Surface &surf, bool x0z) const;
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

class DiodoreMesh : public Mesh {
public:
	String LoadDat(String fileName);
	void SaveDat(String fileName, const Surface &surf);
	
	virtual ~DiodoreMesh() noexcept {}
};

class Wamit : public HydroClass {
public:
	Wamit(const BEM &bem, Hydro *hydro = 0) : HydroClass(bem, hydro) {}
	bool Load(String file, Function <bool(String, int)> Status);
	bool Save(String file, Function <bool(String, int)> Status, bool force_T = false, int qtfHeading = Null);
	bool Save_out(String file);
	virtual ~Wamit() noexcept {}
	
	bool LoadGdfMesh(String file);
	bool LoadDatMesh(String file);
	void SaveGdfMesh(String fileName);
	
	static void Save_hst_static(const MatrixXd &C, String fileName, double rho, double g);
	
	bool Load_frc2(String fileName);
	
protected:
	void ProcessFirstColumnPot(UVector<double> &w, UVector<double> &T);
	bool ProcessFirstColumn1_3(UVector<double> &w, UVector<double> &T);
	
	bool Load_cfg(String fileName);
	int iperin = Null, iperout = Null;
	bool Load_pot(String fileName);
	bool Load_gdf(String fileName);
	
	bool Load_out();							
	void Load_A(FileInLine &in, MatrixXd &A);
	bool Load_Scattering(String fileName);
	bool Load_FK(String fileName);
	
	bool Load_Forces(String fileName, Hydro::Forces &force);
		
	bool Load_1(String fileName);				
	bool Load_3(String fileName);
	bool Load_hst(String fileName);
	bool Load_4(String fileName);
	bool Load_12(String fileName, bool isSum, Function <bool(String, int)> Status);
	bool Load_789(String fileName);
	bool Load_789_0(String fileName, int type, UArray<UArray<UArray<VectorXd>>> &qtf);
	
	void Save_1(String fileName, bool force_T = false);
	void Save_3(String fileName, bool force_T = false);
	void Save_hst(String fileName);
	void Save_4(String fileName, bool force_T = false);
	void Save_12(String fileName, bool isSum, Function <bool(String, int)> Status,
				bool force_T = false, bool force_Deg = true, int qtfHeading = Null);
	void Save_789(String fileName, bool force_T, bool force_Deg);
	void Save_FRC(String fileName);
	void Save_POT(String fileName);
		
	void Save_A(FileOut &out, Function <double(int, int)> fun, const MatrixXd &base, String wavePeriod);
	void Save_AB(FileOut &out, int ifr);
	void Save_Forces(FileOut &out, int ifr);
	void Save_RAO(FileOut &out, int ifr);
};

class HAMS : public Wamit {
public:
	HAMS(const BEM &bem, Hydro *hydro = 0) : Wamit(bem, hydro) {}
	bool Load(String file, Function <bool(String, int)> Status);
	virtual ~HAMS() noexcept {}
	
	bool Load_Settings(String settingsFile);
	bool Load_HydrostaticMesh(String fileName, double rhog);
};

class Fast : public Wamit {
public:
	Fast(const BEM &bem, Hydro *hydro = 0) : Wamit(bem, hydro), WaveNDir(Null), WaveDirRange(Null) {}
	bool Load(String file, Function <bool(String, int)> Status);
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
	Foamm(const BEM &bem, Hydro *hydro = 0) : HydroClass(bem, hydro) {}
	bool Load(String file);
	void Get_Each(int ibody, int idf, int jdf, double from, double to, const UVector<double> &freqs, Function <bool(String, int)> Status, Function <void(String)> FOAMMMessage);
	void Get(const UVector<int> &ibs, const UVector<int> &idfs, const UVector<int> &jdfs,
		const UVector<double> &froms, const UVector<double> &tos, const UVector<UVector<double>> &freqs, 
		Function <bool(String, int)> Status, Function <void(String)> FOAMMMessage);
	virtual ~Foamm() noexcept {}
	
protected:
	bool Load_mat(String fileName, int ib, int jb, bool loadCoeff);
};

class BEMBody : public MoveableAndDeepCopyOption<BEMBody> {
public:
	BEMBody();
	BEMBody(const BEMBody &d, int) : meshFile(d.meshFile), lidFile(d.lidFile),
			ndof(d.ndof), dof(clone(d.dof)), c0(clone(d.c0)), 
			cg(clone(d.cg)), mass(clone(d.mass)), 
			linearDamping(clone(d.linearDamping)), quadraticDamping(clone(d.quadraticDamping)), 
			hydrostaticRestoring(clone(d.hydrostaticRestoring)), externalRestoring(clone(d.externalRestoring)) 
			 {}
	
	String meshFile, lidFile;
	int ndof;
	UVector<bool> dof;
	Vector3d c0;	
	Vector3d cg;
	MatrixXd mass, linearDamping, quadraticDamping, hydrostaticRestoring, externalRestoring;
	
	int GetNDOF() const;
};


class BEMCase {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	BEMCase() {bodies.SetCount(1);}
	BEMCase(const BEMCase &bcase) : h(bcase.h), bodies(clone(bcase.bodies)), 
		Nf(bcase.Nf), minF(bcase.minF), maxF(bcase.maxF), Nh(bcase.Nh),
		minH(bcase.minH), maxH(bcase.maxH), rho(bcase.rho), g(bcase.g),
		xeff(bcase.xeff), yeff(bcase.yeff), irf(bcase.irf), 
		irfStep(bcase.irfStep), irfDuration(bcase.irfDuration),
		showPressure(bcase.showPressure), 
		nFreeX(bcase.nFreeX), nFreeY(bcase.nFreeY), 
		domainX(bcase.domainX), domainY(bcase.domainY),
		nKochin(bcase.nKochin), minK(bcase.minK), maxK(bcase.maxK), solver(bcase.solver) {} 
		
	void Load(String file, const BEM &bem);
	void SaveFolder(String folder, bool bin, int numCases, int numThreads, const BEM &bem, int solver) const;
	UVector<String> Check(int solver) const;
	
	void BeforeSave(String folderBase, int numCases, bool deleteFolder) const;
	
	bool IsDof(int ib, int idf) {return bodies[ib].dof[idf];}
	bool IsNemoh() {return solver <= CAPYTAINE;}
	
	double h = Null;
	
	UVector<BEMBody> bodies;

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
	UVector<String> Check() const;	
	
	virtual ~AQWACase() noexcept {}
};

class HamsCase : public BEMCase {
public:
	bool Load(String fileName);
	void SaveFolder(String folder, bool bin, int numCases, int numThreads, const BEM &bem, int solver) const;
	UVector<String> Check() const;
	
	bool LoadHydrostatic(String fileName);
	
	virtual ~HamsCase() noexcept {}
	
private:
	void SaveFolder0(String folderBase, bool bin, int numCases, const BEM &bem, bool deleteFolder, int numThreads) const;
	static void OutMatrix(FileOut &out, String header, const MatrixXd &mat);
	static void InMatrix(FieldSplit &f, MatrixXd &mat);
		
	void Save_Hydrostatic(String folderInput) const;
	void Save_ControlFile(String folderInput, int _nf, double _minf, double _maxf,
							int numThreads) const;
	void Save_Settings(String folderInput, bool thereIsLid, const BEM &bem) const;
	void Save_Bat(String folder, String batname, String caseFolder, bool bin, String solvName, String meshName) const;
};

class NemohCase : public BEMCase {
public:	
	bool Load(String fileName);
	void SaveFolder(String folder, bool bin, int numCases, int numThreads, const BEM &bem, int solver) const;
	UVector<String> Check() const;
	
	void Save_Cal(String folder, int _nf, double _minf, double _maxf, const UVector<int> &nodes, const UVector<int> &panels, bool isCapy) const;
	
	virtual ~NemohCase() noexcept {}
	
private:
	static int GetNumArgs(const FieldSplit &f);
	void LoadFreeSurface(const FileInLine &in, const FieldSplit &f);
	void LoadKochin(const FileInLine &in, const FieldSplit &f);

	void Save_Id(String folder) const;
	void Save_Bat(String folder, String batname, String caseFolder, bool bin, String preName, String solvName, String postName) const;
	void Save_Mesh_cal(String folder, int ib, String meshFile, Mesh &mesh, int npanels, bool x0z, Vector3d cg, double rho, double g) const;
	void Save_Mesh_bat(String folder, String caseFolder, const UVector<String> &meshes, String meshName, bool bin) const;
	void Save_Input(String folder) const;
	
	void SaveFolder0(String folder, bool bin, int numCases, const BEM &bem, 
					bool deleteFolder, int solver) const;
};

class Nemoh : public HydroClass {
public:
	Nemoh(const BEM &bem, Hydro *hydro = 0) : HydroClass(bem, hydro) {}
	bool Load(String file, double rho = Null);
	void Save(String file);
	virtual ~Nemoh() noexcept {}
	
	bool LoadDatMesh(String file);
	void SaveDatMesh(String file); 
	
	bool Save_KH(String folder) const;
	static bool Save_KH_static(const MatrixXd &C, String fileKH);
	
	void SetFolder(String f)	{folder = f;}
	
	bool Load_Hydrostatics(String subfolder = "Mesh");
	static bool Load_Hydrostatics_static(String folder, int Nb, MatrixXd &cg, MatrixXd &cb, UVector<double> &Vo);
	void Save_Hydrostatics(String subfolder) const;
	static void Save_Hydrostatics_static(String folder, int Nb, const MatrixXd &cg, const MatrixXd &cb, const UVector<double> &Vo);

private:
	String folder;
	BEMCase dcase;
	
	bool Load_Cal(String fileName);
	bool Load_Inf(String fileName);
	bool Load_KH(String subfolder = "Mesh");
	bool Load_Radiation(String fileName);
	bool Load_Excitation(String folder);
	bool Load_Diffraction(String folder);
	bool Load_FroudeKrylov(String folder);
	bool Load_Forces(Hydro::Forces &f, String nfolder, String fileName);
	bool Load_IRF(String fileName);
};

class Aqwa : public HydroClass {
public:
	Aqwa(const BEM &bem, Hydro *hydro = 0) : HydroClass(bem, hydro) {}
	bool Load(String file, double rho = Null);
	bool Save(String file, Function <bool(String, int)> Status);
	virtual ~Aqwa() noexcept {}
	
private:
	bool Load_AH1();
	bool Load_LIS();
	bool Load_QTF();
	void Save_QTF(String file, Function <bool(String, int)> Status);
};


UVector<int> NumSets(int num, int numsets);	


class FieldSplitWamit: public FieldSplit {
public:
	FieldSplitWamit(FileInLine &_in) : FieldSplit(_in) {}
	
	void LoadWamitJoinedFields(String _line);// Trick for "glued" fields in Wamit
};

String FormatWam(double d);

class BEM {
public:
	BEM();
	
	HydroClass &Duplicate(int id);
		
	UArray<HydroClass> hydros;
	UArray<Mesh> surfs;
	
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
	
	UVector<double> headAll;	// Common models data
	UArray<std::complex<double>> headAllMD;
	//UVector<int> orderHeadAll;
	
	int Nb = 0;				
	
	double depth, rho, g, len;
	
	enum DOFType {DOF123, DOFSurgeSway, DOFxyz};
	static DOFType dofType;
	enum HeadingType {HEAD_180_180, HEAD_0_360};
	static HeadingType headingType;
	static const char *strDOFType[];
	static const char *strHeadingType[];
	
	int calcAinf = false, calcAinf_w = false;
	double maxTimeA = Null;
	int numValsA = Null;
	int onlyDiagonal;
	
	String nemohPathPreprocessor, nemohPathSolver, nemohPathPostprocessor, nemohPathNew, nemohPathGREN;
	bool experimental;
	String foammPath;
	String hamsPath, hamsMeshPath;
	int volWarning, volError;
	String csvSeparator;
	
	void LoadBEM(String file, Function <bool(String, int pos)> Status = Null, bool checkDuplicated = false);
	HydroClass &Join(UVector<int> &ids, Function <bool(String, int)> Status = Null);
	void SymmetrizeForces(int id, bool xAxis);
	void Symmetrize(int id);
	void A0(int id);
	void Kirf(int id, double maxT);
	void Ainf(int id);
	void Ainf_w(int id);
	void RAO(int id);
	void OgilvieCompliance(int id, bool zremoval, bool thinremoval, bool decayingTail, bool haskind, UVector<int> &vidof, UVector<int> &vjdof);
	void TranslationTo(int id, double xto, double yto, double zto);
	void DeleteHeadingsFrequencies(int id, const UVector<int> &idFreq, const UVector<int> &idFreqQTF, 
										   const UVector<int> &idHead, const UVector<int> &idHeadMD, const UVector<int> &idHeadQTF);
	void ResetForces(int id, Hydro::FORCE force, bool forceMD, Hydro::FORCE forceQtf);										
	void MultiplyDOF(int id, double factor, const UVector<int> &idDOF, bool a, bool b, bool diag, bool f, bool md, bool qtf);
	void SwapDOF(int id, int ib, int idof1, int idof2);
	
	void RemoveHydro(int id);
		
	void FillFrequencyGapsABForces(int id, bool zero, int maxFreq);
	void FillFrequencyGapsQTF(int id, bool zero, int maxFreq);
	
	void LoadMesh(String file, Function <bool(String, int pos)> Status, bool cleanPanels, bool checkDuplicated);
	void HealingMesh(int id, bool basic, Function <bool(String, int pos)> Status);
	void OrientSurface(int id, Function <bool(String, int)> Status);
	void ImageMesh(int id, int axis);
	void UnderwaterMesh(int id, Function <bool(String, int pos)> Status);
	void RemoveMesh(int id);
	void JoinMesh(int idDest, int idOrig);
	UVector<int> SplitMesh(int id, Function <bool(String, int pos)> Status);
	
	void CopyQTF_MD(int id);
		
	void AddFlatPanel(double x, double y, double z, double size, double panWidthX, double panWidthY);
	void AddRevolution(double x, double y, double z, double size, UVector<Pointf> &vals);
	void AddPolygonalPanel(double x, double y, double z, double size, UVector<Pointf> &vals);
	void AddWaterSurface(int id, char c);
	
	bool LoadSerializeJson();
	bool StoreSerializeJson();
	bool ClearTempFiles();
	static String GetTempFilesFolder() {return AppendFileNameX(GetAppDataFolder(), "BEMRosetta", "Temp");}
	
	void UpdateHeadAll();
	void UpdateHeadAllMD();
	
	const String bemFilesExt = ".1 .2 .3 .hst .4 .12s .12d .frc .pot .out .in .cal .tec .inf .ah1 .lis .qtf .mat .dat .bem .fst";
	const String bstFilesExt = ".in .out .fst .1 .2 .3 .hst .4 .12s .12d .frc .pot .cal .tec .inf .ah1 .lis .qtf .mat .dat .bem";	// Priority
	const UVector<String> bemExtSets = {".1.2.3.hst.4.12s.12d.frc.pot", ".lis.qtf"};
	String bemFilesAst;
	
	int GetBEMExtSet(String file) {
		String ext = GetFileExt(file);
		for (int i = 0; i < bemExtSets.size(); ++i)
			if (bemExtSets[i].Find(ext) >= 0)
				return i;
		return -1;
	}
	
	void Jsonize(JsonIO &json) {
		int idofType, iheadingType;
		if (json.IsLoading()) {
			idofType = 0;
			iheadingType = 0;
			csvSeparator = ";";
		} else {
			idofType = dofType;
			iheadingType = headingType;
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
			("dofType", idofType)
			("headingType", iheadingType)
			("csvSeparator", csvSeparator)
		;
		if (json.IsLoading()) {
			dofType = BEM::DOFType(idofType);
			headingType = BEM::HeadingType(iheadingType);
		}
	}
	
	enum DOF {SURGE = 0, SWAY, HEAVE, ROLL, PITCH, YAW};

	static String StrBDOF(int i, bool abrev) {
		int ib = i/6 + 1;
		int idf = i - (ib - 1)*6;
		return Format("%d%s%s", ib, abrev ? "" : ".", StrDOF(idf, abrev));
	}
	
	static String StrBDOF2(int i, int j, bool abrev) {
		if (i != j) {
			int ib = i/6 + 1;
			int idf = i - (ib - 1)*6;
			int jb = j/6 + 1;
			int jdf = j - (jb - 1)*6;
			if (ib != jb)
				return Format("%d%s%s_%d.%s", ib, abrev ? "" : ".", StrDOF(idf, abrev), jb, StrDOF(jdf, abrev));
			else
				return Format("%d%s%s_%s", ib, abrev ? "" : ".", StrDOF(idf, abrev), StrDOF(jdf, abrev));
		} else
			return StrBDOF(i, abrev);
	}
		
	static String StrBDOFFull(int i) {
		int ib = i/6 + 1;
		int idf = i - (ib - 1)*6;
		return Format("Body #%d. DoF: %s", ib, StrDOF(idf));
	}

	static String StrBDOFFull(int i, int j) {
		if (i != j) {
			int ib = i/6 + 1;
			int idf = i - (ib - 1)*6;
			int jb = j/6 + 1;
			int jdf = j - (jb - 1)*6;
			if (ib != jb)
				return Format("Body #%d, DoF: %s. Body #%d, DoF: %s", ib, StrDOF(idf), jb, StrDOF(jdf));
			else
				return Format("Body #%d. DoF: %s, DoF: %s", ib, StrDOF(idf), StrDOF(jdf));
		} else
			return StrBDOF(i, false);
	}
		
	static int StrDOF_len() {
		auto MaxLen = [&] (const char *str[]) {
			int mx = 0;
			for (int i = 0; i < 6; ++i) {
				if (strlen(str[i]) > mx)
					mx = int(strlen(str[i]));
			}
			return mx;
		};
		if (dofType == DOF123)
			return MaxLen(strDOFnum);
		else if (dofType == DOFSurgeSway)
			return MaxLen(strDOFtext);
		else
			return MaxLen(strDOFxyz);	
	}

	static const char *StrDOF(int i) {
		if (dofType == DOF123)
			return strDOFnum[i];
		else if (dofType == DOFSurgeSway)
			return strDOFtext[i];
		else
			return strDOFxyz[i];
	}
			
	static const char *StrDOFAbrev(int i) {
		if (dofType == DOF123)
			return strDOFnum[i];
		else if (dofType == DOFSurgeSway)
			return strDOFtextAbrev[i];
		else
			return strDOFxyz[i];
	}
	
	static const char *StrDOF(int i, bool abrev) {
		if (abrev)
			return StrDOFAbrev(i);
		else
			return StrDOF(i);
	}
	
	static int DOFStr(const String &str) {
		for (int i = 0; i < 6; ++i)
			if (StrDOF(i) == ToLower(str))
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
			if (StrDOF(i) == ToLower(str))
				return i;
		return -1;
	}
	static const char *strDOFtext[];

private:
	static const char *strDOFnum[];
	static const char *strDOFxyz[];
	static const char *strDOFtextAbrev[];
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
	bool Calc(double x, double y, double rho_water);
	String Test();
	void Jsonize(JsonIO &json);
	
	struct LineType : Moveable<LineType> {
		String name;
		double mass, diameter, bl;
		void Jsonize(JsonIO &json);
	};
	UVector<LineType> lineTypes;
	
	LineType &GetLineType(String name);
	
	struct LineProperty : MoveableAndDeepCopyOption<LineProperty> {
		LineProperty() {}
		LineProperty(const LineProperty &line, int) {
			name = line.name;
			nameType = line.nameType;
			length = line.length;
			from = line.from;
			to = line.to;
			status = line.status;
			lenonfloor = line.lenonfloor;
			theta = line.theta;
			fanchorvessel = line.fanchorvessel;
			fVanchor = line.fVanchor;
			fVvessel = line.fVvessel;
			x = clone(line.x);
			y = clone(line.y);
			z = clone(line.z);
		}
		
		String name, nameType;
		double length;
		String from, to;
		void Jsonize(JsonIO &json);
		
		MooringStatus status;		// Obtained with Calc()
		double lenonfloor;
		UVector<double> x, y, z;
		double theta, fanchorvessel, fVanchor, fVvessel;
	};
	UVector<LineProperty> lineProperties;	

	struct Connection : Moveable<Connection> {
		String name;
		int type;
		double x, y, z;
		void Jsonize(JsonIO &json);
	};
	UVector<Connection> connections;	
	
	Connection &GetConnection(String name);
	
	double depth;
};
	
String FormatDoubleEmpty(double val);
String FormatIntEmpty(int val);


String GetFASTVar(const String &strFile, String varName, String paragraph = "");
void SetFASTVar(String &strFile, String varName, String value, String paragraph = "");
void GetFASTMatrixIds(const String &strFile, String var, int row, int col, int &posIni, int &posEnd);
double GetFASTMatrixVal(const String &strFile, String var, int row, int col);
MatrixXd GetFASTMatrix(const String &strFile, String var, int rows, int cols);
UVector<UVector<String>> GetFASTArray(const String &strFile, String var, String paragraph = "");	

	
#endif
