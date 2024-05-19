// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#ifndef _BEM_Rosetta_BEM_Rosetta_h_
#define _BEM_Rosetta_BEM_Rosetta_h_

#include <Core/Core.h>
#include <SysInfo/SysInfo.h>
#include <ScatterDraw/DataSource.h>
#include <Surface/Surface.h>
#include <STEM4U/Mooring.h>
#include <STEM4U/Utility.h>

using namespace Upp;

enum Symmetry {SYM_NO, SYM_XZ, SYM_YZ, SYM_XZ_YZ, SYM_AXISYMMETRIC};

class BEM;
BEM &Bem();

bool ConsoleMain(const UVector<String>& command, bool gui, Function <bool(String, int pos)> Status);
void SetBuildInfo(String &str);
String GetSystemInfo();

bool PrintStatus(String s, int d);

//class HydroClass; 

class Body : public DeepCopyOption<Body> {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	enum MESH_FMT {WAMIT_GDF,  WAMIT_DAT,  NEMOH_DAT,  NEMOHFS_DAT,   NEMOH_PRE,      AQWA_DAT,  AQWA_LIS, HAMS_PNL,  STL_BIN,     STL_TXT,   EDIT,  MSH_TDYN,   BEM_MESH, DIODORE_DAT,   HYDROSTAR_HST,   ORCA_OWR, MIKE21_GRD, UNKNOWN, NUMMESH};	
	static const char *meshStr[];
	static const bool meshCanSave[];
	static const char *meshExt[];
	
	enum MESH_TYPE {MOVED, UNDERWATER, ALL};
	
	Body() {}
	void Copy(const Body &msh);
	Body(const Body &msh, int) {Copy(msh);}
	
	bool IsEmpty() {return dt.mesh.IsEmpty();}
	
	const char *GetBodyStr() const {
		return meshStr[dt.GetCode()];
	}
	static const char *GetBodyStr(MESH_FMT c) {
		if (c < 0 || c > UNKNOWN)
			return "Unknown";
		return meshStr[c];
	}
	
	static MESH_FMT GetCodeBodyStr(String fmt) {
		fmt = ToLower(Trim(fmt));
		for (int i = 0; i < NUMMESH; ++i)
			if (fmt == ToLower(meshStr[i]))
				return static_cast<MESH_FMT>(i);
		return UNKNOWN;
	}

	static String Load(Body &mesh, String file, double rho, double g, bool cleanPanels, double grid, double eps);
	static String Load(Body &mesh, String file, double rho, double g, bool cleanPanels, double grid, double eps, bool &y0z, bool &x0z);
	static String Load(UArray<Body> &mesh, String file, double rho, double g, bool cleanPanels, double grid, double eps);
	static String Load(UArray<Body> &mesh, String file, double rho, double g, bool cleanPanels, double grid, double eps, bool &y0z, bool &x0z);
	
	String Heal(bool basic, double rho, double g, double grid, double eps, Function <bool(String, int pos)> Status);
	void Orient();
	void Append(const Surface &orig, double rho, double g);
	void Image(int axis);
	void Move(double dx, double dy, double dz, double ax, double ay, double az, 
			  double rho, double g, bool setnewzero);
	void Move(const double *pos, double rho, double g, bool setnewzero);
	void Move(const float *pos, double rho, double g, bool setnewzero);
		
	void AfterLoad(double rho, double g, bool onlyCG, bool isFirstTime, bool massBuoy = true);
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
	
	static void SaveAs(const UArray<Body> &meshes, String fileName, MESH_FMT type, MESH_TYPE meshType, double rho, double g, bool symX, bool symY, int &nNodes, int &nPanels);
	static void SaveAs(const UArray<Body> &meshes, String fileName, MESH_FMT type, MESH_TYPE meshType, double rho, double g, bool symX, bool symY) {
		int nNodes, nPanels;
		SaveAs(meshes, fileName, type, meshType, rho, g, symX, symY, nNodes, nPanels);
	}
	static void SaveAs(const Body &mesh, String fileName, MESH_FMT type, MESH_TYPE meshType, double rho, double g, bool symX, bool symY, int &nNodes, int &nPanels) {
		UArray<Body> meshes;
		meshes.Add(clone(mesh));
		SaveAs(meshes, fileName, type, meshType, rho, g, symX, symY, nNodes, nPanels);
	}
	static void SaveAs(const Body &mesh, String fileName, MESH_FMT type, MESH_TYPE meshType, double rho, double g, bool symX, bool symY) {
		int nNodes, nPanels;
		SaveAs(mesh, fileName, type, meshType, rho, g, symX, symY, nNodes, nPanels);
	}
	
	void SetMass(double m);
	double GetMass() const		{return dt.M(0, 0);}
	
	void Report(double rho) const;
	
	bool IsSymmetricX();
	bool IsSymmetricY();
	
	void Jsonize(JsonIO &json);
		
	// All Body data is here
	class Data {
	public:
		Point3D projectionPos = Null, projectionNeg = Null;
		Pointf cgZ0surface = Null;
		Point3D cb = Null;
		Point3D cg = Null, cg0 = Null, c0 = Null;
		double Vo = Null;   				// Displaced volume
		MatrixXd M,
				 Dlin, Dquad, 
				 C, Cmoor, Cadd,
				 Aadd;
		
		String name;
		String fileName;
		String fileHeader;
		String lidFile;
		
		Surface mesh, under, mesh0;
		
		void SetId(int iid)					{id = iid;}
		int  GetId() const					{return id;}
		
		void	 SetCode(MESH_FMT _code)	{code = _code;}
		MESH_FMT GetCode() const			{return code;}
			
	private:
		MESH_FMT code = UNKNOWN;
		int id = -1;
	};
	Data dt;
	
	static void ResetIdCount()	{idCount = 0;}
	void IncrementIdCount()		{dt.SetId(idCount++);}
	
	static int idCount;
};

class Hydro : public DeepCopyOption<Hydro> {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	enum BEM_FMT {WAMIT, 		  WAMIT_1_3, 					FAST_WAMIT, 				 	HAMS_WAMIT, HAMS, WADAM_WAMIT,   NEMOH,   NEMOHv115, NEMOHv3, SEAFEM_NEMOH,   AQWA,   					  FOAMM,   DIODORE,		BEMROSETTA, 	   ORCAWAVE,   CSV_MAT,    CSV_TABLE,    BEMIOH5,		CAPYTAINE, HYDROSTAR, CAPYNC, 
#ifdef PLATFORM_WIN32	
	ORCAWAVE_OWR, 
#endif
	UNKNOWN, NUMBEM};
	static const char *bemStr[];
	static const bool bemCanSave[];
	static const bool caseCanSave[];
	static const char *bemExt[];
	
	static const char *GetBemStr(BEM_FMT c) {
		if (c < 0 || c > UNKNOWN)
			return "Unknown";
		return bemStr[c];
	}
	
	static BEM_FMT GetCodeBemStr(String fmt) {
		fmt = ToLower(Trim(fmt));
		for (int i = 0; i < NUMBEM; ++i)
			if (fmt == ToLower(bemStr[i]))
				return static_cast<BEM_FMT>(i);
		return UNKNOWN;
	}
	
	//Hydro &hd() 			{return *this;}
	//const Hydro &hd() const {return *this;}
	
	void Copy(const Hydro &hyd);
	Hydro(const Hydro &hyd, int) {Copy(hyd);}
	void SaveAs(String file, Function <bool(String, int)> Status = Null, BEM_FMT type = UNKNOWN, int qtfHeading = Null);
	void Report() const;
	Hydro() {}
	virtual ~Hydro() noexcept {}	
	
	const char *GetCodeStr()	const {
		switch (dt.solver) {
		case WAMIT: 		return t_("Wamit");
		case WAMIT_1_3: 	return t_("Wamit.1.2.3");
		case WADAM_WAMIT: 	return t_("Wadam.Wamit");
		case FAST_WAMIT: 	return t_("Wamit.FAST");
		case NEMOH:			return t_("Nemoh v2");
		case NEMOHv115:		return t_("Nemoh v115");
		case NEMOHv3:		return t_("Nemoh v3");
		case SEAFEM_NEMOH:	return t_("SeaFEM.Nemoh");
		case AQWA:			return t_("AQWA");
		case FOAMM:			return t_("FOAMM");
		case BEMROSETTA:	return t_("BEMRosetta");
		case DIODORE:		return t_("Diodore");
		case ORCAWAVE:		return t_("OrcaWave");
		case CSV_MAT:		return t_("CSV.mat");
		case CSV_TABLE:		return t_("CSV.tab");
		case BEMIOH5:		return t_("BEMIO.h5");
		case HAMS_WAMIT:	return t_("HAMS.1.2.3");
		case HAMS:			return t_("HAMS");
		case CAPYTAINE:		return t_("Capytaine");
		case HYDROSTAR:		return t_("Hydrostar.out");
		case CAPYNC:		return t_("Capytaine.nc");
#ifdef PLATFORM_WIN32	
		case ORCAWAVE_OWR: 	return t_("OrcaWave.owr");
#endif
		case UNKNOWN:		return t_("Unknown");
		case NUMBEM:		NEVER();
		}
		return t_("Unknown");
	}
	
	const char *GetCodeStrAbr() const {
		switch (dt.solver) {
		case WAMIT: 		return t_("Wm.o");
		case WAMIT_1_3: 	return t_("Wm.1");
		case WADAM_WAMIT: 	return t_("WDM");
		case FAST_WAMIT: 	return t_("FST");
		case NEMOH:			return t_("Nmh");
		case NEMOHv115:		return t_("Nmh115");
		case NEMOHv3:		return t_("Nmh3");
		case SEAFEM_NEMOH:	return t_("SFM");
		case AQWA:			return t_("AQW");
		case FOAMM:			return t_("FMM");
		case BEMROSETTA:	return t_("BMR");
		case DIODORE:		return t_("DIO");
		case ORCAWAVE:		return t_("ORC");
		case CSV_MAT:		return t_("CSVm");
		case CSV_TABLE:		return t_("CSVt");
		case BEMIOH5:		return t_("BMh5");
		case HAMS_WAMIT:	return t_("HAMS_W");
		case HAMS:			return t_("HAMS");
		case CAPYTAINE:		return t_("Capy");
		case CAPYNC:		return t_("Capy.nc");
#ifdef PLATFORM_WIN32	
		case ORCAWAVE_OWR: 	return t_("ORC.owr");
#endif
		case HYDROSTAR:		return t_("Hydr");
		case UNKNOWN:		return t_("¿?");
		case NUMBEM:		NEVER();
		}
		return t_("Unknown");
	}
		
	inline bool IsAvailableDOF(int ib, int idf) {
		/*if (dof.IsEmpty())
			return false;*/
		
		int i = ib*6+idf;
		
		if (dt.Ainf.size() > 0)
			if (IsNum(dt.Ainf(i, i)))
				return true;
		
		if (!dt.A.IsEmpty())
			if (dt.A[i][i].size() > 0)
				if (IsNum(dt.A[i][i][0]))
					return true;
		
		return false;			   
	}
	
	bool IsNemoh() {return dt.solver == CAPYTAINE || dt.solver == NEMOH || dt.solver == NEMOHv115 || dt.solver == NEMOHv3;}
	
    double GetMass(int ib) const	{return dt.msh[ib].GetMass();}
    
    int GetHeadId(double hd) const;
    int GetHeadIdMD(const std::complex<double> &h) const;
	
	void BeforeSaveCase(String folderBase, int numCases, bool deleteFolder) const;
		
    /*struct Forces : public DeepCopyOption<Forces> {
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
    };*/
    
    enum FORCE {NONE, ALL, SCATTERING, FK, QTFSUM, QTFDIF};
    
	//typedef UArray<MatrixXcd> Forces;  	 	// [Nh](Nf, 6*Nb) 	
	typedef UArray<UArray<MatrixXcd>> Forces;  	// [Nb][Nh](Nf, 6) 	
  	typedef Forces RAO;
    
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
         
    static void Initialize_QTF(UArray<UArray<UArray<MatrixXcd>>> &qtf, int nb, int nh, int nf) {
        qtf.SetCount(nb);
        for (int ib = 0; ib < nb; ++ib) {
            qtf[ib].SetCount(nh);
        	for (int ih = 0; ih < nh; ++ih) {    
        		qtf[ib][ih].SetCount(6);
        		for (int idf = 0; idf < 6; ++idf) 
        			qtf[ib][ih][idf].setConstant(nf, nf, NaNDouble);
        	}
        }
    }
    
    double GetQTFVal(int ib, int idof, int idh, int ifr1, int ifr2, bool isSum, char what, bool getDim) const;
    MatrixXd GetQTFMat(int ib, int idof, int idh, bool isSum, char what, bool getDim) const;
    
    static void Initialize_MD(UArray<UArray<UArray<VectorXd>>> &md, int nb, int nh, int nf) {
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
    
    void Dimensionalize();
    void Normalize();
    
    static String C_units(int i, int j);
    
    void SetC(int ib, const MatrixXd &K);
	
	String AfterLoad(Function <bool(String, int)> Status = Null);
	
	void Initialize_AB(UArray<UArray<VectorXd>> &a, double val = NaNDouble);
	void Initialize_ABpan(UArray<UArray<UArray<UArray<UArray<double>>>>> &a, double val = NaNDouble);
	
	void Initialize_Forces();
	void Initialize_Forces(Forces &f, int _Nh = -1, double val = NaNDouble);
	void Normalize_Forces(Forces &f);
	void Normalize_RAO(RAO &f);
	void Dimensionalize_Forces(Forces &f);
	void Dimensionalize_RAO(RAO &f);
	void Add_Forces(Forces &to, const Hydro &hy, const Forces &from);
	void Add_RAO(RAO &to, const Hydro &hy, const RAO &from);
	void Symmetrize();
	void GetFexFromFscFfk();
	void GetFscFromFexFfk();
	void GetFfkFromFexFsc();
	void CompleteForces1st();
	
	bool SymmetryRule(int idf6, bool xAxis);
	void Symmetrize_Forces(bool xAxis);
	void Symmetrize_QTF(bool xAxis);
	void Symmetrize_MD(bool xAxis);
	
	void Initialize_PotsRad();
	void Initialize_PotsIncDiff(UArray<UArray<UArray<UArray<std::complex<double>>>>> &pots);
	
	void Initialize_Sts();
	
	void Average(const UArray<Hydro> &hydros, const UVector<int> &ids);
	void Converge(const UArray<Hydro> &hydros, const UVector<int> &ids);
		
	void SortFrequencies();
	void SortHeadings();
	
	void MapNodes(int ib, UVector<Point3D> &points, Tensor<double, 4> &Apan, Tensor<double, 4> &Bpan) const;
	void SaveMap(String fileName, String type, int ifr, bool onlyDiagonal, const UVector<int> &ids, const UVector<Point3D> &points, 
				 const Tensor<double, 4> &Apan, const Tensor<double, 4> &Bpan) const;
	void SaveMap(Grid &g, int ifr, bool onlyDiagonal, const UVector<int> &ids, const UVector<Point3D> &points, 
		const Tensor<double, 4> &Apan, const Tensor<double, 4> &Bpan) const;
			
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
	
	static int GetK_RAO(int i) {
		while (i > 5)
			i -= 6;
		if (i < 3)
			return 0;
		else
			return 1;
	}
	
	static int GetK_C(int i, int j) {
		if (i == 2 && j == 2)
			return 2;	
		else if (i < 3) 
			return 3;
		else
			return 4;
	}
	
	int GetIrregularHead() const;	
	int GetIrregularFreq() const;	
	
	double g_dim()		const;
	double g_ndim()		const;
	double rho_dim()	const;
	double rho_ndim()	const;
	double g_rho_dim()  const;
	double g_rho_ndim() const;
	
	VectorXd Get_w() 							const {return Map<const VectorXd>(dt.w, dt.w.size());}
	VectorXd Get_T()		 					const {
		VectorXd ret(dt.w.size());
		for (int i = 0; i < ret.size(); ++i)
			ret[i] = 2*M_PI/dt.w[i];
		return ret;
	}
	
	double A_dim(int ifr, int idf, int jdf) 	const {return dt.dimen  ? dt.A[idf][jdf][ifr]*rho_dim()/rho_ndim()  : dt.A[idf][jdf][ifr]*(rho_dim()*pow(dt.len, GetK_AB(idf, jdf)));}
	VectorXd A_dim(int idf, int jdf) 			const {return dt.dimen  ? dt.A[idf][jdf]    *(rho_dim()/rho_ndim()) : dt.A[idf][jdf]*     (rho_dim()*pow(dt.len, GetK_AB(idf, jdf)));}
	double A_ndim(int ifr, int idf, int jdf) 	const {return !dt.dimen ? dt.A[idf][jdf][ifr]/**(rho_ndim()/rho_dim())*/ : dt.A[idf][jdf][ifr]/(rho_ndim()*pow(dt.len, GetK_AB(idf, jdf)));}
	VectorXd A_ndim(int idf, int jdf)			const {return !dt.dimen ? dt.A[idf][jdf]/**(rho_ndim()/rho_dim())*/ : dt.A[idf][jdf]*(1/(rho_ndim()*pow(dt.len, GetK_AB(idf, jdf))));}
	double A_(bool ndim, int ifr, int idf, int jdf) const {return ndim ? A_ndim(ifr, idf, jdf) : A_dim(ifr, idf, jdf);}
	MatrixXd A_mat(bool ndim, int ifr, int ib1, int ib2) 	const;
	
	double A_P_dim(int ifr, int idf, int jdf) 	const {return dt.dimen  ? dt.A_P[idf][jdf][ifr]*rho_dim()/rho_ndim()  : dt.A_P[idf][jdf][ifr]*(rho_dim()*pow(dt.len, GetK_AB(idf, jdf)));}
	double A_P_ndim(int ifr, int idf, int jdf) 	const {return !dt.dimen ? dt.A_P[idf][jdf][ifr]/**(rho_ndim()/rho_dim())*/ : dt.A_P[idf][jdf][ifr]/(rho_ndim()*pow(dt.len, GetK_AB(idf, jdf)));}
	double A_P_(bool ndim, int ifr, int idf, int jdf) const {return ndim ? A_P_ndim(ifr, idf, jdf) : A_P_dim(ifr, idf, jdf);}
	
	double A0_dim(int idf, int jdf)   		 	const {return dt.dimen  ? dt.A0(idf, jdf)*rho_dim()/rho_ndim() : dt.A0(idf, jdf)  *(rho_dim()*pow(dt.len, GetK_AB(idf, jdf)));}
	double A0_ndim(int idf, int jdf)  		 	const {return !dt.dimen ? dt.A0(idf, jdf)      : dt.A0(idf, jdf)  /(rho_ndim()*pow(dt.len, GetK_AB(idf, jdf)));}
	double A0_(bool ndim, int idf, int jdf) 	const {return ndim   ? A0_ndim(idf, jdf) : A0_dim(idf, jdf);}
	double Ainf_dim(int idf, int jdf) 		 	const {return dt.dimen  ? dt.Ainf(idf, jdf)*rho_dim()/rho_ndim() : dt.Ainf(idf, jdf)*(rho_dim()*pow(dt.len, GetK_AB(idf, jdf)));}
	MatrixXd Ainf_mat(bool ndim, int ib1, int ib2) const;
	double Ainf_ndim(int idf, int jdf)		 	const {return !dt.dimen ? dt.Ainf(idf, jdf) : dt.Ainf(idf, jdf)/(rho_ndim()*pow(dt.len, GetK_AB(idf, jdf)));}
	double Ainf_(bool ndim, int idf, int jdf) 	const {return ndim   ? Ainf_ndim(idf, jdf) : Ainf_dim(idf, jdf);}
	
	double B_dim(int ifr, int idf, int jdf)  	const {return dt.dimen  ? dt.B[idf][jdf][ifr]*rho_dim()/rho_ndim() : dt.B[idf][jdf][ifr]*(rho_dim()*pow(dt.len, GetK_AB(idf, jdf))*dt.w[ifr]);}
	VectorXd B_dim(int idf, int jdf)  	   		const;
	double B_ndim(int ifr, int idf, int jdf) 	const {return !dt.dimen ? dt.B[idf][jdf][ifr]/**(rho_ndim()/rho_dim())*/ : dt.B[idf][jdf][ifr]/(rho_ndim()*pow(dt.len, GetK_AB(idf, jdf))*dt.w[ifr]);}
	VectorXd B_ndim(int idf, int jdf) 	   		const;
	double B_(bool ndim, int ifr, int idf, int jdf)const {return ndim ? B_ndim(ifr, idf, jdf) : B_dim(ifr, idf, jdf);}	
	MatrixXd B_mat(bool ndim, int ifr, int ib1, int ib2) 	const;
	
	double B_H_dim(int ifr, int idf, int jdf)  	const {return dt.dimen  ? dt.B_H[idf][jdf][ifr]*rho_dim()/rho_ndim() : dt.B_H[idf][jdf][ifr]*(rho_dim()*pow(dt.len, GetK_AB(idf, jdf))*dt.w[ifr]);}
	double B_H_ndim(int ifr, int idf, int jdf) 	const {return !dt.dimen ? dt.B_H[idf][jdf][ifr]/**(rho_ndim()/rho_dim())*/ : dt.B_H[idf][jdf][ifr]/(rho_ndim()*pow(dt.len, GetK_AB(idf, jdf))*dt.w[ifr]);}
	double B_H_(bool ndim, int ifr, int idf, int jdf)const {return ndim ? B_H_ndim(ifr, idf, jdf) : B_H_dim(ifr, idf, jdf);}	
	
	double B_P_dim(int ifr, int idf, int jdf)  	const {return dt.dimen  ? dt.B_P[idf][jdf][ifr]*rho_dim()/rho_ndim() : dt.B_P[idf][jdf][ifr]*(rho_dim()*pow(dt.len, GetK_AB(idf, jdf))*dt.w[ifr]);}
	double B_P_ndim(int ifr, int idf, int jdf) 	const {return !dt.dimen ? dt.B_P[idf][jdf][ifr]/**(rho_ndim()/rho_dim())*/ : dt.B_P[idf][jdf][ifr]/(rho_ndim()*pow(dt.len, GetK_AB(idf, jdf))*dt.w[ifr]);}
	double B_P_(bool ndim, int ifr, int idf, int jdf)const {return ndim ? B_P_ndim(ifr, idf, jdf) : B_P_dim(ifr, idf, jdf);}	
	
	double Kirf_dim(int it, int idf, int jdf)  	   	  const {return dt.dimen ? dt.Kirf[idf][jdf][it]*g_rho_dim()/g_rho_ndim()  : dt.Kirf[idf][jdf][it]*(g_rho_dim()*pow(dt.len, GetK_F(idf)));}
	double Kirf_ndim(int it, int idf, int jdf) 	   	  const {return !dt.dimen ? dt.Kirf[idf][jdf][it] : dt.Kirf[idf][jdf][it]/(g_rho_ndim()*pow(dt.len, GetK_F(idf)));}
	VectorXd Kirf_ndim(int idf, int jdf) 	 		  const {return !dt.dimen ? dt.Kirf[idf][jdf]     : dt.Kirf[idf][jdf]/(g_rho_ndim()*pow(dt.len, GetK_F(idf)));}
	double Kirf_(bool ndim, int it, int idf, int jdf) const {return ndim ? Kirf_ndim(it, idf, jdf) : Kirf_dim(it, idf, jdf);}
	
	double Ainf_w_dim(int ifr, int idf, int jdf)   const {return dt.dimen  ? dt.Ainf_w[idf][jdf][ifr]*rho_dim()/rho_ndim() : dt.Ainf_w[idf][jdf][ifr]*(rho_dim()*pow(dt.len, GetK_AB(idf, jdf)));}
	VectorXd Ainf_w_dim(int idf, int jdf)		   const {return dt.dimen  ? dt.Ainf_w[idf][jdf]    *(rho_dim()/rho_ndim()) : dt.Ainf_w[idf][jdf]*     (rho_dim()*pow(dt.len, GetK_AB(idf, jdf)));}
	double Ainf_w_ndim(int ifr, int idf, int jdf)  const {return !dt.dimen ? dt.Ainf_w[idf][jdf][ifr]/**(rho_ndim()/rho_dim())*/ : dt.Ainf_w[idf][jdf][ifr]/(rho_ndim()*pow(dt.len, GetK_AB(idf, jdf)));}
	VectorXd Ainf_w_ndim(int idf, int jdf)		   const {return !dt.dimen ? dt.Ainf_w[idf][jdf]/**(rho_ndim()/rho_dim())*/ : dt.Ainf_w[idf][jdf]*(1/(rho_ndim()*pow(dt.len, GetK_AB(idf, jdf))));}
	double Ainf_w_(bool ndim, int ifr, int idf, int jdf)const {return ndim   ? Ainf_w_ndim(ifr, idf, jdf) : Ainf_w_dim(ifr, idf, jdf);}
	
	double C_dim(int ib, int idf, int jdf)   	   const {return dt.dimen  ? dt.msh[ib].dt.C(idf, jdf)/**g_rho_dim()/g_rho_ndim()*/  : dt.msh[ib].dt.C(idf, jdf)*(g_rho_dim()*pow(dt.len, GetK_C(idf, jdf)));}
	MatrixXd C_(bool ndim, int ib) 				   const;
	void C_dim();	
	double C_ndim(int ib, int idf, int jdf)  	   const {return !dt.dimen ? dt.msh[ib].dt.C(idf, jdf)  : dt.msh[ib].dt.C(idf, jdf)/(g_rho_ndim()*pow(dt.len, GetK_C(idf, jdf)));}
	double C_(bool ndim, int ib, int idf, int jdf) const {return ndim ? C_ndim(ib, idf, jdf) : C_dim(ib, idf, jdf);}

	double CMoor_dim(int ib, int idf, int jdf)     const {return dt.dimen  ? dt.msh[ib].dt.Cmoor(idf, jdf)/**g_rho_dim()/g_rho_ndim()*/  : dt.msh[ib].dt.Cmoor(idf, jdf)*(g_rho_dim()*pow(dt.len, GetK_C(idf, jdf)));}
	MatrixXd CMoor_(bool ndim, int ib) 			   const;
	void CMoor_dim();	
	double CMoor_ndim(int ib, int idf, int jdf)    const {return !dt.dimen ? dt.msh[ib].dt.Cmoor(idf, jdf)  : dt.msh[ib].dt.Cmoor(idf, jdf)/(g_rho_ndim()*pow(dt.len, GetK_C(idf, jdf)));}
	double CMoor_(bool ndim, int ib, int idf, int jdf)const {return ndim ? CMoor_ndim(ib, idf, jdf) : CMoor_dim(ib, idf, jdf);}

	// MD
	
	double Md_dim(int idf, int ih, int ifr)  const {
		int ib = int(idf/6);
		idf -= 6*ib;
		return dt.dimen ? dt.md[ib][ih][idf](ifr)*g_rho_ndim()/g_rho_dim()  : dt.md[ib][ih][idf](ifr)*(g_rho_dim()*pow(dt.len, GetK_F(idf)));}
	double Md_ndim(int idf, int ih, int ifr) const {
		int ib = int(idf/6);
		idf -= 6*ib;
		return !dt.dimen ? dt.md[ib][ih][idf](ifr) : dt.md[ib][ih][idf](ifr)/(g_rho_ndim()*pow(dt.len, GetK_F(idf)));}
	double Md_(bool ndim, int idf, int ih, int ifr) const {return ndim ? Md_ndim(idf, ih, ifr) : Md_dim(idf, ih, ifr);}	// idf: body, jdf: heading, [Nb][Nh][6](Nf)
	
	VectorXd Md_dof(bool ndim, int _h, int idf) const;
	
	MatrixXd Dlin_dim(int ib) const;
	MatrixXd Dquad_dim(int ib) const;
	
	MatrixXcd QTF_dof(bool ndim, bool isSum, int _h, int idf, int ib) const;
		
	// F
		
	std::complex<double> F_dim(const Forces &f, int _h, int ifr, int idf, int ib)  const {
		return dt.dimen ? f[ib][_h](ifr, idf)*g_rho_ndim()/g_rho_dim()  : f[ib][_h](ifr, idf)*(g_rho_dim()*pow(dt.len, GetK_F(idf)));}
	
	void F_dim(Forces &f);
	
	std::complex<double> F_ndim(const Forces &f, int _h, int ifr, int idf, int ib) const {
		if (!dt.dimen) 
			return f[ib][_h](ifr, idf); 
		else 
			return f[ib][_h](ifr, idf)/(g_rho_ndim()*pow(dt.len, GetK_F(idf)));
	}
	
	std::complex<double> F_(bool ndim, const Forces &f, int _h, int ifr, int idf, int ib) const {
		return ndim ? F_ndim(f, _h, ifr, idf, ib) : F_dim(f, _h, ifr, idf, ib);
	}
	
	template <class T>
	T F_dim(T f, int idf)  	      const {return  dt.dimen ? f*g_rho_ndim()/g_rho_dim() : f*(g_rho_dim()*pow(dt.len, GetK_F(idf)));}
	template <class T>
	T F_ndim(T f, int idf) 	      const {return !dt.dimen ? f : f/(g_rho_ndim()*pow(dt.len, GetK_F(idf)));}
	template <class T>
	T F_(bool ndim, T f, int idf) const {return   ndim ? F_ndim(f, idf) : F_dim(f, idf);}

	VectorXcd F_(bool ndim, const Forces &f, int _h, int ifr, int ib) const;
	VectorXcd F_dof(bool ndim, const Forces &f, int _h, int idf, int ib) const;

	// RAO
	
	std::complex<double> RAO_dim(const RAO &f, int _h, int ifr, int idf, int ib)  const {
		return (g_rho_ndim()/g_rho_dim())*(dt.dimen ? f[ib][_h](ifr, idf)  : f[ib][_h](ifr, idf)*pow(dt.len, GetK_RAO(idf)));}
	
	void RAO_dim(RAO &f);
	
	std::complex<double> RAO_ndim(const RAO &f, int _h, int ifr, int idf, int ib) const {
		return (g_rho_ndim()/g_rho_dim())*(!dt.dimen ? f[ib][_h](ifr, idf) : f[ib][_h](ifr, idf)/pow(dt.len, GetK_RAO(idf)));}
	
	std::complex<double> RAO_(bool ndim, const RAO &f, int _h, int ifr, int idf, int ib) const {
		return ndim ? RAO_ndim(f, _h, ifr, idf, ib) : RAO_dim(f, _h, ifr, idf, ib);
	}
	
	template <class T>
	T RAO_dim(T f, int idf) 	const {
		return (g_rho_ndim()/g_rho_dim())*(dt.dimen  ? f : f/pow(dt.len, GetK_RAO(idf)));
	}
	template <class T>
	T RAO_ndim(T f, int idf) 	const {return (g_rho_ndim()/g_rho_dim())*(!dt.dimen ? f : f*pow(dt.len, GetK_RAO(idf)));}
	template <class T>
	T RAO_(bool ndim, T f, int idf) const {return   ndim ? RAO_ndim(f, idf) : RAO_dim(f, idf);}

	VectorXcd RAO_(bool ndim, const RAO &f, int _h, int ifr, int ib) const;
	VectorXcd RAO_dof(bool ndim, const RAO &f, int _h, int idf, int ib) const;
	VectorXcd RAO_dof(bool ndim, int _h, int idf, int ib) const;
	
	// FOAMM
	
	inline std::complex<double> Z(bool ndim, int ifr, int idf, int jdf) const {
		return std::complex<double>(B_(ndim, ifr, idf, jdf), dt.w[ifr]*(A_(ndim, ifr, idf, jdf) - Ainf_(ndim, idf, jdf))/(!ndim ? 1. : dt.w[ifr]));
	}
	
	std::complex<double> TFS_dim(int ifr, int idf, int jdf) 		const {return dt.dimenSTS  ? dt.sts[idf][jdf].TFS[ifr]*g_rho_dim()/g_rho_ndim() : dt.sts[idf][jdf].TFS[ifr]*(rho_dim()*pow(dt.len, GetK_AB(idf, jdf))*dt.w[ifr]);}
	std::complex<double> TFS_ndim(int ifr, int idf, int jdf) 		const {return !dt.dimenSTS ? dt.sts[idf][jdf].TFS[ifr]/**g_rho_ndim()/g_rho_dim()*/ : dt.sts[idf][jdf].TFS[ifr]/(rho_ndim()*pow(dt.len, GetK_AB(idf, jdf))*dt.w[ifr]);}
	std::complex<double> TFS_(bool ndim, int ifr, int idf, int jdf) const {return ndim ? TFS_ndim(ifr, idf, jdf) : TFS_dim(ifr, idf, jdf);}
	
	double Tdof(int ib, int idf) const;
	double Tdofw(int ib, int idf) const;
	double Theave(int ib) const;
	double Theavew(int ib) const;
	double Troll(int ib) const;
	double Trollw(int ib) const;
	double Tpitch(int ib) const;
	double Tpitchw(int ib) const;
	double GM(int ib, int idf) const;
	double GMroll(int ib) const;
	double GMpitch(int ib) const;
	
	void CheckNaN();
		
	void Jsonize(JsonIO &json);
	
	// All Hydro data is here
	class Data {
	public:
		String file;        			// BEM output file name
		String name;
	    double g = Null;        	   	// gravity
	    double h = Null;    	       	// water depth
	   	double rho = Null;        		// density
	   	double len;						// Length scale
	   	int dimen = Null;				// false if data is dimensionless
	    int Nb = Null;          		// number of bodies
	    int Nf = Null;          		// number of wave frequencies
	    int Nh = Null;          		// number of wave headings
	 	
		UArray<UArray<VectorXd>> A;		// [6*Nb][6*Nb][Nf]	Added mass
		UArray<UArray<VectorXd>> Ainf_w;// [6*Nb][6*Nb][Nf]	Infinite frequency added mass (w)
		UArray<UArray<VectorXd>> A_P;	// [6*Nb][6*Nb][Nf]	Added mass obtained through potentials
	    MatrixXd Ainf;        			// (6*Nb, 6*Nb) 	Infinite frequency added mass
	    MatrixXd A0;        			// (6*Nb, 6*Nb)  	Infinite period added mass
	
	    UArray<UArray<VectorXd>> B; 	// [6*Nb][6*Nb][Nf]	Radiation damping
	    UArray<UArray<VectorXd>> B_H; 	// [6*Nb][6*Nb][Nf]	Radiation damping obtained through Haskind
	    UArray<UArray<VectorXd>> B_P; 	// [6*Nb][6*Nb][Nf]	Radiation damping obtained through potentials
	    
	    UVector<double> head;			// [Nh]             Wave headings (deg)

	    double x_w = Null, y_w = Null;	// 					Wave centre
	    BEM_FMT solver;        			// BEM_FMT			BEM code 
	    
	    UArray<UArray<VectorXd>> Kirf;	// [6*Nb][6*Nb][Nt]	Radiation impulse response function IRF
	    VectorXd Tirf;	  				// [Nt]				Time-window for the calculation of the IRF
	    
	    Forces ex; 						// Excitation
	    Forces sc;			 			// Diffraction scattering
	    Forces fk; 						// Froude-Krylov
	    Forces fk_pot; 					// Froude-Krylov obtained through potentials
	   
	   	RAO rao;
	    
	    String description;
	
	    UArray<UArray<StateSpace>> sts;	// (6*Nb, 6*Nb)		State space data
	    int dimenSTS;					// false if data is dimensionless
	    String stsProcessor;
	    
	    VectorXd  qw;		 			// [Nf]             Wave frequencies
	    VectorXcd qh;					// [Nh]             Wave headings
	    UArray<UArray<UArray<MatrixXcd>>> qtfsum, qtfdif;// [Nb][Nh][6](Nf, Nf)	
	    bool qtfdataFromW = true;
	    int qtftype = 0;				// 7. Control surface, 8. Momentum conservation/Far field, 9. Pressure integration/Near field	
		
		/* The priority is:
											Extension	DOF				Accuracy	Run time	Control surface
		Control surface 					.7			1,2,3,4,5,6		Better		More		Yes
		Pressure integration/Near field		.9			1,2,3,4,5,6		Less		Less		No
		Momentum conservation/Far field		.8			1,2,6									No				*/
		        
	    VectorXcd mdhead;						// [Nh]             Wave headings
		UArray<UArray<UArray<VectorXd>>> md;	// [Nb][Nh][6](Nf)	
		int mdtype = 0;							// 7. Control surface, 8. Momentum conservation/Far field, 9. Pressure integration/Near field		
	    					
	    UVector<double> w;		 				// [Nf]             Wave frequencies
	    int dataFromW = Null;
	    		
	   	UArray<Body> msh;						// [Nb]
	   	
	   	UArray<UArray<UArray<UArray<std::complex<double>>>>> pots_rad;	// [Nb][Np][6][Nf]		Radiation complex potentials
	   	UArray<UArray<UArray<UArray<std::complex<double>>>>> pots_diff;	// [Nb][Np][Nh][Nf]		Diffraction complex potentials
	   	UArray<UArray<UArray<UArray<std::complex<double>>>>> pots_inc;	// [Nb][Np][Nh][Nf]		Incident complex potentials
	   	
	   	Tensor<double, 5> Apan;		// [Nb][Np][6][6][Nf]	Added mass
	   	Tensor<double, 5> Bpan;		// [Nb][Np][6][6][Nf]	Radiation damping
	   	
	   	bool symX = false, symY = false;
		
		void SetId(int _id)			{id = _id;}
		int GetId()	const			{return id;}
	
	private:
		int id = -1;
	};
	Data dt;
	
	static int idCount;
	
	static void ResetIdCount()	{idCount = 0;}
	void IncrementIdCount()		{dt.SetId(idCount++);}

	static const char *strDataToPlot[];
	static String C_units_base(int i, int j);
		
	static void GetOldAB(const UArray<MatrixXd> &oldAB, UArray<UArray<VectorXd>> &AB);
	static void SetOldAB(UArray<MatrixXd> &oldAB, const UArray<UArray<VectorXd>> &AB);
	
	void ResetForces1st(Hydro::FORCE force);
	
	void SaveForce(FileOut &out, Hydro::Forces &f);
	void SaveMD(FileOut &out);
	void SaveC(FileOut &out);
	void SaveM(FileOut &out);
	
public:
	String LoadSerialization(String file);
	void SaveSerialization(String file);
	
	static int LoadHydro(UArray<Hydro> &hydro, String file, Function <bool(String, int)> Status);
	
	void LoadCase(String file, Function <bool(String, int)> Status);
	void SaveFolderCase(String folder, bool bin, int numCases, int numThreads, int solver) const;
	
	void SaveCSVMat(String file);
	void SaveCSVTable(String file);
	
	void SaveDiodoreHDB(String file);
	
	UVector<String> Check(BEM_FMT type) const;
	
	
	enum DataToShow {DATA_A, DATA_B, DATA_AINFW, DATA_KIRF, DATA_FORCE_SC, DATA_FORCE_FK, DATA_FORCE_EX, DATA_RAO, 
					 DATA_STS, DATA_STS2, DATA_MD, DATA_B_H, DATA_A_P, DATA_B_P, DATA_FORCE_FK_P};
	enum DataToPlot {PLOT_A, PLOT_AINF, PLOT_A0, PLOT_B, PLOT_AINFW, PLOT_KIRF, PLOT_FORCE_SC_1, PLOT_FORCE_SC_2,
				 PLOT_FORCE_FK_1, PLOT_FORCE_FK_2, PLOT_FORCE_EX_1, PLOT_FORCE_EX_2, 
				 PLOT_RAO_1, PLOT_RAO_2, PLOT_Z_1, PLOT_Z_2, PLOT_KR_1, PLOT_KR_2, 
				 PLOT_TFS_1, PLOT_TFS_2, PLOT_MD, PLOT_B_H, PLOT_A_P, PLOT_B_P, PLOT_FORCE_FK_1_P, PLOT_FORCE_FK_2_P};
	enum DataMatrix {MAT_K, MAT_A, MAT_DAMP_LIN, MAT_M, MAT_K2, MAT_DAMP_QUAD};
				 
	static const char *StrDataToPlot(DataToPlot dataToPlot) {
		return strDataToPlot[dataToPlot];
	}
	
	bool IsLoadedA	   (int i = 0, int j = 0)const {return dt.A.size() > i && dt.A[i].size() > j && dt.A[i][j].size() > 0 && IsNum(dt.A[i][j][0]);}
	bool IsLoadedAinf_w(int i = 0, int j = 0)const {return dt.Ainf_w.size() > i && dt.Ainf_w[i].size() > j && dt.Ainf_w[i][j].size() > 0 && IsNum(dt.Ainf_w[i][j][0]);}
	bool IsLoadedA_P   (int i = 0, int j = 0)const {return dt.A_P.size() > i && dt.A_P[i].size() > j && dt.A_P[i][j].size() > 0 && IsNum(dt.A_P[i][j][0]);}
	bool IsLoadedAinf  (int i = 0, int j = 0)const {return dt.Ainf.rows() > i && dt.Ainf.cols() > j && IsNum(dt.Ainf(i, j));}
	bool IsLoadedA0	   (int i = 0, int j = 0)const {return dt.A0.rows() > i && dt.A0.cols() > j && IsNum(dt.A0(i, j));}
	bool IsLoadedDlin(int ib = 0, int idf = 0, int jdf = 0)	const {return dt.msh[ib].dt.Dlin.size() > 0 && dt.msh[ib].dt.Dlin.rows() > idf && dt.msh[ib].dt.Dlin.cols() > jdf && IsNum(dt.msh[ib].dt.Dlin(idf, jdf));}
	bool IsLoadedDquad(int ib = 0, int idf = 0, int jdf = 0)const {return dt.msh[ib].dt.Dquad.size() > 0 && dt.msh[ib].dt.Dquad.rows() > idf && dt.msh[ib].dt.Dquad.cols() > jdf && IsNum(dt.msh[ib].dt.Dquad(idf, jdf));}
	bool IsLoadedB	   (int i = 0, int j = 0)const {return dt.B.size() > i && dt.B[i].size() > j && dt.B[i][j].size() > 0 && IsNum(dt.B[i][j][0]);}
	bool IsLoadedB_H   (int i = 0, int j = 0)const {return dt.B_H.size() > i && dt.B_H[i].size() > j && dt.B_H[i][j].size() > 0 && IsNum(dt.B_H[i][j][0]);}
	bool IsLoadedB_P   (int i = 0, int j = 0)const {return dt.B_P.size() > i && dt.B_P[i].size() > j && dt.B_P[i][j].size() > 0 && IsNum(dt.B_P[i][j][0]);}
	bool IsLoadedC(int ib = 0, int idf = 0, int jdf = 0)	const {return dt.msh[ib].dt.C.size() > 0 && dt.msh[ib].dt.C.rows() > idf && dt.msh[ib].dt.C.cols() > jdf && IsNum(dt.msh[ib].dt.C(idf, jdf));}
	bool IsLoadedCMoor(int ib = 0, int idf = 0, int jdf = 0)const {return dt.msh[ib].dt.Cmoor.size() > 0 && dt.msh[ib].dt.Cmoor.rows() > idf && dt.msh[ib].dt.Cmoor.cols() > jdf && IsNum(dt.msh[ib].dt.Cmoor(idf, jdf));}
	bool IsLoadedM(int ib = 0, int idf = 0, int jdf = 0)	const {return dt.msh[ib].dt.M.size() > 0 && dt.msh[ib].dt.M.rows() > idf && dt.msh[ib].dt.M.cols() > jdf && IsNum(dt.msh[ib].dt.M(idf, jdf));}
	
	bool IsLoadedFex(int idf = 0, int ih = 0, int ib = 0)const {return IsLoadedForce(dt.ex, idf, ih, ib);}
	bool IsLoadedFsc(int idf = 0, int ih = 0, int ib = 0)const {return IsLoadedForce(dt.sc, idf, ih, ib);}
	bool IsLoadedFfk(int idf = 0, int ih = 0, int ib = 0)const {return IsLoadedForce(dt.fk, idf, ih, ib);}
	bool IsLoadedFfk_P(int idf = 0, int ih = 0)const {return IsLoadedForce(dt.fk_pot, idf, ih);}
	bool IsLoadedRAO(int idf = 0, int ih = 0)const {return IsLoadedForce(dt.rao,idf, ih);}
	bool IsLoadedForce(const Forces &f, int idf = 0, int ih = 0, int ib = 0)
											 const {return f.size() > ib && f[ib].size() > ih && f[ib][ih].cols() > idf && IsNum(f[ib][ih](0, idf));}
	
	bool IsLoadedStateSpace()	  			 const {return !dt.sts.IsEmpty();}
	bool IsLoadedQTF(bool isSum) 			 const {return isSum ? !dt.qtfsum.IsEmpty() : !dt.qtfdif.IsEmpty();}
	bool IsLoadedMD(int ib = 0, int ih = 0)	 const {return dt.md.size() > ib && dt.md[ib].size() > ih && dt.md[ib][ih].size() == 6 && dt.md[ib][ih][0].size() > 0 && IsNum(dt.md[ib][ih][0](0));}
	static bool IsLoadedMD(const UArray<UArray<UArray<VectorXd>>> &mD, int ib = 0, int ih = 0) {return mD.size() > ib && mD[ib].size() > ih && mD[ib][ih].size() == 6 && mD[ib][ih][0].size() > 0 && IsNum(mD[ib][ih][0](0));}
	bool IsLoadedKirf(int idf=0,int jdf = 0) const {return dt.Kirf.size() > idf && dt.Kirf[idf].size() > jdf && dt.Kirf[idf][jdf].size() > 0 && IsNum(dt.Kirf[idf][jdf][0]);}
	
	bool IsLoadedPotsRad()					 const {return !dt.pots_rad.IsEmpty()  && !dt.pots_rad[0].IsEmpty()  && !dt.pots_rad[0][0].IsEmpty();}
	bool IsLoadedPotsDiff()					 const {return !dt.pots_diff.IsEmpty() && !dt.pots_diff[0].IsEmpty() && !dt.pots_diff[0][0].IsEmpty();}
	bool IsLoadedPotsInc()					 const {return !dt.pots_inc.IsEmpty()  && !dt.pots_inc[0].IsEmpty()  && !dt.pots_inc[0][0].IsEmpty();}
	
	void RemoveThresDOF_A(double thres);
	void RemoveThresDOF_B(double thres);
	void RemoveThresDOF_Force(Forces &f, double thres);
	
	void Compare_rho(Hydro &a);
	void Compare_g(Hydro &a);
	void Compare_h(Hydro &a);
	void Compare_w(Hydro &a);
	void Compare_head(Hydro &a);
	void Compare_Nb(Hydro &a);
	void Compare_A(const UArray<UArray<VectorXd>> &a);
	void Compare_B(const UArray<UArray<VectorXd>> &b);
	void Compare_C(Hydro &a);
	void Compare_cg(Hydro &a);
	void Compare_F(const Forces &a, const Forces &b, String type);
	
	int GetW0();
	void Get3W0(int &id1, int &id2, int &id3);
	void GetA0();
		
	void GetK_IRF(double maxT = 120, int numT = 1000);
	double GetK_IRF_MaxT() const;
	static double GetK_IRF_MaxT(const UVector<double> &w);
	void GetAinf();
	void GetAinf_w();
	void GetRAO();
	void GetB_H(int &num);
	static VectorXcd GetRAO(double w, const MatrixXd &Aw, const MatrixXd &Bw, const VectorXcd &Fwh, 
				const MatrixXd &C, const MatrixXd &M, const MatrixXd &D, const MatrixXd &D2);
	void InitAinf_w();
	void GetOgilvieCompliance(bool zremoval, bool thinremoval, bool decayingTail, UVector<int> &vidof, UVector<int> &vjdof);
	void GetTranslationTo(const MatrixXd &to);
	void GetWaveTo(double xto, double yto);
	String SpreadNegative(Function <bool(String, int)> Status);
	void AddWave(int ib, double dx, double dy);
	
	void DeleteFrequencies(const UVector<int> &idFreq);
	void DeleteFrequenciesQTF(const UVector<int> &idFreqQTF);
	void DeleteHeadings(const UVector<int> &idHead);
	void DeleteHeadingsMD(const UVector<int> &idHead);
	void DeleteHeadingsQTF(const UVector<int> &idHeadQTF);
	void ResetForces(Hydro::FORCE force, bool forceMD, Hydro::FORCE forceQtf);
	void MultiplyDOF(double factor, const UVector<int> &idDOF, bool a, bool b, bool diag, bool f, bool md, bool qtf);
	void SwapDOF(int ib1, int idof1, int ib2, int idof2);
	void SwapDOF(int ib1, int ib2);
	
	void SymmetrizeDOF();
	
	void FillFrequencyGapsABForces(bool zero, int maxFreq);
	void FillFrequencyGapsQTF(bool zero, int maxFreq);
	
	void FillFrequencyGapsABForcesZero();
	void FillFrequencyGapsQTFZero();
	
	void CopyQTF_MD();
	
	void Join(const UVector<Hydro *> &hydrosp);
	
	String S_g()	const {return !IsNum(dt.g)   ? S("-") : Format("%.3f", dt.g);}
	String S_h()	const {return !IsNum(dt.h)   ? S("-") : (dt.h < 0 ? S(t_("INFINITY")) : Format("%.1f", dt.h));}
	String S_rho() 	const {return !IsNum(dt.rho) ? S("-") : Format("%.3f", dt.rho);}
	String S_len() 	const {return !IsNum(dt.len) ? S("-") : Format("%.1f", dt.len);}
	
	static String K_units(bool ndim, int r, int c) {
		if (ndim) {
			if (r < 3 && c < 3)
				return "m²";
			if (r < 3)
				return "m³/rad";
			if (c < 3)
				return "m³";
			return "m⁴/rad";			
		} else {
			if (r < 3 && c < 3)
				return "N/m";
			if (r < 3)
				return "N/rad";
			if (c < 3)
				return "N";
			return "Nm/rad";
		}
	}
	static String Kirf_units(bool ndim, int r, int c) {
		return K_units(ndim, r, c);
	}
	static String B_units(bool ndim, int r, int c) {
		if (ndim) {
			if (r < 3 && c < 3)
				return "m³-rad/s";
			if (r < 3)
				return "m⁴/s";
			if (c < 3)
				return "m⁴/s/rad";
			return "m⁵/s";
		} else {
			if (r < 3 && c < 3)
				return "N/(m/s)";
			if (r < 3)
				return "N/(rad/s)";
			if (c < 3)
				return "N m/(m/s)";
			return "N m/(rad/s)";
		}
	}
	static String D2_units(bool ndim, int r, int c) {
		if (ndim) {
			if (r < 3 && c < 3)
				return "";
			if (r < 3)
				return "";
			if (c < 3)
				return "";
			return "";
		} else {
			if (r < 3 && c < 3)
				return "";
			if (r < 3)
				return "";
			if (c < 3)
				return "";
			return "";
		}
	}
	static String A_units(bool ndim, int r, int c) {
		if (ndim) {
			if (r < 3 && c < 3)
				return "m³";
			if (r < 3)
				return "m⁴/rad";
			if (c < 3)
				return "m⁴";
			return "m⁵/rad";			
		} else {
			if (r < 3 && c < 3)
				return "kg";
			if (r < 3)
				return "kg m/rad";
			if (c < 3)
				return "kg m";
			return "kg m²/rad";
		}
	}
	static String M_units(int r, int c) {
		if (r < 3 && c < 3)
			return "kg";
		else if (r >= 3 && c >= 3)
			return "kg-m²";
		else
			return "kg-m";
	}
	static String F_units(bool ndim, int r) {
		if (ndim) {
			if (r < 3)
				return "m²";
			return "m³";
		} else {
			if (r < 3)
				return "N/m";
			return "Nm/m";
		}
	}
	static String RAO_units(int r) {
		if (r < 3)
			return "m/m";
		return "rad/m";
	}
	static String MD_units(bool ndim, int r) {
		if (ndim) {
			if (r < 3)
				return t_("m");
			return t_("m2");
		} else {
			if (r < 3)
				return t_("N/m2");
			return t_("Nm/m2");
		}
	} 
};

//bool IsNum(const Hydro::Forces &f);

class NemohBody : public Body {
public:
	static String LoadDat(UArray<Body> &mesh, String fileName, bool &x0z);
	static String LoadDatFS(UArray<Body> &mesh, String fileName, bool &x0z);
	static void SaveDat(const UArray<Body> &mesh, String fileName, const Surface &surf, bool x0z, int &npanels);
	static void SavePreBody(String fileName, const Surface &surf);
	void SaveKH(String fileName) const; 
		
	virtual ~NemohBody() noexcept {}

private:
	String LoadDat0(String fileName, bool &x0z);
	static void SaveDat0(String fileName, const Surface &surf, bool x0z, int &npanels);
};

class SalomeBody : public Body {
public:
	static String LoadDat(UArray<Body> &mesh, String fileName);
		
	virtual ~SalomeBody() noexcept {}

private:
	String LoadDat0(String fileName);
};

class HydrostarBody : public Body {
public:
	static String LoadHst(UArray<Body> &mesh, String fileName, bool &y0z, bool &x0z);
		
	virtual ~HydrostarBody() noexcept {}

private:
	static String LoadHst0(String fileName, Body &mesh, Body &damping, bool &y0z, bool &x0z);
};

class HAMSBody : public Body {
public:
	static String LoadPnl(UArray<Body> &mesh, String filefCaseName, bool &y0z, bool &x0z);
	static void SavePnl(String fileName, const Surface &surf, bool y0z, bool x0z);
	
	virtual ~HAMSBody() noexcept {}
};

class WamitBody : public Body {
public:
	static String LoadDat(UArray<Body> &mesh, String fileName);
	static String LoadGdf(UArray<Body> &mesh, String fileName, bool &y0z, bool &x0z);
	static void SaveGdf(String fileName, const Surface &surf, double g, bool y0z, bool x0z);
	void SaveHST(String fileName, double rho, double g) const; 

	virtual ~WamitBody() noexcept {}
};

class AQWABody : public Body {
public:
	static String LoadDat(UArray<Body> &mesh, String fileName, bool &y0z, bool &x0z);
	static String Load_LIS(UArray<Body> &mesh, String fileName, double g, bool &y0z, bool &x0z);
	static void SaveDat(String fileName, const UArray<Body> &meshes, const UArray<Surface> &surf, double rho, double g, bool y0z, bool x0z);
	static void SaveDat(Stream &ret, const UArray<Body> &meshes, const UArray<Surface> &surf, double rho, double g, bool y0z, bool x0z);
	
	virtual ~AQWABody() noexcept {}
};

class ORCABody : public Body {
public:
#ifdef PLATFORM_WIN32
	static String Load_OWR(UArray<Body> &mesh, String fileName, double g, bool &y0z, bool &x0z);
#endif
	virtual ~ORCABody() noexcept {}
};

class DiodoreBody : public Body {
public:
	static String LoadDat(UArray<Body> &mesh, String fileName);
	static void SaveDat(String fileName, const Surface &surf);
	
	virtual ~DiodoreBody() noexcept {}
};

class Wamit : public Hydro {
public:
	Wamit() {}
	String Load(String file, Function <bool(String, int)> Status);
	void Save(String file, Function <bool(String, int)> Status, bool force_T = false, int qtfHeading = Null);
	void Save_out(String file);
	virtual ~Wamit() noexcept {}
	
	bool LoadGdfBody(String file);
	bool LoadDatBody(String file);
	void SaveGdfBody(String fileName);
	
	static void Save_hst_static(const MatrixXd &C, String fileName, double rho, double g);
	
	bool Load_frc2(String fileName);
	void Save_4(String fileName, bool force_T = false);
	
protected:
	void ProcessFirstColumnPot(UVector<double> &w, UVector<double> &T, int iperin);
	bool ProcessFirstColumn1_3(UVector<double> &w, UVector<double> &T, int iperout);
	
	bool Load_cfg(String fileName, int &iperin, int &iperout);
	bool Load_pot(String fileName);
	bool Load_gdf(String fileName);
	bool Load_mmx(String fileName);
	
	bool Load_out(String fileName);							
	void Load_A(FileInLine &in, MatrixXd &A);
	bool Load_Scattering(String fileName, int iperout);
	bool Load_FK(String fileName, int iperout);
	
	bool Load_Forces(String fileName, Hydro::Forces &force, int iperout);
		
	bool Load_1(String fileName, int iperout);				
	bool Load_3(String fileName, int iperout);
	bool Load_hst(String fileName);
	bool Load_4(String fileName, int iperout);
	bool Load_12(String fileName, bool isSum, Function <bool(String, int)> Status);
	bool Load_789(String fileName, int iperout);
	bool Load_789_0(String fileName, int type, UArray<UArray<UArray<VectorXd>>> &qtf, int iperout);
	
	void Save_1(String fileName, bool force_T = false);
	void Save_3(String fileName, bool force_T = false);
	void Save_hst(String fileName);
	void Save_12(String fileName, bool isSum, Function <bool(String, int)> Status,
				bool force_T = false, bool force_Deg = true, int qtfHeading = Null);
	void Save_789(String fileName, bool force_T, bool force_Deg);
	void Save_FRC(String fileName);
	void Save_POT(String fileName);
		
	void Save_A(FileOut &out, Function <double(int, int)> fun, const MatrixXd &base, String wavePeriod);
	void Save_AB(FileOut &out, int ifr);
	void Save_Forces(FileOut &out, int ifr);
	void Save_RAO(FileOut &out, int ifr);
	void Save_MD(FileOut &out, int ifr);
};

class Hams : public Wamit {
public:
	Hams() {}
	String Load(String file, Function <bool(String, int)> Status);
	virtual ~Hams() noexcept {}
	
	bool Load_Settings(String settingsFile);
	bool Load_HydrostaticBody(String fileName, double rhog);
	
	bool Load_In(String fileName);
	void SaveFolder(String folder, bool bin, int numCases, int numThreads, int solver) const;
	UVector<String> Check() const;
	
	bool LoadHydrostatic(String fileName);

private:
	void SaveFolder0(String folderBase, bool bin, int numCases, bool deleteFolder, int numThreads) const;
	static void OutMatrix(FileOut &out, String header, const MatrixXd &mat);
	static void InMatrix(LineParser &f, MatrixXd &mat);
		
	void Save_Hydrostatic(String folderInput) const;
	void Save_ControlFile(String folderInput, int _nf, double _minf, double _maxf,
							int numThreads) const;
	void Save_Settings(String folderInput, bool thereIsLid) const;
	void Save_Bat(String folder, String batname, String caseFolder, bool bin, String solvName, String meshName) const;
};

class Fast : public Wamit {
public:
	Fast() : WaveNDir(Null), WaveDirRange(Null) {}
	String Load(String file, Function <bool(String, int)> Status);
	void Save(String file, Function <bool(String, int)> Status, int qtfHeading = Null);
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

class Foamm : public Hydro {
public:
	Foamm() {}
	String Load(String file);
	void Get_Each(int ibody, int idf, int jdf, double from, double to, const UVector<double> &freqs, Function <bool(String, int)> Status, Function <void(String)> FOAMMMessage);
	void Get(const UVector<int> &ibs, const UVector<int> &idfs, const UVector<int> &jdfs,
		const UVector<double> &froms, const UVector<double> &tos, const UVector<UVector<double>> &freqs, 
		Function <bool(String, int)> Status, Function <void(String)> FOAMMMessage);
	virtual ~Foamm() noexcept {}
	
protected:
	void Load_mat(String fileName, int ib, int jb, bool loadCoeff);
};

class Nemoh : public Hydro {
public:
	Nemoh() {}
	String Load(String file, Function <bool(String, int)> Status = Null, double rho = Null);
	void Save(String file);
	virtual ~Nemoh() noexcept {}
	
	bool LoadDatBody(String file);
	void SaveDatBody(String file); 
	
	void SaveFolder(String folder, bool bin, int numCases, int numThreads, int solver) const;
	
	void Save_Cal(String folder, int _nf, double _minf, double _maxf, const UVector<int> &nodes, const UVector<int> &panels, int solver) const;
	
	bool Save_KH(String folder) const;
	bool Save_Inertia(String folder) const;
	
	bool Load_Hydrostatics(String folder, String subfolder);
	static bool Load_Hydrostatics_static(String folder, int Nb, UArray<Body> &msh);
	void Save_Hydrostatics(String subfolder) const;
	static void Save_Hydrostatics_static(String folder, int Nb, const UArray<Body> &msh);

	static bool Load_6x6(Eigen::MatrixXd &C, String file);
	static bool Save_6x6(const Eigen::MatrixXd &C, String file);
	
	UVector<String> Check() const;

private:
	bool Load_Cal(String fileName);
	bool Load_Inf(String fileName);
	bool Load_KH(String folder, String subfolder);
	bool Load_Radiation(String fileName);
	bool Load_Excitation(String folder);
	bool Load_Diffraction(String folder);
	bool Load_FroudeKrylov(String folder);
	bool Load_Forces(Hydro::Forces &f, String nfolder, String fileName);
	bool Load_IRF(String fileName);
	bool Load_Inertia(String folder, String subfolder);
	bool Load_LinearDamping(String folder, String subfolder);
	bool Load_QTF(String folder, String subfolder, Function <bool(String, int)> Status);
	bool Load_12(String fileName, bool isSum, Function <bool(String, int)> Status);
	
	static int GetNumArgs(const LineParser &f);

	void Save_Id(String folder) const;
	void Save_Bat(String folder, String batname, String caseFolder, bool bin, 
		String preName, String hydroName, String solvName, String postName) const;
	void Save_Body_cal(String folder, int ib, String meshFile, Body &mesh, int npanels, bool x0z, const Point3D &cg, double rho, double g) const;
	void Save_Body_bat(String folder, String caseFolder, const UVector<String> &meshes, String meshName, bool bin) const;
	void Save_Input(String folder, int solver) const;
	
	void SaveFolder0(String folder, bool bin, int numCases,  
					bool deleteFolder, int solver) const;
};

class Aqwa : public Hydro {
public:
	Aqwa() {}
	String Load(String file, Function <bool(String, int)> Status, double rho = Null);
	void Save(String file, Function <bool(String, int)> Status);
	virtual ~Aqwa() noexcept {}
	
	bool Load(String fileName);
	//UVector<String> Check() const;	
	
private:
	bool Load_AH1();
	bool Load_LIS(double &factorMass, Function <bool(String, int)> Status);
	bool Load_QTF(double factorMass);
	void Save_QTF(String file, Function <bool(String, int)> Status);
};

class Diodore : public Hydro {
public:
	Diodore() {}
	String Load(String file, double rho = Null);
	virtual ~Diodore() noexcept {}	
	
private:
	void Load_HDB();
};

class OrcaWave : public Hydro {
public:
	OrcaWave() {}
	String Load(String file, double rho = Null);
	virtual ~OrcaWave() noexcept {}	
	
private:
	void Load_YML_Res();
#ifdef PLATFORM_WIN32
	void Load_OWR();
#endif
};


class OrcaFactors {
public:
	double mass = 1, len = 1, force = 1;
	
	Matrix<double, 6, 6> A, K, B, M;
	Eigen::Vector<double, 6> F, RAO, MD;
	
	void Update() {
		for (int r = 0; r < 6; ++r)	{
			for (int c = 0; c < 6; ++c) {
				A(r, c) = A_(r, c);
				B(r, c) = B_(r, c);
				K(r, c) = K_(r, c);
				M(r, c) = M_(r, c);
			}
			F(r) = F_(r);
			RAO(r) = RAO_(r);
			MD(r) = MD_(r);
		}
	}
	
private:
	double A_(int r, int c) const {
		if (r < 3 && c < 3)
			return mass;
		else if (r >= 3 && c >= 3)
			return len*len*len;
		else
			return mass*len;
	}
	double K_(int r, int c) const {
		if (r < 3 && c < 3)
			return force/len;
		if (r < 3)
			return force;
		if (c < 3)
			return force;
		return force*len;
	}
	double B_(int r, int c) const {
		if (r < 3 && c < 3)
			return force/len;
		if (r < 3)
			return force;
		if (c < 3)
			return force;
		return force*len;
	}
	double M_(int r, int c) const {
		if (r < 3 && c < 3)
			return mass;
		else if (r >= 3 && c >= 3)
			return mass*len*len;
		else
			return mass*len;
	}
	double F_(int r) const {
		if (r < 3) 
			return force/len;
		return force;
	}
	double RAO_(int r) const {
		if (r < 3)
			return 1;
		return 1/len;
	}
	double MD_(int r) const {
		if (r < 3)
			return force/len/len;
		return force/len;
	}
};

class BemioH5 : public Hydro {
public:
	BemioH5() {}
	String Load(String file, double rho = Null);
	void Save(String file);
	virtual ~BemioH5() noexcept {}	
	
private:
	void Load_H5();
};

String CapyNC_Load(const String &file, UArray<Hydro> &hydros, int &num);

class CapyNC : public Hydro {
public:
	CapyNC() {}
	bool Load(String file, double rho = Null);
	virtual ~CapyNC() noexcept {}	
};

UVector<int> NumSets(int num, int numsets);	


class LineParserWamit: public LineParser {
public:
	LineParserWamit(FileInLine &_in) : LineParser(_in) {}
	
	void LoadWamitJoinedFields(String _line);// Trick for "glued" fields in Wamit
};

String FormatWam(double d);

class BEM {
public:
	BEM();
	
	Hydro &Duplicate(int id);
	Hydro &Average(UVector<int> &ids);
		
	UArray<Hydro> hydros;
	UArray<Body> surfs;
	
	int GetHydroId(int id) {
		for (int i = 0; i < hydros.size(); ++i) {
			if (hydros[i].dt.GetId() == id)
				return i;
		}
		return -1;
	}
	int GetBodyId(int id) {	
		for (int i = 0; i < surfs.size(); ++i) {
			if (surfs[i].dt.GetId() == id)
				return i;
		}
		return -1;
	}
		
	static Function <void(String)> Print, PrintWarning, PrintError;	
	
	UVector<double> headAll;	// Common models data
	UArray<std::complex<double>> headAllMD;
	
	int Nb = 0;				
	
	double depth, rho, g, len;
	
	enum DOFType {DOF123, DOFSurgeSway, DOFxyz};
	static DOFType dofType;
	enum HeadingType {HEAD_180_180, HEAD_0_360};
	static HeadingType headingType;
	static const char *strDOFType[];
	static const char *strHeadingType[];
	
	int calcAinf = false, calcAinf_w = false;
	int legend_w_units = true, legend_w_solver = true;
	double maxTimeA = Null;
	int numValsA = Null;
	int onlyDiagonal;
	
	String nemohPath, nemoh115Path, nemoh3Path, nemohPathGREN;
	bool experimental;
	String foammPath;
	String hamsPath, hamsBodyPath;
	int volWarning, volError;
	double roundVal, roundEps;
	String csvSeparator;
	String pythonEnv;
	
	int LoadBEM(String file, Function <bool(String, int pos)> Status = Null, bool checkDuplicated = false);
	Hydro &Join(UVector<int> &ids, Function <bool(String, int)> Status = Null);
	void SymmetrizeForces(int id, bool xAxis);
	void Symmetrize(int id);
	void A0(int id);
	void Kirf(int id, double maxT);
	void Ainf(int id);
	void Ainf_w(int id);
	void RAO(int id);
	void BH(int id, int &num);
	void OgilvieCompliance(int id, bool zremoval, bool thinremoval, bool decayingTail, UVector<int> &vidof, UVector<int> &vjdof);
	void TranslationTo(int id, const MatrixXd &to);
	void WaveTo(int id, double xto, double yto);
	String SpreadNegative(int id, Function <bool(String, int)> Status);
	void DeleteHeadingsFrequencies(int id, const UVector<int> &idFreq, const UVector<int> &idFreqQTF, 
										   const UVector<int> &idHead, const UVector<int> &idHeadMD, const UVector<int> &idHeadQTF);
	void ResetForces(int id, Hydro::FORCE force, bool forceMD, Hydro::FORCE forceQtf);										
	void MultiplyDOF(int id, double factor, const UVector<int> &idDOF, bool a, bool b, bool diag, bool f, bool md, bool qtf);
	void SwapDOF(int id, int ib1, int idof1, int ib2, int idof2);
	
	void RemoveHydro(int id);
		
	void FillFrequencyGapsABForces(int id, bool zero, int maxFreq);
	void FillFrequencyGapsQTF(int id, bool zero, int maxFreq);
	
	void FillFrequencyGapsABForcesZero(int id);
	void FillFrequencyGapsQTFZero(int id);
	
	int LoadBody(String file, Function <bool(String, int pos)> Status, bool cleanPanels, bool checkDuplicated);
	void SaveBody(String fileName, const UVector<int> &ids, Body::MESH_FMT type, Body::MESH_TYPE meshType, bool symX, bool symY);
	void HealingBody(int id, bool basic, Function <bool(String, int pos)> Status);
	void OrientSurface(int id, Function <bool(String, int)> Status);
	void ImageBody(int id, int axis);
	void UnderwaterBody(int id, Function <bool(String, int pos)> Status);
	void RemoveBody(int id);
	void JoinBody(int idDest, int idOrig);
	UVector<int> SplitBody(int id, Function <bool(String, int pos)> Status);
	
	void CopyQTF_MD(int id);
		
	void AddFlatRectangle(double x, double y, double z, double size, double panWidthX, double panWidthY);
	void AddRevolution(double x, double y, double z, double size, UVector<Pointf> &vals);
	void AddPolygonalPanel(double x, double y, double z, double size, UVector<Pointf> &vals);
	void AddWaterSurface(int id, char c);
	
	String LoadSerializeJson();
	bool StoreSerializeJson();
	bool ClearTempFiles();
	static String GetTempFilesFolder() {return AFX(GetAppDataFolder(), "BEMRosetta", "Temp");}
	
	void UpdateHeadAll();
	void UpdateHeadAllMD();
	
	//const String bemFilesExt = ".1 .2 .3 .hst .4 .12s .12d .frc .pot .out .in .cal .tec .inf .ah1 .lis .qtf .mat .dat .bem .fst .yml";
	const String bstFilesExt = ".in .out .fst .1 .2 .3 .hst .4 .12s .12d .frc .pot .mmx .cal .tec .inf .ah1 .lis .qtf .hdb .mat .dat .bem .yml .h5 .nc"	// Priority
#ifdef PLATFORM_WIN32
		" .owr"
#endif
	;
	const UVector<String> bemExtSets = {".1.2.3.hst.4.9.12s.12d.frc.pot.mmx", ".lis.qtf.dat"};	// Any of these files opens all, and it is avoided to load them again
	String bemFilesAst;
	
	int GetBEMExtSet(String file) {
		String ext = ToLower(GetFileExt(file));
		for (int i = 0; i < bemExtSets.size(); ++i)
			if (bemExtSets[i].Find(ext) >= 0)
				return i;
		return -1;
	}
	
	void Jsonize(JsonIO &json) {
		String nemohPathPreprocessor, nemohPathNew;
		int idofType, iheadingType;
		if (json.IsLoading()) {
			idofType = 0;
			iheadingType = 0;
		} else {
			idofType = dofType;
			iheadingType = headingType;
		}
		json
			("depth", depth)
			("rho", rho)
			("g", g)
			("length", len)
			("calcAwinf", calcAinf)
			("calcAwinfw", calcAinf_w)
			("maxTimeA", maxTimeA)
			("numValsA", numValsA)
			("onlyDiagonal", onlyDiagonal)
			("nemohPathPreprocessor", nemohPathPreprocessor)
			("nemohPathGREN", nemohPathGREN)
			("nemohPathNew", nemohPathNew)
			("nemohPath", nemohPath)
			("nemoh115Path", nemoh115Path)
			("nemoh3Path", nemoh3Path)
			("foammPath", foammPath)
			("hamsPath", hamsPath)
			("hamsBodyPath", hamsBodyPath)
			("volWarning", volWarning)
			("volError", volError)
			("roundVal", roundVal)
			("roundEps", roundEps)
			("dofType", idofType)
			("headingType", iheadingType)
			("csvSeparator", csvSeparator)
			("legend_w_units", legend_w_units)
			("legend_w_solver", legend_w_solver)
			("pythonEnv", pythonEnv)
		;
		if (json.IsLoading()) {
			dofType = BEM::DOFType(idofType);
			headingType = BEM::HeadingType(iheadingType);
			if (IsEmpty(nemohPath))
				nemohPath = GetFileFolder(nemohPathPreprocessor);
			if (IsEmpty(nemoh115Path))
				nemoh115Path = GetFileFolder(nemohPathNew);
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
				if ((int)strlen(str[i]) > mx)
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
	static const char *strDOFnum_sub[];
	
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


class Wind {
public:
	String Load(String fileName, String ext = "");
	String Save(String fileName, String ext = "");
	
	const char *GetWindTypeStr() const;
	const String &GetDescription() const {return description;}
	void SetHubHeight(double h)			 {zHub  = float(h);}
	void SetGridHeight(double h)		 {zGrid = float(h);	SetYZ();}
	
	void GetPos(double z, double y, int &idz, int &idy);
	VectorXd Get(int idz, int idy);
	VectorXd GetTime();
	int GetTimeId(double time);
	
	void Report(Grid &grid) const;	

protected:
	String fileName; 
	String description;
	Tensor<double, 4> velocity;		// 4-D array: time, velocity component (1=U, 2=V, 3=W), iy, iz 
	Tensor<double, 3> twrVelocity;	// 3-D array: time, velocity component, iz
	VectorXd yPos;						// horizontal locations y(iy)
	VectorXd zPos;						// vertical locations z(iz)
	//VectorXd zTwr;					// vertical locations of tower points zTwr(iz)
	int nz = Null, ny = Null;		// number of points in the vertical and horizontal directions of the grid
	float dz = Null, dy = Null, dt = Null;		// distance between two points in the vertical [m], horizontal [m], and time [s] dimensions
	float zHub = Null;				// hub height [m]
	float zGrid = Null;				// vertical location of bottom of grid [m above ground level]
	float mffws = Null;				// mean hub-height wind speed
	int nt = Null;					// the number of time steps, INT(4)
	int ntwr = Null;				// the number of tower points, INT(4)
    int fc = -1;						// 1 = 1-component von Karman
										// 2 = 1-component Kaimal
										// 3 = 3-component von Karman
										// 4 = improved von Karman
										// 5 = IEC-2 Kaimal
										// 6 = (not supported)
										// 7 = General Kaimal
										// 8 = Mann model
    
    int turbSimFormat = 7;				// TurbSim format identifier (should = 7 or 8 if periodic), INT(2)
    
    bool LHR = true;		// Default value for Bladed
    
    void SetYZ();
};

class BTSWind : public Wind {
public:
	String LoadBTS(String file);
	String SaveBTS(String file, int fmtSz = -1) const;
	
private:
	void LoadBTSHeader(FileInBinary &file, VectorXf &Vslope, VectorXf &Voffset);
	
	template <class T>
	String LoadBTSBody(FileInBinary &file, const VectorXf &Vslope, const VectorXf &Voffset) {
		int nPts    = ny*nz;
    	int nv      = 3*nPts;               // the size of one time step
    	int nvTwr   = 3*ntwr;
    
        velocity    = Tensor<double, 4>(nt,3,ny,nz);
    	twrVelocity = Tensor<double, 3>(nt,3,ntwr);
  
		Buffer<T> data(nv);
		Buffer<T> datat(nvTwr);
		int fmtSz = sizeof(T);
		for (int it = 0; it < nt; ++it) {
			// get the grid points
			file.Read(data, fmtSz*nv); // read the velocity components for one time step
			
			int ip = 0;
	    	for (int iz = 0; iz < nz; ++iz) {
	    		for (int iy = 0; iy < ny; ++iy) {
	    			for (int k = 0; k < 3; ++k) {
	                    double d = (double(data[ip++]) - Voffset(k))/Vslope(k);
	            		if (d > 200 || d < -200)
	                        return t_("Wrong format");
	                    velocity(it,k,iy,iz) = d;
	                }
	    		}
	    	}
			// get the tower points
			if (ntwr > 0) {
				file.Read(datat, fmtSz*nvTwr);		// read the velocity components for the tower

	            for (int k = 0; k < 3; ++k) {      // scale the data
	                for (int itw = 0; itw < ntwr; ++itw) {
						double d = (double(datat[itw*3 + k]) - Voffset(k))/Vslope(k);
						if (d > 500 || d < -500)
							return t_("Wrong format");
	                	twrVelocity(it,k,itw) = d; 
	                }
	            }
			}
		}
		return "";
	}
	
	void SaveBTSHeader(FileOutBinary &file, VectorXf &Vslope, VectorXf &Voffset, int fmtSz) const;
	
	template <class T>
	void SaveBTSBody(FileOutBinary &file, const VectorXf &Vslope, const VectorXf &Voffset) const {
		int nPts    = ny*nz;
	    int nv      = 3*nPts;               // the size of one time step
	    int nvTwr   = 3*ntwr;
	
		Buffer<T> data(nv);
		Buffer<T> datat(nvTwr);
		int fmtSz = sizeof(T);
		
		for (int it = 0; it < nt; ++it) {
			int ip = 0;
	        for (int iz = 0; iz < nz; ++iz) {
	            for (int iy = 0; iy < ny; ++iy) {
	                for (int k = 0; k < 3; ++k) {
	                    data[ip++] = BetweenVal(T(velocity(it,k,iy,iz)*Vslope(k) + Voffset(k)),	
	                    						std::numeric_limits<T>::lowest(),
	                    						std::numeric_limits<T>::max());
	                }
	            }
	        }			
			file.Write(data.begin(), fmtSz*nv);

			if (ntwr > 0) {
	            for (int k = 0; k < 3; ++k)      // scale the data
	                for (int itw = 0; itw < ntwr; ++itw) 
	                    datat[itw*3 + k] = BetweenVal(T(fround(twrVelocity(it,k,itw)*Vslope(k) + Voffset(k))), 
	                    						 	  std::numeric_limits<T>::lowest(), 
	                    						 	  std::numeric_limits<T>::max());

	            file.Write(datat.begin(), fmtSz*nvTwr);		// read the velocity components for the tower
			}
		}		
	}
};

class WNDWind : public Wind {
public:
	String LoadWND(String file, double _zHub = Null);	
	
private:
	bool LoadSum(String fileName, const UVector<String> &str, UVector<float> &SummVars, float &ZGoffset);
};


class ArrayWind : public Upp::Array<Wind> {
public:
	void Report(Grid &grid);
};


// Compiler options

// Normal
// -Wall -Wextra 

// Extra
// -Wno-unused-parameter -Wno-logical-op-parentheses -Wno-deprecated-copy-with-user-provided-copy -Wno-overloaded-virtual -Wno-missing-braces 

// Full
// -Wshadow	

#endif
