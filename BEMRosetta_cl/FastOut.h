// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#ifndef _BEMRosetta_BEMRosetta_cl_FastOut_h_
#define _BEMRosetta_BEMRosetta_cl_FastOut_h_


//bool FindHydrodynCB(String path, double &ptfmCOBxt, double &ptfmCOByt);
	
class FastOut {
public:
	FastOut();
    virtual ~FastOut() {};
    
	static UVector<String> GetFilesToLoad(String path);
	static String GetFileToLoad(String fileName);
	
	String Load(String fileName, Function <bool(String, int)> Status = Null);
	bool Save(String fileName, Function <bool(String, int)> Status, String type = "", String sep = "", const UVector<int> &ids = UVector<int>());
	
	void AppendLine(int numLine, FastOut &fst);
	
	void Clear();
	bool IsEmpty();
	int GetParameter_throw(String param) const;
	int GetParameterX(String param) const;
	UVector<int> FindParameterMatch(String param) const;
	UVector<String> FindParameterMatchStr(String param) const;
	
	const String &GetParameter(int idparam) const	{return parameters[idparam];}
	const String &GetUnit(int idparam) const		{return units[idparam];}
	int GetParameterCount() const					{return parameters.size();}
	 
	SortedIndex<String> GetParameterList(String filter = ""); 
	SortedIndex<String> GetUnitList(String filter = "");
	SortedVectorMap<String, String> GetList(String filterParam = "", String filterUnits = "");
		
	double GetVal(double time, int idparam) const;
	inline double GetVal(int idtime, int idparam) const	  		{return dataOut[idparam][idtime];}
	inline const UVector<double> &GetUVector(int idparam) const {return dataOut[idparam];}
	inline const UVector<double> &GetUVector(String param) const {
		static UVector<double> bad;
		UVector<int> ids = FindParameterMatch(param);
		if (ids.IsEmpty())
			return bad;
		else
			return GetUVector(ids[0]);
	}
	VectorXd GetVector(int idparam) const {
		const UVector<double> &d = GetUVector(idparam);
		VectorXd ret = Map<const VectorXd>(d, d.size());
		return ret;
	}
	VectorXd GetVector(String param) const {
		const UVector<double> &d = GetUVector(param);
		VectorXd ret = Map<const VectorXd>(d, d.size());
		return ret;
	}
	
	void SetVal(int idparam, double val);
	void SetNextTime(double time);
		
	int GetIdTime(double time) const;
	double GetTimeStart() const;
	double GetTimeEnd()	 const;
	int GetNumData() const;
	bool IsEmpty() const		{return dataOut.IsEmpty();}	
	
	String GetFileName() const	{return fileName;}
	
	int ColFairlead(int i) const{return GetParameter_throw(Format("T[%d]", i-1));}
		
	int GetNumFairlead() const {
		for (int i = 0; i < 100; ++i) {
			if (ColFairlead(i) < 0) 
				return i;
		}
		return Null;
	}
	
	void Serialize(Stream& s) {
        s % dataOut % parameters % units;
        if (s.IsLoading()) {
            parametersd.Clear();
			for (int i = 0; i < parameters.size(); ++i)
				parametersd << ToLower(parameters[i]);
			unitsd.Clear();
			for (int i = 0; i < units.size(); ++i)
				unitsd << ToLower(units[i]);
        }
    }
	
	struct CalcParam {
		CalcParam *Init0(FastOut *_dataFast) {dF = _dataFast; return this;}
		virtual void Init() = 0;
		virtual double Calc(int idt) = 0;
		bool IsEnabled()	{return enabled;}
		String name, units;
		int id = -1;
		
	protected:
		FastOut *dF = nullptr;
		bool enabled = true;
	};
	
	int AddParam(String name, String unit, String description = "") {
		parameters << name;
		units << unit;
		parametersd << ToLower(name);
		unitsd << ToLower(unit);
		dataOut.Add();
		return dataOut.size() - 1;
	}
	
	int GetParamCount() {return parameters.size();}
 
	UVector<String> parameters, units, descriptions;
	UVector<String> parametersd, unitsd;
	UVector<UVector <double> > dataOut;
	UVector<CalcParam *> calcParams;
	
	struct PointParam {
		Point3D pos;
		String name;
		int id = -1;
		
		PointParam(String _name, Point3D _pos, FastOut *_dataFast) {
			name = _name;
			pos = _pos;
			dF = _dataFast;
		}
		
		Point3D Calc(int idtime) {
			Point3D npos;
			TransRot(dF->aff[idtime], pos, npos);	 
			return npos;
		}
		
	protected:
		FastOut *dF = nullptr;	
	}; 
	
	UArray<PointParam> pointParams;
	
	int FindParam(String name) {
		for (int i = 0; i < pointParams.size(); ++i) {
			if (pointParams[i].name == name)
				return i;
		}
		return -1;
	}
	
	UArray<Affine3d> aff;
	double Hx = Null, Hz = Null;			// Position of the hub projection to the blade tip rotation plane. with NacYaw = 0
	
	int idsurge = Null, idsway = Null, idheave = Null, idroll = Null, idpitch = Null, idyaw = Null, idaz = Null, idnacyaw = Null;
	
	double TipRad = Null, OverHang = Null, ShftTilt = Null, Precone = Null, Twr2Shft = Null, TowerHt = Null, baseClearance = Null;
	
	double Hs = Null, Tp = Null, heading = Null;
    
    class id3d : public Moveable<id3d> {
    public:
        int x = Null, y = Null, z = Null;
    };
	UVector<UVector<id3d>> mooringPointIds;
	UVector<int> idPos;
	
private:
	String LoadOut(String fileName, Function <bool(String, int)> Status);
	String LoadOutb(String fileName, Function <bool(String, int)> Status);
	bool SaveOut(String fileName, Function <bool(String, int)> Status, const UVector<int> &ids);
	String LoadCsv(String fileName, Function <bool(String, int)> Status);
	bool SaveCsv(String fileName, Function <bool(String, int)> Status, String sep, const UVector<int> &ids);
	String LoadDb(String fileName, Function <bool(String, int)> Status);
	String Load_LIS(String fileName);
	void AfterLoad();

	String fileName;
	
	int idActualTime = -1;
	
	VectorMap<String, String> syn;		// Synonyms
	
	struct TiltParam : CalcParam {
		TiltParam() {
			name = "PtfmTilt";
			units = "deg";
		}
		virtual void Init() {
			enabled = !(dF->idroll < 0 || dF->idpitch < 0);
		}
		virtual double Calc(int idtime) {
			double roll = dF->GetVal(idtime, dF->idroll);
			double pitch = dF->GetVal(idtime, dF->idpitch);

			if (IsNull(roll) || IsNull(pitch))
				return Null;
							
			return sqrt(roll*roll + pitch*pitch);
		}	
	} ptfmtilt; 

	struct ShiftParam : CalcParam {
		ShiftParam() {
			name = "PtfmShift";
			units = "m";
		}
		virtual void Init() {
			enabled = !(dF->idsurge < 0 || dF->idsway < 0);
		}
		virtual double Calc(int idtime) {
			double surge = dF->GetVal(idtime, dF->idsurge);
			double sway = dF->GetVal(idtime, dF->idsway);
			
			if (IsNull(surge) || IsNull(sway))
				return Null;			
			
			return sqrt(surge*surge + sway*sway);
		}	
	} ptfmshift; 

	bool CalcTipPos(int idBlade, int idtime, double tipdx, double tipdy, double &Tx, double &Ty, double &Tz);
	
	struct BladeTip1xParam : CalcParam {
		BladeTip1xParam() {
			name = "BldTip1x";
			units = "m";
		}
		virtual void Init();
		
		virtual double Calc(int idtime) {
			double tipdx = 0, tipdy = 0;
			if (!IsNull(idTipdx) && !IsNull(idTipdy)) {
				tipdx = Nvl(dF->GetVal(idtime, idTipdx), 0.);
				tipdy = Nvl(dF->GetVal(idtime, idTipdy), 0.);
			}
			if ((isnull = !dF->CalcTipPos(0, idtime, tipdx, tipdy, Tx, Ty, Tz)))
				return Null;
			
			return Tx;
		}	
		double Tx, Ty, Tz;
		int idTipdx = Null, idTipdy = Null;
		bool isnull;
	} bladeTip1xParam; 

	struct BladeTip1yParam : CalcParam {
		BladeTip1yParam() {
			name = "BldTip1y";
			units = "m";
		}
		virtual void Init() {enabled = dF->bladeTip1xParam.IsEnabled();}
		
		virtual double Calc(int idtime) {
			if (dF->bladeTip1xParam.isnull)
				return Null;
			
			return dF->bladeTip1xParam.Ty;
		}	
	} bladeTip1yParam; 


	struct BladeTip1zParam : CalcParam {
		BladeTip1zParam() {
			name = "BldTip1z";
			units = "m";
		}
		virtual void Init() {enabled = dF->bladeTip1xParam.IsEnabled();}
		
		virtual double Calc(int idtime) {
			if (dF->bladeTip1xParam.isnull)
				return Null;
			
			return dF->bladeTip1xParam.Tz;
		}	
	} bladeTip1zParam; 

	struct BladeTip2xParam : CalcParam {
		BladeTip2xParam() {
			name = "BldTip2x";
			units = "m";
		}
		virtual void Init();
		
		virtual double Calc(int idtime) {
			double tipdx = 0, tipdy = 0;
			if (!IsNull(idTipdx) && !IsNull(idTipdy)) {
				tipdx = Nvl(dF->GetVal(idtime, idTipdx), 0.);
				tipdy = Nvl(dF->GetVal(idtime, idTipdy), 0.);
			}
			if ((isnull = !dF->CalcTipPos(1, idtime, tipdx, tipdy, Tx, Ty, Tz)))
				return Null;
			
			return Tx;
		}	
		double Tx, Ty, Tz;
		int idTipdx = Null, idTipdy = Null;
		bool isnull;
	} bladeTip2xParam; 

	struct BladeTip2yParam : CalcParam {
		BladeTip2yParam() {
			name = "BldTip2y";
			units = "m";
		}
		virtual void Init() {enabled = dF->bladeTip2xParam.IsEnabled();}
		
		virtual double Calc(int idtime) {
			if (dF->bladeTip2xParam.isnull)
				return Null;
			
			return dF->bladeTip2xParam.Ty;
		}	
	} bladeTip2yParam; 


	struct BladeTip2zParam : CalcParam {
		BladeTip2zParam() {
			name = "BldTip2z";
			units = "m";
		}
		virtual void Init() {enabled = dF->bladeTip2xParam.IsEnabled();}
		
		virtual double Calc(int idtime) {
			if (dF->bladeTip2xParam.isnull)
				return Null;
			
			return dF->bladeTip2xParam.Tz;
		}	
	} bladeTip2zParam; 
	
	struct BladeTip3xParam : CalcParam {
		BladeTip3xParam() {
			name = "BldTip3x";
			units = "m";
		}
		virtual void Init();
		
		virtual double Calc(int idtime) {
			double tipdx = 0, tipdy = 0;
			if (!IsNull(idTipdx) && !IsNull(idTipdy)) {
				tipdx = Nvl(dF->GetVal(idtime, idTipdx), 0.);
				tipdy = Nvl(dF->GetVal(idtime, idTipdy), 0.);
			}
			if ((isnull = !dF->CalcTipPos(2, idtime, tipdx, tipdy, Tx, Ty, Tz)))
				return Null;
			
			return Tx;
		}	
		double Tx, Ty, Tz;
		int idTipdx = Null, idTipdy = Null;
		bool isnull;
	} bladeTip3xParam; 

	struct BladeTip3yParam : CalcParam {
		BladeTip3yParam() {
			name = "BldTip3y";
			units = "m";
		}
		virtual void Init() {enabled = dF->bladeTip1xParam.IsEnabled();}
		
		virtual double Calc(int idtime) {
			if (dF->bladeTip1xParam.isnull)
				return Null;
			
			return dF->bladeTip3xParam.Ty;
		}	
	} bladeTip3yParam; 


	struct BladeTip3zParam : CalcParam {
		BladeTip3zParam() {
			name = "BldTip3z";
			units = "m";
		}
		virtual void Init() {enabled = dF->bladeTip1xParam.IsEnabled();}
		
		virtual double Calc(int idtime) {
			if (dF->bladeTip3xParam.isnull)
				return Null;
			
			return dF->bladeTip3xParam.Tz;
		}	
	} bladeTip3zParam; 	
				
	struct TwrBsShearParam : CalcParam {
		TwrBsShearParam() {
			name = "TwrBsShear";
			units = "kN";
		}
		virtual void Init() {
			idx = dF->GetParameterX("TwrBsFxt");
			idy = dF->GetParameterX("TwrBsFyt");
			enabled = !(idx < 0 || idy < 0);
		}
		virtual double Calc(int idtime) {
			double fx = dF->GetVal(idtime, idx);
			double fy = dF->GetVal(idtime, idy);
			
			if (IsNull(fx) || IsNull(fy))
				return Null;			
			
			return sqrt(fx*fx + fy*fy);
		}	
		int idx = Null, idy = Null;	
	} twrBsShear;
	
	struct TwrBsBendParam : CalcParam {
		TwrBsBendParam() {
			name = "TwrBsBend";
			units = "kNm";
		}
		virtual void Init() {
			idx = dF->GetParameterX("TwrBsMxt");
			idy = dF->GetParameterX("TwrBsMyt");
			enabled = !(idx < 0 || idy < 0);
		}
		virtual double Calc(int idtime) {
			double mx = dF->GetVal(idtime, idx);
			double my = dF->GetVal(idtime, idy);

			if (IsNull(mx) || IsNull(my))
				return Null;
							
			return sqrt(mx*mx + my*my);
		}	
		int idx = Null, idy = Null;	
	} twrBsBend;	

	struct YawBrShearParam : CalcParam {
		YawBrShearParam() {
			name = "YawBrShear";
			units = "kN";
		}
		virtual void Init() {
			idx = dF->GetParameterX("YawBrFxp");
			idy = dF->GetParameterX("YawBrFyp");
			enabled = !(idx < 0 || idy < 0);
		}
		virtual double Calc(int idtime) {
			double fx = dF->GetVal(idtime, idx);
			double fy = dF->GetVal(idtime, idy);
			
			if (IsNull(fx) || IsNull(fy))
				return Null;
						
			return sqrt(fx*fx + fy*fy);
		}	
		int idx = Null, idy = Null;	
	} yawBrShear;

	struct YawBrBendParam : CalcParam {
		YawBrBendParam() {
			name = "YawBrBend";
			units = "kNm";
		}
		virtual void Init() {
			idx = dF->GetParameterX("YawBrMxp");
			idy = dF->GetParameterX("YawBrMyp");
			enabled = !(idx < 0 || idy < 0);
		}
		virtual double Calc(int idtime) {
			double mx = dF->GetVal(idtime, idx);
			double my = dF->GetVal(idtime, idy);
			
			if (IsNull(mx) || IsNull(my))
				return Null;
						
			return sqrt(mx*mx + my*my);
		}	
		int idx = Null, idy = Null;	
	} yawBrBend;

	struct RootShear1Param : CalcParam {
		RootShear1Param() {
			name = "RootShear1";
			units = "kN";
		}
		virtual void Init() {
			idx = dF->GetParameterX("RootFxc1");
			idy = dF->GetParameterX("RootFyc1");
			enabled = !(idx < 0 || idy < 0);
		}
		virtual double Calc(int idtime) {
			double fx = dF->GetVal(idtime, idx);
			double fy = dF->GetVal(idtime, idy);
			
			return sqrt(fx*fx + fy*fy);
		}	
		int idx = Null, idy = Null;	
	} rootShear1;

	struct RootBend1Param : CalcParam {
		RootBend1Param() {
			name = "RootBend1";
			units = "kNm";
		}
		virtual void Init() {
			idx = dF->GetParameterX("RootMxc1");
			idy = dF->GetParameterX("RootMyc1");
			enabled = !(idx < 0 || idy < 0);
		}
		virtual double Calc(int idtime) {
			double mx = dF->GetVal(idtime, idx);
			double my = dF->GetVal(idtime, idy);
			
			return sqrt(mx*mx + my*my);
		}	
		int idx = Null, idy = Null;	
	} rootBend1;

	struct RootShear2Param : CalcParam {
		RootShear2Param() {
			name = "RootShear2";
			units = "kN";
		}
		virtual void Init() {
			idx = dF->GetParameterX("RootFxc2");
			idy = dF->GetParameterX("RootFyc2");
			enabled = !(idx < 0 || idy < 0);
		}
		virtual double Calc(int idtime) {
			double fx = dF->GetVal(idtime, idx);
			double fy = dF->GetVal(idtime, idy);
			
			return sqrt(fx*fx + fy*fy);
		}	
		int idx = Null, idy = Null;	
	} rootShear2;

	struct RootBend2Param : CalcParam {
		RootBend2Param() {
			name = "RootBend2";
			units = "kNm";
		}
		virtual void Init() {
			idx = dF->GetParameterX("RootMxc2");
			idy = dF->GetParameterX("RootMyc2");
			enabled = !(idx < 0 || idy < 0);
		}
		virtual double Calc(int idtime) {
			double mx = dF->GetVal(idtime, idx);
			double my = dF->GetVal(idtime, idy);
			
			return sqrt(mx*mx + my*my);
		}	
		int idx = Null, idy = Null;	
	} rootBend2;
	
	struct RootShear3Param : CalcParam {
		RootShear3Param() {
			name = "RootShear3";
			units = "kN";
		}
		virtual void Init() {
			idx = dF->GetParameterX("RootFxc3");
			idy = dF->GetParameterX("RootFyc3");
			enabled = !(idx < 0 || idy < 0);
		}
		virtual double Calc(int idtime) {
			double fx = dF->GetVal(idtime, idx);
			double fy = dF->GetVal(idtime, idy);
			
			return sqrt(fx*fx + fy*fy);
		}	
		int idx = Null, idy = Null;	
	} rootShear3;

	struct RootBend3Param : CalcParam {
		RootBend3Param() {
			name = "RootBend3";
			units = "kNm";
		}
		virtual void Init() {
			idx = dF->GetParameterX("RootMxc3");
			idy = dF->GetParameterX("RootMyc3");
			enabled = !(idx < 0 || idy < 0);
		}
		virtual double Calc(int idtime) {
			double mx = dF->GetVal(idtime, idx);
			double my = dF->GetVal(idtime, idy);
			
			return sqrt(mx*mx + my*my);
		}	
		int idx = Null, idy = Null;	
	} rootBend3;
	
	struct NcIMUTAParam : CalcParam {
		NcIMUTAParam() {
			name = "NcIMUTA";
			units = "m/s^2";
		}
		virtual void Init() {
			idx = dF->GetParameterX("NcIMUTAxs");
			idy = dF->GetParameterX("NcIMUTAys");
			idz = dF->GetParameterX("NcIMUTAzs");
			enabled = !(idx < 0 || idy < 0 || idz < 0);
		}
		virtual double Calc(int idtime) {
			double ax = dF->GetVal(idtime, idx);
			double ay = dF->GetVal(idtime, idy);
			double az = dF->GetVal(idtime, idz);
			
			return sqrt(ax*ax + ay*ay + az*az);
		}	
		int idx = Null, idy = Null, idz = Null;	
	} ncIMUTA;
	
	struct Fairten_tParam : CalcParam {
		Fairten_tParam() {
			units = "t";
		}
		virtual ~Fairten_tParam() {}
		
		void Init00(int _id) {
			name = Format("FAIRTEN%d_t", _id);
			id = _id;
		}
		virtual void Init() {
			idFair = dF->GetParameterX(Format("FAIRTEN%d", id));
			enabled = !(idFair < 0);
		}
		virtual double Calc(int idtime) {
			double fairTen = dF->GetVal(idtime, idFair);

			if (IsNull(fairTen))
				return Null;
							
			return fairTen/1000/9.8;
		}	
		int id, idFair;
	};
	UArray<Fairten_tParam> fairTens;	

	struct AnchTenHor : CalcParam {
		AnchTenHor() {
			units = "N";
		}
		virtual ~AnchTenHor() {}
		
		void Init00(int _id) {
			name = Format("ANCHTEN%d_H", _id);
			id = _id;
		}
		virtual void Init() {
			idx = dF->GetParameterX(Format("L%d", id) + "NP0FX");
			idy = dF->GetParameterX(Format("L%d", id) + "NP0FY");
			enabled = idx >= 0  && idy >= 0;
		}
		virtual double Calc(int idtime) {
			double fx = dF->GetVal(idtime, idx);
			double fy = dF->GetVal(idtime, idy);

			if (IsNull(fx) || IsNull(fy))
				return Null;
							
			return sqrt(fx*fx + fy*fy);
		}	
		int id, idx, idy;
	};
	UArray<AnchTenHor> anchTens;

	struct AnchTenV : CalcParam {
		AnchTenV() {
			units = "N";
		}
		virtual ~AnchTenV() {}
		
		void Init00(int _id) {
			name = Format("ANCHTEN%d_V", _id);
			id = _id;
		}
		virtual void Init() {
			idz = dF->GetParameterX(Format("L%d", id) + "NP0FZ");
			enabled = idz >= 0;
		}
		virtual double Calc(int idtime) {
			return dF->GetVal(idtime, idz);
		}	
		int id, idz;
	};
	UArray<AnchTenV> anchTensV;
			
	int id = -1;
	static int staticid;
};


class ParameterMetric : public DeepCopyOption<ParameterMetric> {
public:
	ParameterMetric() {}
	ParameterMetric(const ParameterMetric &d, int) : name(d.name), decimals(d.decimals) {
		metrics = clone(d.metrics);
	}
	String name;
	int decimals;
	UVector<String> metrics;
	
	void Jsonize(JsonIO &json) {
		json
			("name", name)
			("decimals", decimals)
			("metrics", metrics)
		;
	}	
};

struct ParameterMetrics {
	UArray<ParameterMetric> params;	
	
	void Jsonize(JsonIO &json) {
		json
			("params", params)
		;
	}
};

void Calc(const UArray<FastOut> &dataFast, const ParameterMetrics &params, ParameterMetrics &realparams, 
		double start, double end, UVector<UVector<Value>> &table, Function <bool(String, int)> Status = Null);


class FASTCase {
public:
	~FASTCase() {
		//DeleteFolderDeepX(folderCase);
	}
	void Load(String file) {
		String path = GetFileFolder(file);
		
		fast.fileName = file;
		elastodyn.fileName = AFX(path, fast.GetString("EDFile"));
		if (GetFileTitle(elastodyn.fileName).Find("unused") >= 0)
			elastodyn.fileName = "";
		hydrodyn.fileName = AFX(path, fast.GetString("HydroFile"));
		if (GetFileTitle(hydrodyn.fileName).Find("unused") >= 0)
			hydrodyn.fileName = "";
		inflowfile.fileName = AFX(path, fast.GetString("InflowFile"));
		if (GetFileTitle(inflowfile.fileName).Find("unused") >= 0)
			inflowfile.fileName = "";
		subdyn.fileName = AFX(path, fast.GetString("SubFile"));
		if (GetFileTitle(subdyn.fileName).Find("unused") >= 0)
			subdyn.fileName = "";
		
		try {
			dlldat.fileName = GetAbsolutePath(path, hydrodyn.GetString("NLFK_DLL_input"));
		} catch(...) {
		}
		
		elastodyn.IsAvailable();
		int pos;
		String nums = GetFASTVarPos(elastodyn.fileText, "NumTrack", "", pos);
		if (pos >= 0) {
			int num = ScanInt(nums);
			pos = elastodyn.fileText.FindAfter("\n", pos);
			pos = elastodyn.fileText.FindAfter("\n", pos);
			pos = elastodyn.fileText.FindAfter("\n", pos);
			for (int i = 0; i < num; ++i) {
				int npos = elastodyn.fileText.FindAfter("\n", pos);
				String str = elastodyn.fileText.Mid(pos, npos-pos);
				UVector<String> dat = Split(str, " ");
				if (dat.size() < 4)
					continue;
				String name = Trim(dat[0]);
				Point3D point(ScanDouble(dat[1]), ScanDouble(dat[2]), ScanDouble(dat[3]));
				if (IsNull(point))
					continue;
				pointNames << name;
				points << point;
				pos = npos;
			}
		}
	}
	void Save() {
		fast.Save();
		elastodyn.Save();
		hydrodyn.Save();
		dlldat.Save();
	}
	void Setup(String seed, String folder = "");
	void CreateFolderCase(String folder = "");
	
	bool LoadOut() {
		return out.Load(fstFile, Null) > 0;
	}
	virtual bool Postprocess() 		{return false;};
	const String &GetFolderCase()	{return folderCase;}
	
//private:
	class File {
	public:
		String fileName;
		String fileText;
		bool isChanged = false;
		
		void Save() const {
			if (!isChanged || fileName.IsEmpty())
				return;
			if (!SaveFile(fileName, fileText))
				throw Exc(Format(t_("Impossible to save file '%s'"), fileName));
		}
		
		bool IsAvailable() {
			if (fileText.IsEmpty()) 
				fileText = LoadFile(fileName);
			return !fileText.IsEmpty();
		}
		
		String GetString(String var) {
			if (!IsAvailable()) 
				throw Exc(Format(t_("Impossible to read file '%s'"), fileName));
			
			String res;
			UVector<String> vars = Split(var, "/");
			if (vars.size() == 2)
				res = GetFASTVarPos(fileText, vars[1], vars[0], pos);
			else if (vars.size() == 1)		
				res = GetFASTVarPos(fileText, vars[0], Null, pos);
			else
				throw Exc(Format(t_("Wrong variable '%s' in GetString"), var));
			
			if (res == "")
				throw Exc(Format(t_("Unknown variable '%s' in GetString"), var));
			
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
		int GetInt(String var) {
			int ddata = ScanInt(GetString(var));
			if (!IsNum(ddata))
				throw Exc(Format(t_("Wrong variable '%s' in GetInt"), var));
			return ddata;
		}
		bool GetBool(String var) {
			String data = ToLower(GetString(var));
			if (data == "true")
				return true;
			if (data == "false")
				return false;
			int idata = ScanInt(data);
			if (idata == 1)
				return true;
			if (idata == 0)
				return false;
			throw Exc(Format(t_("Wrong variable '%s' in GetBool"), var));
		}
		double GetMatrixVal(String var, int row, int col) {
			if (!IsAvailable()) 
				throw Exc(Format(t_("Impossible to read file '%s'"), fileName));
			return GetFASTMatrixVal(fileText, var, row, col);
		}
		Eigen::MatrixXd GetMatrix(String var, int rows, int cols) {
			if (!IsAvailable()) 
				throw Exc(Format(t_("Impossible to read file '%s'"), fileName));
			return GetFASTMatrix(fileText, var, rows, cols);
		}
		void SetMatrixVal(String var, int row, int col, double val) {
			if (!IsAvailable()) 
				throw Exc(Format(t_("Impossible to read file '%s'"), fileName));
			int posIni, posEnd;
			GetMatrixIds(var, row, col, posIni, posEnd);
			
			int delta = posEnd-posIni-1;
			
			isChanged = true;
			fileText = fileText.Left(posIni) + S(" ") + FDS(val, delta, true) + fileText.Mid(posEnd);
		}
		
		UVector<UVector<String>> GetFASTArray(String var) {
			if (!IsAvailable()) 
				throw Exc(Format(t_("Impossible to read file '%s'"), fileName));
			return ::GetFASTArray(fileText, var);
		}
		
		void SetString(String var, String val) {
			val = S("\"") + val + S("\"");
			SetString0(var, val);
		}
		void SetInt(String var, int val) {
			SetString0(var, FormatInt(val));
		}
		void SetDouble(String var, double val) {
			SetString0(var, FDS(val, 10));
		}
		void SetBool(String var, bool val) {
			SetString0(var, val ? "True" : "False");
		}
		
		int pos = 0;
		
	private:
		void SetString0(String var, String val) {
			if (fileText.IsEmpty()) {
				fileText = LoadFile(fileName);
				if (fileText.IsEmpty())
					throw Exc(Format(t_("Impossible to read file '%s'"), fileName));
			}
			UVector<String> vars = Split(var, "/");
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
			GetFASTMatrixIds(fileText, var, row, col, posIni, posEnd);
		}
	};
	
public:
	File fast, elastodyn, hydrodyn, inflowfile, dlldat, subdyn;
	String folderCase;
	String log;
	String fstFile;
	FastOut out;
	
	UVector<String> pointNames;
	UVector<Point3D> points;
};

class FASTCaseDecay : public FASTCase {
public:
	void Init(BasicBEM::DOF dof, double time, double x, double y, double z, double rx, double ry, double rz);
	virtual bool Postprocess();
	
	BasicBEM::DOF dof;
	//double T, r2;
};

double GetDecayPeriod(FastOut &fst, BasicBEM::DOF dof, double &r2);

double GetRAO(const VectorXd &data, const VectorXd &time, double T, bool onlyFFT, double r2Max);
void GetWaveRegularAmplitude(const FastOut &dataFast, double &T, double &A);
	
#endif
