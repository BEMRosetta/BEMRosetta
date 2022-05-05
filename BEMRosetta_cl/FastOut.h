// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#ifndef _BEMRosetta_BEMRosetta_cl_FastOut_h_
#define _BEMRosetta_BEMRosetta_cl_FastOut_h_


bool FindHydrodyn(String path, double &ptfmCOBxt, double &ptfmCOByt);
	
class FastOut {
public:
	FastOut();
	
	static UVector<String> GetFilesToLoad(String path);
	static String GetFileToLoad(String fileName);
	
	bool Load(String fileName);
	bool Save(String fileName);
	
	void Clear();
	int GetCol(String param) const;
	int FindCol(String param) const;
	UVector<int> FindParameterMatch(String param) const;
	
	const String &GetParameter(int id) const	{return parameters[id];}
	const String &GetUnit(int id) const			{return units[id];}
	int GetParameterCount() const				{return parameters.size();}
	 
	SortedIndex<String> GetParameterList(String filter = ""); 
	SortedIndex<String> GetUnitList(String filter = "");
	SortedVectorMap<String, String> GetList(String filterParam = "", String filterUnits = "");
		
	double GetVal(double time, int idparam) const;
	inline double GetVal(int idtime, int idparam) const	  {return dataOut[idparam][idtime];}
	inline const UVector<double> &GetUVector(int idparam) {return dataOut[idparam];}
	inline const UVector<double> &GetUVector(String param) {
		static UVector<double> bad;
		UVector<int> ids = FindParameterMatch(param);
		if (ids.IsEmpty())
			return bad;
		else
			return GetUVector(ids[0]);
	}
	VectorXd GetVector(String param) {
		const UVector<double> &d = GetUVector(param);
		VectorXd ret = Map<const VectorXd>(d, d.size());
		return ret;
	}
	
	void AddVal(int idparam, double val);
	
	int GetIdTime(double time) const;
	double GetTimeInit() const	{return dataOut[0][0];}
	double GetTimeEnd()	 const	{return dataOut[0][size()-1];}
	int size() const			{return dataOut[0].size();}			
	bool IsEmpty() const		{return dataOut.IsEmpty();}	
	
	String GetLastFile()		{return lastFile;}
	
	int ColFairlead(int i) const{return GetCol(Format("T[%d]", i-1));}
		
	int GetNumFairlead() const {
		for (int i = 0; i < 100; ++i) {
			if (ColFairlead(i) < 0) 
				return i;
		}
		return Null;
	}
	
	void Serialize(Stream& s) {
        s % dataOut % parameters % units;
    }
	
	struct CalcParam {
		void Init0(FastOut *_dataFast) {dataFast = _dataFast;}
		virtual void Init() = 0;
		virtual double Calc(int idt) = 0;
		bool IsEnabled()	{return enabled;}
		
	protected:
		FastOut *dataFast = nullptr;
		bool enabled = true;
	};
	
	struct CalcParams {
		String name, units;
		CalcParam *calc = nullptr;
	};
	
	void AddParam(String name, String units, CalcParam &calc) {
		CalcParams &c = calcParams.Add();
		c.name = name;
		c.units = units;
		c.calc = &calc;
	}
	
	int AddParam(String name, String unit) {
		parameters << name;
		units << unit;
		dataOut.Add();
		return dataOut.size() - 1;
	}
	
	int GetParamCount() {return parameters.size();}
 
	UVector<String> parameters, units;
	UVector<UVector <double> > dataOut;
	UArray<CalcParams> calcParams;

private:
	bool LoadOut(String fileName);
	bool LoadOutb(String fileName);
	bool SaveOut(String fileName);
	void AfterLoad();

	String lastFile;
	Time lastTime;
	
	// shear, bending moment
	
	struct TiltParam : CalcParam {
		virtual void Init() {
			idroll = dataFast->FindCol("PtfmRoll");
			idpitch = dataFast->FindCol("PtfmPitch");
			if (IsNull(idroll) || IsNull(idpitch))
				enabled = false;
		}
		virtual double Calc(int idtime) {
			double roll = dataFast->GetVal(idtime, idroll);
			double pitch = dataFast->GetVal(idtime, idpitch);
			
			return sqrt(roll*roll + pitch*pitch);
		}	
		int idroll = Null, idpitch = Null;	
	} ptfmtilt; 

	struct ShiftParam : CalcParam {
		virtual void Init() {
			idsurge = dataFast->FindCol("PtfmSurge");
			idsway = dataFast->FindCol("PtfmSway");
			if (IsNull(idsurge) || IsNull(idsway))
				enabled = false;
		}
		virtual double Calc(int idtime) {
			double surge = dataFast->GetVal(idtime, idsurge);
			double sway = dataFast->GetVal(idtime, idsway);
			
			return sqrt(surge*surge + sway*sway);
		}	
		int idsurge = Null, idsway = Null;	
	} ptfmshift; 

	struct HeaveCBParam : CalcParam {
		virtual void Init() {
			idheave = dataFast->FindCol("PtfmHeave");
			idpitch = dataFast->FindCol("PtfmPitch");
			idroll = dataFast->FindCol("PtfmRoll");
			idyaw = dataFast->FindCol("PtfmYaw");	
			if (IsNull(idroll) || IsNull(idpitch) || IsNull(idheave) || IsNull(idyaw)) 
				enabled = false;
			else {
				String folder = GetFileFolder(dataFast->GetLastFile());
				if (!FindHydrodyn(folder, ptfmCOBxt, ptfmCOByt)) 
					ptfmCOBxt = ptfmCOByt = Null;
			}
		}
		virtual double Calc(int idtime) {
			double heave = dataFast->GetVal(idtime, idheave);
			if (IsNull(ptfmCOBxt))
				return heave;
	
			double pitch = dataFast->GetVal(idtime, idpitch);
			double roll = dataFast->GetVal(idtime, idroll);
			//double yaw = dataFast->GetVal(idtime, idyaw);
			
			pitch = ToRad(pitch);
			roll = ToRad(roll);
			//yaw = ToRad(yaw);

			return heave - sin(pitch)*ptfmCOBxt + cos(pitch)*sin(roll)*ptfmCOByt;
		}	
		int idheave = Null, idpitch = Null, idroll = Null, idyaw = Null;
		double ptfmCOBxt = Null, ptfmCOByt = Null;	
	} ptfmHeaveCB; 
	
	struct TwrBsShearParam : CalcParam {
		virtual void Init() {
			idx = dataFast->FindCol("TwrBsFxt");
			idy = dataFast->FindCol("TwrBsFyt");
			if (IsNull(idx) || IsNull(idy))
				enabled = false;
		}
		virtual double Calc(int idtime) {
			double fx = dataFast->GetVal(idtime, idx);
			double fy = dataFast->GetVal(idtime, idy);
			
			return sqrt(fx*fx + fy*fy);
		}	
		int idx = Null, idy = Null;	
	} twrBsShear;
	
	struct TwrBsBendParam : CalcParam {
		virtual void Init() {
			idx = dataFast->FindCol("TwrBsMxt");
			idy = dataFast->FindCol("TwrBsMyt");
			if (IsNull(idx) || IsNull(idy))
				enabled = false;
		}
		virtual double Calc(int idtime) {
			double mx = dataFast->GetVal(idtime, idx);
			double my = dataFast->GetVal(idtime, idy);
			
			return sqrt(mx*mx + my*my);
		}	
		int idx = Null, idy = Null;	
	} twrBsBend;	

	struct YawBrShearParam : CalcParam {
		virtual void Init() {
			idx = dataFast->FindCol("YawBrFxn");
			idy = dataFast->FindCol("YawBrFyn");
			if (IsNull(idx) || IsNull(idy))
				enabled = false;
		}
		virtual double Calc(int idtime) {
			double fx = dataFast->GetVal(idtime, idx);
			double fy = dataFast->GetVal(idtime, idy);
			
			return sqrt(fx*fx + fy*fy);
		}	
		int idx = Null, idy = Null;	
	} yawBrShear;

	struct YawBrBendParam : CalcParam {
		virtual void Init() {
			idx = dataFast->FindCol("YawBrMxn");
			idy = dataFast->FindCol("YawBrMyn");
			if (IsNull(idx) || IsNull(idy))
				enabled = false;
		}
		virtual double Calc(int idtime) {
			double mx = dataFast->GetVal(idtime, idx);
			double my = dataFast->GetVal(idtime, idy);
			
			return sqrt(mx*mx + my*my);
		}	
		int idx = Null, idy = Null;	
	} yawBrBend;

	struct RootShear1Param : CalcParam {
		virtual void Init() {
			idx = dataFast->FindCol("RootFxc1");
			idy = dataFast->FindCol("RootFyc1");
			if (IsNull(idx) || IsNull(idy))
				enabled = false;
		}
		virtual double Calc(int idtime) {
			double fx = dataFast->GetVal(idtime, idx);
			double fy = dataFast->GetVal(idtime, idy);
			
			return sqrt(fx*fx + fy*fy);
		}	
		int idx = Null, idy = Null;	
	} rootShear1;

	struct RootBend1Param : CalcParam {
		virtual void Init() {
			idx = dataFast->FindCol("RootMxc1");
			idy = dataFast->FindCol("RootMyc1");
			if (IsNull(idx) || IsNull(idy))
				enabled = false;
		}
		virtual double Calc(int idtime) {
			double mx = dataFast->GetVal(idtime, idx);
			double my = dataFast->GetVal(idtime, idy);
			
			return sqrt(mx*mx + my*my);
		}	
		int idx = Null, idy = Null;	
	} rootBend1;

	struct RootShear2Param : CalcParam {
		virtual void Init() {
			idx = dataFast->FindCol("RootFxc2");
			idy = dataFast->FindCol("RootFyc2");
			if (IsNull(idx) || IsNull(idy))
				enabled = false;
		}
		virtual double Calc(int idtime) {
			double fx = dataFast->GetVal(idtime, idx);
			double fy = dataFast->GetVal(idtime, idy);
			
			return sqrt(fx*fx + fy*fy);
		}	
		int idx = Null, idy = Null;	
	} rootShear2;

	struct RootBend2Param : CalcParam {
		virtual void Init() {
			idx = dataFast->FindCol("RootMxc2");
			idy = dataFast->FindCol("RootMyc2");
			if (IsNull(idx) || IsNull(idy))
				enabled = false;
		}
		virtual double Calc(int idtime) {
			double mx = dataFast->GetVal(idtime, idx);
			double my = dataFast->GetVal(idtime, idy);
			
			return sqrt(mx*mx + my*my);
		}	
		int idx = Null, idy = Null;	
	} rootBend2;
	
	struct RootShear3Param : CalcParam {
		virtual void Init() {
			idx = dataFast->FindCol("RootFxc3");
			idy = dataFast->FindCol("RootFyc3");
			if (IsNull(idx) || IsNull(idy))
				enabled = false;
		}
		virtual double Calc(int idtime) {
			double fx = dataFast->GetVal(idtime, idx);
			double fy = dataFast->GetVal(idtime, idy);
			
			return sqrt(fx*fx + fy*fy);
		}	
		int idx = Null, idy = Null;	
	} rootShear3;

	struct RootBend3Param : CalcParam {
		virtual void Init() {
			idx = dataFast->FindCol("RootMxc3");
			idy = dataFast->FindCol("RootMyc3");
			if (IsNull(idx) || IsNull(idy))
				enabled = false;
		}
		virtual double Calc(int idtime) {
			double mx = dataFast->GetVal(idtime, idx);
			double my = dataFast->GetVal(idtime, idy);
			
			return sqrt(mx*mx + my*my);
		}	
		int idx = Null, idy = Null;	
	} rootBend3;
	
	struct NcIMUTAParam : CalcParam {
		virtual void Init() {
			idx = dataFast->FindCol("NcIMUTAxs");
			idy = dataFast->FindCol("NcIMUTAys");
			idz = dataFast->FindCol("NcIMUTAzs");
			if (IsNull(idx) || IsNull(idy) || IsNull(idz))
				enabled = false;
		}
		virtual double Calc(int idtime) {
			double ax = dataFast->GetVal(idtime, idx);
			double ay = dataFast->GetVal(idtime, idy);
			double az = dataFast->GetVal(idtime, idz);
			
			return sqrt(ax*ax + ay*ay + az*az);
		}	
		int idx = Null, idy = Null, idz = Null;	
	} ncIMUTA;
		
};


class FASTCase {
public:
	~FASTCase() {
		DeleteFolderDeepX(folderCase);
	}
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
	void Setup(String seed);
	void CreateFolderCase();
	
	void LoadOut() {
		out.Load(fstFile);
	}
	virtual void Postprocess() {};
	
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
			UVector<String> vars = Split(var, "/");
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
			if (!IsNum(ddata))
				throw Exc(Format(t_("Wrong variable '%s' in GetDouble"), var));
			return ddata;
		}
		double GetInt(String var) {
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
			return GetFASTMatrixVal(fileText, var, row, col);
		}
		Eigen::MatrixXd GetMatrix(String var, int rows, int cols) {
			if (fileText.IsEmpty()) {
				fileText = LoadFile(fileName);
				if (fileText.IsEmpty())
					throw Exc(Format(t_("Impossible to read file '%s'"), fileName));
			}
			return GetFASTMatrix(fileText, var, rows, cols);
		}
		void SetMatrixVal(String var, int row, int col, double val) {
			int posIni, posEnd;
			GetMatrixIds(var, row, col, posIni, posEnd);
			
			int delta = posEnd-posIni-1;
			
			fileText = fileText.Left(posIni) + S(" ") + FDS(val, delta, true) + fileText.Mid(posEnd);
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
	String fstFile;
	
public:
	File fast, elastodyn, hydrodyn;
	String folderCase;
	String log;
	FastOut out;
};

class FASTCaseDecay : public FASTCase {
public:
	void Init(BEM::DOF dof, double val, double time = 100);
	virtual void Postprocess();
	
	BEM::DOF dof;
	double T, r2;
};

double GetDecayPeriod(FastOut &fst, BEM::DOF dof, double &r2);
	
#endif
