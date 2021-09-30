// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2021, the BEMRosetta author and contributors
#ifndef _BEMRosetta_BEMRosetta_cl_FastOut_h_
#define _BEMRosetta_BEMRosetta_cl_FastOut_h_


bool FindHydrodyn(String path, double &ptfmCOBxt, double &ptfmCOByt);
	
class FastOut {
public:
	FastOut();
	
	static Vector<String> GetFilesToLoad(String path);
	static String GetFileToLoad(String fileName);
	int Load(String fileName);
	
	void Clear();
	int GetCol(String param) const;
	int FindCol(String param) const;
	Upp::Vector<int> FindParameterMatch(String param) const;
	
	const String &GetParameter(int id) const	{return parameters[id];}
	const String &GetUnit(int id) const			{return units[id];}
	int GetParameterCount() const				{return parameters.size();}
	 
	SortedIndex<String> GetParameterList(String filter = ""); 
	SortedIndex<String> GetUnitList(String filter = "");
	SortedVectorMap<String, String> GetList(String filterParam = "", String filterUnits = "");
		
	double GetVal(double time, int idparam) const;
	inline double GetVal(int idtime, int idparam) const	{return dataOut[idparam][idtime];}
	inline const Vector<double> &GetVal(int idparam)	{return dataOut[idparam];}
	inline const Vector<double> &GetVal(String param) {
		static Vector<double> bad;
		Vector<int> ids = FindParameterMatch(param);
		if (ids.IsEmpty())
			return bad;
		else
			return GetVal(ids[0]);
	}
		
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
		FastOut *dataFast = 0;
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

	Upp::Vector<String> parameters, units;
	Upp::Vector<Upp::Vector <double> > dataOut;
	Upp::Array<CalcParams> calcParams;

private:
	bool LoadOut(String fileName);
	bool LoadOutb(String fileName);
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

	
#endif
