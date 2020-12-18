#ifndef _BEMRosetta_BEMRosetta_cl_FastOut_h_
#define _BEMRosetta_BEMRosetta_cl_FastOut_h_


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
		
	double GetVal(double time, int col) const;
	inline double GetVal(int idtime, int col) const	{return dataOut[col][idtime];}
	inline const Vector<double> &GetVal(int col)	{return dataOut[col];}
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
		
	protected:
		FastOut *dataFast = 0;
	};
	
	struct CalcParams {
		String name, units;
		CalcParam *calc;
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
			idroll = dataFast->GetCol("PtfmRoll");
			idpitch = dataFast->GetCol("PtfmPitch");
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
			idsurge = dataFast->GetCol("PtfmSurge");
			idsway = dataFast->GetCol("PtfmSway");
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
			idheave = dataFast->GetCol("PtfmHeave");
			idpitch = dataFast->GetCol("PtfmPitch");
			idroll = dataFast->GetCol("PtfmRoll");
			idyaw = dataFast->GetCol("PtfmYaw");
		}
		virtual double Calc(int idtime) {
			double heave = dataFast->GetVal(idtime, idheave);
			double pitch = dataFast->GetVal(idtime, idpitch);
			double roll = dataFast->GetVal(idtime, idroll);
			double yaw = dataFast->GetVal(idtime, idyaw);
			
			return Null;//FALTA
			
			
		}	
		int idheave = Null, idpitch = Null, idroll = Null, idyaw = Null;	
	} ptfmHeaveCB; 
	
		
	struct YawAccelParam : CalcParam {
		virtual void Init() {
			idax = dataFast->GetCol("YawBrTAxp");
			iday = dataFast->GetCol("YawBrTAyp");
			idaz = dataFast->GetCol("YawBrTAzp");
		}
		virtual double Calc(int idtime) {
			double ax = dataFast->GetVal(idtime, idax);
			double ay = dataFast->GetVal(idtime, iday);
			double az = dataFast->GetVal(idtime, idaz);
			
			return sqrt(ax*ax + ay*ay + az*az);
		}	
		int idax = Null, iday = Null, idaz = Null;	
	} yawaccel;
};


#endif
