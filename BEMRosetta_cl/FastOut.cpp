#include "BEMRosetta.h"

#include <Core/Core.h>
#include <SysInfo/SysInfo.h>
#include <Surface/Surface.h>

using namespace Upp;

#include "FastOut.h"

static int IsTabSpaceRet(int c) {
	if (c == '\t' || c == ' ' || c == '\r' || c == '\n')
		return true; 
	return false;
} 

FastOut::FastOut() {
	ptfmtilt.Init0(this);
	AddParam("PtfmTilt", "deg", ptfmtilt);
	ptfmshift.Init0(this);
	AddParam("PtfmShift", "m", ptfmshift);
	ptfmHeaveCB.Init0(this);
	AddParam("PtfmHeaveCB", "m", ptfmHeaveCB);
	twrBsShear.Init0(this);
	AddParam("TwrBsShear", "kN", twrBsShear);
	twrBsBend.Init0(this);
	AddParam("TwrBsBend", "kN-m", twrBsBend);
	yawBrShear.Init0(this);
	AddParam("YawBrShear", "kN", yawBrShear);
	yawBrBend.Init0(this);
	AddParam("YawBrBend", "kN-m", yawBrBend);	
	rootShear1.Init0(this);
	AddParam("RootShear1", "kN", rootShear1);
	rootShear2.Init0(this);
	AddParam("RootShear2", "kN", rootShear2);
	rootShear3.Init0(this);
	AddParam("RootShear3", "kN", rootShear3);
	rootBend1.Init0(this);
	AddParam("RootBend1", "kN-m", rootBend1);
	rootBend2.Init0(this);
	AddParam("RootBend2", "kN-m", rootBend2);
	rootBend3.Init0(this);
	AddParam("RootBend3", "kN-m", rootBend3);	
	ncIMUTA.Init0(this);
	AddParam("NcIMUTA", "m/s^2", ncIMUTA);
}

Vector<String> FastOut::GetFilesToLoad(String path) {
	Vector<String> ret;
	
	if (TrimBoth(path).IsEmpty())
		return ret;
	
	if (!DirectoryExists(path)) {
		if (FileExists(path))
			ret << GetFileToLoad(path);
		return ret;
	} 
	int64 sz = -1;
	String fileName;
	for (FindFile ff(AppendFileName(path, "*.out*")); ff; ff++) {
		if (ff.IsFile()) { 
			String name = GetFileToLoad(ff.GetPath());
			if (!IsNull(name)) {
				if (GetFileExt(name) == ".outb") {
					ret << name;
					return ret;
				}
				int64 nsz = GetFileLength(name);
				if (nsz > sz) {
					sz = nsz;
					fileName = name;
				}
			}
		}
	}
	if (!fileName.IsEmpty()) {
		ret << fileName;
		return ret;
	}
	for (FindFile ff(AppendFileName(path, "*.*")); ff; ff++) 
		if (ff.IsFolder())
			ret.Append(GetFilesToLoad(ff.GetPath()));

	return ret;
}

String FastOut::GetFileToLoad(String fileName) {
	if (TrimBoth(fileName).IsEmpty())
		return Null;
		
	String strOut = ForceExt(fileName, ".out");
	String strOutB = ForceExt(fileName, ".outb");

	bool exOut = FileExists(strOut);
	bool exOutB = FileExists(strOutB);

	if (exOut && !exOutB) 
		return strOut;
	else if (!exOut && exOutB) 
		return strOutB;
	else if (exOut && exOutB) {
		Time tOut = FileGetTime(strOut);
		Time tOutB = FileGetTime(strOutB);
		if (abs(tOut - tOutB) < 5 || tOutB > tOut) 
			return strOutB;
		else 
			return strOut;
	}
	return Null;
}

int FastOut::Load(String fileName) {	
	Time actualTime = FileGetTime(fileName);
	if (lastFile == fileName && actualTime - lastTime < 5) // Only loads if file is 5 sec older
		return true;
	 
	lastFile = fileName;
	lastTime = actualTime;
	
	String ext = GetFileExt(fileName);
	bool ret = false;
	if (ext == ".out")
		ret = LoadOut(fileName);
	else if (ext == ".outb")
		ret = LoadOutb(fileName);
			
	if (ret)
		AfterLoad();
	return ret;
}
	
bool FastOut::LoadOut(String fileName) {
	String raw = LoadFileBOM(fileName);
	if (raw.IsEmpty()) 
		throw Exc(Format("Problem reading '%s'", fileName)); 
	
	Clear();
	bool begin = false;
	int row = -1;
	int pos = 0, npos = 0;
	int numCol;
	while (npos >= 0) {
		npos = raw.FindAfter("\n", pos);
		String line = raw.Mid(pos, npos-pos);
		Vector<String> fields = Split(line, IsTabSpaceRet, true);
		
		if (!begin) {
			if (!fields.IsEmpty() && fields[0] == "Time") {
				for (int c = 0; c < fields.size(); ++c) 
					parameters << fields[c];
				pos = npos;
				npos = raw.FindAfter("\n", pos);
				line = raw.Mid(pos, npos-pos);
				Vector<String> fields = Split(line, IsTabSpaceRet, true);
				for (int c = 0; c < fields.size(); ++c) 
					units << Replace(Replace(fields[c], "(", ""), ")", "");
				begin = true;
				if (parameters.size() != units.size()) 
					throw Exc("Number of parameters and units do not match");

				numCol = parameters.size();
				dataOut.SetCount(numCol+calcParams.size());
			}
		} else {
			row++;
			if (fields.IsEmpty())
				break;
			if (fields.size() != numCol) 
				throw Exc(Format("Number of values (%d) and parameters (%d) do not match in row %d", fields.size(), numCol, row));

			for (int c = 0; c < fields.size(); ++c) 
				dataOut[c] << ScanDouble(fields[c]);
		}
		pos = npos;
	}
	if (dataOut.size() > 0) {		// Size for calc. fields
		for (int i = numCol; i < dataOut.size(); ++i)
        	dataOut[i].SetCount(dataOut[0].size());
	}
	
	if (dataOut.IsEmpty()) 
		throw Exc(Format("Problem reading '%s'", fileName)); 

	return true;
}

bool FastOut::LoadOutb(String fileName) {
	Clear();
	
	enum FileType {WithTime = 1, WithoutTime, NoCompressWithoutTime, ChanLen_In};

	FileInBinary file(fileName);
	if (!file.IsOpen())
		return false;

	int ChanLen2 = 10;
	int16 FileID = file.ReadB<int16,2>();
	if (FileID == FileType::ChanLen_In) 
		ChanLen2 = file.ReadB<int16,2>();

	int32 NumChans = file.ReadB<int32,4>();
    int32 NumRecs = file.ReadB<int32,4>();

	double TimeScl, TimeOff, TimeOut1, TimeIncr;
    if (FileID == FileType::WithTime) {
        TimeScl = file.ReadB<double,8>(); 
        TimeOff = file.ReadB<double,8>();
    } else {
        TimeOut1 = file.ReadB<double,8>();
        TimeIncr = file.ReadB<double,8>();  
    }

	Buffer<float> ColScl(NumChans), ColOff(NumChans);
    file.ReadB(ColScl, 4*NumChans);	
    file.ReadB(ColOff, 4*NumChans);

	int32 LenDesc = file.ReadB<int32,4>();
    StringBuffer DescStrB(LenDesc);
    file.ReadB(DescStrB, LenDesc);
    String DescStr = DescStrB;
    
    if (FileID == FileType::NoCompressWithoutTime) 
		ChanLen2 = 15;

	parameters.SetCount(NumChans+1); 
	parameters[0] = "Time";
	Buffer<char> name(ChanLen2);
	for (int iChan = 0; iChan < NumChans+1; ++iChan) { 
		file.ReadB(name, ChanLen2); 	
        parameters[iChan] = TrimBoth(String(name, ChanLen2));
    }
    
	units.SetCount(NumChans+1);          		
	units[0] = "s";
	Buffer<char> unit(ChanLen2);
    for (int iChan = 0; iChan < NumChans+1; ++iChan) { 
        file.ReadB(unit, ChanLen2); 			
        units[iChan] = Replace(Replace(TrimBoth(String(unit, ChanLen2)), "(", ""), ")", "");
    }  
    
    int nPts = NumRecs*NumChans;           		   
    dataOut.SetCount(NumChans+1+calcParams.size());
    for (int i = 0; i < dataOut.size(); ++i)
        dataOut[i].SetCount(NumRecs);
    
    Buffer<int32> bufferTime;
    if (FileID == FileType::WithTime) {
        bufferTime.Alloc(NumRecs);
        file.ReadB(bufferTime, 4*NumRecs); 
    }
    
    Buffer<int16> bufferData(nPts);
    file.ReadB(bufferData, 2*nPts); 	
    int ip = 0;
    for (int idt = 0; idt < NumRecs; ++idt) 
    	for (int i = 1; i < NumChans+1; ++i) {
	    	double off = ColOff[i-1];
	    	double scl = ColScl[i-1];
	        dataOut[i][idt] = (bufferData[ip] - off)/scl;
	        ip++;
	}    
    
    if (FileID == FileType::WithTime) {
        for (int idt = 0; idt < NumRecs; ++idt)
            dataOut[0][idt] = ( bufferTime[idt] - TimeOff)/ TimeScl;
    } else {
        for (int idt = 0; idt < NumRecs; ++idt)
            dataOut[0][idt] = TimeOut1 + TimeIncr*idt;
    }
	return true;
}

void FastOut::AfterLoad() {
	for (CalcParams &c : calcParams) {
		if (!c.calc)
			throw Exc("Unexpected error in AfterLoad()");
		c.calc->Init();
		if (c.calc->IsEnabled()) {
			parameters << c.name;
			units << c.units;
			int id = parameters.size()-1;
			for (int idt = 0; idt < dataOut[id].size(); ++idt)
				dataOut[id][idt] = c.calc->Calc(idt);
		}
	}
}

void FastOut::Clear() {
	parameters.Clear();	
	units.Clear();		
	dataOut.Clear();
}

int FastOut::FindCol(String param) const {
	param = ToLower(param);
	for (int c = 0; c < parameters.size(); ++c) {
		if (ToLower(parameters[c]) == param)
			return c;
	}
	return Null;
}

Vector<int> FastOut::FindParameterMatch(String param) const {
	param = ToLower(param);
	Vector<int> ret;
	for (int c = 0; c < parameters.size(); ++c) {
		if (PatternMatch(param, ToLower(parameters[c])))
			ret << c;
	}
	return ret;
}

int FastOut::GetCol(String param) const {
	int ret = FindCol(param);
	if (IsNull(ret))
		throw Exc(Format("Parameter '%s' not found", param));	
	return ret;
}

double FastOut::GetVal(double time, int col) const {
	int idtime = GetIdTime(time);
	if (IsNull(idtime))
		return Null;
	return GetVal(idtime, col);
}

int FastOut::GetIdTime(double time) const {
	if (time < 0)
		return Null;	
	for (int r = 0; r < dataOut[0].size(); ++r) {
		if (dataOut[0][r] >= time)
			return r;
	}
	return Null;
}

SortedIndex<String> FastOut::GetParameterList(String filter) {
	SortedIndex<String> list;
	
	if (filter.IsEmpty()) {
		for (int i = 0; i < GetParameterCount(); ++i) 
			list.FindAdd(GetParameter(i));
	} else {
		filter = "*" + ToLower(filter) + "*";
		for (int i = 0; i < GetParameterCount(); ++i) {
			String str = ToLower(GetParameter(i));
			if (PatternMatch(filter, str))
				list.FindAdd(GetParameter(i));
		}
	}
	//SortIndex(list);
	return list;//.PickKeys();
}

SortedIndex<String> FastOut::GetUnitList(String filter) {
	SortedIndex<String> list;
	
	if (filter.IsEmpty()) {
		for (int i = 0; i < GetParameterCount(); ++i) 
			list.FindAdd(GetUnit(i));
	} else {
		filter = "*" + ToLower(filter) + '*';
		for (int i = 0; i < GetParameterCount(); ++i) {
			String str = ToLower(GetUnit(i));
			if (PatternMatch(filter, str))
				list.FindAdd(GetUnit(i));
		}
	}
	//SortIndex(list);
	return list;//.PickKeys();
}

SortedVectorMap<String, String> FastOut::GetList(String filterParam, String filterUnits) {
	SortedVectorMap<String, String> list;
	
	if (filterParam.IsEmpty() && filterUnits.IsEmpty()) {
		for (int i = 0; i < GetParameterCount(); ++i) 
			list.Add(GetParameter(i), GetUnit(i));
	} else {
		filterParam = "*" + ToLower(filterParam) + '*';
		filterUnits = "*" + ToLower(filterUnits) + '*';
		for (int i = 0; i < GetParameterCount(); ++i) {
			String strParams = ToLower(GetParameter(i));
			String strUnits = ToLower(GetUnit(i));
			if (PatternMatch(filterParam, strParams) && 
				PatternMatch(filterUnits, strUnits))
				list.Add(GetParameter(i), GetUnit(i));
		}
	}
	return list;
}

bool FindHydrodyn(String path, double &ptfmCOBxt, double &ptfmCOByt) {
	for (FindFile ff(AppendFileName(path, "*.dat")); ff; ++ff) {
		if (ff.IsFile()) { 
			String str = LoadFile(ff.GetPath());
			String strx = GetFASTVar(str, "PtfmCOBxt", "");
			if (!IsNull(strx)) {
				ptfmCOBxt = ScanDouble(strx);
				ptfmCOByt = ScanDouble(GetFASTVar(str, "PtfmCOByt", ""));
				return true;	
			}
		} else if (ff.IsFolder()) { 
			if (FindHydrodyn(ff.GetPath(), ptfmCOBxt, ptfmCOByt))
				return true;
		}
	}
	return false;
}