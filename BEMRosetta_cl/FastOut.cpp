// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"

#include <Core/Core.h>
#include <SysInfo/SysInfo.h>
#include <Surface/Surface.h>
#include <ScatterDraw/Equation.h>

using namespace Upp;

#include "FastOut.h"

static int IsTabSpaceRet(int c) {
	if (c == '\t' || c == ' ' || c == '\r' || c == '\n')
		return true; 
	return false;
} 

FastOut::FastOut() {
	ptfmtilt.Init0(this);
	AddParam("+PtfmTilt", "deg", ptfmtilt);
	ptfmshift.Init0(this);
	AddParam("+PtfmShift", "m", ptfmshift);
	ptfmHeaveCB.Init0(this);
	AddParam("+PtfmHeaveCB", "m", ptfmHeaveCB);
	twrBsShear.Init0(this);
	AddParam("+TwrBsShear", "kN", twrBsShear);
	twrBsBend.Init0(this);
	AddParam("+TwrBsBend", "kN-m", twrBsBend);
	yawBrShear.Init0(this);
	AddParam("+YawBrShear", "kN", yawBrShear);
	yawBrBend.Init0(this);
	AddParam("+YawBrBend", "kN-m", yawBrBend);	
	rootShear1.Init0(this);
	AddParam("+RootShear1", "kN", rootShear1);
	rootShear2.Init0(this);
	AddParam("+RootShear2", "kN", rootShear2);
	rootShear3.Init0(this);
	AddParam("+RootShear3", "kN", rootShear3);
	rootBend1.Init0(this);
	AddParam("+RootBend1", "kN-m", rootBend1);
	rootBend2.Init0(this);
	AddParam("+RootBend2", "kN-m", rootBend2);
	rootBend3.Init0(this);
	AddParam("+RootBend3", "kN-m", rootBend3);	
	ncIMUTA.Init0(this);
	AddParam("+NcIMUTA", "m/s^2", ncIMUTA);
}

UVector<String> FastOut::GetFilesToLoad(String path) {
	UVector<String> ret;
	
	if (TrimBoth(path).IsEmpty())
		return ret;
	
	if (!DirectoryExists(path)) {
		if (FileExists(path))
			ret << GetFileToLoad(path);
		return ret;
	} 
	int64 sz = -1;
	String fileName;
	for (FindFile ff(AppendFileNameX(path, "*.out*")); ff; ff++) {
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
	for (FindFile ff(AppendFileNameX(path, "*.*")); ff; ff++) 
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

bool FastOut::Load(String fileName) {
	String ext = GetFileExt(fileName);
	bool ret = false;
	if (ext == ".out")
		ret = LoadOut(fileName);
	else if (ext == ".outb")
		ret = LoadOutb(fileName);
	else {
		fileName = ForceExt(fileName, ".outb");
		if (FileExists(fileName))
			ret = LoadOutb(fileName);
		else {
			fileName = ForceExt(fileName, ".out");
			if (FileExists(fileName))
				ret = LoadOut(fileName);
			return false;
		}
	}
	lastFile = fileName;
	
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
		if (npos < 0) 					// The last row in the file
			npos = raw.GetCount();		
		String line = raw.Mid(pos, npos-pos);
		UVector<String> fields = Split(line, IsTabSpaceRet, true);
		
		if (!begin) {
			if (!fields.IsEmpty() && ToLower(fields[0]) == "time") {
				for (int c = 0; c < fields.size(); ++c) 
					parameters << fields[c];
				pos = npos;
				npos = raw.FindAfter("\n", pos);
				line = raw.Mid(pos, npos-pos);
				UVector<String> fields = Split(line, IsTabSpaceRet, true);
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

			for (int c = 0; c < fields.size(); ++c) 
				dataOut[c] << ScanDouble(fields[c]);
			for (int c = fields.size(); c < numCol; ++c) 
				dataOut[c] << Null;
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

bool FastOut::Save(String fileName, String type, String sep) {
	if (type.IsEmpty())
		type = GetFileExt(fileName);
	if (type == ".out")
		return SaveOut(fileName);
	else if (type == ".csv")
		return SaveCsv(fileName, sep);

	return false;
}
	
bool FastOut::SaveOut(String fileName) {
	String data;
	
	data << "\n\n\n\n\n\n";
	for (int i = 0; i < parameters.size(); ++i) {
		if (i > 0)
			data << "\t";
		data << parameters[i];
	}
	data << "\n";
	for (int i = 0; i < units.size(); ++i) {
		if (i > 0)
			data << "\t";
		data << "(" << units[i] << ")";
	}
	data << "\n";
	for (int idtime = 0; idtime < dataOut[0].size(); ++idtime) {
		if (idtime > 0)
			data << "\n";
		for (int idparam = 0; idparam < dataOut.size(); ++idparam) {
			if (idparam > 0)
				data << "\t";
			if (dataOut[idparam].size() > idtime)
				data << FDS(dataOut[idparam][idtime], 10, true);
		}
	}
	return SaveFile(fileName, data);
}

bool FastOut::SaveCsv(String fileName, String sep) {
	String data;
	
	if (sep.IsEmpty())
		sep = ";";
	for (int i = 0; i < parameters.size(); ++i) {
		if (i > 0)
			data << sep;
		data << parameters[i] << " [" << units[i] << "]";
	}
	data << "\n";
	for (int idtime = 0; idtime < dataOut[0].size(); ++idtime) {
		if (idtime > 0)
			data << "\n";
		for (int idparam = 0; idparam < dataOut.size(); ++idparam) {
			if (idparam > 0)
				data << sep;
			if (dataOut[idparam].size() > idtime)
				data << FDS(dataOut[idparam][idtime], 10, true);
		}
	}
	return SaveFile(fileName, data);
}

bool FastOut::LoadOutb(String fileName) {
	Clear();
	
	enum FILETYPE {WithTime = 1, WithoutTime, NoCompressWithoutTime, ChanLen_In};

	FileInBinary file(fileName);
	if (!file.IsOpen())
		return false;

	int ChanLen2;
	int16 FileType = file.ReadB<int16,2>();
	if (FileType == FILETYPE::ChanLen_In) 
		ChanLen2 = file.ReadB<int16,2>();
	else
		ChanLen2 = 10;

	int32 NumChans = file.ReadB<int32,4>();
    int32 NumRecs = file.ReadB<int32,4>();

	double TimeScl, TimeOff, TimeOut1, TimeIncr;
    if (FileType == FILETYPE::WithTime) {
        TimeScl = file.ReadB<double,8>(); 
        TimeOff = file.ReadB<double,8>();
    } else {
        TimeOut1 = file.ReadB<double,8>();
        TimeIncr = file.ReadB<double,8>();  
    }

	Buffer<float> ChanNames(NumChans);
	Buffer<float> ChanUnits(NumChans);
	
	//Buffer<float> ColMax, ColMin;	
	Buffer<float> ColScl, ColOff;
	Buffer<int32> TmpTimeArray;
	if (FileType != FILETYPE::NoCompressWithoutTime) {
		//ColMax.Alloc(NumChans); 
		//ColMin.Alloc(NumChans);	
		ColOff.Alloc(NumChans);
		ColScl.Alloc(NumChans); 

		if (FileType == FILETYPE::WithTime) 
			TmpTimeArray.Alloc(NumRecs);		
	}

	if (FileType != FILETYPE::NoCompressWithoutTime) {
	    file.ReadB(ColScl, 4*NumChans);	
    	file.ReadB(ColOff, 4*NumChans);
	}

	int32 LenDesc = file.ReadB<int32,4>();
	
    StringBuffer DescStrB(LenDesc);
    file.ReadB(DescStrB, LenDesc);
    String DescStr = DescStrB;

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
    
    // End of header
    
    int nPts = NumRecs*NumChans;           		   
    dataOut.SetCount(NumChans+1+calcParams.size());
    for (int i = 0; i < dataOut.size(); ++i)
        dataOut[i].SetCount(NumRecs);
    
    Buffer<int32> bufferTime;
    if (FileType == FILETYPE::WithTime) {
        bufferTime.Alloc(NumRecs);
        file.ReadB(bufferTime, 4*NumRecs); 
    }
    
    Buffer<int16> bufferData;
    Buffer<double> bufferDataFloat;
    int ip = 0;
    if (FileType == FILETYPE::NoCompressWithoutTime) {
        bufferDataFloat.Alloc(nPts);
        file.ReadB(bufferDataFloat, 8*nPts); 
        for (int idt = 0; idt < NumRecs; ++idt) {
	    	for (int i = 1; i < NumChans+1; ++i) 
		        dataOut[i][idt] = bufferDataFloat[ip++];
        }
    } else {
	    bufferData.Alloc(nPts);
    	file.ReadB(bufferData, 2*nPts); 	
	    for (int idt = 0; idt < NumRecs; ++idt) {
	    	for (int i = 1; i < NumChans+1; ++i) 
		        dataOut[i][idt] = (bufferData[ip++] - ColOff[i-1])/ColScl[i-1];
	    }
    }
    
    if (FileType == FILETYPE::WithTime) {
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

bool FastOut::IsEmpty() {
	return dataOut.IsEmpty();
}

int FastOut::FindCol(String param) const {
	param = ToLower(param);
	for (int c = 0; c < parameters.size(); ++c) {
		if (ToLower(parameters[c]) == param)
			return c;
	}
	return Null;
}

UVector<int> FastOut::FindParameterMatch(String param) const {
	param = ToLower(param);
	UVector<int> ret;
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

void FastOut::SetVal(int idparam, double val) {
	if (idtime == 0)
		return;
		//throw Exc("SetVal idparam == 0 is reserved to time");
	if (idtime < 0)
		throw Exc("SetNextTime has to bec called before SetVal");
	if (idparam < 0)
		return;
	auto &data = dataOut[idparam];
	if (idtime < data.size()) {
		data[idtime] = val;
		return;
	}
	if (data.GetAlloc() == data.size()) 
		data.Reserve(data.size() + 10000);
	data << val;
}

void FastOut::SetNextTime(double time) {
	auto &data = dataOut[0];
	if (idtime >= 0) {
		if (data[idtime] == time)
			return;
		if (data[idtime] > time)
			throw Exc("SetNextTime time is lower than last");
	}
	idtime++;
	if (data.GetAlloc() == data.size()) 
		data.Reserve(data.size() + 10000);
	data << time;
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
	for (FindFile ff(AppendFileNameX(path, "*.dat")); ff; ++ff) {
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

void FASTCase::CreateFolderCase(String folder) {
	Time t = GetSysTime();
	
	std::random_device rd;
	int randomSeed = rd();
	std::default_random_engine re(randomSeed);
	std::mt19937 rng(randomSeed);
	std::uniform_int_distribution<int> gen(0, 999); 
	int rnd = gen(rng);
	
	String cas = Format("%4d%02d%02d_%02d%02d%02d_%03d", t.year, t.month, t.day, t.hour, t.minute, t.second, rnd);
	String base;
	if (folder.IsEmpty())
		base = AppendFileNameX(GetAppDataFolder(), "BEMRosetta", "FASTCases");
	else
		base = folder;
	folderCase = AppendFileNameX(base, cas);
	if (!RealizeDirectory(folderCase))
		throw Exc(Format("Impossible to create folder %s", folderCase));
}
	
void FASTCase::Setup(String seed, String folderCases) {
	CreateFolderCase(folderCases);

	String errorStr;
	DirectoryCopyX(seed, folderCase, false, "*.out*;*.txt;*.dbg", errorStr);
	if (!IsEmpty(errorStr))
		throw Exc(errorStr);
	
	FindFile ffpath(AppendFileNameX(folderCase, "*.fst"));
	if (!ffpath) 
		throw Exc(Format("No .fst file found in folder '%s'", seed));
	
	fstFile = ffpath.GetPath();
	
	Load(fstFile);
}

void FASTCaseDecay::Init(BEM::DOF _dof, double time, double x, double y, double z, double rx, double ry, double rz) {
	dof = _dof;
	
	fast.SetInt("CompInflow", 0);
	fast.SetInt("CompAero", 0);
	fast.SetInt("CompServo", 0);
	fast.SetInt("CompMooring", 0);
	fast.SetInt("OutFileFmt", 2);

	hydrodyn.SetInt("WaveMod", 0);
	hydrodyn.SetBool("WvDiffQTF", false);
	hydrodyn.SetBool("WvSumQTF", false);
	hydrodyn.SetInt("CurrMod", 0);
	hydrodyn.SetInt("PotMod", 1);
	hydrodyn.SetInt("ExctnMod", 0);
	hydrodyn.SetInt("RdtnMod", 1);
	hydrodyn.SetInt("MnDrift", 0);
	hydrodyn.SetInt("NewmanApp", 0);
	hydrodyn.SetInt("DiffQTF", 0);
	hydrodyn.SetInt("SumQTF", 0);

	elastodyn.SetDouble("Gravity", 9.81);
	elastodyn.SetBool("FlapDOF1", false);
	elastodyn.SetBool("FlapDOF2", false);
	elastodyn.SetBool("EdgeDOF", false);
	elastodyn.SetBool("TeetDOF", false);
	elastodyn.SetBool("DrTrDOF", false);
	elastodyn.SetBool("GenDOF", false);
	elastodyn.SetBool("YawDOF", false);
	elastodyn.SetBool("TwFADOF1", false);
	elastodyn.SetBool("TwFADOF2", false);
	elastodyn.SetBool("TwSSDOF1", false);
	elastodyn.SetBool("TwSSDOF2", false);
	elastodyn.SetDouble("RotSpeed", 0);
	elastodyn.SetDouble("NacYaw", 0);
	
	elastodyn.SetBool("PtfmSgDOF", true);
	elastodyn.SetBool("PtfmSwDOF", true);
	elastodyn.SetBool("PtfmHvDOF", true);
	elastodyn.SetBool("PtfmRDOF",  true);
	elastodyn.SetBool("PtfmPDOF",  true);
	elastodyn.SetBool("PtfmYDOF",  true);
	
	elastodyn.SetDouble("PtfmSurge", x);
	elastodyn.SetDouble("PtfmSway", y);
	elastodyn.SetDouble("PtfmHeave", z);
	elastodyn.SetDouble("PtfmRoll", rx);
	elastodyn.SetDouble("PtfmPitch", ry);
	elastodyn.SetDouble("PtfmYaw", rz);
	
	fast.SetDouble("DT", 0.025);
	fast.SetDouble("TMax", time);
}

bool FASTCaseDecay::Postprocess() {
	if (!LoadOut())
		return false;
	//T = GetDecayPeriod(out, dof, r2);
	return true;
}

double GetDecayPeriod(FastOut &fst, BEM::DOF dof, double &r2) {
	r2 = Null;
	
	String param = "ptfm" + S(BEM::strDOFtext[dof]);
	int id = fst.FindCol(param);
	if (IsNull(id)) 
		throw Exc(Format("Param. %s not found in %s", param, fst.GetLastFile()));

	UVector<int> idsx, idsy, idsFixed;
	VectorVectorY<double> vect;
	vect.Init(fst.dataOut, 0, id, idsx, idsy, idsFixed, false);
	
	DampedSinEquation eq;
	ExplicitEquation::FitError err = eq.Fit(vect, r2);
	if (err == ExplicitEquation::NoError)
		return 2*M_PI/abs(eq.GetCoeff(3));
	else
		return Null;
}