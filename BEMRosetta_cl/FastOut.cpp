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
	AddParam(ptfmtilt.Init0(this));
	AddParam(ptfmshift.Init0(this));
	AddParam(ptfmHeaveCB.Init0(this));
	AddParam(twrBsShear.Init0(this));
	AddParam(twrBsBend.Init0(this));
	AddParam(yawBrShear.Init0(this));
	AddParam(yawBrBend.Init0(this));
	AddParam(rootShear1.Init0(this));
	AddParam(rootShear2.Init0(this));
	AddParam(rootShear3.Init0(this));
	AddParam(rootBend1.Init0(this));
	AddParam(rootBend2.Init0(this));
	AddParam(rootBend3.Init0(this));
	AddParam(ncIMUTA.Init0(this));
	for (int i = 0; i < 10; ++i) {
		auto &f = fairTens.Add();
		f.Init00(i+1);
		AddParam(f.Init0(this));	
	}
	
	syn.Add("PtfmRoll", "roll");
	syn.Add("PtfmPitch", "pitch");
	syn.Add("PtfmYaw", "yaw");
	
	for (int i = 0; i < syn.size(); ++i) {
		syn[i] = ToLower(syn[i]);
		syn.SetKey(i, ToLower(syn.GetKey(i)));	
	}
}

UVector<String> FastOut::GetFilesToLoad(String path) {
	UVector<String> ret;
	
	if (TrimBoth(path).IsEmpty())
		return ret;
	
	if (!DirectoryExists(path)) {
		String name = GetFileToLoad(path);
		if (!IsNull(name) && !IsVoid(name)) 
			ret << name;
		return ret;
	} 
	String fileName;
	for (FindFile ff(AppendFileNameX(path, "*.out*")); ff; ff++) {
		if (ff.IsFile()) { 
			String name = GetFileToLoad(ff.GetPath());
			if (!IsNull(name) && !IsVoid(name)) {
				ret << name;
				return ret;
			}
		}
	}
	if (!fileName.IsEmpty()) {
		ret << fileName;
		return ret;
	}
	for (FindFile ff(AppendFileNameX(path, "*.*")); ff; ff++) // Search in inner other folders
		if (ff.IsFolder())
			ret.Append(GetFilesToLoad(ff.GetPath()));

	return ret;
}

String FastOut::GetFileToLoad(String fileName) {
	if (TrimBoth(fileName).IsEmpty())
		return "";
	if (::GetFileName(fileName).Find(".MD") >= 0)
		return String::GetVoid();
	
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
	if (FileExists(fileName))
		return fileName;
				
	return "";
}

bool FastOut::Load(String fileName) {
	String ext = GetFileExt(fileName);
	bool ret = false;
	if (ext == ".out")
		;
	else if (ext == ".outb")
		;
	else if (ext == ".csv")
		;
	else {
		fileName = ForceExt(fileName, ".outb");
		if (FileExists(fileName)) 
			ext = ".outb";
		else {
			fileName = ForceExt(fileName, ".out");	
			if (FileExists(fileName))
				ext = ".out";
			else
				return false;
		}
	}
	fileName = ForceExt(fileName, ext);
	
	if (ext == ".out")
		ret = LoadOut(fileName);
	else if (ext == ".outb")
		ret = LoadOutb(fileName);
	else if (ext == ".csv")
		ret = LoadCsv(fileName);
	
	this->fileName = fileName;
	
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

	if (dataOut.IsEmpty()) 
		throw Exc(Format("Problem reading '%s'", fileName)); 
	
	for (int i = numCol; i < dataOut.size(); ++i)	// Size for calc. fields
    	dataOut[i].SetCount(dataOut[0].size());

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

bool FastOut::LoadCsv(String fileName) {
	Clear();
	
	// Extracts the units from the header
	auto GetUnits = [=](String str, char begin, char end, String &param, String &unit)->bool {
		int idp = str.Find(begin);
		if (idp >= 0) {
			int idep = str.Find(end, idp+1);
			if (idep > 0) {
				unit = str.Mid(idp+1, idep-idp-1);
				param = str.Left(idp) + str.Mid(idep+1);
			} else {
				unit = str.Mid(idp+1);
				param = str.Left(idp);
			}
			return true;
		} else
			return false;
	};	
    			
	int numCol;
	
	char separator = ',';
	char decimalSign = '.';
	int fromRow = 0;
	bool onlyStrings = true;
	if (!ReadCSVFileByLine(fileName, [&](int row, UVector<Value> &result, String &line)-> bool {
		if (row == 0) {
			for (int i = 0; i < result.size(); ++i) {
				String param, unit;
				
				if (!GetUnits(result[i], '(', ')', param, unit))
					if (!GetUnits(result[i], '[', ']', param, unit))
						param = result[i];
				
				parameters << param;
				units << unit;
			}
			numCol = result.size();
			dataOut.SetCount(numCol+calcParams.size());
		} else {
			for (int i = 0; i < result.size(); ++i) 
				dataOut[i] << ScanDouble(AsString(result[i])); 
		}
		return true;
	}, separator, decimalSign, onlyStrings, fromRow))
		throw Exc(Format("Problem reading '%s'", fileName)); 
	
	if (dataOut.IsEmpty()) 
		throw Exc(Format("Problem reading '%s'", fileName)); 
	
	for (int i = numCol; i < dataOut.size(); ++i)	// Size for calc. fields
    	dataOut[i].SetCount(dataOut[0].size());
	
	return true;
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

void FastOut::AfterLoad() {
	parametersd.Clear();
	for (int i = 0; i < parameters.size(); ++i)
		parametersd << ToLower(parameters[i]);
	unitsd.Clear();
	for (int i = 0; i < units.size(); ++i)
		unitsd << ToLower(units[i]);
	
	for (CalcParams &c : calcParams) {
		if (!c.calc)
			throw Exc("Unexpected error in AfterLoad()");
		c.calc->Init();
		if (c.calc->IsEnabled()) {
			parameters << c.name;
			units << c.units;
			parametersd << ToLower(c.name);
			unitsd << ToLower(c.units);
			int id = parameters.size()-1;
			for (int idt = 0; idt < dataOut[id].size(); ++idt)
				dataOut[id][idt] = c.calc->Calc(idt);
		}
	}
}

void FastOut::Clear() {
	parameters.Clear();	
	units.Clear();		
	parametersd.Clear();	
	unitsd.Clear();	
	dataOut.Clear();
}

bool FastOut::IsEmpty() {
	return dataOut.IsEmpty();
}

int FastOut::GetParameterX(String param) const {
	param = ToLower(param);
	int id = Find(parametersd, param);
	if (id < 0) {
		id = syn.Find(param);
		if (id >= 0)
			id = Find(parametersd, syn[id]);
	}
	return id;
}

UVector<int> FastOut::FindParameterMatch(String param) const {
	param = ToLower(param);
	UVector<int> ret;
	for (int c = 0; c < parameters.size(); ++c) {
		if (PatternMatch(param, parametersd[c]))
			ret << c;
	}
	return ret;
}

int FastOut::GetParameter_throw(String param) const {
	int ret = GetParameterX(param);
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
	if (IsNull(time) || time < 0)
		return Null;	
	if (dataOut.size() == 0)
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

void Calc(const UArray<FastOut> &dataFast, const UVector<UVector<String>> &params, double start, bool fromEnd, double end, UVector<UVector<Value>> &table) {
	table.Clear();
	UVector<UVector<double>> fullData(params.size());
	for (const FastOut &fast : dataFast) {
		int idBegin = fast.GetIdTime(start);
		int num = fast.GetNumData();
		if (IsNull(idBegin) || idBegin >= num) 
			throw Exc(t_("Bad start time"));

		int idEnd = fromEnd ? fast.GetIdTime(fast.GetTimeEnd() - end) : fast.GetIdTime(end);
		if (!IsNull(idEnd)) {
			if (idEnd >= num) 
				throw Exc(t_("Bad end time"));
		} else
			idEnd = num-1;
		
		if (idBegin >= idEnd)
			throw Exc(t_("Begin has to be before end time"));
		
		UVector<Value> &t = table.Add();
		
		t << fast.GetFileName();			
		t << fast.GetVal(idBegin, 0);
		t << fast.GetVal(idEnd, 0);
		for (int ip = 0; ip < params.size(); ip++) {
			const UVector<String> &param = params[ip];
			if (param.size() < 3)
				throw Exc(t_("Wrong number of parameters"));
			String strpar = param[0];
			//String format = param[1];
			int id = fast.GetParameterX(strpar);
			if (id < 0) {
				for (int i = 2; i < param.size(); i++) 
					t << "";
			} else {
				const VectorXd data = fast.GetVector(id).segment(idBegin, idEnd - idBegin);
				UVector<double> ndata;
				Copy(data, ndata);
				fullData[ip].Append(ndata);
			
				for (int i = 2; i < param.size(); i++) {
					String stat = param[i];
					double val;
					if (stat == "mean") 
						val = data.mean();
					else if (stat == "min") 
						val = data.minCoeff();
					else if (stat == "max") 
						val = data.maxCoeff();
					else if (stat == "maxval") { 
						double mx = data.maxCoeff();
						double mn = data.minCoeff();
						if (abs(mx) > abs(mn))
							val = mx;
						else
							val = mn;
					} else if (stat == "maxmean") 
						val = data.maxCoeff() - data.mean();
					else if (stat == "minmean") 
						val = data.minCoeff() - data.mean();
					else if (stat == "std") 
						val = sqrt((data.array() - data.mean()).square().sum() / (data.size() - 1));
					else
						throw Exc(Format(t_("Unknown '%s' statistic in parameter '%s'"), stat, strpar));
					
					t << val;//Format("%" + format, val);
				}
			}
		}
	}
	if (dataFast.size() > 1) {
		UVector<Value> &t = table.Add();
		t << t_("Total");
		t << t_("-");
		t << t_("-");
		VectorXd data(table.size()-1), time(table.size()-1);
		int col = 3;
		for (int ip = 0; ip < params.size(); ip++) {
			const UVector<String> &param = params[ip];
			String strpar = param[0];
			String format = param[1];
			for (int i = 2; i < param.size(); i++) {
				for (int row = 0; row < table.size()-1; ++row) {
					data[row] = double(table[row][col]);
					time[row] = double(table[row][2]) - double(table[row][1]);
				}
				col++;
				String stat = param[i];
				double val;
				if (stat == "mean") 
					val = (data.array()*time.array()).sum() / time.sum();
				else if (stat == "min") 
					val = data.minCoeff();
				else if (stat == "max") 
					val = data.maxCoeff();
				else if (stat == "maxval") { 
					double mx = data.maxCoeff();
					double mn = data.minCoeff();
					if (abs(mx) > abs(mn))
						val = mx;
					else
						val = mn;
				} else if (stat == "maxmean") {
					double mean = (data.array()*time.array()).sum() / time.sum();
					val = data.maxCoeff() - mean;
				} else if (stat == "minmean") { 
					double mean = (data.array()*time.array()).sum() / time.sum();
					val = data.minCoeff() - mean;
				} else if (stat == "std") {
					VectorXd d = Map<VectorXd>(fullData[ip], fullData[ip].size());	
					val = sqrt((d.array() - d.mean()).square().sum() / (d.size() - 1));
				} else
					throw Exc(Format(t_("Unknown '%s' statistic in parameter '%s'"), stat, strpar));
				
				t << val;//Format("%" + format, val);
			}
		}
	}
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

double GetDecayPeriod(FastOut &fst, BEM::DOF dof, double &r2, double &damp) {
	r2 = damp = Null;
	
	String param = "ptfm" + S(BEM::strDOFtext[dof]);
	int id = fst.GetParameterX(param);
	if (id < 0) 
		throw Exc(Format("Param. %s not found in %s", param, fst.GetFileName()));

	UVector<int> idsx, idsy, idsFixed;
	VectorVectorY<double> vect;
	vect.Init(fst.dataOut, 0, id, idsx, idsy, idsFixed, false);
	
	DampedSinEquation eq;
	ExplicitEquation::FitError err = eq.Fit(vect, r2);
	if (err == ExplicitEquation::NoError) {
		damp = eq.GetDampingRatio();
		return eq.GetPeriod();
	} else
		return Null;
}