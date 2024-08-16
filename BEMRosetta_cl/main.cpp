// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#if !defined(flagGUI)

#include "BEMRosetta.h"
#include <SysInfo/Crash.h>

#ifdef PLATFORM_WIN32

#if defined(flagBEMR_TEST_DLL) || defined(flagBEMR_TEST_DLL_INTERNAL)	

#ifdef flagBEMR_TEST_DLL_INTERNAL
	#include "export.h"
#endif

CONSOLE_APP_MAIN
{		
#if defined(flagDEBUG) && defined(PLATFORM_WIN32) && !defined(flagBEMR_TEST_DLL)
	GetCrashHandler().Enable();
#endif
	String errorStr;
	try {
#if defined(flagBEMR_TEST_DLL) || defined(flagBEMR_TEST_DLL_INTERNAL)
		const UVector<String>& command = CommandLine();
		
		if (command.GetCount() < 2) 
			throw Exc("Please include in command line binary and BEMRosetta folders");
		
		String binFolder = command[0];
		String bemFolder = command[1];
		String installFolder = AFX(bemFolder, "install");
#endif
#ifdef flagBEMR_TEST_DLL
		FileDelete(AFX(binFolder, "libbemrosetta.exp"));
		FileDelete(AFX(binFolder, "libbemrosetta.lib"));
		
		Dl dll;		
		if (!dll.Load(AFX(binFolder, "libbemrosetta.dll")))
			throw Exc("DLL not found");

		DLLFunction(dll, const char *, DLL_Version, ());
		DLLFunction(dll, void, 		   DLL_ListFunctions, ());
		DLLFunction(dll, const char *, DLL_strListFunctions, ());
		DLLFunction(dll, const char *, DLL_strPythonDeclaration, ());
		DLLFunction(dll, int, 		   DLL_FAST_Load, (const char *));
		DLLFunction(dll, const char *, DLL_FAST_GetParameterName, (int));
		DLLFunction(dll, const char *, DLL_FAST_GetUnitName, (int));
		DLLFunction(dll, int, 		   DLL_FAST_GetParameterId, (const char *name));
		DLLFunction(dll, int, 		   DLL_FAST_GetParameterCount, ());
		DLLFunction(dll, int, 		   DLL_FAST_GetLen, ());		
		DLLFunction(dll, double,	   DLL_FAST_GetTimeStart, ());
		DLLFunction(dll, double,	   DLL_FAST_GetTimeEnd, ());
		DLLFunction(dll, double, 	   DLL_FAST_GetTime, (int idtime));		
		DLLFunction(dll, double, 	   DLL_FAST_GetData, (int idtime, int idparam));
		DLLFunction(dll, double, 	   DLL_FAST_GetAvg,   (int idparam, int idbegin, int idend));
		DLLFunction(dll, double, 	   DLL_FAST_GetArray, (int idparam, int idbegin, int idend, double **, int *));
		
		DLLFunction(dll, double, 	   DemoVectorPy_C, (const double *, int));
		
		DLLFunction(dll, int, 	   	   DLL_FAST_LoadFile, (const char *file));
		DLLFunction(dll, int, 	   	   DLL_FAST_SaveFile, (const char *file));
		DLLFunction(dll, int, 	   	   DLL_FAST_SetVar, (const char *name, const char *paragraph, const char *value));
		DLLFunction(dll, const char *, DLL_FAST_GetVar, (const char *name, const char *paragraph));
#endif

		UVector<double> dat = {1, 2, 3};
		//double res = DemoVectorPy_C(dat, 3);
			
		Cout() << "\nVersion: " << DLL_Version();
		Cout() << "\n\nDLL functions list:\n";
		String strList = DLL_strListFunctions();
		Cout() << strList;
#ifdef flagBEMR_TEST_DLL
		strList = "// BEMRosetta DLL functions list\n\n" + strList;	
		if (!SaveFile(AFX(binFolder, "libbemrosetta.txt"), strList))
			throw Exc(t_("Impossible to save DLL functions list file"));
#endif

		Cout() << "\n\nPython declarations:\n";
		String strPy = DLL_strPythonDeclaration();
		Cout() << strPy;
		
#if defined(flagBEMR_TEST_DLL)
		if (!SaveFile(AFX(binFolder, "libbemrosetta.py"), strPy))
			throw Exc(t_("Impossible to save Python declarations file"));
#elif defined(flagBEMR_TEST_DLL_INTERNAL)
		if (!SaveFile(AFX(GetDesktopFolder(), "libbemrosetta.py"), strPy))
			throw Exc(t_("Impossible to save Python declarations file"));
#endif

#if defined(COMPILER_MSC) && defined(flagBEMR_TEST_DLL)
		String wxs = LoadFile(AFX(installFolder, "BEMRosetta_master.wxs"));
		if (wxs.IsEmpty())
			throw Exc(t_("Installer definition file not found"));
		String strver = LoadFile(AFX(installFolder, "build.txt"));
		if (strver.IsEmpty())
			throw Exc(t_("Version file not found"));
		int ver = ScanInt(strver);
		if (IsNull(ver) || ver < 1)
			throw Exc(t_("Wrong version found"));
		ver++;
		if (!SaveFile(AFX(installFolder, "build.txt"), FormatInt(ver)))
			throw Exc(t_("Impossible to save version file"));
		wxs.Replace("$VERSION$", FormatInt(ver));
		if (!SaveFile(AFX(installFolder, "BEMRosetta.wxs"), wxs))
			throw Exc(t_("Impossible to save installer file"));
#endif

#if defined(flagBEMR_TEST_DLL) || defined(flagBEMR_TEST_DLL_INTERNAL)
		Cout() << "\n\nLoading FAST .out file";
		String outfile = AFX(bemFolder, "examples/fast.out/demo.outb");
		if (!DLL_FAST_Load(outfile))
			throw Exc(Format("Impossible to open file %s", outfile));
		
		Cout() << "\nFAST .out parameters:";
		for (int i = 0; i < DLL_FAST_GetParameterCount(); ++i)
			Cout() << Format(" %s[%s]", DLL_FAST_GetParameterName(i), DLL_FAST_GetUnitName(i));
		
		Cout() << "\nSimulation begins at " << DLL_FAST_GetTimeStart() << " and ends at " << DLL_FAST_GetTimeEnd();
		
		int idptfmheave = DLL_FAST_GetParameterId("ptfmHeave");
		int num;
		double *v;
		DLL_FAST_GetArray(idptfmheave, -1, -1, &v, &num);
		Cout() << "\nRead " << num << " heave values"; 

		double avg = 0, avgv = 0;
		for (int i = 0; i < DLL_FAST_GetLen(); ++i) {
			avg += DLL_FAST_GetData(i, idptfmheave);
			avgv += v[i];
		}
		Cout() << "\nptfmheave_avg = " << avg/DLL_FAST_GetLen();
		Cout() << "\nptfmheave_avg = " << avgv/DLL_FAST_GetLen();
		Cout() << "\nptfmheave_avg = " << DLL_FAST_GetAvg(idptfmheave, -1, -1);
		
		Cout() << "\n\nLoading InflowWind .dat file";
		String datfile = AFX(bemFolder, "examples/fast.out/InflowWind.dat");
		if (!DLL_FAST_LoadFile(datfile))
			throw Exc(Format("Impossible to open file %s", outfile));	
	
		String str;
		str = DLL_FAST_GetVar("WindVziList", "");				
		Cout() << "\nWindVziList: " << str;
		VERIFY(str == "119");
		str = DLL_FAST_GetVar("nx", "");				
		Cout() << "\nnx: " << str;
		VERIFY(str == "64");
		str = DLL_FAST_GetVar("Filename", "");				
		Cout() << "\nFilename: " << str;
		VERIFY(str == "\"unifWind.hh\"");
		str = DLL_FAST_GetVar("Filename", "================== Parameters for Uniform wind file");				
		Cout() << "\nFilename: " << str;
		VERIFY(str == "\"unifWind.hh\"");
		str = DLL_FAST_GetVar("Filename", "================== Parameters for Binary TurbSim");				
		Cout() << "\nFilename: " << str;
		VERIFY(str == "\"TurbSim.bts\"");
		DLL_FAST_SetVar("nx", "", "23");				
		str = DLL_FAST_GetVar("nx", "");				
		Cout() << "\nNew nx: " << str;
		VERIFY(str == "23");
		DLL_FAST_SetVar("Filename", "================== Parameters for Binary TurbSim", "\"New file\"");				
		str = DLL_FAST_GetVar("Filename", "================== Parameters for Binary TurbSim");				
		Cout() << "\nNew Filename: " << str;
		VERIFY(str == "\"New file\"");

	#ifdef flagDEBUG	
		DLL_FAST_SaveFile(AFX(GetDesktopFolder(), "InflowWind_test.dat"));
	#endif
#endif
	} catch (Exc e) {
		errorStr = e;
	} catch(const char *cad) {
		errorStr = cad;
	} catch(const std::string &e) {
		errorStr = e.c_str();	
	} catch (const std::exception &e) {
		errorStr = e.what();
	} catch(...) {
		errorStr = t_("Unknown error");
	}	
	if (!errorStr.IsEmpty()) {
		Cout() << "\n" << Format(t_("Problem found: %s"), errorStr);
		SetExitCode(-1);
	}
#ifdef flagDEBUG
	Cout() << "\nPress a key to end";
	ReadStdIn();
#endif
}

#endif

#if defined(flagBEMR_DLL) || defined(flagBEMR_TEST_DLL_INTERNAL)

#include "FastOut.h"
#include "export.h"
#include "export.brc"

FastOut &DLL_Fastout() {
	static FastOut fast;
	return fast;
}


void DLL_NoPrint() noexcept {
	CoutStreamX::NoPrint();
}
	
const char *DLL_Version() noexcept {
	static String version;
	version << __DATE__ << ", " << __TIME__;
	return version;	
}

const char *DLL_strListFunctions() noexcept {
	static String str;
	
	return str = CleanCFromDeclaration(String(BEMR_DLLexport, BEMR_DLLexport_length));
}

const char *DLL_strPythonDeclaration() noexcept {
	static String str;
	
	return str = GetPythonDeclaration("BEMRosetta", String(BEMR_DLLexport, BEMR_DLLexport_length));	
}

void DLL_ListFunctions() noexcept {
	Cout() << DLL_strListFunctions();	
}

int DLL_FAST_Load(const char *filename) noexcept {
	try {
		String ret = DLL_Fastout().Load(filename);
		if (ret.IsEmpty())
			return 1;
		else {
			CoutX() << Format("Error: %s", ret);
			return 0;
		}
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_Load()";
		return 0;
	}
}

const char *DLL_FAST_GetParameterName(int id) noexcept {
	static String ret;
	try {
		return ret = DLL_Fastout().GetParameter(id);
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetParameterName()";
		return ret = "Error";
	}
}

const char *DLL_FAST_GetUnitName(int id) noexcept {
	static String ret;
	try {
		return ret = DLL_Fastout().GetUnit(id);
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetUnitName()";
		return ret = "Error";
	}
}

int DLL_FAST_GetParameterId(const char *name) noexcept {
	try {
		UVector<int> p = DLL_Fastout().FindParameterMatch(name);
		if (p.IsEmpty())
			return -1;
		else
			return p[0];
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetParameterCount()";
		return Null;
	}
}

int DLL_FAST_GetParameterCount() noexcept {
	try {
		return DLL_Fastout().GetParameterCount();
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetParameterCount()";
		return Null;
	}
}

int DLL_FAST_GetLen() noexcept {
	try {
		return DLL_Fastout().GetNumData();
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetLen()";
		return Null;
	}
}

double DLL_FAST_GetTimeStart() noexcept {
	try {
		return DLL_Fastout().GetTimeStart();
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetTimeStart()";
		return Null;
	}
}

double DLL_FAST_GetTimeEnd() noexcept {
	try {
		return DLL_Fastout().GetTimeEnd();
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetTimeEnd()";
		return Null;
	}
}

double DLL_FAST_GetTime(int idtime) noexcept {
	return DLL_FAST_GetData(idtime, 0);
}

double DLL_FAST_GetData(int idtime, int idparam) noexcept {
	if (idtime < 0) {
		CoutX() << "Error in DLL_FAST_GetData() idtime < 0";
		return Null;
	}
	if (idtime >= DLL_Fastout().GetNumData()) {
		CoutX() << "Error in DLL_FAST_GetData() idtime >= time";
		return Null;
	}
	if (idparam < 0) {
		CoutX() << "Error in DLL_FAST_GetData() idparam < 0";
		return Null;
	}
	if (idparam >= DLL_Fastout().GetParameterCount()) {
		CoutX() << "Error in DLL_FAST_GetData() idparam >= num_params";
		return Null;
	}
		
	return DLL_Fastout().GetVal(idtime, idparam);
}

static void DLL_FAST_GetData(int idparam, int idbegin, int idend, VectorXd &data) {
	if (idparam < 0) 
		throw Exc("idparam < 0");
	if (idparam >= DLL_Fastout().GetParameterCount()) 
		throw Exc("idparam >= num_params");

	if (idbegin < 0)
		idbegin = 0;
	if (idend < 0)
		idend = DLL_Fastout().GetNumData()-1;
	
	if (idbegin > idend) 
		throw Exc("idbegin > idend");
		
	data = DLL_Fastout().GetVector(idparam).segment(idbegin, idend - idbegin + 1);
}

int DLL_FAST_GetArray(int idparam, int idbegin, int idend, double **data, int *num) noexcept {
	static VectorXd v;
	
	try {
		DLL_FAST_GetData(idparam, idbegin, idend, v);
		
		*num = int(v.size());
		*data = v.data();
		return 1;
	} catch (Exc e) {
		CoutX() << Format("Error in DLL_FAST_GetArray(): %s", e);
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetArray()";
	}
	return Null;	
}

double DLL_FAST_GetAvg(int idparam, int idbegin, int idend) noexcept {
	try {
		VectorXd data;
		
		DLL_FAST_GetData(idparam, idbegin, idend, data);
		
		return data.mean();
	} catch (Exc e) {
		CoutX() << Format("Error in DLL_FAST_GetAvg(): %s", e);
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetAvg()";
	}
	return Null;
}

double DLL_FAST_GetMax(int idparam, int idbegin, int idend) noexcept {
	try {
		VectorXd data;
		
		DLL_FAST_GetData(idparam, idbegin, idend, data);
		
		return data.maxCoeff();
	} catch (Exc e) {
		CoutX() << Format("Error in DLL_FAST_GetMax(): %s", e);
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetAvg()";
	}
	return Null;
}

double DLL_FAST_GetMin(int idparam, int idbegin, int idend) noexcept {
	try {
		VectorXd data;
		
		DLL_FAST_GetData(idparam, idbegin, idend, data);
		
		return data.minCoeff();
	} catch (Exc e) {
		CoutX() << Format("Error in DLL_FAST_GetMin(): %s", e);
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetAvg()";
	}
	return Null;
}

int DLL_IsNull(double val) noexcept {return IsNull(val);}

static String fastFileStr;
static String fastFileName;

int DLL_FAST_LoadFile(const char *file) noexcept {
	fastFileStr = LoadFile(file);
	fastFileName = file;
	return !fastFileStr.IsEmpty();
}

int DLL_FAST_SaveFile(const char *file) noexcept {
	bool ret;
	try {
		String sfile(file);
		if (!sfile.IsEmpty()) 
			fastFileName = sfile;
		ret = SaveFile(fastFileName, fastFileStr);
		
	} catch (Exc e) {
		SetConsoleColor(CONSOLE_COLOR::LTYELLOW);
		CoutX() << "\n" << "Error: " << e;
		SetConsoleColor(CONSOLE_COLOR::PREVIOUS);
		return false;
	}
	return ret;	
}

int DLL_FAST_SetVar(const char *name, const char *paragraph, const char *value) noexcept {
	try {
		SetFASTVar(fastFileStr, name, value, paragraph);
	} catch (Exc e) {
		SetConsoleColor(CONSOLE_COLOR::LTYELLOW);
		CoutX() << "\n" << "Error: " << e;
		SetConsoleColor(CONSOLE_COLOR::PREVIOUS);
		return false;
	}
	return true;
}

const char *DLL_FAST_GetVar(const char *name, const char *paragraph) noexcept {
	static String ret;

	try {
		ret = GetFASTVar(fastFileStr, name, paragraph);
	} catch (Exc e) {
		SetConsoleColor(CONSOLE_COLOR::LTYELLOW);
		CoutX() << "\n" << "Error: " << e;
		SetConsoleColor(CONSOLE_COLOR::PREVIOUS);
		return ret = "";
	}
	if (IsVoid(ret))
		return "";
	return ret;
}

double DemoVectorPy_C(const double *v, int num) noexcept {
    double res = 0;
    for (int i = 0; i < num; ++i) 
        res += v[i];
    return res;
}


#endif

#endif

#endif

#if defined(flagBEMR_CL)

	
CONSOLE_APP_MAIN {
	const UVector<String>& command = CommandLine();
	
	if (!ConsoleMain(command, false, PrintStatus))
		SetExitCode(1);
	
#ifdef flagDEBUG
	Cout() << "\nPress enter to end";
	ReadStdIn();
#endif
}
#endif



