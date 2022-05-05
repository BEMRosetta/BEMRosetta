// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#if !defined(flagGUI)

#include "BEMRosetta.h"

#ifdef PLATFORM_WIN32

#if defined(flagBEMR_TEST_DLL) || defined(flagBEMR_TEST_DLL_INTERNAL)	

#ifdef flagBEMR_TEST_DLL_INTERNAL
	#include "export.h"
#endif

CONSOLE_APP_MAIN
{	
	try {
#ifdef flagBEMR_TEST_DLL
		const UVector<String>& command = CommandLine();
		
		if (command.GetCount() < 2) 
			throw Exc("Please include in command line binary and BEMRosetta folders");
		
		String binFolder = command[0];
		String bemFolder = command[1];
		String installFolder = AppendFileNameX(bemFolder, "install");
		
		FileDelete(AppendFileNameX(binFolder, "libbemrosetta.exp"));
		FileDelete(AppendFileNameX(binFolder, "libbemrosetta.lib"));
		
		Dl dll;		
		if (!dll.Load(AppendFileNameX(binFolder, "libbemrosetta.dll")))
			throw Exc("Dll not found");

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
		DLLFunction(dll, double,	   DLL_FAST_GetTimeInit, ());
		DLLFunction(dll, double,	   DLL_FAST_GetTimeEnd, ());
		DLLFunction(dll, double, 	   DLL_FAST_GetTime, (int idtime));		
		DLLFunction(dll, double, 	   DLL_FAST_GetData, (int idtime, int idparam));
		DLLFunction(dll, double, 	   DLL_FAST_GetAvg, (const char *));
		
		DLLFunction(dll, double, 	   DemoVectorPy_C, (const double *, int));
		DLLFunction(dll, double, 	   DemoVectorC_Py, (double **, int *));
#endif

		
		UVector<double> dat = {1, 2, 3};
		double res = DemoVectorPy_C(dat, 3);
		double *v;
		int num;
		DemoVectorC_Py(&v, &num);
		
	
		Cout() << "\nVersion: " << DLL_Version();
		Cout() << "\n\nDLL functions list:\n";
		String strList = DLL_strListFunctions();
		Cout() << strList;
#ifdef flagBEMR_TEST_DLL
		strList = "// BEMRosetta DLL functions list\n\n" + strList;	
		if (!SaveFile(AppendFileNameX(binFolder, "libbemrosetta.txt"), strList))
			throw Exc(t_("Impossible to save DLL functions list file"));
#endif

		Cout() << "\n\nPython declarations:\n";
		String strPy = DLL_strPythonDeclaration();
		Cout() << strPy;
		
#if defined(flagBEMR_TEST_DLL)
		if (!SaveFile(AppendFileNameX(binFolder, "libbemrosetta.py"), strPy))
			throw Exc(t_("Impossible to save Python declarations file"));
#elif defined(flagBEMR_TEST_DLL_INTERNAL)
		if (!SaveFile(AppendFileNameX(GetDesktopFolder(), "libbemrosetta.py"), strPy))
			throw Exc(t_("Impossible to save Python declarations file"));
#endif

#if defined(COMPILER_MSC) && defined(flagBEMR_TEST_DLL)
		String wxs = LoadFile(AppendFileNameX(installFolder, "BEMRosetta_master.wxs"));
		if (wxs.IsEmpty())
			throw Exc(t_("Installer definition file not found"));
		String strver = LoadFile(AppendFileNameX(installFolder, "build.txt"));
		if (strver.IsEmpty())
			throw Exc(t_("Version file not found"));
		int ver = ScanInt(strver);
		if (IsNull(ver) || ver < 1)
			throw Exc(t_("Wrong version found"));
		ver++;
		if (!SaveFile(AppendFileNameX(installFolder, "build.txt"), FormatInt(ver)))
			throw Exc(t_("Impossible to save version file"));
		wxs.Replace("$VERSION$", FormatInt(ver));
		if (!SaveFile(AppendFileNameX(installFolder, "BEMRosetta.wxs"), wxs))
			throw Exc(t_("Impossible to save installer file"));
#endif

#ifdef flagBEMR_TEST_DLL		
		Cout() << "\n\nLoading FAST .out file";
		String outfile = AppendFileNameX(bemFolder, "examples/fast.out/demo.outb");
		if (!DLL_FAST_Load(outfile))
			throw Exc(Format("Impossible to open file %s", outfile));
		
		Cout() << "\nFAST .out parameters:";
		for (int i = 0; i < DLL_FAST_GetParameterCount(); ++i)
			Cout() << Format(" %s[%s]", DLL_FAST_GetParameterName(i), DLL_FAST_GetUnitName(i));
		
		Cout() << "\nSimulation begins at " << DLL_FAST_GetTimeInit() << " and ends at " << DLL_FAST_GetTimeEnd();
		
		int idptfmheave = DLL_FAST_GetParameterId("ptfmheave");
		double avg = 0;
		for (int i = 0; i < DLL_FAST_GetLen(); ++i)
			avg += DLL_FAST_GetData(i, idptfmheave);
		Cout() << "\nptfmheave_avg = " << avg/DLL_FAST_GetLen();
		Cout() << "\nptfmheave_avg = " << DLL_FAST_GetAvg("ptfmheave");
		
		DLLFunction(dll, int, 	   		DLL_FAST_LoadFile, (const char *file));
		DLLFunction(dll, int, 	   		DLL_FAST_SaveFile, (const char *file));
		DLLFunction(dll, int, 	   		DLL_FAST_SetVar, (const char *name, const char *paragraph, const char *value));
		DLLFunction(dll, const char *, 	DLL_FAST_GetVar, (const char *name, const char *paragraph));				

		Cout() << "\n\nLoading InflowWind .dat file";
		String datfile = AppendFileNameX(bemFolder, "examples/fast.out/InflowWind.dat");
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
		DLL_FAST_SaveFile(AppendFileNameX(GetDesktopFolder(), "InflowWind_test.dat"));
	#endif
#endif
	
	} catch (Exc err) {
		Cout() << "\n" << Format(t_("Problem found: %s"), err);
		SetExitCode(-1);
	}
#ifdef flagDEBUG
	Cout() << "\nPress a key to end";
	ReadStdIn();
#endif
}

#endif

#ifndef flagBEMR_TEST_DLL

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
	
	return str = CleanCFromDeclaration(String(DLLexport, DLLexport_length));
}

const char *DLL_strPythonDeclaration() noexcept {
	static String str;
	
	return str = GetPythonDeclaration("BEMRosetta", String(DLLexport, DLLexport_length));	
}

void DLL_ListFunctions() noexcept {
	Cout() << DLL_strListFunctions();	
}

int DLL_FAST_Load(const char *filename) noexcept {
	try {
		return DLL_Fastout().Load(filename);
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
		return DLL_Fastout().size();
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetLen()";
		return Null;
	}
}

double DLL_FAST_GetTimeInit() noexcept {
	try {
		return DLL_Fastout().GetTimeInit();
	} catch (...) {
		CoutX() << "Unknown error in DLL_FAST_GetTimeInit()";
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
		CoutX() << "DLL_FAST_GetData() idtime < 0";
		return Null;
	}
	if (idtime >= DLL_Fastout().size()) {
		CoutX() << "DLL_FAST_GetData() idtime >= time";
		return Null;
	}
	if (idparam < 0) {
		CoutX() << "DLL_FAST_GetData() idparam < 0";
		return Null;
	}
	if (idparam >= DLL_Fastout().GetParameterCount()) {
		CoutX() << "DLL_FAST_GetData() idparam >= num_params";
		return Null;
	}
		
	return DLL_Fastout().GetVal(idtime, idparam);
}

double DLL_FAST_GetAvg(const char *param) noexcept {
	try {
		const UVector<double> &data = DLL_Fastout().GetUVector(param);
		return Eigen::Map<const Eigen::VectorXd>(data, data.size()).mean();
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

int DemoVectorC_Py(double **v, int *num) noexcept {
	static UVector<double> data;
	
	*num = 321;
	data.SetCount(*num);
	for (int i = 0; i < *num; ++i)
		data[i] = i/2.;
	*v = data;
	
	return 1;
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



