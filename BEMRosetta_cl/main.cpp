#if !defined(flagGUI)

#include "BEMRosetta.h"

#ifdef PLATFORM_WIN32

#if defined(flagTEST_DLL)

CONSOLE_APP_MAIN
{
	try {
		const Vector<String>& command = CommandLine();
		
		if (command.GetCount() < 2) 
			throw Exc("Please include in command line binary and BEMRosetta folders");
		
		String binFolder = command[0];
		String bemFolder = command[1];
		 
		Dl dll;		
		if (!dll.Load(AppendFileName(binFolder, "libbemrosetta.dll")))
			throw Exc("Dll not found");

		DLLFunction(dll, const char *, DLL_Version, ());
		DLLFunction(dll, void, 		   DLL_ListFunctions, ());
		DLLFunction(dll, const char *, DLL_strListFunctions, ());
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

		Cout() << "\nVersion: " << DLL_Version();
		Cout() << "\nDLL functions list:\n";
		String strList = DLL_strListFunctions();
		Cout() << strList;
		
		strList = "// DLL functions list\n\n" + strList;
		SaveFile(AppendFileNameX(binFolder, "libbemrosetta.txt"), strList);
		
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
		DLL_FAST_SaveFile(AppendFileNameX(GetDesktopFolder(), "test.dat"));
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

#else

#include "FastOut.h"
#include "export.h"
#include "export.brc"

FastOut &DLL_Fastout() {
	static FastOut fast;
	return fast;
}

extern "C" {
	
const char *DLL_Version() noexcept {
	static String version;
	version << __DATE__ << ", " << __TIME__;
	return version;	
}

const char *DLL_strListFunctions() noexcept {
	static String str;
	
	str = String(DLLexport, DLLexport_length);
	
	str.Replace("	__declspec(dllexport) ", "");
	str.Replace("extern \"C\" {", "");
	str.Replace("};", "");
	str.Replace("\r\n\r\n", "\r\n");
	str.Replace(";", "");
	str.Replace("noexcept", "");
	str.Replace("  ", "");
	str.Replace("\t", "");
	
	return str = Trim(str);	
}

void DLL_ListFunctions() noexcept {
	Cout() << DLL_strListFunctions();	
}

int DLL_FAST_Load(const char *filename) noexcept {
	try {
		return DLL_Fastout().Load(filename);
	} catch (...) {
		Cout() << "Unknown error in DLL_FAST_Load()";
		return 0;
	}
}

const char *DLL_FAST_GetParameterName(int id) noexcept {
	static String ret;
	try {
		return ret = DLL_Fastout().GetParameter(id);
	} catch (...) {
		Cout() << "Unknown error in DLL_FAST_GetParameterName()";
		return ret = "Error";
	}
}

const char *DLL_FAST_GetUnitName(int id) noexcept {
	static String ret;
	try {
		return ret = DLL_Fastout().GetUnit(id);
	} catch (...) {
		Cout() << "Unknown error in DLL_FAST_GetUnitName()";
		return ret = "Error";
	}
}

int DLL_FAST_GetParameterId(const char *name) noexcept {
	try {
		Vector<int> p = DLL_Fastout().FindParameterMatch(name);
		if (p.IsEmpty())
			return -1;
		else
			return p[0];
	} catch (...) {
		Cout() << "Unknown error in DLL_FAST_GetParameterCount()";
		return Null;
	}
}

int DLL_FAST_GetParameterCount() noexcept {
	try {
		return DLL_Fastout().GetParameterCount();
	} catch (...) {
		Cout() << "Unknown error in DLL_FAST_GetParameterCount()";
		return Null;
	}
}

int DLL_FAST_GetLen() noexcept {
	try {
		return DLL_Fastout().size();
	} catch (...) {
		Cout() << "Unknown error in DLL_FAST_GetLen()";
		return Null;
	}
}

double DLL_FAST_GetTimeInit() noexcept {
	try {
		return DLL_Fastout().GetTimeInit();
	} catch (...) {
		Cout() << "Unknown error in DLL_FAST_GetTimeInit()";
		return Null;
	}
}

double DLL_FAST_GetTimeEnd() noexcept {
	try {
		return DLL_Fastout().GetTimeEnd();
	} catch (...) {
		Cout() << "Unknown error in DLL_FAST_GetTimeEnd()";
		return Null;
	}
}

double DLL_FAST_GetTime(int idtime) noexcept {
	return DLL_FAST_GetData(idtime, 0);
}

double DLL_FAST_GetData(int idtime, int idparam) noexcept {
	if (idtime < 0) {
		Cout() << "DLL_FAST_GetData() idtime < 0";
		return Null;
	}
	if (idtime >= DLL_Fastout().size()) {
		Cout() << "DLL_FAST_GetData() idtime >= time";
		return Null;
	}
	if (idparam < 0) {
		Cout() << "DLL_FAST_GetData() idparam < 0";
		return Null;
	}
	if (idparam >= DLL_Fastout().GetParameterCount()) {
		Cout() << "DLL_FAST_GetData() idparam >= num_params";
		return Null;
	}
		
	return DLL_Fastout().GetVal(idtime, idparam);
}

double DLL_FAST_GetAvg(const char *param) noexcept {
	try {
		const Vector<double> &data = DLL_Fastout().GetVal(param);
		return Eigen::Map<const Eigen::VectorXd>(data, data.size()).mean();
	} catch (...) {
		Cout() << "Unknown error in DLL_FAST_GetAvg()";
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
		Cout() << "\n" << "Error: " << e;
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
		Cout() << "\n" << "Error: " << e;
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
		Cout() << "\n" << "Error: " << e;
		SetConsoleColor(CONSOLE_COLOR::PREVIOUS);
		return ret = "";
	}
	if (IsVoid(ret))
		return "";
	return ret;
}

}

#endif

#endif


#if !defined(flagTEST_DLL) && !defined(flagDLL)

CONSOLE_APP_MAIN {
	const Vector<String>& command = CommandLine();
	
	ConsoleMain(command, false);
}

#endif

#endif