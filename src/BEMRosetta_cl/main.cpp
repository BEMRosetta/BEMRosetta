#if !defined(flagGUI)

#include "BEMRosetta.h"

#if defined(flagTEST_DLL)


CONSOLE_APP_MAIN
{
	try {
		Dl dll;		
		if (!dll.Load("libbemrosetta.dll"))
			throw Exc("Dll not found");

		DLLFunction(dll, const char *, DLL_Version, ());
		DLLFunction(dll, void, 		   DLL_ListFunctions, ());
		DLLFunction(dll, int, 		   DLL_FAST_Load, (const char *));
		DLLFunction(dll, const char *, DLL_FAST_GetParameterName, (int));
		DLLFunction(dll, const char *, DLL_FAST_GetUnitName, (int));
		DLLFunction(dll, int, 		   DLL_FAST_GetParameterCount, ());
		DLLFunction(dll, int, 		   DLL_FAST_GetLen, ());		
		DLLFunction(dll, double,	   DLL_FAST_GetTimeInit, ());
		DLLFunction(dll, double,	   DLL_FAST_GetTimeEnd, ());
		DLLFunction(dll, double *, 	   DLL_FAST_GetDataId, (int, int *));		
		DLLFunction(dll, double *, 	   DLL_FAST_GetData, (const char *, int *));
		DLLFunction(dll, double, 	   DLL_FAST_GetAvg, (const char *));

		Cout() << "\nVersion: " << DLL_Version();
		Cout() << "\nDLL libbemrosetta functions list:\n";
		DLL_ListFunctions();
		
		Cout() << "\nLoading FAST .out file";
		if (!DLL_FAST_Load("D:/Desarrollo/BEMRosetta/Material base/FAST/aaaaaaaaa.outb"))
			throw Exc("Impossible to open file");
		
		Cout() << "\nFAST .out parameters:";
		for (int i = 0; i < DLL_FAST_GetParameterCount(); ++i)
			Cout() << Format(" %s[%s]", DLL_FAST_GetParameterName(i), DLL_FAST_GetUnitName(i));
		
		Cout() << "\nSimulation begins at " << DLL_FAST_GetTimeInit() << " and ends at " << DLL_FAST_GetTimeEnd();
		
		int num;
		double *data = DLL_FAST_GetData("ptfmheave", &num);
		Cout() << "\nptfmheave has " << num << " data. ptfmheave[0] = " << data[0];
		Cout() << "\nptfmheave_avg = " << DLL_FAST_GetAvg("ptfmheave");
	
	} catch (Exc err) {
		Cout() << "\n" << Format(t_("Problem found: %s"), err);
	}
	
	Cout() << "\nPress a key to end";
	ReadStdIn();
}

#elif defined (flagDLL)

#include "FastOut.h"
#include "export.h"
#include "export.brc"

const char *DLL_Version() noexcept {
	static String version;
	version << __DATE__ << ", " << __TIME__;
	return version;	
}

void DLL_ListFunctions() noexcept {
	String str(DLLexport, DLLexport_length);
	
	str.Replace("	__declspec(dllexport) ", "");
	str.Replace("extern \"C\" {", "");
	str.Replace("};", "");
	str.Replace("\r\n\r\n", "\r\n");
	str.Replace(";", "");
	str.Replace("noexcept", "");
	str.Replace("  ", "");
	str.Replace("\t", "");
	
	Cout() << Trim(str);	
}

FastOut &DLL_Fastout() {
	static FastOut fast;
	return fast;
}

int DLL_FAST_Load(const char *filename) noexcept {
	try {
		return DLL_Fastout().Load(filename);
	} catch (...) {
		return 0;
	}
}

const char *DLL_FAST_GetParameterName(int id) noexcept {
	static String ret;
	return ret = DLL_Fastout().GetParameter(id);
}

const char *DLL_FAST_GetUnitName(int id) noexcept {
	static String ret;
	return ret = DLL_Fastout().GetUnit(id);
}
	
int DLL_FAST_GetParameterCount() noexcept 	{return DLL_Fastout().GetParameterCount();}
int DLL_FAST_GetLen() noexcept 				{return DLL_Fastout().size();}
double DLL_FAST_GetTimeInit() noexcept		{return DLL_Fastout().GetTimeInit();}
double DLL_FAST_GetTimeEnd() noexcept		{return DLL_Fastout().GetTimeEnd();}

double *DLL_FAST_GetDataId(int id, int *num) noexcept {
	static Vector<double> data;
	data = clone(DLL_Fastout().GetVal(id));
	*num = data.size();
	return data;
}

double *DLL_FAST_GetData(const char *param, int *num) noexcept {
	static Vector<double> data;
	try {
		data = clone(DLL_Fastout().GetVal(param));
	} catch (...) {
	}
	*num = data.size();
	return data;
}

double DLL_FAST_GetAvg(const char *param) noexcept {
	try {
		const Vector<double> &data = DLL_Fastout().GetVal(param);
		return Eigen::Map<const Eigen::VectorXd>(data, data.size()).mean();
	} catch (...) {
	}
	return Null;
}

int DLL_IsNull(double val) noexcept {return IsNull(val);}


#else

CONSOLE_APP_MAIN {
	const Vector<String>& command = CommandLine();
	
	ConsoleMain(command, false);
}

#endif
#endif