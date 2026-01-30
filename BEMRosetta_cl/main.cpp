// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2026, the BEMRosetta author and contributors
#if !defined(flagGUI)

#include "BEMRosetta.h"
#include <SysInfo/Crash.h>
#include "export.h"

#ifdef PLATFORM_WIN32

#if defined(flagBEMR_TEST_DLL) || defined(flagBEMR_TEST_DLL_INTERNAL)	


String GetPythonDeclaration(const String &name, const String &include);
	
String DLL_strPythonDeclaration(String export_h) {
	return GetPythonDeclaration("BEMRosetta", export_h);	
}

String DLL_strListFunctions(String export_h) {
	return CleanCFromDeclaration(export_h);
}


CONSOLE_APP_MAIN
{	
#if defined(flagDEBUG) && defined(PLATFORM_WIN32) && !defined(flagBEMR_TEST_DLL)
	GetCrashHandler().Enable();
#endif
	String errorStr;
	try {
#if defined(flagBEMR_TEST_DLL) || defined(flagBEMR_TEST_DLL_INTERNAL)
		const UVector<String>& command = CommandLine();
		
		if (command.GetCount() < 1) 
			throw Exc("Please include in command line binary and BEMRosetta folders");
		
		String unittestFolder = command[0];
		String binFolder = AFX(unittestFolder, ".\\.test");
		String bemFolder = AFX(unittestFolder, "..");
		String installFolder = AFX(bemFolder, "install");
		String export_h = AFX(bemFolder, "BEMRosetta_cl", "export.h"); 
#endif
#ifdef flagBEMR_TEST_DLL
		FileDelete(AFX(binFolder, "libbemrosetta.exp"));
		FileDelete(AFX(binFolder, "libbemrosetta.lib"));
		
		Dl dll;		
		if (!dll.Load(AFX(binFolder, "libbemrosetta.dll")))
			throw Exc("DLL not found");
		
		DLLFunction(dll, void, 		   DLL_Init, ());
		DLLFunction(dll, const char *, DLL_Version, ());
	
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
		
		DLLFunction(dll, int, 	   	   DLL_FAST_LoadFile, (const char *file));
		DLLFunction(dll, int, 	   	   DLL_FAST_SaveFile, (const char *file));
		DLLFunction(dll, int, 	   	   DLL_FAST_SetVar, (const char *name, const char *paragraph, const char *value));
		DLLFunction(dll, const char *, DLL_FAST_GetVar, (const char *name, const char *paragraph));
#endif

		UVector<double> dat = {1, 2, 3};
		//double res = DLL_DemoVectorPyC(dat, 3);
		
		DLL_Init();	
		Cout() << "\nVersion: " << DLL_Version();
		Cout() << "\n\nDLL functions list:\n";
		String strList = DLL_strListFunctions(LoadFile(export_h));
		Cout() << strList;
#ifdef flagBEMR_TEST_DLL
		strList = "// BEMRosetta DLL functions list\n\n" + strList;	
		if (!SaveFile(AFX(binFolder, "libbemrosetta.txt"), strList))
			throw Exc(t_("Impossible to save DLL functions list file"));
#endif

		Cout() << "\n\nPython declarations:\n";
		String strPy = DLL_strPythonDeclaration(LoadFile(export_h));
		Cout() << strPy;
		
#if defined(flagBEMR_TEST_DLL) || defined(flagBEMR_TEST_DLL_INTERNAL)
		if (!SaveFile(AFX(binFolder, "libbemrosetta.py"), strPy))
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
		DLL_Mesh_Input("../examples/hydrostar/Mesh/Ship.hst");
	    double volx, voly, volz;
	    DLL_Mesh_GetVolume(&volx, &voly, &volz);
	    Cout() << Format("\nVolume            x: %f, y: %f, z: %f", volx, voly, volz);
	    DLL_Mesh_GetUnderwaterVolume(&volx, &voly, &volz);
	    Cout() << Format("\nUnderwater volume x: %f, y: %f, z: %f", volx, voly, volz);
		double *C;
		int dim[2];
		DLL_Mesh_GetHydrostaticStiffness(&C, dim);
		Cout() << "\nStiffness matrix   :\n";
		int ic = 0;
		for (int r = 0; r < dim[0]; ++r) {
			for (int c = 0; c < dim[1]; ++c)
				Cout() << C[ic++] << ",\t";
			Cout() << "\n";
		}
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
		VERIFY(str == "unifWind.hh");
		str = DLL_FAST_GetVar("Filename", "================== Parameters for Uniform wind file");				
		Cout() << "\nFilename: " << str;
		VERIFY(str == "unifWind.hh");
		str = DLL_FAST_GetVar("Filename", "================== Parameters for Binary TurbSim");				
		Cout() << "\nFilename: " << str;
		VERIFY(str == "TurbSim.bts");
		DLL_FAST_SetVar("nx", "", "23");				
		str = DLL_FAST_GetVar("nx", "");				
		Cout() << "\nNew nx: " << str;
		VERIFY(str == "23");
		DLL_FAST_SetVar("Filename", "================== Parameters for Binary TurbSim", "\"New file\"");				
		str = DLL_FAST_GetVar("Filename", "================== Parameters for Binary TurbSim");				
		Cout() << "\nNew Filename: " << str;
		VERIFY(str == "New file");

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

#endif

#endif

#if defined(flagBEMR_CL)

DLL_Data &DLL();

CONSOLE_APP_MAIN {
	const UVector<String>& command = CommandLine();
	
	DLL().Status = PrintStatus;
	if (!DLL().ConsoleMain(command, false))
		SetExitCode(1);
	
#ifdef flagDEBUG
	Cout() << "\nPress enter to end";
	ReadStdIn();
#endif
}
#endif



