// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2026, the BEMRosetta author and contributors
#if !defined(flagGUI)

#include "BEMRosetta.h"
#include <SysInfo/Crash.h>
#include "libbemrosetta.h"

#ifdef PLATFORM_WIN32

#if defined(flagBEMR_TEST_DLL) || defined(flagBEMR_TEST_DLL_INTERNAL)	


String GetPythonDeclaration(const String &name, const String &prefix, const String &include);
	
String BMR_strPythonDeclaration(String export_h) {
	return GetPythonDeclaration("BEMRosetta", "BMR", export_h);	
}

String BMR_strListFunctions(String export_h) {
	return CleanCFromDeclaration(export_h);
}

#if defined(flagBEMR_TEST_DLL_INTERNAL) || defined(flagBEMR_TEST_DLL)
void BEM_Throw();
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
		
		if (command.GetCount() < 1) 
			throw Exc("Please include in command line binary and BEMRosetta folders");
		
		String unittestFolder = command[0];
		String binFolder = AFX(unittestFolder, ".\\.test");
		String bemFolder = AFX(unittestFolder, "..");
		String installFolder = AFX(bemFolder, "install");
		String export_h = AFX(bemFolder, "BEMRosetta_cl", "libbemrosetta.h"); 
#endif
#ifdef flagBEMR_TEST_DLL
		FileDelete(AFX(binFolder, "libbemrosetta.exp"));
		//FileDelete(AFX(binFolder, "libbemrosetta.lib"));
		
		Dl dll;		
		if (!dll.Load(AFX(binFolder, "libbemrosetta.dll")))
			throw Exc("DLL not found");
		
		DLLFunction(dll, void, 		   BMR_Init, ());
		DLLFunction(dll, const char *, BMR_Version, ());
		DLLFunction(dll, const char *, BMR_GetLastError, ());
		
		DLLFunction(dll, void, BMR_Mesh_Load, (const char *file));
    	DLLFunction(dll, void, BMR_Mesh_Translate, (double x, double y, double z));
		DLLFunction(dll, void, BMR_Mesh_C0_Set, (double x, double y, double z));	
		DLLFunction(dll, void, BMR_Mesh_Cg_Set, (double x, double y, double z));
    	DLLFunction(dll, void, BMR_Mesh_Inertia_Set, (const double *data, const int dim[2]));                     
    	DLLFunction(dll, void, BMR_Mesh_LinearDamping_Set, (const double *data, const int dim[2]));
		DLLFunction(dll, void, BMR_Mesh_MooringStiffness_Set, (const double *data, const int dim[2]));
    	DLLFunction(dll, void, BMR_Bem_depth_Set, (double h));
		DLLFunction(dll, void, BMR_Bem_g_Set, (double g));
		DLLFunction(dll, void, BMR_Bem_rho_Set, (double rho));
		DLLFunction(dll, void, BMR_Bem_w_Set, (const double *w, int dim));
		DLLFunction(dll, void, BMR_Bem_headings_Set, (const double *head, int dim));
    	DLLFunction(dll, void, BMR_Bem_Body_LoadMeshFromMesh, (int id));
    	DLLFunction(dll, void, BMR_Bem_SaveCase, (const char *folder, const char *solver, bool x0z, bool y0z, bool irregular, bool autoIrregular, const char *qtfType, bool autoQTF, bool bin, int numCases, int numThreads, bool withPotentials, bool withMesh));
    
    	DLLFunction(dll, void, BMR_Bem_Body_LoadMesh, (const char *file));
		DLLFunction(dll, void, BMR_Bem_Body_Cg_Set, (double x, double y, double z));
		DLLFunction(dll, void, BMR_Bem_Body_Inertia_Set, (const double *data, const int dim[2]));
		DLLFunction(dll, int,  BMR_Bem_Load, (const char *file));
		DLLFunction(dll, void, BMR_Bem_Save, (const char *file));
			
    	DLLFunction(dll, void, BMR_Mesh_Volume_Get, (double *volx, double *voly, double *volz));
    	DLLFunction(dll, void, BMR_Mesh_UnderwaterVolume_Get, (double *volx, double *voly, double *volz));
    	DLLFunction(dll, void, BMR_Mesh_HydrostaticStiffness_Get, (double **data, int dim[2]));
    	
		DLLFunction(dll, int, 		   BMR_FAST_Load, (const char *));
		DLLFunction(dll, const char *, BMR_FAST_GetParameterName, (int));
		DLLFunction(dll, const char *, BMR_FAST_GetUnitName, (int));
		DLLFunction(dll, int, 		   BMR_FAST_GetParameterId, (const char *name));
		DLLFunction(dll, int, 		   BMR_FAST_GetParameterCount, ());
		DLLFunction(dll, int, 		   BMR_FAST_GetLen, ());		
		DLLFunction(dll, double,	   BMR_FAST_GetTimeStart, ());
		DLLFunction(dll, double,	   BMR_FAST_GetTimeEnd, ());
		DLLFunction(dll, double, 	   BMR_FAST_GetTime, (int idtime));		
		DLLFunction(dll, double, 	   BMR_FAST_GetData, (int idtime, int idparam));
		DLLFunction(dll, double, 	   BMR_FAST_GetAvg,   (int idparam, int idbegin, int idend));
		DLLFunction(dll, double, 	   BMR_FAST_GetArray, (int idparam, int idbegin, int idend, double **, int *));
		
		DLLFunction(dll, int, 	   	   BMR_FAST_LoadFile, (const char *file));
		DLLFunction(dll, int, 	   	   BMR_FAST_SaveFile, (const char *file));
		DLLFunction(dll, int, 	   	   BMR_FAST_SetVar, (const char *name, const char *paragraph, const char *value));
		DLLFunction(dll, const char *, BMR_FAST_GetVar, (const char *name, const char *paragraph));
#endif
		BMR_Init();	
		Cout() << "\nVersion: " << BMR_Version();
		Cout() << "\n\nDLL functions list:\n";
		String strList = BMR_strListFunctions(LoadFile(export_h));
		Cout() << strList;
#ifdef flagBEMR_TEST_DLL
		strList = "// BEMRosetta DLL functions list\n\n" + strList;	
		if (!SaveFile(AFX(binFolder, "libbemrosetta.txt"), strList))
			throw Exc(t_("Impossible to save DLL functions list file"));
#endif

		Cout() << "\n\nPython declarations:\n";
		String strPy = BMR_strPythonDeclaration(LoadFile(export_h));
		Cout() << strPy;
		
#if defined(flagBEMR_TEST_DLL) || defined(flagBEMR_TEST_DLL_INTERNAL)
		if (!SaveFile(AFX(binFolder, "libbemrosetta.py"), strPy))
			throw Exc(t_("Impossible to save Python declarations file"));
		if (!FileCopy(export_h, AFX(binFolder, "libbemrosetta.h")))
			throw Exc(t_("Impossible to save C declarations file"));		
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

		{
			printf("\nGenerating new case");
			const char *meshFile = "../examples/capytaine/Orca/Body_1.gdf";	
			BMR_Bem_depth_Set(50);
			BMR_Bem_g_Set(9.81);
			BMR_Bem_rho_Set(1025);
			double w[] = {0.1, 0.5, 1, 1.5, 2};
			BMR_Bem_w_Set(w, sizeof(w)/sizeof(double));
			double head[] = {0, 45, 90};
			BMR_Bem_headings_Set(head, sizeof(head)/sizeof(double));
			BMR_Bem_Body_LoadMesh(meshFile);
			BMR_Bem_Body_Cg_Set(0, 0, 0);
			double M[] = {9.8E5,     0,     0,   0,   0,   0,
						      0, 9.8E5,     0,   0,   0,   0,
						      0,     0, 9.8E5,   0,   0,   0,
						      0,     0,     0, 1E7,   0,   0,
						      0,     0,     0,   0, 1E7,   0,
						      0,     0,     0,   0,   0, 1E7};
			int dim[] = {6, 6};
			BMR_Bem_Body_Inertia_Set(M, dim);
			printf("\nSaving it in Capytaine format");
			BMR_Bem_SaveCase("../unittest/.test/Capy", "Capytaine .py", false, false, true, true, "No", false, false, 1, 4, false, false);
			printf("\nRunning it in Capytaine");
			system("cd /d ..\\unittest\\.test\\Capy && capytaine.bat");
			printf("\nLoading the results in .nc format");
			BMR_Bem_Load("..\\unittest\\.test\\Capy\\capytaine.nc");
			printf("\nSaving the results in .h5 format");
			BMR_Bem_Save("..\\unittest\\.test\\Capy\\capytaine.h5");
		}
	  
		BMR_Mesh_Load("../examples/hydrostar/Mesh/Ship.hst");
	    double volx, voly, volz;
	    BMR_Mesh_Volume_Get(&volx, &voly, &volz);
	    Cout() << F("\nVolume            x: %f, y: %f, z: %f", volx, voly, volz);
	    BMR_Mesh_UnderwaterVolume_Get(&volx, &voly, &volz);
	    Cout() << F("\nUnderwater volume x: %f, y: %f, z: %f", volx, voly, volz);
		double *C;
		int dim[2];
		BMR_Mesh_HydrostaticStiffness_Get(&C, dim);
		Cout() << "\nStiffness matrix   :\n";
		int ic = 0;
		for (int r = 0; r < dim[0]; ++r) {
			for (int c = 0; c < dim[1]; ++c)
				Cout() << C[ic++] << ",\t";
			Cout() << "\n";
		}
		Cout() << "\n\nLoading FAST .out file";
		String outfile = AFX(bemFolder, "examples/fast.out/demo.outb");
		if (!BMR_FAST_Load(outfile))
			throw Exc(F("Impossible to open file %s", outfile));
		
		Cout() << "\nFAST .out parameters:";
		for (int i = 0; i < BMR_FAST_GetParameterCount(); ++i)
			Cout() << F(" %s[%s]", BMR_FAST_GetParameterName(i), BMR_FAST_GetUnitName(i));
		
		Cout() << "\nSimulation begins at " << BMR_FAST_GetTimeStart() << " and ends at " << BMR_FAST_GetTimeEnd();
		
		int idptfmheave = BMR_FAST_GetParameterId("ptfmHeave");
		double *v;
		int num;
		BMR_FAST_GetArray(idptfmheave, -1, -1, &v, &num);
		Cout() << "\nRead " << num << " heave values"; 

		double avg = 0, avgv = 0;
		for (int i = 0; i < BMR_FAST_GetLen(); ++i) {
			avg += BMR_FAST_GetData(i, idptfmheave);
			avgv += v[i];
		}
		Cout() << "\nptfmheave_avg = " << avg/BMR_FAST_GetLen();
		Cout() << "\nptfmheave_avg = " << avgv/BMR_FAST_GetLen();
		Cout() << "\nptfmheave_avg = " << BMR_FAST_GetAvg(idptfmheave, -1, -1);
		
		Cout() << "\n\nLoading InflowWind .dat file";
		String datfile = AFX(bemFolder, "examples/fast.out/InflowWind.dat");
		if (!BMR_FAST_LoadFile(datfile))
			throw Exc(F("Impossible to open file %s", outfile));	
	
		String str;
		str = BMR_FAST_GetVar("WindVziList", "");				
		Cout() << "\nWindVziList: " << str;
		VERIFY(str == "119");
		str = BMR_FAST_GetVar("nx", "");				
		Cout() << "\nnx: " << str;
		VERIFY(str == "64");
		str = BMR_FAST_GetVar("Filename", "");				
		Cout() << "\nFilename: " << str;
		VERIFY(str == "unifWind.hh");
		str = BMR_FAST_GetVar("Filename", "================== Parameters for Uniform wind file");				
		Cout() << "\nFilename: " << str;
		VERIFY(str == "unifWind.hh");
		str = BMR_FAST_GetVar("Filename", "================== Parameters for Binary TurbSim");				
		Cout() << "\nFilename: " << str;
		VERIFY(str == "TurbSim.bts");
		BMR_FAST_SetVar("nx", "", "23");				
		str = BMR_FAST_GetVar("nx", "");				
		Cout() << "\nNew nx: " << str;
		VERIFY(str == "23");
		BMR_FAST_SetVar("Filename", "================== Parameters for Binary TurbSim", "\"New file\"");				
		str = BMR_FAST_GetVar("Filename", "================== Parameters for Binary TurbSim");				
		Cout() << "\nNew Filename: " << str;
		VERIFY(str == "New file");

	#ifdef flagDEBUG	
		BMR_FAST_SaveFile(AFX(unittestFolder, "InflowWind_test.dat"));
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
		Cout() << "\n" << F(t_("Problem found: %s"), errorStr);
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

BMR_Data &BMR();

CONSOLE_APP_MAIN {
	const UVector<String>& command = CommandLine();
	
	BMR().Status = PrintStatus;
	if (!BMR().ConsoleMain(command, false))
		SetExitCode(1);
	
#ifdef flagDEBUG
	Cout() << "\nPress enter to end";
	ReadStdIn();
#endif
}
#endif



