#ifndef flagGUI

#include "BEMRosetta.h"


CONSOLE_APP_MAIN {
	String str = "BEMRosetta Copyright (c) Inaki Zabala\nHydrodynamic coefficients converter for Boundary Element Method solver formats\nVersion beta [Build Info]";
	Hydro::SetBuildInfo(str);
	Cout() << str;
	
	const Vector<String>& command = CommandLine();
	
	String errorStr;
	
	String fileName;
	try {
		Nemoh nem;
		nem.Load(fileName);
		
		Wamit wam;
		wam.Load(fileName);
		
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
	if (!errorStr.IsEmpty())
		Cout() << "\nError: " << errorStr;
	
	
	Cout() << "\nPress enter to end\n";
	ReadStdIn();	
}

#endif