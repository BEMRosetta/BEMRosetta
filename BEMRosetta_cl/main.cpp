#ifndef flagGUI

#include "BEMRosetta.h"


CONSOLE_APP_MAIN {
	String errorStr;
	
	String fileName;
	try {
		Nemoh nem;
		nem.Load(fileName);
		nem.Report();
		Wamit wam;
		wam.Load(fileName);
		wam.Report();
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