#ifndef flagGUI

#include "BEMRosetta.h"


CONSOLE_APP_MAIN {
	String errorStr;
	
	try {
		Nemoh nem;
		nem.Load("C:/Users/0203853/Desktop/Wamit_Nemoh/SAPA/Nemoh/Nemoh.cal");
		nem.Report();
		nem.Load("C:/Users/0203853/Desktop/WEC-Sim-3.1/tutorials/BEMIO/NEMOH/RM3/Nemoh.cal");
		nem.Report();
		nem.Load("C:/Users/0203853/Desktop/Nemoh 1 Prisma/Nemoh.cal");
		nem.Report();
		nem.Load("C:/Users/0203853/Desktop/Nemoh 2 Prismas/Nemoh.cal");
		nem.Report();
		nem.Load("C:/Users/0203853/Desktop/Nemoh 3 Prismas/Nemoh.cal");
		nem.Report();
		Wamit wam;
		wam.Load("C:/Users/0203853/Desktop/WEC-Sim-3.1/tutorials/BEMIO/WAMIT/Cubes/cubes.out");
		wam.Report();
		wam.Load("C:/Users/0203853/Desktop/Wamit_Nemoh/SAPA/Wamit/sapa.out");
		wam.Report();
		wam.Load("C:/Users/0203853/Desktop/Wamit_Nemoh/OSWC/Wamit/oswc.out");
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