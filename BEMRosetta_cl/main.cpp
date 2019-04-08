#ifndef flagGUI

#include "BEMRosetta.h"


void ShowHelp() {
	Cout() << 	"\nUsage: bemrosetta_cl [option] infile outfolder"
				"\nOptions:"
				"\n-h --help     -- Print options"
				"\n-c --convert  -- Convert from input file to output file";
}

CONSOLE_APP_MAIN {
	String str = "BEMRosetta Copyright (c) 2019 Iñaki Zabala\nHydrodynamic coefficients converter for Boundary Element Method solver formats\nVersion beta [Build Info]";
	Hydro::SetBuildInfo(str);
	Cout() << str;
	
	const Vector<String>& command = CommandLine();
	
	String errorStr;
	try {
		if (command.IsEmpty())
			ShowHelp();
		else {
			for (int i = 0; i < command.GetCount(); ++i) {
				if (command[i] == "-h" || command[i] == "--help") {
					ShowHelp();
					break;
				} else if (command[i] == "-c" || command[i] == "--convert") {
					if (command.GetCount() != 3)
						throw Exc("Wrong number of arguments");
					/*Nemoh nem;
					nem.Load(fileName);
					
					Wamit wam;
					wam.Load(fileName);*/
				} else 
					throw Exc(Format("Unknown argument '%s'", command[i]));
			}
		}
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
		Cerr() << "\nError: " << errorStr;
	
	Cout() << "\nPress enter to end\n";
	ReadStdIn();	
}

#endif