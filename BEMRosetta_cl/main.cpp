#ifndef flagGUI

#include "BEMRosetta.h"


void ShowHelp(BEMData &md) {
	Cout() << "\n" << t_("Usage: bemrosetta_cl [options] [-i infile]... [-e outfile]");
	Cout() << "\n";
	Cout() << "\n" << t_("Options:");
	Cout() << "\n" << t_("-h  --help     -- print options");
	Cout() << "\n" << t_("-p  --params   -- set physical parameters:");
	Cout() << "\n" << t_("                 parameter description   units  default value");
	Cout() << "\n" << t_("                    g      gravity       [m/s2]    ") << md.g;
	Cout() << "\n" << t_("                    length length scale  []        ") << md.length;
	Cout() << "\n" << t_("                    rho    water density [Kg/m3]   ") << md.rho;
	Cout() << "\n" << t_("                    depth  water depth   [m]       ") << md.depth;
	Cout() << "\n" << t_("-i  --input    -- load model");
	Cout() << "\n" << t_("-e  --export   -- export from input file to output file");
	Cout() << "\n" << t_("-c  --compare  -- compare input files");
	Cout() << "\n" << t_("-r  --report   -- output last loaded model data");
	Cout() << "\n" << t_("-cl --clear    -- clear loaded model");
	Cout() << "\n";
	Cout() << "\n" << t_("Actions");
	Cout() << "\n" << t_("- are done in sequence: if a physical parameter is changed after export, saved files will not include the change");
	Cout() << "\n" << t_("- can be repeated as desired");
}

void CheckNumArgs(const Vector<String>& command, int i, String param) {
	if (i >= command.GetCount())	
		throw Exc(Format(t_("Missing parameters when reading '%s'"), param));
}

void WamitAdditionalData(BEMData &md, HydroClass &data) {
	if (IsNull(data.hd().g)) 
		data.hd().g = md.g;
	
	if (IsNull(data.hd().len)) 
		data.hd().len = md.length;
	
	if (IsNull(data.hd().rho)) 
		data.hd().rho = md.rho;
	
	if (IsNull(data.hd().h)) 
		data.hd().h = md.depth;
}

CONSOLE_APP_MAIN {
	String str = t_("BEMRosetta Copyright (c) 2019 Iñaki Zabala\nHydrodynamic coefficients converter for Boundary Element Method solver formats\nVersion beta BUILDINFO");
	Hydro::SetBuildInfo(str);
	Cout() << str;
	
	BEMData md;
	
	if (!md.LoadSerializeJson())
		Cout() << "\n" << t_("BEM configuration data are not loaded. Defaults are set");

	const Vector<String>& command = CommandLine();
	
	String errorStr;
	try {
		if (command.IsEmpty()) {
			Cout() << "\n" << t_("Command argument list is empty");
			ShowHelp(md);
		} else {
			for (int i = 0; i < command.GetCount(); i++) {
				if (command[i] == "-h" || command[i] == "--help") {
					ShowHelp(md);
					break;
				} else if (command[i] == "-i" || command[i] == "--input") {
					i++;
					CheckNumArgs(command, i, "--input");
					
					String file = command[i];
					if (!FileExists(file)) 
						throw Exc(Format(t_("File '%s' not found"), file)); 
					
					md.Load(file, STDBACK(WamitAdditionalData));
					Cout() << "\n" << Format(t_("File '%s' loaded"), file);
				} else if (command[i] == "-r" || command[i] == "--report") {
					if (md.hydros.IsEmpty()) 
						throw Exc(t_("No file loaded"));
					int lastId = md.hydros.GetCount() - 1;
					md.hydros[lastId].hd().Report();
				} else if (command[i] == "-cl" || command[i] == "--clear") {
					md.hydros.Clear();
					Cout() << "\n" << t_("Series cleared");
				} else if (command[i] == "-c" || command[i] == "--convert") {
					if (md.hydros.IsEmpty()) 
						throw Exc(t_("No file loaded"));
					i++;
					CheckNumArgs(command, i, "--convert");
					
					String file = command[i];
					
					md.hydros[0].hd().SaveAs(file);
					Cout() << "\n" << Format(t_("File '%s' converted"), file);
				} else if (command[i] == "-p" || command[i] == "--params") {
					i++;
					CheckNumArgs(command, i, "--params");
					
					while (i < command.GetCount()) {
						if (command[i] == "g") {
							i++;
							CheckNumArgs(command, i, "-p g");
							double g = ScanDouble(command[i]);
							if (IsNull(g))
								throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
							md.g = g;
						} else if (command[i] == "length") {
							i++;
							CheckNumArgs(command, i, "-p length");
							double length = ScanDouble(command[i]);
							if (IsNull(length))
								throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
							md.length = length;
						} else if (command[i] == "rho") {
							i++;
							CheckNumArgs(command, i, "-p rho");
							double rho = ScanDouble(command[i]);
							if (IsNull(rho))
								throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
							md.rho = rho;
						} else if (command[i] == "depth") {
							i++;
							CheckNumArgs(command, i, "-p depth");
							double depth = ScanDouble(command[i]);
							if (IsNull(depth))
								throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
							md.depth = depth;
						} else 
							throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
					}
				} else 
					throw Exc(Format(t_("Unknown argument '%s'"), command[i]));
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
		Cerr() << Format("\n%s: %s", t_("Error"), errorStr);
}

#endif