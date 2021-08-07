#include "BEMRosetta.h"



void SetBuildInfo(String &str) {
	String name, mode;
	Time date;
	int version, bits;
	GetCompilerInfo(name, version, date, mode, bits);
	str.Replace("BUILDINFO", Format("%4d%02d%02d%02d, %s, %d bits", 
				date.year, date.month, date.day, date.hour, mode, bits)); 
}

Vector<String> GetCommandLineParams(String str) {
	Vector<String> ret;
	
	bool inQuotes = false, inComment = false;
	String tempstr;
	for (int i = 0; i < str.GetCount(); ++i) {
		const int &c = str[i];
		if (!inComment && c == '\"')
			inQuotes = !inQuotes;
		else if (!inComment && inQuotes)
			tempstr.Cat(c);
		else if (c == '#')
		 	inComment = true;
		else if (c == ' ' || c == '\t' || c == '\r' || c == '\n') {
			if (inComment && c == '\n')
				inComment = false;
			else if (!tempstr.IsEmpty()) {
				ret << tempstr;
				tempstr.Clear();
			}
		} else if (!inComment)
			tempstr.Cat(c);
	}
	if (!tempstr.IsEmpty()) 
		ret << tempstr;

	return ret;
}

void CheckIfAvailableArg(const Vector<String>& command, int i, String param) {
	if (i >= command.size())	
		throw Exc(Format(t_("Missing parameters when reading '%s'"), param));
}

void ShowHelp(BEMData &md) {
	Cout() << "\n" << t_("Usage: bemrosetta_cl [options]");
	Cout() << "\n";
	Cout() << "\n" << t_("Options:");
	Cout() << "\n" << t_("-h  -help             | print options");
	Cout() << "\n" << t_("-paramfile <file> | Params in a file");
	
	Cout() << "\n" << t_("-p  -params           | set physical parameters:");
	Cout() << "\n" << t_("                        parameter description   units   default value");
	Cout() << "\n" << t_("                        g         gravity       [m/s2]  ") << md.g;
	Cout() << "\n" << t_("                        length    length scale  []      ") << md.len;
	Cout() << "\n" << t_("                        rho       water density [kg/m3] ") << md.rho;
	Cout() << "\n" << t_("                        depth     water depth   [m]     ") << md.depth;

	Cout() << "\n" << t_("-bem                  | The next commands are for BEM data");
	Cout() << "\n" << t_("-i  -input <file>     | Load model");
	Cout() << "\n" << t_("-c  -convert <file>   | Export actual model to output file");
	Cout() << "\n" << t_("-setid <id>           | Set the id of the default BEM model");
	Cout() << "\n" << t_("-r  -report           | Output last loaded model data");
	Cout() << "\n" << t_("-cl -clear            | Clear loaded model");
	
	Cout() << "\n" << t_("-mesh                 | The next commands are for mesh data");
	Cout() << "\n" << t_("-i  -input <file>     | Load model");
	Cout() << "\n" << t_("-c  -convert <file>   | Export actual model to output file");
	
	Cout() << "\n" << t_("-getwaterplane        | Extract in new model the waterplane mesh (lid)");
	Cout() << "\n" << t_("-gethull              | Extract in new model the mes underwater hull");
	
	Cout() << "\n" << t_("-setid <id>           | Set the id of the default mesh model");
	Cout() << "\n" << t_("-r  -report           | Output last loaded model data");
	Cout() << "\n" << t_("-cl -clear            | Clear loaded model");
	
	Cout() << "\n";
	Cout() << "\n" << t_("The actions:");
	Cout() << "\n" << t_("- are done in sequence: if a physical parameter is changed after export, saved files will not include the change");
	Cout() << "\n" << t_("- can be repeated as desired");
}

void ConsoleMain(const Vector<String>& _command, bool gui, Function <bool(String, int pos)> Status) {	
	Vector<String> command = clone(_command);
	
	SetConsoleColor(CONSOLE_COLOR::LTCYAN);
	Cout() << "BEMRosetta";
	SetConsoleColor(CONSOLE_COLOR::PREVIOUS);
	
	String str = S(". ") + t_("Copyright (c) 2021. Hydrodynamic coefficients converter for Boundary Element Method solver formats\nVersion beta BUILDINFO");
	SetBuildInfo(str);
	Cout() << str;
	
	SetCurrentDirectory(GetExeFilePath());
	
	BEMData md;
	
	bool firstTime = false;
	if (!md.LoadSerializeJson(firstTime))
		Cout() << "\n" << t_("BEM configuration data has not been loaded. Default values are set");
	
	String errorStr;
	try {
		if (command.IsEmpty()) {
			Cout() << "\n" << t_("Command argument list is empty");
			ShowHelp(md);
		} else {
			String nextcommands = "bem";
			int bemid = -1, meshid = -1;
			for (int i = 0; i < command.size(); i++) {
				String param = ToLower(command[i]);
				if (param == "-h" || param == "-help") {
					ShowHelp(md);
					break;
				} else if (param == "-paramfile") {
					i++;
					CheckIfAvailableArg(command, i, "-paramfile");
					
					String paramfile = command[i];
					String strfile = LoadFile(paramfile);
					if (strfile.IsEmpty())
						throw Exc(Format("-paramfile file '%s' not found", paramfile));
					SetCurrentDirectory(GetFileFolder(paramfile));
					command.Insert(i+1 , pick(GetCommandLineParams(strfile)));
				} else if (param == "-bem") {
					nextcommands = "bem";
				} else if (param == "-mesh") {	
					nextcommands = "mesh";
				} else {
					if (nextcommands == "bem") {
						if (param == "-i" || param == "-input") {
							i++;
							CheckIfAvailableArg(command, i, "--input");
							
							String file = command[i];
							if (!FileExists(file)) 
								throw Exc(Format(t_("File '%s' not found"), file)); 
							
							md.LoadBEM(file, Status, false);
							bemid = md.hydros.size() - 1;
							Cout() << "\n" << Format(t_("File '%s' loaded"), file);
						} else if (param == "-r" || param == "-report") {
							if (md.hydros.IsEmpty()) 
								throw Exc(t_("Report: No file loaded"));
							md.hydros[bemid].hd().Report();
						} else if (param == "-cl" || param == "-clear") {
							md.hydros.Clear();
							bemid = -1;
							Cout() << "\n" << t_("Series cleared");
						} else if (param == "-setid") {
							if (md.hydros.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							i++;
							CheckIfAvailableArg(command, i, "-setid");

							int bemid = ScanInt(command[i]);
							if (IsNull(bemid) || bemid < 0 || bemid > md.hydros.size()-1)
								throw Exc(Format(t_("Invalid id %s"), command[i]));
						} else if (param == "-c" || param == "-convert") {
							if (md.hydros.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							i++;
							CheckIfAvailableArg(command, i, "-convert");
							
							String file = command[i];
							
							if (md.hydros[bemid].hd().SaveAs(file, Status))
								Cout() << "\n" << Format(t_("File '%s' converted"), file);
						} else if (param == "-p" || param == "-params") {
							i++;
							CheckIfAvailableArg(command, i, "--params");
							
							while (i < command.size()) {
								if (command[i] == "g") {
									i++;
									CheckIfAvailableArg(command, i, "-p g");
									double g = ScanDouble(command[i]);
									if (IsNull(g))
										throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
									md.g = g;
								} else if (command[i] == "length") {
									i++;
									CheckIfAvailableArg(command, i, "-p length");
									double len = ScanDouble(command[i]);
									if (IsNull(len))
										throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
									md.len = len;
								} else if (command[i] == "rho") {
									i++;
									CheckIfAvailableArg(command, i, "-p rho");
									double rho = ScanDouble(command[i]);
									if (IsNull(rho))
										throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
									md.rho = rho;
								} else if (command[i] == "depth") {
									i++;
									CheckIfAvailableArg(command, i, "-p depth");
									double depth = ScanDouble(command[i]);
									if (IsNull(depth))
										throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
									md.depth = depth;
								} else 
									throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
							}
						} else 
							throw Exc(Format(t_("Unknown argument '%s'"), command[i]));
					} else if (nextcommands == "mesh") {
						if (param == "-i" || param == "-input") {
							i++;
							CheckIfAvailableArg(command, i, "--input");
							
							String file = command[i];
							if (!FileExists(file)) 
								throw Exc(Format(t_("File '%s' not found"), file)); 
							
							md.LoadMesh(file, Status, false, false);
							meshid = md.surfs.size() - 1;
							Cout() << "\n" << Format(t_("File '%s' loaded"), file);
						} else if (param == "-r" || param == "-report") {
							if (md.surfs.IsEmpty()) 
								throw Exc(t_("Report: No file loaded"));
							md.surfs[meshid].Report(md.rho);
						} else if (param == "-cl" || param == "-clear") {
							md.surfs.Clear();
							meshid = -1;
							Cout() << "\n" << t_("Series cleared");
						} else if (param == "-setid") {
							if (md.surfs.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							i++;
							CheckIfAvailableArg(command, i, "-setid");

							meshid = ScanInt(command[i]);
							if (IsNull(meshid) || meshid < 0 || meshid > md.surfs.size()-1)
								throw Exc(Format(t_("Invalid id %s"), command[i]));
						} else if (param == "-c" || param == "-convert") {
							if (md.surfs.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							i++;
							CheckIfAvailableArg(command, i, "-convert");
							
							String file = command[i];
							
							bool symX = false, symY = false;
							
							if (command.size() > i+1 && !command[i+1].StartsWith("-")) {
								i++;
								String sym = ToLower(command[i]);
								if (sym == "symx")
									symX = true;
								else if (sym == "symy")
									symY = true;
								else if (sym == "symxy")
									symX = symY = true;
								else
									throw Exc(Format(t_("Unknown argument '%s'"), command[i]));
							}
							md.surfs[meshid].SaveAs(file, MeshData::UNKNOWN, md.g, MeshData::ALL, symX, symY);
							Cout() << "\n" << Format(t_("File '%s' converted"), file);
						}  else if (param == "-getwaterplane") {
							md.AddWaterSurface(meshid, 'e');
							meshid = md.surfs.size() - 1;
						} else if (param == "-gethull") {
							md.AddWaterSurface(meshid, 'r');
							meshid = md.surfs.size() - 1;
						} 
					}
				}
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
	if (!errorStr.IsEmpty()) {
		Cerr() << Format("\n%s: %s", t_("Error"), errorStr);
		Cerr() << S("\n\n") + t_("In case of doubt try option -h or --help");
		if (gui)
			Cerr() << S("\n") + t_("or just call command line without arguments to open GUI window");
	}
	Cout() << "\n";
}
