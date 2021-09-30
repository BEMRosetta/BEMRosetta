// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2021, the BEMRosetta author and contributors
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
	
	Cout() << "\n" << t_("-h  -help             # Print options");
	
	Cout() << "\n" << t_("-general              # The next commands are for any data (default)");
	Cout() << "\n" << t_("-paramfile <file>     # Params in a file. New lines are like separators and # indicates a comment");
	Cout() << "\n" << t_("-params <param> <value>   # Set physical parameters:");
	Cout() << "\n" << t_("               g         	# gravity       [m/s2]  ") << md.g;
	Cout() << "\n" << t_("               rho        # water density [kg/m3] ") << md.rho;
	Cout() << "\n" << t_("-echo off/on          # Show text messages");
	Cout() << "\n" << t_("-isEqual \"<value>\"  # Stops if last print is not equal to \"<value>\"");
	
	Cout() << "\n";
	Cout() << "\n" << t_("-bem                  # The next commands are for BEM data");
	Cout() << "\n" << t_("-i  -input <file>     # Load model");
	Cout() << "\n" << t_("-c  -convert <file>   # Export actual model to output file");
	Cout() << "\n" << t_("-setid <id>           # Set the id of the default BEM model");
	Cout() << "\n" << t_("-params <param> <value>   # Set physical parameters:");
	Cout() << "\n" << t_("               length     # length scale  []      ") << md.len;
	Cout() << "\n" << t_("               depth      # water depth   [m]     ") << md.depth;
	Cout() << "\n" << t_("-p  -print <params>   # Prints model data in a row");
	Cout() << "\n" << t_("        <params> nb                  # Number of bodies  []");
	Cout() << "\n" << t_("                 nf                  # Number of frequencies []");
	Cout() << "\n" << t_("                 nh                  # Number of headings    []");
	Cout() << "\n" << t_("                 ainf <dof1> <dof2>  # Ainf(6*Nb, 6*Nb)  [Kg]");
	Cout() << "\n" << t_("-r  -report           # Output last loaded model data");
	Cout() << "\n" << t_("-cl -clear            # Clear loaded model");
	Cout() << "\n";
	Cout() << "\n" << t_("-mesh                 # The next commands are for mesh data");
	Cout() << "\n" << t_("-i  -input <file>     # Load model");
	Cout() << "\n" << t_("-c  -convert <file> <opts>  # Export actual model to output file");
	Cout() << "\n" << t_("        <opts> symx         # - Save only positive X");
	Cout() << "\n" << t_("               symy         # - Save only positive Y");
	Cout() << "\n" << t_("               symx symy    # - Save only positive X and Y");
	Cout() << "\n" << t_("               Wamit.gdf    # - Save in Wamit .gdf format");
	Cout() << "\n" << t_("               Nemoh.dat    # - Save in Nemoh .dat format");
	Cout() << "\n" << t_("               HAMS.pnl     # - Save in HAMS  .pnl format");
	Cout() << "\n" << t_("               STL.Text     # - Save in STL   text format");
	Cout() << "\n" << t_("-t   -translate <x> <y> <z>              # Translate x, y, z [m]");
	Cout() << "\n" << t_("-rot -rotate    <ax> <ay> <az> <cx> <cy> <cz>  # Rotate angle ax, ay, az [deg] around point cx, cy, cz [m]");
	Cout() << "\n" << t_("-cg             <x> <y> <z>       # Sets cg: x, y, z [m]");
	
	Cout() << "\n" << t_("-getwaterplane        # Extract in new model the waterplane mesh (lid)");
	Cout() << "\n" << t_("-gethull              # Extract in new model the mesh underwater hull");
	
	Cout() << "\n" << t_("-setid <id>           # Set the id of the default mesh model");
	Cout() << "\n" << t_("-reset                # Restore the mesh to the initial situation");
	Cout() << "\n" << t_("-r  -report           # Output last loaded model data");
	Cout() << "\n" << t_("-p  -print <params>   # Prints model data in a row");
	Cout() << "\n" << t_("        <params> volume            # volx voly volx [m3]");
	Cout() << "\n" << t_("                 underwatervolume  # volx voly volx [m3]");
	Cout() << "\n" << t_("                 surface           # [m2]");
	Cout() << "\n" << t_("                 underwatersurface # [m2]");
	Cout() << "\n" << t_("                 cb                # cbx cby cbz [m]");
	Cout() << "\n" << t_("                 hydrostiffness");
	Cout() << "\n" << t_("                 >returns K(3,3) [N/m]");
	Cout() << "\n" << t_("                          K(3,4) K(3,5) K(4,3) [N/rad]");
	Cout() << "\n" << t_("                          K(4,4) K(4,5) K(4,6) [Nm/rad]");
	Cout() << "\n" << t_("                          K(5,3) [N/rad]");
	Cout() << "\n" << t_("                          K(5,4) K(5,5) K(5,6) K(6,4) K(6,5) K(6,6) [Nm/rad]");
	Cout() << "\n" << t_("                 inertia <cx> <cy> <cz>   # Inertia tensor around cx, cy, cz [m]");
	Cout() << "\n" << t_("                 >returns Ixx Ixy Ixz Iyx Iyy Iyz Izx Izy Izz [m2]");

	Cout() << "\n" << t_("-cl -clear            # Clear loaded model");
	
	Cout() << "\n";
	Cout() << "\n" << t_("The actions:");
	Cout() << "\n" << t_("- are done in sequence: if a physical parameter is changed after export, saved files will not include the change");
	Cout() << "\n" << t_("- can be repeated as desired");
}

static bool NoPrint(String, int) {return true;}

bool ConsoleMain(const Vector<String>& _command, bool gui, Function <bool(String, int pos)> Status) {	
	Vector<String> command = clone(_command);
	
	SetConsoleColor(CONSOLE_COLOR::LTCYAN);
	Cout() << "BEMRosetta";
	SetConsoleColor(CONSOLE_COLOR::PREVIOUS);
	
	String str = S(". ") + t_("Copyright (c) 2021. Hydrodynamic coefficients converter for Boundary Element Method solver formats\nVersion beta BUILDINFO");
	SetBuildInfo(str);
	Cout() << str;
	
	SetCurrentDirectory(GetExeFilePath());
	
	BEMData md;
	
	bool firstTime;
	if (!md.LoadSerializeJson(firstTime))
		Cout() << "\n" << t_("BEM configuration data has not been loaded. Default values are set");
	
	bool echo = true;
	String lastPrint;
	String errorStr;
	try {
		if (command.IsEmpty()) {
			Cout() << "\n" << t_("Command argument list is empty");
			ShowHelp(md);
		} else {
			String nextcommands = "general";
			int bemid = -1, meshid = -1;
			for (int i = 0; i < command.size(); i++) {
				String param = ToLower(command[i]);
				if (param == "-general") 
					nextcommands = "general";
				else if (param == "-bem") 
					nextcommands = "bem";
				else if (param == "-mesh") 
					nextcommands = "mesh";
				else if (param == "-h" || param == "-help") {
					ShowHelp(md);
					break;
				} else if (param == "-isequal") {
					CheckIfAvailableArg(command, ++i, "-isequal");
					
					String data = Trim(command[i]);
					if (Trim(lastPrint) == data) 
						BEMData::Print("\n" + Format(t_("Last print is equal to \"%s\""), data));
					else
						throw Exc(Format(t_("Last print is not equal to \"%s\""), data));
				} else {
					if (nextcommands == "general") {
					 	if (param == "-paramfile") {
							CheckIfAvailableArg(command, ++i, "-paramfile");
							
							String paramfile = command[i];
							String strfile = LoadFile(paramfile);
							if (strfile.IsEmpty())
								throw Exc(Format("-paramfile file '%s' not found", paramfile));
							SetCurrentDirectory(GetFileFolder(paramfile));
							command.Insert(i+1 , pick(GetCommandLineParams(strfile)));
						} else if (param == "-params") {
							CheckIfAvailableArg(command, i+1, "-params");
							
							while (command.size() > i+1 && !command[i+1].StartsWith("-")) {
								++i;
								if (command[i] == "g") {
									CheckIfAvailableArg(command, ++i, "-p g");
									double g = ScanDouble(command[i]);
									if (IsNull(g))
										throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
									md.g = g;
									BEMData::Print("\n" + Format(t_("Gravity is %f"), g));	
								} else if (command[i] == "rho") {
									CheckIfAvailableArg(command, ++i, "-p rho");
									double rho = ScanDouble(command[i]);
									if (IsNull(rho))
										throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
									md.rho = rho;
									BEMData::Print("\n" + Format(t_("Density is %f"), rho));	
								} else 
									throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
							}
						} else if (param == "-echo") {
							CheckIfAvailableArg(command, ++i, "-echo");
							static Function <void(String)> OldPrint, OldWarning;
							
							String onoff = ToLower(command[i]);
							if (onoff == "on") {
								echo = true;
								if (OldPrint)
									BEMData::Print = OldPrint;
									BEMData::PrintWarning = OldWarning;
							} else if (onoff == "off") {
								echo = false;
								OldPrint = BEMData::Print;
								BEMData::Print.Clear();
								OldWarning = BEMData::PrintWarning;
								BEMData::PrintWarning.Clear();
							} else
								throw Exc(Format(t_("Unknown argument '%s'"), command[i]));
						} else
							throw Exc(Format(t_("Unknown argument '%s'"), command[i]));
					} else if (nextcommands == "bem") {
						if (param == "-i" || param == "-input") {
							CheckIfAvailableArg(command, ++i, "--input");
							
							String file = command[i];
							if (!FileExists(file)) 
								throw Exc(Format(t_("File '%s' not found"), file)); 
							
							md.LoadBEM(file, echo ? Status : NoPrint, false);
							bemid = md.hydros.size() - 1;
							BEMData::Print("\n" + Format(t_("File '%s' loaded"), file));
						} else if (param == "-r" || param == "-report") {
							if (md.hydros.IsEmpty()) 
								throw Exc(t_("Report: No file loaded"));
							md.hydros[bemid].hd().Report();
						} else if (param == "-cl" || param == "-clear") {
							md.hydros.Clear();
							bemid = -1;
							BEMData::Print("\n" + S(t_("BEM data cleared")));	
						} else if (param == "-setid") {
							if (md.hydros.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							
							CheckIfAvailableArg(command, ++i, "-setid");

							int bemid = ScanInt(command[i]);
							if (IsNull(bemid) || bemid < 0 || bemid > md.hydros.size()-1)
								throw Exc(Format(t_("Invalid id %s"), command[i]));
							BEMData::Print("\n" + Format(t_("BEM active model id is %d"), bemid));	
						} else if (param == "-c" || param == "-convert") {
							if (md.hydros.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							
							CheckIfAvailableArg(command, ++i, "-convert");
							
							String file = command[i];
							
							if (md.hydros[bemid].hd().SaveAs(file, echo ? Status : NoPrint)) {
								BEMData::Print("\n" + Format(t_("Model id %d saved as '%s'"), bemid, file));
							}
						} else if (param == "-params") {
							CheckIfAvailableArg(command, i+1, "-params");
							
							while (command.size() > i+1 && !command[i+1].StartsWith("-")) {
								++i;
								if (command[i] == "length") {
									CheckIfAvailableArg(command, ++i, "-p length");
									double len = ScanDouble(command[i]);
									if (IsNull(len))
										throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
									md.len = len;
									BEMData::Print("\n" + Format(t_("length is %f"), len));	
								} else if (command[i] == "depth") {
									CheckIfAvailableArg(command, ++i, "-p depth");
									double depth = ScanDouble(command[i]);
									if (IsNull(depth))
										throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
									md.depth = depth;
									BEMData::Print("\n" + Format(t_("depth is %f"), depth));	
								} else 
									throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
							}
						} else if (param == "-p" || param == "-print") {
							Hydro &data = md.hydros[bemid].hd();
							while (command.size() > i+1 && !command[i+1].StartsWith("-")) {
								i++;
								String param = ToLower(command[i]);
								if (param == "nb") {
									Cout() << "\n";
									BEMData::Print(t_("Nb:") + S(" ")); 
									lastPrint = FormatInt(data.Nb);
									Cout() << lastPrint;
								} else if (param == "nh") {
									Cout() << "\n";
									BEMData::Print(t_("Nh:") + S(" ")); 
									lastPrint = FormatInt(data.Nh);
									Cout() << lastPrint;
								} else if (param == "nf") {
									Cout() << "\n";
									BEMData::Print(t_("Nf:") + S(" ")); 
									lastPrint = FormatInt(data.Nf);
									Cout() << lastPrint;
								} else if (param == "ainf") {
									CheckIfAvailableArg(command, ++i, "Ainf #1 dof");
									int idf = ScanInt(command[i]);
									CheckIfAvailableArg(command, ++i, "Ainf #2 dof");									
									int jdf = ScanInt(command[i]);
									if (idf < 1 || jdf < 1 || idf > 6*data.Nb || jdf > 6*data.Nb)
										throw Exc(Format(t_("Wrong dof in '%s'"), command[i]));
									Cout() << "\n";
									BEMData::Print(Format(t_("Ainf(%d,%d):"), idf, jdf) + " "); 
									lastPrint = Format("%f", data.Awinf_dim(idf-1, jdf-1));
									Cout() << lastPrint;
								}
							}
						} else 
							throw Exc(Format(t_("Unknown argument '%s'"), command[i]));
					} else if (nextcommands == "mesh") {
						if (param == "-i" || param == "-input") {
							CheckIfAvailableArg(command, ++i, "--input");
							
							String file = command[i];
							if (!FileExists(file)) 
								throw Exc(Format(t_("File '%s' not found"), file)); 
							
							md.LoadMesh(file, echo ? Status : NoPrint, false, false);
							meshid = md.surfs.size() - 1;
							BEMData::Print("\n" + Format(t_("File '%s' loaded"), file));
						} else if (param == "-r" || param == "-report") {
							if (md.surfs.IsEmpty()) 
								throw Exc(t_("Report: No file loaded"));
							md.surfs[meshid].Report(md.rho);
						} else if (param == "-cl" || param == "-clear") {
							md.surfs.Clear();
							meshid = -1;
							BEMData::Print("\n" + S(t_("Mesh data cleared")));	
						} else if (param == "-setid") {
							if (md.surfs.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							CheckIfAvailableArg(command, ++i, "-setid");

							meshid = ScanInt(command[i]);
							if (IsNull(meshid) || meshid < 0 || meshid > md.surfs.size()-1)
								throw Exc(Format(t_("Invalid id %s"), command[i]));
							BEMData::Print("\n" + Format(t_("Mesh active model id is %d"), bemid));	
						} else if (param == "-c" || param == "-convert") {
							if (md.surfs.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							CheckIfAvailableArg(command, ++i, "-convert");
							
							String file = command[i];
							
							bool symX = false, symY = false;
							MeshData::MESH_FMT meshFmt = MeshData::UNKNOWN;
							while (command.size() > i+1 && !command[i+1].StartsWith("-")) {
								i++;
								String param = ToLower(command[i]);
								if (param == "symx")
									symX = true;
								else if (param == "symy")
									symY = true;
								else if (param == "symxy")
									symX = symY = true;
								else if (MeshData::GetCodeMeshStr(param) != MeshData::UNKNOWN) {
									meshFmt = MeshData::GetCodeMeshStr(param);
									if (!MeshData::MeshCanSave(meshFmt))
										throw Exc(Format(t_("Saving format '%s' is not implemented"), param));									
								} else
									throw Exc(Format(t_("Unknown argument '%s'"), command[i]));
							}
							md.surfs[meshid].SaveAs(file, meshFmt, md.g, MeshData::ALL, symX, symY);
							BEMData::Print("\n" + Format(t_("Model id %d saved as '%s'"), meshid, file));
						} else if (param == "-t" || param == "-translate") {
							CheckIfAvailableArg(command, ++i, "x");
							double x = ScanDouble(command[i]);
							CheckIfAvailableArg(command, ++i, "y");
							double y = ScanDouble(command[i]);
							CheckIfAvailableArg(command, ++i, "z");
							double z = ScanDouble(command[i]);
							MeshData &data = md.surfs[meshid];
							data.mesh.Translate(x, y, z);
							data.cg.Translate(x, y, z);
							data.AfterLoad(md.rho, md.g, false, false);	
							BEMData::Print("\n" + Format(t_("Mesh id %d translated %f, %f, %f"), meshid, x, y, z)); 
						} else if (param == "-rot" || param == "-rotate") {
							CheckIfAvailableArg(command, ++i, "ax");
							double ax = ScanDouble(command[i]);
							CheckIfAvailableArg(command, ++i, "ay");
							double ay = ScanDouble(command[i]);
							CheckIfAvailableArg(command, ++i, "az");
							double az = ScanDouble(command[i]);
							CheckIfAvailableArg(command, ++i, "cx");
							double cx = ScanDouble(command[i]);
							CheckIfAvailableArg(command, ++i, "cy");
							double cy = ScanDouble(command[i]);
							CheckIfAvailableArg(command, ++i, "cz");
							double cz = ScanDouble(command[i]);
							MeshData &data = md.surfs[meshid];
							data.mesh.Rotate(ax, ay, az, cx, cy, cz);	
							data.cg.Rotate(ax, ay, az, cx, cy, cz);
							data.AfterLoad(md.rho, md.g, false, false);	
							BEMData::Print("\n" + Format(t_("Mesh id %d rotated angles %f, %f, %f around center %f, %f, %f"), meshid, ax, ay, az, cx, cy, cz));
						} else if (param == "-cg") {
							CheckIfAvailableArg(command, ++i, "cgx");
							double x = ScanDouble(command[i]);
							CheckIfAvailableArg(command, ++i, "cgy");
							double y = ScanDouble(command[i]);
							CheckIfAvailableArg(command, ++i, "cgz");
							double z = ScanDouble(command[i]);
							MeshData &data = md.surfs[meshid];
							data.cg.x = x;
							data.cg.y = y;
							data.cg.z = z;
							data.AfterLoad(md.rho, md.g, true, false);
							BEMData::Print("\n" + Format(t_("CG is %f, %f, %f"), x, y, z));
						} else if (param == "-reset") {	
							md.surfs[meshid].Reset(md.rho, md.g);
							BEMData::Print("\n" + Format(t_("Mesh id %d position is reset"), meshid));
						} else if (param == "-getwaterplane") {
							md.AddWaterSurface(meshid, 'e');
							meshid = md.surfs.size() - 1;
							BEMData::Print("\n" + Format(t_("Mesh id %d water plane is got"), meshid));
						} else if (param == "-gethull") {
							md.AddWaterSurface(meshid, 'r');
							meshid = md.surfs.size() - 1;
							BEMData::Print("\n" + Format(t_("Mesh id %d hull is got"), meshid));
						} else if (param == "-p" || param == "-print") {
							MeshData &data = md.surfs[meshid];
							while (command.size() > i+1 && !command[i+1].StartsWith("-")) {
								i++;
								String param = ToLower(command[i]);
								if (param == "volume") {
									Cout() << "\n";
									BEMData::Print(t_("Volume:") + S(" ")); 
									lastPrint = Format("%f %f %f", data.mesh.volumex, data.mesh.volumey, data.mesh.volumez);
									Cout() << lastPrint;
								} else if (param == "underwatervolume") {
									Cout() << "\n";
									BEMData::Print(t_("UnderwaterVolume:") + S(" ")); 
									lastPrint = Format("%f %f %f", data.under.volumex, data.under.volumey, data.under.volumez);
									Cout() << lastPrint;
								} else if (param == "surface") {
									Cout() << "\n";
									BEMData::Print(t_("Surface:") + S(" ")); 
									lastPrint = Format("%f", data.mesh.surface);
									Cout() << lastPrint;
								} else if (param == "underwatersurface") {
									Cout() << "\n";
									BEMData::Print(t_("UnderwaterSurface:") + S(" ")); 
									lastPrint = Format("%f", data.under.surface);
									Cout() << lastPrint;
								} else if (param == "cb") {
									Cout() << "\n";
									BEMData::Print(t_("CB:") + S(" ")); 
									lastPrint = Format("%f %f %f", data.cb.x,data.cb.y, data.cb.z);
									Cout() << lastPrint;
								} else if (param == "hydrostiffness") {
									Cout() << "\n";
									BEMData::Print(t_("HydrostaticStiffness:") + S(" "));
									lastPrint.Clear();
									for (int i = 0; i < 6; ++i) {
										for (int j = 0; j < 6; ++j) {
											if (!Hydro::C_units(i, j).IsEmpty()) 
												lastPrint << data.C(i, j) << " ";
										}
									}
									Cout() << lastPrint;
								} else if (param == "inertia") {
									Cout() << "\n";
									BEMData::Print(t_("Inertia:") + S(" "));
									lastPrint.Clear();
									Eigen::Matrix3d inertia;
									Point3D center;
									CheckIfAvailableArg(command, ++i, "Inertia cx");
									center.x = ScanDouble(command[i]);
									CheckIfAvailableArg(command, ++i, "Inertia cy");
									center.y = ScanDouble(command[i]);
									CheckIfAvailableArg(command, ++i, "Inertia cz");
									center.z = ScanDouble(command[i]);
									data.mesh.GetInertia(inertia, center, true);
									for (int i = 0; i < 3; ++i) 
										for (int j = 0; j < 3; ++j) 
											lastPrint << inertia(i, j) << " ";
									Cout() << lastPrint;
								} else
									throw Exc(Format(t_("Unknown argument '%s'"), command[i]));
							}
						} else 
							throw Exc(Format(t_("Unknown argument '%s'"), command[i]));
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
		Cerr() << S("\n\n") + t_("In case of doubt try option -h or -help") + S("\n");
		if (gui)
			Cerr() << S("\n") + t_("or just call command line without arguments to open GUI window");
		return false;
	}
	Cout() << "\n";
	return true;
}
