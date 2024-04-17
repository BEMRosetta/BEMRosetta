// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2023, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "FastOut.h"
#ifdef PLATFORM_WIN32
#include "orca.h"
#endif
#include <ScatterDraw/ScatterDraw.h>

void SetBuildInfo(String &str) {
	String name, mode;
	Time date;
	int version, bits;
	GetCompilerInfo(name, version, date, mode, bits);
	str.Replace("BUILDINFO", Format("%4d%02d%02d%02d, %s, %d bits", 
				date.year, date.month, date.day, date.hour, mode, bits)); 
}

String GetSystemInfo() {
	String name, mode;
	Time date;
	int version, bits;
	GetCompilerInfo(name, version, date, mode, bits);

	String systemInfo;
	systemInfo << Format(t_("BEMRosetta is at '%s'"), GetExeFilePath());
	systemInfo << "\n" << Format(t_("Build date is %s"), Format(date));
	systemInfo << "\n" << Format(t_("Compiler is %s, version %d, mode %s, %d bits"), name, version, mode, bits);
	
	return systemInfo;
}

UVector<String> GetCommandLineParams(String str) {
	UVector<String> ret;
	
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

void CheckIfAvailableArg(const UVector<String>& command, int i, String param) {
	if (i >= command.size())	
		throw Exc(Format(t_("Missing parameters when reading '%s'"), param));
}

void ShowHelp(BEM &md) {
	Cout() << "\n" << t_("Usage: bemrosetta_cl [options]");
	Cout() << "\n";
	Cout() << "\n" << t_("Options:");
	
	Cout() << "\n" << t_("-h  -help                 # Print options");
	Cout() << "\n" << t_("-echo \"text\"              # Prints \"text\"");
	Cout() << "\n" << t_("-pause");
	Cout() << "\n" << t_("-exit                     # Exits BEMRosetta");
	
	Cout() << "\n" << t_("-general                  # The next commands are for any data (default)");
	Cout() << "\n" << t_("-paramfile <file>         # Params in a file. New lines are like separators and # indicates a comment");
	Cout() << "\n" << t_("-params <param> <value>   # Set physical parameters:");
	Cout() << "\n" << t_("               g          # gravity       [m/s2]  ") << md.g;
	Cout() << "\n" << t_("               rho        # water density [kg/m3] ") << md.rho;
	Cout() << "\n" << t_("-echo off/on              # Show text messages");
	Cout() << "\n" << t_("-csvseparator <sep>       # Sets the separator for .csv files");
	Cout() << "\n" << t_("-isEqual <value>          # Stops if last print is not equal to <value>");
	Cout() << "\n" << t_("-isSimilar <value>        # Stops if <value> is not included in last print");
	
	Cout() << "\n";
	Cout() << "\n" << t_("-bem                      # The next commands are for BEM data");
	Cout() << "\n" << t_("-i  -input <file>         # Load model");
	Cout() << "\n" << t_("-c  -convert <file>       # Export actual model to output file");
	Cout() << "\n" << t_("-convQTFHeads <params>    # Set heading save config. for Wamit .1");
	Cout() << "\n" << t_("               all        # All headings");
	Cout() << "\n" << t_("               allNoCross # All headings but crossed");
	Cout() << "\n" << t_("               <h1> <h2>  # Indicated pair of headings [deg]");
	Cout() << "\n" << t_("-setid <id>               # Set the id of the default BEM model");
	Cout() << "\n" << t_("-params <param> <value>   # Set parameters:");
	Cout() << "\n" << t_("               length     # length scale  []      ") << md.len;
	Cout() << "\n" << t_("               depth      # water depth   [m]     ") << md.depth;
	Cout() << "\n" << t_("-delHead   <h1> ...       # Delete forces for indicated headings [deg]");
	Cout() << "\n" << t_("-delHeadId <h1> ...       # Delete forces for indicated heading ids");
	Cout() << "\n" << t_("-delButHead   <h1> ...    # Delete forces for all but indicated headings [deg]");
	Cout() << "\n" << t_("-delButHeadId <h1> ...    # Delete forces for all but indicated heading ids");
	Cout() << "\n" << t_("-delQTFHead   <h1> <h2> ..# Delete QTF for indicated pairs of headings [deg]");
	Cout() << "\n" << t_("-delQTFHeadId <h1> <h2> ..# Delete QTF for indicated heading ids");
	Cout() << "\n" << t_("-delButQTFHead   <h1> <h2>..# Delete QTF for all but indicated pairs of headings [deg]");
	Cout() << "\n" << t_("-delButQTFHeadId <h1> <h2>..# Delete QTF for all but indicated heading ids");
	Cout() << "\n" << t_("-p  -print <params>       # Prints model data in a row");
	Cout() << "\n" << t_("        <params> nb                  # Number of bodies  []");
	Cout() << "\n" << t_("                 nf                  # Number of frequencies []");
	Cout() << "\n" << t_("                 nh                  # Number of headings    []");
	Cout() << "\n" << t_("                 w                   # List of frequencies   [rad/s]");
	Cout() << "\n" << t_("                 headings            # List of headings      [deg]");
	Cout() << "\n" << t_("                 a <dof1> <dof2>     # List of A(w)(6*Nb, 6*Nb) [Kg]");
	Cout() << "\n" << t_("                 b <dof1> <dof2>     # List of B(w)(6*Nb, 6*Nb)");
	Cout() << "\n" << t_("                 ainf                # Ainf(6*Nb, 6*Nb)  [Kg]");
	Cout() << "\n" << t_("                 Theave <ibody>      # Heave resonance period for body ib [s]");
	Cout() << "\n" << t_("                 Troll <ibody>       # Roll resonance period for body ib [s]");
	Cout() << "\n" << t_("                 Tpitch <ibody>      # Pitch resonance period for body ib [s]");
	Cout() << "\n" << t_("                 GMroll <ibody>      # GM in roll [m]");
	Cout() << "\n" << t_("                 GMpitch <ibody>     # GM in pitch [m]");
	Cout() << "\n" << t_("-r  -report               # Output last loaded model data");
	Cout() << "\n" << t_("-cl -clear                # Clear loaded model");
	Cout() << "\n";
	Cout() << "\n" << t_("-mesh                     # The next commands are for mesh data");
	Cout() << "\n" << t_("-i  -input <file>         # Load model");
	Cout() << "\n" << t_("-c  -convert <file> <opts># Export actual model to output file");
	Cout() << "\n" << t_("        <opts> symx       # - Save only positive X");
	Cout() << "\n" << t_("               symy       # - Save only positive Y");
	Cout() << "\n" << t_("               symx symy  # - Save only positive X and Y");
	Cout() << "\n" << t_("               Wamit.gdf  # - Save in Wamit .gdf format");
	Cout() << "\n" << t_("               Nemoh.dat  # - Save in Nemoh .dat format");
	Cout() << "\n" << t_("               HAMS.pnl   # - Save in HAMS  .pnl format");
	Cout() << "\n" << t_("               STL.Text   # - Save in STL   text format");
	Cout() << "\n" << t_("-t   -translate <x> <y> <z>     # Translate x, y, z [m]");
	Cout() << "\n" << t_("-rot -rotate   <ax> <ay> <az> <cx> <cy> <cz>  # Rotate angle ax, ay, az [deg] around cx, cy, cz [m]");
	Cout() << "\n" << t_("-cg            <x> <y> <z>      # Sets cg: x, y, z [m] cg is the centre of gravity");
	Cout() << "\n" << t_("-c0            <x> <y> <z>      # Sets c0: x, y, z [m] c0 is the centre of motion");
	Cout() << "\n" << t_("-mass          <value>          # Sets the body mass [kg]");
	
	Cout() << "\n" << t_("-getwaterplane                  # Extract in new model the waterplane mesh (lid)");
	Cout() << "\n" << t_("-gethull                        # Extract in new model the mesh underwater hull");
	
	Cout() << "\n" << t_("-setid <id>                     # Set the id of the default mesh model");
	Cout() << "\n" << t_("-reset                          # Restore the mesh to the initial situation");
	Cout() << "\n" << t_("-r  -report                     # Output last loaded model data");
	Cout() << "\n" << t_("-p  -print <params>             # Prints model data in a row");
	Cout() << "\n" << t_("     <params> volume            # volx voly volx [m3]");
	Cout() << "\n" << t_("              underwatervolume  # volx voly volx [m3]");
	Cout() << "\n" << t_("              surface           # [m2]");
	Cout() << "\n" << t_("              underwatersurface # [m2]");
	Cout() << "\n" << t_("              cb                # cbx cby cbz [m]");
	Cout() << "\n" << t_("              cg_vol            # cgx cgy cgz [m]");
	Cout() << "\n" << t_("              cg_surf           # cgx cgy cgz [m]");
	Cout() << "\n" << t_("              hydrostiffness    # Hydrostatic stiffness around the rotation centre");
	Cout() << "\n" << t_("                                # returns K(3,3) [N/m]");
	Cout() << "\n" << t_("                                #         K(3,4) K(3,5) K(4,3) [N/rad]");
	Cout() << "\n" << t_("                                #         K(4,4) K(4,5) K(4,6) [Nm/rad]");
	Cout() << "\n" << t_("                                #         K(5,3) [N/rad]");
	Cout() << "\n" << t_("                                #         K(5,4) K(5,5) K(5,6) K(6,4) K(6,5) K(6,6) [Nm/rad]");
	Cout() << "\n" << t_("              hydrostatic_force # Hydrostatic force around the motion centre");
	Cout() << "\n" << t_("                                # returns Fx, Fy, Fz [N]");
	Cout() << "\n" << t_("                                #         Mx(roll), My(pitch), Mz(yaw) [N·m]");
	Cout() << "\n" << t_("              inertia      <cx> <cy> <cz>");
	Cout() << "\n" << t_("              inertia_vol  <cx> <cy> <cz> # Inertia tensor around cx, cy, cz [m], considering the mass spread in the volume");
	Cout() << "\n" << t_("              inertia_surf <cx> <cy> <cz> # Inertia tensor around cx, cy, cz [m], considering the mass spread in the surface");
	Cout() << "\n" << t_("                                # returns Ixx Ixy Ixz Iyx Iyy Iyz Izx Izy Izz [m2]");
	Cout() << "\n" << t_("              GZ <angle> <from> <to> <delta> # GZ around angle [deg] (0 is around Y axis), from-to-delta [deg]");
	Cout() << "\n" << t_("                                # returns the set of angles [deg] and their gz values [m]");
	Cout() << "\n" << t_("              GM                # returns GMpitch GMroll [m]");

	Cout() << "\n" << t_("-cl -clear                      # Clear loaded model");
	Cout() << "\n";
	Cout() << "\n" << t_("-time                           # The next commands are for time series");
	Cout() << "\n" << t_("-i  -input <file>               # Load file");
	Cout() << "\n" << t_("-c  -convert <file>             # Export actual model to output file");
	Cout() << "\n" << t_("-p  -print <params>             # Prints file data in a row");
	Cout() << "\n" << t_("     <params> list              # Parameter names");
	Cout() << "\n" << t_("              <param> <time>    # Value of <param> in <time>");
	Cout() << "\n" << t_("              <param> avg       # <param> avg");
	Cout() << "\n" << t_("              <param> max       # <param> max");
	Cout() << "\n" << t_("              <param> min       # <param> min");	
	
	Cout() << "\n";
	Cout() << "\n" << t_("-wind                           # The next commands are for wind series");
	Cout() << "\n" << t_("-i  -input <file>               # Load file");
	Cout() << "\n" << t_("-c  -convert <file>             # Export actual model to output file");
	Cout() << "\n" << t_("-setid <id>                     # Set the id of the default BEM model");
	Cout() << "\n" << t_("-params <param> <value>         # Set parameters:");
	Cout() << "\n" << t_("               hubheight        # Hub heignt [m]       ");
	Cout() << "\n" << t_("               gridheight       # Grid base height [m] ");
	Cout() << "\n" << t_("-p  -print <params>             # Prints model data in a row");
	Cout() << "\n" << t_("              time              # Time series [s]");
	Cout() << "\n" << t_("              vel <y> <z> <time># Wind speed norm [m/s] at y, z [m] in <time>");
	Cout() << "\n" << t_("              vel <y> <z> data  # Wind speed norm [m/s] at y, z [m] data series");
	Cout() << "\n" << t_("              vel <y> <z> avg   # Wind speed norm [m/s] at y, z [m] average");
	
	Cout() << "\n" << t_("-r  -report                     # Output loaded model main data");
	Cout() << "\n" << t_("-ra -reportall                  # Output all models main data");
	Cout() << "\n" << t_("-cl -clear                      # Clear loaded model");	
	
#ifdef PLATFORM_WIN32
	Cout() << "\n";
	Cout() << "\n" << t_("-orca                           # The next commands are for OrcaFlex handling (Required to be installed)");
	Cout() << "\n" << t_("-isAvailable                    # Test if OrcaFLEX is installed and the license is available");
	Cout() << "\n" << t_("-rw -runwave <from> <to>        # OrcaWave calculation with <from>, results in <to>");
	Cout() << "\n" << t_("-numtries <num>                 # Number <num> of attempts to connect to the license");
	Cout() << "\n" << t_("-numthread <num>                # Set the number of threads for OrcaWave");
	Cout() << "\n" << t_("-timelog <num>                  # Set minimum clock seconds per each OrcaFlex simulation log");
	Cout() << "\n" << t_("-rf -runflex <from> <to>        # OrcaFlex calculation with <from>, results in <to>");
	Cout() << "\n" << t_("-ls -loadSim <sim>              # Load OrcaFlex simulation in <sim>");
	Cout() << "\n" << t_("-rs -refsystem <cx> <cy> <cz>   # Sets the coordinates of the reference system for outputs");
	Cout() << "\n" << t_("-lp -loadParam <param>          # Load param to be saved");
	Cout() << "\n" << t_("-cp -clearParams                # Clear parameters loaded");
	Cout() << "\n" << t_("-c  -convert <file>             # Export the loaded params of the actual model to output file (csv)");
	Cout() << "\n" << t_("-p  -print <params>             # Prints model data in a row");
	Cout() << "\n" << t_("              numthread         # Number of threads");
	Cout() << "\n" << t_("     <params> list              # Parameter names");
	Cout() << "\n" << t_("              <param> data      # <param> data series");
	Cout() << "\n" << t_("              <param> avg       # <param> avg");
	Cout() << "\n" << t_("              <param> max       # <param> max");
	Cout() << "\n" << t_("              <param> min       # <param> min");	
	Cout() << "\n" << t_("-dll <file/folder>              # File or folder where OrcaFlex .dll is located. Use if detection doesn't work");		
#endif
	Cout() << "\n";
	Cout() << "\n" << t_("The actions:");
	Cout() << "\n" << t_("- are done in sequence: if a physical parameter is changed after export, saved files will not include the change");
	Cout() << "\n" << t_("- can be repeated as desired");
}

static bool NoPrint(String, int) {return true;}

String FileName(String file) {
	file.Replace("*EXEFOLDER*", GetExeFolder());
	return file;
}


#ifdef PLATFORM_WIN32
Function<bool(String, int, const Time &)> Orca::WhenWave = [](String str, int perc, const Time &et)->bool {
	if (IsNull(et))
		BEM::Print("\nCompleted 0%"); 
	else
		BEM::Print(Format("\nCompleted %d%%. Et: %", perc, et)); 
	return 0;
};

Function<bool(String)> Orca::WhenPrint = [](String str)->bool {
	BEM::Print(Format("\n%s", str)); 
	return 0;
};

Time Orca::startCalc = Null, Orca::beginNoLicense = Null, Orca::lastLog;
int64 Orca::noLicenseTime = 0;

#endif

bool ConsoleMain(const UVector<String>& _command, bool gui, Function <bool(String, int pos)> Status) {	
	UVector<String> command = clone(_command);
	
	SetConsoleColor(CONSOLE_COLOR::LTCYAN);
	Cout() << "BEMRosetta";
	SetConsoleColor(CONSOLE_COLOR::PREVIOUS);
	
	String str = S(". ") + t_("Copyright (c) 2024. Hydrodynamic coefficients converter for Boundary Element Method solver formats\nVersion beta BUILDINFO");
	SetBuildInfo(str);
	Cout() << str;
	
	ChangeCurrentDirectory(GetExeFilePath());
	
	BEM bem;
	FastOut fast;
	ArrayWind wind;
	
	bool returnval = true;
	
#ifdef PLATFORM_WIN32
	Orca orca;
	bool dllOrcaLoaded = false;
	UVector<String> paramList;
	UArray<Point3D> centreList;
	Point3D centre = Point3D(0, 0, 0);
	int numOrcaTries = 4;
#endif

	String errorJson = bem.LoadSerializeJson();
	bool firstTime = !errorJson.IsEmpty();
	if (firstTime) 
		Cout() << "\n" << errorJson << "\n" << t_("BEM configuration data is not loaded. Defaults values are set"); 

	
	UVector<String> headParams(2);
	
	bool echo = true;
	String lastPrint;
	String errorStr;
	try {
		if (command.IsEmpty()) {
			Cout() << "\n" << t_("Command argument list is empty");
			ShowHelp(bem);
		} else {
			String nextcommands = "general";
			int bemid = -1, meshid = -1, windid = -1;
			for (int i = 0; i < command.size(); i++) {
				String param = ToLower(command[i]);
				if (param == "-general") 
					nextcommands = "general";
				else if (param == "-bem") 
					nextcommands = "bem";
				else if (param == "-mesh") 
					nextcommands = "mesh";
				else if (param == "-fast" || param == "-time") 
					nextcommands = "time";
				else if (param == "-orca") 
					nextcommands = "orca";
				else if (param == "-wind") 
					nextcommands = "wind";
				else if (param == "-h" || param == "-help") {
					ShowHelp(bem);
					break;
				} else if (param == "-pause") 
					ReadStdIn();
				else if (param == "-exit") 
					break;
				else if (param == "-echo") {
					CheckIfAvailableArg(command, ++i, "-echo");
					
					Cout() << Replace(command[i], "\\n", "\n");
				} else if (param == "-isequal") { 
					CheckIfAvailableArg(command, ++i, "-isequal");
					
					String data = Trim(command[i]);
					if (Trim(lastPrint) == data) 
						BEM::Print("\n" + Format(t_("Last print is equal to \"%s\""), data));
					else
						throw Exc(Format(t_("Last print is not equal to \"%s\""), data));
				} else if (param.StartsWith("-issimilar")) { 
					CheckIfAvailableArg(command, ++i, "-isSimilar");
					
					UVector<String> data = Split(command[i], " ");
					UVector<String> prnt = Split(lastPrint, " ");
					bool isnum = param == "-issimilarnum";
					for (int i = 0; i < min(data.size(), prnt.size()); ++i) {
						if (data[i] != "*") {
							if (isnum) {
								double p = ScanDouble(prnt[i]),
									   d = ScanDouble(data[i]);
								if (!EqualRatio(p, d, 0.001))
									throw Exc(Format(t_("Last print is not equal to \"%s\". \"%s\" != \"%s\""), Trim(command[i]), prnt[i], data[i])); 	
							} else {
								if (prnt[i].Find(data[i]) < 0) 
									throw Exc(Format(t_("Last print is not equal to \"%s\". \"%s\" != \"%s\""), Trim(command[i]), prnt[i], data[i])); 	
							}
						}
					}
					BEM::Print("\n" + Format(t_("Last print is similar to \"%s\""), Trim(command[i])));
				} else {
					if (nextcommands == "general") {
					 	if (param == "-paramfile") {
							CheckIfAvailableArg(command, ++i, "-paramfile");
							
							String paramfile = command[i];
							String strfile = LoadFile(paramfile);
							if (strfile.IsVoid())
								throw Exc(Format("-paramfile file '%s' not found", paramfile));
							ChangeCurrentDirectory(GetFileFolder(paramfile));
							command.Insert(i+1 , pick(GetCommandLineParams(strfile)));
						} else if (param == "-params") {
							CheckIfAvailableArg(command, i+1, "-params");
							
							while (command.size() > i+1 && !command[i+1].StartsWith("-")) {
								++i;
								if (ToLower(command[i]) == "g") {
									CheckIfAvailableArg(command, ++i, "-p g");
									double g = ScanDouble(command[i]);
									if (IsNull(g))
										throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
									bem.g = g;
									BEM::Print("\n" + Format(t_("Gravity is %f"), g));	
								} else if (ToLower(command[i]) == "rho") {
									CheckIfAvailableArg(command, ++i, "-p rho");
									double rho = ScanDouble(command[i]);
									if (IsNull(rho))
										throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
									bem.rho = rho;
									BEM::Print("\n" + Format(t_("Density is %f"), rho));	
								} else 
									throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
							}
						} else if (param == "-echo") {
							CheckIfAvailableArg(command, ++i, "-echo");
							static Function <void(String)> OldPrint, OldWarning;
							
							String onoff = ToLower(command[i]);
							if (onoff == "on") {
								echo = true;
								if (OldPrint) {
									BEM::Print = OldPrint;
									BEM::PrintWarning = OldWarning;
								}
							} else if (onoff == "off") {
								echo = false;
								OldPrint = BEM::Print;
								BEM::Print.Clear();
								OldWarning = BEM::PrintWarning;
								BEM::PrintWarning.Clear();
							} else
								throw Exc(Format(t_("Unknown argument '%s'"), command[i]));
						} else if (param == "-csvseparator") {
							CheckIfAvailableArg(command, ++i, "-csvseparator");
							
							String sep = ToLower(command[i]);				
							
							ScatterDraw::SetDefaultCSVSeparator(sep);
						} else
							throw Exc(Format(t_("Unknown argument '%s'"), command[i]));
					} else if (nextcommands == "bem") {
						if (param == "-i" || param == "-input") {
							CheckIfAvailableArg(command, ++i, "--input");
							
							String file = FileName(command[i]);
							if (!FileExists(file)) 
								throw Exc(Format(t_("File '%s' not found"), file)); 
							
							bem.LoadBEM(file, echo ? Status : NoPrint, false);
							bemid = bem.hydros.size() - 1;
							BEM::Print("\n" + Format(t_("File '%s' loaded"), file));
						} else if (param == "-r" || param == "-report") {
							if (bem.hydros.IsEmpty()) 
								throw Exc(t_("Report: No file loaded"));
							bem.hydros[bemid].hd().Report();
						} else if (param == "-cl" || param == "-clear") {
							bem.hydros.Clear();
							bemid = -1;
							BEM::Print("\n" + S(t_("BEM data cleared")));	
						} else if (param == "-setid") {
							if (bem.hydros.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							
							CheckIfAvailableArg(command, ++i, "-setid");

							bemid = ScanInt(command[i]);
							if (IsNull(bemid) || bemid < 0 || bemid > bem.hydros.size()-1)
								throw Exc(Format(t_("Invalid id %s"), command[i]));
							BEM::Print("\n" + Format(t_("BEM active model id is %d"), bemid));	
						} else if (param == "-c" || param == "-convert") {
							Hydro &hd = bem.hydros[bemid].hd();
							if (bem.hydros.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							
							CheckIfAvailableArg(command, ++i, "-convert");
							
							String file = FileName(command[i]);
							int qtfHeading;
							if (headParams[0] == "all" || IsEmpty(headParams[0]))
								qtfHeading = Null;
							else if (headParams[0] == "allnocross")
								qtfHeading = -1;
							else 
								qtfHeading = FindClosest(hd.qh, std::complex<double>(ScanDouble(headParams[0]), ScanDouble(headParams[1])));
							if (bem.hydros[bemid].hd().SaveAs(file, echo ? Status : NoPrint, Hydro::UNKNOWN, qtfHeading)) 
								BEM::Print("\n" + Format(t_("Model id %d saved as '%s'"), bemid, file));
						} else if (param == "-convqtfheads") {	
							CheckIfAvailableArg(command, ++i, "-convqtfheads");
							if (ToLower(command[i]) == "all") 
								headParams[0] = "all";
							else if (ToLower(command[i]) == "allnocross") 
								headParams[0] = "allnocross";
							else {
								double d;
								d = ScanDouble(command[i]);
								if (IsNull(d))
									throw Exc(Format(t_("Wrong 1st heading '%s'"), command[i]));
								headParams[0] = command[i];
							
								CheckIfAvailableArg(command, ++i, "-convqtfheads 2nd head");
								d = ScanDouble(command[i]);
								if (IsNull(d))
									throw Exc(Format(t_("Wrong 2nd heading '%s'"), command[i]));
								headParams[1] = command[i];
							}
						} else if (param == "-delhead") {	
							Hydro &hd = bem.hydros[bemid].hd();
							UVector<int> ids;
							double head;
							while (i+1 < command.size() && !IsNull(head = ScanDouble(command[i+1]))) {
								int id = FindClosest(hd.head, head);
								FindAdd(ids, id);
								i++;
							}
							hd.DeleteHeadings(ids);	
						} else if (param == "-delheadid") {	
							Hydro &hd = bem.hydros[bemid].hd();
							UVector<int> ids;
							int id;
							while (i+1 < command.size() && !IsNull(id = ScanInt(command[i+1]))) {
								if (id < 0 || id >= hd.head.size())
									throw Exc(Format(t_("Wrong head id '%s'"), command[i+1]));
								FindAdd(ids, id);
								i++;
							}
							hd.DeleteHeadings(ids);	
						} else if (param == "-delbuthead") {
							Hydro &hd = bem.hydros[bemid].hd();	
							UVector<int> ids;
							double head;
							while (i+1 < command.size() && !IsNull(head = ScanDouble(command[i+1]))) {
								int id = FindClosest(hd.head, head);
								FindAdd(ids, id);
								i++;
							}
							UVector<int> idsDel;
							for (int i = 0; i < hd.head.size(); ++i)
								if (Find(ids, i) < 0)
									idsDel << i;
							hd.DeleteHeadings(idsDel);	
						} else if (param == "-delbutheadid") {	
							Hydro &hd = bem.hydros[bemid].hd();
							UVector<int> ids;
							int id;
							while (i+1 < command.size() && !IsNull(id = ScanInt(command[i+1]))) {
								if (id < 0 || id >= hd.head.size())
									throw Exc(Format(t_("Wrong head id '%s'"), command[i+1]));
								FindAdd(ids, id);
							}
							UVector<int> idsDel;
							for (int i = 0; i < hd.head.size(); ++i)
								if (Find(ids, i) < 0)
									idsDel << i;
							hd.DeleteHeadings(idsDel);	
						} else if (param == "-delqtfhead") {
							Hydro &hd = bem.hydros[bemid].hd();
							UVector<int> ids;
							double head1, head2;
							while (i+1 < command.size() && !IsNull(head1 = ScanDouble(command[i+1]))) {
								i++;
								CheckIfAvailableArg(command, ++i, "-delqtfhead 2nd head");
								if (IsNull(head2 = ScanDouble(command[i])))
									throw Exc(Format(t_("Wrong head '%s'"), command[i]));
								int id = FindClosest(hd.qh, std::complex<double>(head1, head2));
								FindAdd(ids, id);
								i++;
							}
							hd.DeleteHeadingsQTF(ids);
						} else if (param == "-delqtfheadid") {	
							Hydro &hd = bem.hydros[bemid].hd();
							UVector<int> ids;
							int id;
							while (i+1 < command.size() && !IsNull(id = ScanInt(command[i+1]))) {
								if (id < 0 || id >= hd.qh.size())
									throw Exc(Format(t_("Wrong head id '%s'"), command[i+1]));
								FindAdd(ids, id);
								i++;
							}
							hd.DeleteHeadingsQTF(ids);	
						} else if (param == "-delbutqtfhead") {
							Hydro &hd = bem.hydros[bemid].hd();
							UVector<int> ids;
							double head1, head2;
							while (i+1 < command.size() && !IsNull(head1 = ScanDouble(command[i+1]))) {
								i++;
								CheckIfAvailableArg(command, ++i, "-delqtfhead 2nd head");
								if (IsNull(head2 = ScanDouble(command[i])))
									throw Exc(Format(t_("Wrong head '%s'"), command[i]));
								int id = FindClosest(hd.qh, std::complex<double>(head1, head2));
								FindAdd(ids, id);
							}
							UVector<int> idsDel;
							for (int i = 0; i < hd.qh.size(); ++i)
								if (Find(ids, i) < 0)
									idsDel << i;
							hd.DeleteHeadingsQTF(idsDel);
						} else if (param == "-delbutqtfheadid") {
							Hydro &hd = bem.hydros[bemid].hd();	
							UVector<int> ids;
							int id;
							while (i+1 < command.size() && !IsNull(id = ScanInt(command[i+1]))) {
								if (id < 0 || id >= hd.qh.size())
									throw Exc(Format(t_("Wrong head id '%s'"), command[i+1]));
								FindAdd(ids, id);
								i++;
							}
							UVector<int> idsDel;
							for (int i = 0; i < hd.qh.size(); ++i)
								if (Find(ids, i) < 0)
									idsDel << i;
							hd.DeleteHeadingsQTF(idsDel);	
						} else if (param == "-params") {
							CheckIfAvailableArg(command, i+1, "-params");
							
							while (command.size() > i+1 && !command[i+1].StartsWith("-")) {
								++i;
								if (ToLower(command[i]) == "length") {
									CheckIfAvailableArg(command, ++i, "-p length");
									double len = ScanDouble(command[i]);
									if (IsNull(len))
										throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
									bem.len = len;
									BEM::Print("\n" + Format(t_("length is %f"), len));	
								} else if (ToLower(command[i]) == "depth") {
									CheckIfAvailableArg(command, ++i, "-p depth");
									double depth = ScanDouble(command[i]);
									if (IsNull(depth))
										throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
									bem.depth = depth;
									BEM::Print("\n" + Format(t_("depth is %f"), depth));	
								} else 
									throw Exc(Format(t_("Wrong command '%s'"), command[i]));
							}
						} else if (param == "-p" || param == "-print") {
							Hydro &data = bem.hydros[bemid].hd();
							while (command.size() > i+1 && !command[i+1].StartsWith("-")) {
								i++;
								String param = ToLower(command[i]);
								if (param == "nb") {
									Cout() << "\n";
									BEM::Print(t_("Nb:") + S(" ")); 
									lastPrint = FormatInt(data.Nb);
									Cout() << lastPrint;
								} else if (param == "nh") {
									Cout() << "\n";
									BEM::Print(t_("Nh:") + S(" ")); 
									lastPrint = FormatInt(data.Nh);
									Cout() << lastPrint;
								} else if (param == "nf") {
									Cout() << "\n";
									BEM::Print(t_("Nf:") + S(" ")); 
									lastPrint = FormatInt(data.Nf);
									Cout() << lastPrint;
								} else if (param == "w" || param == "frequencies") {
									Cout() << "\n";
									BEM::Print(t_("w:") + S(" ")); 
									lastPrint.Clear();
									for (double d : data.w)
										lastPrint << d << " "; 
									Cout() << lastPrint;
								} else if (param == "headings") {
									Cout() << "\n";
									BEM::Print(t_("head:") + S(" ")); 
									lastPrint.Clear();
									for (double d : data.head)
										lastPrint << d << " "; 
									Cout() << lastPrint; 
								} else if (param == "a") {
									CheckIfAvailableArg(command, ++i, "A #dof 1");									
									int idf = ScanInt(command[i]);
									if (idf < 1 || idf > 6*data.Nb)
										throw Exc(Format(t_("Wrong dof in '%s'"), command[i]));
									CheckIfAvailableArg(command, ++i, "A #dof 2");									
									int jdf = ScanInt(command[i]);
									if (jdf < 1 || jdf > 6*data.Nb)
										throw Exc(Format(t_("Wrong dof in '%s'"), command[i]));
									Cout() << "\n";
									BEM::Print(Format(t_("A(%d,%d):"), idf, jdf) + S(" "));
									idf--;	jdf--;
									lastPrint.Clear();
									for (int ifr = 0; ifr < data.Nf; ++ifr) 
										lastPrint << Format("%f ", data.A_dim(ifr, idf, jdf));
									Cout() << lastPrint;
								} else if (param == "ainf") {
									Cout() << "\n";
									BEM::Print(t_("Ainf:") + S(" "));
									lastPrint.Clear();
									for (int idf = 0; idf < 6*data.Nb; ++idf) 
										for (int jdf = 0; jdf < 6*data.Nb; ++jdf) 
											lastPrint << Format("%f ", data.Ainf_dim(idf, jdf));
									Cout() << lastPrint;
								} else if (param == "b") {
									CheckIfAvailableArg(command, ++i, "B #dof 1");									
									int idf = ScanInt(command[i]);
									if (idf < 1 || idf > 6*data.Nb)
										throw Exc(Format(t_("Wrong dof in '%s'"), command[i]));
									idf--;
									CheckIfAvailableArg(command, ++i, "B #dof 2");									
									int jdf = ScanInt(command[i]);
									if (jdf < 1 || jdf > 6*data.Nb)
										throw Exc(Format(t_("Wrong dof in '%s'"), command[i]));
									jdf--;
									Cout() << "\n";
									BEM::Print(t_("B:") + S(" "));
									lastPrint.Clear();
									for (int ifr = 0; ifr < data.Nf; ++ifr) 
										lastPrint << Format("%f ", data.B_dim(ifr, idf, jdf));
									Cout() << lastPrint;
								} else if (param == "theave") {
									CheckIfAvailableArg(command, ++i, "Id. body");
									int ib = ScanInt(command[i]);
									if (ib < 1 || ib > data.Nb)
										throw Exc(Format(t_("Wrong body id. in '%s'"), command[i]));
									Cout() << "\n";
									BEM::Print(Format(t_("Theave(%d):"), ib) + " "); 
									lastPrint = Format("%f", data.Theave(ib-1));
									Cout() << lastPrint;
								} else if (param == "troll") {
									CheckIfAvailableArg(command, ++i, "Id. body");
									int ib = ScanInt(command[i]);
									if (ib < 1 || ib > data.Nb)
										throw Exc(Format(t_("Wrong body id. in '%s'"), command[i]));
									Cout() << "\n";
									BEM::Print(Format(t_("Troll(%d):"), ib) + " "); 
									lastPrint = Format("%f", data.Troll(ib-1));
									Cout() << lastPrint;
								} else if (param == "tpitch") {
									CheckIfAvailableArg(command, ++i, "Id. body");
									int ib = ScanInt(command[i]);
									if (ib < 1 || ib > data.Nb)
										throw Exc(Format(t_("Wrong body id. in '%s'"), command[i]));
									Cout() << "\n";
									BEM::Print(Format(t_("Tpitch(%d):"), ib) + " "); 
									lastPrint = Format("%f", data.Tpitch(ib-1));
									Cout() << lastPrint;
								} else if (param == "gmroll") {
									CheckIfAvailableArg(command, ++i, "Id. body");
									int ib = ScanInt(command[i]);
									if (ib < 1 || ib > data.Nb)
										throw Exc(Format(t_("Wrong body id. in '%s'"), command[i]));
									Cout() << "\n";
									BEM::Print(Format(t_("GMroll(%d):"), ib) + " "); 
									lastPrint = Format("%f", data.GMroll(ib-1));
									Cout() << lastPrint;
								} else if (param == "gmpitch") {
									CheckIfAvailableArg(command, ++i, "Id. body");
									int ib = ScanInt(command[i]);
									if (ib < 1 || ib > data.Nb)
										throw Exc(Format(t_("Wrong body id. in '%s'"), command[i]));
									Cout() << "\n";
									BEM::Print(Format(t_("GMpitch(%d):"), ib) + " "); 
									lastPrint = Format("%f", data.GMpitch(ib-1));
									Cout() << lastPrint;
								} else
									throw Exc(Format(t_("Unknown argument '%s'"), command[i]));
							}
						} else 
							throw Exc(Format(t_("Unknown argument '%s'"), command[i]));
					} else if (nextcommands == "mesh") {
						if (param == "-i" || param == "-input") {
							CheckIfAvailableArg(command, ++i, "--input");
							
							String file = FileName(command[i]);
							if (!FileExists(file)) 
								throw Exc(Format(t_("File '%s' not found"), file)); 
							
							bem.LoadMesh(file, echo ? Status : NoPrint, false, false);		// Doesn't work for multibody .dat
							meshid = bem.surfs.size() - 1;
							BEM::Print("\n" + Format(t_("File '%s' loaded"), file));
						} else if (param == "-r" || param == "-report") {
							if (bem.surfs.IsEmpty()) 
								throw Exc(t_("Report: No file loaded"));
							bem.surfs[meshid].Report(bem.rho);
						} else if (param == "-cl" || param == "-clear") {
							bem.surfs.Clear();
							meshid = -1;
							BEM::Print("\n" + S(t_("Mesh data cleared")));	
						} else if (param == "-setid") {
							if (bem.surfs.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							CheckIfAvailableArg(command, ++i, "-setid");

							meshid = ScanInt(command[i]);
							if (IsNull(meshid) || meshid < 0 || meshid > bem.surfs.size()-1)
								throw Exc(Format(t_("Invalid id %s"), command[i]));
							BEM::Print("\n" + Format(t_("Mesh active model id is %d"), bemid));	
						} else if (param == "-c" || param == "-convert") {
							if (bem.surfs.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							CheckIfAvailableArg(command, ++i, "-convert");
							
							String file = FileName(command[i]);
							
							bool symX = false, symY = false;
							Mesh::MESH_FMT meshFmt = Mesh::UNKNOWN;
							while (command.size() > i+1 && !command[i+1].StartsWith("-")) {
								i++;
								String param = ToLower(command[i]);
								if (param == "symx")
									symX = true;
								else if (param == "symy")
									symY = true;
								else if (param == "symxy")
									symX = symY = true;
								else if (Mesh::GetCodeMeshStr(param) != Mesh::UNKNOWN) {
									meshFmt = Mesh::GetCodeMeshStr(param);
									if (!Mesh::meshCanSave[meshFmt])
										throw Exc(Format(t_("Saving format '%s' is not implemented"), param));									
								} else
									throw Exc(Format(t_("Unknown argument '%s'"), command[i]));
							}
							Mesh::SaveAs(bem.surfs[meshid], file, meshFmt, Mesh::ALL, bem.rho, bem.g, symX, symY);
							BEM::Print("\n" + Format(t_("Model id %d saved as '%s'"), meshid, file));
						} else if (param == "-t" || param == "-translate") {
							CheckIfAvailableArg(command, ++i, "x");
							double x = ScanDouble(command[i]);
							CheckIfAvailableArg(command, ++i, "y");
							double y = ScanDouble(command[i]);
							CheckIfAvailableArg(command, ++i, "z");
							double z = ScanDouble(command[i]);
							Mesh &data = bem.surfs[meshid];
							data.mesh.Translate(x, y, z);
							data.cg.Translate(x, y, z);
							data.AfterLoad(bem.rho, bem.g, false, false);	
							BEM::Print("\n" + Format(t_("Mesh id %d translated %f, %f, %f"), meshid, x, y, z)); 
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
							Mesh &data = bem.surfs[meshid];
							data.mesh.Rotate(ToRad(ax), ToRad(ay), ToRad(az), cx, cy, cz);	
							data.cg.Rotate(ToRad(ax), ToRad(ay), ToRad(az), cx, cy, cz);
							data.AfterLoad(bem.rho, bem.g, false, false);	
							BEM::Print("\n" + Format(t_("Mesh id %d rotated angles %f, %f, %f around centre %f, %f, %f"), meshid, ax, ay, az, cx, cy, cz));
						} else if (param == "-cg") {
							CheckIfAvailableArg(command, ++i, "cgx");
							double x = ScanDouble(command[i]);
							CheckIfAvailableArg(command, ++i, "cgy");
							double y = ScanDouble(command[i]);
							CheckIfAvailableArg(command, ++i, "cgz");
							double z = ScanDouble(command[i]);
							Mesh &data = bem.surfs[meshid];
							data.cg.x = x;
							data.cg.y = y;
							data.cg.z = z;
							data.AfterLoad(bem.rho, bem.g, true, false);
							BEM::Print("\n" + Format(t_("CG is %f, %f, %f"), x, y, z));
						} else if (param == "-c0") {
							CheckIfAvailableArg(command, ++i, "c0x");
							double x = ScanDouble(command[i]);
							CheckIfAvailableArg(command, ++i, "c0y");
							double y = ScanDouble(command[i]);
							CheckIfAvailableArg(command, ++i, "c0z");
							double z = ScanDouble(command[i]);
							Mesh &data = bem.surfs[meshid];
							data.c0.x = x;
							data.c0.y = y;
							data.c0.z = z;
							data.AfterLoad(bem.rho, bem.g, true, false);
							BEM::Print("\n" + Format(t_("CG is %f, %f, %f"), x, y, z));
						} else if (param == "-mass") { 
							CheckIfAvailableArg(command, ++i, "mass");
							double mass = ScanDouble(command[i]);
							Mesh &mesh = bem.surfs[meshid];
							mesh.SetMass(mass);
							mesh.AfterLoad(bem.rho, bem.g, true, false);
							BEM::Print("\n" + Format(t_("Mass is %f"), mass));
						} else if (param == "-reset") {	
							bem.surfs[meshid].Reset(bem.rho, bem.g);
							BEM::Print("\n" + Format(t_("Mesh id %d position is reset"), meshid));
						} else if (param == "-getwaterplane") {
							bem.AddWaterSurface(meshid, 'e');
							meshid = bem.surfs.size() - 1;
							BEM::Print("\n" + Format(t_("Mesh id %d waterplane is got"), meshid));
						} else if (param == "-gethull") {
							bem.AddWaterSurface(meshid, 'r');
							meshid = bem.surfs.size() - 1;
							BEM::Print("\n" + Format(t_("Mesh id %d hull is got"), meshid));
						} else if (param == "-p" || param == "-print") {
							Mesh &data = bem.surfs[meshid];
							while (command.size() > i+1 && !command[i+1].StartsWith("-")) {
								i++;
								String param = ToLower(command[i]);
								if (param == "volume") {
									Cout() << "\n";
									BEM::Print(t_("Volume:") + S(" ")); 
									lastPrint = Format("%f %f %f", data.mesh.volumex, data.mesh.volumey, data.mesh.volumez);
									Cout() << lastPrint;
								} else if (param == "underwatervolume") {
									Cout() << "\n";
									BEM::Print(t_("UnderwaterVolume:") + S(" ")); 
									lastPrint = Format("%f %f %f", data.under.volumex, data.under.volumey, data.under.volumez);
									Cout() << lastPrint;
								} else if (param == "surface") {
									Cout() << "\n";
									BEM::Print(t_("Surface:") + S(" ")); 
									lastPrint = Format("%f", data.mesh.surface);
									Cout() << lastPrint;
								} else if (param == "underwatersurface") {
									Cout() << "\n";
									BEM::Print(t_("UnderwaterSurface:") + S(" ")); 
									lastPrint = Format("%f", data.under.surface);
									Cout() << lastPrint;
								} else if (param == "cb") {
									Cout() << "\n";
									BEM::Print(t_("CB:") + S(" ")); 
									lastPrint = Format("%f %f %f", data.cb.x, data.cb.y, data.cb.z);
									Cout() << lastPrint;
								} else if (param == "cg_vol") {
									Cout() << "\n";
									BEM::Print(t_("CG_vol:") + S(" ")); 
									Point3D cg = data.mesh.GetCentreOfBuoyancy();
									lastPrint = Format("%f %f %f", cg.x, cg.y, cg.z);
									Cout() << lastPrint;
								} else if (param == "cg_surf") {
									Cout() << "\n";
									BEM::Print(t_("CG_surf:") + S(" ")); 
									Point3D cg = data.mesh.GetCentreOfGravity_Surface();
									lastPrint = Format("%f %f %f", cg.x, cg.y, cg.z);
									Cout() << lastPrint;
								} else if (param == "hydrostiffness") {
									Cout() << "\n";
									BEM::Print(t_("HydrostaticStiffness:") + S(" "));
									lastPrint.Clear();
									for (int i = 0; i < 6; ++i) {
										for (int j = 0; j < 6; ++j) {
											if (!Hydro::C_units(i, j).IsEmpty()) 
												lastPrint << data.C(i, j) << " ";
										}
									}
									int idColor = data.under.VolumeMatch(bem.volWarning/100., bem.volError/100.);
									if (idColor == -1)
										lastPrint << ". " << t_("Mesh warning_ Maybe incomplete");
									else if (idColor == -2)
										lastPrint << ". " << t_("Mesh error: Probably incomplete");
									Cout() << lastPrint;
								} else if (param == "hydrostatic_force") {
									Cout() << "\n";
									BEM::Print(t_("HydrostaticForce:") + S(" "));
									lastPrint.Clear();
									Force6D f = data.under.GetHydrostaticForce(data.c0, bem.rho, bem.g);
									for (int i = 0; i < 6; ++i) 				
										lastPrint << f[i] << " ";
									int idColor = data.under.VolumeMatch(bem.volWarning/100., bem.volError/100.);
									if (idColor == -1)
										lastPrint << ". " << t_("Mesh warning: Maybe incomplete");
									else if (idColor == -2)
										lastPrint << ". " << t_("Mesh error: Probably incomplete");
									Cout() << lastPrint;
								} else if (param.StartsWith("inertia")) {
									Cout() << "\n";
									BEM::Print(t_("Inertia:") + S(" "));
									lastPrint.Clear();
									Eigen::Matrix3d inertia;
									Point3D centre;
									CheckIfAvailableArg(command, ++i, "Inertia cx");
									centre.x = ScanDouble(command[i]);
									CheckIfAvailableArg(command, ++i, "Inertia cy");
									centre.y = ScanDouble(command[i]);
									CheckIfAvailableArg(command, ++i, "Inertia cz");
									centre.z = ScanDouble(command[i]);
									bool volume = false;
									if (param == "inertia" || param == "inertia_vol")
										volume = true;
									
									data.mesh.GetInertia33_Radii(inertia, centre, volume, false);
									for (int i = 0; i < 3; ++i) 
										for (int j = 0; j < 3; ++j) 
											lastPrint << inertia(i, j) << " ";
									Cout() << lastPrint;
								} else if (param == "gz") {
									Cout() << "\n";
									BEM::Print(t_("GZ:") + S(" "));
									lastPrint.Clear();
									
									CheckIfAvailableArg(command, ++i, "Angle");
									double angle = ScanDouble(command[i]);
									CheckIfAvailableArg(command, ++i, "From");
									double from = ScanDouble(command[i]);
									CheckIfAvailableArg(command, ++i, "To");
									double to = ScanDouble(command[i]);
									CheckIfAvailableArg(command, ++i, "Delta");
									double delta = ScanDouble(command[i]);
									UVector<double> dataangle, datagz;
									double tolerance = 0.1;
									String errors;
									data.GZ(from, to, delta, angle, bem.rho, bem.g, tolerance, dataangle, datagz, errors);									
									for (int i = 0; i < dataangle.size(); ++i) 
										lastPrint << dataangle[i] << " ";
									lastPrint << "\n";
									for (int i = 0; i < datagz.size(); ++i) 
										lastPrint << datagz[i] << " ";
									Cout() << lastPrint;
								} else if (param == "gm") {
									Cout() << "\n";
									BEM::Print(t_("gm:") + S(" ")); 
									lastPrint = Format("%f %f", data.GMpitch(bem.rho, bem.g), data.GMroll(bem.rho, bem.g));
									Cout() << lastPrint;
								} else
									throw Exc(Format(t_("Unknown argument '%s'"), command[i]));
							}
						} else 
							throw Exc(Format(t_("Unknown argument '%s'"), command[i]));
					} else if (nextcommands == "time") {
						if (param == "-i" || param == "-input") {
							CheckIfAvailableArg(command, ++i, "-input");

							String file = FileName(command[i]);
							if (!FileExists(file)) 
								throw Exc(Format(t_("File '%s' not found"), file)); 
							
							String ret = fast.Load(file);
							if (ret.IsEmpty())
								BEM::Print("\n" + Format(t_("File '%s' loaded"), file));
							else
								BEM::PrintWarning("\n" + Format(t_("Problem loading '%s': %s"), file, ret));
						} else if (param == "-c" || param == "-convert") {
							if (fast.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							
							CheckIfAvailableArg(command, ++i, "-convert");
							
							String file = FileName(command[i]);
							
							if (fast.Save(file, "", ScatterDraw::GetDefaultCSVSeparator())) 
								BEM::Print("\n" + Format(t_("Results saved as '%s'"), file));
						} else if (param == "-p" || param == "-print") {
							while (command.size() > i+1 && !command[i+1].StartsWith("-")) {
								i++;
								String param = ToLower(command[i]);
								if (param == "list") {
									Cout() << "\n";
									BEM::Print(t_("List of parameters:") + S(" ")); 
									lastPrint = "";
									for (int i = 0; i < fast.GetParamCount(); ++i) {
										if (i > 0)
											lastPrint << " ";
										lastPrint << fast.GetParameter(i);
									}
									Cout() << lastPrint;
								} else {
									Cout() << "\n";
									UVector<int> p = fast.FindParameterMatch(param);
									if (p.IsEmpty())
										throw Exc(Format(t_("Parameter '%s' not found"), param));

									int id = p[0];
									i++;
									double time = ScanDouble(command[i]);
									if (!IsNull(time))
										lastPrint = FormatDouble(fast.GetVal(time, id));
									else {
										VectorXd d = fast.GetVector(id);
										CleanNAN_safe(d);		// Just in case. It do happens
										if (command[i] == "avg" || command[i] == "mean") 
											lastPrint = FormatDouble(d.mean());
										else if (command[i] == "max") 
											lastPrint = FormatDouble(d.maxCoeff());
										else if (command[i] == "min") 
											lastPrint = FormatDouble(d.minCoeff());
										else  
											throw Exc(Format(t_("Parameter '%s' not found"), command[i]));
									}
									Cout() << lastPrint;
								}
							}
						} else
							throw Exc(Format(t_("Unknown argument '%s'"), command[i]));
					} else if (nextcommands == "wind") {
						if (param == "-i" || param == "-input") {
							CheckIfAvailableArg(command, ++i, "--input");
							
							String file = FileName(command[i]);
							if (!FileExists(file)) 
								throw Exc(Format(t_("File '%s' not found"), file)); 
							
							Wind &w = wind.Add();
							String ret;
							if (!IsEmpty(ret = w.Load(file)))
								throw Exc(ret);
							
							windid = wind.size() - 1;
							BEM::Print("\n" + Format(t_("File '%s' loaded"), file));
						} else if (param == "-params") {
							CheckIfAvailableArg(command, i+1, "-params");
							
							while (command.size() > i+1 && !command[i+1].StartsWith("-")) {
								++i;
								if (ToLower(command[i]) == "hubheight") {
									CheckIfAvailableArg(command, ++i, "-p hubheight");
									double h = ScanDouble(command[i]);
									if (IsNull(h))
										throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
									wind[windid].SetHubHeight(h);
									BEM::Print("\n" + Format(t_("Hub height is %f"), h));	
								} else if (ToLower(command[i]) == "gridheight") {
									CheckIfAvailableArg(command, ++i, "-p gridheight");
									double h = ScanDouble(command[i]);
									if (IsNull(h))
										throw Exc(Format(t_("Wrong argument '%s'"), command[i]));
									wind[windid].SetGridHeight(h);
									BEM::Print("\n" + Format(t_("Grid height is %f"), h));	
								} else 
									throw Exc(Format(t_("Wrong command '%s'"), command[i]));
							}
						} else if (param == "-r" || param == "-report") {
							if (wind.IsEmpty()) 
								throw Exc(t_("Report: No file loaded"));
							Grid grid;
							wind[windid].Report(grid);
							BEM::Print("\n" + grid.GetString(true, false, " "));
						} else if (param == "-ra" || param == "-reportall") {
							if (wind.IsEmpty()) 
								throw Exc(t_("Report: No file loaded"));
							Grid grid;
							wind.Report(grid);
							BEM::Print("\n" + grid.GetString(true, false, " "));
						} else if (param == "-cl" || param == "-clear") {
							wind.Clear();
							windid = -1;
							BEM::Print("\n" + S(t_("Wind data cleared")));	
						} else if (param == "-setid") {
							if (wind.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							
							CheckIfAvailableArg(command, ++i, "-setid");

							windid = ScanInt(command[i]);
							if (IsNull(windid) || windid < 0 || windid > wind.size()-1)
								throw Exc(Format(t_("Invalid id %s"), command[i]));
							BEM::Print("\n" + Format(t_("Wind active model id is %d"), windid));	
						} else if (param == "-c" || param == "-convert") {
							if (wind.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							
							CheckIfAvailableArg(command, ++i, "-convert");
							
							String file = FileName(command[i]);
							
							String ret;
							if (!IsEmpty(ret = wind[windid].Save(file)))
								throw Exc(ret);
							BEM::Print("\n" + Format(t_("Model id %d saved as '%s'"), windid, file));
						} else if (param == "-p" || param == "-print") {
							Wind &w = wind[windid];
							while (command.size() > i+1 && !command[i+1].StartsWith("-")) {
								i++;
								String param = ToLower(command[i]);
								if (param == "vel") {
									BEM::Print("\n");
									BEM::Print(t_("Velocity:") + S(" "));
									lastPrint.Clear();
									
									CheckIfAvailableArg(command, ++i, "Position Y");
									double ypos = ScanDouble(command[i]);
									CheckIfAvailableArg(command, ++i, "Position Z");
									double zpos = ScanDouble(command[i]);
									int idz, idy;
									w.GetPos(zpos, ypos, idz, idy);
									Eigen::VectorXd v;
									v = w.Get(idz, idy);
									
									CheckIfAvailableArg(command, ++i, "<time>/avg/data");
									String what = command[i];
									if (what == "avg")
										lastPrint << v.mean();
									else if (what == "data") {
										for (int i = 0; i < v.size(); ++i) 
											lastPrint << v(i) << " ";
									} else {
										double time = ScanDouble(what);
										if (IsNull(time))
											throw Exc(t_("Parameter time not found in vel"));
										int id = w.GetTimeId(time);
										lastPrint << v(id);
									}
									BEM::Print(lastPrint);
								} else if (param == "time") {
									BEM::Print("\n");
									BEM::Print(t_("Time:") + S(" "));
									lastPrint.Clear();
									Eigen::VectorXd t;
									t = w.GetTime();
									for (int i = 0; i < t.size(); ++i) 
										lastPrint << t(i) << " ";
									BEM::Print(lastPrint);
								}
							}
						} else 
							throw Exc(Format(t_("Unknown argument '%s'"), command[i]));
					}
#ifdef PLATFORM_WIN32						
					else if (nextcommands == "orca") {
						if (param == "-rw" || param == "-runwave") {
							if (!dllOrcaLoaded) {
								if (orca.FindInit()) 
									dllOrcaLoaded = true;		
							}
							if (!dllOrcaLoaded)
								throw Exc(t_("DLL not found"));
									
							CheckIfAvailableArg(command, ++i, "-runwave from");
							String from = command[i];
							CheckIfAvailableArg(command, ++i, "-runwave to");
							String to = command[i];
							
							int numTry = numOrcaTries;
							String errorStr;
							do {
								try {
									orca.LoadWave(from);
									break;
								} catch (Exc e) {
									errorStr = e;
									if (errorStr.Find("license") < 0)
										break;
									numTry--;
									errorStr.Clear();
								}
								BEM::PrintWarning("\n" + Format(t_("Next try (%d/%d)"), numOrcaTries-numTry+1, numOrcaTries));
							} while (numTry > 0);
							
							if (!IsEmpty(errorStr))
								BEM::PrintWarning("\n" + Format(t_("Problem loading '%s'. %s"), from, errorStr));
							else {
								orca.RunWave();
								orca.SaveWaveResults(to);							
							
								BEM::Print("\n" + Format(t_("Diffraction results saved at '%s'"), to));
							}
						} else if (param == "-rf" || param == "-runflex") {
							if (!dllOrcaLoaded) {
								if (orca.FindInit()) 
									dllOrcaLoaded = true;		
							}
							if (!dllOrcaLoaded)
								throw Exc(t_("DLL not found"));
									
							CheckIfAvailableArg(command, ++i, "-runflex from");
							String from = command[i];
							CheckIfAvailableArg(command, ++i, "-runflex to");
							String to = command[i];
							
							int numTry = numOrcaTries;
							String errorStr;
							do {
								try {
									orca.LoadFlex(from);
									break;
								} catch (Exc e) {
									errorStr = e;
									if (errorStr.Find("license") < 0)
										break;
									numTry--;
									errorStr.Clear();
								}
								BEM::PrintWarning("\n" + Format(t_("Next try (%d/%d)"), numOrcaTries-numTry+1, numOrcaTries));
							} while (numTry > 0);
							
							if (!IsEmpty(errorStr))
								BEM::PrintWarning("\n" + Format(t_("Problem loading '%s'. %s"), from, errorStr));
							else {
								bool saved = false;
								try {
									orca.RunFlex(to, saved);
								} catch (Exc e) {
									Cerr() << "\n" << Format(t_("Error running OrcaFlex: %s"), e);
									returnval = false;
								}
								if (orca.GetModelState() == msSimulationStoppedUnstable)
									Cerr() << "\n" << t_("Simulation aborted: The simulation has become unstable because the solver parameters (time step, iterations num.) or other, have to ve reviewed.");
		
								if (saved)
									BEM::Print("\n" + Format(t_("Simulation results saved at '%s'"), to));
							}
						} else if (param == "-ls" || param == "-loadsim") {
							if (!dllOrcaLoaded) {
								if (orca.FindInit()) 
									dllOrcaLoaded = true;		
							}
							if (!dllOrcaLoaded)
								throw Exc(t_("DLL not found"));
									
							CheckIfAvailableArg(command, ++i, "-loadSim sim");
							String from = command[i];
							
							int numTry = numOrcaTries;
							String errorStr;
							do {
								try {
									orca.LoadFlexSim(from);
									break;
								} catch (Exc e) {
									errorStr = e;
									if (errorStr.Find("license") < 0)
										break;
									numTry--;
									errorStr.Clear();
								}
								BEM::PrintWarning("\n" + Format(t_("Next try (%d/%d)"), numOrcaTries-numTry+1, numOrcaTries));
							} while (numTry > 0);
							
							if (!IsEmpty(errorStr))
								throw Exc(Format(t_("Problem loading '%s'. %s"), from, errorStr));
							else
								BEM::Print("\n" + Format(t_("Simulation '%s' loaded"), from));
						} else if (param == "-numthread") {
							if (!dllOrcaLoaded) {
								if (orca.FindInit()) 
									dllOrcaLoaded = true;		
							}
							if (!dllOrcaLoaded)
								throw Exc(t_("DLL not found"));
							
							CheckIfAvailableArg(command, ++i, "-numthread num");
							int numth = ScanInt(command[i]);
							if (IsNull(numth))
								throw Exc(Format(t_("Invalid thread number %s"), command[i]));							

							orca.SetThreadCount(numth);
							
							BEM::Print("\n" + Format(t_("Thread number set to %d"), numth));
						} else if (param == "-numtries") {
							if (!dllOrcaLoaded) {
								if (orca.FindInit()) 
									dllOrcaLoaded = true;		
							}
							if (!dllOrcaLoaded)
								throw Exc(t_("DLL not found"));
							
							CheckIfAvailableArg(command, ++i, "-numtries num");
							numOrcaTries = ScanInt(command[i]);
							if (IsNull(numOrcaTries))
								throw Exc(Format(t_("Invalid thread number %s"), command[i]));							

							BEM::Print("\n" + Format(t_("Number of OrcaFlex license tries set to %d"), numOrcaTries));
						} else if (param == "-timelog") {
							if (!dllOrcaLoaded) {
								if (orca.FindInit()) 
									dllOrcaLoaded = true;		
							}
							if (!dllOrcaLoaded)
								throw Exc(t_("DLL not found"));
							
							CheckIfAvailableArg(command, ++i, "-timelog num");
							int timeLog = int(StringToSeconds(command[i]));
							if (IsNull(timeLog))
								throw Exc(Format(t_("Invalid time between logs %s"), command[i]));							

							BEM::Print("\n" + Format(t_("OrcaFlex time between logs set to %.1f sec"), timeLog));							
							
							orca.deltaLogSimulation = timeLog;
						} else if (param == "-lp" || param == "-loadparam") {
							CheckIfAvailableArg(command, ++i, "-loadParam");
							
							paramList << command[i];
							centreList << centre;
						} else if (param == "-c" || param == "-convert") {
							if (!dllOrcaLoaded) {
								if (orca.FindInit()) 
									dllOrcaLoaded = true;		
							}
							if (!dllOrcaLoaded)
								throw Exc(t_("DLL not found"));
							
							CheckIfAvailableArg(command, ++i, "-convert");
							
							String file = FileName(command[i]);

							orca.SaveCsv(file, paramList, centreList, ScatterDraw::GetDefaultCSVSeparator());
							BEM::Print("\n" + Format(t_("Results saved as '%s'"), file));
						} else if (param == "-rs" || param == "-refsystem") {	
							CheckIfAvailableArg(command, ++i, "Reference cx");
							centre.x = ScanDouble(command[i]);
							CheckIfAvailableArg(command, ++i, "Reference cy");
							centre.y = ScanDouble(command[i]);
							CheckIfAvailableArg(command, ++i, "Reference cz");
							centre.z = ScanDouble(command[i]);
						} else if (param == "-cp" || param == "-clearparams") {
							paramList.Clear();
							centreList.Clear();
						} else if (param == "-isavailable") {
							if (!dllOrcaLoaded) {
								if (orca.FindInit()) 
									dllOrcaLoaded = true;		
							}
							if (!dllOrcaLoaded) {
								BEM::PrintWarning(S("\n") + t_("OrcaFLEX is not installed"));
								break;
							}
							int numTry = 2;
							do {
								if (orca.IsAvailable())
									break;
								BEM::PrintWarning("\n" + Format(t_("Next try (%d/%d)"), 2-numTry+1, 2));
								numTry--;
							} while (numTry > 0);
							if (numTry == 0) {
								BEM::Print(S("\n") + t_("OrcaFLEX license is not available"));
								break;
							}
							BEM::Print(S("\n") + t_("OrcaFLEX is installed and available"));
						} else if (param == "-p" || param == "-print") {
							if (!dllOrcaLoaded) {
								if (orca.FindInit()) 
									dllOrcaLoaded = true;		
							}
							if (!dllOrcaLoaded)
								throw Exc(t_("DLL not found"));
							
							while (command.size() > i+1 && !command[i+1].StartsWith("-")) {
								i++;
								String param = ToLower(command[i]);
								if (param == "list") {
									Cout() << "\n";
									UVector<String> list = orca.GetFlexSimVariablesList();
									BEM::Print(Format(t_("List of parameters (%d):"), list.size()) + S(" ")); 
									lastPrint = "";
									for (int i = 0; i < list.size(); ++i) {
										if (i > 0)
											lastPrint << ", ";
										lastPrint << list[i];
									}
									Cout() << lastPrint;
								} else if (param == "numthread") {
									int numth = orca.GetThreadCount();
									BEM::Print("\n" + Format("%d", numth));
								} else {
									Cout() << "\n";

									String unit;
									int objType;
									VectorXd d;
									orca.GetFlexSimVar(command[i], centre, unit, objType, d);
									
									i++;
									if (command.size() <= i)
										throw Exc(Format(t_("Last command '%s' is not complete"), param));
									if (command[i] == "data") {
										lastPrint = "";
										for (int i = 0; i < d.size(); ++i) {
											if (i > 0)
												lastPrint << " ";
											lastPrint << d[i];
										}
									} else if (command[i] == "avg" || command[i] == "mean") 
										lastPrint = FormatDouble(d.mean());
									else if (command[i] == "max") 
										lastPrint = FormatDouble(d.maxCoeff());
									else if (command[i] == "min") 
										lastPrint = FormatDouble(d.minCoeff());
									else  
										throw Exc(Format(t_("Parameter '%s' not found in -p"), command[i]));
									Cout() << lastPrint;
								}
							}
						} else if (param == "-dll") {
							CheckIfAvailableArg(command, ++i, "-dll");

							String file = FileName(command[i]);
							if (!FileExists(file)) {
								file = AFX(file, "OrcFxAPI.dll");
								if (!FileExists(file)) 	
									throw Exc(Format(t_("File '%s' not found"), file)); 
							}
							if (!orca.Init(file))
								throw Exc(Format(t_("DLL '%s' not found"), command[i]));
							dllOrcaLoaded = true;
							BEM::Print("\n" + Format(t_("Orca .dll '%s' loaded"), orca.GetDLLPath()));
						} else 
							throw Exc(Format(t_("Unknown argument '%s'"), command[i]));
					}
#endif
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
		Cerr() << "\n" << Format(t_("Error: %s"), errorStr);
		Cerr() << S("\n\n") + t_("In case of doubt try option -h or -help") + S("\n");
		if (gui)
			Cerr() << S("\n") + t_("or just call command line without arguments to open GUI window");
		returnval = false;
	}
	Cout() << "\n";
	return returnval;
}
