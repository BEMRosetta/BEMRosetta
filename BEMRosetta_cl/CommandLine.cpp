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

BEM &Bem() {static BEM bem;		return bem;}

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
	Cout() << "\n" << t_("-params <param> <value/s>       # Set parameters:");
	Cout() << "\n" << t_("        hubheight               # Hub heignt [m]       ");
	Cout() << "\n" << t_("        gridheight              # Grid base height [m] ");
	Cout() << "\n" << t_("        TI <u><v><w>            # Turbulence inten. [\%] (If only u, IEC61400-1 Part 1 is applied: v = 0.8u and w = 0.5u)");
	Cout() << "\n" << t_("        scale <u><v><w>         # Factor to scale velocities []");
	Cout() << "\n" << t_("        powerLaw <pl><zh>       # Power law [] at hubheight [m]");
	Cout() << "\n" << t_("        periodic                # BTS periodic [true/false]");
	Cout() << "\n" << t_("-p  -print <params>             # Prints model data in a row");
	Cout() << "\n" << t_("        time                    # Time series [s]");
	Cout() << "\n" << t_("        vel <y> <z> <time>      # Wind speed norm [m/s] at y, z [m] in <time>");
	Cout() << "\n" << t_("        vel <y> <z> data        # Wind speed norm [m/s] at y, z [m] data series");
	Cout() << "\n" << t_("        vel <y> <z> avg         # Wind speed norm [m/s] at y, z [m] average");
	Cout() << "\n" << t_("        velComp <comp> <y> <z> <time> # Wind speed [m/s] for 0 (u), 1 (v), ... at y, z [m] in <time>");
	Cout() << "\n" << t_("        velComp <comp> <y> <z> data   # Wind speed [m/s] for 0 (u), 1 (v), ... at y, z [m] data series");
	Cout() << "\n" << t_("        velComp <comp> <y> <z> avg    # Wind speed [m/s] for 0 (u), 1 (v), ... at y, z [m] average");
	
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
	
	BEM &bem = Bem();
	FastOut fast;
	ArrayWind wind;
	
	bool returnval = true;
	
#ifdef PLATFORM_WIN32
	Orca orca;
	bool dllOrcaLoaded = false;
	UVector<String> paramList;
	UArray<Point3D> centreList;
	Point3D centreOrca = Point3D(0, 0, 0);
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
			for (int ic = 0; ic < command.size(); ic++) {
				String param = ToLower(command[ic]);
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
				} else if (param == "-pause") {
					BEM::Print(S("\n") + t_("Press Enter to continue..."));
					ReadStdIn();
				} else if (param == "-exit") {
					BEM::Print(S("\n") + t_("Exit by command"));
					break;
				} else if (param == "-echo") {
					CheckIfAvailableArg(command, ++ic, "-echo");
					
					Cout() << Replace(command[ic], "\\n", "\n");
				} else if (param == "-isequal") { 
					CheckIfAvailableArg(command, ++ic, "-isequal");
					
					String data = Trim(command[ic]);
					if (Trim(lastPrint) == data) 
						BEM::Print("\n" + Format(t_("Last print is equal to \"%s\""), data));
					else
						throw Exc(Format(t_("Last print is not equal to \"%s\""), data));
				} else if (param.StartsWith("-issimilar")) { 
					CheckIfAvailableArg(command, ++ic, "-isSimilar");
					
					UVector<String> data = Split(command[ic], " ");
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
					BEM::Print("\n" + Format(t_("Last print is similar to \"%s\""), Trim(command[ic])));
				} else {
					if (nextcommands == "general") {
					 	if (param == "-paramfile") {
							CheckIfAvailableArg(command, ++ic, "-paramfile");
							
							String paramfile = command[ic];
							String strfile = LoadFile(paramfile);
							if (strfile.IsVoid())
								throw Exc(Format("-paramfile file '%s' not found", paramfile));
							ChangeCurrentDirectory(GetFileFolder(paramfile));
							command.Insert(ic+1 , pick(GetCommandLineParams(strfile)));
						} else if (param == "-params") {
							CheckIfAvailableArg(command, ic+1, "-params");
							
							while (command.size() > ic+1 && !command[ic+1].StartsWith("-")) {
								++ic;
								if (ToLower(command[ic]) == "g") {
									CheckIfAvailableArg(command, ++ic, "-p g");
									double g = ScanDouble(command[ic]);
									if (IsNull(g))
										throw Exc(Format(t_("Wrong argument '%s'"), command[ic]));
									bem.g = g;
									BEM::Print("\n" + Format(t_("Gravity is %f"), g));	
								} else if (ToLower(command[ic]) == "rho") {
									CheckIfAvailableArg(command, ++ic, "-p rho");
									double rho = ScanDouble(command[ic]);
									if (IsNull(rho))
										throw Exc(Format(t_("Wrong argument '%s'"), command[ic]));
									bem.rho = rho;
									BEM::Print("\n" + Format(t_("Density is %f"), rho));	
								} else 
									throw Exc(Format(t_("Wrong argument '%s'"), command[ic]));
							}
						} else if (param == "-echo") {
							CheckIfAvailableArg(command, ++ic, "-echo");
							static Function <void(String)> OldPrint, OldWarning;
							
							String onoff = ToLower(command[ic]);
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
								throw Exc(Format(t_("Unknown argument '%s'"), command[ic]));
						} else if (param == "-csvseparator") {
							CheckIfAvailableArg(command, ++ic, "-csvseparator");
							
							String sep = ToLower(command[ic]);				
							
							ScatterDraw::SetDefaultCSVSeparator(sep);
						} else
							throw Exc(Format(t_("Unknown argument '%s'"), command[ic]));
					} else if (nextcommands == "bem") {
						if (param == "-i" || param == "-input") {
							CheckIfAvailableArg(command, ++ic, "--input");
							
							String file = FileName(command[ic]);
							if (!FileExists(file)) 
								throw Exc(Format(t_("File '%s' not found"), file)); 
							
							bem.LoadBEM(file, echo ? Status : NoPrint, false);
							bemid = bem.hydros.size() - 1;
							BEM::Print("\n" + Format(t_("File '%s' loaded"), file));
						} else if (param == "-r" || param == "-report") {
							if (bem.hydros.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							bem.hydros[bemid].Report();
						} else if (param == "-cl" || param == "-clear") {
							bem.hydros.Clear();
							bemid = -1;
							BEM::Print("\n" + S(t_("BEM data cleared")));	
						} else if (param == "-setid") {
							if (bem.hydros.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							
							CheckIfAvailableArg(command, ++ic, "-setid");

							bemid = ScanInt(command[ic]);
							if (IsNull(bemid) || bemid < 0 || bemid > bem.hydros.size()-1)
								throw Exc(Format(t_("Invalid id %s"), command[ic]));
							BEM::Print("\n" + Format(t_("BEM active model id is %d"), bemid));	
						} else if (param == "-c" || param == "-convert") {
							Hydro &hd = bem.hydros[bemid];
							if (bem.hydros.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							
							CheckIfAvailableArg(command, ++ic, "-convert");
							
							String file = FileName(command[ic]);
							int qtfHeading;
							if (headParams[0] == "all" || IsEmpty(headParams[0]))
								qtfHeading = Null;
							else if (headParams[0] == "allnocross")
								qtfHeading = -1;
							else 
								qtfHeading = hd.dt.FindClosestQTFHead(std::complex<double>(ScanDouble(headParams[0]), ScanDouble(headParams[1])));
							bem.hydros[bemid].SaveAs(file, echo ? Status : NoPrint, Hydro::UNKNOWN, qtfHeading);
							BEM::Print("\n" + Format(t_("Model id %d saved as '%s'"), bemid, file));
						} else if (param == "-convqtfheads") {	
							CheckIfAvailableArg(command, ++ic, "-convqtfheads");
							if (ToLower(command[ic]) == "all") 
								headParams[0] = "all";
							else if (ToLower(command[ic]) == "allnocross") 
								headParams[0] = "allnocross";
							else {
								double d;
								d = ScanDouble(command[ic]);
								if (IsNull(d))
									throw Exc(Format(t_("Wrong 1st heading '%s'"), command[ic]));
								headParams[0] = command[ic];
							
								CheckIfAvailableArg(command, ++ic, "-convqtfheads 2nd head");
								d = ScanDouble(command[ic]);
								if (IsNull(d))
									throw Exc(Format(t_("Wrong 2nd heading '%s'"), command[ic]));
								headParams[1] = command[ic];
							}
						} else if (param == "-delhead") {
							if (bem.hydros.IsEmpty()) 
								throw Exc(t_("No file loaded"));	
							Hydro &hd = bem.hydros[bemid];
							UVector<int> ids;
							double head;
							while (ic+1 < command.size() && !IsNull(head = ScanDouble(command[ic+1]))) {
								int id = hd.dt.FindClosestHead(head);
								BEM::Print("\n" + Format(t_("Closest heading to %f is %f"), head, hd.dt.head[id]));
								FindAdd(ids, id);
								ic++;
							}
							hd.DeleteHeadings(ids);	
							BEM::Print("\n" + S(t_("Headings deleted")));
						} else if (param == "-delheadid") {	
							if (bem.hydros.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							Hydro &hd = bem.hydros[bemid];
							UVector<int> ids;
							int id;
							while (ic+1 < command.size() && !IsNull(id = ScanInt(command[ic+1]))) {
								if (id < 0 || id >= hd.dt.head.size())
									throw Exc(Format(t_("Wrong head id '%s'"), command[ic+1]));
								FindAdd(ids, id);
								ic++;
							}
							hd.DeleteHeadings(ids);	
							BEM::Print("\n" + S(t_("Headings deleted")));
						} else if (param == "-delbuthead") {
							if (bem.hydros.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							Hydro &hd = bem.hydros[bemid];	
							UVector<int> ids;
							double head;
							while (ic+1 < command.size() && !IsNull(head = ScanDouble(command[ic+1]))) {
								int id = hd.dt.FindClosestHead(head);
								BEM::Print("\n" + Format(t_("Closest heading to %f is %f"), head, hd.dt.head[id]));
								FindAdd(ids, id);
								ic++;
							}
							UVector<int> idsDel;
							for (int i = 0; i < hd.dt.head.size(); ++i)
								if (Find(ids, i) < 0)
									idsDel << i;
							hd.DeleteHeadings(idsDel);	
							BEM::Print("\n" + S(t_("Headings deleted")));
						} else if (param == "-delbutheadid") {
							if (bem.hydros.IsEmpty()) 
								throw Exc(t_("No file loaded"));	
							Hydro &hd = bem.hydros[bemid];
							UVector<int> ids;
							int id;
							while (ic+1 < command.size() && !IsNull(id = ScanInt(command[ic+1]))) {
								if (id < 0 || id >= hd.dt.head.size())
									throw Exc(Format(t_("Wrong head id '%s'"), command[ic+1]));
								FindAdd(ids, id);
							}
							UVector<int> idsDel;
							for (int i = 0; i < hd.dt.head.size(); ++i)
								if (Find(ids, i) < 0)
									idsDel << i;
							hd.DeleteHeadings(idsDel);
							BEM::Print("\n" + S(t_("Headings deleted")));	
						} else if (param == "-delqtfhead") {
							if (bem.hydros.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							Hydro &hd = bem.hydros[bemid];
							UVector<int> ids;
							double head1, head2;
							while (ic+1 < command.size() && !IsNull(head1 = ScanDouble(command[ic+1]))) {
								ic++;
								CheckIfAvailableArg(command, ++ic, "-delqtfhead 2nd head");
								if (IsNull(head2 = ScanDouble(command[ic])))
									throw Exc(Format(t_("Wrong head '%s'"), command[ic]));
								int id = hd.dt.FindClosestQTFHead(std::complex<double>(head1, head2));
								BEM::Print("\n" + Format(t_("Closest heading to %f:%f is %f:%f"), head1, head2, hd.dt.qhead[id].real(), hd.dt.qhead[id].imag()));
								FindAdd(ids, id);
								ic++;
							}
							hd.DeleteHeadingsQTF(ids);
							BEM::Print("\n" + S(t_("QTF headings deleted")));
						} else if (param == "-delqtfheadid") {	
							if (bem.hydros.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							Hydro &hd = bem.hydros[bemid];
							UVector<int> ids;
							int id;
							while (ic+1 < command.size() && !IsNull(id = ScanInt(command[ic+1]))) {
								if (id < 0 || id >= hd.dt.qhead.size())
									throw Exc(Format(t_("Wrong head id '%s'"), command[ic+1]));
								FindAdd(ids, id);
								ic++;
							}
							hd.DeleteHeadingsQTF(ids);	
							BEM::Print("\n" + S(t_("QTF headings deleted")));
						} else if (param == "-delbutqtfhead") {
							if (bem.hydros.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							Hydro &hd = bem.hydros[bemid];
							UVector<int> ids;
							double head1, head2;
							while (ic+1 < command.size() && !IsNull(head1 = ScanDouble(command[ic+1]))) {
								ic++;
								CheckIfAvailableArg(command, ++ic, "-delqtfhead 2nd head");
								if (IsNull(head2 = ScanDouble(command[ic])))
									throw Exc(Format(t_("Wrong head '%s'"), command[ic]));
								int id = hd.dt.FindClosestQTFHead(std::complex<double>(head1, head2));
								BEM::Print("\n" + Format(t_("Closest heading to %f:%f is %f:%f"), head1, head2, hd.dt.qhead[id].real(), hd.dt.qhead[id].imag()));
								FindAdd(ids, id);
							}
							UVector<int> idsDel;
							for (int i = 0; i < hd.dt.qhead.size(); ++i)
								if (Find(ids, i) < 0)
									idsDel << i;
							hd.DeleteHeadingsQTF(idsDel);
							BEM::Print("\n" + S(t_("QTF headings deleted")));
						} else if (param == "-delbutqtfheadid") {
							if (bem.hydros.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							Hydro &hd = bem.hydros[bemid];	
							UVector<int> ids;
							int id;
							while (ic+1 < command.size() && !IsNull(id = ScanInt(command[ic+1]))) {
								if (id < 0 || id >= hd.dt.qhead.size())
									throw Exc(Format(t_("Wrong head id '%s'"), command[ic+1]));
								FindAdd(ids, id);
								ic++;
							}
							UVector<int> idsDel;
							for (int i = 0; i < hd.dt.qhead.size(); ++i)
								if (Find(ids, i) < 0)
									idsDel << i;
							hd.DeleteHeadingsQTF(idsDel);
							BEM::Print("\n" + S(t_("QTF headings deleted")));	
						} else if (param == "-params") {
							CheckIfAvailableArg(command, ic+1, "-params");
							
							while (command.size() > ic+1 && !command[ic+1].StartsWith("-")) {
								++ic;
								if (ToLower(command[ic]) == "length") {
									CheckIfAvailableArg(command, ++ic, "-p length");
									double len = ScanDouble(command[ic]);
									if (IsNull(len))
										throw Exc(Format(t_("Wrong argument '%s'"), command[ic]));
									bem.len = len;
									BEM::Print("\n" + Format(t_("length is %f"), len));	
								} else if (ToLower(command[ic]) == "depth") {
									CheckIfAvailableArg(command, ++ic, "-p depth");
									double depth = ScanDouble(command[ic]);
									if (IsNull(depth))
										throw Exc(Format(t_("Wrong argument '%s'"), command[ic]));
									bem.depth = depth;
									BEM::Print("\n" + Format(t_("depth is %f"), depth));	
								} else 
									throw Exc(Format(t_("Wrong command '%s'"), command[ic]));
							}
						} else if (param == "-p" || param == "-print") {
							if (bem.hydros.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							Hydro &hy = bem.hydros[bemid];
							while (command.size() > ic+1 && !command[ic+1].StartsWith("-")) {
								ic++;
								String pparam = ToLower(command[ic]);
								if (pparam == "nb") {
									Cout() << "\n";
									BEM::Print(t_("Nb:") + S(" ")); 
									lastPrint = FormatInt(hy.dt.Nb);
									Cout() << lastPrint;
								} else if (pparam == "nh") {
									Cout() << "\n";
									BEM::Print(t_("Nh:") + S(" ")); 
									lastPrint = FormatInt(hy.dt.Nh);
									Cout() << lastPrint;
								} else if (pparam == "nf") {
									Cout() << "\n";
									BEM::Print(t_("Nf:") + S(" ")); 
									lastPrint = FormatInt(hy.dt.Nf);
									Cout() << lastPrint;
								} else if (pparam == "w" || pparam == "frequencies") {
									Cout() << "\n";
									BEM::Print(t_("w:") + S(" ")); 
									lastPrint.Clear();
									for (double d : hy.dt.w)
										lastPrint << d << " "; 
									Cout() << lastPrint;
								} else if (pparam == "headings") {
									Cout() << "\n";
									BEM::Print(t_("head:") + S(" ")); 
									lastPrint.Clear();
									for (double d : hy.dt.head)
										lastPrint << d << " "; 
									Cout() << lastPrint; 
								} else if (pparam == "a") {
									CheckIfAvailableArg(command, ++ic, "A #dof 1");									
									int idf = ScanInt(command[ic]);
									if (idf < 1 || idf > 6*hy.dt.Nb)
										throw Exc(Format(t_("Wrong dof in '%s'"), command[ic]));
									CheckIfAvailableArg(command, ++ic, "A #dof 2");									
									int jdf = ScanInt(command[ic]);
									if (jdf < 1 || jdf > 6*hy.dt.Nb)
										throw Exc(Format(t_("Wrong dof in '%s'"), command[ic]));
									Cout() << "\n";
									BEM::Print(Format(t_("A(%d,%d):"), idf, jdf) + S(" "));
									idf--;	jdf--;
									lastPrint.Clear();
									for (int ifr = 0; ifr < hy.dt.Nf; ++ifr) 
										lastPrint << Format("%f ", hy.A_dim(ifr, idf, jdf));
									Cout() << lastPrint;
								} else if (pparam == "ainf") {
									Cout() << "\n";
									BEM::Print(t_("Ainf:") + S(" "));
									lastPrint.Clear();
									for (int idf = 0; idf < 6*hy.dt.Nb; ++idf) 
										for (int jdf = 0; jdf < 6*hy.dt.Nb; ++jdf) 
											lastPrint << Format("%f ", hy.Ainf_dim(idf, jdf));
									Cout() << lastPrint;
								} else if (pparam == "b") {
									CheckIfAvailableArg(command, ++ic, "B #dof 1");									
									int idf = ScanInt(command[ic]);
									if (idf < 1 || idf > 6*hy.dt.Nb)
										throw Exc(Format(t_("Wrong dof in '%s'"), command[ic]));
									idf--;
									CheckIfAvailableArg(command, ++ic, "B #dof 2");									
									int jdf = ScanInt(command[ic]);
									if (jdf < 1 || jdf > 6*hy.dt.Nb)
										throw Exc(Format(t_("Wrong dof in '%s'"), command[ic]));
									jdf--;
									Cout() << "\n";
									BEM::Print(t_("B:") + S(" "));
									lastPrint.Clear();
									for (int ifr = 0; ifr < hy.dt.Nf; ++ifr) 
										lastPrint << Format("%f ", hy.B_dim(ifr, idf, jdf));
									Cout() << lastPrint;
								} else if (pparam == "theave") {
									CheckIfAvailableArg(command, ++ic, "Id. body");
									int ib = ScanInt(command[ic]);
									if (ib < 1 || ib > hy.dt.Nb)
										throw Exc(Format(t_("Wrong body id. in '%s'"), command[ic]));
									Cout() << "\n";
									BEM::Print(Format(t_("Theave(%d):"), ib) + " "); 
									lastPrint = Format("%f", hy.Theave(ib-1));
									Cout() << lastPrint;
								} else if (pparam == "troll") {
									CheckIfAvailableArg(command, ++ic, "Id. body");
									int ib = ScanInt(command[ic]);
									if (ib < 1 || ib > hy.dt.Nb)
										throw Exc(Format(t_("Wrong body id. in '%s'"), command[ic]));
									Cout() << "\n";
									BEM::Print(Format(t_("Troll(%d):"), ib) + " "); 
									lastPrint = Format("%f", hy.Troll(ib-1));
									Cout() << lastPrint;
								} else if (pparam == "tpitch") {
									CheckIfAvailableArg(command, ++ic, "Id. body");
									int ib = ScanInt(command[ic]);
									if (ib < 1 || ib > hy.dt.Nb)
										throw Exc(Format(t_("Wrong body id. in '%s'"), command[ic]));
									Cout() << "\n";
									BEM::Print(Format(t_("Tpitch(%d):"), ib) + " "); 
									lastPrint = Format("%f", hy.Tpitch(ib-1));
									Cout() << lastPrint;
								} else if (pparam == "gmroll") {
									CheckIfAvailableArg(command, ++ic, "Id. body");
									int ib = ScanInt(command[ic]);
									if (ib < 1 || ib > hy.dt.Nb)
										throw Exc(Format(t_("Wrong body id. in '%s'"), command[ic]));
									Cout() << "\n";
									BEM::Print(Format(t_("GMroll(%d):"), ib) + " "); 
									lastPrint = Format("%f", hy.GMroll(ib-1));
									Cout() << lastPrint;
								} else if (pparam == "gmpitch") {
									CheckIfAvailableArg(command, ++ic, "Id. body");
									int ib = ScanInt(command[ic]);
									if (ib < 1 || ib > hy.dt.Nb)
										throw Exc(Format(t_("Wrong body id. in '%s'"), command[ic]));
									Cout() << "\n";
									BEM::Print(Format(t_("GMpitch(%d):"), ib) + " "); 
									lastPrint = Format("%f", hy.GMpitch(ib-1));
									Cout() << lastPrint;
								} else
									throw Exc(Format(t_("Unknown argument '%s'"), command[ic]));
							}
						} else 
							throw Exc(Format(t_("Unknown argument '%s'"), command[ic]));
					} else if (nextcommands == "mesh") {
						if (param == "-i" || param == "-input") {
							CheckIfAvailableArg(command, ++ic, "--input");
							
							String file = FileName(command[ic]);
							if (!FileExists(file)) 
								throw Exc(Format(t_("File '%s' not found"), file)); 
							
							bem.LoadBody(file, echo ? Status : NoPrint, false, false);		// Doesn't work for multibody .dat
							meshid = bem.surfs.size() - 1;
							BEM::Print("\n" + Format(t_("File '%s' loaded"), file));
						} else if (param == "-r" || param == "-report") {
							if (bem.surfs.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							bem.surfs[meshid].Report(bem.rho);
						} else if (param == "-cl" || param == "-clear") {
							bem.surfs.Clear();
							meshid = -1;
							BEM::Print("\n" + S(t_("Body data cleared")));	
						} else if (param == "-setid") {
							if (bem.surfs.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							CheckIfAvailableArg(command, ++ic, "-setid");

							meshid = ScanInt(command[ic]);
							if (IsNull(meshid) || meshid < 0 || meshid > bem.surfs.size()-1)
								throw Exc(Format(t_("Invalid id %s"), command[ic]));
							BEM::Print("\n" + Format(t_("Body active model id is %d"), bemid));	
						} else if (param == "-c" || param == "-convert") {
							if (bem.surfs.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							CheckIfAvailableArg(command, ++ic, "-convert");
							
							String file = FileName(command[ic]);
							
							bool symX = false, symY = false;
							Body::MESH_FMT meshFmt = Body::UNKNOWN;
							while (command.size() > ic+1 && !command[ic+1].StartsWith("-")) {
								ic++;
								String pparam = ToLower(command[ic]);
								if (pparam == "symx")
									symX = true;
								else if (pparam == "symy")
									symY = true;
								else if (pparam == "symxy")
									symX = symY = true;
								else if (Body::GetCodeBodyStr(pparam) != Body::UNKNOWN) {
									meshFmt = Body::GetCodeBodyStr(pparam);
									if (!Body::meshCanSave[meshFmt])
										throw Exc(Format(t_("Saving format '%s' is not implemented"), pparam));									
								} else
									throw Exc(Format(t_("Unknown argument '%s'"), command[ic]));
							}
							Body::SaveAs(bem.surfs[meshid], file, meshFmt, Body::ALL, bem.rho, bem.g, symX, symY);
							BEM::Print("\n" + Format(t_("Model id %d saved as '%s'"), meshid, file));
						} else if (param == "-t" || param == "-translate") {
							if (bem.surfs.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							CheckIfAvailableArg(command, ++ic, "x");
							double x = ScanDouble(command[ic]);
							CheckIfAvailableArg(command, ++ic, "y");
							double y = ScanDouble(command[ic]);
							CheckIfAvailableArg(command, ++ic, "z");
							double z = ScanDouble(command[ic]);
							Body &msh = bem.surfs[meshid];
							msh.dt.mesh.Translate(x, y, z);
							msh.dt.cg.Translate(x, y, z);
							msh.AfterLoad(bem.rho, bem.g, false, false);	
							BEM::Print("\n" + Format(t_("Body id %d translated %f, %f, %f"), meshid, x, y, z)); 
						} else if (param == "-rot" || param == "-rotate") {
							if (bem.surfs.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							CheckIfAvailableArg(command, ++ic, "ax");
							double ax = ScanDouble(command[ic]);
							CheckIfAvailableArg(command, ++ic, "ay");
							double ay = ScanDouble(command[ic]);
							CheckIfAvailableArg(command, ++ic, "az");
							double az = ScanDouble(command[ic]);
							CheckIfAvailableArg(command, ++ic, "cx");
							double cx = ScanDouble(command[ic]);
							CheckIfAvailableArg(command, ++ic, "cy");
							double cy = ScanDouble(command[ic]);
							CheckIfAvailableArg(command, ++ic, "cz");
							double cz = ScanDouble(command[ic]);
							Body &msh = bem.surfs[meshid];
							msh.dt.mesh.Rotate(ToRad(ax), ToRad(ay), ToRad(az), cx, cy, cz);	
							msh.dt.cg.Rotate(ToRad(ax), ToRad(ay), ToRad(az), cx, cy, cz);
							msh.AfterLoad(bem.rho, bem.g, false, false);	
							BEM::Print("\n" + Format(t_("Body id %d rotated angles %f, %f, %f around centre %f, %f, %f"), meshid, ax, ay, az, cx, cy, cz));
						} else if (param == "-cg") {
							if (bem.surfs.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							CheckIfAvailableArg(command, ++ic, "cgx");
							double x = ScanDouble(command[ic]);
							CheckIfAvailableArg(command, ++ic, "cgy");
							double y = ScanDouble(command[ic]);
							CheckIfAvailableArg(command, ++ic, "cgz");
							double z = ScanDouble(command[ic]);
							Body &msh = bem.surfs[meshid];
							msh.dt.cg.x = x;
							msh.dt.cg.y = y;
							msh.dt.cg.z = z;
							msh.AfterLoad(bem.rho, bem.g, true, false);
							BEM::Print("\n" + Format(t_("CG is %f, %f, %f"), x, y, z));
						} else if (param == "-c0") {
							if (bem.surfs.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							CheckIfAvailableArg(command, ++ic, "c0x");
							double x = ScanDouble(command[ic]);
							CheckIfAvailableArg(command, ++ic, "c0y");
							double y = ScanDouble(command[ic]);
							CheckIfAvailableArg(command, ++ic, "c0z");
							double z = ScanDouble(command[ic]);
							Body &msh = bem.surfs[meshid];
							msh.dt.c0.x = x;
							msh.dt.c0.y = y;
							msh.dt.c0.z = z;
							msh.AfterLoad(bem.rho, bem.g, true, false);
							BEM::Print("\n" + Format(t_("CG is %f, %f, %f"), x, y, z));
						} else if (param == "-mass") { 
							if (bem.surfs.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							CheckIfAvailableArg(command, ++ic, "mass");
							double mass = ScanDouble(command[ic]);
							Body &msh = bem.surfs[meshid];
							msh.SetMass(mass);
							msh.AfterLoad(bem.rho, bem.g, true, false);
							BEM::Print("\n" + Format(t_("Mass is %f"), mass));
						} else if (param == "-reset") {	
							if (bem.surfs.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							bem.surfs[meshid].Reset(bem.rho, bem.g);
							BEM::Print("\n" + Format(t_("Body id %d position is reset"), meshid));
						} else if (param == "-getwaterplane") {
							if (bem.surfs.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							bem.AddWaterSurface(meshid, 'e');
							meshid = bem.surfs.size() - 1;
							BEM::Print("\n" + Format(t_("Body id %d waterplane is got"), meshid));
						} else if (param == "-gethull") {
							if (bem.surfs.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							bem.AddWaterSurface(meshid, 'r');
							meshid = bem.surfs.size() - 1;
							BEM::Print("\n" + Format(t_("Body id %d hull is got"), meshid));
						} else if (param == "-p" || param == "-print") {
							if (bem.surfs.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							Body &msh = bem.surfs[meshid];
							while (command.size() > ic+1 && !command[ic+1].StartsWith("-")) {
								ic++;
								String pparam = ToLower(command[ic]);
								if (pparam == "volume") {
									Cout() << "\n";
									BEM::Print(t_("Volume:") + S(" ")); 
									lastPrint = Format("%f %f %f", msh.dt.mesh.volumex, msh.dt.mesh.volumey, msh.dt.mesh.volumez);
									Cout() << lastPrint;
								} else if (pparam == "underwatervolume") {
									Cout() << "\n";
									BEM::Print(t_("UnderwaterVolume:") + S(" ")); 
									lastPrint = Format("%f %f %f", msh.dt.under.volumex, msh.dt.under.volumey, msh.dt.under.volumez);
									Cout() << lastPrint;
								} else if (pparam == "surface") {
									Cout() << "\n";
									BEM::Print(t_("Surface:") + S(" ")); 
									lastPrint = Format("%f", msh.dt.mesh.surface);
									Cout() << lastPrint;
								} else if (pparam == "underwatersurface") {
									Cout() << "\n";
									BEM::Print(t_("UnderwaterSurface:") + S(" ")); 
									lastPrint = Format("%f", msh.dt.under.surface);
									Cout() << lastPrint;
								} else if (pparam == "cb") {
									Cout() << "\n";
									BEM::Print(t_("CB:") + S(" ")); 
									lastPrint = Format("%f %f %f", msh.dt.cb.x, msh.dt.cb.y, msh.dt.cb.z);
									Cout() << lastPrint;
								} else if (pparam == "cg_vol") {
									Cout() << "\n";
									BEM::Print(t_("CG_vol:") + S(" ")); 
									Point3D cg = msh.dt.mesh.GetCentreOfBuoyancy();
									lastPrint = Format("%f %f %f", cg.x, cg.y, cg.z);
									Cout() << lastPrint;
								} else if (pparam == "cg_surf") {
									Cout() << "\n";
									BEM::Print(t_("CG_surf:") + S(" ")); 
									Point3D cg = msh.dt.mesh.GetCentreOfGravity_Surface();
									lastPrint = Format("%f %f %f", cg.x, cg.y, cg.z);
									Cout() << lastPrint;
								} else if (pparam == "hydrostiffness") {
									Cout() << "\n";
									BEM::Print(t_("HydrostaticStiffness:") + S(" "));
									lastPrint.Clear();
									for (int i = 0; i < 6; ++i) {
										for (int j = 0; j < 6; ++j) {
											if (!Hydro::C_units(i, j).IsEmpty()) 
												lastPrint << msh.dt.C(i, j) << " ";
										}
									}
									int idColor = msh.dt.under.VolumeMatch(bem.volWarning/100., bem.volError/100.);
									if (idColor == -1)
										lastPrint << ". " << t_("Body warning_ Maybe incomplete");
									else if (idColor == -2)
										lastPrint << ". " << t_("Body error: Probably incomplete");
									Cout() << lastPrint;
								} else if (pparam == "hydrostatic_force") {
									Cout() << "\n";
									BEM::Print(t_("HydrostaticForce:") + S(" "));
									lastPrint.Clear();
									Force6D f = msh.dt.under.GetHydrostaticForce(msh.dt.c0, bem.rho, bem.g);
									for (int i = 0; i < 6; ++i) 				
										lastPrint << f[i] << " ";
									int idColor = msh.dt.under.VolumeMatch(bem.volWarning/100., bem.volError/100.);
									if (idColor == -1)
										lastPrint << ". " << t_("Body warning: Maybe incomplete");
									else if (idColor == -2)
										lastPrint << ". " << t_("Body error: Probably incomplete");
									Cout() << lastPrint;
								} else if (pparam.StartsWith("inertia")) {
									Cout() << "\n";
									BEM::Print(t_("Inertia:") + S(" "));
									lastPrint.Clear();
									Eigen::Matrix3d inertia;
									Point3D centre;
									CheckIfAvailableArg(command, ++ic, "Inertia cx");
									centre.x = ScanDouble(command[ic]);
									CheckIfAvailableArg(command, ++ic, "Inertia cy");
									centre.y = ScanDouble(command[ic]);
									CheckIfAvailableArg(command, ++ic, "Inertia cz");
									centre.z = ScanDouble(command[ic]);
									bool volume = false;
									if (pparam == "inertia" || pparam == "inertia_vol")
										volume = true;
									
									msh.dt.mesh.GetInertia33_Radii(inertia, centre, volume, false);
									for (int i = 0; i < 3; ++i) 
										for (int j = 0; j < 3; ++j) 
											lastPrint << inertia(i, j) << " ";
									Cout() << lastPrint;
								} else if (pparam == "gz") {
									Cout() << "\n";
									BEM::Print(t_("GZ:") + S(" "));
									lastPrint.Clear();
									
									CheckIfAvailableArg(command, ++ic, "Angle");
									double angle = ScanDouble(command[ic]);
									CheckIfAvailableArg(command, ++ic, "From");
									double from = ScanDouble(command[ic]);
									CheckIfAvailableArg(command, ++ic, "To");
									double to = ScanDouble(command[ic]);
									CheckIfAvailableArg(command, ++ic, "Delta");
									double delta = ScanDouble(command[ic]);
									UVector<double> dataangle, datagz;
									double tolerance = 0.1;
									String errors;
									msh.GZ(from, to, delta, angle, bem.rho, bem.g, tolerance, dataangle, datagz, errors);									
									for (int i = 0; i < dataangle.size(); ++i) 
										lastPrint << dataangle[i] << " ";
									lastPrint << "\n";
									for (int i = 0; i < datagz.size(); ++i) 
										lastPrint << datagz[i] << " ";
									Cout() << lastPrint;
								} else if (pparam == "gm") {
									Cout() << "\n";
									BEM::Print(t_("gm:") + S(" ")); 
									lastPrint = Format("%f %f", msh.GMpitch(bem.rho, bem.g), msh.GMroll(bem.rho, bem.g));
									Cout() << lastPrint;
								} else
									throw Exc(Format(t_("Unknown argument '%s'"), command[ic]));
							}
						} else 
							throw Exc(Format(t_("Unknown argument '%s'"), command[ic]));
					} else if (nextcommands == "time") {
						if (param == "-i" || param == "-input") {
							CheckIfAvailableArg(command, ++ic, "-input");

							String file = FileName(command[ic]);
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
							
							CheckIfAvailableArg(command, ++ic, "-convert");
							
							String file = FileName(command[ic]);
							
							if (fast.Save(file, "", ScatterDraw::GetDefaultCSVSeparator())) 
								BEM::Print("\n" + Format(t_("Results saved as '%s'"), file));
						} else if (param == "-p" || param == "-print") {
							while (command.size() > ic+1 && !command[ic+1].StartsWith("-")) {
								ic++;
								String pparam = ToLower(command[ic]);
								if (pparam == "list") {
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
									UVector<int> p = fast.FindParameterMatch(pparam);
									if (p.IsEmpty())
										throw Exc(Format(t_("Parameter '%s' not found"), pparam));

									int id = p[0];
									ic++;
									double time = ScanDouble(command[ic]);
									if (!IsNull(time))
										lastPrint = FormatDouble(fast.GetVal(time, id));
									else {
										VectorXd d = fast.GetVector(id);
										CleanNAN_safe(d);		// Just in case. It do happens
										if (command[ic] == "avg" || command[ic] == "mean") 
											lastPrint = FormatDouble(d.mean());
										else if (command[ic] == "max") 
											lastPrint = FormatDouble(d.maxCoeff());
										else if (command[ic] == "min") 
											lastPrint = FormatDouble(d.minCoeff());
										else  
											throw Exc(Format(t_("Parameter '%s' not found"), command[ic]));
									}
									Cout() << lastPrint;
								}
							}
						} else
							throw Exc(Format(t_("Unknown argument '%s'"), command[ic]));
					} else if (nextcommands == "wind") {
						if (param == "-i" || param == "-input") {
							CheckIfAvailableArg(command, ++ic, "--input");
							
							String file = FileName(command[ic]);
							if (!FileExists(file)) 
								throw Exc(Format(t_("File '%s' not found"), file)); 
							
							Wind &w = wind.Add();
							String ret;
							if (!IsEmpty(ret = w.Load(file)))
								throw Exc(ret);
							
							windid = wind.size() - 1;
							BEM::Print("\n" + Format(t_("File '%s' loaded"), file));
						} else if (param == "-params") {
							CheckIfAvailableArg(command, ic+1, "-params");
							
							while (command.size() > ic+1 && !command[ic+1].StartsWith("-")) {
								++ic;
								if (ToLower(command[ic]) == "hubheight") {
									if (wind.IsEmpty()) 
										throw Exc(t_("No file loaded"));
									CheckIfAvailableArg(command, ++ic, "hubheight");
									double h = ScanDouble(command[ic]);
									if (IsNull(h))
										throw Exc(Format(t_("Wrong argument '%s'"), command[ic]));
									wind[windid].SetHubHeight((float)h);
									BEM::Print("\n" + Format(t_("Hub height is %f"), h));	
								} else if (ToLower(command[ic]) == "gridheight") {
									if (wind.IsEmpty()) 
										throw Exc(t_("No file loaded"));
									CheckIfAvailableArg(command, ++ic, "gridheight");
									double h = ScanDouble(command[ic]);
									if (IsNull(h))
										throw Exc(Format(t_("Wrong argument '%s'"), command[ic]));
									wind[windid].SetGridHeight((float)h);
									BEM::Print("\n" + Format(t_("Grid height is %f"), h));	
								} else if (ToLower(command[ic]) == "ti") {	
									if (wind.IsEmpty()) 
										throw Exc(t_("No file loaded"));
									CheckIfAvailableArg(command, ++ic, "TI u");
									double u = ScanDouble(command[ic]);
									if (IsNull(u))
										throw Exc(Format(t_("Wrong argument '%s'"), command[ic]));
									wind[windid].SetTI_u((float)u);
									double v = 0.8*u;
									double w = 0.5*u;
									if (ic+1 < command.size()) {
										double num = ScanDouble(command[ic+1]);
										if (!IsNull(num)) {
											ic++;
											v = num;
											if (ic+1 < command.size()) {		
												num = ScanDouble(command[ic+1]);	
												if (!IsNull(num)) {
													ic++;
													w = num;	
												}
											}
										}
									}
									wind[windid].SetTI_v((float)v);
									wind[windid].SetTI_w((float)w);
									BEM::Print("\n" + Format(t_("Turbulence intensity is u: %f %%, v: %f %%, w: %f %%"), u, v, w));	
								} else if (ToLower(command[ic]) == "scale") {	
									if (wind.IsEmpty()) 
										throw Exc(t_("No file loaded"));
									CheckIfAvailableArg(command, ++ic, "scale u");
									double u = ScanDouble(command[ic]);
									if (IsNull(u))
										throw Exc(Format(t_("Wrong argument '%s'"), command[ic]));
									double v = 1;
									double w = 1;
									if (ic+1 < command.size()) {
										double num = ScanDouble(command[ic+1]);
										if (!IsNull(num)) {
											ic++;
											v = num;
											if (ic+1 < command.size()) {		
												num = ScanDouble(command[ic+1]);	
												if (!IsNull(num)) {
													ic++;
													w = num;	
												}
											}
										}
									}
									wind[windid].SetFactor((float)u, (float)v, (float)w);
									BEM::Print("\n" + Format(t_("Velocities scaled by u: %f %%, v: %f %%, w: %f %%"), u, v, w));	
								} else if (ToLower(command[ic]) == "powerlaw") {	
									if (wind.IsEmpty()) 
										throw Exc(t_("No file loaded"));	
									CheckIfAvailableArg(command, ++ic, "powerLaw");
									double pl = ScanDouble(command[ic]);
									if (IsNull(pl))
										throw Exc(Format(t_("Wrong argument '%s'"), command[ic]));
									CheckIfAvailableArg(command, ++ic, "hub height for power law");
									double zh = ScanDouble(command[ic]);
									if (IsNull(zh))
										throw Exc(Format(t_("Wrong argument '%s'"), command[ic]));
									wind[windid].SetPowerLaw((float)pl, (float)zh);
									BEM::Print("\n" + Format(t_("Power law is %f set at height %f m"), pl, zh));
								} else if (ToLower(command[ic]) == "periodic") {
									if (wind.IsEmpty()) 
										throw Exc(t_("No file loaded"));	
									CheckIfAvailableArg(command, ++ic, "-p periodic");
									String s = ToLower(command[ic]);
									if (s != "true" && s != "false")
										throw Exc(Format(t_("Wrong argument '%s'"), command[ic]));
									wind[windid].SetPeriodic(s == "true" ? true : false);
									BEM::Print("\n" + Format(t_("Set periodic %s"), s));	
								} else 
									throw Exc(Format(t_("Wrong command '%s'"), command[ic]));
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
							
							CheckIfAvailableArg(command, ++ic, "-setid");

							windid = ScanInt(command[ic]);
							if (IsNull(windid) || windid < 0 || windid > wind.size()-1)
								throw Exc(Format(t_("Invalid id %s"), command[ic]));
							BEM::Print("\n" + Format(t_("Wind active model id is %d"), windid));	
						} else if (param == "-c" || param == "-convert") {
							if (wind.IsEmpty()) 
								throw Exc(t_("No file loaded"));
							
							CheckIfAvailableArg(command, ++ic, "-convert");
							
							String file = FileName(command[ic]);
							
							String ret;
							if (!IsEmpty(ret = wind[windid].Save(file)))
								throw Exc(ret);
							BEM::Print("\n" + Format(t_("Model id %d saved as '%s'"), windid, file));
						} else if (param == "-p" || param == "-print") {
							Wind &w = wind[windid];
							while (command.size() > ic+1 && !command[ic+1].StartsWith("-")) {
								ic++;
								String pparam = ToLower(command[ic]);
								if (pparam == "vel") {
									BEM::Print("\n");
									BEM::Print(t_("Velocity:") + S(" "));
									lastPrint.Clear();
									
									CheckIfAvailableArg(command, ++ic, "Position Y");
									double ypos = ScanDouble(command[ic]);
									CheckIfAvailableArg(command, ++ic, "Position Z");
									double zpos = ScanDouble(command[ic]);
									int idz, idy;
									w.GetPos(zpos, ypos, idz, idy);
									Eigen::VectorXd v;
									v = w.GetNorm(idz, idy);
									
									CheckIfAvailableArg(command, ++ic, "<time>/avg/data");
									String what = command[ic];
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
								} else if (pparam == "velcomp") {
									BEM::Print("\n");
									BEM::Print(t_("Velocity component:") + S(" "));
									lastPrint.Clear();
									
									CheckIfAvailableArg(command, ++ic, "Component");
									int icomp = ScanInt(command[ic]);
									CheckIfAvailableArg(command, ++ic, "Position Y");
									double ypos = ScanDouble(command[ic]);
									CheckIfAvailableArg(command, ++ic, "Position Z");
									double zpos = ScanDouble(command[ic]);
									int idz, idy;
									w.GetPos(zpos, ypos, idz, idy);
									Eigen::VectorXd v;
									v = w.Get(icomp, idz, idy);
									
									CheckIfAvailableArg(command, ++ic, "<time>/avg/data");
									String what = command[ic];
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
								} else if (pparam == "time") {
									BEM::Print("\n");
									BEM::Print(t_("Time:") + S(" "));
									lastPrint.Clear();
									Eigen::VectorXd t;
									t = w.GetTime();
									for (int i = 0; i < t.size(); ++i) 
										lastPrint << t(i) << " ";
									BEM::Print(lastPrint);
								} else
									throw Exc(Format(t_("Unknown parameter '%s'"), pparam));
							}
						} else 
							throw Exc(Format(t_("Unknown argument '%s'"), command[ic]));
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
									
							CheckIfAvailableArg(command, ++ic, "-runwave from");
							String from = command[ic];
							CheckIfAvailableArg(command, ++ic, "-runwave to");
							String to = command[ic];
							
							int numTry = numOrcaTries;
							//String errorStr;
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
									
							CheckIfAvailableArg(command, ++ic, "-runflex from");
							String from = command[ic];
							CheckIfAvailableArg(command, ++ic, "-runflex to");
							String to = command[ic];
							
							int numTry = numOrcaTries;
							//String errorStr;
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
									
							CheckIfAvailableArg(command, ++ic, "-loadSim sim");
							String from = command[ic];
							
							int numTry = numOrcaTries;
							//String errorStr;
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
							
							CheckIfAvailableArg(command, ++ic, "-numthread num");
							int numth = ScanInt(command[ic]);
							if (IsNull(numth))
								throw Exc(Format(t_("Invalid thread number %s"), command[ic]));							

							orca.SetThreadCount(numth);
							
							BEM::Print("\n" + Format(t_("Thread number set to %d"), numth));
						} else if (param == "-numtries") {
							if (!dllOrcaLoaded) {
								if (orca.FindInit()) 
									dllOrcaLoaded = true;		
							}
							if (!dllOrcaLoaded)
								throw Exc(t_("DLL not found"));
							
							CheckIfAvailableArg(command, ++ic, "-numtries num");
							numOrcaTries = ScanInt(command[ic]);
							if (IsNull(numOrcaTries))
								throw Exc(Format(t_("Invalid thread number %s"), command[ic]));							

							BEM::Print("\n" + Format(t_("Number of OrcaFlex license tries set to %d"), numOrcaTries));
						} else if (param == "-timelog") {
							if (!dllOrcaLoaded) {
								if (orca.FindInit()) 
									dllOrcaLoaded = true;		
							}
							if (!dllOrcaLoaded)
								throw Exc(t_("DLL not found"));
							
							CheckIfAvailableArg(command, ++ic, "-timelog num");
							int timeLog = int(StringToSeconds(command[ic]));
							if (IsNull(timeLog))
								throw Exc(Format(t_("Invalid time between logs %s"), command[ic]));							

							BEM::Print("\n" + Format(t_("OrcaFlex time between logs set to %.1f sec"), timeLog));							
							
							orca.deltaLogSimulation = timeLog;
						} else if (param == "-lp" || param == "-loadparam") {
							CheckIfAvailableArg(command, ++ic, "-loadParam");
							
							paramList << command[ic];
							centreList << centreOrca;
						} else if (param == "-c" || param == "-convert") {
							if (!dllOrcaLoaded) {
								if (orca.FindInit()) 
									dllOrcaLoaded = true;		
							}
							if (!dllOrcaLoaded)
								throw Exc(t_("DLL not found"));
							
							CheckIfAvailableArg(command, ++ic, "-convert");
							
							String file = FileName(command[ic]);

							orca.SaveCsv(file, paramList, centreList, ScatterDraw::GetDefaultCSVSeparator());
							BEM::Print("\n" + Format(t_("Results saved as '%s'"), file));
						} else if (param == "-rs" || param == "-refsystem") {	
							CheckIfAvailableArg(command, ++ic, "Reference cx");
							centreOrca.x = ScanDouble(command[ic]);
							CheckIfAvailableArg(command, ++ic, "Reference cy");
							centreOrca.y = ScanDouble(command[ic]);
							CheckIfAvailableArg(command, ++ic, "Reference cz");
							centreOrca.z = ScanDouble(command[ic]);
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
							
							while (command.size() > ic+1 && !command[ic+1].StartsWith("-")) {
								ic++;
								String pparam = ToLower(command[ic]);
								if (pparam == "list") {
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
								} else if (pparam == "numthread") {
									int numth = orca.GetThreadCount();
									BEM::Print("\n" + Format("%d", numth));
								} else {
									Cout() << "\n";

									String unit;
									int objType;
									VectorXd d;
									orca.GetFlexSimVar(command[ic], centreOrca, unit, objType, d);
									
									ic++;
									if (command.size() <= ic)
										throw Exc(Format(t_("Last command '%s' is not complete"), pparam));
									if (command[ic] == "data") {
										lastPrint = "";
										for (int i = 0; i < d.size(); ++i) {
											if (i > 0)
												lastPrint << " ";
											lastPrint << d[i];
										}
									} else if (command[ic] == "avg" || command[ic] == "mean") 
										lastPrint = FormatDouble(d.mean());
									else if (command[ic] == "max") 
										lastPrint = FormatDouble(d.maxCoeff());
									else if (command[ic] == "min") 
										lastPrint = FormatDouble(d.minCoeff());
									else  
										throw Exc(Format(t_("Parameter '%s' not found in -p"), command[ic]));
									Cout() << lastPrint;
								}
							}
						} else if (param == "-dll") {
							CheckIfAvailableArg(command, ++ic, "-dll");

							String file = FileName(command[ic]);
							if (!FileExists(file)) {
								file = AFX(file, "OrcFxAPI.dll");
								if (!FileExists(file)) 	
									throw Exc(Format(t_("File '%s' not found"), file)); 
							}
							if (!orca.Init(file))
								throw Exc(Format(t_("DLL '%s' not found"), command[ic]));
							dllOrcaLoaded = true;
							BEM::Print("\n" + Format(t_("Orca .dll '%s' loaded"), orca.GetDLLPath()));
						} else 
							throw Exc(Format(t_("Unknown argument '%s'"), command[ic]));
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
