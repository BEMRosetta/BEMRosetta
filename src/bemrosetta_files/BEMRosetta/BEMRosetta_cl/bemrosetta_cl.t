#ifdef _MSC_VER
#pragma setlocale("C")
#endif
// BEMRosetta.cpp

T_("surge")
esES("")
frFR("")

T_("sway")
esES("")
frFR("")

T_("heave")
esES("")
frFR("")

T_("roll")
esES("")
frFR("")

T_("pitch")
esES("")
frFR("")

T_("yaw")
esES("")
frFR("")

T_("s")
esES("")
frFR("")

T_("w")
esES("")
frFR("")

T_("h")
esES("")
frFR("")

T_("r")
esES("")
frFR("")

T_("p")
esES("")
frFR("")

T_("y")
esES("")
frFR("")

T_("%s is not the same %f<>%f")
esES("")
frFR("")

T_("Density rho")
esES("")
frFR("")

T_("Gravity g")
esES("")
frFR("")

T_("Water depth h")
esES("")
frFR("")

T_("%s is not the same %d<>%d")
esES("")
frFR("")

T_("Number of bodies")
esES("")
frFR("")

T_("Number of frequencies")
esES("")
frFR("")

T_("#%d %s")
esES("")
frFR("")

T_("frequency")
esES("")
frFR("")

T_("Number of headings")
esES("")
frFR("")

T_("%s[%d](%d, %d)")
esES("")
frFR("")

T_("A")
esES("")
frFR("")

T_("B")
esES("")
frFR("")

T_("C")
esES("")
frFR("")

T_("%s(%d, %d)")
esES("")
frFR("")

T_("cg")
esES("")
frFR("")

T_("Conversion to type of file '%s' not supported")
esES("")
frFR("")

T_("Joined files")
esES("")
frFR("")

T_("Water depth does not match between '%s'(%d) and '%s'(%d)")
esES("")
frFR("")

T_("No water depth found in models")
esES("")
frFR("")

T_("Number of bodies does not match between '%s'(%d) and '%s'(%d)")
esES("")
frFR("")

T_("No body found in models")
esES("")
frFR("")

T_("No head found in models")
esES("")
frFR("")

T_("No frequency found in models")
esES("")
frFR("")

T_("%s file '%s'")
esES("")
frFR("")

T_("g [m/s2]: %s, h [m]: %s, rho [kg/m3]: %s, length scale [m]: %s")
esES("")
frFR("")

T_("NONE")
esES("")
frFR("")

T_("delta %s [rad/s]")
esES("")
frFR("")

T_("non constant delta (%s)")
esES("")
frFR("")

T_("%s to %s %s")
esES("")
frFR("")

T_("%s [rad/s]")
esES("")
frFR("")

T_("delta %.1f [deg]")
esES("")
frFR("")

T_("%.1f to %.1f %s")
esES("")
frFR("")

T_("%.1f [deg]")
esES("")
frFR("")

T_("#freqs: %d (%s)")
esES("")
frFR("")

T_("#headings: %d (%s)")
esES("")
frFR("")

T_("#bodies: %d")
esES("")
frFR("")

T_("dof")
esES("")
frFR("")

T_("vol [m3]")
esES("")
frFR("")

T_("Incorrect time for Ainf calculation. Please review it in Options")
esES("")
frFR("")

T_("Incorrect number of time values for Ainf calculation. Please review it "
     "in Options")
esES("")
frFR("")

T_("Obtaining Impulse Response Function")
esES("")
frFR("")

T_("Cancelled by user")
esES("")
frFR("")

T_("Obtaining Infinite-Frequency Added Mass (A_inf)")
esES("")
frFR("")

T_("Loading files")
esES("")
frFR("")

T_("Model '%s' is already loaded")
esES("")
frFR("")

T_("Problem loading '%s'\n%s")
esES("")
frFR("")

T_("Problem loading '%s'")
esES("")
frFR("")

T_("Unknown BEM file extension in '%s'")
esES("")
frFR("")

T_("Problem processing '%s'\n%s")
esES("")
frFR("")

T_("Discarding negligible DOF")
esES("")
frFR("")

T_("Model has more bodies (%d) than previously loaded (%d)")
esES("")
frFR("")

T_("Problem joining models: '%s'\n%s")
esES("")
frFR("")

T_("Loaded mesh '%s'")
esES("")
frFR("")

T_("Model is already loaded")
esES("")
frFR("")

T_("Healing mesh '%s'")
esES("")
frFR("")

T_("Problem loading '%s': %s")
esES("")
frFR("")

T_("The mesh is in good condition")
esES("")
frFR("")

T_("Getting underwater mesh '%s'")
esES("")
frFR("")

T_("Loading '%s'")
esES("")
frFR("")

T_("Error loading '%s'")
esES("")
frFR("")

T_("Saving '%s'")
esES("")
frFR("")

T_("Error saving '%s'")
esES("")
frFR("")

T_("Usage: bemrosetta_cl [options] [-i infile]... [-e outfile]")
esES("")
frFR("")

T_("Options:")
esES("")
frFR("")

T_("-h  --help     -- print options")
esES("")
frFR("")

T_("-p  --params   -- set physical parameters:")
esES("")
frFR("")

T_("                 parameter description   units  default value")
esES("")
frFR("")

T_("                    g      gravity       [m/s2]    ")
esES("")
frFR("")

T_("                    length length scale  []        ")
esES("")
frFR("")

T_("                    rho    water density [Kg/m3]   ")
esES("")
frFR("")

T_("                    depth  water depth   [m]       ")
esES("")
frFR("")

T_("-i  --input    -- load model")
esES("")
frFR("")

T_("-e  --export   -- export from input file to output file")
esES("")
frFR("")

T_("-c  --compare  -- compare input files")
esES("")
frFR("")

T_("-r  --report   -- output last loaded model data")
esES("")
frFR("")

T_("-cl --clear    -- clear loaded model")
esES("")
frFR("")

T_("Actions")
esES("")
frFR("")

T_("- are done in sequence: if a physical parameter is changed after export, "
     "saved files will not include the change")
esES("")
frFR("")

T_("- can be repeated as desired")
esES("")
frFR("")

T_("Missing parameters when reading '%s'")
esES("")
frFR("")

T_("BEMRosetta Copyright (c) 2019 I\303\261aki Zabala\nHydrodynamic coefficients "
     "converter for Boundary Element Method solver formats\nVersion beta BUILDINFO")
esES("")
frFR("")

T_("BEM configuration data are not loaded. Defaults are set")
esES("")
frFR("")

T_("Command argument list is empty")
esES("")
frFR("")

T_("File '%s' not found")
esES("")
frFR("")

T_("File '%s' loaded")
esES("")
frFR("")

T_("No file loaded")
esES("")
frFR("")

T_("Series cleared")
esES("")
frFR("")

T_("File '%s' converted")
esES("")
frFR("")

T_("Wrong argument '%s'")
esES("")
frFR("")

T_("Unknown argument '%s'")
esES("")
frFR("")

T_("Unknown error")
esES("")
frFR("")

T_("Error")
esES("")
frFR("")

T_("In case of doubt try option -h or --help")
esES("")
frFR("")

T_("or just call command line without arguments to open GUI window")
esES("")
frFR("")


// BEMRosetta.h

T_("Wamit")
esES("")
frFR("")

T_("Wamit.1.3")
esES("")
frFR("")

T_("FAST-Wamit")
esES("")
frFR("")

T_("Nemoh")
esES("")
frFR("")

T_("SeaFEM-Nemoh")
esES("")
frFR("")

T_("AQWA")
esES("")
frFR("")

T_("FOAMM")
esES("")
frFR("")

T_("BEMRosetta")
esES("")
frFR("")

T_("Unknown")
esES("")
frFR("")

T_("W.o")
esES("")
frFR("")

T_("W.1")
esES("")
frFR("")

T_("FST")
esES("")
frFR("")

T_("Nmh")
esES("")
frFR("")

T_("SFM")
esES("")
frFR("")

T_("AQW")
esES("")
frFR("")

T_("FMM")
esES("")
frFR("")

T_("BMR")
esES("")
frFR("")

T_("\302\277?")
esES("")
frFR("")

T_("INFINITY")
esES("")
frFR("")

T_("Wamit.gdf")
esES("")
frFR("")

T_("Wamit.dat")
esES("")
frFR("")

T_("Nemoh.dat")
esES("")
frFR("")

T_("Binary.stl")
esES("")
frFR("")

T_("Text.stl")
esES("")
frFR("")

T_("[File: '%s', line: %d]:\n")
esES("")
frFR("")

T_("No data available")
esES("")
frFR("")

T_("Bad %s '%s' in field #%d, line\n'%s'")
esES("")
frFR("")

T_("Field #%d not found in line\n'%s'")
esES("")
frFR("")


// Mesh.cpp

T_("Unknown MESH file extension in '%s'")
esES("")
frFR("")

T_("Model is empty. No panels found")
esES("")
frFR("")

T_("Unknown mesh file type")
esES("")
frFR("")

T_("Limits [m] (%f - %f, %f - %f, %f - %f)")
esES("")
frFR("")

T_("Water-plane area [m2] %f")
esES("")
frFR("")

T_("Surface [m2] %f")
esES("")
frFR("")

T_("Volume [m3] %f")
esES("")
frFR("")

T_("Underwater surface [m2] %f")
esES("")
frFR("")

T_("Underwater volume [m3] %f")
esES("")
frFR("")

T_("Displacement [tm] %f")
esES("")
frFR("")

T_("Center of buoyancy [m] (%f, %f, %f)")
esES("")
frFR("")

T_("Loaded %d panels and %d nodes")
esES("")
frFR("")


// nemoh_mesh.cpp

T_("Impossible to open file '%s'")
esES("")
frFR("")

T_("Format error in Nemoh .dat mesh file")
esES("")
frFR("")

T_("Impossible to open '%s'\n")
esES("")
frFR("")


// wamit_mesh.cpp

T_("'ZONE' field not found")
esES("")
frFR("")

T_("Wrong scale in .gdf file")
esES("")
frFR("")

T_("Number of patches not found in .gdf file")
esES("")
frFR("")

T_("Wrong number of patches in .gdf file")
esES("")
frFR("")

T_("Impossible to open '%s'")
esES("")
frFR("")


// stl.cpp

T_("'solid' field not found")
esES("")
frFR("")

T_("'facet normal' field not found")
esES("")
frFR("")

T_("'outer loop' field not found")
esES("")
frFR("")

T_("Too much vertex in facet")
esES("")
frFR("")

T_("Label '%s' not handled in facet")
esES("")
frFR("")

T_("Binary stl must not begin with 'solid' text")
esES("")
frFR("")


// nemoh.cpp

T_("Hydrostatics file(s) 'Mesh/Hydrostatics*.dat'")
esES("")
frFR("")

T_("Not found")
esES("")
frFR("")

T_("KH file(s) 'Mesh/KH*.dat'")
esES("")
frFR("")

T_("Radiation file 'RadiationCoefficients.tec'")
esES("")
frFR("")

T_("Excitation force file 'ExcitationForce.tec'")
esES("")
frFR("")

T_("Diffraction force file 'DiffractionForce.tec'")
esES("")
frFR("")

T_("Froude Krylov file 'FKForce.tec'")
esES("")
frFR("")

T_("IRF file(s) 'IRF.tec'")
esES("")
frFR("")

T_("Incorrect number of points in x direction %s")
esES("")
frFR("")

T_("Incorrect free surface elevation %s")
esES("")
frFR("")

T_("Incorrect number of Kochin function directions %s")
esES("")
frFR("")

T_("Incorrect Kochin direction %s")
esES("")
frFR("")

T_("Minimum Kochin direction %s has to be lower than maximum direction %s")
esES("")
frFR("")

T_("Incorrect rho %s")
esES("")
frFR("")

T_("Incorrect g %s")
esES("")
frFR("")

T_("Incorrect depth %s")
esES("")
frFR("")

T_("Incorrect number of bodies %s")
esES("")
frFR("")

T_("Incorrect number of points %s")
esES("")
frFR("")

T_("Incorrect number of panels %s")
esES("")
frFR("")

T_("Incorrect DOF number %s in body %d")
esES("")
frFR("")

T_("Incorrect DOF type %d set in body %d")
esES("")
frFR("")

T_("Incorrect number of frequencies %s")
esES("")
frFR("")

T_("Incorrect frequency %s")
esES("")
frFR("")

T_("Minimum frequency %s has to be lower than maximum frequency %s")
esES("")
frFR("")

T_("Incorrect number of headings %s")
esES("")
frFR("")

T_("Incorrect direction %s")
esES("")
frFR("")

T_("Minimum direction %s has to be lower than maximum direction %s")
esES("")
frFR("")

T_("Incorrect IRF step %s")
esES("")
frFR("")

T_("IRF step %s has to be lower than duration %s")
esES("")
frFR("")

T_("Unexpected data %s")
esES("")
frFR("")

T_("Free surface data is already loaded %s")
esES("")
frFR("")

T_("Kochin data is already loaded %s")
esES("")
frFR("")

T_("'empty'")
esES("")
frFR("")

T_("Incorrect min frequency %s")
esES("")
frFR("")

T_("Incorrect min heading %s")
esES("")
frFR("")

T_("Incorrect max heading %s")
esES("")
frFR("")

T_("Minimum heading %s has to be lower than maximum heading %s")
esES("")
frFR("")

T_("Incorrect number of points in x direction %s (0 for no free surface calculation)")
esES("")
frFR("")

T_("Incorrect free surface domain X %s")
esES("")
frFR("")

T_("Incorrect free surface domain Y %s")
esES("")
frFR("")

T_("Impossible to create '%s'")
esES("")
frFR("")

T_("Number of Nemoh cases must be higher than 1 (%d)")
esES("")
frFR("")

T_("Number of Nemoh cases must not be higher than number of frequencies (%d>%d)")
esES("")
frFR("")

T_("Impossible to clean folder '%s'. Maybe it is in use")
esES("")
frFR("")

T_("Problem creating '%s' folder")
esES("")
frFR("")

T_("Problem copying preprocessor file '%s'")
esES("")
frFR("")

T_("Problem copying solver file '%s'")
esES("")
frFR("")

T_("Problem copying postprocessor file '%s'")
esES("")
frFR("")

T_("Problem copying mesh file '%s'")
esES("")
frFR("")

T_("Problem copying gren file '%s'")
esES("")
frFR("")

T_("SeaFEM_Nemoh only allows one body, found %d")
esES("")
frFR("")

T_("DOF does not match in '%s'")
esES("")
frFR("")

T_("DOF does not match in '%s")
esES("")
frFR("")

T_("Number of frequencies higher than the defined in Nemoh.cal file")
esES("")
frFR("")

T_("Number of bodies higher than the defined in Nemoh.cal file")
esES("")
frFR("")


// foamm.cpp

T_("Loading mat file '%s'")
esES("")
frFR("")

T_("Vector w not found")
esES("")
frFR("")

T_("Vector A not found")
esES("")
frFR("")

T_("Vectors w and A size does not match")
esES("")
frFR("")

T_("Vector B not found")
esES("")
frFR("")

T_("Vectors w and B size does not match")
esES("")
frFR("")

T_("Vector Z not found")
esES("")
frFR("")

T_("Vectors w and Z size does not match")
esES("")
frFR("")

T_("Vector TFSResponse not found")
esES("")
frFR("")

T_("Vectors w and TFSResponse size does not match")
esES("")
frFR("")

T_("Matrix A_ss not found")
esES("")
frFR("")

T_("Matrix B_ss not found")
esES("")
frFR("")

T_("Matrix C_ss not found")
esES("")
frFR("")

T_("Matrix Frequencies not found")
esES("")
frFR("")

T_("Matrix FreqRange not found")
esES("")
frFR("")

T_("FOAMM not found. Please set FOAMM path in Options")
esES("")
frFR("")

T_("Processing case %d")
esES("")
frFR("")

T_("Problem creating temporary FOAMM folder '%s'")
esES("")
frFR("")

T_("Problem creating FOAMM file '%s'")
esES("")
frFR("")

T_("Problem writing %s to file '%s'")
esES("")
frFR("")

T_("Problem launching FOAMM from '%s'")
esES("")
frFR("")

T_("Process ended by user")
esES("")
frFR("")


// fast.cpp

T_("File '%s' is not of FAST type")
esES("")
frFR("")

T_("FAST does not support more than one body in file '%s'")
esES("")
frFR("")

T_("No wave headings found in Wamit file")
esES("")
frFR("")

T_("FAST requires simetric wave headings. .3 file headings found from %f to "
     "%f")
esES("")
frFR("")

T_("Wrong format in FAST file '%s'")
esES("")
frFR("")

T_("State Space file '%s'")
esES("")
frFR("")

T_("Number of bodies different to 1 incompatible with FAST")
esES("")
frFR("")

T_("Volume (PtfmVol0) not found in FAST file '%s'")
esES("")
frFR("")

T_("Density (WtrDens) not found in FAST file '%s'")
esES("")
frFR("")

T_("Water depth (WtrDpth) not found in FAST file '%s'")
esES("")
frFR("")

T_("Length scale (WAMITULEN) not found in FAST file '%s'")
esES("")
frFR("")

T_("Number of wave directions (WaveNDir) not found in FAST file '%s'")
esES("")
frFR("")

T_("Range of wave directions (WaveDirRange) not found in FAST file '%s'")
esES("")
frFR("")

T_("Different %s (%f != %f) in FAST file '%s'")
esES("")
frFR("")

T_("volume")
esES("")
frFR("")

T_("density")
esES("")
frFR("")

T_("water depth")
esES("")
frFR("")

T_("length scale")
esES("")
frFR("")

T_("Different %s (%d != %d) in FAST file '%s'")
esES("")
frFR("")

T_("number of wave headings")
esES("")
frFR("")

T_("headings range")
esES("")
frFR("")

T_("Bad format parsing FAST file '%s' for %s")
esES("")
frFR("")

T_("Imposible to save file '%s'")
esES("")
frFR("")

T_(".ss format only allows to save one body. Only first body is saved")
esES("")
frFR("")

T_("BEMRosetta state space matrices obtained with %s")
esES("")
frFR("")

T_("BEMRosetta state space matrices")
esES("")
frFR("")


// wamit.cpp

T_("Loading out file '%s'")
esES("")
frFR("")

T_("Scattering file '%s'")
esES("")
frFR("")

T_("Froude-Krylov file '%s'")
esES("")
frFR("")

T_("Hydrodynamic coefficients A and B .1 file '%s'")
esES("")
frFR("")

T_("Diffraction exciting .3 file '%s'")
esES("")
frFR("")

T_("Hydrostatic restoring file '%s'")
esES("")
frFR("")

T_("RAO file '%s'")
esES("")
frFR("")

T_("Hydrodynamic coefficients A and B file '%s'")
esES("")
frFR("")

T_("Diffraction exciting file '%s'")
esES("")
frFR("")

T_("Water depth has to be positive")
esES("")
frFR("")

T_("Found additional bodies over %d")
esES("")
frFR("")

T_("cg matrix is not dimensioned")
esES("")
frFR("")

T_("Vo matrix is not dimensioned")
esES("")
frFR("")

T_("cb matrix is not dimensioned")
esES("")
frFR("")

T_("C matrix is not dimensioned")
esES("")
frFR("")

T_("Wrong format in Wamit file '%s'")
esES("")
frFR("")

T_("Found additional frequencies over %d")
esES("")
frFR("")

T_("A matrix is not dimensioned")
esES("")
frFR("")

T_("B matrix is not dimensioned")
esES("")
frFR("")

T_("Index (%d, %d) out of bounds")
esES("")
frFR("")

T_("Index [%d](%d, %d) out of bounds")
esES("")
frFR("")

T_("Index (%d) out of bounds")
esES("")
frFR("")

T_("Number of bodies loaded is lower than previous (%d != %d)")
esES("")
frFR("")

T_("Number of frequencies loaded is different than previous (%d != %d)")
esES("")
frFR("")

T_("Frequencies loaded are different than previous\nPrevious: %s\nSeries: "
     "  %s")
esES("")
frFR("")

T_("Periods loaded are different than previous\nPrevious: %s\nSeries:   %s")
esES("")
frFR("")

T_("DOF # does not match (%d, %d)")
esES("")
frFR("")

T_("A[w=inf] is not expected")
esES("")
frFR("")

T_("A[w=0] is not expected")
esES("")
frFR("")

T_("Frequency %f is unknown")
esES("")
frFR("")

T_("Error in file format")
esES("")
frFR("")

T_("Number of headings is different than previous (%d != %d)")
esES("")
frFR("")

T_("Heading %f is unknown")
esES("")
frFR("")

T_("Number of headings loaded is different than previous (%d != %d)")
esES("")
frFR("")

T_("[%s] Periods loaded are different than previous\nPrevious: %s\nSeries: "
     "  %s")
esES("")
frFR("")

T_("No enough data to save (at least 2 frequencies)")
esES("")
frFR("")


// aqwa.cpp

T_("LIS file")
esES("")
frFR("")

T_("AH1 file")
esES("")
frFR("")

T_("Number of headings do not match %d<>%d")
esES("")
frFR("")

T_("Number of frequencies do not match %d<>%d")
esES("")
frFR("")

T_("Unknown body found in COG %d (%d)")
esES("")
frFR("")

T_("Body # does not match in 'HYDSTIFFNESS' %d<>%d")
esES("")
frFR("")

T_("Body # does not match in '%s' %d<>%d")
esES("")
frFR("")

T_("Frequency # does not match in '%s' %d<>%d")
esES("")
frFR("")

T_("Body # does not match in 'FORCERAO' %d<>%d")
esES("")
frFR("")

T_("Heading # does not match in 'FORCERAO' %d<>%d")
esES("")
frFR("")

T_("Frequency # does not match in 'FORCERAO' %d<>%d")
esES("")
frFR("")

T_("Number of bodies not found")
esES("")
frFR("")

T_("Wrong body %d")
esES("")
frFR("")

T_("Number of frequencies not found")
esES("")
frFR("")

T_("Number of headings not found")
esES("")
frFR("")

T_("Format error, '=' not found")
esES("")
frFR("")

T_("= not found")
esES("")
frFR("")

T_("Expected %s data, found '%s'")
esES("")
frFR("")


// Obsolete

T_("Model '%s' already loaded")
esES("")
frFR("")

T_("Unknown file extension in '%s'")
esES("")
frFR("")

T_("Model already loaded")
esES("")
frFR("")

T_("[line %d] No data available")
esES("")
frFR("")

T_("[line %d] Bad %s '%s' in field #%d, line\n'%s'")
esES("")
frFR("")

T_("[line %d] Field #%d not found in line\n'%s'")
esES("")
frFR("")

T_("Unknown file extension")
esES("")
frFR("")

T_("[line %d] 'ZONE' field not found")
esES("")
frFR("")

T_("[line %d] 'solid' field not found")
esES("")
frFR("")

T_("[line %d] 'facet normal' field not found")
esES("")
frFR("")

T_("[line %d] 'outer loop' field not found")
esES("")
frFR("")

T_("[line %d] Too much vertex in facet")
esES("")
frFR("")

T_("[line %d] Too few vertex in facet")
esES("")
frFR("")

T_("[line %d] Label '%s' not handled in facet")
esES("")
frFR("")

T_("[line %d] Incorrect number of points in x direction %s")
esES("")
frFR("")

T_("[line %d] Incorrect free surface elevation %s")
esES("")
frFR("")

T_("[line %d] Incorrect number of Kochin function directions %s")
esES("")
frFR("")

T_("[line %d] Incorrect Kochin direction %s")
esES("")
frFR("")

T_("[line %d] Minimum Kochin direction %s has to be lower than maximum direction "
     "%s")
esES("")
frFR("")

T_("[line %d] Incorrect rho %s")
esES("")
frFR("")

T_("[line %d] Incorrect g %s")
esES("")
frFR("")

T_("[line %d] Incorrect depth %s")
esES("")
frFR("")

T_("[line %d] Incorrect number of bodies %s")
esES("")
frFR("")

T_("[line %d] Incorrect number of points %s")
esES("")
frFR("")

T_("[line %d] Incorrect number of panels %s")
esES("")
frFR("")

T_("[line %d] Incorrect DOF %s in body %d")
esES("")
frFR("")

T_("[line %d] Incorrect DOF type %d set in body %d")
esES("")
frFR("")

T_("[line %d] Incorrect number of frequencies %s")
esES("")
frFR("")

T_("[line %d] Incorrect frequency %s")
esES("")
frFR("")

T_("[line %d] Minimum frequency %s has to be lower than maximum frequency "
     "%s")
esES("")
frFR("")

T_("[line %d] Incorrect number of headings %s")
esES("")
frFR("")

T_("[line %d] Incorrect direction %s")
esES("")
frFR("")

T_("[line %d] Minimum direction %s has to be lower than maximum direction "
     "%s")
esES("")
frFR("")

T_("[line %d] Incorrect IRF step %s")
esES("")
frFR("")

T_("[line %d] IRF step %s has to be lower than duration %s")
esES("")
frFR("")

T_("[line %d] Unexpected data %s")
esES("")
frFR("")

T_("[line %d] Free surface data already loaded %s")
esES("")
frFR("")

T_("[line %d] Kochin data already loaded %s")
esES("")
frFR("")

T_("DOF does not match in %s")
esES("")
frFR("")

T_("[line %d] DOF does not match in %s")
esES("")
frFR("")

T_("[line %d] Number of frequencies higher than the defined in Nemoh.cal file")
esES("")
frFR("")

T_("[line %d] Number of bodies higher than the defined in Nemoh.cal file")
esES("")
frFR("")

T_("Matrix Frequencies_index not found")
esES("")
frFR("")

T_("[line %d] Water depth has to be positive")
esES("")
frFR("")

T_("[line %d] Found additional bodies over %d")
esES("")
frFR("")

T_("[line %d] cg matrix is not dimensioned")
esES("")
frFR("")

T_("[line %d] Vo matrix is not dimensioned")
esES("")
frFR("")

T_("[line %d] cb matrix is not dimensioned")
esES("")
frFR("")

T_("[line %d] C matrix is not dimensioned")
esES("")
frFR("")

T_("[line %d] Wrong format in Wamit file '%s'")
esES("")
frFR("")

T_("[line %d] Found additional frequencies over %d")
esES("")
frFR("")

T_("[line %d] A matrix is not dimensioned")
esES("")
frFR("")

T_("[line %d] B matrix is not dimensioned")
esES("")
frFR("")

T_("[line %d] Index (%d, %d) out of bounds")
esES("")
frFR("")

T_("[line %d] Index [%d](%d, %d) out of bounds")
esES("")
frFR("")

T_("[line %d] Index (%d) out of bounds")
esES("")
frFR("")

T_("[line %d] Number of bodies loaded is lower than previous (%d != %d)")
esES("")
frFR("")

T_("[line %d] Number of frequencies loaded is different than previous (%d "
     "!= %d)")
esES("")
frFR("")

T_("[line %d] Frequencies loaded are different than previous\nPrevious: %s\n"
     "Series:   %s")
esES("")
frFR("")

T_("[line %d] Periods loaded are different than previous\nPrevious: %s\nSeries: "
     "  %s")
esES("")
frFR("")

T_("[line %d] DOF # does not match (%d, %d)")
esES("")
frFR("")

T_("[line %d] A[w=inf] is not expected")
esES("")
frFR("")

T_("[line %d] A[w=0] is not expected")
esES("")
frFR("")

T_("[line %d] Frequency %f is unknown")
esES("")
frFR("")

T_("[line %d] Number of headings is different than previous (%d != %d)")
esES("")
frFR("")

T_("[line %d] Heading %f is unknown")
esES("")
frFR("")

T_("[line %d] Number of headings loaded is different than previous (%d != "
     "%d)")
esES("")
frFR("")

T_("[line %d] Number of headings do not match %d<>%d")
esES("")
frFR("")

T_("[line %d] Number of frequencies do not match %d<>%d")
esES("")
frFR("")

T_("[line %d] Unknown body found in COG %d (%d)")
esES("")
frFR("")

T_("[line %d] Body # does not match in 'HYDSTIFFNESS' %d<>%d")
esES("")
frFR("")

T_("[line %d] Body # does not match in '%s' %d<>%d")
esES("")
frFR("")

T_("[line %d] Frequency # does not match in '%s' %d<>%d")
esES("")
frFR("")

T_("[line %d] Body # does not match in 'FORCERAO' %d<>%d")
esES("")
frFR("")

T_("[line %d] Heading # does not match in 'FORCERAO' %d<>%d")
esES("")
frFR("")

T_("[line %d] Frequency # does not match in 'FORCERAO' %d<>%d")
esES("")
frFR("")

T_("[line %d] Wrong body %d")
esES("")
frFR("")

T_("[line %d] Format error, '=' not found")
esES("")
frFR("")

T_("[line %d] = not found")
esES("")
frFR("")

T_("[line %d] Expected %s data, found '%s'")
esES("")
frFR("")

T_("delta %.1f [\302\272]")
esES("")
frFR("")

T_("%.1f [\302\272]")
esES("")
frFR("")

T_("[%d] No data available")
esES("")
frFR("")

T_("[%d] Bad %s '%s' in field #%d, line\n'%s'")
esES("")
frFR("")

T_("[%d] Field #%d not found in line\n'%s'")
esES("")
frFR("")

T_("Mesh limits (%f - %f, %f - %f, %f - %f)")
esES("")
frFR("")

T_("Mesh water-plane area %f")
esES("")
frFR("")

T_("Mesh surface %f")
esES("")
frFR("")

T_("Mesh volume %f")
esES("")
frFR("")

T_("Mesh underwater surface %f")
esES("")
frFR("")

T_("Mesh underwater volume %f")
esES("")
frFR("")

T_("Center of buoyancy (%f, %f, %f)")
esES("")
frFR("")

T_("[%d] 'ZONE' field not found")
esES("")
frFR("")

T_("[%d] 'solid' field not found")
esES("")
frFR("")

T_("[%d] 'facet normal' field not found")
esES("")
frFR("")

T_("[%d] 'outer loop' field not found")
esES("")
frFR("")

T_("[%d] Too much vertex in facet")
esES("")
frFR("")

T_("[%d] Too few vertex in facet")
esES("")
frFR("")

T_("[%d] Label '%s' not handled in facet")
esES("")
frFR("")

T_("[%d] Incorrect number of points in x direction %s")
esES("")
frFR("")

T_("[%d] Incorrect free surface elevation %s")
esES("")
frFR("")

T_("[%d] Incorrect number of Kochin function directions %s")
esES("")
frFR("")

T_("[%d] Incorrect Kochin direction %s")
esES("")
frFR("")

T_("[%d] Minimum Kochin direction %s has to be lower than maximum direction "
     "%s")
esES("")
frFR("")

T_("[%d] Incorrect rho %s")
esES("")
frFR("")

T_("[%d] Incorrect g %s")
esES("")
frFR("")

T_("[%d] Incorrect depth %s")
esES("")
frFR("")

T_("[%d] Incorrect number of bodies %s")
esES("")
frFR("")

T_("[%d] Incorrect number of points %s")
esES("")
frFR("")

T_("[%d] Incorrect number of panels %s")
esES("")
frFR("")

T_("[%d] Incorrect DOF %s in body %d")
esES("")
frFR("")

T_("[%d] Incorrect DOF type %d set in body %d")
esES("")
frFR("")

T_("[%d] Incorrect number of frequencies %s")
esES("")
frFR("")

T_("[%d] Incorrect frequency %s")
esES("")
frFR("")

T_("[%d] Minimum frequency %s has to be lower than maximum frequency %s")
esES("")
frFR("")

T_("[%d] Incorrect number of directions %s")
esES("")
frFR("")

T_("[%d] Incorrect direction %s")
esES("")
frFR("")

T_("[%d] Minimum direction %s has to be lower than maximum direction %s")
esES("")
frFR("")

T_("[%d] Incorrect IRF step %s")
esES("")
frFR("")

T_("[%d] IRF step %s has to be lower than duration %s")
esES("")
frFR("")

T_("[%d] Unexpected data %s")
esES("")
frFR("")

T_("[%d] Free surface data already loaded %s")
esES("")
frFR("")

T_("[%d] Kochin data already loaded %s")
esES("")
frFR("")

T_("Incorrect number of directions %s")
esES("")
frFR("")

T_("Problem creating %s folder")
esES("")
frFR("")

T_("[%d] DOF does not match in %s")
esES("")
frFR("")

T_("[%d] Number of frequencies higher than the defined in Nemoh.cal file")
esES("")
frFR("")

T_("[%d] Number of bodies higher than the defined in Nemoh.cal file")
esES("")
frFR("")

T_("FOAMM file '%s'")
esES("")
frFR("")

T_("Unknown frequency %f")
esES("")
frFR("")

T_("Unknown heading %f")
esES("")
frFR("")

T_("[%d] Number of headings do not match %d<>%d")
esES("")
frFR("")

T_("[%d] Number of frequencies do not match %d<>%d")
esES("")
frFR("")

T_("[%d] Unknown body found in COG %d (%d)")
esES("")
frFR("")

T_("[%d] Body # does not match in 'HYDSTIFFNESS' %d<>%d")
esES("")
frFR("")

T_("[%d] Body # does not match in '%s' %d<>%d")
esES("")
frFR("")

T_("[%d] Frequency # does not match in '%s' %d<>%d")
esES("")
frFR("")

T_("[%d] Body # does not match in 'FORCERAO' %d<>%d")
esES("")
frFR("")

T_("[%d] Heading # does not match in 'FORCERAO' %d<>%d")
esES("")
frFR("")

T_("[%d] Frequency # does not match in 'FORCERAO' %d<>%d")
esES("")
frFR("")

T_("[%d] Wrong body %d")
esES("")
frFR("")

T_("[%d] Format error, '=' not found")
esES("")
frFR("")

T_("[%d] = not found")
esES("")
frFR("")

T_("[%d] Heading %f not found")
esES("")
frFR("")

T_("[%d] Frequency %f not found")
esES("")
frFR("")

T_("Loading file")
esES("")
frFR("")

T_("Loading mesh '%s'")
esES("")
frFR("")

T_("Wrong format in Nemoh file '%s'")
esES("")
frFR("")

T_("Unknown file format")
esES("")
frFR("")

T_("Healing mesh")
esES("")
frFR("")

T_("Getting mesh normals")
esES("")
frFR("")

T_("delta %.1f [rad/s]")
esES("")
frFR("")

T_("%.1f [rad/s]")
esES("")
frFR("")

T_("No data loaded")
esES("")
frFR("")

T_("Different number of bodies")
esES("")
frFR("")

T_("Different number of wave headings")
esES("")
frFR("")

T_("Wave headings do not match")
esES("")
frFR("")

T_("Wrong nodes found in Nemoh .dat mesh file")
esES("")
frFR("")

T_("Frequencies loaded are different than previous")
esES("")
frFR("")

T_("Periods loaded are different than previous")
esES("")
frFR("")

T_("Wrong nodes found in Wamit .gdf mesh file")
esES("")
frFR("")

T_("Usage: bemrosetta_cl [option] infile outfolder")
esES("")
frFR("")

T_("-h --help     -- print options")
esES("")
frFR("")

T_("-p --params   -- set physical parameters:")
esES("")
frFR("")

T_("                    g      gravity       [m/s2]    9.81")
esES("")
frFR("")

T_("                    length length scale  []        1")
esES("")
frFR("")

T_("                    rho    water density [Kg/m3]   1000")
esES("")
frFR("")

T_("                    depth  depth         [m]       100")
esES("")
frFR("")

T_("-i --input   -- input file")
esES("")
frFR("")

T_("-e --export  -- export from input file to output file")
esES("")
frFR("")

T_("-c --compare -- compare input files")
esES("")
frFR("")

T_("usage: bemrosetta_cl [options] [-i infile]... [-e outfile]")
esES("")
frFR("")

T_("Wrong number of arguments")
esES("")
frFR("")

T_("g [m/s2]: %.3f, h [m]: %s, rho [kg/m3]: %.3f length scale [m]: %.1f")
esES("")
frFR("")

T_("%.3f to %.3f steps %.3f [rad/s]")
esES("")
frFR("")

T_("%.1f to %.1f steps %.1f [\302\272]")
esES("")
frFR("")

T_("#freqs: %d (%s), #headings: %d (%s)")
esES("")
frFR("")

T_("Bad %s '%s' in field #%d, line '%s'")
esES("")
frFR("")

T_("Field #%d not found in line '%s'")
esES("")
frFR("")

T_("g [m/s2]: %.3f, h [m]: %.3f, rho [kg/m3]: %.3f length scale [m]: %.1f")
esES("")
frFR("")
