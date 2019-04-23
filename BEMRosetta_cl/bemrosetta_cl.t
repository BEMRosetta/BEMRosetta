#ifdef _MSC_VER
#pragma setlocale("C")
#endif
// main.cpp

T_("Usage: bemrosetta_cl [option] infile outfolder")

T_("Options:")

T_("-h --help     -- print options")

T_("-p --params   -- set physical parameters:")

T_("                 parameter description   units  default value")

T_("                    g      gravity       [m/s2]    9.81")

T_("                    length length scale  []        1")

T_("                    rho    water density [Kg/m3]   1000")

T_("                    depth  depth         [m]       100")

T_("-i --input   -- input file")

T_("-e --export  -- export from input file to output file")

T_("-c --compare -- compare input files")

T_("usage: bemrosetta_cl [options] [-i infile]... [-e outfile]")

T_("Actions")

T_("- are done in sequence: if a physical parameter is changed after export, "
     "saved files will not include the change")

T_("- can be repeated as desired")

T_("Wrong number of arguments")

T_("BEMRosetta Copyright (c) 2019 I\303\261aki Zabala\nHydrodynamic coefficients "
     "converter for Boundary Element Method solver formats\nVersion beta BUILDINFO")

T_("BEM configuration data are not loaded. Defaults are set")

T_("Wrong argument '%s'")

T_("Unknown argument '%s'")

T_("Unknown error")

T_("Error")


// bemrosetta.cpp

T_("surge")

T_("sway")

T_("heave")

T_("roll")

T_("pitch")

T_("yaw")

T_("s")

T_("w")

T_("h")

T_("r")

T_("p")

T_("y")

T_("%s file '%s'")

T_("g [m/s2]: %.3f, h [m]: %s, rho [kg/m3]: %.3f length scale [m]: %.1f")

T_("INFINITY")

T_("NONE")

T_("%.3f to %.3f steps %.3f [rad/s]")

T_("%.1f to %.1f steps %.1f [\302\272]")

T_("%.1f [\302\272]")

T_("#freqs: %d (%s), #headings: %d (%s)")

T_("#bodies: %d")

T_("dof")

T_("vol [m3]")

T_("No data loaded")

T_("Different number of bodies")

T_("Different number of wave headings")

T_("Wave headings do not match")


// bemrosetta.h

T_("Wamit")

T_("Wamit.1.3")

T_("FAST-Wamit")

T_("Nemoh")

T_("SeaFEM-Nemoh")

T_("Unknown")

T_("Bad %s '%s' in field #%d, line '%s'")

T_("Field #%d not found in line '%s'")


// nemoh.cpp

T_("Loading '%s'")

T_("File '%s' not found")

T_("Hydrostatics file(s) 'Mesh/Hydrostatics*.dat'")

T_("Not found")

T_("KH file(s) 'Mesh/KH*.dat'")

T_("Diffraction force file 'DiffractionForce.tec'")

T_("Froude Krylov file 'FKForce.tec'")

T_("IRF file(s) 'IRF.tec'")

T_("Water depth has to be positive")

T_("Wrong format in Nemoh file '%s'")

T_("SeaFEM_Nemoh only allows one body, found %d")

T_("Impossible to open file '%s'")

T_("Format error in Nemoh .dat mesh file")

T_("Wrong nodes found in Nemoh .dat mesh file")


// wamit.cpp

T_("Loading out file '%s'")

T_("Scattering file '%s'")

T_("Froude-Krylov file '%s'")

T_("Hydrodynamic coefficients A and B .1 file '%s'")

T_("Diffraction exciting .3 file '%s'")

T_("Hydrostatic restoring file '%s'")

T_("RAO file '%s'")

T_("Hydrodynamic coefficients A and B file '%s'")

T_("Diffraction exciting file '%s'")

T_("Wrong format in Wamit file '%s'")

T_("Number of bodies loaded is lower than previous (%d != %d)")

T_("Number of frequencies loaded is different than previous (%d != %d)")

T_("Frequencies loaded are different than previous")

T_("Periods loaded are different than previous")

T_("Error in file format")

T_("Number of headings is different than previous (%d != %d)")

T_("Unknown frequency %f")

T_("Unknown heading %f")

T_("Number of headings loaded is different than previous (%d != %d)")

T_("Impossible to open '%s'")

T_("Wrong number of patches in .gdf file")

T_("Wrong nodes found in Wamit .gdf mesh file")


// fast.cpp

T_("File '%s' is not of FAST type")

T_("FAST does not support more than one body in file '%s'")

T_("No wave headings found in Wamit file")

T_("FAST requires simetric wave headings. .3 file headings found from %f to "
     "%f")

T_("Wrong format in FAST file '%s'")

T_("Saving '%s'")

T_("Different %s (%f != %f) in FAST file '%s'")

T_("volume")

T_("density")

T_("water depth")

T_("length scale")

T_("Different %s (%d != %d) in FAST file '%s'")

T_("number of wave headings")

T_("headings range")

T_("Bad format parsing FAST file '%s' for %s")

T_("Imposible to save file '%s'")


// Obsolete

T_("g [m/s2]: %.3f, h [m]: %.3f, rho [kg/m3]: %.3f length scale [m]: %.1f")
