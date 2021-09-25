#include <Core/Core.h>
#include <BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;


const double rho = 1000, g = 9.81;
const char *bemBenchmarkFolder = "D:/Others/OneDrive - SENER/BEM Benchmark";
BEMData bem;

void Pre(String device, bool isSymYoZ, bool isSymXoZ) {
	Cout() << "\nPreprocessing for BEM benchmark";
	
	String baseFolder = AppendFileNameX(bemBenchmarkFolder, "40_Calculations", device);
	String aqwaFolder = AppendFileNameX(baseFolder, "AQWA");	
	String meshFolder = AppendFileNameX(bemBenchmarkFolder, "30_Mesh files", device);
	
	FindFile ff(AppendFileName(aqwaFolder, "*.dat"));
	BEMCase bcase0;
	bcase0.Load(ff.GetPath());
		
	
	Cout() << "\nGenerating HAMS cases";
	String hamsFolder = AppendFileNameX(baseFolder, "HAMS");			
	DeleteFolderDeepX(hamsFolder); 		DirectoryCreateX(hamsFolder);	
	for (FindFile ff(AppendFileName(meshFolder, "???????.pnl")); ff; ff++) {
		String file = ff.GetPath();
		Cout() << "\n>> " << GetFileName(file);
	
		BEMCase bcase(bcase0);
		bcase.bodies[0].meshFile = file;
		bcase.bodies[0].lidFile = ForceExt(file, "lid.pnl");
		bcase.g = g;
		
		bcase.SaveFolder(AppendFileName(hamsFolder, GetFileTitle(file)), true, 1, 4, bem, BEMCase::HAMS);
	
	}

	Cout() << "\nGenerating Capytaine, Nemoh and Nemohv115 cases";
	String capyFolder = AppendFileNameX(baseFolder, "Capy");			
	DeleteFolderDeep(capyFolder, true); 			
	DirectoryCreateX(capyFolder);
	String nemohFolder = AppendFileNameX(baseFolder, "Nemoh");			
	DeleteFolderDeep(nemohFolder, true); 		
	DirectoryCreateX(nemohFolder);
	String nemohv115Folder = AppendFileNameX(baseFolder, "Nemohv115");	
	DeleteFolderDeep(nemohv115Folder, true); 	
	DirectoryCreateX(nemohv115Folder);
	String strbatCapy;
	for (FindFile ff(AppendFileName(meshFolder, "???????.dat")); ff; ff++) {
		String file = ff.GetPath();
		Cout() << "\n>> " << GetFileName(file);
	
		bool y0z, x0z;
		MeshData data;
		data.Load(file, rho, g, false, y0z, x0z);
		int numPanels = data.under.panels.size();
		int numNodes = data.under.nodes.size();
	
		BEMCase bcase(bcase0);
		bcase.bodies[0].meshFile = file;
		bcase.bodies[0].npanels = numPanels;
		bcase.bodies[0].npoints = numNodes;
		bcase.g = g;
		
		bcase.SaveFolder(AppendFileName(capyFolder, 	 GetFileTitle(file)), false, 1, 4, 	  bem, BEMCase::CAPYTAINE);
		bcase.SaveFolder(AppendFileName(nemohFolder, 	 GetFileTitle(file)), true,  4, Null, bem, BEMCase::NEMOH);
		bcase.SaveFolder(AppendFileName(nemohv115Folder, GetFileTitle(file)), true,  4, Null, bem, BEMCase::NEMOHv115);
		
		strbatCapy << "cd \"" << AppendFileName(capyFolder, 	 GetFileTitle(file)) << "\"\n";
		strbatCapy << "call Capytaine_bat\n";
	}
	SaveFile(AppendFileName(capyFolder, "Capytaine_bat.bat"), strbatCapy);
	

	
	Cout() << "\nLoading AQWA .dat files and saving them in different formats in '30_Mesh files' folder";

	RealizeDirectory(meshFolder);
	MeshData data;
	for (FindFile ff(AppendFileName(aqwaFolder, "*.dat")); ff; ff++) {
		String file = ff.GetPath();
		Cout() << "\n>> " << GetFileName(file);
		
		bool y0z, x0z;
		data.Load(file, rho, g, false, y0z, x0z);
		int numPanels = data.under.panels.size();
		
		String fileOut;
		
		fileOut = Format("%05d__.gdf", numPanels);
		data.SaveAs(AppendFileName(meshFolder, fileOut), MeshData::WAMIT_GDF, g, MeshData::UNDERWATER, false, false);
		
		fileOut = Format("%05d%s%s.gdf", numPanels, isSymYoZ ? "X": "_", isSymXoZ ? "Y": "_");
		data.SaveAs(AppendFileName(meshFolder, fileOut), MeshData::WAMIT_GDF, g, MeshData::UNDERWATER, isSymYoZ, isSymXoZ);

		fileOut = Format("%05d__.dat", numPanels);
		data.SaveAs(AppendFileName(meshFolder, fileOut), MeshData::NEMOH_DAT, g, MeshData::UNDERWATER, false, false);
		
		fileOut = Format("%05d_%s.dat", numPanels, isSymXoZ ? "Y": "_");
		data.SaveAs(AppendFileName(meshFolder, fileOut), MeshData::NEMOH_DAT, g, MeshData::UNDERWATER, isSymYoZ, isSymXoZ);

		fileOut = Format("%05d__.pnl", numPanels);
		data.SaveAs(AppendFileName(meshFolder, fileOut), MeshData::HAMS_PNL, g, MeshData::UNDERWATER, false, false);
		
		MeshData lid;
		lid.mesh.AddWaterSurface(data.mesh, data.under, 'f'); 
		lid.AfterLoad(rho, g, false, false);
	
		fileOut = Format("%05d__lid.pnl", numPanels);
		lid.SaveAs(AppendFileName(meshFolder, fileOut), MeshData::HAMS_PNL, g, MeshData::ALL, false, false);
					
		y0z = x0z = false;
		if (isSymYoZ) {
			y0z = true;
			x0z = false;
		} else if (isSymXoZ) {
			y0z = false;
			x0z = true;
		}
		if (x0z || y0z) {
			fileOut = Format("%05d%s%s.pnl", numPanels, y0z ? "X": "_", x0z ? "Y": "_");
			data.SaveAs(AppendFileName(meshFolder, fileOut), MeshData::HAMS_PNL, g, MeshData::UNDERWATER, y0z, x0z);
		
			fileOut = Format("%05d%s%s%s.pnl", numPanels, y0z ? "X": "_", x0z ? "Y": "_", "lid");
			lid.SaveAs(AppendFileName(meshFolder, fileOut), MeshData::HAMS_PNL, g, MeshData::ALL, y0z, x0z);
		}
		
	}
}

void Post(String device) {
	Cout() << "\nPostprocessing for BEM benchmark"
		   << "\nLoading result files and converted to Wamit .1 format"
		   << "\nDevice to process: " << device;
	
	String infolder  = AppendFileNameX(bemBenchmarkFolder, "40_Calculations", device);
	String outfolder = AppendFileNameX(bemBenchmarkFolder, "50_Results", device);
	
	String inFolder, outFolder;

	Cout() << "\nConvert Wadam files";
	inFolder = AppendFileName(infolder, "Wadam");
	outFolder = AppendFileName(outfolder, "Wadam");
	RealizeDirectory(outFolder);
	for (FindFile ff(AppendFileName(inFolder, "*.out")); ff; ff++) {
		String file = ff.GetPath();
		String name = GetFileTitle(file);
		Cout() << "\n>> " << name;
		
		bem.LoadBEM(file, Null, false);
		bem.hydros[0].hd().SaveAs(AppendFileName(outFolder, ForceExt(name, ".1")), Null, Hydro::WAMIT_1_3);
		bem.hydros.Clear(); 
	}
	DeleteFileDeepWildcardsX(AppendFileName(outFolder, "*.4"));
		
	Cout() << "\nLink number of panels from mesh files with .LIS files";
	inFolder = AppendFileName(infolder, "AQWA");
	outFolder = AppendFileName(outfolder, "AQWA");
	RealizeDirectory(outFolder);
	for (FindFile ff(AppendFileName(inFolder, "*.dat")); ff; ff++) {
		String file = ff.GetPath();
		String name = GetFileTitle(file);
		Cout() << "\n>> " << name;
		
		bool y0z, x0z;
		MeshData data;
		data.Load(file, rho, g, false, y0z, x0z);
		int numPanels = data.under.panels.size();
		
		bem.LoadBEM(ForceExt(file, ".lis"), Null, false);
		bem.hydros[0].hd().SaveAs(AppendFileName(outFolder, Format("%05d__.1", numPanels)), Null, Hydro::WAMIT_1_3);
		bem.hydros.Clear(); 
	}
	DeleteFileDeepWildcardsX(AppendFileName(outFolder, "*.4"));

	
}

	
	
CONSOLE_APP_MAIN
{
	String device = "Corpower";
	String phase = "pre";
	bool isSymYoZ = true, isSymXoZ = true;
	
	bool firstTime;
	if (!bem.LoadSerializeJson(firstTime))
		Cout() << "\n" << t_("BEM configuration data has not been loaded. Default values are set");
	bem.g = g;

	String errorStr;
	try {
		if (phase == "pre")
			Pre(device, isSymYoZ, isSymXoZ);
		else if (phase == "post")
			Post(device);
		else
			throw Exc(Format("Unknown process %s", phase));
		
	} catch (Exc e) {
		errorStr = e;
	} catch(...) {
		errorStr = t_("Unknown error");
	}	
	if (!errorStr.IsEmpty()) {
		Cerr() << Format("\n%s: %s", t_("Error"), errorStr);
	}
#ifdef flagDEBUG
	Cout() << "\nPress return to end";
	ReadStdIn();
#endif
}
