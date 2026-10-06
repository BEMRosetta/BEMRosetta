#include <stdio.h>
#include <stdlib.h> 
#include ".test\\libbemrosetta.h"

void my_error_handler(const char *message, void *dummy) {
	printf("\nError %s", message);
	printf("\nClick Enter to end");
	getchar();	
	exit(-1);
}

int main() {
	
		printf("BEMRosetta C demo\n");

	#ifdef BEMROSETTA_DYNAMIC		// DLL is loaded in runtime
		BEMRosetta bmr = BEMRosetta_Init("libbemrosetta.dll");
	#else							// DLL has to be in the .exe folder or in the PATH
		BEMRosetta bmr = BEMRosetta_Init();
	#endif
	
		bmr.SetErrorHandler(my_error_handler, 0);
		printf("\nBEMRosetta version is %s\n", bmr.Version());

		printf("\n- Mesh handling");
		const char *meshFile = "../examples/capytaine/Orca/Body_1.gdf";
		int idMesh = bmr.Mesh.Load(meshFile);
		printf("\nLoaded mesh '%s'", meshFile);

		double volx, voly, volz;
		bmr.Mesh.UnderwaterVolume.Get(&volx, &voly, &volz);
		printf("\nUnderwater volume %f", (volx + voly + volz)/3);
		double x, y, z;
		bmr.Mesh.Centre.Volume.Get(&x, &y, &z);
		printf("\nCentre of volume is %f, %f, %f", x, y, z);

		printf("\n- BEM case generation");
		bmr.Mesh.Cg.Set(0, 0, -1);
		bmr.Mesh.C0.Set(0, 0, -1);
		double M[] = {9.8E5,     0,     0,   0,   0,   0,
						  0, 9.8E5,     0,   0,   0,   0,
						  0,     0, 9.8E5,   0,   0,   0,
						  0,     0,     0, 1E7,   0,   0,
						  0,     0,     0,   0, 1E7,   0,
						  0,     0,     0,   0,   0, 1E7};
		int dim[] = {6, 6};
		bmr.Mesh.Inertia.Set(M, dim);

		bmr.Bem.New();
		bmr.Bem.depth.Set(50);
		bmr.Bem.g.Set(9.81);
		bmr.Bem.rho.Set(1025);
		double w[] = {0.1, 0.5, 1, 1.5, 2};
		bmr.Bem.w.Set(w, sizeof(w)/sizeof(double));
		double head[] = {0, 45, 90};
		bmr.Bem.headings.Set(head, sizeof(head)/sizeof(double));

		bmr.Bem.LoadMesh(0, idMesh);
		
		printf("\nSaving it in Capytaine format");
		bmr.Bem.SaveCase("../unittest/.test/Capy", "Capytaine .py", false, false, true, true, "No", false, false, 1, 4, false, false);

		printf("\nRunning it in Capytaine");
		system("cd /d ..\\unittest\\.test\\Capy && capytaine.bat");
		
		printf("\nLoading the results in .nc format");
		bmr.Bem.Load("..\\unittest\\.test\\Capy\\capytaine.nc");
		printf("\nSaving the results in .h5 format");
		bmr.Bem.Save("..\\unittest\\.test\\Capy\\capytaine.h5");
	
	printf("\nProgram ended\n");

#ifdef BEMROSETTA_DYNAMIC
    BEMRosetta_Free();
#endif
		
	return 0;
}