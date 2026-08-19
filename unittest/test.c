#include <stdio.h>
#include <stdlib.h> 
#include ".test\\libbemrosetta.h"

int main() {
	printf("BEMRosetta C demo\n");
	
	BMR_Init();
	
	printf("\nBEMRosetta version is %s\n", BMR_Version());
	
	printf("\n- Mesh handling");
	const char *meshFile = "../examples/capytaine/Orca/Body_1.gdf";
	BMR_Mesh_Load(meshFile);
	printf("\nLoaded mesh '%s'", meshFile);

	double volx, voly, volz;
	BMR_Mesh_UnderwaterVolume_Get(&volx, &voly, &volz);
	printf("\nUnderwater volume %f", (volx + voly + volz)/3);
	double x, y, z;
	BMR_Mesh_Centre_Volume_Get(&x, &y, &z);
	printf("\nCentre of volume is %f, %f, %f", x, y, z);

	printf("\n- BEM case generation");
	BMR_Bem_depth_Set(50);
	BMR_Bem_g_Set(9.81);
	BMR_Bem_rho_Set(1025);
	double w[] = {0.1, 0.5, 1, 1.5, 2};
	BMR_Bem_w_Set(w, sizeof(w)/sizeof(double));
	double head[] = {0, 45, 90};
	BMR_Bem_headings_Set(head, sizeof(head)/sizeof(double));
	BMR_Bem_Body_LoadMesh(meshFile);
	BMR_Bem_Body_Cg_Set(0, 0, -1);
	BMR_Bem_Body_C0_Set(0, 0, -1);
	double M[] = {9.8E5,     0,     0,   0,   0,   0,
					  0, 9.8E5,     0,   0,   0,   0,
					  0,     0, 9.8E5,   0,   0,   0,
					  0,     0,     0, 1E7,   0,   0,
					  0,     0,     0,   0, 1E7,   0,
					  0,     0,     0,   0,   0, 1E7};
	int dim[] = {6, 6};
	BMR_Bem_Body_Inertia_Set(M, dim);
	printf("\nSaving it in Capytaine format");
	BMR_Bem_SaveCase("../unittest/.test/Capy", "Capytaine .py", false, false, true, true, "No", false, false, 1, 4, false, false);
	
	printf("\nRunning it in Capytaine");
	system("cd /d ..\\unittest\\.test\\Capy && capytaine.bat");
	
	printf("\n- Loading the results from one format and converting them to other format");
	printf("\nLoading the results in .nc format");
	BMR_Bem_Load("..\\unittest\\.test\\Capy\\capytaine.nc");
	printf("\nSaving the results in .h5 format");
	BMR_Bem_Save("..\\unittest\\.test\\Capy\\capytaine.h5");
			
	printf("\nClick Enter to end");
	getchar();
	return 0;
}