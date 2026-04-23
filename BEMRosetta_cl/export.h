// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2026, the BEMRosetta author and contributors

#ifdef _WIN32
    #define L_EXPORT __declspec(dllexport)
#else
    #define L_EXPORT __attribute__((visibility("default")))
#endif


extern "C" {
	// Libary initialization. Called automatically
	L_EXPORT void BMR_Init() noexcept;
	// Returns the library time and date
	L_EXPORT const char *BMR_Version() noexcept; 		
	
	// Blocks the printing of messages on the screen
	L_EXPORT void BMR_NoPrint() noexcept;	
	// Returns the last error or NULL if there is no error
	L_EXPORT const char *BMR_GetLastError() noexcept;
	// Prints a text
	L_EXPORT void BMR_Echo(const char *str) noexcept;
	// Consider that the Wamit solver to run the cases will be always WamitV6s
	L_EXPORT void BMR_Wamit_V6s_Set(int force) noexcept;
	// Show the calculation dialog when running AQWA
	L_EXPORT void BMR_AQWA_ShowCalculationDialog_Set(int show) noexcept;
		
	// Clear all meshes previously loaded	
	L_EXPORT void BMR_Mesh_Clear() noexcept;
	// Loads a mesh file
	L_EXPORT int BMR_Mesh_Load(const char *file) noexcept;
	// Prints main mesh data
	L_EXPORT void BMR_Mesh_Report() noexcept;
	// Sets the id of the active mesh
	L_EXPORT void BMR_Mesh_Id_Set(int id) noexcept;
	// Gets the id of the active mesh
	L_EXPORT int BMR_Mesh_Id_Get() noexcept;
	// Saves the mesh in the indicated file with the mesh format
	L_EXPORT void BMR_Mesh_Save(const char *file, const char *format, int symX, int symY) noexcept;
	// Translates the mesh
	L_EXPORT void BMR_Mesh_Translate(double x, double y, double z) noexcept;
	// Rotates the mesh. ax,ay,az are the angles in degrees, cx, cy, cz is the centre of rotation
	L_EXPORT void BMR_Mesh_Rotate(double ax, double ay, double az, double cx, double cy, double cz) noexcept;
	// Sets the centre of gravity
	L_EXPORT void BMR_Mesh_Cg_Set(double x, double y, double z) noexcept;
	// Sets the centre of rotation or reference system
	L_EXPORT void BMR_Mesh_C0_Set(double x, double y, double z) noexcept;		
	// Sets the mesh mass	
	L_EXPORT void BMR_Mesh_Mass_Set(double mass) noexcept;	
	// Sets the inertia matrix
	L_EXPORT void BMR_Mesh_Inertia_Set(const double *data, const int dim[2]) noexcept;	
	// Sets the linear damping matrix
	L_EXPORT void BMR_Mesh_LinearDamping_Set(const double *data, const int dim[2]) noexcept;
	// Sets the mooring stiffness matrix
	L_EXPORT void BMR_Mesh_MooringStiffness_Set(const double *data, const int dim[2]) noexcept;
	// Reset the position of the mess to the initial condition	
	L_EXPORT void BMR_Mesh_Reset() noexcept;
	// Duplicates a mesh
	L_EXPORT int BMR_Mesh_Duplicate() noexcept;
	// Extract in new model the waterplane mesh (lid)
	L_EXPORT int BMR_Mesh_GetWaterPlane() noexcept;
	// Extract in new model the mesh underwater hull
	L_EXPORT int BMR_Mesh_GetHull() noexcept;
	// Generates in new model the waterplane lid with mesh size ratio
	L_EXPORT int BMR_Mesh_FillWaterplane(double ratio, int quads) noexcept;
	// Generates in new model a control surface mesh at a distance and with a size ratio
	L_EXPORT int BMR_Mesh_GetControlSurface(double distance, double ratio, int quads) noexcept;
	// Returns an array with the volumes x, y and z
	L_EXPORT void BMR_Mesh_Volume_Get(double *volx, double *voly, double *volz) noexcept;
	// Returns an array with the underwater volumes x, y and z
	L_EXPORT void BMR_Mesh_UnderwaterVolume_Get(double *volx, double *voly, double *volz) noexcept;
	// Returns the body surface
	L_EXPORT double BMR_Mesh_Surface_Get() noexcept;
	// Returns the body wet surface
	L_EXPORT double BMR_Mesh_UnderwaterSurface_Get() noexcept;
	// Gets the centroid of the body
	L_EXPORT void BMR_Mesh_Centre_Volume_Get(double *x, double *y, double *z) noexcept;
	// Gets the centre of gravity of the surface
	L_EXPORT void BMR_Mesh_Centre_Surface_Get(double *x, double *y, double *z) noexcept;
	// Returns the hydrostatic stiffness matrix
	L_EXPORT void BMR_Mesh_HydrostaticStiffness_Get(double **data, int dim[2]) noexcept;
	// Gets the number of panels of the mesh
	L_EXPORT void BMR_Mesh_NumPanels_Get(int *num) noexcept;
	// Gets the envelope around the mesh
	L_EXPORT void BMR_Mesh_VolumeEnvelope_Get(double *minx, double *maxx, double *miny, double *maxy, double *minz, double *maxz) noexcept;
	
	// Clear loaded models
	L_EXPORT void BMR_Bem_Clear() noexcept;
	// Creates a new model model
	L_EXPORT int BMR_Bem_New() noexcept;
	// Sets the id of the active model
	L_EXPORT void BMR_Bem_Id_Set(int id) noexcept;
	// Gets the id of the active model
	L_EXPORT int BMR_Bem_Id_Get() noexcept;
	// Returns the multibody and QTF capabilities of the solver: "Wamit .out", "AQWA .dat",
	// "OrcaWave .yml", "HydroStar .hsg", "HAMS", "HAMS MREL", "Nemoh v3", "Diffrac .xml", "Capytaine .py"
	L_EXPORT void BMR_Bem_Support(const char *solver, int *irregular, int *autoIrregular, int *middle7, int *far8, int *near9, int *autoCS, int *multibody) noexcept;
	// Set the depth
	L_EXPORT void BMR_Bem_depth_Set(double h) noexcept;
	// Sets the gravity
	L_EXPORT void BMR_Bem_g_Set(double g) noexcept;
	// Sets the density
	L_EXPORT void BMR_Bem_rho_Set(double rho) noexcept;
	// Sets the range of frequencies
	L_EXPORT void BMR_Bem_w_Set(const double *w, int dim) noexcept;
	// Sets the range of headings
	L_EXPORT void BMR_Bem_headings_Set(const double *head, int dim) noexcept;
	// Copies a Bem case into other
	L_EXPORT int BMR_Bem_Duplicate() noexcept;
	// Saves the Bem case
	L_EXPORT void BMR_Bem_SaveCase(const char *folder, const char *solver, bool x0z, bool y0z, bool irregular, bool autoIrregular, const char *qtfType, bool autoQTF, bool bin, int numCases, int numThreads, bool withPotentials, bool withMesh) noexcept;
	
	// Sets the id of the active model
	L_EXPORT void BMR_Bem_Body_Id_Set(int id) noexcept;
	// Gets the id of the active model
	L_EXPORT int BMR_Bem_Body_Id_Get() noexcept;
	// Loads a mesh from a file into the actual body
	L_EXPORT void BMR_Bem_Body_LoadMesh(const char *file) noexcept;
	// Loads a lid from a file into the actual body
	L_EXPORT void BMR_Bem_Body_LoadLid(const char *file) noexcept;
	// Loads a control surface from a file into the actual body
	L_EXPORT void BMR_Bem_Body_LoadControlSurface(const char *file) noexcept;
	// Loads a mesh from a loaded mesh into the actual body
	L_EXPORT void BMR_Bem_Body_LoadMeshFromMesh(int id) noexcept;
	// Loads a lid from a loaded mesh into the actual body
	L_EXPORT void BMR_Bem_Body_LoadLidFromMesh(int id) noexcept;
	// Loads a control surface from a loaded mesh into the actual body
	L_EXPORT void BMR_Bem_Body_LoadControlSurfaceFromMesh(int id) noexcept;
	// Sets the centre of reference of the body
	L_EXPORT void BMR_Bem_Body_C0_Set(double x, double y, double z) noexcept;
	// Gets the centre of reference of the body
	L_EXPORT void BMR_Bem_Body_C0_Get(double *x, double *y, double *z) noexcept;
	// Sets the centre of gravity of the body
	L_EXPORT void BMR_Bem_Body_Cg_Set(double x, double y, double z) noexcept;
	// Gets the centre of gravity of the body
	L_EXPORT void BMR_Bem_Body_Cg_Get(double *x, double *y, double *z) noexcept;
	// Sets the 6x6 inertia matrix of the body
	L_EXPORT void BMR_Bem_Body_Inertia_Set(const double *data, const int dim[2]) noexcept;
	// Sets the name of the body
	L_EXPORT void BMR_Bem_Body_Name_Set(const char *name) noexcept;
		
	// Loads a FAST .out or .outb file
	L_EXPORT int BMR_FAST_Load(const char *filename) noexcept;		
	// Returns the parameter name of index id
	L_EXPORT const char *BMR_FAST_GetParameterName(int id) noexcept;
	// Returns the parameter units of index id
	L_EXPORT const char *BMR_FAST_GetUnitName(int id) noexcept;
	// Returns the index id of parameter name
	L_EXPORT int BMR_FAST_GetParameterId(const char *name) noexcept;		
	// Returns the number of parameters
	L_EXPORT int BMR_FAST_GetParameterCount() noexcept;				
	// Returns the number of registers per parameter
	L_EXPORT int BMR_FAST_GetLen() noexcept;			
	// Returns the initial time in seconds
	L_EXPORT double BMR_FAST_GetTimeStart() noexcept;	
	// Returns the end time in seconds
	L_EXPORT double BMR_FAST_GetTimeEnd() noexcept;	
	// Returns the idtime_th time (idtime goes from 0 to BMR_FAST_GetLen())
	L_EXPORT double BMR_FAST_GetTime(int idtime) noexcept;
	// Returns the idtime_th value of parameter idparam (idtime goes from 0 to BMR_FAST_GetLen())
	L_EXPORT double BMR_FAST_GetData(int idtime, int idparam) noexcept;
	// Returns the average value for parameter idparam
	L_EXPORT double BMR_FAST_GetAvg(int idparam, int idbegin, int idend) noexcept;
	// Returns the maximum value for parameter idparam
	L_EXPORT double BMR_FAST_GetMax(int idparam, int idbegin, int idend) noexcept;
	// Returns the minimum value for parameter idparam
	L_EXPORT double BMR_FAST_GetMin(int idparam, int idbegin, int idend) noexcept;
	// Returns an Array of parameter idparam 
	L_EXPORT int BMR_FAST_GetArray(int idparam, int idbegin, int idend, double **data, int *dim) noexcept;
	
	// Open a .dat or .fst FAST file to read or save parameters
	L_EXPORT int BMR_FAST_LoadFile(const char *file) noexcept;
	// Saves the .dat or .fst FAST file opened with FAST_LoadFile() (if file is ""), or to the file indicated in file
	L_EXPORT int BMR_FAST_SaveFile(const char *file) noexcept;
	// Sets the value of a var after paragraph. If paragraph is "", the value is set every time var appears in the file
	L_EXPORT int BMR_FAST_SetVar(const char *name, const char *paragraph, const char *value) noexcept;
	// Reads the value of a var after paragraph. If paragraph is "", it is read the first time var appears in the file
	L_EXPORT const char *BMR_FAST_GetVar(const char *name, const char *paragraph) noexcept;
};

