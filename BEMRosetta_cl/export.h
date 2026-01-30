// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2026, the BEMRosetta author and contributors

#ifdef _WIN32
    #define L_EXPORT __declspec(dllexport)
#else
    #define L_EXPORT __attribute__((visibility("default")))
#endif

extern "C" {
	// Libary initialization. Called automatically
	L_EXPORT void DLL_Init() noexcept;
	// Returns the library time and date
	L_EXPORT const char *DLL_Version() noexcept; 		
	
	// Blocks the printing of messages on the screen
	L_EXPORT void DLL_NoPrint() noexcept;	
	// Returns the last error or NULL if there is no error
	L_EXPORT const char *DLL_GetLastError() noexcept;
	// Prints a text
	L_EXPORT void DLL_Echo(const char *str) noexcept;
	
	// Loads a mesh file
	L_EXPORT void DLL_Mesh_Input(const char *file) noexcept;
	// Prints main mesh data
	L_EXPORT void DLL_Mesh_Report() noexcept;
	// Clear all meshes previously loaded	
	L_EXPORT void DLL_Mesh_Clear() noexcept;
	// Sets the id of the active mesh
	L_EXPORT void DLL_Mesh_SetId(int id) noexcept;
	// Saves the mesh in the indicated file with the mesh format
	L_EXPORT void DLL_Mesh_Convert(const char *file, const char *format, int symX, int symY) noexcept;
	// Translates the mesh
	L_EXPORT void DLL_Mesh_Translate(double x, double y, double z) noexcept;
	// Rotates the mesh. ax,ay,az are the angles in degrees, cx, cy, cz is the centre of rotation
	L_EXPORT void DLL_Mesh_Rotate(double ax, double ay, double az, double cx, double cy, double cz) noexcept;
	// Sets the centre of gravity
	L_EXPORT void DLL_Mesh_Cg(double x, double y, double z) noexcept;
	// Sets the centre of rotation or reference system
	L_EXPORT void DLL_Mesh_C0(double x, double y, double z) noexcept;		
	// Sets the mesh mass	
	L_EXPORT void DLL_Mesh_Mass(double mass) noexcept;		
	// Reset the position of the mess to the initial condition	
	L_EXPORT void DLL_Mesh_Reset() noexcept;
	// Extract in new model the waterplane mesh (lid)
	L_EXPORT void DLL_Mesh_GetWaterPlane() noexcept;
	// Extract in new model the mesh underwater hull
	L_EXPORT void DLL_Mesh_GetHull() noexcept;
	// Generates in new model the waterplane lid with mesh size ratio
	L_EXPORT void DLL_Mesh_FillWaterplane(double ratio) noexcept;
	// Returns an array with the volumes x, y and z
	L_EXPORT void DLL_Mesh_GetVolume(double *volx, double *voly, double *volz) noexcept;
	// Returns an array with the underwater volumes x, y and z
	L_EXPORT void DLL_Mesh_GetUnderwaterVolume(double *volx, double *voly, double *volz) noexcept;
	// Returns the body surface
	L_EXPORT double DLL_Mesh_GetSurface() noexcept;
	// Returns the body wet surface
	L_EXPORT double DLL_Mesh_GetUnderwaterSurface() noexcept;
	// Gets the centre of buoyancy
	L_EXPORT void DLL_Mesh_GetCentreOfGravity(double *x, double *y, double *z) noexcept;
	// Gets the centroid of the body
	L_EXPORT void DLL_Mesh_GetCentreOfGravity_Volume(double *x, double *y, double *z) noexcept;
	// Gets the centre of gravity of the surface
	L_EXPORT void DLL_Mesh_GetCentreOfGravity_Surface(double *x, double *y, double *z) noexcept;
	// Returns the hydrostatic stiffness matrix
	L_EXPORT void DLL_Mesh_GetHydrostaticStiffness(double **data, int dim[2]) noexcept;
		
		
		
	// Loads a FAST .out or .outb file
	L_EXPORT int DLL_FAST_Load(const char *filename) noexcept;		
	// Returns the parameter name of index id
	L_EXPORT const char *DLL_FAST_GetParameterName(int id) noexcept;
	// Returns the parameter units of index id
	L_EXPORT const char *DLL_FAST_GetUnitName(int id) noexcept;
	// Returns the index id of parameter name
	L_EXPORT int DLL_FAST_GetParameterId(const char *name) noexcept;		
	// Returns the number of parameters
	L_EXPORT int DLL_FAST_GetParameterCount() noexcept;				
	// Returns the number of registers per parameter
	L_EXPORT int DLL_FAST_GetLen() noexcept;			
	// Returns the initial time in seconds
	L_EXPORT double DLL_FAST_GetTimeStart() noexcept;	
	// Returns the end time in seconds
	L_EXPORT double DLL_FAST_GetTimeEnd() noexcept;	
	// Returns the idtime_th time (idtime goes from 0 to DLL_FAST_GetLen())
	L_EXPORT double DLL_FAST_GetTime(int idtime) noexcept;
	// Returns the idtime_th value of parameter idparam (idtime goes from 0 to DLL_FAST_GetLen())
	L_EXPORT double DLL_FAST_GetData(int idtime, int idparam) noexcept;
	// Returns the average value for parameter idparam
	L_EXPORT double DLL_FAST_GetAvg(int idparam, int idbegin, int idend) noexcept;
	// Returns the maximum value for parameter idparam
	L_EXPORT double DLL_FAST_GetMax(int idparam, int idbegin, int idend) noexcept;
	// Returns the minimum value for parameter idparam
	L_EXPORT double DLL_FAST_GetMin(int idparam, int idbegin, int idend) noexcept;
	// Returns an Array of parameter idparam 
	L_EXPORT int DLL_FAST_GetArray(int idparam, int idbegin, int idend, double **data, int *dim) noexcept;
	
	// Open a .dat or .fst FAST file to read or save parameters
	L_EXPORT int DLL_FAST_LoadFile(const char *file) noexcept;
	// Saves the .dat or .fst FAST file opened with FAST_LoadFile() (if file is ""), or to the file indicated in file
	L_EXPORT int DLL_FAST_SaveFile(const char *file) noexcept;
	// Sets the value of a var after paragraph. If paragraph is "", the value is set every time var appears in the file
	L_EXPORT int DLL_FAST_SetVar(const char *name, const char *paragraph, const char *value) noexcept;
	// Reads the value of a var after paragraph. If paragraph is "", it is read the first time var appears in the file
	L_EXPORT const char *DLL_FAST_GetVar(const char *name, const char *paragraph) noexcept;
};
