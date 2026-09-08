// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2026, the BEMRosetta author and contributors

#ifdef _WIN32
    #define L_EXPORT __declspec(dllexport)

	#include <windows.h>
#else
    #define L_EXPORT __attribute__((visibility("default")))
    
    #include <dlfcn.h>
#endif

#include <limits.h>
#include <math.h>

#ifdef __cplusplus
	#include <stdexcept>
	#include <string>
	
	inline bool IsError(int x)    {return x == INT_MIN;}
	inline bool IsError(float x)  {return isnan(x);}
	inline bool IsError(double x) {return isnan(x);}
	
	extern "C" {
	#define NOEXCEPT noexcept
	
#else
	#define NOEXCEPT
	#define bool int
	#define true  1
	#define false 0	
	
	#define IsError(x) _Generic((x),  \
	    int:    ((x) == INT_MIN),     \
	    float:  isnan(x),             \
	    double: isnan(x)              \
	)
#endif

#define NullInt		INT_MIN
#define NullFloat	((float)NAN)
#define NullDouble	((double)NAN)


	typedef void (*error_callback_t)(const char *, void *);
	
#ifdef BEMROSETTA_DYNAMIC	
	static error_callback_t (*_BMR_GetErrorHandlerCallback)() = 0;
	static void *(*_BMR_GetErrorHandlerData)() = 0;
#else
	L_EXPORT error_callback_t _BMR_GetErrorHandlerCallback() NOEXCEPT;
	L_EXPORT void *_BMR_GetErrorHandlerData() NOEXCEPT;
#endif

// C DLL functions

	// Sets a callback called when an error is produced
	L_EXPORT void _BMR_SetErrorHandler(void (*message)(const char *, void *), void *data) NOEXCEPT;
 
      
// DLL functions

	// Library initialisation
	L_EXPORT void _BMR_Init() NOEXCEPT;
	// Returns the library time and date
	L_EXPORT const char *_BMR_Version() NOEXCEPT; 		
	// Enables or disables message printing by the functions
	L_EXPORT void _BMR_EnablePrint(bool print) NOEXCEPT;
	
	// Returns the last error or NULL if there is no error
	L_EXPORT const char *_BMR_GetLastError() NOEXCEPT;
	// Clears the last error
	L_EXPORT void _BMR_ClearLastError() NOEXCEPT;
	
	// Consider that the Wamit solver to run the cases will be always WamitV6s
	L_EXPORT void _BMR_Wamit_V6s_Set(int force) NOEXCEPT;
	
	// Show the calculation dialog when running AQWA
	L_EXPORT void _BMR_AQWA_ShowCalculationDialog_Set(int show) NOEXCEPT;
		
	// Clear all meshes previously loaded	
	L_EXPORT void _BMR_Mesh_Clear() NOEXCEPT;
	// Loads a mesh file
	L_EXPORT int _BMR_Mesh_Load(const char *file) NOEXCEPT;
	// Prints main mesh data
	L_EXPORT bool _BMR_Mesh_Report() NOEXCEPT;
	// Sets the id of the active mesh
	L_EXPORT bool _BMR_Mesh_Id_Set(int id) NOEXCEPT;
	// Gets the id of the active mesh
	L_EXPORT int _BMR_Mesh_Id_Get() NOEXCEPT;
	// Saves the mesh in the indicated file with the mesh format
	L_EXPORT bool _BMR_Mesh_Save(const char *file, const char *format, int symX, int symY) NOEXCEPT;
	// Translates the mesh
	L_EXPORT bool _BMR_Mesh_Translate(double x, double y, double z) NOEXCEPT;
	// Rotates the mesh. ax,ay,az are the angles in degrees, cx, cy, cz is the centre of rotation
	L_EXPORT bool _BMR_Mesh_Rotate(double ax, double ay, double az, double cx, double cy, double cz) NOEXCEPT;
	// Sets the centre of gravity
	L_EXPORT bool _BMR_Mesh_Cg_Set(double x, double y, double z) NOEXCEPT;
	// Gets the centre of gravity
	L_EXPORT bool _BMR_Mesh_Cg_Get(double *x, double *y, double *z) NOEXCEPT;
	// Sets the centre of rotation or reference system
	L_EXPORT bool _BMR_Mesh_C0_Set(double x, double y, double z) NOEXCEPT;		
	// Gets the centre of rotation or reference system
	L_EXPORT bool _BMR_Mesh_C0_Get(double *x, double *y, double *z) NOEXCEPT;		
	// Sets the mesh mass	
	L_EXPORT bool _BMR_Mesh_Mass_Set(double mass) NOEXCEPT;	
	// Sets the inertia matrix
	L_EXPORT bool _BMR_Mesh_Inertia_Set(const double *data, const int dim[2]) NOEXCEPT;	
	// Sets the linear damping matrix
	L_EXPORT bool _BMR_Mesh_LinearDamping_Set(const double *data, const int dim[2]) NOEXCEPT;
	// Sets the mooring stiffness matrix
	L_EXPORT bool _BMR_Mesh_MooringStiffness_Set(const double *data, const int dim[2]) NOEXCEPT;
	// Reset the position of the mess to the initial condition	
	L_EXPORT bool _BMR_Mesh_Reset() NOEXCEPT;
	// Duplicates a mesh
	L_EXPORT int _BMR_Mesh_Duplicate() NOEXCEPT;
	// Extract in new model the waterplane mesh (lid)
	L_EXPORT int _BMR_Mesh_GetWaterPlane() NOEXCEPT;
	// Extract in new model the mesh underwater hull
	L_EXPORT int _BMR_Mesh_GetHull() NOEXCEPT;
	// Generates in new model the waterplane lid with mesh size ratio
	L_EXPORT int _BMR_Mesh_FillWaterplane(double ratio, int quads) NOEXCEPT;
	// Generates in new model a control surface mesh at a distance and with a size ratio
	L_EXPORT int _BMR_Mesh_GetControlSurface(double distance, double ratio, int quads, int bottom, int top) NOEXCEPT;
	// Returns an array with the volumes x, y and z
	L_EXPORT bool _BMR_Mesh_Volume_Get(double *volx, double *voly, double *volz) NOEXCEPT;
	// Returns an array with the underwater volumes x, y and z
	L_EXPORT bool _BMR_Mesh_UnderwaterVolume_Get(double *volx, double *voly, double *volz) NOEXCEPT;
	// Returns the body surface
	L_EXPORT double _BMR_Mesh_Surface_Get() NOEXCEPT;
	// Returns the body wet surface
	L_EXPORT double _BMR_Mesh_UnderwaterSurface_Get() NOEXCEPT;
	// Gets the centroid of the body
	L_EXPORT bool _BMR_Mesh_Centre_Volume_Get(double *x, double *y, double *z) NOEXCEPT;
	// Gets the centre of gravity of the surface
	L_EXPORT bool _BMR_Mesh_Centre_Surface_Get(double *x, double *y, double *z) NOEXCEPT;
	// Returns the hydrostatic stiffness matrix
	L_EXPORT bool _BMR_Mesh_HydrostaticStiffness_Get(double **data, int dim[2]) NOEXCEPT;
	// Gets the number of panels of the mesh
	L_EXPORT int _BMR_Mesh_NumPanels_Get() NOEXCEPT;
	// Gets the envelope around the mesh
	L_EXPORT bool _BMR_Mesh_VolumeEnvelope_Get(double *minx, double *maxx, double *miny, double *maxy, double *minz, double *maxz) NOEXCEPT;
	// Sets the name of the body
	L_EXPORT bool _BMR_Mesh_Name_Set(const char *name) NOEXCEPT;
	// Gets the name of the body
	L_EXPORT const char *_BMR_Mesh_Name_Get() NOEXCEPT;
		
	// Clear loaded models
	L_EXPORT void _BMR_Bem_Clear() NOEXCEPT;
	// Creates a new model model
	L_EXPORT int _BMR_Bem_New() NOEXCEPT;
	// Sets the name of the active model
	L_EXPORT bool _BMR_Bem_Name_Set(const char *name) NOEXCEPT;
	// Gets the name of the active model
	L_EXPORT const char *_BMR_Bem_Name_Get() NOEXCEPT;
	// Loads a BEM case
	L_EXPORT int _BMR_Bem_Load(const char *file) NOEXCEPT;
	// Saves a BEM case
	L_EXPORT bool _BMR_Bem_Save(const char *file) NOEXCEPT;
	// Sets the id of the active model
	L_EXPORT bool _BMR_Bem_Id_Set(int id) NOEXCEPT;
	// Gets the id of the active model
	L_EXPORT int _BMR_Bem_Id_Get() NOEXCEPT;
	// Gets the number of bem cases
	L_EXPORT int _BMR_Bem_size() NOEXCEPT;
	// Returns the multibody and QTF capabilities of the solver: "Wamit .out", "AQWA .dat",
	// "OrcaWave .yml", "HydroStar .hsg", "HAMS", "HAMS MREL", "Nemoh v3", "Diffrac .xml", "Capytaine .py"
	L_EXPORT bool _BMR_Bem_Support(const char *solver, int *irregular, int *autoIrregular, int *middle7, int *far8, int *near9, int *autoCS, int *multibody) NOEXCEPT;
	// Sets the wave origin
	L_EXPORT bool _BMR_Bem_WaveOrigin_Set(double x, double y) NOEXCEPT;
	// Set the depth
	L_EXPORT bool _BMR_Bem_depth_Set(double h) NOEXCEPT;
	// Sets the gravity
	L_EXPORT bool _BMR_Bem_g_Set(double g) NOEXCEPT;
	// Sets the density
	L_EXPORT bool _BMR_Bem_rho_Set(double rho) NOEXCEPT;
	// Sets the range of frequencies
	L_EXPORT bool _BMR_Bem_w_Set(const double *w, int dim) NOEXCEPT;
	// Gets the range of frequencies
	L_EXPORT bool _BMR_Bem_w_Get(double **data, int dim[1]) NOEXCEPT;
	// Gets the number of frequencies
	L_EXPORT int _BMR_Bem_w_size() NOEXCEPT;
	// Sets the range of headings
	L_EXPORT bool _BMR_Bem_headings_Set(const double *head, int dim) NOEXCEPT;
	// Gets the range of headings
	L_EXPORT bool _BMR_Bem_headings_Get(double **data, int dim[1]) NOEXCEPT;
	// Gets the number of headings
	L_EXPORT int _BMR_Bem_headings_size() NOEXCEPT;
	// Copies a Bem case into other
	L_EXPORT int _BMR_Bem_Duplicate() NOEXCEPT;
	// Saves the Bem case
	L_EXPORT bool _BMR_Bem_SaveCase(const char *folder, const char *solver, bool x0z, bool y0z, bool irregular, bool autoIrregular, const char *qtfType, bool autoQTF, bool bin, int numCases, int numThreads, bool withPotentials, bool withMesh) NOEXCEPT;
	
	// Gets the incident or Froude-Krylov wave force for all the bodies
	L_EXPORT bool _BMR_Bem_FroudeKrylov_Calc() NOEXCEPT;
	// Fill with zeroes for all the bodies the diffraction forces
	L_EXPORT bool _BMR_Bem_Diffraction_Reset() NOEXCEPT;
	// Fill with zeroes for all the bodies the excitation forces
	L_EXPORT bool _BMR_Bem_Excitation_Reset() NOEXCEPT;
	// Gets the excitation force by summing the Froude-Krylov and diffraction scattering forces
	L_EXPORT bool _BMR_Bem_Excitation_GetFromDiffFK() NOEXCEPT;		// GetFexFromFscFfk
	// Fill with zeroes for all the bodies the added mass
	L_EXPORT bool _BMR_Bem_AddedMass_Reset() NOEXCEPT;
	// Fill with zeroes for all the bodies the radiation damoping
	L_EXPORT bool _BMR_Bem_RadiationDamping_Reset() NOEXCEPT;
	// Sets the added mass of body idBody with the added mass of body idBodyFrom from the Bem case idBemFrom
	L_EXPORT bool _BMR_Bem_AddedMass_Set(int idBodyRow, int idBodyCol, int idBemFrom, int idBodyFromRow, int idBodyFromCol) NOEXCEPT;
	// Sets the radiation damping of body idBody with the added mass of body idBodyFrom from the Bem case idBemFrom
	L_EXPORT bool _BMR_Bem_RadiationDamping_Set(int idBodyRow, int idBodyCol, int idBemFrom, int idBodyFromRow, int idBodyFromCol) NOEXCEPT;
	// Sets the Froude-Krylov force of body idBody with the Froude-Krylov force of body idBodyFrom from the Bem case idBemFrom
	L_EXPORT bool _BMR_Bem_FroudeKrylov_Set(int idBody, int idBemFrom, int idBodyFrom) NOEXCEPT;
	// Sets the diffraction force of body idBody with the diffraction force of body idBodyFrom from the Bem case idBemFrom
	L_EXPORT bool _BMR_Bem_Diffraction_Set(int idBody, int idBemFrom, int idBodyFrom) NOEXCEPT;
	// Sets the excitation force of body idBody with the excitation force of body idBodyFrom from the Bem case idBemFrom
	L_EXPORT bool _BMR_Bem_Excitation_Set(int idBody, int idBemFrom, int idBodyFrom) NOEXCEPT;
	// Maps the potentials on another mesh, creating a new case. rad, diff and inc indicates if the radiation, diffraction, and incident potentials are loaded. relatedToBody indicates if the results are related to the body centre, or to the mesh centre
	L_EXPORT int _BMR_Bem_MapToMesh(int idBody, int idMesh, double tolerance, bool rad, bool diff, bool inc, bool relatedToBody) NOEXCEPT;
	
	// Loads a mesh from a loaded mesh
	L_EXPORT bool _BMR_Bem_Mesh_Load(int idBody, int idMesh) NOEXCEPT;
	// Loads a lid from a loaded mesh
	L_EXPORT bool _BMR_Bem_Lid_Load(int idBody, int idMesh) NOEXCEPT;
	// Loads a control surface from a loaded mesh
	L_EXPORT bool _BMR_Bem_ControlSurface_Load(int idBody, int idMesh) NOEXCEPT;
	// Gets the number of bem cases
	L_EXPORT int _BMR_Bem_Mesh_size() NOEXCEPT;
		
	// Loads a FAST .out or .outb file
	L_EXPORT bool _BMR_FAST_Load(const char *filename) NOEXCEPT;		
	// Returns the parameter name of index id
	L_EXPORT const char *_BMR_FAST_GetParameterName(int id) NOEXCEPT;
	// Returns the parameter units of index id
	L_EXPORT const char *_BMR_FAST_GetUnitName(int id) NOEXCEPT;
	// Returns the index id of parameter name
	L_EXPORT int _BMR_FAST_GetParameterId(const char *name) NOEXCEPT;		
	// Returns the number of parameters
	L_EXPORT int _BMR_FAST_GetParameterCount() NOEXCEPT;				
	// Returns the number of registers per parameter
	L_EXPORT int _BMR_FAST_GetLen() NOEXCEPT;			
	// Returns the initial time in seconds
	L_EXPORT double _BMR_FAST_GetTimeStart() NOEXCEPT;	
	// Returns the end time in seconds
	L_EXPORT double _BMR_FAST_GetTimeEnd() NOEXCEPT;	
	// Returns the idtime_th time (idtime goes from 0 to _BMR_FAST_GetLen())
	L_EXPORT double _BMR_FAST_GetTime(int idtime) NOEXCEPT;
	// Returns the idtime_th value of parameter idparam (idtime goes from 0 to _BMR_FAST_GetLen())
	L_EXPORT double _BMR_FAST_GetData(int idtime, int idparam) NOEXCEPT;
	// Returns the average value for parameter idparam
	L_EXPORT double _BMR_FAST_GetAvg(int idparam, int idbegin, int idend) NOEXCEPT;
	// Returns the maximum value for parameter idparam
	L_EXPORT double _BMR_FAST_GetMax(int idparam, int idbegin, int idend) NOEXCEPT;
	// Returns the minimum value for parameter idparam
	L_EXPORT double _BMR_FAST_GetMin(int idparam, int idbegin, int idend) NOEXCEPT;
	// Returns an Array of parameter idparam 
	L_EXPORT bool _BMR_FAST_GetArray(int idparam, int idbegin, int idend, double **data, int *dim) NOEXCEPT;
	
	// Open a .dat or .fst FAST file to read or save parameters
	L_EXPORT bool _BMR_FAST_LoadFile(const char *file) NOEXCEPT;
	// Saves the .dat or .fst FAST file opened with FAST_LoadFile() (if file is ""), or to the file indicated in file
	L_EXPORT bool _BMR_FAST_SaveFile(const char *file) NOEXCEPT;
	// Sets the value of a var after paragraph. If paragraph is "", the value is set every time var appears in the file
	L_EXPORT bool _BMR_FAST_SetVar(const char *name, const char *paragraph, const char *value) NOEXCEPT;
	// Reads the value of a var after paragraph. If paragraph is "", it is read the first time var appears in the file
	L_EXPORT const char *_BMR_FAST_GetVar(const char *name, const char *paragraph) NOEXCEPT;

// End DLL functions

#ifdef __cplusplus
}
#endif

