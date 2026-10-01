topic "C++ API";
[ $$0,0#00000000000000000000000000000000:Default]
[{_}%EN-US
[A4* C`+`+ API&&]
[A3*;l0; Init()&]
[A3;l0; 
Library initialisation&&]
[A3*;l0; Version()&]
[A3;l0; 
Returns the library time and date&&]
[A3*;l0; EnablePrint(bool print)&]
[A3;l0; 
Enables or disables message printing by the functions&&]
[A3*;l0; GetLastError()&]
[A3;l0; 
Returns the last error or NULL if there is no error&&]
[A3*;l0; ClearLastError()&]
[A3;l0; 
Clears the last error&&]
[A3*;l0; Wamit&]
[A3*;l100; V6s&]
[A3*;l200; Set(int force)&]
[A3;l200; 
Consider that the Wamit solver to run the cases will be always WamitV6s&&]
[A3*;l0; AQWA&]
[A3*;l100; ShowCalculationDialog&]
[A3*;l200; Set(int show)&]
[A3;l200; 
Show the calculation dialog when running AQWA&&]
[A3*;l0; Mesh&]
[A3*;l100; Clear()&]
[A3;l100; 
Clear all meshes previously loaded&&]
[A3*;l100; Load(const char `*file)&]
[A3;l100; 
Loads a mesh file&&]
[A3*;l100; Report()&]
[A3;l100; 
Prints main mesh data&&]
[A3*;l100; Id&]
[A3*;l200; Set(int id)&]
[A3;l200; 
Sets the id of the active mesh&&]
[A3*;l200; Get()&]
[A3;l200; 
Gets the id of the active mesh&&]
[A3*;l100; Save(const char `*file`, const char `*format`, int symX`, int symY)&]
[A3;l100; 
Saves the mesh in the indicated file with the mesh format&&]
[A3*;l100; Translate(double x`, double y`, double z)&]
[A3;l100; 
Translates the mesh&&]
[A3*;l100; Rotate(double ax`, double ay`, double az`, double cx`, double cy`, double cz)&]
[A3;l100; 
Rotates the mesh. ax,ay,az are the angles in degrees, cx, cy, cz is the centre of rotation&&]
[A3*;l100; Cg&]
[A3*;l200; Set(double x`, double y`, double z)&]
[A3;l200; 
Sets the centre of gravity&&]
[A3*;l200; Get(double `*x`, double `*y`, double `*z)&]
[A3;l200; 
Gets the centre of gravity&&]
[A3*;l100; C0&]
[A3*;l200; Set(double x`, double y`, double z)&]
[A3;l200; 
Sets the centre of rotation or reference system&&]
[A3*;l200; Get(double `*x`, double `*y`, double `*z)&]
[A3;l200; 
Gets the centre of rotation or reference system&&]
[A3*;l100; Mass&]
[A3*;l200; Set(double mass)&]
[A3;l200; 
Sets the mesh mass&&]
[A3*;l100; Inertia&]
[A3*;l200; Set(const double `*data`, const int dim`[2`])&]
[A3;l200; 
Sets the inertia matrix&&]
[A3*;l100; LinearDamping&]
[A3*;l200; Set(const double `*data`, const int dim`[2`])&]
[A3;l200; 
Sets the linear damping matrix&&]
[A3*;l100; MooringStiffness&]
[A3*;l200; Set(const double `*data`, const int dim`[2`])&]
[A3;l200; 
Sets the mooring stiffness matrix&&]
[A3*;l100; Reset()&]
[A3;l100; 
Reset the position of the mess to the initial condition&&]
[A3*;l100; Duplicate()&]
[A3;l100; 
Duplicates a mesh&&]
[A3*;l100; GetWaterPlane()&]
[A3;l100; 
Extract in new model the waterplane mesh (lid)&&]
[A3*;l100; GetHull()&]
[A3;l100; 
Extract in new model the mesh underwater hull&&]
[A3*;l100; FillWaterplane(double ratio`, int quads)&]
[A3;l100; 
Generates in new model the waterplane lid with mesh size ratio&&]
[A3*;l100; GetControlSurface(double distance`, double ratio`, int quads`, int bottom`, int top)&]
[A3;l100; 
Generates in new model a control surface mesh at a distance and with a size ratio&&]
[A3*;l100; Volume&]
[A3*;l200; Get(double `*volx`, double `*voly`, double `*volz)&]
[A3;l200; 
Returns an array with the volumes x, y and z&&]
[A3*;l100; UnderwaterVolume&]
[A3*;l200; Get(double `*volx`, double `*voly`, double `*volz)&]
[A3;l200; 
Returns an array with the underwater volumes x, y and z&&]
[A3*;l100; Surface&]
[A3*;l200; Get()&]
[A3;l200; 
Returns the body surface&&]
[A3*;l100; UnderwaterSurface&]
[A3*;l200; Get()&]
[A3;l200; 
Returns the body wet surface&&]
[A3*;l100; Centre&]
[A3*;l200; Volume&]
[A3*;l300; Get(double `*x`, double `*y`, double `*z)&]
[A3;l300; 
Gets the centroid of the body&&]
[A3*;l200; Surface&]
[A3*;l300; Get(double `*x`, double `*y`, double `*z)&]
[A3;l300; 
Gets the centre of gravity of the surface&&]
[A3*;l100; HydrostaticStiffness&]
[A3*;l200; Get(double `*`*data`, int dim`[2`])&]
[A3;l200; 
Returns the hydrostatic stiffness matrix&&]
[A3*;l100; NumPanels&]
[A3*;l200; Get()&]
[A3;l200; 
Gets the number of panels of the mesh&&]
[A3*;l100; VolumeEnvelope&]
[A3*;l200; Get(double `*minx`, double `*maxx`, double `*miny`, double `*maxy`, double `*minz`, double `*maxz)&]
[A3;l200; 
Gets the envelope around the mesh&&]
[A3*;l100; Name&]
[A3*;l200; Set(const char `*name)&]
[A3;l200; 
Sets the name of the body&&]
[A3*;l200; Get()&]
[A3;l200; 
Gets the name of the body&&]
[A3*;l0; Bem&]
[A3*;l100; Clear()&]
[A3;l100; 
Clear loaded models&&]
[A3*;l100; New()&]
[A3;l100; 
Creates a new model model&&]
[A3*;l100; Name&]
[A3*;l200; Set(const char `*name)&]
[A3;l200; 
Sets the name of the active model&&]
[A3*;l200; Get()&]
[A3;l200; 
Gets the name of the active model&&]
[A3*;l100; Load(const char `*file)&]
[A3;l100; 
Loads a BEM case&&]
[A3*;l100; Save(const char `*file)&]
[A3;l100; 
Saves a BEM case&&]
[A3*;l100; Id&]
[A3*;l200; Set(int id)&]
[A3;l200; 
Sets the id of the active model&&]
[A3*;l200; Get()&]
[A3;l200; 
Gets the id of the active model&&]
[A3*;l100; size()&]
[A3;l100; 
Gets the number of bem cases&&]
[A3*;l100; Support(const char `*solver`, int `*irregular`, int `*autoIrregular`, int `*middle7`, int `*far8`, int `*near9`, int `*autoCS`, int `*multibody)&]
[A3;l100; 
Returns the multibody and QTF capabilities of the solver: "Wamit .out", "AQWA .dat","OrcaWave .yml", "HydroStar .hsg", "HAMS", "HAMS MREL", "Nemoh v3", "Diffrac .xml", "Capytaine .py"&&]
[A3*;l100; WaveOrigin&]
[A3*;l200; Set(double x`, double y)&]
[A3;l200; 
Sets the wave origin&&]
[A3*;l100; depth&]
[A3*;l200; Set(double h)&]
[A3;l200; 
Set the depth&&]
[A3*;l100; g&]
[A3*;l200; Set(double g)&]
[A3;l200; 
Sets the gravity&&]
[A3*;l100; rho&]
[A3*;l200; Set(double rho)&]
[A3;l200; 
Sets the density&&]
[A3*;l100; w&]
[A3*;l200; Set(const double `*w`, int dim)&]
[A3;l200; 
Sets the range of frequencies&&]
[A3*;l200; Get(double `*`*data`, int dim`[1`])&]
[A3;l200; 
Gets the range of frequencies&&]
[A3*;l200; size()&]
[A3;l200; 
Gets the number of frequencies&&]
[A3*;l100; headings&]
[A3*;l200; Set(const double `*head`, int dim)&]
[A3;l200; 
Sets the range of headings&&]
[A3*;l200; Get(double `*`*data`, int dim`[1`])&]
[A3;l200; 
Gets the range of headings&&]
[A3*;l200; size()&]
[A3;l200; 
Gets the number of headings&&]
[A3*;l100; Duplicate()&]
[A3;l100; 
Copies a Bem case into other&&]
[A3*;l100; SaveCase(const char `*folder`, const char `*solver`, bool x0z`, bool y0z`, bool irregular`, bool autoIrregular`, const char `*qtfType`, bool autoQTF`, bool bin`, int numCases`, int numThreads`, bool withPotentials`, bool withMesh)&]
[A3;l100; 
Saves the Bem case&&]
[A3*;l100; FroudeKrylov&]
[A3*;l200; Calc()&]
[A3;l200; 
Gets the incident or Froude-Krylov wave force for all the bodies&&]
[A3*;l200; Set(int idBody`, int idBemFrom`, int idBodyFrom)&]
[A3;l200; 
Sets the Froude-Krylov force of body idBody with the Froude-Krylov force of body idBodyFrom from the Bem case idBemFrom&&]
[A3*;l100; Diffraction&]
[A3*;l200; Reset()&]
[A3;l200; 
Fill with zeroes for all the bodies the diffraction forces&&]
[A3*;l200; Set(int idBody`, int idBemFrom`, int idBodyFrom)&]
[A3;l200; 
Sets the diffraction force of body idBody with the diffraction force of body idBodyFrom from the Bem case idBemFrom&&]
[A3*;l100; Excitation&]
[A3*;l200; Reset()&]
[A3;l200; 
Fill with zeroes for all the bodies the excitation forces&&]
[A3*;l200; GetFromDiffFK()&]
[A3;l200; 
Gets the excitation force by summing the Froude-Krylov and diffraction scattering forces&&]
[A3*;l200; Set(int idBody`, int idBemFrom`, int idBodyFrom)&]
[A3;l200; 
Sets the excitation force of body idBody with the excitation force of body idBodyFrom from the Bem case idBemFrom&&]
[A3*;l100; AddedMass&]
[A3*;l200; Reset()&]
[A3;l200; 
Fill with zeroes for all the bodies the added mass&&]
[A3*;l200; Set(int idBodyRow`, int idBodyCol`, int idBemFrom`, int idBodyFromRow`, int idBodyFromCol)&]
[A3;l200; 
Sets the added mass of body idBody with the added mass of body idBodyFrom from the Bem case idBemFrom&&]
[A3*;l100; RadiationDamping&]
[A3*;l200; Reset()&]
[A3;l200; 
Fill with zeroes for all the bodies the radiation damoping&&]
[A3*;l200; Set(int idBodyRow`, int idBodyCol`, int idBemFrom`, int idBodyFromRow`, int idBodyFromCol)&]
[A3;l200; 
Sets the radiation damping of body idBody with the added mass of body idBodyFrom from the Bem case idBemFrom&&]
[A3*;l100; MapToMesh(int idBody`, int idMesh`, double tolerance`, bool rad`, bool diff`, bool inc`, bool relatedToBody)&]
[A3;l100; 
Maps the potentials on another mesh, creating a new case. rad, diff and inc indicates if the radiation, diffraction, and incident potentials are loaded. relatedToBody indicates if the results are related to the body centre, or to the mesh centre&&]
[A3*;l100; Mesh&]
[A3*;l200; Load(int idBody`, int idMesh)&]
[A3;l200; 
Loads a mesh from a loaded mesh&&]
[A3*;l200; size()&]
[A3;l200; 
Gets the number of bem cases&&]
[A3*;l100; Lid&]
[A3*;l200; Load(int idBody`, int idMesh)&]
[A3;l200; 
Loads a lid from a loaded mesh&&]
[A3*;l100; ControlSurface&]
[A3*;l200; Load(int idBody`, int idMesh)&]
[A3;l200; 
Loads a control surface from a loaded mesh&&]
[A3*;l0; FAST&]
[A3*;l100; Load(const char `*filename)&]
[A3;l100; 
Loads a FAST .out or .outb file&&]
[A3*;l100; GetParameterName(int id)&]
[A3;l100; 
Returns the parameter name of index id&&]
[A3*;l100; GetUnitName(int id)&]
[A3;l100; 
Returns the parameter units of index id&&]
[A3*;l100; GetParameterId(const char `*name)&]
[A3;l100; 
Returns the index id of parameter name&&]
[A3*;l100; GetParameterCount()&]
[A3;l100; 
Returns the number of parameters&&]
[A3*;l100; GetLen()&]
[A3;l100; 
Returns the number of registers per parameter&&]
[A3*;l100; GetTimeStart()&]
[A3;l100; 
Returns the initial time in seconds&&]
[A3*;l100; GetTimeEnd()&]
[A3;l100; 
Returns the end time in seconds&&]
[A3*;l100; GetTime(int idtime)&]
[A3;l100; 
Returns the idtime_th time (idtime goes from 0 to _BMR_FAST_GetLen())&&]
[A3*;l100; GetData(int idtime`, int idparam)&]
[A3;l100; 
Returns the idtime_th value of parameter idparam (idtime goes from 0 to _BMR_FAST_GetLen())&&]
[A3*;l100; GetAvg(int idparam`, int idbegin`, int idend)&]
[A3;l100; 
Returns the average value for parameter idparam&&]
[A3*;l100; GetMax(int idparam`, int idbegin`, int idend)&]
[A3;l100; 
Returns the maximum value for parameter idparam&&]
[A3*;l100; GetMin(int idparam`, int idbegin`, int idend)&]
[A3;l100; 
Returns the minimum value for parameter idparam&&]
[A3*;l100; GetArray(int idparam`, int idbegin`, int idend`, double `*`*data`, int `*dim)&]
[A3;l100; 
Returns an Array of parameter idparam&&]
[A3*;l100; LoadFile(const char `*file)&]
[A3;l100; 
Open a .dat or .fst FAST file to read or save parameters&&]
[A3*;l100; SaveFile(const char `*file)&]
[A3;l100; 
Saves the .dat or .fst FAST file opened with FAST_LoadFile() (if file is ""), or to the file indicated in file&&]
[A3*;l100; SetVar(const char `*name`, const char `*paragraph`, const char `*value)&]
[A3;l100; 
Sets the value of a var after paragraph. If paragraph is "", the value is set every time var appears in the file&&]
[A3*;l100; GetVar(const char `*name`, const char `*paragraph)&]
[A3;l100; 
Reads the value of a var after paragraph. If paragraph is "", it is read the first time var appears in the file&&]
