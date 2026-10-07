topic "C++ API";
[ $$0,0#00000000000000000000000000000000:Default]
[{_}%EN-US
[A4* C`+`+ API&&]
[A3*;l0; void [A1*@(150.150.150) .]Init()&]
[A3;l0; 
Library initialisation&&]
[A3*;l0; const char * [A1*@(150.150.150) .]Version()&]
[A3;l0; 
Returns the library time and date&&]
[A3*;l0; void [A1*@(150.150.150) .]EnablePrint(bool print)&]
[A3;l0; 
Enables or disables message printing by the functions&&]
[A3*;l0; const char * [A1*@(150.150.150) .]GetLastError()&]
[A3;l0; 
Returns the last error or NULL if there is no error&&]
[A3*;l0; void [A1*@(150.150.150) .]ClearLastError()&]
[A3;l0; 
Clears the last error&&]
[A3*;l0; [A1*@(150.150.150) .]Wamit&]
[A3*;l150; [A1*@(150.150.150) .Wamit.]V6s&]
[A3*;l300; void [A1*@(150.150.150) .Wamit.V6s.]Set(int force)&]
[A3;l300; 
Consider that the Wamit solver to run the cases will be always WamitV6s&&]
[A3*;l0; [A1*@(150.150.150) .]AQWA&]
[A3*;l150; [A1*@(150.150.150) .AQWA.]ShowCalculationDialog&]
[A3*;l300; void [A1*@(150.150.150) .AQWA.ShowCalculationDialog.]Set(int show)&]
[A3;l300; 
Show the calculation dialog when running AQWA&&]
[A3*;l0; [A1*@(150.150.150) .]Mesh&]
[A3*;l150; void [A1*@(150.150.150) .Mesh.]Clear()&]
[A3;l150; 
Clear all meshes previously loaded&&]
[A3*;l150; int [A1*@(150.150.150) .Mesh.]Load(const char `*file)&]
[A3;l150; 
Loads a mesh file&&]
[A3*;l150; bool [A1*@(150.150.150) .Mesh.]Report()&]
[A3;l150; 
Prints main mesh data&&]
[A3*;l150; [A1*@(150.150.150) .Mesh.]Id&]
[A3*;l300; bool [A1*@(150.150.150) .Mesh.Id.]Set(int id)&]
[A3;l300; 
Sets the id of the active mesh&&]
[A3*;l300; int [A1*@(150.150.150) .Mesh.Id.]Get()&]
[A3;l300; 
Gets the id of the active mesh&&]
[A3*;l150; bool [A1*@(150.150.150) .Mesh.]Save(const char `*file`, const char `*format`, int symX`, int symY)&]
[A3;l150; 
Saves the mesh in the indicated file with the mesh format&&]
[A3*;l150; bool [A1*@(150.150.150) .Mesh.]Translate(double x`, double y`, double z)&]
[A3;l150; 
Translates the mesh&&]
[A3*;l150; bool [A1*@(150.150.150) .Mesh.]Rotate(double ax`, double ay`, double az`, double cx`, double cy`, double cz)&]
[A3;l150; 
Rotates the mesh. ax,ay,az are the angles in degrees, cx, cy, cz is the centre of rotation&&]
[A3*;l150; [A1*@(150.150.150) .Mesh.]Cg&]
[A3*;l300; bool [A1*@(150.150.150) .Mesh.Cg.]Set(double x`, double y`, double z)&]
[A3;l300; 
Sets the centre of gravity&&]
[A3*;l300; bool [A1*@(150.150.150) .Mesh.Cg.]Get(double `*x`, double `*y`, double `*z)&]
[A3;l300; 
Gets the centre of gravity&&]
[A3*;l150; [A1*@(150.150.150) .Mesh.]C0&]
[A3*;l300; bool [A1*@(150.150.150) .Mesh.C0.]Set(double x`, double y`, double z)&]
[A3;l300; 
Sets the centre of rotation or reference system&&]
[A3*;l300; bool [A1*@(150.150.150) .Mesh.C0.]Get(double `*x`, double `*y`, double `*z)&]
[A3;l300; 
Gets the centre of rotation or reference system&&]
[A3*;l150; [A1*@(150.150.150) .Mesh.]Mass&]
[A3*;l300; bool [A1*@(150.150.150) .Mesh.Mass.]Set(double mass)&]
[A3;l300; 
Sets the mesh mass&&]
[A3*;l150; [A1*@(150.150.150) .Mesh.]Inertia&]
[A3*;l300; bool [A1*@(150.150.150) .Mesh.Inertia.]Set(const double `*data`, const int dim`[2`])&]
[A3;l300; 
Sets the inertia matrix&&]
[A3*;l150; [A1*@(150.150.150) .Mesh.]LinearDamping&]
[A3*;l300; bool [A1*@(150.150.150) .Mesh.LinearDamping.]Set(const double `*data`, const int dim`[2`])&]
[A3;l300; 
Sets the linear damping matrix&&]
[A3*;l150; [A1*@(150.150.150) .Mesh.]MooringStiffness&]
[A3*;l300; bool [A1*@(150.150.150) .Mesh.MooringStiffness.]Set(const double `*data`, const int dim`[2`])&]
[A3;l300; 
Sets the mooring stiffness matrix&&]
[A3*;l150; bool [A1*@(150.150.150) .Mesh.]Reset()&]
[A3;l150; 
Reset the position of the mess to the initial condition&&]
[A3*;l150; int [A1*@(150.150.150) .Mesh.]Duplicate()&]
[A3;l150; 
Duplicates a mesh&&]
[A3*;l150; int [A1*@(150.150.150) .Mesh.]GetWaterPlane()&]
[A3;l150; 
Extract in new model the waterplane mesh (lid)&&]
[A3*;l150; int [A1*@(150.150.150) .Mesh.]GetHull()&]
[A3;l150; 
Extract in new model the mesh underwater hull&&]
[A3*;l150; int [A1*@(150.150.150) .Mesh.]FillWaterplane(double ratio`, int quads)&]
[A3;l150; 
Generates in new model the waterplane lid with mesh size ratio&&]
[A3*;l150; int [A1*@(150.150.150) .Mesh.]GetControlSurface(double distance`, double ratio`, int quads`, int bottom`, int top)&]
[A3;l150; 
Generates in new model a control surface mesh at a distance and with a size ratio&&]
[A3*;l150; [A1*@(150.150.150) .Mesh.]Volume&]
[A3*;l300; bool [A1*@(150.150.150) .Mesh.Volume.]Get(double `*volx`, double `*voly`, double `*volz)&]
[A3;l300; 
Returns an array with the volumes x, y and z&&]
[A3*;l150; [A1*@(150.150.150) .Mesh.]UnderwaterVolume&]
[A3*;l300; bool [A1*@(150.150.150) .Mesh.UnderwaterVolume.]Get(double `*volx`, double `*voly`, double `*volz)&]
[A3;l300; 
Returns an array with the underwater volumes x, y and z&&]
[A3*;l150; [A1*@(150.150.150) .Mesh.]Surface&]
[A3*;l300; double [A1*@(150.150.150) .Mesh.Surface.]Get()&]
[A3;l300; 
Returns the body surface&&]
[A3*;l150; [A1*@(150.150.150) .Mesh.]UnderwaterSurface&]
[A3*;l300; double [A1*@(150.150.150) .Mesh.UnderwaterSurface.]Get()&]
[A3;l300; 
Returns the body wet surface&&]
[A3*;l150; [A1*@(150.150.150) .Mesh.]Centre&]
[A3*;l300; [A1*@(150.150.150) .Mesh.Centre.]Volume&]
[A3*;l450; bool [A1*@(150.150.150) .Mesh.Centre.Volume.]Get(double `*x`, double `*y`, double `*z)&]
[A3;l450; 
Gets the centroid of the body&&]
[A3*;l300; [A1*@(150.150.150) .Mesh.Centre.]Surface&]
[A3*;l450; bool [A1*@(150.150.150) .Mesh.Centre.Surface.]Get(double `*x`, double `*y`, double `*z)&]
[A3;l450; 
Gets the centre of gravity of the surface&&]
[A3*;l150; [A1*@(150.150.150) .Mesh.]HydrostaticStiffness&]
[A3*;l300; bool [A1*@(150.150.150) .Mesh.HydrostaticStiffness.]Get(double `*`*data`, int dim`[2`])&]
[A3;l300; 
Returns the hydrostatic stiffness matrix&&]
[A3*;l150; [A1*@(150.150.150) .Mesh.]NumPanels&]
[A3*;l300; int [A1*@(150.150.150) .Mesh.NumPanels.]Get()&]
[A3;l300; 
Gets the number of panels of the mesh&&]
[A3*;l150; [A1*@(150.150.150) .Mesh.]VolumeEnvelope&]
[A3*;l300; bool [A1*@(150.150.150) .Mesh.VolumeEnvelope.]Get(double `*minx`, double `*maxx`, double `*miny`, double `*maxy`, double `*minz`, double `*maxz)&]
[A3;l300; 
Gets the envelope around the mesh&&]
[A3*;l150; [A1*@(150.150.150) .Mesh.]Name&]
[A3*;l300; bool [A1*@(150.150.150) .Mesh.Name.]Set(const char `*name)&]
[A3;l300; 
Sets the name of the body&&]
[A3*;l300; const char * [A1*@(150.150.150) .Mesh.Name.]Get()&]
[A3;l300; 
Gets the name of the body&&]
[A3*;l0; [A1*@(150.150.150) .]Bem&]
[A3*;l150; void [A1*@(150.150.150) .Bem.]Clear()&]
[A3;l150; 
Clear loaded models&&]
[A3*;l150; int [A1*@(150.150.150) .Bem.]New()&]
[A3;l150; 
Creates a new model model&&]
[A3*;l150; [A1*@(150.150.150) .Bem.]Name&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.Name.]Set(const char `*name)&]
[A3;l300; 
Sets the name of the active model&&]
[A3*;l300; const char * [A1*@(150.150.150) .Bem.Name.]Get()&]
[A3;l300; 
Gets the name of the active model&&]
[A3*;l150; int [A1*@(150.150.150) .Bem.]Load(const char `*file)&]
[A3;l150; 
Loads a BEM case&&]
[A3*;l150; bool [A1*@(150.150.150) .Bem.]Save(const char `*file)&]
[A3;l150; 
Saves a BEM case&&]
[A3*;l150; [A1*@(150.150.150) .Bem.]Id&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.Id.]Set(int id)&]
[A3;l300; 
Sets the id of the active model&&]
[A3*;l300; int [A1*@(150.150.150) .Bem.Id.]Get()&]
[A3;l300; 
Gets the id of the active model&&]
[A3*;l150; int [A1*@(150.150.150) .Bem.]size()&]
[A3;l150; 
Gets the number of bem cases&&]
[A3*;l150; bool [A1*@(150.150.150) .Bem.]Support(const char `*solver`, int `*irregular`, int `*autoIrregular`, int `*middle7`, int `*far8`, int `*near9`, int `*autoCS`, int `*multibody)&]
[A3;l150; 
Returns the multibody and QTF capabilities of the solver: "Wamit .out", "AQWA .dat","OrcaWave .yml", "HydroStar .hsg", "HAMS", "HAMS MREL", "Nemoh v3", "Diffrac .xml", "Capytaine .py"&&]
[A3*;l150; [A1*@(150.150.150) .Bem.]WaveOrigin&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.WaveOrigin.]Set(double x`, double y)&]
[A3;l300; 
Sets the wave origin&&]
[A3*;l150; [A1*@(150.150.150) .Bem.]depth&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.depth.]Set(double h)&]
[A3;l300; 
Set the depth&&]
[A3*;l150; [A1*@(150.150.150) .Bem.]g&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.g.]Set(double g)&]
[A3;l300; 
Sets the gravity&&]
[A3*;l150; [A1*@(150.150.150) .Bem.]rho&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.rho.]Set(double rho)&]
[A3;l300; 
Sets the density&&]
[A3*;l150; [A1*@(150.150.150) .Bem.]w&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.w.]Set(const double `*w`, int dim)&]
[A3;l300; 
Sets the range of frequencies&&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.w.]Get(double `*`*data`, int dim`[1`])&]
[A3;l300; 
Gets the range of frequencies&&]
[A3*;l300; int [A1*@(150.150.150) .Bem.w.]size()&]
[A3;l300; 
Gets the number of frequencies&&]
[A3*;l150; [A1*@(150.150.150) .Bem.]headings&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.headings.]Set(const double `*head`, int dim)&]
[A3;l300; 
Sets the range of headings&&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.headings.]Get(double `*`*data`, int dim`[1`])&]
[A3;l300; 
Gets the range of headings&&]
[A3*;l300; int [A1*@(150.150.150) .Bem.headings.]size()&]
[A3;l300; 
Gets the number of headings&&]
[A3*;l150; int [A1*@(150.150.150) .Bem.]Duplicate()&]
[A3;l150; 
Copies a Bem case into other&&]
[A3*;l150; bool [A1*@(150.150.150) .Bem.]SaveCase(const char `*folder`, const char `*solver`, bool x0z`, bool y0z`, bool irregular`, bool autoIrregular`, const char `*qtfType`, bool autoQTF`, bool bin`, int numCases`, int numThreads`, bool withPotentials`, bool withMesh)&]
[A3;l150; 
Saves the Bem case&&]
[A3*;l150; [A1*@(150.150.150) .Bem.]FroudeKrylov&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.FroudeKrylov.]Calc()&]
[A3;l300; 
Gets the incident or Froude-Krylov wave force for all the bodies&&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.FroudeKrylov.]Set(int idBody`, int idBemFrom`, int idBodyFrom)&]
[A3;l300; 
Sets the Froude-Krylov force of body idBody with the Froude-Krylov force of body idBodyFrom from the Bem case idBemFrom&&]
[A3*;l150; [A1*@(150.150.150) .Bem.]Diffraction&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.Diffraction.]Reset()&]
[A3;l300; 
Fill with zeroes for all the bodies the diffraction forces&&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.Diffraction.]Set(int idBody`, int idBemFrom`, int idBodyFrom)&]
[A3;l300; 
Sets the diffraction force of body idBody with the diffraction force of body idBodyFrom from the Bem case idBemFrom&&]
[A3*;l150; [A1*@(150.150.150) .Bem.]Excitation&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.Excitation.]Reset()&]
[A3;l300; 
Fill with zeroes for all the bodies the excitation forces&&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.Excitation.]GetFromDiffFK()&]
[A3;l300; 
Gets the excitation force by summing the Froude-Krylov and diffraction scattering forces&&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.Excitation.]Set(int idBody`, int idBemFrom`, int idBodyFrom)&]
[A3;l300; 
Sets the excitation force of body idBody with the excitation force of body idBodyFrom from the Bem case idBemFrom&&]
[A3*;l150; [A1*@(150.150.150) .Bem.]AddedMass&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.AddedMass.]Reset()&]
[A3;l300; 
Fill with zeroes for all the bodies the added mass&&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.AddedMass.]Set(int idBodyRow`, int idBodyCol`, int idBemFrom`, int idBodyFromRow`, int idBodyFromCol)&]
[A3;l300; 
Sets the added mass of body idBody with the added mass of body idBodyFrom from the Bem case idBemFrom&&]
[A3*;l150; [A1*@(150.150.150) .Bem.]RadiationDamping&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.RadiationDamping.]Reset()&]
[A3;l300; 
Fill with zeroes for all the bodies the radiation damoping&&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.RadiationDamping.]Set(int idBodyRow`, int idBodyCol`, int idBemFrom`, int idBodyFromRow`, int idBodyFromCol)&]
[A3;l300; 
Sets the radiation damping of body idBody with the added mass of body idBodyFrom from the Bem case idBemFrom&&]
[A3*;l150; int [A1*@(150.150.150) .Bem.]MapToMesh(int idBody`, int idMesh`, double tolerance`, bool rad`, bool diff`, bool inc`, bool relatedToBody)&]
[A3;l150; 
Maps the potentials on another mesh, creating a new case. rad, diff and inc indicates if the radiation, diffraction, and incident potentials are loaded. relatedToBody indicates if the results are related to the body centre, or to the mesh centre&&]
[A3*;l150; [A1*@(150.150.150) .Bem.]Mesh&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.Mesh.]Load(int idBody`, int idMesh)&]
[A3;l300; 
Loads a mesh from a loaded mesh&&]
[A3*;l300; int [A1*@(150.150.150) .Bem.Mesh.]size()&]
[A3;l300; 
Gets the number of bem cases&&]
[A3*;l150; [A1*@(150.150.150) .Bem.]Lid&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.Lid.]Load(int idBody`, int idMesh)&]
[A3;l300; 
Loads a lid from a loaded mesh&&]
[A3*;l150; [A1*@(150.150.150) .Bem.]ControlSurface&]
[A3*;l300; bool [A1*@(150.150.150) .Bem.ControlSurface.]Load(int idBody`, int idMesh)&]
[A3;l300; 
Loads a control surface from a loaded mesh&&]
[A3*;l0; [A1*@(150.150.150) .]FAST&]
[A3*;l150; bool [A1*@(150.150.150) .FAST.]Load(const char `*filename)&]
[A3;l150; 
Loads a FAST .out or .outb file&&]
[A3*;l150; const char * [A1*@(150.150.150) .FAST.]GetParameterName(int id)&]
[A3;l150; 
Returns the parameter name of index id&&]
[A3*;l150; const char * [A1*@(150.150.150) .FAST.]GetUnitName(int id)&]
[A3;l150; 
Returns the parameter units of index id&&]
[A3*;l150; int [A1*@(150.150.150) .FAST.]GetParameterId(const char `*name)&]
[A3;l150; 
Returns the index id of parameter name&&]
[A3*;l150; int [A1*@(150.150.150) .FAST.]GetParameterCount()&]
[A3;l150; 
Returns the number of parameters&&]
[A3*;l150; int [A1*@(150.150.150) .FAST.]GetLen()&]
[A3;l150; 
Returns the number of registers per parameter&&]
[A3*;l150; double [A1*@(150.150.150) .FAST.]GetTimeStart()&]
[A3;l150; 
Returns the initial time in seconds&&]
[A3*;l150; double [A1*@(150.150.150) .FAST.]GetTimeEnd()&]
[A3;l150; 
Returns the end time in seconds&&]
[A3*;l150; double [A1*@(150.150.150) .FAST.]GetTime(int idtime)&]
[A3;l150; 
Returns the idtime_th time (idtime goes from 0 to _BMR_FAST_GetLen())&&]
[A3*;l150; double [A1*@(150.150.150) .FAST.]GetData(int idtime`, int idparam)&]
[A3;l150; 
Returns the idtime_th value of parameter idparam (idtime goes from 0 to _BMR_FAST_GetLen())&&]
[A3*;l150; double [A1*@(150.150.150) .FAST.]GetAvg(int idparam`, int idbegin`, int idend)&]
[A3;l150; 
Returns the average value for parameter idparam&&]
[A3*;l150; double [A1*@(150.150.150) .FAST.]GetMax(int idparam`, int idbegin`, int idend)&]
[A3;l150; 
Returns the maximum value for parameter idparam&&]
[A3*;l150; double [A1*@(150.150.150) .FAST.]GetMin(int idparam`, int idbegin`, int idend)&]
[A3;l150; 
Returns the minimum value for parameter idparam&&]
[A3*;l150; bool [A1*@(150.150.150) .FAST.]GetArray(int idparam`, int idbegin`, int idend`, double `*`*data`, int `*dim)&]
[A3;l150; 
Returns an Array of parameter idparam&&]
[A3*;l150; bool [A1*@(150.150.150) .FAST.]LoadFile(const char `*file)&]
[A3;l150; 
Open a .dat or .fst FAST file to read or save parameters&&]
[A3*;l150; bool [A1*@(150.150.150) .FAST.]SaveFile(const char `*file)&]
[A3;l150; 
Saves the .dat or .fst FAST file opened with FAST_LoadFile() (if file is ""), or to the file indicated in file&&]
[A3*;l150; bool [A1*@(150.150.150) .FAST.]SetVar(const char `*name`, const char `*paragraph`, const char `*value)&]
[A3;l150; 
Sets the value of a var after paragraph. If paragraph is "", the value is set every time var appears in the file&&]
[A3*;l150; const char * [A1*@(150.150.150) .FAST.]GetVar(const char `*name`, const char `*paragraph)&]
[A3;l150; 
Reads the value of a var after paragraph. If paragraph is "", it is read the first time var appears in the file&&]
