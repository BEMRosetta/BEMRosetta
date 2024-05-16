// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2023, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#include "FastOut.h"
#include <Eigen/MultiDimMatrixIndex.h>
#ifdef PLATFORM_WIN32
#include "orca.h"

int Orca::deltaLogSimulation = 10;
UVector<int> Orca::objTypes;
UVector<String> Orca::objNames;
UVector<HINSTANCE> Orca::objHandles;
int Orca::actualBlade;
UVector<int> Orca::varIDs, Orca::varBlades;
UVector<String> Orca::varNames, Orca::varFullNames, Orca::varUnits;

VectorMap<int, String> Orca::state = {
	{msReset, "Reset"},
	{msCalculatingStatics, "Calculating Statics"},
	{msInStaticState, "In Static State"},
	{msRunningSimulation, "Running Simulation"},
	{msSimulationStopped, "Simulation finished OK"},
	{msSimulationStoppedUnstable, "Simulation Stopped Unstable"}
};

VectorMap<int, String> Orca::objectTypes = {
	{otNull, "Null"},
	{otGeneral, "General"},
	{otEnvironment, "Environment"},
	{otVessel, "Vessel"},
	{otLine, "Line"},
	{ot6DBuoy, "6DBuoy"},
	{ot3DBuoy, "3DBuoy"},
	{otWinch, "Winch"},
	{otLink , "Link"},
	{otShape , "Shape"},
	{otConstraint, "Constraint"},
	{otTurbine, "Turbine"},
	{otDragChain, "DragChain"},
	{otLineType, "LineType"},
	{otClumpType, "ClumpType"},
	{otWingType, "WingType"},
	{otVesselType, "VesselType"},
	{otDragChainType, "DragChainType"},
	{otFlexJointType, "FlexJointType"},
	{otStiffenerType, "StiffenerType"},
	{otFlexJoint, "FlexJoint"},
	{otAttachedBuoy, "AttachedBuoy"},
	{otFrictionCoefficients, "FrictionCoefficients"},
	{otSolidFrictionCoefficients, "FrictionCoefficients"},
	{otRayleighDampingCoefficients, "RayleighDampingCoefficients"},
	{otWakeModel, "WakeModel"},
	{otPyModel, "PyModel"},
	{otLineContact, "LineContact"},
	{otCodeChecks, "CodeChecks"},
	{otShear7Data, "Shear7Data"},
	{otVIVAData, "VIVAData"},
	{otSupportType, "SupportType"},
	{otMorisonElementType, "MorisonElementType"},
	{otExpansionTable, "ExpansionTable"},
	{otBrowserGroup, "BrowserGroup"},
	{otMultiBodyGroup, "MultiBodyGroup"},
    {otMultipleObjects, "MultipleObjects"},
    {otDragCoefficient, "DragCoefficient"},
    {otAxialStiffness, "AxialStiffness"},
    {otBendingStiffness, "BendingStiffness"},
    {otBendingConnectionStiffness, "BendingConnectionStiffness"},
    {otWingOrientation, "WingOrientation"},
    {otKinematicViscosity, "KinematicViscosity"},
    {otFluidTemperature, "FluidTemperature"},
    {otCurrentSpeed, "CurrentSpeed"},
    {otCurrentDirection, "CurrentDirection"},
    {otExternalFunction, "ExternalFunction"},
    {otHorizontalVariationFactor, "HorizontalVariationFactor"},
    {otLoadForce, "LoadForce"},
    {otLoadMoment, "LoadMoment"},
    {otExpansionFactor, "ExpansionFactor"},
    {otPayoutRate, "PayoutRate"},
    {otWinchPayoutRate, "PayoutRate"},
    {otWinchTension, "WinchTension"},
    {otVerticalVariationFactor, "VerticalVariationFactor"},
    {otTorsionalStiffness, "TorsionalStiffness"},
    {otMinimumBendRadius, "MinimumBendRadius"},
    {otLiftCoefficient, "LiftCoefficient"},
    {otLiftCloseToSeabed, "LiftCloseToSeabed"},
    {otDragCloseToSeabed, "DragCloseToSeabed"},
    {otDragAmplificationFactor, "DragAmplificationFactor"},
    {otLineTypeDiameter, "LineTypeDiameter"},
    {otStressStrainRelationship, "StressStrainRelationship"},
    {otCoatingOrLining, "CoatingOrLining"},
    {otContentsFlowVelocity, "ContentsFlowVelocity"},
    {otAddedMassRateOfChangeCloseToSurface, "AddedMassRateOfChangeCloseToSurface"},
    {otAddedMassCloseToSurface, "AddedMassCloseToSurface"},
    {otContactStiffness, "ContactStiffness"},
    {otSupportsStiffness, "SupportsStiffness"},
    {otConstraintTranslationalStiffness, "ConstraintTranslationalStiffness"},
    {otConstraintRotationalStiffness, "ConstraintRotationalStiffness"},
    {otConstraintTranslationalDamping, "ConstraintTranslationalDamping"},
    {otConstraintRotationalDamping, "ConstraintRotationalDamping"},
    {otAddedMassCloseToSeabed, "AddedMassCloseToSeabed"},
    {otSeabedTangentialResistance, "SeabedTangentialResistance"}
};

Orca::~Orca() {
	if (dll.GetHandle() != 0) {
		int status;
		
		if (wave) {
			DestroyDiffraction(wave, &status);
			if (status != 0)
				throwError("DestroyDiffraction");
		}
		
		if (flex) {
			DestroyModel(flex, &status);
			if (status != 0)
				throwError("DestroyModel");
		}
		FinaliseLibrary(&status);
		if (status != 0)
			throwError("FinaliseLibrary");
	}
}
	
bool Orca::Init(String _dllFile) {
	dllFile = _dllFile;
	if (!dll.Load(dllFile)) {
		BEM::PrintWarning(Format(t_("OrcaFlex DLL is not found in '%s'"), dllFile));
		return false;
	}

	CreateModel  = DLLGetFunction(dll, void, C_CreateModel, (HINSTANCE *handle, HWND hCaller, int *status));
	DestroyModel = DLLGetFunction(dll, void, C_DestroyModel,(HINSTANCE handle, int *status));
	LoadData     = DLLGetFunction(dll, void, C_LoadDataW,   (HINSTANCE handle, LPCWSTR wcs, int *status));
	SaveData     = DLLGetFunction(dll, void, C_SaveDataW,   (HINSTANCE handle, LPCWSTR wcs, int *status));
	
	CreateDiffraction      = DLLGetFunction(dll, void, C_CreateDiffraction, 	  (HINSTANCE *handle, int *status));
	DestroyDiffraction     = DLLGetFunction(dll, void, C_DestroyDiffraction,	  (HINSTANCE handle, int *status));
	LoadDiffractionData    = DLLGetFunction(dll, void, C_LoadDiffractionDataW,    (HINSTANCE handle, LPCWSTR wcs, int *status));
	SaveDiffractionData    = DLLGetFunction(dll, void, C_SaveDiffractionDataW,    (HINSTANCE handle, LPCWSTR wcs, int *status));
	CalculateDiffraction   = DLLGetFunction(dll, void, C_CalculateDiffractionW,	  (HINSTANCE handle, TProgressHandlerProc proc, int *status));
	LoadDiffractionResults = DLLGetFunction(dll, void, C_LoadDiffractionResultsW, (HINSTANCE handle, LPCWSTR wcs, int *status));
	SaveDiffractionResults = DLLGetFunction(dll, void, C_SaveDiffractionResultsW, (HINSTANCE handle, LPCWSTR wcs, int *status));
	CalculateStatics	   = DLLGetFunction(dll, void, C_CalculateStaticsW, 	  (HINSTANCE handle, TProgressHandlerProc proc, int *status));
	RunSimulation		   = DLLGetFunction(dll, void, C_RunSimulation2W, 		  (HINSTANCE handle, TSimulationHandlerProc proc, const TRunSimulationParameters *lpRunSimulationParameters, int *status));
	LoadSimulation		   = DLLGetFunction(dll, void, C_LoadSimulationW,		  (HINSTANCE handle, LPCWSTR wcs, int *status));
	SaveSimulation		   = DLLGetFunction(dll, void, C_SaveSimulationW,		  (HINSTANCE handle, LPCWSTR wcs, int *status));
	GetTimeHistory2		   = DLLGetFunction(dll, void, C_GetTimeHistory2W,		  (HINSTANCE handle, void *nil, const TPeriod *period, int varID, double *lpValues, int *status));
	GetVarID			   = DLLGetFunction(dll, void, C_GetVarIDW,				  (HINSTANCE handle, LPCWSTR wcs, int *lpVarID, int *status));
	GetNumOfSamples		   = DLLGetFunction(dll, int,  C_GetNumOfSamples, 		  (HINSTANCE handle, const TPeriod *period, int *status));
	EnumerateObjects	   = DLLGetFunction(dll, void, C_EnumerateObjectsW, 	  (HINSTANCE handle, TEnumerateObjectsProc proc, int *lpNumOfObjects, int *status));
	EnumerateVars2		   = DLLGetFunction(dll, void, C_EnumerateVars2W, 		  (HINSTANCE handle, const TObjectExtra2 *objectextra, int ResultType, 
																						TEnumerateVarsProc EnumerateVarsProc, int *lpNumberOfVars, int *status));
	ObjectCalled		   = DLLGetFunction(dll, void, C_ObjectCalledW,			  (HINSTANCE handle, LPCWSTR lpObjectName, TObjectInfo *lpObjectInfo, int *status));
	CGetModelState		   = DLLGetFunction(dll, void, C_GetModelState,		  	  (HINSTANCE handle, int *lpModelState, int *status));

	GetDataType_	   	   = DLLGetFunction(dll, void, C_GetDataTypeW,		  	  (HINSTANCE handle, LPCWSTR name, int *lpType, int *status));
	GetDataRowCount_	   = DLLGetFunction(dll, void, C_GetDataRowCountW,		  (HINSTANCE handle, LPCWSTR name, int *lpCount, int *status));
	GetDataInteger		   = DLLGetFunction(dll, void, C_GetDataIntegerW,		  (HINSTANCE handle, LPCWSTR name, int index, int *lpData, int *status));
	GetDataDouble		   = DLLGetFunction(dll, void, C_GetDataDoubleW,		  (HINSTANCE handle, LPCWSTR name, int index, double *lpData, int *status));
	GetDataString		   = DLLGetFunction(dll, int,  C_GetDataStringW,		  (HINSTANCE handle, LPCWSTR name, int index, LPWSTR lpData, int *status));
	
	SetModelThreadCount = DLLGetFunction(dll, void, C_SetModelThreadCount, (HINSTANCE handle, int threadCount, int *status));
	GetModelThreadCount = DLLGetFunction(dll, int, C_GetModelThreadCount,  (HINSTANCE handle, int *status));
	
	GetDiffractionOutput0 = DLLGetFunction(dll, void, C_GetDiffractionOutput, (HINSTANCE handle, int OutputType, int *lpOutputSize, void *lpOutput, int *lpStatus));

	RegisterLicenceNotFoundHandler = DLLGetFunction(dll, void, C_RegisterLicenceNotFoundHandler, 	  (TLicenceNotFoundHandlerProc Handler, int *lpStatus));
	
	GetLastErrorString = DLLGetFunction(dll, int,  C_GetLastErrorStringW, (LPCWSTR wcs));
	FinaliseLibrary    = DLLGetFunction(dll, void, C_FinaliseLibrary,     (int *status));
	
	return true;
}

int Orca::GetDiffractionOutput(HINSTANCE handle, int OutputType, int *lpOutputSize, void *lpOutput) {
	int lpStatus;
	GetDiffractionOutput0(handle, OutputType, lpOutputSize, lpOutput, &lpStatus);
	return lpStatus;
}

void __stdcall Orca::DiffractionHandlerProc(HINSTANCE handle, LPCWSTR lpProgress, BOOL *lpCancel) {
	static int lastPerc = -1;
	String msg = WideToString(lpProgress);
	int perc = -1;
	int idp = msg.Find("%");
	if (idp >= 0) {
		int idbp = msg.ReverseFind(" ", idp);
		if (idbp >= 0) 
			perc = ScanInt(msg.Mid(idbp+1, idp-idbp));
	}
	if (lastPerc == perc) 
		return;
	
	lastPerc = perc;
	
	Time et;
	if (perc == 0)
		et = Null;
	else {
		int64 sec = GetSysTime() - startCalc - noLicenseTime;
		int64 estsec = int64(100*sec/double(perc));
		et = GetSysTime() + estsec;
	}
	*lpCancel = WhenWave(msg, perc, et);
}

void __stdcall Orca::StaticsHandlerProc(HINSTANCE handle, LPCWSTR lpProgress, BOOL *lpCancel) {
	String msg = WideToString(lpProgress);
	*lpCancel = WhenPrint(msg);
}

void __stdcall Orca::LicenceNotFoundHandler(int action, BOOL *lpAttemptReconnection, INT_PTR *lpData) {
	switch (action) {
	case lrBegin:
		beginNoLicense = GetSysTime();
		Sleep(60*1000); 
		*lpAttemptReconnection = TRUE;
		WhenPrint("License lost. Attemping reconnection in a minute");
		*lpData = 1;
		return;
	case lrContinue:
		*lpAttemptReconnection = TRUE;	//*lpData < 10; 
		if (*lpAttemptReconnection) {
			Sleep(60*1000);
			WhenPrint(Format("License lost for %d min. Attemping reconnection in a minute", *lpData));
			(*lpData)++; 
		}
		return;
	case lrEnd:
		noLicenseTime += (GetSysTime() - beginNoLicense);
	}
}

void __stdcall Orca::SimulationHandlerProc(HINSTANCE handle, double simulationTime, double simulationStart, double simulationStop, BOOL *lpCancel) {
	Time tm = GetSysTime();
	if (IsNull(startCalc))		// Time starts here. Statics calculation delay is discarded
		startCalc = tm;
	else if (tm - lastLog < deltaLogSimulation)
		return;
	lastLog = tm;
	double elapsed  = simulationTime - simulationStart,
		   total    = simulationStop - simulationStart;
	int64  elapsedT = tm - startCalc - noLicenseTime,
		   pending  = int64(elapsedT*(total/elapsed - 1));		// elapsedT*total/elapsed - elapsedT
	*lpCancel = WhenPrint(Format("Elap/Total:%.1f/%.0f ET:%s Clk/Sim:%.1f", elapsed, total,
										   					  SecondsToString(double(pending), 0, false, false, true, false, true), elapsedT/elapsed));
}

void __stdcall Orca::EnumerateObjectsProc(HINSTANCE handle, const TObjectInfo *info) {
	objTypes << info->ObjectType;
	WString str(info->ObjectName);
	objNames << str.ToString();
	objHandles << info->ObjectHandle;
}

void __stdcall Orca::EnumerateVarsProc(const TVarInfo *lpVarInfo) {
	int id = Find(varIDs, lpVarInfo->VarID);
	if (id < 0 || (actualBlade > 0 && varBlades[id] != actualBlade)) {
		varIDs << lpVarInfo->VarID;
		varBlades << actualBlade;
		varNames << WideToString(lpVarInfo->lpVarName);
		varFullNames << WideToString(lpVarInfo->lpFullName);
		varUnits << WideToString(lpVarInfo->lpVarUnits);
	}
}
			
bool Orca::FindInit() {
	UArray<SoftwareDetails> orcadata = GetSoftwareDetails("*OrcaFlex*");	// Get installed versions
	if (orcadata.IsEmpty()) {
		BEM::PrintWarning("OrcaFlex is not installed");
		return false;
	}
	
	int iversion = 0;
	UVector<int> version = orcadata[0].GetVersion();
	for (int i = 1; i < orcadata.size(); ++i) {							// Get the newest version from the installed
		UVector<int> each = orcadata[i].GetVersion();
		if (SoftwareDetails::IsHigherVersion(each, version)) {			// Version are numbers separated by .
			iversion = i;
			version = pick(each);
		}
	}
	String arch;
#ifdef CPU_64
	arch = "Win64";
#else
	arch = "Win32";
#endif
	String path = AFX(orcadata[iversion].path, "OrcFxAPI", arch, "OrcFxAPI.dll");	// Assembles the path
	
	return Init(path);
}
		
int Orca::GetDataType(HINSTANCE handle, const wchar_t *name) {
	int status;
	int ret;
	
	GetDataType_(handle, name, &ret, &status);
	if (status != 0)
		throwError(Format("Load GetDataType %s", name));	
	return ret;
}

int Orca::GetDataRowCount(HINSTANCE handle, const wchar_t *name) {
	int status;
	int ret;
	
	GetDataRowCount_(handle, name, &ret, &status);
	if (status != 0)
		throwError(Format("Load GetDataRowCount %s", name));	
	return ret;
}

int Orca::GetInt(HINSTANCE handle, const wchar_t *name, int id) {
	int status;
	int ret;
	
	GetDataInteger(handle, name, id, &ret, &status);
	if (status != 0)
		throwError(Format("Load GetDataInteger %s", name));	
	return ret;
}

double Orca::GetDouble(HINSTANCE handle, const wchar_t *name, int id) {
	int status;
	double ret;
	
	GetDataDouble(handle, name, id, &ret, &status);
	if (status != 0)
		throwError(Format("Load GetDataDouble %s", name));	
	return ret;
}
	
String Orca::GetString(HINSTANCE handle, const wchar_t *name, int id) {
	int status;
		
	int len = GetDataString(wave, name, -1, NULL, &status);
	if (status != 0)
		throwError(Format("Load GetDataString %s", name));	
	Buffer<wchar_t> rw(len);
	LPWSTR wcs = (LPWSTR)rw.begin();
	GetDataString(wave, name, -1, wcs, &status);
	if (status != 0)
		throwError(Format("Load GetDataString 2 %s", name));	
	return WideToString(wcs, len);
}
		
void Orca::LoadParameters(Hydro &hy) {
	int sz;
	
	OrcaFactors factor;	
	
	factor.len = FactorLen(GetString(wave, L"LengthUnits"));
	factor.mass = FactorMass(GetString(wave, L"MassUnits"));
	factor.force = FactorForce(GetString(wave, L"ForceUnits"));

	factor.Update();

	hy.dt.symX = hy.dt.symY = false;
	
	hy.dt.g = GetDouble(wave, L"g")*factor.len;

	hy.dt.h = GetDouble(wave, L"WaterDepth");
	if (hy.dt.h == 1e307)
		hy.dt.h = -1;
	else
		hy.dt.h *= factor.len;
	
	hy.dt.rho = GetDouble(wave, L"WaterDensity")*factor.mass/factor.len/factor.len/factor.len;
			
	if (GetDiffractionOutput(wave, dotAngularFrequencies, &sz, NULL))
		throwError("Load dotAngularFrequencies");	
	
	hy.dt.Nf = sz/sizeof(double);
	hy.dt.w.SetCount(hy.dt.Nf);
	if (GetDiffractionOutput(wave, dotAngularFrequencies, &sz, hy.dt.w.begin()))
		throwError("Load dotAngularFrequencies 2");	
		
	/*hy.dt.dt.T.SetCount(hy.dt.dt.Nf);
	for (int i = 0; i < hy.dt.dt.Nf; ++i)	
		hy.dt.dt.T[i] = 2*M_PI/hy.dt.dt.w[i];	*/
	
	if (GetDiffractionOutput(wave, dotHeadings, &sz, NULL))
		throwError("Load dotHeadings");	
	
	hy.dt.Nh = sz/sizeof(double);
	hy.dt.head.SetCount(hy.dt.Nh);
	if (GetDiffractionOutput(wave, dotHeadings, &sz, hy.dt.head.begin()))
		throwError("Load dotHeadings 2");	
	
	if (GetDiffractionOutput(wave, dotHydrostaticResults, &sz, NULL))
		throwError("Load dotHydrostaticResults");			
	
	hy.dt.Nb = GetInt(wave, L"NumberOfIncludedBodies");
	
	if (sz/hy.dt.Nb != sizeof(TDiffractionBodyHydrostaticInfo))
		throwError("Incompatible OrcaFlex version. TDiffractionBodyHydrostaticInfo size does not match");			
	
	hy.dt.Nb = sz/sizeof(TDiffractionBodyHydrostaticInfo);
	Buffer<TDiffractionBodyHydrostaticInfo> bodies(hy.dt.Nb);
	if (GetDiffractionOutput(wave, dotHydrostaticResults, &sz, bodies.begin()))
		throwError("Load dotHydrostaticResults 2");
	
	hy.dt.msh.SetCount(hy.dt.Nb);
	//hy.dt.dt.Vo.SetCount(hy.dt.dt.Nb);
	//hy.dt.dt.cb.resize(3, hy.dt.dt.Nb);
	//hy.dt.dt.cg.resize(3, hy.dt.dt.Nb);
	//hy.dt.dt.C.SetCount(hy.dt.dt.Nb);
	for (int ib = 0; ib < hy.dt.Nb; ++ib) 
		hy.dt.msh[ib].dt.C.resize(6, 6);
	//hy.dt.dt.M.SetCount(hy.dt.dt.Nb);
	for (int ib = 0; ib < hy.dt.Nb; ++ib) 
		hy.dt.msh[ib].dt.M.resize(6, 6);
	
	//hy.dt.dt.c0.setConstant(3, hy.dt.dt.Nb, 0);
	//hy.dt.dt.names.SetCount(hy.dt.dt.Nb);
	
	for (int ib = 0; ib < hy.dt.Nb; ++ib) {
		const TDiffractionBodyHydrostaticInfo &b = bodies[ib];
		
		hy.dt.msh[ib].dt.Vo = b.Volume*factor.len*factor.len*factor.len;
		for (int idf = 0; idf < 3; ++idf) {
			hy.dt.msh[ib].dt.cb[idf] = b.CentreOfBuoyancy[idf]*factor.len;
			hy.dt.msh[ib].dt.cg[idf] = b.CentreOfMass[idf]*factor.len;
		}
		for (int r = 0; r < 6; ++r) {
			for (int c = 0; c < 6; ++c) {
				hy.dt.msh[ib].dt.C(r, c) = b.RestoringMatrix[r][c]*factor.K(r, c);
				hy.dt.msh[ib].dt.M(r, c) = b.InertiaMatrix[r][c]*factor.M(r, c);
			}
		}
		hy.dt.msh[ib].dt.name = FormatInt(ib+1);
	}
	
	auto LoadAB = [&](UArray<UArray<VectorXd>> &ab, int type, const char *stype, const Matrix<double, 6, 6> &factor) {
		hy.Initialize_AB(ab);
	
		if (GetDiffractionOutput(wave, type, &sz, NULL))
			throwError("Load dotAddedMass_Radiation");	
		
		if (sz/sizeof(double) != 6*hy.dt.Nb*6*hy.dt.Nb*hy.dt.Nf)		
			throw Exc(Format("Wrong %s size (%d <> %d)", stype, int(sz/sizeof(double)), 6*hy.dt.Nb*6*hy.dt.Nb*hy.dt.Nf));
		
		MultiDimMatrixRowMajor<double> a(hy.dt.Nf, 6*hy.dt.Nb, 6*hy.dt.Nb);
		if (GetDiffractionOutput(wave, type, &sz, a.begin()))
			throwError("Load dotAddedMass_Radiation 2");		
		
		for (int r = 0; r < 6*hy.dt.Nb; ++r)
			for (int c = 0; c < 6*hy.dt.Nb; ++c)
				for (int ifr = 0; ifr < hy.dt.Nf; ++ifr)
		 			ab[r][c][ifr] = a(ifr, r, c)*factor(r%6, c%6);
	};
	
	LoadAB(hy.dt.A, dotAddedMass, "added mass", factor.A);
	LoadAB(hy.dt.B, dotDamping, "radiation damping", factor.B);
	
	if (GetDiffractionOutput(wave, dotInfiniteFrequencyAddedMass, &sz, NULL))
		throwError("Load dotInfiniteFrequencyAddedMass");	
	
	if (sz/sizeof(double) != 6*hy.dt.Nb*6*hy.dt.Nb)		
		throw Exc(Format("Wrong %s size (%d <> %d)", "infinite frequency added mass", int(sz/sizeof(double)), 6*hy.dt.Nb*6*hy.dt.Nb));
	
	MultiDimMatrixRowMajor<double> a(6*hy.dt.Nb, 6*hy.dt.Nb);
	if (GetDiffractionOutput(wave, dotInfiniteFrequencyAddedMass, &sz, a.begin()))
		throwError("Load dotInfiniteFrequencyAddedMass 2");		
	
	hy.dt.Ainf.resize(6*hy.dt.Nb, 6*hy.dt.Nb);
	for (int r = 0; r < 6*hy.dt.Nb; ++r)
		for (int c = 0; c < 6*hy.dt.Nb; ++c)
	 		hy.dt.Ainf(r, c) = a(r, c)*factor.A(r%6, c%6);	
	
	
	auto LoadF = [&](Hydro::Forces &f, int type, const char *stype, const Eigen::Vector<double, 6> &factor) {
		hy.Initialize_Forces(f);
		
		if (GetDiffractionOutput(wave, type, &sz, NULL))
			throwError("Load dotLoadRAOs");	
		
		if (sz/sizeof(TComplex) != 6*hy.dt.Nb*hy.dt.Nf*hy.dt.Nh)		
			throw Exc(Format("Wrong %s size (%d <> %d)", stype, int(sz/sizeof(TComplex)), 6*hy.dt.Nb*6*hy.dt.Nb*hy.dt.Nf));
		
		MultiDimMatrixRowMajor<TComplex> a(hy.dt.Nh, hy.dt.Nf, 6*hy.dt.Nb);
		if (GetDiffractionOutput(wave, type, &sz, a.begin()))
			throwError("Load dotLoadRAOs 2");		
		
		for (int r = 0; r < 6*hy.dt.Nb; ++r) {
			for (int ih = 0; ih < hy.dt.Nh; ++ih) {
				for (int ifr = 0; ifr < hy.dt.Nf; ++ifr) {
					f.force[ih](ifr, r).real(a(ih, ifr, r).Re*factor(r%6));
					f.force[ih](ifr, r).imag(a(ih, ifr, r).Im*factor(r%6));
				}
			}
		}
	};
	
	LoadF(hy.dt.ex, dotLoadRAOsDiffraction, "diffraction force", factor.F);
	
	LoadF(hy.dt.rao, dotDisplacementRAOs, "RAO", factor.RAO);
	

	if (GetDiffractionOutput(wave, dotQTFAngularFrequencies, &sz, NULL))
		throwError("Load dotAngularFrequencies");	

	Buffer<double> qfreq(sz/sizeof(double));
	if (GetDiffractionOutput(wave, dotQTFAngularFrequencies, &sz, qfreq))
		throwError("Load dotQTFAngularFrequencies 2");	
	
	sz /= sizeof(double);
	int Nqw = sz/3;
	UVector<double> qw;
	for (int i = 0; i < sz; i += 3) 
		FindAdd(qw, qfreq[i]);
	Sort(qw);
	Copy(qw, hy.dt.qw);
	
	if (GetDiffractionOutput(wave, dotMeanDriftHeadingPairs, &sz, NULL))
		throwError("Load dotMeanDriftHeadingPairs");	
	
	Buffer<double> qmh(sz/sizeof(double));
	if (GetDiffractionOutput(wave, dotMeanDriftHeadingPairs, &sz, qmh.begin()))
		throwError("Load dotMeanDriftHeadingPairs 2");	

	sz /= (2*sizeof(double));
	hy.dt.mdhead.resize(sz);
	for (int i = 0; i < sz; i++) 
		hy.dt.mdhead[i] = std::complex<double>(qmh[2*i], qmh[2*i+1]);
	
	auto LoadMD = [&](int type, const char *stype)->bool {
		if (GetDiffractionOutput(wave, type, &sz, NULL))
			throwError(Format("Load %s", stype));	
	
		if (sz == 0)
			return false;
		
		if (sz/sizeof(TComplex) != 6*hy.dt.Nb*hy.dt.Nf*hy.dt.mdhead.size())		
			throw Exc(Format("Wrong %s size (%d <> %d)", stype, int(sz/sizeof(TComplex)), 6*hy.dt.Nb*hy.dt.Nf*hy.dt.mdhead.size()));
		
		MultiDimMatrixRowMajor<TComplex> md((int)hy.dt.mdhead.size(), hy.dt.Nf, 6*hy.dt.Nb);
		if (GetDiffractionOutput(wave, type, &sz, md.begin()))
			throwError(Format("Load %s 2", stype));		
		
		Hydro::Initialize_MD(hy.dt.md, hy.dt.Nb, int(hy.dt.mdhead.size()), hy.dt.Nf);
		
		int Nb = hy.dt.Nb;
		if (type == dotMeanDriftLoadControlSurface)
			hy.dt.mdtype = 7;
		else if (type == dotMeanDriftLoadPressureIntegration)
			hy.dt.mdtype = 9;
		else {
			hy.dt.mdtype = 8;
			Nb = 1;			// Momentum Conservation handles only 1 body
		}
		
		for (int ib = 0; ib < Nb; ++ib) 
			for (int ih = 0; ih < hy.dt.mdhead.size(); ++ih) 
				for (int idf = 0; idf < 6; ++idf) 
					for (int ifr = 0; ifr < hy.dt.Nf; ++ifr) {
						double re = md(ih, ifr, 6*ib+idf).Re;
						double im = md(ih, ifr, 6*ib+idf).Im;
						hy.dt.md[ib][ih][idf](ifr) = Sign(re)*sqrt(sqr(re) + sqr(im))*factor.MD(idf);
					}
					
		return true;
	};
	
	if (!LoadMD(dotMeanDriftLoadControlSurface, "Mean Drift Control Surface"))
		if (!LoadMD(dotMeanDriftLoadPressureIntegration, "Mean Drift Pressure Integration"))
			LoadMD(dotMeanDriftLoadMomentumConservation, "Mean Drift Momentum Conservation");

	if (GetDiffractionOutput(wave, dotQTFHeadingPairs, &sz, NULL))
		throwError("Load dotQTFHeadingPairs");	
				
	Buffer<double> qh(sz/sizeof(double));
	if (GetDiffractionOutput(wave, dotQTFHeadingPairs, &sz, qh.begin()))
		throwError("Load dotQTFHeadingPairs 2");	
	
	sz /= (2*sizeof(double));
	hy.dt.qh.resize(sz);
	for (int i = 0; i < sz; i++) {
		hy.dt.qh[i].real(qh[2*i]);
		hy.dt.qh[i].imag(qh[2*i+1]);
	}
	
	auto LoadQTF = [&](int type, const char *stype)->bool {
		MultiDimMatrixRowMajor<TComplex> qtf, qtfDirect;
		{
			if (GetDiffractionOutput(wave, type, &sz, NULL))
				throwError(Format("Load %s", stype));	
		
			if (sz == 0)
				return false;
			
			if (sz/sizeof(TComplex) != 6*hy.dt.Nb*Nqw*hy.dt.qh.size())		
				throw Exc(Format("Wrong %s size (%d <> %d)", stype, int(sz/sizeof(TComplex)), 6*hy.dt.Nb*Nqw*hy.dt.qh.size()));
			
			qtf.Resize((int)hy.dt.qh.size(), Nqw, 6*hy.dt.Nb);
			if (GetDiffractionOutput(wave, type, &sz, qtf.begin()))
				throwError(Format("Load %s 2", stype));		
		}
		{
			if (GetDiffractionOutput(wave, dotDirectPotentialLoad, &sz, NULL))
				throwError(Format("Load %s", stype));	
		
			if (sz == 0)
				return false;
			
			if (sz/sizeof(TComplex) != 6*hy.dt.Nb*Nqw*hy.dt.qh.size())		
				throw Exc(Format("Wrong %s size (%d <> %d)", stype, int(sz/sizeof(TComplex)), 6*hy.dt.Nb*Nqw*hy.dt.qh.size()));
			
			qtfDirect.Resize((int)hy.dt.qh.size(), Nqw, 6*hy.dt.Nb);
			if (GetDiffractionOutput(wave, dotDirectPotentialLoad, &sz, qtfDirect.begin()))
				throwError(Format("Load %s 2", stype));		
		}
		int Nb = hy.dt.Nb;
		
		for (int ib = 0; ib < Nb; ++ib) 
			for (int ih = 0; ih < hy.dt.qh.size(); ++ih) 
				for (int idf = 0; idf < 6; ++idf) 
					for (int ifr = 0, iNqw = 0; ifr < 3*Nqw; ifr += 3, iNqw++) {
						qtf(ih, iNqw, 6*ib + idf).Re += qtfDirect(ih, iNqw, 6*ib + idf).Re;
						qtf(ih, iNqw, 6*ib + idf).Im += qtfDirect(ih, iNqw, 6*ib + idf).Im;
					}
					
		Hydro::Initialize_QTF(hy.dt.qtfsum, hy.dt.Nb, int(hy.dt.qh.size()), int(hy.dt.qw.size()));
		Hydro::Initialize_QTF(hy.dt.qtfdif, hy.dt.Nb, int(hy.dt.qh.size()), int(hy.dt.qw.size()));
		
		if (type == dotQuadraticLoadFromControlSurface)
			hy.dt.qtftype = 7;
		else if (type == dotQuadraticLoadFromPressureIntegration)
			hy.dt.qtftype = 9;
	
		for (int ib = 0; ib < Nb; ++ib) 
			for (int ih = 0; ih < hy.dt.qh.size(); ++ih) 
				for (int idf = 0; idf < 6; ++idf) 
					for (int ifr = 0, iNqw = 0; ifr < 3*Nqw; ifr += 3, iNqw++) {
						double freq1 = qfreq[ifr];
						double freq2 = qfreq[ifr+1];
						double freq12= qfreq[ifr+2];
						int ifr1 = Find(hy.dt.qw, freq1);
						int ifr2 = Find(hy.dt.qw, freq2);
						if (abs(freq1 + freq2 - freq12) < 0.001)
							hy.dt.qtfsum[ib][ih][idf](ifr1, ifr2) = std::complex<double>(qtf(ih, iNqw, 6*ib+idf).Re, 
																					  qtf(ih, iNqw, 6*ib+idf).Im)*factor.F(idf);
						else 
							hy.dt.qtfdif[ib][ih][idf](ifr1, ifr2) = std::complex<double>(qtf(ih, iNqw, 6*ib+idf).Re, 
																					  qtf(ih, iNqw, 6*ib+idf).Im)*factor.F(idf);
					}

			
		return true;
	};
		
	if (!LoadQTF(dotQuadraticLoadFromControlSurface, "QTF Control Surface"))
		LoadQTF(dotQuadraticLoadFromPressureIntegration, "QTF Pressure Integration");

	int Np;
	if (GetDiffractionOutput(wave, dotPanelCount, &sz, NULL))
		throwError("Load dotPanelCount");
	if (GetDiffractionOutput(wave, dotPanelCount, &sz, &Np))
		throwError("Load dotPanelCount 2");
	
	if (GetDiffractionOutput(wave, dotPanelGeometry, &sz, NULL))
		throwError("Load dotPanelGeometry");			
		
	if (sz/Np != sizeof(TDiffractionPanelGeometry))
		throwError("Incompatible OrcaFlex version. TDiffractionPanelGeometry size does not match");			
	
	Buffer<TDiffractionPanelGeometry> panels(Np);
	if (GetDiffractionOutput(wave, dotPanelGeometry, &sz, panels.begin()))
		throwError("Load dotPanelGeometry 2");
	
	hy.dt.msh.SetCount(hy.dt.Nb);
	for (int ib = 0; ib < hy.dt.Nb; ++ib) {
		hy.dt.msh[ib].dt.SetCode(Mesh::ORCA_OWR);
		hy.dt.msh[ib].dt.c0 = Point3D(0, 0, 0);
		//hy.dt.msh[ib].cg = Point3D(hy.dt.cg(0, ib), hy.dt.cg(1, ib), hy.dt.cg(2, ib));
	}
	for (int ip = 0; ip < Np; ++ip) {
		const TDiffractionPanelGeometry &pan = panels[ip];
		int ib = pan.ObjectId;
		
		if (ib < 0)		// Artificial damping lid
			continue;
		
		Panel &p = hy.dt.msh[ib].dt.mesh.panels.Add();
		
		for (int i = 0; i < 4; ++i) {
			const TVector &v0 = pan.Vertices[i];
			
			if (i == 3 && std::isnan<double>(v0[0]))
				p.id[i] = p.id[0];
			else {
				Point3D pnt(v0[0]*factor.len, v0[1]*factor.len, v0[2]*factor.len);
				p.id[i] = FindAdd(hy.dt.msh[ib].dt.mesh.nodes, pnt);
			}
		}
	}
	
	if (GetDiffractionOutput(wave, dotPanelPressureRadiation, &sz, NULL))
		throwError("Load dotPanelPressureRadiation");			
	
	if (sz > 0) {
		if (sz/sizeof(TComplex) != 6*hy.dt.Nb*hy.dt.Nf*Np)
			throw Exc(Format("Wrong %s size (%d <> %d)", "dotPanelPressureRadiation", int(sz/sizeof(TComplex)), 6*hy.dt.Nb*hy.dt.Nf*Np));
					
		MultiDimMatrixRowMajor<TComplex> pres(6*hy.dt.Nb, hy.dt.Nf, Np);
		if (GetDiffractionOutput(wave, dotPanelPressureRadiation, &sz, pres.begin()))
			throwError("Load dotPanelPressureRadiation 2");
		
		hy.Initialize_PotsRad();
		for (int ib = 0; ib < hy.dt.Nb; ++ib)
			for (int idf = 0; idf < 6; ++idf)
				for (int ifr = 0; ifr < hy.dt.Nf; ++ifr)
					for (int ip = 0; ip < Np; ++ip) {
						const TComplex &c = pres(idf + 6*ib, ifr, ip);
						hy.dt.pots_rad[ib][ip][idf][ifr] = std::complex<double>(-c.Re, c.Im)/hy.dt.rho;
							//std::complex<double>(c.Im, -c.Re)/hy.dt.rho/hy.dt.w[ifr]; // press = i w rho phi(potential)
					}
	}
}
	
#endif