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

	hy.symX = hy.symY = false;
	
	hy.g = GetDouble(wave, L"g")*factor.len;

	hy.h = GetDouble(wave, L"WaterDepth");
	if (hy.h == 1e307)
		hy.h = -1;
	else
		hy.h *= factor.len;
	
	hy.rho = GetDouble(wave, L"WaterDensity")*factor.mass/factor.len/factor.len/factor.len;
			
	if (GetDiffractionOutput(wave, dotAngularFrequencies, &sz, NULL))
		throwError("Load dotAngularFrequencies");	
	
	hy.Nf = sz/sizeof(double);
	hy.w.SetCount(hy.Nf);
	if (GetDiffractionOutput(wave, dotAngularFrequencies, &sz, hy.w.begin()))
		throwError("Load dotAngularFrequencies 2");	
		
	hy.T.SetCount(hy.Nf);
	for (int i = 0; i < hy.Nf; ++i)	
		hy.T[i] = 2*M_PI/hy.w[i];	
	
	if (GetDiffractionOutput(wave, dotHeadings, &sz, NULL))
		throwError("Load dotHeadings");	
	
	hy.Nh = sz/sizeof(double);
	hy.head.SetCount(hy.Nh);
	if (GetDiffractionOutput(wave, dotHeadings, &sz, hy.head.begin()))
		throwError("Load dotHeadings 2");	
	
	if (GetDiffractionOutput(wave, dotHydrostaticResults, &sz, NULL))
		throwError("Load dotHydrostaticResults");			
	
	hy.Nb = GetInt(wave, L"NumberOfIncludedBodies");
	
	if (sz/hy.Nb != sizeof(TDiffractionBodyHydrostaticInfo))
		throwError("Incompatible OrcaFlex version. TDiffractionBodyHydrostaticInfo size does not match");			
	
	hy.Nb = sz/sizeof(TDiffractionBodyHydrostaticInfo);
	Buffer<TDiffractionBodyHydrostaticInfo> bodies(hy.Nb);
	if (GetDiffractionOutput(wave, dotHydrostaticResults, &sz, bodies.begin()))
		throwError("Load dotHydrostaticResults 2");
	
	hy.Vo.SetCount(hy.Nb);
	hy.cb.resize(3, hy.Nb);
	hy.cg.resize(3, hy.Nb);
	hy.C.SetCount(hy.Nb);
	for (int ib = 0; ib < hy.Nb; ++ib) 
		hy.C[ib].resize(6, 6);
	hy.M.SetCount(hy.Nb);
	for (int ib = 0; ib < hy.Nb; ++ib) 
		hy.M[ib].resize(6, 6);
	
	hy.c0.setConstant(3, hy.Nb, 0);
	hy.names.SetCount(hy.Nb);
	
	for (int ib = 0; ib < hy.Nb; ++ib) {
		const TDiffractionBodyHydrostaticInfo &b = bodies[ib];
		
		hy.Vo[ib] = b.Volume*factor.len*factor.len*factor.len;
		for (int idf = 0; idf < 3; ++idf) {
			hy.cb(idf, ib) = b.CentreOfBuoyancy[idf]*factor.len;
			hy.cg(idf, ib) = b.CentreOfMass[idf]*factor.len;
		}
		for (int r = 0; r < 6; ++r) {
			for (int c = 0; c < 6; ++c) {
				hy.C[ib](r, c) = b.RestoringMatrix[r][c]*factor.K(r, c);
				hy.M[ib](r, c) = b.InertiaMatrix[r][c]*factor.M(r, c);
			}
		}
		hy.names[ib] = FormatInt(ib+1);
	}
	
	auto LoadAB = [&](UArray<UArray<VectorXd>> &ab, int type, const char *stype, const Matrix<double, 6, 6> &factor) {
		hy.Initialize_AB(ab);
	
		if (GetDiffractionOutput(wave, type, &sz, NULL))
			throwError("Load dotAddedMass_Radiation");	
		
		if (sz/sizeof(double) != 6*hy.Nb*6*hy.Nb*hy.Nf)		
			throw Exc(Format("Wrong %s size (%d <> %d)", stype, int(sz/sizeof(double)), 6*hy.Nb*6*hy.Nb*hy.Nf));
		
		MultiDimMatrixRowMajor<double> a(hy.Nf, 6*hy.Nb, 6*hy.Nb);
		if (GetDiffractionOutput(wave, type, &sz, a.begin()))
			throwError("Load dotAddedMass_Radiation 2");		
		
		for (int r = 0; r < 6*hy.Nb; ++r)
			for (int c = 0; c < 6*hy.Nb; ++c)
				for (int ifr = 0; ifr < hy.Nf; ++ifr)
		 			ab[r][c][ifr] = a(ifr, r, c)*factor(r%6, c%6);
	};
	
	LoadAB(hy.A, dotAddedMass, "added mass", factor.A);
	LoadAB(hy.B, dotDamping, "radiation damping", factor.B);
	
	if (GetDiffractionOutput(wave, dotInfiniteFrequencyAddedMass, &sz, NULL))
		throwError("Load dotInfiniteFrequencyAddedMass");	
	
	if (sz/sizeof(double) != 6*hy.Nb*6*hy.Nb)		
		throw Exc(Format("Wrong %s size (%d <> %d)", "infinite frequency added mass", int(sz/sizeof(double)), 6*hy.Nb*6*hy.Nb));
	
	MultiDimMatrixRowMajor<double> a(6*hy.Nb, 6*hy.Nb);
	if (GetDiffractionOutput(wave, dotInfiniteFrequencyAddedMass, &sz, a.begin()))
		throwError("Load dotInfiniteFrequencyAddedMass 2");		
	
	hy.Ainf.resize(6*hy.Nb, 6*hy.Nb);
	for (int r = 0; r < 6*hy.Nb; ++r)
		for (int c = 0; c < 6*hy.Nb; ++c)
	 		hy.Ainf(r, c) = a(r, c)*factor.A(r%6, c%6);	
	
	
	auto LoadF = [&](Hydro::Forces &f, int type, const char *stype, const Eigen::Vector<double, 6> &factor) {
		hy.Initialize_Forces(f);
		
		if (GetDiffractionOutput(wave, type, &sz, NULL))
			throwError("Load dotLoadRAOs");	
		
		if (sz/sizeof(TComplex) != 6*hy.Nb*hy.Nf*hy.Nh)		
			throw Exc(Format("Wrong %s size (%d <> %d)", stype, int(sz/sizeof(TComplex)), 6*hy.Nb*6*hy.Nb*hy.Nf));
		
		MultiDimMatrixRowMajor<TComplex> a(hy.Nh, hy.Nf, 6*hy.Nb);
		if (GetDiffractionOutput(wave, type, &sz, a.begin()))
			throwError("Load dotLoadRAOs 2");		
		
		for (int r = 0; r < 6*hy.Nb; ++r) {
			for (int ih = 0; ih < hy.Nh; ++ih) {
				for (int ifr = 0; ifr < hy.Nf; ++ifr) {
					f.force[ih](ifr, r).real(a(ih, ifr, r).Re*factor(r%6));
					f.force[ih](ifr, r).imag(a(ih, ifr, r).Im*factor(r%6));
				}
			}
		}
	};
	
	LoadF(hy.ex, dotLoadRAOsDiffraction, "diffraction force", factor.F);
	
	LoadF(hy.rao, dotDisplacementRAOs, "RAO", factor.RAO);
	

	if (GetDiffractionOutput(wave, dotQTFAngularFrequencies, &sz, NULL))
		throwError("Load dotAngularFrequencies");	

	Buffer<double> qfreq(sz/sizeof(double));
	if (GetDiffractionOutput(wave, dotQTFAngularFrequencies, &sz, qfreq))
		throwError("Load dotQTFAngularFrequencies 2");	
	
	int Nqw = sz /= sizeof(double);
	UVector<double> qw;
	for (int i = 0; i < sz; i += 3) 
		FindAdd(qw, qfreq[i]);
	Copy(qw, hy.qw);
	
	if (GetDiffractionOutput(wave, dotMeanDriftHeadingPairs, &sz, NULL))
		throwError("Load dotMeanDriftHeadingPairs");	
	
	Buffer<double> qmh(sz/sizeof(double));
	if (GetDiffractionOutput(wave, dotMeanDriftHeadingPairs, &sz, qmh.begin()))
		throwError("Load dotMeanDriftHeadingPairs 2");	

	sz /= (2*sizeof(double));
	hy.mdhead.resize(sz);
	for (int i = 0; i < sz; i++) {
		hy.mdhead[i].real(qmh[2*i]);
		hy.mdhead[i].imag(qmh[2*i+1]);
	}
	
	auto LoadMD = [&](int type, const char *stype)->bool {
		if (GetDiffractionOutput(wave, type, &sz, NULL))
			throwError(Format("Load %s", stype));	
	
		if (sz == 0)
			return false;
		
		if (sz/sizeof(TComplex) != 6*hy.Nb*hy.Nf*hy.mdhead.size())		
			throw Exc(Format("Wrong %s size (%d <> %d)", stype, int(sz/sizeof(TComplex)), 6*hy.Nb*hy.Nf*hy.mdhead.size()));
		
		MultiDimMatrixRowMajor<TComplex> md(hy.mdhead.size(), hy.Nf, 6*hy.Nb);
		if (GetDiffractionOutput(wave, type, &sz, md.begin()))
			throwError(Format("Load %s 2", stype));		
		
		Hydro::Initialize_MD(hy.md, hy.Nb, int(hy.mdhead.size()), hy.Nf);
		
		int Nb = hy.Nb;
		if (type == dotMeanDriftLoadControlSurface)
			hy.mdtype = 7;
		else if (type == dotMeanDriftLoadPressureIntegration)
			hy.mdtype = 9;
		else {
			hy.mdtype = 8;
			Nb = 1;			// Momentum Conservation handles only 1 body
		}
		
		for (int ib = 0; ib < Nb; ++ib) 
			for (int ih = 0; ih < hy.mdhead.size(); ++ih) 
				for (int idf = 0; idf < 6; ++idf) 
					for (int ifr = 0; ifr < hy.Nf; ++ifr) {
						double re = md(ih, ifr, 6*ib+idf).Re;
						double im = md(ih, ifr, 6*ib+idf).Im;
						hy.md[ib][ih][idf](ifr) = Sign(re)*sqrt(sqr(re) + sqr(im))*factor.MD(idf);
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
	hy.qh.resize(sz);
	for (int i = 0; i < sz; i++) {
		hy.qh[i].real(qh[2*i]);
		hy.qh[i].imag(qh[2*i+1]);
	}
	
	auto LoadQTF = [&](int type, const char *stype)->bool {
		if (GetDiffractionOutput(wave, type, &sz, NULL))
			throwError(Format("Load %s", stype));	
	
		if (sz == 0)
			return false;
		
		if (sz/sizeof(TComplex) != 6*hy.Nb*(Nqw/3)*hy.qh.size())		
			throw Exc(Format("Wrong %s size (%d <> %d)", stype, int(sz/sizeof(TComplex)), 6*hy.Nb*(Nqw/3)*hy.qh.size()));
		
		MultiDimMatrixRowMajor<TComplex> qtf(hy.qh.size(), Nqw, 6*hy.Nb);
		if (GetDiffractionOutput(wave, type, &sz, qtf.begin()))
			throwError(Format("Load %s 2", stype));		
		
		Hydro::Initialize_QTF(hy.qtfsum, hy.Nb, int(hy.qh.size()), hy.qw.size());
		Hydro::Initialize_QTF(hy.qtfdif, hy.Nb, int(hy.qh.size()), hy.qw.size());
		
		int Nb = hy.Nb;
		if (type == dotQuadraticLoadFromControlSurface)
			hy.qtftype = 7;
		else if (type == dotQuadraticLoadFromPressureIntegration)
			hy.qtftype = 9;
	
		
		for (int ib = 0; ib < Nb; ++ib) 
			for (int ih = 0; ih < hy.qh.size(); ++ih) 
				for (int idf = 0; idf < 6; ++idf) 
					for (int ifr = 0; ifr < Nqw; ifr += 3) {
						int ifr1 = Find(hy.qw, qfreq[ifr]);
						if (ifr1 < 0)
							throw Exc(Format("Frequency %.3f not found", qfreq[ifr]));
						int ifr2 = Find(hy.qw, qfreq[ifr+1]);
						if (ifr2 < 0)
							throw Exc(Format("Frequency %.3f not found", qfreq[ifr]));
						if (abs(qfreq[ifr] + qfreq[ifr+1] - qfreq[ifr+2]) < 0.001)
							hy.qtfsum[ib][ih][idf](ifr1, ifr2) = std::complex<double>(qtf(ih, ifr, 6*ib+idf).Re, qtf(ih, ifr, 6*ib+idf).Im)*factor.F(idf);
						else
							hy.qtfdif[ib][ih][idf](ifr1, ifr2) = std::complex<double>(qtf(ih, ifr, 6*ib+idf).Re, qtf(ih, ifr, 6*ib+idf).Im)*factor.F(idf);
					}

			
		return true;
	};
		
	//if (!LoadQTF(dotQuadraticLoadFromControlSurface, "QTF Control Surface"))
	//	LoadQTF(dotQuadraticLoadFromPressureIntegration, "QTF Pressure Integration");

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
	
	hy.meshes.SetCount(hy.Nb);
	for (int ib = 0; ib < hy.Nb; ++ib) {
		hy.meshes[ib].SetCode(Mesh::ORCA_OWR);
		hy.meshes[ib].c0 = Point3D(hy.c0(0, ib), hy.c0(1, ib), hy.c0(2, ib));
		hy.meshes[ib].cg = Point3D(hy.cg(0, ib), hy.cg(1, ib), hy.cg(2, ib));
	}
	for (int ip = 0; ip < Np; ++ip) {
		const TDiffractionPanelGeometry &pan = panels[ip];
		int ib = pan.ObjectId;
		
		if (ib < 0)		// Artificial damping lid
			continue;
		
		Panel &p = hy.meshes[ib].mesh.panels.Add();
		
		for (int i = 0; i < 4; ++i) {
			const TVector &v0 = pan.Vertices[i];
			
			if (i == 3 && std::isnan<double>(v0[0]))
				p.id[i] = p.id[0];
			else {
				Point3D pnt(v0[0]*factor.len, v0[1]*factor.len, v0[2]*factor.len);
				p.id[i] = FindAdd(hy.meshes[ib].mesh.nodes, pnt);
			}
		}
	}
	
	if (GetDiffractionOutput(wave, dotPanelPressureRadiation, &sz, NULL))
		throwError("Load dotPanelPressureRadiation");			
	
	if (sz > 0) {
		if (sz/sizeof(TComplex) != 6*hy.Nb*hy.Nf*Np)
			throw Exc(Format("Wrong %s size (%d <> %d)", "dotPanelPressureRadiation", int(sz/sizeof(TComplex)), 6*hy.Nb*hy.Nf*Np));
					
		MultiDimMatrixRowMajor<TComplex> pres(6*hy.Nb, hy.Nf, Np);
		if (GetDiffractionOutput(wave, dotPanelPressureRadiation, &sz, pres.begin()))
			throwError("Load dotPanelPressureRadiation 2");
		
		hy.Initialize_Pots();
		for (int ib = 0; ib < hy.Nb; ++ib)
			for (int idf = 0; idf < 6; ++idf)
				for (int ifr = 0; ifr < hy.Nf; ++ifr)
					for (int ip = 0; ip < Np; ++ip) {
						const TComplex &c = pres(idf + 6*ib, ifr, ip);
						hy.pots[ib][ip][idf][ifr] = std::complex<double>(-c.Im, c.Re)/hy.rho/hy.w[ifr]; // press = i w rho phi(potential)
					}
	}
}
	
#endif