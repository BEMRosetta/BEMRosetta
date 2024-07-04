// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#ifdef PLATFORM_WIN32
#include "orca.h"
#endif

String OrcaWave::Load(String _file, double) {
	dt.file = _file;
	dt.name = GetFileTitle(dt.file);
	dt.dimen = true;
	dt.len = 1;
	dt.Nb = Null;
	
	try {
		BEM::Print("\n\n" + Format(t_("Loading '%s'"), dt.file));

		BEM::Print("\n- " + S(t_("YML file")));

#ifdef PLATFORM_WIN32
		if (GetFileExt(dt.file) == ".owr")
			Load_OWR();
		else
#endif
			Load_OF_YML();

		if (IsNull(dt.Nb))
			return t_("No data found");
	
		/*dof.Clear();	dof.SetCount(Nb, 0);
		for (int i = 0; i < Nb; ++i)
			dof[i] = 6;*/
	} catch (Exc e) {
		//BEM::PrintError(Format("\n%s: %s", t_("Error"), e));
		//lastError = e;
		return e;
	}
	return String();
}

#ifdef PLATFORM_WIN32
void OrcaWave::Load_OWR() {
	String fileName = ForceExtSafer(dt.file, ".owr");
	
	Orca orca;
	
	orca.LoadWaveResults(fileName);
	
	dt.solver = Hydro::ORCAWAVE_OWR;
	
	//dt.dataFromW = true;
	dt.dimen = true;
		
	dt.x_w = dt.y_w = 0;
	
	dt.Nb = dt.Nf = dt.Nh = Null;
	
	orca.LoadParameters(*this);
	
}
#endif

void OrcaWave::Load_OF_YML() {
	String fileName = ForceExtSafer(dt.file, ".yml");
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		throw Exc(in.Str() + "\n" + t_("File not found or blocked"));
	
	dt.solver = Hydro::ORCAFLEX_YML;
	
	//dt.dataFromW = true;
	char dataFrom = 'r';
	//bool rad_s = true;
	dt.dimen = true;
	
	dt.x_w = dt.y_w = 0;
	
	dt.Nb = dt.Nf = dt.Nh = Null;

	UArray<Point3D> c0s;	

	auto Origin = [&] (int ib, const UVector<String> &norig) {
		UVector<String> snorig(norig.size());
		for (int i = 0; i < norig.size(); ++i)
			snorig[i] = Trim(norig[i]);
		if (IsNull(c0s[ib])) 
			c0s[ib].Set(ScanDouble(snorig[0]), ScanDouble(snorig[1]), ScanDouble(snorig[2]));
		else if (snorig[0] == "~" || snorig[1] == "~" || snorig[2] == "~")
			;
		else if (c0s[ib].x != ScanDouble(snorig[0]) || c0s[ib].y != ScanDouble(snorig[1]) || c0s[ib].z != ScanDouble(snorig[2]))
			throw Exc(in.Str() + "\n"  + Format(t_("RAOOrigin for body %d (%s, %s, %s) diferent than the previously set (%f, %f, %f)"), 
							ib, snorig[0], snorig[1], snorig[2], c0s[ib].x, c0s[ib].y, c0s[ib].z));
	};
	auto Phase = [&] (const UVector<String> &norig) {
		UVector<String> snorig(norig.size());
		for (int i = 0; i < norig.size(); ++i)
			snorig[i] = Trim(norig[i]);
		if (snorig[0] == "~" || snorig[1] == "~" || snorig[2] == "~")
			;
		else if (!IsEqualRange<UVector<String>>({"0","0","0"}, snorig))
			throw Exc(in.Str() + "\n"  + Format(t_("Only PhaseOrigin:[0,0,0] is supported. Read (%s, %s, %s)"), snorig[0], snorig[1], snorig[2]));
	};
		
	OrcaFactors factor;	
	
	int ib = -1;
	dt.Nb = 0; 
	
	dt.g = 9.80665;		// Default value used when SI units
	
	YmlParser fy(in);

	FileInLine::Pos fpos = in.GetPos();

	while(fy.GetLine()) {
		if (fy.FirstIs("LoadRAOCalculationMethod")) 
			throw Exc(t_("This .yml is an OrcaWave case"));
		else if (fy.FirstIs("General")) {
			if (fy.FirstIs("UnitsSystem")) {
				if (fy.GetVal() == "SI") {
					factor.mass = factor.force = 1000;
					factor.len = 1;	
				} else if (fy.GetVal() == "User") 
					;
				else
					throw Exc(in.Str() + "\n" + Format(t_("Only SI and User units are supported. Read '%s'"), fy.GetVal()));
			} else if (fy.FirstIs("LengthUnits")) 
				factor.len = FactorLen(fy.GetVal());
			else if (fy.FirstIs("MassUnits")) 
				factor.mass = FactorMass(fy.GetVal());
			else if (fy.FirstIs("ForceUnits")) 
				factor.force = FactorForce(fy.GetVal());
			else if (fy.FirstIs("g")) 
				dt.g = ScanDouble(fy.GetVal())*factor.len;
		} else if (fy.FirstIs("Environment")) {
			if (fy.FirstIs("WaterSurfaceZ") && fy.GetVal() != "0") 
				throw Exc(in.Str() + "\n" + Format(t_("Only WaterSurfaceZ 0 is supported. Read '%s'"), fy.GetVal()));
		} else if (fy.FirstIs("VesselTypes")) {
			if (fy.FirstIs("Name")) {
				if (fy.Index() != dt.Nb)
					throw Exc(in.Str() + "\n" + t_("Failed body count"));
				dt.msh.SetCount(dt.Nb+1);
				dt.msh[dt.Nb].dt.name = fy.GetVal();
				dt.Nb++;
				c0s << Null;
			} else if (fy.FirstIs("WavesReferredToBy") && fy.Index() == 0) {		// Only for the first body
				String val = fy.GetVal();
				if (val.Find("(rad/s)") >= 0)
					//dt.dataFromW = true;
					dataFrom = 'r';
				else if (val.Find("(s)") >= 0)
					//dt.dataFromW = false;
					dataFrom = 's';
				else if (val.Find("(Hz)") >= 0)
					dataFrom = 'h';
				else
					throw Exc(in.Str() + "\n"  + Format(t_("Unknown data in WavesReferredToBy: %s"), val));
				
				/*if (val.Find("(rad/s)") >= 0)
					rad_s = true;
				else if (val.Find("(Hz)") >= 0)	 
					rad_s = false;
				else if (val.Find("(s)") >= 0)	 	
					;
				else
					throw Exc(in.Str() + "\n"  + Format(t_("Unknown data in WavesReferredToBy: %s"), val));*/
			} else if (fy.FirstIs("SurgePositive")) {	
				if (fy.GetVal() != "forward")
					throw Exc(in.Str() + "\n"  + Format(t_("Only SurgePositive: 'forward' is supported. Read '%s'"), fy.GetVal()));
			} else if (fy.FirstIs("SwayPositive")) {	
				if (fy.GetVal() != "port")
					throw Exc(in.Str() + "\n"  + Format(t_("Only SwayPositive: 'port' is supported. Read '%s'"), fy.GetVal()));
			} else if (fy.FirstIs("HeavePositive")) {	
				if (fy.GetVal() != "up")
					throw Exc(in.Str() + "\n"  + Format(t_("Only HeavePositive: 'up' is supported. Read '%s'"), fy.GetVal()));
			} else if (fy.FirstIs("RollPositiveStarboard")) {	
				if (fy.GetVal() != "down")
					throw Exc(in.Str() + "\n"  + Format(t_("Only RollPositiveStarboard: 'down' is supported. Read '%s'"), fy.GetVal()));
			} else if (fy.FirstIs("PitchPositiveBowe")) {	
				if (fy.GetVal() != "down")
					throw Exc(in.Str() + "\n"  + Format(t_("Only PitchPositiveBowe: 'down' is supported. Read '%s'"), fy.GetVal()));
			} else if (fy.FirstIs("YawPositiveBow")) {	
				if (fy.GetVal() != "port")
					throw Exc(in.Str() + "\n"  + Format(t_("Only YawPositiveBow: 'port' is supported. Read '%s'"), fy.GetVal()));
			} else if (fy.FirstIs("QTFConventionsRotationOrder")) {	
				if (fy.GetVal() != "RzRyRx")
					throw Exc(in.Str() + "\n"  + Format(t_("Only QTFConventionsRotationOrder: 'RzRyRx is supported. Read '%s'"), fy.GetVal()));
			/*} else if (fy.FirstIs("QTFConventionsRotationAxes")) {	
				if (fy.GetVal() != "original")
					throw Exc(in.Str() + "\n"  + Format(t_("Only QTFConventionsRotationAxes: 'original' is supported. Read '%s'"), fy.GetVal()));
			} else if (fy.FirstIs("QTFConventionsFrameOfReference")) {	
				if (fy.GetVal() != "earth")
					throw Exc(in.Str() + "\n"  + Format(t_("Only QTFConventionsFrameOfReference: 'earth' is supported. Read '%s'"), fy.GetVal()));
					*/
			} else if (fy.FirstIs("Draughts")) {
				ib = fy.GetIndex()[1];
				if (fy.FirstIs("DisplacementRAOs")) {
					if (fy.FirstIs("RAOOrigin")) 
						Origin(ib, fy.GetVector());
					else if (fy.FirstIs("PhaseOrigin")) 
						Phase(fy.GetVector());
					else if (fy.FirstIs("RAOs")) {
						if (fy.FirstIs("RAODirection") && fy.GetIndex()[1] == 0) {		// Only for the first body
							FindAdd(dt.head, ScanDouble(fy.GetVal()));
						}
					}
				} else if (fy.FirstIs("LoadRAOs")) {
					if (fy.FirstIs("RAOOrigin")) 
						Origin(ib, fy.GetVector());	
					else if (fy.FirstIs("PhaseOrigin")) 
						Phase(fy.GetVector());
					else if (fy.FirstIs("RAOs")) {
						if (fy.FirstIs("RAODirection") && fy.GetIndex()[1] == 0) {		// Only for the first body
							FindAdd(dt.head, ScanDouble(fy.GetVal()));
						}
					}
				} else if (fy.FirstIs("WaveDrift")) {
					if (fy.FirstIs("RAOOrigin")) {
						Origin(ib, fy.GetVector());	
					} else if (fy.FirstIs("RAOs")) {
						if (fy.FirstIs("RAODirection") && fy.GetIndex()[1] == 0) {		// Only for the first body
							FindAdd(dt.head, ScanDouble(fy.GetVal()));
						}
					}
				} else if (fy.FirstIs("SumFrequencyQTFs")) {
					if (fy.FirstIs("RAOOrigin")) 
						Origin(ib, fy.GetVector());	
					else if (fy.FirstIs("PhaseOrigin")) 
						Phase(fy.GetVector());
				} else if (fy.FirstIs("OtherDampingOrigin")) {
					if (!IsEqualRange<UVector<double>>({0,0,0}, fy.GetVectorDouble()))
						throw Exc(in.Str() + "\n"  + Format(t_("Only OtherDampingOrigin:[0,0,0] is supported. Read '%s'"), fy.StrVar()));
				} else if (fy.FirstIs("ReferenceOrigin")) {
					if (!IsEqualRange<UVector<double>>({0,0,0}, fy.GetVectorDouble()))
						throw Exc(in.Str() + "\n"  + Format(t_("Only ReferenceOrigin:[0,0,0] is supported. Read '%s'"), fy.StrVar()));
				} else if (fy.FirstIs("ReferenceOriginDatumPosition")) {
					if (!IsEqualRange<UVector<double>>({0,0,0}, fy.GetVectorDouble()))
						throw Exc(in.Str() + "\n"  + Format(t_("Only ReferenceOriginDatumPosition:[0,0,0] is supported. Read '%s'"), fy.StrVar()));
				} else if (fy.FirstIs("FrequencyDependentAddedMassAndDamping")) {
					if (fy.FirstIs("AMDPeriodOrFrequency")) {
						if (((dataFrom == 'r' || dataFrom == 'h') && fy.GetVal() != "Infinity") ||
							 (dataFrom == 's' && fy.GetVal() != "0")) {
							if (dt.w.size() != fy.Index() - 1)		// -1 because Infinity is the first
								throw Exc(in.Str() + "\n" + t_("Failed frequencies count"));			
							dt.w << ScanDouble(fy.GetVal());
						}
					}
				}
			}
		} else if (fy.FirstIs("MultibodyGroups")) {
			if (fy.FirstIs("Bodies")) {
				if (fy.FirstIs("Name")) 
					ib = fy.Index();
			} else if (fy.FirstIs("MultibodyAddedMassAndDamping")) {
				if (fy.FirstIs("AMDPeriodOrFrequency")) {
					if (fy.GetVal() != "Infinity") {
						if (dt.w.size() != fy.Index() - 1)
							throw Exc(in.Str() + "\n" + t_("Failed frequencies count"));			
						dt.w << ScanDouble(fy.GetVal());
					}
				}
			}
		}
	}
	
	factor.Update();
	
	if (dt.Nb == 0)
		throw Exc(S("\n") + t_("No body found"));
	
	dt.msh.SetCount(dt.Nb);
	//Vo.SetCount(Nb, NaNDouble);
	//cg.setConstant(3, Nb, NaNDouble);
	
	//c0.resize(3, Nb);
	for (int iib = 0; iib < dt.Nb; ++iib)				
		for (int idf = 0; idf < 3; ++idf)
			dt.msh[iib].dt.c0[idf] = c0s[iib][idf];
	
	//cb.setConstant(3, Nb, NaNDouble);
	//C.SetCount(Nb);
	for (int iib = 0; iib < dt.Nb; ++iib) 
		dt.msh[iib].dt.C.setConstant(6, 6, 0);
	//M.SetCount(Nb);
	for (int iib = 0; iib < dt.Nb; ++iib) 
		dt.msh[iib].dt.M.setConstant(6, 6, 0);

	dt.Nf = dt.w.size();
	dt.Nh = dt.head.size();
	
	dt.Ainf.setConstant(dt.Nb*6, dt.Nb*6, NaNDouble);
	
	Initialize_AB(dt.A);
	Initialize_AB(dt.B);
		
	if (dataFrom == 'h') {
		for (int i = 0; i < dt.w.size(); ++i)
				dt.w[i] *= 2*M_PI;	
	} else if (dataFrom == 's'){
		for (int i = 0; i < dt.w.size(); ++i)
			dt.w[i] = 2*M_PI/dt.w[i];	
	}	

	Initialize_Forces(dt.ex);
	Initialize_Forces(dt.rao);

	int rrow = -1, ccol = -1;
	int idf = -1;
	bool diffFullQTF = false;
	
	in.SeekPos(fpos);
	
	while(fy.GetLine()) {
		if (fy.FirstIs("Environment")) {
			if (fy.FirstIs("Density")) 
				dt.rho = ScanDouble(fy.GetVal())*factor.mass/factor.len/factor.len/factor.len;		// In kg/m3
			else if (fy.FirstIs("WaterDepth")) { 
				String sh = fy.GetVal(); 
				if (ToLower(sh) == "infinite")
					dt.h = -1;
				else	
					dt.h = ScanDouble(sh)*factor.len;
			}
		} else if (fy.FirstIs("VesselTypes")) {
			if (fy.FirstIs("Name")) 
				ib = fy.GetIndex()[1];
			else if (fy.FirstIs("Draughts")) {
				Eigen::MatrixXd &inertia = dt.msh[ib].dt.M;
				if (fy.FirstIs("Mass")) 
					inertia(0, 0) = inertia(1, 1) = inertia(2, 2) = ScanDouble(fy.GetVal())*factor.mass;
				else if (fy.FirstMatch("MomentOfInertiaTensor*")) {			// Referred to cg
					UVector<UVector<double>> mat = fy.GetMatrixDouble();
					
					inertia(3, 3) = mat[0][0]*factor.M(3, 3);		// In kg
					inertia(3, 4) = mat[0][1]*factor.M(3, 4);
					inertia(3, 5) = mat[0][2]*factor.M(3, 5);
					inertia(4, 3) = mat[1][0]*factor.M(4, 3);
					inertia(4, 4) = mat[1][1]*factor.M(4, 4);
					inertia(4, 5) = mat[1][2]*factor.M(4, 5);
					inertia(5, 3) = mat[2][0]*factor.M(5, 3);
					inertia(5, 4) = mat[2][1]*factor.M(5, 4);
					inertia(5, 5) = mat[2][2]*factor.M(5, 5);
				} else if (fy.FirstMatch("HydrostaticStiffnessz*")) {
					UVector<UVector<double>> mat = fy.GetMatrixDouble();
					
					for (int r = 0; r < 3; ++r)				// Only heave, roll, pitch
						for (int c = 0; c < 3; ++c)
							dt.msh[ib].dt.C(r+2, c+2) = mat[r][c]*factor.K(r+2, c+2);
				} else if (fy.FirstIs("CentreOfMass")) {
					UVector<double> line = fy.GetVectorDouble();
					
					dt.msh[ib].dt.cg.x = line[0]*factor.len;
					dt.msh[ib].dt.cg.y = line[1]*factor.len;
					dt.msh[ib].dt.cg.z = line[2]*factor.len;
				} else if (fy.FirstIs("CentreOfBuoyancy")) {
					UVector<double> line = fy.GetVectorDouble();
					
					dt.msh[ib].dt.cb.x = line[0]*factor.len;
					dt.msh[ib].dt.cb.y = line[1]*factor.len;
					dt.msh[ib].dt.cb.z = line[2]*factor.len;
				} else if (fy.FirstIs("DisplacedVolume")) 
					dt.msh[ib].dt.Vo = ScanDouble(fy.GetVal());
				else if (fy.FirstIs("DisplacementRAOs")) {
					if (fy.FirstIs("RAOs")) {
						if (fy.FirstMatch("RAOPeriodOrFrequency*")) {	
							int idh = fy.Index();
							if (idh < 0 || idh >= dt.head.size())
								throw Exc(in.Str() + "\n" + t_("Wrong heading"));
								
							UVector<UVector<double>> mat = fy.GetMatrixDouble();
							
							if (mat.size() != dt.Nf || mat[0].size() != 13)
								throw Exc(in.Str() + "\n"  + Format(t_("Wrong number of numbers in DisplacementRAOs matrix"), fy.GetText()));

							for (int ifr = 0; ifr < dt.Nf; ++ifr) 
								for (int idof = 0; idof < 6; ++idof) 
									dt.rao[ib][idh](ifr, idof) = std::polar<double>(mat[ifr][1 + 2*idof]*factor.RAO(idof), ToRad(mat[ifr][1 + 2*idof + 1]));
						}
					}
				}  else if (fy.FirstIs("LoadRAOs")) {
					if (fy.FirstIs("RAOs")) {
						if (fy.FirstMatch("RAOPeriodOrFrequency*")) {	
							int idh = fy.Index();
							if (idh < 0 || idh >= dt.head.size())
								throw Exc(in.Str() + "\n" + t_("Wrong heading"));
								
							UVector<UVector<double>> mat = fy.GetMatrixDouble();
							
							if (mat.size() != dt.Nf || mat[0].size() != 13)
								throw Exc(in.Str() + "\n"  + Format(t_("Wrong number of numbers in LoadRAOs matrix"), fy.GetText()));

							for (int ifr = 0; ifr < dt.Nf; ++ifr) 
								for (int idof = 0; idof < 6; ++idof) 
									dt.ex[ib][idh](ifr, idof) = std::polar<double>(mat[ifr][1 + 2*idof]*factor.F(idof), ToRad(mat[ifr][1 + 2*idof + 1]));
						}
					}
				} else if (fy.FirstIs("WaveDriftQTFMethod"))
					diffFullQTF = fy.GetVal() == "Full QTFs";
				else if (!diffFullQTF && fy.FirstIs("WaveDrift")) {
					if (fy.FirstIs("RAOs")) {
						if (fy.FirstMatch("RAOPeriodOrFrequency*")) {	
							if (!IsLoadedMD()) {
								dt.mdhead.resize(dt.head.size());
								for (int ih = 0; ih < dt.head.size(); ++ih)
									dt.mdhead[ih] = std::complex<double>(dt.head[ih], dt.head[ih]);
								Hydro::Initialize_MD(dt.md, dt.Nb, int(dt.mdhead.size()), dt.Nf);
							}
							
							int idh = fy.Index();
							if (idh < 0 || idh >= dt.head.size())
								throw Exc(in.Str() + "\n" + t_("Wrong heading"));
								
							UVector<UVector<double>> mat = fy.GetMatrixDouble();
							
							if (mat.size() != dt.Nf || mat[0].size() != 7)
								throw Exc(in.Str() + "\n"  + Format(t_("Wrong number of numbers in WaveDrift matrix"), fy.GetText()));

							for (int ifr = 0; ifr < dt.Nf; ++ifr) 
								for (int idof = 0; idof < 6; ++idof) 
									dt.md[ib][idh][idof](ifr) = mat[ifr][1 + idof]*factor.MD(idof);
						}
					}
				} else if (fy.FirstIs("SumFrequencyQTFs") || (diffFullQTF && fy.FirstIs("WaveDrift"))) {
					if (fy.FirstMatch("RAOPeriodOrFrequency*")) {
						UVector<UVector<double>> mat = fy.GetMatrixDouble();
						
						UArray<UArray<UArray<MatrixXcd>>> &q = diffFullQTF ? dt.qtfdif : dt.qtfsum;
						double phmult = !diffFullQTF ? 1 : -1;		// Difference is conjugate-symmetric
						
						if (!IsLoadedQTF(!diffFullQTF)) {		// Gets frequencies and headings
							::Copy(dt.w, dt.qw);
							
							UArray<std::complex<double>> qh;
							for (int row = 0; row < mat.size(); ++row) {
								if (mat[row].size() != 16)
									throw Exc(in.Str() + "\n"  + t_("Wrong data in list"));
								double h1 = mat[row][2],
									   h2 = mat[row][3];
									   
								FindAddDelta(qh, std::complex<double>(h1, h2), 0.0001);	
							}
							::Copy(qh, dt.qhead);
	
	
							Hydro::Initialize_QTF(q, dt.Nb, int(qh.size()), dt.Nf);
							dt.mdtype = dt.qtftype = 9;
						}
						diffFullQTF = false;
						
						for (int row = 0; row < mat.size(); ++row) {
							if (mat[row].size() != 16)
								throw Exc(in.Str() + "\n"  + t_("Wrong data in list"));
							double w1 = mat[row][0],
								   w2 = mat[row][1],
								   h1 = mat[row][2],
								   h2 = mat[row][3];
								   
							if (dataFrom == 'h') {
								w1 *= 2*M_PI;
								w2 *= 2*M_PI;
							} else if (dataFrom == 's') {
								w1 = 2*M_PI/w1;
								w2 = 2*M_PI/w2;
							}
							int ifr1 = FindDelta(dt.qw, w1, 0.0001),
								ifr2 = FindDelta(dt.qw, w2, 0.0001),
								ih = FindDelta(dt.qhead, std::complex<double>(h1, h2), 0.0001);	
							if (ifr1 < 0)
								throw Exc(in.Str() + "\n"  + Format(t_("Wrong frequency '%d' in QTF"), w1));
							if (ifr2 < 0)
								throw Exc(in.Str() + "\n"  + Format(t_("Wrong frequency '%d' in QTF"), w2));
							if (ih < 0)
								throw Exc(in.Str() + "\n"  + Format(t_("Wrong head (%d,%d) in QTF"), h1, h2));
							for (int iidf = 0; iidf < 6; ++iidf) {
								double mag = mat[row][4+iidf*2]*factor.F(iidf);
								double ph  = ToRad(mat[row][4+iidf*2+1]);
								q[ib][ih][iidf](ifr1, ifr2) = std::polar<double>(mag, ph);
								q[ib][ih][iidf](ifr2, ifr1) = std::polar<double>(mag, ph*phmult);	// Diagonal
							}
						}
					}
				} else if (fy.FirstIs("FrequencyDependentAddedMassAndDamping")) {
					if (fy.FirstMatch("AddedMassMatrixX*")) {
						idf = fy.Index()-1;
						if (idf < -1 || idf >= dt.w.size())		// Infinity is the first
							throw Exc(in.Str() + "\n" + t_("Wrong frequency"));			
						
						UVector<UVector<double>> mat = fy.GetMatrixDouble();
						
						if (idf == -1) {
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) 
									dt.Ainf(r, c) = mat[r][c]*factor.A(r, c); 
							}
						} else {
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) 
									dt.A[r][c](idf) = mat[r][c]*factor.A(r, c); 
							}
						}
					} else if (fy.FirstMatch("DampingX*")) {
						idf = fy.Index()-1;
						if (idf < -1 || idf >= dt.w.size())		// Infinity is the first
							throw Exc(in.Str() + "\n" + t_("Wrong frequency"));			
						//idf--;
						
						if (idf >= 0) {
							UVector<UVector<double>> mat = fy.GetMatrixDouble();
							
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) 
									dt.B[r][c](idf) = mat[r][c]*factor.B(r, c); 
							}
						}
					}
				}
			}
		} else if (fy.FirstIs("MultibodyGroups")) {
			if (fy.FirstIs("Bodies")) {
				if (fy.FirstIs("Name")) {
					ib = fy.Index();
					dt.msh[ib].dt.name = fy.GetVal();
				} else if (fy.FirstIs("DisplacedVolume")) 
					dt.msh[ib].dt.Vo = ScanDouble(fy.GetVal())*factor.len*factor.len*factor.len;
				else if (fy.FirstIs("CentreOfBuoyancy")) {
					UVector<double> line = fy.GetVectorDouble();
					
					dt.msh[ib].dt.cb.x = line[0]*factor.len;
					dt.msh[ib].dt.cb.y = line[1]*factor.len;
					dt.msh[ib].dt.cb.z = line[2]*factor.len;
				} else if (fy.FirstMatch("HydrostaticStiffnessz*")) {
					UVector<UVector<double>> mat = fy.GetMatrixDouble();
					
					for (int r = 0; r < 3; ++r)				// Only heave, roll, pitch
						for (int c = 0; c < 3; ++c)
							dt.msh[ib].dt.C(r+2, c+2) = mat[r][c]*factor.K(r, c);
				}
			} else if (fy.FirstIs("MultibodyAddedMassAndDamping")) {
				if (fy.FirstIs("AMDPeriodOrFrequency")) {
					idf = fy.Index()-1;
					if (idf < -1 || idf >= dt.w.size())		// Infinity is the first
						throw Exc(in.Str() + "\n" + t_("Wrong frequency"));	
						/*
					if (fy.GetVal() == "Infinity")
						idf = -1;
					else {
						idf = Find(w, ScanDouble(fy.GetVal()));
						if (idf < 0)
							throw Exc(in.Str() + "\n" + t_("Wrong frequency"));	
					}*/
				} else if (fy.FirstIs("Matrices")) {
					if (fy.FirstIs("Row")) 
						rrow = ScanInt(fy.GetVal())-1;
					else if (fy.FirstIs("Column")) 
						ccol = ScanInt(fy.GetVal())-1;
					else if (fy.FirstMatch("AddedMassX*")) {
						UVector<UVector<double>> mat = fy.GetMatrixDouble();
						
						if (idf == -1) {
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) {
									dt.Ainf(r+rrow*6, c+ccol*6) = mat[r][c]*factor.A(r, c); 
									if (rrow != ccol)		// Fill the symmetric
										dt.Ainf(r+ccol*6, c+rrow*6) = dt.Ainf(r+rrow*6, c+ccol*6); 
								}
							}
						} else {
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) {
									dt.A[r+rrow*6][c+ccol*6](idf) = mat[r][c]*factor.A(r, c); 
									if (rrow != ccol)		// Fill the symmetric
										dt.A[r+ccol*6][c+rrow*6](idf) = dt.A[r+rrow*6][c+ccol*6](idf); 
								}
							}
						}
					} else if (fy.FirstMatch("DampingX*")) {
						if (idf >= 0) {
							UVector<UVector<double>> mat = fy.GetMatrixDouble();
							
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) {
									dt.B[r+rrow*6][c+ccol*6](idf) = mat[r][c]*factor.B(r, c); 
									if (rrow != ccol)		// Fill the symmetric
										dt.B[r+ccol*6][c+rrow*6](idf) = dt.B[r+rrow*6][c+ccol*6](idf);
								}
							}
						}
					}
				}
			}
		}
	}
	
	if (dt.Nb == 0)
		throw Exc(t_("Incorrect .yml format"));
	
	// Inertia matrices have to be translated from cg to c0
	for (int iib = 0; iib < dt.Nb; ++iib) //{
		//Point3D cg(cg.col(ib));
		//Point3D c0(c0.col(ib));
		Surface::TranslateInertia66(dt.msh[iib].dt.M, dt.msh[iib].dt.cg, dt.msh[iib].dt.cg, dt.msh[iib].dt.c0);
	//}
}

void OrcaWave::SaveFolder_OW_YML(String folder, bool bin, int numThreads, bool withPotentials, bool withMesh, bool withQTF, bool x0z, bool y0z) const {
	String exeName = "bemrosetta_cl";
	if (bin) {
		String exe = GetExeFilePath();		
		String ext = GetFileExt(exe);
		String exe_cl = AFX(GetFileFolder(GetExeFilePath()), exeName + ext);
		if (FileExists(exe_cl))
			FileCopy(exe_cl, AFX(folder, exeName + ext));
		else {
			exeName = "bemrosetta";
			FileCopy(exe, AFX(folder, exeName + ext));
		}
	}
	
	String name = GetFileTitle(folder); 
	String fileYaml = AFX(folder, name + ".wave.yml");
	String fileBat  = AFX(folder, name + ".bat");
	FileOut bat;
	if (!bat.Open(fileBat))
		throw Exc(Format(t_("Impossible to open file '%s'"), fileBat));
	bat << Format("%s -orca -numtries 5 -numthread %d -rw \"%s\" \"%s\"", exeName, numThreads, fileYaml, AFX(folder, name + ".flex.yml"));
	
	FileOut	out;
	if (!out.Open(fileYaml))
		throw Exc(Format(t_("Impossible to open file '%s'"), fileYaml));
	
	for (int ib = 0; ib < dt.Nb; ++ib) {
		const Body &b = dt.msh[ib];
		Body::SaveAs(b, AFX(folder, Format("Body_%d.gdf", ib+1)), Body::WAMIT_GDF, Body::UNDERWATER, dt.rho, dt.g, y0z, x0z);
	}
	
	String machineName, domain, ip4, ip6;
	GetNetworkInfo(machineName, domain, ip4, ip6);
	Time t = GetSysTime(); 
	
	out << 	"\%YAML 1.1\n"
			"# Type: Diffraction\n"
			"# Program: OrcaWave 11.4c\n"
			"# File: " << fileYaml << "\n"
			"# Created: " << Format("%02d:%02d", t.hour, t.minute) << " on " << Format("%02d/%02d/%04d", t.day, t.month, t.year) << "\n"
			"# User: " << GetUserName() << "\n"
			"# Machine: " << machineName << "\n"
			"---\n";
			
	out << 	"# Model\n"
			"UnitsSystem: User\n"
			"LengthUnits: m\n"
			"MassUnits: kg\n"
			"ForceUnits: N\n"
			"g: " << Format("%.5f", dt.g) << "\n";
	
	out << 	"# Calculation & output\n"
			"SolveType: " << (withQTF ? "Full QTF calculation" : "Potential and source formulations") << "\n"
			"LoadRAOCalculationMethod: Diffraction\n"
			"QuadraticLoadPressureIntegration: No\n"
			"QuadraticLoadControlSurface: " << (withQTF ? "Yes" : "No") << "\n"
			"QuadraticLoadMomentumConservation: No\n";
	out << 	(withQTF ? "PreferredQuadraticLoadCalculationMethod: Control surface\n" : "");
	out << 	"HasResonanceDampingLid: No\n"
			"LengthTolerance: 100e-9\n"
			"WaterlineZTolerance: 1e-6\n"
			"WaterlineGapTolerance: 1e-6\n"
			"DivideNonPlanarPanels: Yes\n"
			"LinearSolverMethod: Direct LU\n"
			"OutputPanelPressures: " << (withPotentials ? "Yes" : "No") << "\n"
			"OutputPanelVelocities: No\n"
			"OutputBodyWireFrames: " << (withMesh ? "Yes" : "No") << "\n"
			"OutputIntermediateResults: " << (withPotentials ? "Yes" : "No") << "\n"
			"ValidatePanelArrangement: No\n"
			"BodyVolumeWarningLevel: 1e-12\n"
			"PanelAspectRatioWarningLevel: 25\n"
			"PanelsPerWavelengthWarningLevel: 5\n";
	
	out << 	"# Environment\n"
			"WaterDepth: " << (dt.h > 0 ? Format("%.2f", dt.h) : "Infinity") << "\n"
			"WaterDensity: " << Format("%.2f", dt.rho) << "\n"
			"WavesReferredToBy: frequency (rad/s)\n"
			"HasWaveSpectrumForDragLinearisation: No\n"
			"MorisonFluidVelocity: Undisturbed incident wave\n"
			"PeriodOrFrequency:\n";
	for (int ifr = 0; ifr < dt.Nf; ++ifr)
		out << "  - " << Format("%.5f", dt.w[ifr]) << "\n";
	out << 	"WaveHeading:\n";
	for (int ih = 0; ih < dt.Nh; ++ih)
		out << "  - " << Format("%.1f", dt.head[ih]) << "\n";
	
	if (withQTF) 
		out << 	"QTFMinCrossingAngle: 0\n"
				"QTFMaxCrossingAngle: 180\n"
				"QTFMinPeriodOrFrequency: 0\n"
				"QTFMaxPeriodOrFrequency: Infinity\n"
				"QTFFrequencyTypes: Both\n";
	
	out << 	"# Bodies\n"
			"Bodies:\n";
	for (int ib = 0; ib < dt.Nb; ++ib) {
		const Body &b = dt.msh[ib];
		const Body::Data &d = b.dt;
		
		double panelSize = 0;	///////////////////////////////////////////////////////////////////////
		double separationFromBody = 0;		//////////////////////////////////////////////////////
		
		out <<	"  - BodyName: " << d.name << "\n"
				"    BodyMeshPosition: [0, 0, 0]\n"
				"    BodyMeshAttitude: [0, 0, 0]\n"
				"    BodyIncludedInAnalysis: Yes\n"
				"    BodyMeshFileName: " << Format(".\\Body_%d.gdf", ib+1) << "\n"
				"    BodyMeshFormat: Wamit gdf\n"
				"    BodyMeshLengthUnits: m\n"
				"    BodyMeshSymmetry: ";
		if (y0z) {	
			if (x0z)
				out << 	"xz and yz planes";
			else
				out << 	"yz plane";
		} else if (x0z)
			out << 	"xz plane";
		else
			out << 	"None";
		out << 	"\n"
				"    BodyMeshDipolePanels:\n"
				"    BodyAddInteriorSurfacePanels: Yes\n"
				"    BodyInteriorSurfacePanelMethod: Triangulation method\n";
				
		if (withQTF)
			out <<	"    BodyControlSurfaceType: Automatically generated\n"
					"    BodyControlSurfacePanelSize: " << panelSize << "\n"
					"    BodyControlSurfaceSeparationFromBody: " << separationFromBody << "\n"
					"    BodyControlSurfaceIncludeFreeSurface: Yes\n";		
    	out <<	"    BodyOrcaFlexImportSymmetry: Use global mesh symmetry\n"
    			"    BodyOrcaFlexImportLength: ~\n"
    			"    BodyHydrostaticIntegralMethod: Standard\n"
    			"    BodyHydrostaticStiffnessMethod: Displacement\n"
    			"    BodyInertiaSpecifiedBy: Matrix (for a general body)\n"
    			"    BodyCentreOfMass: [" << Format("%f, %f, %f", d.cg.x, d.cg.y, d.cg.z) << "]\n"
    			"    BodyMass: " << b.GetMass() << "\n"
    			"    BodyInertiaTensorRx, BodyInertiaTensorRy, BodyInertiaTensorRz:\n";
		if (d.M.size() > 0)
			out << 	"      - [" << Format("%f, %f, %f", d.M(3, 3), d.M(3, 4), d.M(3, 5)) << "]\n"
					"      - [" << Format("%f, %f, %f", d.M(4, 3), d.M(4, 4), d.M(4, 5)) << "]\n"
					"      - [" << Format("%f, %f, %f", d.M(5, 3), d.M(5, 4), d.M(5, 5)) << "]\n";
		else
			out <<	"      - [0.001, 0, 	0]\n"
					"      - [0, 	 0.001,	0]\n"
					"      - [0, 	 0, 	0.001]\n";
    	out << 	"    BodyInertiaTensorOriginType: Body origin\n"
    			"    BodyExternalStiffnessMatrixx, BodyExternalStiffnessMatrixy, BodyExternalStiffnessMatrixz, BodyExternalStiffnessMatrixRx, BodyExternalStiffnessMatrixRy, BodyExternalStiffnessMatrixRz:\n";
    	if (d.Cadd.size() > 0) {
    		for (int r = 0; r < 6; ++r) {
    			out << "      - [";
    			for (int c = 0; c < 6; ++c) {
    				if (c > 0)
    					out << ", ";
    				out << d.Cadd(r, c);
    			}
    			out << "]\n";
    		}
		} else {
			for (int r = 0; r < 6; ++r) 
				out << "      - [0, 0, 0, 0, 0, 0]\n";
		}
    	out <<	"    BodyExternalStiffnessMatrixOriginType: Body origin\n"
    			"    BodyExternalDampingMatrixx, BodyExternalDampingMatrixy, BodyExternalDampingMatrixz, BodyExternalDampingMatrixRx, BodyExternalDampingMatrixRy, BodyExternalDampingMatrixRz:\n";
    	if (d.Dlin.size() > 0) {
    		for (int r = 0; r < 6; ++r) {
    			out << "      - [";
    			for (int c = 0; c < 6; ++c) {
    				if (c > 0)
    					out << ", ";
    				out << d.Dlin(r, c);
    			}
    			out << "]\n";
    		}
		} else {
			for (int r = 0; r < 6; ++r) 
				out << "      - [0, 0, 0, 0, 0, 0]\n";
		}		
    	out <<	"    BodyExternalDampingMatrixOriginType: Body origin\n"
    			"    BodyConnectionParent: Free\n"
    			"    BodyIncreaseRollDampingToTarget: No\n"
    			"    BodyFixedDOFx: No\n"
    			"    BodyFixedDOFy: No\n"
    			"    BodyFixedDOFz: No\n"
    			"    BodyFixedDOFRx: No\n"
    			"    BodyFixedDOFRy: No\n"
    			"    BodyFixedDOFRz: No\n";
	}
	out <<	"# Field points\n"
			"DetectAndSkipFieldPointsInsideBodies: Yes\n";
	if (withQTF) 
		out << 	"# QTFs\n"
				"QTFCalculationMethod: Both\n"
				"PreferredQTFCalculationMethod: Direct method\n"
				"FreeSurfacePanelledZoneType: Defined by mesh file\n"
				"FreeSurfacePanelledZoneMeshFileName: n"
				"FreeSurfacePanelledZoneMeshFormat: Wamit gdf\n"
				"FreeSurfacePanelledZoneMeshLengthUnits: m\n";
	out << 	"...\n";
}