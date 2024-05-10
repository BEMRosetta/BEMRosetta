// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#ifdef PLATFORM_WIN32
#include "orca.h"
#endif

String OrcaWave::Load(String file, double) {
	hd().file = file;
	hd().name = GetFileTitle(file);
	hd().dimen = true;
	hd().len = 1;
	hd().solver = Hydro::ORCAWAVE;
	hd().Nb = Null;
	
	try {
		BEM::Print("\n\n" + Format(t_("Loading '%s'"), file));

		BEM::Print("\n- " + S(t_("YML file")));

#ifdef PLATFORM_WIN32
		if (GetFileExt(file) == ".owr")
			Load_OWR();
		else
#endif
			Load_YML_Res();

		if (IsNull(hd().Nb))
			return t_("No data found");
	
		/*hd().dof.Clear();	hd().dof.SetCount(hd().Nb, 0);
		for (int i = 0; i < hd().Nb; ++i)
			hd().dof[i] = 6;*/
	} catch (Exc e) {
		//BEM::PrintError(Format("\n%s: %s", t_("Error"), e));
		//hd().lastError = e;
		return e;
	}
	return String();
}

#ifdef PLATFORM_WIN32
void OrcaWave::Load_OWR() {
	String fileName = ForceExtSafer(hd().file, ".owr");
	
	Orca orca;
	
	orca.LoadWaveResults(fileName);
	
	hd().dataFromW = true;
	hd().dimen = true;
		
	hd().x_w = hd().y_w = 0;
	
	hd().Nb = hd().Nf = hd().Nh = Null;
	
	orca.LoadParameters(hd());
	
}
#endif

void OrcaWave::Load_YML_Res() {
	String fileName = ForceExtSafer(hd().file, ".yml");
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		throw Exc(in.Str() + "\n" + t_("File not found or blocked"));
	
	hd().dataFromW = true;
	bool rad_s = true;
	hd().dimen = true;
	
	hd().x_w = hd().y_w = 0;
	
	hd().Nb = hd().Nf = hd().Nh = Null;

	UArray<Point3D> c0s;	

	auto Origin =[&] (int ib, const UVector<String> &norig) {
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
	auto Phase =[&] (const UVector<String> &norig) {
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
	int Nb = 0; 
	
	hd().g = 9.80665;		// Default value used when SI units
	
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
				hd().g = ScanDouble(fy.GetVal())*factor.len;
		} else if (fy.FirstIs("Environment")) {
			if (fy.FirstIs("WaterSurfaceZ") && fy.GetVal() != "0") 
				throw Exc(in.Str() + "\n" + Format(t_("Only WaterSurfaceZ 0 is supported. Read '%s'"), fy.GetVal()));
		} else if (fy.FirstIs("VesselTypes")) {
			if (fy.FirstIs("Name")) {
				if (fy.Index() != Nb)
					throw Exc(in.Str() + "\n" + t_("Failed body count"));
				hd().msh.SetCount(Nb+1);
				hd().msh[Nb].name = fy.GetVal();
				Nb++;
				c0s << Null;
			} else if (fy.FirstIs("WavesReferredToBy") && fy.Index() == 0) {		// Only for the first body
				String val = fy.GetVal();
				if (val.Find("frequency") >= 0)
					hd().dataFromW = true;
				else if (val.Find("period") >= 0)
					hd().dataFromW = false;
				else
					throw Exc(in.Str() + "\n"  + Format(t_("Unknown data in WavesReferredToBy: %s"), val));
				
				if (val.Find("(rad/s)") >= 0)
					rad_s = true;
				else if (val.Find("(Hz)") >= 0)	 
					rad_s = false;
				else if (val.Find("(s)") >= 0)	 	
					;
				else
					throw Exc(in.Str() + "\n"  + Format(t_("Unknown data in WavesReferredToBy: %s"), val));
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
							FindAdd(hd().head, ScanDouble(fy.GetVal()));
						}
					}
				} else if (fy.FirstIs("LoadRAOs")) {
					if (fy.FirstIs("RAOOrigin")) 
						Origin(ib, fy.GetVector());	
					else if (fy.FirstIs("PhaseOrigin")) 
						Phase(fy.GetVector());
					else if (fy.FirstIs("RAOs")) {
						if (fy.FirstIs("RAODirection") && fy.GetIndex()[1] == 0) {		// Only for the first body
							FindAdd(hd().head, ScanDouble(fy.GetVal()));
						}
					}
				} else if (fy.FirstIs("WaveDrift")) {
					if (fy.FirstIs("RAOOrigin")) {
						Origin(ib, fy.GetVector());	
					} else if (fy.FirstIs("RAOs")) {
						if (fy.FirstIs("RAODirection") && fy.GetIndex()[1] == 0) {		// Only for the first body
							FindAdd(hd().head, ScanDouble(fy.GetVal()));
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
						if (fy.GetVal() != "Infinity") {
							if (hd().w.size() != fy.Index() - 1)		// -1 because Infinity is the first
								throw Exc(in.Str() + "\n" + t_("Failed frequencies count"));			
							hd().w << ScanDouble(fy.GetVal());
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
						if (hd().w.size() != fy.Index() - 1)
							throw Exc(in.Str() + "\n" + t_("Failed frequencies count"));			
						hd().w << ScanDouble(fy.GetVal());
					}
				}
			}
		}
	}
	
	factor.Update();
	
	if (Nb == 0)
		throw Exc(S("\n") + t_("No body found"));
	
	hd().Nb = Nb;
	hd().msh.SetCount(hd().Nb);
	//hd().Vo.SetCount(hd().Nb, NaNDouble);
	//hd().cg.setConstant(3, hd().Nb, NaNDouble);
	
	//hd().c0.resize(3, hd().Nb);
	for (int ib = 0; ib < Nb; ++ib)				
		for (int idf = 0; idf < 3; ++idf)
			hd().msh[ib].c0[idf] = c0s[ib][idf];
	
	//hd().cb.setConstant(3, hd().Nb, NaNDouble);
	//hd().C.SetCount(hd().Nb);
	for (int ib = 0; ib < hd().Nb; ++ib) 
		hd().msh[ib].C.setConstant(6, 6, 0);
	//hd().M.SetCount(hd().Nb);
	for (int ib = 0; ib < hd().Nb; ++ib) 
		hd().msh[ib].M.setConstant(6, 6, 0);

	hd().Nf = hd().w.size();
	hd().Nh = hd().head.size();
	
	hd().Ainf.setConstant(hd().Nb*6, hd().Nb*6, NaNDouble);
	
	hd().Initialize_AB(hd().A);
	hd().Initialize_AB(hd().B);
		
	if (hd().dataFromW) {
		if (!rad_s) {
			for (int i = 0; i < hd().w.size(); ++i)
				hd().w[i] *= 2*M_PI;	
		}
		/*hd().T.SetCount(hd().w.size());
		for (int i = 0; i < hd().w.size(); ++i)
			hd().T[i] = 2*M_PI/hd().w[i];	*/
	} else {
		//hd().T = pick(hd().w);
		//hd().w.SetCount(hd().T.size());
		for (int i = 0; i < hd().w.size(); ++i)
			hd().w[i] = 2*M_PI/hd().w[i];	
	}	

	hd().Initialize_Forces(hd().ex);
	hd().Initialize_Forces(hd().rao);

	int row = -1, col = -1;
	int idf = -1;
	bool diffFullQTF = false;
	
	in.SeekPos(fpos);
	

	while(fy.GetLine()) {
		if (fy.FirstIs("Environment")) {
			if (fy.FirstIs("Density")) 
				hd().rho = ScanDouble(fy.GetVal())*factor.mass/factor.len/factor.len/factor.len;		// In kg/m3
			else if (fy.FirstIs("WaterDepth")) { 
				String h = fy.GetVal(); 
				if (ToLower(h) == "infinite")
					hd().h = -1;
				else	
					hd().h = ScanDouble(h)*factor.len;
			}
		} else if (fy.FirstIs("VesselTypes")) {
			if (fy.FirstIs("Name")) 
				ib = fy.GetIndex()[1];
			else if (fy.FirstIs("Draughts")) {
				Eigen::MatrixXd &inertia = hd().msh[ib].M;
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
							hd().msh[ib].C(r+2, c+2) = mat[r][c]*factor.K(r+2, c+2);
				} else if (fy.FirstIs("CentreOfMass")) {
					UVector<double> line = fy.GetVectorDouble();
					
					hd().msh[ib].cg.x = line[0]*factor.len;
					hd().msh[ib].cg.y = line[1]*factor.len;
					hd().msh[ib].cg.z = line[2]*factor.len;
				} else if (fy.FirstIs("CentreOfBuoyancy")) {
					UVector<double> line = fy.GetVectorDouble();
					
					hd().msh[ib].cb.x = line[0]*factor.len;
					hd().msh[ib].cb.y = line[1]*factor.len;
					hd().msh[ib].cb.z = line[2]*factor.len;
				} else if (fy.FirstIs("DisplacedVolume")) 
					hd().msh[ib].Vo = ScanDouble(fy.GetVal());
				else if (fy.FirstIs("DisplacementRAOs")) {
					if (fy.FirstIs("RAOs")) {
						if (fy.FirstMatch("RAOPeriodOrFrequency*")) {	
							int idh = fy.Index();
							if (idh < 0 || idh >= hd().head.size())
								throw Exc(in.Str() + "\n" + t_("Wrong heading"));
								
							UVector<UVector<double>> mat = fy.GetMatrixDouble();
							
							if (mat.size() != hd().Nf || mat[0].size() != 13)
								throw Exc(in.Str() + "\n"  + Format(t_("Wrong number of numbers in DisplacementRAOs matrix"), fy.GetText()));

							for (int ifr = 0; ifr < hd().Nf; ++ifr) 
								for (int idof = 0; idof < 6; ++idof) 
									hd().rao.force[idh](ifr, idof+6*ib) = std::polar<double>(mat[ifr][1 + 2*idof]*factor.RAO(idof), ToRad(mat[ifr][1 + 2*idof + 1]));
						}
					}
				}  else if (fy.FirstIs("LoadRAOs")) {
					if (fy.FirstIs("RAOs")) {
						if (fy.FirstMatch("RAOPeriodOrFrequency*")) {	
							int idh = fy.Index();
							if (idh < 0 || idh >= hd().head.size())
								throw Exc(in.Str() + "\n" + t_("Wrong heading"));
								
							UVector<UVector<double>> mat = fy.GetMatrixDouble();
							
							if (mat.size() != hd().Nf || mat[0].size() != 13)
								throw Exc(in.Str() + "\n"  + Format(t_("Wrong number of numbers in LoadRAOs matrix"), fy.GetText()));

							for (int ifr = 0; ifr < hd().Nf; ++ifr) 
								for (int idof = 0; idof < 6; ++idof) 
									hd().ex.force[idh](ifr, idof+6*ib) = std::polar<double>(mat[ifr][1 + 2*idof]*factor.F(idof), ToRad(mat[ifr][1 + 2*idof + 1]));
						}
					}
				} else if (fy.FirstIs("WaveDriftQTFMethod"))
					diffFullQTF = fy.GetVal() == "Full QTFs";
				else if (!diffFullQTF && fy.FirstIs("WaveDrift")) {
					if (fy.FirstIs("RAOs")) {
						if (fy.FirstMatch("RAOPeriodOrFrequency*")) {	
							if (!hd().IsLoadedMD()) {
								hd().mdhead.resize(hd().head.size());
								for (int ih = 0; ih < hd().head.size(); ++ih)
									hd().mdhead[ih] = std::complex<double>(hd().head[ih], hd().head[ih]);
								Hydro::Initialize_MD(hd().md, hd().Nb, int(hd().mdhead.size()), hd().Nf);
							}
							
							int idh = fy.Index();
							if (idh < 0 || idh >= hd().head.size())
								throw Exc(in.Str() + "\n" + t_("Wrong heading"));
								
							UVector<UVector<double>> mat = fy.GetMatrixDouble();
							
							if (mat.size() != hd().Nf || mat[0].size() != 7)
								throw Exc(in.Str() + "\n"  + Format(t_("Wrong number of numbers in WaveDrift matrix"), fy.GetText()));

							for (int ifr = 0; ifr < hd().Nf; ++ifr) 
								for (int idof = 0; idof < 6; ++idof) 
									hd().md[ib][idh][idof](ifr) = mat[ifr][1 + idof]*factor.MD(idof);
						}
					}
				} else if (fy.FirstIs("SumFrequencyQTFs") || (diffFullQTF && fy.FirstIs("WaveDrift"))) {
					if (fy.FirstMatch("RAOPeriodOrFrequency*")) {
						UVector<UVector<double>> mat = fy.GetMatrixDouble();
						
						UArray<UArray<UArray<MatrixXcd>>> &q = diffFullQTF ? hd().qtfdif : hd().qtfsum;
						double phmult = !diffFullQTF ? 1 : -1;		// Difference is conjugate-symmetric
						
						if (!hd().IsLoadedQTF(!diffFullQTF)) {		// Gets frequencies and headings
							::Copy(hd().w, hd().qw);
							
							UArray<std::complex<double>> qh;
							for (int row = 0; row < mat.size(); ++row) {
								if (mat[row].size() != 16)
									throw Exc(in.Str() + "\n"  + t_("Wrong data in list"));
								double h1 = mat[row][2],
									   h2 = mat[row][3];
									   
								FindAddDelta(qh, std::complex<double>(h1, h2), 0.0001);	
							}
							::Copy(qh, hd().qh);
	
	
							Hydro::Initialize_QTF(q, hd().Nb, int(hd().qh.size()), hd().Nf);
							hd().mdtype = hd().qtftype = 9;
						}
						diffFullQTF = false;
						
						for (int row = 0; row < mat.size(); ++row) {
							if (mat[row].size() != 16)
								throw Exc(in.Str() + "\n"  + t_("Wrong data in list"));
							double w1 = mat[row][0],
								   w2 = mat[row][1],
								   h1 = mat[row][2],
								   h2 = mat[row][3];
								   
							if (hd().dataFromW) {
								if (!rad_s) {
									w1 *= 2*M_PI;
									w2 *= 2*M_PI;
								}
							} else {
								w1 = 2*M_PI/w1;
								w2 = 2*M_PI/w2;
							}
							int ifr1 = FindDelta(hd().qw, w1, 0.0001),
								ifr2 = FindDelta(hd().qw, w2, 0.0001),
								ih = FindDelta(hd().qh, std::complex<double>(h1, h2), 0.0001);	
							if (ifr1 < 0)
								throw Exc(in.Str() + "\n"  + Format(t_("Wrong frequency '%d' in QTF"), w1));
							if (ifr2 < 0)
								throw Exc(in.Str() + "\n"  + Format(t_("Wrong frequency '%d' in QTF"), w2));
							if (ih < 0)
								throw Exc(in.Str() + "\n"  + Format(t_("Wrong head (%d,%d) in QTF"), h1, h2));
							for (int idf = 0; idf < 6; ++idf) {
								double mag = mat[row][4+idf*2]*factor.F(idf);
								double ph  = ToRad(mat[row][4+idf*2+1]);
								q[ib][ih][idf](ifr1, ifr2) = std::polar<double>(mag, ph);
								q[ib][ih][idf](ifr2, ifr1) = std::polar<double>(mag, ph*phmult);	// Diagonal
							}
						}
					}
				} else if (fy.FirstIs("FrequencyDependentAddedMassAndDamping")) {
					if (fy.FirstMatch("AddedMassMatrixX*")) {
						idf = fy.Index()-1;
						if (idf < -1 || idf >= hd().w.size())		// Infinity is the first
							throw Exc(in.Str() + "\n" + t_("Wrong frequency"));			
						
						UVector<UVector<double>> mat = fy.GetMatrixDouble();
						
						if (idf == -1) {
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) 
									hd().Ainf(r, c) = mat[r][c]*factor.A(r, c); 
							}
						} else {
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) 
									hd().A[r][c](idf) = mat[r][c]*factor.A(r, c); 
							}
						}
					} else if (fy.FirstMatch("DampingX*")) {
						idf = fy.Index()-1;
						if (idf < -1 || idf >= hd().w.size())		// Infinity is the first
							throw Exc(in.Str() + "\n" + t_("Wrong frequency"));			
						//idf--;
						
						if (idf >= 0) {
							UVector<UVector<double>> mat = fy.GetMatrixDouble();
							
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) 
									hd().B[r][c](idf) = mat[r][c]*factor.B(r, c); 
							}
						}
					}
				}
			}
		} else if (fy.FirstIs("MultibodyGroups")) {
			if (fy.FirstIs("Bodies")) {
				if (fy.FirstIs("Name")) {
					ib = fy.Index();
					hd().msh[ib].name = fy.GetVal();
				} else if (fy.FirstIs("DisplacedVolume")) 
					hd().msh[ib].Vo = ScanDouble(fy.GetVal())*factor.len*factor.len*factor.len;
				else if (fy.FirstIs("CentreOfBuoyancy")) {
					UVector<double> line = fy.GetVectorDouble();
					
					hd().msh[ib].cb.x = line[0]*factor.len;
					hd().msh[ib].cb.y = line[1]*factor.len;
					hd().msh[ib].cb.z = line[2]*factor.len;
				} else if (fy.FirstMatch("HydrostaticStiffnessz*")) {
					UVector<UVector<double>> mat = fy.GetMatrixDouble();
					
					for (int r = 0; r < 3; ++r)				// Only heave, roll, pitch
						for (int c = 0; c < 3; ++c)
							hd().msh[ib].C(r+2, c+2) = mat[r][c]*factor.K(r, c);
				}
			} else if (fy.FirstIs("MultibodyAddedMassAndDamping")) {
				if (fy.FirstIs("AMDPeriodOrFrequency")) {
					idf = fy.Index()-1;
					if (idf < -1 || idf >= hd().w.size())		// Infinity is the first
						throw Exc(in.Str() + "\n" + t_("Wrong frequency"));	
						/*
					if (fy.GetVal() == "Infinity")
						idf = -1;
					else {
						idf = Find(hd().w, ScanDouble(fy.GetVal()));
						if (idf < 0)
							throw Exc(in.Str() + "\n" + t_("Wrong frequency"));	
					}*/
				} else if (fy.FirstIs("Matrices")) {
					if (fy.FirstIs("Row")) 
						row = ScanInt(fy.GetVal())-1;
					else if (fy.FirstIs("Column")) 
						col = ScanInt(fy.GetVal())-1;
					else if (fy.FirstMatch("AddedMassX*")) {
						UVector<UVector<double>> mat = fy.GetMatrixDouble();
						
						if (idf == -1) {
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) {
									hd().Ainf(r+row*6, c+col*6) = mat[r][c]*factor.A(r, c); 
									if (row != col)		// Fill the symmetric
										hd().Ainf(r+col*6, c+row*6) = hd().Ainf(r+row*6, c+col*6); 
								}
							}
						} else {
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) {
									hd().A[r+row*6][c+col*6](idf) = mat[r][c]*factor.A(r, c); 
									if (row != col)		// Fill the symmetric
										hd().A[r+col*6][c+row*6](idf) = hd().A[r+row*6][c+col*6](idf); 
								}
							}
						}
					} else if (fy.FirstMatch("DampingX*")) {
						if (idf >= 0) {
							UVector<UVector<double>> mat = fy.GetMatrixDouble();
							
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) {
									hd().B[r+row*6][c+col*6](idf) = mat[r][c]*factor.B(r, c); 
									if (row != col)		// Fill the symmetric
										hd().B[r+col*6][c+row*6](idf) = hd().B[r+row*6][c+col*6](idf);
								}
							}
						}
					}
				}
			}
		}
	}
	
	if (hd().Nb == 0)
		throw Exc(t_("Incorrect .yml format"));
	
	// Inertia matrices have to be translated from cg to c0
	for (int ib = 0; ib < hd().Nb; ++ib) //{
		//Point3D cg(hd().cg.col(ib));
		//Point3D c0(hd().c0.col(ib));
		Surface::TranslateInertia66(hd().msh[ib].M, hd().msh[ib].cg, hd().msh[ib].cg, hd().msh[ib].c0);
	//}
}