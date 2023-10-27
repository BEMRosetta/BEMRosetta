// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"


bool OrcaWave::Load(String file, double) {
	hd().file = file;
	hd().name = GetFileTitle(file);
	hd().dimen = true;
	hd().len = 1;
	hd().code = Hydro::ORCAWAVE;
	hd().Nb = Null;
	
	try {
		BEM::Print("\n\n" + Format(t_("Loading '%s'"), file));

		BEM::Print("\n- " + S(t_("YML file")));
		if (!Load_YML_Res()) 
			BEM::PrintWarning(S(": ** YML file ") + t_("Not found") + "**");
		
		if (IsNull(hd().Nb))
			return false;
	
		hd().dof.Clear();	hd().dof.SetCount(hd().Nb, 0);
		for (int i = 0; i < hd().Nb; ++i)
			hd().dof[i] = 6;
	} catch (Exc e) {
		BEM::PrintError(Format("\n%s: %s", t_("Error"), e));
		//hd().lastError = e;
		return false;
	}
	
	return true;
}

bool OrcaWave::Load_YML_Res() {
	String fileName = ForceExt(hd().file, ".yml");
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		return false;
	
	hd().dataFromW = true;
	bool rad_s = true;
	hd().dimen = true;
	
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
		
	double factorMass = 1, factorLen = 1, factorForce = 1;	
	
	int ib = -1;
	int Nb = 0; 

	YmlParser fy(in);

	FileInLine::Pos fpos = in.GetPos();

	while(fy.GetLine()) {
		if (fy.FirstIs("LoadRAOCalculationMethod")) 
			throw Exc(t_("This .yml is an OrcaWave case"));
		else if (fy.FirstIs("General")) {
			if (fy.FirstIs("UnitsSystem")) {
				if (fy.GetVal() == "SI") {
					factorMass = factorForce = 1000;
					factorLen = 1;	
				} else if (fy.GetVal() == "User") 
					;
				else
					throw Exc(in.Str() + "\n" + Format(t_("Only SI and User units are supported. Read '%s'"), fy.GetVal()));
			} else if (fy.FirstIs("LengthUnits")) {
				if (fy.GetVal() == "m") 
					factorLen = 1;	
				else if (fy.GetVal() == "mm") 
					factorLen = 1E-3;
				else if (fy.GetVal() == "cm") 
					factorLen = 1E-2;
				else if (fy.GetVal() == "km") 
					factorLen = 1000;
				else
					throw Exc(in.Str() + "\n" + Format(t_("This length unit is not supported. Read '%s'"), fy.GetVal()));
			} else if (fy.FirstIs("MassUnits")) {
				if (fy.GetVal() == "kg") 
					factorMass = 1;	
				else if (fy.GetVal() == "te") 
					factorMass = 1000;
				else
					throw Exc(in.Str() + "\n" + Format(t_("This mass unit is not supported. Read '%s'"), fy.GetVal())); 
			} else if (fy.FirstIs("ForceUnits")) {
				if (fy.GetVal() == "N") 
					factorForce = 1;	
				else if (fy.GetVal() == "kN") 
					factorForce = 1000;
				else if (fy.GetVal() == "MN") 
					factorForce = 1E6;
				else
					throw Exc(in.Str() + "\n" + Format(t_("This force unit is not supported. Read '%s'"), fy.GetVal())); 
			} else if (fy.FirstIs("g")) 
				hd().g = ScanDouble(fy.GetVal())*factorLen;
		} else if (fy.FirstIs("Environment")) {
			if (fy.FirstIs("WaterSurfaceZ") && fy.GetVal() != "0") 
				throw Exc(in.Str() + "\n" + Format(t_("Only WaterSurfaceZ 0 is supported. Read '%s'"), fy.GetVal()));
		} else if (fy.FirstIs("VesselTypes")) {
			if (fy.FirstIs("Name")) {
				if (fy.Index() != Nb)
					throw Exc(in.Str() + "\n" + t_("Failed body count"));
				hd().names << fy.GetVal();
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
	
	auto factorA = [&](int r, int c)->double {
		if (r < 3 && c < 3)
			return factorMass;
		else if (r >= 3 && c >= 3)
			return factorMass*factorLen*factorLen;
		else
			return factorMass*factorLen;
	};
	auto factorK = [&](int r, int c)->double {
		if (r < 3 && c < 3)
			return factorForce/factorLen;
		if (r < 3)
			return factorForce;
		if (c < 3)
			return factorForce;
		return factorForce*factorLen;
	};
	auto factorB = [&](int r, int c)->double {
		if (r < 3 && c < 3)
			return factorForce/factorLen;
		if (r < 3)
			return factorForce;
		if (c < 3)
			return factorForce;
		return factorForce*factorLen;
	};
	auto factorM = [&](int r, int c)->double {
		if (r < 3 && c < 3)
			return factorMass;
		else if (r >= 3 && c >= 3)
			return factorMass*factorLen*factorLen;
		else
			return factorMass*factorLen;
	};
	auto factorF = [&](int r)->double {
		if (r < 3)
			return factorForce/factorLen;
		return factorForce;
	};
	auto factorRAO = [&](int r)->double {
		if (r < 3)
			return 1;
		return 1/factorLen;
	};
	auto factorMD = [&](int r)->double {
		if (r < 3)
			return factorForce/factorLen/factorLen;
		return factorForce/factorLen;
	}; 
	
	
	if (Nb == 0)
		throw Exc(S("\n") + t_("No body found"));
	
	hd().Nb = Nb;
	hd().Vo.SetCount(hd().Nb, NaNDouble);
	hd().cg.setConstant(3, hd().Nb, NaNDouble);
	
	hd().c0.resize(3, hd().Nb);
	for (int ib = 0; ib < Nb; ++ib)				
		for (int idf = 0; idf < 3; ++idf)
			hd().c0(idf, ib) = c0s[ib][idf];
	
	hd().cb.setConstant(3, hd().Nb, NaNDouble);
	hd().C.SetCount(hd().Nb);
	for (int ib = 0; ib < hd().Nb; ++ib) 
		hd().C[ib].setConstant(6, 6, 0);
	hd().M.SetCount(hd().Nb);
	for (int ib = 0; ib < hd().Nb; ++ib) 
		hd().M[ib].setConstant(6, 6, 0);

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
		hd().T.SetCount(hd().w.size());
		for (int i = 0; i < hd().w.size(); ++i)
			hd().T[i] = 2*M_PI/hd().w[i];	
	} else {
		hd().T = pick(hd().w);
		hd().w.SetCount(hd().T.size());
		for (int i = 0; i < hd().w.size(); ++i)
			hd().w[i] = 2*M_PI/hd().T[i];	
	}	

	hd().Initialize_Forces(hd().ex);
	hd().Initialize_Forces(hd().rao);

	int row = -1, col = -1;
	int idf = -1;
	bool diffFullQTF = false;
	
	in.SeekPos(fpos);
	
	auto FillInertia = [](Eigen::MatrixXd &inertia, double cgx, double cgy, double cgz) {
		double &mass = inertia(0, 0);
		inertia(1, 5) = inertia(5, 1) =  cgx*mass;
		inertia(2, 4) = inertia(4, 2) = -cgx*mass;
		inertia(2, 3) = inertia(3, 2) =  cgy*mass;
		inertia(0, 5) = inertia(5, 0) = -cgy*mass;
		inertia(0, 4) = inertia(4, 0) =  cgz*mass;
		inertia(1, 3) = inertia(3, 1) = -cgz*mass;
	};
		
	while(fy.GetLine()) {
		if (fy.FirstIs("Environment")) {
			if (fy.FirstIs("Density")) 
				hd().rho = ScanDouble(fy.GetVal())*factorMass/factorLen/factorLen/factorLen;		// In kg/m3
			else if (fy.FirstIs("WaterDepth")) { 
				String h = fy.GetVal(); 
				if (ToLower(h) == "infinite")
					hd().h = -1;
				else	
					hd().h = ScanDouble(h)*factorLen;
			}
		} else if (fy.FirstIs("VesselTypes")) {
			if (fy.FirstIs("Name")) 
				ib = fy.GetIndex()[1];
			else if (fy.FirstIs("Draughts")) {
				Eigen::MatrixXd &inertia = hd().M[ib];
				if (fy.FirstIs("Mass")) 
					inertia(0, 0) = inertia(1, 1) = inertia(2, 2) = ScanDouble(fy.GetVal())*factorMass;
				else if (fy.FirstMatch("MomentOfInertiaTensor*")) {
					UVector<UVector<double>> mat = fy.GetMatrixDouble();
					
					inertia(3, 3) = mat[0][0]*factorM(3, 3);		// In kg
					inertia(3, 4) = mat[0][1]*factorM(3, 4);
					inertia(3, 5) = mat[0][2]*factorM(3, 5);
					inertia(4, 3) = mat[1][0]*factorM(4, 3);
					inertia(4, 4) = mat[1][1]*factorM(4, 4);
					inertia(4, 5) = mat[1][2]*factorM(4, 5);
					inertia(5, 3) = mat[2][0]*factorM(5, 3);
					inertia(5, 4) = mat[2][1]*factorM(5, 4);
					inertia(5, 5) = mat[2][2]*factorM(5, 5);
					
					if (IsNum(hd().cg(0, ib)))
						FillInertia(inertia, hd().cg(0, ib), hd().cg(1, ib), hd().cg(2, ib));
				} else if (fy.FirstMatch("HydrostaticStiffnessz*")) {
					UVector<UVector<double>> mat = fy.GetMatrixDouble();
					
					for (int r = 0; r < 3; ++r)				// Only heave, roll, pitch
						for (int c = 0; c < 3; ++c)
							hd().C[ib](r+2, c+2) = mat[r][c]*factorK(r+2, c+2);
				} else if (fy.FirstIs("CentreOfMass")) {
					UVector<double> line = fy.GetVectorDouble();
					
					hd().cg(0, ib) = line[0]*factorLen;
					hd().cg(1, ib) = line[1]*factorLen;
					hd().cg(2, ib) = line[2]*factorLen;
					
					if (inertia(0, 0) > 0)
						FillInertia(inertia, hd().cg(0, ib), hd().cg(1, ib), hd().cg(2, ib));
				} else if (fy.FirstIs("CentreOfBuoyancy")) {
					UVector<double> line = fy.GetVectorDouble();
					
					hd().cb(0, ib) = line[0]*factorLen;
					hd().cb(1, ib) = line[1]*factorLen;
					hd().cb(2, ib) = line[2]*factorLen;
				} else if (fy.FirstIs("DisplacedVolume")) 
					hd().Vo[ib] = ScanDouble(fy.GetVal());
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
									hd().rao.force[idh](ifr, idof+6*ib) = std::polar<double>(mat[ifr][1 + 2*idof]*factorRAO(idof), ToRad(mat[ifr][1 + 2*idof + 1]));
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
									hd().ex.force[idh](ifr, idof+6*ib) = std::polar<double>(mat[ifr][1 + 2*idof]*factorF(idof), ToRad(mat[ifr][1 + 2*idof + 1]));
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
									hd().md[ib][idh][idof](ifr) = mat[ifr][1 + idof]*factorMD(idof);
						}
					}
				} else if (fy.FirstIs("SumFrequencyQTFs") || (diffFullQTF && fy.FirstIs("WaveDrift"))) {
					if (fy.FirstMatch("RAOPeriodOrFrequency*")) {
						UVector<UVector<double>> mat = fy.GetMatrixDouble();
						
						UArray<UArray<UArray<MatrixXcd>>> &q = diffFullQTF ? hd().qtfdif : hd().qtfsum;
						double phmult = !diffFullQTF ? 1 : -1;		// Difference is conjugate-symmetric
						
						if (!hd().IsLoadedQTF(!diffFullQTF)) {		// Gets frequencies and headings
							hd().qw.resize(hd().Nf);
							for (int iw = 0; iw < hd().Nf; ++iw) 
								hd().qw[iw] = hd().w[iw];
							
							UArray<std::complex<double>> qh;
							for (int row = 0; row < mat.size(); ++row) {
								if (mat[row].size() != 16)
									throw Exc(in.Str() + "\n"  + t_("Wrong data in list"));
								double h1 = mat[row][2],
									   h2 = mat[row][3];
									   
								FindAddDelta(qh, std::complex<double>(h1, h2), 0.0001);	
							}
							Copy(qh, hd().qh);
	
	
							Hydro::Initialize_QTF(q, hd().Nb, int(hd().qh.size()), hd().Nf);
							hd().mdtype = 9;
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
								double mag = mat[row][4+idf*2]*factorF(idf);
								double ph  = ToRad(mat[row][4+idf*2+1]);
								q[ib][ih][idf](ifr1, ifr2) = std::polar<double>(mag, ph);
								q[ib][ih][idf](ifr2, ifr1) = std::polar<double>(mag, ph*phmult);
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
									hd().Ainf(r, c) = mat[r][c]*factorA(r, c); 
							}
						} else {
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) 
									hd().A[r][c](idf) = mat[r][c]*factorA(r, c); 
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
									hd().B[r][c](idf) = mat[r][c]*factorB(r, c); 
							}
						}
					}
				}
			}
		} else if (fy.FirstIs("MultibodyGroups")) {
			if (fy.FirstIs("Bodies")) {
				if (fy.FirstIs("Name")) {
					ib = fy.Index();
					hd().names[ib] = fy.GetVal();
				} else if (fy.FirstIs("DisplacedVolume")) 
					hd().Vo[ib] = ScanDouble(fy.GetVal())*factorLen*factorLen*factorLen;
				else if (fy.FirstIs("CentreOfBuoyancy")) {
					UVector<double> line = fy.GetVectorDouble();
					
					hd().cb(0, ib) = line[0]*factorLen;
					hd().cb(1, ib) = line[1]*factorLen;
					hd().cb(2, ib) = line[2]*factorLen;
				} else if (fy.FirstMatch("HydrostaticStiffnessz*")) {
					UVector<UVector<double>> mat = fy.GetMatrixDouble();
					
					for (int r = 0; r < 3; ++r)				// Only heave, roll, pitch
						for (int c = 0; c < 3; ++c)
							hd().C[ib](r+2, c+2) = mat[r][c]*factorK(r, c);
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
								for (int c = 0; c < 6; ++c) 
									hd().Ainf(r+row*6, c+col*6) = mat[r][c]*factorA(r, c); 
							}
						} else {
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) 
									hd().A[r+row*6][c+col*6](idf) = mat[r][c]*factorA(r, c); 
							}
						}
					} else if (fy.FirstMatch("DampingX*")) {
						if (idf >= 0) {
							UVector<UVector<double>> mat = fy.GetMatrixDouble();
							
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) 
									hd().B[r+row*6][c+col*6](idf) = mat[r][c]*factorB(r, c); 
							}
						}
					}
				}
			}
		}
	}

	if (hd().Nb == 0)
		throw Exc(t_("Incorrect .yml format"));
	
	return true;	
}