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
		hd().lastError = e;
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
	
	int Nb = 0; 

	YmlParser fy(in);

	FileInLine::Pos fpos = in.GetPos();

	while(fy.GetLine()) {
		if (fy.FirstIs("LoadRAOCalculationMethod")) 
			throw Exc(t_("This .yml is an OrcaWave case"));
		else if (fy.FirstIs("General")) {
			if (fy.FirstIs("UnitsSystem") && fy.GetVal() != "SI") 
				throw Exc(in.Str() + "\n" + Format(t_("Only SI units are supported. Read '%s'"), fy.GetVal()));
		} else if (fy.FirstIs("Environment")) {
			if (fy.FirstIs("WaterSurfaceZ") && fy.GetVal() != "0") 
				throw Exc(in.Str() + "\n" + Format(t_("Only WaterSurfaceZ 0 is supported. Read '%s'"), fy.GetVal()));
		} else if (fy.FirstIs("VesselTypes")) {
			if (fy.FirstIs("Name")) {
				if (fy.Index() != Nb)
					throw Exc(in.Str() + "\n" + t_("Failed body count"));
				hd().names << fy.GetVal();
				Nb++;
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
				if (fy.FirstIs("DisplacementRAOs")) {
					if (fy.FirstIs("RAOOrigin")) {
						if (!IsEqualRange<UVector<double>>({0,0,0}, fy.GetVectorDouble()))
							throw Exc(in.Str() + "\n"  + Format(t_("Only RAOOrigin:[0,0,0] is supported. Read '%s'"), fy.StrVar()));	
					} else if (fy.FirstIs("PhaseOrigin")) {
						if (!IsEqualRange<UVector<double>>({0,0,0}, fy.GetVectorDouble()))
							throw Exc(in.Str() + "\n"  + Format(t_("Only PhaseOrigin:[0,0,0] is supported. Read '%s'"), fy.StrVar()));
					} else if (fy.FirstIs("RAOs")) {
						if (fy.FirstIs("RAODirection") && fy.GetIndex()[1] == 0) {		// Only for the first body
							if (hd().head.size() != fy.Index())
								throw Exc(in.Str() + "\n" + t_("Failed headings count"));			
							hd().head << ScanDouble(fy.GetVal());
						}
					}
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
			if (fy.FirstIs("MultibodyAddedMassAndDamping")) {
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
	
	if (Nb == 0)
		throw Exc(S("\n") + t_("No body found"));
	
	hd().Nb = Nb;
	hd().Vo.SetCount(hd().Nb, NaNDouble);
	hd().cg.setConstant(3, hd().Nb, NaNDouble);
	hd().c0.setConstant(3, hd().Nb, 0);				// OrcaWave reference is 0,0,0
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
	hd().A.SetCount(6*hd().Nb);
	hd().B.SetCount(6*hd().Nb);
	for (int i = 0; i < 6*hd().Nb; ++i) {
		hd().A[i].SetCount(6*hd().Nb);
		hd().B[i].SetCount(6*hd().Nb);
		for (int j = 0; j < 6*hd().Nb; ++j) {
			hd().A[i][j].setConstant(hd().Nf, NaNDouble);// In Wamit, unloaded DOFs are considered negligible	
			hd().B[i][j].setConstant(hd().Nf, NaNDouble);	
		}
	}
		
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

	int ib = -1, row = -1, col = -1;
	int idf = -1;
	
	in.SeekPos(fpos);
		
	while(fy.GetLine()) {
		if (fy.FirstIs("Environment")) {
			if (fy.FirstIs("Density")) 
				hd().rho = ScanDouble(fy.GetVal())*1000;		// In kg/m3
			else if (fy.FirstIs("WaterDepth")) { 
				hd().rho = ScanDouble(fy.GetVal())*1000;		// In kg/m3
				String h = fy.GetVal(); 
				if (ToLower(h) == "infinite")
					hd().h = -1;
				else	
					hd().h = ScanDouble(h);
			}
		} else if (fy.FirstIs("VesselTypes")) {
			if (fy.FirstIs("Draughts")) {
				ib = fy.GetIndex()[1];
				Eigen::MatrixXd &inertia = hd().M[ib];
				if (fy.FirstIs("Mass")) 
					inertia(0, 0) = inertia(1, 1) = inertia(2, 2) = ScanDouble(fy.GetVal())*1000;
				else if (fy.FirstMatch("MomentOfInertiaTensor*")) {
					UVector<UVector<double>> mat = fy.GetMatrixDouble();
					
					inertia(3, 3) = mat[0][0]*1000;				// In kg
					inertia(3, 4) = mat[0][1]*1000;
					inertia(3, 5) = mat[0][2]*1000;
					inertia(4, 3) = mat[1][0]*1000;
					inertia(4, 4) = mat[1][1]*1000;
					inertia(4, 5) = mat[1][2]*1000;
					inertia(5, 3) = mat[2][0]*1000;
					inertia(5, 4) = mat[2][1]*1000;
					inertia(5, 5) = mat[2][2]*1000;
				} else if (fy.FirstMatch("HydrostaticStiffnessz*")) {
					UVector<UVector<double>> mat = fy.GetMatrixDouble();
					
					for (int r = 0; r < 3; ++r)				// Only heave, roll, pitch
						for (int c = 0; c < 3; ++c)
							hd().C[ib](r+2, c+2) = mat[r][c]*1000;
				} else if (fy.FirstIs("CentreOfMass")) {
					UVector<double> line = fy.GetVectorDouble();
					
					hd().cg(0, ib) = line[0];
					hd().cg(1, ib) = line[1];
					hd().cg(2, ib) = line[2];
				} else if (fy.FirstIs("CentreOfBuoyancy")) {
					UVector<double> line = fy.GetVectorDouble();
					
					hd().cb(0, ib) = line[0];
					hd().cb(1, ib) = line[1];
					hd().cb(2, ib) = line[2];
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
									hd().rao.force[idh](ifr, idof+6*ib) = std::polar<double>(mat[ifr][1 + 2*idof], ToRad(mat[ifr][1 + 2*idof + 1]));
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
									hd().ex.force[idh](ifr, idof+6*ib) = std::polar<double>(mat[ifr][1 + 2*idof]*10, ToRad(mat[ifr][1 + 2*idof + 1]));
						}
					}
				} else if (fy.FirstIs("WaveDrift")) {
					if (fy.FirstIs("RAOs")) {
						if (fy.FirstMatch("RAOPeriodOrFrequency*")) {	
							if (hd().md.size() == 0) {
								hd().mdhead.resize(hd().head.size());
								for (int ih = 0; ih < hd().head.size(); ++ih)
									hd().mdhead[ih] = std::complex<double>(hd().head[ih], hd().head[ih]);
								Hydro::InitMD(hd().md, hd().Nb, int(hd().mdhead.size()), hd().Nf);
							}
							
							int idh = fy.Index();
							if (idh < 0 || idh >= hd().head.size())
								throw Exc(in.Str() + "\n" + t_("Wrong heading"));
								
							UVector<UVector<double>> mat = fy.GetMatrixDouble();
							
							if (mat.size() != hd().Nf || mat[0].size() != 7)
								throw Exc(in.Str() + "\n"  + Format(t_("Wrong number of numbers in WaveDrift matrix"), fy.GetText()));

							for (int ifr = 0; ifr < hd().Nf; ++ifr) 
								for (int idof = 0; idof < 6; ++idof) 
									hd().md[ib][idh][idof](ifr) = mat[ifr][1 + idof]*10;
						}
					}
				} else if (fy.FirstIs("FrequencyDependentAddedMassAndDamping")) {
					if (fy.FirstMatch("AddedMassMatrixX*")) {
						idf = fy.Index()-1;
						if (idf < -1 || idf >= hd().w.size())		// Infinity is the first
							throw Exc(in.Str() + "\n" + t_("Wrong frequency"));			
						//idf--;
						
						UVector<UVector<double>> mat = fy.GetMatrixDouble();
						
						if (idf == -1) {
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) 
									hd().Ainf(r, c) = mat[r][c]*100000; 
							}
						} else {
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) 
									hd().A[r][c](idf) = mat[r][c]*100000; 
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
									hd().B[r][c](idf) = mat[r][c]*100000; 
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
					hd().Vo[ib] = ScanDouble(fy.GetVal());
				else if (fy.FirstIs("CentreOfBuoyancy")) {
					UVector<double> line = fy.GetVectorDouble();
					
					hd().cb(0, ib) = line[0];
					hd().cb(1, ib) = line[1];
					hd().cb(2, ib) = line[2];
				} else if (fy.FirstMatch("HydrostaticStiffnessz*")) {
					UVector<UVector<double>> mat = fy.GetMatrixDouble();
					
					for (int r = 0; r < 3; ++r)				// Only heave, roll, pitch
						for (int c = 0; c < 3; ++c)
							hd().C[ib](r+2, c+2) = mat[r][c]*1000;
				}
			} else if (fy.FirstIs("MultibodyAddedMassAndDamping")) {
				if (fy.FirstIs("AMDPeriodOrFrequency")) {
					if (fy.GetVal() == "Infinity")
						idf = -1;
					else {
						idf = Find(hd().w, ScanDouble(fy.GetVal()));
						if (idf < 0)
							throw Exc(in.Str() + "\n" + t_("Wrong frequency"));	
					}
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
									hd().Ainf(r+row*6, c+col*6) = mat[r][c]*100000; 
							}
						} else {
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) 
									hd().A[r+row*6][c+col*6](idf) = mat[r][c]*100000; 
							}
						}
					} else if (fy.FirstMatch("DampingX*")) {
						if (idf >= 0) {
							UVector<UVector<double>> mat = fy.GetMatrixDouble();
							
							for (int r = 0; r < 6; ++r) {
								for (int c = 0; c < 6; ++c) 
									hd().B[r+row*6][c+col*6](idf) = mat[r][c]*100000; 
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