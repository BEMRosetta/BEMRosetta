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
		if (!Load_YML()) 
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

bool OrcaWave::Load_YML() {
	
	auto GetVector = [](FieldSplit &f, String str = "")->UVector<double> {
		UVector<String> list;

		while (true) {
			if (list.IsEmpty()) {
				int id = str.Find("[");
				if (id >= 0) {
					str = str.Mid(id+1);
					list = Split(str, ",");
					str.Clear();
				}
			} else {
				list.Append(Split(str, ","));
				str.Clear();
			}
			if (!list.IsEmpty() && Last(list).EndsWith("]")) 
				break;
			str << Trim(f.GetLine());
		}
		
		UVector<double> ret;
		for (String &s : list)
			ret << ScanDouble(s);
		return ret;
	};
	
	String fileName = ForceExt(hd().file, ".yml");
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		return false;
	
	String line; 
	FieldSplit f(in);
	f.IsSeparator = IsTabSpace;
	
	hd().dataFromW = true;
	bool rad_s = true;
	hd().dimen = true;
	
	hd().Nb = hd().Nf = hd().Nh = Null;
	
	int Nb = 0; 

	FileInLine::Pos fpos = in.GetPos();
	
	bool invessel = false;
	while(!in.IsEof()) {
		f.GetLine_discard_empty();
		if (f.GetText(0) == "LoadRAOCalculationMethod:") 		// Only in OrcaWave cases
			throw Exc(t_("This .yml is an OrcaWave case"));
		else if (f.GetText(0) == "UnitsSystem:") {
			if (f.GetText(1) != "SI")
				throw Exc(in.Str() + "\n"  + Format(t_("Only SI units are supported %s"), f.GetText()));
		} else if (f.GetText(0) == "WaterSurfaceZ:") {
			if (f.GetText(1) != "0")
				throw Exc(in.Str() + "\n"  + Format(t_("Only WaterSurfaceZ 0 is supported %s"), f.GetText()));
		} else if (f.GetText(0) == "WavesReferredToBy:") {
			if (f.GetText(1) == "frequency")
				hd().dataFromW = true;
			else if (f.GetText(1) == "period")
				hd().dataFromW = false;
			else
				throw Exc(in.Str() + "\n"  + Format(t_("Unknown data in WavesReferredToBy: %s"), f.GetText()));
			if (f.GetText(2) == "(rad/s)")
				rad_s = true;
			else if (f.GetText(2) == "(Hz)")	 
				rad_s = false;
			else if (f.GetText(2) == "(s)")	 	
				;
			else
				throw Exc(in.Str() + "\n"  + Format(t_("Unknown data in WavesReferredToBy: %s"), f.GetText()));
		} else if (f.GetText(0) == "SurgePositive:") {
			if (f.GetText(1) != "forward")
				throw Exc(in.Str() + "\n"  + Format(t_("Only SurgePositive: forward is supported %s"), f.GetText()));
		} else if (f.GetText(0) == "SwayPositive:") {
			if (f.GetText(1) != "port")
				throw Exc(in.Str() + "\n"  + Format(t_("Only SwayPositive: port is supported %s"), f.GetText()));
		} else if (f.GetText(0) == "HeavePositive:") {
			if (f.GetText(1) != "up")
				throw Exc(in.Str() + "\n"  + Format(t_("Only HeavePositive: up is supported %s"), f.GetText()));
		} else if (f.GetText(0) == "RollPositiveStarboard:") {
			if (f.GetText(1) != "down")
				throw Exc(in.Str() + "\n"  + Format(t_("Only RollPositiveStarboard: down is supported %s"), f.GetText()));
		} else if (f.GetText(0) == "PitchPositiveBowe:") {
			if (f.GetText(1) != "down")
				throw Exc(in.Str() + "\n"  + Format(t_("Only PitchPositiveBow: down is supported %s"), f.GetText()));
		} else if (f.GetText(0) == "YawPositiveBow:") {
			if (f.GetText(1) != "port")
				throw Exc(in.Str() + "\n"  + Format(t_("Only YawPositiveBow: port is supported %s"), f.GetText()));
		} else if (f.GetText(0) == "QTFConventionsRotationOrder:") {
			if (f.GetText(1) != "RzRyRx")
				throw Exc(in.Str() + "\n"  + Format(t_("Only QTFConventionsRotationOrder: RzRyRx is supported %s"), f.GetText()));
		} else if (f.GetText(0) == "QTFConventionsRotationAxes:") {
			if (f.GetText(1) != "original")
				throw Exc(in.Str() + "\n"  + Format(t_("Only QTFConventionsRotationAxes: original is supported %s"), f.GetText()));
		} else if (f.GetText(0) == "QTFConventionsFrameOfReference:") {
			if (f.GetText(1) != "earth")
				throw Exc(in.Str() + "\n"  + Format(t_("Only QTFConventionsFrameOfReference: earth is supported %s"), f.GetText()));
		} else if (f.GetText(0) == "RAOOrigin:") {
			String line = f.GetText();
			line.Replace(" ", "");
			line.Replace("\t", "");
			if (line != "RAOOrigin:[0,0,0]")
				throw Exc(in.Str() + "\n"  + Format(t_("Only RAOOrigin:[0,0,0] is supported %s"), f.GetText()));
		} else if (f.GetText(0) == "PhaseOrigin:") {
			String line = f.GetText();
			line.Replace(" ", "");
			line.Replace("\t", "");
			if (line != "PhaseOrigin:[0,0,0]" && line != "PhaseOrigin:[~,~,~]")
				throw Exc(in.Str() + "\n"  + Format(t_("Only PhaseOrigin:[0,0,0] is supported %s"), f.GetText()));
		} else if (f.GetText(0) == "OtherDampingOrigin:") {
			String line = f.GetText();
			line.Replace(" ", "");
			line.Replace("\t", "");
			if (line != "OtherDampingOrigin:[0,0,0]")
				throw Exc(in.Str() + "\n"  + Format(t_("Only OtherDampingOrigin:[0,0,0] is supported %s"), f.GetText()));
		} else if (f.GetText(0) == "ReferenceOrigin:") {
			String line = f.GetText();
			line.Replace(" ", "");
			line.Replace("\t", "");
			if (line != "ReferenceOrigin:[0,0,0]")
				throw Exc(in.Str() + "\n"  + Format(t_("Only ReferenceOrigin:[0,0,0] is supported %s"), f.GetText()));
		} else if (f.GetText(0) == "ReferenceOriginDatumPosition:") {
			String line = f.GetText();
			line.Replace(" ", "");
			line.Replace("\t", "");
			if (line != "ReferenceOriginDatumPosition:[0,0,0]")
				throw Exc(in.Str() + "\n"  + Format(t_("Only ReferenceOriginDatumPosition:[0,0,0] is supported %s"), f.GetText()));
		} else if (f.GetText(0) == "VesselTypes:")
			invessel = true;
		else if (invessel && f.GetText().StartsWith("      - Name:"))
			Nb++;
		else if (f.GetText(0) == "-" && f.GetText(1) == "RAODirection:")
			FindAdd(hd().head, f.GetDouble(2));
		else if (f.GetText(0) == "-" && f.GetText(1) == "AMDPeriodOrFrequency:") {
			if (f.GetText(2) != "Infinity")
				hd().w << f.GetDouble(2);
		} else if (f.GetText(0) == "Vessels:")
			break;
	}
						
	if (Nb == 0)
		throw Exc(S("\n") + t_("No body found"));
	
	hd().Nb = Nb;
	hd().names.SetCount(hd().Nb);
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

	in.SeekPos(fpos);

	char inforce = '\0';
	invessel = false;
	int ib = -1;
	int idFreq = -2;		// First is infinity
	int idHead = -1;
	
	double mass = Null;
	while(!f.IsEof()) {
		f.GetLine_discard_empty();
		if (f.size() == 0)
			;
		else if (f.GetText(0) == "Density:") 
			hd().rho = f.GetDouble(1)*1000;		// In kg/m3
		else if (f.GetText(0) == "WaterDepth:") {
			String h = f.GetText(1); 
			if (ToLower(h) == "infinite")
				hd().h = -1;
			else	
				hd().h = ScanDouble(h);
		} else if (f.GetText(0) == "VesselTypes:")
			invessel = true;
		else if (invessel && f.GetText().StartsWith("      - Name:")) {
			ib++;	
			hd().names[ib] = f.GetText(2);		
		} else if (f.GetText(0) == "Mass:") 
			mass = f.GetDouble(1)*1000;
		else if (invessel && f.GetText().StartsWith("        MomentOfInertiaTensor")) {	
			Eigen::MatrixXd &inertia = hd().M[ib];
			
			inertia(0, 0) = inertia(1, 1) = inertia(2, 2) = mass;
			
			UVector<double> line;
			
			line = GetVector(f);			
			inertia(3, 3) = line[0]*1000;
			inertia(3, 4) = line[1]*1000;
			inertia(3, 5) = line[2]*1000;
			line = GetVector(f);			
			inertia(4, 3) = line[0]*1000;
			inertia(4, 4) = line[1]*1000;
			inertia(4, 5) = line[2]*1000;
			line = GetVector(f);			
			inertia(5, 3) = line[0]*1000;
			inertia(5, 4) = line[1]*1000;
			inertia(5, 5) = line[2]*1000;
		} else if (invessel && f.GetText().StartsWith("        CentreOfMass:")) {
			Eigen::MatrixXd &inertia = hd().M[ib];
			
			UVector<double> line = GetVector(f, f.GetText().Mid(21));
			
			hd().cg(0, ib) = line[0];
			hd().cg(1, ib) = line[1];
			hd().cg(2, ib) = line[2];

			double cx = mass*hd().cg(0, ib);
			double cy = mass*hd().cg(1, ib);
			double cz = mass*hd().cg(2, ib);
			inertia(1, 5) = inertia(5, 1) =  cx;
			inertia(2, 4) = inertia(4, 2) = -cx;
			inertia(0, 5) = inertia(5, 0) = -cy;
			inertia(2, 3) = inertia(3, 2) =  cy;
			inertia(0, 4) = inertia(4, 0) =  cz;
			inertia(1, 3) = inertia(3, 1) = -cz;
		} else if (invessel && f.GetText().StartsWith("          - AMDPeriodOrFrequency:")) 
			idFreq++; 	
		else if (invessel && f.GetText().StartsWith("            AddedMassMatrixX")) {
			if (idFreq == -1) {
				for (int r = 0; r < 6; ++r) {
					UVector<double> line = GetVector(f); 
					for (int c = 0; c < 6; ++c) 
						hd().Ainf(r, c) = line[c]*1000; 
				}
			} else {
				for (int r = 0; r < 6; ++r) {
					UVector<double> line = GetVector(f); 
					for (int c = 0; c < 6; ++c) 
						hd().A[r][c](idFreq) = line[c]*1000; 
				}
			}
		} else if (invessel && f.GetText().StartsWith("            DampingX")) {
			if (idFreq >= 0) {
				for (int r = 0; r < 6; ++r) {
					UVector<double> line = GetVector(f); 
					for (int c = 0; c < 6; ++c) 
						hd().B[r][c](idFreq) = line[c]*1000; 
				}
			}
		} else if (f.GetText().StartsWith("        DisplacementRAOs:")) {
			inforce = 'r';		
			idHead = -1;
		} else if (f.GetText().StartsWith("        LoadRAOs:")) {
			inforce = 'f';		
			idHead = -1;
		} else if (f.GetText().StartsWith("        WaveDrift:")) {
			inforce = 'd';		
			idHead = -1;
			if (hd().md.size() == 0) {
				hd().mdhead.resize(hd().head.size());
				for (int ih = 0; ih < hd().head.size(); ++ih)
					hd().mdhead[ih] = std::complex<double>(hd().head[ih], hd().head[ih]);
				Hydro::InitMD(hd().md, hd().Nb, int(hd().mdhead.size()), hd().Nf);
			}
		} else if (f.GetText().StartsWith("              RAOPeriodOrFrequency")) {
			idHead++;
			for (int ifr = 0; ifr < hd().Nf; ++ifr) {
				UVector<double> line = GetVector(f); 
				
				if (inforce == 'f') {
					if (line.size() != 13)
						throw Exc(in.Str() + "\n"  + Format(t_("Wrong number of numbers in line %s"), f.GetText()));
					for (int idof = 0; idof < 6; ++idof) 
						hd().ex .force[idHead](ifr, idof + 6*ib) = std::polar<double>(line[1 + 2*idof]*1000, ToRad(line[1 + 2*idof + 1]));
				} else if (inforce == 'r') {
					if (line.size() != 13)
						throw Exc(in.Str() + "\n"  + Format(t_("Wrong number of numbers in line %s"), f.GetText()));
					for (int idof = 0; idof < 6; ++idof) 
						hd().rao.force[idHead](ifr, idof + 6*ib) = std::polar<double>(line[1 + 2*idof], ToRad(line[1 + 2*idof + 1]));
				} else if (inforce == 'd') {
					if (line.size() != 7)
						throw Exc(in.Str() + "\n"  + Format(t_("Wrong number of numbers in line %s"), f.GetText()));
					for (int idof = 0; idof < 6; ++idof) 
						hd().md[ib][idHead][idof](ifr) = line[1 + idof]*1000;
				} else
					throw Exc(in.Str() + "\n"  + Format(t_("Wrong RAOPeriodOrFrequency status %s"), f.GetText()));
			}
		}
	}
	if (hd().Nb == 0)
		throw Exc(t_("Incorrect .yml format"));
	
	return true;	
}