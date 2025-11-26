// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2024, the BEMRosetta author and contributors
#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <SurfaceCanvas/SurfaceCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <TabBar/TabBar.h>
#include <DropGrid/DropGrid.h>

#include <BEMRosetta_cl/BEMRosetta.h>
#include <STEM4U/SeaWaves.h>

using namespace Upp;

#include "main.h"

void MainTools::Init() {
	Add(tab.SizePos());
	
	waveInfo.Init();
	tab.Add(waveInfo.SizePos(), t_("Waves calculator"));
	mooringInfo.Init();
	tab.Add(mooringInfo.SizePos(), t_("Chain & Rope calculator"));
#ifdef PLATFORM_WIN32
	orcaLicense.Init();
	tab.Add(orcaLicense.SizePos(), t_("OrcaWave license"));
#endif 
};


bool Tools_WaveInfo::Init() {
	CtrlLayout(*this);
	
	gridResultsT.MultiSelect().SetLineCy(EditField::GetStdHeight());
	gridResultsT.WhenBar = THISBACK(OnArrayBar);
	gridResultsT.AddColumn(t_("Parameter"));
	gridResultsT.AddColumn(t_("Value"));
	gridResultsT.Add(t_("Celerity (C) [m/s]"));
	gridResultsT.Add(t_("Conditions"));
	gridResultsT.Add(t_("Wave number (k) [rad/m]"));
	gridResultsT.Add(t_("Wave length (L) [m]"));
	gridResultsT.Add(t_("Group Celerity [m/s]"));
	
	gridResultsH.MultiSelect().SetLineCy(EditField::GetStdHeight());
	gridResultsH.WhenBar = THISBACK(OnArrayBar);
	gridResultsH.AddColumn(t_("Parameter"));
	gridResultsH.AddColumn(t_("Value"));
	gridResultsH.Add(t_("Power flux [kW/m]"));
	
	opInfinite.WhenAction = [&]() {editDepth.Enable(!opInfinite);	OnCalc();};		
	editDepth.WhenAction = THISBACK(OnCalc);
	editDepth <<= 50;
	
	opTw = prevTw = 0;
	opTw.WhenAction = [&]() {
		double val;
		if (prevTw == 0)
			val = editTw;
		else if (prevTw == 1)
			val = 2*M_PI/editTw;
		else
			val = 1/editTw;
		
		if (opTw == 0)
			editTw <<= val;
		else if (opTw == 1)
			editTw <<= 2*M_PI/val;
		else
			editTw <<= 1/val;
			
		prevTw = opTw;
	};
	editTw <<= 12;
	editTw.WhenAction = THISBACK(OnCalc);
	editTw.Pattern("%.6g");
	editH <<= 2;
	editH.WhenAction = THISBACK(OnCalc);
	OnCalc();
	
	editDelta <<= 1;
	editDelta.WhenAction = THISBACK(OnDelta);
	
	return true;
}

void Tools_WaveInfo::OnCalc() {
	double g = Bem().g;
	double rho = Bem().rho;
	
	if ((!opInfinite && editDepth < 0) || editTw < 0) {
		for (int r = 0; r < gridResultsT.GetCount(); ++r)
			gridResultsT.Set(r, 1, Null);		
		for (int r = 0; r < gridResultsH.GetCount(); ++r)
			gridResultsH.Set(r, 1, Null);
		return;
	}
	double depth = opInfinite ? -1 : double(editDepth);
	double T;
	if (opTw == 0)
		T = editTw;
	else if (opTw == 1)
		T = 2*M_PI/editTw;
	else
		T = 1/editTw;
	
	gridResultsT.Set(0, 1, FDS(SeaWaves::Celerity(T, depth, g), 8));
	String cond;
	switch (SeaWaves::GetSeaType(T, depth, g)) {
	case SeaWaves::SHALLOW: 		cond = t_("Shallow");		break;
	case SeaWaves::INTERMEDIATE: 	cond = t_("Intermediate");	break;
	default:						cond = t_("Deep");
	}
	gridResultsT.Set(1, 1, cond);
	gridResultsT.Set(2, 1, FDS(SeaWaves::WaveNumber(T, depth, g, true), 8));
	gridResultsT.Set(3, 1, FDS(SeaWaves::WaveLength(T, depth, g), 8));
	gridResultsT.Set(4, 1, FDS(SeaWaves::GroupCelerity(T, depth, g), 8));	
	
	gridResultsH.Set(0, 1, FDS(SeaWaves::Power(T, editH, depth, g, rho), 8));
}

void Tools_WaveInfo::OnDelta() {
	if (editDelta < 0) 
		return;
	editDepth.SetInc(editDelta);	
	editTw.SetInc(editDelta);
	editH.SetInc(editDelta);
}

void Tools_WaveInfo::OnArrayBar(Bar &menu) {
	menu.Add(t_("Select all"), Null, THISBACK(ArraySelect)).Key(K_CTRL_A).Help(t_("Select all rows"));
	menu.Add(t_("Copy"), THISBACK(ArrayCopy)).Key(K_CTRL_C).Help(t_("Copy selected rows"));
}

void Tools_WaveInfo::ArrayCopy() {
	WaitCursor waitcursor;
	gridResultsT.SetClipboard(true, true); 
	gridResultsH.SetClipboard(true, true); 
}

void Tools_WaveInfo::ArraySelect() {
	WaitCursor waitcursor;
	gridResultsT.Select(0, gridResultsT.GetCount(), true);
	gridResultsH.Select(0, gridResultsH.GetCount(), true);
}


bool Tools_MooringInfo::Init() {
	CtrlLayout(*this);
	
	opStud = 0;
	opStud.WhenAction = [&] {OnCalc();};
	
	opForce = 2;
	opForce.WhenAction = [&] {OnCalc();};
	
	opType = 0;
	opType.WhenAction = [&] {OnCalc();};
	
	gridResults.MultiSelect().SetLineCy(EditField::GetStdHeight());
	gridResults.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, gridResults, true);};
	gridResults.AddColumn(t_("Parameter"), 35);
	gridResults.AddColumn(t_("Value"), 10);
	gridResults.AddColumn(t_("Description"), 40);
	gridResults.Add(t_("Outer diameter [m]"), 										"", t_("The diameter of a circular line with equivalent cross section"));
	gridResults.Add(t_("Outer contact diameter [m]"),								"", t_("Chain link width"));
	gridResults.Add(t_("Mass/length [kg/m]"), 										"", t_("Dry weight per unit of length"));	
	gridResults.Add(t_("Proof Load [N]"));
	gridResults.Add(t_("Minimum Breaking Load MBL [N]"));
	gridResults.Add(t_("Longitudinal/Axial/Tangential stiffness EA [N]"), 			"", t_("Product of Young's modulus (E) and cross-sectional area (A)"));
	gridResults.Add(t_("Longitudinal/Axial/Tangential drag coefficient"));
	gridResults.Add(t_("Longitudinal/Axial/Tangential drag coefficient CdAx"));
	gridResults.Add(t_("Longitudinal/Axial/Tangential drag diameter [m]"));
	gridResults.Add(t_("Longitudinal/Axial/Tangential added mass coefficient CaAx"));
	gridResults.Add(t_("Transverse/Normal drag coefficient"));
	gridResults.Add(t_("Transverse/Normal drag coefficient Cd"));
	gridResults.Add(t_("Transverse/Normal drag diameter [m]"));
	gridResults.Add(t_("Transverse/Normal added mass coefficient Ca"));
	
	gridGrade.SetLineCy(EditField::GetStdHeight());	
	gridGrade.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, gridGrade, true);};
	gridGrade.AddColumn(t_("Grade"));
	UVector<String> grade = {"API-2F", "R3", "R3S", "R4", "R4S", "R5"};
	for (const String &d : grade)
		gridGrade.Add(d);
	gridGrade.WhenSel =[&] {
		int row = gridGrade.GetCursor();
		if (row < 0)
			return;
		OnCalc();
	};
	gridGrade.SetCursor(3);	
	
	/*gridCatalogue.SetLineCy(EditField::GetStdHeight());	
	gridCatalogue.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, gridCatalogue, true);};
	gridCatalogue.AddColumn(t_("Diameter [mm]"));
	UVector<int> cat = {54, 56, 58, 60, 62, 64, 66, 68, 70, 73, 76, 78, 81, 84, 87, 90, 92, 95, 97, 100, 102, 105, 107, 111, 114, 117, 120, 122, 
						124, 127, 130, 132, 137, 142, 147, 152, 157, 162, 165, 168, 171, 175, 178, 180, 185, 188, 191, 194, 197, 200, 205, 220};
	for (int d : cat)
		gridCatalogue.Add(d);
	gridCatalogue.WhenSel =[&] {
		int row = gridCatalogue.GetCursor();
		if (row < 0)
			return;
		editDiam <<= gridCatalogue.Get(row, 0);
		OnCalc();
	};
	gridCatalogue.SetCursor(38);*/
	
	editDiam.WhenAction = [&] {OnCalc();};
	editDiam <<= 160;
	
	OnCalc();

	return true;
}

void Tools_MooringInfo::OnCalc() {
	enum Type {CHAIN, WIRE, POLYESTER, NYLON, HPME};
	
	opStud.Enable(opType == CHAIN);
	gridGrade.Enable(opType == CHAIN);
	
	double g = Bem().g;
	double rho = Bem().rho;
	
	double d = ~editDiam;			// Catalogue diameter
	if (IsNull(d))
		return;
	
	String uForce;
	double fcoeff;
	if (opForce == 0) {
		uForce = "N";
		fcoeff = 1;
	} else if (opForce == 1) {
		uForce = "kN";
		fcoeff = 1/1000.;
	} else {
		uForce = "t";
		fcoeff = 1/1000./9.8;
	}
	gridResults.Set(3, 0, Format(t_("Proof Load [%s]"), uForce));
	gridResults.Set(4, 0, Format(t_("Minimum Breaking Load MBL [%s]"), uForce));
	gridResults.Set(5, 0, Format(t_("Axial stiffness EA [%s]"), uForce));
	
	double od = Null, ocd = Null, ml = Null, pl = Null, mbl = Null, 
		   as = Null, adc = Null, adcm = Null, add = Null, aam = Null, 
		   			  ndc = Null, ndcm = Null, ndd = Null, nam = Null;
	String 			                                    smbl,
					  sadc, sadcm, 						   saam,
					  sndc, sndcm, 						   snam;
	if (opType == CHAIN) {
		if (opStud == 0)
			od = 1.8*d/1000;
		else
			od = 1.89*d/1000;
		if (opStud == 0)
			ocd = 3.35*d/1000;
		else
			ocd = 3.6*d/1000;
		if (opStud == 0)
			ml = 19.9*d*d/1000;
		else
			ml = 21.9*d*d/1000;
		int grow = gridGrade.GetCursor();
		if (grow >= 0) {
			String grade = gridGrade.Get(grow, 0);
			VectorMap<String, UVector<double>> gradeVals;
			// 				   PL Studless, PL Studlink and BL
			gradeVals.Add("API-2F", {0.0140, 0.0140, 0.0211});
			gradeVals.Add("R3",  	{0.0156, 0.0156, 0.0223});		
			gradeVals.Add("R3S", 	{0.0174, 0.018,  0.0249});
			gradeVals.Add("R4",  	{0.0192, 0.0216, 0.0274});
			gradeVals.Add("R4S", 	{0.0213, 0.024,  0.0304});
			gradeVals.Add("R5",  	{0.0223, 0.0251, 0.032});
			
			const UVector<double> &coeff = gradeVals.Get(grade);
			
			pl = fcoeff*coeff[opStud]*d*d*(44 - 0.08*d)*1000;	
			smbl = t_("DNV-OS-E302");	
			mbl = fcoeff*coeff[2]*d*d*(44 - 0.08*d)*1000;	
			
			if (opStud == 0)
				as = fcoeff*0.854*d*d*100000;
			else
				as = fcoeff*1.01*d*d*100000;
			
			sadc = t_("DNV-OS-E301 (OrcaFlex)");
			if (opStud == 0)			
				adc = 1.15;
			else
				adc = 1.4;	
			
			sadcm = t_("DNV-OS-E301 (MoorDyn)");
			adcm = adc*(d/1000)/od;
			add = d/M_PI/1000;	// d / PI
			
			saam = t_("BV NR493 R04 Jul 2021");
			aam = 0.5;
			
			sndc = t_("DNV-OS-E301 (OrcaFlex)");
			if (opStud == 0)
				ndc = 2.4;
			else
				ndc = 2.6;
			
			sndcm = t_("DNV-OS-E301 (MoorDyn)");
			ndcm = ndc*(d/1000)/od;
			ndd = d/1000;
			snam = t_("BV NR493 R04 Jul 2021");
			nam = 1;
		} 
	} else if (opType == WIRE) {
		od = ocd = 1.18*d/1000;
		ml = 5293*d*d/1000000;
		smbl = t_("The IEA Wind Task 49 Reference Floating Wind Array Design Basis");
		mbl = fcoeff*1022*d*d;
		as = fcoeff*97.1*d*d*1000;
		
		saam = t_("BV NR493 R04 Jul 2021");
		aam = 0;
		snam = t_("BV NR493 R04 Jul 2021");
		nam = 1;
	} else if (opType == POLYESTER) {
		od = ocd = 0.79*d/1000;
		ml = 679*d*d/1000000;
		smbl = t_("The IEA Wind Task 49 Reference Floating Wind Array Design Basis");
		mbl = fcoeff*308*d*d;
		as = fcoeff*(6795.572491625*d*d);		// (4359599.684*d - 600097037.135);

		saam = t_("BV NR493 R04 Jul 2021");
		aam = 0.15;
		snam = t_("BV NR493 R04 Jul 2021");
		nam = 1.1;
	} else if (opType == NYLON) {
		od = ocd = 0.81*d/1000;
		ml = 585*d*d/1000000;
		smbl = t_("The IEA Wind Task 49 Reference Floating Wind Array Design Basis");
		mbl = fcoeff*(207*d*d + 230*d*d*d/1000);
		saam = t_("BV NR493 R04 Jul 2021");
		aam = 0.15;
		snam = t_("BV NR493 R04 Jul 2021");
		nam = 1.1;
	} else if (opType == HPME) {
		od = ocd = 0.80*d/1000;
		ml = 496*d*d/1000000;
		smbl = t_("The IEA Wind Task 49 Reference Floating Wind Array Design Basis");
		mbl = fcoeff*(580*d*d + 651*d*d*d/1000);
		saam = t_("BV NR493 R04 Jul 2021");
		aam = 0.15;
		snam = t_("BV NR493 R04 Jul 2021");
		nam = 1.1;
	}
	int row = 0;	
	gridResults.Set(row++, 1, Nvl2(od,  FDS(od, 6),  S("-")));										// Outer diameter [m]: The diameter of a circular line with equivalent cross section
	gridResults.Set(row++, 1, Nvl2(ocd, FDS(ocd, 6), S("-")));										// Outer contact diameter [m]: Chain link width
	gridResults.Set(row++, 1, Nvl2(ml,  FDS(ml, 9),  S("-")));										// Mass/length [kg/m]: Dry weight per unit of length
	gridResults.Set(row++, 1, Nvl2(pl,  FDS(pl, 8),  S("-")));										// Proof Load [N]
	gridResults.Set(row  , 1, Nvl2(mbl, FDS(mbl, 8), S("-")));	gridResults.Set(row++, 2, smbl);	// Minimum Breaking Load MBL [N]
	gridResults.Set(row++, 1, Nvl2(as,  FDS(as, 8),  S("-")));										// Longitudinal/Axial/Tangential stiffness EA [N]: Product of Young's modulus (E) and cross-sectional area (A)
	gridResults.Set(row  , 1, Nvl2(adc, FDS(adc, 4), S("-")));	gridResults.Set(row++, 2, sadc);	// Longitudinal/Axial/Tangential drag coefficient
	gridResults.Set(row  , 1, Nvl2(adcm,FDS(adcm,6), S("-")));	gridResults.Set(row++, 2, sadcm);	// Longitudinal/Axial/Tangential drag coefficient CdAx
	gridResults.Set(row++, 1, Nvl2(add, FDS(add, 8), S("-")));										// Longitudinal/Axial/Tangential drag diameter [m]
	gridResults.Set(row  , 1, Nvl2(aam, FDS(aam, 4), S("-")));	gridResults.Set(row++, 2, saam);	// Longitudinal/Axial/Tangential added mass coefficient CaAx
	gridResults.Set(row  , 1, Nvl2(ndc, FDS(ndc, 4), S("-")));	gridResults.Set(row++, 2, sndc);	// Transverse/Normal drag coefficient 	
	gridResults.Set(row  , 1, Nvl2(ndcm,FDS(ndcm,6), S("-")));	gridResults.Set(row++, 2, sndcm);	// Transverse/Normal drag coefficient Cd
	gridResults.Set(row++, 1, Nvl2(ndd, FDS(ndd, 4), S("-")));										// Transverse/Normal drag diameter [m]
	gridResults.Set(row  , 1, Nvl2(nam, FDS(nam, 4), S("-")));	gridResults.Set(row++, 2, snam);	// Transverse/Normal added mass coefficient Ca
}

#ifdef PLATFORM_WIN32

void Tools_OrcaLicense::Init() {
	CtrlLayout(*this);
	
	butCapture.WhenAction = [&]() {
		int aval = orca.IsAvailable();
		
		if (aval == Orca::NOT_INSTALLED) {
			PromptOK(t_("OrcaWave is not installed"));
			return;
		} else if (aval ==	Orca::NOT_AVAILABLE) 
			PromptOK(t_("OrcaWave is not available"));
		
		String label = butCapture.GetLabel();
		
		if (label == t_("Check license")) {
			butCapture.SetLabel(t_("Stop checking"));
			labCapture.SetText(t_("Checking if an OrcaWave license is available ..."));
			SetTimeCallback(int(-60*1000), [&]() {CheckAvailable();}, 12);
			CheckAvailable();
		} else if (label == t_("Stop checking")) {
			butCapture.SetLabel(t_("Check license"));
			labCapture.SetText("");
			KillTimeCallback(12);
		}
	};
}

void Tools_OrcaLicense::CheckAvailable() {
	labCapture.SetText(t_("Checking now ..."));
	Ctrl::ProcessEvents();
	if (orca.IsAvailable() == Orca::AVAILABLE) {
		KillTimeCallback(12);
		butCapture.SetLabel(t_("Check license"));
		labCapture.SetText(Format(t_("OrcaWave is now available (%`)"), GetSysTime()));
		PromptOK(t_("OrcaWave is now available"));
	} else // Orca::NOT_AVAILABLE) 
		labCapture.SetText(t_("Checking if an OrcaWave license is available ..."));
}


#endif