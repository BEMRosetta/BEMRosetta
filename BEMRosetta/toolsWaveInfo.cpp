// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2024, the BEMRosetta author and contributors
#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas.h>
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
	tab.Add(waveInfo.SizePos(), t_("Waves"));
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
	gridResultsT.Set(2, 1, FDS(SeaWaves::WaveNumber(T, depth, true), 8));
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


#ifdef PLATFORM_WIN32

void Tools_OrcaLicense::Init() {
	CtrlLayout(*this);
	
	butCapture.WhenAction = [&]() {
		if (!orca.IsLoaded()) {
			if (!orca.FindInit()) 
				PromptOK(t_("OrcaWave is not available"));
		}
		
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
	if (orca.IsAvailable()) {
		KillTimeCallback(12);
		butCapture.SetLabel(t_("Check license"));
		labCapture.SetText(Format(t_("OrcaWave is now available (%`)"), GetSysTime()));
		PromptOK(t_("OrcaWave is now available"));
	} else
		labCapture.SetText(t_("Checking if an OrcaWave license is available ..."));
}


#endif