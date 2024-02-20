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
	
	editDepth <<= 50;
	editDepth.WhenAction = THISBACK(OnCalc);
	editT <<= 12;
	editT.WhenAction = THISBACK(OnCalc);
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
	
	if (editDepth < 0 || editT < 0) {
		for (int r = 0; r < gridResultsT.GetCount(); ++r)
			gridResultsT.Set(r, 1, Null);		
		for (int r = 0; r < gridResultsH.GetCount(); ++r)
			gridResultsH.Set(r, 1, Null);
		return;
	}
	gridResultsT.Set(0, 1, FDS(SeaWaves::Celerity(editT, editDepth, g), 8));
	String cond;
	switch (SeaWaves::GetSeaType(editT, editDepth, g)) {
	case SeaWaves::SHALLOW: 		cond = t_("Shallow");		break;
	case SeaWaves::INTERMEDIATE: 	cond = t_("Intermediate");	break;
	default:						cond = t_("Deep");
	}
	gridResultsT.Set(1, 1, cond);
	gridResultsT.Set(2, 1, FDS(SeaWaves::WaveNumber(editT, editDepth, true), 8));
	gridResultsT.Set(3, 1, FDS(SeaWaves::WaveLength(editT, editDepth, g), 8));
	gridResultsT.Set(4, 1, FDS(SeaWaves::GroupCelerity(editT, editDepth, g), 8));	
	
	gridResultsH.Set(0, 1, FDS(SeaWaves::Power(editT, editH, editDepth, g, rho), 8));
}

void Tools_WaveInfo::OnDelta() {
	if (editDelta < 0) 
		return;
	editDepth.SetInc(editDelta);	
	editT.SetInc(editDelta);
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
