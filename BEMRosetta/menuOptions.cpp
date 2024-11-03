// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2023, the BEMRosetta author and contributors
#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <TabBar/TabBar.h>
#include <DropGrid/DropGrid.h>

#include <BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;

#include "main.h"

using namespace Eigen;

void MenuOptions::Init(BEM &_bem) {
	CtrlLayout(*this);
	
	bem = &_bem;
	butSave <<= butSave2 <<= THISBACK(OnSave);
	
	g.isMouseEnter = rho.isMouseEnter = dofType.isMouseEnter = headingType.isMouseEnter = false;
	
	for (int i = 0; BasicBEM::strDOFType[i][0] != '\0'; ++i)
		dofType.Add(BasicBEM::strDOFType[i]);
	for (int i = 0; BasicBEM::strHeadingType[i][0] != '\0'; ++i)
		headingType.Add(BasicBEM::strHeadingType[i]);
}

void MenuOptions::InitSerialize(bool ret, bool &openOptions) {
	;
}

void MenuOptions::Load() {
	g <<= bem->g;
	rho <<= bem->rho;
	len <<= bem->len;
	depth <<= bem->depth;
	//discardNegDOF <<= bem->discardNegDOF;
	//thres <<= bem->thres;
	calcAinf <<= bem->calcAinf;
	calcAinf_w <<= bem->calcAinf_w;
	maxTimeA <<= bem->maxTimeA;
	numValsA <<= bem->numValsA;	
	onlyDiagonal <<= bem->onlyDiagonal;
	nemoh3Path <<= bem->nemoh3Path;
	nemoh115Path <<= bem->nemoh115Path;
	nemohPath <<= bem->nemohPath;
	nemohPathGREN <<= bem->nemohPathGREN;
	foammPath <<= bem->foammPath;
	hamsPath <<= bem->hamsPath;
	hamsBodyPath <<= bem->hamsBodyPath;
	aqwaPath <<= bem->aqwaPath;
	wamitPath <<= bem->wamitPath;
	volWarning <<= bem->volWarning;
	volError <<= bem->volError;
	roundVal <<= bem->roundVal;
	roundEps <<= bem->roundEps;
	csvSeparator <<= bem->csvSeparator;
	legend_w_units <<= bem->legend_w_units;
	legend_w_solver <<= bem->legend_w_solver;
	pythonEnv <<= bem->pythonEnv;
	zeroIfEmpty <<= bem->zeroIfEmpty;
	
	dofType.SetIndex(bem->dofType);
	headingType.SetIndex(bem->headingType);
}

void MenuOptions::OnSave() {
	String warning;
	
	if (Trim(~nemohPath) != "" && !DirectoryExists(~nemohPath))
		warning << Format("'%s' file doesn't exist\n", ~nemohPath);
	if (Trim(~nemoh115Path) != "" && !DirectoryExists(~nemoh115Path))
		warning << Format("'%s' file doesn't exist\n", ~nemoh115Path);
	if (Trim(~nemoh3Path) != "" && !DirectoryExists(~nemoh3Path))
		warning << Format("'%s' file doesn't exist\n", ~nemoh3Path);
	if (Trim(~nemohPathGREN) != "" && !FileExists(~nemohPathGREN))
		warning << Format("'%s' file doesn't exist\n", ~nemohPathGREN);
	
	if (warning.IsEmpty() || 
		ErrorOKCancel(DeQtfLf(Format(t_("Some errors found:\n%s Do you wish to save them?"), warning)))) {
		bem->g = ~g;
		bem->rho = ~rho;
		bem->len = ~len;
		bem->depth = ~depth;
		//bem->discardNegDOF = ~discardNegDOF;
		//bem->thres = ~thres;
		bem->calcAinf = ~calcAinf;
		bem->calcAinf_w = ~calcAinf_w;
		bem->maxTimeA = ~maxTimeA;
		bem->numValsA = ~numValsA;	
		bem->onlyDiagonal = ~onlyDiagonal;
		bem->nemohPath = ~nemohPath;
		bem->nemoh115Path = ~nemoh115Path;
		bem->nemoh3Path = ~nemoh3Path;	
		bem->nemohPathGREN = ~nemohPathGREN;
		bem->foammPath = ~foammPath;
		bem->hamsPath = ~hamsPath;
		bem->hamsBodyPath = ~hamsBodyPath;
		bem->aqwaPath = ~aqwaPath;
		bem->wamitPath = ~wamitPath;
		bem->volWarning = ~volWarning;
		bem->volError = ~volError;
		bem->roundVal = ~roundVal;
		bem->roundEps = ~roundEps;
		bem->csvSeparator = ~csvSeparator;
		ScatterDraw::SetDefaultCSVSeparator(~csvSeparator);
		bem->legend_w_units = ~legend_w_units;
		bem->legend_w_solver = ~legend_w_solver;
		bem->pythonEnv = ~pythonEnv;
		bem->zeroIfEmpty = ~zeroIfEmpty;
		
		bem->dofType = BasicBEM::DOFType(dofType.GetIndex());
		bem->headingType = BasicBEM::HeadingType(headingType.GetIndex());
		bem->UpdateHeadAll();
		bem->UpdateHeadAllMD();
		Ma().OptionsUpdated(rho, g, bem->dofType, bem->headingType);
	}
	Ma().SetLastTab();
}

bool MenuOptions::IsChanged() {
	if (!EqualDecimals(bem->g, ~g, 8)) 
		return true;
	if (!EqualDecimals(bem->rho, ~rho, 8))
		return true;
	if (!EqualDecimals(bem->len, ~len, 8))
		return true;
	if (!EqualDecimals(bem->depth, ~depth, 8))
		return true;
	//if (bem->discardNegDOF != ~discardNegDOF)
	//	return true;
	//if (!EqualDecimals(bem->thres, ~thres, 8)) 
	//	return true;
	if (bem->calcAinf != ~calcAinf)
		return true;
	if (bem->calcAinf_w != ~calcAinf_w)
		return true;
	if (!EqualDecimals(bem->maxTimeA, ~maxTimeA, 8))
		return true;
	if (bem->numValsA != ~numValsA)
		return true;
	if (bem->onlyDiagonal != ~onlyDiagonal)
		return true;
	if (bem->nemohPath != ~nemohPath)
		return true;
	if (bem->nemoh115Path != ~nemoh115Path)
		return true;
	if (bem->nemoh3Path != ~nemoh3Path)
		return true;
	if (bem->nemohPathGREN != ~nemohPathGREN)
		return true;
	if (bem->foammPath != ~foammPath)
		return true;
	if (bem->hamsPath != ~hamsPath)
		return true;
	if (bem->hamsBodyPath != ~hamsBodyPath)
		return true;
	if (bem->aqwaPath != ~aqwaPath)
		return true;
	if (bem->wamitPath != ~wamitPath)
		return true;
	if (bem->volWarning != ~volWarning)
		return true;
	if (bem->volError != ~volError)
		return true;
	if (bem->roundVal != ~roundVal)
		return true;
	if (bem->roundEps != ~roundEps)
		return true;
	if (bem->dofType != dofType.GetIndex())
		return true;
	if (bem->headingType != headingType.GetIndex())
		return true;
	if (bem->csvSeparator != ~csvSeparator)
		return true;
	if (bem->legend_w_solver != ~legend_w_solver)
		return true;
	if (bem->legend_w_units != ~legend_w_units)
		return true;
	if (bem->pythonEnv != ~pythonEnv)
		return true;
	if (bem->zeroIfEmpty != ~zeroIfEmpty)
		return true;
				
	return false;
}