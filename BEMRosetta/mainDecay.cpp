// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2021, the BEMRosetta author and contributors
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


void MainDecay::Init() {
	CtrlLayout(*this);	
	
	files.Init();
	staticfiles.Add(files.SizePos());

	const String moorFiles = ".json";
	String moorFilesAst = clone(moorFiles);
	moorFilesAst.Replace(".", "*.");
	file.Type(Format("All supported mooring files (%s)", moorFiles), moorFilesAst);
	file.AllFilesType();
	file.BrowseRightWidth(40).UseOpenFolder(true).BrowseOpenFolderWidth(10);
	butLoad.Tip(t_("Loads mooring file")).WhenAction = [&] {OnLoad();};
	butSave.Tip(t_("Saves mooring file")).WhenAction = [&] {OnSave();};
}


bool MainDecay::OnLoad() {
	try {

	} catch (const Exc &e) {
		Exclamation(DeQtfLf(e));
		return false;
	}
	
	return true;
}

bool MainDecay::OnSave() {
	try {

	} catch (const Exc &e) {
		Exclamation(DeQtfLf(e));
		return false;
	}
	
	return true;
}

void MainDecay::OnUpdate() {
	Exclamation("Pendiente");	
}
