// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U/Controls4U.h>
#include <Scatter/ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas/GLCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <TabBar/TabBar.h>
#include <DropGrid/DropGrid.h>

#include <BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;

#include "main.h"


void MainDecay_Files::Init() {
	CtrlLayout(*this);
	
	parray = &array;
	
	array.WhenSel = THISBACK(ArrayOnCursor);

	butAdd <<= THISBACK(ArrayOnAdd);
	butDuplicate <<= THISBACK(ArrayOnDuplicate);
	butRemove <<= THISBACK(ArrayOnRemove);
	
	InitArray();
	
	edFile.WhenAction 	= [&] {ArrayUpdateCursor();};
	edrow.WhenAction 	= [&] {ArrayUpdateCursor();};
	edcolt.WhenAction	= [&] {ArrayUpdateCursor();};
	edcolz.WhenAction	= [&] {ArrayUpdateCursor();};
}
	
void MainDecay_Files::InitArray() {
	array.Reset();
	array.SetLineCy(EditField::GetStdHeight());
	array.AddColumn(t_("File"), 40);
	array.AddColumn(t_("#Row"), 10);
	array.AddColumn(t_("#Col t"), 10);
	array.AddColumn(t_("#Col z"), 10);
}

bool MainDecay_Files::ArrayUpdateCursor() {
	int id = array.GetCursor();
	if (id < 0) {
		if (array.GetCount() == 0) {
			InitArray();
			array.Add();
			id = 0;
		} else
			id = array.GetCount()-1;
	}	
	
	array.Set(id, 0, ~edFile);
	array.Set(id, 1, ~edrow);
	array.Set(id, 2, ~edcolt);
	array.Set(id, 3, ~edcolz);
	
	return true;
}

void MainDecay_Files::ArrayOnCursor() {
	int id = array.GetCursor();
	if (id < 0)
		return;
	
	edFile <<= array.Get(id, 0);
	edrow  <<= array.Get(id, 1);
	edcolt <<= array.Get(id, 2);
	edcolz <<= array.Get(id, 3);
}

void MainDecay_Files::ArrayClear() {
	edFile <<= Null;
	edrow  <<= Null;
	edcolt <<= Null;
	edcolz <<= Null;
}

void MainDecay_Files::Load() {
	InitArray();
/*	const auto &mooring = *pmooring;
	
	for (int id = 0; id < mooring.lineTypes.size(); ++id) {
		const auto &val = mooring.lineTypes[id];
		array.Add(val.name, val.mass, val.diameter);
	}
	if (mooring.lineTypes.size() > 0)
		array.SetCursor(0);*/
}

void MainDecay_Files::Save() {
/*	auto &mooring = *pmooring;
	
	mooring.lineTypes.SetCount(array.GetCount());
	for (int id = 0; id < mooring.lineTypes.size(); ++id) {
		mooring.lineTypes[id].name = array.Get(id, 0);
		mooring.lineTypes[id].mass = array.Get(id, 1);
		mooring.lineTypes[id].diameter = array.Get(id, 2);
	}*/
}	
