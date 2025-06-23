// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <SurfaceCanvas/SurfaceCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <TabBar/TabBar.h>
#include <DropGrid/DropGrid.h>

#include <BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;

#include "main.h"


void MainMoor_LinesTypes::Init(Mooring &mooring) {
	CtrlLayout(*this);
	
	parray = &array;
	pmooring = &mooring;
	
	array.WhenSel = THISBACK(ArrayOnCursor);

	butAdd <<= THISBACK(ArrayOnAdd);
	butDuplicate <<= THISBACK(ArrayOnDuplicate);
	butRemove <<= THISBACK(ArrayOnRemove);
	
	InitArray();
	
	edName.WhenAction 		= [&] {ArrayUpdateCursor();};
	edMass.WhenAction 		= [&] {ArrayUpdateCursor();};
	edBL.WhenAction 		= [&] {ArrayUpdateCursor();};
	edDiameter.WhenAction	= [&] {ArrayUpdateCursor();};
}
	
void MainMoor_LinesTypes::InitArray() {
	array.Reset();
	array.SetLineCy(EditField::GetStdHeight()).MultiSelect();
	array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, array);};
	
	array.AddColumn(t_("Name"), 40);
	array.AddColumn(t_("Mass"), 40);
	array.AddColumn(t_("Diameter"), 40);
	array.AddColumn(t_("BL"), 40);
}

bool MainMoor_LinesTypes::ArrayUpdateCursor() {
	int id = array.GetCursor();
	if (id < 0) {
		if (array.GetCount() == 0) {
			InitArray();
			array.Add();
			id = 0;
		} else
			id = array.GetCount()-1;
	}	
	
	array.Set(id, 0, ~edName);
	array.Set(id, 1, ~edMass);
	array.Set(id, 2, ~edDiameter);
	array.Set(id, 3, ~edBL);
	
	return true;
}

void MainMoor_LinesTypes::ArrayOnCursor() {
	int id = array.GetCursor();
	if (id < 0)
		return;
	
	edName <<= array.Get(id, 0);
	edMass <<= array.Get(id, 1);
	edDiameter <<= array.Get(id, 2);
	edBL <<= array.Get(id, 3);
}

void MainMoor_LinesTypes::ArrayClear() {
	edName <<= Null;
	edMass <<= Null;
	edDiameter <<= Null;
	edBL <<= Null;
}

void MainMoor_LinesTypes::Load() {
	InitArray();
	const auto &mooring = *pmooring;
	
	for (int id = 0; id < mooring.lineTypes.size(); ++id) {
		const auto &val = mooring.lineTypes[id];
		array.Add(val.name, val.mass, val.diameter, val.bl);
	}
	if (mooring.lineTypes.size() > 0)
		array.SetCursor(0);
}

void MainMoor_LinesTypes::Save() {
	auto &mooring = *pmooring;
	
	mooring.lineTypes.SetCount(array.GetCount());
	for (int id = 0; id < mooring.lineTypes.size(); ++id) {
		mooring.lineTypes[id].name = array.Get(id, 0);
		mooring.lineTypes[id].mass = array.Get(id, 1);
		mooring.lineTypes[id].diameter = array.Get(id, 2);
		mooring.lineTypes[id].bl = array.Get(id, 3);
	}
}	
