// SPDX-License-Identifier: GPL-3.0-or-later
#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <TabBar/TabBar.h>

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
	edDiameter.WhenAction	= [&] {ArrayUpdateCursor();};
}
	
void MainMoor_LinesTypes::InitArray() {
	array.Reset();
	array.SetLineCy(EditField::GetStdHeight());
	array.AddColumn(t_("Name"), 40);
	array.AddColumn(t_("Mass"), 40);
	array.AddColumn(t_("Diameter"), 40);
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
	
	return true;
}

void MainMoor_LinesTypes::ArrayOnCursor() {
	int id = array.GetCursor();
	if (id < 0)
		return;
	
	edName <<= array.Get(id, 0);
	edMass <<= array.Get(id, 1);
	edDiameter <<= array.Get(id, 2);
}

void MainMoor_LinesTypes::ArrayClear() {
	edName <<= Null;
	edMass <<= Null;
	edDiameter <<= Null;
}

void MainMoor_LinesTypes::Load() {
	InitArray();
	const auto &mooring = *pmooring;
	
	for (int id = 0; id < mooring.lineTypes.size(); ++id) {
		const auto &val = mooring.lineTypes[id];
		array.Add(val.name, val.mass, val.diameter);
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
	}
}	
