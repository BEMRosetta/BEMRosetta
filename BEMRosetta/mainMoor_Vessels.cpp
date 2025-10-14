// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2025, the BEMRosetta author and contributors
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


void MainMoor_Vessels::Init(Mooring &mooring) {
	CtrlLayout(*this);
	
	parray = &array;
	pmooring = &mooring;
	
	array.WhenSel = THISBACK(ArrayOnCursor);

	butAdd <<= THISBACK(ArrayOnAdd);
	butDuplicate <<= THISBACK(ArrayOnDuplicate);
	butRemove <<= THISBACK(ArrayOnRemove);
	
	Load();
	
	edName.WhenAction	= [&] {ArrayUpdateCursor();};
}
	
void MainMoor_Vessels::InitArray() {
	array.Reset();
	array.SetLineCy(EditField::GetStdHeight()).MultiSelect();
	array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, array, true, true, [&]{ArrayUpdateCursor();});};
	
	array.AddColumn(t_("Name"), 40);
}

bool MainMoor_Vessels::ArrayUpdateCursor() {
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
	
	GetParentCtrl<MainMoor>(this).OnUpdate();
	
	return true;
}

void MainMoor_Vessels::ArrayOnCursor() {
	int id = array.GetCursor();
	if (id < 0)
		return;
	
	edName <<= array.Get(id, 0);
}

void MainMoor_Vessels::ArrayClear() {
	edName <<= Null;
}

void MainMoor_Vessels::Load() {
	InitArray();
	const auto &mooring = *pmooring;
	
	for (int id = 0; id < mooring.vessels.size(); ++id) {
		const String &con = mooring.vessels[id].name;
		array.Add(con);
	}
	if (mooring.connections.size() > 0)
		array.SetCursor(0);
}

void MainMoor_Vessels::Save() {
	auto &mooring = *pmooring;
	
	mooring.vessels.SetCount(array.GetCount());
	for (int id = 0; id < mooring.vessels.size(); ++id) 
		mooring.vessels[id].name = array.Get(id, 0);
}