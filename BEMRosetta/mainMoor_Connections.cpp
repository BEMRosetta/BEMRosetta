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


void MainMoor_Connections::Init(Mooring &mooring) {
	CtrlLayout(*this);
	
	parray = &array;
	pmooring = &mooring;
	
	array.WhenSel = THISBACK(ArrayOnCursor);

	butAdd <<= THISBACK(ArrayOnAdd);
	butDuplicate <<= THISBACK(ArrayOnDuplicate);
	butRemove.WhenAction = [&]{ArrayOnRemove();	GetDefinedParent<MainMoor>(this).OnUpdate();};
	
	InitArray();
	
	edName.WhenAction	= [&] {ArrayUpdateCursor();};
	dropVessel.WhenAction=[&] {ArrayUpdateCursor();};
	edx.WhenAction 		= [&] {ArrayUpdateCursor();};
	edy.WhenAction 		= [&] {ArrayUpdateCursor();};
	edz.WhenAction 		= [&] {ArrayUpdateCursor();};
}
	
void MainMoor_Connections::InitArray() {
	array.Reset();
	array.SetLineCy(EditField::GetStdHeight()).MultiSelect();
	array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, array, true, true, [&]{ArrayUpdateCursor();});};
	
	array.AddColumn(t_("Name"), 40);
	array.AddColumn(t_("Vessel"), 30);
	array.AddColumn(t_("x [m]"), 40);
	array.AddColumn(t_("y [m]"), 40);
	array.AddColumn(t_("z [m]"), 40);
	
}

bool MainMoor_Connections::ArrayUpdateCursor() {
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
	array.Set(id, 1, ~dropVessel);
	array.Set(id, 2, ~edx);
	array.Set(id, 3, ~edy);
	array.Set(id, 4, ~edz);
	
	GetDefinedParent<MainMoor>(this).OnUpdate();
	
	return true;
}

void MainMoor_Connections::ArrayOnCursor() {
	int id = array.GetCursor();
	if (id < 0)
		return;
	
	edName <<= array.Get(id, 0);
	dropVessel <<= array.Get(id, 1);
	edx <<= ScanDouble(array.Get(id, 2).ToString());
	edy <<= ScanDouble(array.Get(id, 3).ToString());
	edz <<= ScanDouble(array.Get(id, 4).ToString());
}

void MainMoor_Connections::ArrayClear() {
	edName <<= Null;
	dropVessel <<= Null;
	edx <<= Null;
	edy <<= Null;
	edz <<= Null;
}

void MainMoor_Connections::Load() {
	InitArray();
	const auto &mooring = *pmooring;
	
	for (int id = 0; id < mooring.connections.size(); ++id) {
		const auto &con = mooring.connections[id];
		array.Add(con.name, con.where, con.x, con.y, con.z);
	}
	LoadDrop();
	if (mooring.connections.size() > 0)
		array.SetCursor(0);
}

void MainMoor_Connections::Save() {
	auto &mooring = *pmooring;
	
	mooring.connections.SetCount(array.GetCount());
	for (int id = 0; id < mooring.connections.size(); ++id) {
		mooring.connections[id].name = array.Get(id, 0);
		mooring.connections[id].where = array.Get(id, 1);
		mooring.connections[id].x = ScanDouble(array.Get(id, 2).ToString());
		mooring.connections[id].y = ScanDouble(array.Get(id, 3).ToString());
		mooring.connections[id].z = ScanDouble(array.Get(id, 4).ToString());
	}
}

void MainMoor_Connections::LoadDrop() {
	auto &mooring = *pmooring;
	
	String strVessel = ~dropVessel;
	dropVessel.Clear();
	for (int i = 0; i < mooring.vessels.size(); ++i) 
		dropVessel.Add(mooring.vessels[i].name);
	dropVessel.Add("fixed");
	dropVessel <<= strVessel;
}