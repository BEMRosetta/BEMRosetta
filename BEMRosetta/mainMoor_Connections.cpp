// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
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


void MainMoor_Connections::Init(Mooring &mooring) {
	CtrlLayout(*this);
	
	parray = &array;
	pmooring = &mooring;
	
	array.WhenSel = THISBACK(ArrayOnCursor);

	butAdd <<= THISBACK(ArrayOnAdd);
	butDuplicate <<= THISBACK(ArrayOnDuplicate);
	butRemove <<= THISBACK(ArrayOnRemove);
	
	InitArray();
	
	edName.WhenAction	= [&] {ArrayUpdateCursor();};
	opVessel.WhenAction	= [&] {ArrayUpdateCursor();};
	edx.WhenAction 		= [&] {ArrayUpdateCursor();};
	edy.WhenAction 		= [&] {ArrayUpdateCursor();};
	edz.WhenAction 		= [&] {ArrayUpdateCursor();};
}
	
void MainMoor_Connections::InitArray() {
	array.Reset();
	array.SetLineCy(EditField::GetStdHeight());
	array.AddColumn(t_("Name"), 40);
	array.AddColumn(t_("Vessel"), 30).With([](One<Ctrl>& x) {x.Create<Option>().SetReadOnly();});
	array.AddColumn(t_("x"), 40);
	array.AddColumn(t_("y"), 40);
	array.AddColumn(t_("z"), 40);
	
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
	array.Set(id, 1, ~opVessel);
	array.Set(id, 2, ~edx);
	array.Set(id, 3, ~edy);
	array.Set(id, 4, ~edz);
	
	return true;
}

void MainMoor_Connections::ArrayOnCursor() {
	int id = array.GetCursor();
	if (id < 0)
		return;
	
	edName <<= array.Get(id, 0);
	opVessel <<= array.Get(id, 1);
	edx <<= array.Get(id, 2);
	edy <<= array.Get(id, 3);
	edz <<= array.Get(id, 4);
}

void MainMoor_Connections::ArrayClear() {
	edName <<= Null;
	opVessel <<= Null;
	edx <<= Null;
	edy <<= Null;
	edz <<= Null;
}

void MainMoor_Connections::Load() {
	InitArray();
	const auto &mooring = *pmooring;
	
	for (int id = 0; id < mooring.connections.size(); ++id) {
		const auto &con = mooring.connections[id];
		array.Add(con.name, con.vessel, con.x, con.y, con.z);
	}
	if (mooring.connections.size() > 0)
		array.SetCursor(0);
}

void MainMoor_Connections::Save() {
	auto &mooring = *pmooring;
	
	mooring.connections.SetCount(array.GetCount());
	for (int id = 0; id < mooring.connections.size(); ++id) {
		mooring.connections[id].name = array.Get(id, 0);
		mooring.connections[id].vessel = array.Get(id, 1);
		mooring.connections[id].x = array.Get(id, 2);
		mooring.connections[id].y = array.Get(id, 3);
		mooring.connections[id].z = array.Get(id, 4);
	}
}