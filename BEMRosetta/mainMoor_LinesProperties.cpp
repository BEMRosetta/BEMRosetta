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


void MainMoor_LineProperties::Init(Mooring &mooring) {
	CtrlLayout(*this);
	
	parray = &array;
	pmooring = &mooring;
	
	array.WhenSel = THISBACK(ArrayOnCursor);

	butAdd <<= THISBACK(ArrayOnAdd);
	butDuplicate <<= THISBACK(ArrayOnDuplicate);
	butRemove <<= THISBACK(ArrayOnRemove);
	
	InitArray();
	
	edName.WhenAction 		= [&] {ArrayUpdateCursor();};
	dropLineType.WhenAction	= [&] {ArrayUpdateCursor();};
	edLength.WhenAction		= [&] {ArrayUpdateCursor();};
	dropFrom.WhenAction		= [&] {ArrayUpdateCursor();};
	dropTo.WhenAction		= [&] {ArrayUpdateCursor();};
}
	
void MainMoor_LineProperties::InitArray() {
	array.Reset();
	array.SetLineCy(EditField::GetStdHeight());
	array.AddColumn(t_("Name"), 40);
	array.AddColumn(t_("Line type"), 40);
	array.AddColumn(t_("Length"), 40);
	array.AddColumn(t_("From"), 40);
	array.AddColumn(t_("To"), 40);
}

bool MainMoor_LineProperties::ArrayUpdateCursor() {
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
	array.Set(id, 1, ~dropLineType);
	array.Set(id, 2, ~edLength);
	array.Set(id, 3, ~dropFrom);
	array.Set(id, 4, ~dropTo);

	auto &mooring = *pmooring;
	Mooring temp = clone(mooring);
	
	Save();
	String error = mooring.Test();
	if (!error.IsEmpty()) {
		Exclamation(error);
		mooring = clone(temp);
		Load();
		return false;
	}
	
	return true;
}

void MainMoor_LineProperties::ArrayOnCursor() {
	int id = array.GetCursor();
	if (id < 0)
		return;
	
	edName 		 <<= array.Get(id, 0);
	dropLineType <<= array.Get(id, 1);
	edLength 	 <<= array.Get(id, 2);
	dropFrom 	 <<= array.Get(id, 3);
	dropTo		 <<= array.Get(id, 4);
}

void MainMoor_LineProperties::ArrayClear() {
	edName <<= Null;
	dropLineType <<= Null;
	edLength <<= Null;
	dropFrom <<= Null;
	dropTo <<= Null;
}

void MainMoor_LineProperties::Load() {
	InitArray();
	const auto &mooring = *pmooring;
	
	for (int id = 0; id < mooring.lineProperties.size(); ++id) {
		const auto &val = mooring.lineProperties[id];
		array.Add(val.name, val.nameType, val.length, val.from, val.to);
	}
	LoadDrop();
	if (mooring.lineProperties.size() > 0)
		array.SetCursor(0);
}

void MainMoor_LineProperties::Save() {
	auto &mooring = *pmooring;
	
	mooring.lineProperties.SetCount(array.GetCount());
	for (int id = 0; id < mooring.lineProperties.size(); ++id) {
		mooring.lineProperties[id].name 	= array.Get(id, 0);
		mooring.lineProperties[id].nameType = array.Get(id, 1);
		mooring.lineProperties[id].length 	= array.Get(id, 2);
		mooring.lineProperties[id].from 	= array.Get(id, 3);
		mooring.lineProperties[id].to 		= array.Get(id, 4);
	}
}

void MainMoor_LineProperties::LoadDrop() {
	auto &mooring = *pmooring;
	
	String strFrom = ~dropFrom;
	String strTo   = ~dropTo;
	dropFrom.Clear();
	dropTo.Clear();
	for (int i = 0; i < mooring.connections.size(); ++i) {
		dropFrom.Add(mooring.connections[i].name);
		dropTo.  Add(mooring.connections[i].name);
	}
	dropFrom <<= strFrom;
	dropTo   <<= strTo;
	
	String strLineType = ~dropLineType;
	dropLineType.Clear();
	for (int i = 0; i < mooring.lineTypes.size(); ++i) 
		dropLineType.Add(mooring.lineTypes[i].name);
	dropLineType <<= strLineType;
	
}
