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


void MainMoor_LineProperties::Init(Mooring &mooring) {
	CtrlLayout(*this);
	
	parray = &array;
	pmooring = &mooring;
	
	array.WhenSel = THISBACK(ArrayOnCursor);

	butAdd <<= THISBACK(ArrayOnAdd);
	butDuplicate <<= THISBACK(ArrayOnDuplicate);
	butRemove.WhenAction = [&]{ArrayOnRemove();	GetParentCtrl<MainMoor>(this).OnUpdate();};
	
	InitArray();
	
	edName.WhenAction 		= [&] {ArrayUpdateCursor();};
	dropLineType.WhenAction	= [&] {ArrayUpdateCursor();};
	edLength.WhenAction		= [&] {ArrayUpdateCursor();};
	edNumSeg.WhenAction		= [&] {ArrayUpdateCursor();};
	dropFrom.WhenAction		= [&] {ArrayUpdateCursor();};
	dropTo.WhenAction		= [&] {ArrayUpdateCursor();};
	butTaut.WhenAction 		= [&] {
		String sfrom = ~dropFrom;
		Point3D pfrom = Null;
		for (const Mooring::Connection &c : mooring.connections)
			if (c.name == sfrom) {
				pfrom = Point3D(c.x, c.y, c.z);
				break;
			}
		if (IsNull(pfrom))
			return;
		String sto = ~dropTo;
		Point3D pto = Null;
		for (const Mooring::Connection &c : mooring.connections)
			if (c.name == sto) {
				pto = Point3D(c.x, c.y, c.z);
				break;
			}
		if (IsNull(pto))
			return;

		double elong = ~GetParentCtrl<MainMoor>(this).right.edElong;
		
		edLength <<= Distance(pfrom, pto)*(100 + elong)/100.;
		ArrayUpdateCursor();
	};
}
	
void MainMoor_LineProperties::InitArray() {
	array.Reset();
	array.SetLineCy(EditField::GetStdHeight()).MultiSelect();
	array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, array, true, true, [&]{ArrayUpdateCursor();});};
	
	array.AddColumn(t_("Name"), 40);
	array.AddColumn(t_("Line type"), 40);
	array.AddColumn(t_("From"), 40);
	array.AddColumn(t_("To"), 40);
	array.AddColumn(t_("Length [m]"), 40);
	array.AddColumn(t_("Num. Segs."), 40);
	array.AddColumn(t_("From x"), 10);
	array.AddColumn(t_("From y"), 10);
	array.AddColumn(t_("From z"), 10);
	array.AddColumn(t_("To x"), 10);
	array.AddColumn(t_("To y"), 10);
	array.AddColumn(t_("To z"), 10);
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
	
	int col = 0;		
	array.Set(id, col++, ~edName);
	array.Set(id, col++, ~dropLineType);
	array.Set(id, col++, ~dropFrom);
	array.Set(id, col++, ~dropTo);
	array.Set(id, col++, ~edLength);
	array.Set(id, col++, ~edNumSeg);

	Mooring &mooring = *pmooring;

	const Mooring::Connection *pfrom = mooring.GetConnectionP(~dropFrom);
	if (pfrom) {
		array.Set(id, col++, pfrom->x);
		array.Set(id, col++, pfrom->y);
		array.Set(id, col++, pfrom->z);
	} else {
		array.Set(id, col++, "-");
		array.Set(id, col++, "-");
		array.Set(id, col++, "-");
	}
	const Mooring::Connection *pto = mooring.GetConnectionP(~dropTo);
	if (pfrom) {
		array.Set(id, col++, pto->x);
		array.Set(id, col++, pto->y);
		array.Set(id, col++, pto->z);
	} else {
		array.Set(id, col++, "-");
		array.Set(id, col++, "-");
		array.Set(id, col++, "-");
	}
	
	Save();
	String error = mooring.Test();
	if (!error.IsEmpty()) {
		BEM::PrintError(error);
		Load();
		return false;
	}
	
	GetParentCtrl<MainMoor>(this).OnUpdate();
	
	return true;
}

void MainMoor_LineProperties::ArrayOnCursor() {
	int id = array.GetCursor();
	if (id < 0)
		return;
	
	int col = 0;
	edName 		 <<= array.Get(id, col++);
	dropLineType <<= array.Get(id, col++);
	dropFrom 	 <<= array.Get(id, col++);
	dropTo		 <<= array.Get(id, col++);
	edLength 	 <<= ScanDouble(array.Get(id, col++).ToString());
	edNumSeg 	 <<= ScanInt(array.Get(id, col++).ToString());
}

void MainMoor_LineProperties::ArrayClear() {
	edName <<= Null;
	dropLineType <<= Null;
	edLength <<= Null;
	edNumSeg <<= Null;
	dropFrom <<= Null;
	dropTo <<= Null;
}

void MainMoor_LineProperties::Load() {
	InitArray();
	const auto &mooring = *pmooring;
	
	for (int id = 0; id < mooring.lineProperties.size(); ++id) {
		const auto &val = mooring.lineProperties[id];
		String fx, fy, fz, tx, ty, tz;
		const Mooring::Connection *pfrom = mooring.GetConnectionP(val.from);
		if (pfrom) {
			fx = FormatDouble(pfrom->x);
			fy = FormatDouble(pfrom->y);
			fz = FormatDouble(pfrom->z);
		} else
			fx = fy = fz = "-";
		const Mooring::Connection *pto = mooring.GetConnectionP(val.to);
		if (pto) {
			tx = FormatDouble(pto->x);
			ty = FormatDouble(pto->y);
			tz = FormatDouble(pto->z);
		} else
			tx = ty = tz = "-";	
		array.Add(val.name, val.nameType, val.from, val.to, val.length, val.numseg, fx, fy, fz, tx, ty, tz);
	}
	LoadDrop();
	if (mooring.lineProperties.size() > 0)
		array.SetCursor(0);
}

void MainMoor_LineProperties::Save() {
	auto &mooring = *pmooring;
	
	mooring.lineProperties.SetCount(array.GetCount());
	for (int id = 0; id < mooring.lineProperties.size(); ++id) {
		int col = 0;
		mooring.lineProperties[id].name 	= array.Get(id, col++);
		mooring.lineProperties[id].nameType = array.Get(id, col++);
		mooring.lineProperties[id].from 	= array.Get(id, col++);
		mooring.lineProperties[id].to 		= array.Get(id, col++);
		mooring.lineProperties[id].length 	= ScanDouble(array.Get(id, col++).ToString());
		mooring.lineProperties[id].numseg 	= ScanInt(array.Get(id, col++).ToString());
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

void MainMoor_LineProperties::Renumber() {
	auto &mooring = *pmooring;
			
	for (int i = 0; i < mooring.lineProperties.size(); ++i)	
		mooring.lineProperties[i].name = FormatInt(i+1);
	Load();
}

void MainMoor_LineProperties::DeleteUnused() {
	auto &mooring = *pmooring;
	
	for (int i = mooring.lineProperties.size()-1; i >= 0; --i)	{
		auto &line = mooring.lineProperties[i];
		bool foundFrom = false, foundTo = false;
		for (const auto &con : mooring.connections) {
			if (con.name == line.from)
				foundFrom = true;
			else if (con.name == line.to)
				foundTo = true;
			if (foundFrom && foundTo)
				break;
		}
		if (!foundFrom || !foundTo)
			mooring.lineProperties.Remove(i);
	}
	Load();	
}

