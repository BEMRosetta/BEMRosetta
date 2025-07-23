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


MainViewData::MainViewData() {
	butClose << [=]   {OnClose();};
	butClear << [=]   {OnClear();};
	butRemove << [=]  {OnRemove();};
	butExtract << [=] {OnExtract();};
}

void MainViewData::OnClose() {
	MainBody &mainBody = GetDefinedParent<MainBody>(this);
	mainBody.CloseSplitter();	
}

void MainViewData::OnClear() {
	int id = tab.Get();
	if (id < 0)
		return; 
	models[id].arrayFacetsAll2.array.ClearSelection(true);
}

void MainViewData::OnRemove() {
	int id = tab.Get();
	if (id < 0)
		return;

	UVector<int> list;
	ArrayCtrl &array = models[id].arrayFacetsAll2.array;
	for (int i = 0; i < array.GetCount(); ++i) {
		if (array.IsSelected(i))
			list << i;
	}
	Body &surf = Bem().surfs[id];
	surf.RemovePanels(list, Bem().rho, Bem().g);
	models[id].arrayFacetsAll2.array.ClearSelection(true);
	
	MainBody &mainBody = GetDefinedParent<MainBody>(this);
	mainBody.UpdateLast(id);	
}

void MainViewData::OnExtract() {
	int id = tab.Get();
	if (id < 0)
		return;

	UVector<int> list;
	ArrayCtrl &array = models[id].arrayFacetsAll2.array;
	for (int i = 0; i < array.GetCount(); ++i) {
		if (array.IsSelected(i))
			list << i;
	}
	Body &surf = Bem().surfs[id];
	MainBody &mainBody = GetDefinedParent<MainBody>(this);
		
	try {
		Bem().AddPanels(surf, list);
		
		Body &msh = Last(Bem().surfs);
		msh.dt.name = t_("Extracted");
		msh.dt.fileName =  "";
		
		msh.AfterLoad(Bem().rho, Bem().g, false, true);
		
		msh.Report(Bem().rho);
		
		mainBody.AddRow(msh);
		mainBody.After();
		OnAddedModel(mainBody.mainView);
		mainBody.OnOpt();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
	
	surf.RemovePanels(list, Bem().rho, Bem().g);
	models[id].arrayFacetsAll2.array.ClearSelection(true);
}

void MainViewData::Init() {
	CtrlLayout(*this);
}

void MainViewData::SelBody(int id, int row, bool select) {
	tab.Set(id);
	models[id].tab.Set(0);
	ArrayCtrl &array = models[id].arrayFacetsAll2.array;
	bool isSelected = array.IsSelected(row);
	if (select && !isSelected)
		array.Select(row, true);	// Select the panel
	else if (!select && isSelected)
		array.Select(row, false);	// Deselect the panel
}

void MainViewData::OnAddedModel(MainView &mainView) {
	int idx = Bem().surfs.size()-1;
	
	MainViewDataEach &model = models.Add();
	model.Init(Bem().surfs[idx], mainView);
	tab.Add(model.SizePos(), Bem().surfs[idx].dt.name);
}

void MainViewData::OnRefresh() {
	for (int i = 0; i < models.size(); ++i)
		models[i].OnRefresh();
}

void MainViewData::Clear() {
	models.Clear();
	tab.Reset();
}

void MainViewData::ReLoad(MainView &mainView) {
	Clear();
	
	for (int i = 0; i < Bem().surfs.size(); ++i)  {
		MainViewDataEach &model = models.Add();
		model.Init(Bem().surfs[i], mainView);
		tab.Add(model.SizePos(), Bem().surfs[i].dt.name);		
	}
}

void MainViewDataEach::DataSourceFacets::Init(Body &_mesh, int _col, bool _all) {
	pmesh = &_mesh;	
	col = _col; 
	all = _all;
}

Value MainViewDataEach::DataSourceFacets::Format(const Value& q) const {
	ASSERT(pmesh);
	int iq = q;
	if (col < 0)
		return iq + 1;
	else {
		if (all) {
			if (iq >= pmesh->dt.mesh.panels.size())
				return Null;
			if (col == 3 && pmesh->dt.mesh.panels[iq].IsTriangle())
				return "-";
			else
				return pmesh->dt.mesh.panels[iq].id[col]+1;
		} else {
			if (iq >= pmesh->dt.under.panels.size())
				return Null;
			if (col == 3 && pmesh->dt.under.panels[iq].IsTriangle())
				return "-";
			else
				return pmesh->dt.under.panels[iq].id[col]+1;
		}
	}
}

void MainViewDataEach::DataSourceNodes::Init(Body &_mesh, int _xyz, int _origMovedUnder) {
	pmesh = &_mesh;	
	xyz = _xyz;
	origMovedUnder = _origMovedUnder;
}

Value MainViewDataEach::DataSourceNodes::Format(const Value& q) const {
	ASSERT(pmesh);
	int iq = q;
	if (origMovedUnder == 0 && pmesh->dt.mesh.nodes.size() <= iq)
		return Null;
	if (origMovedUnder == 1 && pmesh->dt.under.nodes.size() <= iq)
		return Null;
	
	const Point3D &p = origMovedUnder == 0 ? pmesh->dt.mesh.nodes[iq] : pmesh->dt.under.nodes[iq];
	if (xyz == -1)
		return iq + 1;		// id
	else if (xyz == 0)
		return p.x;			// x
	else if (xyz == 1)
		return p.y;			// y
	else
		return p.z;			// z
}

void MainViewDataEach::UpdateStatus() {
	MainBody &mainBody = GetDefinedParent<MainBody>(this);
	
	bool show = mainBody.GetShowBody();
	selectedPanels = ArrayCtrlSelectedGet(arrayFacetsAll2.array);
	selectedNodes = ArrayCtrlSelectedGet(arrayNodesMoved.array);

	int numPanels = selectedPanels.size();
	int numNodes  = selectedNodes.size();
	String strPanels = numPanels > 0 ? FormatInt(numPanels) : S(t_("no"));
	String strNodes  = numNodes > 0  ? FormatInt(numNodes)  : S(t_("no"));
	
	if (numPanels + numNodes > 0) {
		if (!show)
			status.Set(Format(t_("%s is hidden in Plot menu so selection will not be shown"), t_("Mesh")));
		else
			status.Set(Format(t_("Selected %s panels and %s nodes"), strPanels, strNodes));
	} else
		status.Set("");
}

void MainViewDataEach::Init(Body &msh, MainView &mainView) {
	CtrlLayout(arrayFacetsAll2);
	//CtrlLayout(arrayFacetsUnder);
	CtrlLayout(arrayNodesMoved);
	//CtrlLayout(arrayNodesUnder);
	
	arrayFacetsAll2.title.SetText(t_("Facet node ids"));
	arrayFacetsAll2.array.NoFocusSetCursor().MultiSelect().SetLineCy(EditField::GetStdHeight()).HeaderObject();
	arrayFacetsAll2.array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, arrayFacetsAll2.array);};
	dataSourceFacetsAll.SetCount(5);
	dataSourceFacetsAll[0].Init(msh, -1, true);
	arrayFacetsAll2.array.AddRowNumColumn(t_("#panel"), 60).SetConvert(dataSourceFacetsAll[0]);
	for (int c = 0; c < 4; ++c) {
		dataSourceFacetsAll[c+1].Init(msh, c, true);
		arrayFacetsAll2.array.AddRowNumColumn(Format(t_("#%d"), c+1), 60).SetConvert(dataSourceFacetsAll[c+1]);
	}
	arrayFacetsAll2.array.WhenSel = [&] {
		MainBody &mainBody = GetDefinedParent<MainBody>(this);
		UpdateStatus();
		msh.dt.mesh.SelPanels(selectedPanels);	
		//arrayFacetsUnder.array.ClearSelection();	
		mainView.FullRefresh(mainBody);		
		//lastSel = 0;
	};
		
	/*arrayFacetsUnder.title.SetText(t_("Facet node ids"));
	arrayFacetsUnder.array.MultiSelect().SetLineCy(EditField::GetStdHeight()).HeaderObject();
	arrayFacetsUnder.array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, arrayFacetsUnder.array);};
	dataSourceFacetsUnder.SetCount(5);
	dataSourceFacetsUnder[0].Init(msh, -1, true);
	arrayFacetsUnder.array.AddRowNumColumn(t_("#panel"), 60).SetConvert(dataSourceFacetsUnder[0]);
	for (int c = 0; c < 4; ++c) {
		dataSourceFacetsUnder[c+1].Init(msh, c, false);
		arrayFacetsUnder.array.AddRowNumColumn(Format(t_("#%d"), c+1), 60).SetConvert(dataSourceFacetsUnder[c+1]);
	}
	arrayFacetsUnder.array.WhenSel = [&] {
		MainBody &mainBody = GetDefinedParent<MainBody>(this);
		UpdateStatus(true);
		msh.dt.under.SelPanels(selectedPanels);
		arrayFacetsAll2.array.ClearSelection();	
		mainView.FullRefresh(mainBody);			
		//lastSel = 1;	
	};*/
	
	const char *xyz[] = {"x", "y", "z"};

	arrayNodesMoved.title.SetText(t_("Node coordinates"));
	arrayNodesMoved.array.NoFocusSetCursor().MultiSelect().SetLineCy(EditField::GetStdHeight()).HeaderObject();
	arrayNodesMoved.array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, arrayNodesMoved.array);};
	dataSourceNodesMoved.SetCount(4);
	dataSourceNodesMoved[0].Init(msh, -1, 0);
	arrayNodesMoved.array.AddRowNumColumn(t_("#node"), 60).SetConvert(dataSourceNodesMoved[0]);
	for (int c = 0; c < 3; ++c) {
		dataSourceNodesMoved[c+1].Init(msh, c, 0);
		arrayNodesMoved.array.AddRowNumColumn(Format(t_("%s"), xyz[c]), 80).SetConvert(dataSourceNodesMoved[c+1]);
	}
	arrayNodesMoved.array.WhenSel = [&] {
		MainBody &mainBody = GetDefinedParent<MainBody>(this);
		UpdateStatus();
		msh.dt.mesh.SelNodes(selectedNodes);
		//arrayNodesUnder.array.ClearSelection();
		mainView.FullRefresh(mainBody);			
		//lastSel = 2;
	};
	/*	
	arrayNodesUnder.title.SetText(t_("Node coordinates"));
	arrayNodesUnder.array.MultiSelect().SetLineCy(EditField::GetStdHeight()).HeaderObject();
	arrayNodesUnder.array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, arrayNodesUnder.array);};
	dataSourceNodesUnder.SetCount(4);
	dataSourceNodesUnder[0].Init(msh, -1, 1);
	arrayNodesUnder.array.AddRowNumColumn(t_("#node"), 60).SetConvert(dataSourceNodesUnder[0]);
	for (int c = 0; c < 3; ++c) {
		dataSourceNodesUnder[c+1].Init(msh, c, 1);
		arrayNodesUnder.array.AddRowNumColumn(Format(t_("%s"), xyz[c]), 80).SetConvert(dataSourceNodesUnder[c+1]);
	}
	arrayNodesUnder.array.WhenSel = [&] {
		MainBody &mainBody = GetDefinedParent<MainBody>(this);
		UpdateStatus(true);
		msh.dt.under.SelNodes(selectedNodes);	
		arrayNodesMoved.array.ClearSelection();
		mainView.FullRefresh(mainBody);			
		//lastSel = 3;
	};
	*/
	moved.Horz(arrayFacetsAll2.SizePos(), arrayNodesMoved.SizePos());
	//movedUnder.Horz(arrayFacetsUnder.SizePos(), arrayNodesUnder.SizePos());	
	  					 
	tab.Add(moved.SizePos(), t_("All mesh"));
	//tab.Add(movedUnder.SizePos(), t_("Only underwater"));
	Add(tab.SizePos());	
	AddFrame(status);
	
	OnRefresh();
	//timeCallback.Set(-1000, THISBACK(OnTimer));
}

void MainViewDataEach::OnRefresh() {
	const Body &msh = dataSourceFacetsAll[0].GetBody();
	int num;
	
	num = msh.dt.mesh.panels.size();
	arrayFacetsAll2.array.GoBegin();
	arrayFacetsAll2.array.Clear();
	arrayFacetsAll2.array.ClearSelection();
	arrayFacetsAll2.array.SetVirtualCount(num);
	arrayFacetsAll2.array.Refresh();
	arrayFacetsAll2.numRows.SetText(FormatInt(num));
		
	/*num = msh.dt.under.panels.size();
	arrayFacetsUnder.array.Clear();
	arrayFacetsUnder.array.ClearSelection();
	arrayFacetsUnder.array.SetVirtualCount(num);
	arrayFacetsUnder.array.Refresh();
	arrayFacetsUnder.numRows.SetText(FormatInt(num));*/
	
	num = msh.dt.mesh.nodes.size();
	arrayNodesMoved.array.Clear();
	arrayNodesMoved.array.ClearSelection();
	arrayNodesMoved.array.SetVirtualCount(num);
	arrayNodesMoved.array.Refresh();
	arrayNodesMoved.numRows.SetText(FormatInt(num));
	
	/*num = msh.dt.under.nodes.size();
	arrayNodesUnder.array.Clear();
	arrayNodesUnder.array.ClearSelection();
	arrayNodesUnder.array.SetVirtualCount(num);
	arrayNodesUnder.array.Refresh();
	arrayNodesUnder.numRows.SetText(FormatInt(num));*/
}

/*
void MainViewDataEach::OnTimer() {
	switch (lastSel) {
	case 0:	arrayFacetsAll2.array.WhenSel();	break;
	case 1:	arrayFacetsUnder.array.WhenSel();	break;
	case 2:	arrayNodesMoved.array.WhenSel();	break;
	case 3:	arrayNodesUnder.array.WhenSel();	break;
	}
}
*/
