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

void MainViewData::Init() {
	Add(tab.SizePos());	
}

void MainViewData::OnAddedModel(MainView &mainView) {
	int id = Bem().surfs.size()-1;
	
	MainViewDataEach &model = models.Add();
	model.Init(Bem().surfs[id], mainView);
	tab.Add(model.SizePos(), Bem().surfs[id].name);
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
		tab.Add(model.SizePos(), Bem().surfs[i].name);		
	}
}

void MainViewDataEach::DataSourceFacets::Init(Mesh &_mesh, int _col, bool _all) {
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
			if (iq >= pmesh->mesh.panels.size())
				return Null;
			if (col == 3 && pmesh->mesh.panels[iq].IsTriangle())
				return "-";
			else
				return pmesh->mesh.panels[iq].id[col]+1;
		} else {
			if (iq >= pmesh->under.panels.size())
				return Null;
			if (col == 3 && pmesh->under.panels[iq].IsTriangle())
				return "-";
			else
				return pmesh->under.panels[iq].id[col]+1;
		}
	}
}

void MainViewDataEach::DataSourceNodes::Init(Mesh &_mesh, int _xyz, int _origMovedUnder) {
	pmesh = &_mesh;	
	xyz = _xyz;
	origMovedUnder = _origMovedUnder;
}

Value MainViewDataEach::DataSourceNodes::Format(const Value& q) const {
	ASSERT(pmesh);
	int iq = q;
	if (origMovedUnder == 0 && pmesh->mesh.nodes.size() <= iq)
		return Null;
	if (origMovedUnder == 1 && pmesh->under.nodes.size() <= iq)
		return Null;
	
	const Point3D &p = origMovedUnder == 0 ? pmesh->mesh.nodes[iq] : pmesh->under.nodes[iq];
	if (xyz == -1)
		return iq + 1;
	else if (xyz == 0)
		return p.x;
	else if (xyz == 1)
		return p.y;
	else
		return p.z;
}

void MainViewDataEach::UpdateStatus(bool under) {
	MainMesh &mainMesh = GetDefinedParent<MainMesh>(this);
	
	bool show;
	if (!under) {
		show = mainMesh.GetShowMesh();
		selectedPanels = ArrayCtrlSelectedGet(arrayFacetsAll2.array);
		selectedNodes = ArrayCtrlSelectedGet(arrayNodesMoved.array);
	} else {
		show = mainMesh.GetShowUnderwater();
		selectedPanels = ArrayCtrlSelectedGet(arrayFacetsUnder.array);
		selectedNodes = ArrayCtrlSelectedGet(arrayNodesUnder.array);
	}
	int numPanels = selectedPanels.size();
	int numNodes  = selectedNodes.size();
	String strPanels = numPanels > 0 ? FormatInt(numPanels) : S(t_("no"));
	String strNodes  = numNodes > 0  ? FormatInt(numNodes)  : S(t_("no"));
	
	if (numPanels + numNodes > 0) {
		if (!show)
			status.Set(Format(t_("%s is hidden in Plot menu so selection will not be shown"), 
						under ? t_("Underwater mesh") : t_("Mesh")));
		else
			status.Set(Format(t_("Selected %s panels and %s nodes"), strPanels, strNodes));
	} else
		status.Set("");
}

void MainViewDataEach::Init(Mesh &_mesh, MainView &mainView) {
	CtrlLayout(arrayFacetsAll2);
	CtrlLayout(arrayFacetsUnder);
	CtrlLayout(arrayNodesMoved);
	CtrlLayout(arrayNodesUnder);
	
	arrayFacetsAll2.title.SetText(t_("Facet node ids"));
	arrayFacetsAll2.array.MultiSelect().SetLineCy(EditField::GetStdHeight()).HeaderObject();
	arrayFacetsAll2.array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, arrayFacetsAll2.array);};
	dataSourceFacetsAll.SetCount(5);
	dataSourceFacetsAll[0].Init(_mesh, -1, true);
	arrayFacetsAll2.array.AddRowNumColumn(t_("#panel"), 60).SetConvert(dataSourceFacetsAll[0]);
	for (int c = 0; c < 4; ++c) {
		dataSourceFacetsAll[c+1].Init(_mesh, c, true);
		arrayFacetsAll2.array.AddRowNumColumn(Format(t_("#%d"), c+1), 60).SetConvert(dataSourceFacetsAll[c+1]);
	}
	arrayFacetsAll2.array.WhenSel = [&] {
		UpdateStatus(false);
		_mesh.mesh.SelPanels(selectedPanels);	
		arrayFacetsUnder.array.ClearSelection();	
		mainView.gl.Refresh();		
		lastSel = 0;
	};
		
	arrayFacetsUnder.title.SetText(t_("Facet node ids"));
	arrayFacetsUnder.array.MultiSelect().SetLineCy(EditField::GetStdHeight()).HeaderObject();
	arrayFacetsUnder.array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, arrayFacetsUnder.array);};
	dataSourceFacetsUnder.SetCount(5);
	dataSourceFacetsUnder[0].Init(_mesh, -1, true);
	arrayFacetsUnder.array.AddRowNumColumn(t_("#panel"), 60).SetConvert(dataSourceFacetsUnder[0]);
	for (int c = 0; c < 4; ++c) {
		dataSourceFacetsUnder[c+1].Init(_mesh, c, false);
		arrayFacetsUnder.array.AddRowNumColumn(Format(t_("#%d"), c+1), 60).SetConvert(dataSourceFacetsUnder[c+1]);
	}
	arrayFacetsUnder.array.WhenSel = [&] {
		UpdateStatus(true);
		_mesh.under.SelPanels(selectedPanels);
		arrayFacetsAll2.array.ClearSelection();	
		mainView.gl.Refresh();	
		lastSel = 1;	
	};
	
	const char *xyz[] = {"x", "y", "z"};

	arrayNodesMoved.title.SetText(t_("Node coordinates"));
	arrayNodesMoved.array.MultiSelect().SetLineCy(EditField::GetStdHeight()).HeaderObject();
	arrayNodesMoved.array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, arrayNodesMoved.array);};
	dataSourceNodesMoved.SetCount(4);
	dataSourceNodesMoved[0].Init(_mesh, -1, 1);
	arrayNodesMoved.array.AddRowNumColumn(t_("#node"), 60).SetConvert(dataSourceNodesMoved[0]);
	for (int c = 0; c < 3; ++c) {
		dataSourceNodesMoved[c+1].Init(_mesh, c, 1);
		arrayNodesMoved.array.AddRowNumColumn(Format(t_("%s"), xyz[c]), 80).SetConvert(dataSourceNodesMoved[c+1]);
	}
	arrayNodesMoved.array.WhenSel = [&] {
		UpdateStatus(false);
		_mesh.mesh.SelNodes(selectedNodes);
		arrayNodesUnder.array.ClearSelection();
		mainView.gl.Refresh();	
		lastSel = 2;
	};
		
	arrayNodesUnder.title.SetText(t_("Node coordinates"));
	arrayNodesUnder.array.MultiSelect().SetLineCy(EditField::GetStdHeight()).HeaderObject();
	arrayNodesUnder.array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, arrayNodesUnder.array);};
	dataSourceNodesMoved.SetCount(4);
	dataSourceNodesMoved[0].Init(_mesh, -1, 2);
	arrayNodesUnder.array.AddRowNumColumn(t_("#node"), 60).SetConvert(dataSourceNodesMoved[0]);
	for (int c = 0; c < 3; ++c) {
		dataSourceNodesMoved[c+1].Init(_mesh, c, 2);
		arrayNodesUnder.array.AddRowNumColumn(Format(t_("%s"), xyz[c]), 80).SetConvert(dataSourceNodesMoved[c+1]);
	}
	arrayNodesUnder.array.WhenSel = [&] {
		UpdateStatus(true);
		_mesh.under.SelNodes(selectedNodes);	
		arrayNodesMoved.array.ClearSelection();
		mainView.gl.Refresh();	
		lastSel = 3;
	};
	
	moved.Horz(arrayFacetsAll2.SizePos(), arrayNodesMoved.SizePos());
	movedUnder.Horz(arrayFacetsUnder.SizePos(), arrayNodesUnder.SizePos());	
	  					 
	tab.Add(moved.SizePos(), t_("All mesh"));
	tab.Add(movedUnder.SizePos(), t_("Only underwater"));
	Add(tab.SizePos());	
	AddFrame(status);
	
	OnRefresh();
	timeCallback.Set(-1000, THISBACK(OnTimer));
}

void MainViewDataEach::OnRefresh() {
	const Mesh &mesh = dataSourceFacetsAll[0].GetMesh();
	int num;
	
	num = mesh.mesh.panels.size();
	arrayFacetsAll2.array.GoBegin();
	arrayFacetsAll2.array.Clear();
	arrayFacetsAll2.array.ClearSelection();
	arrayFacetsAll2.array.SetVirtualCount(num);
	arrayFacetsAll2.array.Refresh();
	arrayFacetsAll2.numRows.SetText(FormatInt(num));
		
	num = mesh.under.panels.size();
	arrayFacetsUnder.array.Clear();
	arrayFacetsUnder.array.ClearSelection();
	arrayFacetsUnder.array.SetVirtualCount(num);
	arrayFacetsUnder.array.Refresh();
	arrayFacetsUnder.numRows.SetText(FormatInt(num));
	
	num = mesh.mesh.nodes.size();
	arrayNodesMoved.array.Clear();
	arrayNodesMoved.array.ClearSelection();
	arrayNodesMoved.array.SetVirtualCount(num);
	arrayNodesMoved.array.Refresh();
	arrayNodesMoved.numRows.SetText(FormatInt(num));
	
	num = mesh.under.nodes.	size();
	arrayNodesUnder.array.Clear();
	arrayNodesUnder.array.ClearSelection();
	arrayNodesUnder.array.SetVirtualCount(num);
	arrayNodesUnder.array.Refresh();
	arrayNodesUnder.numRows.SetText(FormatInt(num));
}

void MainViewDataEach::OnTimer() {
	switch (lastSel) {
	case 0:	arrayFacetsAll2.array.WhenSel();	break;
	case 1:	arrayFacetsUnder.array.WhenSel();	break;
	case 2:	arrayNodesMoved.array.WhenSel();	break;
	case 3:	arrayNodesUnder.array.WhenSel();	break;
	}
}

void MainMeshW::Init(MainMesh &_mesh, const Image &icon, const Image &largeIcon, Function <void()> _WhenClose) {
	WhenClose = _WhenClose;
	LoadFromJson(mesh, StoreAsJson(_mesh));
	mesh.Init();
	Add(mesh.SizePos());
	Title(t_("BEMRosetta Mesh Viewer")).Sizeable().Zoomable().Icon(icon, largeIcon);
}


void VideoCtrl::Init() {
	CtrlLayout(*this);
	
	butPlay <<= THISBACK(OnPlay);
	butStop <<= THISBACK(OnStop);
	butStop.Disable();
	butRecord <<= THISBACK(OnRecord);
	
	if (IsNull(deltaT)) 
		deltaT <<= 0.5;
	
	editTime << "0";
}

void VideoCtrl::OnPlay() {
	butPlay.Disable();
	butStop.Enable();
	butRecord.Disable();
	
	if (playing) 
		return;	
	if(ExistsTimeCallback(1))
		return;
	playing = true;
	SetTimeCallback(-200, THISBACK(TimerFun), 1);
	time.Reset();
	
}

void VideoCtrl::OnStop() {
	butPlay.Enable();        
    butStop.Disable();
    butRecord.Enable();
	
	KillTimeCallback(1);
		
    playing = false; 
}

void VideoCtrl::OnRecord() {
	recording = !recording;
	
	butPlay.Enable(!recording);
	deltaT.Enable(!recording);
    butStop.Disable();
    
	if (!recording) {
		butRecord.SetLabel(t_("Record"));
		KillTimeCallback(1);
	} else {
		video.Create<BasicVideoSequence>();
		dT = deltaT;
		butRecord.SetLabel(t_("..."));
		SetTimeCallback(-200, THISBACK(TimerFun), 1);
		time.Reset();
	}
}

void VideoCtrl::TimerFun() {
	if (!playing && !recording)
		return;

	editTime <<= SecondsToString(time.Seconds(), 1);
}

