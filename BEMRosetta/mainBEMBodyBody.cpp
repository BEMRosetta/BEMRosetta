// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2024, the BEMRosetta author and contributors
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


void BodyBody::Init() {
	CtrlLayout(*this);
	
	nodes.MultiSelect().SetLineCy(EditField::GetStdHeight()).HeaderObject();
	nodes.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, nodes, true);};
	
	panels.MultiSelect().SetLineCy(EditField::GetStdHeight()).HeaderObject().Absolute();
	panels.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, panels, true);};
	
	dropFreq.WhenAction = [&] {
		int ifr = dropFreq.GetData();
		for (auto &d : grd.dataSourcePanels) 
			d.ifr = ifr;
	
		Refresh();
	};
	
	dropHead.WhenAction = [&] {
		int ih = dropHead.GetData();
		for (auto &d : grd.dataSourcePanels) 
			d.ih = ih;
		
		grd.UpdatePanelHeaders();
		ArrayCtrlVirtual_UpdateHeaders(panels, grd.grdPanels);
		
		Refresh();
	};
	
	opPotPress <<= 0;
	opPotPress.WhenAction = [&] {
		for (auto &d : grd.dataSourcePanels) 
			d.pot = opPotPress == 0;
		
		grd.UpdatePanelHeaders();
		ArrayCtrlVirtual_UpdateHeaders(panels, grd.grdPanels);
		
		Refresh();
	};
}

void BodyBody::Load(int idx, int ib) {
	const Hydro &hy = Bem().hydros[idx];
	
	int nNodes, nPanels;
	grd.Load(idx, ib, nNodes, nPanels);
	
	numNodes <<= nNodes;
	ArrayCtrlVirtual(nodes, grd.grdNodes);
	
	numPanels <<= nPanels;	
	ArrayCtrlVirtual(panels, grd.grdPanels);
	
	double freq = dropFreq.GetValue();
	dropFreq.Clear();
	bool showFreq = hy.IsLoadedPotsRad() || hy.IsLoadedPotsDif() || hy.IsLoadedPotsInc() || hy.IsLoadedPotsIncBMR();
	labFreq.Enable (showFreq);
	dropFreq.Enable(showFreq);
	
	if (showFreq) {
		for (int ifr = 0; ifr < hy.dt.Nf; ++ifr)
			dropFreq.Add(ifr, hy.dt.w[ifr]);
		if (!IsNull(freq))
			dropFreq.SetValue(freq);
		else
			dropFreq.SetData(0);
	}
	
	double head = dropHead.GetValue();
	dropHead.Clear();
	bool showHead = hy.IsLoadedPotsDif() || hy.IsLoadedPotsInc() || hy.IsLoadedPotsIncBMR();
	labHead.Enable (showHead);
	dropHead.Enable(showHead);
	
	if (showHead) {
		for (int ih = 0; ih < hy.dt.Nh; ++ih)
			dropHead.Add(ih, hy.dt.head[ih]);
		if (!IsNull(head))
			dropHead.SetValue(head);
		else
			dropHead.SetData(0);
	}	
}
	
void MainBodyTable::Init() {
	Add(tab.SizePos());
}

bool MainBodyTable::Load() {
	try {
		MainBEM &mbm = GetDefinedParent<MainBEM>(this);
		int idx = ArrayModel_IndexHydro(mbm.listLoaded);
	
		tab.Reset();
		bodies.Clear();
		
		UArray<Hydro> &hydros = Bem().hydros; 
		if (hydros.IsEmpty() || idx < 0) 
			return false;
		
		const Hydro &hy = hydros[idx];
		if (hy.dt.msh.IsEmpty() || hy.dt.msh[0].dt.mesh.panels.IsEmpty())
			return false;
		
		for (int ib = 0; ib < hy.dt.Nb; ++ib) {
			BodyBody &b = bodies.Add();
			b.Init();
			b.Load(idx, ib);
			tab.Add(b.SizePos(), Format("%d. %s", ib+1, hy.dt.msh[ib].dt.name));
		}
		return true;
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
		return false;
	}
}