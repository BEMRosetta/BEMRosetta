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

#define IMAGECLASS ImgMoor
#define IMAGEFILE <BEMRosetta/main.iml>
#include <Draw/iml.h>

#include "main.h"


void MainMoor::Init() {
	Add(splitter.SizePos());

	scat.Add(scatLateral, 1, 0).Add(scatUp, 0, 0);
	scat.WhenHeights =[&](int width, int height, UVector<int> &heights) {
		heights[0] = max(0, height - width);
		heights[1] = width;
	};
	splitter.Horz(scat.SizePos(), right.SizePos());
	splitter.SetPos(3000, 0);
	
	CtrlLayout(right);	
	
	right.edVessX <<= 0;
	right.edVessY <<= 0;
	right.edDepth <<= 100;
	
	right.edVessX.WhenAction = THISBACK(OnUpdate);
	right.edVessY.WhenAction = THISBACK(OnUpdate);
	
	lineTypes.Init(mooring);
	lineProperties.Init(mooring);
	lineConnections.Init(mooring);

	right.tab.Add(lineTypes.SizePos(), t_("Types"));
	right.tab.Add(lineConnections.SizePos(), t_("Connections"));
	right.tab.Add(lineProperties.SizePos(), t_("Properties"));
	right.tab.WhenSet = [&] {
		lineTypes.Save();
		lineProperties.Save();
		lineConnections.Save();
		if (right.tab.IsAt(lineProperties)) 
			lineProperties.LoadDrop();	
		
	};
	
	const String moorFiles = ".json";
	String moorFilesAst = clone(moorFiles);
	moorFilesAst.Replace(".", "*.");
	right.fileMoor.Type(Format("All supported mooring files (%s)", moorFiles), moorFilesAst);
	right.fileMoor.AllFilesType();
	right.fileMoor.BrowseRightWidth(40).UseOpenFolder(true).BrowseOpenFolderWidth(10);
	right.butLoad.Tip(t_("Loads mooring file")).WhenAction = [&] {OnLoad();};
	right.butSave.Tip(t_("Saves mooring file")).WhenAction = [&] {OnSave();};
	//butSave.Enable(false);
	right.butUpdate.Tip(t_("Update mooring"))  .WhenAction = [&] {OnUpdate();};
	
	right.results.Reset();
	right.results.SetLineCy(EditField::GetStdHeight()).MultiSelect();
	right.results.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, right.results, false);};
	right.results.AddColumn("", 5);
	right.results.AddColumn("", 20);
	right.results.AddColumn("", 10);
	right.results.AddColumn("", 10);
	right.results.AddColumn("", 10);
	right.results.AddColumn("", 10);
	right.results.AddColumn("", 10);
	right.results.AddColumn("", 10);
	right.results.AddColumn("", 10);
	right.results.AddColumn("", 10);
	right.results.AddColumn("", 10);
	
	scatLateral.SetTitle(t_("Side view")).SetTitleFont(Arial(12)).SetMargin(70, 25, 30, 50);
	scatUp.     SetTitle(t_("Top view")). SetTitleFont(Arial(12)).SetMargin(70, 25, 30, 50);
}


bool MainMoor::OnLoad() {
	try {
		if (!mooring.Load(~right.fileMoor)) {
			BEM::PrintError(Format("Problem loading %s file", DeQtf(~right.fileMoor)));
			return false;
		}
		lineTypes.Load();
		lineProperties.Load();
		lineConnections.Load();
		right.edDepth <<= mooring.depth;
	} catch (const Exc &e) {
		BEM::PrintError(DeQtfLf(e));
		return false;
	}
	//butSave.Enable(true);
	
	OnUpdate();
	
	return true;
}

bool MainMoor::OnSave() {
	try {
		lineTypes.Save();
		lineProperties.Save();
		lineConnections.Save();
		mooring.depth = ~right.edDepth;
		if (!mooring.Save(~right.fileMoor)) {
			BEM::PrintError(Format(t_("Problem loading %s file"), DeQtf(~right.fileMoor)));
			return false;
		}
	} catch (const Exc &e) {
		BEM::PrintError(DeQtfLf(e));
		return false;
	}
	
	return true;
}

void MainMoor::OnUpdate() {
	try {
		right.results.Clear();
		scatLateral.RemoveAllSeries();
		scatUp.RemoveAllSeries();
		
		lineTypes.Save();
		lineProperties.Save();
		lineConnections.Save();
		mooring.depth = ~right.edDepth;
		
		px.Clear();	py.Clear(); pz.Clear();
		for (auto &conn : mooring.connections) {
			px << conn.x + (conn.type == 'v' ? double(~right.edVessX) : 0.);
			py << conn.y + (conn.type == 'v' ? double(~right.edVessY) : 0.);
			pz << conn.z + mooring.depth;
		}
		
		scatLateral.AddSeries(px, pz).NoPlot().MarkStyle<CircleMarkPlot>().NoSeriesLegend();
		scatLateral.ZoomToFit(true, true);
		
		scatUp.AddSeries(px, py).NoPlot().MarkStyle<CircleMarkPlot>().NoSeriesLegend();
		scatUp.ZoomToFit(true, true);
		

		if (!mooring.Calc(double(~right.edVessX), double(~right.edVessY), Bem().rho))
			return;
		
		right.results.Add(t_("Line"), t_("Status"), t_("Ffairlead [t]"), t_("Angle [º]"), 
				t_("Fh vessel [t]"), t_("Fv vessel [t]"), t_("Fv anchor [t]"), t_("Len Seabed [m]"), t_("Footprint [m]"),
				t_("FootprintX [m]"), t_("FootprintY [m]"));
	
		for (int i = 0; i < mooring.lineProperties.size(); ++i) {
			const auto &line = mooring.lineProperties[i];
			if (line.status != MooringStatus::BROKEN && line.status != MooringStatus::BL_EXCEDEED) {
				double lenx = abs(line.x[0] - line.x.Top());
				double leny = abs(line.y[0] - line.y.Top());
				double len = sqrt(sqr(lenx) + sqr(leny));
				right.results.Add(line.name, InitCaps(MooringStatusStr(line.status)), 
							FormatF(sqrt(sqr(line.fVvessel) + sqr(line.fanchorvessel))/1000./9.8, 1),
							FormatF(ToDeg(atan2(line.fVvessel, line.fanchorvessel)), 1),
							FormatF(line.fanchorvessel/1000./9.8, 1), 
							FormatF(line.fVvessel/1000./9.8, 1), 
							FormatF(line.fVanchor/1000./9.8, 1), 
							FormatF(line.lenonfloor, 2), 
							FormatF(len, 2),
							FormatF(lenx, 2),
							FormatF(leny, 2));
			} else
				right.results.Add(line.name, InitCaps(MooringStatusStr(line.status)));
		}
		
		for (int i = 0; i < mooring.lineProperties.size(); ++i) {
			auto &line = mooring.lineProperties[i];
			scatLateral.AddSeries(line.x, line.z).NoMark().Legend(line.name).Units("m", "m").Stroke(1).NoDash();
			if (line.status == MooringStatus::BROKEN || line.status == MooringStatus::BL_EXCEDEED) 
				scatLateral.Dash(ScatterDraw::LINE_DASHED);
		}
		scatLateral.ZoomToFit(true, true);
		
		for (int i = 0; i < mooring.lineProperties.size(); ++i) {
			auto &line = mooring.lineProperties[i];
			scatUp.AddSeries(line.x, line.y).NoMark().Legend(line.name).Units("m", "m").Stroke(1).NoDash();
			if (line.status == MooringStatus::BROKEN || line.status == MooringStatus::BL_EXCEDEED) 
				scatUp.Dash(ScatterDraw::LINE_DASHED);
		}
		scatUp.ZoomToFit(true, true);
	} catch (const Exc &e) {	
		BEM::PrintError(DeQtfLf(e));
		return;
	}
}
