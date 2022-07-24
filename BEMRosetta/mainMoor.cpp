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

#define IMAGECLASS ImgMoor
#define IMAGEFILE <BEMRosetta/main.iml>
#include <Draw/iml.h>

#include "main.h"


void MainMoor::Init() {
	CtrlLayout(*this);	

	edVessX <<= 0;
	edVessY <<= 0;
	edDepth <<= 100;
	
	edVessX.WhenAction = THISBACK(OnUpdate);
	edVessY.WhenAction = THISBACK(OnUpdate);
	
	lineTypes.Init(mooring);
	lineProperties.Init(mooring);
	lineConnections.Init(mooring);

	tab.Add(lineTypes.SizePos(), t_("Types"));
	tab.Add(lineConnections.SizePos(), t_("Connections"));
	tab.Add(lineProperties.SizePos(), t_("Properties"));
	tab.WhenSet = [&] {
		lineTypes.Save();
		lineProperties.Save();
		lineConnections.Save();
		if (tab.IsAt(lineProperties)) 
			lineProperties.LoadDrop();	
		
	};
	
	const String moorFiles = ".json";
	String moorFilesAst = clone(moorFiles);
	moorFilesAst.Replace(".", "*.");
	file.Type(Format("All supported mooring files (%s)", moorFiles), moorFilesAst);
	file.AllFilesType();
	file.BrowseRightWidth(40).UseOpenFolder(true).BrowseOpenFolderWidth(10);
	butLoad.Tip(t_("Loads mooring file")).WhenAction = [&] {OnLoad();};
	butSave.Tip(t_("Saves mooring file")).WhenAction = [&] {OnSave();};
	butSave.Enable(false);
	butUpdate.Tip(t_("Update mooring"))  .WhenAction = [&] {OnUpdate();};
	
	results.Reset();
	results.SetLineCy(EditField::GetStdHeight());
	results.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, results);};
	results.AddColumn(t_("Line"), 8);
	results.AddColumn(t_("Status"), 10);
	results.AddColumn(t_("Ffair [t]"), 10);
	results.AddColumn(t_("Angle [º]"), 10);
	results.AddColumn(t_("Fhor [t]"), 10);
	results.AddColumn(t_("Fvanchor [t]"), 10);
	results.AddColumn(t_("Fvvessel [t]"), 10);
	results.AddColumn(t_("Lenseabed [m]"), 10);
	results.AddColumn(t_("Footprint [m]"), 10);
	results.AddColumn(t_("FootprintX [m]"), 10);
	results.AddColumn(t_("FootprintY [m]"), 10);
}


bool MainMoor::OnLoad() {
	try {
		if (!mooring.Load(~file)) {
			Exclamation(Format("Problem loading %s file", DeQtf(~file)));
			return false;
		}
		lineTypes.Load();
		lineProperties.Load();
		lineConnections.Load();
		edDepth <<= mooring.depth;
	} catch (const Exc &e) {
		Exclamation(DeQtfLf(e));
		return false;
	}
	butSave.Enable(true);
	
	OnUpdate();
	
	return true;
}

bool MainMoor::OnSave() {
	try {
		lineTypes.Save();
		lineProperties.Save();
		lineConnections.Save();
		mooring.depth = ~edDepth;
		if (!mooring.Save(~file)) {
			Exclamation(Format(t_("Problem loading %s file"), DeQtf(~file)));
			return false;
		}
	} catch (const Exc &e) {
		Exclamation(DeQtfLf(e));
		return false;
	}
	
	return true;
}

void MainMoor::OnUpdate() {
	try {
		results.Clear();
		scatLateral.RemoveAllSeries().SetSequentialXAll();
		scatUp.RemoveAllSeries().SetSequentialXAll();
		
		if (!mooring.Calc(-double(~edVessX), -double(~edVessY), Bem().rho))
			return;
		
		for (int i = 0; i < mooring.lineProperties.size(); ++i) {
			const auto &line = mooring.lineProperties[i];
			double lenx = abs(line.x[0] - line.x.Top());
			double leny = abs(line.y[0] - line.y.Top());
			double len = sqrt(sqr(lenx) + sqr(leny));
			results.Add(line.name, InitCaps(MooringStatusStr(line.status)), 
						FormatF(sqrt(sqr(line.fVvessel) + sqr(line.fanchorvessel))/1000./9.8, 0),
						FormatF(ToDeg(atan2(line.fVvessel, line.fanchorvessel)), 0),
						FormatF(line.fanchorvessel/1000./9.8, 0), FormatF(line.fVanchor/1000./9.8, 0), 
						FormatF(line.fVvessel/1000./9.8, 0), FormatF(line.lenonfloor, 1), len, lenx, leny);
		}
		
		for (int i = 0; i < mooring.lineProperties.size(); ++i) {
			auto &line = mooring.lineProperties[i];
			scatLateral.AddSeries(line.x, line.z).NoMark().Legend(line.name).Units("m", "m").Stroke(1);
		}
		scatLateral.ZoomToFit(true, true);
		
		for (int i = 0; i < mooring.lineProperties.size(); ++i) {
			auto &line = mooring.lineProperties[i];
			scatUp.AddSeries(line.x, line.y).NoMark().Legend(line.name).Units("m", "m").Stroke(1);
		}
		scatUp.ZoomToFit(true, true);
	} catch (const Exc &e) {	
		Exclamation(DeQtfLf(e));
		return;
	}
}
