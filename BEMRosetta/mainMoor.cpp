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

	EM().Log("Initialising mooring...");

	splitter.Horz(left.SizePos(), right.SizePos());
	splitter.SetPos(5000, 0);

	left.view.SetShowMesh(SurfaceView::SHOW_MESH_FACES);
	left.view.SetBackgroundColor(::Color(220, 220, 230)).SetLineThickness(1);
	left.view.SetShowColor(SurfaceView::SHOW_DARKER)
		.SetLightDir(Point3D(1, 0, 1));
	left.view.SetCentre(Point3D(0, 0, 0)).SetRotationXYZ();
	left.view.SetSort(false);		// No hidden items
	
	CtrlLayout(right);	
	
	EM().Log("Initialising mooring... 1");
	
	right.edDepth <<= 100;
	right.edDepth.WhenAction = [&] {OnUpdate();};
	left.edNumSeg <<= 50;
	left.edNumSeg.WhenAction = [&] {OnUpdate();};
	right.edDt <<= 0.001;
	right.edElong <<= 0;
	
	left.arrayPos.Reset();
	
	right.showMin.WhenAction = [&] {OnUpdate();};
	
	EM().Log("Initialising mooring... lineTypes");
	lineTypes.Init(mooring);
	EM().Log("Initialising mooring... lineProperties");
	lineProperties.Init(mooring);
	EM().Log("Initialising mooring... lineConnections");
	lineConnections.Init(mooring);
	EM().Log("Initialising mooring... lineVessels");
	lineVessels.Init(mooring);

	right.tab.Add(lineTypes.SizePos(), t_("Line Types"));
	right.tab.Add(lineConnections.SizePos(), t_("Connections"));
	right.tab.Add(lineProperties.SizePos(), t_("Line Properties"));
	right.tab.Add(lineVessels.SizePos(), t_("Vessels"));
	right.tab.WhenSet = [&] {
		static bool wasAtVessel = false;
		
		lineTypes.Save();
		lineProperties.Save();
		lineConnections.Save();
		lineVessels.Save();
		if (right.tab.IsAt(lineProperties)) {
			right.butRenumber.Tip(t_("Renumber lines consecutively")).Show();
			right.butDelete.Tip(t_("Remove lines connected to removed points")).Show();
			lineProperties.LoadDrop();
		} else if (right.tab.IsAt(lineConnections)) {
			right.butRenumber.Tip(t_("Renumber connection points consecutively")).Show();
			right.butDelete.Tip(t_("Remove unused connections points")).Show();
			lineConnections.LoadDrop();	
			OnUpdate();
		} else {
			right.butRenumber.Hide();
			right.butDelete.Hide();
		}
		
		if (wasAtVessel) {
			LoadVesselPositionArray();		// Modifying vessels forces update
			OnUpdate();
		}
		wasAtVessel = right.tab.IsAt(lineVessels);
	};
	EM().Log("Initialising mooring... 2");
	
	const String moorFiles = ".json";
	String moorFilesAst = clone(moorFiles);
	moorFilesAst.Replace(".", "*.");
	right.fileMoor.Type(Format("All supported mooring files (%s)", moorFiles), moorFilesAst);
	right.fileMoor.AllFilesType();
	right.fileMoor.BrowseRightWidth(40).UseOpenFolder(true).BrowseOpenFolderWidth(10);
	right.fileMoor.								WhenChange = [&] {OnLoad(); return true;}; 
	right.butLoad.Tip(t_("Loads mooring file")).WhenAction = [&] {OnLoad();};
	right.butSave.Tip(t_("Saves mooring file")).WhenAction = [&] {OnSave();};
	//right.butUpdate.Tip(t_("Update mooring"))  .WhenAction = [&] {OnUpdate();};
	
	EM().Log("Initialising mooring... 3");
	
	right.arrayresults.Reset();
	right.arrayresults.SetLineCy(EditField::GetStdHeight()).MultiSelect();
	right.arrayresults.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, right.arrayresults, false);};
	right.arrayresults.AddColumn("", 5);
	right.arrayresults.AddColumn("", 20);
	right.arrayresults.AddColumn("", 10);
	right.arrayresults.AddColumn("", 10);
	right.arrayresults.AddColumn("", 10);
	right.arrayresults.AddColumn("", 10);
	right.arrayresults.AddColumn("", 10);
	right.arrayresults.AddColumn("", 10);
	right.arrayresults.AddColumn("", 10);
	right.arrayresults.AddColumn("", 10);
	right.arrayresults.AddColumn("", 10);
	right.arrayresults.AddColumn("", 10);
#ifndef flagDEBUG
	right.arrayresults.Hide();
	right.labresults.Hide();
#endif
	EM().Log("Initialising mooring... 3.1");
	right.opMooring <<= 0;
	right.opMooring.Tip(t_("Outputs the mooring lines intermediate positions"));
	right.edDt.     Tip(t_("Time step to use in mooring integration"));
	EM().Log("Initialising mooring... 3.2");
	right.dropExport.Add("json").Add("MoorDyn");
	right.dropExport.WhenAction = [&] {
		right.opMooring.Enable(right.dropExport.GetData() == "MoorDyn");
		right.edDt.     Enable(right.dropExport.GetData() == "MoorDyn");
	};
	EM().Log(Format("Initialising mooring... 3.3 '%d'", dropExportId));
	if (dropExportId < 0 || dropExportId > 1)
		dropExportId = 0;
	EM().Log(Format("Initialising mooring... 3.3 '%d'", dropExportId));
	right.dropExport.SetIndex(dropExportId);
	right.dropExport.WhenAction();
	
	EM().Log("Initialising mooring... 4");
	
	right.butRenumber.WhenAction = [&] {Renumber(right.tab.IsAt(lineProperties));		FullRefresh(false);};
	right.butRenumber.Hide();
	right.butDelete.WhenAction   = [&] {DeleteUnused(right.tab.IsAt(lineProperties));	FullRefresh(false);};
	right.butDelete.Hide();
	
	left.opConnection.WhenAction = [&]{FullRefresh(false);};
	left.opLines.WhenAction = [&]{FullRefresh(false);};
		
	CtrlLayout(left);
	
	EM().Log("Initialising mooring... 5");
	
	LoadVesselPositionArray();
}

void MainMoor::Renumber(bool lines) {
	if (lines)
		lineProperties.Renumber();
	else
		lineConnections.Renumber();
}

void MainMoor::DeleteUnused(bool lines) {
	if (lines)
		lineProperties.DeleteUnused();
	else
		lineConnections.DeleteUnused();
}

void MainMoor::LoadVesselPositionArray() {
	left.arrayPos.Reset();
	left.arrayPos.SetLineCy(EditField::GetStdHeight());
	left.arrayPos.AddColumn(t_("Vessel"));
	left.arrayPos.AddColumn(t_("x")).Ctrls([this](int i, One<Ctrl>& ctrl) {
		ctrl.Create<EditDoubleSpin>().SetInc(1.).SetFrame(NullFrame()).WhenAction = [&] {OnUpdate();};
	});
	left.arrayPos.AddColumn(t_("y")).Ctrls([this](int i, One<Ctrl>& ctrl) {
		ctrl.Create<EditDoubleSpin>().SetInc(1.).SetFrame(NullFrame()).WhenAction = [&] {OnUpdate();};
	});
	left.arrayPos.HeaderTab(1).SetMargin(0);
	left.arrayPos.HeaderTab(2).SetMargin(0);
#ifdef flagDEBUG
	left.arrayPos.AddColumn(t_("Fx")).SetFormat("%.1f");
	left.arrayPos.AddColumn(t_("Fy")).SetFormat("%.1f");
#endif
	for (const Mooring::Vessel &v : mooring.vessels)
		left.arrayPos.Add(v.name, 0., 0., 0., 0.);
}

bool MainMoor::OnLoad() {
	try {
		if (!mooring.Load(~right.fileMoor)) {
			if (!mooring.LoadMoordyn(~right.fileMoor)) {
				BEM::PrintError(Format("Problem loading %s file", ~right.fileMoor));
				return false;
			}
		}
		lineTypes.Load();
		lineProperties.Load();
		lineConnections.Load();
		lineVessels.Load();
		right.edDepth <<= mooring.depth;
		right.edDt <<= mooring.dtM;
		
		LoadVesselPositionArray();
	} catch (const Exc &e) {
		BEM::PrintError(e);
		return false;
	}
	OnUpdate(true);
	
	return true;
}

bool MainMoor::OnSave() {
	try {
		String format = right.dropExport.GetValue();
		
		FileSel fs;
		
		String fileName = ~right.fileMoor;
		String ext;
		if (format == "json")
			ext = ".json";
		else
			ext = ".dat";
		fs.Type(Format(t_("%s file"), format), "*" + ext);
		fs.ActiveType(0);
		fs.Set(ForceExtSafer(fileName, ext));
		fs.ActiveDir(GetFileDirectory(fileName));
		
		if (!fs.ExecuteSaveAs(Format(t_("Save mooring file as %s"), format)))
			return false;
		
		fileName = ~fs;
		
		lineTypes.Save();
		lineProperties.Save();
		lineConnections.Save();
		lineVessels.Save();
		mooring.depth = ~right.edDepth;
		mooring.dtM = ~right.edDt;
		
		bool ret;
		if (format == "json")
			ret = mooring.Save(fileName);
		else
			ret = mooring.SaveMoordyn(fileName, right.opMooring, right.opFairleads, right.opAnchors);
		if (!ret) {
			BEM::PrintError(Format(t_("Problem saving %s file"), fileName));
			return false;
		}
	} catch (const Exc &e) {
		BEM::PrintError(e);
		return false;
	}
	
	return true;
}

void MainMoor::LoadDragDrop() {
	GuiLock __;
	
	Sort(filesToDrop);
	for (int i = filesToDrop.size()-1; i > 0; --i)
		if (ToLower(GetFileTitle(filesToDrop[i])) == ToLower(GetFileTitle(filesToDrop[i-1])))
			filesToDrop.Remove(i);
	
	if (filesToDrop.IsEmpty())
		return;

	bool followWithErrors = false;
	right.fileMoor <<= filesToDrop[0];
	Status(Format(t_("Loading '%s'"), filesToDrop[0]));
	OnLoad();
}
	
void MainMoor::DragAndDrop(Point , PasteClip& d) {
	GuiLock __;
	if (IsDragAndDropSource())
		return;
	if (AcceptFiles(d)) {
		filesToDrop = GetFiles(d);
		timerDrop.Set(0, [=] {LoadDragDrop();});
		return;
	}
	timerDrop.Kill();
}

bool MainMoor::Key(dword key, int ) {
	GuiLock __;
	if (key == K_CTRL_V) {
		filesToDrop = GetFiles(Ctrl::Clipboard());
		timerDrop.Set(0, [=] {LoadDragDrop();});
		return true;
	}
	return false;
}

void MainMoor::OnUpdate(bool fit) {
	try {
		right.arrayresults.Clear();
		
		lineTypes.Save();
		lineProperties.Save();
		lineConnections.Save();
		lineVessels.Save();
		mooring.depth = ~right.edDepth;
		
		for (int r = 0; r < left.arrayPos.GetCount(); ++r) {
			String name = left.arrayPos.Get(r, 0);
			int id = mooring.FindVessel(name);
			if (id < 0) {
				Exclamation("Unknown vessel");
				return;
			}
			mooring.vessels[id].dx = ScanDouble(left.arrayPos.Get(r, 1).ToString());
			mooring.vessels[id].dy = ScanDouble(left.arrayPos.Get(r, 2).ToString());
		}
		
		try {
			if (!mooring.Calc(Bem().rho, ~left.edNumSeg)) {
				Status(t_("Problem in line calculation"));
				return;
			}
		} catch(Exc e) {
			Status(e);
			return;	
		}
		right.arrayresults.Add(t_("Line"), t_("Status"), t_("Ffairlead [t]"), t_("Angle_v [º]"), 
				t_("Fh vessel [t]"), t_("Fv vessel [t]"), t_("Fv anchor [t]"), t_("Len Seabed [m]"), 
				t_("Ratio Seabed [\%]"), t_("Footprint [m]"), t_("FootprintX [m]"), t_("FootprintY [m]"));
	
		UVector<Pointf> forcevessel(mooring.vessels.size(), Pointf(0, 0));
		
		for (int il = 0; il < mooring.lineProperties.size(); ++il) {
			const auto &line = mooring.lineProperties[il];
			bool broken = line.status == MooringStatus::BROKEN || line.status == MooringStatus::BL_EXCEDEED || line.status == MooringStatus::TAUT || line.status == MooringStatus::CALCULATION_PROBLEM;
			if (!broken) {
				double lenx = abs(line.x[0] - line.x.Top());
				double leny = abs(line.y[0] - line.y.Top());
				double len = sqrt(sqr(lenx) + sqr(leny));
				right.arrayresults.Add(	line.name, 
									InitCaps(MooringStatusStr(line.status)), 
									FormatF(sqrt(sqr(line.fVvessel) + sqr(line.fanchorvessel))/1000./9.8, 1),
									FormatF(ToDeg(atan2(line.fVvessel, line.fanchorvessel)), 1),
									FormatF(line.fanchorvessel/1000./9.8, 1), 
									FormatF(line.fVvessel/1000./9.8, 1), 
									FormatF(line.fVanchor/1000./9.8, 1), 
									FormatF(line.lenonfloor, 2), 
									FormatF(100*line.lenonfloor/line.length, 1) + "%", 
									FormatF(len, 2),
									FormatF(lenx, 2),
									FormatF(leny, 2));
			} else
				right.arrayresults.Add(line.name, InitCaps(MooringStatusStr(line.status)));
						
			int id;
			for (int ic = 0; ic < mooring.connections.size(); ++ic) {
				if (mooring.connections[ic].name == line.from) {
					id = mooring.FindVessel(mooring.connections[ic].where);
					if (id >= 0 && !IsNull(forcevessel[id])) {
						if (broken)
							forcevessel[id] = Null;
						else
							forcevessel[id] += Pointf(line.fanchorvessel*cos(line.theta)/1000./9.8, line.fanchorvessel*sin(line.theta)/1000./9.8);
					}
				}
				if (mooring.connections[ic].name == line.to) {
					id = mooring.FindVessel(mooring.connections[ic].where);
					if (id >= 0 && !IsNull(forcevessel[id])) {
						if (broken)
							forcevessel[id] = Null;
						else
							forcevessel[id] += Pointf(line.fanchorvessel*cos(line.theta)/1000./9.8, line.fanchorvessel*sin(line.theta)/1000./9.8);
					}
				}
			}
		}
		
		for (int r = 0; r < mooring.vessels.size(); ++r) {
			if (!IsNull(forcevessel[r])) {
				left.arrayPos.Set(r, 3, forcevessel[r].x);
				left.arrayPos.Set(r, 4, forcevessel[r].y);
			} else {
				left.arrayPos.Set(r, 3, 0);
				left.arrayPos.Set(r, 4, 0);
			}
		}
		if (!mooring.FindClosest(cl))
			right.minDistance <<= "";
		else {
			const Mooring::LineProperty &prop1 = mooring.lineProperties[cl.line1],
							   			&prop2 = mooring.lineProperties[cl.line2];
			right.minDistance <<= Format("%.1f m. From line %d at (%.1f, %.1f, %.1f) to line %d at (%.1f, %.1f, %.1f)", cl.distance,
										cl.line1+1, prop1.x[cl.point1], prop1.y[cl.point1], prop1.z[cl.point1], 
										cl.line2+1, prop2.x[cl.point2], prop2.y[cl.point2], prop2.z[cl.point2]);
		}
		FullRefresh(fit);		
	} catch (const Exc &e) {	
		BEM::PrintError(e);
		return;
	}
}

void MainMoor::FullRefresh(bool fit) {
	left.view.Clear();
	
	// No sort to hide, painted in order
	// Shadow
	for (const auto &line : mooring.lineProperties) {
		if (line.status == MooringStatus::TAUT) {
			UVector<double> z(line.x.size(), -mooring.depth);
			left.view.PaintLines(line.x, line.y, z, LtGray(), 4);
		} else {
			UVector<double> z(line.x.size(), -mooring.depth);
			UVector<double> x = clone(line.x);
			UVector<double> y = clone(line.y);
			for (int i = 0; i < x.size(); ++i) {
				if (line.z[i] <= -mooring.depth + 0.1)
					x[i] = y[i] = Null;
			}
			left.view.PaintLines(x, y, z, LtGray(), 4);
		}
	}
	// Lines
	for (const auto &line : mooring.lineProperties) {
		if (line.x.size() < 2)
			continue;
		
		::Color c = line.status == MooringStatus::BROKEN || line.status == MooringStatus::BL_EXCEDEED ? LtRed() : LtBlue();
		left.view.PaintLines(line.x, line.y, line.z, c, 1);
		if (left.opLines) {
			double x, y, z;
			if (line.x.size() == 2) {
				x = line.x[0] + 0.33*(line.x[1] - line.x[0]);		// 1/3 of the first connection	
				y = line.y[0] + 0.33*(line.y[1] - line.y[0]);
				z = line.z[0] + 0.33*(line.z[1] - line.z[0]);
			} else {
				int id = int(line.x.size()*1./3.);
				x = line.x[id];
				y = line.y[id];
				z = line.z[id];
			}
			left.view.PaintText(x, y, z, Format("[3@(67.120.120) %s]", line.name));
		}
	}
	// Points
	for (const Mooring::Connection &con : mooring.connections) {
		double dx = 0, dy = 0;
		::Color col;
		int id = mooring.FindVessel(con.where);
		if (id >= 0) {
			dx = mooring.vessels[id].dx;
			dy = mooring.vessels[id].dy;
			if (con.z <= -mooring.depth + 0.1)
				col = ::Color(153, 76, 0);		// Anchor
			else
				col = ::Color(67, 120, 120);
		} else
			col = ::Color(153, 76, 0);			// Anchor
		left.view.PaintCube(con.x + dx, con.y + dy, con.z, 2, col);
		if (left.opConnection)
			left.view.PaintText(con.x + dx, con.y + dy, con.z, Format("[3@(153.76.0) %s]", con.name));
	}
	// Closest
	if (right.showMin && !IsNull(cl.distance)) {
		const Mooring::LineProperty &prop1 = mooring.lineProperties[cl.line1],
						   			&prop2 = mooring.lineProperties[cl.line2];
		left.view.PaintLine(prop1.x[cl.point1], prop1.y[cl.point1], prop1.z[cl.point1],
					   prop2.x[cl.point2], prop2.y[cl.point2], prop2.z[cl.point2], LtRed(), 10);
	}
	
	if (fit) {
		left.view.SetRotation(Value3D(ToRad(-45), 0, ToRad(45)));
		left.view.ZoomToFit();
	}
	left.view.Render();
	Refresh();
}