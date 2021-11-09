// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2021, the BEMRosetta author and contributors
#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <TabBar/TabBar.h>

#include <BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;

#include "main.h"


void MainBEM::Init() {
	CtrlLayout(*this);
	
	CtrlLayout(menuOpen);
	menuOpen.file.WhenChange = THISBACK(OnLoad);
	menuOpen.file.BrowseRightWidth(40).UseOpenFolder().BrowseOpenFolderWidth(10);
	menuOpen.butLoad << [&] {menuOpen.file.DoGo();};
	
	ArrayModel_Init(listLoaded).MultiSelect();
	listLoaded.WhenSel = [&] {
		OnMenuConvertArraySel();
		menuFOAMM.OnCursor();
		mainTab.GetItem(mainTab.Find(mainQTF)).Enable(mainQTF.Load());
		UpdateButtons();
	};
	listLoaded.WhenBar = [&](Bar &menu) {
		listLoaded.StdBar(menu);
		menu.Add(listLoaded.GetCount() > 0, t_("Open file folder"), Null, [&]{
			LaunchWebBrowser(GetFileFolder(ArrayModel_GetFileName(listLoaded)));}).Help(t_("Opens file explorer in the file folder"));
		menu.Add(listLoaded.GetCount() > 0, t_("Deselect all"), Null, [&]{listLoaded.ClearSelection();})
			.Help(t_("Deselect all table rows"));
	};

	menuOpen.butRemove.Disable();	
	menuOpen.butRemove <<= THISBACK(OnRemove);
	menuOpen.butRemoveSelected.Disable();	
	menuOpen.butRemoveSelected <<= THISBACK1(OnRemoveSelected, false);
	menuOpen.butJoin.Disable();	
	menuOpen.butJoin <<= THISBACK(OnJoin);
	menuOpen.butDuplicate.Disable();	
	menuOpen.butDuplicate <<= THISBACK(OnDuplicate);
	menuOpen.butDescription.Disable();
	menuOpen.butDescription <<= THISBACK(OnDescription);
	
	CtrlLayout(menuProcess);
	menuProcess.butSymX.Disable();	
	menuProcess.butSymX <<= THISBACK1(OnSymmetrize, true);
	menuProcess.butSymY.Disable();
	menuProcess.butSymY <<= THISBACK1(OnSymmetrize, false);
	menuProcess.butA0.Disable();	
	menuProcess.butA0 <<= THISBACK1(OnKirfAinf, Hydro::PLOT_A0);
	menuProcess.butAinf.Disable();	
	menuProcess.butAinf <<= THISBACK1(OnKirfAinf, Hydro::PLOT_AINF);
	menuProcess.butKirf.Disable();	
	menuProcess.butKirf <<= THISBACK1(OnKirfAinf, Hydro::PLOT_K);
	
	CtrlLayout(menuAdvanced);
	menuAdvanced.butAinfw <<= THISBACK1(OnKirfAinf, Hydro::PLOT_AINFW);
	menuAdvanced.butOgilvie <<= THISBACK(OnOgilvie);
	menuAdvanced.butAinfw.Disable();	
	menuAdvanced.butOgilvie.Disable();
	
	CtrlLayout(menuConvert);
	menuConvert.file.WhenChange = THISBACK(OnConvert);
	menuConvert.file.BrowseRightWidth(40).UseOpenFolder(true).BrowseOpenFolderWidth(10);
	menuConvert.file.SelLoad(false);
	menuConvert.butConvert << [&] {menuConvert.file.DoGo();};
	menuConvert.opt 	   << [&] {OnOpt();};
	
	OnOpt();
	
	CtrlLayout(menuPlot);
	menuPlot.butZoomToFit << [&] {GetSelScatter().ZoomToFit(true, true);};
	menuPlot.autoFit	  << [&] {
		LoadSelTab(Bem());
		menuPlot.fromY0.Enable(~menuPlot.autoFit);
	};
	menuPlot.fromY0 	<< [&] {LoadSelTab(Bem());};
	menuPlot.opwT 		<< [&] {LoadSelTab(Bem());};
	menuPlot.showPoints << [&] {LoadSelTab(Bem());};
	menuPlot.showNdim 	<< [&] {LoadSelTab(Bem());};
	
	OnOpt();
	
	menuFOAMM.Init(*this, mainSetupFOAMM);
	
	OnOpt();
		
	menuTab.Add(menuOpen.SizePos(), 	t_("Load"));
	menuTab.Add(menuProcess.SizePos(), 	t_("Process")).Disable();
	menuTab.Add(menuAdvanced.SizePos(), t_("Advanced")).Disable();
	menuTab.Add(menuConvert.SizePos(), 	t_("Save as")).Disable();
	menuTab.Add(menuPlot.SizePos(), 	t_("Plot")).Disable();
	menuTab.Add(menuFOAMM.SizePos(), 	t_("FOAMM State Space")).Disable();
	
	menuTab.WhenSet = [&] {
		LOGTAB(menuTab);
		bool setupfoamm = false;
		if (menuTab.IsAt(menuFOAMM)) {
			setupfoamm = true;
			if (!FileExists(Bem().foammPath))
				Status(t_("FOAMM not found. Please set FOAMM path in Options"), 10000);	
		} else if (menuTab.IsAt(menuPlot) || menuTab.IsAt(menuOpen) || 
				   menuTab.IsAt(menuProcess) || menuTab.IsAt(menuAdvanced)) 
			setupfoamm = true;
		
		if (!setupfoamm) 
			mainTab.Set(0);
		
		if (menuTab.IsAt(menuFOAMM)) 
			mainTab.Set(mainSetupFOAMM);
		else if (menuTab.IsAt(menuConvert)) 
			listLoaded.WhenSel(); 
		
		ShowMenuPlotItems();
	};
	
	mainTab.WhenSet = [&] {
		LOGTAB(mainTab);
		Upp::Vector<int> ids = ArrayModel_IdsHydro(listLoaded);
		bool plot = true, convertProcess = true, ismenuFOAMM = false;

		if (ids.IsEmpty())
			plot = convertProcess = false;
		else if (mainTab.IsAt(mainStiffness)) {
			plot = false;
			mainStiffness.Load(Bem().hydros, ids);
		} else if (mainTab.IsAt(mainA))
			mainA.Load(Bem(), ids);
		else if (mainTab.IsAt(mainB))
			mainB.Load(Bem(), ids);
		else if (mainTab.IsAt(mainK))
			mainK.Load(Bem(), ids);
		else if (mainTab.IsAt(mainAinfw))
			mainAinfw.Load(Bem(), ids);
		else if (mainTab.IsAt(mainForceSC))
			mainForceSC.Load(Bem(), ids);
		else if (mainTab.IsAt(mainForceFK))
			mainForceFK.Load(Bem(), ids);
		else if (mainTab.IsAt(mainForceEX))
			mainForceEX.Load(Bem(), ids);
		else if (mainTab.IsAt(mainRAO))
			mainRAO.Load(Bem(), ids);
		else if (mainTab.IsAt(mainStateSpace)) {
			mainStateSpace.Load(Bem(), ids);
			ismenuFOAMM = true;
		} else if (mainTab.IsAt(mainSetupFOAMM)) {
			menuTab.Set(menuFOAMM);
			ismenuFOAMM = true;
			if (mainSetupFOAMM.arrayCases.GetCount() == 0 && ids.size() == 1) 
				listLoaded.SetCursor(0);
			menuFOAMM.OnCursor();
		} else if (mainTab.IsAt(mainQTF))
			mainQTF.Load();
		else if (menuTab.IsAt(menuFOAMM)) 
			;
		else 
			plot = false;
		
		TabCtrl::Item& tabMenuPlot = menuTab.GetItem(menuTab.Find(menuPlot));
		tabMenuPlot.Enable(plot);
		TabCtrl::Item& tabMenuProcess = menuTab.GetItem(menuTab.Find(menuProcess));
		tabMenuProcess.Enable(convertProcess);
		TabCtrl::Item& tabMenuAdvanced = menuTab.GetItem(menuTab.Find(menuAdvanced));
		tabMenuAdvanced.Enable(convertProcess);
		TabCtrl::Item& tabMenuConvert = menuTab.GetItem(menuTab.Find(menuConvert));
		tabMenuConvert.Enable(convertProcess);
		if (plot) {
			tabMenuPlot.Text(t_("Plot"));
			tabMenuProcess.Text(t_("Process"));
			tabMenuAdvanced.Text(t_("Advanced"));
		} else {
			tabMenuPlot.Text("");
			tabMenuProcess.Text("");
			tabMenuAdvanced.Text("");
		}
		
		if (convertProcess) {
			tabMenuProcess.Text(t_("Process"));
			tabMenuAdvanced.Text(t_("Advanced"));
			tabMenuConvert.Text(t_("Save as"));
		} else {
			tabMenuProcess.Text("");
			tabMenuConvert.Text("");
			tabMenuAdvanced.Text("");
		}
		TabCtrl::Item& tabMenuFOAMM = menuTab.GetItem(menuTab.Find(menuFOAMM));
		tabMenuFOAMM.Enable(convertProcess);
		if (convertProcess) 
			tabMenuFOAMM.Text(t_("FOAMM State Space"));
		else 
			tabMenuFOAMM.Text("");
		
		if (!ismenuFOAMM && plot && (menuTab.IsAt(menuFOAMM) || menuTab.IsAt(mainSetupFOAMM))) 
			menuTab.Set(menuPlot);
		
		ShowMenuPlotItems();
	};
	mainTab.WhenSet();
	
	mainSummary.Init();
	mainTab.Add(mainSummary.SizePos(), t_("Summary"));
	
	mainArrange.Init();
	mainTab.Add(mainArrange.SizePos(), t_("Arrange DOF")).Disable();
	
	mainStiffness.Init();
	mainTab.Add(mainStiffness.SizePos(), t_("K")).Disable();
	
	mainA.Init(Hydro::DATA_A);
	mainTab.Add(mainA.SizePos(), t_("A")).Disable();
	
	mainB.Init(Hydro::DATA_B);
	mainTab.Add(mainB.SizePos(), t_("B")).Disable();
	
	mainK.Init(Hydro::DATA_K);
	mainTab.Add(mainK.SizePos(), t_("Kirf")).Disable();
	
	mainForceEX.Init(Hydro::DATA_FORCE_EX);
	mainTab.Add(mainForceEX.SizePos(), t_("Fex")).Disable();
	
	mainForceSC.Init(Hydro::DATA_FORCE_SC);
	mainTab.Add(mainForceSC.SizePos(), t_("Fsc")).Disable();
	
	mainForceFK.Init(Hydro::DATA_FORCE_FK);
	mainTab.Add(mainForceFK.SizePos(), t_("Ffk")).Disable();
	
	mainRAO.Init(Hydro::DATA_RAO);
	mainTab.Add(mainRAO.SizePos(), t_("RAO")).Disable();

	mainAinfw.Init(Hydro::DATA_AINFW);
	mainTab.Add(mainAinfw.SizePos(), t_("A∞(ω)")).Disable();
	
	mainSetupFOAMM.Init();
	mainTab.Add(mainSetupFOAMM.SizePos(), t_("Setup FOAMM")).Disable();

	mainStateSpace.Init();
	mainTab.Add(mainStateSpace.SizePos(), t_("State Space")).Disable();

	mainQTF.Init();
	mainTab.Add(mainQTF.SizePos(), t_("QTF")).Disable();
}

void MainBEM::ShowMenuPlotItems() {
	menuPlot.showNdim.Enable();

	bool show = true, showwT = true;
	if (mainTab.IsAt(mainSetupFOAMM)) 
		showwT = false;
	else if (mainTab.IsAt(mainK)) 
		showwT = false;
	else if (mainTab.IsAt(mainQTF))
		show = false;
	
	menuPlot.opwT.Enable(showwT);
	menuPlot.butZoomToFit.Enable(show);
	menuPlot.autoFit.Enable(show);
	menuPlot.fromY0.Enable(show);
	menuPlot.showPoints.Enable(show);
}
	
void MainBEM::OnMenuConvertArraySel() {
	int id = ArrayModel_IdHydro(listLoaded);
	if (id < 0)
		return;
	
	String file = ~menuConvert.file;
	String folder = GetFileFolder(file);
	String ext = ToLower(GetFileExt(file));
	String fileName = GetFileTitle(ArrayModel_GetFileName(listLoaded));
	if (fileName.IsEmpty())
		fileName = ArrayModel_GetTitle(listLoaded);
	file = AppendFileNameX(folder, fileName + ext);
	menuConvert.file <<= file;
}

void MainBEM::InitSerialize(bool ret) {
	if (!ret || IsNull(menuPlot.autoFit)) 
		menuPlot.autoFit = true;
	
	if (!ret || IsNull(menuPlot.fromY0)) 
		menuPlot.fromY0 = false;
	
	if (!ret || IsNull(menuPlot.opwT)) 
		menuPlot.opwT = 0;

	if (!ret || IsNull(menuPlot.showPoints)) 
		menuPlot.showPoints = true;
	
	if (!ret || IsNull(menuPlot.showNdim)) 
		menuPlot.showNdim = false;

	if (!ret || IsNull(menuConvert.opt)) 
		menuConvert.opt = 0;
}

void MainBEM::LoadSelTab(BEMData &bem) {
	Upp::Vector<int> ids = ArrayModel_IdsHydro(listLoaded);
	int id = mainTab.Get();
	if (id == mainTab.Find(mainStateSpace))
		mainStateSpace.Load(bem, ids);
	else if (id == mainTab.Find(mainStiffness))
		mainStiffness.Load(bem.hydros, ids);
	else if (id == mainTab.Find(mainSummary) || id == mainTab.Find(mainArrange))
		;
	else if (id == mainTab.Find(mainSetupFOAMM))
		mainSetupFOAMM.Load();
	else if (id == mainTab.Find(mainQTF))
		mainQTF.Load();
	else 
		GetSelABForce().Load(bem, ids);
}

MainABForce &MainBEM::GetSelABForce() {
	int id = mainTab.Get();
	Ctrl *ctrl = mainTab.GetItem(id).GetSlave();
	if (!ctrl)
		throw Exc(t_("Object not found in GetSelABForce()"));
	if (typeid(MainABForce) != typeid(*ctrl))
		throw Exc(t_("Unexpected type in GetSelABForce()"));
	return *(static_cast<MainABForce*>(ctrl));
}

MainStateSpace &MainBEM::GetSelStateSpace() {
	int id = mainTab.Get();
	Ctrl *ctrl = mainTab.GetItem(id).GetSlave();
	if (!ctrl)
		throw Exc(t_("Object not found in GetSelStateSpace()"));
	if (typeid(MainStateSpace) != typeid(*ctrl))
		throw Exc(t_("Unexpected type in GetSelStateSpace()"));
	return *(static_cast<MainStateSpace*>(ctrl));
}

ScatterCtrl &MainBEM::GetSelScatter() {
	int id = mainTab.Get();
	if (id == mainTab.Find(mainStateSpace)) {
		TabCtrl &tab = GetSelStateSpace().tab;
		Ctrl *ctrl = tab.GetItem(tab.Get()).GetSlave();
		if (!ctrl)
			throw Exc(t_("Object not found in GetSelScatter(1)"));
		if (typeid(MainStateSpacePlot) != typeid(*ctrl))
			throw Exc(t_("Unexpected type in GetSelScatter(1)"));		
		MainStateSpacePlot *mainStateSpacePlot = static_cast<MainStateSpacePlot*>(ctrl);
		return mainStateSpacePlot->mainPlot.scatt;
	} else {
		TabCtrl &tab = GetSelABForce().tab;
		Ctrl *ctrl = tab.GetItem(tab.Get()).GetSlave();
		if (!ctrl)
			throw Exc(t_("Object not found in GetSelScatter(2)"));
		if (typeid(MainPlot) != typeid(*ctrl))
			throw Exc(t_("Unexpected type in GetSelScatter(2)"));		
		MainPlot *mainPlot = static_cast<MainPlot*>(ctrl);
		return mainPlot->scatt;
	}
}

void MainBEM::OnOpt() {
	menuOpen.file.ClearTypes();

	menuOpen.file.Type(Format(t_("All supported BEM files (%s)"), Bem().bemFilesExt), Bem().bemFilesAst);
	menuOpen.file.AllFilesType();
	String extOpen = ToLower(GetFileExt(menuOpen.file.GetData().ToString()));
	if (extOpen.IsEmpty())
		menuOpen.file.ActiveType(0);
	else if (Bem().bemFilesExt.Find(extOpen) >= 0)
		menuOpen.file.ActiveType(0);
	else
		menuOpen.file.ActiveType(1);
	
	menuConvert.file.ClearTypes();
	switch (menuConvert.opt) {
	case 0:	menuConvert.file <<= ForceExtSafe(~menuConvert.file, ".out"); 	
			menuConvert.file.Type(t_("Wamit .out file"), "*.out");
			break;	
	case 1:	menuConvert.file <<= ForceExtSafe(~menuConvert.file, ".1"); 	
			menuConvert.file.Type(t_("Wamit .1.2.3.4.hst.12s.12d file"), "*.1 *.3 *.hst *.4 *.12s *.12d");
			break;
	case 2:	menuConvert.file <<= ForceExtSafe(~menuConvert.file, ".dat"); 
			menuConvert.file.Type(t_("FAST HydroDyn file"), "*.dat");
			break;
	case 3:	menuConvert.file <<= ForceExtSafe(~menuConvert.file, ".bem"); 
			menuConvert.file.Type(t_("BEMRosetta file"), "*.bem");
			break;
	default:menuConvert.file.Type(t_("All converted files"), "*.1 *.3 *.hst *.4 *.12s *.12d *.dat *.bem *.out");
			break;
	}
	String extConv = ToLower(GetFileExt(menuConvert.file.GetData().ToString()));
	if (extConv.IsEmpty())
		extConv = "-";
	if (String(".1 .2 .3 .hst .4 .12s .12d").Find(extConv) >= 0)
		menuConvert.file.ActiveType(0);
	else if (String(".dat").Find(extConv) >= 0)
		menuConvert.file.ActiveType(1);
	else
		menuConvert.file.ActiveType(2);
}

bool MainBEM::OnLoad() {
	String file = ~menuOpen.file;
	return OnLoadFile(file);
}

bool MainBEM::OnLoadFile(String file) {
	GuiLock __;
	try {
		Progress progress(t_("Loading BEM files..."), 100); 
		
		for (int i = 0; i < Bem().hydros.size(); ++i) {
			if (ForceExt(Bem().hydros[i].hd().file, ".") == ForceExt(file, ".")) {
				if (!PromptYesNo(t_("Model is already loaded") + S("&") + t_("Do you wish to open it anyway?")))
					return false;
				break;
			}
		}
		
		WaitCursor wait;
		
		Bem().LoadBEM(file, [&](String str, int _pos) {
			progress.SetText(str); 
			progress.SetPos(_pos); 
			return !progress.Canceled();
		}, false);
		
		int id = Bem().hydros.size()-1;
		HydroClass &data = Bem().hydros[id];
		
		data.hd().Report();
		mainSummary.Report(data.hd(), id);
		if (data.hd().Nf < 0)
			return false;
		
		ArrayModel_Add(listLoaded, data.hd().GetCodeStr(), data.hd().name, data.hd().file, data.hd().GetId());
		
		UpdateButtons();

		Upp::Vector<int> ids = ArrayModel_IdsHydro(listLoaded);
		
		mainArrange.Load(Bem().hydros, ids);	
		mainTab.GetItem(mainTab.Find(mainArrange)).Enable(true);	
		mainTab.GetItem(mainTab.Find(mainStiffness)).Enable(mainStiffness.Load(Bem().hydros, ids));
		mainTab.GetItem(mainTab.Find(mainA)).Enable(mainA.Load(Bem(), ids));	
		mainTab.GetItem(mainTab.Find(mainB)).Enable(mainB.Load(Bem(), ids));
		mainTab.GetItem(mainTab.Find(mainK)).Enable(mainK.Load(Bem(), ids));
		// if (Bem().experimental)
			mainTab.GetItem(mainTab.Find(mainAinfw)).Enable(mainAinfw.Load(Bem(), ids));
		mainTab.GetItem(mainTab.Find(mainForceSC)).Enable(mainForceSC.Load(Bem(), ids));
		mainTab.GetItem(mainTab.Find(mainForceFK)).Enable(mainForceFK.Load(Bem(), ids));
		mainTab.GetItem(mainTab.Find(mainForceEX)).Enable(mainForceEX.Load(Bem(), ids));
		mainTab.GetItem(mainTab.Find(mainRAO)).Enable(mainRAO.Load(Bem(), ids));
		mainTab.GetItem(mainTab.Find(mainSetupFOAMM)).Enable(data.hd().IsLoadedB());
		mainTab.GetItem(mainTab.Find(mainQTF)).Enable(mainQTF.Load());
		if (data.hd().IsLoadedStateSpace())
			mainTab.GetItem(mainTab.Find(mainStateSpace)).Enable(true);
		
		mainTab.WhenSet();
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
		return false;
	}
	return true;
}

void MainBEM::OnRemove() {
	OnRemoveSelected(true);
}

void MainBEM::OnRemoveSelected(bool all) {
	bool selected = false;
	for (int row = listLoaded.GetCount()-1; row >= 0; --row) {
		if (all || listLoaded.IsSelected(row)) {
			int id = ArrayModel_IdHydro(listLoaded, row);
			Bem().hydros.Remove(id);
			mainArrange.Remove(row);
			listLoaded.Remove(row);
			selected = true;
		}
	}	// Only one available => directly selected
	if (!selected && listLoaded.GetCount() == 1) {	
		int id = ArrayModel_IdHydro(listLoaded, 0);
		Bem().hydros.Remove(id);
		mainArrange.Remove(0);
		listLoaded.Remove(0);
		selected = true;
	}		
	if (!selected) {
		Exclamation(t_("No model selected"));
		return;
	}
 	mainSummary.Clear();
	for (int i = 0; i < Bem().hydros.size(); ++i)
		mainSummary.Report(Bem().hydros[i].hd(), i);
	
	UpdateButtons();
	
	Upp::Vector<int> ids = ArrayModel_IdsHydro(listLoaded);
	
	mainArrange.Load(Bem().hydros, ids);	
	mainTab.GetItem(mainTab.Find(mainArrange)).Enable(ids.size() > 0);	
	mainTab.GetItem(mainTab.Find(mainStiffness)).Enable(mainStiffness.Load(Bem().hydros, ids));
	mainTab.GetItem(mainTab.Find(mainA)).Enable(mainA.Load(Bem(), ids));	
	mainTab.GetItem(mainTab.Find(mainB)).Enable(mainB.Load(Bem(), ids));
	mainTab.GetItem(mainTab.Find(mainK)).Enable(mainK.Load(Bem(), ids));
	// if (Bem().experimental)
		mainTab.GetItem(mainTab.Find(mainAinfw)).Enable(mainAinfw.Load(Bem(), ids));
	mainTab.GetItem(mainTab.Find(mainForceSC)).Enable(mainForceSC.Load(Bem(), ids));
	mainTab.GetItem(mainTab.Find(mainForceFK)).Enable(mainForceFK.Load(Bem(), ids));
	mainTab.GetItem(mainTab.Find(mainForceEX)).Enable(mainForceEX.Load(Bem(), ids));
	mainTab.GetItem(mainTab.Find(mainRAO)).Enable(mainRAO.Load(Bem(), ids));
	mainTab.GetItem(mainTab.Find(mainQTF)).Enable(mainQTF.Load());
	mainTab.GetItem(mainTab.Find(mainSetupFOAMM)).Enable(ids.size() > 0);
	
	mainTab.WhenSet();
}

void MainBEM::UpdateButtons() {
	int numrow = listLoaded.GetCount();
	int numsel = ArrayCtrlSelectedGetCount(listLoaded);
	menuOpen.butRemove.Enable(numrow > 0);
	menuOpen.butRemoveSelected.Enable(numsel > 0);
	menuOpen.butJoin.Enable(numsel > 1);
	menuOpen.butDuplicate.Enable(numsel == 1);
	menuOpen.butDescription.Enable(numsel == 1 || numrow == 1);
	menuProcess.butSymX.Enable	(numsel == 1 || numrow == 1);
	menuProcess.butSymY.Enable	(numsel == 1 || numrow == 1);
	menuProcess.butKirf.Enable	(numsel == 1 || numrow == 1);
	menuProcess.butA0.Enable  	(numsel == 1 || numrow == 1);
	menuProcess.butAinf.Enable	(numsel == 1 || numrow == 1);
	menuAdvanced.butAinfw.Enable  (numsel == 1 || numrow == 1);
	menuAdvanced.butOgilvie.Enable(numsel == 1 || numrow == 1);
	menuConvert.butConvert.Enable(numsel == 1 || numrow == 1);
}

void MainBEM::OnJoin() {
	Upp::Vector<int> idsjoin, rowsJoin;
	for (int row = listLoaded.GetCount()-1; row >= 0; --row) {
		if (listLoaded.IsSelected(row)) {
			rowsJoin << row;
			int id = ArrayModel_IdHydro(listLoaded, row);
			idsjoin << id;
		}
	}
	if (idsjoin.IsEmpty()) {
		Exclamation(t_("No model selected"));
		return;
	}
	if (idsjoin.size() == 1) {
		Exclamation(t_("Please select more than one model"));
		return;
	}
	try {
		WaitCursor wait;
		Progress progress(t_("Joining selected BEM files..."), 100); 
		
		HydroClass &data = Bem().Join(idsjoin, [&](String str, int _pos) {
			progress.SetText(str); 
			progress.SetPos(_pos); 
			return !progress.Canceled();
		});
		
		mainSummary.Clear();
		mainArrange.Clear();
		
		ArrayModel_Add(listLoaded, data.hd().GetCodeStr(), data.hd().name, data.hd().file, data.hd().GetId());
		ArrayModel_RowsHydroDel(listLoaded, rowsJoin);
	
		Upp::Vector<int> ids = ArrayModel_IdsHydro(listLoaded);
	
		for (int id = 0; id < Bem().hydros.size(); ++id) {
			const Hydro &data = Bem().hydros[id].hd();
			mainSummary.Report(data, id);
			mainArrange.Load(Bem().hydros, ids);
		}
			
		if (Bem().hydros.size() > 0) 
			listLoaded.SetCursor(0);

		mainArrange.Load(Bem().hydros, ids);	
		mainTab.GetItem(mainTab.Find(mainArrange)).Enable(ids.size() > 0);		
		mainTab.GetItem(mainTab.Find(mainStiffness)).Enable(mainStiffness.Load(Bem().hydros, ids));
		mainTab.GetItem(mainTab.Find(mainA)).Enable(mainA.Load(Bem(), ids));	
		mainTab.GetItem(mainTab.Find(mainB)).Enable(mainB.Load(Bem(), ids));
		mainTab.GetItem(mainTab.Find(mainK)).Enable(mainK.Load(Bem(), ids));
		// if (Bem().experimental)
			mainTab.GetItem(mainTab.Find(mainAinfw)).Enable(mainAinfw.Load(Bem(), ids));
		mainTab.GetItem(mainTab.Find(mainForceSC)).Enable(mainForceSC.Load(Bem(), ids));
		mainTab.GetItem(mainTab.Find(mainForceFK)).Enable(mainForceFK.Load(Bem(), ids));
		mainTab.GetItem(mainTab.Find(mainForceEX)).Enable(mainForceEX.Load(Bem(), ids));
		mainTab.GetItem(mainTab.Find(mainRAO)).Enable(mainRAO.Load(Bem(), ids));
		mainTab.GetItem(mainTab.Find(mainQTF)).Enable(mainQTF.Load());
		mainTab.GetItem(mainTab.Find(mainSetupFOAMM)).Enable(true);
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}
		
	UpdateButtons();
}

void MainBEM::OnDuplicate() {
	try {
		int id = GetIdOneSelected();
		if (id < 0) 
			return;

		HydroClass &data = Bem().Duplicate(id);
		
		mainSummary.Clear();
		
		ArrayModel_Add(listLoaded, data.hd().GetCodeStr(), data.hd().name, data.hd().file, data.hd().GetId());
	
		for (int i = 0; i < Bem().hydros.size(); ++i)
			mainSummary.Report(Bem().hydros[i].hd(), i);
		
		Upp::Vector<int> ids = ArrayModel_IdsHydro(listLoaded);
		
		mainTab.GetItem(mainTab.Find(mainForceSC)).Enable(mainForceSC.Load(Bem(), ids));
		mainTab.GetItem(mainTab.Find(mainForceFK)).Enable(mainForceFK.Load(Bem(), ids));
		mainTab.GetItem(mainTab.Find(mainForceEX)).Enable(mainForceEX.Load(Bem(), ids));
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}
}

void MainBEM::OnSymmetrize(bool xAxis) {
	try {
		int id = GetIdOneSelected();
		if (id < 0) 
			return;

		WaitCursor wait;

		Progress progress(t_("Symmetrizing forces and RAOs in selected BEM file..."), 100); 
		
		Bem().Symmetrize(id, xAxis);
		
		mainSummary.Clear();
		for (int i = 0; i < Bem().hydros.size(); ++i)
			mainSummary.Report(Bem().hydros[i].hd(), i);
		
		Upp::Vector<int> ids = ArrayModel_IdsHydro(listLoaded);
		
		mainTab.GetItem(mainTab.Find(mainForceSC)).Enable(mainForceSC.Load(Bem(), ids));
		mainTab.GetItem(mainTab.Find(mainForceFK)).Enable(mainForceFK.Load(Bem(), ids));
		mainTab.GetItem(mainTab.Find(mainForceEX)).Enable(mainForceEX.Load(Bem(), ids));
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}
}

void MainBEM::OnKirfAinf(Hydro::DataToPlot param) {
	try {
		int id = GetIdOneSelected();
		if (id < 0) 
			return;
			
		Progress progress(Format(t_("Calculating %s in selected BEM file..."), Hydro::StrDataToPlot(param)), 100); 
		
		double maxT = Null;
		
		if (param == Hydro::PLOT_K || 
		   ((param == Hydro::PLOT_AINF || param == Hydro::PLOT_AINFW) 
		   	 && !Bem().hydros[id].hd().IsLoadedKirf())) {
		 	maxT = Bem().hydros[id].hd().GetK_IRF_MaxT();
			if (maxT < 0)
				maxT = Bem().maxTimeA;
			else if (Bem().maxTimeA > maxT) {
				if (!PromptYesNo(Format(t_("Defined time for Kirf calculation (%.1f) may be longer than advised (%.1f). Do you wish to used advised time?"), Bem().maxTimeA, maxT)))
					maxT = Bem().maxTimeA;
			} else
				maxT = Bem().maxTimeA;
		}
		   
		WaitCursor wait;

		if (param == Hydro::PLOT_A0)
			Bem().A0(id);
		else if (param == Hydro::PLOT_K)
			Bem().Kirf(id, maxT);
		else if (param == Hydro::PLOT_AINF) {
			if (!Bem().hydros[id].hd().IsLoadedKirf())
				Bem().Kirf(id, maxT);
			Bem().Ainf(id);
		} else if (param == Hydro::PLOT_AINFW) {
			if (!Bem().hydros[id].hd().IsLoadedKirf())
				Bem().Kirf(id, maxT);
			if (!Bem().hydros[id].hd().IsLoadedAinf())
				Bem().Ainf(id);
			Bem().Ainf_w(id);
		}
		
		mainSummary.Clear();
		for (int i = 0; i < Bem().hydros.size(); ++i)
			mainSummary.Report(Bem().hydros[i].hd(), i);
		
		Upp::Vector<int> ids = ArrayModel_IdsHydro(listLoaded);
		
		mainTab.GetItem(mainTab.Find(mainA)).Enable(mainA.Load(Bem(), ids));	
		mainTab.GetItem(mainTab.Find(mainK)).Enable(mainK.Load(Bem(), ids));
		// if (Bem().experimental)
			mainTab.GetItem(mainTab.Find(mainAinfw)).Enable(mainAinfw.Load(Bem(), ids));
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}
}

void MainBEM::OnOgilvie() {
	try {
		int id = GetIdOneSelected();
		if (id < 0) 
			return;
			
		Progress progress(t_("Calculating Ogilvie compliance in selected BEM file..."), 100); 
		
		WaitCursor wait;
		
		double maxT = Null;
		
		Bem().OgilvieCompliance(id, ~menuAdvanced.opZremoval, ~menuAdvanced.opThinremoval, ~menuAdvanced.opDecayingTail);
				
		mainSummary.Clear();
		for (int i = 0; i < Bem().hydros.size(); ++i)
			mainSummary.Report(Bem().hydros[i].hd(), i);
		
		Upp::Vector<int> ids = ArrayModel_IdsHydro(listLoaded);
		
		mainTab.GetItem(mainTab.Find(mainA)).Enable(mainA.Load(Bem(), ids));
		mainTab.GetItem(mainTab.Find(mainB)).Enable(mainB.Load(Bem(), ids));	
		mainTab.GetItem(mainTab.Find(mainK)).Enable(mainK.Load(Bem(), ids));
		// if (Bem().experimental)
			mainTab.GetItem(mainTab.Find(mainAinfw)).Enable(mainAinfw.Load(Bem(), ids));
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}	
}

void MainBEM::OnDescription() {
	int id = GetIdOneSelected();
	if (id < 0) 
		return;

	WithDescription<TopWindow> w;
	CtrlLayout(w);
	w.Title(t_("Enter model name and description"));
	w.butOK 	  << [&] {w.Close();};
	w.name  	  <<= Bem().hydros[id].hd().name;
	w.description <<= Bem().hydros[id].hd().description;
	
	w.Execute();
	
	Bem().hydros[id].hd().name = ~w.name;
	Bem().hydros[id].hd().description = ~w.description;
	
	ArrayModel_Change(listLoaded, Bem().GetHydroId(id), Null, ~w.name, Null);
	
	mainSummary.Clear();
	for (int i = 0; i < Bem().hydros.size(); ++i)
		mainSummary.Report(Bem().hydros[i].hd(), i);
}

int MainBEM::AskQtfHeading(const Hydro &hydro) {
	int numH = hydro.qtfhead.size();
	if (numH <= 1)
		return Null;
	
	WithSaveQTF<TopWindow> dialog;
	CtrlLayout(dialog);
	
	dialog.Title(t_("Please choose the QTF headings to save"));
	
	dialog.swHeadings <<= 0;
	dialog.dropHeadings.Enable(false);
	
	int id0 = Null;
	for (int i = 0; i < numH; ++i) {
		dialog.dropHeadings.Add(hydro.qtfhead[i]);
		if (abs(hydro.qtfhead[i]) < 0.001)
			id0 = i;
	}
	if (!IsNull(id0))
		dialog.dropHeadings.SetIndex(id0);
	else
		dialog.dropHeadings.SetIndex(0);
	
	dialog.swHeadings << [&] {
		dialog.dropHeadings.Enable(dialog.swHeadings == 2);
	};
	
	bool cancel = true;
	dialog.ok		<< [&] {cancel = false;	dialog.Close();};
	dialog.cancel	<< [&] {dialog.Close();}; 
	dialog.Execute();
	if (cancel) 
		throw Exc(t_("Cancelled by user"));
	
	if (dialog.swHeadings == 0)
		return -1;
	else if (dialog.swHeadings == 1)
		return -2;
	else
		return dialog.dropHeadings.GetIndex();
}

bool MainBEM::OnConvert() {
	String file = ~menuConvert.file;
	
	try {
		int id = GetIdOneSelected();
		if (id < 0) 
			return false;
		
		int qtfHeading = Null;
		
		Hydro::BEM_SOFT type;	
		switch (menuConvert.opt) {
		case 0:	type = Hydro::WAMIT;	
				break;
		case 1:	type = Hydro::WAMIT_1_3;	
				qtfHeading = AskQtfHeading(Bem().hydros[id].hd());
				break;
		case 2:	type = Hydro::FAST_WAMIT;	
				qtfHeading = AskQtfHeading(Bem().hydros[id].hd());	
				break;
		case 3:	type = Hydro::BEMROSETTA;	
				break;
		case 4:	type = Hydro::UNKNOWN;		
				break;
		default:throw Exc(t_("Unknown type in OnConvert()"));
		}
		
		Progress progress(t_("Saving BEM files..."), 100); 
		progress.Granularity(1000);
		
		Bem().hydros[id].hd().SaveAs(file, [&](String str, int _pos) {
			if (!IsEmpty(str))
				progress.SetText(str); 
			if (_pos >= 0)
				progress.SetPos(_pos); 
			return !progress.Canceled();}, type, qtfHeading);
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
		return false;
	}
	return true;
}

int MainBEM::GetIdOneSelected() {
	int id = -1;
	for (int row = 0; row < listLoaded.GetCount(); ++row) {
		if (listLoaded.IsSelected(row)) {
			id = ArrayModel_IdHydro(listLoaded, row);
			break;
		}
	}	// Only one available => directly selected
	if (id < 0 && listLoaded.GetCount() == 1)
		id = ArrayModel_IdHydro(listLoaded, 0);
	if (id < 0) {
		Exclamation(t_("No model selected"));
		return -1;
	}
	if (ArrayCtrlSelectedGetCount(listLoaded) > 1) {
		Exclamation(t_("Please select just one model"));
		return -1;
	}
	return id;
}

void MainBEM::Jsonize(JsonIO &json) {
	if (json.IsLoading()) {
		menuPlot.autoFit = Null;
		menuPlot.fromY0 = Null;
		menuPlot.opwT = Null;
		menuPlot.showPoints = Null;
		menuPlot.showNdim = Null;
		menuConvert.opt = Null;
	}
	json
		("menuOpen_file", menuOpen.file)
		("menuConvert_file", menuConvert.file)
		("menuConvert_opt", menuConvert.opt)
		("menuPlot_autoFit", menuPlot.autoFit)
		("menuPlot_fromY0", menuPlot.fromY0)
		("menuPlot_opwT", menuPlot.opwT)
		("menuPlot_showPoints", menuPlot.showPoints)
		("menuPlot_showNdim", menuPlot.showNdim)
		("menuProcess_opZremoval", menuAdvanced.opZremoval)
		("menuProcess_opThinremoval", menuAdvanced.opThinremoval)
		("menuProcess_opDecayingTail", menuAdvanced.opDecayingTail)
		("mainStiffness", mainStiffness)
	;
}

String MainBEM::BEMFile(String fileFolder) const {
	if (DirectoryExists(fileFolder)) {
		int bestipos = INT_MAX;
		for (FindFile ff(AppendFileNameX(fileFolder, "*.*")); ff; ff++) {
			if (ff.IsFile()) {
				String name = ToLower(ff.GetName());
				if (GetFileExt(name) == ".bem")
					return ff.GetPath();
				if (name == "nemoh.cal")
					return ff.GetPath();
				int ipos = Bem().bemFilesExt.Find(GetFileExt(name));
 				if (ipos >= 0 && ipos < bestipos) {
					fileFolder = ff.GetPath();
					bestipos = ipos;	// It takes the file with most probable extension
				}
			} else if (ff.IsFolder()) {
				if (ff.GetName() == "Output") {
					String hamsFolder = AppendFileNameX(fileFolder, "Output", "Wamit_format");
					if (DirectoryExists(hamsFolder))
						return BEMFile(hamsFolder);
				}
			}
		}
	}
	return fileFolder;
}

void MainBEM::LoadDragDrop() {
	GuiLock __;
	
	bool followWithErrors = false;
	for (int i = 0; i < filesToDrop.size(); ++i) {
		String file = BEMFile(filesToDrop[i]);
		menuOpen.file <<= file;
		Status(Format(t_("Loading '%s'"), file));
		if (!OnLoad() && !followWithErrors && filesToDrop.size() - i > 1) {
			if (!PromptYesNo(Format(t_("Do you wish to load the pending %d files?"), filesToDrop.size() - i - 1)))
				return;
			followWithErrors = true;
		}
		ProcessEvents();
	}
	timerDrop.Kill();
}

void MainBEM::DragAndDrop(Point , PasteClip& d) {
	GuiLock __;
	if (IsDragAndDropSource())
		return;
	if (AcceptFiles(d)) {
		filesToDrop = GetFiles(d); 
		timerDrop.Set(0, [=] {LoadDragDrop();});
	}
}

bool MainBEM::Key(dword key, int ) {
	GuiLock __;
	if (key == K_CTRL_V) {
		filesToDrop = GetFiles(Ctrl::Clipboard());
		timerDrop.Set(0, [=] {LoadDragDrop();});
		return true;
	}
	return false;
}

void MainSummary::Init() {
	CtrlLayout(*this);
	array.SetLineCy(EditField::GetStdHeight()).MultiSelect();
	array.WhenBar = [&](Bar &menu) {	ArrayCtrlWhenBar(menu, array);};
}

void MainSummary::Clear() {
	array.Reset();
	array.SetLineCy(EditField::GetStdHeight()).MultiSelect();;
}

void MainSummaryCoeff::Report(const Hydro &data, int id) {
	if (array.GetColumnCount() == 0)
		array.AddColumn("Param");
	if (id >= array.GetColumnCount()-1)
		array.AddColumn(Format("#%d %s", id+1, data.name));
	int row = 0;
	int col = id + 1;
	
	array.Set(row, 0, t_("File"));				array.Set(row++, col, data.file);
	array.Set(row, 0, t_("Name"));				array.Set(row++, col, data.name);
	array.Set(row, 0, t_("Description"));		array.Set(row++, col, data.description);
	array.Set(row, 0, t_("Software"));			array.Set(row++, col, data.GetCodeStr());
	array.Set(row, 0, t_("g [m/s2]"));			array.Set(row++, col, data.S_g());
	array.Set(row, 0, t_("rho [kg/m3]"));		array.Set(row++, col, data.S_rho());
	array.Set(row, 0, t_("h (water depth) [m]"));array.Set(row++,col, data.S_h());
	array.Set(row, 0, t_("length scale [m]"));	array.Set(row++, col, data.S_len());
	
	array.Set(row, 0, t_("#frequencies"));		array.Set(row++, col, Nvl(data.Nf, 0)); 
	if (!data.w.IsEmpty()) {
		array.Set(row, 0, t_("freq_0 [rad/s]"));	array.Set(row++, col, data.w[0]);
	} else {
		array.Set(row, 0, t_("freq_0 [rad/s]"));	array.Set(row++, col, "-");
	}
	if (data.w.size() > 1) {
		array.Set(row, 0, t_("freq_end [rad/s]"));	array.Set(row++, col, data.w[data.w.size()-1]);
		if (data.GetIrregularFreq() < 0) { 
			array.Set(row, 0, t_("freq_delta [rad/s]"));array.Set(row++, col, data.w[1] - data.w[0]);
		} else {
			String strHead;
			for (int i = 0; i < data.w.size(); ++i) {
				if (i > 0)
					strHead << ", ";
				strHead << data.w[i];
			}
			array.Set(row, 0, t_("freq_delta [rad/s]"));
			array.Set(row++, col, Format(t_("Non constant delta (%s)"), strHead));
		}
	} else {
		array.Set(row, 0, t_("freq_end [rad/s]"));	array.Set(row++, col, "-");
		array.Set(row, 0, t_("freq_delta [rad/s]"));array.Set(row++, col, "-");
	}
	
	array.Set(row, 0, t_("#headings"));			array.Set(row++, col, Nvl(data.Nh, 0));
	if (!data.head.IsEmpty()) {
		array.Set(row, 0, t_("head_0 [º]"));	array.Set(row++, col, data.head[0]);
	} else {
		array.Set(row, 0, t_("head_0 [º]"));	array.Set(row++, col, "-");
	}
	if (data.head.size() > 1) {
		array.Set(row, 0, t_("head_end [º]"));	array.Set(row++, col, data.head[data.head.size()-1]);
		if (data.GetIrregularHead() < 0) { 
			array.Set(row, 0, t_("head_delta [º]"));array.Set(row++, col, data.head[1] - data.head[0]);
		} else {
			String strHead;
			for (int i = 0; i < data.head.size(); ++i) {
				if (i > 0)
					strHead << ", ";
				strHead << data.head[i];
			}
			array.Set(row, 0, t_("head_delta [º]"));
			array.Set(row++, col, Format(t_("Non constant delta (%s)"), strHead));
		}
	} else {
		array.Set(row, 0, t_("head_end [º]"));	array.Set(row++, col, "-");
		array.Set(row, 0, t_("head_delta [º]"));array.Set(row++, col, "-");
	}
	
	array.Set(row, 0, t_("A0 available"));		array.Set(row++, col, data.IsLoadedA0()   ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("A∞ available"));		array.Set(row++, col, data.IsLoadedAinf() ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("A available"));		array.Set(row++, col, data.IsLoadedA() 	   ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("B available"));		array.Set(row++, col, data.IsLoadedB() 	   ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("K available"));		array.Set(row++, col, data.IsLoadedC() 	   ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("Inertia available"));	array.Set(row++, col, data.IsLoadedM() 	   ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("Fex available"));		array.Set(row++, col, data.IsLoadedFex()   ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("Fsc available"));		array.Set(row++, col, data.IsLoadedFsc()   ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("Ffk available"));		array.Set(row++, col, data.IsLoadedFfk()   ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("RAO available"));		array.Set(row++, col, data.IsLoadedRAO()   ? t_("Yes") : t_("No"));
	
	array.Set(row, 0, t_("#bodies"));			array.Set(row++, col, data.Nb);
	for (int ib = 0; ib < data.Nb; ++ib) {
		String sib = Format("#%d", ib+1);
		if (data.names.size() > ib) {
			sib += " " + data.names[ib];
			array.Set(row, 0, sib + " " + t_("Name"));		array.Set(row++, col, data.names[ib]);
		} else {
			array.Set(row, 0, sib + " " + t_("Name"));		array.Set(row++, col, "-");
		}
		array.Set(row, 0, sib + " " + t_("#dof"));
		if (data.dof.size() > ib) 
			array.Set(row++, col, data.dof[ib]);
		else 
			array.Set(row++, col, "-");
		
		array.Set(row, 0, sib + " " + t_("Vsub [m3]"));
		if (data.Vo.size() > ib && !IsNull(data.Vo[ib])) 
			array.Set(row++, col, FormatDoubleSize(data.Vo[ib], 10, false));
		else 
			array.Set(row++, col, "-");
		
		array.Set(row, 0, sib + " " + t_("Cg [m]"));
		if (data.cg.size() > 3*ib && !IsNull(data.cg(0, ib))) 
			array.Set(row++, col, Format(t_("%s, %s, %s"),
									FormatDoubleSize(data.cg(0, ib), 10, false),
									FormatDoubleSize(data.cg(1, ib), 10, false),
									FormatDoubleSize(data.cg(2, ib), 10, false)));
		else
			array.Set(row++, col, "-");

		array.Set(row, 0, sib + " " + t_("Cb [m]"));
		if (data.cb.size() > 3*ib && !IsNull(data.cb(0, ib))) 
			array.Set(row++, col, Format(t_("%s, %s, %s"),
									FormatDoubleSize(data.cb(0, ib), 10, false),
									FormatDoubleSize(data.cb(1, ib), 10, false),
									FormatDoubleSize(data.cb(2, ib), 10, false)));
		else
			array.Set(row++, col, "-");
		
		array.Set(row, 0, sib + " " + t_("C0 [m]"));
		if (data.c0.size() > 3*ib && !IsNull(data.c0(0, ib))) 
			array.Set(row++, col, Format(t_("%s, %s, %s"),
									FormatDoubleSize(data.c0(0, ib), 10, false),
									FormatDoubleSize(data.c0(1, ib), 10, false),
									FormatDoubleSize(data.c0(2, ib), 10, false)));
		else
			array.Set(row++, col, "-");
		
		array.Set(row, 0, sib + " " + t_("Water plane area [m2]"));
		if (data.C.size() > ib && data.C[ib].size() > 0) {
			double wPlaneArea = data.C_ndim(ib, 2, 2);
			array.Set(row++, col, FormatDoubleSize(wPlaneArea, 10, false));		
			for (int i = 0; i < 6; ++i) {
				for (int j = 0; j < 6; ++j) {
					if (!Hydro::C_units(i, j).IsEmpty()) {
						array.Set(row, 0, sib + " " + Format(t_("K(%d,%d) [%s]"), i+1, j+1, Hydro::C_units(i, j)));	array.Set(row++, col, FormatDoubleSize(data.C_dim(ib, i, j), 10, false));		
					}
				}
			}
		} else {
			array.Set(row++, col, "-");		
			for (int i = 0; i < 6; ++i) {
				for (int j = 0; j < 6; ++j) {
					if (!Hydro::C_units(i, j).IsEmpty()) {
						array.Set(row, 0, sib + " " + Format(t_("K(%d,%d) [%s]"), i+1, j+1, Hydro::C_units(i, j)));	array.Set(row++, col, "-");		
					}
				}
			}
		}
		array.Set(row, 0, sib + " " + t_("Mass [kg]"));
		if (data.M.size() > ib && data.M[ib].size() > 0) {
			double mass = data.M[ib](0, 0); 
			array.Set(row++, col, FormatDoubleSize(mass, 10, false));	
			for (int i = 3; i < 6; ++i) {
				for (int j = 3; j < 6; ++j) {
					array.Set(row, 0, sib + " " + Format(t_("Inertia(%d,%d) [kg/m2]"), i+1-3, j+1-3));	
					array.Set(row++, col, FormatDoubleSize(data.M[ib](i, j), 10, false));		
				}
			}
		} else {
			array.Set(row++, col, "-");	
			for (int i = 3; i < 6; ++i) {
				for (int j = 3; j < 6; ++j) {
					array.Set(row, 0, sib + " " + Format(t_("Inertia(%d,%d) [kg/m2]"), i+1-3, j+1-3));	
					array.Set(row++, col, "-");		
				}
			}
		}
	}	
}

void MainOutput::Init() {
	CtrlLayout(*this);
	cout.SetReadOnly();
	Print(t_("BEMRosetta\nHydrodynamic coefficients viewer and converter for Boundary Element Method solver formats"));
}

void MainOutput::Print(String str) {
	cout.Append(str);
	cout.ScrollEnd();
}

void MainArrange::Init() {
	CtrlLayout(*this);
}

void MainArrange::Clear() {
	tab.Reset();
	arrangeDOF.Clear();
}

void MainArrange::Load(Upp::Array<HydroClass> &hydros, const Upp::Vector<int> &ids) {
	MainArrange::Clear();
	for (int i = 0; i < ids.size(); ++i) {
		int icase = ids[i];
		ArrangeDOF &arr = arrangeDOF.Add();
		arr.Init(hydros[icase].hd());
		tab.Add(arr.SizePos(), Format("[%s] %s", hydros[icase].hd().GetCodeStr(), hydros[icase].hd().name));
	}
}

void MainArrange::Remove(int c) {
	tab.Remove(c);
	arrangeDOF.Remove(c);
}

	
void MainQTF::Init() {
	CtrlLayout(*this);
	
	listCases.Reset();
	listCases.SetLineCy(EditField::GetStdHeight());
	listCases.AddColumn(t_("# Body"), 15);
	listCases.AddColumn(t_("Heading 1"), 20);
	listCases.AddColumn(t_("Heading 2"), 20);
	
	opShow.Clear();
	opShow.Add(MAGNITUDE, t_("Magnitude")).Add(PHASE, t_("Phase")).Add(REAL, t_("Real")).Add(IMAGINARY, t_("Imaginary"));
	opShow.SetIndex(0);

	opDOF.Clear();
	for (int i = 0; i < 6; ++i)
		opDOF.Add(i, InitCaps(Hydro::StrDOF_base(i)));
	opDOF.SetIndex(0);
	
	opQTF.Clear();
	opQTF.Add(FSUM, t_("Sum")).Add(FDIFFERENCE, t_("Difference"));
	opQTF.SetIndex(0);
	
	opShow << [&] {listCases.WhenSel();};
	opDOF  << [&] {listCases.WhenSel();};
	opQTF  << [&] {listCases.WhenSel();};
	
	listCases.WhenSel = [&] {
		if (idHydro < 0)
			return;
		
		int id = listCases.GetCursor();
		if (id < 0)
			return;

		try {
			WaitCursor wait;
			
			MainBEM &mbm = GetDefinedParent<MainBEM>(this);
			
			bool ndim = mbm.menuPlot.showNdim;
			bool show_w = mbm.menuPlot.opwT == 0;
		
			const Hydro &hd = Bem().hydros[idHydro].hd();
			
			int ib = int(listCases.Get(id, 0))-1;
			int ih1 = hd.GetQTFHeadId(listCases.Get(id, 1));
			int ih2 = hd.GetQTFHeadId(listCases.Get(id, 2));
			int idof = opDOF.GetIndex();
			int qtfNf = hd.qtfw.size();
			
			listQTF.Reset();
			listQTF.SetLineCy(EditField::GetStdHeight()).HeaderObject().Absolute();
			listQTF.MultiSelect();
			listQTF.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, listQTF);};
	
			listQTF.AddColumn(show_w ? t_("ω [rad/s]") : t_("T [s]"), 60);
			for (int c = 0; c < qtfNf; ++c)
				listQTF.AddColumn(FormatDoubleSize(show_w ? hd.qtfw[c] : hd.qtfT[c], 8), 90);
			for (int r = 0; r < qtfNf; ++r)
				listQTF.Add(FormatDoubleSize(show_w ? hd.qtfw[r] : hd.qtfT[r], 8));
			
			const Upp::Array<Hydro::QTF> &qtfList = opQTF.GetData() == FSUM ? hd.qtfsum : hd.qtfdif;
			
			double mn = DBL_MAX, mx = DBL_MIN;
			for (int ifr1 = 0; ifr1 < qtfNf; ++ifr1) {
				for (int ifr2 = 0; ifr2 < qtfNf; ++ifr2) {
					int idq = 0;
					if ((idq = hd.GetQTFId(idq, qtfList, hd.qtfCases, ib, ih1, ih2, ifr1, ifr2)) < 0)
						continue;
					double val;
					switch(int(opShow.GetData())) {
					case MAGNITUDE:
						val = hd.F_(ndim, qtfList[idq].fma[idof], idof);	break;
					case REAL:
						val = hd.F_(ndim, qtfList[idq].fre[idof], idof);	break;
					case IMAGINARY:
						val = hd.F_(ndim, qtfList[idq].fim[idof], idof);	break;
					}
					mn = min(mn, val);
					mx = max(mx, val);
				}
			}
			for (int ifr1 = 0; ifr1 < qtfNf; ++ifr1) {
				for (int ifr2 = 0; ifr2 < qtfNf; ++ifr2) {
					int idq = 0;
					if ((idq = hd.GetQTFId(idq, qtfList, hd.qtfCases, ib, ih1, ih2, ifr1, ifr2)) < 0)
						listQTF.Set(ifr2, 1+ifr1, "-");
					else {
						if (PHASE == opShow.GetData()) 
							listQTF.Set(ifr2, 1+ifr1, FormatDoubleSize(qtfList[idq].fph[idof], 10, false));
						else {
							double val;
							switch(int(opShow.GetData())) {
							case MAGNITUDE:
								val = hd.F_(ndim, qtfList[idq].fma[idof], idof);	break;
							case REAL:
								val = hd.F_(ndim, qtfList[idq].fre[idof], idof);	break;
							case IMAGINARY:
								val = hd.F_(ndim, qtfList[idq].fim[idof], idof);	break;
							}
							
							::Color backColor = GetRainbowColor((val - mn)/(mx - mn), White(), LtBlue(), 0);
							::Color color = Black();
							if (Grayscale(backColor) < 150)
								color = White();
							
							String str = FormatDoubleSize(val, 10, false);
							
							listQTF.Set(ifr2, 1+ifr1, AttrText(str).Center().Ink(color).Paper(backColor));
						}
					}
				}
			}
		} catch (Exc e) {
			Exclamation(DeQtfLf(e));
		}
	};
}
	
bool MainQTF::Load() {
	listCases.Clear();
	listQTF.Reset();
	
	try {
		MainBEM &mbm = GetDefinedParent<MainBEM>(this);
		
		idHydro = -1;
		for (int row = 0; row < mbm.listLoaded.GetCount(); ++row) {
			if (mbm.listLoaded.IsSelected(row)) {
				idHydro = ArrayModel_IdHydro(mbm.listLoaded, row);
				break;
			}
		}	// Only one available => directly selected
		if (idHydro < 0 && mbm.listLoaded.GetCount() == 1)
			idHydro = ArrayModel_IdHydro(mbm.listLoaded, 0);
		if (idHydro < 0) 
			return false;
		
		if (ArrayCtrlSelectedGetCount(mbm.listLoaded) > 1) 
			return false;
		
		const Hydro &hd = Bem().hydros[idHydro].hd();
		
		if (hd.qtfsum.IsEmpty() && hd.qtfdif.IsEmpty())
			return false;
		
		for (int i = 0; i < hd.qtfCases.ib.size(); ++i) 
			listCases.Add(hd.qtfCases.ib[i]+1, hd.qtfhead[hd.qtfCases.ih1[i]], hd.qtfhead[hd.qtfCases.ih2[i]]);

		if (listCases.GetCount() > 0)
			listCases.SetCursor(0);
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}				
	return true;		
}

void MainBEMW::Init(MainBEM &_bem, const Image &icon, const Image &largeIcon) {
	LoadFromJson(bem, StoreAsJson(_bem));
	bem.Init();
	Add(bem.SizePos());
	Title(t_("BEMRosetta BEM Coefficients Processing")).Sizeable().Zoomable().Icon(icon, largeIcon);
}