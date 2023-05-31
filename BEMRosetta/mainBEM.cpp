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


void MainBEM::Init() {
	CtrlLayout(*this);
	
	CtrlLayout(menuOpen);
	menuOpen.file.WhenChange = THISBACK(OnLoad);
	menuOpen.file.BrowseRightWidth(40).UseOpenFolder().BrowseOpenFolderWidth(10);
	menuOpen.butLoad << [&] {menuOpen.file.DoGo();};
	
	ArrayModel_Init(listLoaded).MultiSelect();
	listLoaded.WhenSel = [&] {
		OnMenuAdvancedArraySel();
		menuFOAMM.OnCursor();
		//mainTab.GetItem(mainTab.Find(mainQTF)).Enable(mainQTF.Load());
		if (mainTab.Find(mainQTF) == mainTab.Get())
			mainQTF.Load();
		UpdateButtons();
	};
	listLoaded.WhenBar = [&](Bar &menu) {
		listLoaded.StdBar(menu);
		menu.Add(listLoaded.GetCount() > 0, t_("Open file folder"), Null, [&]{
			LaunchWebBrowser(GetFileFolder(ArrayModel_GetFileName(listLoaded)));}).Help(t_("Opens file explorer in the file folder"));
		menu.Add(listLoaded.GetCount() > 0, t_("Remove"), Null, [&]{
			OnRemoveSelected(false);}).Help(t_("Remove model"));	
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
	menuOpen.butExport <<= THISBACK(OnConvert);
	menuOpen.butExport.Tip(t_("Exports data file"));
	for (int i = 0; i < Hydro::GetBemStrCount(); ++i)
		if (Hydro::bemCanSave[i])
			menuOpen.dropExport.Add(Hydro::GetBemStr(static_cast<Hydro::BEM_FMT>(i)));
	menuOpen.dropExport.SetIndex(dropExportId);

	CtrlLayout(menuProcess);
	menuProcess.butSymX.Disable();	
	menuProcess.butSymX <<= THISBACK1(OnSymmetrizeForces, true);
	menuProcess.butSymY.Disable();
	menuProcess.butSymY <<= THISBACK1(OnSymmetrizeForces, false);
	menuProcess.butA0.Disable();	
	menuProcess.butA0 <<= THISBACK1(OnKirfAinf, Hydro::PLOT_A0);
	menuProcess.butAinf.Disable();	
	menuProcess.butAinf <<= THISBACK1(OnKirfAinf, Hydro::PLOT_AINF);
	menuProcess.butKirf.Disable();	
	menuProcess.butKirf <<= THISBACK1(OnKirfAinf, Hydro::PLOT_KIRF);
	menuProcess.butRAO.Disable();	
	menuProcess.butRAO <<= THISBACK(OnRAO);
	menuProcess.butSymmetrize <<= THISBACK(OnSymmetrize);
	
	menuProcess.butABForces << THISBACK(OnABForces);
	menuProcess.butQTF << THISBACK(OnQTF);

	menuProcess.butABForcesZero << THISBACK(OnABForcesZero);
	menuProcess.butQTFZero << THISBACK(OnQTFZero);
	
	menuProcess.opFill.Tip(t_("Fills with zeroes or with interpolated values"));

	auto DropDOF = [&](DropList &drop1, DropList &drop2)->bool {
		if (drop1.GetIndex() < 1 || drop2.GetIndex() < 1)
			return false;
		if (drop1.GetIndex() == drop2.GetIndex())
			return false;
		return true;
	};
	
	menuProcess.dropDOF1.WhenAction = [&] {
		menuProcess.butSwapDOF.Show(DropDOF(menuProcess.dropDOF1, menuProcess.dropDOF2));
	};
	menuProcess.dropDOF1.Add(t_("Choose"));
		for (int i = 0; i < 6; ++i)
			menuProcess.dropDOF1.Add(BEM::StrDOF(i));
	menuProcess.dropDOF1.SetIndex(0);
	
	menuProcess.dropDOF2.WhenAction = [&] {
		menuProcess.butSwapDOF.Show(DropDOF(menuProcess.dropDOF1, menuProcess.dropDOF2));
	};
	menuProcess.dropDOF2.Add(t_("Choose"));
		for (int i = 0; i < 6; ++i)
			menuProcess.dropDOF2.Add(BEM::StrDOF(i));
	menuProcess.dropDOF2.SetIndex(0);
	
	menuProcess.butSwapDOF << THISBACK(OnSwapDOF);
	menuProcess.butSwapDOF.Hide();
	
	if (IsNull(menuProcess.opFill))
		menuProcess.opFill <<= 0;
	if (IsNull(menuProcess.maxFreq))
		menuProcess.maxFreq <<= 150;
	
	menuProcess.butQTF_MD << THISBACK(OnQTF_MD);
	
	CtrlLayout(menuProcess2);
	
	auto DropChecked = [&](DropGrid &drop)->bool {
		for (int i = 0; i < drop.GetCount(); ++i) {
			if (drop.GetList().Get(i, 0) == true) 
				return true;
		}
		return false;
	};
	
	menuProcess2.factor <<= 0;
	menuProcess2.dropDOF.AddColumn("", 20);
	menuProcess2.dropDOF.AddColumn("", 50);
	menuProcess2.dropDOF.GetList().GetColumn(0).Option();
	menuProcess2.dropDOF.Width(100);
	menuProcess2.dropDOF.GetList().Sorting(false);
	menuProcess2.dropDOF.OnFocus = [&] {
		menuProcess2.dropDOF.DropGrid::GotFocus();
		bool show = DropChecked(menuProcess2.dropDOF) && 
					(menuProcess2.opA || menuProcess2.opAd || 
					 menuProcess2.opB || menuProcess2.opBd || 
					 menuProcess2.opF || menuProcess2.opMD ||
					 menuProcess2.opQTF);
		menuProcess2.butResetDOF.Show(show);
		menuProcess2.butResetDOF0.Show(show);
	};
	menuProcess2.opA.WhenAction = menuProcess2.opAd.WhenAction =
		menuProcess2.opB.WhenAction = menuProcess2.opBd.WhenAction =
		menuProcess2.opF.WhenAction = menuProcess2.opMD.WhenAction =
		menuProcess2.opQTF.WhenAction = menuProcess2.dropDOF.OnFocus;
	
	
	menuProcess2.butResetDOF << THISBACK1(OnMultiplyDOF, false);
	menuProcess2.butResetDOF.Hide();
	menuProcess2.butResetDOF0 << THISBACK1(OnMultiplyDOF, true);
	menuProcess2.butResetDOF0.Hide();
	
	menuProcess2.dropFreq.Tip(t_("Selects the frequencies to be removed from the first order hydrodynamic coefficients and mean drift."));		
	menuProcess2.dropFreq.AddColumn("", 20);
	menuProcess2.dropFreq.AddColumn("", 50);
	menuProcess2.dropFreq.GetList().GetColumn(0).Option();
	menuProcess2.dropFreq.Width(100);
	menuProcess2.dropFreq.GetList().Sorting(false);
	menuProcess2.dropFreq.OnFocus = [&] {
		menuProcess2.dropFreq.DropGrid::GotFocus();
		menuProcess2.butRemoveFreq.Show(DropChecked(menuProcess2.dropFreq) || 
										DropChecked(menuProcess2.dropFreqQTF));
	};
	menuProcess2.dropFreq.Tip(t_("Selects the frequencies to be removed from the QTF."));
	menuProcess2.dropFreqQTF.AddColumn("", 20);
	menuProcess2.dropFreqQTF.AddColumn("", 50);
	menuProcess2.dropFreqQTF.GetList().GetColumn(0).Option();
	menuProcess2.dropFreqQTF.Width(100);
	menuProcess2.dropFreqQTF.GetList().Sorting(false);
	menuProcess2.butRemoveFreq << THISBACK(OnDeleteHeadingsFrequencies);
	menuProcess2.dropFreqQTF.OnFocus = [&] {
		menuProcess2.dropFreqQTF.DropGrid::GotFocus();
		menuProcess2.butRemoveFreq.Show(DropChecked(menuProcess2.dropFreq) || 
										DropChecked(menuProcess2.dropFreqQTF));
	};
	menuProcess2.butRemoveFreq.Hide();
	
	menuProcess2.dropHead.AddColumn("", 20);
	menuProcess2.dropHead.AddColumn(t_("Heading [º]"), 50);
	menuProcess2.dropHead.GetList().GetColumn(0).Option();
	menuProcess2.dropHead.Width(100);
	menuProcess2.dropHead.GetList().Sorting(false);
	menuProcess2.dropHead.OnFocus = [&] {
		menuProcess2.dropHead.DropGrid::GotFocus();
		menuProcess2.butRemoveHead.Show(DropChecked(menuProcess2.dropHead) || 
										DropChecked(menuProcess2.dropHeadMD) || 
										DropChecked(menuProcess2.dropHeadQTF));
	};
	menuProcess2.dropHeadMD.AddColumn("", 10);
	menuProcess2.dropHeadMD.AddColumn(t_("Heading [º]"), 50);
	menuProcess2.dropHeadMD.GetList().GetColumn(0).Option();
	menuProcess2.dropHeadMD.Width(150);
	menuProcess2.dropHeadMD.GetList().Sorting(false);
	menuProcess2.dropHeadMD.OnFocus = [&] {
		menuProcess2.dropHeadMD.DropGrid::GotFocus();
		menuProcess2.butRemoveHead.Show(DropChecked(menuProcess2.dropHead) || 
										DropChecked(menuProcess2.dropHeadMD) || 
										DropChecked(menuProcess2.dropHeadQTF));
	};
	menuProcess2.dropHeadQTF.AddColumn("", 10);
	menuProcess2.dropHeadQTF.AddColumn(t_("Heading [º]"), 50);
	menuProcess2.dropHeadQTF.GetList().GetColumn(0).Option();
	menuProcess2.dropHeadQTF.Width(150);
	menuProcess2.dropHeadQTF.GetList().Sorting(false);
	menuProcess2.dropHeadQTF.OnFocus = [&] {
		menuProcess2.dropHeadQTF.DropGrid::GotFocus();
		menuProcess2.butRemoveHead.Show(DropChecked(menuProcess2.dropHead) || 
										DropChecked(menuProcess2.dropHeadMD) || 
										DropChecked(menuProcess2.dropHeadQTF));
	};
	menuProcess2.butRemoveHead << THISBACK(OnDeleteHeadingsFrequencies);
	menuProcess2.butRemoveHead.Hide();
	
	menuProcess2.dropForce.AddColumn("", 20);
	menuProcess2.dropForce.AddColumn(t_("Force type"), 50);
	menuProcess2.dropForce.GetList().GetColumn(0).Option();
	menuProcess2.dropForce.Width(120);
	menuProcess2.dropForce.GetList().Sorting(false);
	menuProcess2.dropForce.OnFocus = [&] {
		menuProcess2.dropForce.DropGrid::GotFocus();
		menuProcess2.butRemoveForces.Show(DropChecked(menuProcess2.dropForce) || 
										  DropChecked(menuProcess2.dropForceMD) ||
										  DropChecked(menuProcess2.dropForceQTF));
	};
	menuProcess2.dropForceMD.AddColumn("", 20);
	menuProcess2.dropForceMD.AddColumn(t_("Force type"), 50);
	menuProcess2.dropForceMD.GetList().GetColumn(0).Option();
	menuProcess2.dropForceMD.Width(100);
	menuProcess2.dropForceMD.GetList().Sorting(false);
	menuProcess2.dropForceMD.OnFocus = [&] {
		menuProcess2.dropForceMD.DropGrid::GotFocus();
		menuProcess2.butRemoveForces.Show(DropChecked(menuProcess2.dropForce) || 
										  DropChecked(menuProcess2.dropForceMD) ||
										  DropChecked(menuProcess2.dropForceQTF));
	};
	menuProcess2.dropForceQTF.AddColumn("", 20);
	menuProcess2.dropForceQTF.AddColumn(t_("Force type"), 50);
	menuProcess2.dropForceQTF.GetList().GetColumn(0).Option();
	menuProcess2.dropForceQTF.Width(100);
	menuProcess2.dropForceQTF.GetList().Sorting(false);
	menuProcess2.dropForceQTF.OnFocus = [&] {
		menuProcess2.dropForceQTF.DropGrid::GotFocus();
		menuProcess2.butRemoveForces.Show(DropChecked(menuProcess2.dropForce) || 
										  DropChecked(menuProcess2.dropForceMD) ||
										  DropChecked(menuProcess2.dropForceQTF));
	};
	menuProcess2.butRemoveForces << THISBACK(OnResetForces);
	menuProcess2.butRemoveForces.Hide();
	
	CtrlLayout(menuAdvanced);
	menuAdvanced.butAinfw <<= THISBACK1(OnKirfAinf, Hydro::PLOT_AINFW);
	menuAdvanced.butOgilvie <<= THISBACK(OnOgilvie);
	menuAdvanced.butAinfw.Disable();	
	menuAdvanced.butOgilvie.Disable();
	menuAdvanced.opDecayingTail.Disable();
	menuAdvanced.opThinremoval.Disable();
	menuAdvanced.opZremoval.Disable();
	menuAdvanced.butConvergence <<= THISBACK(OnConvergence);
	menuAdvanced.butConvergence.Disable();
	menuAdvanced.butAverage <<= THISBACK(OnAverage);
	menuAdvanced.butAverage.Disable();
	//menuAdvanced.opHaskind.Disable();
	menuAdvanced.opHaskind.Hide();		// Not ready
	menuAdvanced.butUpdateCrot << THISBACK(OnUpdateCrot);

	OnOpt();
	
	CtrlLayout(menuPlot);
	menuPlot.butZoomToFit << [&] {GetSelScatter().ZoomToFit(true, true);};
	menuPlot.autoFit	  << [&] {
		LoadSelTab(Bem());
		menuPlot.fromY0.Enable(~menuPlot.autoFit);
	};
	menuPlot.fromY0 	<< [&] {LoadSelTab(Bem());};
	menuPlot.opwT 		<< [&] {LoadSelTab(Bem());};
	menuPlot.opMP 		<< [&] {LoadSelTab(Bem());};
	menuPlot.showPoints << [&] {LoadSelTab(Bem());};
	menuPlot.showNdim 	<< [&] {LoadSelTab(Bem());};
	
	menuPlot.head1st.NoHeader().MultiSelect();
	menuPlot.head1st.AddColumn("");
	menuPlot.headMD.NoHeader().MultiSelect();
	menuPlot.headMD.AddColumn("");
	menuPlot.headMD.AddColumn("");
	//menuPlot.headQTF.NoHeader().MultiSelect();
	//menuPlot.headQTF.AddColumn("");
	
	OnOpt();
	
	menuFOAMM.Init(*this, mainSetupFOAMM);
	
	OnOpt();
		
	menuTab.Add(menuOpen.SizePos(), 	t_("Load"));
	menuTab.Add(menuPlot.SizePos(), 	t_("Plot")).Disable();
	menuTab.Add(menuProcess.SizePos(), 	t_("Process")).Disable();
	menuTab.Add(menuProcess2.SizePos(), t_("Remove & Mult")).Disable();
	menuTab.Add(menuAdvanced.SizePos(), t_("Advanced")).Disable();
	menuTab.Add(menuFOAMM.SizePos(), 	t_("FOAMM State Space")).Disable();
	
	menuTab.WhenSet = [&] {
		LOGTAB(menuTab);
		bool setupfoamm = false;
		if (menuTab.IsAt(menuFOAMM)) {
			setupfoamm = true;
			if (!FileExists(Bem().foammPath))
				Status(t_("FOAMM not found. Please set FOAMM path in Options"), 10000);	
		} else if (menuTab.IsAt(menuPlot) || menuTab.IsAt(menuOpen) || 
				   menuTab.IsAt(menuProcess) || menuTab.IsAt(menuProcess2) ||
				   menuTab.IsAt(menuAdvanced)) 
			setupfoamm = true;
		
		if (!setupfoamm) 
			mainTab.Set(0);
		
		if (menuTab.IsAt(menuFOAMM)) 
			mainTab.Set(mainSetupFOAMM);
		
		ShowMenuPlotItems();
	};
	
	mainTab.WhenSet = [&] {
		LOGTAB(mainTab);
		UVector<int> ids = ArrayModel_IdsHydro(listLoaded);
		bool plot = true, convertProcess = true, ismenuFOAMM = false;
		int is = -1;			// 0: 1st, 1: QTF, 2: MD

		if (ids.IsEmpty())
			plot = convertProcess = false;
		else if (mainTab.IsAt(mainMatrixK)) 
			mainMatrixK.Load(Bem().hydros, ids, ~menuPlot.showNdim);
		else if (mainTab.IsAt(mainMatrixA))
			mainMatrixA.Load(Bem().hydros, ids, ~menuPlot.showNdim);
		else if (mainTab.IsAt(mainMatrixM)) {
			plot = false;
			mainMatrixM.Load(Bem().hydros, ids, false);
		} else if (mainTab.IsAt(mainMatrixDlin)) {
			plot = false;
			mainMatrixDlin.Load(Bem().hydros, ids, false);
		} else if (mainTab.IsAt(mainA))
			mainA.Load(Bem(), ids);
		else if (mainTab.IsAt(mainB))
			mainB.Load(Bem(), ids);
		else if (mainTab.IsAt(mainMD)) {
			mainMD.Load(Bem(), ids, menuPlot.headMD.GetCursor());
			is = 2;
		} else if (mainTab.IsAt(mainK))
			mainK.Load(Bem(), ids);
		else if (mainTab.IsAt(mainAinfw))
			mainAinfw.Load(Bem(), ids);
		else if (mainTab.IsAt(mainForceSC)) {
			is = 0;
			mainForceSC.Load(Bem(), ids, menuPlot.head1st.GetCursor());
		} else if (mainTab.IsAt(mainForceFK)) {
			is = 0;
			mainForceFK.Load(Bem(), ids, menuPlot.head1st.GetCursor());
		} else if (mainTab.IsAt(mainForceEX)) {
			is = 0;
			mainForceEX.Load(Bem(), ids, menuPlot.head1st.GetCursor());
		} else if (mainTab.IsAt(mainRAO))
			mainRAO.Load(Bem(), ids, menuPlot.head1st.GetCursor());
		else if (mainTab.IsAt(mainStateSpace)) {
			mainStateSpace.Load(Bem(), ids);
			ismenuFOAMM = true;
		} else if (mainTab.IsAt(mainSetupFOAMM)) {
			menuTab.Set(menuFOAMM);
			ismenuFOAMM = true;
			if (mainSetupFOAMM.arrayCases.GetCount() == 0 && ids.size() == 1) 
				listLoaded.SetCursor(0);
			menuFOAMM.OnCursor();
		} else if (mainTab.IsAt(mainQTF)) {
			mainQTF.Load();
			is = 1;
		} else if (menuTab.IsAt(menuFOAMM)) 
			;
		else 
			plot = false;
		
		if (is != 1)
			mainQTF.Unload();
		
		menuPlot.labHead1st.Show(is == 0);	menuPlot.head1st.Show(is == 0);
		menuPlot.labHeadQTF.Show(is == 1);	menuPlot.headQTF.Show(is == 1);
		menuPlot.labHeadMD.Show(is == 2);	menuPlot.headMD.Show(is == 2);
			
		TabCtrl::Item& tabMenuPlot = menuTab.GetItem(menuTab.Find(menuPlot));
		tabMenuPlot.Enable(plot);
		TabCtrl::Item& tabMenuProcess = menuTab.GetItem(menuTab.Find(menuProcess));
		tabMenuProcess.Enable(convertProcess);
		TabCtrl::Item& tabMenuProcess2 = menuTab.GetItem(menuTab.Find(menuProcess2));
		tabMenuProcess2.Enable(convertProcess);
		TabCtrl::Item& tabMenuAdvanced = menuTab.GetItem(menuTab.Find(menuAdvanced));
		tabMenuAdvanced.Enable(convertProcess);

		if (plot) {
			tabMenuPlot.Text(t_("Plot"));
			tabMenuProcess.Text(t_("Process"));
			tabMenuProcess2.Text(t_("Remove"));
			tabMenuAdvanced.Text(t_("Advanced"));
		} else {
			tabMenuPlot.Text("");
			tabMenuProcess.Text("");
			tabMenuProcess2.Text("");
			tabMenuAdvanced.Text("");
		}
		
		if (convertProcess) {
			tabMenuProcess.Text(t_("Process"));
			tabMenuProcess2.Text(t_("Remove"));
			tabMenuAdvanced.Text(t_("Advanced"));
		} else {
			tabMenuProcess.Text("");
			tabMenuProcess2.Text("");
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
	
	mainMatrixK.Init(Hydro::MAT_K);
	mainTab.Add(mainMatrixK.SizePos(), t_("K")).Disable();
		
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

	mainMatrixA.Init(Hydro::MAT_A);
	mainTab.Add(mainMatrixA.SizePos(), t_("A∞")).Disable();
	
	mainAinfw.Init(Hydro::DATA_AINFW);
	mainTab.Add(mainAinfw.SizePos(), t_("A∞(ω)")).Disable();

	mainMatrixM.Init(Hydro::MAT_M);
	mainTab.Add(mainMatrixM.SizePos(), t_("M")).Disable();
	
	mainMatrixDlin.Init(Hydro::MAT_DAMP_LIN);
	mainTab.Add(mainMatrixDlin.SizePos(), t_("Lin. Damp.")).Disable();
	
	mainMD.Init(Hydro::DATA_MD);
	mainTab.Add(mainMD.SizePos(), t_("Mean Drift")).Disable();
		
	mainQTF.Init(*this);
	mainTab.Add(mainQTF.SizePos(), t_("QTF")).Disable();
	
	mainSetupFOAMM.Init();
	mainTab.Add(mainSetupFOAMM.SizePos(), t_("FOAMM")).Disable();

	mainStateSpace.Init();
	mainTab.Add(mainStateSpace.SizePos(), t_("State Space")).Disable();
	
	UpdateButtons();
	saveFolder = GetDesktopFolder();
}

void MainBEM::ShowMenuPlotItems() {
	menuPlot.showNdim.Enable();

	bool show = true, showwT = true, showComplex = false, showDim = true;
	if (mainTab.IsAt(mainSetupFOAMM)) {
		showwT = false;
		showComplex = true;
	} else if (mainTab.IsAt(mainK)) 
		showwT = false;
	else if (mainTab.IsAt(mainQTF)) {
		show = true;
		showComplex = true;
	} else if (mainTab.IsAt(mainForceSC) || mainTab.IsAt(mainForceFK) || mainTab.IsAt(mainForceEX) || 
			   mainTab.IsAt(mainRAO) || mainTab.IsAt(mainStateSpace)) {
		showComplex = true;
		if (mainTab.IsAt(mainRAO))
			showDim = false;
	}
		
	menuPlot.showNdim.Enable(showDim);
	menuPlot.opwT.Enable(showwT);
	menuPlot.opMP.Enable(showComplex);
	menuPlot.butZoomToFit.Enable(show);
	menuPlot.autoFit.Enable(show);
	menuPlot.fromY0.Enable(show);
	menuPlot.showPoints.Enable(show);
}

void MainBEM::OnMenuAdvancedArraySel() {
	int id = ArrayModel_IdHydro(listLoaded);
	if (id < 0)
		return;
	
	Hydro &data = Bem().hydros[id].hd();
	menuAdvanced.x_0.Enable(data.c0.size() == 3);	
	menuAdvanced.y_0.Enable(data.c0.size() == 3);	
	menuAdvanced.z_0.Enable(data.c0.size() == 3);	
		
	if (data.c0.size() == 3) {
		menuAdvanced.x_0 <<= data.c0(0);
		menuAdvanced.y_0 <<= data.c0(1);
		menuAdvanced.z_0 <<= data.c0(2);
	} else {
		menuAdvanced.x_0 <<= Null;
		menuAdvanced.y_0 <<= Null;
		menuAdvanced.z_0 <<= Null;
	}
}

void MainBEM::InitSerialize(bool ret) {
	if (!ret || IsNull(menuPlot.autoFit)) 
		menuPlot.autoFit = true;
	
	if (!ret || IsNull(menuPlot.fromY0)) 
		menuPlot.fromY0 = false;
	
	if (!ret || IsNull(menuPlot.opwT)) 
		menuPlot.opwT = 0;

	if (!ret || IsNull(menuPlot.opMP)) 
		menuPlot.opMP = 0;
	
	if (!ret || IsNull(menuPlot.showPoints)) 
		menuPlot.showPoints = true;
	
	if (!ret || IsNull(menuPlot.showNdim)) 
		menuPlot.showNdim = false;
}

void MainBEM::LoadSelTab(BEM &bem) {
	UVector<int> ids = ArrayModel_IdsHydro(listLoaded);
	int id = mainTab.Get();
	if (id == mainTab.Find(mainStateSpace))
		mainStateSpace.Load(bem, ids);
	else if (id == mainTab.Find(mainMatrixK))
		mainMatrixK.Load(bem.hydros, ids, ~menuPlot.showNdim);
	else if (id == mainTab.Find(mainMatrixA))
		mainMatrixA.Load(bem.hydros, ids, ~menuPlot.showNdim);
	else if (id == mainTab.Find(mainMatrixM))
		mainMatrixM.Load(bem.hydros, ids, false);
	else if (id == mainTab.Find(mainMatrixDlin))
		mainMatrixDlin.Load(bem.hydros, ids, false);
	else if (id == mainTab.Find(mainSummary))
		;
	else if (id == mainTab.Find(mainSetupFOAMM))
		mainSetupFOAMM.Load();
	else if (id == mainTab.Find(mainQTF))
		mainQTF.Load();
	else 
		GetSelABForce().Load(bem, ids);
	
	UpdateButtons();
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

	menuOpen.file.Type(Format(t_("All supported BEM files (%s)"), Bem().bstFilesExt/*bemFilesExt*/), Bem().bemFilesAst);
	menuOpen.file.AllFilesType();
	String extOpen = ToLower(GetFileExt(menuOpen.file.GetData().ToString()));
	if (extOpen.IsEmpty())
		menuOpen.file.ActiveType(0);
	else if (Bem().bstFilesExt/*bemFilesExt*/.Find(extOpen) >= 0)
		menuOpen.file.ActiveType(0);
	else
		menuOpen.file.ActiveType(1);
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
			if (ForceExt(Bem().hydros[i].hd().file, ".") == ForceExt(file, ".") &&
				(Bem().GetBEMExtSet(file) < 0 || 
				 Bem().GetBEMExtSet(file) == Bem().GetBEMExtSet(Bem().hydros[i].hd().file))) {
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
		Hydro &data = Bem().hydros[id].hd();
		
		data.Report();
		mainSummary.Report(data, id);
		if (data.Nf < 0)
			return false;
		
		ArrayModel_Add(listLoaded, data.GetCodeStr(), data.name, data.file, data.GetId());
		
		UpdateButtons();

		AfterBEM();
		
		// Sets headings to 0
		{
			int rowid = Null;
			double distanceTo0 = std::numeric_limits<double>::max();
			for (int r = 0; r < menuPlot.head1st.GetCount(); ++r) {
				double dist = abs(double(menuPlot.head1st.Get(r, 0)));
				if (distanceTo0 > dist) {
					rowid = r;
					distanceTo0 = dist;
				}
			}
			if (!IsNull(rowid))
				menuPlot.head1st.SetCursor(rowid);
		}
		{
			int rowid = Null;
			double distanceTo0 = std::numeric_limits<double>::max();
			for (int r = 0; r < menuPlot.headMD.GetCount(); ++r) {
				double dist = abs(double(menuPlot.headMD.Get(r, 0))) + abs(double(menuPlot.headMD.Get(r, 1)));
				if (distanceTo0 > dist) {
					rowid = r;
					distanceTo0 = dist;
				}
			}
			if (!IsNull(rowid))
				menuPlot.headMD.SetCursor(rowid);
		}
		
		mainTab.WhenSet();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
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
			Bem().RemoveHydro(id);
			listLoaded.Remove(row);
			selected = true;
		}
	}	// Only one available => directly selected
	if (!selected && listLoaded.GetCount() == 1) {	
		int id = ArrayModel_IdHydro(listLoaded, 0);
		Bem().RemoveHydro(id);
		listLoaded.Remove(0);
		selected = true;
	}		
	if (!selected) {
		BEM::PrintError(t_("No model selected"));
		return;
	}
	UpdateButtons();
	
	AfterBEM();
	
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
	menuOpen.dropExport.Enable(numsel == 1);
	menuOpen.butExport.Enable(numsel == 1);
	menuProcess.butSymX.		Enable(numsel == 1 || numrow == 1);
	menuProcess.butSymY.		Enable(numsel == 1 || numrow == 1);
	menuProcess.butKirf.		Enable(numsel == 1 || numrow == 1);
	menuProcess.butA0.			Enable(numsel == 1 || numrow == 1);
	menuProcess.butAinf.		Enable(numsel == 1 || numrow == 1);
	menuProcess.butRAO.			Enable(numsel == 1 || numrow == 1);
	menuAdvanced.butAinfw.		Enable(numsel == 1 || numrow == 1);
	menuAdvanced.butOgilvie.	Enable(numsel == 1 || numrow == 1);
	menuAdvanced.butConvergence.Enable(numsel >= 3);
	menuAdvanced.butAverage.Enable(numsel >= 2);
	menuAdvanced.opDecayingTail.Enable(numsel == 1 || numrow == 1);
	menuAdvanced.opThinremoval. Enable(numsel == 1 || numrow == 1);
	menuAdvanced.opZremoval.	Enable(numsel == 1 || numrow == 1);
	menuAdvanced.opHaskind.		Enable(numsel == 1 || numrow == 1);

	{	
		int row = menuPlot.head1st.GetCursor();
		menuPlot.head1st.Clear();
		for (int ih = 0; ih < Bem().headAll.size(); ++ih)
			menuPlot.head1st.Add(Bem().headAll[ih/*Bem().orderHeadAll[i]*/]);
		if (row >= 0)
			menuPlot.head1st.SetCursor(row);
		
		menuPlot.head1st.WhenCursor = [&] {
			int row = menuPlot.head1st.GetCursor();
			if (row < 0)
				return;
			
			UVector<int> ids = ArrayModel_IdsHydro(listLoaded);
			
			mainForceSC.Load(Bem(), ids, row);
			mainForceFK.Load(Bem(), ids, row);
			mainForceEX.Load(Bem(), ids, row);
			mainRAO.Load(Bem(), ids, row);
		};
	}
	{
		int row = menuPlot.headMD.GetCursor();
		menuPlot.headMD.Clear();
		for (int ih = 0; ih < Bem().headAllMD.size(); ++ih)
			menuPlot.headMD.Add(Bem().headAllMD[ih].real(), Bem().headAllMD[ih].imag());
		if (row >= 0)
			menuPlot.headMD.SetCursor(row);
		
		menuPlot.headMD.WhenCursor = [&] {
			int row = menuPlot.headMD.GetCursor();
			if (row < 0)
				return;
			
			UVector<int> ids = ArrayModel_IdsHydro(listLoaded);
		
			mainMD.Load(Bem(), ids, row);
		};		
	}
	
	bool show_w = menuPlot.opwT == 0;
	bool show_ma_ph = menuPlot.opMP == 0;
	
	menuProcess.dropBody.Clear();
	menuProcess2.dropFreq.GetList().GetColumn(1).Name(show_w ? t_("ω [rad/s]") : t_("T [s]"));
	menuProcess2.dropFreqQTF.GetList().GetColumn(1).Name(show_w ? t_("ω [rad/s]") : t_("T [s]"));
	menuProcess2.dropFreq.Clear();
	menuProcess2.dropFreqQTF.Clear();
	menuProcess2.dropHead.Clear();
	menuProcess2.dropHeadMD.Clear();
	menuProcess2.dropHeadQTF.Clear();
	menuProcess2.dropDOF.Clear();
	menuProcess2.dropForce.Clear();
	menuProcess2.dropForceMD.Clear();
	menuProcess2.dropForceQTF.Clear();
	menuProcess2.opA <<= false;
	menuProcess2.opAd <<= false;
	menuProcess2.opB <<= false;
	menuProcess2.opBd <<= false;
	menuProcess2.opF <<= false;
	menuProcess2.opQTF <<= false;
	
	int id = GetIdOneSelected(false);
	if (id >= 0) { 
		Hydro &data = Bem().hydros[id].hd();
		for (int i = 0; i < data.Nb; ++i)
			menuProcess.dropBody.Add(i+1);
		menuProcess.dropBody.SetIndex(0);
		for (int i = 0; i < data.w.size(); ++i)
			menuProcess2.dropFreq.Add(false, show_w ? data.w[i] : data.T[i]);
		for (int i = 0; i < data.qw.size(); ++i)
			menuProcess2.dropFreqQTF.Add(false, show_w ? data.qw[i] : 2*M_PI/data.qw[i]);
		for (int i = 0; i < data.head.size(); ++i)
			menuProcess2.dropHead.Add(false, data.head[i]);
		for (int i = 0; i < data.mdhead.size(); ++i)
			menuProcess2.dropHeadMD.Add(false, Format("%.1f-%.1f", data.mdhead[i].real(), data.mdhead[i].imag()));
		for (int i = 0; i < data.qh.size(); ++i)
			menuProcess2.dropHeadQTF.Add(false, Format("%.1f-%.1f", data.qh[i].real(), data.qh[i].imag()));
		for (int i = 0; i < 6; ++i)
			menuProcess2.dropDOF.Add(false, BEM::StrDOF(i));
		
		if (data.IsLoadedFex())
			menuProcess2.dropForce.Add(false, t_("All"));
		if (data.IsLoadedFsc())
			menuProcess2.dropForce.Add(false, t_("Diffraction"));
		if (data.IsLoadedFfk())
			menuProcess2.dropForce.Add(false, t_("Froude-Krylov"));

		if (data.IsLoadedMD())
			menuProcess2.dropForceMD.Add(false, t_("All"));

		if (data.IsLoadedQTF(true) || data.IsLoadedQTF(false))
			menuProcess2.dropForceQTF.Add(false, t_("All"));
		if (data.IsLoadedQTF(true))
			menuProcess2.dropForceQTF.Add(false, t_("Summation"));
		if (data.IsLoadedQTF(false))
			menuProcess2.dropForceQTF.Add(false, t_("Difference"));
	}
}

void MainBEM::OnJoin() {
	UVector<int> idsjoin, rowsJoin;
	for (int row = listLoaded.GetCount()-1; row >= 0; --row) {
		if (listLoaded.IsSelected(row)) {
			rowsJoin << row;
			int id = ArrayModel_IdHydro(listLoaded, row);
			idsjoin << id;
		}
	}
	if (idsjoin.IsEmpty()) {
		BEM::PrintError(t_("No model selected"));
		return;
	}
	if (idsjoin.size() == 1) {
		BEM::PrintError(t_("Please select more than one model"));
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
		
		ArrayModel_Add(listLoaded, data.hd().GetCodeStr(), data.hd().name, data.hd().file, data.hd().GetId());
		ArrayModel_RowsHydroDel(listLoaded, rowsJoin);
	
		for (int id = 0; id < Bem().hydros.size(); ++id) {
			const Hydro &data = Bem().hydros[id].hd();
			mainSummary.Report(data, id);
		}
			
		if (Bem().hydros.size() > 0) 
			listLoaded.SetCursor(0);

		AfterBEM();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
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

		AfterBEM();
		
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
}

void MainBEM::OnSymmetrizeForces(bool xAxis) {
	try {
		int id = GetIdOneSelected();
		if (id < 0) 
			return;

		WaitCursor wait;

		Progress progress(t_("Symmetrizing forces and RAOs in selected BEM file..."), 100); 
		
		Bem().SymmetrizeForces(id, xAxis);
		
		AfterBEM();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
}

void MainBEM::OnKirfAinf(Hydro::DataToPlot param) {
	try {
		int id = GetIdOneSelected();
		if (id < 0) 
			return;
			
		Progress progress(Format(t_("Calculating %s in selected BEM file..."), Hydro::StrDataToPlot(param)), 100); 
		
		double maxT = Null;
		
		if (param == Hydro::PLOT_KIRF || 
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
		else if (param == Hydro::PLOT_KIRF)
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
		
		AfterBEM();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
}

void MainBEM::OnRAO() {
	try {
		int id = GetIdOneSelected();
		if (id < 0) 
			return;
			
		Progress progress(t_("Calculating RAO in selected BEM file..."), 100); 
		
		WaitCursor wait;

		Bem().RAO(id);
		
		AfterBEM();	
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
}	

void MainBEM::OnSymmetrize() {
	try {
		int id = GetIdOneSelected();
		if (id < 0) 
			return;
			
		WaitCursor wait;
		
		Bem().Symmetrize(id);
				
		AfterBEM();

	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MainBEM::OnOgilvie() {
	String str;
	try {
		int id = GetIdOneSelected();
		if (id < 0) 
			return;
			
		Progress progress(t_("Calculating Ogilvie compliance in selected BEM file..."), 100); 
		
		WaitCursor wait;
		
		//double maxT = Null;
		UVector<int> vidof, vjdof;
		Bem().OgilvieCompliance(id, ~menuAdvanced.opZremoval, ~menuAdvanced.opThinremoval, ~menuAdvanced.opDecayingTail, ~menuAdvanced.opHaskind, vidof, vjdof);
				
		AfterBEM();
		
		if (vidof.size() > 0) {
			str = "The degrees of freedom healed are:\n";
			for (int i = 0; i < vidof.size(); ++i) {
				str << Format("[%d, %d] ", vidof[i]+1, vjdof[i]+1);
				if (i > 0 && vidof[i-1] != vidof[i])
					str << "\n";
			}
		} else
			str = "No degrees of freedom has been healed";
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
	PromptOK(DeQtfLf(str));
}

void MainBEM::OnAverage() {
	try {
		UVector<int> ids;
		for (int row = 0; row < listLoaded.GetCount(); ++row) {
			if (listLoaded.IsSelected(row)) {
				int id = ArrayModel_IdHydro(listLoaded, row);
				ids << id;
			}
		}		
		if (ids.size() < 2) {
			BEM::PrintError(t_("Not enough models selected"));
			return;
		}
			
		Progress progress(t_("Calculating average..."), 100); 
		
		WaitCursor wait;
		
		HydroClass &data = Bem().Average(ids);
		
		mainSummary.Clear();
		
		ArrayModel_Add(listLoaded, data.hd().GetCodeStr(), data.hd().name, data.hd().file, data.hd().GetId());
				
		AfterBEM();

	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MainBEM::OnConvergence() {
	String str;
	try {
		UVector<int> ids;
		for (int row = 0; row < listLoaded.GetCount(); ++row) {
			if (listLoaded.IsSelected(row)) {
				int id = ArrayModel_IdHydro(listLoaded, row);
				ids << id;
			}
		}		
		if (ids.size() < 3) {
			BEM::PrintError(t_("Not enough models selected"));
			return;
		}
			
		Progress progress(t_("Calculating asympthotic convergence..."), 100); 
		
		WaitCursor wait;
		
		UVector<int> vidof, vjdof;
		
		//Bem().OgilvieCompliance(id, ~menuAdvanced.opZremoval, ~menuAdvanced.opThinremoval, ~menuAdvanced.opDecayingTail, ~menuAdvanced.opHaskind, vidof, vjdof);
				
		AfterBEM();
		
		if (vidof.size() > 0) {
			str = "The degrees of freedom healed are:\n";
			for (int i = 0; i < vidof.size(); ++i) {
				str << Format("[%d, %d] ", vidof[i]+1, vjdof[i]+1);
				if (i > 0 && vidof[i-1] != vidof[i])
					str << "\n";
			}
		} else
			str = "No degrees of freedom has been healed";
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
	PromptOK(DeQtfLf(str));
	


	
}

void MainBEM::OnResetForces() {
	try {
		int id = GetIdOneSelected();
		if (id < 0) 
			return;
		
		
		Hydro::FORCE force = Hydro::NONE;
		for (int i = 0; i < menuProcess2.dropForce.GetCount(); ++i) {
			if (menuProcess2.dropForce.GetList().Get(i, 0) == true) {
				String val = menuProcess2.dropForce.GetList().Get(i, 1);
				if (val == t_("All")) {
					force = Hydro::ALL;
					break;
				} else if (val == t_("Diffraction")) {
					if (force == Hydro::NONE)
						force = Hydro::SCATTERING;
					else {
						force = Hydro::ALL;
						break;
					}
				} else if (val == t_("Froude-Krylov")) {
					if (force == Hydro::NONE)
						force = Hydro::FK;
					else {
						force = Hydro::ALL;
						break;
					}
				}
			}
		}

		bool forceMD = false;
		for (int i = 0; i < menuProcess2.dropForceMD.GetCount(); ++i) {
			if (menuProcess2.dropForceMD.GetList().Get(i, 0) == true) {
				String val = menuProcess2.dropForceMD.GetList().Get(i, 1);
				if (val == t_("All")) {
					forceMD = true;
					break;
				}
			}
		}
		
		Hydro::FORCE forceQtf = Hydro::NONE;
		for (int i = 0; i < menuProcess2.dropForceQTF.GetCount(); ++i) {
			if (menuProcess2.dropForceQTF.GetList().Get(i, 0) == true) {
				String val = menuProcess2.dropForceQTF.GetList().Get(i, 1);
				if (val == t_("All")) {
					forceQtf = Hydro::ALL;
					break;
				} else if (val == t_("Summation")) {
					if (forceQtf == Hydro::NONE)
						forceQtf = Hydro::QTFSUM;
					else {
						forceQtf = Hydro::ALL;
						break;
					}
				} else if (val == t_("Difference")) {
					if (forceQtf == Hydro::NONE)
						forceQtf = Hydro::QTFDIF;
					else {
						forceQtf = Hydro::ALL;
						break;
					}
				}
			}
		}
					
		WaitCursor wait;
		
		Bem().ResetForces(id, force, forceMD, forceQtf);
				
		mainSummary.Clear();
		for (int i = 0; i < Bem().hydros.size(); ++i)
			mainSummary.Report(Bem().hydros[i].hd(), i);
		
		UVector<int> ids = ArrayModel_IdsHydro(listLoaded);
		
		mainTab.GetItem(mainTab.Find(mainForceSC)).Enable(mainForceSC.Load(Bem(), ids, menuPlot.head1st.GetCursor()));
		mainTab.GetItem(mainTab.Find(mainForceFK)).Enable(mainForceFK.Load(Bem(), ids, menuPlot.head1st.GetCursor()));
		mainTab.GetItem(mainTab.Find(mainForceEX)).Enable(mainForceEX.Load(Bem(), ids, menuPlot.head1st.GetCursor()));
		mainTab.GetItem(mainTab.Find(mainRAO)).Enable(mainRAO.Load(Bem(), ids));
		mainTab.GetItem(mainTab.Find(mainQTF)).Enable(mainQTF.Load());
	
		LoadSelTab(Bem());
		
		if (menuProcess2.dropForce.GetCount() > 0)
			menuProcess2.dropForce.SetIndex(0);
		if (menuProcess2.dropForceQTF.GetCount() > 0)
			menuProcess2.dropForceQTF.SetIndex(0);
		menuProcess2.dropForce.OnFocus();
		menuProcess2.dropForceQTF.OnFocus();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MainBEM::OnMultiplyDOF(bool isReset) {
	try {
		int id = GetIdOneSelected();
		if (id < 0) 
			return;
		
		if (isReset)
			menuProcess2.factor <<= 0;
			
		double factor = ~menuProcess2.factor;
		if (IsNull(factor))
			factor = 0;
				
		UVector<int> idDOF;
		for (int i = 0; i < 6; ++i)
			if (menuProcess2.dropDOF.GetList().Get(i, 0) == true)
				idDOF << i;
			
		WaitCursor wait;
		
		bool a = menuProcess2.opA || menuProcess2.opAd;
		bool b = menuProcess2.opB || menuProcess2.opBd;
		bool diag = menuProcess2.opAd || menuProcess2.opBd;
		
		Bem().MultiplyDOF(id, factor, idDOF, a, b, diag, menuProcess2.opF, menuProcess2.opMD, menuProcess2.opQTF);
				
		AfterBEM();	
		LoadSelTab(Bem());
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MainBEM::OnSwapDOF() {
	try {
		int id = GetIdOneSelected();
		if (id < 0) 
			return;
		
		WaitCursor wait;
		
		Bem().SwapDOF(id, menuProcess.dropBody.GetIndex(), 
				menuProcess.dropDOF1.GetIndex()-1, menuProcess.dropDOF2.GetIndex()-1);
				
		mainSummary.Clear();
		for (int i = 0; i < Bem().hydros.size(); ++i)
			mainSummary.Report(Bem().hydros[i].hd(), i);
		
		UVector<int> ids = ArrayModel_IdsHydro(listLoaded);
		
		mainTab.GetItem(mainTab.Find(mainA)).Enable(mainA.Load(Bem(), ids));
		mainTab.GetItem(mainTab.Find(mainB)).Enable(mainB.Load(Bem(), ids));
		mainTab.GetItem(mainTab.Find(mainMD)).Enable(mainMD.Load(Bem(), ids, menuPlot.headMD.GetCursor()));
		mainTab.GetItem(mainTab.Find(mainK)).Enable(mainK.Load(Bem(), ids));
		mainTab.GetItem(mainTab.Find(mainAinfw)).Enable(mainAinfw.Load(Bem(), ids));
		mainTab.GetItem(mainTab.Find(mainForceSC)).Enable(mainForceSC.Load(Bem(), ids, menuPlot.head1st.GetCursor()));
		mainTab.GetItem(mainTab.Find(mainForceFK)).Enable(mainForceFK.Load(Bem(), ids, menuPlot.head1st.GetCursor()));
		mainTab.GetItem(mainTab.Find(mainForceEX)).Enable(mainForceEX.Load(Bem(), ids, menuPlot.head1st.GetCursor()));
		mainTab.GetItem(mainTab.Find(mainRAO)).Enable(mainRAO.Load(Bem(), ids, menuPlot.head1st.GetCursor()));	
		mainTab.GetItem(mainTab.Find(mainQTF)).Enable(mainQTF.Load());
	
		LoadSelTab(Bem());
		
		if (menuProcess.dropDOF1.GetCount() > 0)
			menuProcess.dropDOF1.SetIndex(0);
		if (menuProcess.dropDOF2.GetCount() > 0)
			menuProcess.dropDOF2.SetIndex(0);
		menuProcess.dropDOF1.WhenAction();
		menuProcess.dropDOF2.WhenAction();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MainBEM::OnABForces() {
	try {
		int id = GetIdOneSelected();
		if (id < 0) 
			return;
		
		WaitCursor wait;
		
		Bem().FillFrequencyGapsABForces(id, ~menuProcess.opFill == 0, ~menuProcess.maxFreq);
				
		AfterBEM();
		LoadSelTab(Bem());
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MainBEM::OnQTF() {
	try {
		int id = GetIdOneSelected();
		if (id < 0) 
			return;
		
		WaitCursor wait;
		
		Bem().FillFrequencyGapsQTF(id, ~menuProcess.opFill == 0, ~menuProcess.maxFreq);
				
		AfterBEM();
		LoadSelTab(Bem());
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MainBEM::OnABForcesZero() {
	try {
		int id = GetIdOneSelected();
		if (id < 0) 
			return;
		
		WaitCursor wait;
		
		Bem().FillFrequencyGapsABForcesZero(id);
		
		AfterBEM();
		LoadSelTab(Bem());
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MainBEM::OnQTFZero() {
	try {
		int id = GetIdOneSelected();
		if (id < 0) 
			return;
		
		WaitCursor wait;
		
		Bem().FillFrequencyGapsQTFZero(id);
		
		AfterBEM();
		LoadSelTab(Bem());
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MainBEM::OnQTF_MD() {
	try {
		int id = GetIdOneSelected();
		if (id < 0) 
			return;
		
		WaitCursor wait;
		
		Bem().CopyQTF_MD(id);
		
		Bem().UpdateHeadAllMD();
		
		AfterBEM();
		LoadSelTab(Bem());
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MainBEM::OnDeleteHeadingsFrequencies() {
	try {
		int id = GetIdOneSelected();
		if (id < 0) 
			return;
		
		UVector<int> idFreq, idFreqQTF, idHead, idHeadMD, idHeadQTF;
		for (int i = 0; i < menuProcess2.dropFreq.GetList().GetRowCount(); ++i)
			if (menuProcess2.dropFreq.GetList().Get(i, 0) == true)
				idFreq << i;
		for (int i = 0; i < menuProcess2.dropFreqQTF.GetList().GetRowCount(); ++i)
			if (menuProcess2.dropFreqQTF.GetList().Get(i, 0) == true)
				idFreqQTF << i;
		for (int i = 0; i < menuProcess2.dropHead.GetList().GetRowCount(); ++i)
			if (menuProcess2.dropHead.GetList().Get(i, 0) == true)
				idHead << i;
		for (int i = 0; i < menuProcess2.dropHeadMD.GetList().GetRowCount(); ++i)
			if (menuProcess2.dropHeadMD.GetList().Get(i, 0) == true)
				idHeadMD << i;
		for (int i = 0; i < menuProcess2.dropHeadQTF.GetList().GetRowCount(); ++i)
			if (menuProcess2.dropHeadQTF.GetList().Get(i, 0) == true)
				idHeadQTF << i;
		
		WaitCursor wait;
		
		Bem().DeleteHeadingsFrequencies(id, idFreq, idFreqQTF, idHead, idHeadMD, idHeadQTF);
		
		AfterBEM();	
		LoadSelTab(Bem());
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MainBEM::OnUpdateCrot() {
	try {
		int id = GetIdOneSelected();
		if (id < 0) 
			return;
		
		WaitCursor wait;
		
		Bem().TranslationTo(id, ~menuAdvanced.x_0, ~menuAdvanced.y_0, ~menuAdvanced.z_0);
				
		AfterBEM();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MainBEM::AfterBEM() {
	mainSummary.Clear();
	for (int id = 0; id < Bem().hydros.size(); ++id) {
		Hydro &data = Bem().hydros[id].hd();
		mainSummary.Report(data, id);
	}
	
	UVector<int> ids = ArrayModel_IdsHydro(listLoaded);

	Progress progress(t_("Processing loaded data..."), 14);
	int pos = 0;
	mainTab.GetItem(mainTab.Find(mainMatrixA)).Enable(mainMatrixA.Load(Bem().hydros, ids, ~menuPlot.showNdim));		progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainMatrixM)).Enable(mainMatrixM.Load(Bem().hydros, ids, false));					progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainMatrixK)).Enable(mainMatrixK.Load(Bem().hydros, ids, ~menuPlot.showNdim));		progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainMatrixDlin)).Enable(mainMatrixDlin.Load(Bem().hydros, ids, false));			progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainA)).Enable(mainA.Load(Bem(), ids));											progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainAinfw)).Enable(mainAinfw.Load(Bem(), ids));									progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainB)).Enable(mainB.Load(Bem(), ids));											progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainK)).Enable(mainK.Load(Bem(), ids));											progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainMD)).Enable(mainMD.Load(Bem(), ids, menuPlot.headMD.GetCursor()));				progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainForceSC)).Enable(mainForceSC.Load(Bem(), ids, menuPlot.head1st.GetCursor()));	progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainForceFK)).Enable(mainForceFK.Load(Bem(), ids, menuPlot.head1st.GetCursor()));	progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainForceEX)).Enable(mainForceEX.Load(Bem(), ids, menuPlot.head1st.GetCursor()));	progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainRAO)).Enable(mainRAO.Load(Bem(), ids, menuPlot.head1st.GetCursor()));			progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainQTF)).Enable(mainQTF.Load());													progress.SetPos(pos++);

	bool isLoadedSS = false;
	for (int id = 0; id < Bem().hydros.size(); ++id) {
		Hydro &data = Bem().hydros[id].hd();
		if (data.IsLoadedStateSpace()) {
			isLoadedSS = true;
			break;
		}
	}	
	mainTab.GetItem(mainTab.Find(mainStateSpace)).Enable(isLoadedSS);
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
	int numH = int(hydro.qh.size());
	if (numH <= 1)
		return Null;
	
	WithSaveQTF<TopWindow> dialog;
	CtrlLayout(dialog);
	
	dialog.Title(t_("Please choose the QTF headings to save"));
	
	dialog.swHeadings <<= 0;
	dialog.dropHeadings.Enable(false);
	
	int id0 = Null;
	for (int i = 0; i < numH; ++i) {
		dialog.dropHeadings.Add(Format("%f-%f", hydro.qh[i].real(), hydro.qh[i].imag()));
		if (abs(hydro.qh[i].real()) < 0.001)
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
		return Null;
	else if (dialog.swHeadings == 1)
		return -1;
	else
		return dialog.dropHeadings.GetIndex();
}

void MainBEM::OnConvert() {
	GuiLock __;
	
	try {
		int id = GetIdOneSelected();
		if (id < 0) 
			return;
		
		Status(t_("Saving BEM data"));
		String fileType = ~menuOpen.dropExport;
		Hydro::BEM_FMT type = Hydro::GetCodeBemStr(fileType);
		String ext = Replace(Hydro::bemExt[type], "*", "");
		
		FileSel fs;
		
		for (int i = 0; i < Hydro::GetBemStrCount(); ++i)
			if (Hydro::bemCanSave[i] && (i == type || i == Hydro::UNKNOWN)) 
				fs.Type(Hydro::GetBemStr(static_cast<Hydro::BEM_FMT>(i)), Hydro::bemExt[i]);
		
		fs.ActiveType(0);
		fs.Set(ForceExt(~menuOpen.file, ext));
		fs.ActiveDir(saveFolder);
		
		if (!fs.ExecuteSaveAs(Format(t_("Save BEM data as %s"), fileType)))
			return;
		
		String fileName = ~fs;
		
		int qtfHeading = Null;
		
		if (type == Hydro::WAMIT_1_3 || type == Hydro::FAST_WAMIT || (type == Hydro::UNKNOWN && GetFileExt(fileName) == ".1"))	
			qtfHeading = AskQtfHeading(Bem().hydros[id].hd());	
		
		Progress progress(t_("Saving BEM files..."), 100); 
		progress.Granularity(1000);
		
		WaitCursor wait;
		
		Bem().hydros[id].hd().SaveAs(fileName, [&](String str, int _pos) {
			if (!IsEmpty(str))
				progress.SetText(str); 
			if (_pos >= 0)
				progress.SetPos(_pos); 
			return !progress.Canceled();}, type, qtfHeading);
			
		saveFolder = GetFileFolder(~fs);
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
}

int MainBEM::GetIdOneSelected(bool complain) {
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
		if (complain)
			BEM::PrintError(t_("No model selected"));
		return -1;
	}
	if (ArrayCtrlSelectedGetCount(listLoaded) > 1) {
		if (complain)
			BEM::PrintError(t_("Please select just one model"));
		return -1;
	}
	return id;
}

void MainBEM::Jsonize(JsonIO &json) {
	if (json.IsLoading()) {
		menuPlot.autoFit = Null;
		menuPlot.fromY0 = Null;
		menuPlot.opwT = Null;
		menuPlot.opMP = Null;
		menuPlot.showPoints = Null;
		menuPlot.showNdim = Null;
		dropExportId = 2;
		menuProcess.opFill = Null;
		menuProcess.maxFreq = Null;
	} else {
		dropExportId = menuOpen.dropExport.GetIndex();
	}
	json
		("menuOpen_file", menuOpen.file)
		("menuOpen_saveFolder", saveFolder)
		("menuOpen_dropExport", dropExportId)
		("menuPlot_autoFit", menuPlot.autoFit)
		("menuPlot_fromY0", menuPlot.fromY0)
		("menuPlot_opwT", menuPlot.opwT)
		("menuPlot_opMP", menuPlot.opMP)
		("menuPlot_showPoints", menuPlot.showPoints)
		("menuPlot_showNdim", menuPlot.showNdim)
		("menuProcess_opZremoval", menuAdvanced.opZremoval)
		("menuProcess_opThinremoval", menuAdvanced.opThinremoval)
		("menuProcess_opDecayingTail", menuAdvanced.opDecayingTail)
		("menuProcess_opFill", menuProcess.opFill) 
		("menuProcess_maxFreq", menuProcess.maxFreq)
		("opHaskind", menuAdvanced.opHaskind)
		("mainStiffness", mainMatrixK)
		("mainMatrixA", mainMatrixA)
		("mainMatrixDlin", mainMatrixDlin)
		("mainMatrixM", mainMatrixM)
	;
}

String MainBEM::BEMFile(String fileFolder) const {
	if (ToLower(fileFolder) == "id.dat")		// This a Nemoh file
		return "";
	if (DirectoryExists(fileFolder)) {
		int bestipos = INT_MAX;
		for (FindFile ff(AppendFileNameX(fileFolder, "*.*")); ff; ff++) {
			if (ff.IsFile()) {
				String name = ToLower(ff.GetName());
				if (GetFileExt(name) == ".bem")
					return ff.GetPath();
				if (name == "nemoh.cal")
					return ff.GetPath();
				if (ff.GetName() == "ControlFile.in") 
					return ff.GetPath();
				int ipos = Bem().bstFilesExt.Find(GetFileExt(name));
 				if (ipos >= 0 && ipos < bestipos) {
					fileFolder = ff.GetPath();
					bestipos = ipos;	// It takes the file with most probable extension
				}
			} else if (ff.IsFolder() && ff.GetName() == "Input") {
				String controlFile = AppendFileNameX(ff.GetPath(), "ControlFile.in");
				if (FileExists(controlFile))
					return controlFile;
			}
		}
	}
	return fileFolder;
}

void MainBEM::LoadDragDrop() {
	GuiLock __;
	
	if (filesToDrop.size() == 1 && DirectoryExists(First(filesToDrop))) {	// If a folder is dropped, all the files are included
		String folder = First(filesToDrop);
		filesToDrop.Clear();
		for (FindFile ff(AppendFileNameX(folder, "*.*")); ff; ff++) {
			if (ff.IsFile())
				filesToDrop << ff.GetPath();
		}
	}
	for (int i = filesToDrop.size()-1; i >= 0; --i) {						// Remove files with unknown extensions
		if (Bem().bstFilesExt.Find(ToLower(GetFileExt(filesToDrop[i]))) < 0)
			filesToDrop.Remove(i);
	}
	
	Sort(filesToDrop);
	
	UVector<int> sets(filesToDrop.size());									// Remove files in a set
	for (int i = 0; i < filesToDrop.size(); ++i)
		sets[i] = Bem().GetBEMExtSet(filesToDrop[i]);
	
	for (int i = filesToDrop.size()-1; i > 0; --i) {
		for (int j = 0; j < i; ++j) {
			if (sets[i] >= 0 && 
				sets[i] == sets[j] && 
				ToLower(GetFileFolder(filesToDrop[i])) == ToLower(GetFileFolder(filesToDrop[j])) && 
				ToLower(GetFileTitle (filesToDrop[i])) == ToLower(GetFileTitle (filesToDrop[j]))) {		// Removes files that are loaded in a set, like .lis .qtf, or .1 .3 .hst
				filesToDrop.Remove(i);
				break;
			}
		}
	}
	if (filesToDrop.IsEmpty()) {
		BEM::PrintError(t_("Impossible to load files"));
		timerDrop.Kill();
		return;	
	}
	bool followWithErrors = false;
	for (int i = 0; i < filesToDrop.size(); ++i) {
		String file = BEMFile(filesToDrop[i]);
		if (!file.IsEmpty()) {
			menuOpen.file <<= file;
			Status(Format(t_("Loading '%s'"), file));
			if (!OnLoad() && !followWithErrors && filesToDrop.size() - i > 1) {
				if (!PromptYesNo(Format(t_("Do you wish to load the pending %d files?"), filesToDrop.size() - i - 1))) {
					timerDrop.Kill();
					return;
				}
				followWithErrors = true;
			}
			ProcessEvents();
		}
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
	
	Upp::Color lightRed = Upp::Color(255, 165, 158);								
	
	array.Set(row, 0, t_("File"));				array.Set(row++, col, data.file);
	array.Set(row, 0, t_("Name"));				array.Set(row++, col, data.name);
	array.Set(row, 0, t_("Description"));		array.Set(row++, col, data.description);
	array.Set(row, 0, t_("Software"));			array.Set(row++, col, data.GetCodeStr());
	array.Set(row, 0, t_("g [m/s2]"));			array.Set(row++, col, data.S_g());
	array.Set(row, 0, t_("rho [kg/m3]"));		array.Set(row++, col, data.S_rho());
	array.Set(row, 0, t_("h (water depth) [m]"));array.Set(row++,col, data.S_h());
	array.Set(row, 0, t_("length scale [m]"));	array.Set(row++, col, data.S_len());
	
	array.Set(row, 0, t_("#frequencies"));		array.Set(row++, col, Nvl2(data.Nf, 0)); 
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
	
	array.Set(row, 0, t_("#1st order headings"));			array.Set(row++, col, Nvl2(data.Nh, 0));
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
	array.Set(row, 0, t_("A available"));		array.Set(row++, col, data.IsLoadedA() 	  ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("B available"));		array.Set(row++, col, data.IsLoadedB() 	  ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("K available"));		array.Set(row++, col, data.IsLoadedC() 	  ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("Inertia available"));	array.Set(row++, col, data.IsLoadedM() 	  ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("Fex available"));		array.Set(row++, col, data.IsLoadedFex()  ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("Fsc available"));		array.Set(row++, col, data.IsLoadedFsc()  ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("Ffk available"));		array.Set(row++, col, data.IsLoadedFfk()  ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("RAO available"));		array.Set(row++, col, data.IsLoadedRAO()  ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("Linear damping available"));	array.Set(row++, col, data.IsLoadedDlin()   ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("Mean Drift available"));		array.Set(row++, col, data.IsLoadedMD() 	? t_("Yes") : t_("No"));
	
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
		if (data.Vo.size() > ib && IsNum(data.Vo[ib])) {
			
			array.Set(row++, col, AttrText(FDS(data.Vo[ib], 10, false)).Paper(data.Vo[ib] < 0 ? lightRed : Null));
		} else 
			array.Set(row++, col, "-");
		
		array.Set(row, 0, sib + " " + t_("Cg [m]"));
		if (data.cg.size() > 3*ib && IsNum(data.cg(0, ib))) 
			array.Set(row++, col, Format(t_("%s, %s, %s"),
									FDS(data.cg(0, ib), 10, false),
									FDS(data.cg(1, ib), 10, false),
									FDS(data.cg(2, ib), 10, false)));
		else
			array.Set(row++, col, "-");

		array.Set(row, 0, sib + " " + t_("Cb [m]"));
		if (data.cb.size() > 3*ib && IsNum(data.cb(0, ib))) 
			array.Set(row++, col, Format(t_("%s, %s, %s"),
									FDS(data.cb(0, ib), 10, false),
									FDS(data.cb(1, ib), 10, false),
									FDS(data.cb(2, ib), 10, false)));
		else
			array.Set(row++, col, "-");
		
		array.Set(row, 0, sib + " " + t_("C0 [m]"));
		if (data.c0.size() > 3*ib && IsNum(data.c0(0, ib))) 
			array.Set(row++, col, Format(t_("%s, %s, %s"),
									FDS(data.c0(0, ib), 10, false),
									FDS(data.c0(1, ib), 10, false),
									FDS(data.c0(2, ib), 10, false)));
		else
			array.Set(row++, col, "-");
		
		array.Set(row, 0, sib + " " + t_("Waterplane area [m2]"));
		if (data.C.size() > ib && data.C[ib].size() > 0) {
			double wPlaneArea = data.C_dim(ib, 2, 2);
			array.Set(row++, col, FDS(wPlaneArea, 10, false));		
			for (int i = 0; i < 6; ++i) {
				for (int j = 0; j < 6; ++j) {
					if (!Hydro::C_units(i, j).IsEmpty()) {
						array.Set(row, 0, sib + " " + Format(t_("K(%d,%d) [%s]"), i+1, j+1, Hydro::C_units(i, j)));	array.Set(row++, col, FDS(data.C_dim(ib, i, j), 10, false));		
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
		array.Set(row,   0, sib + " " + t_("Theave [s]"));
		array.Set(row+1, 0, sib + " " + t_("Troll [s]"));
		array.Set(row+2, 0, sib + " " + t_("Tpitch [s]"));
		if (IsNum(data.rho) && IsNum(data.g) && 
			data.M.size() > ib && data.M[ib].size() > 0 && 
			data.C.size() > ib && data.C[ib].size() > 0) {
			array.Set(row++, col, FDS(data.Theave(ib), 5, false));
			array.Set(row++, col, FDS(data.Troll(ib), 5, false));
			array.Set(row++, col, FDS(data.Tpitch(ib), 5, false));
		} else {
			array.Set(row++, col, "-");	
			array.Set(row++, col, "-");	
			array.Set(row++, col, "-");	
		}
		array.Set(row,   0, sib + " " + t_("GMroll [m]"));
		array.Set(row+1, 0, sib + " " + t_("GMpitch [m]"));
		if (IsNum(data.rho) && IsNum(data.g) && data.IsLoadedC()) {
			array.Set(row++, col, FDS(data.GMroll(ib), 5, false));
			array.Set(row++, col, FDS(data.GMpitch(ib), 5, false));
		} else {
			array.Set(row++, col, "-");	
			array.Set(row++, col, "-");	
		}
	}	
}

void MainOutput::Init() {
	CtrlLayout(*this);
	cout.SetReadOnly();
	Print(t_("BEMRosetta\nHydrodynamic coefficients viewer and converter for Boundary Element Method solver formats\n"));
}

void MainOutput::Print(String str) {
	cout.Append(str);
	cout.ScrollEnd();
}

void MainBEMW::Init(MainBEM &_bem, const Image &icon, const Image &largeIcon, Function <void()> _WhenClose) {
	WhenClose = _WhenClose;
	LoadFromJson(bem, StoreAsJson(_bem));
	bem.Init();
	Add(bem.SizePos());
	Title(t_("BEMRosetta BEM Coefficients Processing")).Sizeable().Zoomable().Icon(icon, largeIcon);
}