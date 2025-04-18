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
	MainBEMBody::Init();
	
	CtrlLayout(menuOpen);
	menuOpen.file.WhenChange = THISBACK(OnLoad);
	menuOpen.file.BrowseRightWidth(40).UseOpenFolder().BrowseOpenFolderWidth(10);
	menuOpen.butLoad << [&] {menuOpen.file.DoGo();};
	
	ArrayModel_Init(listLoaded).MultiSelect();
	listLoaded.WhenSel = [&] {
		OnMenuAdvancedArraySel(true);
		menuFOAMM.OnCursor();
		//mainTab.GetItem(mainTab.Find(mainQTF)).Enable(mainQTF.Load());
		if (mainTab.Find(mainQTF) == mainTab.Get())
			mainQTF.Load();
		else if (mainTab.Find(mainBody) == mainTab.Get())
			mainBody.Load();
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
	menuOpen.butRename.Disable();
	menuOpen.butRename <<= THISBACK(OnDescription);
	menuOpen.butSolve.Disable();
	menuOpen.butSolve <<= THISBACK(OnSolve);
	menuOpen.butExport <<= THISBACK(OnConvert);
	menuOpen.butExport.Tip(t_("Exports data file"));
	for (int i = 0; i < Hydro::NUMBEM; ++i)
		if (Hydro::bemCanSave[i])
			menuOpen.dropExport.Add(Hydro::GetBemStr(static_cast<Hydro::BEM_FMT>(i)));
	menuOpen.dropExport.SetIndex(dropExportId);

	CtrlLayout(menuProcess);
	menuProcess.butSym.Disable();	
	menuProcess.butSym <<= THISBACK(OnSymmetrizeForces);
	
	menuProcess.dropSym.Add(BasicBEM::SYM_NO, "No symmetry").Add(BasicBEM::SYM_XZ, "XZ").Add(BasicBEM::SYM_YZ, "YZ").Add(BasicBEM::SYM_XZ_YZ, "XZ+YZ");//.Add(SYM_AXISYMMETRIC, "Axisymmetric");
	menuProcess.dropSym.SetIndex(0);
	menuProcess.dropSym.WhenAction = [&]() {UpdateButtons();};
	
	menuProcess.butA0.Disable();	
	menuProcess.butA0 <<= THISBACK1(OnKirfAinf, Hydro::PLOT_A0);
	menuProcess.butAinf.Disable();	
	menuProcess.butAinf <<= THISBACK1(OnKirfAinf, Hydro::PLOT_AINF);
	menuProcess.butKirf.Disable();	
	menuProcess.butKirf <<= THISBACK1(OnKirfAinf, Hydro::PLOT_KIRF);
	menuProcess.butRAO.Disable();	
	menuProcess.critDamp.Disable();	
	menuProcess.critDamp = 0.05;
	menuProcess.butRAO <<= THISBACK(OnRAO);
	menuProcess.butSymmetrize <<= THISBACK(OnSymmetrize);
	
	menuProcess.butABForces << THISBACK(OnABForces);
	menuProcess.butQTF << THISBACK(OnQTF);

	menuProcess.butABForcesZero << THISBACK(OnABForcesZero);
	menuProcess.butQTFZero << THISBACK(OnQTFZero);
	
	menuProcess.opFill.Tip(t_("Fills with zeroes or with interpolated values"));

	auto DropDOF = [&](DropList &b1, DropList &dof1, DropList &b2, DropList &dof2)->bool {
		if (dof1.GetIndex() != dof2.GetIndex() && (dof1.GetIndex() == 0 || dof2.GetIndex() == 0)) 
			return false;
		if (b1.GetIndex() == b2.GetIndex()) {
			if (dof1.GetIndex() != 0) {
				if (dof1.GetIndex() == dof2.GetIndex())
					return false;
			} else 
				return false;
		}
		return true;
	};
	
	menuProcess.dropDOF1.WhenAction = [&] {
		menuProcess.butSwapDOF.Show(DropDOF(menuProcess.dropBody1, menuProcess.dropDOF1, menuProcess.dropBody2, menuProcess.dropDOF2));
	};
	menuProcess.dropDOF1.Add(t_("All"));
	for (int i = 0; i < 6; ++i)
		menuProcess.dropDOF1.Add(BEM::StrDOF(i));
	menuProcess.dropDOF1.SetIndex(0);
	
	menuProcess.dropDOF2.WhenAction = menuProcess.dropDOF1.WhenAction;
	menuProcess.dropDOF2.Add(t_("All"));
	for (int i = 0; i < 6; ++i)
		menuProcess.dropDOF2.Add(BEM::StrDOF(i));
	menuProcess.dropDOF2.SetIndex(0);
	
	menuProcess.dropBody1.WhenAction = menuProcess.dropBody2.WhenAction = menuProcess.dropDOF1.WhenAction;
	
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
					 menuProcess2.opStiffness || menuProcess2.opQTF);
		menuProcess2.butResetDOF.Show(show);
		menuProcess2.butResetDOF0.Show(show);
	};
	menuProcess2.opA.WhenAction = menuProcess2.opAd.WhenAction =
		menuProcess2.opB.WhenAction = menuProcess2.opBd.WhenAction =
		menuProcess2.opF.WhenAction = menuProcess2.opMD.WhenAction =
		menuProcess2.opQTF.WhenAction = menuProcess2.opStiffness.WhenAction = menuProcess2.dropDOF.OnFocus;
	
	
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
	menuAdvanced.butAinfw.Disable();
	menuAdvanced.butBH.WhenAction = menuAdvanced.numBH.WhenEnter = [&]() {OnBH(int(menuAdvanced.numBH));};
	menuAdvanced.butBH.Disable();
	menuAdvanced.butOgilvie <<= THISBACK(OnOgilvie);
	menuAdvanced.butOgilvie.Disable();
	menuAdvanced.opDecayingTail.Disable();
	menuAdvanced.opThinremoval.Disable();
	menuAdvanced.opZremoval.Disable();
	menuAdvanced.butConvergence <<= THISBACK(OnConvergence);
	menuAdvanced.butConvergence.Disable();
	menuAdvanced.butAverage <<= THISBACK(OnAverage);
	menuAdvanced.butAverage.Disable();
	//menuAdvanced.butUpdateCrot << THISBACK(OnUpdateCrot);
	menuAdvanced.butUpdateCrot.SetCtrl(menuAdvancedReference).Tip(t_("Click to change the centre of rotation"));
	menuAdvanced.c_array.AddColumn(t_("Body"));
	menuAdvanced.c_array.AddColumn(t_("x"));
	menuAdvanced.c_array.AddColumn(t_("y"));
	menuAdvanced.c_array.AddColumn(t_("z"));
	menuAdvanced.c_array.WhenLeftDouble = [&]() {menuAdvanced.butUpdateCrot.WhenAction();};
	menuAdvanced.butUpdateCwave << THISBACK(OnUpdateCwave);
	menuAdvanced.x_w.WhenEnter << THISBACK(OnUpdateCwave);
	menuAdvanced.y_w.WhenEnter << THISBACK(OnUpdateCwave);
	
	OnOpt();
	
	CtrlLayout(menuPlot);
	menuPlot.butZoomToFit << [&] {GetSelScatter().ZoomToFit(true, true);};
	menuPlot.autoFit	  << [&] {
		LoadSelTab(Bem());
		menuPlot.fromY0.Enable(~menuPlot.autoFit);
	};
	menuPlot.fromY0 	<< [&]{LoadSelTab(Bem());};
	menuPlot.opwT 		<< [&]{mainQTF.pf = Null;	LoadSelTab(Bem());};
	menuPlot.opMP 		<< [&]{LoadSelTab(Bem());};
	menuPlot.showPoints << [&]{LoadSelTab(Bem());};
	menuPlot.showNdim 	<< [&]{LoadSelTab(Bem());};
	menuPlot.opAinf 	<< [&]{LoadSelTab(Bem());};
	menuPlot.opA0 		<< [&]{LoadSelTab(Bem());};
	menuPlot.opB 		<< [&]{LoadSelTab(Bem());};
	menuPlot.opApot 	<< [&]{LoadSelTab(Bem());};
	menuPlot.opBhask 	<< [&]{LoadSelTab(Bem());};
	menuPlot.opBpot 	<< [&]{LoadSelTab(Bem());};
	menuPlot.opFfkpot 	<< [&]{LoadSelTab(Bem());};
	menuPlot.opFscpot 	<< [&]{LoadSelTab(Bem());};
		
	menuPlot.head1st.NoHeader().MultiSelect();
	menuPlot.head1st.AddColumn("");
	menuPlot.headMD.NoHeader().MultiSelect();
	menuPlot.headMD.AddColumn("");
	menuPlot.headMD.AddColumn("");
	menuPlot.butList.SetCtrl(menuPlotList).Tip(t_("Click to select headings"));
	menuPlotList.Init(*this, menuPlot.head1st, menuPlot.headMD, menuPlot.headQTF);
	
	CtrlLayout(menuBody);
	menuBody.butSpreadNegative 	<<= THISBACK(OnSpreadNegative);
	menuBody.butMapNodes 		<<= THISBACK(OnMapNodes);
	menuBody.butMapMeshes 		<<= THISBACK(OnMapMeshes);
	menuBody.butFill1st 		<<= THISBACK(OnFill1st);
	
	OnOpt();
	
	menuFOAMM.Init(*this, mainSetupFOAMM);
	
	OnOpt();
		
	menuTab.Add(menuOpen.SizePos(), 	t_("File"));
	menuTab.Add(menuPlot.SizePos(), 	t_("Plot")).Disable();
	menuTab.Add(menuProcess.SizePos(), 	t_("Process")).Disable();
	menuTab.Add(menuProcess2.SizePos(), t_("Remove & Mult")).Disable();
	menuTab.Add(menuAdvanced.SizePos(), t_("Advanced")).Disable();
	menuTab.Add(menuFOAMM.SizePos(), 	t_("FOAMM State Space")).Disable();
	menuTab.Add(menuBody.SizePos(), 	t_("Body")).Disable();
	
	menuTab.WhenSet = [&] {
		LOGTAB(menuTab);
		bool setupfoamm = false;
		if (menuTab.IsAt(menuFOAMM)) {
			setupfoamm = true;
			if (!FileExists(Bem().foammPath))
				Status(t_("FOAMM not found. Please set FOAMM path in Options"), 10000);	
		} else if (menuTab.IsAt(menuPlot) || menuTab.IsAt(menuOpen) || 
				   menuTab.IsAt(menuProcess) || menuTab.IsAt(menuProcess2) ||
				   menuTab.IsAt(menuAdvanced) || menuTab.IsAt(menuBody)) 
			setupfoamm = true;
		
		if (!setupfoamm) 
			mainTab.Set(0);
		
		if (menuTab.IsAt(menuFOAMM)) 
			mainTab.Set(mainSetupFOAMM);
		else if (menuTab.IsAt(menuBody)) 
			mainTab.Set(mainBody);
		
		ShowMenuPlotItems();
	};
	
	mainTab.WhenSet = [&] {
		LOGTAB(mainTab);
		UVector<int> idxs = ArrayModel_IndexsHydro(listLoaded);
		bool plot = true, convertProcess = true, ismenuFOAMM = false, ismesh = false;
		int is = -1;			// 0: 1st, 1: QTF, 2: MD

		if (idxs.IsEmpty())
			plot = convertProcess = false;
		else if (mainTab.IsAt(mainMatrixK)) 
			mainMatrixK.Load(Bem().hydros, idxs, ~menuPlot.showNdim);
		else if (mainTab.IsAt(mainMatrixK2)) 
			mainMatrixK2.Load(Bem().hydros, idxs, ~menuPlot.showNdim);
		else if (mainTab.IsAt(mainMatrixA))
			mainMatrixA.Load(Bem().hydros, idxs, ~menuPlot.showNdim);
		else if (mainTab.IsAt(mainMatrixM)) {
			plot = false;
			mainMatrixM.Load(Bem().hydros, idxs, false);
		} else if (mainTab.IsAt(mainMatrixDlin)) {
			plot = false;
			mainMatrixDlin.Load(Bem().hydros, idxs, false);
		} else if (mainTab.IsAt(mainMatrixDquad)) {
			plot = false;
			mainMatrixDquad.Load(Bem().hydros, idxs, false);
		} else if (mainTab.IsAt(mainA))
			mainA.Load(idxs);
		else if (mainTab.IsAt(mainB))
			mainB.Load(idxs);
		else if (mainTab.IsAt(mainMD)) {
			mainMD.Load(idxs, menuPlot.headMD.GetCursor());
			is = 2;
		} else if (mainTab.IsAt(mainK))
			mainK.Load(idxs);
		else if (mainTab.IsAt(mainAinfw))
			mainAinfw.Load(idxs);
		else if (mainTab.IsAt(mainForceSC)) {
			is = 0;
			mainForceSC.Load(idxs, menuPlot.head1st.GetCursor());
		} else if (mainTab.IsAt(mainForceFK)) {
			is = 0;
			mainForceFK.Load(idxs, menuPlot.head1st.GetCursor());
		} else if (mainTab.IsAt(mainForceEX)) {
			is = 0;
			mainForceEX.Load(idxs, menuPlot.head1st.GetCursor());
		} else if (mainTab.IsAt(mainRAO)) {
			is = 0;
			mainRAO.Load(idxs, menuPlot.head1st.GetCursor());
		} else if (mainTab.IsAt(mainStateSpace)) {
			mainStateSpace.Load(idxs);
			ismenuFOAMM = true;
		} else if (mainTab.IsAt(mainSetupFOAMM)) {
			menuTab.Set(menuFOAMM);
			ismenuFOAMM = true;
			if (mainSetupFOAMM.arrayCases.GetCount() == 0 && idxs.size() == 1) 
				listLoaded.SetCursor(0);
			menuFOAMM.OnCursor();
		} else if (mainTab.IsAt(mainQTF)) {
			mainQTF.Load();
			is = 1;
		} else if (menuTab.IsAt(menuFOAMM)) 
			;
		else if (mainTab.IsAt(mainBody)) {
			mainBody.Load();	
			plot = false;
		} else {
			plot = false;
		}
		
		for (int idx : idxs) {
			const Hydro &hy = Bem().hydros[idx];
			if (!hy.dt.msh.IsEmpty()) {
				ismesh = true;
				break;
			}
		}
			
		if (is != 1)
			mainQTF.Unload();
		
		menuPlot.labHead1st.Show(is == 0);	menuPlot.head1st.Show(is == 0);
		menuPlot.labHeadQTF.Show(is == 1);	menuPlot.headQTF.Show(is == 1);
		menuPlot.labHeadMD.Show(is == 2);	menuPlot.headMD.Show(is == 2);

		menuPlot.butList.Show(is >= 0);
		menuPlotList.head1st.Show(is == 0);
		menuPlotList.headQTF.Show(is == 1);
		menuPlotList.headMD.Show(is == 2);
		if (is == 0)
			menuPlotList.label.SetText(t_("Headings 1st"));
		else if (is == 1)
			menuPlotList.label.SetText(t_("Headings QTF"));
		else if (is == 2)
			menuPlotList.label.SetText(t_("Head. Mean Drift"));
			
		TabCtrl::Item& tabMenuPlot = menuTab.GetItem(menuTab.Find(menuPlot));
		tabMenuPlot.Enable(plot);
		TabCtrl::Item& tabMenuProcess = menuTab.GetItem(menuTab.Find(menuProcess));
		tabMenuProcess.Enable(convertProcess);
		TabCtrl::Item& tabMenuProcess2 = menuTab.GetItem(menuTab.Find(menuProcess2));
		tabMenuProcess2.Enable(convertProcess);
		TabCtrl::Item& tabMenuAdvanced = menuTab.GetItem(menuTab.Find(menuAdvanced));
		tabMenuAdvanced.Enable(convertProcess);
		TabCtrl::Item& tabMenuBody = menuTab.GetItem(menuTab.Find(menuBody));
		tabMenuBody.Enable(convertProcess);
		
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
			if (menuTab.IsAt(menuPlot) || menuTab.IsAt(menuProcess) || menuTab.IsAt(menuProcess2) || menuTab.IsAt(menuAdvanced))
				menuTab.Set(0);
			else if (menuTab.IsAt(menuBody))
				menuTab.Set(menuBody);
		}
		
		if (convertProcess) {
			tabMenuProcess.Text(t_("Process"));
			tabMenuProcess2.Text(t_("Remove"));
			tabMenuAdvanced.Text(t_("Advanced"));
		} else {
			tabMenuProcess.Text("");
			tabMenuProcess2.Text("");
			tabMenuAdvanced.Text("");
			if (menuTab.IsAt(menuProcess) || menuTab.IsAt(menuProcess2) || menuTab.IsAt(menuAdvanced))
				menuTab.Set(0);
		}
		
		if (ismesh) {
			tabMenuBody.Text(t_("Body"));
			tabMenuBody.Enable();
		} else {
			tabMenuBody.Text("");
			tabMenuBody.Disable();
			if (menuTab.IsAt(menuBody))
				menuTab.Set(0);
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
	mainTab.Add(mainMatrixK.SizePos(), t_("StiffHydr")).Disable();
		
	mainA.Init(Hydro::DATA_A);
	mainTab.Add(mainA.SizePos(), t_("A")).Disable();
	
	mainB.Init(Hydro::DATA_B);
	mainTab.Add(mainB.SizePos(), t_("B")).Disable();
	
	mainK.Init(Hydro::DATA_KIRF);
	mainTab.Add(mainK.SizePos(), t_("IRF")).Disable();
	
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

	mainMatrixM.Init(Hydro::MAT_M);
	mainTab.Add(mainMatrixM.SizePos(), t_("M")).Disable();
	
	mainMatrixDlin.Init(Hydro::MAT_DAMP_LIN);
	mainTab.Add(mainMatrixDlin.SizePos(), t_("DampLin")).Disable();
	
	mainMatrixDquad.Init(Hydro::MAT_DAMP_QUAD);
	mainTab.Add(mainMatrixDquad.SizePos(), t_("DampQuad")).Disable();

	mainMatrixK2.Init(Hydro::MAT_KMOOR);
	mainTab.Add(mainMatrixK2.SizePos(), t_("StiffMoor")).Disable();
	
	mainMD.Init(Hydro::DATA_MD);
	mainTab.Add(mainMD.SizePos(), t_("Mean Drift")).Disable();
		
	mainQTF.Init(*this);
	mainTab.Add(mainQTF.SizePos(), t_("QTF")).Disable();
	
	mainAinfw.Init(Hydro::DATA_AINFW);
	mainTab.Add(mainAinfw.SizePos(), t_("A∞(ω)")).Disable();
	
	mainSetupFOAMM.Init();
	mainTab.Add(mainSetupFOAMM.SizePos(), t_("FOAMM")).Disable();

	mainStateSpace.Init();
	mainTab.Add(mainStateSpace.SizePos(), t_("State Space")).Disable();
	
	mainBody.Init();
	mainTab.Add(mainBody.SizePos(), t_("Mesh & Potentials")).Disable();
	
	UpdateButtons();
	saveFolder = GetDesktopFolder();
}


void MainBEM::ShowMenuPlotItems() {
	menuPlot.showNdim.Enable();

	bool showScatter = true, showwT = true, showComplex = false, showDim = true, 
		 showA = false, showB = false, showFfk = false, showFsc = false;
		 
	if (mainTab.IsAt(mainSetupFOAMM)) {
		showwT = false;
		showComplex = true;
	} else if (mainTab.IsAt(mainK)) {
		showwT = false;
	} else if (mainTab.IsAt(mainMatrixK)) {
		showwT = false;
		showScatter = false;
	} else if (mainTab.IsAt(mainQTF)) {
		showComplex = true;
	} else if (mainTab.IsAt(mainForceSC) || mainTab.IsAt(mainForceFK) || mainTab.IsAt(mainForceEX) || 
			   mainTab.IsAt(mainRAO) || mainTab.IsAt(mainStateSpace)) {
		showComplex = true;
		if (mainTab.IsAt(mainForceFK))
			showFfk = true;
		else if (mainTab.IsAt(mainForceSC))
			showFsc = true;
	} else if (mainTab.IsAt(mainA))
		showA = true;
	else if (mainTab.IsAt(mainB))
		showB = true;
	else if (mainTab.IsAt(mainBody)) {
		showScatter = false;
	}
		
	menuPlot.showNdim.Enable(showDim);
	menuPlot.opwT.Enable(showwT);
	menuPlot.opMP.Show(showComplex);
	menuPlot.labMP.Show(showComplex || showA || showB);

	menuPlot.butZoomToFit.Enable(showScatter);
	menuPlot.autoFit.Enable(showScatter);
	menuPlot.fromY0.Enable(showScatter);
	menuPlot.showPoints.Enable(showScatter);

	menuPlot.opAinf.Show(showA);
	menuPlot.opA0.Show(showA);
	menuPlot.opB.Show(showA);
	menuPlot.opApot.Show(showA);
	
	menuPlot.opBhask.Show(showB);
	menuPlot.opBpot.Show(showB);
	
	menuPlot.opFfkpot.Show(showFfk);
	menuPlot.opFscpot.Show(showFsc);
}

void MainBEM::OnMenuAdvancedArraySel(bool updateBH) {
	int idx = ArrayModel_IndexHydro(listLoaded);
	if (idx < 0)
		return;
	
	Hydro &hy = Bem().hydros[idx];
	
	menuAdvanced.x_w = hy.dt.x_w;
	menuAdvanced.y_w = hy.dt.y_w;
	
	menuAdvanced.numBH.Enable(hy.dt.Nh > 0);
	if (updateBH)
		menuAdvanced.numBH = hy.dt.Nh;
	menuAdvanced.labNum.SetText(Format("/\%d", hy.dt.Nh));
	
	menuAdvanced.c_array.Clear();
	for (int ib = 0; ib < hy.dt.Nb; ++ib)
		menuAdvanced.c_array.Add		 (Format("%d.%s", ib+1, hy.dt.msh[ib].dt.name), hy.dt.msh[ib].dt.c0.x, hy.dt.msh[ib].dt.c0.y, hy.dt.msh[ib].dt.c0.z);			
	
	menuAdvancedReference.c_array.Clear();
	for (int ib = 0; ib < hy.dt.Nb; ++ib)
		menuAdvancedReference.c_array.Add(Format("%d.%s", ib+1, hy.dt.msh[ib].dt.name), hy.dt.msh[ib].dt.c0.x, hy.dt.msh[ib].dt.c0.y, hy.dt.msh[ib].dt.c0.z);
	menuAdvancedReference.Init(*this, idx);
		
	menuAdvanced.labelBodyAxis.SetLabel(Format(t_("Body Axis (%d)"), hy.dt.Nb)); 
	
	menuAdvanced.butUpdateCrot.Enable(hy.dt.Nb > 0);
	menuAdvanced.c_array.Enable(hy.dt.Nb > 0);
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
	
	if (!ret || IsNull(menuPlot.opAinf)) 
		menuPlot.opAinf = true;	
	if (!ret || IsNull(menuPlot.opA0)) 
		menuPlot.opA0 = true;
	if (!ret || IsNull(menuPlot.opApot)) 
		menuPlot.opApot = true;
	if (!ret || IsNull(menuPlot.opB)) 
		menuPlot.opB = false;
	
	if (!ret || IsNull(menuPlot.opBhask)) 
		menuPlot.opBhask = true;
	if (!ret || IsNull(menuPlot.opBpot)) 
		menuPlot.opBpot = true;
	
	if (!ret || IsNull(menuPlot.opFfkpot)) 
		menuPlot.opFfkpot = true;
	
	if (!ret || IsNull(menuPlot.opFscpot)) 
		menuPlot.opFscpot = true;
	
	if (!ret || IsNull(menuPlot.showPoints)) 
		menuPlot.showPoints = true;
	
	if (!ret || IsNull(menuPlot.showNdim)) 
		menuPlot.showNdim = false;
}

void MainBEM::LoadSelTab(BEM &bem) {
	UVector<int> idxs = ArrayModel_IndexsHydro(listLoaded);
	int id = mainTab.Get();
	if (id == mainTab.Find(mainStateSpace))
		mainStateSpace.Load(idxs);
	else if (id == mainTab.Find(mainMatrixK))
		mainMatrixK.Load(bem.hydros, idxs, ~menuPlot.showNdim);
	else if (id == mainTab.Find(mainMatrixK2))
		mainMatrixK2.Load(bem.hydros, idxs, ~menuPlot.showNdim);
	else if (id == mainTab.Find(mainMatrixA))
		mainMatrixA.Load(bem.hydros, idxs, ~menuPlot.showNdim);
	else if (id == mainTab.Find(mainMatrixM))
		mainMatrixM.Load(bem.hydros, idxs, false);
	else if (id == mainTab.Find(mainMatrixDlin))
		mainMatrixDlin.Load(bem.hydros, idxs, false);
	else if (id == mainTab.Find(mainMatrixDquad))
		mainMatrixDquad.Load(bem.hydros, idxs, false);
	else if (id == mainTab.Find(mainSummary))
		;
	else if (id == mainTab.Find(mainSetupFOAMM))
		mainSetupFOAMM.Load();
	else if (id == mainTab.Find(mainQTF))
		mainQTF.Load();
	else if (id == mainTab.Find(mainBody))
		mainBody.Load();
	else 
		GetSelABForce().Load(idxs);
	
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
			if (ForceExt(Bem().hydros[i].dt.file, ".") == ForceExt(file, ".") &&
				(Bem().GetBEMExtSet(file) < 0 || 
				 Bem().GetBEMExtSet(file) == Bem().GetBEMExtSet(Bem().hydros[i].dt.file))) {
				if (!PromptYesNo(t_("Model is already loaded") + S("&") + t_("Do you wish to open it anyway?")))
					return false;
				break;
			}
		}
		
		WaitCursor wait;
		
		int num = Bem().LoadBEM(file, [&](String str, int _pos) {
			progress.SetText(str); 
			progress.SetPos(_pos); 
			return !progress.Canceled();
		}, false);
		
		//int id = Bem().hydros.size()-1;
		for (int idx = Bem().hydros.size() - num; idx < Bem().hydros.size(); ++idx) {
			Hydro &hy = Bem().hydros[idx];
		
			hy.Report();
			mainSummary.Report(hy, idx);
			if (hy.dt.Nf < 0)
				return false;
			
			ArrayModel_Add(listLoaded, hy.GetCodeStr(), hy.dt.name, hy.dt.file, hy.dt.GetId());
		}
		
		//UpdateButtons();

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
			
			menuPlotList.Set1st();
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
			
			menuPlotList.SetMD();
		}
		
		mainTab.WhenSet();
	} catch (Exc e) {
		if (!e.IsEmpty())
			BEM::PrintError(t_("Error: ") + e);
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
			int id = ArrayModel_IndexHydro(listLoaded, row);
			Bem().RemoveHydro(id);
			listLoaded.Remove(row);
			selected = true;
		}
	}	// Only one available => directly selected
	if (!selected && listLoaded.GetCount() == 1) {	
		int id = ArrayModel_IndexHydro(listLoaded, 0);
		Bem().RemoveHydro(id);
		listLoaded.Remove(0);
		selected = true;
	}		
	if (!selected) {
		BEM::PrintError(t_("No model selected"));
		return;
	}
	//UpdateButtons();
	
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
	menuOpen.butRename.Enable(numsel == 1 || numrow == 1);
	menuOpen.butSolve.Enable(numsel == 1 || numrow == 1);
	menuOpen.dropExport.Enable(numsel == 1);
	menuOpen.butExport.Enable(numsel == 1);

	menuProcess.butSym.			Enable((numsel == 1 || numrow == 1) && (int)menuProcess.dropSym.GetData() > 0);
	
	menuProcess.butKirf.		Enable(numsel == 1 || numrow == 1);
	menuProcess.butA0.			Enable(numsel == 1 || numrow == 1);
	menuProcess.butAinf.		Enable(numsel == 1 || numrow == 1);
	menuProcess.butRAO.			Enable(numsel == 1 || numrow == 1);
	menuProcess.critDamp.		Enable(numsel == 1 || numrow == 1);
	menuAdvanced.butAinfw.		Enable(numsel == 1 || numrow == 1);
	menuAdvanced.butBH.			Enable(numsel >= 1);
	menuAdvanced.butOgilvie.	Enable(numsel == 1 || numrow == 1);
	menuAdvanced.butConvergence.Enable(numsel >= 3);
	menuAdvanced.butAverage.	Enable(numsel >= 2);
	menuAdvanced.opDecayingTail.Enable(numsel == 1 || numrow == 1);
	menuAdvanced.opThinremoval. Enable(numsel == 1 || numrow == 1);
	menuAdvanced.opZremoval.	Enable(numsel == 1 || numrow == 1);

	{	
		int row = menuPlot.head1st.GetCursor();
		double h;
		if (row >= 0)
			h = menuPlot.head1st.Get(row, 0);
		
		menuPlot.head1st.Clear();
		for (int ih = 0; ih < Bem().headAll.size(); ++ih) {
			if (Bem().headAll[ih] == h)
				row = ih;
			menuPlot.head1st.Add(Bem().headAll[ih]);
		}
		if (row >= 0)
			menuPlot.head1st.SetCursor(row);
		
		menuPlot.head1st.WhenCursor = [&] {
			int row = menuPlot.head1st.GetCursor();
			if (row < 0)
				return;
			
			menuPlotList.head1st.SetCursor(row);
			
			double h = menuPlot.head1st.Get(row, 0);
			int ih = FindClosest(Bem().headAll, h); 
			
			UVector<int> idxs = ArrayModel_IndexsHydro(listLoaded);
			
			mainForceSC.Load(idxs, ih);
			mainForceFK.Load(idxs, ih);
			mainForceEX.Load(idxs, ih);
			mainRAO.Load(idxs, ih);
		};
		menuPlot.head1st.WhenLeftDouble = [&] {menuPlot.butList.WhenAction();};
	}
	{
		int row = menuPlot.headMD.GetCursor();
		std::complex<double> h;
		if (row >= 0) {
			h.real(menuPlot.headMD.Get(row, 0));
			h.imag(menuPlot.headMD.Get(row, 1));
		}
		menuPlot.headMD.Clear();
		for (int ih = 0; ih < Bem().headAllMD.size(); ++ih) {
			if (h == Bem().headAllMD[ih])
				row = ih;
			menuPlot.headMD.Add(Bem().headAllMD[ih].real(), Bem().headAllMD[ih].imag());
		}
		if (row >= 0)
			menuPlot.headMD.SetCursor(row);
		
		menuPlot.headMD.WhenCursor = [&] {
			int row = menuPlot.headMD.GetCursor();
			if (row < 0)
				return;
			
			menuPlotList.headMD.SetCursor(row);
			
			std::complex<double> h(menuPlot.headMD.Get(row, 0), menuPlot.headMD.Get(row, 1));
			int ih = FindClosest(Bem().headAllMD, h); 
			
			UVector<int> idxs = ArrayModel_IndexsHydro(listLoaded);
		
			mainMD.Load(idxs, ih);
		};	
		menuPlot.headMD.WhenLeftDouble = [&] {menuPlot.butList.WhenAction();};	
	}
	
	bool show_w = menuPlot.opwT == 0;
	//bool show_ma_ph = menuPlot.opMP == 0;
	
	menuProcess.dropBody1.Clear();
	menuProcess.dropBody2.Clear();
	
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
	
	int idx = GetIndexOneSelected(false);
	if (idx >= 0) { 
		Hydro &hy = Bem().hydros[idx];
		
		for (int i = 0; i < hy.dt.Nb; ++i)
			menuProcess.dropBody1.Add(i+1);
		menuProcess.dropBody1.SetIndex(0);
		
		if (hy.dt.Nb > 1) {
			for (int i = 0; i < hy.dt.Nb; ++i)
				menuProcess.dropBody2.Add(i+1);
			menuProcess.dropBody2.SetIndex(0);
		}
		menuProcess.dropBody2.Enable(hy.dt.Nb > 1);
		
		VectorXd TT	= hy.Get_T();			
		for (int i = 0; i < hy.dt.w.size(); ++i)
			menuProcess2.dropFreq.Add(false, show_w ? hy.dt.w[i] : TT[i]);
		for (int i = 0; i < hy.dt.qw.size(); ++i)
			menuProcess2.dropFreqQTF.Add(false, show_w ? hy.dt.qw[i] : 2*M_PI/hy.dt.qw[i]);
		for (int i = 0; i < hy.dt.head.size(); ++i)
			menuProcess2.dropHead.Add(false, hy.dt.head[i]);
		for (int i = 0; i < hy.dt.mdhead.size(); ++i)
			menuProcess2.dropHeadMD.Add(false, Format("%.1f-%.1f", hy.dt.mdhead[i].real(), hy.dt.mdhead[i].imag()));
		for (int i = 0; i < hy.dt.qhead.size(); ++i)
			menuProcess2.dropHeadQTF.Add(false, Format("%.1f-%.1f", hy.dt.qhead[i].real(), hy.dt.qhead[i].imag()));
		for (int i = 0; i < 6; ++i)
			menuProcess2.dropDOF.Add(false, BEM::StrDOF(i));
		
		if (hy.IsLoadedFex())
			menuProcess2.dropForce.Add(false, t_("All"));
		if (hy.IsLoadedFsc())
			menuProcess2.dropForce.Add(false, t_("Diffraction"));
		if (hy.IsLoadedFfk())
			menuProcess2.dropForce.Add(false, t_("Froude-Krylov"));

		if (hy.IsLoadedMD())
			menuProcess2.dropForceMD.Add(false, t_("All"));

		if (hy.IsLoadedQTF(true) || hy.IsLoadedQTF(false))
			menuProcess2.dropForceQTF.Add(false, t_("All"));
		if (hy.IsLoadedQTF(true))
			menuProcess2.dropForceQTF.Add(false, t_("Summation"));
		if (hy.IsLoadedQTF(false))
			menuProcess2.dropForceQTF.Add(false, t_("Difference"));
	}
	OnMenuAdvancedArraySel(false);
}

void MainBEM::OnJoin() {
	UVector<int> idsjoin, rowsJoin;
	for (int row = listLoaded.GetCount()-1; row >= 0; --row) {
		if (listLoaded.IsSelected(row)) {
			rowsJoin << row;
			int idx = ArrayModel_IndexHydro(listLoaded, row);
			idsjoin << idx;
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
		
		Hydro &hy = Bem().Join(idsjoin, [&](String str, int _pos) {
			progress.SetText(str); 
			progress.SetPos(_pos); 
			return !progress.Canceled();
		});
		
		mainSummary.Clear();
		
		ArrayModel_Add(listLoaded, hy.GetCodeStr(), hy.dt.name, hy.dt.file, hy.dt.GetId());
		ArrayModel_RowsHydroDel(listLoaded, rowsJoin);
	
		for (int idx = 0; idx < Bem().hydros.size(); ++idx) {
			const Hydro &hhy = Bem().hydros[idx];
			mainSummary.Report(hhy, idx);
		}
			
		if (Bem().hydros.size() > 0) 
			listLoaded.SetCursor(0);

	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
	AfterBEM();		
}

void MainBEM::OnDuplicate() {
	try {
		int idx = GetIndexOneSelected();
		if (idx < 0) 
			return;

		Hydro &hy = Bem().Duplicate(idx);
		
		mainSummary.Clear();
		
		ArrayModel_Add(listLoaded, hy.GetCodeStr(), hy.dt.name, hy.dt.file, hy.dt.GetId());

		AfterBEM();
		
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
}

void MainBEM::OnSymmetrizeForces() {
	try {
		int idx = GetIndexOneSelected();
		if (idx < 0) 
			return;

		WaitCursor wait;

		Progress progress(t_("Symmetrizing forces and RAOs in selected BEM file..."), 100); 
		
		if (~menuProcess.dropSym == BasicBEM::SYM_YZ || ~menuProcess.dropSym == BasicBEM::SYM_XZ_YZ)
			Bem().SymmetrizeForces(idx, false);
		if (~menuProcess.dropSym == BasicBEM::SYM_XZ || ~menuProcess.dropSym == BasicBEM::SYM_XZ_YZ)
			Bem().SymmetrizeForces(idx, true);
		
		AfterBEM();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
}

void MainBEM::OnKirfAinf(Hydro::DataToPlot param) {
	try {
		int idx = GetIndexOneSelected();
		if (idx < 0) 
			return;
			
		Progress progress(Format(t_("Calculating %s in selected BEM file..."), Hydro::StrDataToPlot(param)), 100); 
		
		double maxT = Null;
		
		if (param == Hydro::PLOT_KIRF || 
		   ((param == Hydro::PLOT_AINF || param == Hydro::PLOT_AINFW) 
		   	 && !Bem().hydros[idx].IsLoadedKirf())) {
		 	maxT = Bem().hydros[idx].GetK_IRF_MaxT();
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
			Bem().A0(idx);
		else if (param == Hydro::PLOT_KIRF)
			Bem().Kirf(idx, maxT);
		else if (param == Hydro::PLOT_AINF) {
			if (!Bem().hydros[idx].IsLoadedKirf())
				Bem().Kirf(idx, maxT);
			Bem().Ainf(idx);
		} else if (param == Hydro::PLOT_AINFW) {
			if (!Bem().hydros[idx].IsLoadedKirf())
				Bem().Kirf(idx, maxT);
			if (!Bem().hydros[idx].IsLoadedAinf())
				Bem().Ainf(idx);
			Bem().Ainf_w(idx);
		}
		
		AfterBEM();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
}

void MainBEM::OnBH(int num) {
	try {
		int idx = GetIndexOneSelected();
		if (idx < 0) 
			return;

		WaitCursor wait;
		
		Bem().BH(idx, num);
			
		menuAdvanced.numBH <<= num;
			
		AfterBEM();
		
		mainTab.Set(mainTab.Find(mainB));
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
}

void MainBEM::OnRAO() {
	try {
		int idx = GetIndexOneSelected();
		if (idx < 0) 
			return;

		double critDamp = ~menuProcess.critDamp;
		if (IsNull(critDamp))
			critDamp = 0;
			
		Progress progress(t_("Calculating RAO in selected BEM file..."), 100); 
		
		WaitCursor wait;

		Bem().RAO(idx, critDamp);
		
		AfterBEM();	
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
}	

void MainBEM::OnSymmetrize() {
	try {
		int idx = GetIndexOneSelected();
		if (idx < 0) 
			return;
			
		WaitCursor wait;
		
		Bem().Symmetrize(idx);
				
		AfterBEM();

	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MainBEM::OnOgilvie() {
	String str;
	try {
		int idx = GetIndexOneSelected();
		if (idx < 0) 
			return;
			
		Progress progress(t_("Calculating Ogilvie compliance in selected BEM file..."), 100); 
		
		WaitCursor wait;
		
		//double maxT = Null;
		UVector<int> vidof, vjdof;
		Bem().OgilvieCompliance(idx, ~menuAdvanced.opZremoval, ~menuAdvanced.opThinremoval, ~menuAdvanced.opDecayingTail, vidof, vjdof);
				
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
				int id = ArrayModel_IndexHydro(listLoaded, row);
				ids << id;
			}
		}		
		if (ids.size() < 2) {
			BEM::PrintError(t_("Not enough models selected"));
			return;
		}
			
		Progress progress(t_("Calculating average..."), 100); 
		
		WaitCursor wait;
		
		Hydro &hy = Bem().Average(ids);
		
		mainSummary.Clear();
		
		ArrayModel_Add(listLoaded, hy.GetCodeStr(), hy.dt.name, hy.dt.file, hy.dt.GetId());
				
		AfterBEM();

	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MainBEM::OnConvergence() {
	String str;
	try {
		UVector<int> idxs;
		for (int row = 0; row < listLoaded.GetCount(); ++row) {
			if (listLoaded.IsSelected(row)) {
				int idx = ArrayModel_IndexHydro(listLoaded, row);
				idxs << idx;
			}
		}		
		if (idxs.size() < 3) {
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
		int idx = GetIndexOneSelected();
		if (idx < 0) 
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
		
		Bem().ResetForces(idx, force, forceMD, forceQtf);
				
		mainSummary.Clear();
		for (int i = 0; i < Bem().hydros.size(); ++i)
			mainSummary.Report(Bem().hydros[i], i);
		
		UVector<int> idxs = ArrayModel_IndexsHydro(listLoaded);
		
		mainTab.GetItem(mainTab.Find(mainForceSC)).Enable(mainForceSC.Load(idxs, menuPlot.head1st.GetCursor()));
		mainTab.GetItem(mainTab.Find(mainForceFK)).Enable(mainForceFK.Load(idxs, menuPlot.head1st.GetCursor()));
		mainTab.GetItem(mainTab.Find(mainForceEX)).Enable(mainForceEX.Load(idxs, menuPlot.head1st.GetCursor()));
		mainTab.GetItem(mainTab.Find(mainRAO)).Enable(mainRAO.Load(idxs));
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
		int idx = GetIndexOneSelected();
		if (idx < 0) 
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
		
		Bem().MultiplyDOF(idx, factor, idDOF, a, b, diag, menuProcess2.opF, menuProcess2.opMD, menuProcess2.opQTF, menuProcess2.opStiffness);
				
		AfterBEM();	
		LoadSelTab(Bem());
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MainBEM::OnSwapDOF() {
	try {
		int idx = GetIndexOneSelected();
		if (idx < 0) 
			return;
		
		WaitCursor wait;
		
		Bem().SwapDOF(idx, menuProcess.dropBody1.GetIndex(), menuProcess.dropDOF1.GetIndex() - 1, 
						  menuProcess.dropBody2.GetIndex(), menuProcess.dropDOF2.GetIndex() - 1);
				
		mainSummary.Clear();
		for (int i = 0; i < Bem().hydros.size(); ++i)
			mainSummary.Report(Bem().hydros[i], i);
		
		UVector<int> idxs = ArrayModel_IndexsHydro(listLoaded);
		
		mainTab.GetItem(mainTab.Find(mainA)).Enable(mainA.Load(idxs));
		mainTab.GetItem(mainTab.Find(mainB)).Enable(mainB.Load(idxs));
		mainTab.GetItem(mainTab.Find(mainMD)).Enable(mainMD.Load(idxs, menuPlot.headMD.GetCursor()));
		mainTab.GetItem(mainTab.Find(mainK)).Enable(mainK.Load(idxs));
		mainTab.GetItem(mainTab.Find(mainAinfw)).Enable(mainAinfw.Load(idxs));
		mainTab.GetItem(mainTab.Find(mainForceSC)).Enable(mainForceSC.Load(idxs, menuPlot.head1st.GetCursor()));
		mainTab.GetItem(mainTab.Find(mainForceFK)).Enable(mainForceFK.Load(idxs, menuPlot.head1st.GetCursor()));
		mainTab.GetItem(mainTab.Find(mainForceEX)).Enable(mainForceEX.Load(idxs, menuPlot.head1st.GetCursor()));
		mainTab.GetItem(mainTab.Find(mainRAO)).Enable(mainRAO.Load(idxs, menuPlot.head1st.GetCursor()));	
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
		int idx = GetIndexOneSelected();
		if (idx < 0) 
			return;
		
		WaitCursor wait;
		
		Bem().FillFrequencyGapsABForces(idx, ~menuProcess.opFill == 0, ~menuProcess.maxFreq);
				
		AfterBEM();
		LoadSelTab(Bem());
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MainBEM::OnQTF() {
	try {
		int idx = GetIndexOneSelected();
		if (idx < 0) 
			return;
		
		WaitCursor wait;
		
		Bem().FillFrequencyGapsQTF(idx, ~menuProcess.opFill == 0, ~menuProcess.maxFreq);
				
		AfterBEM();
		LoadSelTab(Bem());
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MainBEM::OnABForcesZero() {
	try {
		int idx = GetIndexOneSelected();
		if (idx < 0) 
			return;
		
		WaitCursor wait;
		
		Bem().FillFrequencyGapsABForcesZero(idx);
		
		AfterBEM();
		LoadSelTab(Bem());
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MainBEM::OnQTFZero() {
	try {
		int id = GetIndexOneSelected();
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
		int idx = GetIndexOneSelected();
		if (idx < 0) 
			return;
		
		WaitCursor wait;
		
		Bem().CopyQTF_MD(idx);
		
		Bem().UpdateHeadAllMD();
		
		AfterBEM();
		LoadSelTab(Bem());
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MainBEM::OnDeleteHeadingsFrequencies() {
	try {
		int id = GetIndexOneSelected();
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

void MainBEM::OnUpdateCwave() {
	try {
		int idx = GetIndexOneSelected();
		if (idx < 0) 
			return;
		
		WaitCursor wait;
		
		Bem().WaveTo(idx, ~menuAdvanced.x_w, ~menuAdvanced.y_w);
				
		AfterBEM();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MainBEM::OnSpreadNegative() {
	try {
		int id = GetIndexOneSelected();
		if (id < 0) 
			return;
		
		Progress progress(t_("Spreading negative values..."), 100); 
		
		String errors = Bem().SpreadNegative(id, [&](String str, int _pos) {
			progress.SetText(str); 
			progress.SetPos(_pos); 
			return !progress.Canceled();
		});
		
		progress.Close();
		
		if (!errors.IsEmpty())
			Exclamation(t_("Some negative panels cannot be spread in:&") + DeQtfLf(errors));
				
		AfterBEM();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}


void MainBEM::OnMapNodes() {
	try {
		int idx = GetIndexOneSelected();
		if (idx < 0) {
			Exclamation(t_("Please select a file"));
			return;
		}
		
		const Hydro &hy = Bem().hydros[idx];
		if (hy.dt.msh.IsEmpty() || hy.dt.msh[0].dt.mesh.panels.IsEmpty()) {
			Exclamation(t_("No mesh is available"));
			return;
		}
		
		int ib = mainBody.GetIb();
		
		if (!hy.IsLoadedPotsRad(ib)) {
			Exclamation(Format(t_("No radiation potentials/pressures are available for body %d"), ib+1));
			return;
		}
		
		mapNodes.Init(idx, ib);
		mapNodes.Execute();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
}

void MapNodes::Init(int _idx, int _ib) {
	idx = _idx;
	ib = _ib;
	
	CtrlLayout(*this);
	Title(Format(t_("Paste points and map to them the %d%s body mesh properties"), ib+1, Ordinal(ib+1)));
	
	butClose <<= THISBACK(OnClose);
	butPaste <<= THISBACK(OnPasteNodes);
	butMap <<= THISBACK(OnMapNodes);
	butExport <<= THISBACK(OnExport);
	
	const Hydro &hy = Bem().hydros[idx];
	
	dropFreq.Clear();
	for (int ifr = 0; ifr < hy.dt.Nf; ++ifr)
		dropFreq.Add(ifr, hy.dt.w[ifr]);
	dropFreq.SetIndex(0);
	
	dropExport.Clear();
	dropExport.Add(".csv").Add(".xlsx");
	dropExport.SetIndex(1);
	
	dropFreq.WhenAction = [&] {RefreshTable();};
	
	RefreshTable();
}

void MapNodes::RefreshTable() {
	arrayNodes.Reset();
	arrayNodes.SetLineCy(EditField::GetStdHeight()).MultiSelect().HeaderObject().Absolute();
	arrayNodes.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, arrayNodes, true);};
	
	const Hydro &hy = Bem().hydros[idx];
	Grid g;
	hy.SaveMap(g, dropFreq.GetData(), Bem().onlyDiagonal, ids, points, A_pan, B_pan);
	ArrayCtrlFill(arrayNodes, g, false);
	
	numNodes <<= arrayNodes.GetCount();
}	
		
void MapNodes::OnPasteNodes() {
	String str = ReadClipboardText();
	
	CSVParameters par;
	StringStream sin(str);
	if (!GuessCSVStream(sin, true, par)) {
		Exclamation(t_("Problem pasting nodes coordinates"));
		return;
	}
	if (par.parameters.size() < 3 || par.parameters.size() > 4) {
		Exclamation(t_("Incorrect number of columns"));
		return;
	}
	
	WaitCursor wait;
	
	ids.Clear();
	points.Clear();
	A_pan = Tensor<double, 4>();
	B_pan = Tensor<double, 4>();	
	
	sin.Seek(par.beginData);

	const char *endptr;	
	for (int row = 0; !sin.IsEof(); ++row) {
		UVector<String> data = Split(sin.GetLine(), par.separator, par.repetition);
		if (data.size() >= 3) {
			Point3D &p = points.Add();
			if (par.parameters.size() == 3) {
				ids << row+1;
				for (int c = 0; c < min(3, data.size()); ++c) 
					p[c] = ScanDouble(data[c], &endptr, par.decimalSign == ',');
			} else {
				ids << ScanInt(data[0]);
				for (int c = 0; c < min(3, data.size()); ++c) 
					p[c] = ScanDouble(data[c+1], &endptr, par.decimalSign == ',');
			}
		}
	}
	
	RefreshTable();
}

void MapNodes::OnMapNodes() {
	const Hydro &hy = Bem().hydros[idx];
		
	hy.MapNodes(ib, points, A_pan, B_pan);
	
	RefreshTable();
}

void MapNodes::OnExport(){
	String fileType = ~dropExport;
	
	FileSel fs;
	
	if (fileType == ".csv")
		fs.Type("Comma-separated values", "*.csv");
	else if (fileType == ".xlsx")
		fs.Type("Spreadsheet", "*.xlsx");
	fs.Type("Any file", "*.*");
	fs.ActiveType(0);
	fs.ActiveDir(GetDesktopFolder());
	
	if (!fs.ExecuteSaveAs(Format(t_("Save nodes data as %s"), fileType)))
		return;
	
	WaitCursor wait;
	
	String fileName = ~fs;

	int freqId = Null;
	if (opSaveAll)	
 		freqId = dropFreq.GetData();
	
	const Hydro &hy = Bem().hydros[idx];
	
	hy.SaveMap(fileName, fileType, freqId, Bem().onlyDiagonal, ids, points, A_pan, B_pan);
}

void MainBEM::OnMapMeshes() {
	try {
		int idx = GetIndexOneSelected();
		if (idx < 0) {
			Exclamation(t_("Please select a file"));
			return;
		}
		
		const Hydro &hy = Bem().hydros[idx];
		if (hy.dt.msh.IsEmpty() || hy.dt.msh[0].dt.mesh.panels.IsEmpty()) {
			Exclamation(t_("No mesh is available"));
			return;
		}
		
		int ib = mainBody.GetIb();
		
		if (!hy.IsLoadedPotsRad(ib)) {
			Exclamation(Format(t_("No radiation potentials/pressures are available for body %d"), ib+1));
			return;
		}
		
		mapMeshes.Init(idx, ib);
		mapMeshes.Execute();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
}

void MainBEM::OnFill1st() {
	try {
		int idx = GetIndexOneSelected();
		if (idx < 0) {
			Exclamation(t_("Please select a file"));
			return;
		}
		
		Hydro &hy = Bem().hydros[idx];
		hy.FillWithPotentials();
		
		AfterBEM();
		LoadSelTab(Bem());
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MapMeshes::Init(int _idx, int _ib) {
	idx = _idx;
	ib = _ib;
	
	CtrlLayout(*this);
	Title(Format(t_("Select meshes and map the %d%s body properties"), ib+1, Ordinal(ib+1)));
	text.Background(Null);
	int dots = 7*StdFont().GetHeight();
	text.SetQTF(Format(t_("[A+%d This option allows the hydrodynamic coefficients to be divided into a series of sectional bodies defined by their meshes.&"
				"These meshes are obtained by dividing the mesh of the calculated body into parts.&"
				"Therefore, the meshes of each sectional body do not have to be closed underwater.]"), dots));
	
	butClose <<= THISBACK(OnClose);
	butMapMeshes <<= THISBACK(OnMapMeshes);
	opOneMany = 0;
	
	ArrayModel_Init(listLoaded, true, {10, 10, 100}); 
	for (int i = 0; i < Bem().surfs.size(); ++i) {
		Body &msh = Bem().surfs[i];
		ArrayModel_Add(listLoaded, msh.GetBodyStr(), msh.dt.name, msh.dt.fileName, msh.dt.GetId(),
					optionsPlot, Null);
	}
}

void MapMeshes::OnMapMeshes() {
	try {
		UVector<int> idmeshes;
		for (int row = 0; row < listLoaded.GetCount(); ++row) {
			if (ArrayModel_IsVisible(listLoaded, row)) {
				int idxx = ArrayModel_IndexBody(listLoaded, row);
				if (idxx < 0)
					throw Exc(t_("Unexpected problem in OnMapMeshes()"));
				idmeshes << idxx;
			}
		}
		WaitCursor wait;
		
		int idFrom = Bem().hydros.size();
		Bem().MapMeshes(idx, ib, idmeshes, int(~opOneMany) == 0);
		
		for (int i = idFrom; i < Bem().hydros.size(); ++i) {
			const Hydro &hy = Bem().hydros[i];
			ArrayModel_Add(Ma().mainBEM.listLoaded, hy.GetCodeStr(), hy.dt.name, hy.dt.file, hy.dt.GetId());
		}
		
		Ma().mainBEM.mainSummary.Clear();
		for (int i = 0; i < Bem().hydros.size(); ++i) {
			const Hydro &hy = Bem().hydros[i];
			Ma().mainBEM.mainSummary.Report(hy, i);
		}
		
		Ma().mainBEM.AfterBEM();	
		Ma().mainBEM.LoadSelTab(Bem());
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
	Close();
}
	
MenuAdvancedReference::MenuAdvancedReference() {
	CtrlLayout(*this);
	
	c_array.Editing().Clipboard().Sorting(false);
	c_array.AddColumn(t_("Body"));
	c_array.AddColumn(t_("x")).Edit(edit[0]);
	c_array.AddColumn(t_("y")).Edit(edit[1]);
	c_array.AddColumn(t_("z")).Edit(edit[2]);
	
	cancel << [&] {Close();}; 
	ok 	   << [&] {
		try {
			Hydro &hy = Bem().hydros[idx];
			
			if (c_array.GetRowCount() != hy.dt.Nb)
				throw Exc(t_("Number of centres is different than the number of bodies"));
		
			MatrixXd to(3, hy.dt.Nb);
			for (int ib = 0; ib < hy.dt.Nb; ++ib) 
				for (int idf = 0; idf < 3; ++idf) 
					to(idf, ib) = c_array.Get(ib, idf+1);
			
			bool justSetCentres	= false;
			for (int ib = 0; ib < hy.dt.Nb; ++ib) {
				if (IsNull(hy.dt.msh[ib].dt.c0)) {
					PromptOK(t_("No translation will be done. New centres will be set."));
					justSetCentres = true;
					break;
				}
			}
			if (justSetCentres) {
				for (int ib = 0; ib < hy.dt.Nb; ++ib) 
					for (int idf = 0; idf < 3; ++idf) 
						hy.dt.msh[ib].dt.c0[idf] = c_array.Get(ib, idf+1);
			} else {
				WaitCursor wait;
		
				Progress progress(t_("Translating data..."), 100);
				
				Bem().TranslationTo(idx, to, [&](String str, int _pos) {
					progress.SetText(str); 
					progress.SetPos(_pos); 
					return !progress.Canceled();
				});
			}
			mbem->AfterBEM();
			
			Close();
		} catch (Exc e) {
			BEM::PrintError(DeQtfLf(e));
		}
	};
}

MenuPlotList::MenuPlotList() {
	CtrlLayout(*this);
}

void MainBEM::AfterBEM() {
	mainSummary.Clear();
	for (int idx = 0; idx < Bem().hydros.size(); ++idx) {
		Hydro &hy = Bem().hydros[idx];
		mainSummary.Report(hy, idx);
	}
	
	UVector<Point3D> c0;
	UVector<Pointf> c;
	for (int idx = 0; idx < Bem().hydros.size(); ++idx) {
		for (int ib = 0; ib < Bem().hydros[idx].dt.msh.size(); ++ib) 
			c0 << Bem().hydros[idx].dt.msh[ib].dt.c0;
		c << Pointf(Bem().hydros[idx].dt.x_w, Bem().hydros[idx].dt.y_w);
	}
	bool sameSystem = true;
	if (c.size() > 0) {
		for (int i = 1; i < c.size(); ++i) {
			if (!IsNull(c[i]) &&
				(!EqualDecimals(c[i].x, c[i-1].x, 3) || 
				 !EqualDecimals(c[i].y, c[i-1].y, 3))) {
				sameSystem = false;
				break;
			}
		}
	}
	if (sameSystem && c0.size() > 1) {
		for (int i = 1; i < c0.size(); ++i) {
			if (!IsNull(c0[i]) && 
				(!EqualDecimals(c0[i].x, c0[i-1].x, 3) || 
				 !EqualDecimals(c0[i].y, c0[i-1].y, 3) || 
				 !EqualDecimals(c0[i].z, c0[i-1].z, 3))) {
				sameSystem = false;
				break;
			}
		}
	}
	errorMsg.Show(!sameSystem);
	if (!sameSystem) 
		errorMsg.SetLabel(t_("Some bodies global or body axis mismatch"));
	
	UVector<int> idxs = ArrayModel_IndexsHydro(listLoaded);

	Progress progress(t_("Processing loaded data..."), 18);
	int pos = 0;
	mainTab.GetItem(mainTab.Find(mainMatrixA)).Enable(mainMatrixA.Load(Bem().hydros, idxs, ~menuPlot.showNdim));	progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainMatrixM)).Enable(mainMatrixM.Load(Bem().hydros, idxs, false));					progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainMatrixK)).Enable(mainMatrixK.Load(Bem().hydros, idxs, ~menuPlot.showNdim));	progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainMatrixK2)).Enable(mainMatrixK2.Load(Bem().hydros, idxs, ~menuPlot.showNdim));	progress.SetPos(pos++);	
	mainTab.GetItem(mainTab.Find(mainMatrixDlin)).Enable(mainMatrixDlin.Load(Bem().hydros, idxs, false));			progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainMatrixDquad)).Enable(mainMatrixDquad.Load(Bem().hydros, idxs, false));			progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainA)).Enable(mainA.Load(idxs));													progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainAinfw)).Enable(mainAinfw.Load(idxs));											progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainB)).Enable(mainB.Load(idxs));													progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainK)).Enable(mainK.Load(idxs));													progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainMD)).Enable(mainMD.Load(idxs, menuPlot.headMD.GetCursor()));					progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainForceSC)).Enable(mainForceSC.Load(idxs, menuPlot.head1st.GetCursor()));		progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainForceFK)).Enable(mainForceFK.Load(idxs, menuPlot.head1st.GetCursor()));		progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainForceEX)).Enable(mainForceEX.Load(idxs, menuPlot.head1st.GetCursor()));		progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainRAO)).Enable(mainRAO.Load(idxs, menuPlot.head1st.GetCursor()));				progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainQTF)).Enable(mainQTF.Load());													progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainSetupFOAMM)).Enable(/*data.IsLoadedB() && */idxs.size() > 0);					progress.SetPos(pos++);
	mainTab.GetItem(mainTab.Find(mainBody)).Enable(mainBody.Load());												progress.SetPos(pos++);
	
	bool isLoadedSS = false;
	for (int idx = 0; idx < Bem().hydros.size(); ++idx) {
		Hydro &hy = Bem().hydros[idx];
		if (hy.IsLoadedStateSpace()) {
			isLoadedSS = true;
			break;
		}
	}	
	mainTab.GetItem(mainTab.Find(mainStateSpace)).Enable(isLoadedSS);
	
	int id = mainTab.Get();
	if (id >= 0 && mainTab.GetItem(id).IsEnabled() == false)
		mainTab.Set(0);
	
	UpdateButtons();
}

void MainBEM::OnDescription() {
	try {
		int idx = GetIndexOneSelected();
		if (idx < 0) 
			return;
	
		WithDescription<TopWindow> w;
		CtrlLayout(w);
		w.Title(t_("Enter model name and description"));
		w.butOK 	  << [&] {w.Close();};
		w.name  	  <<= Bem().hydros[idx].dt.name;
		w.description <<= Bem().hydros[idx].dt.description;
		
		w.Execute();
		
		Bem().hydros[idx].dt.name = ~w.name;
		Bem().hydros[idx].dt.description = ~w.description;
		
		ArrayModel_Change(listLoaded, Bem().GetHydroIndex(idx), Null, ~w.name, Null);
		
		mainSummary.Clear();
		for (int i = 0; i < Bem().hydros.size(); ++i)
			mainSummary.Report(Bem().hydros[i], i);
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
}

void MainBEM::OnSolve() {
	try {
		int idx = GetIndexOneSelected();
		if (idx < 0) 
			return;
	
		Ma().GetMainSolver().loadFrom <<= Bem().hydros[idx].dt.file;
		Ma().tab.Set(Ma().GetMainSolver());
		Ma().GetMainSolver().butLoad.WhenAction();

	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

int MainBEM::AskQtfHeading(const Hydro &hy) {
	int numH = int(hy.dt.qhead.size());
	if (numH <= 1)
		return Null;
	
	WithSaveQTF<TopWindow> dialog;
	CtrlLayout(dialog);
	
	dialog.Title(t_("Please choose the QTF headings to save"));
	
	dialog.swHeadings <<= 0;
	dialog.dropHeadings.Enable(false);
	
	int id0 = Null;
	for (int i = 0; i < numH; ++i) {
		dialog.dropHeadings.Add(Format("%.3f-%.3f", hy.dt.qhead[i].real(), hy.dt.qhead[i].imag()));
		if (abs(hy.dt.qhead[i].real()) < 0.001)
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
		int idx = GetIndexOneSelected();
		if (idx < 0) 
			return;
		
		Status(t_("Saving BEM data"));
		String fileType = ~menuOpen.dropExport;
		Hydro::BEM_FMT type = Hydro::GetCodeBemStr(fileType);
		String ext = Replace(Hydro::bemExt[type], "*", "");
		
		FileSel fs;
		
		for (int i = 0; i < Hydro::NUMBEM; ++i)
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
			qtfHeading = AskQtfHeading(Bem().hydros[idx]);	
		
		Progress progress(t_("Saving BEM files..."), 100); 
		progress.Granularity(1000);
		
		WaitCursor wait;
		
		Bem().hydros[idx].SaveAs(fileName, [&](String str, int _pos) {
			if (!IsEmpty(str))
				progress.SetText(str); 
			if (_pos >= 0)
				progress.SetPos(_pos); 
			return !progress.Canceled();}, type, qtfHeading, mainBody.GetIb());
			
		saveFolder = GetFileFolder(~fs);
	} catch (Exc e) {
		BEM::PrintError(e);
	}
}

int MainBEM::GetIndexOneSelected(bool complain) {
	int id = -1;
	for (int row = 0; row < listLoaded.GetCount(); ++row) {
		if (listLoaded.IsSelected(row)) {
			id = ArrayModel_IndexHydro(listLoaded, row);
			break;
		}
	}	// Only one available => directly selected
	if (id < 0 && listLoaded.GetCount() == 1)
		id = ArrayModel_IndexHydro(listLoaded, 0);
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
/*
UVector<int> MainBEM::GetIdsSelected(bool complain) {
	UVector<int> ret = ArrayCtrlSelectedGet(listLoaded);

	if (ret.IsEmpty() && listLoaded.GetCount() == 1)
		ret << ArrayModel_IndexHydro(listLoaded, 0);
	if (ret.IsEmpty()) {
		if (complain)
			BEM::PrintError(t_("No model selected"));
		return ret;
	}
	return ret;
}	*/

void MainBEM::Jsonize(JsonIO &json) {
	if (json.IsLoading()) {
		menuPlot.autoFit = Null;
		menuPlot.fromY0 = Null;
		menuPlot.opwT = Null;
		menuPlot.opMP = Null;
		menuPlot.showPoints = Null;
		menuPlot.showNdim = Null;
		menuPlot.opAinf = Null;
		menuPlot.opA0 = Null;
		menuPlot.opApot = Null;
		menuPlot.opB = Null;
		menuPlot.opBhask = Null;
		menuPlot.opBpot = Null;
		menuPlot.opFfkpot = Null;
		menuPlot.opFscpot = Null;	
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
		("mainStiffness", mainMatrixK)
		("mainStiffness2", mainMatrixK2)
		("mainMatrixA", mainMatrixA)
		("mainMatrixDlin", mainMatrixDlin)
		("mainMatrixDquad", mainMatrixDquad)
		("mainMatrixM", mainMatrixM)
		("menuPlot_opAinf", menuPlot.opAinf)
		("menuPlot_opA0", menuPlot.opA0)
		("menuPlot_opApot", menuPlot.opApot)
		("menuPlot_opB", menuPlot.opB)
		("menuPlot_opBhask", menuPlot.opBhask)
		("menuPlot_opBpot", menuPlot.opBpot)
		("menuPlot_opFfkpot", menuPlot.opFfkpot)
		("menuPlot_opFscpot", menuPlot.opFscpot)
	;
	if (json.IsLoading()) {
		if (IsNull(dropExportId) || dropExportId < 0)
			dropExportId = 0;
	}
}

String MainBEM::BEMFile(String fileFolder) const {
	if (ToLower(fileFolder) == "id.dat")		// This a Nemoh file
		return "";
	if (DirectoryExists(fileFolder)) {
		int bestipos = INT_MAX;
		for (FindFile ff(AFX(fileFolder, "*.*")); ff; ff++) {
			if (ff.IsFile()) {
				String name = ToLower(ff.GetName());
				if (GetFileExt(name) == ".bemr")
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
				String controlFile = AFX(ff.GetPath(), "ControlFile.in");
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
		for (FindFile ff(AFX(folder, "*.*")); ff; ff++) {
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
		BEM::PrintError(t_("Unknown file types"));
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

void MainSummaryCoeff::Report(const Hydro &hy, int id) {
	if (array.GetColumnCount() == 0)
		array.AddColumn("Param");
	if (id >= array.GetColumnCount()-1)
		array.AddColumn(Format("#%d %s", id+1, hy.dt.name));
	int row = 0;
	int col = id + 1;
	
	Upp::Color lightRed = Upp::Color(255, 165, 158);								
	
	const MainBEM &bm = GetDefinedParent<MainBEM>(this);
	Upp::Color color = ArrayModel_GetColor(bm.listLoaded, id);
	::Color textColor = Black();
	if (Grayscale(color) < 150)
		textColor = White();
						
	array.Set(row, 0, t_("File"));				array.Set(row++, col, AttrText(hy.dt.file).Paper(color).Ink(textColor).Bold()); 
	array.Set(row, 0, t_("Name"));				array.Set(row++, col, hy.dt.name);
	array.Set(row, 0, t_("Description"));		array.Set(row++, col, hy.dt.description);
	array.Set(row, 0, t_("Software"));			array.Set(row++, col, hy.GetCodeStr());
	array.Set(row, 0, t_("g [m/s2]"));			array.Set(row++, col, hy.S_g());
	array.Set(row, 0, t_("rho [kg/m³]"));		array.Set(row++, col, hy.S_rho());
	array.Set(row, 0, t_("h (water depth) [m]"));array.Set(row++,col, hy.S_h());
	array.Set(row, 0, t_("length scale [m]"));	array.Set(row++, col, hy.S_len());
	
	array.Set(row, 0, t_("#frequencies"));		array.Set(row++, col, Nvl2(hy.dt.Nf, 0)); 
	if (!hy.dt.w.IsEmpty()) {
		array.Set(row, 0, t_("freq 0 [rad/s]"));	array.Set(row++, col, hy.dt.w[0]);
	} else {
		array.Set(row, 0, t_("freq 0 [rad/s]"));	array.Set(row++, col, "-");
	}
	if (hy.dt.w.size() > 1) {
		array.Set(row, 0, t_("freq end [rad/s]"));	array.Set(row++, col, hy.dt.w[hy.dt.w.size()-1]);
		double avg;
		if (hy.GetIrregularFreq(avg) < 0) { 
			array.Set(row, 0, t_("freq delta [rad/s]"));	array.Set(row++, col, avg);
		} else {
			String strHead;
			for (int i = 0; i < hy.dt.w.size(); ++i) {
				if (i > 0)
					strHead << ", ";
				strHead << hy.dt.w[i];
			}
			array.Set(row, 0, t_("freq delta [rad/s]"));
			array.Set(row++, col, Format(t_("Non constant delta (%s)"), strHead));
		}
	} else {
		array.Set(row, 0, t_("freq end [rad/s]"));	array.Set(row++, col, "-");
		array.Set(row, 0, t_("freq delta [rad/s]"));array.Set(row++, col, "-");
	}
	
	array.Set(row, 0, t_("#1st order headings"));			array.Set(row++, col, Nvl2(hy.dt.Nh, 0));
	if (!hy.dt.head.IsEmpty()) {
		array.Set(row, 0, t_("head 0 [º]"));	array.Set(row++, col, hy.dt.head[0]);
	} else {
		array.Set(row, 0, t_("head 0 [º]"));	array.Set(row++, col, "-");
	}
	if (hy.dt.head.size() > 1) {
		array.Set(row, 0, t_("head end [º]"));	array.Set(row++, col, hy.dt.head[hy.dt.head.size()-1]);
		double avg;
		if (hy.GetIrregularHead(avg) < 0) { 
			array.Set(row, 0, t_("head delta [º]"));	array.Set(row++, col, avg);
		} else {
			String strHead;
			for (int i = 0; i < hy.dt.head.size(); ++i) {
				if (i > 0)
					strHead << ", ";
				strHead << hy.dt.head[i];
			}
			array.Set(row, 0, t_("head delta [º]"));
			array.Set(row++, col, Format(t_("Non constant delta (%s)"), strHead));
		}
	} else {
		array.Set(row, 0, t_("head end [º]"));	array.Set(row++, col, "-");
		array.Set(row, 0, t_("head delta [º]"));array.Set(row++, col, "-");
	}
	
	array.Set(row, 0, t_("A0 available"));		array.Set(row++, col, hy.IsLoadedA0()   ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("A∞ available"));		array.Set(row++, col, hy.IsLoadedAinf() ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("A available"));		array.Set(row++, col, hy.IsLoadedA() 	  ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("B available"));		array.Set(row++, col, hy.IsLoadedB() 	  ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("Hydro stiff available"));		array.Set(row++, col, hy.IsLoadedC() 	  ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("Moor stiff available"));		array.Set(row++, col, hy.IsLoadedCMoor()? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("Inertia available"));	array.Set(row++, col, hy.IsLoadedM() 	  ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("Fex available"));		array.Set(row++, col, hy.IsLoadedFex()  ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("Fsc available"));		array.Set(row++, col, hy.IsLoadedFsc()  ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("Ffk available"));		array.Set(row++, col, hy.IsLoadedFfk()  ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("RAO available"));		array.Set(row++, col, hy.IsLoadedRAO()  ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("Linear damping available"));		array.Set(row++, col, hy.IsLoadedDlin()   ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("Quadratic damping available"));	array.Set(row++, col, hy.IsLoadedDquad()  ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("Mean Drift available"));		array.Set(row++, col, hy.IsLoadedMD() 	? t_("Yes") : t_("No"));
	
	array.Set(row, 0, t_("#bodies"));			array.Set(row++, col, hy.dt.Nb);
	for (int ib = 0; ib < hy.dt.Nb; ++ib) {
		String sib = Format("#%d", ib+1);
		//sib += " " + data.msh[ib].name;
		array.Set(row, 0, sib + " " + t_("Name"));		array.Set(row++, col, hy.dt.msh[ib].dt.name);
		/*array.Set(row, 0, sib + " " + t_("#dof"));
		if (data.dof.size() > ib) 
			array.Set(row++, col, data.dof[ib]);
		else 
			array.Set(row++, col, "-");*/
		
		array.Set(row, 0, sib + " " + t_("Vsub [m³]"));
		if (/*data.Vo.size() > ib && */IsNum(hy.dt.msh[ib].dt.Vo)) {
			
			array.Set(row++, col, AttrText(FDS(hy.dt.msh[ib].dt.Vo, 10, false)).Paper(hy.dt.msh[ib].dt.Vo < 0 ? lightRed : Null));
		} else 
			array.Set(row++, col, "-");
		
		array.Set(row, 0, sib + " " + t_("Cg [m]"));
		if (IsNum(hy.dt.msh[ib].dt.cg)) 
			array.Set(row++, col, Format(t_("%s, %s, %s"),
									FDS(hy.dt.msh[ib].dt.cg.x, 10, false),
									FDS(hy.dt.msh[ib].dt.cg.y, 10, false),
									FDS(hy.dt.msh[ib].dt.cg.z, 10, false)));
		else
			array.Set(row++, col, "-");

		array.Set(row, 0, sib + " " + t_("Cb [m]"));
		if (IsNum(hy.dt.msh[ib].dt.cb)) 
			array.Set(row++, col, Format(t_("%s, %s, %s"),
									FDS(hy.dt.msh[ib].dt.cb.x, 10, false),
									FDS(hy.dt.msh[ib].dt.cb.y, 10, false),
									FDS(hy.dt.msh[ib].dt.cb.z, 10, false)));
		else
			array.Set(row++, col, "-");
		
		array.Set(row, 0, sib + " " + t_("C0 [m]"));
		if (IsNum(hy.dt.msh[ib].dt.c0)) 
			array.Set(row++, col, Format(t_("%s, %s, %s"),
									FDS(hy.dt.msh[ib].dt.c0.x, 10, false),
									FDS(hy.dt.msh[ib].dt.c0.y, 10, false),
									FDS(hy.dt.msh[ib].dt.c0.z, 10, false)));
		else
			array.Set(row++, col, "-");
		
		array.Set(row, 0, sib + " " + t_("Waterplane area [m²]"));
		if (/*data.C.size() > ib && */hy.dt.msh[ib].dt.C.size() > 0) {
			double wPlaneArea = hy.C_dim(ib, 2, 2);
			array.Set(row++, col, FDS(wPlaneArea, 10, false));		
			for (int i = 0; i < 6; ++i) {
				for (int j = 0; j < 6; ++j) {
					if (!Hydro::C_units(i, j).IsEmpty()) {
						array.Set(row, 0, sib + " " + Format(t_("K(%d,%d) [%s]"), i+1, j+1, Hydro::C_units(i, j)));	array.Set(row++, col, FDS(hy.C_dim(ib, i, j), 10, false));		
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
		//array.Set(row,   0, sib + " " + t_("Theave(∞) [s]"));
		array.Set(row, 0, sib + " " + t_("Theave(ω) [s]"));
		//array.Set(row+2, 0, sib + " " + t_("Troll(∞)  [s]"));
		array.Set(row+1, 0, sib + " " + t_("Troll(ω)  [s]"));
		//array.Set(row+4, 0, sib + " " + t_("Tpitch(∞) [s]"));
		array.Set(row+2, 0, sib + " " + t_("Tpitch(ω) [s]"));
		if (/*IsNum(data.rho) && IsNum(data.g) &&*/ 
			/*data.M.size() > ib && */hy.dt.msh[ib].dt.M.size() > 0 && 
			/*data.C.size() > ib && */hy.dt.msh[ib].dt.C.size() > 0) {
			//array.Set(row++, col, FDS(data.Theave (ib), 5, false, "-"));
			array.Set(row++, col, FDS(hy.Tdof(ib, 2), 5, false, FDS(hy.Tdof_inf(ib, 2), 5, false, "-") + S(" (∞)")));
			//array.Set(row++, col, FDS(data.Troll  (ib), 5, false, "-"));
			array.Set(row++, col, FDS(hy.Tdof(ib, 3), 5, false, FDS(hy.Tdof_inf(ib, 3), 5, false, "-") + S(" (∞)")));
			//array.Set(row++, col, FDS(data.Tpitch (ib), 5, false, "-"));
			array.Set(row++, col, FDS(hy.Tdof(ib, 4), 5, false, FDS(hy.Tdof_inf(ib, 4), 5, false, "-") + S(" (∞)")));
		} else {
			array.Set(row++, col, "-");	
			array.Set(row++, col, "-");	
			array.Set(row++, col, "-");	
			//array.Set(row++, col, "-");	
			//array.Set(row++, col, "-");	
			//array.Set(row++, col, "-");	
		}
		array.Set(row,   0, sib + " " + t_("GMroll  [m]"));
		array.Set(row+1, 0, sib + " " + t_("GMpitch [m]"));
		if (/*IsNum(data.rho) && IsNum(data.g) && */hy.IsLoadedC()) {
			array.Set(row++, col, FDS(hy.GMroll(ib), 5, false, "-"));
			array.Set(row++, col, FDS(hy.GMpitch(ib), 5, false, "-"));
		} else {
			array.Set(row++, col, "-");	
			array.Set(row++, col, "-");	
		}
	}	
}

void MainOutput::Init() {
	CtrlLayout(*this);
	cout.SetReadOnly();
	Print(t_("BEMRosetta\nHydrodynamic solvers viewer and converter\n"));
}

void MainOutput::Print(String str) {
	cout.Append(str);
	cout.ScrollEnd();
}
