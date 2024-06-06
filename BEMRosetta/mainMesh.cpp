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


void MainBody::Init() {
	MainBEMBody::Init();
	
	OnOpt();
	
	CtrlLayout(menuOpen);
	menuOpen.file.WhenChange = THISBACK(OnLoad);
	menuOpen.file.BrowseRightWidth(40).UseOpenFolder(true).BrowseOpenFolderWidth(10);
	menuOpen.butLoad.WhenAction = [&] {menuOpen.file.DoGo();};

	ArrayModel_Init(listLoaded, true).MultiSelect();
	listLoaded.WhenSel = [&] {
		OnMenuOpenArraySel();
		OnMenuProcessArraySel();
		OnMenuMoveArraySel();
		OnMenuAdvancedArraySel();
		LoadSelTab(Bem());
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
	menuOpen.butRemove.WhenAction = THISBACK(OnRemove);
	menuOpen.butRemoveSelected.Disable();	
	menuOpen.butRemoveSelected.WhenAction = THISBACK1(OnRemoveSelected, false);
	menuOpen.butJoin.Disable();	
	menuOpen.butJoin.WhenAction = THISBACK(OnJoin);
	menuOpen.butSplit.Disable();	
	menuOpen.butSplit.WhenAction = THISBACK(OnSplit);
	menuOpen.butExport <<= THISBACK(OnConvertBody);
	for (int i = 0; i < Body::NUMMESH; ++i)
		if (Body::meshCanSave[i])
			menuOpen.dropExport.Add(Body::GetBodyStr(static_cast<Body::MESH_FMT>(i)));
	menuOpen.dropExport.SetIndex(dropExportId);
	menuOpen.dropExport <<= THISBACK(OnOpt);
	
	OnOpt();
	
	CtrlLayout(menuPlot);
	menuPlot.showBody.Tip(t_("Shows all the loaded mesh")).WhenAction 	    	= [&] {LoadSelTab(Bem());};
	menuPlot.showNormals.Tip(t_("Shows panel normals as arrows")).WhenAction 	= [&] {LoadSelTab(Bem());};
	menuPlot.showWaterLevel.Tip(t_("Shows the cut of the mesh with the water line")).WhenAction  = [&] {LoadSelTab(Bem());};
	menuPlot.showSkewed.Tip(t_("Shows skewed panels")).WhenAction      			= [&] {LoadSelTab(Bem());};
	menuPlot.showFissure.Tip(t_("Shows fissures in the hull")).WhenAction     	= [&] {LoadSelTab(Bem());};
	menuPlot.showMultiPan.Tip(t_("Shows wrong panels")).WhenAction    			= [&] {LoadSelTab(Bem());};
	menuPlot.showAxis.Tip(t_("Shows system axis")).WhenAction  					= [&] {mainView.gl.Refresh();};
	menuPlot.showLimits.Tip(t_("Shows boundaris of the geometry")).WhenAction 	= [&] {mainView.gl.Refresh();};
	menuPlot.showCb.Tip(t_("Shows the centre of buoyancy")).WhenAction  		= [&] {mainView.gl.Refresh();};
	menuPlot.showCg.Tip(t_("Shows the centre of gravity")).WhenAction  			= [&] {mainView.gl.Refresh();};
	menuPlot.showCr.Tip(t_("Shows the centre of motion")).WhenAction  		    = [&] {mainView.gl.Refresh();};
	menuPlot.showSel.Tip(t_("Shows volume around selected object")).WhenAction  = [&] {mainView.gl.Refresh();};
	menuPlot.showCr.Tip(t_("Shows the lines")).WhenAction  						= [&] {mainView.gl.Refresh();};	
	menuPlot.showUnderwater.Tip(t_("Shows nderwater mesh")).WhenAction  		= [&] {mainView.gl.Refresh();};
	menuPlot.butXYZ.Tip(t_("Orients the camera as isometric")).WhenAction  		= [&] {mainView.gl.View(true, true, true);};
	menuPlot.butXoY.Tip(t_("Orients the camera through Z axis")).WhenAction  	= [&] {mainView.gl.View(true, true, false);};	
	menuPlot.butYoZ.Tip(t_("Orients the camera through X axis")).WhenAction  	= [&] {mainView.gl.View(false, true, true);};
	menuPlot.butXoZ.Tip(t_("Orients the camera through Y axis")).WhenAction  	= [&] {mainView.gl.View(true, false, true);};
	menuPlot.butFit.Tip(t_("Zooms the camera to fit the bodies")).WhenAction	= [&] {mainView.gl.ZoomToFit();};
	menuPlot.showBodyData.Tip(t_("Shows a list of panels and nodes")).WhenAction= [&] {splitterAll.SetButton(0);};
	menuPlot.showBodyData.Tip(t_("Controls for 3D playing")).WhenAction			= [&] {splitterVideo.SetButton(0);};
	menuPlot.backColor.Tip(t_("Sets the background color")).WhenAction      	= [&] {LoadSelTab(Bem());};
	menuPlot.lineThickness.Tip(t_("Sets the thickness of the mesh wireframe")).WhenAction = [&] {LoadSelTab(Bem());};
	
	styleRed = styleGreen = styleBlue = Button::StyleNormal();
	styleRed.textcolor[0] = styleRed.textcolor[1] = styleRed.textcolor[2] = LtRed();
	styleGreen.textcolor[0] = styleGreen.textcolor[1] = styleGreen.textcolor[2] = Green();
	styleBlue.textcolor[0] = styleBlue.textcolor[1] = styleBlue.textcolor[2] = LtBlue();
	menuPlot.butYoZ.SetStyle(styleRed);
	menuPlot.butXoZ.SetStyle(styleGreen);
	menuPlot.butXoY.SetStyle(styleBlue);
	
	OnOpt();
	
	CtrlLayout(menuProcess);
	
	menuProcess.x_g <<= Null;
	menuProcess.x_g.WhenEnter = THISBACK2(OnUpdate, NONE, true);
	menuProcess.y_g <<= Null;
	menuProcess.y_g.WhenEnter = THISBACK2(OnUpdate, NONE, true);
	menuProcess.z_g <<= Null;
	menuProcess.z_g.WhenEnter = THISBACK2(OnUpdate, NONE, true);
	menuProcess.x_0 <<= 0;
	menuProcess.x_0.WhenEnter = THISBACK2(OnUpdate, NONE, true);
	menuProcess.y_0 <<= 0;
	menuProcess.y_0.WhenEnter = THISBACK2(OnUpdate, NONE, true);
	menuProcess.z_0 <<= 0;
	menuProcess.z_0.WhenEnter = THISBACK2(OnUpdate, NONE, true);
	menuProcess.mass <<= 0;
	menuProcess.mass.WhenEnter = THISBACK2(OnUpdate, NONE, true);
	menuProcess.butUpdateCg  << THISBACK2(OnUpdate, NONE, true);
	menuProcess.butUpdateCg.Tip(t_("Sets the centre of gravity"));
	menuProcess.x_g.WhenEnter << THISBACK2(OnUpdate, NONE, true);
	menuProcess.y_g.WhenEnter << THISBACK2(OnUpdate, NONE, true);
	menuProcess.z_g.WhenEnter << THISBACK2(OnUpdate, NONE, true);
	
	menuProcess.butCgtoC0.WhenAction = [&] {
		menuProcess.x_0 <<= ~menuProcess.x_g;
		menuProcess.y_0 <<= ~menuProcess.y_g;
		menuProcess.z_0 <<= ~menuProcess.z_g;
		OnUpdate(NONE, true);
	};
	menuProcess.butCgtoC0.Tip(t_("Sets the centre of motion with the centre of gravity"));
	menuProcess.butC0toCg.WhenAction = [&] {
		menuProcess.x_g <<= ~menuProcess.x_0;
		menuProcess.y_g <<= ~menuProcess.y_0;
		menuProcess.z_g <<= ~menuProcess.z_0;
		OnUpdate(NONE, true);
	};
	menuProcess.butC0toCg.Tip(t_("Sets the centre of gravity with the centre of motion"));
	menuProcess.butUpdateCrot  <<= THISBACK2(OnUpdate, NONE, true);
	menuProcess.butUpdateCrot.Tip(t_("Sets the centre of the body axis"));	
	menuProcess.x_0.WhenEnter << THISBACK2(OnUpdate, NONE, true);
	menuProcess.y_0.WhenEnter << THISBACK2(OnUpdate, NONE, true);
	menuProcess.z_0.WhenEnter << THISBACK2(OnUpdate, NONE, true);	
	
	menuProcess.butUpdateMassVol  <<= THISBACK(OnUpdateMass);
	menuProcess.butUpdateMassVol.Tip(t_("Sets mass from inmersed volume"));
	menuProcess.butUpdateMass  <<= THISBACK2(OnUpdate, NONE, true);
	menuProcess.butUpdateMass.Tip(t_("Sets mass"));
	
	menuProcess.butImageX <<= THISBACK1(OnImage, 0);
	menuProcess.butImageX.Tip(t_("Mirrors the mesh in X axis"));
	menuProcess.butImageY <<= THISBACK1(OnImage, 1);
	menuProcess.butImageY.Tip(t_("Mirrors the mesh in Y axis"));
	menuProcess.butImageZ <<= THISBACK1(OnImage, 2);
	menuProcess.butImageZ.Tip(t_("Mirrors the mesh in Z axis"));
	
	menuProcess.butSimplify <<= THISBACK1(OnHealing, true);
	menuProcess.butSimplify.Tip(t_("Simplify mesh removing duplicated elements"));
	menuProcess.butFullHealing <<= THISBACK1(OnHealing, false);
	menuProcess.butFullHealing.Tip(t_("Tries to fix problems in mesh"));
	menuProcess.butOrientSurface <<= THISBACK(OnOrientSurface);
	menuProcess.butOrientSurface.Tip(t_("Orient all face normals to one side. Set show normals to see the results"));
	
	menuProcess.butWaterFill  <<= THISBACK1(OnAddWaterSurface, 'f');
	menuProcess.butWaterFill.Tip(t_("Generates waterplane mesh based on how the hull crosses the waterplane"));
	menuProcess.butWaterPlane <<= THISBACK1(OnAddWaterSurface, 'e');
	menuProcess.butWaterPlane.Tip(t_("Extracts waterplane from a mesh that already includes it"));
	menuProcess.butHull		  <<= THISBACK1(OnAddWaterSurface, 'r');
	menuProcess.butHull.Tip(t_("Extracts underwater (wet) hull from a mesh"));
	
	menuProcess.rat_x <<= 1;
	menuProcess.rat_y <<= 1;
	menuProcess.rat_z <<= 1;
	menuProcess.butScale <<= THISBACK(OnScale);
	menuProcess.butScale.Tip(t_("Scales the mesh"));
	
	menuProcess.butInertia <<= THISBACK(OnInertia);
	
	CtrlLayout(menuMove);	
	
	menuMove.butReset <<= THISBACK(OnReset);
	menuMove.butReset.Tip(t_("Translates the mesh"));
	
	menuMove.butUpdateCrot  <<= THISBACK2(OnUpdate, NONE, false);
	menuMove.butUpdateCrot.Tip(t_("Sets the centre of the body axis"));	
	
	menuMove.t_x <<= 0;
	menuMove.t_x.WhenEnter = THISBACK2(OnUpdate, ROTATE, false);
	menuMove.t_y <<= 0;
	menuMove.t_y.WhenEnter = THISBACK2(OnUpdate, ROTATE, false);
	menuMove.t_z <<= 0;
	menuMove.t_z.WhenEnter = THISBACK2(OnUpdate, ROTATE, false);
	menuMove.a_x <<= 0;
	menuMove.a_x.WhenEnter = THISBACK2(OnUpdate, ROTATE, false);
	menuMove.a_y <<= 0;
	menuMove.a_y.WhenEnter = THISBACK2(OnUpdate, ROTATE, false);
	menuMove.a_z <<= 0;
	menuMove.a_z.WhenEnter = THISBACK2(OnUpdate, ROTATE, false);
	menuMove.butUpdatePos <<= THISBACK2(OnUpdate, MOVE, false);
	menuMove.butUpdatePos.Tip(t_("Translates the mesh"));
	menuMove.butUpdateAng <<= THISBACK2(OnUpdate, ROTATE, false);
	menuMove.butUpdateAng.Tip(t_("Rotates the mesh"));	
	
	menuMove.butArchimede <<= THISBACK(OnArchimede);
	menuMove.butArchimede.Tip(t_("Lets the body fall to rest"));	
	
	menuMove.opZArchimede.WhenAction = [&] {menuMove.t_z.Enable(!menuMove.opZArchimede);};
	
	CtrlLayout(menuEdit);
	menuEdit.edit_x <<= 0;
	menuEdit.edit_y <<= 0;
	menuEdit.edit_z <<= 0;
	menuEdit.edit_size <<= 1;
	
	menuEdit.butPanel <<= THISBACK(OnAddPanel);
	menuEdit.panWidthX.WhenEnter << THISBACK(OnAddPanel);
	menuEdit.panWidthY.WhenEnter << THISBACK(OnAddPanel);
	
	menuEdit.revolutionList.AddColumn("H").Ctrls<EditString>().HeaderTab().SetMargin(-2);
	menuEdit.revolutionList.AddColumn("V").Ctrls<EditString>().HeaderTab().SetMargin(-2);
	menuEdit.revolutionList.AddColumn("");
	menuEdit.revolutionList.ColumnWidths("10 10 2");
	menuEdit.revolutionList.NoHeader().MultiSelect();
	menuEdit.revolutionList.SetLineCy(int(EditField::GetStdHeight()*2/3));
	menuEdit.revolutionList.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, menuEdit.revolutionList, true, true);};
	
	menuEdit.butRevolution << THISBACK(OnAddRevolution);
	menuEdit.revolutionList.WhenLeftDouble << THISBACK(OnAddRevolution);
	
	menuEdit.polynomialList.AddColumn("H").Ctrls<EditString>().HeaderTab().SetMargin(-2);
	menuEdit.polynomialList.AddColumn("V").Ctrls<EditString>().HeaderTab().SetMargin(-2);
	menuEdit.polynomialList.AddColumn("");
	menuEdit.polynomialList.ColumnWidths("10 10 2");
	menuEdit.polynomialList.NoHeader().MultiSelect();
	menuEdit.polynomialList.SetLineCy(int(EditField::GetStdHeight()*2/3));
	menuEdit.polynomialList.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, menuEdit.polynomialList, true, true);};
	
	menuEdit.butPolynomial << THISBACK(OnAddPolygonalPanel);
	menuEdit.polynomialList.WhenLeftDouble << THISBACK(OnAddPolygonalPanel);
		
	menuTab.Add(menuOpen.SizePos(),    	t_("Load"));
	menuTab.Add(menuPlot.SizePos(),    	t_("Plot")).Disable();
	menuTab.Add(menuMove.SizePos(), 	t_("Move")).Disable();
	menuTab.Add(menuProcess.SizePos(), 	t_("Process")).Disable();	
	menuTab.Add(menuEdit.SizePos(), 	t_("Edit"));
		
	mainViewData.Init();		// Un-comment this to view video controls
	//splitterVideo.Vert(mainView.SizePos(), videoCtrl.SizePos());
	//splitterVideo.SetPositions(8500, 9900).SetInitialPositionId(1).SetButtonNumber(1).SetButtonWidth(20);
	//splitterAll.Horz(splitterVideo.SizePos(), mainViewData.SizePos());
	splitterAll.Horz(mainView.SizePos(), mainViewData.SizePos());
	splitterAll.SetPositions(6000, 9900).SetInitialPositionId(1).SetButtonNumber(1).SetButtonWidth(20);
	splitterAll.Tip(t_("")).WhenAction = [&] {mainView.SetPaintSelect(splitterAll.GetPos() < 9900);};
	mainTab.Add(splitterAll.SizePos(), t_("View"));
	mainView.Init(*this);
	
	String bitmapFolder = AFX(GetDesktopFolder(), "BEMRosetta Body Images");
	int idBitmapFolder = 0;
	
	videoCtrl.Init([&](UVector<int> &ids)->int {
			ids = ArrayModel_IdsBody(listLoaded);
			int num = ArrayCtrlSelectedGetCount(listLoaded);
			if (num > 1) {
				BEM::PrintError(t_("Please select just one model"));
				return -1;
			}
			int id;
			if (num == 0 && listLoaded.GetCount() == 1)
				id = ArrayModel_IdBody(listLoaded, 0);
			else {
			 	id = ArrayModel_IdBody(listLoaded);
				if (id < 0) {
					BEM::PrintError(t_("Please select a model to process"));
					return -1;
				}
			}
			return id;
		}, [&](int id, const UVector<int> &ids, const Point3D &pos, const Point3D &angle, const Point3D &c0, bool full, bool saveBitmap) {
			Body &msh = Bem().surfs[id];			
			
			msh.dt.cg.TransRot(pos.x, pos.y, pos.z, ToRad(angle.x), ToRad(angle.y), ToRad(angle.z), c0.x, c0.y, c0.z);
			msh.dt.mesh.TransRot(pos.x, pos.y, pos.z, ToRad(angle.x), ToRad(angle.y), ToRad(angle.z), c0.x, c0.y, c0.z);
			
			menuProcess.x_g <<= msh.dt.cg.x;
			menuProcess.y_g <<= msh.dt.cg.y;
			menuProcess.z_g <<= msh.dt.cg.z;
			
			msh.AfterLoad(Bem().rho, Bem().g, false, false);
			
			if (full)
				mainStiffness.Load(Bem().surfs, ids);
			mainView.CalcEnvelope();
			if (full)
				mainSummary.Report(Bem().surfs, id);
			
			mainView.gl.Refresh();
			if (saveBitmap) {
				RealizeDirectory(bitmapFolder);
				DeleteFileDeepWildcardsX(bitmapFolder);
				mainView.gl.SaveToFile(AFX(bitmapFolder, Format("Image%4d", idBitmapFolder++)));
			}
			if (full)
				mainViewData.OnRefresh();
		});
	
	
	mainSummary.Init();
	mainTab.Add(mainSummary.SizePos(), t_("Summary"));

	mainStiffness.Init(Hydro::MAT_K);
	mainTab.Add(mainStiffness.SizePos(), t_("Hydrostatic Stiffness"));
	
	mainStiffness.opMassBuoy.WhenAction = THISBACK2(OnUpdate, NONE, true);
	
	mainStiffness2.Init(Hydro::MAT_K2);
	mainTab.Add(mainStiffness2.SizePos(), t_("Mooring Stiffness"));
	
	mainGZ.Init();
	mainTab.Add(mainGZ.SizePos(), t_("Stability GZ"));
			
	mainTab.WhenSet = [&] {
		LOGTAB(mainTab);
		UVector<int> ids = ArrayModel_IdsBody(listLoaded);
		bool plot = true, move = false, convertProcess = true;
		if (Bem().surfs.IsEmpty()) 
			plot = convertProcess = false;
		else if (mainTab.IsAt(splitterAll)) 
			;
		else if (mainTab.IsAt(mainStiffness)) {
			plot = false;
			move = true;
			mainStiffness.Load(Bem().surfs, ids);
		} else 
			plot = false;
		
		TabCtrl::Item& tabMenuPlot = menuTab.GetItem(menuTab.Find(menuPlot));
		tabMenuPlot.Enable(plot);
		TabCtrl::Item& tabMenuMove = menuTab.GetItem(menuTab.Find(menuMove));
		tabMenuMove.Enable(convertProcess);
		TabCtrl::Item& tabMenuProcess = menuTab.GetItem(menuTab.Find(menuProcess));
		tabMenuProcess.Enable(convertProcess);
		if (plot) {
			tabMenuPlot.Text(t_("Plot"));
			menuTab.Set(menuPlot);
		} else {
			tabMenuPlot.Text("");
			if (move)
				menuTab.Set(menuMove);
			else
				menuTab.Set(menuOpen);
		}
		if (convertProcess) {
			tabMenuProcess.Text(t_("Process"));
			tabMenuMove.Text(t_("Move"));
		} else {
			tabMenuProcess.Text("");
			tabMenuMove.Text("");
		}
	};
	mainTab.WhenSet();
	
	menuTab.WhenSet = [&] {
	LOGTAB(menuTab);
	
	};
	menuTab.WhenSet();
	
	UpdateButtons();
	saveFolder = GetDesktopFolder();
}

void MainBody::OnMenuOpenArraySel() {
	int id = ArrayModel_IdBody(listLoaded);
	if (id < 0)
		return;
	
	Body::MESH_FMT type = Body::GetCodeBodyStr(~menuOpen.dropExport);
	menuOpen.symX <<= ((type == Body::WAMIT_GDF || type == Body::AQWA_DAT) && Bem().surfs[id].IsSymmetricX());
	menuOpen.symY <<= Bem().surfs[id].IsSymmetricY();
}

void MainBody::OnMenuProcessArraySel() {
	int id = ArrayModel_IdBody(listLoaded);
	if (id < 0)
		return;
	
	Body &msh = Bem().surfs[id];
	if (!IsNull(msh.dt.cg)) {
		menuProcess.x_g <<= msh.dt.cg.x;
		menuProcess.y_g <<= msh.dt.cg.y;
		menuProcess.z_g <<= msh.dt.cg.z;
	} else
		menuProcess.x_g <<= menuProcess.y_g <<= menuProcess.z_g <<= Null;
	
	menuProcess.x_0 <<= msh.dt.c0.x;
	menuProcess.y_0 <<= msh.dt.c0.y;
	menuProcess.z_0 <<= msh.dt.c0.z;
	menuProcess.mass <<= msh.GetMass();
}

void MainBody::OnMenuMoveArraySel() {
	int id = ArrayModel_IdBody(listLoaded);
	if (id < 0)
		return;
	
	Body &msh = Bem().surfs[id];
	menuMove.x_0 <<= msh.dt.c0.x;
	menuMove.y_0 <<= msh.dt.c0.y;
	menuMove.z_0 <<= msh.dt.c0.z;	
}

void MainBody::OnMenuAdvancedArraySel() {
	int id = ArrayModel_IdBody(listLoaded);
	if (id < 0)
		return;
}

void MainBody::OnArraySel() {
	OnMenuOpenArraySel();
	OnMenuProcessArraySel();
	OnMenuMoveArraySel();
	OnMenuAdvancedArraySel();
}
	
void MainBody::InitSerialize(bool ret) {
	if (!ret || IsNull(menuPlot.showBody)) 
		menuPlot.showBody = true;	
	if (!ret || IsNull(menuPlot.showNormals)) 
		menuPlot.showNormals = true;	
	if (!ret || IsNull(menuPlot.showWaterLevel)) 
		menuPlot.showWaterLevel = true;	
	if (!ret || IsNull(menuPlot.showSkewed)) 
		menuPlot.showSkewed = false;
	if (!ret || IsNull(menuPlot.showFissure)) 
		menuPlot.showFissure = false;
	if (!ret || IsNull(menuPlot.showMultiPan)) 
		menuPlot.showMultiPan = false;
	if (!ret || IsNull(menuPlot.showAxis)) 
		menuPlot.showAxis = true;
	if (!ret || IsNull(menuPlot.showLimits)) 
		menuPlot.showLimits = false;
	if (!ret || IsNull(menuPlot.showCg)) 
		menuPlot.showCg = true;
	if (!ret || IsNull(menuPlot.showCb)) 
		menuPlot.showCb = true;
	if (!ret || IsNull(menuPlot.showCr)) 
		menuPlot.showCr = true;
	if (!ret || IsNull(menuPlot.showLines)) 
		menuPlot.showLines = true;
	if (!ret || IsNull(menuPlot.showSel)) 
		menuPlot.showSel = true;	
	if (!ret || IsNull(menuPlot.lineThickness)) 
		menuPlot.lineThickness <<= 1;
	if (!ret || IsNull(menuPlot.backColor)) 
		menuPlot.backColor <<= White();
				
	if (!ret || IsNull(menuOpen.optBodyType)) 
		menuOpen.optBodyType = 0;
	if (!ret || IsNull(menuOpen.opClean)) 
		menuOpen.opClean = false;
}

void MainBody::LoadSelTab(BEM &bem) {
	const UVector<int> &ids = ArrayModel_IdsBody(listLoaded);
	if (mainTab.Get() == mainTab.Find(mainStiffness))
		mainStiffness.Load(bem.surfs, ids);
	else if (mainTab.Get() == mainTab.Find(mainView))
		mainView.gl.Refresh();
}

void MainBody::OnOpt() {
	menuOpen.file.ClearTypes(); 

	String meshFiles = Body::GetMeshExt();	//".gdf .dat .stl .pnl .msh .mesh .hst .grd .nc";
	
	String meshFilesAst = clone(meshFiles);
	meshFiles.Replace("*.", ".");
	menuOpen.file.Type(Format("All supported mesh files (%s)", meshFiles), meshFilesAst);
	menuOpen.file.AllFilesType();
	String extView = ToLower(GetFileExt(menuOpen.file.GetData().ToString()));
	if (extView.IsEmpty())
		menuOpen.file.ActiveType(0);
	else if (meshFiles.Find(extView) >= 0)
		menuOpen.file.ActiveType(0);
	else
		menuOpen.file.ActiveType(1);
	
	menuOpen.symX.Disable();
	menuOpen.symY.Disable();
	Body::MESH_FMT type = Body::GetCodeBodyStr(~menuOpen.dropExport);
		
	switch (type) {
	case Body::WAMIT_GDF:	menuOpen.symX.Enable();
							menuOpen.symY.Enable();
							break;
	case Body::HAMS_PNL:	menuOpen.symX.Enable();
							menuOpen.symY.Enable();
							break;
	case Body::AQWA_DAT:	menuOpen.symX.Enable();
							menuOpen.symY.Enable();
							break;
	case Body::NEMOH_DAT:	menuOpen.symY.Enable();
							break;	
	default:				break;		
	}
}

void MainBody::AfterAdd(String file, int num) {
	mainTab.Set(mainView);
	for (int id = Bem().surfs.size() - num; id < Bem().surfs.size(); ++id) {
		Body &surf = Bem().surfs[id];
		
		surf.Report(Bem().rho);
		
		AddRow(surf);
	}

	mainView.CalcEnvelope();
	
	mainView.gl.Enable();
	mainView.gl.ZoomToFit();
	
	After();	
	OnArraySel();	
}

bool MainBody::OnLoad() {
	GuiLock __;
	
	String file = ~menuOpen.file;
		
	try {
		Progress progress(t_("Loading mesh file..."), 100); 
		
		UVector<int> ids = ArrayModel_IdsBody(listLoaded);
		for (int i = 0; i < ids.size(); ++i) {
			if (Bem().surfs[ids[i]].dt.fileName == file) {
				if (!PromptYesNo(t_("Model is already loaded") + S("&") + t_("Do you wish to open it anyway?")))
					return false;
				break;
			}
		}
		mainView.gl.Disable();
		
		WaitCursor waitcursor;

		int num = Bem().LoadBody(file, [&](String str, int _pos) {
			progress.SetText(str); 
			progress.SetPos(_pos); 
			return !progress.Canceled();
		}, ~menuOpen.opClean, false);
		
		//int id = Bem().surfs.size()-1;
		for (int id = Bem().surfs.size() - num; id < Bem().surfs.size(); ++id) {
			Body &msh = Bem().surfs[id];
			
			if (!msh.dt.mesh.IsEmpty() && ~menuMove.opZArchimede) {
				double dz = 0.1;
				Surface under;
				if (msh.dt.mesh.TranslateArchimede(msh.GetMass(), Bem().rho, dz, under)) {
					msh.dt.cg.Translate(0, 0, dz);
					videoCtrl.AddReg(Point3D(0, 0, dz));
				} else
					BEM::PrintError(t_("Problem readjusting the Z value to comply with displacement"));
				
				Ma().Status(Format(t_("Loaded '%s', and translated vertically %f m to comply with displacement"), file, dz));
			} else
				Ma().Status(Format(t_("Loaded '%s'"), file));			
		}
		
		AfterAdd(file, num);
		
		mainViewData.OnAddedModel(mainView);
		OnOpt();
	} catch (Exc e) {
		mainView.gl.Enable();
		BEM::PrintError(e);
		return false;
	}
	return true;
}

void MainBody::OnConvertBody() {
	GuiLock __;
	
	try {
		String fileType = ~menuOpen.dropExport;
		Body::MESH_FMT type = Body::GetCodeBodyStr(fileType);
		String ext = Replace(Body::meshExt[type], "*", "");
		
		UVector<int> sel = ArrayCtrlSelectedGet(listLoaded);
		
		if (sel.size() > 1 && type != Body::AQWA_DAT) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		if (sel.size() == 0) {
			if (listLoaded.GetCount() == 1)
				sel << ArrayModel_IdBody(listLoaded, 0);
			else {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		}
		
		Status(t_("Saving mesh data"));
		
		FileSel fs;
		
		for (int i = 0; i < Body::NUMMESH; ++i)
			if (Body::meshCanSave[i] && (i == type || i == Body::UNKNOWN)) 
				fs.Type(Body::GetBodyStr(static_cast<Body::MESH_FMT>(i)), Body::meshExt[i]);
		
		fs.ActiveType(0);
		fs.Set(ForceExtSafer(~menuOpen.file, ext));
		fs.ActiveDir(saveFolder);
		
		if (!fs.ExecuteSaveAs(Format(t_("Save mesh file as %s"), fileType)))
			return;
		
		String fileName = ~fs;
		
		Progress progress(t_("Saving mesh file..."), 100); 
		progress.Granularity(1000);
		
		WaitCursor waitcursor;
		
		Bem().SaveBody(fileName, sel, type, 
							   static_cast<Body::MESH_TYPE>(int(~menuOpen.optBodyType)),
							   ~menuOpen.symX, ~menuOpen.symY);	
							   
		saveFolder = GetFileFolder(~fs);
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
}

void MainBody::OnReset() {
	GuiLock __;
	
	try {
		UVector<int> ids = ArrayModel_IdsBody(listLoaded);
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		int id;
		if (num == 0 && listLoaded.GetCount() == 1)
			id = ArrayModel_IdBody(listLoaded, 0);
		else {
		 	id = ArrayModel_IdBody(listLoaded);
			if (id < 0) {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		}
				
		Body &msh = Bem().surfs[id];
		msh.Reset(Bem().rho, Bem().g);

	 	mainStiffness.Load(Bem().surfs, ids);
		mainView.CalcEnvelope();
		mainSummary.Report(Bem().surfs, id);
		
		mainView.gl.Refresh();
		mainViewData.OnRefresh();
		
		//menuMove.pos_x <<= data.mesh.GetPos().x;
		//menuMove.pos_y <<= data.mesh.GetPos().y;
		//menuMove.pos_z <<= data.mesh.GetPos().z;
		//menuMove.ang_x <<= data.mesh.GetAngle().x;
		//menuMove.ang_y <<= data.mesh.GetAngle().y;
		//menuMove.ang_z <<= data.mesh.GetAngle().z;
		menuProcess.x_g <<= msh.dt.cg.x; 
		menuProcess.y_g <<= msh.dt.cg.y;
		menuProcess.z_g <<= msh.dt.cg.z;
		
		Ma().Status(t_("Model oriented on the initial layout"));
		
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
}

void MainBody::OnUpdateMass() {
	GuiLock __;
	
	try {
		UVector<int> ids = ArrayModel_IdsBody(listLoaded);
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		int id;
		if (num == 0 && listLoaded.GetCount() == 1)
			id = ArrayModel_IdBody(listLoaded, 0);
		else {
		 	id = ArrayModel_IdBody(listLoaded);
			if (id < 0) {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		}
		
		Body &msh = Bem().surfs[id];
		
		double mass = msh.dt.under.volume*Bem().rho;
		menuProcess.mass <<= mass;

		OnUpdate(NONE, true);
		
		Ma().Status(Format(t_("Mass updated to %f kg"), mass));
		
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MainBody::OnScale() {
	GuiLock __;
	
	try {
		UVector<int> ids = ArrayModel_IdsBody(listLoaded);
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		int id;
		if (num == 0 && listLoaded.GetCount() == 1)
			id = ArrayModel_IdBody(listLoaded, 0);
		else {
		 	id = ArrayModel_IdBody(listLoaded);
			if (id < 0) {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		}
				
		Body &msh = Bem().surfs[id];
	
		WaitCursor wait;
		
		double rx = double(~menuProcess.rat_x) - 1;
		double ry = double(~menuProcess.rat_y) - 1;
		double rz = double(~menuProcess.rat_z) - 1;
		
		msh.dt.mesh.Scale(rx, ry, rz, msh.dt.c0);
		msh.dt.cg.Translate(rx*(msh.dt.cg.x - msh.dt.c0.x), ry*(msh.dt.cg.y - msh.dt.c0.y),
						    rz*(msh.dt.cg.z - msh.dt.c0.z));
		
		Ma().Status(Format(t_("Model scaled %f, %f, %f"), rx, ry, rz));
				
		msh.AfterLoad(Bem().rho, Bem().g, NONE, false);
			
		mainStiffness.Load(Bem().surfs, ids);
		mainView.CalcEnvelope();
		mainSummary.Report(Bem().surfs, id);
			
		mainView.gl.Refresh();
		mainViewData.OnRefresh();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}		
}

void MainBody::OnInertia() {
	GuiLock __;
	
	try {
		UVector<int> ids = ArrayModel_IdsBody(listLoaded);
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		int id;
		if (num == 0 && listLoaded.GetCount() == 1)
			id = ArrayModel_IdBody(listLoaded, 0);
		else {
		 	id = ArrayModel_IdBody(listLoaded);
			if (id < 0) {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		}
				
		Body &mesh = Bem().surfs[id];
	
		WithMenuBodyProcessInertia<TopWindow> dialog;
		CtrlLayout(dialog);
		
		double volume = mesh.dt.mesh.volume;
		dialog.volume = volume;
		
		dialog.arrayVol.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, dialog.arrayVol);};
		dialog.arrayVol.SetLineCy(EditField::GetStdHeight());
	
		dialog.arraySurf.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, dialog.arraySurf);};
		dialog.arraySurf.SetLineCy(EditField::GetStdHeight());

		auto Action = [&](bool vol, ArrayCtrl &array, const Point3D &cg) {
			Point3D c0(~dialog.x_0, ~dialog.y_0, ~dialog.z_0);
			
			Matrix3d inertia3;
			mesh.dt.mesh.GetInertia33(inertia3, c0, vol, false);
			MatrixXd inertia6;
			mesh.dt.mesh.GetInertia66(inertia6, inertia3, cg, c0, false);
			if (dialog.opMass == 2)
				inertia6 *= volume;
			else if (dialog.opMass < 2)
				inertia6 *= double(~dialog.mass);
			
			array.Reset();
			array.SetLineCy(EditField::GetStdHeight()).MultiSelect();
			if (dialog.opMass < 3) {
				for (int i = 0; i < 6; ++i)
					array.AddColumn(BEM::StrDOF(i));
	
				for (int r = 0; r < 6; ++r)		
					for (int c = 0; c < 6; ++c)
						array.Set(r, c, FDS(inertia6(r, c), 8, false));
			} else {
				for (int i = 3; i < 6; ++i)
					array.AddColumn(BEM::StrDOF(i));
	
				for (int r = 0; r < 3; ++r)	{	
					for (int c = 0; c < 3; ++c) {
						double val = inertia3(r, c);
						int sign = Sign(val);
						array.Set(r, c, FDS(sign*sqrt(abs(val)), 8, false));
					}
				}
			}
		};
		auto ActionV = [&]() {
			Action(true, dialog.arrayVol, Point3D(~dialog.x_g_v, ~dialog.y_g_v, ~dialog.z_g_v));
		};
		auto ActionS = [&](){
			Action(false, dialog.arraySurf, Point3D(~dialog.x_g_s, ~dialog.y_g_s, ~dialog.z_g_s));
		};
				
		dialog.Title(t_("Mass matrices obtained from mesh"));
		dialog.opC0 = 0;
		dialog.opMass = 0;

		dialog.x_0.WhenEnter = [&]{ActionV();ActionS();};
		dialog.y_0.WhenEnter = [&]{ActionV();ActionS();};
		dialog.z_0.WhenEnter = [&]{ActionV();ActionS();};
		
		dialog.mass.WhenEnter = [&]{ActionV();ActionS();};

		auto opC0_WhenAction = [&](bool action) {
			dialog.x_0.Enable(dialog.opC0 != 0);
			dialog.y_0.Enable(dialog.opC0 != 0);
			dialog.z_0.Enable(dialog.opC0 != 0);
			if (dialog.opC0 == 0) {
				dialog.x_0 = mesh.dt.c0.x;
				dialog.y_0 = mesh.dt.c0.y;
				dialog.z_0 = mesh.dt.c0.z;
			}
			if (action) {
				ActionV();
				ActionS();
			}
		};
		auto opCG_v_WhenAction = [&](bool action) {
			Point3D c = mesh.dt.mesh.GetCentreOfBuoyancy();
			dialog.x_g_v = c.x;
			dialog.y_g_v = c.y;
			dialog.z_g_v = c.z;
			if (action)
				ActionV();
		};
		auto opCG_s_WhenAction = [&](bool action) {
			Point3D c = mesh.dt.mesh.GetCentreOfGravity_Surface();
			dialog.x_g_s = c.x;
			dialog.y_g_s = c.y;
			dialog.z_g_s = c.z;
			if (action)
				ActionS();
		};
		auto opMass_WhenAction = [&](bool action) {
			dialog.mass.Enable(dialog.opMass == 1);
			if (dialog.opMass == 0) 
				dialog.mass = mesh.GetMass();
			else if (dialog.opMass > 1)
				dialog.mass = Null;
			else {
				if (IsNull(dialog.mass))
					dialog.mass = mesh.GetMass();
			}
			String str;
			if (dialog.opMass == 3) 
				str = t_("Radii of gyration");
			else if (dialog.opMass == 2)
				str = t_("Volume moments of inertia");
			else				
				str = t_("Moments of inertia");
			
			dialog.labInertiaV.SetLabel(str);
			dialog.labInertiaS.SetLabel(str);
				
			if (action) {
				ActionV();
				ActionS();
			}
		};
				
		dialog.opC0.WhenAction   = [&]{opC0_WhenAction(true);};
		dialog.opMass.WhenAction = [&]{opMass_WhenAction(true);};
		
		opC0_WhenAction(false);
		opCG_v_WhenAction(false);
		opCG_s_WhenAction(false);
		opMass_WhenAction(false);
		
		ActionV();
		ActionS();
		
		dialog.update	<< [&] {ActionV();ActionS();};		
		dialog.ok		<< [&] {dialog.Close();};
		dialog.Execute();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}		
}
	
void MainBody::OnUpdate(Action action, bool fromMenuProcess) {
	GuiLock __;
	
	try {
		if (fromMenuProcess) {
			menuMove.x_0 <<= ~menuProcess.x_0;
			menuMove.y_0 <<= ~menuProcess.y_0;
			menuMove.z_0 <<= ~menuProcess.z_0;
		} else {
			menuProcess.x_0 <<= ~menuMove.x_0;
			menuProcess.y_0 <<= ~menuMove.y_0;
			menuProcess.z_0 <<= ~menuMove.z_0;
		}
		
		UVector<int> ids = ArrayModel_IdsBody(listLoaded);
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		int id;
		if (num == 0 && listLoaded.GetCount() == 1)
			id = ArrayModel_IdBody(listLoaded, 0);
		else {
		 	id = ArrayModel_IdBody(listLoaded);
			if (id < 0) {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		}
				
		Body &msh = Bem().surfs[id];

		double mass = ~menuProcess.mass;
		double x_g = ~menuProcess.x_g;
		double y_g = ~menuProcess.y_g;
		double z_g = ~menuProcess.z_g;
		double x_0 = ~menuProcess.x_0;
		double y_0 = ~menuProcess.y_0;
		double z_0 = ~menuProcess.z_0;
		double t_x = ~menuMove.t_x;
		double t_y = ~menuMove.t_y;
		double t_z = ~menuMove.t_z;
		double a_x = ~menuMove.a_x;
		double a_y = ~menuMove.a_y;
		double a_z = ~menuMove.a_z;

		if (action == NONE && (IsNull(mass) || IsNull(x_g) || IsNull(y_g) || IsNull(z_g))) {
			BEM::PrintError(t_("Please fill CG data"));
			return;
		}
		if (action == NONE && (IsNull(mass) || IsNull(x_0) || IsNull(y_0) || IsNull(z_0))) {
			BEM::PrintError(t_("Please fill centre of rotation data"));
			return;
		}
				
		if (action == MOVE && (IsNull(t_x) || IsNull(t_y) || IsNull(t_z))) {
			BEM::PrintError(t_("Please fill translation data"));
			return;
		}
		if (action == ROTATE && (IsNull(a_x) || IsNull(a_y) || IsNull(a_z))) {
			BEM::PrintError(t_("Please fill rotation data"));
			return;
		}
		
		WaitCursor wait;

		if (action == MOVE) {
			msh.dt.cg.Translate(t_x, t_y, t_z);
			msh.dt.mesh.Translate(t_x, t_y, t_z);
			videoCtrl.AddReg(Point3D(t_x, t_y, t_z));
			
			Ma().Status(Format(t_("Model moved %f, %f, %f"), t_x, t_y, t_z));
			
		} else if (action == ROTATE) {
			msh.dt.cg.Rotate(ToRad(a_x), ToRad(a_y), ToRad(a_z), x_0, y_0, z_0);
			msh.dt.mesh.Rotate(ToRad(a_x), ToRad(a_y), ToRad(a_z), x_0, y_0, z_0);
			videoCtrl.AddReg(Point3D(a_x, a_y, a_z), Point3D(x_0, y_0, z_0));
			
			if (~menuMove.opZArchimede) {
				double dz = 0;
				Surface under;
				msh.dt.mesh.TranslateArchimede(msh.GetMass(), Bem().rho, dz, under);
				if (!IsNull(dz)) {
					msh.dt.cg.Translate(0, 0, dz);
					videoCtrl.AddReg(Point3D(0, 0, dz));
				} else
					BEM::PrintError(t_("Problem readjusting the Z value to comply with displacement"));
				
				Ma().Status(Format(t_("Model rotated %f, %f, %f deg. around %f, %f, %f, and translated vertically %f m to comply with displacement"), a_x, a_y, a_z, x_0, y_0, z_0, dz));
			} else
				Ma().Status(Format(t_("Model rotated %f, %f, %f around %f, %f, %f"), a_x, a_y, a_z, x_0, y_0, z_0));
		} else if (action == NONE) {
			msh.SetMass(mass);
			msh.dt.cg.Set(x_g, y_g, z_g);
			msh.dt.c0.Set(x_0, y_0, z_0);
		}
		
		menuProcess.x_g <<= msh.dt.cg.x;
		menuProcess.y_g <<= msh.dt.cg.y;
		menuProcess.z_g <<= msh.dt.cg.z;
		
		//menuMove.pos_x <<= data.mesh.GetPos().x;
		//menuMove.pos_y <<= data.mesh.GetPos().y;
		//menuMove.pos_z <<= data.mesh.GetPos().z;
		//menuMove.ang_x <<= data.mesh.GetAngle().x;
		//menuMove.ang_y <<= data.mesh.GetAngle().y;
		//menuMove.ang_z <<= data.mesh.GetAngle().z;
		
		msh.AfterLoad(Bem().rho, Bem().g, action == NONE, false, mainStiffness.opMassBuoy);
		
	 	mainStiffness.Load(Bem().surfs, ids);
		mainView.CalcEnvelope();
		mainSummary.Report(Bem().surfs, id);
		
		mainView.gl.Refresh();
		mainViewData.OnRefresh();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
}

void MainBody::OnArchimede() {
	GuiLock __;
	
	try {
		UVector<int> ids = ArrayModel_IdsBody(listLoaded);
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		int id;
		if (num == 0 && listLoaded.GetCount() == 1)
			id = ArrayModel_IdBody(listLoaded, 0);
		else {
		 	id = ArrayModel_IdBody(listLoaded);
			if (id < 0) {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		}
		
		Body &msh = Bem().surfs[id];
		
		double dz = 0.5, droll = 0.5, dpitch = 0.5;
		Surface under;
		if (!msh.dt.mesh.Archimede(msh.GetMass(), msh.dt.cg, msh.dt.c0, Bem().rho, Bem().g, dz, droll, dpitch, under))
			BEM::PrintError(t_("Problem readjusting the Z, roll and pitch values to comply with buoyancy"));
		
		menuProcess.x_g <<= msh.dt.cg.x;
		menuProcess.y_g <<= msh.dt.cg.y;
		menuProcess.z_g <<= msh.dt.cg.z;
				
		//menuMove.pos_x <<= data.mesh.GetPos().x;
		//menuMove.pos_y <<= data.mesh.GetPos().y;
		//menuMove.pos_z <<= data.mesh.GetPos().z;
		//const Point3D &ang = data.mesh.GetAngle();
		//menuMove.ang_x <<= data.mesh.GetAngle().x;
		//menuMove.ang_y <<= data.mesh.GetAngle().y;
		//menuMove.ang_z <<= data.mesh.GetAngle().z;
		
		msh.AfterLoad(Bem().rho, Bem().g, false, false);
		
	 	mainStiffness.Load(Bem().surfs, ids);
		mainView.CalcEnvelope();
		mainSummary.Report(Bem().surfs, id);
		
		mainView.gl.Refresh();
		mainViewData.OnRefresh();
		
		Ma().Status(Format(t_("Model rotated %f, %f, 0 deg and translated vertically %f m to comply with displacement"), droll, dpitch, dz));
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}			
}			
	

void MainBody::OnAddPanel() {
	GuiLock __;
	
	if (IsNull(menuEdit.edit_x) || IsNull(menuEdit.edit_y) || IsNull(menuEdit.edit_z)) {
		BEM::PrintError(t_("Panel position has to be set"));
		return;
	}
	if (IsNull(menuEdit.edit_size) || double(~menuEdit.edit_size) <= 0) {
		BEM::PrintError(t_("Wrong mesh size"));
		return;
	}
	if (IsNull(menuEdit.panWidthX) || IsNull(menuEdit.panWidthY)) {
		BEM::PrintError(t_("Panel width and height has to be set"));
		return;
	}
	
	WaitCursor waitcursor;
	mainView.gl.Disable();
	try {
		Bem().AddFlatRectangle(~menuEdit.edit_x, ~menuEdit.edit_y, ~menuEdit.edit_z, ~menuEdit.edit_size, 
							 ~menuEdit.panWidthX, ~menuEdit.panWidthY);
		
		Body &msh = Bem().surfs[Bem().surfs.size()-1];
		msh.dt.name = t_("Panel");
		msh.dt.fileName =  "";
		
		msh.AfterLoad(Bem().rho, Bem().g, false, true);
		
		msh.Report(Bem().rho);
		AddRow(msh);
		After();
		mainViewData.OnAddedModel(mainView);
		OnOpt();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
	mainView.gl.Enable();
}

void MainBody::OnAddRevolution() {
	GuiLock __;
	
	UVector<Pointf> vals;
	for (int r = 0; r < menuEdit.revolutionList.GetCount(); ++r) {
		Pointf &val = vals.Add();
		val.x = ScanDouble(AsString(menuEdit.revolutionList.Get(r, 0)));
		if (IsNull(val.x)) {
			BEM::PrintError(Format(t_("Incorrect data in row %d, col %d"), r, 0));
			return;
		}
		val.y = ScanDouble(AsString(menuEdit.revolutionList.Get(r, 1)));
		if (IsNull(val.x)) {
			BEM::PrintError(Format(t_("Incorrect data in row %d, col %d"), r, 1));
			return;
		}
	}
	if (vals.size() < 2) {
		BEM::PrintError(t_("Unsufficient value number in list"));
		return;
	}
	
	WaitCursor waitcursor;
	mainView.gl.Disable();
	try {
		Bem().AddRevolution(~menuEdit.edit_x, ~menuEdit.edit_y, ~menuEdit.edit_z, ~menuEdit.edit_size, vals);
		
		Body &msh = Bem().surfs[Bem().surfs.size()-1];
		msh.dt.name = t_("Revolution");
		msh.dt.fileName =  "";
		
		msh.AfterLoad(Bem().rho, Bem().g, false, true);
		
		msh.Report(Bem().rho);
		AddRow(msh);
		After();
		mainViewData.OnAddedModel(mainView);
		OnOpt();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
	mainView.gl.Enable();
}

void MainBody::OnAddPolygonalPanel() {
	GuiLock __;
	
	UVector<Pointf> vals;
	for (int r = 0; r < menuEdit.polynomialList.GetCount(); ++r) {
		Pointf &val = vals.Add();
		val.x = ScanDouble(AsString(menuEdit.polynomialList.Get(r, 0)));
		if (IsNull(val.x)) {
			BEM::PrintError(Format(t_("Incorrect data in row %d, col %d"), r, 0));
			return;
		}
		val.y = ScanDouble(AsString(menuEdit.polynomialList.Get(r, 1)));
		if (IsNull(val.x)) {
			BEM::PrintError(Format(t_("Incorrect data in row %d, col %d"), r, 1));
			return;
		}
	}
	if (vals.size() < 3) {
		BEM::PrintError(t_("Unsufficient value number in list"));
		return;
	}
	
	WaitCursor waitcursor;
	mainView.gl.Disable();
	try {
		Bem().AddPolygonalPanel(~menuEdit.edit_x, ~menuEdit.edit_y, ~menuEdit.edit_z, ~menuEdit.edit_size, vals);
		
		Body &msh = Bem().surfs[Bem().surfs.size()-1];
		msh.dt.name = t_("Polynomial");
		msh.dt.fileName =  "";
		
		msh.AfterLoad(Bem().rho, Bem().g, false, true);
		
		msh.Report(Bem().rho);
		AddRow(msh);
		After();
		mainViewData.OnAddedModel(mainView);
		OnOpt();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
	mainView.gl.Enable();
}

void MainBody::OnAddWaterSurface(char c) {
	GuiLock __;

	try {
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		int id;
		if (num == 0 && listLoaded.GetCount() == 1)
			id = ArrayModel_IdBody(listLoaded, 0);
		else {
		 	id = ArrayModel_IdBody(listLoaded);
			if (id < 0) {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		}
			
		WaitCursor waitcursor;
		mainView.gl.Disable();
	
		Bem().AddWaterSurface(id, c);
		
		Body &nw = Last(Bem().surfs);
		AddRow(nw);
		After();
		mainViewData.OnAddedModel(mainView);
		OnOpt();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
	mainView.gl.Enable();
}

void MainBody::OnHealing(bool basic) {
	GuiLock __;
	
	try {
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		int id;
		if (num == 0 && listLoaded.GetCount() == 1)
			id = ArrayModel_IdBody(listLoaded, 0);
		else {
		 	id = ArrayModel_IdBody(listLoaded);
			if (id < 0) {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		}
		
		WaitCursor waitcursor;
		Progress progress(t_("Healing mesh file..."), 100); 
		mainView.gl.Disable();
		
		Bem().HealingBody(id, basic, [&](String str, int _pos) {progress.SetText(str); progress.SetPos(_pos); return progress.Canceled();});
		
		UVector<int> ids = ArrayModel_IdsBody(listLoaded);
	 	mainStiffness.Load(Bem().surfs, ids);
		mainView.CalcEnvelope();
		mainSummary.Report(Bem().surfs, id);
		
		mainView.gl.Refresh();
		mainViewData.OnRefresh();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
	mainView.gl.Enable();
}

void MainBody::OnOrientSurface() {
	GuiLock __;
	
	try {
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		int id;
		if (num == 0 && listLoaded.GetCount() == 1)
			id = ArrayModel_IdBody(listLoaded, 0);
		else {
		 	id = ArrayModel_IdBody(listLoaded);
			if (id < 0) {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		}
		
		WaitCursor waitcursor;
		Progress progress(t_("Orienting mesh surface..."), 100); 
		mainView.gl.Disable();
		
		Bem().OrientSurface(id, [&](String str, int _pos) {progress.SetText(str); progress.SetPos(_pos); return progress.Canceled();});
		
		Bem().surfs[id].AfterLoad(Bem().rho, Bem().g, false, false);
		
		UVector<int> ids = ArrayModel_IdsBody(listLoaded);
	 	mainStiffness.Load(Bem().surfs, ids);
		mainView.CalcEnvelope();
		mainSummary.Report(Bem().surfs, id);
		
		mainView.gl.Refresh();
		mainViewData.OnRefresh();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
	mainView.gl.Enable();
}

void MainBody::OnImage(int axis) {
	GuiLock __;
	
	String saxis = (axis == 0) ? "X" : ((axis == 1) ? "Y" : "Z");

	try {
		UVector<int> ids = ArrayModel_IdsBody(listLoaded);
		int id = ArrayModel_IdBody(listLoaded);
		if (id < 0) {
			BEM::PrintError(t_("Please select a model to process"));
			return;
		}
		
		WaitCursor waitcursor;
		
		Body &msh = Bem().surfs[id];

		msh.SetMass(~menuProcess.mass);
		if (axis == 0)
			msh.dt.cg.x = -msh.dt.cg.x;
		else if (axis == 1)
			msh.dt.cg.y = -msh.dt.cg.y;
		else
			msh.dt.cg.z = -msh.dt.cg.z;
		
		msh.dt.mesh.Image(axis);
	
		msh.AfterLoad(Bem().rho, Bem().g, false, false);
		
	 	mainStiffness.Load(Bem().surfs, ids);
		mainView.CalcEnvelope();
		mainSummary.Report(Bem().surfs, id);
		
		mainView.gl.Refresh();
		mainViewData.OnRefresh();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MainBody::OnRemove() {
	WaitCursor waitcursor;

	OnRemoveSelected(true);
}

void MainBody::OnRemoveSelected(bool all) {	
	bool selected = false;
	
	UVector<int> sel = ArrayCtrlSelectedGet(listLoaded);
	
	for (int r = listLoaded.GetCount()-1; r >= 0; --r) {
		if (all || Find(sel, r) >= 0) {
			int id = ArrayModel_IdBody(listLoaded, r);
			Bem().RemoveBody(id);
			listLoaded.Remove(r);
			selected = true;
		}
	}	// Only one available => directly selected
	if (!selected && listLoaded.GetCount() == 1) {
		int id = ArrayModel_IdBody(listLoaded, 0);
		Bem().RemoveBody(id);
		listLoaded.Remove(0);
		selected = true;		
	}	
	if (!selected) {
		BEM::PrintError(t_("No model selected"));
		return;
	}

	UVector<int> ids = ArrayModel_IdsBody(listLoaded);
	mainStiffness.Load(Bem().surfs, ids);
	mainViewData.ReLoad(mainView);
	
	mainGZ.Clear(false);
	
	After();
}

void MainBody::OnJoin() {
	GuiLock __;
	
	try {	
		//bool selected = false;
		int idDest = Null;
		for (int r = 0; r < listLoaded.GetCount(); ++r) {
			if (listLoaded.IsSelected(r)) {
				if (IsNull(idDest))
					idDest = ArrayModel_IdBody(listLoaded, r);
				else
					idDest = min(idDest, ArrayModel_IdBody(listLoaded, r));
			}
		}
		if (IsNull(idDest)) {
			BEM::PrintError(t_("No model joined"));
			return;
		}
		
		WaitCursor waitcursor;
		
		for (int r = listLoaded.GetCount()-1; r >= 0; --r) {
			if (listLoaded.IsSelected(r)) {
				int id = ArrayModel_IdBody(listLoaded, r);
				if (idDest != id) {
					Bem().JoinBody(idDest, id);
					RemoveRow(r);
					//selected = true;
				}
			}
		}	
	
		UVector<int> ids = ArrayModel_IdsBody(listLoaded);
		mainStiffness.Load(Bem().surfs, ids);
		mainViewData.ReLoad(mainView);
		
		After();
	} catch (Exc e) {
		mainView.gl.Enable();
		BEM::PrintError(DeQtfLf(e));
	}
}

void MainBody::AddRow(const Body &msh) {
	ArrayModel_Add(listLoaded, msh.GetBodyStr(), msh.dt.name, msh.dt.fileName, msh.dt.GetId(),
					optionsPlot, [&] {mainView.gl.Refresh();});
}
		
void MainBody::RemoveRow(int row) {
	listLoaded.Remove(row);
}	
	
void MainBody::OnSplit() {
	GuiLock __;
	
	String file = ~menuOpen.file;
		
	try {
		Progress progress(t_("Splitting mesh file..."), 100); 
		
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num == 0 && listLoaded.GetCount() != 1) {
			BEM::PrintError(t_("No model selected"));
			return;
		}
		if (num > 1) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		WaitCursor waitcursor;
				
		UVector<int> idsmesh;
		int row = -1;
		for (row = listLoaded.GetCount()-1; row >= 0; --row) {
			if (listLoaded.IsSelected(row)) 
				break;
		}	// Only one available => directly selected
		if (row < 0 && listLoaded.GetCount() == 1) 
			row = 0;
	
		if (idsmesh.size() == 1) {
			BEM::PrintError(t_("The mesh is monolithic so it cannot be automatically split"));
			return;
		}
		int id = ArrayModel_IdBody(listLoaded, row);
		String fileName = Bem().surfs[id].dt.fileName;
		String name = Bem().surfs[id].dt.name;
		idsmesh = Bem().SplitBody(id, [&](String str, int _pos) {
			progress.SetText(str); 
			progress.SetPos(_pos); 
			progress.Refresh();
			return true;
		});	
		
		RemoveRow(row);
		
		for (int i = 0; i < idsmesh.size(); ++i) {
			int idm = idsmesh[i];
			Body &msh = Bem().surfs[idm];
			
			mainTab.Set(mainView);
			
			msh.Report(Bem().rho);
	
			msh.dt.fileName = fileName;
			msh.dt.name = name + S("-") + FormatInt(i+1);		
			AddRow(msh);
		}
		
		mainViewData.ReLoad(mainView);
				
		After();
	} catch (Exc e) {
		mainView.gl.Enable();
		BEM::PrintError(DeQtfLf(e));
	}
}

void MainBody::UpdateButtons() {
	int numrow = listLoaded.GetCount();
	int numsel = ArrayCtrlSelectedGetCount(listLoaded);
	menuOpen.butRemove.Enable(numrow > 0);
	menuOpen.butRemoveSelected.Enable(numsel > 0);
	menuOpen.butJoin.Enable(numsel > 1);
	menuOpen.butSplit.Enable(numsel == 1 || numrow == 1);
	menuOpen.dropExport.Enable(numsel >= 1);
	menuOpen.butExport.Enable(numsel >= 1);
	menuOpen.optBodyType.Enable(numsel == 1);
	menuOpen.symX.Enable(numsel == 1);
	menuOpen.symY.Enable(numsel == 1);
	menuOpen.labBody.Enable(numsel == 1);
	menuOpen.labSymmetry.Enable(numsel == 1);
	
	menuProcess.butUpdateCg.Enable(numsel == 1 || numrow == 1);
	
	menuMove.butUpdatePos.Enable(numsel == 1 || numrow == 1);
	menuMove.butUpdateAng.Enable(numsel == 1 || numrow == 1);
	menuMove.butReset.Enable(numsel == 1 || numrow == 1);
	
	menuProcess.butImageX.Enable(numsel == 1 || numrow == 1);
	menuProcess.butImageY.Enable(numsel == 1 || numrow == 1);
	menuProcess.butImageZ.Enable(numsel == 1 || numrow == 1);
	menuProcess.butSimplify.Enable(numsel == 1 || numrow == 1);
	menuProcess.butFullHealing.Enable(numsel == 1 || numrow == 1);
	menuProcess.butOrientSurface.Enable(numsel == 1 || numrow == 1);
	menuProcess.butWaterFill.Enable(numsel == 1 || numrow == 1);
	menuProcess.butWaterPlane.Enable(numsel == 1 || numrow == 1);
	menuProcess.butHull.Enable(numsel == 1 || numrow == 1);
}

void MainBody::After() {
	UpdateButtons();

	mainView.CalcEnvelope();	

	mainSummary.Clear();
	for (int row = 0; row < listLoaded.GetCount(); ++row) {
		int id = ArrayModel_IdBody(listLoaded, row);
		mainSummary.Report(Bem().surfs, id);
	}		

	mainTab.WhenSet();
	
	mainView.gl.Enable();
	mainView.gl.Refresh();
}
	
void MainBody::Jsonize(JsonIO &json) {
	if (json.IsLoading()) {
		menuPlot.showBody = Null;	
		menuPlot.showNormals = Null;	
		menuPlot.showWaterLevel = Null;	
		menuPlot.showSkewed = Null;
		menuPlot.showFissure = Null;
		menuPlot.showMultiPan = Null;
		menuPlot.showAxis = Null;
		menuPlot.showLimits = Null;
		menuPlot.showCg = Null;
		menuPlot.showCb = Null;
		menuPlot.showCr = Null;
		menuPlot.showLines = Null;
		menuPlot.showSel = Null;
		menuOpen.optBodyType = Null;
		menuOpen.opClean = false;
		menuMove.opZArchimede = Null;
		dropExportId = 2;
	} else
		dropExportId = menuOpen.dropExport.GetIndex();
	json
		("menuOpen_file", menuOpen.file)
		("menuOpen_saveFolder", saveFolder)
		("menuOpen_dropExport", dropExportId)
		("menuOpen_optBodyType", menuOpen.optBodyType)
		("menuOpen_opClean", menuOpen.opClean)
		("menuPlot_showBody", menuPlot.showBody)		
		("menuPlot_showNormals", menuPlot.showNormals)	
		("menuPlot_showSkewed", menuPlot.showSkewed)	
		("menuPlot_showFissure", menuPlot.showFissure)	
		("menuPlot_showMultiPan", menuPlot.showMultiPan)	
		("menuPlot_showAxis", menuPlot.showAxis)
		("menuPlot_showLimits", menuPlot.showLimits)
		("menuPlot_showCg", menuPlot.showCg)
		("menuPlot_showCb", menuPlot.showCb)
		("menuPlot_showCr", menuPlot.showCr)
		("menuPlot_showSel", menuPlot.showSel)
		("menuPlot_showLines", menuPlot.showLines)
		("menuPlot_showUnderwater", menuPlot.showUnderwater)
		("menuPlot_showWaterLevel", menuPlot.showWaterLevel)
		("menuPlot_backColor", menuPlot.backColor)
		("menuPlot_lineThickness", menuPlot.lineThickness)	
		("mainStiffness", mainStiffness)
		("menuMove_opZArchimede", menuMove.opZArchimede)
		("mainGZ", mainGZ)
		("videoCtrl", videoCtrl)
	;
	if (json.IsLoading()) {
		if (IsNull(dropExportId) || dropExportId < 0)
			dropExportId = 0;
	}	
}

void MainSummaryBody::Report(const UArray<Body> &surfs, int id) {
	const Body &msh = surfs[id];
	String name = msh.dt.name;
	
	if (array.GetColumnCount() == 0)
		array.AddColumn("Param");
	if (id >= array.GetColumnCount()-1)
		array.AddColumn(Format("#%d %s", id+1, name));
	int row = 0;
	int col = id + 1;

	int idColor;
	Upp::Color backColorBody;
	idColor = msh.dt.mesh.VolumeMatch(Bem().volWarning/100., Bem().volError/100.);
	if (idColor == -2 || msh.dt.mesh.surface < 0)
		backColorBody = Upp::Color(255, 165, 158);	// Light red
	else if (idColor == -1)
		backColorBody = Upp::Color(255, 255, 150);	// Light yellow
	
	Upp::Color backColorUnder;
	idColor = msh.dt.under.VolumeMatch(Bem().volWarning/100., Bem().volError/100.);
	if (idColor == -2)
		backColorUnder = Upp::Color(255, 165, 158);
	else if (idColor == -1)
		backColorUnder = Upp::Color(255, 255, 150);
									
	bool healing = msh.dt.mesh.healing;

	const MainBody &mn = GetDefinedParent<MainBody>(this);
	Upp::Color color = ArrayModel_GetColor(mn.listLoaded, id);//GetColorId(id);
	::Color textColor = Black();
	if (Grayscale(color) < 150)
		textColor = White();
														
	array.Set(row, 0, t_("File"));				array.Set(row++, col, AttrText(msh.dt.fileName).Paper(color).Ink(textColor).Bold()); 
	array.Set(row, 0, t_("Name"));				array.Set(row++, col, name + (healing ? (S(" ") + t_("(healed)")) : ""));
	array.Set(row, 0, t_("Format"));			array.Set(row++, col, msh.GetBodyStr());	
	
	array.Set(row, 0, t_("# Panels"));			array.Set(row++, col, msh.dt.mesh.panels.size());
	array.Set(row, 0, t_("# Nodes"));			array.Set(row++, col, msh.dt.mesh.nodes.size());

	array.Set(row, 0, t_("Surface [m²]"));		array.Set(row++, col, FDS(msh.dt.mesh.surface, 8, false));
	array.Set(row, 0, t_("Volume [m3] Vavg (Vx,Vy,Vz)"));		  array.Set(row++, col, AttrText(Format(t_("%s (%s, %s, %s)"), 
														FDS(msh.dt.mesh.volume,  10, false),
														FDS(msh.dt.mesh.volumex, 10, false),
														FDS(msh.dt.mesh.volumey, 10, false),
														FDS(msh.dt.mesh.volumez, 10, false))).Paper(backColorBody));
	
	array.Set(row, 0, t_("Wetted surface [m²]"));array.Set(row++, col, FDS(msh.dt.under.surface, 10, false));
	array.Set(row, 0, t_("Immersed volume [m3] Vavg (Vx,Vy,Vz)")); array.Set(row++, col, AttrText(Format(t_("%s (%s, %s, %s)"), 
														FDS(msh.dt.under.volume,  10, false),
														FDS(msh.dt.under.volumex, 10, false),
														FDS(msh.dt.under.volumey, 10, false),
														FDS(msh.dt.under.volumez, 10, false))).Paper(backColorUnder));
	array.Set(row, 0, t_("Displacement [kg]")); array.Set(row++, col, FDS(msh.dt.under.volume*Bem().rho, 10, false));
	array.Set(row, 0, t_("Cg [m]"));			
	if (!IsNull(msh.dt.cg))
		array.Set(row++, col, Format(t_("%s, %s, %s"),
														FDS(msh.dt.cg.x, 10, false),			
														FDS(msh.dt.cg.y, 10, false),
														FDS(msh.dt.cg.z, 10, false)));
	else
		array.Set(row++, col, "-");
			
	array.Set(row, 0, t_("Cb [m]"));
	if (!IsNull(msh.dt.cb))	
		array.Set(row++, col, Format(t_("%s, %s, %s"),  FDS(msh.dt.cb.x, 10, false),			
														FDS(msh.dt.cb.y, 10, false),
														FDS(msh.dt.cb.z, 10, false)));
	else 
		array.Set(row++, col, "-");
	
	array.Set(row, 0, t_("C0 [m]"));			array.Set(row++, col, Format(t_("%s, %s, %s"),
														FDS(msh.dt.c0.x, 10, false),			
														FDS(msh.dt.c0.y, 10, false),
														FDS(msh.dt.c0.z, 10, false)));

	array.Set(row, 0, t_("GMroll [m]"));
	double gmpitch = msh.GMpitch(Bem().rho, Bem().g); 	
	if (IsNum(gmpitch))
 		array.Set(row++, col, FDS(gmpitch, 5, false));
	else
		array.Set(row++, col, "-");
	array.Set(row, 0, t_("GMpitch [m]")); 		
	double gmroll = msh.GMroll(Bem().rho, Bem().g); 	
	if (IsNum(gmroll))
		array.Set(row++, col, FDS(gmroll, 5, false));
	else
		array.Set(row++, col, "-");
	
	array.Set(row, 0, t_("GZ [m]"));
	if (!IsNull(msh.dt.cb) && !IsNull(msh.dt.cg)) {
		Direction3D cgcb = msh.dt.cg - msh.dt.cb;
		double gz = sqrt(sqr(cgcb.x) + sqr(cgcb.y));
		array.Set(row++, col, AttrText(FormatDouble(gz, 4)));	
	} else
		array.Set(row++, col, "-");
												
	array.Set(row, 0, t_("Surface projection Z-axis (waterplane area) [m²]"));	
												array.Set(row++, col, Format(t_("%s - %s = %s"),
														FDS(-msh.dt.projectionPos.z, 10, false),
														FDS(msh.dt.projectionNeg.z,  10, false),
														FDS(msh.dt.projectionPos.z + msh.dt.projectionNeg.z, 10, false)));
	
	array.Set(row, 0, t_("Waterplane geometric centre (centre of flotation) [m]"));
	if (!IsNull(msh.dt.cgZ0surface)) 
		array.Set(row++, col, Format(t_("%s, %s"), FDS(msh.dt.cgZ0surface.x, 10, false),			
												   FDS(msh.dt.cgZ0surface.y, 10, false)));
	else
		array.Set(row++, col, "-");
	array.Set(row, 0, t_("Surface projection X-axis [m²]"));	
												array.Set(row++, col, Format(t_("%s - %s = %s"),
														FDS(-msh.dt.projectionPos.x, 10, false),
														FDS(msh.dt.projectionNeg.x,  10, false),
														FDS(msh.dt.projectionPos.x + msh.dt.projectionNeg.x, 10, false)));
	
	array.Set(row, 0, t_("Surface projection Y-axis [m²]"));	
												array.Set(row++, col, Format(t_("%s - %s = %s"),
														FDS(-msh.dt.projectionPos.y, 10, false),
														FDS(msh.dt.projectionNeg.y,  10, false),
														FDS(msh.dt.projectionPos.y + msh.dt.projectionNeg.y, 10, false)));
	
	array.Set(row, 0, t_("Dimensions [m]"));	array.Set(row++, col, Format(t_("From (%s, %s, %s) to (%s, %s, %s)"),
														FDS(msh.dt.mesh.env.minX, 10, false),
														FDS(msh.dt.mesh.env.minY, 10, false),
														FDS(msh.dt.mesh.env.minZ, 10, false),
														FDS(msh.dt.mesh.env.maxX, 10, false),
														FDS(msh.dt.mesh.env.maxY, 10, false),
														FDS(msh.dt.mesh.env.maxZ, 10, false)));

	//Force6D f = data.under.GetHydrostaticForce(data.c0, Bem().rho, Bem().g);	
	Force6D fcb = msh.dt.under.GetHydrostaticForceCB(msh.dt.c0, msh.dt.cb, Bem().rho, Bem().g);	
	
	array.Set(row, 0, t_("Hydrostatic forces [N]"));
	
	if (!IsNull(msh.dt.cb))
		array.Set(row++, col, AttrText(Format(t_("%s, %s, %s"),
														FDS(fcb[0], 10, false),
														FDS(fcb[1], 10, false),
														FDS(fcb[2], 10, false))).Paper(backColorUnder));
	else
		array.Set(row++, col, "-");

	array.Set(row, 0, t_("Hydrostatic moments [N·m]"));
	if (!IsNull(msh.dt.cb))
		array.Set(row++, col, AttrText(Format(t_("%s, %s, %s"),
														FDS(fcb[3], 10, false),
														FDS(fcb[4], 10, false),
														FDS(fcb[5], 10, false))));							
	else
		array.Set(row++, col, "-");
		
	array.Set(row, 0, t_("Mass [kg]"));			array.Set(row++, col, FormatF(msh.GetMass(), 1));
	
	Force6D fcg = Surface::GetMassForce(msh.dt.c0, msh.dt.cg, msh.GetMass(), Bem().g);	
	
	array.Set(row, 0, t_("Mass moments [N·m]"));
	if (!IsNull(fcg)) 
		array.Set(row++, col, Format(t_("%s, %s, %s"),
														FDS(fcg[3], 10, false),
														FDS(fcg[4], 10, false),
														FDS(fcg[5], 10, false)));
	else
		array.Set(row++, col, "-");
	
	array.Set(row, 0, t_("Mass+Hydrostatics forces [N]"));
	if (!IsNull(fcg) && !IsNull(fcb))
		array.Set(row++, col, AttrText(Format(t_("%s, %s, %s"),
														"0",
														"0",
														FDS(fcg[2]+fcb[2], 10, false))).Paper(backColorUnder));
	else
		array.Set(row++, col, "-");
	
	array.Set(row, 0, t_("Mass+Hydrostatics moments [N·m]"));
	if (!IsNull(fcb))
		array.Set(row++, col, AttrText(Format(t_("%s, %s, %s"),
														FDS(fcg[3]+fcb[3], 10, false),
														FDS(fcg[4]+fcb[4], 10, false),
														"0")));		
	else
		array.Set(row++, col, "-");
														
	array.Set(row++, 0, t_("Stiffness Matrix"));	
	if (msh.dt.C.size() > 0) {
		for (int i = 0; i < 6; ++i) {
			for (int j = 0; j < 6; ++j) {
				if (!Hydro::C_units(i, j).IsEmpty()) {
					array.Set(row, 0, Format(t_("K(%d,%d) [%s]"), i+1, j+1, Hydro::C_units(i, j)));	array.Set(row++, col, Format("%12E", msh.dt.C(i, j)));		
				}
			}
		}
	} else {
		for (int i = 0; i < 6; ++i) {
			for (int j = 0; j < 6; ++j) {
				if (!Hydro::C_units(i, j).IsEmpty()) {
					array.Set(row, 0, Format(t_("K(%d,%d) [%s]"), i+1, j+1, Hydro::C_units(i, j)));	array.Set(row++, col, "-");		
				}
			}
		}
	}
	
	array.Set(row, 0, t_("Healing"));			array.Set(row++, col, healing ? t_("Yes") : t_("No"));

	array.Set(row, 0, t_("# Segments"));		array.Set(row++, col, !healing ? Null : msh.dt.mesh.segments.size());
	array.Set(row, 0, t_("# Seg Waterplane"));	array.Set(row++, col, !healing ? Null : msh.dt.mesh.segWaterlevel.size());
	array.Set(row, 0, t_("# Seg leak"));		array.Set(row++, col, !healing ? Null : msh.dt.mesh.segTo1panel.size());
	array.Set(row, 0, t_("# Seg 3 panels"));	array.Set(row++, col, !healing ? Null : msh.dt.mesh.segTo3panel.size());

//	array.Set(row, 0, t_("# Panels off"));		array.Set(row++, col, !healing ? Null : data.mesh.numUnprocessed);

	array.Set(row, 0, t_("# Triangles"));		array.Set(row++, col, !healing ? Null : msh.dt.mesh.numTriangles);
	array.Set(row, 0, t_("# BiQuads"));			array.Set(row++, col, !healing ? Null : msh.dt.mesh.numBiQuads);
	array.Set(row, 0, t_("# MonoQuads"));		array.Set(row++, col, !healing ? Null : msh.dt.mesh.numMonoQuads);
	array.Set(row, 0, t_("# Dup panels"));		array.Set(row++, col, !healing ? Null : msh.dt.mesh.numDupPan);
	array.Set(row, 0, t_("# Dup nodes"));		array.Set(row++, col, !healing ? Null : msh.dt.mesh.numDupP);
	array.Set(row, 0, t_("# Skewed pan"));		array.Set(row++, col, !healing ? Null : msh.dt.mesh.numSkewed);
	array.Set(row, 0, t_("Avg. side [m]"));		array.Set(row++, col, !healing ? Null : msh.dt.mesh.GetAvgLenSegment());
}

void MainView::Init(MainBody &parent) {
	CtrlLayout(*this);
	main = &parent;
	
	gl.SetEnv(env);
	gl.WhenPaint = THISBACK(OnPaint);	
}
	
void MainView::OnPaint() {
	gl.SetLineThickness(~GetMenuPlot().lineThickness);
	gl.SetBackgroundColor(~GetMenuPlot().backColor);
	
	if (~GetMenuPlot().showAxis && GetMain().listLoaded.GetCount() > 0) 
		gl.PaintAxis(0, 0, 0, env.LenRef()/4.);	
	
	if (~GetMenuPlot().showLimits && Bem().surfs.size() > 0) 
		gl.PaintCuboid(Point3D(env.maxX, env.maxY, env.maxZ), Point3D(env.minX, env.minY, env.minZ), Gray());
	
	for (int row = 0; row < GetMain().listLoaded.GetCount(); ++row) {
		if (ArrayModel_IsVisible(GetMain().listLoaded, row)) {
			int id = ArrayModel_IdBody(GetMain().listLoaded, row);
			if (id < 0)
				throw Exc(t_("Unexpected problem in OnPaint()"));
			
			double len = env.LenRef()/10;
			bool showNormals = ~GetMenuPlot().showNormals && ~GetMenuPlot().showBody;
			bool showNormalsUnderwater = ~GetMenuPlot().showNormals && ~GetMenuPlot().showUnderwater;
			if (~GetMenuPlot().showNormals && !~GetMenuPlot().showBody && !~GetMenuPlot().showUnderwater)
				showNormals = true;
			
			const Upp::Color &color = ArrayModel_GetColor(GetMain().listLoaded, row);
			const Body &msh = Bem().surfs[id];
			
			gl.PaintSurface(msh.dt.mesh, color, ~GetMenuPlot().showBody, 	
				showNormals);
				
			gl.PaintSurface(msh.dt.under, color, ~GetMenuPlot().showUnderwater, 
				showNormalsUnderwater);
			
			if (~GetMenuPlot().showBody)
				gl.PaintSegments(msh.dt.mesh, color);
			
			if (~GetMenuPlot().showSkewed)
				gl.PaintSegments(msh.dt.mesh.skewed, LtRed());
			if (~GetMenuPlot().showFissure)
				gl.PaintSegments(msh.dt.mesh.segTo1panel, LtRed());
			if (~GetMenuPlot().showWaterLevel)
				gl.PaintSegments(msh.dt.under.segWaterlevel, LtBlue());
			if (~GetMenuPlot().showMultiPan)
				gl.PaintSegments(msh.dt.mesh.segTo3panel, Black());
			
			if (~GetMenuPlot().showCb) {
				gl.PaintDoubleAxis(msh.dt.cb, len, LtBlue());
				gl.PaintCube(msh.dt.cb, len/10, LtBlue());
			}
			if (~GetMenuPlot().showCg) {
				gl.PaintDoubleAxis(msh.dt.cg, len, Black());
				gl.PaintCube(msh.dt.cg, len/5, Black());
			}
			if (~GetMenuPlot().showCr) {
				gl.PaintDoubleAxis(msh.dt.c0, len*10, Cyan());
				gl.PaintCube(msh.dt.c0, len/20, Gray());
			}
			if (~GetMenuPlot().showLines) 
				gl.PaintLines(msh.dt.mesh.lines, color);
			
			if (paintSelect) {
				if (~GetMenuPlot().showBody) {
					const UVector<int> &nod = msh.dt.mesh.GetSelNodes();
					for (int in = 0; in < nod.size(); ++in)
						gl.PaintCube(msh.dt.mesh.nodes[nod[in]], len/20, LtBlue());
					const UVector<int> &pan = msh.dt.mesh.GetSelPanels();
					const UVector<Point3D> &nodes = msh.dt.mesh.nodes;
					for (int ip = 0; ip < pan.size(); ++ip) {
						const Panel &panel = msh.dt.mesh.panels[pan[ip]];
						gl.PaintQuad(nodes[panel.id[0]], nodes[panel.id[1]], nodes[panel.id[2]], nodes[panel.id[3]], LtRed(), .2);
					}
				}
				if (~GetMenuPlot().showUnderwater) {
					const UVector<int> &nod = msh.dt.under.GetSelNodes();
					for (int in = 0; in < nod.size(); ++in)
						gl.PaintCube(msh.dt.under.nodes[nod[in]], len/20, LtBlue());
					const UVector<int> &pan = msh.dt.under.GetSelPanels();
					const UVector<Point3D> &nodes = msh.dt.under.nodes;
					for (int ip = 0; ip < pan.size(); ++ip) {
						const Panel &panel = msh.dt.under.panels[pan[ip]];
						gl.PaintQuad(nodes[panel.id[0]], nodes[panel.id[1]], nodes[panel.id[2]], nodes[panel.id[3]], LtRed(), .2);
					}
				}
			}
			if (~GetMenuPlot().showSel && ArrayModel_IsSelected(GetMain().listLoaded, row)) { 
				double minX = msh.dt.mesh.env.minX*.99; double maxX = msh.dt.mesh.env.maxX*1.01;
				double minY = msh.dt.mesh.env.minY*.99; double maxY = msh.dt.mesh.env.maxY*1.01;
				double minZ = msh.dt.mesh.env.minZ*.99; double maxZ = msh.dt.mesh.env.maxZ*1.01;
				
				gl.PaintCuboid(Point3D(maxX, maxY, maxZ), Point3D(minX, minY, minZ), color);
				gl.PaintCube(Point3D(maxX, maxY, minZ), len/10, color);
				gl.PaintCube(Point3D(maxX, minY, minZ), len/10, color);
				gl.PaintCube(Point3D(minX, maxY, minZ), len/10, color);
				gl.PaintCube(Point3D(minX, minY, minZ), len/10, color);
				gl.PaintCube(Point3D(maxX, maxY, maxZ), len/10, color);
				gl.PaintCube(Point3D(maxX, minY, maxZ), len/10, color);
				gl.PaintCube(Point3D(minX, maxY, maxZ), len/10, color);
				gl.PaintCube(Point3D(minX, minY, maxZ), len/10, color);
			}
		}
	}
}

const WithMenuBodyPlot<StaticRect> &MainView::GetMenuPlot() const {
	return main->menuPlot;
}

void MainView::CalcEnvelope() {
	env.Reset();
	//for (int i = 0; i < Bem().surfs.size(); ++i)
	for (int row = 0; row < GetMain().listLoaded.GetCount(); ++row) {
		int id = ArrayModel_IdBody(GetMain().listLoaded, row);
		if (id < 0)
			throw Exc("Unexpected problem in CalcEnvelope()");
		env.MixEnvelope(Bem().surfs[id].dt.mesh.env);
	}
}

void MainBody::LoadDragDrop() {
	GuiLock __;
	
	Sort(filesToDrop);
	for (int i = filesToDrop.size()-1; i > 0; --i)
		if (ToLower(GetFileTitle(filesToDrop[i])) == ToLower(GetFileTitle(filesToDrop[i-1])))
			filesToDrop.Remove(i);
		
	bool followWithErrors = false;
	for (int i = 0; i < filesToDrop.size(); ++i) {
		menuOpen.file <<= filesToDrop[i];
		Status(Format(t_("Loading '%s'"), filesToDrop[i]));
		if (!OnLoad() && !followWithErrors && filesToDrop.size() - i > 1) {
			if (!PromptYesNo(Format(t_("Do you wish to load the pending %d files?"), filesToDrop.size() - i - 1)))
				return;
			followWithErrors = true;
		}
		ProcessEvents();
	}
}
	
void MainBody::DragAndDrop(Point , PasteClip& d) {
	GuiLock __;
	if (IsDragAndDropSource())
		return;
	if (AcceptFiles(d)) {
		filesToDrop = GetFiles(d);
		timerDrop.Set(0, [=] {LoadDragDrop();});
	}
}

bool MainBody::Key(dword key, int ) {
	GuiLock __;
	if (key == K_CTRL_V) {
		filesToDrop = GetFiles(Ctrl::Clipboard());
		timerDrop.Set(0, [=] {LoadDragDrop();});
		return true;
	}
	return false;
}

void MainGZ::Init() {
	CtrlLayout(*this);
	
	if (IsNull(edFrom))
		edFrom <<= -45;	
	if (IsNull(edTo))	
		edTo <<= 45;
	if (IsNull(edDelta))
		edDelta <<= 1;
	if (IsNull(edAngleFrom))
		edAngleFrom <<= 0;
	if (IsNull(edAngleTo))
		edAngleTo <<= 0;
	if (IsNull(edAngleDelta))
		edAngleDelta <<= 30;
	if (IsNull(edTolerance))
		edTolerance <<= 5;
	
	butUpdate.WhenAction = THISBACK(OnUpdate);
	
	splitter.Vert(scatter.SizePos(), array.SizePos());
	splitter.SetPos(7000, 0);
	
	scatter.ShowAllMenus().SetSequentialXAll().SetFastViewX();
	scatter.SetLabelX(t_("Angle [º]"));
	scatter.SetTitle(t_("GZ around Y axis"));
	scatter.SetLegendFillColor(White());
	scatter.SetMargin(140, 140, 50, 80);
	
	Clear(true);
	
	edAngleDelta.WhenAction = [&] {edAngleTo.Enable(!(IsNull(edAngleDelta) || double(~edAngleDelta) == 0));};
}

void MainGZ::OnUpdate() {
	try {
		MainBody &mm = GetDefinedParent<MainBody>(this);
		
		idOpened = ArrayModel_IdBody(mm.listLoaded);
		if (idOpened < 0) {
			BEM::PrintError(t_("Please select a model to process"));
			return;
		}
		
		Clear(true);
		
		double angleFrom = ~edAngleFrom;
		if (IsNull(angleFrom)) {
			BEM::PrintError(t_("Wrong angle from"));
			return;
		}
			
		
		if (double(~edFrom) > double(~edTo)) {
			BEM::PrintError(t_("Wrong Y angle range"));
			return;
		}
		double angleTo, angleDelta;
		if (edAngleTo.IsEnabled()) {
			angleTo = double(~edAngleTo);
			angleDelta =  double(~edAngleDelta);
			if (angleDelta == 0) {
				angleTo = angleFrom;
				angleDelta = 1;
			}
		} else {
			angleTo = angleFrom;
			angleDelta = 1;
		}
		if (angleFrom > angleTo) {
			BEM::PrintError(t_("Wrong Z angle range"));
			return;
		}
		Body &msh = Bem().surfs[idOpened];	
	
		int numAngle = 1 + int((angleTo - angleFrom)/angleDelta);
		
		Progress progress(t_("GZ calculation..."), 100*numAngle); 
		
		datagz.Clear();
		dataMoment.Clear();
		
		scatter.SetTitle(Format(t_("GZ around Y axis at (%.2f, %.2f, %.2f)"), msh.dt.c0.x, msh.dt.c0.y, msh.dt.c0.z));
		
		String errors;
		int iangle = 0;
		for (double angle = double(~edAngleFrom); angle <= angleTo; angle += angleDelta, iangle++) {
			UVector<double> &dgz = datagz.Add();
			UVector<double> &dMoment = dataMoment.Add();
			UVector<double> vol, disp, wett, wplane, draft;
			UVector<Point3D> cb, cg;
			
			msh.GZ(~edFrom, ~edTo, ~edDelta, angle, Bem().rho, Bem().g, double(~edTolerance)/100., 
				[&](String, int pos)->bool {
					progress.SetPos(pos + 100*iangle);
					return !progress.Canceled();
				}, dangle, dgz, dMoment, vol, disp, wett, wplane, draft, cb, cg, errors);
			
			Upp::Color color = ScatterDraw::GetNewColor(iangle);
			scatter.AddSeries(dangle, dgz).NoMark().Legend(Format(t_("GZ %.1f"), angle))
				   .Units(t_("m"), t_("sec")).Stroke(2, color).Dash(LINE_SOLID);	
			scatter.AddSeries(dangle, dMoment).NoMark()./*MarkStyle<CircleMarkPlot>().*/MarkWidth(8).MarkColor(color)
				   .Legend(Format(t_("Healing lever %.1f"), angle))
				   .Units(t_("N·m"), t_("sec")).NoPlot().SetDataSecondaryY();
			scatter.ZoomToFit(true, true);
			double mx = ceil(scatter.GetYMax());
			double mn = floor(scatter.GetYMin());
			scatter.SetXYMin(Null, mn);
			scatter.SetRange(Null, mx-mn);
			int sz = int((mx-mn)/5);
			scatter.SetMajorUnits(Null, sz > 0 ? 1 : 0.5);
			
			for (int i = 0; i < dangle.size(); ++i) {
				array.AddColumn(Format("%.1f %.1f", angle, dangle[i]), 60);
				int row = 0, col = array.GetColumnCount()-1;
				array.Set(row++, col, angle);
				array.Set(row++, col, dangle[i]);
				int numdec = 8;
				array.Set(row++, col, FDS(dgz[i], numdec));				
				array.Set(row++, col, FDS(dMoment[i], numdec));
				array.Set(row++, col, FDS(disp[i], numdec));
				array.Set(row++, col, FDS(vol[i], numdec));
				array.Set(row++, col, FDS(wett[i], numdec));
				array.Set(row++, col, FDS(wplane[i], numdec));
				array.Set(row++, col, FDS(draft[i], numdec));
				array.Set(row++, col, FDS(cb[i].x, numdec));
				array.Set(row++, col, FDS(cb[i].y, numdec));
				array.Set(row++, col, FDS(cb[i].z, numdec));
				array.Set(row++, col, FDS(cg[i].x, numdec));
				array.Set(row++, col, FDS(cg[i].y, numdec));
				array.Set(row++, col, FDS(cg[i].z, numdec));
			}
		}
		
		if (~opShowEnvelope) {
			mingz.Clear();
			mingz.SetCount(dangle.size(), 0);
			for (int i = 0; i < dangle.size(); ++i) {
				for (iangle = 0; iangle < datagz.size(); ++iangle) {
					double d = datagz[iangle][i];
					if (dangle[i] < 0) {
						if (d > 0) {
							mingz[i] = 0;
							break;
						} else if (mingz[i] == 0) 
							mingz[i] = d;
						else 
							mingz[i] = mingz[i] < d ? d : mingz[i];
					} else {
						if (d < 0) {
							mingz[i] = 0;
							break;
						} else if (mingz[i] == 0) 
							mingz[i] = d;
						else 
							mingz[i] = mingz[i] > d ? d : mingz[i];
					}
				}
			}
			scatter.AddSeries(dangle, mingz).NoMark().Legend(t_("Envelope"))
					   .Units(t_("m"), t_("sec")).Stroke(4, Black()).Dash(LINE_SOLID);
		}
		if (!errors.IsEmpty()) 
			BEM::PrintError(DeQtfLf("Errors found in mesh:\n" + errors));
	} catch(Exc e) {
		BEM::PrintError(DeQtfLf(e));
		Clear(true);
	}
}

void MainGZ::Clear(bool force) {
	if (!force) {
		MainBody &mm = GetDefinedParent<MainBody>(this);
		
		UVector<int> ids = ArrayModel_IdsBody(mm.listLoaded);
		if (Find(ids, idOpened) >= 0)
			return;
	}
	scatter.RemoveAllSeries();
	
	array.Reset();
	array.SetLineCy(EditField::GetStdHeight());
	array.HeaderObject().Absolute();
	array.MultiSelect();
	array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, array);};
	array.AddColumn(t_("Parameter"), 100);
	int row = 0;
	array.Set(row++, 0, t_("ZoY plane [º]"));
	array.Set(row++, 0, t_("Angle [º]"));
	array.Set(row++, 0, t_("GZ [m]"));
	array.Set(row++, 0, t_("Heeling lever [N·m]"));
	array.Set(row++, 0, t_("Displacement [kg]"));
	array.Set(row++, 0, t_("Sub. volume [m3]"));
	array.Set(row++, 0, t_("Wetted area [m²]"));
	array.Set(row++, 0, t_("Waterpl. area [m²]"));
	array.Set(row++, 0, t_("Draft [m]"));
	array.Set(row++, 0, t_("Cb_x [m]"));
	array.Set(row++, 0, t_("Cb_y [m]"));
	array.Set(row++, 0, t_("Cb_z [m]"));
	array.Set(row++, 0, t_("Cg_x [m]"));
	array.Set(row++, 0, t_("Cg_y [m]"));
	array.Set(row++, 0, t_("Cg_z [m]"));
}
		
void MainGZ::Jsonize(JsonIO &json) {
	if (json.IsLoading()) {
		edFrom <<= Null;	
		edTo <<= Null;
		edDelta <<= Null;
		edAngleFrom <<= Null;
		edAngleTo <<= Null;
		edAngleDelta <<= Null;
		edTolerance <<= Null;
	}
	json
		("edFrom", edFrom)
		("edTo", edTo)
		("edDelta", edDelta)
		("edAngleFrom", edAngleFrom)
		("edAngleTo", edAngleTo)
		("edAngleDelta", edAngleDelta)
		("opShowEnvelope", opShowEnvelope)
		("edTolerance", edTolerance)
	;
}
