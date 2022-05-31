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


void MainMesh::Init() {
	CtrlLayout(*this);
	
	OnOpt();
	
	CtrlLayout(menuOpen);
	menuOpen.file.WhenChange = THISBACK(OnLoad);
	menuOpen.file.BrowseRightWidth(40).UseOpenFolder(true).BrowseOpenFolderWidth(10);
	menuOpen.butLoad.Tip(t_("Loads mesh file")).WhenAction = [&] {menuOpen.file.DoGo();};

	ArrayModel_Init(listLoaded, true).MultiSelect();
	listLoaded.WhenSel = [&] {
		OnMenuConvertArraySel();
		OnMenuProcessArraySel();
		OnMenuAdvancedArraySel();
		LoadSelTab(Bem());
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
	
	menuOpen.butRemove.Tip(t_("Removes all loaded files")).Disable();	
	menuOpen.butRemove.WhenAction = THISBACK(OnRemove);
	menuOpen.butRemoveSelected.Tip(t_("Removes selected files")).Disable();	
	menuOpen.butRemoveSelected.WhenAction = THISBACK1(OnRemoveSelected, false);
	menuOpen.butJoin.Tip(t_("Joins selected meshes")).Disable();	
	menuOpen.butJoin.WhenAction = THISBACK(OnJoin);
	menuOpen.butSplit.Tip(t_("Splits mesh in parts (if parts are not joined together)")).Disable();	
	menuOpen.butSplit.WhenAction = THISBACK(OnSplit);
	menuOpen.opClean.Tip(t_("Cleans duplicated panels when loading (it may be slow!)"));
	
	CtrlLayout(menuConvert);
	menuConvert.file.WhenChange = THISBACK(OnConvertMesh);
	menuConvert.file.BrowseRightWidth(40).UseOpenFolder(true).BrowseOpenFolderWidth(10);
	menuConvert.file.SelLoad(false);
	menuConvert.butConvert.Tip(t_("Saves mesh file")).WhenAction = [&] {menuConvert.file.DoGo();};

	menuConvert.opt.WhenAction = [&] {OnOpt();};

	OnOpt();
	
	CtrlLayout(menuPlot);
	menuPlot.showMesh.Tip(t_("Shows all the loaded mesh")).WhenAction 	    	= [&] {LoadSelTab(Bem());};
	menuPlot.showNormals.Tip(t_("Shows panel normals as arrows")).WhenAction 	= [&] {LoadSelTab(Bem());};
	menuPlot.showWaterLevel.Tip(t_("Shows the cut of the mesh with the water line")).WhenAction  = [&] {LoadSelTab(Bem());};
	menuPlot.showSkewed.Tip(t_("Shows skewed panels")).WhenAction      			= [&] {LoadSelTab(Bem());};
	menuPlot.showFissure.Tip(t_("Shows fissures in the hull")).WhenAction     	= [&] {LoadSelTab(Bem());};
	menuPlot.showMultiPan.Tip(t_("Shows wrong panels")).WhenAction    			= [&] {LoadSelTab(Bem());};
	menuPlot.showAxis.Tip(t_("Shows system axis")).WhenAction  					= [&] {mainView.gl.Refresh();};
	menuPlot.showLimits.Tip(t_("Shows boundaris of the geometry")).WhenAction 	= [&] {mainView.gl.Refresh();};
	menuPlot.showCb.Tip(t_("Shows the centre of buoyancy")).WhenAction  		= [&] {mainView.gl.Refresh();};
	menuPlot.showCg.Tip(t_("Shows the centre of gravity")).WhenAction  			= [&] {mainView.gl.Refresh();};
	menuPlot.showCr.Tip(t_("Shows the centre of rotation")).WhenAction  		= [&] {mainView.gl.Refresh();};
	menuPlot.showSel.Tip(t_("Shows volume around selected object")).WhenAction  = [&] {mainView.gl.Refresh();};	
	menuPlot.showUnderwater.Tip(t_("Shows nderwater mesh")).WhenAction  		= [&] {mainView.gl.Refresh();};
	menuPlot.butXYZ.Tip(t_("Orients the camera as isometric")).WhenAction  		= [&] {mainView.gl.View(true, true, true);};
	menuPlot.butXoY.Tip(t_("Orients the camera through Z axis")).WhenAction  	= [&] {mainView.gl.View(true, true, false);};	
	menuPlot.butYoZ.Tip(t_("Orients the camera through X axis")).WhenAction  	= [&] {mainView.gl.View(false, true, true);};
	menuPlot.butXoZ.Tip(t_("Orients the camera through Y axis")).WhenAction  	= [&] {mainView.gl.View(true, false, true);};
	menuPlot.butFit.Tip(t_("Zooms the camera to fit the bodies")).WhenAction	= [&] {mainView.gl.ZoomToFit();};
	menuPlot.showMeshData.Tip(t_("Shows a list of panels and nodes")).WhenAction= [&] {splitterAll.SetButton(0);};
	menuPlot.showMeshData.Tip(t_("Controls for 3D playing")).WhenAction			= [&] {splitterVideo.SetButton(0);};
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
	
	menuProcess.x_g <<= 0;
	menuProcess.x_g.WhenEnter = THISBACK2(OnUpdate, NONE, true);
	menuProcess.y_g <<= 0;
	menuProcess.y_g.WhenEnter = THISBACK2(OnUpdate, NONE, true);
	menuProcess.z_g <<= 0;
	menuProcess.z_g.WhenEnter = THISBACK2(OnUpdate, NONE, true);
	menuProcess.x_0 <<= 0;
	menuProcess.x_0.WhenEnter = THISBACK2(OnUpdate, NONE, true);
	menuProcess.y_0 <<= 0;
	menuProcess.y_0.WhenEnter = THISBACK2(OnUpdate, NONE, true);
	menuProcess.z_0 <<= 0;
	menuProcess.z_0.WhenEnter = THISBACK2(OnUpdate, NONE, true);
	menuProcess.mass <<= 0;
	menuProcess.mass.WhenEnter = THISBACK2(OnUpdate, NONE, true);
	menuProcess.butUpdateCg  <<= THISBACK2(OnUpdate, NONE, true);
	menuProcess.butUpdateCg.Tip(t_("Sets the centre of gravity"));
	menuProcess.butCgtoC0.WhenAction = [&] {
		menuProcess.x_0 <<= ~menuProcess.x_g;
		menuProcess.y_0 <<= ~menuProcess.y_g;
		menuProcess.z_0 <<= ~menuProcess.z_g;
		OnUpdate(NONE, true);
	};
	menuProcess.butCgtoC0.Tip(t_("Sets the centre of rotation with the centre of gravity"));
	menuProcess.butC0toCg.WhenAction = [&] {
		menuProcess.x_g <<= ~menuProcess.x_0;
		menuProcess.y_g <<= ~menuProcess.y_0;
		menuProcess.z_g <<= ~menuProcess.z_0;
		OnUpdate(NONE, true);
	};
	menuProcess.butC0toCg.Tip(t_("Sets the centre of gravity with the centre of rotation"));
	menuProcess.butUpdateCrot  <<= THISBACK2(OnUpdate, NONE, true);
	menuProcess.butUpdateCrot.Tip(t_("Sets the centre of rotation"));	
	
	menuProcess.butUpdateMass  <<= THISBACK(OnUpdateMass);
	menuProcess.butUpdateMass.Tip(t_("Sets mass from inmersed volume"));

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
	
	CtrlLayout(menuMove);	
	
	menuMove.butReset <<= THISBACK(OnReset);
	menuMove.butReset.Tip(t_("Translates the mesh"));
	
	menuMove.butUpdateCrot  <<= THISBACK2(OnUpdate, NONE, false);
	menuMove.butUpdateCrot.Tip(t_("Sets the centre of rotation"));	
	
	
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
	
	//menuMove.pos_x <<= 0;
	//menuMove.pos_y <<= 0;
	//menuMove.pos_z <<= 0;
	//menuMove.ang_x <<= 0;
	//menuMove.ang_y <<= 0;
	//menuMove.ang_z <<= 0;
	
	menuMove.opZArchimede.WhenAction = [&] {menuMove.t_z.Enable(!menuMove.opZArchimede);};
	
	CtrlLayout(menuEdit);
	menuEdit.edit_x <<= 0;
	menuEdit.edit_y <<= 0;
	menuEdit.edit_z <<= 0;
	menuEdit.edit_size <<= 1;
	
	menuEdit.butPanel <<= THISBACK(OnAddPanel);
	
	menuEdit.revolutionList.AddColumn("H").Ctrls<EditString>().HeaderTab().SetMargin(-2);
	menuEdit.revolutionList.AddColumn("V").Ctrls<EditString>().HeaderTab().SetMargin(-2);
	menuEdit.revolutionList.AddColumn("");
	menuEdit.revolutionList.ColumnWidths("10 10 2");
	menuEdit.revolutionList.NoHeader().MultiSelect();
	menuEdit.revolutionList.SetLineCy(int(EditField::GetStdHeight()*2/3));
	menuEdit.revolutionList.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, menuEdit.revolutionList, true, true);};
	
	menuEdit.butRevolution <<= THISBACK(OnAddRevolution);

	menuEdit.polynomialList.AddColumn("H").Ctrls<EditString>().HeaderTab().SetMargin(-2);
	menuEdit.polynomialList.AddColumn("V").Ctrls<EditString>().HeaderTab().SetMargin(-2);
	menuEdit.polynomialList.AddColumn("");
	menuEdit.polynomialList.ColumnWidths("10 10 2");
	menuEdit.polynomialList.NoHeader().MultiSelect();
	menuEdit.polynomialList.SetLineCy(int(EditField::GetStdHeight()*2/3));
	menuEdit.polynomialList.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, menuEdit.polynomialList, true, true);};
	
	menuEdit.butPolynomial    <<= THISBACK(OnAddPolygonalPanel);
		
	menuTab.Add(menuOpen.SizePos(),    	t_("Load"));
	menuTab.Add(menuPlot.SizePos(),    	t_("Plot")).Disable();
	menuTab.Add(menuMove.SizePos(), 	t_("Move")).Disable();
	menuTab.Add(menuProcess.SizePos(), 	t_("Process")).Disable();	
	menuTab.Add(menuEdit.SizePos(), 	t_("Edit"));
	menuTab.Add(menuConvert.SizePos(), 	t_("Save as")).Disable();
		
	mainViewData.Init();
	splitterVideo.Vert(mainView.SizePos(), videoCtrl.SizePos());
	splitterVideo.SetPositions(8500, 9900).SetInitialPositionId(1).SetButtonNumber(1).SetButtonWidth(20);
	splitterAll.Horz(splitterVideo.SizePos(), mainViewData.SizePos());
	splitterAll.SetPositions(6000, 9900).SetInitialPositionId(1).SetButtonNumber(1).SetButtonWidth(20);
	splitterAll.Tip(t_("")).WhenAction = [&] {mainView.SetPaintSelect(splitterAll.GetPos() < 9900);};
	mainTab.Add(splitterAll.SizePos(), t_("View"));
	mainView.Init();
	
	String bitmapFolder = AppendFileNameX(GetDesktopFolder(), "BEMRosetta Mesh Images");
	int idBitmapFolder = 0;
	
	videoCtrl.Init([&](UVector<int> &ids)->int {
			ids = ArrayModel_IdsMesh(listLoaded);
			int num = ArrayCtrlSelectedGetCount(listLoaded);
			if (num > 1) {
				Exclamation(t_("Please select just one model"));
				return -1;
			}
			int id;
			if (num == 0 && listLoaded.GetCount() == 1)
				id = ArrayModel_IdMesh(listLoaded, 0);
			else {
			 	id = ArrayModel_IdMesh(listLoaded);
				if (id < 0) {
					Exclamation(t_("Please select a model to process"));
					return -1;
				}
			}
			return id;
		}, [&](int id, const UVector<int> &ids, const Point3D &pos, const Point3D &angle, const Point3D &c0, bool full, bool saveBitmap) {
			Mesh &data = Bem().surfs[id];			
			
			data.cg.TransRot(pos.x, pos.y, pos.z, ToRad(angle.x), ToRad(angle.y), ToRad(angle.z), c0.x, c0.y, c0.z);
			data.mesh.TransRot(pos.x, pos.y, pos.z, ToRad(angle.x), ToRad(angle.y), ToRad(angle.z), c0.x, c0.y, c0.z);
			
			menuProcess.x_g <<= data.cg.x;
			menuProcess.y_g <<= data.cg.y;
			menuProcess.z_g <<= data.cg.z;
			
			data.AfterLoad(Bem().rho, Bem().g, false, false);
			
			if (full)
				mainStiffness.Load(Bem().surfs, ids);
			mainView.CalcEnvelope();
			if (full)
				mainSummary.Report(Bem().surfs, id);
			
			mainView.gl.Refresh();
			if (saveBitmap) {
				RealizeDirectory(bitmapFolder);
				DeleteFileDeepWildcardsX(bitmapFolder);
				mainView.gl.SaveToFile(AppendFileNameX(bitmapFolder, Format("Image%4d", idBitmapFolder++)));
			}
			if (full)
				mainViewData.OnRefresh();
		});
	
	
	mainSummary.Init();
	mainTab.Add(mainSummary.SizePos(), t_("Summary"));

	mainStiffness.Init(Hydro::MAT_K);
	mainTab.Add(mainStiffness.SizePos(), t_("K Stiffness Matrix"));
	
	mainGZ.Init();
	mainTab.Add(mainGZ.SizePos(), t_("GZ"));
			
	mainTab.WhenSet = [&] {
		LOGTAB(mainTab);
		UVector<int> ids = ArrayModel_IdsMesh(listLoaded);
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
		TabCtrl::Item& tabMenuConvert = menuTab.GetItem(menuTab.Find(menuConvert));
		tabMenuConvert.Enable(convertProcess);
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
			tabMenuConvert.Text(t_("Save as"));
			tabMenuProcess.Text(t_("Process"));
			tabMenuMove.Text(t_("Move"));
		} else {
			tabMenuConvert.Text("");
			tabMenuProcess.Text("");
			tabMenuMove.Text("");
		}
	};
	mainTab.WhenSet();
	
	menuTab.WhenSet = [&] {
	LOGTAB(menuTab);
		if (menuTab.IsAt(menuConvert)) 
			listLoaded.WhenSel();
	};
	menuTab.WhenSet();
}

void MainMesh::OnMenuProcessArraySel() {
	int id = ArrayModel_IdMesh(listLoaded);
	if (id < 0)
		return;
	
	Mesh &data = Bem().surfs[id];
	menuProcess.x_g <<= data.cg.x;
	menuProcess.y_g <<= data.cg.y;
	menuProcess.z_g <<= data.cg.z;
	menuProcess.x_0 <<= data.c0.x;
	menuProcess.y_0 <<= data.c0.y;
	menuProcess.z_0 <<= data.c0.z;
	menuProcess.mass <<= data.mass;
}

void MainMesh::OnMenuMoveArraySel() {
	int id = ArrayModel_IdMesh(listLoaded);
	if (id < 0)
		return;
	
	Mesh &data = Bem().surfs[id];
	menuMove.x_0 <<= data.c0.x;
	menuMove.y_0 <<= data.c0.y;
	menuMove.z_0 <<= data.c0.z;	
}

void MainMesh::OnMenuAdvancedArraySel() {
	int id = ArrayModel_IdMesh(listLoaded);
	if (id < 0)
		return;
	
	//Mesh &data = Bem().surfs[id];
	
}

void MainMesh::OnArraySel() {
	OnMenuConvertArraySel();
	OnMenuProcessArraySel();
	OnMenuMoveArraySel();
	OnMenuAdvancedArraySel();
}
	
void MainMesh::OnMenuConvertArraySel() {
	int id = ArrayModel_IdMesh(listLoaded);
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

	menuConvert.symX <<= (ext == ".gdf" && Bem().surfs[id].IsSymmetricX());
	menuConvert.symY <<= Bem().surfs[id].IsSymmetricY();
	
	UpdateButtons();
}
	
void MainMesh::InitSerialize(bool ret) {
	if (!ret || IsNull(menuPlot.showMesh)) 
		menuPlot.showMesh = true;	
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
	if (!ret || IsNull(menuPlot.showSel)) 
		menuPlot.showSel = true;	
	if (!ret || IsNull(menuPlot.lineThickness)) 
		menuPlot.lineThickness <<= 1;
	if (!ret || IsNull(menuPlot.backColor)) 
		menuPlot.backColor <<= White();
				
	if (!ret || IsNull(menuConvert.opt)) 
		menuConvert.opt = 0;
	if (!ret || IsNull(menuConvert.optMeshType)) 
		menuConvert.optMeshType = 0;
}

void MainMesh::LoadSelTab(BEM &bem) {
	const UVector<int> &ids = ArrayModel_IdsMesh(listLoaded);
	if (mainTab.Get() == mainTab.Find(mainStiffness))
		mainStiffness.Load(bem.surfs, ids);
	else if (mainTab.Get() == mainTab.Find(mainView))
		mainView.gl.Refresh();
}

void MainMesh::OnOpt() {
	menuOpen.file.ClearTypes(); 

	const String meshFiles = ".gdf .dat .stl .pnl .msh .mesh";
	String meshFilesAst = clone(meshFiles);
	meshFilesAst.Replace(".", "*.");
	menuOpen.file.Type(Format("All supported mesh files (%s)", meshFiles), meshFilesAst);
	menuOpen.file.AllFilesType();
	String extView = ToLower(GetFileExt(menuOpen.file.GetData().ToString()));
	if (extView.IsEmpty())
		menuOpen.file.ActiveType(0);
	else if (meshFiles.Find(extView) >= 0)
		menuOpen.file.ActiveType(0);
	else
		menuOpen.file.ActiveType(1);
	
	menuConvert.file.ClearTypes();
	menuConvert.symX.Disable();
	switch (menuConvert.opt) {
	case 0:	menuConvert.file <<= ForceExtSafe(~menuConvert.file, ".gdf"); 	
			menuConvert.file.Type("Wamit .gdf file", "*.gdf");
			menuConvert.symX.Enable();
			break;
	case 1:	menuConvert.file <<= ForceExtSafe(~menuConvert.file, ".dat"); 	
			menuConvert.file.Type("Nemoh .dat file", "*.dat");
			break;
	case 2:	menuConvert.file <<= ForceExtSafe(~menuConvert.file, "."); 	
			menuConvert.file.Type("Nemoh pre mesh file", "*.");
			break;
	case 3:	menuConvert.file <<= ForceExtSafe(~menuConvert.file, ".pnl"); 	
			menuConvert.file.Type("HAMS .pnl file", "*.pnl");
			menuConvert.symX.Enable();
			break;
	case 4:	menuConvert.file <<= ForceExtSafe(~menuConvert.file, ".dat"); 	
			menuConvert.file.Type("Diodore .dat file", "*.dat");
			menuConvert.symX.Enable();
			break;			
	case 5:	menuConvert.file <<= ForceExtSafe(~menuConvert.file, ".stl"); 	
			menuConvert.file.Type("STL binary .stl file", "*.stl");
			break;
	case 6:	menuConvert.file <<= ForceExtSafe(~menuConvert.file, ".stl"); 	
			menuConvert.file.Type("STL text .stl file", "*.stl");
			break;
	default:menuConvert.file.Type("All converted files", "*.gdf *.dat *.pnl *.stl");
			menuConvert.symX.Enable();
			break;
	}
	String extConvmesh = ToLower(GetFileExt(menuConvert.file.GetData().ToString()));
	if (extConvmesh.IsEmpty())
		menuConvert.file.ActiveType(0);
	else if (String(".gdf").Find(extConvmesh) >= 0)
		menuConvert.file.ActiveType(0);
	else if (String(".dat").Find(extConvmesh) >= 0)
		menuConvert.file.ActiveType(1);
	else if (String(".pnl").Find(extConvmesh) >= 0)
		menuConvert.file.ActiveType(2);
	else if (String(".stl").Find(extConvmesh) >= 0)
		menuConvert.file.ActiveType(3);
}

void MainMesh::AfterAdd(String file) {
	int id = Bem().surfs.size() - 1;
	Mesh &surf = Bem().surfs[id];
		
	mainTab.Set(mainView);
	
	surf.Report(Bem().rho);
	
	AddRow(surf);

	mainView.CalcEnvelope();
	
	mainView.gl.Enable();
	mainView.gl.ZoomToFit();
	
	After();	
	OnArraySel();	
}

bool MainMesh::OnLoad() {
	GuiLock __;
	
	String file = ~menuOpen.file;
		
	try {
		Progress progress(t_("Loading mesh file..."), 100); 
		
		UVector<int> ids = ArrayModel_IdsMesh(listLoaded);
		for (int i = 0; i < ids.size(); ++i) {
			if (Bem().surfs[ids[i]].fileName == file) {
				if (!PromptYesNo(t_("Model is already loaded") + S("&") + t_("Do you wish to open it anyway?")))
					return false;
				break;
			}
		}
		mainView.gl.Disable();
		
		WaitCursor waitcursor;

		Bem().LoadMesh(file, [&](String str, int _pos) {progress.SetText(str); progress.SetPos(_pos); return progress.Canceled();}, ~menuOpen.opClean, false);
		
		int id = Bem().surfs.size()-1;
		Mesh &data = Bem().surfs[id];
		
		if (!data.mesh.IsEmpty() && ~menuMove.opZArchimede) {
			double dz = 0.1;
			Surface under;
			if (data.mesh.TranslateArchimede(data.mass, Bem().rho, dz, under)) {
				data.cg.Translate(0, 0, dz);
				videoCtrl.AddReg(Point3D(0, 0, dz));
				//menuMove.pos_z <<= data.mesh.GetPos().z;
			} else
				Exclamation(t_("Problem readjusting the Z value to comply with displacement"));
			
			ma().Status(Format(t_("Loaded '%s', and translated vertically %f m to comply with displacement"), file, dz));
		} else
			ma().Status(Format(t_("Loaded '%s'"), file));			
		
		AfterAdd(file);
		
		mainViewData.OnAddedModel(mainView);
		
	} catch (Exc e) {
		mainView.gl.Enable();
		Exclamation(DeQtfLf(e));
		return false;
	}
	return true;
}

bool MainMesh::OnConvertMesh() {
	GuiLock __;
	
	try {
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			Exclamation(t_("Please select just one model"));
			return false;
		}
		int id;
		if (num == 0 && listLoaded.GetCount() == 1)
			id = ArrayModel_IdMesh(listLoaded, 0);
		else {
		 	id = ArrayModel_IdMesh(listLoaded);
			if (id < 0) {
				Exclamation(t_("Please select a model to process"));
				return false;
			}
		}
		
		Mesh::MESH_FMT type;	
		switch (menuConvert.opt) {
		case 0:	type = Mesh::WAMIT_GDF;	break;
		case 1:	type = Mesh::NEMOH_DAT;	break;
		case 2:	type = Mesh::NEMOH_PRE;	break;
		case 3:	type = Mesh::HAMS_PNL;	break;
		case 4:	type = Mesh::DIODORE_DAT;break;
		case 5:	type = Mesh::STL_BIN;	break;
		case 6:	type = Mesh::STL_TXT;	break;
		case 7:	type = Mesh::UNKNOWN;	break;
		default: throw Exc(t_("Unknown type in OnConvert()"));
		}
		
		WaitCursor waitcursor;
		Progress progress(t_("Saving mesh file..."), 100); 
		
		Bem().surfs[id].SaveAs(~menuConvert.file, type, Bem().g, 
							   static_cast<Mesh::MESH_TYPE>(int(~menuConvert.optMeshType)),
							   ~menuConvert.symX, ~menuConvert.symY);	
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
		return false;
	}
	return true;
}

void MainMesh::OnReset() {
	GuiLock __;
	
	try {
		UVector<int> ids = ArrayModel_IdsMesh(listLoaded);
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			Exclamation(t_("Please select just one model"));
			return;
		}
		int id;
		if (num == 0 && listLoaded.GetCount() == 1)
			id = ArrayModel_IdMesh(listLoaded, 0);
		else {
		 	id = ArrayModel_IdMesh(listLoaded);
			if (id < 0) {
				Exclamation(t_("Please select a model to process"));
				return;
			}
		}
				
		Mesh &data = Bem().surfs[id];
		data.Reset(Bem().rho, Bem().g);

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
		menuProcess.x_g <<= data.cg.x; 
		menuProcess.y_g <<= data.cg.y;
		menuProcess.z_g <<= data.cg.z;
		
		ma().Status(t_("Model oriented on the initial layout"));
		
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}
}

void MainMesh::OnUpdateMass() {
	GuiLock __;
	
	try {
		UVector<int> ids = ArrayModel_IdsMesh(listLoaded);
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			Exclamation(t_("Please select just one model"));
			return;
		}
		int id;
		if (num == 0 && listLoaded.GetCount() == 1)
			id = ArrayModel_IdMesh(listLoaded, 0);
		else {
		 	id = ArrayModel_IdMesh(listLoaded);
			if (id < 0) {
				Exclamation(t_("Please select a model to process"));
				return;
			}
		}
		
		Mesh &data = Bem().surfs[id];
		
		double mass = data.under.volume*Bem().rho;
		menuProcess.mass <<= mass;

		OnUpdate(NONE, true);
		
		ma().Status(Format(t_("Mass updated to %f kg"), mass));
		
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}	
}

void MainMesh::OnUpdate(Action action, bool fromMenuProcess) {
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
		
		UVector<int> ids = ArrayModel_IdsMesh(listLoaded);
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			Exclamation(t_("Please select just one model"));
			return;
		}
		int id;
		if (num == 0 && listLoaded.GetCount() == 1)
			id = ArrayModel_IdMesh(listLoaded, 0);
		else {
		 	id = ArrayModel_IdMesh(listLoaded);
			if (id < 0) {
				Exclamation(t_("Please select a model to process"));
				return;
			}
		}
				
		Mesh &data = Bem().surfs[id];

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
			Exclamation(t_("Please fill CG data"));
			return;
		}
		if (action == NONE && (IsNull(mass) || IsNull(x_0) || IsNull(y_0) || IsNull(z_0))) {
			Exclamation(t_("Please fill centre of rotation data"));
			return;
		}
				
		if (action == MOVE && (IsNull(t_x) || IsNull(t_y) || IsNull(t_z))) {
			Exclamation(t_("Please fill translation data"));
			return;
		}
		if (action == ROTATE && (IsNull(a_x) || IsNull(a_y) || IsNull(a_z))) {
			Exclamation(t_("Please fill rotation data"));
			return;
		}
		
		WaitCursor wait;

		if (action == MOVE) {
			data.cg.Translate(t_x, t_y, t_z);
			data.mesh.Translate(t_x, t_y, t_z);
			videoCtrl.AddReg(Point3D(t_x, t_y, t_z));
			
			ma().Status(Format(t_("Model moved %f, %f, %f"), t_x, t_y, t_z));
			
		} else if (action == ROTATE) {
			data.cg.Rotate(ToRad(a_x), ToRad(a_y), ToRad(a_z), x_0, y_0, z_0);
			data.mesh.Rotate(ToRad(a_x), ToRad(a_y), ToRad(a_z), x_0, y_0, z_0);
			videoCtrl.AddReg(Point3D(a_x, a_y, a_z), Point3D(x_0, y_0, z_0));
			
			if (~menuMove.opZArchimede) {
				double dz = 0;
				Surface under;
				data.mesh.TranslateArchimede(data.mass, Bem().rho, dz, under);
				if (!IsNull(dz)) {
					data.cg.Translate(0, 0, dz);
					videoCtrl.AddReg(Point3D(0, 0, dz));
				} else
					Exclamation(t_("Problem readjusting the Z value to comply with displacement"));
				
				ma().Status(Format(t_("Model rotated %f, %f, %f deg. around %f, %f, %f, and translated vertically %f m to comply with displacement"), a_x, a_y, a_z, x_0, y_0, z_0, dz));
			} else
				ma().Status(Format(t_("Model rotated %f, %f, %f around %f, %f, %f"), a_x, a_y, a_z, x_0, y_0, z_0));
		} else if (action == NONE) {
			data.mass = mass;
			data.cg.Set(x_g, y_g, z_g);
			data.c0.Set(x_0, y_0, z_0);
		}
		
		menuProcess.x_g <<= data.cg.x;
		menuProcess.y_g <<= data.cg.y;
		menuProcess.z_g <<= data.cg.z;
		
		//menuMove.pos_x <<= data.mesh.GetPos().x;
		//menuMove.pos_y <<= data.mesh.GetPos().y;
		//menuMove.pos_z <<= data.mesh.GetPos().z;
		//menuMove.ang_x <<= data.mesh.GetAngle().x;
		//menuMove.ang_y <<= data.mesh.GetAngle().y;
		//menuMove.ang_z <<= data.mesh.GetAngle().z;
		
		data.AfterLoad(Bem().rho, Bem().g, action == NONE, false);
		
	 	mainStiffness.Load(Bem().surfs, ids);
		mainView.CalcEnvelope();
		mainSummary.Report(Bem().surfs, id);
		
		mainView.gl.Refresh();
		mainViewData.OnRefresh();
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}
}

void MainMesh::OnArchimede() {
	GuiLock __;
	
	try {
		UVector<int> ids = ArrayModel_IdsMesh(listLoaded);
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			Exclamation(t_("Please select just one model"));
			return;
		}
		int id;
		if (num == 0 && listLoaded.GetCount() == 1)
			id = ArrayModel_IdMesh(listLoaded, 0);
		else {
		 	id = ArrayModel_IdMesh(listLoaded);
			if (id < 0) {
				Exclamation(t_("Please select a model to process"));
				return;
			}
		}
		
		Mesh &data = Bem().surfs[id];
		
		double dz = 0.5, droll = 0.5, dpitch = 0.5;
		Surface under;
		if (!data.mesh.Archimede(data.mass, data.cg, data.c0, Bem().rho, Bem().g, dz, droll, dpitch, under))
			Exclamation(t_("Problem readjusting the Z, roll and pitch values to comply with buoyancy"));
		
		menuProcess.x_g <<= data.cg.x;
		menuProcess.y_g <<= data.cg.y;
		menuProcess.z_g <<= data.cg.z;
				
		//menuMove.pos_x <<= data.mesh.GetPos().x;
		//menuMove.pos_y <<= data.mesh.GetPos().y;
		//menuMove.pos_z <<= data.mesh.GetPos().z;
		//const Point3D &ang = data.mesh.GetAngle();
		//menuMove.ang_x <<= data.mesh.GetAngle().x;
		//menuMove.ang_y <<= data.mesh.GetAngle().y;
		//menuMove.ang_z <<= data.mesh.GetAngle().z;
		
		data.AfterLoad(Bem().rho, Bem().g, false, false);
		
	 	mainStiffness.Load(Bem().surfs, ids);
		mainView.CalcEnvelope();
		mainSummary.Report(Bem().surfs, id);
		
		mainView.gl.Refresh();
		mainViewData.OnRefresh();
		
		ma().Status(Format(t_("Model rotated %f, %f, 0 deg and translated vertically %f m to comply with displacement"), droll, dpitch, dz));
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}			
}			
	

void MainMesh::OnAddPanel() {
	GuiLock __;
	
	if (IsNull(menuEdit.edit_x) || IsNull(menuEdit.edit_y) || IsNull(menuEdit.edit_z)) {
		Exclamation(t_("Panel position has to be set"));
		return;
	}
	if (IsNull(menuEdit.edit_size) || double(~menuEdit.edit_size) <= 0) {
		Exclamation(t_("Wrong mesh size"));
		return;
	}
	if (IsNull(menuEdit.panWidthX) || IsNull(menuEdit.panWidthY)) {
		Exclamation(t_("Panel width and height has to be set"));
		return;
	}
	
	WaitCursor waitcursor;
	mainView.gl.Disable();
	try {
		Bem().AddFlatPanel(~menuEdit.edit_x, ~menuEdit.edit_y, ~menuEdit.edit_z, ~menuEdit.edit_size, 
							 ~menuEdit.panWidthX, ~menuEdit.panWidthY);
		
		Mesh &surf = Bem().surfs[Bem().surfs.size()-1];
		surf.name = t_("Panel");
		surf.fileName =  "";
		
		surf.AfterLoad(Bem().rho, Bem().g, false, true);
		
		surf.Report(Bem().rho);
		AddRow(surf);
		After();
		mainViewData.OnAddedModel(mainView);
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}
	mainView.gl.Enable();
}

void MainMesh::OnAddRevolution() {
	GuiLock __;
	
	UVector<Pointf> vals;
	for (int r = 0; r < menuEdit.revolutionList.GetCount(); ++r) {
		Pointf &val = vals.Add();
		val.x = ScanDouble(AsString(menuEdit.revolutionList.Get(r, 0)));
		if (IsNull(val.x)) {
			Exclamation(Format(t_("Incorrect data in row %d, col %d"), r, 0));
			return;
		}
		val.y = ScanDouble(AsString(menuEdit.revolutionList.Get(r, 1)));
		if (IsNull(val.x)) {
			Exclamation(Format(t_("Incorrect data in row %d, col %d"), r, 1));
			return;
		}
	}
	if (vals.size() < 2) {
		Exclamation(t_("Unsufficient value number in list"));
		return;
	}
	
	WaitCursor waitcursor;
	mainView.gl.Disable();
	try {
		Bem().AddRevolution(~menuEdit.edit_x, ~menuEdit.edit_y, ~menuEdit.edit_z, ~menuEdit.edit_size, vals);
		
		Mesh &surf = Bem().surfs[Bem().surfs.size()-1];
		surf.name = t_("Revolution");
		surf.fileName =  "";
		
		surf.AfterLoad(Bem().rho, Bem().g, false, true);
		
		surf.Report(Bem().rho);
		AddRow(surf);
		After();
		mainViewData.OnAddedModel(mainView);
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}
	mainView.gl.Enable();
}

void MainMesh::OnAddPolygonalPanel() {
	GuiLock __;
	
	UVector<Pointf> vals;
	for (int r = 0; r < menuEdit.polynomialList.GetCount(); ++r) {
		Pointf &val = vals.Add();
		val.x = ScanDouble(AsString(menuEdit.polynomialList.Get(r, 0)));
		if (IsNull(val.x)) {
			Exclamation(Format(t_("Incorrect data in row %d, col %d"), r, 0));
			return;
		}
		val.y = ScanDouble(AsString(menuEdit.polynomialList.Get(r, 1)));
		if (IsNull(val.x)) {
			Exclamation(Format(t_("Incorrect data in row %d, col %d"), r, 1));
			return;
		}
	}
	if (vals.size() < 3) {
		Exclamation(t_("Unsufficient value number in list"));
		return;
	}
	
	WaitCursor waitcursor;
	mainView.gl.Disable();
	try {
		Bem().AddPolygonalPanel(~menuEdit.edit_x, ~menuEdit.edit_y, ~menuEdit.edit_z, ~menuEdit.edit_size, vals);
		
		Mesh &surf = Bem().surfs[Bem().surfs.size()-1];
		surf.name = t_("Polynomial");
		surf.fileName =  "";
		
		surf.AfterLoad(Bem().rho, Bem().g, false, true);
		
		surf.Report(Bem().rho);
		AddRow(surf);
		After();
		mainViewData.OnAddedModel(mainView);
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}
	mainView.gl.Enable();
}

void MainMesh::OnAddWaterSurface(char c) {
	GuiLock __;

	try {
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			Exclamation(t_("Please select just one model"));
			return;
		}
		int id;
		if (num == 0 && listLoaded.GetCount() == 1)
			id = ArrayModel_IdMesh(listLoaded, 0);
		else {
		 	id = ArrayModel_IdMesh(listLoaded);
			if (id < 0) {
				Exclamation(t_("Please select a model to process"));
				return;
			}
		}
			
		WaitCursor waitcursor;
		mainView.gl.Disable();
	
		Bem().AddWaterSurface(id, c);
		
		Mesh &surf = Bem().surfs[Bem().surfs.size()-1];
		AddRow(surf);
		After();
		mainViewData.OnAddedModel(mainView);
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}
	mainView.gl.Enable();
}

void MainMesh::OnHealing(bool basic) {
	GuiLock __;
	
	try {
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			Exclamation(t_("Please select just one model"));
			return;
		}
		int id;
		if (num == 0 && listLoaded.GetCount() == 1)
			id = ArrayModel_IdMesh(listLoaded, 0);
		else {
		 	id = ArrayModel_IdMesh(listLoaded);
			if (id < 0) {
				Exclamation(t_("Please select a model to process"));
				return;
			}
		}
		
		WaitCursor waitcursor;
		Progress progress(t_("Healing mesh file..."), 100); 
		mainView.gl.Disable();
		
		Bem().HealingMesh(id, basic, [&](String str, int _pos) {progress.SetText(str); progress.SetPos(_pos); return progress.Canceled();});
		
		UVector<int> ids = ArrayModel_IdsMesh(listLoaded);
	 	mainStiffness.Load(Bem().surfs, ids);
		mainView.CalcEnvelope();
		mainSummary.Report(Bem().surfs, id);
		
		mainView.gl.Refresh();
		mainViewData.OnRefresh();
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}	
	mainView.gl.Enable();
}

void MainMesh::OnOrientSurface() {
	GuiLock __;
	
	try {
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			Exclamation(t_("Please select just one model"));
			return;
		}
		int id;
		if (num == 0 && listLoaded.GetCount() == 1)
			id = ArrayModel_IdMesh(listLoaded, 0);
		else {
		 	id = ArrayModel_IdMesh(listLoaded);
			if (id < 0) {
				Exclamation(t_("Please select a model to process"));
				return;
			}
		}
		
		WaitCursor waitcursor;
		Progress progress(t_("Orienting mesh surface..."), 100); 
		mainView.gl.Disable();
		
		Bem().OrientSurface(id, [&](String str, int _pos) {progress.SetText(str); progress.SetPos(_pos); return progress.Canceled();});
		
		Bem().surfs[id].AfterLoad(Bem().rho, Bem().g, false, false);
		
		UVector<int> ids = ArrayModel_IdsMesh(listLoaded);
	 	mainStiffness.Load(Bem().surfs, ids);
		mainView.CalcEnvelope();
		mainSummary.Report(Bem().surfs, id);
		
		mainView.gl.Refresh();
		mainViewData.OnRefresh();
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}	
	mainView.gl.Enable();
}

void MainMesh::OnImage(int axis) {
	GuiLock __;
	
	String saxis = (axis == 0) ? "X" : ((axis == 1) ? "Y" : "Z");

	try {
		UVector<int> ids = ArrayModel_IdsMesh(listLoaded);
		int id = ArrayModel_IdMesh(listLoaded);
		if (id < 0) {
			Exclamation(t_("Please select a model to process"));
			return;
		}
		
		WaitCursor waitcursor;
		
		Mesh &data = Bem().surfs[id];

		data.mass = ~menuProcess.mass;
		if (axis == 0)
			data.cg.x = -data.cg.x;
		else if (axis == 1)
			data.cg.y = -data.cg.y;
		else
			data.cg.z = -data.cg.z;
		
		data.mesh.Image(axis);
	
		data.AfterLoad(Bem().rho, Bem().g, false, false);
		
	 	mainStiffness.Load(Bem().surfs, ids);
		mainView.CalcEnvelope();
		mainSummary.Report(Bem().surfs, id);
		
		mainView.gl.Refresh();
		mainViewData.OnRefresh();
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}	
}

void MainMesh::OnRemove() {
	WaitCursor waitcursor;

	OnRemoveSelected(true);
}

void MainMesh::OnRemoveSelected(bool all) {	
	bool selected = false;
	
	UVector<int> sel = ArrayCtrlSelectedGet(listLoaded);
	
	for (int r = listLoaded.GetCount()-1; r >= 0; --r) {
		if (all || Find(sel, r) >= 0) {
			int id = ArrayModel_IdMesh(listLoaded, r);
			Bem().RemoveMesh(id);
			listLoaded.Remove(r);
			selected = true;
		}
	}	// Only one available => directly selected
	if (!selected && listLoaded.GetCount() == 1) {
		int id = ArrayModel_IdMesh(listLoaded, 0);
		Bem().RemoveMesh(id);
		listLoaded.Remove(0);
		selected = true;		
	}	
	if (!selected) {
		Exclamation(t_("No model selected"));
		return;
	}

	UVector<int> ids = ArrayModel_IdsMesh(listLoaded);
	mainStiffness.Load(Bem().surfs, ids);
	mainViewData.ReLoad(mainView);
	
	mainGZ.Clear(false);
	
	After();
}

void MainMesh::OnJoin() {
	GuiLock __;
	
	try {	
		bool selected = false;
		int idDest = Null;
		for (int r = 0; r < listLoaded.GetCount(); ++r) {
			if (listLoaded.IsSelected(r)) {
				if (IsNull(idDest))
					idDest = ArrayModel_IdMesh(listLoaded, r);
				else
					idDest = min(idDest, ArrayModel_IdMesh(listLoaded, r));
			}
		}
		if (IsNull(idDest)) {
			Exclamation(t_("No model joined"));
			return;
		}
		
		WaitCursor waitcursor;
		
		for (int r = listLoaded.GetCount()-1; r >= 0; --r) {
			if (listLoaded.IsSelected(r)) {
				int id = ArrayModel_IdMesh(listLoaded, r);
				if (idDest != id) {
					Bem().JoinMesh(idDest, id);
					RemoveRow(r);
					selected = true;
				}
			}
		}	
	
		UVector<int> ids = ArrayModel_IdsMesh(listLoaded);
		mainStiffness.Load(Bem().surfs, ids);
		mainViewData.ReLoad(mainView);
		
		After();
	} catch (Exc e) {
		mainView.gl.Enable();
		Exclamation(DeQtfLf(e));
	}
}

void MainMesh::AddRow(const Mesh &surf) {
	ArrayModel_Add(listLoaded, surf.GetCodeMeshStr(), surf.name, surf.fileName, surf.GetId(),
					optionsPlot, [&] {mainView.gl.Refresh();});
}
		
void MainMesh::RemoveRow(int row) {
	listLoaded.Remove(row);
}	
	
void MainMesh::OnSplit() {
	GuiLock __;
	
	String file = ~menuOpen.file;
		
	try {
		Progress progress(t_("Splitting mesh file..."), 100); 
		
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num == 0 && listLoaded.GetCount() != 1) {
			Exclamation(t_("No model selected"));
			return;
		}
		if (num > 1) {
			Exclamation(t_("Please select just one model"));
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
			Exclamation(t_("The mesh is monolithic so it cannot be automatically split"));
			return;
		}
		int id = ArrayModel_IdMesh(listLoaded, row);
		String fileName = Bem().surfs[id].fileName;
		String name = Bem().surfs[id].name;
		idsmesh = Bem().SplitMesh(id, [&](String str, int _pos) {
			progress.SetText(str); 
			progress.SetPos(_pos); 
			progress.Refresh();
			return true;
		});	
		
		RemoveRow(row);
		
		for (int i = 0; i < idsmesh.size(); ++i) {
			int id = idsmesh[i];
			Mesh &surf = Bem().surfs[id];
			
			mainTab.Set(mainView);
			
			surf.Report(Bem().rho);
	
			surf.fileName = fileName;
			surf.name = name + S("-") + FormatInt(i+1);		
			AddRow(surf);
		}
		
		mainViewData.ReLoad(mainView);
				
		After();
	} catch (Exc e) {
		mainView.gl.Enable();
		Exclamation(DeQtfLf(e));
	}
}

void MainMesh::UpdateButtons() {
	int numrow = listLoaded.GetCount();
	int numsel = ArrayCtrlSelectedGetCount(listLoaded);
	menuOpen.butRemove.Enable(numrow > 0);
	menuOpen.butRemoveSelected.Enable(numsel > 0);
	menuOpen.butJoin.Enable(numsel > 1);
	menuOpen.butSplit.Enable(numsel == 1 || numrow == 1);
	menuConvert.butConvert.Enable(numsel == 1 || numrow == 1);
	
	menuProcess.butUpdateCg.Enable(numsel == 1 || numrow == 1);
	menuMove.butUpdatePos.Enable(numsel == 1 || numrow == 1);
	menuMove.butUpdateAng.Enable(numsel == 1 || numrow == 1);
	menuProcess.butImageX.Enable(numsel == 1 || numrow == 1);
	menuProcess.butImageY.Enable(numsel == 1 || numrow == 1);
	menuProcess.butImageZ.Enable(numsel == 1 || numrow == 1);
	menuProcess.butSimplify.Enable(numsel == 1 || numrow == 1);
	menuProcess.butFullHealing.Enable(numsel == 1 || numrow == 1);
	menuProcess.butOrientSurface.Enable(numsel == 1 || numrow == 1);
	menuProcess.butWaterFill.Enable(numsel == 1 || numrow == 1);
	menuProcess.butWaterPlane.Enable(numsel == 1 || numrow == 1);
	menuProcess.butHull.Enable(numsel == 1 || numrow == 1);
	menuMove.butReset.Enable(numsel == 1 || numrow == 1);
}

void MainMesh::After() {
	UpdateButtons();

	mainView.CalcEnvelope();	

	mainSummary.Clear();
	for (int row = 0; row < listLoaded.GetCount(); ++row) {
		int id = ArrayModel_IdMesh(listLoaded, row);
		mainSummary.Report(Bem().surfs, id);
	}		

	mainTab.WhenSet();
	
	mainView.gl.Enable();
	mainView.gl.Refresh();
}
	
void MainMesh::Jsonize(JsonIO &json) {
	if (json.IsLoading()) {
		menuPlot.showMesh = Null;	
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
		menuPlot.showSel = Null;
		menuConvert.opt = Null;
		menuConvert.optMeshType = Null;
		menuMove.opZArchimede = Null;
	}
	json
		("menuOpen_file", menuOpen.file)
		("menuConvert_file", menuConvert.file)	
		("menuConvert_opt", menuConvert.opt)
		("menuConvert_optMesh", menuConvert.optMeshType)
		("menuPlot_showMesh", menuPlot.showMesh)		
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
		("menuPlot_showUnderwater", menuPlot.showUnderwater)
		("menuPlot_showWaterLevel", menuPlot.showWaterLevel)
		("menuPlot_backColor", menuPlot.backColor)
		("menuPlot_lineThickness", menuPlot.lineThickness)	
		("mainStiffness", mainStiffness)
		("menuMove_opZArchimede", menuMove.opZArchimede)
		("mainGZ", mainGZ)
		("videoCtrl", videoCtrl)
	;
}

void MainSummaryMesh::Report(const UArray<Mesh> &surfs, int id) {
	const Mesh &data = surfs[id];
	String name = data.name;
	
	if (array.GetColumnCount() == 0)
		array.AddColumn("Param");
	if (id >= array.GetColumnCount()-1)
		array.AddColumn(Format("#%d %s", id+1, name));
	int row = 0;
	int col = id + 1;

	int idColor;
	Upp::Color backColorMesh;
	idColor = data.mesh.VolumeMatch(Bem().volWarning/100., Bem().volError/100.);
	if (idColor == -2 || data.mesh.surface < 0)
		backColorMesh = Upp::Color(255, 165, 158);	// Light red
	else if (idColor == -1)
		backColorMesh = Upp::Color(255, 255, 150);	// Light yellow
	
	Upp::Color backColorUnder;
	idColor = data.under.VolumeMatch(Bem().volWarning/100., Bem().volError/100.);
	if (idColor == -2)
		backColorUnder = Upp::Color(255, 165, 158);
	else if (idColor == -1)
		backColorUnder = Upp::Color(255, 255, 150);
									
	bool healing = data.mesh.healing;
	
	array.Set(row, 0, t_("File"));				array.Set(row++, col, data.fileName);
	array.Set(row, 0, t_("Name"));				array.Set(row++, col, name + (healing ? (S(" ") + t_("(healed)")) : ""));
	array.Set(row, 0, t_("Format"));			array.Set(row++, col, data.GetCodeMeshStr());	
	
	array.Set(row, 0, t_("# Panels"));			array.Set(row++, col, data.mesh.panels.size());
	array.Set(row, 0, t_("# Nodes"));			array.Set(row++, col, data.mesh.nodes.size());

	array.Set(row, 0, t_("Surface [m2]"));		array.Set(row++, col, FDS(data.mesh.surface, 8, false));
	array.Set(row, 0, t_("Volume [m3] Vavg (Vx,Vy,Vz)"));		  array.Set(row++, col, AttrText(Format(t_("%s (%s, %s, %s)"), 
														FDS(data.mesh.volume,  10, false),
														FDS(data.mesh.volumex, 10, false),
														FDS(data.mesh.volumey, 10, false),
														FDS(data.mesh.volumez, 10, false))).Paper(backColorMesh));
	
	array.Set(row, 0, t_("Wetted surface [m2]"));array.Set(row++, col, FDS(data.under.surface, 10, false));
	array.Set(row, 0, t_("Immersed volume [m3] Vavg (Vx,Vy,Vz)")); array.Set(row++, col, AttrText(Format(t_("%s (%s, %s, %s)"), 
														FDS(data.under.volume,  10, false),
														FDS(data.under.volumex, 10, false),
														FDS(data.under.volumey, 10, false),
														FDS(data.under.volumez, 10, false))).Paper(backColorUnder));
	array.Set(row, 0, t_("Displacement [kg]")); array.Set(row++, col, FDS(data.under.volume*Bem().rho, 10, false));
	array.Set(row, 0, t_("Cg [m]"));			array.Set(row++, col, Format(t_("%s, %s, %s"),
														FDS(data.cg.x, 10, false),			
														FDS(data.cg.y, 10, false),
														FDS(data.cg.z, 10, false)));
	array.Set(row, 0, t_("Cb [m]"));
	if (!IsNull(data.cb))	
		array.Set(row++, col, Format(t_("%s, %s, %s"),  FDS(data.cb.x, 10, false),			
														FDS(data.cb.y, 10, false),
														FDS(data.cb.z, 10, false)));
	else 
		array.Set(row++, col, "-");
	array.Set(row, 0, t_("C0 [m]"));			array.Set(row++, col, Format(t_("%s, %s, %s"),
														FDS(data.c0.x, 10, false),			
														FDS(data.c0.y, 10, false),
														FDS(data.c0.z, 10, false)));

	array.Set(row, 0, t_("GMroll [m]"));
	double gmpitch = data.GMpitch(Bem().rho, Bem().g); 	
	if (IsNum(gmpitch))
 		array.Set(row++, col, FDS(gmpitch, 5, false));
	else
		array.Set(row++, col, "-");
	array.Set(row, 0, t_("GMpitch [m]")); 		
	double gmroll = data.GMroll(Bem().rho, Bem().g); 	
	if (IsNum(gmroll))
		array.Set(row++, col, FDS(gmroll, 5, false));
	else
		array.Set(row++, col, "-");
	array.Set(row, 0, t_("GZ [m]"));
	if (!IsNull(data.cb)) {
		Direction3D cgcb = data.cg - data.cb;
		double gz = sqrt(sqr(cgcb.x) + sqr(cgcb.y));
		array.Set(row++, col, AttrText(FormatDouble(gz, 4)));	
	} else
		array.Set(row++, col, "-");
												
	array.Set(row, 0, t_("Surface projection Z-axis (Waterplane Area) [m2]"));	
												array.Set(row++, col, Format(t_("%s - %s = %s"),
														FDS(-data.zProjectionPos, 10, false),
														FDS(data.zProjectionNeg,  10, false),
														FDS(data.zProjectionPos+data.zProjectionNeg, 10, false)));
	
	array.Set(row, 0, t_("Waterplane geometric centre [m]"));
	if (!IsNull(data.cgZ0surface)) 
		array.Set(row++, col, Format(t_("%s, %s"), FDS(data.cgZ0surface.x, 10, false),			
												   FDS(data.cgZ0surface.y, 10, false)));
	else
		array.Set(row++, col, "-");
	array.Set(row, 0, t_("Surface projection X-axis [m2]"));	
												array.Set(row++, col, Format(t_("%s - %s = %s"),
														FDS(-data.xProjectionPos, 10, false),
														FDS(data.xProjectionNeg,  10, false),
														FDS(data.xProjectionPos+data.xProjectionNeg, 10, false)));
	
	array.Set(row, 0, t_("Surface projection Y-axis [m2]"));	
												array.Set(row++, col, Format(t_("%s - %s = %s"),
														FDS(-data.yProjectionPos, 10, false),
														FDS(data.yProjectionNeg,  10, false),
														FDS(data.yProjectionPos+data.yProjectionNeg, 10, false)));
	
	array.Set(row, 0, t_("Dimensions [m]"));	array.Set(row++, col, Format(t_("From (%s, %s, %s) to (%s, %s, %s)"),
														FDS(data.mesh.env.minX, 10, false),
														FDS(data.mesh.env.minY, 10, false),
														FDS(data.mesh.env.minZ, 10, false),
														FDS(data.mesh.env.maxX, 10, false),
														FDS(data.mesh.env.maxY, 10, false),
														FDS(data.mesh.env.maxZ, 10, false)));

	//Force6D f = data.under.GetHydrostaticForce(data.c0, Bem().rho, Bem().g);	
	Force6D fcb;
	if (!IsNull(data.cb))
		fcb = data.under.GetHydrostaticForceCB(data.c0, data.cb, Bem().rho, Bem().g);	
	
	array.Set(row, 0, t_("Hydrostatic forces [N]"));
	
	if (!IsNull(data.cb))
		array.Set(row++, col, AttrText(Format(t_("%s, %s, %s"),
														FDS(fcb[0], 10, false),
														FDS(fcb[1], 10, false),
														FDS(fcb[2], 10, false))).Paper(backColorUnder));
	else
		array.Set(row++, col, "-");

	array.Set(row, 0, t_("Hydrostatic moments [N·m]"));
	if (!IsNull(data.cb))
		array.Set(row++, col, AttrText(Format(t_("%s, %s, %s"),
														FDS(fcb[3], 10, false),
														FDS(fcb[4], 10, false),
														FDS(fcb[5], 10, false))));							
	else
		array.Set(row++, col, "-");
		
	array.Set(row, 0, t_("Mass [kg]"));			array.Set(row++, col, FormatF(data.mass, 1));
	
	Force6D fcg = Surface::GetMassForce(data.c0, data.cg, data.mass, Bem().g);	
	
	array.Set(row, 0, t_("Mass moments [N·m]"));array.Set(row++, col, Format(t_("%s, %s, %s"),
														FDS(fcg[3], 10, false),
														FDS(fcg[4], 10, false),
														FDS(fcg[5], 10, false)));

	array.Set(row, 0, t_("Mass+Hydrostatics forces [N]"));
	if (!IsNull(data.cb))
		array.Set(row++, col, AttrText(Format(t_("%s, %s, %s"),
														"0",
														"0",
														FDS(fcg[2]+fcb[2], 10, false))).Paper(backColorUnder));
	else
		array.Set(row++, col, "-");
	
	array.Set(row, 0, t_("Mass+Hydrostatics moments [N·m]"));
	if (!IsNull(data.cb))
		array.Set(row++, col, AttrText(Format(t_("%s, %s, %s"),
														FDS(fcg[3]+fcb[3], 10, false),
														FDS(fcg[4]+fcb[4], 10, false),
														"0")));		
	else
		array.Set(row++, col, "-");
														
	array.Set(row++, 0, t_("Stiffness Matrix"));	
	if (data.C.size() > 0) {
		for (int i = 0; i < 6; ++i) {
			for (int j = 0; j < 6; ++j) {
				if (!Hydro::C_units(i, j).IsEmpty()) {
					array.Set(row, 0, Format(t_("K(%d,%d) [%s]"), i+1, j+1, Hydro::C_units(i, j)));	array.Set(row++, col, Format("%12E", data.C(i, j)));		
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

	array.Set(row, 0, t_("# Segments"));		array.Set(row++, col, !healing ? Null : data.mesh.segments.size());
	array.Set(row, 0, t_("# Seg Waterplane"));	array.Set(row++, col, !healing ? Null : data.mesh.segWaterlevel.size());
	array.Set(row, 0, t_("# Seg leak"));		array.Set(row++, col, !healing ? Null : data.mesh.segTo1panel.size());
	array.Set(row, 0, t_("# Seg 3 panels"));	array.Set(row++, col, !healing ? Null : data.mesh.segTo3panel.size());

	array.Set(row, 0, t_("# Panels off"));		array.Set(row++, col, !healing ? Null : data.mesh.numUnprocessed);

	array.Set(row, 0, t_("# Triangles"));		array.Set(row++, col, !healing ? Null : data.mesh.numTriangles);
	array.Set(row, 0, t_("# BiQuads"));			array.Set(row++, col, !healing ? Null : data.mesh.numBiQuads);
	array.Set(row, 0, t_("# MonoQuads"));		array.Set(row++, col, !healing ? Null : data.mesh.numMonoQuads);
	array.Set(row, 0, t_("# Dup panels"));		array.Set(row++, col, !healing ? Null : data.mesh.numDupPan);
	array.Set(row, 0, t_("# Dup nodes"));		array.Set(row++, col, !healing ? Null : data.mesh.numDupP);
	array.Set(row, 0, t_("# Skewed pan"));		array.Set(row++, col, !healing ? Null : data.mesh.numSkewed);
}

void MainView::Init() {
	CtrlLayout(*this);
	main = &GetDefinedParent<MainMesh>(this);
	
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
			int id = ArrayModel_IdMesh(GetMain().listLoaded, row);
			if (id < 0)
				throw Exc(t_("Unexpected problem in OnPaint()"));
			
			double len = env.LenRef()/10;
			bool showNormals = ~GetMenuPlot().showNormals && ~GetMenuPlot().showMesh;
			bool showNormalsUnderwater = ~GetMenuPlot().showNormals && ~GetMenuPlot().showUnderwater;
			if (~GetMenuPlot().showNormals && !~GetMenuPlot().showMesh && !~GetMenuPlot().showUnderwater)
				showNormals = true;
			
			const Upp::Color &color = ArrayModel_GetColor(GetMain().listLoaded, row);
			const Mesh &mesh = Bem().surfs[id];
			
			gl.PaintSurface(mesh.mesh, color, ~GetMenuPlot().showMesh, 	
				showNormals);
				
			gl.PaintSurface(mesh.under, color, ~GetMenuPlot().showUnderwater, 
				showNormalsUnderwater);
			
			if (~GetMenuPlot().showSkewed)
				gl.PaintSegments(mesh.mesh.skewed, LtRed());
			if (~GetMenuPlot().showFissure)
				gl.PaintSegments(mesh.mesh.segTo1panel, LtRed());
			if (~GetMenuPlot().showWaterLevel)
				gl.PaintSegments(mesh.under.segWaterlevel, LtBlue());
			if (~GetMenuPlot().showMultiPan)
				gl.PaintSegments(mesh.mesh.segTo3panel, Black());
			
			if (~GetMenuPlot().showCb) {
				gl.PaintDoubleAxis(mesh.cb, len, LtBlue());
				gl.PaintCube(mesh.cb, len/10, LtBlue());
			}
			if (~GetMenuPlot().showCg) {
				gl.PaintDoubleAxis(mesh.cg, len, Black());
				gl.PaintCube(mesh.cg, len/5, Black());
			}
			if (~GetMenuPlot().showCr) {
				gl.PaintDoubleAxis(mesh.c0, len*10, Cyan());
				gl.PaintCube(mesh.c0, len/20, Gray());
			}
			if (paintSelect) {
				if (~GetMenuPlot().showMesh) {
					const UVector<int> &nod = mesh.mesh.GetSelNodes();
					for (int in = 0; in < nod.size(); ++in)
						gl.PaintCube(mesh.mesh.nodes[nod[in]], len/20, LtBlue());
					const UVector<int> &pan = mesh.mesh.GetSelPanels();
					const UVector<Point3D> &nodes = mesh.mesh.nodes;
					for (int ip = 0; ip < pan.size(); ++ip) {
						const Panel &panel = mesh.mesh.panels[pan[ip]];
						gl.PaintQuad(nodes[panel.id[0]], nodes[panel.id[1]], nodes[panel.id[2]], nodes[panel.id[3]], LtRed(), .2);
					}
				}
				if (~GetMenuPlot().showUnderwater) {
					const UVector<int> &nod = mesh.under.GetSelNodes();
					for (int in = 0; in < nod.size(); ++in)
						gl.PaintCube(mesh.under.nodes[nod[in]], len/20, LtBlue());
					const UVector<int> &pan = mesh.under.GetSelPanels();
					const UVector<Point3D> &nodes = mesh.under.nodes;
					for (int ip = 0; ip < pan.size(); ++ip) {
						const Panel &panel = mesh.under.panels[pan[ip]];
						gl.PaintQuad(nodes[panel.id[0]], nodes[panel.id[1]], nodes[panel.id[2]], nodes[panel.id[3]], LtRed(), .2);
					}
				}
			}
			if (~GetMenuPlot().showSel && ArrayModel_IsSelected(GetMain().listLoaded, row)) { 
				double minX = mesh.mesh.env.minX*.99; double maxX = mesh.mesh.env.maxX*1.01;
				double minY = mesh.mesh.env.minY*.99; double maxY = mesh.mesh.env.maxY*1.01;
				double minZ = mesh.mesh.env.minZ*.99; double maxZ = mesh.mesh.env.maxZ*1.01;
				
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

const WithMenuMeshPlot<StaticRect> &MainView::GetMenuPlot() const {
	return main->menuPlot;
}

void MainView::CalcEnvelope() {
	env.Reset();
	//for (int i = 0; i < Bem().surfs.size(); ++i)
	for (int row = 0; row < GetMain().listLoaded.GetCount(); ++row) {
		int id = ArrayModel_IdMesh(GetMain().listLoaded, row);
		if (id < 0)
			throw Exc("Unexpected problem in CalcEnvelope()");
		env.MixEnvelope(Bem().surfs[id].mesh.env);
	}
}

void MainMesh::LoadDragDrop() {
	GuiLock __;
	
	Sort(filesToDrop);
	for (int i = filesToDrop.size()-1; i > 0; --i)
		if (GetFileTitle(filesToDrop[i]) == GetFileTitle(filesToDrop[i-1]))
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
	
void MainMesh::DragAndDrop(Point , PasteClip& d) {
	GuiLock __;
	if (IsDragAndDropSource())
		return;
	if (AcceptFiles(d)) {
		filesToDrop = GetFiles(d);
		timerDrop.Set(0, [=] {LoadDragDrop();});
	}
}

bool MainMesh::Key(dword key, int ) {
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
		MainMesh &mm = GetDefinedParent<MainMesh>(this);
		
		idOpened = ArrayModel_IdMesh(mm.listLoaded);
		if (idOpened < 0) {
			Exclamation(t_("Please select a model to process"));
			return;
		}
		
		Clear(true);
		
		if (double(~edFrom) > double(~edTo)) {
			Exclamation(t_("Wrong Y angle range"));
			return;
		}
		double angleTo, angleDelta;
		if (edAngleTo.IsEnabled()) {
			angleTo = double(~edAngleTo);
			angleDelta =  double(~edAngleDelta);
			if (angleDelta == 0) {
				angleTo = double(~edAngleFrom);
				angleDelta = 1;
			}
		} else {
			angleTo = double(~edAngleFrom);
			angleDelta = 1;
		}
		if (double(~edAngleFrom) > angleTo) {
			Exclamation(t_("Wrong Z angle range"));
			return;
		}
		Mesh &mesh = Bem().surfs[idOpened];	
	
		int numAngle = 1 + int((angleTo - double(~edAngleFrom))/angleDelta);
		
		Progress progress(t_("GZ calculation..."), 100*numAngle); 
		
		datagz.Clear();
		dataMoment.Clear();
		
		scatter.SetTitle(Format(t_("GZ around Y axis at (%.2f, %.2f, %.2f)"), mesh.c0.x, mesh.c0.y, mesh.c0.z));
		
		String errors;
		int iangle = 0;
		for (double angle = double(~edAngleFrom); angle <= angleTo; angle += angleDelta, iangle++) {
			UVector<double> &dgz = datagz.Add();
			UVector<double> &dMoment = dataMoment.Add();
			UVector<double> vol, disp, wett, wplane, draft;
			UVector<Point3D> cb, cg;
			
			mesh.GZ(~edFrom, ~edTo, ~edDelta, angle, Bem().rho, Bem().g, double(~edTolerance)/100., 
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
				for (int iangle = 0; iangle < datagz.size(); ++iangle) {
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
			Exclamation(DeQtfLf("Errors found in mesh:\n" + errors));
	} catch(Exc e) {
		Exclamation(DeQtfLf(e));
		Clear(true);
	}
}

void MainGZ::Clear(bool force) {
	if (!force) {
		MainMesh &mm = GetDefinedParent<MainMesh>(this);
		
		UVector<int> ids = ArrayModel_IdsMesh(mm.listLoaded);
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
	array.Set(row++, 0, t_("Wetted area [m2]"));
	array.Set(row++, 0, t_("Waterpl. area [m2]"));
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
