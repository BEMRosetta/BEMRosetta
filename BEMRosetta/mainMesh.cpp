// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2025, the BEMRosetta author and contributors
#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <SurfaceCanvas/SurfaceCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <TabBar/TabBar.h>
#include <DropGrid/DropGrid.h>

#include <BEMRosetta_cl/BEMRosetta.h>

#define IMAGECLASS ImgM
#define IMAGEFILE <BEMRosetta/main.iml>
#include <Draw/iml.h>


using namespace Upp;

#include "main.h"

	
	
void MainBody::Init() {
	MainBEMBody::Init();
	
	OnOpt();
	
	CtrlLayout(menuOpen);
	menuOpen.file.WhenChange = THISBACK(OnLoad);
	menuOpen.file.BrowseRightWidth(40).UseOpenFolder(true).BrowseOpenFolderWidth(10);
	menuOpen.butLoad.WhenAction = [&] {menuOpen.file.DoGo();};
	menuOpen.rectangleButtons.SetBackground(SColorFace());
	
	ArrayModel_Init(listLoaded, true).MultiSelect();
	listLoaded.WhenSel = [&] {
		OnArraySel();
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
	menuProcess.butJoin.Disable();	
	menuProcess.butJoin.WhenAction = THISBACK(OnJoin);
	menuOpen.butDuplicate.Disable();	
	menuOpen.butDuplicate.WhenAction = THISBACK(OnDuplicate);
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
	menuPlot.showAxis.Tip(t_("Shows system axis")).WhenAction  					= [&] {mainView.FullRefresh(*this);};
	menuPlot.showLimits.Tip(t_("Shows boundaries of the geometry")).WhenAction 	= [&] {mainView.FullRefresh(*this);};
	menuPlot.showCb.Tip(t_("Shows the centre of buoyancy")).WhenAction  		= [&] {mainView.FullRefresh(*this);};
	menuPlot.showCg.Tip(t_("Shows the centre of gravity")).WhenAction  			= [&] {mainView.FullRefresh(*this);};
	menuPlot.showCr.Tip(t_("Shows the centre of motion")).WhenAction  		    = [&] {mainView.FullRefresh(*this);};
	menuPlot.showSel.Tip(t_("Shows volume around selected object")).WhenAction  = [&] {mainView.FullRefresh(*this);};
	menuPlot.showPoints.Tip(t_("Shows the points")).WhenAction  				= [&] {mainView.FullRefresh(*this);};
	menuPlot.showLines.Tip(t_("Shows the lines")).WhenAction  					= [&] {mainView.FullRefresh(*this);};
	menuPlot.showUnderwater.Tip(t_("Shows underwater mesh")).WhenAction  		= [&] {
		menuPlot.showBody.Enable(!menuPlot.showUnderwater);
		mainView.FullRefresh(*this);
	};
	menuPlot.butXYZ.Tip(t_("Orients the camera as isometric")).WhenAction  		= [&] {mainView.surf.SetRotation(Value3D(-M_PI/4, 0, M_PI/4));};
	menuPlot.butYoZ.Tip(t_("Orients the camera through X axis")).WhenAction  	= [&] {mainView.surf.SetRotation(Value3D(M_PI/2, M_PI, M_PI/2));};
	menuPlot.butXoZ.Tip(t_("Orients the camera through Y axis")).WhenAction  	= [&] {mainView.surf.SetRotation(Value3D(M_PI/2, M_PI, 0));};
	menuPlot.butXoY.Tip(t_("Orients the camera through Z axis")).WhenAction  	= [&] {mainView.surf.SetRotation(Value3D(0, 0, 0));};	
	menuPlot.butFit.Tip(t_("Zooms the camera to fit the bodies")).WhenAction	= [&] {mainView.surf.ZoomToFit();};
	menuPlot.showBodyData.Tip(t_("Shows a list of panels and nodes")).WhenAction= [&] {splitterAll.SetButton(0);};
	
	menuPlot.lineThickness.SetInc(1).Min(1).Tip(t_("Sets the thickness of the mesh wireframe")).WhenAction = [&] {LoadSelTab(Bem());};
	
	menuPlot.backColor.Tip(t_("Sets the background color")).WhenAction      	= [&] {mainView.RenderRefresh(*this);};
	menuPlot.lightX 	<< [=]{mainView.RenderRefresh(*this);};					menuPlot.lightX = 1;
	menuPlot.lightY 	<< [=]{mainView.RenderRefresh(*this);};					menuPlot.lightY = 1;
	menuPlot.lightZ 	<< [=]{mainView.RenderRefresh(*this);};					menuPlot.lightZ = 0;
			
	styleRed = styleGreen = styleBlue = Button::StyleNormal();
	styleRed.textcolor[0] = styleRed.textcolor[1] = styleRed.textcolor[2] = LtRed();
	styleGreen.textcolor[0] = styleGreen.textcolor[1] = styleGreen.textcolor[2] = Green();
	styleBlue.textcolor[0] = styleBlue.textcolor[1] = styleBlue.textcolor[2] = LtBlue();
	menuPlot.butYoZ.SetStyle(styleRed);
	menuPlot.butXoZ.SetStyle(styleGreen);
	menuPlot.butXoY.SetStyle(styleBlue);
	
	menuPlot.opShowColor.Add(SurfaceView::SHOW_DARKER, "SHOW_DARKER").Add(SurfaceView::SHOW_BRIGHTER, "SHOW_BRIGHTER")
				 .Add(SurfaceView::SHOW_FLAT, "SHOW_FLAT").SetData(0);
	menuPlot.opShowColor << [&] {mainView.RenderRefresh(*this);};
	
	menuPlot.opShowMesh.Add(SurfaceView::SHOW_MESH, "SHOW_MESH").Add(SurfaceView::SHOW_VISIBLE_MESH, "SHOW_VISIBLE_MESH")
				.Add(SurfaceView::SHOW_FACES, "SHOW_FACES").Add(SurfaceView::SHOW_MESH_FACES, "SHOW_MESH_FACES").SetData(3);
	menuPlot.opShowMesh << [&] {
		OnOpt();
		mainView.ViewRefresh(*this);
	};
	
	menuPlot.cx 		<< [=]{mainView.FullRefresh(*this);};					menuPlot.cx <<= 0;
	menuPlot.cy 		<< [=]{mainView.FullRefresh(*this);};					menuPlot.cy <<= 0;
	menuPlot.cz 		<< [=]{mainView.FullRefresh(*this);};					menuPlot.cz <<= 0;	
	
	CtrlLayout(menuAnimation);
	menuAnimation.edTime.Pattern("%.2f");
	menuAnimation.edTime 	<< [=]{
		if (!IsNull(menuAnimation.edTime)) {
			playing = true;
			menuAnimation.slider <<= ~menuAnimation.edTime;
			double tprev = tmGetTimeX();
			mainView.FullRefresh(*this);
			Ctrl::ProcessEvents();	
			double deltat = tmGetTimeX() - tprev;
			menuAnimation.labtime <<= Format("%.2f", 1./deltat);
			playing = false;
		}
	};				
	if (!IsNull(Bem().fast.GetTimeStart()))
		menuAnimation.edTime <<= Bem().fast.GetTimeStart();	

	menuAnimation.slider.Step(1, true);
	menuAnimation.slider.SetMajorTicksSize(4).SetThickness(2);
	if (!IsNull(Bem().fast.GetTimeStart()) && !IsNull(Bem().fast.GetTimeEnd()))
		menuAnimation.slider.MinMax((int)Bem().fast.GetTimeStart(), (int)Bem().fast.GetTimeEnd()).Step(1, true);	
	menuAnimation.slider <<= Bem().fast.GetTimeStart();	
	menuAnimation.slider.WhenAction = [&] {
		playing = true;
		menuAnimation.edTime <<= ~menuAnimation.slider;
		double tprev = tmGetTimeX();
		mainView.FullRefresh(*this);
		Ctrl::ProcessEvents();	
		double deltat = tmGetTimeX() - tprev;
		menuAnimation.labtime <<= Format("%.2f", 1./deltat);
		playing = false;
	};
	
	menuAnimation.butInit.SetImage(ImgM::resultset_first()).Tip(t_("Go to animation start")).WhenAction = [&]{
		playing = false;
		menuAnimation.edTime <<= Bem().fast.GetTimeStart();
		menuAnimation.slider <<= Bem().fast.GetTimeStart();	
	};
	menuAnimation.butPlayStep	.SetImage(ImgM::resultset_nextpause()).Tip(t_("Play step by step"));
	menuAnimation.butPlayStep.WhenPush   = [&]{PlayStep();};
	menuAnimation.butPlayStep.WhenRepeat = [&]{PlayStep();};
	menuAnimation.butPlay		.SetImage(ImgM::resultset_next())	  .Tip(t_("Play animation"))		.WhenAction = [&]{Play(true);};
	menuAnimation.butPlayReverse.SetImage(ImgM::resultset_previous()) .Tip(t_("Play animation reverse")).WhenAction = [&]{Play(false);};
	menuAnimation.butEnd		.SetImage(ImgM::resultset_last())	  .Tip(t_("Go to animation end"))   .WhenAction = [&]{
		playing = false;
		menuAnimation.edTime <<= Bem().fast.GetTimeEnd();
		menuAnimation.slider <<= Bem().fast.GetTimeEnd();	
	};
	
	menuAnimation.arraySpeed.AddColumn(t_("Speed"), 10);
	menuAnimation.arraySpeed.Add("x1");
	menuAnimation.arraySpeed.Add("x2");
	menuAnimation.arraySpeed.Add("x4");
	menuAnimation.arraySpeed.Add("x10");
	menuAnimation.arraySpeed.Add("x50");
	menuAnimation.arraySpeed.SetCursor(0);
	
	menuAnimation.butSetT0 << THISBACK(OnSetT0);
	
	mainView.surf.SetCanSelect();
	
	mainView.surf.WhenSelect << [=](int idBody, int idSubBody, bool select) {
		//if (splitterAll.GetPositionId() == 0 && idBody < 0)
		//	splitterAll.SetButton(0);	// Closes the right table		
		//else 
		if (splitterAll.GetPositionId() == 1 && idBody >= 0)
			splitterAll.SetButton(0);		// Opens the right table
		if (idBody >= 0 && idBody < 10000) 
			mainViewData.SelBody(idBody%10000, idSubBody, select);	// Selects the tabs in the right table. Select the panel
		else if (idBody >= 10000)
			Exclamation(t_("Only main body mesh can be edited"));
		//else
		//mainView.FullRefresh(*this);								// Updates view
	};
	
	mainView.FullRefresh(*this);
	
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
	menuProcess.butUpdateCg   << THISBACK2(OnUpdate, NONE, true);
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
	
	menuProcess.meshRatio <<= 1;
	menuProcess.meshRatio.Tip(t_("Ratio of generated mesh size to current mesh size"));
	menuProcess.meshRatio.Tip(t_("Mix triangles and quads, or set only quads"));
	menuProcess.lambda = 0.3;
	menuProcess.mu = -0.2;
	menuProcess.iterations = 10;
	menuProcess.butSmooth <<= THISBACK(OnSmooth);
	menuProcess.butSmooth.Tip(t_("Smooths the mesh using the Taubin method"));
	
	menuProcess.butInertia.SetCtrl(menuProcessInertia).Tip(t_("Click to get the inertia matrix"));
	//menuProcess.butInertia.SetAutoOpen();
	
	CtrlLayout(menuMove);	
	
	menuMove.butReset <<= THISBACK(OnReset);
	menuMove.butReset.Tip(t_("Translates the mesh to the original position"));
	
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
	menuMove.butArchimede.Tip(t_("Let the body fall to rest"));	
	
	menuMove.butPCA <<= THISBACK(OnPCA);
	menuMove.butPCA.Tip(t_("Orient the device to the principal axis"));	
	
	menuMove.opZArchimede.WhenAction = [&] {OnOpt();};
	menuMove.backArch.SetBackground(SColorFace());
	
	CtrlLayout(menuEdit);
	menuEdit.edit_x <<= 0;
	menuEdit.edit_y <<= 0;
	menuEdit.edit_z <<= 0;
	menuEdit.edit_size <<= 1;
	
	menuEdit.edit_cx <<= 0;
	menuEdit.edit_cy <<= 0;
	menuEdit.edit_cz <<= 0;
	menuEdit.opClose <<= true;
	menuEdit.butExtrude <<= THISBACK(OnExtrude);
	
	menuEdit.butPanel <<= THISBACK(OnAddPanel);
	menuEdit.panWidthX.WhenEnter << THISBACK(OnAddPanel);
	menuEdit.panWidthY.WhenEnter << THISBACK(OnAddPanel);
	
	menuEdit.butRevolution.SetCtrl(dialogRevolution).Tip(t_("Generates a revolution mesh"));
	menuEdit.butPolygon.SetCtrl(dialogPolygon).Tip(t_("Generates a polygonal flat panel"));
	
	menuEdit.butExtract <<= THISBACK(OnExtract);
	menuEdit.cutXYZ = 2;
	menuEdit.cutPosNeg = 1;
	
	menuEdit.rat_x <<= 1;
	menuEdit.rat_y <<= 1;
	menuEdit.rat_z <<= 1;
	menuEdit.butScale <<= THISBACK(OnScale);
	menuEdit.butScale.Tip(t_("Scales the mesh"));
	
	menuEdit.butImageX <<= THISBACK1(OnImage, 0);
	menuEdit.butImageX.Tip(t_("Mirrors the mesh in X axis"));
	menuEdit.butImageY <<= THISBACK1(OnImage, 1);
	menuEdit.butImageY.Tip(t_("Mirrors the mesh in Y axis"));
	menuEdit.butImageZ <<= THISBACK1(OnImage, 2);
	menuEdit.butImageZ.Tip(t_("Mirrors the mesh in Z axis"));
	
	dialogRevolution.butOK << THISBACK(OnAddRevolution);
	dialogPolygon.butOK << THISBACK(OnAddPolygonalPanel);
	
	CtrlLayout(menuStability);
	menuStability.heelingMoment <<= 0;
	
	struct MyConvert : public Convert {
		virtual Value Format(const Value& q) const	{return FDS(q, 7);}	
		virtual Value Scan(const Value& text) const	{return ScanDouble(text.ToString());}
	};

	menuStability.heelingMoment.SetConvert(Single<MyConvert>());
	menuStability.butPointsA.SetCtrl(dialogPointsA).Tip(t_("Unprotected points (lead to progressive flooding)"));
	dialogPointsA.WhenClose = [&] {
		int idx = ArrayModel_IndexBody(listLoaded);
		if (idx < 0) 
			return;
		
		Bem().surfs[idx].Reset(Bem().rho, Bem().g);
		dialogPointsA.FromGrid(Bem().surfs[idx].cdt.controlPointsA);
		Bem().surfs[idx].cdt.controlPointsA0 = clone(Bem().surfs[idx].cdt.controlPointsA);
		
		menuStability.labPointsA.SetText(Format(t_("%d points"), dialogPointsA.grid.GetRowCount()));
		menuStability.labPointsA.SetFont(menuStability.labPointsA.GetFont().Bold(dialogPointsA.grid.GetRowCount() > 0));
	};
	menuStability.butPointsB.SetCtrl(dialogPointsB).Tip(t_("Weathertight openings (lead to progressive flooding)"));
	dialogPointsB.WhenClose = [&] {
		int idx = ArrayModel_IndexBody(listLoaded);
		if (idx < 0) 
			return;
		
		Bem().surfs[idx].Reset(Bem().rho, Bem().g);
		dialogPointsB.FromGrid(Bem().surfs[idx].cdt.controlPointsB);
		Bem().surfs[idx].cdt.controlPointsB0 = clone(Bem().surfs[idx].cdt.controlPointsB);
		
		menuStability.labPointsB.SetText(Format(t_("%d points"), dialogPointsB.grid.GetRowCount()));
		menuStability.labPointsB.SetFont(menuStability.labPointsB.GetFont().Bold(dialogPointsB.grid.GetRowCount() > 0));
	};
			
	menuStability.butPointsC.SetCtrl(dialogPointsC).Tip(t_("Watertight openings (do not lead to progressive flooding)"));
	dialogPointsC.WhenClose = [&] {
		int idx = ArrayModel_IndexBody(listLoaded);
		if (idx < 0) 
			return;
		
		Bem().surfs[idx].Reset(Bem().rho, Bem().g);
		dialogPointsC.FromGrid(Bem().surfs[idx].cdt.controlPointsC);
		Bem().surfs[idx].cdt.controlPointsC0 = clone(Bem().surfs[idx].cdt.controlPointsC);
		
		menuStability.labPointsC.SetText(Format(t_("%d points"), dialogPointsC.grid.GetRowCount()));
		menuStability.labPointsC.SetFont(menuStability.labPointsC.GetFont().Bold(dialogPointsC.grid.GetRowCount() > 0));
	};
	
	menuStability.butLoads.SetCtrl(dialogLoads).Tip(t_("Loads"));
	dialogLoads.WhenClose = [&] {
		int idx = ArrayModel_IndexBody(listLoaded);
		if (idx < 0) 
			return;
		
		Bem().surfs[idx].Reset(Bem().rho, Bem().g);
		dialogLoads.FromGrid(Bem().surfs[idx].cdt.controlLoads);
		Bem().surfs[idx].cdt.controlLoads0 = clone(Bem().surfs[idx].cdt.controlLoads);
		
		menuStability.labLoads.SetText(Format(t_("%d loads"), dialogLoads.grid.GetRowCount()));
		menuStability.labLoads.SetFont(menuStability.labLoads.GetFont().Bold(dialogLoads.grid.GetRowCount() > 0));
	};
	
	menuStability.butDamage.SetCtrl(dialogDamage).Tip(t_("Damage"));
	dialogDamage.WhenClose = [&] {
		int idx = ArrayModel_IndexBody(listLoaded);
		if (idx < 0) 
			return;
		int num = 0;
		Bem().surfs[idx].cdt.damagedBodies.Clear();
		for (int r = 0; r < dialogDamage.grid.GetRowCount(); ++r) {
			int iddx = Bem().GetBodyIndex(dialogDamage.grid.Get(r, 0));
			if (iddx != idx && dialogDamage.grid.Get(r, 1) == true) {
				Bem().surfs[idx].cdt.damagedBodies << &(Bem().surfs[iddx]);
				num++; 
			}
		}
		menuStability.labDamage.SetText(Format(t_("%d damage"), num));
		menuStability.labDamage.SetFont(menuStability.labDamage.GetFont().Bold(num > 0));
	};
	
	menuStability.butArchimede <<= THISBACK(OnArchimede);
	menuStability.butArchimede.Tip(t_("Let the body fall to rest"));	
	
	menuStability.file.WhenChange = [&]() {menuStability.butLoad.WhenAction(); return true;};
	menuStability.file.BrowseRightWidth(40).UseOpenFolder().BrowseOpenFolderWidth(10);
	menuStability.file.ActiveDir(saveFolder);
	menuStability.file.Type(t_("Openings and loads"), "*.json");
	menuStability.file.ActiveType(0);
	menuStability.butLoad.WhenAction = [&]() {
		int idx = ArrayModel_IndexBody(listLoaded);
		if (idx < 0)
			return;
		if (!LoadFromJsonFile(Bem().surfs[idx].cdt, ~menuStability.file))
			Exclamation(t_("Impossible to load file"));
		OnMenuStabilityArraySel();
		Bem().surfs[idx].Reset(Bem().rho, Bem().g);
	};
	
	menuStability.butSave.WhenAction = [&]() {
		int idx = ArrayModel_IndexBody(listLoaded);
		if (idx < 0)
			return;
		
		FileSel fs;
		fs.Type(t_("Openings and loads"), "*.json");
		fs.ActiveDir(saveFolder);
		fs.ActiveType(0);
		fs.Set(ForceExtSafer(~menuStability.file, ".json"));
		
		if (!fs.ExecuteSaveAs(t_("Save openings and loads data"))) 
			Exclamation(t_("Cancelled by the user"));
		
		String fileName = ~fs;
				
		if (!StoreAsJsonFile(Bem().surfs[idx].cdt, fileName, true))
			Exclamation(t_("Impossible to save file"));
	};
	
	menuTab.Add(menuOpen.SizePos(),    	t_("File"));
	menuTab.Add(menuPlot.SizePos(),    	t_("Plot")).Disable();
	menuTab.Add(menuMove.SizePos(), 	t_("Move")).Disable();
	menuTab.Add(menuProcess.SizePos(), 	t_("Process")).Disable();
	menuTab.Add(menuStability.SizePos(),t_("Stability")).Disable();
	menuTab.Add(menuEdit.SizePos(), 	t_("Edit"));
	menuTab.Add(menuAnimation.SizePos(),t_("Animation")).Disable();
	
	menuStability.butClear.WhenAction = [&]() {
		int idx = ArrayModel_IndexBody(listLoaded);
		if (idx < 0)
			return;
		
		Bem().surfs[idx].cdt.Reset();
		OnMenuStabilityArraySel();
		Bem().surfs[idx].Reset(Bem().rho, Bem().g);
	};
	
	mainViewData.Init();
	splitterAll.Horz(mainView.SizePos(), mainViewData.SizePos());
	splitterAll.SetPositions(6000, 9900).SetInitialPositionId(1).SetButtonNumber(1).SetButtonWidth(20);
	splitterAll.Tip(t_("")).WhenAction = [&] {mainView.SetPaintSelect(splitterAll.GetPos() < 9900);};
	mainTab.Add(splitterAll.SizePos(), t_("View"));
	
	String bitmapFolder = AFX(GetDesktopFolder(), "BEMRosetta Body Images");
	int idBitmapFolder = 0;
	
	mainSummary.Init();
	mainTab.Add(mainSummary.SizePos(), t_("Summary"));

	mainStiffness.Init(Hydro::MAT_K);
	mainTab.Add(mainStiffness.SizePos(), t_("Hydrostatic Stiffness")).Disable();
	
	mainStiffness.opMassBuoy.Tip(t_("Obtain stiffness matrix including the effect of mass and buoyancy"));
	mainStiffness.opMassBuoy.WhenAction = THISBACK2(OnUpdate, NONE, true);
	mainStiffness.opMassBuoy.Hide();
	
	mainM.Init(Hydro::MAT_M);
	mainTab.Add(mainM.SizePos(), t_("Inertia")).Disable();
	mainM.label.SetText(t_("Inertia Matrices"));
	mainM.opMassBuoy.Hide();
	
	mainStiffness2.Init(Hydro::MAT_KMOOR);
	mainTab.Add(mainStiffness2.SizePos(), t_("Mooring Stiffness")).Disable();
	mainStiffness2.label.SetText(t_("Mooring Stiffness Matrices"));
	mainStiffness2.opMassBuoy.Hide();
	
	mainGZ.Init();
	mainTab.Add(mainGZ.SizePos(), t_("Stability GZ")).Disable();
			
	mainTab.WhenSet = [&] {
		LOGTAB(mainTab);
		UVector<int> idxs = ArrayModel_IndexsBody(listLoaded);
		bool plot = true, move = false, convertProcess = true, stability = false;
		if (Bem().surfs.IsEmpty()) 
			plot = convertProcess = false;
		else if (mainTab.IsAt(splitterAll)) 
			;
		else if (mainTab.IsAt(mainM)) {
			plot = false;
			move = true;
			mainM.Load(Bem().surfs, idxs, false);
		} else if (mainTab.IsAt(mainStiffness)) {
			plot = false;
			move = true;
			mainStiffness.Load(Bem().surfs, idxs, true);
		} else if (mainTab.IsAt(mainStiffness2)) {
			plot = false;
			move = true;
			mainStiffness2.Load(Bem().surfs, idxs, true);
		} else if (mainTab.IsAt(mainGZ)) 
			stability = true;
		else 
			plot = false;
		
		TabCtrl::Item& tabMenuPlot = menuTab.GetItem(menuTab.Find(menuPlot));
		tabMenuPlot.Enable(plot);
		TabCtrl::Item& tabMenuAnimation = menuTab.GetItem(menuTab.Find(menuAnimation));
		tabMenuAnimation.Enable(plot && !Bem().fast.IsEmpty());
		TabCtrl::Item& tabMenuMove = menuTab.GetItem(menuTab.Find(menuMove));
		tabMenuMove.Enable(convertProcess);
		TabCtrl::Item& tabMenuProcess = menuTab.GetItem(menuTab.Find(menuProcess));
		tabMenuProcess.Enable(convertProcess);
		TabCtrl::Item& tabMenuStability = menuTab.GetItem(menuTab.Find(menuStability));
		tabMenuStability.Enable(convertProcess);
		if (stability)
			menuTab.Set(menuStability);
		else if (plot) {
			tabMenuPlot.Text(t_("Plot"));
			menuTab.Set(menuPlot);
		} else {
			tabMenuPlot.Text("");
			if (move)
				menuTab.Set(menuMove);
			else
				menuTab.Set(menuOpen);
		}
		
		if (plot && !Bem().fast.IsEmpty())
			tabMenuAnimation.Text(t_("Animation"));
		else
			tabMenuAnimation.Text(t_(""));
		
		if (convertProcess) {
			tabMenuProcess.Text(t_("Process"));
			tabMenuMove.Text(t_("Move"));
			tabMenuStability.Text(t_("Stability"));
		} else {
			tabMenuProcess.Text("");
			tabMenuMove.Text("");
			tabMenuStability.Text("");
		}
	};
	mainTab.WhenSet();
	
	menuTab.WhenSet = [&] {
		LOGTAB(menuTab);
		
		if (playing) {
			menuTab.Set(menuAnimation);
			return;
		}
		if (menuTab.IsAt(menuStability))
			mainTab.Set(mainGZ);
	};
	menuTab.WhenSet();
	
	UpdateButtons();
	saveFolder = GetDesktopFolder();
}

void MainBody::PlayEnableCtrls(bool enable, int forward) {
	menuAnimation.butInit.Enable(enable);
	menuAnimation.butEnd.Enable(enable);
	menuAnimation.butPlayStep.Enable(enable);
	menuAnimation.arraySpeed.Enable(enable);
	menuAnimation.edTime.Enable(enable);
	menuAnimation.slider.Enable(enable);
	menuAnimation.butPlay.Enable(enable);
	menuAnimation.butPlayReverse.Enable(enable);
	menuAnimation.butPlayStep.Enable(enable);
	if (forward == 1)
		menuAnimation.butPlay.Enable(true);
	else if (forward == 0)
		menuAnimation.butPlayStep.Enable(true);
	else
		menuAnimation.butPlayReverse.Enable(true);
}

void MainBody::PlayStep() {
	playing = true;
	String sspeed = menuAnimation.arraySpeed.Get(menuAnimation.arraySpeed.GetCursor(), 0);
	double speed = ScanDouble(sspeed.Mid(1));
	
	double tm = ~menuAnimation.edTime;
	PlayEnableCtrls(false, 0);
	double tprev = tmGetTimeX();

	mainView.FullRefresh(*this);
	Ctrl::ProcessEvents();	
	double t = tmGetTimeX();
	double deltat = t - tprev;

	tm += deltat*speed;
	if (tm >= Bem().fast.GetTimeEnd())
		tm = Bem().fast.GetTimeEnd();

	menuAnimation.labtime <<= Format("%.2f", 1./deltat);
	menuAnimation.edTime <<= tm;
	menuAnimation.slider <<= tm;
	
	PlayEnableCtrls(true, 0);	
	playing = false;
}

void MainBody::Play(bool forward) {
	if (playing) {
		playing = false;
		if (forward)
			menuAnimation.butPlay		.SetImage(ImgM::resultset_next())	 .Tip(t_("Play animation"));
		else
			menuAnimation.butPlayReverse.SetImage(ImgM::resultset_previous()).Tip(t_("Play animation reverse"));
	} else {
		playing = true;
		
		if (forward)
			menuAnimation.butPlay		.SetImage(ImgM::resultset_pause()).Tip(t_("Pause animation"));
		else
			menuAnimation.butPlayReverse.SetImage(ImgM::resultset_pause()).Tip(t_("Pause animation"));
		
		String sspeed = menuAnimation.arraySpeed.Get(menuAnimation.arraySpeed.GetCursor(), 0);
		double speed = ScanDouble(sspeed.Mid(1));
		
		double tm = ~menuAnimation.edTime;
		PlayEnableCtrls(false, forward ? 1 : -1);
		double tprev = tmGetTimeX();
		long count = 0;
		while (playing) {
			mainView.FullRefresh(*this);
			Ctrl::ProcessEvents();	
			double t = tmGetTimeX();
			double deltat = t - tprev;
			if (forward) {
				tm += deltat*speed;
				if (tm >= Bem().fast.GetTimeEnd()) {
					tm = Bem().fast.GetTimeEnd();
					Play(forward);
				}
			} else {
				tm -= deltat*speed;
				if (tm <= Bem().fast.GetTimeStart()) {
					tm = Bem().fast.GetTimeStart();
					Play(forward);
				}
			}
			if (!(count%20))
				menuAnimation.labtime <<= Format("%.2f", 1./deltat);
			menuAnimation.edTime <<= tm;
			menuAnimation.slider <<= tm;
			tprev = t;
			count++;
		}
		PlayEnableCtrls(true, forward ? 1 : -1);
	}
}			
	
void MainBody::OnMenuOpenArraySel() {
	int idx = ArrayModel_IndexBody(listLoaded);
	if (idx < 0)
		return;
	
	Body::MESH_FMT type = Body::GetCodeBodyStr(~menuOpen.dropExport);
	menuOpen.symX <<= ((type == Body::WAMIT_GDF || type == Body::AQWA_DAT) && Bem().surfs[idx].IsSymmetricX());
	menuOpen.symY <<= Bem().surfs[idx].IsSymmetricY();
	
	dialogDamage.SelectId(idx);
}

void MainBody::OnMenuProcessArraySel() {
	int idx = ArrayModel_IndexBody(listLoaded);
	if (idx < 0)
		return;
	
	Body &msh = Bem().surfs[idx];
	if (!IsNull(msh.dt.cg)) {
		menuProcess.x_g <<= msh.dt.cg.x;
		menuProcess.y_g <<= msh.dt.cg.y;
		menuProcess.y_g <<= msh.dt.cg.y;
		menuProcess.z_g <<= msh.dt.cg.z;
	} else
		menuProcess.x_g <<= menuProcess.y_g <<= menuProcess.z_g <<= Null;
	
	menuProcess.x_0 <<= msh.dt.c0.x;
	menuProcess.y_0 <<= msh.dt.c0.y;
	menuProcess.z_0 <<= msh.dt.c0.z;
	menuProcess.mass <<= msh.GetMass();
	
	menuProcessInertia.Init(*this, idx);
}

void MainBody::OnMenuMoveArraySel() {
	int id = ArrayModel_IndexBody(listLoaded);
	if (id < 0)
		return;
}

void MainBody::OnMenuAdvancedArraySel() {
	int id = ArrayModel_IndexBody(listLoaded);
	if (id < 0)
		return;
}

void MainBody::OnMenuStabilityArraySel() {
	int idx = ArrayModel_IndexBody(listLoaded);
	if (idx < 0)
		return;
	
	dialogPointsA.ToGrid(Bem().surfs[idx].cdt.controlPointsA0);
	dialogPointsB.ToGrid(Bem().surfs[idx].cdt.controlPointsB0);
	dialogPointsC.ToGrid(Bem().surfs[idx].cdt.controlPointsC0);
	dialogLoads.ToGrid(Bem().surfs[idx].cdt.controlLoads0);
	
	menuStability.labPointsA.SetText(Format(t_("%d points"), dialogPointsA.grid.GetRowCount()));
	menuStability.labPointsA.SetFont(menuStability.labPointsA.GetFont().Bold(dialogPointsA.grid.GetRowCount() > 0));
	menuStability.labPointsB.SetText(Format(t_("%d points"), dialogPointsB.grid.GetRowCount()));
	menuStability.labPointsB.SetFont(menuStability.labPointsB.GetFont().Bold(dialogPointsB.grid.GetRowCount() > 0));
	menuStability.labPointsC.SetText(Format(t_("%d points"), dialogPointsC.grid.GetRowCount()));
	menuStability.labPointsC.SetFont(menuStability.labPointsC.GetFont().Bold(dialogPointsC.grid.GetRowCount() > 0));
	menuStability.labLoads.SetText(Format(t_("%d loads"), dialogLoads.grid.GetRowCount()));	
	menuStability.labLoads.SetFont(menuStability.labLoads.GetFont().Bold(dialogLoads.grid.GetRowCount() > 0));
}

void MainBody::OnArraySel() {
	OnMenuOpenArraySel();
	OnMenuProcessArraySel();
	OnMenuMoveArraySel();
	OnMenuAdvancedArraySel();
	OnMenuStabilityArraySel();
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
	if (!ret || IsNull(menuPlot.showPoints)) 
		menuPlot.showPoints = true;
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
	const UVector<int> &idxs = ArrayModel_IndexsBody(listLoaded);
	if (mainTab.Get() == mainTab.Find(mainStiffness))
		mainStiffness.Load(bem.surfs, idxs, true);
	else //if (mainTab.Get() == mainTab.Find(mainView))
		mainView.FullRefresh(*this);
}

void MainBody::OnOpt() {
	menuOpen.file.ClearTypes(); 

	String meshFiles = Body::GetMeshExt();	//".gdf .dat .stl .pnl .msh .mesh .hst .grd .obj .nc";
	
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
	case Body::MIKE21_GRD:	menuOpen.symX.Enable();
							menuOpen.symY.Enable();
							break;		
	case Body::BEM_MESH:	menuOpen.symX.Disable();
							menuOpen.symY.Disable();
	default:				break;		
	}
	
	menuOpen.optBodyType.Enable(type != Body::BEM_MESH);
	
	menuPlot.opShowColor.Enable(menuPlot.opShowMesh.GetIndex() != SurfaceView::SHOW_MESH && 
								menuPlot.opShowMesh.GetIndex() != SurfaceView::SHOW_VISIBLE_MESH);
								
	menuAnimation.Enable(!Bem().fast.IsEmpty());							
	
	menuMove.t_z.Enable(!menuMove.opZArchimede);
	
	if (!Bem().fast.IsEmpty()) {
		double tm = ~menuAnimation.edTime;
		if (IsNull(tm) || tm < Bem().fast.GetTimeStart() || tm > Bem().fast.GetTimeEnd())
			menuAnimation.edTime <<= tm = Bem().fast.GetTimeStart();
		menuAnimation.slider.MinMax((int)Bem().fast.GetTimeStart(), (int)Bem().fast.GetTimeEnd());
		menuAnimation.slider.SetMajorTicks().SetMinorTicks(menuAnimation.slider.GetMajorTicks()/4);
		menuAnimation.slider <<= tm;
	}
}

void MainBody::AfterAdd(String file, int num, bool firstLoad) {
	mainTab.Set(mainView);
	for (int idx = Bem().surfs.size() - num; idx < Bem().surfs.size(); ++idx) {
		Body &surf = Bem().surfs[idx];
		
		surf.Report(Bem().rho);
		
		AddRow(surf);
		
		UpdateLast(idx);
	}
	
	mainView.surf.Enable();
	mainView.surf.Render();
	if (firstLoad)
		mainView.surf.SetRotationXYZ();
	
	After();	
	mainView.surf.ZoomToFit();
	
	OnArraySel();	
}

bool MainBody::OnLoad() {
	GuiLock __;
	
	String file = ~menuOpen.file;
		
	try {
		Progress progress(t_("Loading mesh file..."), 100); 
		
		UVector<int> idxs = ArrayModel_IndexsBody(listLoaded);
		for (int i = 0; i < idxs.size(); ++i) {
			if (Bem().surfs[idxs[i]].dt.fileName == file) {
				if (!PromptYesNo(t_("Model is already loaded") + S("&") + t_("Do you wish to open it anyway?")))
					return false;
				break;
			}
		}
		mainView.surf.Disable();

		bool firstLoad = Bem().surfs.IsEmpty();
		
		WaitCursor waitcursor;

		int num = Bem().LoadBody(file, [&](String str, int _pos) {
			progress.SetText(str); 
			progress.SetPos(_pos); 
			return !progress.Canceled();
		}, ~menuOpen.opClean, false, idxs);
		
		menuAnimation.animationFile <<= Bem().fast.GetFileName();
		
		for (int idx = Bem().surfs.size() - num; idx < Bem().surfs.size(); ++idx) {
			Body &msh = Bem().surfs[idx];
			
			if (!msh.dt.mesh.IsEmpty() && ~menuMove.opZArchimede && msh.GetMass_all() > 0) {
				double dz = 0.1;
				if (msh.GetMass_all() == 0)
					throw Exc(t_("Set mass before fitting buoyancy"));
				if (!msh.TranslateArchimede(Bem().rho, 0.05, dz)) 
					Exclamation(t_("Problem readjusting the Z value to comply with displacement (Archimede).&Mesh loaded as-is"));
					
				msh.AfterLoad(Bem().rho, Bem().g, false, false, true);
				
				Ma().Status(Format(t_("Loaded '%s', and translated vertically %f m to comply with displacement"), file, dz));
			} else
				Ma().Status(Format(t_("Loaded '%s'"), file));			
		}
		
		AfterAdd(file, num, firstLoad);
		
		if (num > 0)
			mainViewData.OnAddedModel(mainView);
		OnOpt();
	} catch (Exc e) {
		mainView.surf.Enable();
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
				sel << ArrayModel_IndexBody(listLoaded, 0);
			else {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		} else {
			for (int i = 0; i < sel.size(); ++i)
				sel[i] = ArrayModel_IndexBody(listLoaded, sel[i]);//Bem().GetBodyIndex(i);
		}
		
		if (type == Body::AQWA_DAT) {
			bool mismatch = false;
			for (int i = 0; i < sel.size(); ++i) {
				if (Bem().surfs[sel[i]].dt.c0 != Bem().surfs[sel[i]].dt.cg) {
					mismatch = true;
					break;
				}
			}
			if (mismatch && !PromptOKCancel(t_("In some body the centre of gravity is not the same that the body axis.&Do you wish to continue?")))
				return;
		}
		
		Status(t_("Saving mesh data"));
		
		FileSel fs;
		
		for (int i = 0; i < Body::NUMMESH; ++i)
			if (Body::meshCanSave[i] && (i == type || i == Body::UNKNOWN)) 
				fs.Type(Body::GetBodyStr(static_cast<Body::MESH_FMT>(i)), Body::meshExt[i]);
		
		
		String fileName;
		int id = ArrayModel_IndexBody(listLoaded);
		if (id < 0)
			fileName = ~menuOpen.file;
		else
			fileName = ArrayModel_GetFileName(listLoaded, id);
		
		fs.ActiveType(0);
		fs.Set(ForceExtSafer(fileName, ext));
		fs.ActiveDir(saveFolder);
		
		if (!fs.ExecuteSaveAs(Format(t_("Save mesh file as %s"), fileType)))
			return;
		
		fileName = ~fs;
		
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
		UVector<int> idxs = ArrayModel_IndexsBody(listLoaded);
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		int id;
		if (num == 0 && listLoaded.GetCount() == 1)
			id = ArrayModel_IndexBody(listLoaded, 0);
		else {
		 	id = ArrayModel_IndexBody(listLoaded);
			if (id < 0) {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		}
				
		Body &msh = Bem().surfs[id];
		msh.Reset(Bem().rho, Bem().g);

		UpdateLast(id);
		
		mainView.FullRefresh(*this);
		mainViewData.OnRefresh();

		menuProcess.x_g <<= msh.dt.cg.x; 
		menuProcess.y_g <<= msh.dt.cg.y;
		menuProcess.z_g <<= msh.dt.cg.z;
		
		Ma().Status(t_("Model oriented on the initial layout"));
		
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
}

void MainBody::OnSetT0() {
	GuiLock __;
	
	try {
		UVector<int> idxs = ArrayModel_IndexsBody(listLoaded);
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		int id;
		if (num == 0 && listLoaded.GetCount() == 1)
			id = ArrayModel_IndexBody(listLoaded, 0);
		else {
		 	id = ArrayModel_IndexBody(listLoaded);
			if (id < 0) {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		}
				
		Body &msh = Bem().surfs[id];
		msh.SetT0(Bem().rho, Bem().g);

		UpdateLast(id);
		
		mainView.FullRefresh(*this);
		mainViewData.OnRefresh();

		menuProcess.x_g <<= msh.dt.cg.x; 
		menuProcess.y_g <<= msh.dt.cg.y;
		menuProcess.z_g <<= msh.dt.cg.z;
		
		Ma().Status(t_("Model is set the initial layout"));
		
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
}

void MainBody::OnUpdateMass() {
	GuiLock __;
	
	try {
		UVector<int> ids = ArrayModel_IndexsBody(listLoaded);
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		int idx;
		if (num == 0 && listLoaded.GetCount() == 1)
			idx = ArrayModel_IndexBody(listLoaded, 0);
		else {
		 	idx = ArrayModel_IndexBody(listLoaded);
			if (idx < 0) {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		}
		
		Body &msh = Bem().surfs[idx];
		
		if (msh.dt.under.volume <= 0) {
			BEM::PrintError(t_("No body volume is under the water"));
			return;
		}
		if (msh.dt.under.VolumeMatch(0.05, 0.05) < 0) {
			BEM::PrintError(t_("Underwater mesh have a problem. Please review normals"));
			return;
		}
		
		double mass = msh.dt.under.volume*Bem().rho;
		menuProcess.mass <<= mass;
		menuProcessInertia.mass <<= mass;

		OnUpdate(NONE, true);
		
		Ma().Status(Format(t_("Mass updated to %f kg"), mass));
		
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MainBody::OnScale() {
	GuiLock __;
	
	try {
		UVector<int> ids = ArrayModel_IndexsBody(listLoaded);
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		int id;
		if (num == 0 && listLoaded.GetCount() == 1)
			id = ArrayModel_IndexBody(listLoaded, 0);
		else {
		 	id = ArrayModel_IndexBody(listLoaded);
			if (id < 0) {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		}
				
		Body &msh = Bem().surfs[id];
	
		WaitCursor wait;
		
		double rx = double(~menuEdit.rat_x) - 1;
		double ry = double(~menuEdit.rat_y) - 1;
		double rz = double(~menuEdit.rat_z) - 1;
		
		msh.dt.mesh.Scale(rx, ry, rz, msh.dt.c0);
		if (!IsNull(msh.dt.cg))
			msh.dt.cg.Translate(rx*(msh.dt.cg.x - msh.dt.c0.x), ry*(msh.dt.cg.y - msh.dt.c0.y),
							    rz*(msh.dt.cg.z - msh.dt.c0.z));
		
		Ma().Status(Format(t_("Model scaled %f, %f, %f"), rx, ry, rz));
				
		msh.AfterLoad(Bem().rho, Bem().g, NONE, false);
			
		UpdateLast(id);
			
		mainView.FullRefresh(*this);
		mainViewData.OnRefresh();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}		
}

void MainBody::OnSmooth() {
	GuiLock __;
	
	try {
		UVector<int> ids = ArrayModel_IndexsBody(listLoaded);
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		int id;
		if (num == 0 && listLoaded.GetCount() == 1)
			id = ArrayModel_IndexBody(listLoaded, 0);
		else {
		 	id = ArrayModel_IndexBody(listLoaded);
			if (id < 0) {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		}
				
		Body &msh = Bem().surfs[id];
	
		WaitCursor wait;
		
		double lambda = double(~menuProcess.lambda);
		double mu = double(~menuProcess.mu);
		
		if (IsNull(lambda) && IsNull(mu)) {
			BEM::PrintError(t_("Please fill the parameters"));
			return;
		}
		if (!IsNull(mu) && mu > 0) {
			BEM::PrintError(t_("Mu (Laplacian backward ratio) has to be negative"));
			return;
		}
		if (!IsNull(mu) && lambda < abs(mu)) {
			BEM::PrintError(t_("Lambda (Laplacian forward ratio) has to be higher than Mu (Laplacian backward ratio)"));
			return;
		}
		int iterations = int(~menuProcess.iterations);
		
		msh.dt.mesh.SmoothTaubin(lambda, mu, iterations, 0);		// 0 weightless, 1, weighted, 2, implicit (to be implemented)
		
		Ma().Status(t_("Model smoothed"));
				
		msh.AfterLoad(Bem().rho, Bem().g, NONE, false);
			
		UpdateLast(id);
			
		mainView.FullRefresh(*this);
		mainViewData.OnRefresh();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}		
}

void MainBody::OnUpdate(Action action, bool fromMenuProcess) {
	GuiLock __;
	
	try {
		UVector<int> idxs = ArrayModel_IndexsBody(listLoaded);
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		int idx;
		if (num == 0 && listLoaded.GetCount() == 1)
			idx = ArrayModel_IndexBody(listLoaded, 0);
		else {
		 	idx = ArrayModel_IndexBody(listLoaded);
			if (idx < 0) {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		}
				
		Body &msh = Bem().surfs[idx];

		double mass = ~menuProcess.mass;
		menuProcessInertia.mass <<= mass;
		
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

		msh.SetMass(mass);

		if (action == NONE && (IsNull(mass) || mass <= 0)) {
			BEM::PrintError(t_("Please set the mass"));
			return;
		}
		if (action == NONE && (IsNull(x_g) || IsNull(y_g) || IsNull(z_g))) {
			BEM::PrintError(t_("Please set the CG"));
			return;
		}
		if (action == NONE && (IsNull(x_0) || IsNull(y_0) || IsNull(z_0))) {
			BEM::PrintError(t_("Please set centre of rotation centre"));
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

		if (~menuMove.opZArchimede) 
			if (msh.GetMass_all() == 0)
				throw Exc(t_("Set mass before fitting buoyancy"));
				
		if (action == MOVE) {
			msh.Translate(t_x, t_y, t_z);
			
			if (~menuMove.opZArchimede) {
				double dz = 0;
				if (!msh.TranslateArchimede(Bem().rho, 0.05, dz))
					BEM::PrintError(t_("Problem readjusting the Z value to comply with displacement (Archimede)"));
				else
					Ma().Status(Format(t_("Model moved %f, %f, %f"), t_x, t_y, t_z + dz));
			} else
				Ma().Status(Format(t_("Model moved %f, %f, %f"), t_x, t_y, t_z));
			
		} else if (action == ROTATE) {
			msh.Rotate(ToRad(a_x), ToRad(a_y), ToRad(a_z), x_0, y_0, z_0);
			
			if (~menuMove.opZArchimede) {
				double dz = 0;
				if (!msh.TranslateArchimede(Bem().rho, 0.05, dz))
					BEM::PrintError(t_("Problem readjusting the Z value to comply with displacement (Archimede)"));
				else
					Ma().Status(Format(t_("Model rotated %f, %f, %f deg. around %f, %f, %f, and translated vertically %f m to comply with displacement"), a_x, a_y, a_z, x_0, y_0, z_0, dz));
			} else
				Ma().Status(Format(t_("Model rotated %f, %f, %f around %f, %f, %f"), a_x, a_y, a_z, x_0, y_0, z_0));
		} else if (action == NONE) {
			Point3D ncg(x_g, y_g, z_g);
			if (msh.dt.cg != ncg) {
				msh.dt.cg.Set(x_g, y_g, z_g);
				msh.dt.cg0 = clone(msh.dt.cg);
				BEM::PrintError(t_("Warning: Changing the position of the cg may have altered the inertia matrix."));
			}
			Point3D nc0(x_0, y_0, z_0);
			if (msh.dt.c0 != nc0) {
				Surface::TranslateInertia66(msh.dt.M, Point3D(x_g, y_g, z_g), msh.dt.c0, nc0);
				msh.dt.c0.Set(nc0);
			}
		}
		
		menuProcess.x_g <<= msh.dt.cg.x;
		menuProcess.y_g <<= msh.dt.cg.y;
		menuProcess.z_g <<= msh.dt.cg.z;
		
		msh.AfterLoad(Bem().rho, Bem().g, action == NONE, false, mainStiffness.opMassBuoy);
		
		UpdateLast(idx);
		
		mainView.FullRefresh(*this);
		mainViewData.OnRefresh();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
}

void MainBody::UpdateLast(int id) {
	mainTab.GetItem(mainTab.Find(mainStiffness)).Enable(mainStiffness.Load(Bem().surfs, ArrayModel_IndexsBody(listLoaded), true));
	mainTab.GetItem(mainTab.Find(mainM)).Enable(mainM.Load(Bem().surfs, ArrayModel_IndexsBody(listLoaded), false));
	mainTab.GetItem(mainTab.Find(mainStiffness2)).Enable(mainStiffness2.Load(Bem().surfs, ArrayModel_IndexsBody(listLoaded), true));
	mainTab.GetItem(mainTab.Find(mainGZ)).Enable(Bem().surfs.size() > 0);

	mainSummary.Report(Bem().surfs, id);
	
	Bem().surfs[id].dt.mesh.ClearSelPanels();
}

void MainBody::OnArchimede() {
	GuiLock __;
	
	try {
		UVector<int> idxs = ArrayModel_IndexsBody(listLoaded);
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		int idx;
		if (num == 0 && listLoaded.GetCount() == 1)
			idx = ArrayModel_IndexBody(listLoaded, 0);
		else {
		 	idx = ArrayModel_IndexBody(listLoaded);
			if (idx < 0) {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		}
		
		Body &msh = Bem().surfs[idx];
		
		if (msh.GetMass_all() == 0)
			throw Exc(t_("Set mass before fitting buoyancy"));
		if (IsNull(msh.dt.cg))
			throw Exc(t_("Set cog before fitting buoyancy"));

		WaitCursor waitcursor;
		
		double roll, pitch, dz;
		if (!msh.Archimede(Bem().rho, Bem().g, 0.05, roll, pitch, dz))
			BEM::PrintError(t_("An equilibrium position has not been found.\nApparently this device is unstable."));
		
		menuProcess.x_g <<= msh.dt.cg.x;
		menuProcess.y_g <<= msh.dt.cg.y;
		menuProcess.z_g <<= msh.dt.cg.z;
		
		msh.AfterLoad(Bem().rho, Bem().g, false, false);
		
	 	/*mainStiffness.Load(Bem().surfs, ids);
		mainView.CalcEnvelope();
		mainSummary.Report(Bem().surfs, id);*/
		UpdateLast(idx);
		
		mainView.FullRefresh(*this);
		mainViewData.OnRefresh();
		
		Ma().Status(Format(t_("Model rotated %f, %f deg. around %f, %f, %f, and translated vertically %f m to comply with displacement"), 
						ToDeg(roll), ToDeg(pitch), msh.dt.cg.x, msh.dt.cg.y, msh.dt.cg.z, dz));
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}			
}			

void MainBody::OnPCA() {
	GuiLock __;
	
	try {
		UVector<int> idxs = ArrayModel_IndexsBody(listLoaded);
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		int idx;
		if (num == 0 && listLoaded.GetCount() == 1)
			idx = ArrayModel_IndexBody(listLoaded, 0);
		else {
		 	idx = ArrayModel_IndexBody(listLoaded);
			if (idx < 0) {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		}
		
		Body &msh = Bem().surfs[idx];
		
		double yaw;
		msh.PCA(yaw);
		
		menuProcess.x_g <<= msh.dt.cg.x;
		menuProcess.y_g <<= msh.dt.cg.y;
		menuProcess.z_g <<= msh.dt.cg.z;
		
		msh.AfterLoad(Bem().rho, Bem().g, false, false);
		
		UpdateLast(idx);
		
		mainView.FullRefresh(*this);
		mainViewData.OnRefresh();
		
		Ma().Status(Format(t_("Model rotated in yaw %f deg. around %f, %f, %f to be oriented to main axis"), 
							ToDeg(yaw), msh.dt.cg.x, msh.dt.cg.y, msh.dt.cg.z));
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
	mainView.surf.Disable();
	try {
		double x = ~menuEdit.edit_x;
		double y = ~menuEdit.edit_y;
		double z = ~menuEdit.edit_z;
		double widthX = ~menuEdit.panWidthX;
		double widthY = ~menuEdit.panWidthY;
		Bem().AddFlatRectangle(x - widthX/2., y - widthY/2., z, ~menuEdit.edit_size, widthX, widthY);
		
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
	mainView.surf.Enable();
}


void MainBody::OnExtract() {
	GuiLock __;
	
	WaitCursor waitcursor;
	mainView.surf.Disable();
	try {
		UVector<int> idxs = ArrayModel_IndexsBody(listLoaded);
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		int idx;
		if (num == 0 && listLoaded.GetCount() == 1)
			idx = ArrayModel_IndexBody(listLoaded, 0);
		else {
		 	idx = ArrayModel_IndexBody(listLoaded);
			if (idx < 0) {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		}
		
		Bem().Extract(idx, ~menuEdit.cutXYZ, ~menuEdit.cutPosNeg);
		
		Body &msh = Bem().surfs[Bem().surfs.size()-1];
		msh.dt.name = t_("Extracted");
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
	mainView.surf.Enable();
}

DropCtrlDialogRevolution::DropCtrlDialogRevolution() {
	CtrlLayout(*this);
		
	list.Appending().Removing().Editing().Sorting(false).MultiSelect().Clipboard().ExtraPaste().SelectRow(false);
	list.SetToolBar();
	list.AddColumn(t_("y"), 15).Edit(y);
	list.AddColumn(t_("z"), 15).Edit(z);
	list.Add();
	
	scatter.SetLabelX("y").SetLabelY("z").SetMargin(50, 20, 20, 50).ShowAllMenus();
	scatter.AddSeries(points);
	
	SetDeactivate(false);
	
	butCancel.WhenAction = [&]() {Close();};
	y.WhenAction = z.WhenAction = list.WhenPaste = list.WhenRemoveRow = THISBACK(UpdatePlot);
}
	
void DropCtrlDialogRevolution::UpdatePlot() {
	points.Clear();
	for (int r = 0; r < list.GetRowCount(); ++r) {
		double yy = list.Get(r, 0);
		double zz = list.Get(r, 1);
		if (!IsNull(yy) && !IsNull(zz))
			points << Pointf(yy, zz);
	}
	scatter.ZoomToFit(true, true, .2);
}

DropCtrlDialogPolygon::DropCtrlDialogPolygon() {
	CtrlLayout(*this);
		
	list.Appending().Removing().Editing().Sorting(false).MultiSelect().Clipboard().ExtraPaste().SelectRow(false);
	list.SetToolBar();
	list.AddColumn(t_("x"), 15).Edit(x);
	list.AddColumn(t_("y"), 15).Edit(y);
	list.Add();

	scatter.SetLabelX("x").SetLabelY("y").SetMargin(50, 20, 20, 50).ShowAllMenus();
	scatter.AddSeries(points);
		
	SetDeactivate(false);
	
	butCancel.WhenAction = [&]() {Close();};
	
	x.WhenAction = y.WhenAction = list.WhenPaste = list.WhenRemoveRow = THISBACK(UpdatePlot);
}
	
void DropCtrlDialogPolygon::UpdatePlot() {
	points.Clear();
	for (int r = 0; r < list.GetRowCount(); ++r) {
		double xx = list.Get(r, 0);
		double yy = list.Get(r, 1);
		if (!IsNull(xx) && !IsNull(yy))
			points << Pointf(xx, yy);
	}
	scatter.ZoomToFit(true, true, .2);
}

void MainBody::OnAddRevolution() {
	GuiLock __;
	
	dialogRevolution.Close();
	
	UVector<Pointf> vals;
	for (int r = 0; r < dialogRevolution.list.GetCount(); ++r) {
		Pointf &val = vals.Add();
		val.x = ScanDouble(AsString(dialogRevolution.list.Get(r, 0)));
		if (IsNull(val.x)) {
			BEM::PrintError(Format(t_("Incorrect data in row %d, col %d"), r, 0));
			return;
		}
		val.y = ScanDouble(AsString(dialogRevolution.list.Get(r, 1)));
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
	mainView.surf.Disable();
	try {
		Bem().AddRevolution(~menuEdit.edit_x, ~menuEdit.edit_y, ~menuEdit.edit_z, ~menuEdit.edit_size, vals);
		
		Body &msh = Last(Bem().surfs);
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
	mainView.surf.Enable();
}

void MainBody::OnAddPolygonalPanel() {
	GuiLock __;
	
	dialogPolygon.Close();
	
	UVector<Pointf> vals;
	for (int r = 0; r < dialogPolygon.list.GetCount(); ++r) {
		Pointf &val = vals.Add();
		val.x = ScanDouble(AsString(dialogPolygon.list.Get(r, 0)));
		if (IsNull(val.x)) {
			BEM::PrintError(Format(t_("Incorrect data in row %d, col %d"), r, 0));
			return;
		}
		val.y = ScanDouble(AsString(dialogPolygon.list.Get(r, 1)));
		if (IsNull(val.y)) {
			BEM::PrintError(Format(t_("Incorrect data in row %d, col %d"), r, 1));
			return;
		}
	}
	if (vals.size() < 3) {
		BEM::PrintError(t_("Unsufficient value number in list"));
		return;
	}
	
	WaitCursor waitcursor;
	mainView.surf.Disable();
	try {
		Bem().AddPolygonalPanel(~menuEdit.edit_x, ~menuEdit.edit_y, ~menuEdit.edit_z, ~menuEdit.edit_size, vals, false);
		
		Body &msh = Last(Bem().surfs);
		msh.dt.name = t_("Polygon");
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
	mainView.surf.Enable();
}

void MainBody::OnExtrude() {
	GuiLock __;

	try {
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		int idx;
		if (num == 0 && listLoaded.GetCount() == 1)
			idx = ArrayModel_IndexBody(listLoaded, 0);
		else {
		 	idx = ArrayModel_IndexBody(listLoaded);
			if (idx < 0) {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		}

		WaitCursor waitcursor;
		mainView.surf.Disable();
		
		Bem().Extrude(idx, ~menuEdit.edit_cx, ~menuEdit.edit_cy, ~menuEdit.edit_cz, ~menuEdit.opClose);
	
		Body &msh = Bem().surfs[idx];
		
		msh.AfterLoad(Bem().rho, Bem().g, false, false, true, true);
		
		msh.Report(Bem().rho);
		After();
		mainViewData.OnAddedModel(mainView);
		OnOpt();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
	mainView.surf.Enable();	
}

void MainBody::OnAddWaterSurface(char c) {
	GuiLock __;

	try {
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		int idx;
		if (num == 0 && listLoaded.GetCount() == 1)
			idx = ArrayModel_IndexBody(listLoaded, 0);
		else {
		 	idx = ArrayModel_IndexBody(listLoaded);
			if (idx < 0) {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		}
			
		WaitCursor waitcursor;
		mainView.surf.Disable();
	
		Bem().AddWaterSurface(idx, c, menuProcess.meshRatio, menuProcess.opQuad);
		
		Body &nw = Last(Bem().surfs);
		AddRow(nw);
		After();
		mainViewData.OnAddedModel(mainView);
		OnOpt();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
	mainView.surf.Enable();
}

void MainBody::OnHealing(bool basic) {
	GuiLock __;
	
	try {
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		int idx;
		if (num == 0 && listLoaded.GetCount() == 1)
			idx = ArrayModel_IndexBody(listLoaded, 0);
		else {
		 	idx = ArrayModel_IndexBody(listLoaded);
			if (idx < 0) {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		}
		
		WaitCursor waitcursor;
		Progress progress(t_("Healing mesh file..."), 100); 
		mainView.surf.Disable();
		
		Bem().HealingBody(idx, basic, [&](String str, int _pos) {progress.SetText(str); progress.SetPos(_pos); return progress.Canceled();});
		
		UpdateLast(idx);
		
		mainView.FullRefresh(*this);
		mainViewData.OnRefresh();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
	mainView.surf.Enable();
}

void MainBody::OnOrientSurface() {
	GuiLock __;
	
	try {
		int num = ArrayCtrlSelectedGetCount(listLoaded);
		if (num > 1) {
			BEM::PrintError(t_("Please select just one model"));
			return;
		}
		int idx;
		if (num == 0 && listLoaded.GetCount() == 1)
			idx = ArrayModel_IndexBody(listLoaded, 0);
		else {
		 	idx = ArrayModel_IndexBody(listLoaded);
			if (idx < 0) {
				BEM::PrintError(t_("Please select a model to process"));
				return;
			}
		}
		
		WaitCursor waitcursor;
		Progress progress(t_("Orienting mesh surface..."), 100); 
		mainView.surf.Disable();
		
		Bem().OrientSurface(idx, [&](String str, int _pos) {progress.SetText(str); progress.SetPos(_pos); return progress.Canceled();});
		
		Bem().surfs[idx].AfterLoad(Bem().rho, Bem().g, false, false);
		
		UpdateLast(idx);
		
		mainView.FullRefresh(*this);
		mainViewData.OnRefresh();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
	mainView.surf.Enable();
}

void MainBody::OnImage(int axis) {
	GuiLock __;
	
	String saxis = (axis == 0) ? "X" : ((axis == 1) ? "Y" : "Z");

	try {
		int idx = ArrayModel_IndexBody(listLoaded);
		if (idx < 0) {
			BEM::PrintError(t_("Please select a model to process"));
			return;
		}
		
		WaitCursor waitcursor;
		
		Body &msh = Bem().surfs[idx];

		msh.SetMass(~menuProcess.mass);
		if (axis == 0)
			msh.dt.cg.x = -msh.dt.cg.x;
		else if (axis == 1)
			msh.dt.cg.y = -msh.dt.cg.y;
		else
			msh.dt.cg.z = -msh.dt.cg.z;
		
		msh.dt.mesh.Mirror(axis);
	
		msh.AfterLoad(Bem().rho, Bem().g, false, false);
		
		UpdateLast(idx);
		
		mainView.FullRefresh(*this);
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
	
	mainView.Disable();
	
	UVector<int> sel = ArrayCtrlSelectedGet(listLoaded);
	
	for (int r = listLoaded.GetCount()-1; r >= 0; --r) {
		if (all || Find(sel, r) >= 0) {
			int idx = ArrayModel_IndexBody(listLoaded, r);
			Bem().RemoveBody(idx);
			listLoaded.Remove(r);
			dialogDamage.RemoveId(idx);
			selected = true;
		}
	}	// Only one available => directly selected
	if (!selected && listLoaded.GetCount() == 1) {
		int idx = ArrayModel_IndexBody(listLoaded, 0);
		Bem().RemoveBody(idx);
		listLoaded.Remove(0);
		dialogDamage.RemoveId(idx);
		selected = true;		
	}	
	mainView.Enable();
		
	if (!selected) {
		BEM::PrintError(t_("No model selected"));
		return;
	}

	UVector<int> idxs = ArrayModel_IndexsBody(listLoaded);
	mainStiffness.Load(Bem().surfs, idxs, true);
	mainViewData.ReLoad(mainView);
	
	mainGZ.ClearX(false);
	
	After();
}

void MainBody::OnJoin() {
	GuiLock __;
	
	try {	
		UVector<int> idsJoin;
		for (int r = 0; r < listLoaded.GetCount(); ++r) {
			if (listLoaded.IsSelected(r)) {
				int id = ArrayModel_IndexBody(listLoaded, r);
				idsJoin << id;
			}
		}
		if (idsJoin.IsEmpty()) {
			BEM::PrintError(t_("No model joined"));
			return;
		}
		
		WaitCursor waitcursor;
		
		Bem().JoinBody(idsJoin);
	
		Body &msh = Last(Bem().surfs);
		msh.dt.name = t_("Joined");
		msh.AfterLoad(Bem().rho, Bem().g, false, true);
		
		msh.Report(Bem().rho);
		AddRow(msh);
		After();
		mainViewData.OnAddedModel(mainView);
		OnOpt();
		
	} catch (Exc e) {
		mainView.surf.Enable();
		BEM::PrintError(DeQtfLf(e));
	}
}

void MainBody::OnDuplicate() {
	GuiLock __;
	
	try {	
		int idx = ArrayModel_IndexBody(listLoaded);
		if (idx < 0) {
			BEM::PrintError(t_("Please select a model to duplicate"));
			return;
		}
		
		WaitCursor waitcursor;
		
		Bem().DuplicateBody(idx);
	
		Body &msh = Last(Bem().surfs);
		msh.dt.name = t_("Duplicated");
		msh.AfterLoad(Bem().rho, Bem().g, false, true);
		
		msh.Report(Bem().rho);
		AddRow(msh);
		After();
		mainViewData.OnAddedModel(mainView);
		OnOpt();
		
	} catch (Exc e) {
		mainView.surf.Enable();
		BEM::PrintError(DeQtfLf(e));
	}
}

void MainBody::AddRow(const Body &msh) {
	ArrayModel_Add(listLoaded, msh.GetBodyStr(), msh.dt.name, msh.dt.fileName, msh.dt.GetId(),
					optionsPlot, [&] {mainView.FullRefresh(*this);});
	dialogDamage.AddId(msh.dt.GetId(), msh.dt.name, msh.dt.fileName);
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
				
		UVector<int> idxsmesh;
		int row = -1;
		for (row = listLoaded.GetCount()-1; row >= 0; --row) {
			if (listLoaded.IsSelected(row)) 
				break;
		}	// Only one available => directly selected
		if (row < 0 && listLoaded.GetCount() == 1) 
			row = 0;
	
		if (idxsmesh.size() == 1) {
			BEM::PrintError(t_("The mesh is monolithic so it cannot be automatically split"));
			return;
		}
		int idx = ArrayModel_IndexBody(listLoaded, row);
		String fileName = Bem().surfs[idx].dt.fileName;
		String name = Bem().surfs[idx].dt.name;
		idxsmesh = Bem().SplitBody(idx, [&](String str, int _pos) {
			progress.SetText(str); 
			progress.SetPos(_pos); 
			progress.Refresh();
			return true;
		});	
		
		RemoveRow(row);
		
		for (int i = 0; i < idxsmesh.size(); ++i) {
			int idxm = idxsmesh[i];
			Body &msh = Bem().surfs[idxm];
			
			mainTab.Set(mainView);
			
			msh.Report(Bem().rho);
	
			msh.dt.fileName = fileName;
			msh.dt.name = name + S("-") + FormatInt(i+1);		
			AddRow(msh);
		}
		
		mainViewData.ReLoad(mainView);
				
		After();
	} catch (Exc e) {
		mainView.surf.Enable();
		BEM::PrintError(DeQtfLf(e));
	}
}

void MainBody::UpdateButtons() {
	int numrow = listLoaded.GetCount();
	int numsel = ArrayCtrlSelectedGetCount(listLoaded);
	menuOpen.butRemove.Enable(numrow > 0);
	menuOpen.butRemoveSelected.Enable(numsel > 0);
	menuOpen.butJoin.Enable(numsel > 1);
	menuProcess.butJoin.Enable(numsel > 1);
	menuOpen.butDuplicate.Enable(numsel == 1);
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
	
	menuEdit.butImageX.Enable(numsel == 1 || numrow == 1);
	menuEdit.butImageY.Enable(numsel == 1 || numrow == 1);
	menuEdit.butImageZ.Enable(numsel == 1 || numrow == 1);
	
	menuEdit.butExtract.Enable(numsel == 1 || numrow == 1);
	menuEdit.cutXYZ.Enable(numsel == 1 || numrow == 1);
	menuEdit.cutPosNeg.Enable(numsel == 1 || numrow == 1);
	
	menuEdit.butScale.Enable(numsel == 1 || numrow == 1);
	menuEdit.rat_x.Enable(numsel == 1 || numrow == 1);
	menuEdit.rat_y.Enable(numsel == 1 || numrow == 1);
	menuEdit.rat_z.Enable(numsel == 1 || numrow == 1);
	
	menuProcess.butSimplify.Enable(numsel == 1 || numrow == 1);
	menuProcess.butFullHealing.Enable(numsel == 1 || numrow == 1);
	menuProcess.butOrientSurface.Enable(numsel == 1 || numrow == 1);
	menuProcess.butWaterFill.Enable(numsel == 1 || numrow == 1);
	menuProcess.butWaterPlane.Enable(numsel == 1 || numrow == 1);
	menuProcess.butHull.Enable(numsel == 1 || numrow == 1);
	
	menuMove.opZArchimede.WhenAction();
}

void MainBody::After() {
	UpdateButtons();

	mainSummary.Clear();
	for (int row = 0; row < listLoaded.GetCount(); ++row) {
		int idx = ArrayModel_IndexBody(listLoaded, row);
		if (idx < 0)
			throw Exc("Unexpected problem in After()");
		mainSummary.Report(Bem().surfs, idx);
	}		

	mainTab.WhenSet();
	
	mainView.surf.Enable();
	mainView.FullRefresh(*this);
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
		menuPlot.showPoints = Null;
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
		("menuPlot_showPoints", menuPlot.showPoints)
		("menuPlot_showUnderwater", menuPlot.showUnderwater)
		("menuPlot_showWaterLevel", menuPlot.showWaterLevel)
		("menuPlot_backColor", menuPlot.backColor)
		("menuPlot_lineThickness", menuPlot.lineThickness)	
		("mainStiffness", mainStiffness)
		("menuMove_opZArchimede", menuMove.opZArchimede)
		("mainGZ", mainGZ)
		("menuStability_file", menuStability.file)
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

	const MainBody &mn = GetParentCtrl<MainBody>(this);
	Upp::Color color = ArrayModel_GetColor(mn.listLoaded, id);
	::Color textColor = Black();
	if (Grayscale(color) < 150)
		textColor = White();
														
	array.Set(row, 0, t_("File"));				array.Set(row++, col, AttrText(msh.dt.fileName).Paper(color).Ink(textColor).Bold()); 
	array.Set(row, 0, t_("Name"));				array.Set(row++, col, name + (healing ? (S(" ") + t_("(healed)")) : ""));
	array.Set(row, 0, t_("Format"));			array.Set(row++, col, msh.GetBodyStr());	
	
	array.Set(row, 0, t_("# Panels"));			array.Set(row++, col, msh.dt.mesh.panels.size());
	array.Set(row, 0, t_("# Nodes"));			array.Set(row++, col, msh.dt.mesh.nodes.size());

	array.Set(row, 0, t_("Surface [m²]"));		array.Set(row++, col, FDS(msh.dt.mesh.surface, 8, false));
	array.Set(row, 0, t_("Volume [m³] Vavg (Vx,Vy,Vz)"));		  array.Set(row++, col, AttrText(Format(t_("%s (%s, %s, %s)"), 
														FDS(msh.dt.mesh.volume,  10, false),
														FDS(msh.dt.mesh.volumex, 10, false),
														FDS(msh.dt.mesh.volumey, 10, false),
														FDS(msh.dt.mesh.volumez, 10, false))).Paper(backColorBody));
	
	array.Set(row, 0, t_("Wetted surface [m²]"));array.Set(row++, col, FDS(msh.dt.under.surface, 10, false));
	array.Set(row, 0, t_("Immersed volume [m³] Vavg (Vx,Vy,Vz)")); array.Set(row++, col, AttrText(Format(t_("%s (%s, %s, %s)"), 
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
		Vector3D cgcb = msh.dt.cg - msh.dt.cb;
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
	Force6D fcb = Surface::GetHydrostaticForceCB(msh.dt.c0, msh.dt.cb, msh.dt.under.volume, Bem().rho, Bem().g);	
	
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

	array.Set(row, 0, t_("# Triangles"));		array.Set(row++, col, !healing ? Null : msh.dt.mesh.numTriangles);
	array.Set(row, 0, t_("# BiQuads"));			array.Set(row++, col, !healing ? Null : msh.dt.mesh.numBiQuads);
	array.Set(row, 0, t_("# MonoQuads"));		array.Set(row++, col, !healing ? Null : msh.dt.mesh.numMonoQuads);
	array.Set(row, 0, t_("# Dup panels"));		array.Set(row++, col, !healing ? Null : msh.dt.mesh.numDupPan);
	array.Set(row, 0, t_("# Dup nodes"));		array.Set(row++, col, !healing ? Null : msh.dt.mesh.numDupP);
	array.Set(row, 0, t_("# Skewed pan"));		array.Set(row++, col, !healing ? Null : msh.dt.mesh.numSkewed);
	array.Set(row, 0, t_("Avg. side [m]"));		array.Set(row++, col, !healing ? Null : msh.dt.mesh.GetAvgLenSegment());
}

void MainView::ViewRefresh(MainBody &mainBody) {
	if (!IsEnabled())
		return;
	
	WithMenuBodyPlot<StaticRect> &menuPlot = mainBody.menuPlot;
	
	surf.SetShowMesh((SurfaceView::ShowMesh)(int)menuPlot.opShowMesh.GetData());
		
	Refresh();	
}

void MainView::RenderRefresh(MainBody &mainBody) {
	if (!IsEnabled())
		return;
 	
 	WithMenuBodyPlot<StaticRect> &menuPlot = mainBody.menuPlot;
 
	surf.SetBackgroundColor(~menuPlot.backColor).SetLineThickness(~menuPlot.lineThickness);
	surf.SetShowColor((SurfaceView::ShowColor)(int)menuPlot.opShowColor.GetData())
		.SetLightDir(Point3D(menuPlot.lightX, menuPlot.lightY, menuPlot.lightZ));
	surf.SetCentre(Point3D(~menuPlot.cx, ~menuPlot.cy, ~menuPlot.cz));
	
	surf.Render();
	ViewRefresh(mainBody);
}
		
void MainView::FullRefresh(MainBody &mainBody) {
	if (!IsEnabled())
		return;
	
	VolumeEnvelope env;
	for (int row = 0; row < mainBody.listLoaded.GetCount(); ++row) {
		if (ArrayModel_IsVisible(mainBody.listLoaded, row)) {
			int idx = ArrayModel_IndexBody(mainBody.listLoaded, row);
			if (idx >= 0) {
				const Body &msh = Bem().surfs[idx];
				env.MixEnvelope(msh.dt.mesh.env);
			}
		}
	}
	surf.Clear();
	
	if (env.IsNull()) {
		RenderRefresh(mainBody);
		return;
	}
	
	WithMenuBodyPlot<StaticRect> &menuPlot = mainBody.menuPlot;
	WithMenuBodyAnimation<StaticRect> &menuAnimation = mainBody.menuAnimation;
	
	if (menuPlot.showAxis && mainBody.listLoaded.GetCount() > 0) 
		surf.PaintAxis(0, 0, 0, env.LenRef()/4., -10);	
	
	if (menuPlot.showLimits && Bem().surfs.size() > 0) {
		surf.PaintGrid(Point3D(env.maxX, env.maxY, env.minZ), Point3D(env.minX, env.minY, env.minZ), WhiteGray(), SurfaceView::UNDER_ALL);
		surf.PaintCuboid(Point3D(env.maxX, env.maxY, env.maxZ), Point3D(env.minX, env.minY, env.minZ), WhiteGray(), -25);
	}
	
	double actualTime = ~menuAnimation.edTime;
	UVector<double> pos(6);
	int idtime;
	if (Bem().fast.IsEmpty())
		idtime = Null;
	else {
		idtime = Bem().fast.GetIdTime(actualTime);
		if (IsNull(idtime)) {
			if (!IsNull(actualTime) && actualTime >= 0) {
				menuAnimation.edTime <<= prevTime;
				idtime = Bem().fast.GetIdTime(prevTime);
				if (IsNull(idtime))
					idtime = 0;		// Probably reloaded shorter .out
			} else
				idtime = 0;
		} else
			prevTime = actualTime;
		
		for (int i = 0; i < 3; ++i)
			pos[i] = Bem().fast.GetVal(idtime, Bem().fast.idPos[i]);		
		for (int i = 3; i < 6; ++i)
			pos[i] = ToRad(Bem().fast.GetVal(idtime, Bem().fast.idPos[i]));		
	}
	
	for (int row = 0; row < mainBody.listLoaded.GetCount(); ++row) {
		if (ArrayModel_IsVisible(mainBody.listLoaded, row)) {
			int idx = ArrayModel_IndexBody(mainBody.listLoaded, row);
			if (idx < 0)
				throw Exc(t_("Unexpected problem in OnPaint()"));
			
			double len = env.LenRef()/10;
			bool showNormals = menuPlot.showNormals && menuPlot.showBody;
			bool showNormalsUnderwater = menuPlot.showNormals && menuPlot.showUnderwater;
			if (menuPlot.showNormals && !menuPlot.showBody && !menuPlot.showUnderwater)
				showNormals = true;
			
			const Upp::Color &color = ArrayModel_GetColor(mainBody.listLoaded, row);
			Body &msh = Bem().surfs[idx];


			if (mainBody.IsPlaying() && (menuPlot.showBody || menuPlot.showUnderwater) && msh.dt.GetCode() != Body::MOORING_MESH && !IsNull(idtime))
				msh.Move(pos, Bem().rho, Bem().g, false);
			
			if (menuPlot.showBody && !menuPlot.showUnderwater)
				surf.PaintSurface(msh.dt.mesh, color, color, 1, idx, showNormals, msh.dt.mesh.avgFacetSideLen);
			
			if (menuPlot.showUnderwater)
				surf.PaintSurface(msh.dt.under, color, color, 1, idx+10000, showNormalsUnderwater, msh.dt.mesh.avgFacetSideLen);
			
			if (menuPlot.showBody)
				surf.PaintSegments(msh.dt.mesh, color);
			
			if (menuPlot.showSkewed)
				surf.PaintSegments(msh.dt.mesh.skewed, LtRed());
			if (menuPlot.showFissure)
				surf.PaintSegments(msh.dt.mesh.segTo1panel, LtRed());
			if (menuPlot.showWaterLevel)
				surf.PaintSegments(msh.dt.under.segWaterlevel, LtBlue(), SurfaceView::OVER_ALL);
			if (menuPlot.showMultiPan)
				surf.PaintSegments(msh.dt.mesh.segTo3panel, Black());
			
			if (menuPlot.showCb && !IsNull(msh.dt.cb)) {
				surf.PaintDoubleAxis(msh.dt.cb, len, LtBlue(), SurfaceView::OVER_ALL);
				surf.PaintCube(msh.dt.cb, len/10, LtBlue());
				if (msh.cdt.damagedBodies.size() > 0) {
					Point3D cball = msh.GetCB_all();
					surf.PaintDoubleAxis(cball, len, LtBlue(), SurfaceView::OVER_ALL);
					surf.PaintCube(cball, len/5, LtBlue(), SurfaceView::OVER_ALL);
				}
			}
			if (menuPlot.showCg && !IsNull(msh.dt.cg)) {
				surf.PaintDoubleAxis(msh.dt.cg, len, Black(), SurfaceView::OVER_ALL);
				surf.PaintCube(msh.dt.cg, len/5, Black());
				if (msh.cdt.controlLoads.size() > 0) {
					Point3D cgall = msh.GetCG_all();
					surf.PaintDoubleAxis(cgall, len, Black(), SurfaceView::OVER_ALL);
					surf.PaintCube(cgall, len/3, Black(), SurfaceView::OVER_ALL);
				}
			}
			if (menuPlot.showCr) {
				surf.PaintDoubleAxis(msh.dt.c0, len*10, Cyan(), SurfaceView::OVER_ALL);
				surf.PaintCube(msh.dt.c0, len/20, Gray(), SurfaceView::OVER_ALL);
			}
			if (menuPlot.showLines && !IsNull(idtime) && msh.dt.GetCode() == Body::MOORING_MESH) {
				MeshBody::SetLineTime(msh, idtime);
				for (const auto &l : msh.dt.mesh.lines) {
					if (l.toPlot.IsEmpty()) 
						surf.PaintLines(l.lines, color, 2);
					else
						surf.PaintLines(l.toPlot, color, 2);
				}
			}
					
			if (menuPlot.showPoints) {
				for (int i = 0; i < msh.cdt.controlPointsA.size(); ++i)
					surf.PaintCube(msh.cdt.controlPointsA[i].p, len/5, msh.cdt.controlPointsA[i].p.z > 0 ? Red() : LtBlue());
				for (int i = 0; i < msh.cdt.controlPointsB.size(); ++i)
					surf.PaintCube(msh.cdt.controlPointsB[i].p, len/5, msh.cdt.controlPointsB[i].p.z > 0 ? Red() : LtBlue());
				for (int i = 0; i < msh.cdt.controlPointsC.size(); ++i)
					surf.PaintCube(msh.cdt.controlPointsC[i].p, len/5, msh.cdt.controlPointsC[i].p.z > 0 ? Red() : LtBlue());
				for (int i = 0; i < msh.cdt.controlLoads.size(); ++i)
					if (msh.cdt.controlLoads[i].loaded)
						surf.PaintCube(msh.cdt.controlLoads[i].p, len/5, msh.cdt.controlLoads[i].p.z > 0 ? LtRed() : LtBlue());
			}
			
			if (paintSelect) {
				if (menuPlot.showBody) {
					const UVector<int> &nod = msh.dt.mesh.GetSelNodes();
					for (int in = 0; in < nod.size(); ++in)
						surf.PaintCube(msh.dt.mesh.nodes[nod[in]], len/20, LtBlue());
					const UVector<int> &pan = msh.dt.mesh.GetSelPanels();
					const UVector<Point3D> &nodes = msh.dt.mesh.nodes;
					for (int ip = 0; ip < pan.size(); ++ip) {
						const Panel &panel = msh.dt.mesh.panels[pan[ip]];
						surf.PaintMesh(nodes[panel.id[0]], nodes[panel.id[1]], nodes[panel.id[2]], nodes[panel.id[3]], LtRed());//, .2);
					}
				}
				if (menuPlot.showUnderwater) {
					const UVector<int> &nod = msh.dt.under.GetSelNodes();
					for (int in = 0; in < nod.size(); ++in)
						surf.PaintCube(msh.dt.under.nodes[nod[in]], len/20, LtBlue());
					const UVector<int> &pan = msh.dt.under.GetSelPanels();
					const UVector<Point3D> &nodes = msh.dt.under.nodes;
					for (int ip = 0; ip < pan.size(); ++ip) {
						const Panel &panel = msh.dt.under.panels[pan[ip]];
						surf.PaintMesh(nodes[panel.id[0]], nodes[panel.id[1]], nodes[panel.id[2]], nodes[panel.id[3]], LtRed());//, .2);
					}
				}
			}
			if (menuPlot.showSel && ArrayModel_IsSelected(mainBody.listLoaded, row)) { 
				double minX = msh.dt.mesh.env.minX*.99; double maxX = msh.dt.mesh.env.maxX*1.01;
				double minY = msh.dt.mesh.env.minY*.99; double maxY = msh.dt.mesh.env.maxY*1.01;
				double minZ = msh.dt.mesh.env.minZ*.99; double maxZ = msh.dt.mesh.env.maxZ*1.01;
				
				surf.PaintCuboid(Point3D(maxX, maxY, maxZ), Point3D(minX, minY, minZ), color);
				surf.PaintCube(Point3D(maxX, maxY, minZ), len/10, color);
				surf.PaintCube(Point3D(maxX, minY, minZ), len/10, color);
				surf.PaintCube(Point3D(minX, maxY, minZ), len/10, color);
				surf.PaintCube(Point3D(minX, minY, minZ), len/10, color);
				surf.PaintCube(Point3D(maxX, maxY, maxZ), len/10, color);
				surf.PaintCube(Point3D(maxX, minY, maxZ), len/10, color);
				surf.PaintCube(Point3D(minX, maxY, maxZ), len/10, color);
				surf.PaintCube(Point3D(minX, minY, maxZ), len/10, color);
			}
		}
	}
	RenderRefresh(mainBody);
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
		return;
	}
	timerDrop.Kill();
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
	scatter.SetLabelX(t_("Angle [deg]")).SetLabelY(t_("GZ [m]")).SetLabelY2(t_("Heeling lever [N·m]"));
	scatter.SetTitle(t_("GZ around Y axis"));
	scatter.SetLegendFillColor(White());
	scatter.SetMargin(140, 140, 50, 80);
	
	opShowEnvelope 	<< [&]	{OnUpdate();};
	opShowPOIs 		<< [&]	{OnUpdate();};
	
	ClearX();
	
	edAngleDelta.WhenAction = [&] {edAngleTo.Enable(!(IsNull(edAngleDelta) || double(~edAngleDelta) == 0));};
}

void MainGZ::OnUpdate() {
	try {
		MainBody &mm = GetParentCtrl<MainBody>(this);
		
		idxOpened = ArrayModel_IndexBody(mm.listLoaded);
		if (idxOpened < 0) {
			BEM::PrintError(t_("Please select a model to process"));
			return;
		}
		
		ClearX();
		
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
		
		Body &msh = Bem().surfs[idxOpened];	
	
		if (IsNull(msh.dt.cg)) {
			BEM::PrintError(t_("cg is undefined"));
			return;
		}
		
		int numAngle = 1 + int((angleTo - angleFrom)/angleDelta);
		
		Progress progress(t_("GZ calculation..."), 100*numAngle); 
		
		datagz.Clear();
		dataMoment.Clear();
		datazA.Clear();
		datazB.Clear();
		datazC.Clear();
	
		moment.Clear();
	
		scatter.SetTitle(Format(t_("GZ around Y axis at (%.2f, %.2f, %.2f)"), msh.dt.c0.x, msh.dt.c0.y, msh.dt.c0.z));
		scatter.SetMargin(100, 130, 30, 70);
		
		String errors;
		int iangle = 0;
		for (double angle = double(~edAngleFrom); angle <= angleTo; angle += angleDelta, iangle++) {
			UVector<double> &dgz = datagz.Add();
			UVector<double> &dMoment = dataMoment.Add();
			UArray<UVector<double>> &zA = datazA.Add();
			UArray<UVector<double>> &zB = datazB.Add();
			UArray<UVector<double>> &zC = datazC.Add();

			UVector<double> vol, disp, wett, wplane, draft;
			UVector<Point3D> cb, cg;
			
			msh.GZ(~edFrom, ~edTo, ~edDelta, angle, Bem().rho, Bem().g, double(~edTolerance)/100., 
				[&](String, int pos)->bool {
					progress.SetPos(pos + 100*iangle);
					return !progress.Canceled();
				}, dangle, dgz, dMoment, vol, disp, wett, wplane, draft, cb, cg, errors, zA, zB, zC);
			
			Upp::Color color = ScatterDraw::GetNewColor(iangle);
			
			if (moment.IsEmpty()) {
				double mom = mm.menuStability.heelingMoment;
				if (!IsNull(mom) && mom > 0) {
					for (double a : dangle)
						moment << mom*cos(ToRad(a));
					scatter.AddSeries(dangle, moment).NoMark().Legend(t_("Moment"))
						   .Stroke(3, Cyan()).Dash(ScatterDraw::LINE_SOLID).SetDataSecondaryY();							
				}
			}
			scatter.AddSeries(dangle, dgz).NoMark().Legend(Format(t_("GZ %.1f"), angle))
				   .Stroke(2, color).Dash(ScatterDraw::LINE_SOLID);	
			scatter.AddSeries(dangle, dMoment).NoMark()./*MarkStyle<CircleMarkPlot>().*/MarkWidth(8).MarkColor(color)
				   .Legend("").NoPlot().SetDataSecondaryY();
			if (~opShowPOIs) {
				for (int j = 0; j < zA.size(); ++j)
					scatter.AddSeries(dangle, zA[j]).NoMark().Legend(Format(t_("Unpr.%s.z %.1f"), msh.cdt.controlPointsA[j].name, angle))
					   			.Stroke(1, color).Dash(ScatterDraw::LINE_DASHED);	
				for (int j = 0; j < zB.size(); ++j)
					scatter.AddSeries(dangle, zB[j]).NoMark().Legend(Format(t_("WeaT.%s.z %.1f"), msh.cdt.controlPointsB[j].name, angle))
					   			.Stroke(1, color).Dash(ScatterDraw::LINE_DOTTED);	
				for (int j = 0; j < zC.size(); ++j)
					scatter.AddSeries(dangle, zC[j]).NoMark().Legend(Format(t_("WatT.%s.z %.1f"), msh.cdt.controlPointsC[j].name, angle))
					   			.Stroke(1, color).Dash(ScatterDraw::LINE_DOTTED_FINER);	
			}
			scatter.ZoomToFit(true, true);/*
			double mx = ceil(scatter.GetYMax());
			double mn = floor(scatter.GetYMin());
			scatter.SetXYMin(Null, mn);
			scatter.SetRange(Null, mx-mn);
			int sz = int((mx-mn)/5);
			scatter.SetMajorUnits(Null, sz > 0 ? 1 : 0.5);*/
			
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
				for (int j = 0; j < zA.size(); ++j)
					array.Set(row++, col, FDS(zA[j][i], numdec));
				for (int j = 0; j < zB.size(); ++j)
					array.Set(row++, col, FDS(zB[j][i], numdec));
				for (int j = 0; j < zC.size(); ++j)
					array.Set(row++, col, FDS(zC[j][i], numdec));
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
					   .Units(t_("m"), t_("sec")).Stroke(4, Black()).Dash(ScatterDraw::LINE_SOLID);
		}
		if (!errors.IsEmpty()) 
			BEM::PrintError(DeQtfLf("Errors found in mesh:\n" + errors));
	} catch(Exc e) {
		BEM::PrintError(DeQtfLf(e));
		ClearX();
	}
}

void MainGZ::ClearX(bool all) {
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
	array.Set(row++, 0, t_("Sub. volume [m³]"));
	array.Set(row++, 0, t_("Wetted area [m²]"));
	array.Set(row++, 0, t_("Waterpl. area [m²]"));
	array.Set(row++, 0, t_("Draft [m]"));
	array.Set(row++, 0, t_("Cb_x [m]"));
	array.Set(row++, 0, t_("Cb_y [m]"));
	array.Set(row++, 0, t_("Cb_z [m]"));
	array.Set(row++, 0, t_("Cg_x [m]"));
	array.Set(row++, 0, t_("Cg_y [m]"));
	array.Set(row++, 0, t_("Cg_z [m]"));
	
	if (idxOpened >= 0 && all) {
		Body &msh = Bem().surfs[idxOpened];
		
		for (int j = 0; j < msh.cdt.controlPointsA.size(); ++j)
			array.Set(row++, 0, Format(t_("Unpr.%s.z [m]"), msh.cdt.controlPointsA[j].name));
		for (int j = 0; j < msh.cdt.controlPointsB.size(); ++j)
			array.Set(row++, 0, Format(t_("WeaT.%s.z [m]"), msh.cdt.controlPointsB[j].name));
		for (int j = 0; j < msh.cdt.controlPointsC.size(); ++j)
			array.Set(row++, 0, Format(t_("WatT.%s.z [m]"), msh.cdt.controlPointsC[j].name));
	}
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
		("opShowPOIs", opShowPOIs)
	;
}
