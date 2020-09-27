#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <TabBar/TabBar.h>

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
		LoadSelTab(Bem());
	};
	listLoaded.WhenBar = [&](Bar &menu) {
		listLoaded.StdBar(menu);
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
	menuConvert.butConvert.WhenAction = [&] {menuConvert.file.DoGo();};

	menuConvert.opt.WhenAction = [&] {OnOpt();};

	OnOpt();
	
	CtrlLayout(menuPlot);
	menuPlot.showMesh.WhenAction 	    = [&] {LoadSelTab(Bem());};
	menuPlot.showNormals.WhenAction     = [&] {LoadSelTab(Bem());};
	menuPlot.showWaterLevel.WhenAction  = [&] {LoadSelTab(Bem());};
	menuPlot.showSkewed.WhenAction      = [&] {LoadSelTab(Bem());};
	menuPlot.showFissure.WhenAction     = [&] {LoadSelTab(Bem());};
	menuPlot.showMultiPan.WhenAction    = [&] {LoadSelTab(Bem());};
	menuPlot.showAxis.WhenAction  		= [&] {mainView.gl.Refresh();};
	menuPlot.showLimits.WhenAction 		= [&] {mainView.gl.Refresh();};
	menuPlot.showCb.WhenAction  		= [&] {mainView.gl.Refresh();};
	menuPlot.showCg.WhenAction  		= [&] {mainView.gl.Refresh();};
	menuPlot.showUnderwater.WhenAction  = [&] {mainView.gl.Refresh();};
	menuPlot.butXYZ.WhenAction  		= [&] {mainView.gl.View(true, true, true);};
	menuPlot.butXoY.WhenAction  		= [&] {mainView.gl.View(true, true, false);};	
	menuPlot.butYoZ.WhenAction  		= [&] {mainView.gl.View(false, true, true);};
	menuPlot.butXoZ.WhenAction  		= [&] {mainView.gl.View(true, false, true);};
	menuPlot.butFit.WhenAction			= [&] {mainView.gl.ZoomToFit();};
	menuPlot.showMeshData.WhenAction	= [&] {mainVAll.SetButton(0);};
	
	styleRed = styleGreen = styleBlue = Button::StyleNormal();
	styleRed.textcolor[0] = styleRed.textcolor[1] = styleRed.textcolor[2] = LtRed();
	styleGreen.textcolor[0] = styleGreen.textcolor[1] = styleGreen.textcolor[2] = Green();
	styleBlue.textcolor[0] = styleBlue.textcolor[1] = styleBlue.textcolor[2] = LtBlue();
	menuPlot.butYoZ.SetStyle(styleRed);
	menuPlot.butXoZ.SetStyle(styleGreen);
	menuPlot.butXoY.SetStyle(styleBlue);
	
	OnOpt();
	
	CtrlLayout(menuProcess);
	menuProcess.cg_x <<= 0;
	menuProcess.cg_x.WhenEnter = THISBACK1(OnUpdate, NONE);
	menuProcess.cg_y <<= 0;
	menuProcess.cg_y.WhenEnter = THISBACK1(OnUpdate, NONE);
	menuProcess.cg_z <<= 0;
	menuProcess.cg_z.WhenEnter = THISBACK1(OnUpdate, NONE);
	menuProcess.mass <<= 0;
	menuProcess.mass.WhenEnter = THISBACK1(OnUpdate, NONE);
	menuProcess.butUpdateCg  <<= THISBACK1(OnUpdate, NONE);
	
	menuProcess.t_x <<= 0;
	menuProcess.t_x.WhenEnter = THISBACK1(OnUpdate, ROTATE);
	menuProcess.t_y <<= 0;
	menuProcess.t_y.WhenEnter = THISBACK1(OnUpdate, ROTATE);
	menuProcess.t_z <<= 0;
	menuProcess.t_z.WhenEnter = THISBACK1(OnUpdate, ROTATE);
	menuProcess.a_x <<= 0;
	menuProcess.a_x.WhenEnter = THISBACK1(OnUpdate, ROTATE);
	menuProcess.a_y <<= 0;
	menuProcess.a_y.WhenEnter = THISBACK1(OnUpdate, ROTATE);
	menuProcess.a_z <<= 0;
	menuProcess.a_z.WhenEnter = THISBACK1(OnUpdate, ROTATE);
	menuProcess.c_x <<= 0;
	menuProcess.c_x.WhenEnter = THISBACK1(OnUpdate, ROTATE);
	menuProcess.c_y <<= 0;
	menuProcess.c_y.WhenEnter = THISBACK1(OnUpdate, ROTATE);
	menuProcess.c_z <<= 0;
	menuProcess.c_z.WhenEnter = THISBACK1(OnUpdate, ROTATE);
	menuProcess.butUpdatePos <<= THISBACK1(OnUpdate, MOVE);
	menuProcess.butUpdateAng <<= THISBACK1(OnUpdate, ROTATE);
	menuProcess.butImageX <<= THISBACK1(OnImage, 0);
	menuProcess.butImageY <<= THISBACK1(OnImage, 1);
	menuProcess.butImageZ <<= THISBACK1(OnImage, 2);
	
	menuProcess.butBasicHealing <<= THISBACK1(OnHealing, true);
	menuProcess.butFullHealing <<= THISBACK1(OnHealing, false);
	menuProcess.butOrientSurface <<= THISBACK(OnOrientSurface);
	
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
	
	menuEdit.butPolynomial <<= THISBACK(OnAddPolygonalPanel);
		
	menuTab.Add(menuOpen.SizePos(),    	t_("Load"));
	menuTab.Add(menuPlot.SizePos(),    	t_("Plot")).Disable();
	menuTab.Add(menuProcess.SizePos(), 	t_("Process")).Disable();
	menuTab.Add(menuEdit.SizePos(), 	t_("Edit"));
	menuTab.Add(menuConvert.SizePos(), 	t_("Save as")).Disable();
		
	mainViewData.Init();
	mainVAll.Horz(mainView, mainViewData);
	mainVAll.SetPositions(6000, 9970).SetInitialPositionId(1).SetButtonNumber(1);
	mainVAll.WhenAction = [&] {mainView.SetPaintSelect(mainVAll.GetPos() < 9950);};
	mainTab.Add(mainVAll.SizePos(), t_("View"));
	mainView.Init();
	
	mainSummary.Init();
	mainTab.Add(mainSummary.SizePos(), t_("Summary"));

	mainStiffness.Init();
	mainTab.Add(mainStiffness.SizePos(), t_("K Stiffness Matrix"));
			
	mainTab.WhenSet = [&] {
		LOGTAB(mainTab);
		Upp::Vector<int> ids = ArrayModel_IdsMesh(listLoaded);
		bool plot = true, convertProcess = true;
		if (Bem().surfs.IsEmpty()) 
			plot = convertProcess = false;
		else if (mainTab.IsAt(mainVAll)) 
			;
		else if (mainTab.IsAt(mainStiffness)) {
			plot = false;
			mainStiffness.Load(Bem().surfs, ids);
		} else 
			plot = false;
		
		TabCtrl::Item& tabMenuPlot = menuTab.GetItem(menuTab.Find(menuPlot));
		tabMenuPlot.Enable(plot);
		TabCtrl::Item& tabMenuConvert = menuTab.GetItem(menuTab.Find(menuConvert));
		tabMenuConvert.Enable(convertProcess);
		TabCtrl::Item& tabMenuProcess = menuTab.GetItem(menuTab.Find(menuProcess));
		tabMenuProcess.Enable(convertProcess);
		if (plot) {
			tabMenuPlot.Text(t_("Plot"));
			menuTab.Set(menuPlot);
		} else {
			tabMenuPlot.Text("");
			menuTab.Set(menuOpen);
		}
		if (convertProcess) {
			tabMenuConvert.Text(t_("Save as"));
			tabMenuProcess.Text(t_("Process"));
		} else {
			tabMenuConvert.Text("");
			tabMenuProcess.Text("");
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
	
	MeshData &data = Bem().surfs[id];
	menuProcess.cg_x <<= data.cg.x;
	menuProcess.cg_y <<= data.cg.y;
	menuProcess.cg_z <<= data.cg.z;
	menuProcess.mass <<= data.mass;
}

void MainMesh::OnArraySel() {
	OnMenuConvertArraySel();
	OnMenuProcessArraySel();
}

void MainMesh::OnMenuConvertArraySel() {
	int id = ArrayModel_IdMesh(listLoaded);
	if (id < 0)
		return;
	
	String file = ~menuConvert.file;
	String folder = GetFileFolder(file);
	String ext = GetFileExt(file);
	String fileName = GetFileTitle(ArrayModel_GetFileName(listLoaded));
	file = AppendFileName(folder, fileName + ext);
	
	menuConvert.file <<= file;
	menuConvert.symX <<= Bem().surfs[id].IsSymmetricX();
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
			
	if (!ret || IsNull(menuConvert.opt)) 
		menuConvert.opt = 0;
	if (!ret || IsNull(menuConvert.optMeshType)) 
		menuConvert.optMeshType = 0;
}

void MainMesh::LoadSelTab(BEMData &bem) {
	const Upp::Vector<int> &ids = ArrayModel_IdsMesh(listLoaded);
	if (mainTab.Get() == mainTab.Find(mainStiffness))
		mainStiffness.Load(bem.surfs, ids);
	else 
		mainView.gl.Refresh();
}

void MainMesh::OnOpt() {
	menuOpen.file.ClearTypes(); 

	const String meshFiles = ".gdf .dat .stl";
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
	switch (menuConvert.opt) {
	case 0:	menuConvert.file <<= ForceExtSafe(~menuConvert.file, ".gdf"); 	
			menuConvert.file.Type("Wamit .gdf file", "*.gdf");
			break;
	case 1:	menuConvert.file <<= ForceExtSafe(~menuConvert.file, ".dat"); 	
			menuConvert.file.Type("Nemoh .dat file", "*.dat");
			break;
	case 2:	menuConvert.file <<= ForceExtSafe(~menuConvert.file, "."); 	
			menuConvert.file.Type("Nemoh pre mesh file", "*.");
			break;
	case 3:	menuConvert.file <<= ForceExtSafe(~menuConvert.file, ".stl"); 	
			menuConvert.file.Type("STL binary .stl file", "*.stl");
			break;
	case 4:	menuConvert.file <<= ForceExtSafe(~menuConvert.file, ".stl"); 	
			menuConvert.file.Type("STL text .stl file", "*.stl");
			break;
	default:menuConvert.file.Type("All converted files", "*.gdf *.dat *.stl");
			break;
	}
	String extConvmesh = ToLower(GetFileExt(menuConvert.file.GetData().ToString()));
	if (extConvmesh.IsEmpty())
		menuConvert.file.ActiveType(0);
	else if (String(".gdf").Find(extConvmesh) >= 0)
		menuConvert.file.ActiveType(0);
	else if (String(".dat").Find(extConvmesh) >= 0)
		menuConvert.file.ActiveType(1);
	else if (String(".stl").Find(extConvmesh) >= 0)
		menuConvert.file.ActiveType(2);
}

void MainMesh::AfterAdd(String file) {
	int id = Bem().surfs.GetCount() - 1;
	MeshData &surf = Bem().surfs[id];
		
	mainTab.Set(mainView);
	
	surf.Report(Bem().rho);
	
	AddRow(surf);

	mainView.CalcEnvelope();
	
	mainView.gl.Enable();
	mainView.gl.ZoomToFit();
	
	After();		
}

bool MainMesh::OnLoad() {
	GuiLock __;
	
	String file = ~menuOpen.file;
		
	try {
		Progress progress(t_("Loading mesh file..."), 100); 
		
		Upp::Vector<int> ids = ArrayModel_IdsMesh(listLoaded);
		for (int i = 0; i < ids.GetCount(); ++i) {
			if (Bem().surfs[ids[i]].fileName == file) {
				if (!PromptYesNo(t_("Model is already loaded") + S("&") + t_("Do you wish to open it anyway?")))
					return false;
				break;
			}
		}
		mainView.gl.Disable();
		
		WaitCursor waitcursor;

		Bem().LoadMesh(file, [&](String str, int _pos) {progress.SetText(str); progress.SetPos(_pos);}, ~menuOpen.opClean, false);
		
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
		
		MeshData::MESH_FMT type;	
		switch (menuConvert.opt) {
		case 0:	type = MeshData::WAMIT_GDF;	break;
		case 1:	type = MeshData::NEMOH_DAT;	break;
		case 2:	type = MeshData::NEMOH_PRE;	break;
		case 3:	type = MeshData::STL_BIN;	break;
		case 4:	type = MeshData::STL_TXT;	break;
		case 5:	type = MeshData::UNKNOWN;	break;
		default: throw Exc(t_("Unknown type in OnConvert()"));
		}
		
		WaitCursor waitcursor;
		Progress progress(t_("Saving mesh file..."), 100); 
		
		Bem().surfs[id].SaveAs(~menuConvert.file, type, Bem().g, 
							   static_cast<MeshData::MESH_TYPE>(int(~menuConvert.optMeshType)),
							   ~menuConvert.symX, ~menuConvert.symY);	
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
		return false;
	}
	return true;
}

void MainMesh::OnUpdate(Action action) {
	GuiLock __;
	
	try {
		Upp::Vector<int> ids = ArrayModel_IdsMesh(listLoaded);
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
				
		MeshData &data = Bem().surfs[id];

		double mass = ~menuProcess.mass;
		double cg_x = ~menuProcess.cg_x;
		double cg_y = ~menuProcess.cg_y;
		double cg_z = ~menuProcess.cg_z;
		double t_x = ~menuProcess.t_x;
		double t_y = ~menuProcess.t_y;
		double t_z = ~menuProcess.t_z;
		double a_x = ~menuProcess.a_x;
		double a_y = ~menuProcess.a_y;
		double a_z = ~menuProcess.a_z;
		double c_x = ~menuProcess.c_x;
		double c_y = ~menuProcess.c_y;
		double c_z = ~menuProcess.c_z;

		if (action == NONE && (IsNull(mass) || IsNull(cg_x) || IsNull(cg_y) || IsNull(cg_z))) {
			Exclamation(t_("Please fill CG data"));
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
		if (action == ROTATE && (IsNull(c_x) || IsNull(c_y) || IsNull(c_z))) {
			Exclamation(t_("Please fill center of rotation data"));
			return;
		}
		
		WaitCursor wait;

		if (action == MOVE) {
			data.cg.Translate(t_x, t_y, t_z);
			data.mesh.Translate(t_x, t_y, t_z);
		} else if (action == ROTATE) {
			data.cg.Rotate(a_x, a_y, a_z, c_x, c_y, c_z);
			data.mesh.Rotate(a_x, a_y, a_z, c_x, c_y, c_z);
		} else if (action == NONE) {
			data.mass = mass;
			data.cg.Set(cg_x, cg_y, cg_z);
		}
		
		data.AfterLoad(Bem().rho, Bem().g, action == NONE);
		
	 	mainStiffness.Load(Bem().surfs, ids);
		mainView.CalcEnvelope();
		mainSummary.Report(Bem().surfs, id);
		
		mainView.gl.Refresh();
		mainViewData.OnRefresh();
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
		
		MeshData &surf = Bem().surfs[Bem().surfs.GetCount()-1];
		surf.name = t_("Panel");
		surf.fileName =  "";
		
		surf.AfterLoad(Bem().rho, Bem().g, false);
		
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
	
	Upp::Vector<Pointf> vals;
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
	if (vals.GetCount() < 2) {
		Exclamation(t_("Unsufficient value number in list"));
		return;
	}
	
	WaitCursor waitcursor;
	mainView.gl.Disable();
	try {
		Bem().AddRevolution(~menuEdit.edit_x, ~menuEdit.edit_y, ~menuEdit.edit_z, ~menuEdit.edit_size, vals);
		
		MeshData &surf = Bem().surfs[Bem().surfs.GetCount()-1];
		surf.name = t_("Revolution");
		surf.fileName =  "";
		
		surf.AfterLoad(Bem().rho, Bem().g, false);
		
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
	
	Upp::Vector<Pointf> vals;
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
	if (vals.GetCount() < 3) {
		Exclamation(t_("Unsufficient value number in list"));
		return;
	}
	
	WaitCursor waitcursor;
	mainView.gl.Disable();
	try {
		Bem().AddPolygonalPanel(~menuEdit.edit_x, ~menuEdit.edit_y, ~menuEdit.edit_z, ~menuEdit.edit_size, vals);
		
		MeshData &surf = Bem().surfs[Bem().surfs.GetCount()-1];
		surf.name = t_("Polynomial");
		surf.fileName =  "";
		
		surf.AfterLoad(Bem().rho, Bem().g, false);
		
		surf.Report(Bem().rho);
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
		
		Bem().HealingMesh(id, basic, [&](String str, int _pos) {progress.SetText(str); progress.SetPos(_pos);});
		
		Bem().surfs[id].AfterLoad(Bem().rho, Bem().g, false);
		
		Upp::Vector<int> ids = ArrayModel_IdsMesh(listLoaded);
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
		
		Bem().OrientSurface(id, [&](String str, int _pos) {progress.SetText(str); progress.SetPos(_pos);});
		
		Bem().surfs[id].AfterLoad(Bem().rho, Bem().g, false);
		
		Upp::Vector<int> ids = ArrayModel_IdsMesh(listLoaded);
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
		Upp::Vector<int> ids = ArrayModel_IdsMesh(listLoaded);
		int id = ArrayModel_IdMesh(listLoaded);
		if (id < 0) {
			Exclamation(t_("Please select a model to process"));
			return;
		}
		
		WaitCursor waitcursor;
		
		MeshData &data = Bem().surfs[id];

		data.mass = ~menuProcess.mass;
		if (axis == 0)
			data.cg.x = -data.cg.x;
		else if (axis == 1)
			data.cg.y = -data.cg.y;
		else
			data.cg.z = -data.cg.z;
		
		data.mesh.Image(axis);
	
		data.AfterLoad(Bem().rho, Bem().g, false);
		
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
	
	mainSummary.Clear();
	listLoaded.Clear();
	menuOpen.butRemove.Disable();
	menuOpen.butRemoveSelected.Disable();
	menuOpen.butJoin.Disable();
	menuOpen.butSplit.Disable();
	optionsPlot.Clear();
	
	mainStiffness.Clear();
	mainViewData.Clear();
	
	Bem().surfs.Clear();
	mainView.env.Reset();
	
	mainTab.WhenSet();
	mainView.gl.Refresh();
}

void MainMesh::OnRemoveSelected(bool all) {	
	bool selected = false;
	
	Upp::Vector<int> sel = ArrayCtrlSelectedGet(listLoaded);
	
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

	Upp::Vector<int> ids = ArrayModel_IdsMesh(listLoaded);
	mainStiffness.Load(Bem().surfs, ids);
	mainViewData.ReLoad(mainView);
	
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
	
		Upp::Vector<int> ids = ArrayModel_IdsMesh(listLoaded);
		mainStiffness.Load(Bem().surfs, ids);
		mainViewData.ReLoad(mainView);
		
		After();
	} catch (Exc e) {
		mainView.gl.Enable();
		Exclamation(DeQtfLf(e));
	}
}

void MainMesh::AddRow(const MeshData &surf) {
	ArrayModel_Add(listLoaded, surf.GetCodeStr(), surf.name, surf.fileName, surf.GetId(),
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
				
		Upp::Vector<int> idsmesh;
		int row = -1;
		for (row = listLoaded.GetCount()-1; row >= 0; --row) {
			if (listLoaded.IsSelected(row)) 
				break;
		}	// Only one available => directly selected
		if (row < 0 && listLoaded.GetCount() == 1) 
			row = 0;
	
		if (idsmesh.GetCount() == 1) {
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
		});	
		
		RemoveRow(row);
		
		for (int i = 0; i < idsmesh.GetCount(); ++i) {
			int id = idsmesh[i];
			MeshData &surf = Bem().surfs[id];
			
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
	menuProcess.butUpdatePos.Enable(numsel == 1 || numrow == 1);
	menuProcess.butUpdateAng.Enable(numsel == 1 || numrow == 1);
	menuProcess.butImageX.Enable(numsel == 1 || numrow == 1);
	menuProcess.butImageY.Enable(numsel == 1 || numrow == 1);
	menuProcess.butImageZ.Enable(numsel == 1 || numrow == 1);
	menuProcess.butBasicHealing.Enable(numsel == 1 || numrow == 1);
	menuProcess.butFullHealing.Enable(numsel == 1 || numrow == 1);
	menuProcess.butOrientSurface.Enable(numsel == 1 || numrow == 1);
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
		menuConvert.opt = Null;
		menuConvert.optMeshType = Null;
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
		("menuPlot_showUnderwater", menuPlot.showUnderwater)
		("menuPlot_showWaterLevel", menuPlot.showWaterLevel)	
	;
}

void MainSummaryMesh::Report(const Upp::Array<MeshData> &surfs, int id) {
	const MeshData &data = surfs[id];
	String name = data.name;
	
	if (array.GetColumnCount() == 0)
		array.AddColumn("Param");
	if (id >= array.GetColumnCount()-1)
		array.AddColumn(Format("#%d %s", id+1, name));
	int row = 0;
	int col = id + 1;
	
	bool healing = data.mesh.healing;
	
	array.Set(row, 0, t_("File"));				array.Set(row++, col, data.fileName);
	array.Set(row, 0, t_("Name"));				array.Set(row++, col, name + (healing ? (S(" ") + t_("(healed)")) : ""));
	array.Set(row, 0, t_("Format"));			array.Set(row++, col, data.GetCodeStr());	
	
	array.Set(row, 0, t_("# Panels"));			array.Set(row++, col, data.mesh.panels.GetCount());
	array.Set(row, 0, t_("# Nodes"));			array.Set(row++, col, data.mesh.nodes.GetCount());

	array.Set(row, 0, t_("Surface [m2]"));		array.Set(row++, col, FormatDouble(data.mesh.surface, 6, FD_EXP));
	array.Set(row, 0, t_("Volume [m3]"));		array.Set(row++, col, Format(t_("%s (%s, %s, %s)"), 
														FormatDouble(data.mesh.volume,  6, FD_EXP),
														FormatDouble(data.mesh.volumex, 6, FD_EXP),
														FormatDouble(data.mesh.volumey, 6, FD_EXP),
														FormatDouble(data.mesh.volumez, 6, FD_EXP)));
	
	array.Set(row, 0, t_("Immersed surface [m2]"));array.Set(row++, col, FormatDouble(data.under.surface, 6, FD_EXP));
	array.Set(row, 0, t_("Immersed volume [m3]")); array.Set(row++, col, Format(t_("%s (%s, %s, %s)"), 
														FormatDouble(data.under.volume,  6, FD_EXP),
														FormatDouble(data.under.volumex, 6, FD_EXP),
														FormatDouble(data.under.volumey, 6, FD_EXP),
														FormatDouble(data.under.volumez, 6, FD_EXP)));
	array.Set(row, 0, t_("Displacement [Kg]")); array.Set(row++, col, FormatDouble(data.under.volume*Bem().rho, 9, FD_EXP));
	array.Set(row, 0, t_("Cg [m]"));			array.Set(row++, col, Format(t_("%s, %s, %s"),
														FormatDouble(data.cg.x, 3, FD_EXP),			
														FormatDouble(data.cg.y, 3, FD_EXP),
														FormatDouble(data.cg.z, 3, FD_EXP)));
	array.Set(row, 0, t_("Cb [m]"));			array.Set(row++, col, Format(t_("%s, %s, %s"),
														FormatDouble(data.cb.x, 3, FD_EXP),			
														FormatDouble(data.cb.y, 3, FD_EXP),
														FormatDouble(data.cb.z, 3, FD_EXP)));
	array.Set(row, 0, t_("Water Plane Area [m2]"));	array.Set(row++, col, FormatDouble(data.waterPlaneArea, 6, FD_EXP));
	
	array.Set(row, 0, t_("Dimensions [m]"));	array.Set(row++, col, Format(t_("From (%s, %s, %s) to (%s, %s, %s)"),
														FormatDouble(data.mesh.env.minX, 3, FD_EXP),
														FormatDouble(data.mesh.env.minY, 3, FD_EXP),
														FormatDouble(data.mesh.env.minZ, 3, FD_EXP),
														FormatDouble(data.mesh.env.maxX, 3, FD_EXP),
														FormatDouble(data.mesh.env.maxY, 3, FD_EXP),
														FormatDouble(data.mesh.env.maxZ, 3, FD_EXP)));

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

	array.Set(row, 0, t_("# Segments"));		array.Set(row++, col, !healing ? Null : data.mesh.segments.GetCount());
	array.Set(row, 0, t_("# Seg Water Plane"));	array.Set(row++, col, !healing ? Null : data.mesh.segWaterlevel.GetCount());
	array.Set(row, 0, t_("# Seg leak"));		array.Set(row++, col, !healing ? Null : data.mesh.segTo1panel.GetCount());
	array.Set(row, 0, t_("# Seg 3 panels"));	array.Set(row++, col, !healing ? Null : data.mesh.segTo3panel.GetCount());

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
	if (~GetMenuPlot().showAxis) 
		gl.PaintAxis(0, 0, 0, env.LenRef()/4.);	
	
	if (~GetMenuPlot().showLimits) 
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
			const MeshData &mesh = Bem().surfs[id];
			
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
				gl.PaintCube(mesh.cb, len/10, LtGray());
			}
			if (~GetMenuPlot().showCg) {
				gl.PaintDoubleAxis(mesh.cg, len, Black());
				gl.PaintCube(mesh.cg, len/10, LtGray());
			}
			if (paintSelect) {
				if (~GetMenuPlot().showMesh) {
					const Upp::Vector<int> &nod = mesh.mesh.GetSelNodes();
					for (int in = 0; in < nod.GetCount(); ++in)
						gl.PaintCube(mesh.mesh.nodes[nod[in]], len/20, LtBlue());
					const Upp::Vector<int> &pan = mesh.mesh.GetSelPanels();
					const Upp::Vector<Point3D> &nodes = mesh.mesh.nodes;
					for (int ip = 0; ip < pan.GetCount(); ++ip) {
						const Panel &panel = mesh.mesh.panels[pan[ip]];
						gl.PaintQuad(nodes[panel.id[0]], nodes[panel.id[1]], nodes[panel.id[2]], nodes[panel.id[3]], LtRed(), .2);
					}
				}
				if (~GetMenuPlot().showUnderwater) {
					const Upp::Vector<int> &nod = mesh.under.GetSelNodes();
					for (int in = 0; in < nod.GetCount(); ++in)
						gl.PaintCube(mesh.under.nodes[nod[in]], len/20, LtBlue());
					const Upp::Vector<int> &pan = mesh.under.GetSelPanels();
					const Upp::Vector<Point3D> &nodes = mesh.under.nodes;
					for (int ip = 0; ip < pan.GetCount(); ++ip) {
						const Panel &panel = mesh.under.panels[pan[ip]];
						gl.PaintQuad(nodes[panel.id[0]], nodes[panel.id[1]], nodes[panel.id[2]], nodes[panel.id[3]], LtRed(), .2);
					}
				}
			}
			if (ArrayModel_IsSelected(GetMain().listLoaded, row)) { 
				double minX = mesh.mesh.env.minX; double maxX = mesh.mesh.env.maxX;
				double minY = mesh.mesh.env.minY; double maxY = mesh.mesh.env.maxY;
				double minZ = mesh.mesh.env.minZ; double maxZ = mesh.mesh.env.maxZ;
				
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
	//for (int i = 0; i < Bem().surfs.GetCount(); ++i)
	for (int row = 0; row < GetMain().listLoaded.GetCount(); ++row) {
		int id = ArrayModel_IdMesh(GetMain().listLoaded, row);
		if (id < 0)
			throw Exc("Unexpected problem in CalcEnvelope()");
		env.MixEnvelope(Bem().surfs[id].mesh.env);
	}
}

void MainMesh::LoadDragDrop(const Upp::Vector<String> &files) {
	bool followWithErrors = false;
	for (int i = 0; i < files.GetCount(); ++i) {
		menuOpen.file <<= files[i];
		Status(Format(t_("Loading '%s'"), files[i]));
		if (!OnLoad() && !followWithErrors && files.GetCount() - i > 1) {
			if (!PromptYesNo(Format(t_("Do you wish to load the pending %d files?"), files.GetCount() - i - 1)))
				return;
			followWithErrors = true;
		}
		ProcessEvents();
	}
}
	
void MainMesh::DragAndDrop(Point , PasteClip& d) {
	if (IsDragAndDropSource())
		return;
	if (AcceptFiles(d)) {
		Upp::Vector<String> files = GetFiles(d);
		LoadDragDrop(files);
	}
}

bool MainMesh::Key(dword key, int ) {
	if (key == K_CTRL_V) {
		Upp::Vector<String> files = GetFiles(Ctrl::Clipboard());
		LoadDragDrop(files);
		return true;
	}
	return false;
}

void MainViewData::Init() {
	Add(tab.SizePos());	
}

void MainViewData::OnAddedModel(MainView &mainView) {
	int id = Bem().surfs.GetCount()-1;
	
	MainViewDataEach &model = models.Add();
	model.Init(Bem().surfs[id], mainView);
	tab.Add(model.SizePos(), Bem().surfs[id].name);
}

void MainViewData::OnRefresh() {
	for (int i = 0; i < models.GetCount(); ++i)
		models[i].OnRefresh();
}

void MainViewData::Clear() {
	models.Clear();
	tab.Reset();
}

void MainViewData::ReLoad(MainView &mainView) {
	Clear();
	
	for (int i = 0; i < Bem().surfs.GetCount(); ++i)  {
		MainViewDataEach &model = models.Add();
		model.Init(Bem().surfs[i], mainView);
		tab.Add(model.SizePos(), Bem().surfs[i].name);		
	}
}

void MainViewDataEach::DataSourceFacets::Init(MeshData &_mesh, int _col, bool _all) {
	pmesh = &_mesh;	
	col = _col; 
	all = _all;
}

Value MainViewDataEach::DataSourceFacets::Format(const Value& q) const {
	ASSERT(pmesh);
	int iq = q;
	if (col < 0)
		return iq + 1;
	else {
		if (all) {
			if (iq >= pmesh->mesh.panels.GetCount())
				return Null;
			if (col == 3 && pmesh->mesh.panels[iq].IsTriangle())
				return "-";
			else
				return pmesh->mesh.panels[iq].id[col]+1;
		} else {
			if (iq >= pmesh->under.panels.GetCount())
				return Null;
			if (col == 3 && pmesh->under.panels[iq].IsTriangle())
				return "-";
			else
				return pmesh->under.panels[iq].id[col]+1;
		}
	}
}

void MainViewDataEach::DataSourceNodes::Init(MeshData &_mesh, int _xyz, int _origMovedUnder) {
	pmesh = &_mesh;	
	xyz = _xyz;
	origMovedUnder = _origMovedUnder;
}

Value MainViewDataEach::DataSourceNodes::Format(const Value& q) const {
	ASSERT(pmesh);
	int iq = q;
	if (origMovedUnder == 0 && pmesh->mesh.nodes.GetCount() <= iq)
		return Null;
	if (origMovedUnder == 1 && pmesh->under.nodes.GetCount() <= iq)
		return Null;
	
	const Point3D &p = origMovedUnder == 0 ? pmesh->mesh.nodes[iq] : pmesh->under.nodes[iq];
	if (xyz == -1)
		return iq + 1;
	else if (xyz == 0)
		return p.x;
	else if (xyz == 1)
		return p.y;
	else
		return p.z;
}

void MainViewDataEach::UpdateStatus(bool under) {
	MainMesh &mainMesh = GetDefinedParent<MainMesh>(this);
	
	bool show;
	if (!under) {
		show = mainMesh.GetShowMesh();
		selectedPanels = ArrayCtrlSelectedGet(arrayFacetsAll2.array);
		selectedNodes = ArrayCtrlSelectedGet(arrayNodesMoved.array);
	} else {
		show = mainMesh.GetShowUnderwater();
		selectedPanels = ArrayCtrlSelectedGet(arrayFacetsUnder.array);
		selectedNodes = ArrayCtrlSelectedGet(arrayNodesUnder.array);
	}
	int numPanels = selectedPanels.GetCount();
	int numNodes  = selectedNodes.GetCount();
	String strPanels = numPanels > 0 ? FormatInt(numPanels) : S(t_("no"));
	String strNodes  = numNodes > 0  ? FormatInt(numNodes)  : S(t_("no"));
	
	if (numPanels + numNodes > 0) {
		if (!show)
			status.Set(Format(t_("%s is hidden in Plot menu so selection will not be shown"), 
						under ? t_("Underwater mesh") : t_("Mesh")));
		else
			status.Set(Format(t_("Selected %s panels and %s nodes"), strPanels, strNodes));
	} else
		status.Set("");
}

void MainViewDataEach::Init(MeshData &_mesh, MainView &mainView) {
	CtrlLayout(arrayFacetsAll2);
	CtrlLayout(arrayFacetsUnder);
	CtrlLayout(arrayNodesMoved);
	CtrlLayout(arrayNodesUnder);
	
	arrayFacetsAll2.title.SetText(t_("Facet node ids"));
	arrayFacetsAll2.array.MultiSelect().SetLineCy(EditField::GetStdHeight()).HeaderObject();
	arrayFacetsAll2.array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, arrayFacetsAll2.array);};
	dataSourceFacetsAll.SetCount(5);
	dataSourceFacetsAll[0].Init(_mesh, -1, true);
	arrayFacetsAll2.array.AddRowNumColumn(t_("#panel"), 60).SetConvert(dataSourceFacetsAll[0]);
	for (int c = 0; c < 4; ++c) {
		dataSourceFacetsAll[c+1].Init(_mesh, c, true);
		arrayFacetsAll2.array.AddRowNumColumn(Format(t_("#%d"), c+1), 60).SetConvert(dataSourceFacetsAll[c+1]);
	}
	arrayFacetsAll2.array.WhenSel = [&] {
		UpdateStatus(false);
		_mesh.mesh.SelPanels(selectedPanels);	
		arrayFacetsUnder.array.ClearSelection();	
		mainView.gl.Refresh();		
		lastSel = 0;
	};
		
	arrayFacetsUnder.title.SetText(t_("Facet node ids"));
	arrayFacetsUnder.array.MultiSelect().SetLineCy(EditField::GetStdHeight()).HeaderObject();
	arrayFacetsUnder.array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, arrayFacetsUnder.array);};
	dataSourceFacetsUnder.SetCount(5);
	dataSourceFacetsUnder[0].Init(_mesh, -1, true);
	arrayFacetsUnder.array.AddRowNumColumn(t_("#panel"), 60).SetConvert(dataSourceFacetsUnder[0]);
	for (int c = 0; c < 4; ++c) {
		dataSourceFacetsUnder[c+1].Init(_mesh, c, false);
		arrayFacetsUnder.array.AddRowNumColumn(Format(t_("#%d"), c+1), 60).SetConvert(dataSourceFacetsUnder[c+1]);
	}
	arrayFacetsUnder.array.WhenSel = [&] {
		UpdateStatus(true);
		_mesh.under.SelPanels(selectedPanels);
		arrayFacetsAll2.array.ClearSelection();	
		mainView.gl.Refresh();	
		lastSel = 1;	
	};
	
	const char *xyz[] = {"x", "y", "z"};

	arrayNodesMoved.title.SetText(t_("Node coordinates"));
	arrayNodesMoved.array.MultiSelect().SetLineCy(EditField::GetStdHeight()).HeaderObject();
	arrayNodesMoved.array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, arrayNodesMoved.array);};
	dataSourceNodesMoved.SetCount(4);
	dataSourceNodesMoved[0].Init(_mesh, -1, 1);
	arrayNodesMoved.array.AddRowNumColumn(t_("#node"), 60).SetConvert(dataSourceNodesMoved[0]);
	for (int c = 0; c < 3; ++c) {
		dataSourceNodesMoved[c+1].Init(_mesh, c, 1);
		arrayNodesMoved.array.AddRowNumColumn(Format(t_("%s"), xyz[c]), 80).SetConvert(dataSourceNodesMoved[c+1]);
	}
	arrayNodesMoved.array.WhenSel = [&] {
		UpdateStatus(false);
		_mesh.mesh.SelNodes(selectedNodes);
		arrayNodesUnder.array.ClearSelection();
		mainView.gl.Refresh();	
		lastSel = 2;
	};
		
	arrayNodesUnder.title.SetText(t_("Node coordinates"));
	arrayNodesUnder.array.MultiSelect().SetLineCy(EditField::GetStdHeight()).HeaderObject();
	arrayNodesUnder.array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, arrayNodesUnder.array);};
	dataSourceNodesMoved.SetCount(4);
	dataSourceNodesMoved[0].Init(_mesh, -1, 2);
	arrayNodesUnder.array.AddRowNumColumn(t_("#node"), 60).SetConvert(dataSourceNodesMoved[0]);
	for (int c = 0; c < 3; ++c) {
		dataSourceNodesMoved[c+1].Init(_mesh, c, 2);
		arrayNodesUnder.array.AddRowNumColumn(Format(t_("%s"), xyz[c]), 80).SetConvert(dataSourceNodesMoved[c+1]);
	}
	arrayNodesUnder.array.WhenSel = [&] {
		UpdateStatus(true);
		_mesh.under.SelNodes(selectedNodes);	
		arrayNodesMoved.array.ClearSelection();
		mainView.gl.Refresh();	
		lastSel = 3;
	};
	
	moved.Horz(arrayFacetsAll2.SizePos(), arrayNodesMoved.SizePos());
	movedUnder.Horz(arrayFacetsUnder.SizePos(), arrayNodesUnder.SizePos());	
	  					 
	tab.Add(moved.SizePos(), t_("All mesh"));
	tab.Add(movedUnder.SizePos(), t_("Only underwater"));
	Add(tab.SizePos());	
	AddFrame(status);
	
	OnRefresh();
	timeCallback.Set(-1000, THISBACK(OnTimer));
}

void MainViewDataEach::OnRefresh() {
	const MeshData &mesh = dataSourceFacetsAll[0].GetMesh();
	int num;
	
	num = mesh.mesh.panels.GetCount();
	arrayFacetsAll2.array.GoBegin();
	arrayFacetsAll2.array.Clear();
	arrayFacetsAll2.array.ClearSelection();
	arrayFacetsAll2.array.SetVirtualCount(num);
	arrayFacetsAll2.array.Refresh();
	arrayFacetsAll2.numRows.SetText(FormatInt(num));
		
	num = mesh.under.panels.GetCount();
	arrayFacetsUnder.array.Clear();
	arrayFacetsUnder.array.ClearSelection();
	arrayFacetsUnder.array.SetVirtualCount(num);
	arrayFacetsUnder.array.Refresh();
	arrayFacetsUnder.numRows.SetText(FormatInt(num));
	
	num = mesh.mesh.nodes.GetCount();
	arrayNodesMoved.array.Clear();
	arrayNodesMoved.array.ClearSelection();
	arrayNodesMoved.array.SetVirtualCount(num);
	arrayNodesMoved.array.Refresh();
	arrayNodesMoved.numRows.SetText(FormatInt(num));
	
	num = mesh.under.nodes.	GetCount();
	arrayNodesUnder.array.Clear();
	arrayNodesUnder.array.ClearSelection();
	arrayNodesUnder.array.SetVirtualCount(num);
	arrayNodesUnder.array.Refresh();
	arrayNodesUnder.numRows.SetText(FormatInt(num));
}

void MainViewDataEach::OnTimer() {
	switch (lastSel) {
	case 0:	arrayFacetsAll2.array.WhenSel();	break;
	case 1:	arrayFacetsUnder.array.WhenSel();	break;
	case 2:	arrayNodesMoved.array.WhenSel();	break;
	case 3:	arrayNodesUnder.array.WhenSel();	break;
	}
}

void MainMeshW::Init(MainMesh &_mesh, const Image &icon, const Image &largeIcon) {
	LoadFromJson(mesh, StoreAsJson(_mesh));
	mesh.Init();
	Add(mesh.SizePos());
	Title(t_("BEMRosetta Mesh Viewer")).Sizeable().Zoomable().Icon(icon, largeIcon);
}