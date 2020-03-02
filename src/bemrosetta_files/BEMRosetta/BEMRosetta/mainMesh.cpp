#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <CtrlScroll/CtrlScroll.h>

#include <BEMRosetta/BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;

#include "main.h"


void MainMesh::Init() {
	CtrlLayout(*this);
	
	OnOpt();
	
	CtrlLayout(menuOpen);
	menuOpen.file.WhenChange = THISBACK(OnLoad);
	menuOpen.file.BrowseRightWidth(40).UseOpenFolder(true).BrowseOpenFolderWidth(10);
	menuOpen.butLoad.Tip(t_("Loads mesh file")).WhenAction = [&] {menuOpen.file.DoGo();};

	ArrayModel_Init(menuOpen.arrayModel).MultiSelect(); 
	
	menuOpen.butRemove.Tip(t_("Removes all loaded files")).Disable();	
	menuOpen.butRemove.WhenAction = THISBACK(OnRemove);
	menuOpen.butRemoveSelected.Tip(t_("Removes selected files")).Disable();	
	menuOpen.butRemoveSelected.WhenAction = THISBACK1(OnRemoveSelected, false);
	menuOpen.butJoin.Tip(t_("Join selected meshes")).Disable();	
	//menuOpen.butRemoveSelected.WhenAction = THISBACK1(OnRemoveSelected, false);
	menuOpen.butSplit.Tip(t_("Split mesh in parts")).Disable();	
	//menuOpen.butRemoveSelected.WhenAction = THISBACK1(OnRemoveSelected, false);
	
	CtrlLayout(menuConvert);
	menuConvert.file.WhenChange = THISBACK(OnConvertMesh);
	menuConvert.file.BrowseRightWidth(40).UseOpenFolder(true).BrowseOpenFolderWidth(10);
	menuConvert.butLoad.WhenAction = [&] {OnConvertMesh();};

	ArrayModel_Init(menuConvert.arrayModel);
	
	menuConvert.arrayModel.WhenCursor = THISBACK(OnMenuConvertArraySel);
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
	
	ArrayModel_Init(menuPlot.arrayModel, true);
	
	OnOpt();
	
	CtrlLayout(menuProcess);
	menuProcess.cg_x <<= 0;
	menuProcess.cg_x.WhenEnter = THISBACK1(OnUpdate, false);
	menuProcess.cg_y <<= 0;
	menuProcess.cg_y.WhenEnter = THISBACK1(OnUpdate, false);
	menuProcess.cg_z <<= 0;
	menuProcess.cg_z.WhenEnter = THISBACK1(OnUpdate, false);
	menuProcess.mass <<= 0;
	menuProcess.mass.WhenEnter = THISBACK1(OnUpdate, false);
	menuProcess.butUpdateCg <<= THISBACK1(OnUpdate, false);
	
	menuProcess.t_x <<= 0;
	menuProcess.t_x.WhenEnter = THISBACK1(OnUpdate, false);
	menuProcess.t_y <<= 0;
	menuProcess.t_y.WhenEnter = THISBACK1(OnUpdate, false);
	menuProcess.t_z <<= 0;
	menuProcess.t_z.WhenEnter = THISBACK1(OnUpdate, false);
	menuProcess.a_x <<= 0;
	menuProcess.a_x.WhenEnter = THISBACK1(OnUpdate, false);
	menuProcess.a_y <<= 0;
	menuProcess.a_y.WhenEnter = THISBACK1(OnUpdate, false);
	menuProcess.a_z <<= 0;
	menuProcess.a_z.WhenEnter = THISBACK1(OnUpdate, false);
	menuProcess.c_x <<= 0;
	menuProcess.c_x.WhenEnter = THISBACK1(OnUpdate, false);
	menuProcess.c_y <<= 0;
	menuProcess.c_y.WhenEnter = THISBACK1(OnUpdate, false);
	menuProcess.c_z <<= 0;
	menuProcess.c_z.WhenEnter = THISBACK1(OnUpdate, false);
	menuProcess.moveType <<= 0;
	menuProcess.moveType.Transparent(false);
	menuProcess.butUpdatePos <<= THISBACK1(OnUpdate, false);
	menuProcess.butUpdateAng <<= THISBACK1(OnUpdate, false);
	menuProcess.butHealing <<= THISBACK(OnHealing);
	menuProcess.butImageX <<= THISBACK1(OnImage, 0);
	menuProcess.butImageY <<= THISBACK1(OnImage, 1);
	menuProcess.butImageZ <<= THISBACK1(OnImage, 2);
	
	ArrayModel_Init(menuProcess.arrayModel); 
	
	menuProcess.arrayModel.WhenCursor = THISBACK(OnMenuProcessArraySel);

	menuTab.Add(menuOpen.SizePos(),    	t_("Load"));
	menuTab.Add(menuPlot.SizePos(),    	t_("Plot")).Disable();
	menuTab.Add(menuProcess.SizePos(),t_("Process")).Disable();
	menuTab.Add(menuConvert.SizePos(), 	t_("Save as")).Disable();
		
	mainView.Init(menuPlot, menuOpen.arrayModel);
	mainViewData.Init();
	mainVAll.Horz(mainView, mainViewData);
	mainVAll.SetPositions(6000, 9970).SetInitialPositionId(1).SetButtonNumber(1);
	mainVAll.WhenAction = [&] {mainView.SetPaintSelect(mainVAll.GetPos() < 9950);};
	mainTab.Add(mainVAll.SizePos(), t_("View"));
	
	mainSummary.Init();
	mainTab.Add(mainSummary.SizePos(), t_("Summary"));

	mainStiffness.Init();
	mainTab.Add(mainStiffness.SizePos(), t_("K Stiffness Matrix"));
			
	mainTab.WhenSet = [&] {
		LOGTAB(mainTab);
		Vector<int> ids = ArrayModel_IdsMesh(menuOpen.arrayModel);
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
		TabCtrl::Item& tabmenuProcess = menuTab.GetItem(menuTab.Find(menuProcess));
		tabmenuProcess.Enable(convertProcess);
		if (plot) {
			tabMenuPlot.Text(t_("Plot"));
			menuTab.Set(menuPlot);
		} else {
			tabMenuPlot.Text("");
			menuTab.Set(menuOpen);
		}
		if (convertProcess) {
			tabMenuConvert.Text(t_("Convert"));
			tabmenuProcess.Text(t_("Process"));
		} else {
			tabMenuConvert.Text("");
			tabmenuProcess.Text("");
		}
	};
	mainTab.WhenSet();
	
	menuTab.WhenSet = [&] {
	LOGTAB(menuTab);
		if (menuTab.IsAt(menuConvert)) 
			menuConvert.arrayModel.WhenCursor();
	};
	menuTab.WhenSet();
}

void MainMesh::OnMenuProcessArraySel() {
	int id = ArrayModel_IdMesh(menuProcess.arrayModel);
	if (id < 0)
		return;
	
	MeshData &data = Bem().surfs[id];
	menuProcess.cg_x <<= data.cg.x;
	menuProcess.cg_y <<= data.cg.y;
	menuProcess.cg_z <<= data.cg.z;
	menuProcess.mass <<= data.mass;
	menuProcess.t_x  <<= data.mesh.x;
	menuProcess.t_y  <<= data.mesh.y;
	menuProcess.t_z  <<= data.mesh.z;
	menuProcess.a_x  <<= data.mesh.a_x;
	menuProcess.a_y  <<= data.mesh.a_y;
	menuProcess.a_z  <<= data.mesh.a_z;
	menuProcess.c_x  <<= data.mesh.c_x;
	menuProcess.c_y  <<= data.mesh.c_y;
	menuProcess.c_z  <<= data.mesh.c_z;
	menuProcess.moveType <<= 0;
}

void MainMesh::OnMenuConvertArraySel() {
	int id = ArrayModel_IdMesh(menuConvert.arrayModel);
	if (id < 0)
		return;
	
	String file = ~menuConvert.file;
	String folder = GetFileFolder(file);
	String ext = GetFileExt(file);
	String fileName = GetFileTitle(ArrayModel_GetFileName(menuConvert.arrayModel));
	file = AppendFileName(folder, fileName + ext);
	menuConvert.file <<= file;
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
	if (!ret || IsNull(menuConvert.optMesh)) 
		menuConvert.optMesh = 0;
}

void MainMesh::LoadSelTab(BEMData &bem) {
	const Vector<int> &ids = ArrayModel_IdsMesh(menuOpen.arrayModel);
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

void MainMesh::AfterLoad(String file) {
	int id = Bem().surfs.GetCount() - 1;
	MeshData &surf = Bem().surfs[id];
	
	ArrayModel_Add(menuOpen.arrayModel, surf.GetCodeStr(), GetFileTitle(file), surf.fileName, surf.GetId());
	
	mainView.CalcEnvelope();
	mainView.gl.ZoomToFit();
	mainTab.Set(mainView);
	mainTab.WhenSet();
	
	surf.Report(Bem().rho);
	mainSummary.Report(Bem().surfs, id);		
	
	menuOpen.butRemove.Enable();
	menuOpen.butRemoveSelected.Enable();
	ArrayModel_Add(menuPlot.arrayModel, surf.GetCodeStr(), GetFileTitle(file), surf.fileName, surf.GetId(),
					optionsPlot, [&] {mainView.gl.Refresh();});
	if (ArrayModel_IdMesh(menuPlot.arrayModel) < 0)
		menuPlot.arrayModel.SetCursor(0);
	ArrayModel_Add(menuConvert.arrayModel, surf.GetCodeStr(), GetFileTitle(file), surf.fileName, surf.GetId());
	if (ArrayModel_IdMesh(menuConvert.arrayModel) < 0)
		menuConvert.arrayModel.SetCursor(0);
	ArrayModel_Add(menuProcess.arrayModel, surf.GetCodeStr(), GetFileTitle(file), surf.fileName, surf.GetId());
	if (ArrayModel_IdMesh(menuProcess.arrayModel) < 0)
		menuProcess.arrayModel.SetCursor(0);		
}

bool MainMesh::OnLoad() {
	GuiLock __;
	
	String file = ~menuOpen.file;
		
	try {
		WaitCursor waitcursor;
		Progress progress(t_("Loading mesh file..."), 100); 
		
		Vector<int> ids = ArrayModel_IdsMesh(menuOpen.arrayModel);
		for (int i = 0; i < ids.GetCount(); ++i) {
			if (Bem().surfs[ids[i]].fileName == file) {
				if (!PromptYesNo(t_("Model is already loaded") + S("&") + t_("Do you wish to open it anyway?")))
					return false;
				break;
			}
		}
		mainView.gl.Disable();
		
		Bem().LoadMesh(file, [&](String str, int _pos) {progress.SetText(str); progress.SetPos(_pos);}, false);
		
		AfterLoad(file);
		
		mainView.gl.Enable();
		mainView.gl.Refresh();
		
		mainViewData.OnAddedModel(mainView);
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
		return false;
	}
	return true;
}

bool MainMesh::OnConvertMesh() {
	try {
		int id = ArrayModel_IdMesh(menuConvert.arrayModel);
		if (id < 0) {
			Exclamation(t_("Please select a model to process"));
			return false;
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
								static_cast<MeshData::MESH_TYPE>(int(~menuConvert.optMesh)));	
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
		return false;
	}
	return true;
}

void MainMesh::OnUpdate(bool forceMoved) {
	try {
		Vector<int> ids = ArrayModel_IdsMesh(menuOpen.arrayModel);
		int id = ArrayModel_IdMesh(menuProcess.arrayModel);
		if (id < 0) {
			Exclamation(t_("Please select a model to process"));
			return;
		}
		
		WaitCursor wait;
		
		MeshData &data = Bem().surfs[id];

		data.mass = ~menuProcess.mass;
		data.cg0.Set(~menuProcess.cg_x, ~menuProcess.cg_y, ~menuProcess.cg_z);
		
		double t_x, t_y, t_z, a_x, a_y, a_z, c_x, c_y, c_z;
		if (~menuProcess.moveType == 0) {
			t_x = ~menuProcess.t_x;
			t_y = ~menuProcess.t_y;
			t_z = ~menuProcess.t_z;
			a_x = ~menuProcess.a_x;
			a_y = ~menuProcess.a_y;
			a_z = ~menuProcess.a_z;
		} else {
			t_x = data.mesh.x   + double(~menuProcess.t_x);
			t_y = data.mesh.y   + double(~menuProcess.t_y);
			t_z = data.mesh.z   + double(~menuProcess.t_z);
			a_x = data.mesh.a_x + double(~menuProcess.a_x);
			a_y = data.mesh.a_y + double(~menuProcess.a_y);
			a_z = data.mesh.a_z + double(~menuProcess.a_z);
		}
		c_x = ~menuProcess.c_x;
		c_y = ~menuProcess.c_y;
		c_z = ~menuProcess.c_z;
		
		if (IsNull(t_x) || IsNull(t_y) || IsNull(t_z)) {
			Exclamation(t_("Please fill translation data"));
			return;
		}
		if (IsNull(a_x) || IsNull(a_y) || IsNull(a_z)) {
			Exclamation(t_("Please fill rotation data"));
			return;
		}
		if (IsNull(c_x) || IsNull(c_y) || IsNull(c_z)) {
			Exclamation(t_("Please fill center of rotation data"));
			return;
		}
		
		data.cg.MoveTo(data.cg0, t_x, t_y, t_z, a_x, a_y, a_z, c_x, c_y, c_z);

		bool isMoved = forceMoved || data.mesh.IsMoved(t_x, t_y, t_z, a_x, a_y, a_z);
						
		if (isMoved)
			data.mesh.MoveTo(t_x, t_y, t_z, a_x, a_y, a_z, c_x, c_y, c_z);
	
		data.AfterLoad(Bem().rho, Bem().g, !isMoved);
		
	 	mainStiffness.Load(Bem().surfs, ids);
		mainView.CalcEnvelope();
		mainSummary.Report(Bem().surfs, id);
		
		mainView.gl.Refresh();
		mainViewData.OnRefresh();
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}
}

void MainMesh::OnHealing() {
	try {
		int id = ArrayModel_IdMesh(menuProcess.arrayModel);
		if (id < 0) {
			Exclamation(t_("Please select a model to process"));
			return;
		}
		
		WaitCursor waitcursor;
		Progress progress(t_("Healing mesh file..."), 100); 
		mainView.gl.Disable();
		
		Bem().HealingMesh(id, [&](String str, int _pos) {progress.SetText(str); progress.SetPos(_pos);});
		
		OnUpdate(true);
		mainView.gl.Enable();
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}	
}

void MainMesh::OnImage(int axis) {
	String saxis = (axis == 0) ? "X" : ((axis == 1) ? "Y" : "Z");

	try {
		Vector<int> ids = ArrayModel_IdsMesh(menuOpen.arrayModel);
		int id = ArrayModel_IdMesh(menuProcess.arrayModel);
		if (id < 0) {
			Exclamation(t_("Please select a model to process"));
			return;
		}
		
		WaitCursor waitcursor;
		
		MeshData &data = Bem().surfs[id];

		data.mass = ~menuProcess.mass;
		if (axis == 0)
			data.cg0.x = -data.cg0.x;
		else if (axis == 1)
			data.cg0.y = -data.cg0.y;
		else
			data.cg0.z = -data.cg0.z;
		
		data.cg = data.cg0;

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
	menuOpen.arrayModel.Clear();
	menuOpen.butRemove.Disable();
	menuOpen.butRemoveSelected.Disable();
	menuPlot.arrayModel.Clear();
	optionsPlot.Clear();
	menuConvert.arrayModel.Clear();
	menuProcess.arrayModel.Clear();
	
	mainStiffness.Clear();
	mainViewData.Clear();
	
	Bem().surfs.Clear();
	mainView.env.Reset();
	
	mainTab.WhenSet();
	mainView.gl.Refresh();
}

void MainMesh::OnRemoveSelected(bool all) {	
	bool selected = false;
	for (int r = menuOpen.arrayModel.GetCount()-1; r >= 0; --r) {
		if (all || menuOpen.arrayModel.IsSelected(r)) {
			int id = ArrayModel_IdMesh(menuOpen.arrayModel, r);
			Bem().surfs.Remove(id);
			menuOpen.arrayModel.Remove(r);
			menuPlot.arrayModel.Remove(r);
			menuConvert.arrayModel.Remove(r);
			menuProcess.arrayModel.Remove(r);
			selected = true;
		}
	}	// Only one available => directly selected
	if (!selected && menuOpen.arrayModel.GetCount() == 1) {
		int id = ArrayModel_IdMesh(menuOpen.arrayModel, 0);
		Bem().surfs.Remove(id);
		menuOpen.arrayModel.Remove(0);
		menuPlot.arrayModel.Remove(0);
		menuConvert.arrayModel.Remove(0);
		menuProcess.arrayModel.Remove(0);
		selected = true;		
	}	
	if (!selected) {
		Exclamation(t_("No model selected"));
		return;
	}
 	mainSummary.Clear();
	//for (int i = 0; i < Bem().surfs.GetCount(); ++i) 
	for (int row = 0; row < menuOpen.arrayModel.GetCount(); ++row) {
		int id = ArrayModel_IdMesh(menuOpen.arrayModel, row);
		mainSummary.Report(Bem().surfs, id);
	}
	Vector<int> ids = ArrayModel_IdsMesh(menuOpen.arrayModel);
	mainStiffness.Load(Bem().surfs, ids);
	mainViewData.Clear();//OnRefresh();
	
	int numrow = menuOpen.arrayModel.GetCount();
	menuOpen.butRemove.Enable(numrow > 0);
	menuOpen.butRemoveSelected.Enable(numrow > 0);
	
	mainView.CalcEnvelope();
	
	mainTab.WhenSet();
	mainView.gl.Refresh();	
}

void MainMesh::Jsonize(JsonIO &json) {
	json
		("menuOpen_file", menuOpen.file)
		("menuConvert_file", menuConvert.file)	
		("menuConvert_opt", menuConvert.opt)
		("menuConvert_optMesh", menuConvert.optMesh)
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
	String name = GetFileTitle(data.fileName);
	
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

	array.Set(row, 0, t_("Surface [m2]"));		array.Set(row++, col, FormatDouble(data.mesh.surface, 3, FD_EXP));
	array.Set(row, 0, t_("Volume [m3]"));		array.Set(row++, col, Format(t_("%s (%s, %s, %s)"), 
														FormatDouble(data.mesh.volume,  1, FD_EXP),
														FormatDouble(data.mesh.volumex, 1, FD_EXP),
														FormatDouble(data.mesh.volumey, 1, FD_EXP),
														FormatDouble(data.mesh.volumez, 1, FD_EXP)));
	
	array.Set(row, 0, t_("Immersed surface [m2]"));array.Set(row++, col, FormatDouble(data.under.surface, 3, FD_EXP));
	array.Set(row, 0, t_("Immersed volume [m3]")); array.Set(row++, col, Format(t_("%s (%s, %s, %s)"), 
														FormatDouble(data.under.volume,  1, FD_EXP),
														FormatDouble(data.under.volumex, 1, FD_EXP),
														FormatDouble(data.under.volumey, 1, FD_EXP),
														FormatDouble(data.under.volumez, 1, FD_EXP)));
	array.Set(row, 0, t_("Displacement [tm]")); array.Set(row++, col, FormatDouble(data.under.volume*Bem().rho/1000, 3, FD_EXP));
	array.Set(row, 0, t_("Cg [m]"));			array.Set(row++, col, Format(t_("%s, %s, %s"),
														FormatDouble(data.cg.x, 3, FD_EXP),			
														FormatDouble(data.cg.y, 3, FD_EXP),
														FormatDouble(data.cg.z, 3, FD_EXP)));
	array.Set(row, 0, t_("Cb [m]"));			array.Set(row++, col, Format(t_("%s, %s, %s"),
														FormatDouble(data.cb.x, 3, FD_EXP),			
														FormatDouble(data.cb.y, 3, FD_EXP),
														FormatDouble(data.cb.z, 3, FD_EXP)));
	array.Set(row, 0, t_("Water Plane Area [m2]"));	array.Set(row++, col, FormatDouble(data.waterPlaneArea, 3, FD_EXP));
	
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

void MainView::Init(const WithMenuPlotMesh<StaticRect> &_menuPlot, const ArrayCtrl &_array) {
	CtrlLayout(*this);
	menuPlot = &_menuPlot;
	arrayModel = &_array;
	
	gl.SetEnv(env);
	gl.WhenPaint = THISBACK(OnPaint);	
}
	
void MainView::OnPaint() {
	if (~GetMenuPlot().showAxis) 
		gl.PaintAxis(0, 0, 0, env.LenRef()/4.);	
	
	if (~GetMenuPlot().showLimits) 
		gl.PaintCuboid(Point3D(env.maxX, env.maxY, env.maxZ), Point3D(env.minX, env.minY, env.minZ), Gray());
	
	for (int row = 0; row < arrayModel->GetCount(); ++row) {
		if (ArrayModel_IsVisible(GetMenuPlot().arrayModel, row)) {
			int id = ArrayModel_IdMesh(*arrayModel, row);
			
			double len = env.LenRef()/10;
			bool showNormals = ~GetMenuPlot().showNormals && ~GetMenuPlot().showMesh;
			bool showNormalsUnderwater = ~GetMenuPlot().showNormals && ~GetMenuPlot().showUnderwater;
			if (~GetMenuPlot().showNormals && !~GetMenuPlot().showMesh && !~GetMenuPlot().showUnderwater)
				showNormals = true;
			
			const Upp::Color &color = ArrayModel_GetColor(GetMenuPlot().arrayModel, row);
			const MeshData &mesh = Bem().surfs[id];
			gl.PaintSurface(mesh.mesh, color, ~GetMenuPlot().showMesh, 	
				showNormals, false, ~GetMenuPlot().showSkewed, 
				~GetMenuPlot().showFissure, ~GetMenuPlot().showMultiPan);
			gl.PaintSurface(mesh.under, color, ~GetMenuPlot().showUnderwater, 
				showNormalsUnderwater, ~GetMenuPlot().showWaterLevel, false, 
				false, false);
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
					const Vector<int> &nod = mesh.mesh.GetSelNodes();
					for (int in = 0; in < nod.GetCount(); ++in)
						gl.PaintCube(mesh.mesh.nodes[nod[in]], len/20, LtBlue());
					const Vector<int> &pan = mesh.mesh.GetSelPanels();
					const Vector<Point3D> &nodes = mesh.mesh.nodes;
					for (int ip = 0; ip < pan.GetCount(); ++ip) {
						const Panel &panel = mesh.mesh.panels[pan[ip]];
						gl.PaintQuad(nodes[panel.id[0]], nodes[panel.id[1]], nodes[panel.id[2]], nodes[panel.id[3]], LtRed(), .2);
					}
				}
				if (~GetMenuPlot().showUnderwater) {
					const Vector<int> &nod = mesh.under.GetSelNodes();
					for (int in = 0; in < nod.GetCount(); ++in)
						gl.PaintCube(mesh.under.nodes[nod[in]], len/20, LtBlue());
					const Vector<int> &pan = mesh.under.GetSelPanels();
					const Vector<Point3D> &nodes = mesh.under.nodes;
					for (int ip = 0; ip < pan.GetCount(); ++ip) {
						const Panel &panel = mesh.under.panels[pan[ip]];
						gl.PaintQuad(nodes[panel.id[0]], nodes[panel.id[1]], nodes[panel.id[2]], nodes[panel.id[3]], LtRed(), .2);
					}
				}
			}
		}
	}
}

void MainView::CalcEnvelope() {
	env.Reset();
	//for (int i = 0; i < Bem().surfs.GetCount(); ++i)
	for (int row = 0; row < arrayModel->GetCount(); ++row) {
		int id = ArrayModel_IdMesh(*arrayModel, row);
		env.MixEnvelope(Bem().surfs[id].mesh.env);
	}
}

void MainMesh::DragAndDrop(Point , PasteClip& d) {
	if (IsDragAndDropSource())
		return;
	if (AcceptFiles(d)) {
		Vector<String> files = GetFiles(d);
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
}

bool MainMesh::Key(dword key, int ) {
	if (key == K_CTRL_V) {
		Vector<String> files = GetFiles(Ctrl::Clipboard());
		bool followWithErrors = false;
		for (int i = 0; i < files.GetCount(); ++i) {
			menuOpen.file <<= files[i];
			Status(Format(t_("Loading '%s'"), files[i]));
			if (!OnLoad() && !followWithErrors && files.GetCount() - i > 1) {
				if (!PromptYesNo(Format(t_("Do you wish to load the pending %d files?"), files.GetCount() - i - 1)))
					return true;
				followWithErrors = true;
			}
			ProcessEvents();
		}
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
	tab.Add(model.SizePos(), GetFileTitle(Bem().surfs[id].fileName));
}

void MainViewData::OnRefresh() {
	for (int i = 0; i < models.GetCount(); ++i)
		models[i].OnRefresh();
}

void MainViewData::Clear() {
	models.Clear();
	tab.Reset();
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
	if (origMovedUnder == 0 && pmesh->mesh.nodes0.GetCount() <= iq)
		return Null;
	if (origMovedUnder == 1 && pmesh->mesh.nodes.GetCount() <= iq)
		return Null;
	if (origMovedUnder == 2 && pmesh->under.nodes.GetCount() <= iq)
		return Null;
	
	const Point3D &p = origMovedUnder == 0 ? pmesh->mesh.nodes0[iq] : (origMovedUnder == 1 ? pmesh->mesh.nodes[iq] : pmesh->under.nodes[iq]);
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
		selectedPanels = ArrayCtrlGetSelected(arrayFacetsAll2.array);
		selectedNodes = ArrayCtrlGetSelected(arrayNodesMoved.array);
	} else {
		show = mainMesh.GetShowUnderwater();
		selectedPanels = ArrayCtrlGetSelected(arrayFacetsUnder.array);
		selectedNodes = ArrayCtrlGetSelected(arrayNodesUnder.array);
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
	arrayFacetsAll2.array.WhenCursor = [&] {
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
	arrayFacetsUnder.array.WhenCursor = [&] {
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
	arrayNodesMoved.array.WhenCursor = [&] {
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
	arrayNodesUnder.array.WhenCursor = [&] {
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
	
	num = mesh.under.nodes.GetCount();
	arrayNodesUnder.array.Clear();
	arrayNodesUnder.array.ClearSelection();
	arrayNodesUnder.array.SetVirtualCount(num);
	arrayNodesUnder.array.Refresh();
	arrayNodesUnder.numRows.SetText(FormatInt(num));
}

void MainViewDataEach::OnTimer() {
	switch (lastSel) {
	case 0:	arrayFacetsAll2.array.WhenCursor();		break;
	case 1:	arrayFacetsUnder.array.WhenCursor();	break;
	case 2:	arrayNodesMoved.array.WhenCursor();		break;
	case 3:	arrayNodesUnder.array.WhenCursor();		break;
	}
}
