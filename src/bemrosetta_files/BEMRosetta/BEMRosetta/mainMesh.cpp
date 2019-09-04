#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas.h>

#include <BEMRosetta/BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;

#include "main.h"


void MainMesh::Init() {
	CtrlLayout(*this);
	
	OnOpt();
	
	CtrlLayout(menuOpen);
	menuOpen.file.WhenChange = THISBACK(OnView);
	menuOpen.file.BrowseRightWidth(40).UseOpenFolder(true).BrowseOpenFolderWidth(10);
	menuOpen.butLoad.WhenAction = [&] {OnView();};

	menuOpen.arrayModel.NoHeader().NoVertGrid().AutoHideSb();
	menuOpen.arrayModel.AddColumn("", 20);	
	menuOpen.arrayModel.AddColumn("", 20); 
	
	menuOpen.butRemove.Disable();	
	menuOpen.butRemove.WhenAction = [&] {
		WaitCursor waitcursor;
		
		mainSummary.Clear();
		menuOpen.arrayModel.Clear();
		menuOpen.butRemove.Disable();
		menuConvert.arrayModel.Clear();
		menuStability.arrayModel.Clear();
		
		mainStiffness.Clear();
		mainViewData.Clear();
		
		ma().bem.surfs.Clear();
		mainView.env.Reset();
		
		mainTab.WhenSet();
		mainView.gl.Refresh();
	};
	
	CtrlLayout(menuConvert);
	menuConvert.file.WhenChange = THISBACK(OnConvertMesh);
	menuConvert.file.BrowseRightWidth(40).UseOpenFolder(true).BrowseOpenFolderWidth(10);
	menuConvert.butLoad.WhenAction = [&] {OnConvertMesh();};

	menuConvert.arrayModel.NoHeader().NoVertGrid().AutoHideSb();
	menuConvert.arrayModel.AddColumn("", 20);	
	menuConvert.arrayModel.AddColumn("", 20);
	
	menuConvert.opt.WhenAction = [&] {OnOpt();};

	OnOpt();
	
	CtrlLayout(menuPlot);
	menuPlot.showMesh.WhenAction 	    = [&] {LoadSelTab(ma().bem);};
	menuPlot.showNormals.WhenAction     = [&] {LoadSelTab(ma().bem);};
	menuPlot.showWaterLevel.WhenAction  = [&] {LoadSelTab(ma().bem);};
	menuPlot.showSkewed.WhenAction      = [&] {LoadSelTab(ma().bem);};
	menuPlot.showFissure.WhenAction     = [&] {LoadSelTab(ma().bem);};
	menuPlot.showMultiPan.WhenAction    = [&] {LoadSelTab(ma().bem);};
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
		
	OnOpt();
	
	CtrlLayout(menuStability);
	menuStability.cg_x <<= 0;
	menuStability.cg_x.WhenEnter = THISBACK1(OnUpdate, false);
	menuStability.cg_y <<= 0;
	menuStability.cg_y.WhenEnter = THISBACK1(OnUpdate, false);
	menuStability.cg_z <<= 0;
	menuStability.cg_z.WhenEnter = THISBACK1(OnUpdate, false);
	menuStability.butUpdateCg <<= THISBACK1(OnUpdate, false);
	menuStability.t_x <<= 0;
	menuStability.t_x.WhenEnter = THISBACK1(OnUpdate, false);
	menuStability.t_y <<= 0;
	menuStability.t_y.WhenEnter = THISBACK1(OnUpdate, false);
	menuStability.t_z <<= 0;
	menuStability.t_z.WhenEnter = THISBACK1(OnUpdate, false);
	menuStability.a_x <<= 0;
	menuStability.a_x.WhenEnter = THISBACK1(OnUpdate, false);
	menuStability.a_y <<= 0;
	menuStability.a_y.WhenEnter = THISBACK1(OnUpdate, false);
	menuStability.a_z <<= 0;
	menuStability.a_z.WhenEnter = THISBACK1(OnUpdate, false);
	menuStability.c_x <<= 0;
	menuStability.c_x.WhenEnter = THISBACK1(OnUpdate, false);
	menuStability.c_y <<= 0;
	menuStability.c_y.WhenEnter = THISBACK1(OnUpdate, false);
	menuStability.c_z <<= 0;
	menuStability.c_z.WhenEnter = THISBACK1(OnUpdate, false);
	menuStability.moveType <<= 0;
	menuStability.moveType.Transparent(false);
	menuStability.butUpdatePos <<= THISBACK1(OnUpdate, false);
	menuStability.butHealing <<= THISBACK(OnHealing);
	
	menuStability.arrayModel.NoHeader().NoVertGrid().AutoHideSb();
	menuStability.arrayModel.AddColumn("", 20);	
	menuStability.arrayModel.AddColumn("", 20); 
	
	menuStability.arrayModel.WhenSel = THISBACK(OnMenuConvertArraySel);

	menuTab.Add(menuOpen.SizePos(),    	t_("Open"));
	menuTab.Add(menuPlot.SizePos(),    	t_("Plot")).Disable();
	menuTab.Add(menuStability.SizePos(),t_("Process")).Disable();
	menuTab.Add(menuConvert.SizePos(), 	t_("Convert")).Disable();
		
	mainView.Init(menuPlot);
	mainTab.Add(mainView.SizePos(), t_("View"));
	
	mainSummary.Init();
	mainTab.Add(mainSummary.SizePos(), t_("Summary"));

	mainStiffness.Init();
	mainTab.Add(mainStiffness.SizePos(), t_("K Stiffness Matrix"));
	
	mainViewData.Init();
	mainTab.Add(mainViewData.SizePos(), t_("Data"));
		
	mainTab.WhenSet = [&] {
		bool plot = true, convertProcess = true;
		if (ma().bem.surfs.IsEmpty()) {
			plot = false;
			convertProcess = false;
		} else if (mainTab.IsAt(mainView)) 
			;
		else if (mainTab.IsAt(mainStiffness)) {
			plot = false;
			mainStiffness.Load(ma().bem.surfs);
		} else 
			plot = false;
		
		TabCtrl::Item& tabMenuPlot = menuTab.GetItem(menuTab.Find(menuPlot));
		tabMenuPlot.Enable(plot);
		TabCtrl::Item& tabMenuConvert = menuTab.GetItem(menuTab.Find(menuConvert));
		tabMenuConvert.Enable(convertProcess);
		TabCtrl::Item& tabMenuStability = menuTab.GetItem(menuTab.Find(menuStability));
		tabMenuStability.Enable(convertProcess);
		if (plot) {
			tabMenuPlot.Text(t_("Plot"));
			menuTab.Set(menuPlot);
		} else {
			tabMenuPlot.Text("");
			menuTab.Set(menuOpen);
		}
		if (convertProcess) {
			tabMenuConvert.Text(t_("Convert"));
			tabMenuStability.Text(t_("Process"));
		} else {
			tabMenuConvert.Text("");
			tabMenuStability.Text("");
		}
	};
	mainTab.WhenSet();
}

void MainMesh::OnMenuConvertArraySel() {
	int id = menuStability.arrayModel.GetCursor();
	if (id < 0)
		return;
	MeshData &data = ma().bem.surfs[id];
	menuStability.cg_x <<= data.cg.x;
	menuStability.cg_y <<= data.cg.y;
	menuStability.cg_z <<= data.cg.z;
	menuStability.t_x  <<= data.mesh.x;
	menuStability.t_y  <<= data.mesh.y;
	menuStability.t_z  <<= data.mesh.z;
	menuStability.a_x  <<= data.mesh.a_x;
	menuStability.a_y  <<= data.mesh.a_y;
	menuStability.a_z  <<= data.mesh.a_z;
	menuStability.c_x  <<= data.mesh.c_x;
	menuStability.c_y  <<= data.mesh.c_y;
	menuStability.c_z  <<= data.mesh.c_z;
	menuStability.moveType = 0;
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
	if (!ret || IsNull(menuConvert.optPosition)) 
		menuConvert.optPosition = 0;
}

void MainMesh::LoadSelTab(BEMData &bem) {
	if (mainTab.Get() == mainTab.Find(mainStiffness))
		mainStiffness.Load(bem.surfs);
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
	case 2:	menuConvert.file <<= ForceExtSafe(~menuConvert.file, ".stl"); 	
			menuConvert.file.Type("STL binary .stl file", "*.stl");
			break;
	case 3:	menuConvert.file <<= ForceExtSafe(~menuConvert.file, ".stl"); 	
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
	mainView.CalcEnvelope();
	mainView.gl.ZoomToFit();
	mainTab.Set(mainView);
	mainTab.WhenSet();
	
	int id = ma().bem.surfs.GetCount() - 1;
	MeshData &surf = ma().bem.surfs[id];
	
	surf.Report();
	mainSummary.Report(surf, id);		
	
	menuOpen.arrayModel.Add(surf.GetCodeStr(), GetFileTitle(file));
	menuOpen.butRemove.Enable();
	menuConvert.arrayModel.Add(surf.GetCodeStr(), GetFileTitle(file));
	if (menuConvert.arrayModel.GetCursor() < 0)
		menuConvert.arrayModel.SetCursor(0);
	menuStability.arrayModel.Add(surf.GetCodeStr(), GetFileTitle(file));
	if (menuStability.arrayModel.GetCursor() < 0)
		menuStability.arrayModel.SetCursor(0);		
}

bool MainMesh::OnView() {
	String file = ~menuOpen.file;
	
	try {
		WaitCursor waitcursor;
		Progress progress(t_("Loading mesh file..."), 100); 
		mainView.gl.Disable();
		
		ma().bem.LoadMesh(file, [&](String str, int _pos) {progress.SetText(str); progress.SetPos(_pos);});
		
		AfterLoad(file);
		
		mainView.gl.Enable();
		mainView.gl.Refresh();
		
		mainViewData.OnAddedModel();
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
		return false;
	}
	return true;
}

bool MainMesh::OnConvertMesh() {
	try {
		int id = menuConvert.arrayModel.GetCursor();
		if (id < 0) {
			Exclamation(t_("Please select a model to process"));
			return false;
		}
		MeshData::MESH_FMT type;	
		switch (menuConvert.opt) {
		case 0:	type = MeshData::WAMIT_GDF;	break;
		case 1:	type = MeshData::NEMOH_DAT;	break;
		case 2:	type = MeshData::STL_BIN;	break;
		case 3:	type = MeshData::STL_TXT;	break;
		case 4:	type = MeshData::UNKNOWN;	break;
		default: throw Exc(t_("Unknown type in OnConvert()"));
		}
		
		WaitCursor waitcursor;
		Progress progress(t_("Saving mesh file..."), 100); 
		
		bool meshAll = ~menuConvert.optMesh == 0;
		bool positionOriginal = ~menuConvert.optPosition == 0;
		
		ma().bem.surfs[id].SaveAs(~menuConvert.file, type, ma().bem.g, meshAll, positionOriginal);	
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
		return false;
	}
	return true;
}

void MainMesh::OnUpdate(bool forceMoved) {
	try {
		int id = menuStability.arrayModel.GetCursor();
		if (id < 0) {
			Exclamation(t_("Please select a model to process"));
			return;
		}
		WaitCursor wait;
		
		MeshData &data = ma().bem.surfs[id];

		data.cg0.Set(~menuStability.cg_x, ~menuStability.cg_y, ~menuStability.cg_z);
		
		double t_x, t_y, t_z, a_x, a_y, a_z, c_x, c_y, c_z;
		if (~menuStability.moveType == 0) {
			t_x = ~menuStability.t_x;
			t_y = ~menuStability.t_y;
			t_z = ~menuStability.t_z;
			a_x = ~menuStability.a_x;
			a_y = ~menuStability.a_y;
			a_z = ~menuStability.a_z;
		} else {
			t_x = data.mesh.x   + static_cast<double>(~menuStability.t_x);
			t_y = data.mesh.y   + static_cast<double>(~menuStability.t_y);
			t_z = data.mesh.z   + static_cast<double>(~menuStability.t_z);
			a_x = data.mesh.a_x + static_cast<double>(~menuStability.a_x);
			a_y = data.mesh.a_y + static_cast<double>(~menuStability.a_y);
			a_z = data.mesh.a_z + static_cast<double>(~menuStability.a_z);
		}
		c_x = ~menuStability.c_x;
		c_y = ~menuStability.c_y;
		c_z = ~menuStability.c_z;
				
		data.cg.MoveTo(data.cg0, t_x, t_y, t_z, a_x, a_y, a_z, c_x, c_y, c_z);

		bool isMoved = forceMoved || data.mesh.IsMoved(t_x, t_y, t_z, a_x, a_y, a_z);
						
		if (isMoved)
			data.mesh.MoveTo(t_x, t_y, t_z, a_x, a_y, a_z, c_x, c_y, c_z);
	
		data.AfterLoad(ma().bem.rho, ma().bem.g, !isMoved);
		
	 	mainStiffness.Load(ma().bem.surfs);
		mainSummary.Report(data, id);
		mainView.CalcEnvelope();
		mainView.gl.Refresh();
		mainViewData.OnRefresh();
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}
}

void MainMesh::OnHealing() {
	try {
		int id = menuStability.arrayModel.GetCursor();
		if (id < 0) {
			Exclamation(t_("Please select a model to process"));
			return;
		}
		WaitCursor waitcursor;
		Progress progress(t_("Healing mesh file..."), 100); 
		mainView.gl.Disable();
		
		ma().bem.HealingMesh(id, [&](String str, int _pos) {progress.SetText(str); progress.SetPos(_pos);});
		
		OnUpdate(true);
		mainView.gl.Enable();
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}	
}

void MainMesh::Jsonize(JsonIO &json) {
	json
		("menuOpen_file", menuOpen.file)
		("menuConvert_file", menuConvert.file)
		("menuConvert_opt", menuConvert.opt)
		("menuConvert_optMesh", menuConvert.optMesh)
		("menuConvert_optPosition", menuConvert.optPosition)
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

void MainSummaryMesh::Report(const MeshData &data, int id) {
	String name = GetFileTitle(data.file);
	
	if (array.GetColumnCount() == 0)
		array.AddColumn("Param");
	if (id >= array.GetColumnCount()-1)
		array.AddColumn(Format("#%d %s", id+1, name));
	int row = 0;
	int col = id + 1;
	
	bool healing = data.mesh.healing;
	
	array.Set(row, 0, t_("File"));				array.Set(row++, col, data.file);
	array.Set(row, 0, t_("Name"));				array.Set(row++, col, name + (healing ? (x_(" ") + t_("(healed)")) : ""));
	array.Set(row, 0, t_("Format"));			array.Set(row++, col, data.GetCodeStr());	
	
	array.Set(row, 0, t_("# Panels"));			array.Set(row++, col, data.mesh.panels.GetCount());
	array.Set(row, 0, t_("# Nodes"));			array.Set(row++, col, data.mesh.nodes.GetCount());

	array.Set(row, 0, t_("Surface [m2]"));		array.Set(row++, col, FormatDouble(data.mesh.surface, 3, FD_EXP));
	array.Set(row, 0, t_("Volume [m3]"));		array.Set(row++, col, FormatDouble(data.mesh.volume, 3, FD_EXP));
	
	array.Set(row, 0, t_("Immersed surface [m2]"));array.Set(row++, col, FormatDouble(data.under.surface, 3, FD_EXP));
	array.Set(row, 0, t_("Immersed volume [m3]")); array.Set(row++, col, FormatDouble(data.under.volume, 3, FD_EXP));
	array.Set(row, 0, t_("Cb(x) [m]"));			array.Set(row++, col, FormatDouble(data.cb.x, 3, FD_EXP));
	array.Set(row, 0, t_("Cb(y) [m]"));			array.Set(row++, col, FormatDouble(data.cb.y, 3, FD_EXP));
	array.Set(row, 0, t_("Cb(z) [m]"));			array.Set(row++, col, FormatDouble(data.cb.z, 3, FD_EXP));
	
	array.Set(row, 0, t_("Water Plane Area [m2]"));	array.Set(row++, col, FormatDouble(data.waterPlaneArea, 3, FD_EXP));
	
	array.Set(row, 0, t_("Min X [m]"));			array.Set(row++, col, FormatDouble(data.mesh.env.minX, 3, FD_EXP));	
	array.Set(row, 0, t_("Max X [m]"));			array.Set(row++, col, FormatDouble(data.mesh.env.maxX, 3, FD_EXP));
	array.Set(row, 0, t_("Min Y [m]"));			array.Set(row++, col, FormatDouble(data.mesh.env.minY, 3, FD_EXP));
	array.Set(row, 0, t_("Max Y [m]"));			array.Set(row++, col, FormatDouble(data.mesh.env.maxY, 3, FD_EXP));
	array.Set(row, 0, t_("Min Z [m]"));			array.Set(row++, col, FormatDouble(data.mesh.env.minZ, 3, FD_EXP));
	array.Set(row, 0, t_("Max Z [m]"));			array.Set(row++, col, FormatDouble(data.mesh.env.maxZ, 3, FD_EXP));

	array.Set(row++, 0, t_("Stiffness Matrix"));	
	if (data.c.size() > 0) {
		for (int i = 0; i < 6; ++i) {
			for (int j = 0; j < 6; ++j) {
				if (!Hydro::C_units(i, j).IsEmpty()) {
					array.Set(row, 0, Format(t_("K(%d,%d) [%s]"), i+1, j+1, Hydro::C_units(i, j)));	array.Set(row++, col, Format("%12E", data.c(i, j)));		
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

void MainView::Init(const WithMenuPlotMesh<StaticRect> &_menuPlot) {
	CtrlLayout(*this);
	menuPlot = &_menuPlot;
	
	gl.SetEnv(env);
	gl.WhenPaint = THISBACK(OnPaint);	
}
	
static Color GetColor(int id) {
	static Color col[10] = {Color(93,185,86), Color(180,96,189), Color(176,177,67), Color(106,126,205),
							Color(219,142,68), Color(79,186,172), Color(203,83,67), Color(83,133,68),
							Color(199,89,128), Color(147,119,58)};
	int mod = id % 10;
	return col[mod];
}

void MainView::OnPaint() {
	if (~GetMenuPlot().showAxis) 
		gl.PaintAxis(0, 0, 0, env.Max()/4.);	
	
	if (~GetMenuPlot().showLimits) 
		gl.PaintCuboid(Point3D(env.maxX, env.maxY, env.maxZ), Point3D(env.minX, env.minY, env.minZ), Gray());
	
	for (int i = 0; i < ma().bem.surfs.GetCount(); ++i) {
		bool showNormals = ~GetMenuPlot().showNormals && ~GetMenuPlot().showMesh;
		bool showNormalsUnderwater = ~GetMenuPlot().showNormals && ~GetMenuPlot().showUnderwater;
		if (~GetMenuPlot().showNormals && !~GetMenuPlot().showMesh && !~GetMenuPlot().showUnderwater)
			showNormals = true;
		gl.PaintSurface(ma().bem.surfs[i].mesh, GetColor(i), ~GetMenuPlot().showMesh, 	
			showNormals, false, ~GetMenuPlot().showSkewed, 
			~GetMenuPlot().showFissure, ~GetMenuPlot().showMultiPan);
		gl.PaintSurface(ma().bem.surfs[i].under, GetColor(i), ~GetMenuPlot().showUnderwater, 
			showNormalsUnderwater, ~GetMenuPlot().showWaterLevel, false, 
			false, false);
		if (~GetMenuPlot().showCb) {
			double len = env.Max()/20;
			gl.PaintDoubleAxis(ma().bem.surfs[i].cb, len, LtBlue());
			gl.PaintCube(ma().bem.surfs[i].cb, len/10, LtGray());
		}
		if (~GetMenuPlot().showCg) {
			double len = env.Max()/20;
			gl.PaintDoubleAxis(ma().bem.surfs[i].cg, len, Black());
			gl.PaintCube(ma().bem.surfs[i].cg, len/10, LtGray());
		}
	}
}

void MainView::CalcEnvelope() {
	env.Reset();
	for (int i = 0; i < ma().bem.surfs.GetCount(); ++i)
		env.MixEnvelope(ma().bem.surfs[i].mesh.env);
}

void MainViewData::Init() {
	Add(tab.SizePos());	
}

void MainViewData::OnAddedModel() {
	int id = ma().bem.surfs.GetCount()-1;
	
	MainViewDataEach &model = models.Add();
	model.Init(ma().bem.surfs[id]);
	tab.Add(model.SizePos(), GetFileTitle(ma().bem.surfs[id].file));
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
	if (all) {
		if (col == 3 && pmesh->mesh.panels[q].IsTriangle())
			return "-";
		else
			return pmesh->mesh.panels[q].id[col]+1;
	} else {
		if (col == 3 && pmesh->under.panels[q].IsTriangle())
			return "-";
		else
			return pmesh->under.panels[q].id[col]+1;
	}
}

void MainViewDataEach::DataSourceNodes::Init(MeshData &_mesh, int _xyz, int _origMovedUnder) {
	pmesh = &_mesh;	
	xyz = _xyz;
	origMovedUnder = _origMovedUnder;
}

Value MainViewDataEach::DataSourceNodes::Format(const Value& q) const {
	ASSERT(pmesh);
	const Point3D &p = origMovedUnder == 0 ? pmesh->mesh.nodes0[q] : (origMovedUnder == 1 ? pmesh->mesh.nodes[q] : pmesh->under.nodes[q]);
	if (xyz == -1)
		return static_cast<int>(q) + 1;
	else if (xyz == 0)
		return p.x;
	else if (xyz == 1)
		return p.y;
	else
		return p.z;
}

void MainViewDataEach::Init(MeshData &_mesh) {
	CtrlLayout(arrayFacetsAll);
	CtrlLayout(arrayFacetsAll2);
	CtrlLayout(arrayFacetsUnder);
	CtrlLayout(arrayNodesOrig);
	CtrlLayout(arrayNodesMoved);
	CtrlLayout(arrayNodesUnder);
	
	arrayFacetsAll.title.SetText(t_("Facet nodes ids"));
	arrayFacetsAll.array.MultiSelect().SetLineCy(EditField::GetStdHeight()).HeaderObject().Absolute();
	arrayFacetsAll.array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, arrayFacetsAll.array);};
	dataSourceFacetsAll.SetCount(4);
	for (int c = 0; c < 4; ++c) {
		dataSourceFacetsAll[c].Init(_mesh, c, true);
		arrayFacetsAll.array.AddRowNumColumn(Format(t_("#%d"), c+1), 60).SetConvert(dataSourceFacetsAll[c]);
	}
	
	arrayFacetsAll2.title.SetText(t_("Facet node ids"));
	arrayFacetsAll2.array.MultiSelect().SetLineCy(EditField::GetStdHeight()).HeaderObject().Absolute();
	arrayFacetsAll2.array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, arrayFacetsAll2.array);};
	dataSourceFacetsAll.SetCount(4);
	for (int c = 0; c < 4; ++c) {
		dataSourceFacetsAll[c].Init(_mesh, c, true);
		arrayFacetsAll2.array.AddRowNumColumn(Format(t_("#%d"), c+1), 60).SetConvert(dataSourceFacetsAll[c]);
	}
	
	arrayFacetsUnder.title.SetText(t_("Facet nodes ids"));
	arrayFacetsUnder.array.MultiSelect().SetLineCy(EditField::GetStdHeight()).HeaderObject().Absolute();
	arrayFacetsUnder.array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, arrayFacetsUnder.array);};
	dataSourceFacetsUnder.SetCount(4);
	for (int c = 0; c < 4; ++c) {
		dataSourceFacetsUnder[c].Init(_mesh, c, false);
		arrayFacetsUnder.array.AddRowNumColumn(Format(t_("#%d"), c+1), 60).SetConvert(dataSourceFacetsUnder[c]);
	}

	const char *xyz[] = {"x", "y", "z"};

	arrayNodesOrig.title.SetText(t_("Original node coordinates"));
	arrayNodesOrig.array.MultiSelect().SetLineCy(EditField::GetStdHeight()).HeaderObject().Absolute();
	arrayNodesOrig.array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, arrayNodesOrig.array);};
	dataSourceNodesOrig.SetCount(4);
	dataSourceNodesOrig[0].Init(_mesh, -1, 0);
	arrayNodesOrig.array.AddRowNumColumn(t_("#node"), 60).SetConvert(dataSourceNodesOrig[0]);
	for (int c = 0; c < 3; ++c) {
		dataSourceNodesOrig[c+1].Init(_mesh, c, 0);
		arrayNodesOrig.array.AddRowNumColumn(Format(t_("%s"), xyz[c]), 80).SetConvert(dataSourceNodesOrig[c+1]);
	}
	
	arrayNodesMoved.title.SetText(t_("Moved node coordinates"));
	arrayNodesMoved.array.MultiSelect().SetLineCy(EditField::GetStdHeight()).HeaderObject().Absolute();
	arrayNodesMoved.array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, arrayNodesMoved.array);};
	dataSourceNodesMoved.SetCount(4);
	dataSourceNodesMoved[0].Init(_mesh, -1, 1);
	arrayNodesMoved.array.AddRowNumColumn(t_("#node"), 60).SetConvert(dataSourceNodesMoved[0]);
	for (int c = 0; c < 3; ++c) {
		dataSourceNodesMoved[c+1].Init(_mesh, c, 1);
		arrayNodesMoved.array.AddRowNumColumn(Format(t_("%s"), xyz[c]), 80).SetConvert(dataSourceNodesMoved[c+1]);
	}
	
	arrayNodesUnder.title.SetText(t_("Node coordinates"));
	arrayNodesUnder.array.MultiSelect().SetLineCy(EditField::GetStdHeight()).HeaderObject().Absolute();
	arrayNodesUnder.array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, arrayNodesUnder.array);};
	dataSourceNodesMoved.SetCount(4);
	dataSourceNodesMoved[0].Init(_mesh, -1, 2);
	arrayNodesUnder.array.AddRowNumColumn(t_("#node"), 60).SetConvert(dataSourceNodesMoved[0]);
	for (int c = 0; c < 3; ++c) {
		dataSourceNodesMoved[c+1].Init(_mesh, c, 2);
		arrayNodesUnder.array.AddRowNumColumn(Format(t_("%s"), xyz[c]), 80).SetConvert(dataSourceNodesMoved[c+1]);
	}
	
	orig.Horz(arrayFacetsAll.SizePos(), arrayNodesOrig.SizePos());
	moved.Horz(arrayFacetsAll2.SizePos(), arrayNodesMoved.SizePos());
	movedUnder.Horz(arrayFacetsUnder.SizePos(), arrayNodesUnder.SizePos());	
	  					 
	tab.Add(orig.SizePos(), t_("Original"));
	tab.Add(moved.SizePos(), t_("Moved"));
	tab.Add(movedUnder.SizePos(), t_("Moved underwater"));
	Add(tab.SizePos());	
	
	OnRefresh();
}

void MainViewDataEach::OnRefresh() {
	const MeshData &mesh = dataSourceFacetsAll[0].GetMesh();
	int num;
	
	num = mesh.mesh.panels.GetCount();
	arrayFacetsAll.array.SetVirtualCount(num);
	arrayFacetsAll.array.Refresh();
	arrayFacetsAll.numRows.SetText(FormatInt(num));
	arrayFacetsAll2.array.SetVirtualCount(num);
	arrayFacetsAll2.array.Refresh();
	arrayFacetsAll2.numRows.SetText(FormatInt(num));
		
	num = mesh.under.panels.GetCount();
	arrayFacetsUnder.array.SetVirtualCount(num);
	arrayFacetsUnder.array.Refresh();
	arrayFacetsUnder.numRows.SetText(FormatInt(num));
	
	num = mesh.mesh.nodes0.GetCount();
	arrayNodesOrig.array.SetVirtualCount(num);
	arrayNodesOrig.array.Refresh();
	arrayNodesOrig.numRows.SetText(FormatInt(num));
	
	num = mesh.mesh.nodes.GetCount();
	arrayNodesMoved.array.SetVirtualCount(num);
	arrayNodesMoved.array.Refresh();
	arrayNodesMoved.numRows.SetText(FormatInt(num));
	
	num = mesh.under.nodes.GetCount();
	arrayNodesUnder.array.SetVirtualCount(num);
	arrayNodesUnder.array.Refresh();
	arrayNodesUnder.numRows.SetText(FormatInt(num));
}

