#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <CtrlScroll/CtrlScroll.h>

#include <BEMRosetta/BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;

#define IMAGECLASS Img2
#define IMAGEFILE <BEMRosetta/BEMRosetta/main.iml>
#include <Draw/iml.h>

#include "main.h"
#include "clip.brc"

void MainBEM::Init() {
	CtrlLayout(*this);
	
	mbm(this);
	
	CtrlLayout(menuOpen);
	menuOpen.file.WhenChange = THISBACK(OnLoad);
	menuOpen.file.BrowseRightWidth(40).UseOpenFolder().BrowseOpenFolderWidth(10);
	menuOpen.butLoad.WhenAction = [&] {menuOpen.file.DoGo();};
	
	menuOpen.arrayModel.NoHeader().NoVertGrid().AutoHideSb().MultiSelect();
	menuOpen.arrayModel.AddColumn("", 20);	
	menuOpen.arrayModel.AddColumn("", 20);
	
	menuOpen.butRemove.Disable();	
	menuOpen.butRemove.WhenAction = THISBACK(OnRemove);
	menuOpen.butRemoveSelected.Disable();	
	menuOpen.butRemoveSelected.WhenAction = THISBACK(OnRemoveSelected);
	menuOpen.butJoin.Disable();	
	menuOpen.butJoin.WhenAction = THISBACK(OnJoin);
	
	CtrlLayout(menuConvert);
	menuConvert.file.WhenChange = THISBACK(OnConvert);
	menuConvert.file.BrowseRightWidth(40).UseOpenFolder(true).BrowseOpenFolderWidth(10);
	menuConvert.butLoad.WhenAction = [&] {OnConvert();};

	menuConvert.arrayModel.NoHeader().NoVertGrid().AutoHideSb();
	menuConvert.arrayModel.AddColumn("", 20);	
	menuConvert.arrayModel.AddColumn("", 20);
	
	menuConvert.opt.WhenAction = [&] {OnOpt();};
	
	OnOpt();
	
	CtrlLayout(menuPlot);
	menuPlot.butZoomToFit.WhenAction = [&] {GetSelScatter().ZoomToFit(true, true);};
	menuPlot.autoFit.WhenAction 	 = [&] {LoadSelTab(ma().bem);};
	menuPlot.opwT.WhenAction 	 	 = [&] {LoadSelTab(ma().bem);};
	menuPlot.showPoints.WhenAction 	 = [&] {LoadSelTab(ma().bem);};
	menuPlot.showPhase.WhenAction 	 = [&] {LoadSelTab(ma().bem);};
	menuPlot.showNdim.WhenAction 	 = [&] {LoadSelTab(ma().bem);};
	
	OnOpt();
	
	menuFOAMM.Init(*this, ma().bem);
	
	OnOpt();
		
	menuTab.Add(menuOpen.SizePos(), 	t_("Load"));
	menuTab.Add(menuConvert.SizePos(), 	t_("Save as")).Disable();
	menuTab.Add(menuPlot.SizePos(), 	t_("Plot")).Disable();
	if (ma().bem.experimentalFOAMM) 
		menuTab.Add(menuFOAMM.SizePos(), t_("FOAMM State Space")).Disable();
	
	menuTab.WhenSet = [&] {
		if (menuTab.IsAt(menuFOAMM)) {
			mainTab.Hide();
			TabCtrl::Item& tabMenuPlot = menuTab.GetItem(menuTab.Find(menuPlot));
			tabMenuPlot.Enable(true);
			tabMenuPlot.Text(t_("Plot"));
		} else
		 	mainTab.Show();
	};

	mainTab.WhenSet = [&] {
		bool plot = true, convertProcess = true;
		if (ma().bem.hydros.IsEmpty())
			plot = convertProcess = false;
		else if (mainTab.IsAt(mainStiffness)) {
			plot = false;
			mainStiffness.Load(ma().bem.hydros);
		} else if (mainTab.IsAt(mainA)) {
			mainA.Load(ma().bem);
			menuPlot.showPhase.Enable(false);
		} else if (mainTab.IsAt(mainB)) {
			mainB.Load(ma().bem);
			menuPlot.showPhase.Enable(false);
		} else if (mainTab.IsAt(mainForceSC)) {
			mainForceSC.Load(ma().bem);
			menuPlot.showPhase.Enable(true);
		} else if (mainTab.IsAt(mainForceFK)) {
			mainForceFK.Load(ma().bem);
			menuPlot.showPhase.Enable(true);
		} else if (mainTab.IsAt(mainForceEX)) {
			mainForceEX.Load(ma().bem);
			menuPlot.showPhase.Enable(true);
		} else if (mainTab.IsAt(mainRAO)) {
			mainRAO.Load(ma().bem);
			menuPlot.showPhase.Enable(true);
		} else if (mainTab.IsAt(mainStateSpace)) {
			mainStateSpace.Load(ma().bem);
			menuPlot.showPhase.Enable(true);
		} else if (menuTab.IsAt(menuFOAMM)) 
			menuPlot.showPhase.Enable(false);
		else {
			plot = false;
			menuPlot.showPhase.Enable(false);
		}
		TabCtrl::Item& tabMenuPlot = menuTab.GetItem(menuTab.Find(menuPlot));
		tabMenuPlot.Enable(plot);
		TabCtrl::Item& tabMenuConvert = menuTab.GetItem(menuTab.Find(menuConvert));
		tabMenuConvert.Enable(convertProcess);
		if (plot) {
			tabMenuPlot.Text(t_("Plot"));
			menuTab.Set(menuPlot);
		} else {
			tabMenuPlot.Text("");
			menuTab.Set(menuOpen);
		}		
		if (convertProcess) {
			tabMenuConvert.Text(t_("Save as"));
		} else {
			tabMenuConvert.Text("");
		}
		
		if (ma().bem.experimentalFOAMM) {
			TabCtrl::Item& tabMenuFOAMM = menuTab.GetItem(menuTab.Find(menuFOAMM));
			tabMenuFOAMM.Enable(convertProcess);
			if (convertProcess) 
				tabMenuFOAMM.Text(t_("FOAMM State Space"));
			else 
				tabMenuFOAMM.Text("");
		}
	};
	mainTab.WhenSet();
	
	mainSummary.Init();
	mainTab.Add(mainSummary.SizePos(), t_("Summary"));
	
	mainArrange.Init();
	mainTab.Add(mainArrange.SizePos(), t_("Arrange DOF")).Disable();
	
	mainStiffness.Init();
	mainTab.Add(mainStiffness.SizePos(), t_("K")).Disable();
	
	mainA.Init(DATA_A);
	mainTab.Add(mainA.SizePos(), t_("A")).Disable();
	
	mainB.Init(DATA_B);
	mainTab.Add(mainB.SizePos(), t_("B")).Disable();

	mainForceEX.Init(DATA_FORCE_EX);
	mainTab.Add(mainForceEX.SizePos(), t_("Fex")).Disable();
	
	mainForceSC.Init(DATA_FORCE_SC);
	mainTab.Add(mainForceSC.SizePos(), t_("Fsc")).Disable();
	
	mainForceFK.Init(DATA_FORCE_FK);
	mainTab.Add(mainForceFK.SizePos(), t_("Ffk")).Disable();
	
	mainRAO.Init(DATA_RAO);
	mainTab.Add(mainRAO.SizePos(), t_("RAO")).Disable();

	mainStateSpace.Init();
	mainTab.Add(mainStateSpace.SizePos(), t_("State Space")).Disable();
}

void MainBEM::InitSerialize(bool ret) {
	if (!ret || IsNull(menuPlot.autoFit)) 
		menuPlot.autoFit = true;
	
	if (!ret || IsNull(menuPlot.opwT)) 
		menuPlot.opwT = 0;

	if (!ret || IsNull(menuPlot.showPoints)) 
		menuPlot.showPoints = true;
	
	if (!ret || IsNull(menuPlot.showPhase)) 
		menuPlot.showPhase = true;

	if (!ret || IsNull(menuPlot.showNdim)) 
		menuPlot.showNdim = false;

	if (!ret || IsNull(menuConvert.opt)) 
		menuConvert.opt = 0;
}

void MainBEM::LoadSelTab(BEMData &bem) {
	int id = mainTab.Get();
	if (id == mainTab.Find(mainStateSpace))
		mainStateSpace.Load(bem);
	else if (id == mainTab.Find(mainStiffness))
		mainStiffness.Load(bem.hydros);
	else if (id == mainTab.Find(mainSummary) || id == mainTab.Find(mainArrange))
		;
	else 
		GetSelABForce().Load(bem);
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
		return mainStateSpacePlot->scatter;
	} else {
		TabCtrl &tab = GetSelABForce().tab;
		Ctrl *ctrl = tab.GetItem(tab.Get()).GetSlave();
		if (!ctrl)
			throw Exc(t_("Object not found in GetSelScatter(2)"));
		if (typeid(MainPlot) != typeid(*ctrl))
			throw Exc(t_("Unexpected type in GetSelScatter(2)"));		
		MainPlot *mainPlot = static_cast<MainPlot*>(ctrl);
		return mainPlot->scatter;
	}
}

void MainBEM::OnOpt() {
	menuOpen.file.ClearTypes();

	menuOpen.file.Type(Format(t_("All supported BEM files (%s)"), ma().bem.bemFilesExt), ma().bem.bemFilesAst);
	menuOpen.file.AllFilesType();
	String extOpen = ToLower(GetFileExt(menuOpen.file.GetData().ToString()));
	if (extOpen.IsEmpty())
		menuOpen.file.ActiveType(0);
	else if (ma().bem.bemFilesExt.Find(extOpen) >= 0)
		menuOpen.file.ActiveType(0);
	else
		menuOpen.file.ActiveType(1);
	
	menuConvert.file.ClearTypes();
	switch (menuConvert.opt) {
	case 0:	menuConvert.file <<= ForceExtSafe(~menuConvert.file, ".1"); 	
			menuConvert.file.Type(t_("Wamit .1.3.hst file"), "*.1 *.3 *.hst");
			break;
	case 1:	menuConvert.file <<= ForceExtSafe(~menuConvert.file, ".dat"); 
			menuConvert.file.Type(t_("FAST HydroDyn file"), "*.dat");
			break;
	default:menuConvert.file.Type(t_("All converted files"), "*.1 *.3 *.hst *.dat");
			break;
	}
	String extConv = ToLower(GetFileExt(menuConvert.file.GetData().ToString()));
	if (extConv.IsEmpty())
		extConv = "-";
	if (String(".1 .3 .hst").Find(extConv) >= 0)
		menuConvert.file.ActiveType(0);
	else if (String(".dat").Find(extConv) >= 0)
		menuConvert.file.ActiveType(1);
	else
		menuConvert.file.ActiveType(2);
	
	/*menuFOAMM.file.ClearTypes();
	menuFOAMM.file.Type(t_("Maynooth COER FOAMM file *.mat"), "*.mat");
	menuFOAMM.file.AllFilesType();
	String extFOAMM = ToLower(GetFileExt(menuFOAMM.file.GetData().ToString()));
	if (extFOAMM.IsEmpty())
		menuFOAMM.file.ActiveType(0);
	else if (extFOAMM == ".mat")
		menuFOAMM.file.ActiveType(0);
	else
		menuFOAMM.file.ActiveType(1);*/	
}

bool MainBEM::OnLoad() {
	String file = ~menuOpen.file;
	return OnLoadFile(file);
}

bool MainBEM::OnLoadFile(String file) {
	try {
		WaitCursor wait;
		Progress progress(t_("Loading BEM files..."), 100); 
		
		ma().bem.Load(file, [&](String str, int _pos) {
			progress.SetText(str); 
			progress.SetPos(_pos); 
			return !progress.Canceled();
		});
		int id = ma().bem.hydros.GetCount()-1;
		HydroClass &data = ma().bem.hydros[id];
		
		data.hd().Report();
		mainSummary.Report(data.hd(), id);
		if (data.hd().Nf < 0)
			return false;
		
		mainArrange.Load(ma().bem.hydros);
		
		menuOpen.arrayModel.Add(data.hd().GetCodeStr(), data.hd().name);
		menuOpen.butRemove.Enable();
		menuOpen.butRemoveSelected.Enable();
		if (menuOpen.arrayModel.GetCount() > 1)
			menuOpen.butJoin.Enable();
		menuConvert.arrayModel.Add(data.hd().GetCodeStr(), data.hd().name);
		if (menuConvert.arrayModel.GetCursor() < 0)
			menuConvert.arrayModel.SetCursor(0);
		menuFOAMM.arrayModel.Add(data.hd().GetCodeStr(), data.hd().name);
		if (menuFOAMM.arrayModel.GetCursor() < 0)
			menuFOAMM.arrayModel.SetCursor(0);
		mainTab.GetItem(mainTab.Find(mainArrange)).Enable(true);	
		mainTab.GetItem(mainTab.Find(mainStiffness)).Enable(true);
		mainTab.GetItem(mainTab.Find(mainA)).Enable(mainA.Load(ma().bem));	
		mainTab.GetItem(mainTab.Find(mainB)).Enable(mainB.Load(ma().bem));
		mainTab.GetItem(mainTab.Find(mainForceSC)).Enable(mainForceSC.Load(ma().bem));
		mainTab.GetItem(mainTab.Find(mainForceFK)).Enable(mainForceFK.Load(ma().bem));
		mainTab.GetItem(mainTab.Find(mainForceEX)).Enable(mainForceEX.Load(ma().bem));
		mainTab.GetItem(mainTab.Find(mainRAO)).Enable(mainRAO.Load(ma().bem));
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
	ma().bem.hydros.Clear();
	
	mainSummary.Clear();
	menuOpen.arrayModel.Clear();
	menuOpen.butRemove.Disable();
	menuOpen.butRemoveSelected.Disable();
	menuOpen.butJoin.Disable();
	menuConvert.arrayModel.Clear();
	menuFOAMM.arrayModel.Clear();
	
	mainArrange.Clear();	mainTab.GetItem(mainTab.Find(mainArrange)).Disable();
	mainStiffness.Clear();	mainTab.GetItem(mainTab.Find(mainStiffness)).Disable();
	mainA.Clear();			mainTab.GetItem(mainTab.Find(mainA)).Disable();
	mainB.Clear();			mainTab.GetItem(mainTab.Find(mainB)).Disable();
	mainForceSC.Clear();	mainTab.GetItem(mainTab.Find(mainForceSC)).Disable();
	mainForceFK.Clear();	mainTab.GetItem(mainTab.Find(mainForceFK)).Disable();
	mainForceEX.Clear();	mainTab.GetItem(mainTab.Find(mainForceEX)).Disable();
	mainRAO.Clear();		mainTab.GetItem(mainTab.Find(mainRAO)).Disable();
	
	mainTab.WhenSet();
}

void MainBEM::OnRemoveSelected() {
	bool selected = false;
	for (int r = menuOpen.arrayModel.GetCount()-1; r >= 0; --r) {
		if (menuOpen.arrayModel.IsSelected(r)) {
			ma().bem.hydros.Remove(r);
			mainArrange.Remove(r);
			menuOpen.arrayModel.Remove(r);
			menuConvert.arrayModel.Remove(r);
			menuFOAMM.arrayModel.Remove(r);
			selected = true;
		}
	}
	if (!selected) {
		Exclamation(t_("No model selected"));
		return;
	}
 	mainSummary.Clear();
	for (int i = 0; i < ma().bem.hydros.GetCount(); ++i)
		mainSummary.Report(ma().bem.hydros[i].hd(), i);
	
	int numrow = menuOpen.arrayModel.GetCount();
	menuOpen.butRemove.Enable(numrow > 0);
	menuOpen.butRemoveSelected.Enable(numrow > 0);
	menuOpen.butJoin.Enable(numrow > 1);
	
	mainArrange.Load(ma().bem.hydros);	mainTab.GetItem(mainTab.Find(mainArrange)).Enable(ma().bem.hydros.GetCount() > 0);	
	mainTab.GetItem(mainTab.Find(mainStiffness)).Enable(true);
	mainTab.GetItem(mainTab.Find(mainA)).Enable(mainA.Load(ma().bem));	
	mainTab.GetItem(mainTab.Find(mainB)).Enable(mainB.Load(ma().bem));
	mainTab.GetItem(mainTab.Find(mainForceSC)).Enable(mainForceSC.Load(ma().bem));
	mainTab.GetItem(mainTab.Find(mainForceFK)).Enable(mainForceFK.Load(ma().bem));
	mainTab.GetItem(mainTab.Find(mainForceEX)).Enable(mainForceEX.Load(ma().bem));
	mainTab.GetItem(mainTab.Find(mainRAO)).Enable(mainRAO.Load(ma().bem));
	mainTab.GetItem(mainTab.Find(mainStateSpace)).Enable(mainStateSpace.Load(ma().bem));
	
	mainTab.WhenSet();
}

void MainBEM::OnJoin() {
	Vector<int> ids;
	bool selected = false;
	for (int r = menuOpen.arrayModel.GetCount()-1; r >= 0; --r) {
		if (menuOpen.arrayModel.IsSelected(r)) {
			ids << r;
			selected = true;
		}
	}
	if (!selected) {
		Exclamation(t_("No model selected"));
		return;
	}
	try {
		WaitCursor wait;
		Progress progress(t_("Joining selected BEM files..."), 100); 
		
		ma().bem.Join(ids, [&](String str, int _pos) {
			progress.SetText(str); 
			progress.SetPos(_pos); 
			return !progress.Canceled();
		});
		
		mainSummary.Clear();
		mainArrange.Clear();
		menuOpen.arrayModel.Clear();
		menuConvert.arrayModel.Clear();
		menuFOAMM.arrayModel.Clear();
		for (int i = 0; i < ma().bem.hydros.GetCount(); ++i) {
			const Hydro &data = ma().bem.hydros[i].hd();
			mainSummary.Report(data, i);
			mainArrange.Load(ma().bem.hydros);
			menuOpen.arrayModel.Add(data.GetCodeStr(), data.name);
			menuConvert.arrayModel.Add(data.GetCodeStr(), data.name);
			menuFOAMM.arrayModel.Add(data.GetCodeStr(), data.name);
		}
		if (ma().bem.hydros.GetCount() > 0) {
			menuOpen.arrayModel.SetCursor(0);
			menuConvert.arrayModel.SetCursor(0);
			menuFOAMM.arrayModel.SetCursor(0);
		}
		
		mainArrange.Load(ma().bem.hydros);	mainTab.GetItem(mainTab.Find(mainArrange)).Enable(ma().bem.hydros.GetCount() > 0);		
		mainTab.GetItem(mainTab.Find(mainStiffness)).Enable(true);
		mainTab.GetItem(mainTab.Find(mainA)).Enable(mainA.Load(ma().bem));	
		mainTab.GetItem(mainTab.Find(mainB)).Enable(mainB.Load(ma().bem));
		mainTab.GetItem(mainTab.Find(mainForceSC)).Enable(mainForceSC.Load(ma().bem));
		mainTab.GetItem(mainTab.Find(mainForceFK)).Enable(mainForceFK.Load(ma().bem));
		mainTab.GetItem(mainTab.Find(mainForceEX)).Enable(mainForceEX.Load(ma().bem));
		mainTab.GetItem(mainTab.Find(mainRAO)).Enable(mainRAO.Load(ma().bem));
		mainTab.GetItem(mainTab.Find(mainStateSpace)).Enable(mainStateSpace.Load(ma().bem));
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}
		
	int numrow = menuOpen.arrayModel.GetCount();
	menuOpen.butJoin.Enable(numrow > 1);
}

bool MainBEM::OnConvert() {
	String file = ~menuConvert.file;
	
	try {
		int id = menuConvert.arrayModel.GetCursor();
		if (id < 0) {
			Exclamation(t_("Please select a model to export"));
			return false;
		}
		Hydro::BEM_SOFT type;	
		switch (menuConvert.opt) {
		case 0:	type = Hydro::WAMIT_1_3;	break;
		case 1:	type = Hydro::FAST_WAMIT;	break;
		case 2:	type = Hydro::UNKNOWN;		break;
		default: throw Exc(t_("Unknown type in OnConvert()"));
		}
		ma().bem.hydros[id].hd().SaveAs(file, type);	
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
		return false;
	}
	return true;
}

void MainBEM::Jsonize(JsonIO &json) {
	json
		("menuOpen_file", menuOpen.file)
		("menuConvert_file", menuConvert.file)
		("menuConvert_opt", menuConvert.opt)
		("menuPlot_autoFit", menuPlot.autoFit)
		("menuPlot_opwT", menuPlot.opwT)
		("menuPlot_showPoints", menuPlot.showPoints)
		("menuPlot_showPhase", menuPlot.showPhase)
		("menuPlot_showNdim", menuPlot.showNdim)
		//("menuFOAMM_file", menuFOAMM.file)
	;
}

String MainBEM::BEMFile(String fileFolder) const {
	if (DirectoryExists(fileFolder)) {
		int bestipos = INT_MAX;
		FindFile ff(AppendFileName(fileFolder, "*.*"));
		while (ff) {
			if (ff.IsFile()) {
				int ipos = ma().bem.bemFilesExt.Find(GetFileExt(ff.GetName()));
 				if (ipos >= 0 && ipos < bestipos) {
					fileFolder = ff.GetPath();
					bestipos = ipos;	// It takes the file with most probable extension
				}
			}
			ff.Next();
		}
	}
	return fileFolder;
}

void MainBEM::DragAndDrop(Point , PasteClip& d) {
	if (IsDragAndDropSource())
		return;
	if (AcceptFiles(d)) {
		Vector<String> files = GetFiles(d);
		for (int i = 0; i < files.GetCount(); ++i) {
			String file = BEMFile(files[i]);
			menuOpen.file <<= file;
			OnLoad();
		}
	}
}

bool MainBEM::Key(dword key, int ) {
	if (key == K_CTRL_V) {
		Vector<String> files = GetFiles(Ctrl::Clipboard());
		for (int i = 0; i < files.GetCount(); ++i) {
			String file = BEMFile(files[i]);
			menuOpen.file <<= file;
			OnLoad();
		}
		return true;
	}
	return false;
}

void MainSummary::Init() {
	CtrlLayout(*this);
	array.SetLineCy(EditField::GetStdHeight()).MultiSelect();
	array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, array);};
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
	array.Set(row, 0, t_("Soft"));				array.Set(row++, col, data.GetCodeStr());
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
	if (data.w.GetCount() > 1) {
		array.Set(row, 0, t_("freq_end [rad/s]"));	array.Set(row++, col, data.w[data.w.GetCount()-1]);
		if (data.GetIrregularFreq() < 0) { 
			array.Set(row, 0, t_("freq_delta [rad/s]"));array.Set(row++, col, data.w[1] - data.w[0]);
		} else {
			String strHead;
			for (int i = 0; i < data.w.GetCount(); ++i) {
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
	if (data.head.GetCount() > 1) {
		array.Set(row, 0, t_("head_end [º]"));	array.Set(row++, col, data.head[data.head.GetCount()-1]);
		if (data.GetIrregularHead() < 0) { 
			array.Set(row, 0, t_("head_delta [º]"));array.Set(row++, col, data.head[1] - data.head[0]);
		} else {
			String strHead;
			for (int i = 0; i < data.head.GetCount(); ++i) {
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
	
	array.Set(row, 0, t_("A0 available"));		array.Set(row++, col, data.IsLoadedAw0()   ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("Ainf available"));	array.Set(row++, col, data.IsLoadedAwinf() ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("A available"));		array.Set(row++, col, data.IsLoadedA() 	   ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("B available"));		array.Set(row++, col, data.IsLoadedB() 	   ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("K available"));		array.Set(row++, col, data.IsLoadedC() 	   ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("Fex available"));		array.Set(row++, col, data.IsLoadedFex()   ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("Fsc available"));		array.Set(row++, col, data.IsLoadedFsc()   ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("Ffk available"));		array.Set(row++, col, data.IsLoadedFfk()   ? t_("Yes") : t_("No"));
	array.Set(row, 0, t_("RAO available"));		array.Set(row++, col, data.IsLoadedRAO()   ? t_("Yes") : t_("No"));
	
	array.Set(row, 0, t_("#bodies"));			array.Set(row++, col, data.Nb);
	for (int ib = 0; ib < data.Nb; ++ib) {
		String sib = Format("#%d", ib+1);
		if (data.names.GetCount() > ib) {
			sib += " " + data.names[ib];
			array.Set(row, 0, sib + " " + t_("Name"));		array.Set(row++, col, data.names[ib]);
		} else {
			array.Set(row, 0, sib + " " + t_("Name"));		array.Set(row++, col, "-");
		}
		array.Set(row, 0, sib + " " + t_("#dof"));
		if (data.dof.GetCount() > ib) 
			array.Set(row++, col, data.dof[ib]);
		else 
			array.Set(row++, col, "-");
		
		array.Set(row, 0, sib + " " + t_("Vsub [m3]"));
		if (data.Vo.size() > ib) 
			array.Set(row++, col, FormatDouble(data.Vo[ib], 3, FD_EXP));
		else 
			array.Set(row++, col, "-");
		
		if (data.cg.size() > 3*ib) {
			array.Set(row, 0, sib + " " + t_("Cg(x) [m]"));	array.Set(row++, col, FormatDouble(data.cg(0, ib), 3, FD_EXP));
			array.Set(row, 0, sib + " " + t_("Cg(y) [m]"));	array.Set(row++, col, FormatDouble(data.cg(1, ib), 3, FD_EXP));
			array.Set(row, 0, sib + " " + t_("Cg(z) [m]"));	array.Set(row++, col, FormatDouble(data.cg(2, ib), 3, FD_EXP));
		} else {
			array.Set(row, 0, sib + " " + t_("Cg(x) [m]"));	array.Set(row++, col, "-");
			array.Set(row, 0, sib + " " + t_("Cg(y) [m]"));	array.Set(row++, col, "-");
			array.Set(row, 0, sib + " " + t_("Cg(z) [m]"));	array.Set(row++, col, "-");
		}
		if (data.cb.size() > 3*ib) {
			array.Set(row, 0, sib + " " + t_("Cb(x) [m]"));	array.Set(row++, col, FormatDouble(data.cb(0, ib), 3, FD_EXP));
			array.Set(row, 0, sib + " " + t_("Cb(y) [m]"));	array.Set(row++, col, FormatDouble(data.cb(1, ib), 3, FD_EXP));
			array.Set(row, 0, sib + " " + t_("Cb(z) [m]"));	array.Set(row++, col, FormatDouble(data.cb(2, ib), 3, FD_EXP));
		} else {
			array.Set(row, 0, sib + " " + t_("Cb(x) [m]"));	array.Set(row++, col, "-");
			array.Set(row, 0, sib + " " + t_("Cb(y) [m]"));	array.Set(row++, col, "-");
			array.Set(row, 0, sib + " " + t_("Cb(z) [m]"));	array.Set(row++, col, "-");
		}
		if (data.C.GetCount() > ib) {
			for (int i = 0; i < 6; ++i) {
				for (int j = 0; j < 6; ++j) {
					if (!Hydro::C_units(i, j).IsEmpty()) {
						array.Set(row, 0, sib + " " + Format(t_("K(%d,%d) [%s]"), i+1, j+1, Hydro::C_units(i, j)));	array.Set(row++, col, FormatDouble(data.C_dim(ib, i, j), 3, FD_EXP));		
					}
				}
			}
		} else {
			for (int i = 0; i < 6; ++i) {
				for (int j = 0; j < 6; ++j) {
					if (!Hydro::C_units(i, j).IsEmpty()) {
						array.Set(row, 0, sib + " " + Format(t_("K(%d,%d) [%s]"), i+1, j+1, Hydro::C_units(i, j)));	array.Set(row++, col, "-");		
					}
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

void MainArrange::Load(Upp::Array<HydroClass> &hydro) {
	MainArrange::Clear();
	for (int i = 0; i < hydro.GetCount(); ++i) {
		ArrangeDOF &arr = arrangeDOF.Add();
		arr.Init(hydro[i].hd());
		tab.Add(arr.SizePos(), Format("[%s] %s", hydro[i].hd().GetCodeStr(), hydro[i].hd().name));
	}
}

void MainArrange::Remove(int c) {
	tab.Remove(c);
	arrangeDOF.Remove(c);
}

MainBEM &mbm(MainBEM *m) {
	static MainBEM *mp = 0;
	if (m)
		mp = m;
	return *mp;
}

void MenuFOAMM::Init(MainBEM &mainBEM, BEMData &bem) {
	CtrlLayout(*this);
	
	butLoad.WhenAction 	= [&] {
		if (OnFOAMM()) {
			int id = mainBEM.mainTab.Find(mainBEM.mainStateSpace);
			mainBEM.mainTab.GetItem(id).Enable(true);
			int idPlot = mainBEM.menuTab.Find(mainBEM.menuPlot);
			mainBEM.menuTab.GetItem(idPlot).Enable(true);
			mainBEM.mainTab.Set(0);
			mainBEM.mainTab.Set(id);
		}
	};
	
	arrayModel.NoHeader().NoVertGrid().AutoHideSb();
	arrayModel.AddColumn("", 20);	
	arrayModel.AddColumn("", 20);	
	arrayModel.WhenSel = [&] {WhenSelArrayModel(bem);};
	
	arrayCases.SetLineCy(EditField::GetStdHeight());
	arrayCases.AddColumn(t_("Selected"), 30);
	arrayCases.AddColumn(t_("Body"), 20);
	arrayCases.AddColumn(t_("Row"), 20);	
	arrayCases.AddColumn(t_("Column"), 20);
	arrayCases.AddColumn(t_("From (rad/s)"), 40);
	arrayCases.AddColumn(t_("To (rad/s)"), 40);
	arrayCases.AddColumn(t_("Frequencies (rad/s)"), 60);
	arrayCases.WhenSel = [&] {WhenSelArrayCases();};
	
	fromFreq.WhenAction = THISBACK(WhenArrayCases);
	fromFreq.SetRectangle(rectFrom, THISBACK(WhenFocus));
	toFreq.WhenAction   = THISBACK(WhenArrayCases);
	toFreq.SetRectangle(rectTo, THISBACK(WhenFocus));
	
	select.SetRectangle(rectSelected, THISBACK(WhenFocus));
	addFreq.SetRectangle(rectSelected, THISBACK(WhenFocus));
	removeFreq.SetRectangle(rectSelected, THISBACK(WhenFocus));
	select.WhenAction = THISBACK(WhenArrayCases);
	addFreq.WhenAction 	= [&] {
			select.Add();
			WhenArrayCases();			
		};
	addFreq.SetRectangle(rectSelected, THISBACK(WhenFocus));
	removeFreq.WhenAction = [&] {
			select.Remove();
			WhenArrayCases();
		};
	removeFreq.SetRectangle(rectSelected, THISBACK(WhenFocus));
	
	plotsReal.Init();
	plotsImag.Init();
	splitter.Vert(plotsReal.SizePos(), plotsImag.SizePos());
	
	plotsReal.scatter.WhenPainter = THISBACK(OnPainter);
	plotsReal.scatter.WhenLeftDown = THISBACK(OnLeftDown);
	
	foammLogo.Set(Img2::FOAMM());
	foammWorking.LoadBuffer(String(animatedStar, animatedStar_length));
	foammWorking.Hide();
	status.Hide();
	progress.Hide();
	butCancel.Hide();
	butCancel.WhenAction = [&] {isCancelled = true;};
}

void MenuFOAMM::WhenFocus(StaticRectangle *rect) {
	rectFrom.Hide();
	rectTo.Hide();
	rectSelected.Hide();
	if (rect)
		rect->Show();
	rectActual = rect;
}

void MenuFOAMM::OnPainter(Painter &w) {
	int plotW = plotsReal.scatter.GetPlotWidth(), plotH = plotsReal.scatter.GetPlotHeight();
	
	for (int i = 0; i < select.GetCount(); ++i) {
		double xFreq = plotsReal.scatter.GetPosX(ScanDouble(AsString(select.Get(i))));
		DrawLineOpa(w, xFreq, 0, xFreq, plotH, 1, 1, 2, LtBlue(), "2 2");
	}
	if (!IsNull(fromFreq)) {
		double xFrom = plotsReal.scatter.GetPosX(~fromFreq);
		FillRectangleOpa(w, 0, 0, xFrom, plotH, 0.5, Null, LtBlue());
	}
	if (!IsNull(toFreq)) {
		double xTo = plotsReal.scatter.GetPosX(~toFreq);
		FillRectangleOpa(w, xTo, 0, plotW, plotH, 0.5, Null, LtBlue());
	}
}

void MenuFOAMM::OnLeftDown(Point p) {
	double freq = plotsReal.scatter.GetRealPosX(p.x);

	if (!plotsReal.scatter.IsEmpty()) {
		DataSource &data = plotsReal.scatter.GetDataSource(0);
		freq = data.CloserX(freq);
	}
	
	if (rectActual == &rectFrom) {	
		fromFreq <<= freq;
		fromFreq.WhenAction();
	} else if (rectActual == &rectTo) {	
		toFreq <<= freq;
		toFreq.WhenAction();
	} else if (rectActual == &rectSelected) {	
		select.Set(freq);
		select.WhenAction();
	}
}
	
void MenuFOAMM::WhenSelArrayModel(BEMData &bem) {
	arrayCases.Clear();
	options.Clear();
	
	int id = arrayModel.GetCursor();
	if (id < 0)
		return;
	
	ASSERT(id < ma().bem.hydros.GetCount());
	
	const Hydro &hydro = ma().bem.hydros[id].hd();
	
	for (int ib = 0; ib < hydro.Nb; ++ib) {
		for (int idof = 0; idof < 6; ++idof) {
			for (int jdof = 0; jdof < 6; ++jdof) {
				if (!bem.onlyDiagonal || idof == jdof) {
					int idf = ib*6 + idof;
					int jdf = ib*6 + jdof;
	
					if (hydro.IsLoadedA() && hydro.IsLoadedB() && !IsNull(hydro.A[0](idf, jdf)) && !IsNull(hydro.B[0](idf, jdf))) {
						arrayCases.Add(false, ib+1, Hydro::StrDOFAbrev_base(idof), Hydro::StrDOFAbrev_base(jdof));
						int row = arrayCases.GetCount()-1;
						arrayCases.SetCtrl(row, 0, options.Add());
						options.Top() << [=] {options[row].SetFocus();};
					}
				}
			}
		}
	}
	if (arrayCases.GetCount() > 0)
		arrayCases.SetCursor(0);
}

void MenuFOAMM::WhenSelArrayCases() {
	try {
		int row = arrayCases.GetCursor();
		if (row < 0)
			return;
	
		bool opChoose = arrayCases.Get(row, 0);
		fromFreq <<= arrayCases.Get(row, 4);
		toFreq   <<= arrayCases.Get(row, 5);
	
		String freqs = arrayCases.Get(row, 6);
		Vector<String> afreqs = Split(freqs, ';');
		select.Clear();
		for (int i = 0; i < afreqs.GetCount(); ++i)
			select.Add(ScanDouble(afreqs[i]));	
		
		if (opChoose)
			ma().Status(Check(~fromFreq, ~toFreq, ~freqs));
		
		int id = arrayModel.GetCursor();
		if (id < 0)
			return;
		
		const Hydro &hydro = ma().bem.hydros[id].hd();
		
		int ib = int(arrayCases.Get(row, 1)) - 1;
		int idof = Hydro::DOFStrAbrev(arrayCases.Get(row, 2));
		int jdof = Hydro::DOFStrAbrev(arrayCases.Get(row, 3));
		
		plotsReal.Init(idof + 6*ib, jdof + 6*ib, DATA_STS_MA);
		plotsReal.Load(hydro);

		plotsImag.Init(idof + 6*ib, jdof + 6*ib, DATA_STS_PH);
		plotsImag.Load(hydro);


	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
		return;
	}
}

void MenuFOAMM::WhenArrayCases() {
	int row = arrayCases.GetCursor();
	if (row < 0)
		return;
	
	arrayCases.Set(row, 4, ~fromFreq);
	arrayCases.Set(row, 5, ~toFreq);
	
	String freqs;
	for (int i = 0; i < select.GetCount(); ++i) {
		if (!freqs.IsEmpty())
			freqs << ";";
		double freq = select.Get(i);
		if (!plotsReal.scatter.IsEmpty()) {
			DataSource &data = plotsReal.scatter.GetDataSource(0);
			freq = data.CloserX(freq);
		}
		freqs << freq;
	}
	arrayCases.Set(row, 6, freqs);
	
	//if (~opChoose)
		ma().Status(Check(~fromFreq, ~toFreq, ~freqs));
	
	plotsReal.scatter.Refresh();
}

String MenuFOAMM::Check(double fromFreq, double toFreq, String freqs) {
	if (IsNull(fromFreq))
		return t_("From: frequency is empty");

	if (IsNull(toFreq))
		return t_("To: frequency is empty");

	if (toFreq <= fromFreq) 
		return t_("From: frequency has to be lower than To: frequency");
			
	Vector<String> afreqs = Split(freqs, ';');
	if (afreqs.IsEmpty()) 
		return t_("No frequency has been selected");
	
	Vector<double> unique;
	for (int i = 0; i < afreqs.GetCount(); ++i) {
		double freq = ScanDouble(afreqs[i]);
		if (freq < fromFreq || freq > toFreq) 
			return t_("Selected frequencies have to be between lower and higher limits");
		FindAddRatio(unique, freq, 0.001);
	}
	if (unique.GetCount() != afreqs.GetCount()) 
		return t_("Some selected frequencies are repeated");
	
	return String("");
}

bool MenuFOAMM::OnFOAMM() {
	Vector<int> ibs, idofs, jdofs;
	Vector<double> froms, tos;
	Vector<Vector<double>> freqs;
	String ret;
	
	try {
		int id = arrayModel.GetCursor();
		if (id < 0)
			return false;
		for (int row = 0; row < arrayCases.GetCount(); ++row) {
			bool proc = arrayCases.Get(row, 0);
			if (proc) {
				int ib = int(arrayCases.Get(row, 1))-1;
				ibs << ib;
				String sidf = arrayCases.Get(row, 2);
				idofs << Hydro::DOFStrAbrev(sidf);
				String sjdf = arrayCases.Get(row, 3);
				jdofs << Hydro::DOFStrAbrev(sjdf);
				double from = arrayCases.Get(row, 4);
				double to = arrayCases.Get(row, 5);
				String strfreqs = arrayCases.Get(row, 6);
				String err = Check(fromFreq, toFreq, strfreqs);
				if (!err.IsEmpty()) {
					Exclamation(Format(t_("Problem in body %d (%s, %s): %s"), ib+1, sidf, sjdf, err));
					return false;		
				}
				froms << from;
				tos << to;
				Vector<double> &f = freqs.Add();
				Vector<String> fs = Split(strfreqs, ';');
				for (int i = 0; i < fs.GetCount(); ++i)
					f << ScanDouble(fs[i]);
			}
		}
		if (ibs.IsEmpty()) {
			Exclamation(t_("No case has been selected"));
			return false;			
		}
		foammWorking.Show();
		foammWorking.Play();
		status.Show();
		progress.Show();
		butCancel.Show();
		butLoad.Disable();
		WaitCursor wait;
		isCancelled = false;
		status.SetText(t_("Starts processing"));
		Foamm &foamm = static_cast<Foamm&>(ma().bem.hydros[id]);
		foamm.Get(ibs, idofs, jdofs, froms, tos, freqs,
			[&](String str, int pos)->bool {
				if (!str.IsEmpty())
					status.SetText(str);	
				if (IsNull(pos))
					; 
				else
					progress.Set(pos, 100);
				ProcessEvents(); 
				return isCancelled;
			}, 
			[&](String str) {
				if (!str.IsEmpty()) {
					str.Replace("\r", "");
					str.Replace("\n\n", "\n");
					Exclamation(t_("FOAMM message:&") + DeQtfLf(str));
				}
				ProcessEvents(); 
			});
	} catch (Exc e) {
		ret = DeQtfLf(e);
	}
	foammWorking.Hide();
	foammWorking.Stop();
	status.Hide();
	progress.Hide();
	butCancel.Hide();
	butLoad.Enable();
	if (!ret.IsEmpty()) {
		Exclamation(ret);
		return false;
	}
	return true;
}
