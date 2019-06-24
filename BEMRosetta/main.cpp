#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas.h>

#include <BEMRosetta/BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;

#define IMAGECLASS Img
#define IMAGEFILE <BEMRosetta/BEMRosetta/main.iml>
#include <Draw/iml.h>

#define TOPICFILE <BEMRosetta/BEMRosetta/main.tpp/all.i>
#include <Core/topic_group.h>

#include "main.h"

void Main::Init() {
	CtrlLayout(*this, "BEMRosetta");
	Sizeable().Zoomable().SetMinSize(Size(800, 600));
	Icon(Img::Rosetta64());
	LargeIcon(Img::Rosetta256());
	ma(this);
	
	bool firstTime = false;
	if (!bem.LoadSerializeJson()) {
		firstTime = true;
		Cout() << "\n" << t_("BEM configuration data is not loaded. Defaults are set");
	}
	if (!LoadSerializeJson()) {
		firstTime = true;
		Cout() << "\n" << t_("Configuration data is not loaded. Defaults are set");
	}
	
	CtrlLayout(menuOpen);
	menuOpen.file <<= THISBACK(OnLoad);
	menuOpen.file.BrowseRightWidth(40).UseOpenFolder().BrowseOpenFolderWidth(10);
	menuOpen.butLoad <<= THISBACK(OnLoad);
	
	//menuOpen.arrayModel.SetLineCy(EditField::GetStdHeight());
	menuOpen.arrayModel.NoHeader().NoVertGrid().AutoHideSb();
	menuOpen.arrayModel.AddColumn("", 20);	
	menuOpen.arrayModel.AddColumn("", 20);
	
	menuOpen.butRemove.Disable();	
	menuOpen.butRemove.WhenAction = [&] {
		bem.hydros.Clear();
		mainSummary.Clear();
		menuOpen.arrayModel.Clear();
		menuOpen.butRemove.Disable();
		menuConvert.arrayModel.Clear();	
		mainArrange.Clear();	mainTab.GetItem(mainTab.Find(mainArrange)).Disable();
		mainA.Clear();			mainTab.GetItem(mainTab.Find(mainA)).Disable();
		mainB.Clear();			mainTab.GetItem(mainTab.Find(mainB)).Disable();
		mainForceSC.Clear();	mainTab.GetItem(mainTab.Find(mainForceSC)).Disable();
		mainForceFK.Clear();	mainTab.GetItem(mainTab.Find(mainForceFK)).Disable();
		mainForceEX.Clear();	mainTab.GetItem(mainTab.Find(mainForceEX)).Disable();
		mainRAO.Clear();		mainTab.GetItem(mainTab.Find(mainRAO)).Disable();
	};
	
	CtrlLayout(menuConvert);
	menuConvert.file <<= THISBACK(OnConvert);
	menuConvert.file.BrowseRightWidth(40).UseOpenFolder(true).BrowseOpenFolderWidth(10);
	menuConvert.butLoad <<= THISBACK(OnConvert);

	//menuConvert.arrayModel.SetLineCy(EditField::GetStdHeight());
	menuConvert.arrayModel.NoHeader().NoVertGrid().AutoHideSb();
	menuConvert.arrayModel.AddColumn("", 20);	
	menuConvert.arrayModel.AddColumn("", 20);
	
	menuConvert.opt.WhenAction = [&] {OnOpt();};
	
	OnOpt();
	
	CtrlLayout(menuPlot);
	menuPlot.butZoomToFit.WhenAction = [&] {GetSelPlot().scatter.ZoomToFit(true, true);};
	menuPlot.autoFit.WhenAction 	 = [&] {LoadSelTab(bem.hydros);};
	menuPlot.opwT.WhenAction 	 	 = [&] {LoadSelTab(bem.hydros);};
	menuPlot.showPoints.WhenAction 	 = [&] {LoadSelTab(bem.hydros);};
	menuPlot.showPhase.WhenAction 	 = [&] {LoadSelTab(bem.hydros);};
	menuPlot.showAdim.WhenAction 	 = [&] {LoadSelTab(bem.hydros);};
	
	CtrlLayout(menuView);
	menuView.file <<= THISBACK(OnView);
	menuView.file.BrowseRightWidth(40).UseOpenFolder(true).BrowseOpenFolderWidth(10);
	menuView.butLoad <<= THISBACK(OnView);

	menuView.butRemove.WhenAction = [&] {
		WaitCursor waitcursor;
		
		surfs.Clear();
	};
	
	menuOptions.Init(bem);
	menuOptions.Load();
	
	menuAbout.Init();
	
	menuTab.Add(menuOpen.SizePos(), 	t_("Open"));
	menuTab.Add(menuConvert.SizePos(), 	t_("Convert"));
	menuTab.Add(menuPlot.SizePos(), 	t_("Plot")).Disable();
	menuTab.Add(menuView.SizePos(), 	t_("Mesh"));
	menuTab.Add(menuOptions.SizePos(), 	t_("Options"));
	menuTab.Add(menuAbout.SizePos(), 	t_("About"));
	
	menuTab.WhenSet = [&] {
		if (menuTab.IsAt(menuAbout)) 
			mainTab.Hide();
		else if (menuTab.IsAt(menuOptions)) {
			menuOptions.Load();
			mainTab.Show();
		} else {
			if (menuTab.IsAt(menuView)) { 
				mainTab.Set(mainView);
				menuTab.Set(menuView);
			}
			mainTab.Show();
		}
		if (!menuTab.IsAt(menuOptions) && menuOptions.IsChanged()) {
			if (PromptYesNo(t_("Options have changed&Do you want to save them?")))
				menuOptions.OnSave();
			else
				menuOptions.Load();
		}
	};	
	
	mainTab.WhenSet = [&] {
		bool plot = true;
		if (mainTab.IsAt(mainA)) {
			mainA.Load(bem.hydros);
			menuPlot.showPhase.Enable(false);
		} else if (mainTab.IsAt(mainB)) {
			mainB.Load(bem.hydros);
			menuPlot.showPhase.Enable(false);
		} else if (mainTab.IsAt(mainForceSC)) {
			mainForceSC.Load(bem.hydros);
			menuPlot.showPhase.Enable(true);
		} else if (mainTab.IsAt(mainForceFK)) {
			mainForceFK.Load(bem.hydros);
			menuPlot.showPhase.Enable(true);
		} else if (mainTab.IsAt(mainForceEX)) {
			mainForceEX.Load(bem.hydros);
			menuPlot.showPhase.Enable(true);
		} else if (mainTab.IsAt(mainRAO)) {
			mainRAO.Load(bem.hydros);
			menuPlot.showPhase.Enable(true);
		} else if (mainTab.IsAt(mainView)) { 
			plot = false;
			menuTab.Set(menuView);	
			mainTab.Set(mainView);	 
		} else if (mainTab.IsAt(mainStateSpace)) {
			mainStateSpace.Load(bem.hydros);
			menuPlot.showPhase.Enable(true);
		} else {
			plot = false;
			menuPlot.showPhase.Enable(false);
		}
		TabCtrl::Item& plotIt = menuTab.GetItem(menuTab.Find(menuPlot));
		plotIt.Enable(plot);
		if (plot) {
			plotIt.Text("Plot");
			menuTab.Set(menuPlot);
		} else
			plotIt.Text("");
	};
	
	mainSummary.Init();
	mainTab.Add(mainSummary.SizePos(), t_("Summary"));
	
	mainArrange.Init();
	mainTab.Add(mainArrange.SizePos(), t_("Arrange DOF")).Disable();
	
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
		
	mainView.Init();
	mainTab.Add(mainView.SizePos(), t_("View"));
		
	mainOutput.Init();
	mainTab.Add(mainOutput.SizePos(), t_("Output"));
	
	Hydro::Print 	  	= [this](String s) {printf("%s", ~s); mainOutput.Print(s);};
	Hydro::PrintWarning = [this](String s) {printf("%s", ~s); mainOutput.Print(s);};
	Hydro::PrintError 	= [this](String s) {printf("%s", ~s); mainOutput.Print(s); mainTab.Set(mainOutput);};
	
	if (firstTime)
		menuTab.Set(menuAbout);
}

void Main::LoadSelTab(Upp::Array<HydroClass> &hydros) {
	if (mainTab.Get() == mainTab.Find(mainStateSpace))
		mainStateSpace.Load(hydros);
	else 
		GetSelTab().Load(hydros);
}

MainABForce &Main::GetSelTab() {
	return *(static_cast<MainABForce*>(mainTab.GetItem(mainTab.Get()).GetSlave()));
}

MainPlot &Main::GetSelPlot() {
	TabCtrl &tab = GetSelTab().tab;
	return *(static_cast<MainPlot*>(tab.GetItem(tab.Get()).GetSlave()));
}

static String ForceExtSafe(String fileName, String ext) {
	if (fileName.IsEmpty())
		return String();
	return ForceExt(fileName, ext);
}

void Main::OnOpt() {
	menuOpen.file.ClearTypes();
	
//	menuOpen.file.Type("Wamit .1.3.4.hst file", "*.1 *.3 *.4 *.hst");
//	menuOpen.file.Type("Wamit .out file", "*.out");	
//	menuOpen.file.Type("FAST HydroDyn file", "*.dat");	
//	menuOpen.file.Type("Nemoh .cal file", "*.cal");	
//	menuOpen.file.Type("SeaFEM .flavia.inf file", "*.inf");
//	menuOpen.file.Type("AQWA .AH1.LIS file", "*.ah1 *.lis");
//	menuOpen.file.Type("FOAMM .mat file", "*.mat");
	menuOpen.file.Type("All supported BEM files", "*.1 *.3 *.hst *.4 *.out *.dat *.cal *.inf *.ah1 *.lis *.mat");
	menuOpen.file.AllFilesType();
	String extOpen = ToLower(GetFileExt(menuOpen.file.GetData().ToString()));
	if (extOpen.IsEmpty())
		extOpen = "-";
	if (String(".1 .3 .4 .hst").Find(extOpen) >= 0)
		menuOpen.file.ActiveType(0);
	else if (String(".out").Find(extOpen) >= 0)
		menuOpen.file.ActiveType(1);
	else if (String(".dat").Find(extOpen) >= 0)
		menuOpen.file.ActiveType(2);
	else if (String(".cal").Find(extOpen) >= 0)
		menuOpen.file.ActiveType(3);
	else if (String(".inf").Find(extOpen) >= 0)
		menuOpen.file.ActiveType(4);
	else if (String(".ah1 .lis").Find(extOpen) >= 0)
		menuOpen.file.ActiveType(5);
	else if (String(".mat").Find(extOpen) >= 0)
		menuOpen.file.ActiveType(6);
	else
		menuOpen.file.ActiveType(7);
	
	menuConvert.file.ClearTypes();
	switch (menuConvert.opt) {
	case 0:	menuConvert.file <<= ForceExtSafe(~menuConvert.file, ".1"); 	
			menuConvert.file.Type("Wamit .1.3.hst file", "*.1 *.3 *.hst");
			break;
	case 1:	menuConvert.file <<= ForceExtSafe(~menuConvert.file, ".dat"); 
			menuConvert.file.Type("FAST HydroDyn file", "*.dat");
			break;
	default:menuConvert.file.Type("All converted files", "*.1 *.3 *.hst *.dat");
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
	
	menuView.file.ClearTypes(); 

	menuView.file.Type("Wamit .gdf file", "*.gdf");	
	menuView.file.Type("Nemoh .dat and Wamit panel .dat file", "*.dat");
	menuView.file.Type("All supported mesh files", "*.gdf *.dat");
	menuView.file.AllFilesType();
	String extView = ToLower(GetFileExt(menuView.file.GetData().ToString()));
	if (extView.IsEmpty())
		extView = "-";
	if (String(".gdf").Find(extView) >= 0)
		menuView.file.ActiveType(0);
	else if (String(".dat").Find(extView) >= 0)
		menuView.file.ActiveType(1);
	else
		menuView.file.ActiveType(2);
}

void Main::OnLoad() {
	String file = ~menuOpen.file;
	
	try {
		WaitCursor wait;
		Progress progress(t_("Loading BEM files..."), 100); 
		
		bem.Load(file, [&](String str, int pos) {progress.SetText(str); progress.SetPos(pos);});
	
		String strError;
		if (!HydroClass::MatchCoeffStructure(bem.hydros, strError)) {
			int num = bem.hydros.GetCount()-1;
			bem.hydros.SetCount(num);
			//menuOpen.arrayModel.SetCount(num);
			//menuConvert.arrayModel.SetCount(num);
			Exclamation(Format(t_("Model '%s' does not match with the one previously loaded: %s"), DeQtf(file), strError));
			return;
		}
		
		int id = bem.hydros.GetCount()-1;
		HydroClass &data = bem.hydros[id];
		
		data.hd().Report();
		mainSummary.Report(data.hd(), id);
		if (data.hd().Nf < 0)
			return;
		
		mainArrange.Load(bem.hydros);
		
		menuOpen.arrayModel.Add(data.hd().GetCodeStr(), data.hd().name);
		menuOpen.butRemove.Enable();
		menuConvert.arrayModel.Add(data.hd().GetCodeStr(), data.hd().name);
		if (menuConvert.arrayModel.GetCursor() < 0)
			menuConvert.arrayModel.SetCursor(0);
		mainTab.GetItem(mainTab.Find(mainArrange)).Enable(true);	
		mainTab.GetItem(mainTab.Find(mainA)).Enable(mainA.Load(bem.hydros));	
		mainTab.GetItem(mainTab.Find(mainB)).Enable(mainB.Load(bem.hydros));
		mainTab.GetItem(mainTab.Find(mainForceSC)).Enable(mainForceSC.Load(bem.hydros));
		mainTab.GetItem(mainTab.Find(mainForceFK)).Enable(mainForceFK.Load(bem.hydros));
		mainTab.GetItem(mainTab.Find(mainForceEX)).Enable(mainForceEX.Load(bem.hydros));
		mainTab.GetItem(mainTab.Find(mainRAO)).Enable(mainRAO.Load(bem.hydros));
		mainTab.GetItem(mainTab.Find(mainView)).Enable(true);
		if (data.hd().IsLoadedStateSpace())
			mainTab.GetItem(mainTab.Find(mainStateSpace)).Enable(true);
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}
}
/*
void Main::WindowAdditionalData(BEMData &bem, HydroClass &data) {
	if (!IsNull(data.hd().g) && !IsNull(data.hd().len) && !IsNull(data.hd().rho) 
		&& !IsNull(data.hd().h))
		return;
		
	WithAdditionalDataLoad<TopWindow> dialog;
	CtrlLayoutOK(dialog, t_("Set pending data"));
	if (IsNull(data.hd().g)) 
		data.hd().g = bem.g;
	else
		dialog.g.Disable();
	
	if (IsNull(data.hd().len)) 
		data.hd().len = bem.length;
	else
		dialog.len.Disable();
	
	if (IsNull(data.hd().rho)) 
		data.hd().rho = bem.rho;
	else
		dialog.rho.Disable();
	
	if (IsNull(data.hd().h)) 
		data.hd().h = bem.depth;
	else
		dialog.h.Disable();
	if (data.hd().h == -1)
		dialog.h <<= "INFINITY";
	else
		dialog.h <<= FormatDouble(data.hd().h);
	
	CtrlRetriever rf;
	rf
		(dialog.g, data.hd().g)
		(dialog.len, data.hd().len)
		(dialog.rho, data.hd().rho)
		//(dialog.h, data.hd().h)
	;
	dialog.Execute();
	rf.Retrieve();		
	
	if (~dialog.h == "INFINITY")
		data.hd().h = -1;
	else
		data.hd().h = ScanDouble(dialog.h.GetData().ToString());
}
*/		
void Main::OnConvert() {
	try {
		int id = menuConvert.arrayModel.GetCursor();
		if (id < 0) {
			Exclamation(t_("Please select a model to export"));
			return;
		}
		Hydro::BEM_SOFT type;	
		switch (menuConvert.opt) {
		case 0:	type = Hydro::WAMIT_1_3;	break;
		case 1:	type = Hydro::FAST_WAMIT;	break;
		case 2:	type = Hydro::UNKNOWN;		break;
		default: throw Exc(t_("Unknown type in OnConvert()"));
		}
		bem.hydros[id].hd().SaveAs(~menuConvert.file, type);	
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
	}
}

void Main::OnView() {
	WaitCursor waitcursor;
	
	String file = ~menuView.file;
	for (int i = 0; i < surfs.GetCount(); ++i) {
		if (surfs[i].mh().file == file) {
			Exclamation(t_("Model already loaded"));
			return;
		}
	}
	String ext = ToLower(GetFileExt(file));
	if (ext == ".dat") {
		Nemoh &data = surfs.Create<Nemoh>(bem);
		if (!data.LoadDatMesh(file)) {
			surfs.SetCount(surfs.GetCount()-1);
			Wamit &data = surfs.Create<Wamit>(bem);
			if (!data.LoadDatMesh(file)) {
				Exclamation(DeQtfLf(Format(t_("Problem loading '%s'") + x_("\n%s"), file, data.mh().GetLastError())));	
				surfs.SetCount(surfs.GetCount()-1);
				return;
			}		
		} 
	} else if (ext == ".gdf") {
		Wamit &data = surfs.Create<Wamit>(bem);
		if (!data.LoadGdfMesh(file)) {
			Exclamation(DeQtfLf(Format(t_("Problem loading '%s'") + x_("\n%s"), file, data.mh().GetLastError())));	
			surfs.SetCount(surfs.GetCount()-1);
			return;
		}
	} else {
		Exclamation(DeQtfLf(Format(t_("Problem loading '%s'") + x_("\n%s"), file, t_("Unknown file format"))));	
		return;
	}
	mainView.CalcEnvelope();
	mainView.ZoomToFit();
	mainTab.Set(mainView);
}

Main::~Main() {
	if (!closed)
		Close();
}

void Main::Close(bool store) {
	if (store) {
		bem.StoreSerializeJson();
		StoreSerializeJson();
	}
	Thread::ShutdownThreads(); 
	RejectBreak(IDOK);		// Empty EditStrings does not disturb
	TopWindow::Close();
	closed = true;
}

void Main::Jsonize(JsonIO &json) {
	json
		("menuOpen_file", menuOpen.file)
		("menuConvert_file", menuConvert.file)
		("menuConvert_opt", menuConvert.opt)
		("menuPlot_autoFit", menuPlot.autoFit)
		("menuPlot_opwT", menuPlot.opwT)
		("menuPlot_showPoints", menuPlot.showPoints)
		("menuPlot_showPhase", menuPlot.showPhase)
		("menuPlot_showAdim", menuPlot.showAdim)
		("menuView_file", menuView.file)
	;
}

void MenuOptions::Init(BEMData &bem) {
	CtrlLayout(*this);
	
	this->bem = &bem;
	butSave <<= THISBACK(OnSave);
}

void MenuOptions::Load() {
	g <<= bem->g;
	rho <<= bem->rho;
	length <<= bem->length;
	depth <<= bem->depth;
	discardNegDOF <<= bem->discardNegDOF;
	thres <<= bem->thres;
	calcAwinf <<= bem->calcAwinf;
	maxTimeA <<= bem->maxTimeA;
	numValsA <<= bem->numValsA;	
}

void MenuOptions::OnSave() {
	bem->g = ~g;
	bem->rho = ~rho;
	bem->length = ~length;
	bem->depth = ~depth;
	bem->discardNegDOF = ~discardNegDOF;
	bem->thres = ~thres;
	bem->calcAwinf = ~calcAwinf;
	bem->maxTimeA = ~maxTimeA;
	bem->numValsA = ~numValsA;	
}

bool MenuOptions::IsChanged() {
	if ((bem->g != ~g) || (bem->rho != ~rho) || (bem->length != ~length) || (bem->depth != ~depth) || 
		(bem->discardNegDOF != ~discardNegDOF) || (bem->thres != ~thres) || 
		(bem->calcAwinf != ~calcAwinf) || (bem->maxTimeA != ~maxTimeA) || (bem->numValsA != ~numValsA))
		return true;
	
	return false;
}

void MenuAbout::Init() {
	CtrlLayout(*this);
	
	String qtf = GetTopic(x_("BEMRosetta/BEMRosetta/main/About$en-us")); 
	Hydro::SetBuildInfo(qtf);
	info.SetQTF(qtf);
}
	
void MainSummary::Init() {
	CtrlLayout(*this);
	array.SetLineCy(EditField::GetStdHeight()).MultiSelect();
	array.WhenBar = THISBACK(OnArrayBar);
}

void MainSummary::Clear() {
	array.Reset();
	array.SetLineCy(EditField::GetStdHeight()).MultiSelect();;
}

void MainSummary::Report(Hydro &data, int id) {
	if (array.GetColumnCount() == 0)
		array.AddColumn("Param");
	if (id >= array.GetColumnCount()-1)
		array.AddColumn(Format("#%d %s", id+1, data.name));
	int row = 0;
	int col = id + 1;
	
	array.Set(row, 0, t_("File"));				array.Set(row++, col, data.file);
	array.Set(row, 0, t_("Name"));				array.Set(row++, col, data.name);
	array.Set(row, 0, t_("Soft"));				array.Set(row++, col, data.GetCodeStr());
	String sg = IsNull(data.g) ? x_("unknown") : Format("%.3f", data.g);
	array.Set(row, 0, t_("g [m/s2]"));			array.Set(row++, col, sg);
	String srho = IsNull(data.rho) ? x_("unknown") : Format("%.3f", data.rho);
	array.Set(row, 0, t_("rho [kg/m3]"));		array.Set(row++, col, srho);
	array.Set(row, 0, t_("h (water depth) [m]"));array.Set(row++, col,data.h < 0 ? x_(t_("INFINITY")) : FormatDouble(data.h));
	String slen = IsNull(data.len) ? x_("unknown") : Format("%.1f", data.len);
	array.Set(row, 0, t_("length scale [m]"));	array.Set(row++, col, slen);
	
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
	array.Set(row, 0, t_("C available"));		array.Set(row++, col, data.IsLoadedC() 	   ? t_("Yes") : t_("No"));
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
						array.Set(row, 0, sib + " " + Format(t_("C(%d,%d) [%s]"), i+1, j+1, Hydro::C_units(i, j)));	array.Set(row++, col, FormatDouble(data.C_dim(ib, i, j), 3, FD_EXP));		
					}
				}
			}
		} else {
			for (int i = 0; i < 6; ++i) {
				for (int j = 0; j < 6; ++j) {
					if (!Hydro::C_units(i, j).IsEmpty()) {
						array.Set(row, 0, sib + " " + Format(t_("C(%d,%d) [%s]"), i+1, j+1, Hydro::C_units(i, j)));	array.Set(row++, col, "-");		
					}
				}
			}
		}
	}	
}

void MainSummary::OnArrayBar(Bar &menu) {
	menu.Add(t_("Select all"), Null, THISBACK(ArraySelect)).Key(K_CTRL_A)
								.Help(t_("Select all rows"));
								
	int count = array.GetSelectCount();
	if (count == 0)
		menu.Add(t_("No row selected"), Null, Null).Enable(false).Bold(true);
	else {
		menu.Add(Format(t_("Selected %d rows"), count), Null, Null).Enable(false).Bold(true);
		menu.Add(t_("Copy"), ScatterImgP::Copy(), THISBACK(ArrayCopy)).Key(K_CTRL_C)
									.Help(t_("Copy selected rows to clipboard"));
	}
}

void MainSummary::ArrayCopy() {
	array.SetClipboard(true, true);
}

void MainSummary::ArraySelect() {
	array.Select(0, array.GetCount(), true);
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
	tab.Reset();
	arrangeDOF.Clear();
	for (int i = 0; i < hydro.GetCount(); ++i) {
		ArrangeDOF &arr = arrangeDOF.Add();
		arr.Init(hydro[i].hd());
		tab.Add(arr.SizePos(), Format("[%s] %s", hydro[i].hd().GetCodeStr(), hydro[i].hd().name));
	}
}

void MainABForce::Init(DataToShow dataToShow) {
	CtrlLayout(*this);
	
	this->dataToShow = dataToShow;
	
	selTab = 0;
	isFilling = false;
	tab.WhenSet = [&] {
		if (!isFilling)
			selTab = tab.Get();
	};
}

void MainABForce::Clear() {
	tab.Reset();
	selTab = 0;
}

bool MainABForce::Load(Upp::Array<HydroClass> &hydro) {
	if (hydro.IsEmpty())
		return false;
	isFilling = true;
	tab.Reset();
	String format;
	switch (dataToShow) {
	case DATA_A:		format = t_("A%s%s");		break;		
	case DATA_B:		format = t_("B%s%s");		break;
	case DATA_FORCE_SC:	format = t_("Fsc%s%.1fº");	break;
	case DATA_FORCE_FK:	format = t_("Ffk%s%.1fº");	break;
	case DATA_FORCE_EX:	format = t_("Fex%s%.1fº");	break;
	case DATA_RAO:		format = t_("RAO%s%.1fº");	break;
	}
	int sdof = 6*hydro[0].hd().Nb;
	if (dataToShow == DATA_A || dataToShow == DATA_B) {
		plots.SetCount(sdof);
		for (int i = 0; i < sdof; ++i) {
			plots[i].SetCount(sdof);
			for (int j = 0; j < sdof; ++j) {
				plots[i][j].Init(i, j, Null, dataToShow);
				if (plots[i][j].Load(hydro))
					tab.Add(plots[i][j].SizePos(), Format(format, Hydro::StrDOFAbrev(i), Hydro::StrDOFAbrev(j)));
			}
		}
	} else {
		int Nh = hydro[0].hd().Nh;
		if (Nh < 0)
			return false;
		plots.SetCount(Nh);
		for (int ih = 0; ih < Nh; ++ih) 
			plots[ih].SetCount(sdof);
		for (int i = 0; i < sdof; ++i) {
			for (int ih = 0; ih < Nh; ++ih) {
				plots[ih][i].Init(i, ih, hydro[0].hd().head[ih], dataToShow);
				if (plots[ih][i].Load(hydro))
					tab.Add(plots[ih][i].SizePos(), Format(format, Hydro::StrDOFAbrev(i), hydro[0].hd().head[ih]));
			}
		}
	}
	
	isFilling = false;
	if (tab.GetCount() == 0)
		return false;
	else if (tab.GetCount() > selTab)	
		tab.Set(selTab);
	return true;
}

String FormatDOF(int i, int j) {
	if (i != j)
		return Format("%s_%s", Hydro::StrDOF(i), Hydro::StrDOF(j));
	else
		return Hydro::StrDOF(i);
}

void MainPlot::Init(int i, int j_h, double h, DataToShow dataToShow) {
	CtrlLayout(*this);
	
	this->i = i;
	this->j_h = j_h;
	this->dataToShow = dataToShow;
	scatter.ShowAllMenus();
	String title, labelY, labelY2;
	switch (dataToShow) {
	case DATA_A:		title = Format(t_("Added mass %s"), FormatDOF(i, j_h));		
						labelY = t_("Added mass");				
						break;		
	case DATA_B:		title = Format(t_("Radiation damping %s"), FormatDOF(i, j_h));
						labelY = t_("Radiation damping");		
						break;
	case DATA_FORCE_SC:	title = Format(t_("Diffraction scattering force %s heading %.1fº"), Hydro::StrDOF(i), h);
						labelY = t_("Diffraction scattering force");
						labelY2 = t_("Diffraction scattering force phase [º]");	
						break;
	case DATA_FORCE_FK:	title = Format(t_("Froude-Krylov force %s heading %.1fº"), Hydro::StrDOF(i), h);
						labelY = t_("Froude-Krylov force");		
						labelY2 = t_("Froude-Krylov force phase [º]");			
						break;
	case DATA_FORCE_EX:	title = Format(t_("Excitation Force %s heading %.1fº"), Hydro::StrDOF(i), h);
						labelY = t_("Excitation force");		
						labelY2 = t_("Excitation force phase [º]");				
						break;
	case DATA_RAO:		title = Format(t_("Response Amplitude Operator %s heading %.1fº"), Hydro::StrDOF(i), h);
						labelY = t_("RAO []");		
						labelY2 = t_("RAO phase [º]");							
						break;
	}
	scatter.SetTitle(title);
	scatter.SetLabelY(labelY);
	if (!labelY2.IsEmpty()) {
		scatter.SetDrawY2Reticle(true);
		scatter.SetPlotAreaRightMargin(80);
		scatter.SetLabelY2(labelY2);
	}
}

bool MainPlot::Load(Upp::Array<HydroClass> &hydro) {
	scatter.RemoveAllSeries();
	ABF_source.SetCount(hydro.GetCount());
	ABF_source2.SetCount(hydro.GetCount());
	Ainf_source.SetCount(hydro.GetCount());
	
	bool adim = ma().menuPlot.showAdim;
	int markW = ma().menuPlot.showPoints ? 10 : 0;
	bool show_w = ma().menuPlot.opwT == 0;
	if (show_w)
		scatter.SetLabelX(t_("w [rad/s]"));
	else
		scatter.SetLabelX(t_("T [s]"));
	
	bool loaded = false;
	for (int id = 0; id < hydro.GetCount(); ++id) {
		Hydro &hy = hydro[id].hd();
		String nameType = Format("%s(%s)", hy.name, hy.GetCodeStrAbr());
		if (dataToShow == DATA_A) {
			Upp::Color acolor = Null;
			if (hy.IsLoadedA()) {
				if (ABF_source[id].Init(hy, i, j_h, PLOT_A, show_w, adim)) {
					loaded = true;
					scatter.AddSeries(ABF_source[id]).Legend(Format(t_("A_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>();
					if (!adim)
						scatter.Units("Ns2/m");
					double dummy;
					scatter.GetStroke(scatter.GetCount()-1, dummy, acolor);
				}
			}
			if (hy.IsLoadedAwinf()) {
				if (Ainf_source[id].Init(hy, i, j_h, PLOT_AINF, show_w, adim)) {
					loaded = true;
					scatter.AddSeries(Ainf_source[id]).Legend(Format(t_("Ainf_%s"), nameType)).Dash(LINE_DOTTED).Stroke(2, acolor).NoMark();
					if (!adim)
						scatter.Units("Ns2/m");
				}
			}
		} else if (dataToShow == DATA_B && hy.IsLoadedB()) {
			if (ABF_source[id].Init(hy, i, j_h, PLOT_B, show_w, adim)) {
				loaded = true;
				scatter.AddSeries(ABF_source[id]).Legend(Format(t_("B_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>();
				if (!adim)
					scatter.Units("Ns/m");
			}
		} else if (dataToShow == DATA_FORCE_SC && hy.IsLoadedFsc()) {
			if (ABF_source[id].Init(hy, i, j_h, PLOT_FORCE_SC_MA, show_w, adim)) {
				loaded = true;
				scatter.AddSeries(ABF_source[id]).Legend(Format(t_("Fsc_ma_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>();
				if (!adim)
					scatter.Units("N");
				if (ABF_source2[id].Init(hy, i, j_h, PLOT_FORCE_SC_PH, show_w, adim)) {
					loaded = true;
					if (ma().menuPlot.showPhase)
						scatter.AddSeries(ABF_source2[id]).Legend(Format(t_("Fsc_ph_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetDataSecondaryY();
				}
			}
		} else if (dataToShow == DATA_FORCE_FK && hy.IsLoadedFfk()) {
			if (ABF_source[id].Init(hy, i, j_h, PLOT_FORCE_FK_MA, show_w, adim)) {
				loaded = true;
				scatter.AddSeries(ABF_source[id]).Legend(Format(t_("Ffk_ma_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>();
				if (!adim)
					scatter.Units("N");
				if (ABF_source2[id].Init(hy, i, j_h, PLOT_FORCE_FK_PH, show_w, adim)) {
					loaded = true;
					if (ma().menuPlot.showPhase)
						scatter.AddSeries(ABF_source2[id]).Legend(Format(t_("Ffk_ph_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetDataSecondaryY();
				}
			}
		} else if (dataToShow == DATA_FORCE_EX && hy.IsLoadedFex()) {
			if (ABF_source[id].Init(hy, i, j_h, PLOT_FORCE_EX_MA, show_w, adim)) {
				loaded = true;
				scatter.AddSeries(ABF_source[id]).Legend(Format(t_("Fex_ma_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>();
				if (!adim)
					scatter.Units("N");
				if (ABF_source2[id].Init(hy, i, j_h, PLOT_FORCE_EX_PH, show_w, adim)) {
					loaded = true;
					if (ma().menuPlot.showPhase)
						scatter.AddSeries(ABF_source2[id]).Legend(Format(t_("Fex_ph_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetDataSecondaryY();
				}
			}
		} else if (dataToShow == DATA_RAO && hy.IsLoadedRAO()) {
			if (ABF_source[id].Init(hy, i, j_h, PLOT_RAO_MA, show_w, adim)) {
				loaded = true;
				scatter.AddSeries(ABF_source[id]).Legend(Format(t_("RAO_ma_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>();
				if (ABF_source2[id].Init(hy, i, j_h, PLOT_RAO_PH, show_w, adim)) {
					loaded = true;
					if (ma().menuPlot.showPhase)
						scatter.AddSeries(ABF_source2[id]).Legend(Format(t_("RAO_ph_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetDataSecondaryY();
				}
			}
		}
	}
	if (ma().menuPlot.autoFit)
		scatter.ZoomToFit(true, true);
	return loaded;
}

void MainStateSpace::Init() {
	CtrlLayout(*this);
	
	scatterZ.ShowAllMenus();
	scatterZ.SetTitle(t_("Frequency response"));
	scatterZ.SetDrawY2Reticle(true);
}

bool MainStateSpace::Load(Upp::Array<HydroClass> &hydro) {
	scatterZ.RemoveAllSeries();
	Z_source.SetCount(hydro.GetCount());
	Z_source2.SetCount(hydro.GetCount());
	TFS_source.SetCount(hydro.GetCount());
	TFS_source2.SetCount(hydro.GetCount());
	
	bool adim = ma().menuPlot.showAdim;
	int markW = ma().menuPlot.showPoints ? 10 : 0;
	bool show_w = ma().menuPlot.opwT == 0;
	if (show_w) 
		scatterZ.SetLabelX(t_("w [rad/s]"));
	else 
		scatterZ.SetLabelX(t_("T [s]"));
	
	bool loaded = false;
	for (int id = 0; id < hydro.GetCount(); ++id) {	
		if (hydro[id].hd().IsLoadedStateSpace()) {
			if (Z_source[id].Init(hydro[id].hd(), 0, 0, PLOT_Z_MA, show_w, adim)) {
				loaded = true;
				scatterZ.AddSeries(Z_source[id]).Legend(Format(t_("Z Magnitude %s"), hydro[id].hd().name)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().Units("dB");
				if (Z_source2[id].Init(hydro[id].hd(), 0, 0, PLOT_Z_PH, show_w, adim)) {
					loaded = true;
					if (ma().menuPlot.showPhase)
						scatterZ.AddSeries(Z_source2[id]).Legend(Format(t_("Z Phase %s"), hydro[id].hd().name)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().Units("rad").SetDataSecondaryY();
				}
			}
			if (TFS_source[id].Init(hydro[id].hd(), 0, 0, PLOT_TFS_MA, show_w, adim)) {
				loaded = true;
				scatterZ.AddSeries(TFS_source[id]).Legend(Format(t_("TFSResponse Magnitude %s"), hydro[id].hd().name)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().Units("dB");
				if (TFS_source2[id].Init(hydro[id].hd(), 0, 0, PLOT_TFS_PH, show_w, adim)) {
					loaded = true;
					if (ma().menuPlot.showPhase)
						scatterZ.AddSeries(TFS_source2[id]).Legend(Format(t_("TFSResponse Phase %s"), hydro[id].hd().name)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().Units("rad").SetDataSecondaryY();
				}
			}
		}
	}
	if (ma().menuPlot.autoFit) 
		scatterZ.ZoomToFit(true, true);
	return loaded;
}

void MainView::Init() {
	CtrlLayout(*this);
	
	gl.WhenPaint = THISBACK(OnPaint);	
}


static Color GetColor(int id) {
	static Color col[10] = {Color(93,185,86), Color(180,96,189), Color(176,177,67), Color(106,126,205),
							Color(219,142,68), Color(79,186,172), Color(203,83,67), Color(83,133,68),
							Color(199,89,128), Color(147,119,58)};
	int mod = id % 10;
	if (id < 10)
		return col[id];
	else
		return col[mod];
}

void MainView::OnPaint() {
	double mx = max(max(maxX, maxY), maxZ)/4.;
	gl.PaintAxis(mx, mx, mx);	
	
	for (int i = 0; i < ma().surfs.GetCount(); ++i)
		gl.PaintSurface(ma().surfs[i].mh(), GetColor(i));
}

void MainView::CalcEnvelope() {
	env.Reset();
	for (int i = 0; i < ma().surfs.GetCount(); ++i)
		env.MixEnvelope(ma().surfs[i].mh().env);
	
	maxX = env.maxX;
	if (env.minX < 0)
		maxX = max(maxX, -env.minX);
	
	maxY = env.maxY;
	if (env.minY < 0)
		maxY = max(maxY, -env.minY);
		
	maxZ = env.maxZ;
	if (env.minZ < 0)
		maxZ = max(maxZ, -env.minZ);
}

void MainView::ZoomToFit() {
	//double mx = max(max(env.maxX, env.maxY), env.maxZ);
	
	//gl.SetZoomFactor(mx*1.25);
}
			
Main &ma(Main *m) {
	static Main *mp = 0;
	if (m)
		mp = m;
	return *mp;
}


GUI_APP_MAIN {
	ConsoleOutput console;
	
	Ctrl::SetAppName(t_("Hydrodynamic coefficents viewer and converter"));
	Ctrl::GlobalBackPaint();
	
	String errorStr;
	try {
		Main main;
		
		main.Init();
		main.Run();
		
		main.Close(true);
	} catch (Exc e) {
		errorStr = e;
	} catch(const char *cad) {
		errorStr = cad;
	} catch(const std::string &e) {
		errorStr = e.c_str();	
	} catch (const std::exception &e) {
		errorStr = e.what();
	} catch(...) {
		errorStr = t_("Unknown error");
	}	
	if (!errorStr.IsEmpty())
		Exclamation(t_("Internal error:") + x_("&") + DeQtf(errorStr) + x_("&") + t_("Program ended"));
}

