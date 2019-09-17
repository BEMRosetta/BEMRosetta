#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas.h>

#include <BEMRosetta/BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;

#include "main.h"


void MainBEM::Init() {
	CtrlLayout(*this);
	
	mbm(this);
	
	CtrlLayout(menuOpen);
	menuOpen.file.WhenChange = THISBACK(OnLoad);
	menuOpen.file.BrowseRightWidth(40).UseOpenFolder().BrowseOpenFolderWidth(10);
	menuOpen.butLoad.WhenAction = [&] {OnLoad();};
	
	menuOpen.arrayModel.NoHeader().NoVertGrid().AutoHideSb();
	menuOpen.arrayModel.AddColumn("", 20);	
	menuOpen.arrayModel.AddColumn("", 20);
	
	menuOpen.butRemove.Disable();	
	menuOpen.butRemove.WhenAction = [&] {
		ma().bem.hydros.Clear();
		
		mainSummary.Clear();
		menuOpen.arrayModel.Clear();
		menuOpen.butRemove.Disable();
		menuConvert.arrayModel.Clear();	
		
		mainArrange.Clear();	mainTab.GetItem(mainTab.Find(mainArrange)).Disable();
		mainStiffness.Clear();	mainTab.GetItem(mainTab.Find(mainStiffness)).Disable();
		mainA.Clear();			mainTab.GetItem(mainTab.Find(mainA)).Disable();
		mainB.Clear();			mainTab.GetItem(mainTab.Find(mainB)).Disable();
		mainForceSC.Clear();	mainTab.GetItem(mainTab.Find(mainForceSC)).Disable();
		mainForceFK.Clear();	mainTab.GetItem(mainTab.Find(mainForceFK)).Disable();
		mainForceEX.Clear();	mainTab.GetItem(mainTab.Find(mainForceEX)).Disable();
		mainRAO.Clear();		mainTab.GetItem(mainTab.Find(mainRAO)).Disable();
	};
	
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
	menuPlot.butZoomToFit.WhenAction = [&] {GetSelPlot().scatter.ZoomToFit(true, true);};
	menuPlot.autoFit.WhenAction 	 = [&] {LoadSelTab(ma().bem);};
	menuPlot.opwT.WhenAction 	 	 = [&] {LoadSelTab(ma().bem);};
	menuPlot.showPoints.WhenAction 	 = [&] {LoadSelTab(ma().bem);};
	menuPlot.showPhase.WhenAction 	 = [&] {LoadSelTab(ma().bem);};
	menuPlot.showNdim.WhenAction 	 = [&] {LoadSelTab(ma().bem);};
	
	OnOpt();
	
	CtrlLayout(menuFOAMM);
	menuFOAMM.file.WhenChange = THISBACK(OnFOAMM);
	menuFOAMM.file.BrowseRightWidth(40).UseOpenFolder(true).BrowseOpenFolderWidth(10);
	cancelFOAMM = false;
	menuFOAMM.progress.Show(false);
	menuFOAMM.butCancel.Show(false);
	menuFOAMM.butLoad.WhenAction 	= [&] {OnFOAMM();};
	menuFOAMM.butCancel.WhenAction 	= [&] {cancelFOAMM = true;};
	
	menuFOAMM.arrayModel.NoHeader().NoVertGrid().AutoHideSb();
	menuFOAMM.arrayModel.AddColumn("", 20);	
	menuFOAMM.arrayModel.AddColumn("", 20);	
	
	OnOpt();
		
	menuTab.Add(menuOpen.SizePos(), 	t_("Open"));
	menuTab.Add(menuConvert.SizePos(), 	t_("Convert"));
	menuTab.Add(menuPlot.SizePos(), 	t_("Plot")).Disable();
	if (ma().bem.experimentalFOAMM) 
		menuTab.Add(menuFOAMM.SizePos(), t_("State Space"));
	
	mainTab.WhenSet = [&] {
		bool plot = true;
		if (ma().bem.hydros.IsEmpty())
			plot = false;
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
		} else {
			plot = false;
			menuPlot.showPhase.Enable(false);
		}
		TabCtrl::Item& plotIt = menuTab.GetItem(menuTab.Find(menuPlot));
		plotIt.Enable(plot);
		if (plot) {
			plotIt.Text(t_("Plot"));
			menuTab.Set(menuPlot);
		} else {
			plotIt.Text("");
			menuTab.Set(menuOpen);
		}
	};
	
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
	if (mainTab.Get() == mainTab.Find(mainStateSpace))
		mainStateSpace.Load(bem);
	else if (mainTab.Get() == mainTab.Find(mainStiffness))
		mainStiffness.Load(bem.hydros);
	else 
		GetSelTab().Load(bem);
}

MainABForce &MainBEM::GetSelTab() {
	return *(static_cast<MainABForce*>(mainTab.GetItem(mainTab.Get()).GetSlave()));
}

MainPlot &MainBEM::GetSelPlot() {
	TabCtrl &tab = GetSelTab().tab;
	return *(static_cast<MainPlot*>(tab.GetItem(tab.Get()).GetSlave()));
}

void MainBEM::OnOpt() {
	menuOpen.file.ClearTypes();
	const String bemFiles = ".1 .3 .hst .4 .out .dat .cal .inf .ah1 .lis .mat";
	String bemFilesAst = clone(bemFiles);
	bemFilesAst.Replace(".", "*.");
	menuOpen.file.Type(Format("All supported BEM files (%s)", bemFiles), bemFilesAst);
	menuOpen.file.AllFilesType();
	String extOpen = ToLower(GetFileExt(menuOpen.file.GetData().ToString()));
	if (extOpen.IsEmpty())
		menuOpen.file.ActiveType(0);
	else if (bemFiles.Find(extOpen) >= 0)
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
	
	menuFOAMM.file.ClearTypes();
	menuFOAMM.file.Type(t_("Maynooth COER FOAMM file *.mat"), "*.mat");
	menuFOAMM.file.AllFilesType();
	String extFOAMM = ToLower(GetFileExt(menuFOAMM.file.GetData().ToString()));
	if (extFOAMM.IsEmpty())
		menuFOAMM.file.ActiveType(0);
	else if (extFOAMM == ".mat")
		menuFOAMM.file.ActiveType(0);
	else
		menuFOAMM.file.ActiveType(1);	
}

bool MainBEM::OnLoad() {
	String file = ~menuOpen.file;
	
	try {
		WaitCursor wait;
		Progress progress(t_("Loading BEM files..."), 100); 
		
		ma().bem.Load(file, [&](String str, int _pos) {progress.SetText(str); progress.SetPos(_pos);});
		
		int id = ma().bem.hydros.GetCount()-1;
		HydroClass &data = ma().bem.hydros[id];
		
		data.hd().Report();
		mainSummary.Report(data.hd(), id);
		if (data.hd().Nf < 0)
			return false;
		
		mainArrange.Load(ma().bem.hydros);
		
		menuOpen.arrayModel.Add(data.hd().GetCodeStr(), data.hd().name);
		menuOpen.butRemove.Enable();
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
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
		return false;
	}
	return true;
}
	
bool MainBEM::OnConvert() {
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
		ma().bem.hydros[id].hd().SaveAs(~menuConvert.file, type);	
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
		return false;
	}
	return true;
}

bool MainBEM::OnFOAMM() {
	try {
		int id = menuFOAMM.arrayModel.GetCursor();
		if (id < 0) {
			Exclamation(t_("Please select a model to get State Space"));
			return false;
		}
		menuFOAMM.progress.Show(true);
		menuFOAMM.butCancel.Show(true);
		ma().bem.hydros[id].hd().GetFOAMM(~menuFOAMM.file, [&](String str) {
				menuFOAMM.progress++; 
				if (!str.IsEmpty()) {
					//menuFOAMM.progress.Show(false);
					//menuFOAMM.butCancel.Show(false);
					str.Replace("\r", "");
					str.Replace("\n\n", "\n");
					Exclamation (DeQtfLf(str));
					//return true;
				}
				ProcessEvents(); 
				return cancelFOAMM;
			});
		menuFOAMM.progress.Show(false);
		menuFOAMM.butCancel.Show(false);
		cancelFOAMM = false;	
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
		menuFOAMM.progress.Show(false);
		menuFOAMM.butCancel.Show(false);
		cancelFOAMM = false;
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
	;
}

void MainBEM::DragAndDrop(Point p, PasteClip& d) {
	if (IsDragAndDropSource())
		return;
	if (AcceptFiles(d)) {
		Vector<String> files = GetFiles(d);
		for (int i = 0; i < files.GetCount(); ++i) {
			menuOpen.file <<= files[i];
			OnLoad();
		}
	}
}

bool MainBEM::Key(dword key, int count) {
	if (key == K_CTRL_V) {
		Vector<String> files = GetFiles(Ctrl::Clipboard());
		for (int i = 0; i < files.GetCount(); ++i) {
			menuOpen.file <<= files[i];
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
	tab.Reset();
	arrangeDOF.Clear();
	for (int i = 0; i < hydro.GetCount(); ++i) {
		ArrangeDOF &arr = arrangeDOF.Add();
		arr.Init(hydro[i].hd());
		tab.Add(arr.SizePos(), Format("[%s] %s", hydro[i].hd().GetCodeStr(), hydro[i].hd().name));
	}
}

void MainABForce::Init(DataToShow _dataToShow) {
	CtrlLayout(*this);
	
	dataToShow = _dataToShow;
	
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

bool MainABForce::Load(BEMData &bem) {
	try {
		Upp::Array<HydroClass> &hydro = bem.hydros; 
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
		int sdof = 6*bem.Nb;
		if (dataToShow == DATA_A || dataToShow == DATA_B) {
			plots.SetCount(sdof);
			for (int i = 0; i < sdof; ++i) {
				plots[i].SetCount(sdof);
				for (int j = 0; j < sdof; ++j) {
					if (!bem.onlyDiagonal || i == j) {
						plots[i][j].Init(i, j, dataToShow);
						if (plots[i][j].Load(hydro)) {
							if (i != j)
								tab.Add(plots[i][j].SizePos(), Format(format, Hydro::StrDOFAbrev(i), Hydro::StrDOFAbrev(j)));
							else
								tab.Add(plots[i][j].SizePos(), Format(format, Hydro::StrDOF(i), ""));
						}
					}
				}
			}
		} else {
			int Nh = bem.head.GetCount();
			if (Nh < 0)
				return false;
			plots.SetCount(Nh);
			for (int ih = 0; ih < Nh; ++ih) 
				plots[ih].SetCount(sdof);
			for (int i = 0; i < sdof; ++i) {
				for (int ih = 0; ih < Nh; ++ih) {
					plots[ih][i].Init(i, bem.head[ih], dataToShow);
					if (plots[ih][i].Load(hydro))
						tab.Add(plots[ih][i].SizePos(), Format(format, Hydro::StrDOFAbrev(i), bem.head[ih]));
				}
			}
		}
		
		isFilling = false;
		if (tab.GetCount() == 0)
			return false;
		else if (tab.GetCount() > selTab)	
			tab.Set(selTab);
		return true;
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
		return false;
	}
}

String FormatDOF(int i, int j) {
	if (i != j)
		return Format("%s_%s", Hydro::StrDOF(i), Hydro::StrDOF(j));
	else
		return Hydro::StrDOF(i);
}

void MainPlot::Init(int _idof, double jdof_ih, DataToShow _dataToShow) {
	CtrlLayout(*this);
	
	idof = _idof;
	jdof = int(jdof_ih);
	heading = jdof_ih;
	dataToShow = _dataToShow;
	scatter.ShowAllMenus();
	String title, labelY, labelY2;
	switch (dataToShow) {
	case DATA_A:		title = Format(t_("Added mass %s"), FormatDOF(idof, jdof));		
						labelY = t_("Added mass");				
						break;		
	case DATA_B:		title = Format(t_("Radiation damping %s"), FormatDOF(idof, jdof));
						labelY = t_("Radiation damping");		
						break;
	case DATA_FORCE_SC:	title = Format(t_("Diffraction scattering force %s heading %.1fº"), Hydro::StrDOF(idof), heading);
						labelY = t_("Diffraction scattering force");
						labelY2 = t_("Diffraction scattering force phase [º]");	
						break;
	case DATA_FORCE_FK:	title = Format(t_("Froude-Krylov force %s heading %.1fº"), Hydro::StrDOF(idof), heading);
						labelY = t_("Froude-Krylov force");		
						labelY2 = t_("Froude-Krylov force phase [º]");			
						break;
	case DATA_FORCE_EX:	title = Format(t_("Excitation Force %s heading %.1fº"), Hydro::StrDOF(idof), heading);
						labelY = t_("Excitation force");		
						labelY2 = t_("Excitation force phase [º]");				
						break;
	case DATA_RAO:		title = Format(t_("Response Amplitude Operator %s heading %.1fº"), Hydro::StrDOF(idof), heading);
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
	
	dim = !mbm().menuPlot.showNdim;
	markW = mbm().menuPlot.showPoints ? 10 : 0;
	show_w = mbm().menuPlot.opwT == 0;
	if (show_w)
		scatter.SetLabelX(t_("w [rad/s]"));
	else
		scatter.SetLabelX(t_("T [s]"));
	
	bool loaded = false;
	for (int id = 0; id < hydro.GetCount(); ++id) {
		Hydro &hy = hydro[id].hd();
		int ih = -1;
		if (dataToShow != DATA_A && dataToShow != DATA_B) 
			ih = hy.GetHeadId(heading);
		String nameType = Format("%s(%s)", hy.name, hy.GetCodeStrAbr());
		if (dataToShow == DATA_A) {
			Upp::Color acolor = Null;
			if (hy.IsLoadedA()) {
				if (ABF_source[id].Init(hy, idof, jdof, PLOT_A, show_w, !dim)) {
					loaded = true;
					scatter.AddSeries(ABF_source[id]).Legend(Format(t_("A_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>();
					if (dim)
						scatter.Units("Ns2/m");
					double dummy;
					scatter.GetStroke(scatter.GetCount()-1, dummy, acolor);
				}
			}
			if (hy.IsLoadedAwinf()) {
				if (Ainf_source[id].Init(hy, idof, jdof, PLOT_AINF, show_w, !dim)) {
					loaded = true;
					scatter.AddSeries(Ainf_source[id]).Legend(Format(t_("Ainf_%s"), nameType)).Dash(LINE_DOTTED).Stroke(2, acolor).NoMark();
					if (dim)
						scatter.Units("Ns2/m");
				}
			}
		} else if (dataToShow == DATA_B && hy.IsLoadedB()) {
			if (ABF_source[id].Init(hy, idof, jdof, PLOT_B, show_w, !dim)) {
				loaded = true;
				scatter.AddSeries(ABF_source[id]).Legend(Format(t_("B_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>();
				if (dim)
					scatter.Units("Ns/m");
			}
		} else if (dataToShow == DATA_FORCE_SC && hy.IsLoadedFsc() && ih >= 0) {
			if (ABF_source[id].Init(hy, idof, ih, PLOT_FORCE_SC_MA, show_w, !dim)) {
				loaded = true;
				scatter.AddSeries(ABF_source[id]).Legend(Format(t_("Fsc_ma_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>();
				if (dim)
					scatter.Units("N");
				if (ABF_source2[id].Init(hy, idof, ih, PLOT_FORCE_SC_PH, show_w, !dim)) {
					loaded = true;
					if (mbm().menuPlot.showPhase)
						scatter.AddSeries(ABF_source2[id]).Legend(Format(t_("Fsc_ph_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetDataSecondaryY();
				}
			}
		} else if (dataToShow == DATA_FORCE_FK && hy.IsLoadedFfk() && ih >= 0) {
			if (ABF_source[id].Init(hy, idof, ih, PLOT_FORCE_FK_MA, show_w, !dim)) {
				loaded = true;
				scatter.AddSeries(ABF_source[id]).Legend(Format(t_("Ffk_ma_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>();
				if (dim)
					scatter.Units("N");
				if (ABF_source2[id].Init(hy, idof, ih, PLOT_FORCE_FK_PH, show_w, !dim)) {
					loaded = true;
					if (mbm().menuPlot.showPhase)
						scatter.AddSeries(ABF_source2[id]).Legend(Format(t_("Ffk_ph_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetDataSecondaryY();
				}
			}
		} else if (dataToShow == DATA_FORCE_EX && hy.IsLoadedFex() && ih >= 0) {
			if (ABF_source[id].Init(hy, idof, ih, PLOT_FORCE_EX_MA, show_w, !dim)) {
				loaded = true;
				scatter.AddSeries(ABF_source[id]).Legend(Format(t_("Fex_ma_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>();
				if (dim)
					scatter.Units("N");
				if (ABF_source2[id].Init(hy, idof, ih, PLOT_FORCE_EX_PH, show_w, !dim)) {
					loaded = true;
					if (mbm().menuPlot.showPhase)
						scatter.AddSeries(ABF_source2[id]).Legend(Format(t_("Fex_ph_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetDataSecondaryY();
				}
			}
		} else if (dataToShow == DATA_RAO && hy.IsLoadedRAO() && ih >= 0) {
			if (ABF_source[id].Init(hy, idof, ih, PLOT_RAO_MA, show_w, !dim)) {
				loaded = true;
				scatter.AddSeries(ABF_source[id]).Legend(Format(t_("RAO_ma_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>();
				if (ABF_source2[id].Init(hy, idof, ih, PLOT_RAO_PH, show_w, !dim)) {
					loaded = true;
					if (mbm().menuPlot.showPhase)
						scatter.AddSeries(ABF_source2[id]).Legend(Format(t_("RAO_ph_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetDataSecondaryY();
				}
			}
		}
	}
	if (mbm().menuPlot.autoFit)
		scatter.ZoomToFit(true, true);
	return loaded;
}

void MainStateSpace::Init(ArrayCtrl &array) {
	array.Reset();
	array.NoHeader().SetLineCy(EditField::GetStdHeight()).HeaderObject().Absolute();
	array.MultiSelect().SpanWideCells();
	array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, array);};
}

void MainStateSpace::Init() {
	scatter.ShowAllMenus();
	scatter.SetTitle(t_("Frequency response")).SetTitleFont(SansSerif(12));
	scatter.SetPlotAreaLeftMargin(70);
	
	splitter.Horz(tab.SizePos(), scatter.SizePos());
	Add(splitter.SizePos());
}

bool MainStateSpace::Load(BEMData &bem) {
	Upp::Array<HydroClass> &hydros = bem.hydros;
	int hnum = hydros.GetCount();
	
	scatter.RemoveAllSeries();
	Z_source.SetCount(hnum);
	Z_source2.SetCount(hnum);
	TFS_source.SetCount(hnum);
	TFS_source2.SetCount(hnum);
	
	bool dim = !mbm().menuPlot.showNdim;
	int markW = mbm().menuPlot.showPoints ? 10 : 0;
	bool show_w = mbm().menuPlot.opwT == 0;
	if (show_w) 
		scatter.SetLabelX(t_("w [rad/s]"));
	else 
		scatter.SetLabelX(t_("T [s]"));
	
	scatter.SetDrawY2Reticle(mbm().menuPlot.showPhase);
	scatter.SetPlotAreaRightMargin(mbm().menuPlot.showPhase ? 50 : 20);
	
	bool loaded = false;
	for (int id = 0; id < hydros.GetCount(); ++id) {
		Hydro &hydro = hydros[id].hd();	
		if (hydro.IsLoadedStateSpace()) {
			if (Z_source[id].Init(hydro, 0, 0, PLOT_Z_MA, show_w, !dim)) {
				loaded = true;
				scatter.AddSeries(Z_source[id]).Legend(Format(t_("Z Magnitude %s"), hydro.name)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().Units("dB");
				if (Z_source2[id].Init(hydro, 0, 0, PLOT_Z_PH, show_w, !dim)) {
					loaded = true;
					if (mbm().menuPlot.showPhase)
						scatter.AddSeries(Z_source2[id]).Legend(Format(t_("Z Phase %s"), hydro.name)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetDataSecondaryY().Units("rad");
				}
			}
			if (TFS_source[id].Init(hydro, 0, 0, PLOT_TFS_MA, show_w, !dim)) {
				loaded = true;
				scatter.AddSeries(TFS_source[id]).Legend(Format(t_("TFSResponse Magnitude %s"), hydro.name)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().Units("dB");
				if (TFS_source2[id].Init(hydro, 0, 0, PLOT_TFS_PH, show_w, !dim)) {
					loaded = true;
					if (mbm().menuPlot.showPhase)
						scatter.AddSeries(TFS_source2[id]).Legend(Format(t_("TFSResponse Phase %s"), hydro.name)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetDataSecondaryY().Units("rad");
				}
			}
		}
	}
	if (mbm().menuPlot.autoFit) 
		scatter.ZoomToFit(true, true);
	
	tab.Reset();
	arrays.Clear();
	
	for (int id = 0; id < hydros.GetCount(); ++id) {
		Hydro &hydro = hydros[id].hd();
		int row = 0;
		if (hydro.A_ss.size() > 0 || hydro.B_ss.size() > 0 || hydro.C_ss.size() > 0) {
			loaded = true;
			ArrayCtrl &array = arrays.Add();
			Init(array);
			tab.Add(array.SizePos(), hydro.name);
			if (hydro.A_ss.size() > 0) {
				if (hydro.A_ss.cols() > array.GetColumnCount()) {
					int ncols = static_cast<int>(hydro.A_ss.cols()) - array.GetColumnCount();
					for (int i = 0; i < ncols; ++i)
						array.AddColumn("", 80);
				}
				array.Set(row++, 0, AttrText(t_("A_ss")).Bold());
				for (int r = 0; r < hydro.A_ss.rows(); ++r)	{		
					for (int c = 0; c < hydro.A_ss.cols(); ++c)
						array.Set(row + r, c, hydro.A_ss(r, c));
				}
				row += static_cast<int>(hydro.A_ss.rows());
			}
			if (hydro.B_ss.size() > 0) {
				array.Set(row++, 0, AttrText(t_("B_ss")).Bold());
				for (int r = 0; r < hydro.B_ss.size(); ++r)		
					array.Set(row, r, hydro.B_ss(r));
				row++;
			}
			if (hydro.C_ss.size() > 0) {
				array.Set(row++, 0, AttrText(t_("C_ss")).Bold());
				for (int c = 0; c < hydro.C_ss.size(); ++c)			
					array.Set(row, c, hydro.C_ss(c));
				row++;
			}
			if (hydro.ssFrequencies.size() > 0) {
				array.Set(row++, 0, AttrText(t_("Frequencies")).Bold());
				for (int c = 0; c < hydro.ssFrequencies.size(); ++c)			
					array.Set(row, c, hydro.ssFrequencies[c]);
				row++;
			}
			if (hydro.ssFreqRange.size() > 0) {
				array.Set(row++, 0, AttrText(t_("FreqRange")).Bold());
				for (int c = 0; c < hydro.ssFreqRange.size(); ++c)			
					array.Set(row, c, hydro.ssFreqRange[c]);
				row++;
			}
			if (hydro.ssFrequencies_index.size() > 0) {
				array.Set(row++, 0, AttrText(t_("Frequencies_index")).Bold());
				for (int c = 0; c < hydro.ssFrequencies_index.size(); ++c)			
					array.Set(row, c, hydro.ssFrequencies_index[c]);
				row++;
			}						
			if (!IsNull(hydro.ssMAE)) {
				array.Set(row++, 0, AttrText(t_("MAE")).Bold());
				array.Set(row++, 0, hydro.ssMAE);
			}
		}
	}
	return loaded;
}

MainBEM &mbm(MainBEM *m) {
	static MainBEM *mp = 0;
	if (m)
		mp = m;
	return *mp;
}