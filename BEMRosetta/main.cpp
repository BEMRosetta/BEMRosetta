#include "main.h"

#define IMAGECLASS Img
#define IMAGEFILE <BEMRosetta/BEMRosetta/main.iml>
#include <Draw/iml_source.h>

#define TOPICFILE <BEMRosetta/BEMRosetta/main.tpp/all.i>
#include <Core/topic_group.h>

void Main::Init() {
	CtrlLayout(*this, "BEMRosetta");
	Sizeable().Zoomable().SetMinSize(Size(800, 600));
	Icon(Img::Rosetta64());
	LargeIcon(Img::Rosetta256());
	ma(this);
	
	int ret = LoadSerializeJson();
	if (ret == 0)
		throw Exc(t_("Configuration file corrupted"));
	
	CtrlLayout(menuOpen);
	menuOpen.file.WhenChange = THISBACK(OnLoad);
	menuOpen.file.BrowseRightWidth(40);
	menuOpen.butLoad <<= THISBACK(OnLoad);
	if (IsNull(menuOpen.optLoadIn))
		menuOpen.optLoadIn = 1;
	if (IsNull(menuOpen.opt))
		menuOpen.opt = 4;
	menuOpen.opt <<= THISBACK(OnOpt);
	menuOpen.butRemove.WhenAction = [&] {
		hydro.Clear();
		mainSummary.Clear();
		menuConvert.arrayModel.Clear();	
		mainArrange.Clear();	mainTab.GetItem(mainTab.Find(mainArrange)).Disable();
		mainA.Clear();			mainTab.GetItem(mainTab.Find(mainA)).Disable();
		mainB.Clear();			mainTab.GetItem(mainTab.Find(mainB)).Disable();
		mainForceSC.Clear();	mainTab.GetItem(mainTab.Find(mainForceSC)).Disable();
		mainForceFK.Clear();	mainTab.GetItem(mainTab.Find(mainForceFK)).Disable();
		mainForceEX.Clear();	mainTab.GetItem(mainTab.Find(mainForceEX)).Disable();
	};
	
	CtrlLayout(menuConvert);
	menuConvert.file.WhenChange = THISBACK(OnConvert);
	menuConvert.file.BrowseRightWidth(40);
	menuConvert.butLoad <<= THISBACK(OnConvert);
	if (IsNull(menuConvert.opt))
		menuConvert.opt = 0;
	menuConvert.opt <<= THISBACK(OnOpt);
	menuConvert.arrayModel.SetLineCy(EditField::GetStdHeight());
	menuConvert.arrayModel.NoHeader().NoVertGrid().AutoHideSb();
	menuConvert.arrayModel.AddColumn("", 20);	
	menuConvert.arrayModel.AddColumn("", 20);
	
	OnOpt();
	
	CtrlLayout(menuPlot);
	menuPlot.butZoomToFit.WhenAction = [&] {GetSelPlot().scatter.ZoomToFit(true, true);};
	menuPlot.autoFit.WhenAction 	 = [&] {GetSelTab().Load(hydro);};
	menuPlot.showPoints.WhenAction 	 = [&] {GetSelTab().Load(hydro);};
	
	menuTab.Add(menuOpen.SizePos(), "Open");
	menuTab.Add(menuConvert.SizePos(), "Convert");
	menuTab.Add(menuPlot.SizePos(), "Plot").Disable();
	
	menuAbout.Init();
	menuTab.Add(menuAbout.SizePos(), "About");
	
	menuTab.WhenSet = [&] {
		if (menuTab.IsAt(menuAbout)) 
			mainTab.Hide();
		else 
			mainTab.Show();
	};	
	
	if (ret == 1)
		menuTab.Set(menuAbout);
	
	mainTab.WhenSet = [&] {
		bool plot = true;
		if (mainTab.IsAt(mainA)) 
			mainA.Load(hydro);
		else if (mainTab.IsAt(mainB)) 
			mainB.Load(hydro);
		else if (mainTab.IsAt(mainForceSC)) 
			mainForceSC.Load(hydro);
		else if (mainTab.IsAt(mainForceFK)) 
			mainForceFK.Load(hydro);
		else if (mainTab.IsAt(mainForceEX)) 
			mainForceEX.Load(hydro);
		else 
			plot = false;
		TabCtrl::Item& plotIt = menuTab.GetItem(menuTab.Find(menuPlot));
		plotIt.Enable(plot);
		if (plot) {
			plotIt.Text("Plot");
			menuTab.Set(menuPlot);
		} else
			plotIt.Text("");
	};
	
	mainSummary.Init();
	mainTab.Add(mainSummary.SizePos(), "Summary");
	
	mainArrange.Init();
	mainTab.Add(mainArrange.SizePos(), "Arrange DOF").Disable();
	
	mainA.Init(DATA_A);
	mainTab.Add(mainA.SizePos(), "A").Disable();
	
	mainB.Init(DATA_B);
	mainTab.Add(mainB.SizePos(), "B").Disable();

	mainForceEX.Init(DATA_FORCE_EX);
	mainTab.Add(mainForceEX.SizePos(), "Fex").Disable();
	
	mainForceSC.Init(DATA_FORCE_SC);
	mainTab.Add(mainForceSC.SizePos(), "Fsc").Disable();
	
	mainForceFK.Init(DATA_FORCE_FK);
	mainTab.Add(mainForceFK.SizePos(), "Ffk").Disable();
		
	mainOutput.Init();
	mainTab.Add(mainOutput.SizePos(), "Output");
	
	Hydro::Print = [this](String s) {mainOutput.Print(s);};
	Hydro::PrintError = [this](String s) {mainOutput.Print(s); mainTab.Set(mainOutput);};
}

MainABForce &Main::GetSelTab() {
	return *(static_cast<MainABForce*>(mainTab.GetItem(mainTab.Get()).GetSlave()));
}

MainPlot &Main::GetSelPlot() {
	TabCtrl &tab = GetSelTab().tab;
	return *(static_cast<MainPlot*>(tab.GetItem(tab.Get()).GetSlave()));
}

void Main::OnOpt() {
	menuOpen.file.ClearTypes();
	switch (menuOpen.opt) {
	case 0:	menuOpen.file.Type("Wamit .1.3.hst file", "*.1 *.3 *.hst");		menuOpen.file = ForceExt(~menuOpen.file, ".1");		break;
	case 1:	menuOpen.file.Type("Wamit .out file", "*.out");					menuOpen.file = ForceExt(~menuOpen.file, ".out"); 	break;
	case 2:	menuOpen.file.Type("FAST HydroDyn file", "*.dat");				menuOpen.file = ForceExt(~menuOpen.file, ".dat"); 	break;
	case 3:	menuOpen.file.Type("Nemoh .cal file", "*.cal");					menuOpen.file = ForceExt(~menuOpen.file, ".cal"); 	break;
	case 4:	menuOpen.file.Type("SeaFEM .flavia.inf file", "*.inf");			menuOpen.file = ForceExt(~menuOpen.file, ".inf"); 	break;
	case 5:	menuOpen.file.Type("Wamit, Nemoh and FAST files", "*.out *.1 *.3 *.hst *.cal *.dat");	break;
	}
	menuOpen.file.AllFilesType();
	
	menuConvert.file.ClearTypes();
	switch (menuConvert.opt) {
	case 0:	menuConvert.file.Type("Wamit .1.3.hst file", "*.1 *.3 *.hst");	menuConvert.file = ForceExt(~menuConvert.file, ".1"); 	break;
	case 1:	menuConvert.file.Type("Wamit .out file", "*.out");				menuConvert.file = ForceExt(~menuConvert.file, ".out"); break;
	case 2:	menuConvert.file.Type("FAST HydroDyn file", "*.dat");			menuConvert.file = ForceExt(~menuConvert.file, ".dat"); break;
	}
	menuConvert.file.AllFilesType();
}

void Main::OnLoad() {
	String file = ~menuOpen.file;
	for (int i = 0; i < hydro.GetCount(); ++i) {
		if (hydro[i].file == file) {
			Exclamation("Model already loaded");
			return;
		}
	}
	String ext = GetFileExt(file);
	if (ext == ".cal") {
		Nemoh &data = hydro.Create<Nemoh>();
		if (!data.Load(file)) {
			Exclamation(DeQtfLf(Format("Problem loading '%s'\n%s", file, data.GetLastError())));	
			hydro.SetCount(hydro.GetCount()-1);
			return;
		}
	} else if (ext == ".inf") {
		Nemoh &data = hydro.Create<Nemoh>();
		if (!data.Load(file)) {
			Exclamation(DeQtfLf(Format("Problem loading '%s'\n%s", file, data.GetLastError())));	
			hydro.SetCount(hydro.GetCount()-1);
			return;
		}
	} else if (ext == ".out") {
		Wamit &data = hydro.Create<Wamit>();
		if (!data.Load(file, Null)) {
			Exclamation(DeQtfLf(Format("Problem loading '%s'\n%s", file, data.GetLastError())));		
			hydro.SetCount(hydro.GetCount()-1);
			return;
		}
		if (IsNull(data.rho)) {
			WithWamitLoad<TopWindow> dialog;
			CtrlLayoutOK(dialog, "Set Wamit data");
			dialog.rho <<= 1000;
			dialog.Execute();
			data.rho = ~dialog.rho;
		}
		
		data.Dimensionalize();
		
	} else if (ext == ".dat") {
		Fast &data = hydro.Create<Fast>();
		if (!data.Load(file, Null)) {
			Exclamation(DeQtfLf(Format("Problem loading '%s'\n%s", file, data.GetLastError())));		
			hydro.SetCount(hydro.GetCount()-1);
			return;
		}
		if (IsNull(data.g)) {
			WithFastLoad<TopWindow> dialog;
			CtrlLayoutOK(dialog, "Set Fast-Wamit data");
			dialog.g <<= 9.81;
			dialog.Execute();
			data.g = ~dialog.g;
		}
		
		data.Dimensionalize();
		
	} else if (ext == ".1" || ext == ".3" || ext == ".hst") {
		Fast &data = hydro.Create<Fast>();
		if (!data.Load(file, Null)) {
			Exclamation(DeQtfLf(Format("Problem loading '%s'\n%s", file, data.GetLastError())));		
			hydro.SetCount(hydro.GetCount()-1);
			return;
		}
		WithWamit13Load<TopWindow> dialog;
		CtrlLayoutOK(dialog, "Set Wamit data");
		data.g = 9.81;
		data.len = 1;
		data.rho = 1000;
		data.h = 100;
		CtrlRetriever rf;
		rf
			(dialog.g, data.g)
			(dialog.len, data.len)
			(dialog.rho, data.rho)
			(dialog.h, data.h)
		;
		dialog.Execute();
		rf.Retrieve();

		data.Dimensionalize();
		
	} else {
		Exclamation(Format("Unknown file extension in '%s'", DeQtf(file)));
		return;
	}
	if (menuOpen.optLoadIn == 0) {
		hydro.Remove(0, hydro.GetCount()-1);
		mainSummary.array.Reset();
	}
	String strError;
	if (!Hydro::MatchCoeffStructure(hydro, strError)) {
		hydro.SetCount(hydro.GetCount()-1);
		Exclamation(Format("Model '%s' does not match with previous: %s", DeQtf(file), strError));
	}
	
	int id = hydro.GetCount()-1;
	Hydro &data = hydro[id];
	data.Report();
	mainSummary.Report(data, id);
	if (data.Nh < 0 || data.Nf < 0)
		return;
	mainArrange.Load(hydro);
	menuConvert.arrayModel.Add(data.GetCodeStr(), data.name);
	if (menuConvert.arrayModel.GetCursor() < 0)
		menuConvert.arrayModel.SetCursor(0);
	mainTab.GetItem(mainTab.Find(mainArrange)).Enable(true);	
	mainTab.GetItem(mainTab.Find(mainA)).Enable(mainA.Load(hydro));	
	mainTab.GetItem(mainTab.Find(mainB)).Enable(mainB.Load(hydro));
	mainTab.GetItem(mainTab.Find(mainForceSC)).Enable(mainForceSC.Load(hydro));
	mainTab.GetItem(mainTab.Find(mainForceFK)).Enable(mainForceFK.Load(hydro));
	mainTab.GetItem(mainTab.Find(mainForceEX)).Enable(mainForceEX.Load(hydro));
}

void Main::OnConvert() {
	int id = menuConvert.arrayModel.GetCursor();
	if (id < 0) {
		Exclamation("Please select a model to export");
		return;
	}
	Hydro::BEM_SOFT type;	
	switch (menuConvert.opt) {
	case 0:	type = Hydro::WAMIT_1_3;	break;
	case 1:	type = Hydro::WAMIT;		break;
	case 2:	type = Hydro::FAST_WAMIT;	break;
	case 3:	type = Hydro::NEMOH;		break;
	default: throw Exc("Unknown type in OnConvert()");
	}
	hydro[id].Save(~menuConvert.file, type);	
}

Main::~Main() {
	if (!closed)
		Close();
}

void Main::Close(bool store) {
	if (store) 
		StoreSerializeJson();
	
	Thread::ShutdownThreads(); 
	RejectBreak(IDOK);		// Empty EditStrings does not disturb
	TopWindow::Close();
	closed = true;
}

void Main::Jsonize(JsonIO &json) {
	json
		("menuOpen_file", menuOpen.file)
		("menuOpen_opt", menuOpen.opt)
		("menuOpen_optLoadIn", menuOpen.optLoadIn)
		("menuConvert_file", menuConvert.file)
		("menuConvert_opt", menuConvert.opt)
		("menuPlot_autoFit", menuPlot.autoFit)
		("menuPlot_showPoints", menuPlot.showPoints)
	;
}

void MenuAbout::Init() {
	CtrlLayout(*this);
	
	String qtf = GetTopic(String("BEMRosetta/main/About$en-us")); 
	HelpHandler(qtf);
	info.SetQTF(qtf);
}

void MenuAbout::HelpHandler(String &str) {
	String name, mode;
	Time date;
	int version, bits;
	GetCompilerInfo(name, version, date, mode, bits);
	str.Replace(DeQtf("[Build Info]"), Format("%4d%02d%02d%02d, %s, %d bits", 
						date.year, date.month, date.day, date.hour, mode, bits)); 
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
	array.Set(row, 0, "File");				array.Set(row++, col, data.file);
	array.Set(row, 0, "Name");				array.Set(row++, col, data.name);
	array.Set(row, 0, "Soft");				array.Set(row++, col, data.GetCodeStr());
	array.Set(row, 0, "g [m/s2]");			array.Set(row++, col, data.g);
	array.Set(row, 0, "rho [kg/m3]");		array.Set(row++, col, data.rho);
	array.Set(row, 0, "h [m]");				array.Set(row++, col, data.h == INFINITY ? "INFINITY" : FormatDouble(data.h));
	array.Set(row, 0, "length scale [m]");	array.Set(row++, col, data.len);
	array.Set(row, 0, "#frequencies");		array.Set(row++, col, Nvl(data.Nf, 0));
	if (!data.w.IsEmpty()) {
		array.Set(row, 0, "freq_0 [rad/s]");	array.Set(row++, col, data.w[0]);
		array.Set(row, 0, "freq_end [rad/s]");	array.Set(row++, col, data.w[data.w.GetCount()-1]);
		array.Set(row, 0, "freq_delta [rad/s]");array.Set(row++, col, data.w[1]-data.w[0]);
	} else {
		array.Set(row, 0, "freq_0 [rad/s]");	array.Set(row++, col, "-");
		array.Set(row, 0, "freq_end [rad/s]");	array.Set(row++, col, "-");
		array.Set(row, 0, "freq_delta [rad/s]");array.Set(row++, col, "-");
	}
	array.Set(row, 0, "#headings");			array.Set(row++, col, Nvl(data.Nh, 0));
	if (!data.head.IsEmpty()) {
		array.Set(row, 0, "head_0 [º]");		array.Set(row++, col, data.head[0]);
	} else {
		array.Set(row, 0, "head_0 [º]");		array.Set(row++, col, "-");
	}
	if (data.head.GetCount() > 1) {
		array.Set(row, 0, "head_end [º]");	array.Set(row++, col, data.head[data.head.GetCount()-1]);
		array.Set(row, 0, "head_delta [º]");array.Set(row++, col, data.head[1] - data.head[0]);
	} else {
		array.Set(row, 0, "head_end [º]");	array.Set(row++, col, "-");
		array.Set(row, 0, "head_delta [º]");array.Set(row++, col, "-");
	}
	array.Set(row, 0, "A0 available");		array.Set(row++, col, data.IsLoadedAw0() ? "Yes" : "No");
	array.Set(row, 0, "Ainf available");	array.Set(row++, col, data.IsLoadedAwinf() ? "Yes" : "No");
	array.Set(row, 0, "A available");		array.Set(row++, col, data.IsLoadedA() ? "Yes" : "No");
	array.Set(row, 0, "B available");		array.Set(row++, col, data.IsLoadedB() ? "Yes" : "No");
	array.Set(row, 0, "C available");		array.Set(row++, col, data.IsLoadedC() ? "Yes" : "No");
	array.Set(row, 0, "Fex available");		array.Set(row++, col, data.IsLoadedFex() ? "Yes" : "No");
	array.Set(row, 0, "Fsc available");		array.Set(row++, col, data.IsLoadedFsc() ? "Yes" : "No");
	array.Set(row, 0, "Ffk available");		array.Set(row++, col, data.IsLoadedFfk() ? "Yes" : "No");
	array.Set(row, 0, "#bodies");			array.Set(row++, col, data.Nb);
	for (int ib = 0; ib < data.Nb; ++ib) {
		String sib = Format("#%d", ib+1);
		if (data.names.GetCount() > ib) {
			sib += " " + data.names[ib];
			array.Set(row, 0, sib + " Name");		array.Set(row++, col, data.names[ib]);
		} else {
			array.Set(row, 0, sib + " Name");		array.Set(row++, col, "-");
		}
		array.Set(row, 0, sib + " #dof");
		if (data.dof.GetCount() > ib) 
			array.Set(row++, col, data.dof[ib]);
		else 
			array.Set(row++, col, "-");
		
		array.Set(row, 0, sib + " Vsub [m3]");
		if (data.Vo.size() > ib) 
			array.Set(row++, col, data.Vo[ib]);
		else 
			array.Set(row++, col, "-");
		
		if (data.cg.size() > 3*ib) {
			array.Set(row, 0, sib + " Cg(x) [m]");	array.Set(row++, col, data.cg(0, ib));
			array.Set(row, 0, sib + " Cg(y) [m]");	array.Set(row++, col, data.cg(1, ib));
			array.Set(row, 0, sib + " Cg(z) [m]");	array.Set(row++, col, data.cg(2, ib));
		} else {
			array.Set(row, 0, sib + " Cg(x) [m]");	array.Set(row++, col, "-");
			array.Set(row, 0, sib + " Cg(y) [m]");	array.Set(row++, col, "-");
			array.Set(row, 0, sib + " Cg(z) [m]");	array.Set(row++, col, "-");
		}
		if (data.cb.size() > 3*ib) {
			array.Set(row, 0, sib + " Cb(x) [m]");	array.Set(row++, col, data.cb(0, ib));
			array.Set(row, 0, sib + " Cb(y) [m]");	array.Set(row++, col, data.cb(1, ib));
			array.Set(row, 0, sib + " Cb(z) [m]");	array.Set(row++, col, data.cb(2, ib));
		} else {
			array.Set(row, 0, sib + " Cb(x) [m]");	array.Set(row++, col, "-");
			array.Set(row, 0, sib + " Cb(y) [m]");	array.Set(row++, col, "-");
			array.Set(row, 0, sib + " Cb(z) [m]");	array.Set(row++, col, "-");
		}
		if (data.C.GetCount() > ib) {
			for (int i = 0; i < 6; ++i) {
				for (int j = 0; j < 6; ++j) {
					if (!Hydro::C_units(i, j).IsEmpty()) {
						array.Set(row, 0, sib + Format(" C(%d,%d) [%s]", i+1, j+1, Hydro::C_units(i, j)));	array.Set(row++, col, data.C[ib](i, j));		
					}
				}
			}
		} else {
			for (int i = 0; i < 6; ++i) {
				for (int j = 0; j < 6; ++j) {
					if (!Hydro::C_units(i, j).IsEmpty()) {
						array.Set(row, 0, sib + Format(" C(%d,%d) [%s]", i+1, j+1, Hydro::C_units(i, j)));	array.Set(row++, col, "-");		
					}
				}
			}
		}
	}	
}

void MainSummary::OnArrayBar(Bar &menu) 
{
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
	Print("BEMRosetta\nHydrodynamic coefficients viewer and converter for Boundary Element Method solver formats");
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

void MainArrange::Load(Upp::Array<Hydro> &hydro) {
	tab.Reset();
	arrangeDOF.Clear();
	for (int i = 0; i < hydro.GetCount(); ++i) {
		ArrangeDOF &arr = arrangeDOF.Add();
		arr.Init(hydro[i]);
		tab.Add(arr.SizePos(), Format("[%s] %s", hydro[i].GetCodeStr(), hydro[i].name));
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

bool MainABForce::Load(Upp::Array<Hydro> &hydro) {
	if (hydro.IsEmpty())
		return false;
	isFilling = true;
	tab.Reset();
	String format;
	switch (dataToShow) {
	case DATA_A:		format = "A%s%s";		break;		
	case DATA_B:		format = "B%s%s";		break;
	case DATA_FORCE_SC:	format = "Fsc%s%.1fº";	break;
	case DATA_FORCE_FK:	format = "Ffk%s%.1fº";	break;
	case DATA_FORCE_EX:	format = "Fex%s%.1fº";	break;
	}
	int sdof = 6*hydro[0].Nb;
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
		int Nh = hydro[0].Nh;
		if (Nh < 0)
			return false;
		plots.SetCount(Nh);
		for (int ih = 0; ih < Nh; ++ih) 
			plots[ih].SetCount(sdof);
		for (int i = 0; i < sdof; ++i) {
			for (int ih = 0; ih < Nh; ++ih) {
				plots[ih][i].Init(i, ih, hydro[0].head[ih], dataToShow);
				if (plots[ih][i].Load(hydro))
					tab.Add(plots[ih][i].SizePos(), Format(format, Hydro::StrDOFAbrev(i), hydro[0].head[ih]));
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

void MainPlot::Init(int i, int j_h, double h, DataToShow dataToShow) {
	CtrlLayout(*this);
	
	this->i = i;
	this->j_h = j_h;
	this->dataToShow = dataToShow;
	scatter.ShowAllMenus();
	String title, labelY;
	switch (dataToShow) {
	case DATA_A:		title = Format("Added mass A%s_%s", Hydro::StrDOF(i), Hydro::StrDOF(j_h));		
						labelY = "Added mass [Ns2/m]";				break;		
	case DATA_B:		title = Format("Radiation damping B%s_%s", Hydro::StrDOF(i), Hydro::StrDOF(j_h));
						labelY = "Radiation damping [Ns/m]";		break;
	case DATA_FORCE_SC:	title = Format("Diffraction scattering force %s heading %.1fº", Hydro::StrDOF(i), h);
						labelY = "Diffraction scattering force [N]";break;
	case DATA_FORCE_FK:	title = Format("Froude-Krylov force %s heading %.1fº", Hydro::StrDOF(i), h);
						labelY = "Froude-Krylov force [N]";			break;
	case DATA_FORCE_EX:	title = Format("Excitation Force %s heading %.1fº", Hydro::StrDOF(i), h);
						labelY = "Excitation Force [N]";			break;
						labelY = "Excitation Force [N]";			break;
	}
	scatter.SetTitle(title);
	scatter.SetLabelY(labelY);
	scatter.SetLabelX("w [rad/s]");
}

bool MainPlot::Load(Upp::Array<Hydro> &hydro) {
	scatter.RemoveAllSeries();
	ABF_source.SetCount(hydro.GetCount());
	Ainf_source.SetCount(hydro.GetCount());
	int markW = ma().menuPlot.showPoints ? 10 : 0;
	bool loaded = false;
	for (int id = 0; id < hydro.GetCount(); ++id) {	
		if (dataToShow == DATA_A) {
			Upp::Color acolor = Null;
			if (hydro[id].IsLoadedA()) {
				if (ABF_source[id].Init(hydro[id], i, j_h, PLOT_A)) {
					loaded = true;
					scatter.AddSeries(ABF_source[id]).Legend("A_" + hydro[id].name).SetMarkWidth(markW);
					double dummy;
					scatter.GetStroke(scatter.GetCount()-1, dummy, acolor);
				}
			}
			if (hydro[id].IsLoadedAwinf()) {
				if (Ainf_source[id].Init(hydro[id], i, j_h, PLOT_AINF)) {
					loaded = true;
					scatter.AddSeries(Ainf_source[id]).Legend("Ainf_" + hydro[id].name).Dash(LINE_DOTTED).Stroke(2, acolor).SetMarkWidth(markW);
				}
			}
		} else if (dataToShow == DATA_B && hydro[id].IsLoadedB()) {
			if (ABF_source[id].Init(hydro[id], i, j_h, PLOT_B)) {
				loaded = true;
				scatter.AddSeries(ABF_source[id]).Legend("B_" + hydro[id].name).SetMarkWidth(markW);
			}
		} else if (dataToShow == DATA_FORCE_SC && hydro[id].IsLoadedFsc()) {
			if (ABF_source[id].Init(hydro[id], i, j_h, PLOT_FORCE_SC_MA)) {
				loaded = true;
				scatter.AddSeries(ABF_source[id]).Legend("Fsc_ma_" + hydro[id].name).SetMarkWidth(markW);
			}
		} else if (dataToShow == DATA_FORCE_FK && hydro[id].IsLoadedFfk()) {
			if (ABF_source[id].Init(hydro[id], i, j_h, PLOT_FORCE_FK_MA)) {
				loaded = true;
				scatter.AddSeries(ABF_source[id]).Legend("Ffk_ma_" + hydro[id].name).SetMarkWidth(markW);
			}
		} else if (dataToShow == DATA_FORCE_EX && hydro[id].IsLoadedFex()) {
			if (ABF_source[id].Init(hydro[id], i, j_h, PLOT_FORCE_EX_MA)) {
				loaded = true;
				scatter.AddSeries(ABF_source[id]).Legend("Fex_ma_" + hydro[id].name).SetMarkWidth(markW);
			}
		}
	}
	if (ma().menuPlot.autoFit)
		scatter.ZoomToFit(true, true);
	return loaded;
}

Main &ma(Main *m) {
	static Main *mp = 0;
	if (m)
		mp = m;
	return *mp;
}

void Assert(const char *s) {
	NotPanic();
	throw Exc(s);
}
	
GUI_APP_MAIN {
	Ctrl::SetAppName(t_("Hydrodynamic coefficents viewer and converter"));
	
	String errorStr;
	try {
		SetAssertFailedHook(Assert);
		
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
		Exclamation(t_("Internal error:") + String("&") + DeQtf(errorStr) + String("&") + t_("Program ended"));
}

