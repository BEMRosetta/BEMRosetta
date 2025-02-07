// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include <CtrlLib/CtrlLib.h>

#define IMAGECLASS Img
#define IMAGEFILE <BEMRosetta/main.iml>
#include <Draw/iml.h>

#define TOPICFILE <BEMRosetta/main.tpp/all.i>
#include <Core/topic_group.h>

#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <TabBar/TabBar.h>
#include <DropGrid/DropGrid.h>


using namespace Upp;

#include <BEMRosetta_cl/BEMRosetta.h>
#include <BEMRosetta_cl/functions.h>

#include "main.h"


using namespace Eigen;

FreqSelector::FreqSelector() {
	add.SetImage(Img::add());
	add << [&] {AddField();};
	Add(add.LeftPos(0, 19).TopPos(3, 19));
};

void Main::Init() {
	BEM::Print 	  	  = [this](String s) {printf("%s", ~s); mainOutput.Print(s);};
	BEM::PrintWarning = [this](String s) {printf("%s", ~s); mainOutput.Print(s); /*Status(s);*/};
	BEM::PrintError   = [this](String s) {printf("%s", ~s); /*tab.Set(mainOutput); */Status(s);	Exclamation(DeQtfLf(s));};
	
	LOG("Init");
	Sizeable().Zoomable().SetMinSize(Size(800, 600));
	Icon(Img::Rosetta64());
	LargeIcon(Img::Rosetta256());
	Ma(this);
	CtrlLayout(*this);
	
	String name, mode;
	Time date;
	int version, bits;
	GetCompilerInfo(name, version, date, mode, bits);
	
	SetLanguage(GetSystemLNG());
	String sdate = Format("%mon %d", date.month, date.year);
	
	SetLanguage(LNG_('E', 'N', 'U', 'S'));
	
	String title;
	if (parameter == "")
		title = "BEMRosetta";
	else if (parameter == "bem")
		title = "BEMRosetta BEM models viewer";
	else if (parameter == "mesh")
		title = "BEMRosetta mesh files viewer";
	else if (parameter == "time")
		title = "BEMRosetta Time domain results viewer";
	else
		throw Exc(Format(t_("Unknown -gui parameter %s"), parameter));
		
	Title(title + " " + sdate + (Bem().experimental ? " EXPERIMENTAL" : ""));

	rectangleButtons.SetBackground(SColorFace());
	
	String errorJson;
	bool firstTime;
	
	mainSolver.InitBeforeSerialize();
	
	errorJson = Bem().LoadSerializeJson();
	firstTime = !errorJson.IsEmpty();
	if (firstTime) {
		String str = errorJson + "\n" + t_("BEM config. data is not loaded. Defaults values are set"); 
		if (errorJson != t_("First time")) 
			Exclamation(DeQtfLf(str + "\n" + t_("Please send the above file to the authors.\nCheck it beforehand in case protected data is included.")));
		LOG(str);
	}
	
	LOG("BEM configuration loaded");
	
	if (!Bem().ClearTempFiles()) 
		Cout() << "\n" << t_("BEM temporary files folder cannot be created");
	
	bool openOptions = false;
	errorJson = LoadSerializeJson(firstTime, openOptions);
	if (!errorJson.IsEmpty()) {
		String str = errorJson + "\n" + t_("Config. data is not loaded. Defaults values are set"); 
		if (errorJson != t_("First time")) 
			Exclamation(DeQtfLf(str + "\n" + t_("Please send the above file to the authors.\nCheck it beforehand in case protected data is included.")));
		LOG(str);
	}

	LOG("Configuration loaded");
	
	if (parameter.IsEmpty() || parameter == "mesh") {
		mainBody.Init();			LOG("Init Body");
		tab.Add(mainBody.SizePos(), t_("Vessel Mesh"));
	}
	if (parameter.IsEmpty()) {
		mainSolver.Init();		LOG("Init Nemoh");
		//tab.Add(mainNemohScroll.AddPaneV(mainSolver).SizePos(), tabTexts[TAB_NEMOH]);
		tab.Add(mainSolver.SizePos(), t_("BEM Solver"));
	}
	
	if (parameter.IsEmpty() || parameter == "bem") {
		mainBEM.Init();				LOG("Init BEM");
		tab.Add(mainBEM.SizePos(),  t_("Hydro Coeff"));
	}
	
	if (parameter.IsEmpty() || parameter == "time") {
		mainFAST.Init(GetBEMRosettaDataFolder(), bar);	LOG("Init FAST");
		tab.Add(mainFAST.SizePos(), t_("Time series"));
	}
	
	if (parameter.IsEmpty()) {
		mainTools.Init();
		tab.Add(mainTools.SizePos(), t_("Tools"));
	}
	
	if (parameter.IsEmpty()) {
		mainMoor.Init();	LOG("Init Moor");
#ifdef flagDEBUG
			tab.Add(mainMoor.SizePos(), t_("Mooring"));
#else
		if (Bem().experimental)
			tab.Add(mainMoor.SizePos(), t_("Mooring"));
#endif
	}
	
	if (parameter.IsEmpty()) {
		mainDecay.Init();	LOG("Init Decay");
		if (false/*Bem().experimental*/)
			tab.Add(mainDecay.SizePos(), t_("Decay"));
	}
	
	if (parameter.IsEmpty())
		tab.Add().Disable();
	
	mainOutput.Init();			LOG("Init Output");
	tab.Add(mainOutput.SizePos(), t_("Output"));	
	
	menuOptions.Init(Bem());		LOG("Init Options");
	menuOptions.Load();			LOG("Init Options.Load");
	tab.Add(menuOptionsScroll.AddPaneV(menuOptions).SizePos(), t_("Options"));
	
	menuAbout.Init();			LOG("Init About");
	tab.Add(menuAbout.SizePos(), t_("About"));
	
	editrho.OnLeftDown = [&](Point, dword) {tab.Set(menuOptionsScroll); menuOptions.rho.SetFocus(); menuOptions.rho.Underline(1);};
	editrho.SetReadOnly();
	editrho <<= Bem().rho;
	
	editg.OnLeftDown = [&](Point, dword) {tab.Set(menuOptionsScroll); menuOptions.g.SetFocus(); menuOptions.g.Underline(1);};
	editg.SetReadOnly();
	editg <<= Bem().g;
	
	editdofType.OnLeftDown = [&](Point, dword) {tab.Set(menuOptionsScroll); menuOptions.dofType.SetFocus(); menuOptions.dofType.Underline(1);};
	editdofType.SetReadOnly();
	editdofType <<= BasicBEM::strDOFType[Bem().dofType];
	
	editHeadingType.OnLeftDown = [&](Point, dword) {tab.Set(menuOptionsScroll); menuOptions.headingType.SetFocus(); menuOptions.headingType.Underline(1);};
	editHeadingType.SetReadOnly();
	editHeadingType <<= BasicBEM::strHeadingType[Bem().headingType];
	
	butWindow.SetImage(Img::application_double()).Tip(t_("Open new window"));
	butWindow.Hide();

	butWindow << [&] {
		if (tab.IsAt(mainBody)) 
			LaunchCommand(Format("\"%s\" -gui mesh", GetExeFilePath()));
		else if (tab.IsAt(mainBEM)) 
			LaunchCommand(Format("\"%s\" -gui bem", GetExeFilePath()));
		else if (tab.IsAt(mainFAST)) 
			LaunchCommand(Format("\"%s\" -gui time", GetExeFilePath()));
	};

	tab.WhenSet = [&] {
		LOGTAB(tab);
		if (tab.IsAt(menuOptionsScroll)) 
			menuOptions.Load();
		else if (tab.IsAt(mainSolver))			// mainNemohScroll)) 
			mainSolver.Load();
		if (!tab.IsAt(menuOptionsScroll) && menuOptions.IsChanged()) {
			if (PromptYesNo(t_("Options have changed&Do you want to save them?")))
				menuOptions.OnSave();
			else
				menuOptions.Load();
		}
		
		if (tab.IsAt(mainBody)) {
			butWindow.Show(true);
			lastTabS = tab.GetItem(tab.Get()).GetText();
		} else if (tab.IsAt(mainBEM)) {
			butWindow.Show(true);
			lastTabS = tab.GetItem(tab.Get()).GetText();
		} else 	if (tab.IsAt(mainSolver)) {		//	mainNemohScroll)) {
			butWindow.Show(false);
			lastTabS = tab.GetItem(tab.Get()).GetText();
		} else 	if (tab.IsAt(mainMoor)) {
			butWindow.Show(false);
			lastTabS = tab.GetItem(tab.Get()).GetText();
		} else 	if (tab.IsAt(mainDecay)) {
			butWindow.Show(false);
			lastTabS = tab.GetItem(tab.Get()).GetText();
		} else if (tab.IsAt(mainFAST)) {
			butWindow.Show(true);
			lastTabS = tab.GetItem(tab.Get()).GetText();
		} else 
			butWindow.Show(false);
			
		if (tab.IsAt(mainBody)) 
			mainBody.mainTab.WhenSet();
		else if (tab.IsAt(mainBEM)) 
			mainBEM.mainTab.WhenSet();
	};	
	
	if (firstTime) {
		tab.Set(menuAbout);
	} else if (openOptions)
		tab.Set(menuOptionsScroll);
	else 
		SetLastTab();
	
	tab.WhenSet();
	
	AddFrame(bar);
}

void Main::OptionsUpdated(double rho, double g, int dofType, int headingType) {
	mainBEM.OnOpt();
	mainBody.OnOpt();
	
	editg <<= g;
	editrho <<= rho;
	editdofType <<= BasicBEM::strDOFType[dofType];
	editHeadingType <<= BasicBEM::strHeadingType[headingType];
	
	mainBEM.UpdateButtons();
}

String Main::LoadSerializeJson(bool &firstTime, bool &openOptions) {
	String ret;
	String folder = GetBEMRosettaDataFolder();
	if (!DirectoryCreateX(folder))
		ret = Format(t_("Impossible to create folder '%s' to store configuration file"), folder);
	else {
		String fileName = AFX(folder, "config.cf");
		if (!FileExists(fileName)) 
			ret = t_("First time");
		else {
			String jsonText = LoadFile(fileName);
			if (jsonText.IsEmpty())
				ret = Format(t_("Configuration file '%s' is empty"), fileName);
			else {
				ret = LoadFromJsonError(*this, jsonText);
				if (!ret.IsEmpty())
					ret = Format(t_("Problem loading configuration file '%s': %s"), fileName, ret);
			}
			if (!ret.IsEmpty()) {
				DirectoryCreateX(AFX(folder, "Errors"));
				FileCopy(fileName, AFX(folder, "Errors", "config.cf"));
			}
		}
	}
	
	firstTime = !ret.IsEmpty();
	
	mainBody.InitSerialize(!firstTime);
	mainSolver.InitAfterSerialize(!firstTime);
	mainBEM.InitSerialize(!firstTime);
	menuOptions.InitSerialize(!firstTime, openOptions);
	
	if (firstTime)
		tab.Set(menuAbout);
	
	return ret;
}

bool Main::StoreSerializeJson() {
	String folder = GetBEMRosettaDataFolder();
	if (!DirectoryCreateX(folder))
		return false;
	String fileName = AFX(folder, "config.cf");
	return StoreAsJsonFile(*this, fileName, true);
}

Main::~Main() noexcept {
	if (!closed)
		CloseMain(false);
}

void Main::Close() {
	CloseMain(!closed);
}

void Main::CloseMain(bool store) {
	if (store) {
		Bem().StoreSerializeJson();
		StoreSerializeJson();
	}
	Thread::ShutdownThreads(); 
	RejectBreak(IDOK);		// Empty EditStrings does not disturb
	try {
		TopWindow::Close();
	} catch(...) {
	}
	closed = true;
}

void Main::Jsonize(JsonIO &json) {
	json
		("mesh", mainBody)
		("nemoh", mainSolver)
		("bem", mainBEM)
		("mooring", mainMoor)
		("decay", mainDecay)
		("lastTabS", lastTabS)
	;
}


void MenuAbout::Init() {
	CtrlLayout(*this);
	
	String qtf = GetTopic(S("BEMRosetta/main/About$en-us")); 
	SetBuildInfo(qtf);
	qtf.Replace("SYSTEMINFO", DeQtf(GetSystemInfo()));
	info.SetQTF(qtf);
}

			
Main &Ma(Main *m) {
	static Main *mp = 0;
	if (m)
		mp = m;
	if (!mp)
		throw Exc(t_("Main is not initialized"));	
	return *mp;
}

void Status(String str, int time)	{Ma().Status(str, time);}

String TabText(const TabCtrl &tab) {
	int id = tab.Get();
	if (id < 0)
		return String();
	return tab.GetItem(id).GetText();
}

const Color &GetColorId(int id) {
	static Color col[10] = {Color(93,185,86), Color(180,96,189), Color(176,177,67), Color(106,126,205),
							Color(219,142,68), Color(79,186,172), Color(203,83,67), Color(83,133,68),
							Color(199,89,128), Color(147,119,58)};
	int mod = id % 10;
	return col[mod];
}

struct RectDisplay : public Display {
	virtual void Paint(Draw& w, const Rect& r, const Value& q,
	                   Color , Color , dword ) const {
		w.DrawRect(r, q);
	}
};

ArrayCtrl &ArrayModel_Init(ArrayCtrl &array, bool option, const UVector<int> &w) {
	static const UVector<int> ww = {30, 30, 30};
	const UVector<int> *pww = &w;
	if (w.IsEmpty())
		pww = &ww;
	
	array.Reset();
	array.NoHeader().NoVertGrid().AutoHideSb();
	array.HeaderObject().HideTab(array.AddColumn().HeaderTab().GetIndex());
	array.AddColumn("", 5).SetDisplay(Single<RectDisplay>());
	if (option)
		array.AddColumn("", 10);
	else
		array.AddColumn("", 0);
	array.AddColumn("", (*pww)[0]);	
	array.AddColumn("", (*pww)[1]);
	array.AddColumn("", (*pww)[2]);
	array.HeaderTab(1).SetMargin(-2);
	if (option)
		array.HeaderTab(2).SetMargin(-2);
	return array;
}

void ArrayModel_Add(ArrayCtrl &array, String codeStr, String title, String fileName, int id, 
					UArray<Option> &option, Function <void()>OnPush) {
	if (id < 0)
		throw Exc(t_("Wrong body id added"));
	array.Add(id, GetColorId(id), true, codeStr, title, fileName);
	int row = array.GetCount()-1;
	Option & opt = option.Add();
	array.SetCtrl(row, 2, opt);
	opt << OnPush;
	array.SetCursor(array.GetCount());
}

void ArrayModel_Add(ArrayCtrl &array, String codeStr, String title, String fileName, int id) {
 	array.Add(id, GetColorId(id), Null, codeStr, title, fileName);
 	array.SetCursor(array.GetCount());
}

void ArrayModel_Change(ArrayCtrl &array, int id, String codeStr, String title, String fileName) {
	for (int row = 0; row < array.GetCount(); ++row) {
		if (array.Get(row, 0) == id) {
			if (!IsNull(codeStr))
				array.Set(row, 3, codeStr);
			if (!IsNull(title))
				array.Set(row, 4, title);
			if (!IsNull(fileName))
				array.Set(row, 5, fileName);
			return;
		}
	}
	throw Exc(Format(t_("Id %d not found in ArrayModel_Change()"), id));
}

int ArrayModel_Id(const ArrayCtrl &array) {
	int row = array.GetCursor();
	if (row < 0)
		return -1;
	return array.Get(row, 0);
}

int ArrayModel_Id(const ArrayCtrl &array, int row) {
	return array.Get(row, 0);
}

int ArrayModel_IndexBody(const ArrayCtrl &array) {
	int id;
	if (array.GetCount() == 1)
		id = ArrayModel_Id(array, 0);
	else {
		id = ArrayModel_Id(array);
		if (id < 0)
			return -1;
	}
	return Bem().GetBodyIndex(id);
}

int ArrayModel_IndexBody(const ArrayCtrl &array, int row) {
	int id = ArrayModel_Id(array, row);
	if (id < 0)
		return -1;
	return Bem().GetBodyIndex(id);
}

int ArrayModel_IndexHydro(const ArrayCtrl &array) {
	int id = ArrayModel_Id(array);
	if (id < 0)
		return -1;
	return Bem().GetHydroIndex(id);
}

int ArrayModel_IndexHydro(const ArrayCtrl &array, int row) {
	int id = ArrayModel_Id(array, row);
	if (id < 0)
		return -1;
	return Bem().GetHydroIndex(id);
}

UVector<int> ArrayModel_IndexsHydro(const ArrayCtrl &array) {		
	UVector<int> ids;
	for (int row = 0; row < array.GetCount(); ++row) {
		int id = ArrayModel_IndexHydro(array, row);		
		if (id >= 0)
			ids << id;
	}
	return ids;
}

UVector<int> ArrayModel_IndexsBody(const ArrayCtrl &array) {		
	UVector<int> ids;
	for (int row = 0; row < array.GetCount(); ++row) {
		int id = ArrayModel_IndexBody(array, row);		
		if (id >= 0)
			ids << id;
	}
	return ids;
}

void ArrayModel_IdsHydroDel(ArrayCtrl &array, const UVector<int> &idxs) {		
	for (int row = array.GetCount() - 1; row >= 0 ; --row) {
		int idxrow = ArrayModel_IndexHydro(array, row);
		for (auto idx : idxs) {
			if (idxrow == idx) {
				array.Remove(row);
				break;
			}
		}
	}
}

void ArrayModel_RowsHydroDel(ArrayCtrl &array, const UVector<int> &rows) {		
	for (int row = array.GetCount() - 1; row >= 0 ; --row) {
		for (auto rw : rows) {
			if (row == rw) {
				array.Remove(row);
				break;
			}
		}
	}
}

bool ArrayModel_IsVisible(const ArrayCtrl &array, int row) {
	return array.Get(row, 2);
}

bool ArrayModel_IsSelected(const ArrayCtrl &array, int row) {
	return array.IsSelected(row);
}

const Color& ArrayModel_GetColor(const ArrayCtrl &array, int row) {
	return GetColorId(array.Get(row, 0));
}

String ArrayModel_GetTitle(ArrayCtrl &array, int row) {
	if (row < 0) 
		row = array.GetCursor();
	if (row < 0)
		return String();
 	return array.Get(row, 4);
}

String ArrayModel_GetFileName(ArrayCtrl &array, int row) {
	if (row < 0) 
		row = array.GetCursor();
	if (row < 0)
		return String();
	return array.Get(row, 5);
}


GUI_APP_MAIN {
#if defined(PLATFORM_WIN32) 
	GetCrashHandler().Enable();
	#ifndef flagDEBUG
	if (EM().Init("BEMRosetta", "BEMRosetta", EM().DefaultExitError, Null))
		return;
	InstallPanicMessageBox([](const char *title, const char *text) {
		EM().Log(Format("%s: %s", title, text));
		throw Exc(text);
	});
	#endif
#endif
	
	const UVector<String>& command = CommandLine();

	String errorStr;

	if (!command.IsEmpty() && command[0] != "-gui") {
		try {
			ConsoleOutput con(true);
			
			ConsoleMain(command, true, PrintStatus);
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
		if (!errorStr.IsEmpty()) {
			Cout() << "\n" << Format(t_("Problem found: %s"), errorStr);
			SetExitCode(-1);
		}
		return;
	}
	
	Ctrl::SetAppName(t_("Hydrodynamic coefficients viewer and converter"));
	Ctrl::GlobalBackPaint();

	try {
		Main main;
		if (command.size() >= 1) {
			if (command[0] == "-gui") {
				if (command.size() < 2)
					throw Exc("-gui option requires to indicate window to show");
				main.parameter = ToLower(command[1]);
			} else
				throw Exc(Format("Unknown command %s", command[0]));
		}
		main.Init();
		main.OpenMain();
		
		Ctrl::EventLoop();
		main.Close();
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
		Exclamation(t_("Internal error:") + S("&") + DeQtf(errorStr) + S("&") + t_("Program ended"));
}
