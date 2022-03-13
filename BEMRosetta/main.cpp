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

#include "main.h"


using namespace Eigen;

FreqSelector::FreqSelector() {
	add.SetImage(Img::add());
	add << [&] {AddField();};
	Add(add.LeftPos(0, 19).TopPos(3, 19));
};

static String GetBEMRosettaDataFolder() {
	return AppendFileNameX(GetAppDataFolder(), "BEMRosetta");
}

void Main::Init() {
	LOG("Init");
	Sizeable().Zoomable().SetMinSize(Size(800, 600));
	Icon(Img::Rosetta64());
	LargeIcon(Img::Rosetta256());
	ma(this);
	CtrlLayout(*this);
	
	String name, mode;
	Time date;
	int version, bits;
	GetCompilerInfo(name, version, date, mode, bits);
	
	SetLanguage(GetSystemLNG());
	String sdate = Format("%mon %d", date.month, date.year);
	
	SetLanguage(LNG_('E', 'N', 'U', 'S'));
	
	Title(S("BEMRosetta") + " " + sdate + (Bem().experimental ? " EXPERIMENTAL" : ""));

	rectangleButtons.SetBackground(SColorFace());
	
	tabTexts << t_("Mesh Handling") << t_("BEM Solver") << t_("Hydrodynamic Coefficients") 
			 << t_("Mooring") << t_("Decay") << t_("FAST .out+b Reader");
		
	bool firstTime = !bem.LoadSerializeJson();
	if (firstTime) 
		Cout() << "\n" << t_("BEM configuration data is not loaded. Defaults values are set");
	
	LOG("BEM configuration loaded");
	
	if (IsNull(lastTab))
		lastTab = 0;
	
	if (!bem.ClearTempFiles()) 
		Cout() << "\n" << t_("BEM temporary files folder cannot be created");
	
	bool openOptions = false;
	if (!LoadSerializeJson(firstTime, openOptions)) 
		Cout() << "\n" << t_("Configuration data is not loaded. Defaults are set");
	
	LOG("Configuration loaded");
	
	if (menuOptions.showTabMesh) {
		mainMesh.Init();			LOG("Init Mesh");
		tab.Add(mainMesh.SizePos(), tabTexts[TAB_MESH]);
	}
	if (menuOptions.showTabNemoh) {
		mainSolver.Init(bem);		LOG("Init Nemoh");
		//tab.Add(mainNemohScroll.AddPaneV(mainSolver).SizePos(), tabTexts[TAB_NEMOH]);
		tab.Add(mainSolver.SizePos(), tabTexts[TAB_NEMOH]);
	}
	if (menuOptions.showTabCoeff) {
		mainBEM.Init();				LOG("Init BEM");
		tab.Add(mainBEM.SizePos(),  tabTexts[TAB_COEFF]);
	}
	if (menuOptions.showTabMoor) {
		tab.Add().Disable();
		mainMoor.Init();	LOG("Init Moor");
		if (false/*Bem().experimental*/)
			tab.Add(mainMoor.SizePos(), tabTexts[TAB_MOOR]);
	}
	if (menuOptions.showTabDecay) {
		if (!menuOptions.showTabMoor)
			tab.Add().Disable();
		mainDecay.Init();	LOG("Init Decay");
		if (false/*Bem().experimental*/)
			tab.Add(mainDecay.SizePos(), tabTexts[TAB_DECAY]);
	}
	if (menuOptions.showTabFAST) {
		if (!menuOptions.showTabMoor && !menuOptions.showTabDecay)	
			tab.Add().Disable();
		mainFAST.Init(GetBEMRosettaDataFolder(), bar);	LOG("Init FAST");
		tab.Add(mainFAST.SizePos(), tabTexts[TAB_FAST]);
	}
	
	tab.Add().Disable();
	mainOutput.Init();			LOG("Init Output");
	tab.Add(mainOutput.SizePos(), t_("Output"));	
	menuOptions.Init(bem);		LOG("Init Options");
	menuOptions.Load();			LOG("Init Options.Load");
	tab.Add(menuOptionsScroll.AddPaneV(menuOptions).SizePos(),t_("Options"));
	tab.Add().Disable();
	menuAbout.Init();			LOG("Init About");
	tab.Add(menuAbout.SizePos(), t_("About"));
	
	editrho.OnLeftDown = [&](Point, dword) {tab.Set(menuOptionsScroll); menuOptions.rho.SetFocus(); menuOptions.rho.Underline(1);};
	editrho.SetReadOnly();
	editrho <<= bem.rho;
	
	editg.OnLeftDown = [&](Point, dword) {tab.Set(menuOptionsScroll); menuOptions.g.SetFocus(); menuOptions.g.Underline(1);};
	editg.SetReadOnly();
	editg <<= bem.g;
	
	editdofType.OnLeftDown = [&](Point, dword) {tab.Set(menuOptionsScroll); menuOptions.dofType.SetFocus(); menuOptions.dofType.Underline(1);};
	editdofType.SetReadOnly();
	editdofType <<= BEM::strDOFType[bem.dofType];
	
	editHeadingType.OnLeftDown = [&](Point, dword) {tab.Set(menuOptionsScroll); menuOptions.headingType.SetFocus(); menuOptions.headingType.Underline(1);};
	editHeadingType.SetReadOnly();
	editHeadingType <<= BEM::strHeadingType[bem.headingType];
	
	butWindow.SetImage(Img::application_double()).SetLabel(t_("New window")).Tip(t_("Open new window"));
	butWindow.Hide();
	
	butWindow << [&] {
		if (tab.IsAt(mainMesh)) {
			AddWindow();
			MainMeshW *mainMeshW = new MainMeshW();
			mainMeshW->Init(mainMesh, Img::Rosetta64(), Img::Rosetta256(), [&]() {DeleteWindow();});
			mainMeshW->OpenMain();
		} else if (tab.IsAt(mainBEM)) {
			AddWindow();
			MainBEMW *mainBEMW = new MainBEMW();
			mainBEMW->Init(mainBEM, Img::Rosetta64(), Img::Rosetta256(), [&]() {DeleteWindow();});
			mainBEMW->OpenMain();
		} else if (tab.IsAt(mainFAST)) {
			AddWindow();
			MainFASTW *mainFASTW = new MainFASTW();
			mainFASTW->Init(GetBEMRosettaDataFolder(), Img::Rosetta64(), Img::Rosetta256(), bar, [&]() {DeleteWindow();});
			mainFASTW->OpenMain();
		}
	};
	
	tab.WhenSet = [&] {
		LOGTAB(tab);
		if (tab.IsAt(menuOptionsScroll)) 
			menuOptions.Load();
		else if (tab.IsAt(mainSolver))			// mainNemohScroll)) 
			mainSolver.Load(bem);
		if (!tab.IsAt(menuOptionsScroll) && menuOptions.IsChanged()) {
			if (PromptYesNo(t_("Options have changed&Do you want to save them?")))
				menuOptions.OnSave();
			else
				menuOptions.Load();
		}
		
		if (tab.IsAt(mainMesh)) {
#ifdef PLATFORM_POSIX
			butWindow.Show(false);
#else
			butWindow.Show(true);
#endif
			lastTab = ~tab;
		} else if (tab.IsAt(mainBEM)) {
			butWindow.Show(true);
			lastTab = ~tab;
		} else 	if (tab.IsAt(mainSolver)) {		//	mainNemohScroll)) {
			butWindow.Show(false);
			lastTab = ~tab;
		} else 	if (tab.IsAt(mainMoor)) {
			butWindow.Show(false);
			lastTab = ~tab;
		} else 	if (tab.IsAt(mainDecay)) {
			butWindow.Show(false);
			lastTab = ~tab;
		} else if (tab.IsAt(mainFAST)) {
			butWindow.Show(true);
			lastTab = ~tab;
		} else 
			butWindow.Show(false);
			
		if (tab.IsAt(mainMesh)) 
			mainMesh.mainTab.WhenSet();
		else if (tab.IsAt(mainBEM)) 
			mainBEM.mainTab.WhenSet();
	};	
	
	if (firstTime) {
		lastTab = 0;
		tab.Set(menuAbout);
	} else if (openOptions)
		tab.Set(menuOptionsScroll);
	else 
		SetLastTab();
	
	tab.WhenSet();
	
	AddFrame(bar);
	
	BEM::Print 	  	  = [this](String s) {printf("%s", ~s); mainOutput.Print(s);};
	BEM::PrintWarning = [this](String s) {printf("%s", ~s); mainOutput.Print(s); /*Status(s);*/};
	BEM::PrintError   = [this](String s) {printf("%s", ~s); mainOutput.Print(s); tab.Set(mainOutput); Status(s);};
}

void Main::OptionsUpdated(double rho, double g, int dofType, int headingType) {
	mainBEM.OnOpt();
	mainMesh.OnOpt();
	
	editg <<= g;
	editrho <<= rho;
	editdofType <<= BEM::strDOFType[dofType];
	editHeadingType <<= BEM::strHeadingType[headingType];
}

bool Main::LoadSerializeJson(bool &firstTime, bool &openOptions) {
	bool ret;
	String folder = GetBEMRosettaDataFolder();
	if (!DirectoryCreateX(folder))
		ret = false;
	else {
		String fileName = AppendFileNameX(folder, "config.cf");
		if (!FileExists(fileName)) 
			ret = false;
		else {
			String jsonText = LoadFile(fileName);
			if (jsonText.IsEmpty())
				ret = false;
			else {
				if (!LoadFromJson(*this, jsonText))
					ret = false;
				else
					ret = true;
			}
		}
	}
	
	mainMesh.InitSerialize(ret);
	mainSolver.InitSerialize(ret);
	mainBEM.InitSerialize(ret);
	menuOptions.InitSerialize(ret, openOptions);
	
	if (!ret)
		tab.Set(menuAbout);
	
	firstTime = !ret;
	
	return ret;
}

bool Main::StoreSerializeJson() {
	String folder = GetBEMRosettaDataFolder();
	if (!DirectoryCreateX(folder))
		return 0;
	String fileName = AppendFileNameX(folder, "config.cf");
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
	if (numWindows > 0) {
		if (!PromptOKCancel(t_("There windows opened.&Do you want to close anyway?")))	
			return;
	}
	if (store) {
		bem.StoreSerializeJson();
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
		("mesh", mainMesh)
		("nemoh", mainSolver)
		("bem", mainBEM)
		("mooring", mainMoor)
		("decay", mainDecay)
		("menuOptions", menuOptions)
		("lastTab", lastTab)
	;
}

void MenuOptions::Init(BEM &_bem) {
	CtrlLayout(*this);
	
	bem = &_bem;
	butSave  <<= THISBACK(OnSave);
	butSave2 <<= THISBACK(OnSave);
	
	g.isMouseEnter = rho.isMouseEnter = dofType.isMouseEnter = headingType.isMouseEnter = false;
	
	arrayShown.AddColumn("");
	arrayShown.Add();	arrayShown.CreateCtrl<Option>(Main::TAB_MESH,  0, false).SetLabel(ma().tabTexts[Main::TAB_MESH]);
	arrayShown.Add();	arrayShown.CreateCtrl<Option>(Main::TAB_NEMOH, 0, false).SetLabel(ma().tabTexts[Main::TAB_NEMOH]);
	arrayShown.Add();	arrayShown.CreateCtrl<Option>(Main::TAB_COEFF, 0, false).SetLabel(ma().tabTexts[Main::TAB_COEFF]);
	arrayShown.Add();	arrayShown.CreateCtrl<Option>(Main::TAB_MOOR,  0, false).SetLabel(ma().tabTexts[Main::TAB_MOOR]).Show(Bem().experimental);
	arrayShown.Add();	arrayShown.CreateCtrl<Option>(Main::TAB_DECAY, 0, false).SetLabel(ma().tabTexts[Main::TAB_DECAY]).Show(Bem().experimental);
	arrayShown.Add();	arrayShown.CreateCtrl<Option>(Main::TAB_FAST,  0, false).SetLabel(ma().tabTexts[Main::TAB_FAST]);
	
	for (int i = 0; BEM::strDOFType[i][0] != '\0'; ++i)
		dofType.Add(BEM::strDOFType[i]);
	for (int i = 0; BEM::strHeadingType[i][0] != '\0'; ++i)
		headingType.Add(BEM::strHeadingType[i]);
}

void MenuOptions::InitSerialize(bool ret, bool &openOptions) {
	if (!ret || IsNull(showTabMesh)) 
		showTabMesh = true;
	if (!ret || IsNull(showTabNemoh)) 
		showTabNemoh = true;
	if (!ret || IsNull(showTabCoeff)) 
		showTabCoeff = true;
	if (!ret || IsNull(showTabMoor)) 
		showTabMoor = true;
	if (!ret || IsNull(showTabDecay)) 
		showTabDecay = true;
	if (!ret || IsNull(showTabFAST)) 
		showTabFAST = true;
	if (!showTabMesh && !showTabNemoh && !showTabCoeff && !showTabFAST && !showTabMoor && !showTabDecay) {
		openOptions = true;
		Exclamation(t_("[* No tab selected to be shown]&To show them, choose them in [* 'Options/General/Tabs shown']"));
	}
}

void MenuOptions::Load() {
	g <<= bem->g;
	rho <<= bem->rho;
	len <<= bem->len;
	depth <<= bem->depth;
	//discardNegDOF <<= bem->discardNegDOF;
	//thres <<= bem->thres;
	calcAinf <<= bem->calcAinf;
	calcAinf_w <<= bem->calcAinf_w;
	maxTimeA <<= bem->maxTimeA;
	numValsA <<= bem->numValsA;	
	onlyDiagonal <<= bem->onlyDiagonal;
	nemohPathPreprocessor <<= bem->nemohPathPreprocessor;
	nemohPathSolver <<= bem->nemohPathSolver;
	nemohPathPostprocessor <<= bem->nemohPathPostprocessor;
	nemohPathNew <<= bem->nemohPathNew;
	nemohPathGREN <<= bem->nemohPathGREN;
	foammPath <<= bem->foammPath;
	hamsPath <<= bem->hamsPath;
	hamsMeshPath <<= bem->hamsMeshPath;
	volWarning <<= bem->volWarning;
	volError <<= bem->volError;
	
	arrayShown.GetCtrl(Main::TAB_MESH,  0)->SetData(showTabMesh);
	arrayShown.GetCtrl(Main::TAB_NEMOH, 0)->SetData(showTabNemoh);
	arrayShown.GetCtrl(Main::TAB_COEFF, 0)->SetData(showTabCoeff);
	arrayShown.GetCtrl(Main::TAB_MOOR,  0)->SetData(showTabMoor);
	arrayShown.GetCtrl(Main::TAB_DECAY, 0)->SetData(showTabDecay);
	arrayShown.GetCtrl(Main::TAB_FAST,  0)->SetData(showTabFAST);
	
	dofType.SetIndex(bem->dofType);
	headingType.SetIndex(bem->headingType);
}

void MenuOptions::OnSave() {
	bem->g = ~g;
	bem->rho = ~rho;
	bem->len = ~len;
	bem->depth = ~depth;
	//bem->discardNegDOF = ~discardNegDOF;
	//bem->thres = ~thres;
	bem->calcAinf = ~calcAinf;
	bem->calcAinf_w = ~calcAinf_w;
	bem->maxTimeA = ~maxTimeA;
	bem->numValsA = ~numValsA;	
	bem->onlyDiagonal = ~onlyDiagonal;
	bem->nemohPathPreprocessor = ~nemohPathPreprocessor;
	bem->nemohPathSolver = ~nemohPathSolver;
	bem->nemohPathPostprocessor = ~nemohPathPostprocessor;	
	bem->nemohPathNew = ~nemohPathNew;
	bem->nemohPathGREN = ~nemohPathGREN;
	bem->foammPath = ~foammPath;
	bem->hamsPath = ~hamsPath;
	bem->hamsMeshPath = ~hamsMeshPath;
	bem->volWarning = ~volWarning;
	bem->volError = ~volError;
	
	showTabMesh  = arrayShown.GetCtrl(Main::TAB_MESH,  0)->GetData();
	showTabNemoh = arrayShown.GetCtrl(Main::TAB_NEMOH, 0)->GetData();
	showTabCoeff = arrayShown.GetCtrl(Main::TAB_COEFF, 0)->GetData();
	showTabMoor  = arrayShown.GetCtrl(Main::TAB_MOOR,  0)->GetData();
	showTabDecay = arrayShown.GetCtrl(Main::TAB_DECAY, 0)->GetData();
	showTabFAST  = arrayShown.GetCtrl(Main::TAB_FAST,  0)->GetData();
	
	bem->dofType = BEM::DOFType(dofType.GetIndex());
	bem->headingType = BEM::HeadingType(headingType.GetIndex());
	bem->UpdateHeadAll();
	ma().OptionsUpdated(rho, g, bem->dofType, bem->headingType);
	
	ma().SetLastTab();
}

bool MenuOptions::IsChanged() {
	if (!EqualDecimals(bem->g, ~g, 8)) 
		return true;
	if (!EqualDecimals(bem->rho, ~rho, 8))
		return true;
	if (!EqualDecimals(bem->len, ~len, 8))
		return true;
	if (!EqualDecimals(bem->depth, ~depth, 8))
		return true;
	//if (bem->discardNegDOF != ~discardNegDOF)
	//	return true;
	//if (!EqualDecimals(bem->thres, ~thres, 8)) 
	//	return true;
	if (bem->calcAinf != ~calcAinf)
		return true;
	if (bem->calcAinf_w != ~calcAinf_w)
		return true;
	if (!EqualDecimals(bem->maxTimeA, ~maxTimeA, 8))
		return true;
	if (bem->numValsA != ~numValsA)
		return true;
	if (bem->onlyDiagonal != ~onlyDiagonal)
		return true;
	if (bem->nemohPathPreprocessor != ~nemohPathPreprocessor)
		return true;
	if (bem->nemohPathSolver != ~nemohPathSolver)
		return true;
	if (bem->nemohPathPostprocessor != ~nemohPathPostprocessor)
		return true;
	if (bem->nemohPathNew != ~nemohPathNew)
		return true;
	if (bem->nemohPathGREN != ~nemohPathGREN)
		return true;
	if (bem->foammPath != ~foammPath)
		return true;
	if (bem->hamsPath != ~hamsPath)
		return true;
	if (bem->hamsMeshPath != ~hamsMeshPath)
		return true;
	if (bem->volWarning != ~volWarning)
		return true;
	if (bem->volError != ~volError)
		return true;
	if (bem->dofType != dofType.GetIndex())
		return true;
	if (bem->headingType != headingType.GetIndex())
		return true;
	
	if (showTabMesh  != arrayShown.GetCtrl(Main::TAB_MESH,  0)->GetData())
		return true;
	if (showTabNemoh != arrayShown.GetCtrl(Main::TAB_NEMOH, 0)->GetData())
		return true;
	if (showTabCoeff != arrayShown.GetCtrl(Main::TAB_COEFF, 0)->GetData())
		return true;
	if (showTabMoor  != arrayShown.GetCtrl(Main::TAB_MOOR,  0)->GetData())
		return true;
	if (showTabDecay != arrayShown.GetCtrl(Main::TAB_DECAY, 0)->GetData())
		return true;
	if (showTabFAST  != arrayShown.GetCtrl(Main::TAB_FAST,  0)->GetData())
		return true;

	return false;
}

void MenuAbout::Init() {
	CtrlLayout(*this);
	
	String qtf = GetTopic(S("BEMRosetta/main/About$en-us")); 
	SetBuildInfo(qtf);
	qtf.Replace("SYSTEMINFO", DeQtf(GetSystemInfo()));
	info.SetQTF(qtf);
}

			
Main &ma(Main *m) {
	static Main *mp = 0;
	if (m)
		mp = m;
	if (!mp)
		throw Exc(t_("Main is not initialized"));	
	return *mp;
}

BEM &Bem()							{return ma().bem;}
void Status(String str, int time)	{ma().Status(str, time);}

void OnPanic(const char *title, const char *text) {
	throw Exc(Format(t_("Error type 1 %s: %s"), title, text));	
}

void OnAssert(const char *text) {
	throw Exc(Format(t_("Error type 2: %s"), text));	
}


GUI_APP_MAIN {
	InstallPanicMessageBox(OnPanic);
	//SetAssertFailedHook(OnAssert);
	
	const Upp::Vector<String>& command = CommandLine();
	
	if (!command.IsEmpty()) {
		ConsoleOutput con(true);
		
		ConsoleMain(command, true, PrintStatus);
		return;
	}
	
	Ctrl::SetAppName(t_("Hydrodynamic coefficients viewer and converter"));
	Ctrl::GlobalBackPaint();
	
	String errorStr;
	Main *main = new Main();
	try {
		main->Init();
		main->OpenMain();
		
		Ctrl::EventLoop(main);
		
		main->Close();		
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
	delete(main);

	if (!errorStr.IsEmpty())
		Exclamation(t_("Internal error:") + S("&") + DeQtf(errorStr) + S("&") + t_("Program ended"));
}

String ForceExtSafe(String fileName, String ext) {
	if (fileName.IsEmpty())
		return String();
	return ForceExt(fileName, ext);
}

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

ArrayCtrl &ArrayModel_Init(ArrayCtrl &array, bool option) {
	array.NoHeader().NoVertGrid().AutoHideSb();
	array.HeaderObject().HideTab(array.AddColumn().HeaderTab().GetIndex());
	array.AddColumn("", 5).SetDisplay(Single<RectDisplay>());
	if (option)
		array.AddColumn("", 10);
	else
		array.AddColumn("", 0);
	array.AddColumn("", 30);	
	array.AddColumn("", 30);
	array.AddColumn("", 30);
	array.HeaderTab(1).SetMargin(-2);
	if (option)
		array.HeaderTab(2).SetMargin(-2);
	return array;
}

void ArrayModel_Add(ArrayCtrl &array, String codeStr, String title, String fileName, int id, 
					Upp::Array<Option> &option, Function <void()>OnPush) {
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

int ArrayModel_IdMesh(const ArrayCtrl &array) {
	int id;
	if (array.GetCount() == 1)
		id = ArrayModel_Id(array, 0);
	else {
		id = ArrayModel_Id(array);
		if (id < 0)
			return -1;
	}
	return Bem().GetMeshId(id);
}

int ArrayModel_IdMesh(const ArrayCtrl &array, int row) {
	int id = ArrayModel_Id(array, row);
	if (id < 0)
		return -1;
	return Bem().GetMeshId(id);
}

int ArrayModel_IdHydro(const ArrayCtrl &array) {
	int id = ArrayModel_Id(array);
	if (id < 0)
		return -1;
	return Bem().GetHydroId(id);
}

int ArrayModel_IdHydro(const ArrayCtrl &array, int row) {
	int id = ArrayModel_Id(array, row);
	if (id < 0)
		return -1;
	return Bem().GetHydroId(id);
}

Upp::Vector<int> ArrayModel_IdsHydro(const ArrayCtrl &array) {		
	Upp::Vector<int> ids;
	for (int row = 0; row < array.GetCount(); ++row) {
		int id = ArrayModel_IdHydro(array, row);		
		if (id >= 0)
			ids << id;
	}
	return ids;
}

Upp::Vector<int> ArrayModel_IdsMesh(const ArrayCtrl &array) {		
	Upp::Vector<int> ids;
	for (int row = 0; row < array.GetCount(); ++row) 
		ids << ArrayModel_IdMesh(array, row);		
	return ids;
}

void ArrayModel_IdsHydroDel(ArrayCtrl &array, const Upp::Vector<int> &ids) {		
	for (int row = array.GetCount() - 1; row >= 0 ; --row) {
		int idrow = ArrayModel_IdHydro(array, row);
		for (auto id : ids) {
			if (idrow == id) {
				array.Remove(row);
				break;
			}
		}
	}
}

void ArrayModel_RowsHydroDel(ArrayCtrl &array, const Upp::Vector<int> &rows) {		
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
