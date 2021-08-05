#include <CtrlLib/CtrlLib.h>

using namespace Upp;

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

#include <BEMRosetta_cl/BEMRosetta.h>

#include "main.h"


using namespace Eigen;

FreqSelector::FreqSelector() {
	add.SetImage(Img::add());
	add << [&] {AddField();};
	Add(add.LeftPos(0, 19).TopPos(3, 19));
};

static String GetBEMRosettaDataFolder() {
	return AppendFileName(GetAppDataFolder(), "BEMRosetta");
}

void Main::Init() {
	LOG("Init");
	Sizeable().Zoomable().SetMinSize(Size(800, 600));
	Icon(Img::Rosetta64());
	LargeIcon(Img::Rosetta256());
	ma(this);
	
	Title(S("BEMRosetta") + (Bem().experimental ? " EXPERIMENTAL" : ""));

	tabTexts << t_("Mesh Handling") << t_("BEM Solver") << t_("Hydrodynamic Coefficients") 
			 << t_("Mooring") << t_("Decay") << t_("FAST .out Reader");
		
	bool firstTime = false, openOptions = false;
	if (!bem.LoadSerializeJson(firstTime)) 
		Cout() << "\n" << t_("BEM configuration data is not loaded. Defaults are set");
	
	LOG("BEM configuration loaded");
	
	if (IsNull(lastTab))
		lastTab = 0;
	
	if (!bem.ClearTempFiles()) 
		Cout() << "\n" << t_("BEM temporary files folder cannot be created");
	
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
		if (Bem().experimental)
			tab.Add(mainMoor.SizePos(), tabTexts[TAB_MOOR]);
	}
	if (menuOptions.showTabDecay) {
		if (!menuOptions.showTabMoor)
			tab.Add().Disable();
		mainDecay.Init();	LOG("Init Decay");
		if (Bem().experimental)
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
	tab.Add(menuAbout.SizePos(),  t_("About"));
	
	Add(tab.SizePos());	

	Add(labrho.SetLabel(t_("rho [kg/m3]:")).RightPosZ(128, 80).TopPosZ(1, 22));
	Add(editrho.SetReadOnly().RightPosZ(98, 40).TopPosZ(2, 20));
	editrho <<= bem.rho;
	
	Add(labg.SetLabel(t_("Gravity [m/s2]:")).RightPosZ(252, 80).TopPosZ(1, 22));
	Add(editg.SetReadOnly().RightPosZ(218, 40).TopPosZ(2, 20));
	editg <<= bem.g;
	
	butWindow.SetImage(Img::application_double()).SetLabel(t_("New window")).Tip(t_("Open new window"));
	Add(butWindow.RightPosZ(2, 90).TopPosZ(0, 22));
	butWindow.Hide();
	
	butWindow << [&] {
		if (tab.IsAt(mainMesh)) {
			MainMeshW *mainMeshW = new MainMeshW();
			mainMeshW->Init(mainMesh, Img::Rosetta64(), Img::Rosetta256());
			mainMeshW->OpenMain();
		} else if (tab.IsAt(mainBEM)) {
			MainBEMW *mainBEMW = new MainBEMW();
			mainBEMW->Init(mainBEM, Img::Rosetta64(), Img::Rosetta256());
			mainBEMW->OpenMain();
		} else if (tab.IsAt(mainFAST)) {
			MainFASTW *mainFASTW = new MainFASTW();
			mainFASTW->Init(GetBEMRosettaDataFolder(), Img::Rosetta64(), Img::Rosetta256(), bar);
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
		} else {
			butWindow.Show(false);
			lastTab = 0;
		}
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
		tab.Set(lastTab);
	
	tab.WhenSet();
	
	AddFrame(bar);
	
	BEMData::Print 	  	  = [this](String s) {printf("%s", ~s); mainOutput.Print(s);};
	BEMData::PrintWarning = [this](String s) {printf("%s", ~s); mainOutput.Print(s); Status(s);};
	BEMData::PrintError   = [this](String s) {printf("%s", ~s); mainOutput.Print(s); tab.Set(mainOutput); Status(s);};
}

void Main::OptionsUpdated(double rho, double g) {
	mainBEM.OnOpt();
	mainMesh.OnOpt();
	
	editg <<= g;
	editrho <<= rho;
}

bool Main::LoadSerializeJson(bool &firstTime, bool &openOptions) {
	bool ret;
	String folder = GetBEMRosettaDataFolder();
	if (!DirectoryCreateX(folder))
		ret = false;
	else {
		String fileName = AppendFileName(folder, "config.cf");
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
	String fileName = AppendFileName(folder, "config.cf");
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

void MenuOptions::Init(BEMData &_bem) {
	CtrlLayout(*this);
	
	bem = &_bem;
	butSave  <<= THISBACK(OnSave);
	butSave2 <<= THISBACK(OnSave);
	
	arrayShown.AddColumn("");
	arrayShown.Add();	arrayShown.CreateCtrl<Option>(Main::TAB_MESH,  0, false).SetLabel(ma().tabTexts[Main::TAB_MESH]);
	arrayShown.Add();	arrayShown.CreateCtrl<Option>(Main::TAB_NEMOH, 0, false).SetLabel(ma().tabTexts[Main::TAB_NEMOH]);
	arrayShown.Add();	arrayShown.CreateCtrl<Option>(Main::TAB_COEFF, 0, false).SetLabel(ma().tabTexts[Main::TAB_COEFF]);
	arrayShown.Add();	arrayShown.CreateCtrl<Option>(Main::TAB_MOOR,  0, false).SetLabel(ma().tabTexts[Main::TAB_MOOR]).Show(Bem().experimental);
	arrayShown.Add();	arrayShown.CreateCtrl<Option>(Main::TAB_DECAY, 0, false).SetLabel(ma().tabTexts[Main::TAB_DECAY]).Show(Bem().experimental);
	arrayShown.Add();	arrayShown.CreateCtrl<Option>(Main::TAB_FAST,  0, false).SetLabel(ma().tabTexts[Main::TAB_FAST]);
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
	calcAwinf <<= bem->calcAwinf;
	calcAwinfw <<= bem->calcAwinfw;
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
	
	arrayShown.GetCtrl(Main::TAB_MESH,  0)->SetData(showTabMesh);
	arrayShown.GetCtrl(Main::TAB_NEMOH, 0)->SetData(showTabNemoh);
	arrayShown.GetCtrl(Main::TAB_COEFF, 0)->SetData(showTabCoeff);
	arrayShown.GetCtrl(Main::TAB_MOOR,  0)->SetData(showTabMoor);
	arrayShown.GetCtrl(Main::TAB_DECAY, 0)->SetData(showTabDecay);
	arrayShown.GetCtrl(Main::TAB_FAST,  0)->SetData(showTabFAST);
}

void MenuOptions::OnSave() {
	bem->g = ~g;
	bem->rho = ~rho;
	bem->len = ~len;
	bem->depth = ~depth;
	//bem->discardNegDOF = ~discardNegDOF;
	//bem->thres = ~thres;
	bem->calcAwinf = ~calcAwinf;
	bem->calcAwinfw = ~calcAwinfw;
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
	
	showTabMesh  = arrayShown.GetCtrl(Main::TAB_MESH,  0)->GetData();
	showTabNemoh = arrayShown.GetCtrl(Main::TAB_NEMOH, 0)->GetData();
	showTabCoeff = arrayShown.GetCtrl(Main::TAB_COEFF, 0)->GetData();
	showTabMoor  = arrayShown.GetCtrl(Main::TAB_MOOR,  0)->GetData();
	showTabDecay = arrayShown.GetCtrl(Main::TAB_DECAY, 0)->GetData();
	showTabFAST  = arrayShown.GetCtrl(Main::TAB_FAST,  0)->GetData();
	
	ma().OptionsUpdated(rho, g);
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
	if (bem->calcAwinf != ~calcAwinf)
		return true;
	if (bem->calcAwinfw != ~calcAwinfw)
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
	info.SetQTF(qtf);
}


void MainStiffness::Init() {
	CtrlLayout(*this);
}

void MainStiffness::Clear() {
	array.Reset();
	array.NoHeader().SetLineCy(EditField::GetStdHeight()).HeaderObject().Absolute();
	array.MultiSelect().SpanWideCells();
	array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, array);};
}

void MainStiffness::AddPrepare(int &row0, int &col0, String name, int icase, String bodyName, int ibody, int idc) {
	row0 = ibody*9 + 1;
	col0 = icase*8;
	
	array.Set(0, col0, AttrText(name).Bold());
	array.Set(2, col0, AttrText(" ").Paper(GetColorId(idc)));
	
	while (array.GetColumnCount() < col0 + 7) {
		if (icase > 0) {
			array.AddColumn("", 10);
			int icol = array.GetColumnCount() - 1;
			array.HeaderObject().Tab(icol).SetMargin(0);
		}
		array.AddColumn("", 20);
		for (int i = 0; i < 2; ++i)
			array.AddColumn("", 20);	
		for (int i = 0; i < 4; ++i)
			array.AddColumn("", 70);
	}
	array.Set(row0, col0, AttrText(Format(t_("#%d body. %s"), ibody + 1, bodyName)).Bold());
	for (int i = 0; i < 6; ++i) {
		array.Set(row0 + 1, 	col0 + i + 1, AttrText(FormatInt(i + 1)).Bold().Align(ALIGN_CENTER));
		array.Set(row0 + i + 2, col0, 	  	  AttrText(FormatInt(i + 1)).Bold().Align(ALIGN_CENTER));
	}
	for (int i = 1; i <= ibody; ++i) 
		array.SetLineCy(i*9, 10);	
}

void MainStiffness::Add(const MeshData &mesh, int icase, bool button) {
	String name = mesh.fileName;
	const MatrixXd &K = mesh.C;
	int idc = mesh.GetId();
	
	int row0, col0;
	AddPrepare(row0, col0, name, icase, "", 0, idc);
	
	if (K.size() == 0)
		return;
	for (int r = 0; r < 6; ++r) {
		for (int c = 0; c < 6; ++c)
			array.Set(row0 + r + 2, col0 + c + 1, AttrText(FormatDouble(K(r, c), 7, FD_EXP|FD_CAP_E|FD_REL)).Align(ALIGN_RIGHT));
	}
	if (button) {
		array.CreateCtrl<Button>(row0, col0+5, false).SetLabel(t_("Save")).Tip(t_("Saves to Wamit .hst stiffness matrix format"))
			<< [&] {
				FileSel fs;
				fs.Type(t_("Wamit stiffness matrix format"), "*.hst");
				if (fs.ExecuteSaveAs(t_("Save to Wamit .hst stiffness matrix format"))) 
					mesh.SaveHST(~fs, Bem().rho, Bem().g);
			};
	}
	if (button && Bem().hydros.size() > 0) {
		array.CreateCtrl<Button>(row0, col0+6, false).SetLabel(t_("Copy")).Tip(t_("Copies matrix and paste it in selected BEM Coefficients file and body"))
			<< [=] {
				WithBEMList<TopWindow> w;
				CtrlLayout(w);
				w.Title(t_("Copies matrix and paste it in selected BEM Coefficients file and body"));
				w.array.SetLineCy(EditField::GetStdHeight());
				w.array.AddColumn(t_("File"), 20);
				w.array.AddColumn(t_("Body"), 10);
				w.array.HeaderObject().HideTab(w.array.AddColumn().HeaderTab().GetIndex());
				w.array.HeaderObject().HideTab(w.array.AddColumn().HeaderTab().GetIndex());
				for (int f = 0; f < Bem().hydros.size(); ++f) {
					const Hydro &hy = Bem().hydros[f].hd();
					for (int ib = 0; ib < hy.Nb; ++ib)
						w.array.Add(hy.name, hy.names[ib].IsEmpty() ? AsString(ib+1) : hy.names[ib], f, ib);
				}
				bool cancel = true;
				w.butSelect << [&] {cancel = false;	w.Close();};
				w.butCancel << [&] {w.Close();}; 
				w.Execute();
				if (!cancel) {
					int id = w.array.GetCursor();
					if (id < 0)
						return; 
					int f = w.array.Get(id, 2);
					int ib = w.array.Get(id, 3);
					Bem().hydros[f].hd().SetC(ib, K);
				}
			 };
	}
}
	
void MainStiffness::Add(String name, int icase, String bodyName, int ibody, const Hydro &hydro, int idc) {
	int row0, col0;
	AddPrepare(row0, col0, name, icase, bodyName, ibody, idc);

	if (hydro.C.IsEmpty() || hydro.C[ibody].size() == 0)
		return;

	for (int r = 0; r < 6; ++r) {
		for (int c = 0; c < 6; ++c)
			array.Set(row0 + r + 2, col0 + c + 1, AttrText(FormatDouble(hydro.C_dim(ibody, r, c), 7, FD_EXP|FD_CAP_E|FD_REL)).Align(ALIGN_RIGHT));
	}
}

bool MainStiffness::Load(Upp::Array<HydroClass> &hydros, const Upp::Vector<int> &ids) {
	Clear();
	
	for (int i = 0; i < ids.size(); ++i) {
		int isurf = ids[i];
		Hydro &hydro = hydros[isurf].hd();
		for (int ibody = 0; ibody < hydro.Nb; ++ibody) 
			Add(hydro.name, i, hydro.names[ibody], ibody, hydro, hydro.GetId());
	}
	return true;
}	

void MainStiffness::Load(Upp::Array<MeshData> &surfs, const Upp::Vector<int> &ids) {
	Clear();

	for (int i = 0; i < ids.size(); ++i) {
		int isurf = ids[i];
		if (isurf >= 0)	
			Add(surfs[isurf], i, true);
	}
}
			
Main &ma(Main *m) {
	static Main *mp = 0;
	if (m)
		mp = m;
	if (!mp)
		throw Exc(t_("Main is not initialized"));	
	return *mp;
}

BEMData &Bem()						{return ma().bem;}
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
		
		ConsoleMain(command, true);
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
	                   Color ink, Color paper, dword style) const {
		w.DrawRect(r, q);
	}
};

ArrayCtrl &ArrayModel_Init(ArrayCtrl &array, bool option) {
	array.NoHeader().NoVertGrid().AutoHideSb();
	array.HeaderObject().HideTab(array.AddColumn().HeaderTab().GetIndex());
	array.AddColumn("", 5).SetDisplay(Single<RectDisplay>());
	if (option)
		array.AddColumn("", 10);
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
 	array.Add(id, GetColorId(id), codeStr, title, fileName);
 	array.SetCursor(array.GetCount());
}

void ArrayModel_Change(ArrayCtrl &array, int id, String codeStr, String title, String fileName) {
	for (int row = 0; row < array.GetCount(); ++row) {
		if (array.Get(row, 0) == id) {
			if (!IsNull(codeStr))
				array.Set(row, 2, codeStr);
			if (!IsNull(title))
				array.Set(row, 3, title);
			if (!IsNull(fileName))
				array.Set(row, 4, fileName);
			return;
		}
	}
	throw Exc(t_("Id not found in ArrayModel_Change()"));
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

String ArrayModel_GetFileName(ArrayCtrl &array, int row) {
	if (row < 0) 
		row = array.GetCursor();
	if (row < 0)
		return String();
	return array.Get(row, 5);
}

String ArrayModel_GetTitle(ArrayCtrl &array, int row) {
	if (row < 0) 
		row = array.GetCursor();
	if (row < 0)
		return String();
 	return array.Get(row, 4);
}