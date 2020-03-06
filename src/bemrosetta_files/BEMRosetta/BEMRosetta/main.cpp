#include <CtrlLib/CtrlLib.h>

using namespace Upp;

#define IMAGECLASS Img
#define IMAGEFILE <BEMRosetta/BEMRosetta/main.iml>
#include <Draw/iml.h>

#define TOPICFILE <BEMRosetta/BEMRosetta/main.tpp/all.i>
#include <Core/topic_group.h>

#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <CtrlScroll/CtrlScroll.h>

#include <BEMRosetta/BEMRosetta_cl/BEMRosetta.h>

#include "main.h"

using namespace Eigen;

FreqSelector::FreqSelector() {
	add.SetImage(Img::add());
	add.WhenAction = [&] {AddField();};
	Add(add.LeftPos(0, 19).TopPos(3, 19));
};
	
void Main::Init() {
	LOG("Init");
	Title("BEMRosetta");
	Sizeable().Zoomable().SetMinSize(Size(800, 600));
	Icon(Img::Rosetta64());
	LargeIcon(Img::Rosetta256());
	ma(this);

	bool firstTime = false;
	if (!bem.LoadSerializeJson()) {
		firstTime = true;
		Cout() << "\n" << t_("BEM configuration data is not loaded. Defaults are set");
	}
	LOG("BEM configuration loaded");
	
	if (IsNull(lastTab))
		lastTab = 0;
	
	if (!bem.ClearTempFiles()) 
		Cout() << "\n" << t_("BEM temporary files folder cannot be created");
	if (!LoadSerializeJson()) {
		firstTime = true;
		Cout() << "\n" << t_("Configuration data is not loaded. Defaults are set");
	}
	LOG("Configuration loaded");
	
	mainMesh.Init();			LOG("Init Mesh");
	mainNemoh.Init(bem);		LOG("Init Nemoh");
	mainBEM.Init();				LOG("Init BEM");
	mainOutput.Init();			LOG("Init Output");
	menuOptions.Init(bem);		LOG("Init Options");
	menuOptions.Load();			LOG("Init Options.Load");
	menuAbout.Init();			LOG("Init About");
	
	tab.Add(mainMesh.SizePos(),   t_("Mesh"));
	tab.Add(mainNemoh.SizePos(),  t_("Nemoh"));
	tab.Add(mainBEM.SizePos(),    t_("Coefficients"));
	tab.Add().Disable();
	tab.Add(mainOutput.SizePos(), t_("Output"));
	tab.Add(menuOptions.SizePos(),t_("Options"));
	tab.Add().Disable();
	tab.Add(menuAbout.SizePos(),  t_("About"));
	
	Add(tab.SizePos());	

	Add(labrho.SetLabel(t_("rho [Kg/m3]:")).RightPosZ(58, 80).TopPosZ(0, 19));
	Add(editrho.SetReadOnly().RightPosZ(28, 40).TopPosZ(1, 19));
	editrho <<= bem.rho;
	
	Add(labg.SetLabel(t_("Gravity [m/s2]:")).RightPosZ(182, 80).TopPosZ(0, 19));
	Add(editg.SetReadOnly().RightPosZ(148, 40).TopPosZ(1, 19));
	editg <<= bem.g;
	
	butWindow.SetImage(Img::application_double());
	Add(butWindow.RightPosZ(0, 20).TopPosZ(0, 20));
	butWindow.Hide();
	
	butWindow.WhenAction = [&] {
		if (tab.IsAt(mainMesh)) {
			MainMeshW *mainMeshW = new MainMeshW();
			mainMeshW->Init(mainMesh);
			mainMeshW->OpenMain();
		} else if (tab.IsAt(mainBEM)) {
			MainBEMW *mainBEMW = new MainBEMW();
			mainBEMW->Init(mainBEM);
			mainBEMW->OpenMain();
		}
	};
	
	tab.WhenSet = [&] {
		if (tab.Get() < 0)
			return;
		
		LOGTAB(tab);
		if (tab.IsAt(menuOptions)) 
			menuOptions.Load();
		else if (tab.IsAt(mainNemoh)) 
			mainNemoh.Load(bem);
		if (!tab.IsAt(menuOptions) && menuOptions.IsChanged()) {
			if (PromptYesNo(t_("Options have changed&Do you want to save them?")))
				menuOptions.OnSave();
			else
				menuOptions.Load();
		}
		
		if (tab.IsAt(mainMesh) || tab.IsAt(mainBEM)) {
			butWindow.Show(true);
			lastTab = ~tab;
		} else 	if (tab.IsAt(mainNemoh)) {
			butWindow.Show(false);
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
	
	tab.Set(-1);
	if (firstTime) {
		lastTab = 0;
		tab.Set(menuAbout);
	} else
		tab.Set(lastTab);
	
	AddFrame(bar);
	
	BEMData::Print 	  	  = [this](String s) {printf("%s", ~s); mainOutput.Print(s);};
	BEMData::PrintWarning = [this](String s) {printf("%s", ~s); mainOutput.Print(s);};
	BEMData::PrintError   = [this](String s) {printf("%s", ~s); mainOutput.Print(s); tab.Set(mainOutput);};
}

void Main::OptionsUpdated(double rho, double g) {
	mainBEM.OnOpt();
	mainMesh.OnOpt();
	
	editg <<= g;
	editrho <<= rho;
}

bool Main::LoadSerializeJson() {
	bool ret;
	String folder = AppendFileName(GetAppDataFolder(), "BEMRosetta");
	DirectoryCreate(folder);
	if (!DirectoryExists(folder))
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
	mainNemoh.InitSerialize(ret);
	mainBEM.InitSerialize(ret);
	
	if (!ret)
		tab.Set(menuAbout);
	
	return ret;
}

bool Main::StoreSerializeJson() {
	String folder = AppendFileName(GetAppDataFolder(), "BEMRosetta");
	DirectoryCreate(folder);
	if (!DirectoryExists(folder))
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
		("nemoh", mainNemoh)
		("bem", mainBEM)
		("lastTab", lastTab)
	;
}

void MenuOptions::Init(BEMData &_bem) {
	CtrlLayout(*this);
	
	bem = &_bem;
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
	onlyDiagonal <<= bem->onlyDiagonal;
	nemohPathPreprocessor <<= bem->nemohPathPreprocessor;
	nemohPathSolver <<= bem->nemohPathSolver;
	nemohPathPostprocessor <<= bem->nemohPathPostprocessor;
	nemohPathNew <<= bem->nemohPathNew;
	nemohPathGREN <<= bem->nemohPathGREN;
	//experimental <<= bem->experimental;
	foammPath <<= bem->foammPath;
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
	bem->onlyDiagonal = ~onlyDiagonal;
	bem->nemohPathPreprocessor = ~nemohPathPreprocessor;
	bem->nemohPathSolver = ~nemohPathSolver;
	bem->nemohPathPostprocessor = ~nemohPathPostprocessor;	
	bem->nemohPathNew = ~nemohPathNew;
	bem->nemohPathGREN = ~nemohPathGREN;
	//bem->experimental = ~experimental;	
	bem->foammPath = ~foammPath;
	
	ma().OptionsUpdated(rho, g);
}

bool MenuOptions::IsChanged() {
	if (TruncDecimals(bem->g, 8) !=  TruncDecimals(double(~g), 8)) 
		return true;
	if (TruncDecimals(bem->rho, 8) !=  TruncDecimals(double(~rho), 8))
		return true;
	if (bem->length != ~length)
		return true;
	if (bem->depth != ~depth)
		return true;
	if (bem->discardNegDOF != ~discardNegDOF)
		return true;
	if (bem->thres != ~thres) 
		return true;
	if (bem->calcAwinf != ~calcAwinf)
		return true;
	if (bem->maxTimeA != ~maxTimeA)
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
	//if (bem->experimental != ~experimental)
	//	return true;
	if (bem->foammPath != ~foammPath)
		return true;
	
	return false;
}

void MenuAbout::Init() {
	CtrlLayout(*this);
	
	String qtf = GetTopic(S("BEMRosetta/BEMRosetta/main/About$en-us")); 
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
		array.CreateCtrl<Button>(row0, col0+5, false).SetLabel(t_("Save")).Tip(t_("Save to Wamit .hst stiffness matrix format"))
			.WhenAction = [&] {
				FileSel fs;
				fs.Type(t_("Wamit stiffness matrix format"), "*.hst");
				if (fs.ExecuteSaveAs(t_("Save to Wamit .hst stiffness matrix format"))) 
					mesh.SaveHST(~fs, Bem().rho, Bem().g);
			};
	}
	if (button && Bem().hydros.GetCount() > 0) {
		array.CreateCtrl<Button>(row0, col0+6, false).SetLabel(t_("Copy")).Tip(t_("Copy matrix and paste it in selected BEM Coefficients file and body"))
			.WhenAction = [=] {
				WithBEMList<TopWindow> w;
				CtrlLayout(w);
				w.Title(t_("Copy matrix and paste it in selected BEM Coefficients file and body"));
				w.array.SetLineCy(EditField::GetStdHeight());
				w.array.AddColumn(t_("File"), 20);
				w.array.AddColumn(t_("Body"), 10);
				w.array.HeaderObject().HideTab(w.array.AddColumn().HeaderTab().GetIndex());
				w.array.HeaderObject().HideTab(w.array.AddColumn().HeaderTab().GetIndex());
				for (int f = 0; f < Bem().hydros.GetCount(); ++f) {
					const Hydro &hy = Bem().hydros[f].hd();
					for (int ib = 0; ib < hy.Nb; ++ib)
						w.array.Add(hy.name, hy.names[ib].IsEmpty() ? AsString(ib+1) : hy.names[ib], f, ib);
				}
				bool cancel = true;
				w.butSelect.WhenAction = [&] {cancel = false;	w.Close();};
				w.butCancel.WhenAction = [&] {w.Close();}; 
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

bool MainStiffness::Load(Upp::Array<HydroClass> &hydros, const Vector<int> &ids) {
	Clear();
	
	for (int i = 0; i < ids.GetCount(); ++i) {
		int isurf = ids[i];
		Hydro &hydro = hydros[isurf].hd();
		for (int ibody = 0; ibody < hydro.Nb; ++ibody) 
			Add(hydros[isurf].hd().name, i, hydros[isurf].hd().names[ibody], ibody, hydro, hydro.GetId());
	}
	return true;
}	

void MainStiffness::Load(Upp::Array<MeshData> &surfs, const Vector<int> &ids) {
	Clear();

	for (int i = 0; i < ids.GetCount(); ++i) {
		int isurf = ids[i];	
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
	
	ConsoleOutput console;
	
	const Vector<String>& command = CommandLine();
	
	if (!command.IsEmpty()) {
		ConsoleMain(command, true);
		return;
	}
	
	Ctrl::SetAppName(t_("Hydrodynamic coefficents viewer and converter"));
	Ctrl::GlobalBackPaint();
	
	String errorStr;
	try {
		Main main;
		
		main.Init();
		main.Run();
		
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
	opt.WhenAction = OnPush;
}

void ArrayModel_Add(ArrayCtrl &array, String codeStr, String title, String fileName, int id) {
 	array.Add(id, GetColorId(id), codeStr, title, fileName);
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
	int id = ArrayModel_Id(array);
	if (id < 0)
		return -1;
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

Vector<int> ArrayModel_IdsHydro(const ArrayCtrl &array) {		
	Vector<int> ids;
	for (int row = 0; row < array.GetCount(); ++row) 
		ids << ArrayModel_IdHydro(array, row);		
	return ids;
}

Vector<int> ArrayModel_IdsMesh(const ArrayCtrl &array) {		
	Vector<int> ids;
	for (int row = 0; row < array.GetCount(); ++row) 
		ids << ArrayModel_IdMesh(array, row);		
	return ids;
}

void ArrayModel_IdsHydroDel(ArrayCtrl &array, const Vector<int> &ids) {		
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

void ArrayModel_RowsHydroDel(ArrayCtrl &array, const Vector<int> &rows) {		
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

const Color& ArrayModel_GetColor(const ArrayCtrl &array, int row) {
	return GetColorId(array.Get(row, 0));
}

String ArrayModel_GetFileName(ArrayCtrl &array, int row) {
	if (row < 0) 
		row = array.GetCursor();
	if (row < 0)
		return String();
	if (array.GetCtrl(row, 2))
		return array.Get(row, 5);
 	return array.Get(row, 4);
}