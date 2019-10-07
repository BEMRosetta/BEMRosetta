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

#include <BEMRosetta/BEMRosetta_cl/BEMRosetta.h>

#include "main.h"

void Main::Init() {
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
	if (!bem.ClearTempFiles()) 
		Cout() << "\n" << t_("BEM temporary files folder cannot be created");
	if (!LoadSerializeJson()) {
		firstTime = true;
		Cout() << "\n" << t_("Configuration data is not loaded. Defaults are set");
	}
	
	mainMesh.Init();
	mainNemoh.Init(bem);
	mainBEM.Init();
	mainOutput.Init();
	menuOptions.Init(bem);
	menuOptions.Load();
	menuAbout.Init();
	if (bem.experimental) {
		tab.Add(mainMesh.SizePos(),   t_("Mesh"));
		tab.Add(mainNemoh.SizePos(),  t_("Nemoh"));
	}
	tab.Add(mainBEM.SizePos(),    t_("Coefficients"));
	tab.Add().Disable();
	tab.Add(mainOutput.SizePos(), t_("Output"));
	tab.Add(menuOptions.SizePos(),t_("Options"));
	tab.Add().Disable();
	tab.Add(menuAbout.SizePos(),  t_("About"));
	
	Add(tab.SizePos());	
		
	tab.WhenSet = [&] {
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
	};	
	
	if (firstTime)
		tab.Set(menuAbout);
	
	AddFrame(bar);
	
	BEMData::Print 	  	  = [this](String s) {printf("%s", ~s); mainOutput.Print(s);};
	BEMData::PrintWarning = [this](String s) {printf("%s", ~s); mainOutput.Print(s);};
	BEMData::PrintError   = [this](String s) {printf("%s", ~s); mainOutput.Print(s); tab.Set(mainOutput);};
}

void Main::OptionsUpdated() {
	mainBEM.OnOpt();
	mainMesh.OnOpt();
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

Main::~Main() {
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
	experimental <<= bem->experimental;
	experimentalFOAMM <<= bem->experimentalFOAMM;
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
	bem->experimental = ~experimental;	
	bem->experimentalFOAMM = ~experimentalFOAMM;
	bem->foammPath = ~foammPath;
	
	ma().OptionsUpdated();
}

bool MenuOptions::IsChanged() {
	if (TruncDecimals(bem->g, 8) !=  TruncDecimals(static_cast<double>(~g), 8)) 
		return true;
	if (TruncDecimals(bem->rho, 8) !=  TruncDecimals(static_cast<double>(~rho), 8))
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
	if (bem->experimental != ~experimental)
		return true;
	if (bem->experimentalFOAMM != ~experimentalFOAMM)
		return true;
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

void MainStiffness::AddPrepare(int &row0, int &col0, String name, int icase, String bodyName, int ibody) {
	row0 = ibody*9 + 1;
	col0 = icase*8;
	
	array.Set(0, col0, AttrText(name).Bold());
	
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

void MainStiffness::Add(String name, int icase, const MatrixXd &K) {
	int row0, col0;
	AddPrepare(row0, col0, name, icase, "", 0);
	
	if (K.size() == 0)
		return;
	for (int r = 0; r < 6; ++r) {
		for (int c = 0; c < 6; ++c)
			array.Set(row0 + r + 2, col0 + c + 1, AttrText(FormatDouble(K(r, c), 7, FD_EXP|FD_CAP_E|FD_REL)).Align(ALIGN_RIGHT));
	}	
}
	
void MainStiffness::Add(String name, int icase, String bodyName, int ibody, const Hydro &hydro) {
	int row0, col0;
	AddPrepare(row0, col0, name, icase, bodyName, ibody);

	if (hydro.C.IsEmpty())
		return;

	for (int r = 0; r < 6; ++r) {
		for (int c = 0; c < 6; ++c)
			array.Set(row0 + r + 2, col0 + c + 1, AttrText(FormatDouble(hydro.C_dim(ibody, r, c), 7, FD_EXP|FD_CAP_E|FD_REL)).Align(ALIGN_RIGHT));
	}
}

void MainStiffness::Load(Upp::Array<HydroClass> &hydros) {
	Clear();
	
	for (int icase = 0; icase < hydros.GetCount(); ++icase) {
		Hydro &hydro = hydros[icase].hd();
		for (int ibody = 0; ibody < hydro.Nb; ++ibody) 
			Add(hydros[icase].hd().name, icase, hydros[icase].hd().names[ibody], ibody, hydro);
	}
}	

void MainStiffness::Load(Upp::Array<MeshData> &surfs) {
	Clear();

	for (int icase = 0; icase < surfs.GetCount(); ++icase) 
		Add(GetFileTitle(surfs[icase].file), icase, surfs[icase].c);
}
			
Main &ma(Main *m) {
	static Main *mp = 0;
	if (m)
		mp = m;
	return *mp;
}

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

