#ifndef _BEM_Rosetta_GUI_BEM_Rosetta_GUI_h
#define _BEM_Rosetta_GUI_BEM_Rosetta_GUI_h

#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas.h>

#include <BEMRosetta/BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;

class RichTextView2 : public RichTextView {
public:
	RichTextView2() {
		zoomlevel = 5;
		Background(Color(245, 245, 245));
	}
	virtual void Layout() {
		RichTextView::Layout();
		PageWidth(int(zoomlevel*GetSize().cx));
	}
	double zoomlevel;
};


#define LAYOUTFILE <BEMRosetta/BEMRosetta/main.lay>
#include <CtrlCore/lay.h>

#define IMAGECLASS Img
#define IMAGEFILE <BEMRosetta/BEMRosetta/main.iml>
#include <Draw/iml_header.h>

#include "arrange.h"

enum DataToShow {DATA_A, DATA_B, DATA_FORCE_SC, DATA_FORCE_FK, DATA_FORCE_EX};
enum DataToPlot {PLOT_A, PLOT_AINF, PLOT_B, PLOT_FORCE_SC_MA, PLOT_FORCE_SC_PH,
				PLOT_FORCE_FK_MA, PLOT_FORCE_FK_PH, PLOT_FORCE_EX_MA, PLOT_FORCE_EX_PH};
    
class HydroSource : public DataSource {
public:
	HydroSource() : data(0) {}
	HydroSource(Hydro &data, int i, int j_h, DataToPlot dataToPlot)	{Init(data, i, j_h, dataToPlot);}
	bool Init(Hydro *data, int i, int j_h, DataToPlot dataToPlot) 	{return Init(*data, i, j_h, dataToPlot);}
	bool Init(Hydro &data, int i, int j_h, DataToPlot dataToPlot) {
		this->data = &data;
		this->i = data.dofOrder[i];
		if (dataToPlot == PLOT_A || dataToPlot == PLOT_AINF || dataToPlot == PLOT_B)
			this->j_h = data.dofOrder[j_h];
		else
			this->j_h = j_h;
		this->dataToPlot = dataToPlot;
		if (IsNaN(y(0)))
			return false;
		return true;
	}
	virtual inline double y(int64 id) {
		switch (dataToPlot) {
		case PLOT_A:			return data->A[int(id)](i, j_h);
		case PLOT_AINF:			return data->Awinf(i, j_h);
		case PLOT_B:			return data->B[int(id)](i, j_h);
		case PLOT_FORCE_SC_MA:	return data->sc.ma[j_h](int(id), i);
		case PLOT_FORCE_SC_PH:	return data->sc.ph[j_h](int(id), i);
		case PLOT_FORCE_FK_MA:	return data->fk.ma[j_h](int(id), i);
		case PLOT_FORCE_FK_PH:	return data->fk.ph[j_h](int(id), i);
		case PLOT_FORCE_EX_MA:	return data->ex.ma[j_h](int(id), i);
		case PLOT_FORCE_EX_PH:	return data->ex.ph[j_h](int(id), i);
		default:				NEVER();	return Null;
		}
	}
	virtual inline double x(int64 id) 	{return data->w[int(id)];}
	virtual int64 GetCount()		  	{return data->Nf;}
	
private:
	Hydro *data;
	int i, j_h;
	DataToPlot dataToPlot;
};

class MenuAbout : public WithMenuAbout<StaticRect> {
public:
	typedef MenuAbout CLASSNAME;
	void Init();
	void HelpHandler(String &str);
};

class MainSummary : public WithMainSummary<StaticRect> {
public:
	typedef MainSummary CLASSNAME;
	void Init();
	void Clear();
	void Report(Hydro &data, int id);
	void OnArrayBar(Bar &menu); 
	void ArrayCopy();
	void ArraySelect();
};

class MainOutput : public WithMainOutput<StaticRect> {
public:
	typedef MainOutput CLASSNAME;
	void Init();
	void Print(String str);
};

class MainArrange : public WithMainArrange<StaticRect> {
public:
	typedef MainArrange CLASSNAME;
	void Init();
	void Clear();
	void Load(Upp::Array<HydroClass> &hydro);
	
private:
	Upp::Array<ArrangeDOF> arrangeDOF;
};

class MainView : public WithMainView<StaticRect> {
public:
	typedef MainView CLASSNAME;
	
	MainView() : maxX(Null), maxY(Null), maxZ(Null) {}
	void Init();
	
	void CalcEnvelope();
		
	VolumeEnvelope env;
	double maxX, maxY, maxZ;
	
	void ZoomToFit();
	
private:
	void OnPaint();
};


class MainPlot : public WithMainPlot<StaticRect> {
public:
	typedef MainPlot CLASSNAME;
	void Init(int i, int j_h, double h, DataToShow dataToShow);
	bool Load(Upp::Array<HydroClass> &hydro);
	
	Upp::Array<HydroSource> ABF_source, ABF_source2;
	Upp::Array<HydroSource> Ainf_source;
	int i, j_h;
	DataToShow dataToShow;
};

class MainABForce : public WithMainABForce<StaticRect> {
public:
	typedef MainABForce CLASSNAME;
	void Init(DataToShow dataToShow);
	void Clear();
	bool Load(Upp::Array<HydroClass> &hydro);
	
	Upp::Array<Upp::Array<MainPlot>> plots;

private:
	int selTab;
	bool isFilling;
	DataToShow dataToShow;
};

class Main : public WithMain<TopWindow> {
public:
	typedef Main CLASSNAME;
	
	Main() : closed(false) {}
	virtual ~Main();
	void Init();
	void Close(bool store = false);
	
	void OnLoad();
	void OnConvert();
	void OnView();
	void OnOpt();
	
	int LoadSerializeJson() {
		String folder = AppendFileName(GetAppDataFolder(), "BEMRosetta");
		DirectoryCreate(folder);
		if (!DirectoryExists(folder))
			return 0;
		String fileName = AppendFileName(folder, "config.cf");
		if (!FileExists(fileName)) 
			return 1;
		String jsonText = LoadFile(fileName);
		if (jsonText.IsEmpty())
			return 0;
		if (!LoadFromJson(*this, jsonText))
			return 0;
		return 2;
	}
	bool StoreSerializeJson() {
		String folder = AppendFileName(GetAppDataFolder(), "BEMRosetta");
		DirectoryCreate(folder);
		if (!DirectoryExists(folder))
			return 0;
		String fileName = AppendFileName(folder, "config.cf");
		return StoreAsJsonFile(*this, fileName, true);
	}
	void Jsonize(JsonIO &json);
		
	WithMenuOpen<StaticRect> menuOpen;
	WithMenuConvert<StaticRect> menuConvert;
	WithMenuPlot<StaticRect> menuPlot;
	WithMenuView<StaticRect> menuView;
	MenuAbout menuAbout;
	
	MainSummary mainSummary;
	MainArrange mainArrange;
	MainABForce mainA;
	MainABForce mainB;
	MainABForce mainForceSC, mainForceFK, mainForceEX;
	MainView mainView;
	MainOutput mainOutput;
	
	Upp::Array<HydroClass> hydros;
	Upp::Array<MeshClass> surfs;
	
private:
	MainPlot &GetSelPlot();
	MainABForce &GetSelTab();
		
	bool closed;
};

Main &ma(Main *m = 0);

#endif
