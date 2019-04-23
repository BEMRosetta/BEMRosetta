#ifndef _BEM_Rosetta_GUI_BEM_Rosetta_GUI_h
#define _BEM_Rosetta_GUI_BEM_Rosetta_GUI_h

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


#include "arrange.h"

enum DataToShow {DATA_A, DATA_B, DATA_FORCE_SC, DATA_FORCE_FK, DATA_FORCE_EX, DATA_RAO};
enum DataToPlot {PLOT_A, PLOT_AINF, PLOT_B, PLOT_FORCE_SC_MA, PLOT_FORCE_SC_PH,
				 PLOT_FORCE_FK_MA, PLOT_FORCE_FK_PH, PLOT_FORCE_EX_MA, PLOT_FORCE_EX_PH, PLOT_RAO_MA, PLOT_RAO_PH};
    
class HydroSource : public DataSource {
public:
	HydroSource() : data(0) {}
	HydroSource(Hydro &data, int i, int j_h, DataToPlot dataToPlot, bool show_w) {
		Init(data, i, j_h, dataToPlot, show_w);
	}
	bool Init(Hydro *data, int i, int j_h, DataToPlot dataToPlot, bool show_w) 	{
		return Init(*data, i, j_h, dataToPlot, show_w);
	}
	bool Init(Hydro &data, int i, int j_h, DataToPlot dataToPlot, bool show_w) {
		this->data = &data;
		this->i = data.dofOrder[i];
		if (dataToPlot == PLOT_A || dataToPlot == PLOT_AINF || dataToPlot == PLOT_B)
			j_h = data.dofOrder[j_h];
		this->j_h = j_h;
		this->dataToPlot = dataToPlot;
		this->show_w = show_w;
		if (IsNull(y(0)))
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
		case PLOT_RAO_MA:		return data->rao.ma[j_h](int(id), i);
		case PLOT_RAO_PH:		return data->rao.ph[j_h](int(id), i);
		default:				NEVER();	return Null;
		}
	}
	virtual inline double x(int64 id) 	{
		if (show_w)
			return data->w[int(id)];
		else
			return data->T[int(id)];
	}
	virtual int64 GetCount()		  	{return data->Nf;}
	
private:
	Hydro *data;
	int i, j_h;
	DataToPlot dataToPlot;
	bool show_w;
};

class MenuOptions : public WithMenuOptions<StaticRect> {
public:
	typedef MenuOptions CLASSNAME;
	MenuOptions() : md(0) {}
	void Init(BEMData &md);
	void Load();
	void OnSave();
	bool IsChanged();
	
private:
	BEMData *md;
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

typedef class MainABForce MainRAO;

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
	
	void WindowWamitAdditionalData(BEMData &md, HydroClass &data);
		
	bool LoadSerializeJson() {
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
		
		if (!ret || IsNull(menuOpen.optLoadIn)) 
			menuOpen.optLoadIn = 1;

		if (!ret || IsNull(menuPlot.autoFit)) 
			menuPlot.autoFit = true;
		
		if (!ret || IsNull(menuPlot.opwT)) 
			menuPlot.opwT = 0;
	
		if (!ret || IsNull(menuPlot.showPoints)) 
			menuPlot.showPoints = true;
		
		if (!ret || IsNull(menuPlot.showPhase)) 
			menuPlot.showPhase = true;

		if (!ret || IsNull(menuView.optLoadIn)) 
			menuView.optLoadIn = 1;
		
		if (!ret || IsNull(menuConvert.opt)) 
			menuConvert.opt = 0;
		
		if (!ret)
			menuTab.Set(menuAbout);
		
		return ret;
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
	MenuOptions menuOptions;
	MenuAbout menuAbout;
	
	MainSummary mainSummary;
	MainArrange mainArrange;
	MainABForce mainA;
	MainABForce mainB;
	MainABForce mainForceSC, mainForceFK, mainForceEX;
	MainRAO mainRAO;
	MainView mainView;
	MainOutput mainOutput;
	
	BEMData md;
	Upp::Array<MeshClass> surfs;
	
private:
	MainPlot &GetSelPlot();
	MainABForce &GetSelTab();
		
	bool closed;
};

Main &ma(Main *m = 0);

#endif
