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

#include <ScatterDraw/Unpedantic.h>
#define LAYOUTFILE <BEMRosetta/BEMRosetta/main.lay>
#include <CtrlCore/lay.h>
#include <ScatterDraw/Pedantic.h>

#include "arrange.h"

enum DataToShow {DATA_A, DATA_B, DATA_FORCE_SC, DATA_FORCE_FK, DATA_FORCE_EX, DATA_RAO};
enum DataToPlot {PLOT_A, PLOT_AINF, PLOT_B, PLOT_FORCE_SC_MA, PLOT_FORCE_SC_PH,
				 PLOT_FORCE_FK_MA, PLOT_FORCE_FK_PH, PLOT_FORCE_EX_MA, PLOT_FORCE_EX_PH, 
				 PLOT_RAO_MA, PLOT_RAO_PH, PLOT_Z_MA, PLOT_Z_PH, PLOT_TFS_MA, PLOT_TFS_PH};

template <class T>
T magnitude(const std::complex<T> &val) {
	return sqrt(pow2(val.real()) + pow2(val.imag()));
}

template <class T>
T phase(const std::complex<T> &val) {
	return atan2(val.imag(), val.real());
}

String ForceExtSafe(String fileName, String ext);
   
class HydroSource : public DataSource {
public:
	HydroSource() : data(0) {}
	HydroSource(Hydro &_data, int i, int _j_h, DataToPlot _dataToPlot, bool _show_w, bool _ndim) {
		Init(_data, i, _j_h, _dataToPlot, _show_w, _ndim);
	}
	bool Init(Hydro *_data, int i, int _j_h, DataToPlot _dataToPlot, bool _show_w, bool _ndim) 	{
		return Init(*_data, i, _j_h, _dataToPlot, _show_w, _ndim);
	}
	bool Init(Hydro &_data, int _idof, int _j_h, DataToPlot _dataToPlot, bool _show_w, bool _ndim) {
		data = &_data;
		dataToPlot = _dataToPlot;
		if (_idof >= _data.dofOrder.GetCount())
			return false;
		idof = _data.dofOrder[_idof];
		if (dataToPlot == PLOT_A || dataToPlot == PLOT_AINF || dataToPlot == PLOT_B)
			_j_h = _data.dofOrder[_j_h];
		j_h = _j_h;
		show_w = _show_w;
		ndim = _ndim;
		if (IsNullData())
			return false;
		return true;
	}
	inline bool IsNullData() {
		ASSERT(data != 0);
		switch (dataToPlot) {
		case PLOT_A:			return IsNull(data->A[0](idof, j_h));
		case PLOT_AINF:			return IsNull(data->Awinf(idof, j_h));
		case PLOT_B:			return IsNull(data->B[0](idof, j_h));
		case PLOT_FORCE_SC_MA:	return IsNull(data->sc.ma[j_h](0, idof));
		case PLOT_FORCE_SC_PH:	return IsNull(data->sc.ph[j_h](0, idof));
		case PLOT_FORCE_FK_MA:	return IsNull(data->fk.ma[j_h](0, idof));
		case PLOT_FORCE_FK_PH:	return IsNull(data->fk.ph[j_h](0, idof));
		case PLOT_FORCE_EX_MA:	return IsNull(data->ex.ma[j_h](0, idof));
		case PLOT_FORCE_EX_PH:	return IsNull(data->ex.ph[j_h](0, idof));
		case PLOT_RAO_MA:		return IsNull(data->rao.ma[j_h](0, idof));
		case PLOT_RAO_PH:		return IsNull(data->rao.ph[j_h](0, idof));
		case PLOT_Z_MA:			return IsNull(magnitude(data->Z[0]));
		case PLOT_Z_PH:			return IsNull(phase(data->Z[0]));
		case PLOT_TFS_MA:		return IsNull(magnitude(data->TFSResponse[0]));
		case PLOT_TFS_PH:		return IsNull(phase(data->TFSResponse[0]));
		default:				NEVER();	return true;
		}
	}
	virtual inline double y(int64 id) {
		ASSERT(data != 0);
		switch (dataToPlot) {
		case PLOT_A:			return data->A_(ndim, int(id), idof, j_h);
		case PLOT_AINF:			return data->Awinf_(ndim, idof, j_h);
		case PLOT_B:			return data->B_(ndim, int(id), idof, j_h);
		case PLOT_FORCE_SC_MA:	return data->F_ma_(ndim, data->sc, j_h, int(id), idof);
		case PLOT_FORCE_SC_PH:	return data->sc.ph[j_h](int(id), idof);
		case PLOT_FORCE_FK_MA:	return data->F_ma_(ndim, data->fk, j_h, int(id), idof);
		case PLOT_FORCE_FK_PH:	return data->fk.ph[j_h](int(id), idof);
		case PLOT_FORCE_EX_MA:	return data->F_ma_(ndim, data->ex, j_h, int(id), idof);
		case PLOT_FORCE_EX_PH:	return data->ex.ph[j_h](int(id), idof);
		case PLOT_RAO_MA:		return data->F_ma_(ndim, data->rao, j_h, int(id), idof);
		case PLOT_RAO_PH:		return data->rao.ph[j_h](int(id), idof);
		case PLOT_Z_MA:			return magnitude(data->Z[int(id)]);
		case PLOT_Z_PH:			return phase(data->Z[int(id)]);
		case PLOT_TFS_MA:		return magnitude(data->TFSResponse[int(id)]);
		case PLOT_TFS_PH:		return phase(data->TFSResponse[int(id)]);
		default:				NEVER();	return Null;
		}
	}
	virtual inline double x(int64 id) 	{
		ASSERT(data != 0);
		if (show_w)
			return data->w[int(id)];
		else
			return data->T[int(id)];
	}
	virtual int64 GetCount()		  	{ASSERT(data != 0); return data->Nf;}
	
private:
	Hydro *data;
	int idof, j_h;
	DataToPlot dataToPlot;
	bool show_w, ndim;
};


class MenuOptions : public WithMenuOptions<StaticRect> {
public:
	typedef MenuOptions CLASSNAME;
	MenuOptions() : bem(0) {}
	void Init(BEMData &bem);
	void Load();
	void OnSave();
	bool IsChanged();
	
private:
	BEMData *bem;
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
};

class MainSummaryCoeff : public MainSummary {
public:
	typedef MainSummaryCoeff CLASSNAME;

	void Report(const Hydro &data, int id);
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
	
	MainView() {}
	void Init(const WithMenuPlotMesh<StaticRect> &menuPlot);
	void CalcEnvelope();
	void OnPaint();
	const WithMenuPlotMesh<StaticRect> &GetMenuPlot() const {return *menuPlot;}
	void SetPaintSelect(bool _paintSelect)					{paintSelect = _paintSelect;}
	
	VolumeEnvelope env;
	
private:
	const WithMenuPlotMesh<StaticRect> *menuPlot = 0;
	bool paintSelect = true;
};


class MainViewDataEach : public StaticRect {
public:
	MainViewDataEach() {}
	void Init(MeshData &_mesh, MainView &mainView);
	void OnRefresh();
	
	TabCtrl tab;
	Splitter orig, moved, movedUnder;
	WithMainPlotList<StaticRect> arrayFacetsAll, arrayNodesOrig,
								 arrayFacetsAll2, arrayNodesMoved,
								 arrayFacetsUnder, arrayNodesUnder;
	
	class DataSourceFacets : public Convert {
	public:
		DataSourceFacets() : pmesh(0), col(0), all(true) {}
		void Init(MeshData &_mesh, int _col, bool _all);
		Value Format(const Value& q) const;
		inline const MeshData &GetMesh()	{return *pmesh;}
		
	private:
		MeshData *pmesh;
		int col;
		bool all;
	};
	class DataSourceNodes : public Convert {
	public:
		DataSourceNodes() : pmesh(0), xyz(0), origMovedUnder(0) {}
		void Init(MeshData &_mesh, int _xyz, int _origMovedUnder);
		Value Format(const Value& q) const;
		
	private:
		MeshData *pmesh;
		int xyz;
		int origMovedUnder;
	};
	
	Upp::Array<DataSourceFacets> dataSourceFacetsAll, dataSourceFacetsUnder;
	Upp::Array<DataSourceNodes> dataSourceNodesOrig, dataSourceNodesMoved;
};

class MainViewData : public StaticRect {
public:
	typedef MainViewData CLASSNAME;
	
	void Init();
	void OnAddedModel(MainView &mainView);
	void OnRefresh();
	void Clear();
	
private:
	TabCtrl tab;
	Upp::Array<MainViewDataEach> models;
};

class MainSummaryMesh : public MainSummary {
public:
	typedef MainSummaryCoeff CLASSNAME;

	void Report(const MeshData &surf, int id);
};

class MainStiffness : public WithMainStiffness<StaticRect> {
public:
	typedef MainStiffness CLASSNAME;
	
	void Init();
	void Clear();
	void Load(Upp::Array<HydroClass> &hydros);
	void Load(Upp::Array<MeshData> &surfs);
	
private:
	void AddPrepare(int &row0, int &icol0, String name, int icase, String bodyName, int ibody);
	void Add(String name, int icase, const MatrixXd &K);
	void Add(String name, int icase, String bodyName, int ibody, const Hydro &hydro);
	
};

class MainPlot : public WithMainPlot<StaticRect> {
public:
	typedef MainPlot CLASSNAME;
	void Init(int idof, double jdof_ih, DataToShow dataToShow);
	bool Load(Upp::Array<HydroClass> &hydro);
	
	Upp::Array<HydroSource> ABF_source, ABF_source2;
	Upp::Array<HydroSource> Ainf_source;
	int idof, jdof;
	double heading;
	DataToShow dataToShow;
	
	bool dim;
	int markW;
	bool show_w;
	bool showPhase;		
};

class MainABForce : public WithMainABForce<StaticRect> {
public:
	typedef MainABForce CLASSNAME;
	void Init(DataToShow dataToShow);
	void Clear();
	bool Load(BEMData &bem);
	
	Upp::Array<Upp::Array<MainPlot>> plots;

private:
	int selTab;
	bool isFilling;
	DataToShow dataToShow;
};

class MainStateSpace : public StaticRect {
public:
	typedef MainStateSpace CLASSNAME;
	void Init();
	void Init(ArrayCtrl &array);
	bool Load(BEMData &bem);
	
	Splitter splitter;
	TabCtrl tab;
	ScatterCtrl scatter;
	
	Upp::Array<HydroSource> Z_source, Z_source2, TFS_source, TFS_source2;
	Upp::Array<ArrayCtrl> arrays;
};

typedef class MainABForce MainRAO;

class MainMesh : public WithMain<StaticRect> {
public:
	typedef MainMesh CLASSNAME;
	
	void Init();
	void InitSerialize(bool ret);
	
	void AfterLoad(String file);
	bool OnLoad();
	bool OnConvertMesh();
	void OnUpdate(bool forceMoved);
	void OnHealing();
	void OnOpt();
	void OnMenuConvertArraySel() ;
	
	void LoadSelTab(BEMData &bem);
		
	void Jsonize(JsonIO &json);
		
	WithMenuMesh<StaticRect> menuOpen;
	WithMenuConvertMesh<StaticRect> menuConvert;
	WithMenuPlotMesh<StaticRect> menuPlot;
	WithMenuMeshStability<StaticRect> menuStability;

private:	
	MainView mainView;
	MainViewData mainViewData;
	SplitterButton mainVAll;
	MainSummaryMesh mainSummary;
	MainStiffness mainStiffness;

	virtual void DragAndDrop(Point p, PasteClip& d);
	virtual bool Key(dword key, int count);
};

class MainNemoh : public WithNemoh<StaticRect> {
public:
	typedef MainNemoh CLASSNAME;

	void Init(const BEMData &bem);
	void InitSerialize(bool ret);
	
	void Load(const BEMData &bem);
	void Load(const NemohCal &data);
	void Save(NemohCal &data);
		
	void Jsonize(JsonIO &json);
	
private:
	bool OnLoad();
	bool OnSave(const BEMData &bem);
	void OnCursor();
	void arrayOnCursor();
	void arrayUpdateCursor();
	void arrayClear();
	void arrayOnAdd();
	void arrayOnDuplicate();
	void arrayOnRemove();
	void InitArray();
	
	virtual void DragAndDrop(Point p, PasteClip& d);
	virtual bool Key(dword key, int count);
};

class MainBEM : public WithMain<StaticRect> {
public:
	typedef MainBEM CLASSNAME;
	
	void Init();
	void InitSerialize(bool ret);

	bool OnLoad();
	bool OnConvert();
	void OnOpt();
	bool OnFOAMM();
		
	void Jsonize(JsonIO &json);
		
	WithMenuOpen<StaticRect> menuOpen;
	WithMenuConvert<StaticRect> menuConvert;
	WithMenuPlot<StaticRect> menuPlot;
	WithMenuStateSpace<StaticRect> menuFOAMM;
	
	MainSummaryCoeff mainSummary;
	MainArrange mainArrange;
	MainABForce mainA;
	MainABForce mainB;
	MainABForce mainForceSC, mainForceFK, mainForceEX;
	MainRAO mainRAO;
	MainStateSpace mainStateSpace;
	MainStiffness mainStiffness;
		
private:
	MainPlot &GetSelPlot();
	MainABForce &GetSelTab();
	void LoadSelTab(BEMData &bem);
	
	virtual void DragAndDrop(Point p, PasteClip& d);
	virtual bool Key(dword key, int count);
	bool cancelFOAMM;
};

class Main : public TopWindow {
public:
	typedef Main CLASSNAME;
	
	Main() : closed(false) {}
	virtual ~Main();
	virtual void Close();
	void CloseMain(bool store);

	void Init();

	void OptionsUpdated();

	bool LoadSerializeJson();
	bool StoreSerializeJson();
	
	void Jsonize(JsonIO &json);

	BEMData bem;
		
private:
	TabCtrl tab;
	
	MainMesh mainMesh;
	MainNemoh mainNemoh;
	MainBEM mainBEM;
	MainOutput mainOutput;
	
	MenuOptions menuOptions;
	MenuAbout menuAbout;
	
	bool closed;
};

Main &ma(Main *m = 0);
MainBEM &mbm(MainBEM *m = 0);

#endif
