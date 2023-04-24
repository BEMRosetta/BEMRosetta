// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#ifndef _BEM_Rosetta_GUI_BEM_Rosetta_GUI_h
#define _BEM_Rosetta_GUI_BEM_Rosetta_GUI_h

#include <BEMRosetta_cl/FastOut.h>
#include "auxiliar.h"
#include "FastScatter.h"



class FreqSelector : public StaticRect {
public:
	typedef FreqSelector CLASSNAME;
	
	FreqSelector();
	
	void Init(Function <void()>_WhenAction, RectEnterSet &_frameSet) {
		OnAction = _WhenAction;
		frameSet = &_frameSet;
	}
	void Clear() {
		for (int i = 0; i < edits.size(); ++i)
			frameSet->Remove(edits[i].GetRectEnter());
		edits.Clear();
		Rect rc = add.GetRect();
		rc.right = rc.Width();
		rc.left = 0;
		add.SetRect(rc);
		pos = 0;
		Layout();
	}
	int size() 			{return edits.size();}
	double Get(int id) 		{return ~edits[id];}
	bool IsSelected(int id)	{return edits[id].IsShownFrame();}
	
	int GetSelected() {
		for (int i = 0; i < edits.size(); ++i)
			if(IsSelected(i))
				return i;
		return -1;
	}
		
	void Set(int id, double val) {
		edits[id] <<= val;
		edits[id].WhenAction();
	} 
	
	void AddField(double val = Null) {
		UnderlineCtrl<EditDouble> &edit = edits.Add();
		edit.WhenAction = OnAction;
		edit <<= val;
		frameSet->Add(edit.GetRectEnter());
		Add(edit.LeftPos(pos*fWidth + fThick, fWidth-2*fThick).TopPos(3, 21));
		Rect rc = add.GetRect();
		rc.right = (pos+1)*fWidth + rc.Width();
		rc.left = (pos+1)*fWidth;
		add.SetRect(rc);
		edit.SetFocus();
		pos++;
		Layout();
	}
	void CheckFocus(double freq) {
		for (int i = 0; i < edits.size(); ++i) {
			if (abs(double(~edits[i]) - freq) < 0.000001) {
				frameSet->OnEnter(edits[i].GetRectEnter());
				break;
			}
		}
	}

private:
	UArray<UnderlineCtrl<EditDouble>> edits;
	Button add;
	int pos = 0;
	int fWidth = 40, fThick = 3;
	
	Function <void()> OnAction;
	RectEnterSet *frameSet = 0;
};


String TabText(const TabCtrl &tab);

#define LOGTAB(tab) LOG(Format("Tab %s, Set %s, file: %s, line: %d", #tab, TabText(tab), __FILE__, __LINE__))

#include "premain.h"

#define LAYOUTFILE <BEMRosetta/main.lay>
#include <CtrlCore/lay.h>


class CompareParameters : public WithCompareParameters<StaticRect> {
public:
	void Init(ScatterDraw& scatter, SplitterButton &psplitter);
	void Init(Hydro::DataToShow data);
	void Load();
	~CompareParameters();
	
private:
	ScatterDraw *pscatter = nullptr;
	SplitterButton *psplitter = nullptr;
	Hydro::DataToShow dataToShow;
	static UVector<CompareParameters *> plist;
};


String ForceExtSafe(String fileName, String ext);
   
class HydroSource : public DataSource {
public:
	HydroSource() {}
	HydroSource(Hydro &_data, int i, int _j_h, Hydro::DataToPlot _dataToPlot, bool _show_w, bool _ndim, bool _show_ma_ph) {
		Init(_data, i, _j_h, _dataToPlot, _show_w, _ndim, _show_ma_ph);
	}
	bool Init(Hydro *_data, int i, int _j_dof, Hydro::DataToPlot _dataToPlot, bool _show_w, bool _ndim, bool _show_ma_ph) 	{
		return Init(*_data, i, _j_dof, _dataToPlot, _show_w, _ndim, _show_ma_ph);
	}
	bool Init(const Hydro &_data, int _idf, int _j_dof, Hydro::DataToPlot _dataToPlot, bool _show_w, bool _ndim, bool _show_ma_ph) {
		data = &_data;
		dataToPlot = _dataToPlot;
		idf = _idf;
		jdf = _j_dof;
		show_w = _show_w;
		show_ma_ph = _show_ma_ph;
		ndim = _ndim;
		if (IsNullData())
			return false;
		return true;
	}
	inline bool IsNullData() {
		ASSERT(data != 0);
		switch (dataToPlot) {
		case Hydro::PLOT_A:			return !data->IsLoadedA(idf, jdf);
		case Hydro::PLOT_AINF:		return !data->IsLoadedAinf(idf, jdf);
		case Hydro::PLOT_A0:		return !data->IsLoadedA0(idf, jdf);
		case Hydro::PLOT_B:			return !data->IsLoadedB(idf, jdf);
		case Hydro::PLOT_MD:		return !data->IsLoadedMD(int(idf/6), jdf);
		case Hydro::PLOT_KIRF:		return !data->IsLoadedKirf(idf, jdf);
		case Hydro::PLOT_AINFW:		return !data->IsLoadedAinf_w(idf, jdf);
		case Hydro::PLOT_FORCE_SC_1:	
		case Hydro::PLOT_FORCE_SC_2:return !data->IsLoadedFsc(idf, jdf);		// jdf: heading, idf: body
		case Hydro::PLOT_FORCE_FK_1:	
		case Hydro::PLOT_FORCE_FK_2:return !data->IsLoadedFfk(idf, jdf);
		case Hydro::PLOT_FORCE_EX_1:
		case Hydro::PLOT_FORCE_EX_2:return !data->IsLoadedFex(idf, jdf);
		case Hydro::PLOT_RAO_1:
		case Hydro::PLOT_RAO_2:		return !data->IsLoadedRAO(idf, jdf);
		case Hydro::PLOT_TFS_1:
		case Hydro::PLOT_TFS_2:		return !data->sts[idf][jdf].TFS.IsEmpty();
		case Hydro::PLOT_Z_1:
		case Hydro::PLOT_Z_2:		return !data->IsLoadedA(idf, jdf) || !data->IsLoadedAinf(idf, jdf) || !data->IsLoadedB(idf, jdf);
		default: 		NEVER();	return true;
		}
	}
	virtual inline double y(int64 id) {
		ASSERT(data != 0);
		switch (dataToPlot) {
		case Hydro::PLOT_A:			return data->A_(ndim, int(id), idf, jdf);
		case Hydro::PLOT_AINF:		return data->Ainf_(ndim, idf, jdf);
		case Hydro::PLOT_A0:		return data->A0_(ndim, idf, jdf);
		case Hydro::PLOT_B:			return data->B_(ndim, int(id), idf, jdf);
		case Hydro::PLOT_KIRF:		return data->Kirf_(ndim, int(id), idf, jdf);
		case Hydro::PLOT_AINFW:		return data->Ainf_w_(ndim, int(id), idf, jdf);
		case Hydro::PLOT_MD:		return data->Md_(ndim, idf, jdf, int(id));		// idf: body, jdf: heading, [Nb][Nh][6](Nf)
		case Hydro::PLOT_FORCE_SC_1:return show_ma_ph ? abs (data->F_(ndim, data->sc, jdf, int(id), idf)) : 
														real(data->F_(ndim, data->sc, jdf, int(id), idf));
		case Hydro::PLOT_FORCE_SC_2:return show_ma_ph ? arg (data->F_(ndim, data->sc, jdf, int(id), idf)) : 
														imag(data->F_(ndim, data->sc, jdf, int(id), idf));
		case Hydro::PLOT_FORCE_FK_1:return show_ma_ph ? abs (data->F_(ndim, data->fk, jdf, int(id), idf)) : 
														real(data->F_(ndim, data->fk, jdf, int(id), idf));
		case Hydro::PLOT_FORCE_FK_2:return show_ma_ph ? arg (data->F_(ndim, data->fk, jdf, int(id), idf)) : 
														imag(data->F_(ndim, data->fk, jdf, int(id), idf));
		case Hydro::PLOT_FORCE_EX_1:return show_ma_ph ? abs (data->F_(ndim, data->ex, jdf, int(id), idf)) : 
														real(data->F_(ndim, data->ex, jdf, int(id), idf));
		case Hydro::PLOT_FORCE_EX_2:return show_ma_ph ? arg (data->F_(ndim, data->ex, jdf, int(id), idf)) : 
														imag(data->F_(ndim, data->ex, jdf, int(id), idf));
		case Hydro::PLOT_RAO_1:		return show_ma_ph ? abs (data->R(jdf, int(id), idf)) : 
														real(data->R(jdf, int(id), idf));
		case Hydro::PLOT_RAO_2:		return show_ma_ph ? arg (data->R(jdf, int(id), idf)) : 
														imag(data->R(jdf, int(id), idf));
		case Hydro::PLOT_TFS_1:		return show_ma_ph ? abs (data->TFS_(ndim, int(id), idf, jdf)) : 
														real(data->TFS_(ndim, int(id), idf, jdf));
		case Hydro::PLOT_TFS_2:		return show_ma_ph ? arg (data->TFS_(ndim, int(id), idf, jdf)) : 
														imag(data->TFS_(ndim, int(id), idf, jdf));
		case Hydro::PLOT_Z_1:		return show_ma_ph ? abs (data->Z(ndim, int(id), idf, jdf)) : 
														real(data->Z(ndim, int(id), idf, jdf));
		case Hydro::PLOT_Z_2:		return show_ma_ph ? arg (data->Z(ndim, int(id), idf, jdf)) : 
														imag(data->Z(ndim, int(id), idf, jdf));
		default:					NEVER();	return Null;
		}
	}
	virtual inline double x(int64 id) {
		ASSERT(data != 0);
		
		if ((dataToPlot == Hydro::PLOT_AINF || dataToPlot == Hydro::PLOT_A0) && id == 1 && data->Nf > 0)
			id = data->Nf - 1;
		
		if (dataToPlot == Hydro::PLOT_KIRF)
			return data->Tirf[int(id)];
		
		if (show_w) {
			if (dataToPlot == Hydro::PLOT_A0)
				return 0;
			else if (dataToPlot == Hydro::PLOT_AINF && data->w.size() == 0)
				return double(id);
			return data->w[int(id)];
		} else {
			if (dataToPlot == Hydro::PLOT_AINF)
				return 0;
			else if (dataToPlot == Hydro::PLOT_A0 && data->w.size() == 0)
				return double(id);
			return data->T[int(id)];
		}
	}
	virtual int64 GetCount() const {
		ASSERT(data != 0); 
		
		if (dataToPlot == Hydro::PLOT_AINF)
			return show_w  ? 2 : 1;		
		if (dataToPlot == Hydro::PLOT_A0)
			return !show_w ? 2 : 1;
		if (dataToPlot == Hydro::PLOT_KIRF)
			return data->Tirf.size();
		return data->Nf;
	}
	
private:
	const Hydro *data = nullptr;
	int idf = -1, jdf = -1;
	Hydro::DataToPlot dataToPlot;
	bool show_w = false, show_ma_ph = true, ndim = 0;
};


class MenuOptions : public WithMenuOptions<StaticRect> {
public:
	typedef MenuOptions CLASSNAME;
	
	MenuOptions() {}
	void Init(BEM &bem);
	void Load();
	void OnSave();
	bool IsChanged();
	void InitSerialize(bool ret, bool &openOptions);
	
private:
	BEM *bem = nullptr;
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

class MainMesh;
	
class MainView : public WithMainView<StaticRect> {
public:
	typedef MainView CLASSNAME;
	
	MainView() {}
	void Init();
	void CalcEnvelope();
	void OnPaint();
	const WithMenuMeshPlot<StaticRect> &GetMenuPlot() const;
	const MainMesh &GetMain() const 						{return *main;}			
	void SetPaintSelect(bool _paintSelect)					{paintSelect = _paintSelect;}
	
	VolumeEnvelope env;
	
private:
	const MainMesh *main = nullptr;
	bool paintSelect = true;
};

class VideoSequence {
public:
	struct VideoRecord {
		Point3D pos, angle, c0;
		double time;
		void Jsonize(JsonIO &json) {
			json
				("pos",   pos)
				("angle", angle)
				("c0", 	  c0)
				("time",  time)
			;
		}
	};
	
	virtual ~VideoSequence() = 0; 
	
	virtual bool Load(String fileName) = 0;
	virtual bool Save(String fileName) const = 0;
	virtual void Clear() = 0;
	virtual int GetRecord(VideoRecord &record, double time) = 0;
	virtual int size() = 0;
	template <class T>  bool Is() const	{return dynamic_cast<const T*>(this);}
	virtual bool HasTime() = 0;
	virtual void SetDeltaTime(double dt) = 0;
};

class BasicVideoSequence : public VideoSequence {
public:
	virtual ~BasicVideoSequence() {}
	
	virtual bool Load(String fileName) {
		return LoadFromJsonFile(*this, fileName);
	}
	virtual bool Save(String fileName) const {
		return StoreAsJsonFile(*this, fileName, true);
	}
	virtual void Clear() {
		reg.Clear();
	}
	virtual int GetRecord(VideoRecord &record, double time) {
		int askedId = min(reg.size()-1, int(time/deltaT));
		
		if (askedId > lastId) {
			lastId++;
			record = clone(reg[lastId]);
			record.time = askedId*deltaT;
			if (lastId == reg.size()-1) {
				lastId = -1;
				return 2;		// Last pos
			}
			return 1;
		} else 
			return 0;			// Pos already processed
	}
	virtual int size() {return reg.size();}
	void Add(const Point3D &pos, const Point3D &angle, const Point3D &c0) {
		auto &rg = reg.Add();
		rg.pos = clone(pos);
		rg.angle = clone(angle);
		rg.c0 = clone(c0);
		rg.time = 0;
	}
	virtual bool HasTime() {return false;};
	virtual void SetDeltaTime(double dt)	{deltaT = dt;}
	
	void Jsonize(JsonIO &json) {
		json
			("reg", reg)
		;
	}
	
private:
	UArray<VideoRecord> reg;
	double deltaT;
	int lastId = -1;
};

class VideoCtrl : public WithVideoCtrl<StaticRect> {
public:
	typedef VideoCtrl CLASSNAME;
	
	enum VideoType {UNKNOWN, JSON, OpenFAST, CSV};
	
	VideoCtrl() {}
	void Init(Function <int(UVector<int> &ids)> _GetMeshId, Function <void(int id, const UVector<int> &ids, const Point3D &pos, const Point3D &angle, const Point3D &c0, bool full, bool saveBitmap)> _Action);
	
	~VideoCtrl() {
		video.Clear();
	}
	
	void OnPlay();
	void OnRecord();
	void OnRecordClear();
	void OnLoad();
	void OnSave();
	
	VideoType GetVideoType(String name);
	
	bool IsBasicOpened() {return video && video->Is<BasicVideoSequence>();}
	
	void AddReg(const Point3D &pos, const Point3D &angle, const Point3D &c0) {
		if (!IsBasicOpened())
			return;
		static_cast<BasicVideoSequence&>(*video).Add(pos, angle, c0);
		numFrames <<= static_cast<BasicVideoSequence&>(*video).size();
	}
	void AddReg(const Point3D &angle, const Point3D &c0) {
		if (!IsBasicOpened())
			return;
		static_cast<BasicVideoSequence&>(*video).Add(Point3D(0, 0, 0), angle, c0);
		numFrames <<= static_cast<BasicVideoSequence&>(*video).size();
	}
	void AddReg(const Point3D &pos) {
		if (!IsBasicOpened())
			return;
		static_cast<BasicVideoSequence&>(*video).Add(pos, Point3D(0, 0, 0), Point3D(0, 0, 0));	
		numFrames <<= static_cast<BasicVideoSequence&>(*video).size();
	}	
	void ClearReg() {
		if (!IsBasicOpened())
			return;
		static_cast<BasicVideoSequence&>(*video).Clear();
		numFrames <<= 0;
	}
	void Jsonize(JsonIO &json) {
		json
			("deltaT", deltaT)
			("editFile", editFile)
			("opSaveBitmap", opSaveBitmap)
		;
	}
	
private:	
	void TimerFun();
	bool playing = false, recording = false;
	RealTimeStop time;
	double dT;
	
	UVector<int> ids;
	int meshId = -1;
	
	Function <int(UVector<int> &ids)> GetMeshId;
	Function<void(int meshid, const UVector<int> &ids, const Point3D &pos, const Point3D &angle, const Point3D &c0, bool full, bool saveBitmap)> Action;
	
	One<VideoSequence> video;
};

class MainViewDataEach : public StaticRect {
public:
	typedef MainViewDataEach CLASSNAME;
	
	MainViewDataEach() {}
	void Init(Mesh &_mesh, MainView &mainView);
	void OnRefresh();
	
	TabCtrl tab;
	Splitter moved, movedUnder;
	WithMainPlotList<StaticRect> arrayFacetsAll2, arrayNodesMoved,
								 arrayFacetsUnder, arrayNodesUnder;
	StatusBar status;
	
	class DataSourceFacets : public Convert {
	public:
		DataSourceFacets() : pmesh(0), col(0), all(true) {}
		void Init(Mesh &_mesh, int _col, bool _all);
		Value Format(const Value& q) const;
		inline const Mesh &GetMesh()	{return *pmesh;}
		
	private:
		Mesh *pmesh;
		int col;
		bool all;
	};
	class DataSourceNodes : public Convert {
	public:
		DataSourceNodes() : pmesh(0), xyz(0), origMovedUnder(0) {}
		void Init(Mesh &_mesh, int _xyz, int _origMovedUnder);
		Value Format(const Value& q) const;
		
	private:
		Mesh *pmesh;
		int xyz;
		int origMovedUnder;
	};
	
	UArray<DataSourceFacets> dataSourceFacetsAll, dataSourceFacetsUnder;
	UArray<DataSourceNodes> dataSourceNodesMoved, dataSourceNodesUnder;

private:
	void OnTimer();
	TimeCallback timeCallback;
	void UpdateStatus(bool under);
	UVector<int> selectedPanels, selectedNodes;  
	int lastSel = -1;
};

class MainViewData : public StaticRect {
public:
	typedef MainViewData CLASSNAME;
	
	void Init();
	void OnAddedModel(MainView &mainView);
	void OnRefresh();
	void Clear();
	void ReLoad(MainView &mainView);
	
private:
	TabCtrl tab;
	UArray<MainViewDataEach> models;
};

class MainSummaryMesh : public MainSummary {
public:
	typedef MainSummaryCoeff CLASSNAME;

	void Report(const UArray<Mesh> &surfs, int id);
};

class MainMatrixKA : public WithMainMatrixKA<StaticRect> {
public:
	typedef MainMatrixKA CLASSNAME;
	
	void Init(Hydro::DataMatrix what);
	void OnOp(int id);
	void Clear();
	bool Load(UArray<HydroClass> &hydros, const UVector<int> &ids, bool ndim);
	void Load(UArray<Mesh> &surfs, const UVector<int> &ids);
	
	void Jsonize(JsonIO &json) {
		bool opdigits = ~opDigits;
		if (json.IsLoading()) 
			opdigits = true;
		else
			opdigits = ~opDigits;
		json
			("opEmptyZero", opEmptyZero)
			("numDigits", numDigits)
			("numDecimals", numDecimals)
			("expRatio", expRatio)
			("opDigits", opdigits)
			("hstFolder", hstFolder)
		;
		if (json.IsLoading()) {
			opDigits <<= opdigits;
			opDecimals <<= !opdigits;
			if (IsNull(opEmptyZero)) 
				opEmptyZero <<= false;
			if (IsNull(numDigits)) 
				numDigits <<= 6;
			if (IsNull(numDecimals)) 
				numDecimals <<= 0;
			if (IsNull(expRatio)) 
				expRatio <<= 5;
			if (IsNull(hstFolder))
				hstFolder = GetDesktopFolder();
		}
	}
	
private:
	void AddPrepare(int &row0, int &icol0, String name, int icase, String bodyName, int ibody, int idc);
	void Add(const Mesh &mesh, int icase, bool button);
	void Add(String name, int icase, String bodyName, int ibody, const Hydro &hydro, int idc, bool ndim);
	void PrintData();
		
	UArray<Eigen::MatrixXd> data;
	UArray<int> row0s, col0s;
	
	Hydro::DataMatrix what;
	bool Ndim;
	
	String hstFolder;
};

class MainGZ : public WithMainGZ<StaticRect> {
public:
	typedef MainGZ CLASSNAME;

	void Init();
	void Clear(bool force);
	
	void OnUpdate();
	void Jsonize(JsonIO &json);
	
private:
	UArray<UVector<double>> datagz, dataMoment;
	UVector<double> dangle, mingz;
	int idOpened = Null;
	
	ScatterCtrl scatter;
	ArrayCtrl array;
};

class MainBEM;

class MainPlot : public StaticRect {
public:
	typedef MainPlot CLASSNAME;
	
	~MainPlot()	{Clear();}
		
	void Init(bool vert);
	void Init(int idf, double jdf_ih, Hydro::DataToShow dataToShow, double _heading1 = Null);
	bool Load(const UArray<HydroClass> &hydro, const MainBEM &mbm, const UVector<int> &ids);
	bool Load(const Hydro &hy, const MainBEM &mbm);
	void LoadEach(const Hydro &hy, int id, bool &loaded, int idc = -1);
	void Clear();
	void RefreshScatter()	{scatt.Refresh();	scatP.Refresh();}
	
	UArray<HydroSource> ABFZ_source, ABFZ_source2;
	UArray<HydroSource> Ainf_source, A0_source;	
	UArray<HydroSource> TFS_source, TFS_source2;
		
	int plot_idf, plot_jdf;
	double heading0, heading1;
	Hydro::DataToShow dataToShow;
	
	bool dim;
	int markW;
	bool show_w, show_ma_ph;
	
	ScatterCtrl scatt, scatP;
	Splitter splitter;
	CompareParameters compare;
	SplitterButton splitCompare;

private:
	bool isInit = false;	
};

class MainABForce : public StaticRect {
public:
	typedef MainABForce CLASSNAME;
	
	void Init(Hydro::DataToShow dataToShow);
	void Clear();
	bool Load(BEM &bem, const UVector<int> &ids, int ih = Null);
	
	void UpdateHead(BEM &bem, int ih);
	void UpdateHeadMD(BEM &bem, int ih);
		
	TabCtrl tab;
	UArray<UArray<MainPlot>> plots;

private:
	int selTab;
	bool isFilling;
	Hydro::DataToShow dataToShow;
};

class MainStateSpacePlot : public StaticRect {
public:
	typedef MainStateSpacePlot CLASSNAME;
	
	void Init(int _idf, int _jdf);
	bool Load(UArray<HydroClass> &hydro, const UVector<int> &ids, const MainBEM &mbm);
	void InitArray(ArrayCtrl &array);
	
	MainPlot mainPlot;
	Splitter splitterTab;
	TabCtrl tab;
	UArray<ArrayCtrl> arrays;
};

class MainStateSpace : public WithMainStateSpace<StaticRect> {
public:
	typedef MainStateSpace CLASSNAME;
	
	void Init();
	void Clear();
	void Init(ArrayCtrl &array);
	bool Load(BEM &bem, const UVector<int> &ids);
	
	UArray<UArray<MainStateSpacePlot>> plots;
	
private:
	int selTab;
	bool isFilling;
};

typedef class MainABForce MainRAO;

class MainMesh : public WithMainBEMMesh<StaticRect> {
public:
	typedef MainMesh CLASSNAME;
	
	void Init();
	void InitSerialize(bool ret);
	
	enum Action {NONE, MOVE, ROTATE};
	
	void AfterAdd(String file, int num);
	void After();
	bool OnLoad();
	void OnRemove();
	void OnReset();
	void OnRemoveSelected(bool all);
	void OnJoin();
	void OnSplit();
	void OnConvertMesh();
	void OnUpdate(Action action, bool fromMenuProcess);
	void OnScale();
	void OnTranslateArchimede(bool fromMenuProcess);
	void OnArchimede();
	void OnUpdateMass();
	void OnHealing(bool basic);
	void OnOrientSurface();
	void OnImage(int axis);
	void OnOpt();
	void OnArraySel();
	void OnMenuOpenArraySel();
	void OnMenuProcessArraySel();
	void OnMenuAdvancedArraySel();
	void OnMenuMoveArraySel();
	void OnAddPanel();
	void OnAddRevolution();
	void OnAddPolygonalPanel();
	void OnAddWaterSurface(char c);
	void UpdateButtons();
			
	void AddRow(const Mesh &surf);
	void RemoveRow(int row);
	
	void LoadSelTab(BEM &bem);
		
	void Jsonize(JsonIO &json);
		
	WithMenuMesh<StaticRect> menuOpen;
	WithMenuMeshPlot<StaticRect> menuPlot;
	WithMenuMeshProcess<StaticRect> menuProcess;
	WithMenuMeshMove<StaticRect> menuMove;
	WithMenuMeshEdit<StaticRect> menuEdit;
	
	bool GetShowMesh()			{return menuPlot.showMesh;}
	bool GetShowUnderwater()	{return menuPlot.showUnderwater;}

private:	
	MainView mainView;
	VideoCtrl videoCtrl;
	MainViewData mainViewData;
	SplitterButton splitterAll, splitterVideo;
	MainSummaryMesh mainSummary;
	MainMatrixKA mainStiffness;
	MainGZ mainGZ;

	UArray<Option> optionsPlot;
	
	virtual void DragAndDrop(Point p, PasteClip& d);
	virtual bool Key(dword key, int count);
	void LoadDragDrop();
	
	Button::Style styleRed, styleGreen, styleBlue;
	
	TimeCallback timerDrop;
	UVector<String> filesToDrop;
	
	String saveFolder;
	int dropExportId;
};

class MainMeshW : public TopWindow {
public:
	typedef MainMeshW CLASSNAME;
	
	void Init(MainMesh &_mesh, const Image &icon, const Image &largeIcon, Function <void()> _WhenClose);
	virtual void Close() {
		WhenClose();
		delete this;
	}
	
	MainMesh mesh;

private:
	Function <void()> WhenClose;
};

class MainSolver : public WithMainSolver<StaticRect> {
public:
	typedef MainSolver CLASSNAME;

	void Init(const BEM &bem);
	void InitSerialize(bool ret);
	
	void Load(String file, const BEM &bem);
	void Load(const BEM &bem);
	
	bool Save(BEMCase &data, bool isNemoh);
	
	void Jsonize(JsonIO &json);
	
	WithMainSolver_Bodies<StaticRect> bodies;
	WithMainSolver_Nemoh<StaticRect> nemoh;
	WithMainSolver_HAMS<StaticRect> hams;
	
	CtrlScroll nemohScroll;
	//CtrlScroll hamsScroll;
	
private:
	bool OnLoad(const BEM &bem);
	bool OnSave(const BEM &bem);
	void OnCursor();
	void arrayOnCursor();
	bool ArrayUpdateCursor();
	void arrayClear();
	void arrayOnAdd();
	void arrayOnDuplicate();
	void arrayOnRemove();
	void InitArray(bool isNemoh);
	
	void InitGrid(GridCtrl &grid, EditDouble edit[]);
	void LoadMatrix(GridCtrl &grid, const Eigen::MatrixXd &mat);
	
	EditDouble editMass[6], editLinear[6], editQuadratic[6], editInternal[6], editExternal[6], editAdd[6];
	int dropSolverVal = 0;
};

class ArrayFields {
public:
	//void Init(Mooring &mooring) = 0;
	virtual void Load() = 0;
	virtual void Save() = 0;

private:
	virtual void InitArray() = 0;	
	virtual bool ArrayUpdateCursor() = 0;
	virtual void ArrayOnCursor() = 0;
	virtual void ArrayClear() = 0;
	
protected:	
	ArrayCtrl *parray = nullptr;

	void ArrayOnAdd() {
		ArrayCtrl &array = *parray;
		if (array.GetCount() == 0)
			InitArray();
		array.Add();
		array.SetCursor(array.GetCount()-1);	
		ArrayClear();
		ArrayUpdateCursor();
	}	
	void ArrayOnDuplicate() {
		ArrayCtrl &array = *parray;
		if (array.GetCount() == 0) {
			Exclamation(t_("No row available to duplicate"));
			return;
		}
		int id = array.GetCursor();
		if (id < 0) {
			Exclamation(t_("Please select the row to duplicate"));
			return;
		}
		int nr = id + 1;
		array.Insert(nr);
		for (int c = 0; c < array.GetColumnCount(); ++c)
			array.Set(nr, c, array.Get(id, c));
		array.Disable();
		array.SetCursor(nr);
		array.Enable();
		ArrayOnCursor();
	}
	
	void ArrayOnRemove() {
		ArrayCtrl &array = *parray;
		if (array.GetCount() == 0) {
			Exclamation(t_("No row available to remove"));
			return;
		}
		int id = array.GetCursor();
		if (id < 0) {
			Exclamation(t_("Please select the row to remove"));
			return;
		}
		array.Remove(id);
		if (id >= array.GetCount()) {
			id = array.GetCount()-1;
			if (id < 0) {
				ArrayClear();
				ArrayUpdateCursor();
				return;
			}
		} 
		array.SetCursor(id);
	}

};

class MainMoor_LinesTypes : public WithMainMoor_LinesTypes<StaticRect>, public ArrayFields {
public:
	typedef MainMoor_LinesTypes CLASSNAME;
	
	void Init(Mooring &mooring);
	virtual void Load();
	virtual void Save();
	
private:
	Mooring *pmooring = nullptr;
	
	virtual void InitArray();	
	virtual bool ArrayUpdateCursor();
	virtual void ArrayOnCursor();
	virtual void ArrayClear();
};

class MainMoor_LineProperties : public WithMainMoor_LineProperties<StaticRect>, public ArrayFields  {
public:
	typedef MainMoor_LineProperties CLASSNAME;
	
	void Init(Mooring &mooring);
	virtual void Load();
	void LoadDrop();
	virtual void Save();
	
private:
	Mooring *pmooring = nullptr;
	
	virtual void InitArray();	
	virtual bool ArrayUpdateCursor();
	virtual void ArrayOnCursor();
	virtual void ArrayClear();
};

class MainMoor_Connections : public WithMainMoor_Connections<StaticRect>, public ArrayFields  {
public:
	typedef MainMoor_Connections CLASSNAME;
	
	void Init(Mooring &mooring);
	virtual void Load();
	virtual void Save();
	
private:
	Mooring *pmooring = nullptr;
	
	virtual void InitArray();	
	virtual bool ArrayUpdateCursor();
	virtual void ArrayOnCursor();
	virtual void ArrayClear();
};

class MainMoor : public StaticRect {
public:
	typedef MainMoor CLASSNAME;

	void Init();
	
	Splitter splitter;
	Box scat;
	ScatterCtrl scatLateral, scatUp;
	WithMainMoorRight<StaticRect> right;
	MainMoor_LinesTypes lineTypes;
	MainMoor_LineProperties lineProperties;
	MainMoor_Connections lineConnections;
	
	void Jsonize(JsonIO &json) {
		if (json.IsLoading())
			right.file <<= "";
		json
			("file", right.file)
		;
	}
	
private:
	Mooring mooring;
	
	UVector<double> px, py, pz;
	
	bool OnLoad();
	bool OnSave();
	void OnUpdate();
};

class MainDecay_Files : public WithMainDecay_Files<StaticRect>, public ArrayFields {
public:
	typedef MainDecay_Files CLASSNAME;
	
	virtual void Init();
	virtual void Load();
	virtual void Save();
	
private:
	virtual void InitArray();	
	virtual bool ArrayUpdateCursor();
	virtual void ArrayOnCursor();
	virtual void ArrayClear();
};

class MainDecay : public WithMainDecay<StaticRect> {
public:
	typedef MainDecay CLASSNAME;

	void Init();
	
	void Jsonize(JsonIO &json) {
		if (json.IsLoading())
			file <<= "";
		json
			("file", file)
		;
	}
	
private:
	MainDecay_Files files;
	
	bool OnLoad();
	bool OnSave();
	void OnUpdate();
};

class MainSetupFOAMM : public WithMainStateSpaceSetup<StaticRect> {
public:
	typedef MainSetupFOAMM CLASSNAME;
	
	void Init();
	
	void Load()				{WhenSelArrayCases();}
	void WhenSelArrayModel(int id, BEM &bem);
	void WhenSelArrayCases();
	void WhenArrayCases();
	void WhenArrayFreq();
	String Check(double fromFreq, double toFreq, String freqs);
	
	void WhenFocus();
	void OnPainter(Painter &w, ScatterCtrl *scat);
	void OnMouse(Point p, dword, ScatterCtrl::MouseAction action, ScatterCtrl *scat);
	//void OnMove(Point p);
	
	bool Get(UVector<int> &ibs, UVector<int> &idfs, UVector<int> &jdfs,
		UVector<double> &froms, UVector<double> &tos, UVector<UVector<double>> &freqs); 
	void Clear();
	
	MainPlot plots;
	
	FreqSelector selector;

private:
	UArray<Option> options;
	int id = -1;
	RectEnterSet frameSet;
};

class QTFTabDof : public StaticRect {
public:
	typedef QTFTabDof CLASSNAME;

	void Init(int posSplitter, int ib, int idof);
	void Load(const Hydro &hd, int ib, int ih, int idof, bool ndim, bool show_w, bool show_ma_ph, bool isSum, bool opBilinear, bool showPoints, bool fromY0, bool autoFit, int posSplitter);
	
	Splitter splitter;
	
	int type = 0;
	bool isSum;
	int ib, ih, idof;
	bool ndim, show_w;
	bool showPoints, fromY0, autoFit;
	std::complex<double> head;
	
	struct Data {
		ArrayCtrl array;
		Box sc;
		ScatterCtrl surf, scatter;
		UVector<double> zData, xAxis;
		TableDataVector dataSurf;
		UArray<UArray<Pointf>> dataPlot;
		bool isUp;
		bool show_ma_ph;
		String labelY, units, ma_ph;
	} up, down;
	
	double GetData(const Hydro &hd, const Data &data, int idh, int ifr1, int ifr2);
	double qwT(const Hydro &hd, int id);
	
private:
	Box leftsplit, rightsplit;
	Pointf pf = Null;
	
	void UpdateArray(const Hydro &hd, bool show_ma_ph, Data &data, bool opBilinear);
	void OnClick(Point p, int idof, ScatterCtrl::MouseAction action);
	void DoClick(Data &up, int idof);
		
	void OnPainter(Painter &w)		{OnPaint(w);}
	void OnDraw(Draw &w)			{OnPaint(w);}
	
	template <class T>
	void OnPaint(T& w) {
		if (!up.surf.IsSurf())
			return;
		ScatterCtrl &s = up.surf;
		double mn = s.GetSurfMinX(), mx = s.GetSurfMaxX();
		if (type == 0)
			DrawLineOpa(w, s.GetPosX(mn), s.GetPosY(mn), s.GetPosX(mx), s.GetPosY(mx), 1, 1, 3, LtRed(), "2 2");
		else if (type == 1)
			DrawLineOpa(w, s.GetPosX(mn), s.GetPosY(mx), s.GetPosX(mx), s.GetPosY(mn), 1, 1, 3, LtRed(), "2 2");
		else {
			if (IsNull(pf))
				return;
			int px = fround(up.surf.GetPosX(pf.x)), py = fround(up.surf.GetPosY(pf.y));
			if (type == 2)
				DrawLineOpa(w, s.GetPosX(mn), py, 			 s.GetPosX(mx), py,    		   1, 1, 3, LtRed(), "2 2");	
			else
				DrawLineOpa(w, px, 			  s.GetPosY(mn), px,    		s.GetPosY(mx), 1, 1, 3, LtRed(), "2 2");	
		}
	}
};

class MainQTF : public WithMainQTF<StaticRect> {
public:
	typedef MainQTF CLASSNAME;
	
	void Init(MainBEM &parent);	
	bool Load();
	void Unload(int idf = -1);
	void OnHeadingsSel(ArrayCtrl *listHead);
	void OnSurf();
	
	enum Mag {MAGNITUDE, PHASE, REAL, IMAGINARY};
	enum Show {FSUM, FDIFFERENCE};
	
private:
	MainBEM *_mbm = nullptr;
	UArray<QTFTabDof> dof;
	bool isLoading = false;
	
	int idof = -1, ib = -1;
	int posSplitter = 3000;
	std::complex<double> head;
	bool isSumm = true;
};

class MenuFOAMM : public WithMenuStateSpace<StaticRect> {
public:
	typedef MenuFOAMM CLASSNAME;
	
	void Init(MainBEM &_mainBEM, MainSetupFOAMM &_setup);
	bool OnFOAMM();
	void OnCursor();
	
	void Clear();

private:
	bool isCancelled;
	MainSetupFOAMM *setup = nullptr;
};

class MainBEM : public WithMainBEMMesh<StaticRect> {
public:
	typedef MainBEM CLASSNAME;
	
	void Init();
	void InitSerialize(bool ret);

	bool OnLoad();
	bool OnLoadFile(String file);
	void OnConvert();
	void OnOpt();
	void OnRemove();
	void OnRemoveSelected(bool all);
	void OnJoin();
	void OnSymmetrizeForces(bool xAxis);
	void OnSymmetrize();
	void OnDuplicate();
	void OnKirfAinf(Hydro::DataToPlot param);
	void OnRAO();
	void OnOgilvie();
	void OnConvergence();
	void OnUpdateCrot();
	void OnDeleteHeadingsFrequencies();
	void OnResetForces();
	void OnMultiplyDOF(bool isReset);
	void OnSwapDOF();
	void OnDescription();
	void OnMenuAdvancedArraySel();
	void OnSelListLoaded();
	void UpdateButtons();
	void ShowMenuPlotItems();
	void OnABForces();
	void OnQTF();
	void OnABForcesZero();
	void OnQTFZero();
	void OnQTF_MD();
	
	void AfterBEM();
		
	int GetIdOneSelected(bool complain = true);
		
	void Jsonize(JsonIO &json);
		
	WithMenuOpen<StaticRect> menuOpen;
	WithMenuProcess<StaticRect> menuProcess;
	WithMenuProcess2<StaticRect> menuProcess2;
	WithMenuAdvanced<StaticRect> menuAdvanced;
	WithMenuPlot<StaticRect> menuPlot;
	MenuFOAMM menuFOAMM;
	
	MainSummaryCoeff mainSummary;
	MainABForce mainA;
	MainABForce mainB;
	MainABForce mainAinfw;
	MainABForce mainK;
	MainABForce mainForceSC, mainForceFK, mainForceEX;
	MainABForce mainMD;
	MainRAO mainRAO;
	MainStateSpace mainStateSpace;
	MainMatrixKA mainMatrixK;
	MainMatrixKA mainMatrixA;
	MainMatrixKA mainMatrixDlin;
	MainMatrixKA mainMatrixM;
	MainSetupFOAMM mainSetupFOAMM;
	MainQTF mainQTF;
		
private:
	ScatterCtrl &GetSelScatter();
	MainABForce &GetSelABForce();
	MainStateSpace &GetSelStateSpace();
	void LoadSelTab(BEM &bem);

	int AskQtfHeading(const Hydro &hydro);
		
	virtual void DragAndDrop(Point p, PasteClip& d);
	virtual bool Key(dword key, int count);
	void LoadDragDrop();
	String BEMFile(String fileFolder) const;
	
	TimeCallback timerDrop;
	UVector<String> filesToDrop;
	
	String saveFolder;
	int dropExportId;
};

class MainBEMW : public TopWindow {
public:
	typedef MainBEMW CLASSNAME;
	
	void Init(MainBEM &_bem, const Image &icon, const Image &largeIcon, Function <void()> _WhenClose);
	virtual void Close() {
		WhenClose();
		delete this;
	}
	
	MainBEM bem;
	
private:
	Function <void()> WhenClose;
};


class Main : public WithMain<TopWindow> {
public:
	typedef Main CLASSNAME;
	
	Main() : closed(false) {}
	virtual ~Main() noexcept;
	virtual void Close();
	void CloseMain(bool store);

	void Init();

	void OptionsUpdated(double rho, double g, int dofType, int headingType);

	bool LoadSerializeJson(bool &firstTime, bool &openOptions);
	bool StoreSerializeJson();
	
	void Jsonize(JsonIO &json);

	BEM bem;
	
	void Status(String str = String(), int time = 6000)	{
		if (!str.IsEmpty()) {
			bar.Temporary(str, time);
			BEM::Print("\n" + str);
		} else
			bar.EndTemporary();
		ProcessEvents();
	}
	
	void SetLastTab()	{tab.Set(lastTab);};
	
	void AddWindow()	{numWindows++;}
	void DeleteWindow()	{numWindows--;}
	
private:
	int numWindows = 0;
	
	int lastTab;
	
	MainSolver mainSolver;
	MainBEM mainBEM;
	MainOutput mainOutput;
	MainMesh mainMesh;
	FastScatterTabs mainFAST;
	MainMoor mainMoor;
	MainDecay mainDecay;
	
	MenuOptions menuOptions;
	CtrlScroll menuOptionsScroll;
	MenuAbout menuAbout;
	
	bool closed;
	
	StatusBar bar;
};

ArrayCtrl &ArrayModel_Init(ArrayCtrl &array, bool push = false);
void ArrayModel_Add(ArrayCtrl &array, String codeStr, String title, String fileName, int id);
void ArrayModel_Add(ArrayCtrl &array, String codeStr, String title, String fileName, int id, 
					UArray<Option> &option, Function <void()>OnPush);
int ArrayModel_Id(const ArrayCtrl &array);
int ArrayModel_Id(const ArrayCtrl &array, int row);
int ArrayModel_IdMesh(const ArrayCtrl &array);
int ArrayModel_IdMesh(const ArrayCtrl &array, int row);
int ArrayModel_IdHydro(const ArrayCtrl &array);
int ArrayModel_IdHydro(const ArrayCtrl &array, int row);
UVector<int> ArrayModel_IdsHydro(const ArrayCtrl &array);		
UVector<int> ArrayModel_IdsMesh(const ArrayCtrl &array);
void ArrayModel_IdsHydroDel(ArrayCtrl &array, const UVector<int> &ids);
void ArrayModel_RowsHydroDel(ArrayCtrl &array, const UVector<int> &ids);
bool ArrayModel_IsVisible(const ArrayCtrl &array, int row);
bool ArrayModel_IsSelected(const ArrayCtrl &array, int row);
const Color& ArrayModel_GetColor(const ArrayCtrl &array, int row);
String ArrayModel_GetFileName(ArrayCtrl &array, int row = -1);
String ArrayModel_GetTitle(ArrayCtrl &array, int row = -1);

void ArrayModel_Change(ArrayCtrl &array, int id, String codeStr, String title, String fileName);
		
Main &ma(Main *m = 0);
BEM &Bem();
void Status(String str = String(), int time = 2000);

	
#endif
