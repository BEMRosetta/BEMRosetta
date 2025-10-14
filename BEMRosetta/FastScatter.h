// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#ifndef _BEMRosetta_BEMRosetta_FastScatter_h_
#define _BEMRosetta_BEMRosetta_FastScatter_h_

template <class T>
class WithClick : public T {
private:
	virtual Image CursorImage(Point, dword) 	{return Image::Hand();}
	virtual void LeftDown(Point, dword) 		{T::WhenAction();}	
};

#define LAYOUTFILE <BEMRosetta/FastScatter.lay>
#include <CtrlCore/lay.h>

class FastScatter;

class ScatterLeft : public StaticRect {
public: 
	typedef ScatterLeft CLASSNAME;
		
	ScatterLeft() {
		Add(box.SizePos());
	}
	void ClearScatter() {
		scatter.Clear();
		box.Clear();
	}
	void AddScatter() {
		ScatterCtrl &s = scatter.Add();
		s.ShowAllMenus().SetMode(ScatterDraw::MD_DRAW);
		if (scatter.size() > 1)
			s.LinkedWith(scatter[0]);
		box.Add(s, box.GetRows(), 0);
	}
	int AddAll(bool addScatter) {
		dataFast.Add();
		dataSource.Add();
		if (addScatter || scatter.size() == 0) 
			AddScatter();
		return dataFast.size()-1;
	}
	void EnableX(bool en = true) {
		for (ScatterCtrl &s : scatter)
			s.Enable(en);
		Enable(en);
	}
	int Find(String fileName) {
		for (int iff = 0; iff < dataFast.size(); ++iff) {
			String file = dataFast[iff].GetFileName();
			if (fileName == file || 
				(ForceExt(fileName, "") == ForceExt(file, "") && 
				 PatternMatch(".out*", GetFileExt(fileName)) && 
				 PatternMatch(".out*", GetFileExt(file)))) 
				return iff;
		}
		return -1;
	}
	
	class DataSource : public Convert {
	public:
		DataSource() : datafast(0), col(0) {}
		void Init(FastOut &_datafast, int _col)	{datafast = &_datafast;	col = _col;};	
		Value Format(const Value& q) const {
			ASSERT(datafast);
			return datafast->dataOut[col][q];
		}
	private:
		FastOut *datafast;
		int col;
	};
	
	int scatterSize() {return scatter.size();}
	
	UArray<ScatterCtrl> scatter;
	UArray<UArray<DataSource>> dataSource;
	UArray<FastOut> dataFast;
	
private:
	Box box;
};

class FastScatterBase : public WithScatterBase<StaticRect> {
public:
	typedef FastScatterBase CLASSNAME;
	
	void Init(FastScatter *parent, Function <bool(String)> OnFile, Function <void(String)> OnCopyTabs, StatusBar &statusBar);
	void Clear();
	
	WithSearchColumn<StaticRect> leftSearch, rightSearch;
	
	void LoadParams();
	void SaveParams();
	
	void OnLoad();
	bool OnLoad(String filename);
	
	void SelPaste(String str = "");
	
	bool IsEmpty()		{return left.dataFast.IsEmpty();}
	
	ScatterLeft left;
	double timeStart = Null, timeEnd = Null;
	
private:
	FastScatter *fastScatter = nullptr;
	
	bool OnLoad0(String filename);
	
	void UpdateButtons(bool on);
	void OnSaveAs();
	void OnFilter(bool show);
	void ShowSelected(bool zoomtofit);
	bool AddParameter(String param, ArrayCtrl *parray);
	void WhenArrayLeftDouble(ArrayCtrl *parray);
	
	void OnDropInsert(int line, PasteClip& d, ArrayCtrl &array);
	//void OnDrop(PasteClip& d, ArrayCtrl &array);
	void OnDrag(ArrayCtrl &array, bool remove);
	
	String SelectedStr();
	void SelCopy();
	void SelCopyTabs();
	
	void OnTimer();
	TimeCallback timer;
	TimeStop timeStop;
	
	StatusBar *statusBar = nullptr;
	
	WithScatterRightTop<StaticRect> rightT;
	WithScatterRightBottom<StaticRect> rightB;
	WithScatterRightBottomUp<StaticRect> rightUp;
	WithScatterRightBottomDown<StaticRect> rightDown;
	Box right;
	
	RectEnterSet frameSet;
	
	Function <bool(String)> WhenFile;
	Function <void(String)> WhenCopyTabs;
	
	struct Params {
		UVector<String> left, right;
		void Jsonize(JsonIO &json) {
			json
				("left", left)
				("right", right)
			;
		}
		void Get(const ArrayCtrl &aleft, const ArrayCtrl &aright);
		void Set(FastScatterBase &parent, ArrayCtrl &aleft, ArrayCtrl &aright) const;
	};
	
	String saveFolder;
	
	UArray<EditString> edits;
	UArray<EditDouble> editNums;
	
	bool disableShowSelected = false;
};

class FastScatter : public StaticRect {
public:
	typedef FastScatter CLASSNAME;

	void Init(Function <bool(String)> OnFile, Function <void(String)> OnCopyTabs, StatusBar &statusBar);
	void OnCalc();
	bool OnLoadCompare();
	void OnSaveCompare();
	
	WithCompare<StaticRect> compare;
	FastScatterBase fscbase;
	SplitterButton splitCompare;
	
	ParameterMetrics params;
	
private:
	UArray<EditString> edits;
	
	void ParamsToGrid();
	void GridToParams();
};

class FastScatterTabs : public StaticRect {
public:
	typedef FastScatterTabs CLASSNAME;

	virtual ~FastScatterTabs();
	void Init(String appDataFolder, StatusBar &statusBar);
	void OnTab();
	void AddFile(String filename);
	void OnCloseTab(Value key);
	
	void Jsonize(JsonIO &json) {
		if (json.IsStoring()) {
			if (history.size() > 100)			// Saves 100 max
				history.Trim(100);
			if (historyMetrics.size() > 100)
				historyMetrics.Trim(100);
		}
		json
			("history", history)
			("historyMetrics", historyMetrics)
		;
	}
	
	TabBar tabBar;
	UArray<FastScatter> tabScatters;
	Upp::Index<int> tabKeys;
	
	bool loadingDragDrop = false;
	
private:
	int tabCounter = 0;
	
	Upp::Index<String> history;
	Upp::Index<String> historyMetrics;
	
	virtual void DragAndDrop(Point p, PasteClip& d);
	virtual bool Key(dword key, int count);
	bool LoadDragDrop(const UVector<String> &files);
	
	void AddHistory(String filename);
	void AddHistoryMetrics(String filename);
		
	TabBar::Style styleTab;
	
	StatusBar *statusBar = nullptr;
	
	friend FastScatter;
};


#endif
