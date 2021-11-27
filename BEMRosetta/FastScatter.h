// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2021, the BEMRosetta author and contributors
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

class FastScatter : public WithFastScatter<StaticRect> {
public:
	typedef FastScatter CLASSNAME;
	
	void Init(Function <void(String)> OnFile, Function <void(String)> OnCopyTabs, StatusBar &statusBar);
	void Clear();
	
	WithSearchColumn<StaticRect> leftSearch, rightSearch;
	
	void LoadParams();
	void SaveParams();
	
	void SelPaste(String str = "");
	
private:
	bool OnLoad(bool justUpdate = false);
	void OnFilter(bool show);
	void ShowSelected();
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
	
	FastOut datafast;
	
	class DataSource : public Convert {
	public:
		DataSource() : datafast(0), col(0) {}
		void Init(FastOut &_datafast, int _col)	{datafast = &_datafast;	col = _col;};	
		Value Format(const Value& q) const;
	private:
		FastOut *datafast;
		int col;
	};
	Upp::Array<DataSource> dataSource;
	
	WithFastScatterLeft<StaticRect> left;
	WithFastScatterRight<StaticRect> right;
	
	RectEnterSet frameSet;
	
	Function <void(String)> WhenFile;
	Function <void(String)> WhenCopyTabs;
	
	struct Params {
		Vector<String> left, right;
		void Jsonize(JsonIO &json) {
			json
				("left", left)
				("right", right)
			;
		}
		void Get(const ArrayCtrl &aleft, const ArrayCtrl &aright);
		void Set(ArrayCtrl &aleft, ArrayCtrl &aright) const;
	};
};
	
class FastScatterTabs : public StaticRect {
public:
	typedef FastScatterTabs CLASSNAME;

	virtual ~FastScatterTabs();
	void Init(String appDataFolder, StatusBar &statusBar);
	void OnTab();
	void AddTab(String filename);
	void OnCloseTab(Value key);
	
	void Jsonize(JsonIO &json) {
		json
			("history", history)
		;
	}
	
private:
	TabBar tabBar;
	Upp::Array<FastScatter> tabScatters;
	Upp::Index<int> tabKeys;
	Vector<String> fileNames;
	int tabCounter = 0;
	
	Index<String> history;
	
	virtual void DragAndDrop(Point p, PasteClip& d);
	virtual bool Key(dword key, int count);
	bool LoadDragDrop(const Upp::Vector<String> &files);
	
	void AddHistory(String filename);
		
	TabBar::Style styleTab;
	
	StatusBar *statusBar = nullptr;
};

class MainFASTW : public TopWindow {
public:
	typedef MainFASTW CLASSNAME;
	
	void Init(String appDataFolder, const Image &icon, const Image &largeIcon, StatusBar &statusBar);
	virtual void Close() 	{delete this;}
	
	FastScatterTabs fast;
};

#endif
