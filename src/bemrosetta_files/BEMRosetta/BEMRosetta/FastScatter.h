#ifndef _BEMRosetta_BEMRosetta_FastScatter_h_
#define _BEMRosetta_BEMRosetta_FastScatter_h_

template <class T>
class WithClick : public T {
private:
	virtual Image CursorImage(Point, dword) 	{return Image::Hand();}
	virtual void LeftDown(Point, dword) 		{T::WhenAction();}	
};

#define LAYOUTFILE <BEMRosetta/BEMRosetta/FastScatter.lay>
#include <CtrlCore/lay.h>

class FastScatter : public WithFastScatter<StaticRect> {
public:
	typedef FastScatter CLASSNAME;
	
	void Init();
		
private:
	bool OnLoad();
	void OnFilter();
	void ShowSelected();
		
	void OnTimer();
	TimeCallback timer;
	
	StatusBar statusBar;
	
	FastOut datafast;
	
	/*class DataSource : public Convert {
	public:
		DataSource() : datafast(0), col(0) {}
		void Init(FastOut &datafast, int col)	{this->datafast = &datafast;	this->col = col;};	
		Value Format(const Value& q) const;
	private:
		FastOut *datafast;
		int col;
	};
	Upp::Array<DataSource> dataSource;*/
	
	WithFastScatterLeft<StaticRect> left;
	WithFastScatterRight<StaticRect> right;
	WithSearchColumn<StaticRect> leftSearch, rightSearch;
};
	
class FastScatterTabs : public StaticRect {
public:
	typedef FastScatterTabs CLASSNAME;

	void Init();
	void OnTab();
	void AddTab(String filename);
	void OnCloseTab(Value key);
	
private:
	TabBar tabBar;
	Upp::Array<FastScatter> tabScatters;
	Upp::Index<int> tabKeys;
	int tabCounter = 0;
	
	virtual void DragAndDrop(Point p, PasteClip& d);
	virtual bool Key(dword key, int count);
	bool LoadDragDrop(const Vector<String> &files);
};


#endif
