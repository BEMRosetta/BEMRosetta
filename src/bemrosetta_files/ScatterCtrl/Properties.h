#ifndef _ScatterCtrl_Properties_h_
#define _ScatterCtrl_Properties_h_

#define  LAYOUTFILE <ScatterCtrl/ScatterCtrl.lay>
#include <CtrlCore/lay.h>

#define IMAGECLASS ScatterImgP
#define IMAGEFILE <ScatterCtrl/ScatterCtrl.iml>
#include <Draw/iml_header.h>

class FontSelect : public WithFontSelector<TopWindow> {
typedef FontSelect CLASSNAME;
public:
	Event<Font> WhenAction;

	FontSelect() {
		CtrlLayoutExit(*this, t_("Font"));
		FrameLess(true);
		
		face.WhenAction = THISBACK(Select);
		height.WhenAction = THISBACK(Select);
		bold.WhenAction = THISBACK(Select);
		italic.WhenAction = THISBACK(Select);
		naa.WhenAction = THISBACK(Select);
		face.Clear();
		Upp::Index<String> fni;
		for(int i = 0; i < Font::GetFaceCount(); i++) {
			if(Font::GetFaceInfo(i)) {
				String n = Font::GetFaceName(i);
				if(fni.Find(n) < 0) {
					fni.Add(n);
					face.Add(i, n);
				}
			}
		}
		face.SetIndex(0);
		height.ClearList();
		for(int i = 6; i < 64; i++)
			height.Add(i);
		Select();
	}
	void Set(Font f) {
		int fi = f.GetFace();
		if(!face.HasKey(fi)) {
			fi = face.FindValue(f.GetFaceName());
			if(fi < 0)
				fi = Font::COURIER;
			else
				fi = face.GetKey(fi);
		}
		face.SetData(fi);
		Select();
		height.SetData(f.GetHeight());
		for(int i = 0; i < height.GetCount(); i++) {
			int q = height.GetKey(i);
			if(f.GetHeight() <= q) {
				height.SetData(q);
				break;
			}
		}
		bold = f.IsBold();
		italic = f.IsItalic();
		naa = f.IsNonAntiAliased();
	}
	Font Get() {
		Font f(face.GetData(), height.GetData());
		f.Bold(bold);
		f.Italic(italic);
		f.NonAntiAliased(naa);
		return f;
	}
	void Execute(Ctrl &_parent) {
		Open(this);
		Rect rec = GetRect();
		Size sz = rec.GetSize();
		rec.left = _parent.GetScreenRect().left;
		rec.top = _parent.GetScreenRect().bottom;
		rec.right = rec.left + sz.cx;
		rec.bottom = rec.top + sz.cy;
		SetRect(rec);
		Run();
	}
	
private:
	void Select() {WhenAction(Get());}
};	


class MeasuresTab : public WithMeasures<StaticRect> {
public:
	typedef MeasuresTab CLASSNAME;
	
	void Init(ScatterCtrl &scatter);
	void Change();

private:
	ScatterCtrl *pscatter;
};

class TextsTab : public WithTexts<StaticRect> {
public:
	typedef TextsTab CLASSNAME;
	
	void Init(ScatterCtrl &scatter);
	void DoShowText();
	
private:
	ScatterCtrl *pscatter;
	
	void Change();
	void OnFontTitle();
	void OnChangeFontTitle(Font font);
	void OnFontLabel();
	void OnChangeFontLabel(Font font);
};

class LegendTab : public WithLegend<StaticRect> {
public:
	typedef LegendTab CLASSNAME;
	
	void Init(ScatterCtrl &scatter);

private:
	ScatterCtrl *pscatter;
	
	void Change();
	void ChangeAnchor(Option *op);
};

class SeriesTab : public Splitter {
public:
	typedef SeriesTab CLASSNAME;
	
	SeriesTab() : dashCount(DashStyle::GetCount()) {}
	virtual ~SeriesTab() {
		DashStyle::UnregisterFrom(dashCount);
	}
	void Init(ScatterCtrl& scatter);
	
private:
	ScatterCtrl *pscatter;
	int dashCount;
	
	void Change();
	void ChangeMark();
	void UpdateFields();
	void OnMoveUp(); 
	void OnMoveDown();
	void OnDelete();
	
	void Init0();
	
	WithSeriesLeft<StaticRect> left;
	WithSeriesRight<StaticRect> right;
};

class GeneralTab : public WithGeneral<StaticRect> {
public:
	typedef GeneralTab CLASSNAME;
	
	void Init(ScatterCtrl &scatter);

private:
	ScatterCtrl *pscatter;
	
	void Change();
	void ChangeAll();
};

class DataDlg : public WithData<TopWindow> {
public:
	typedef DataDlg CLASSNAME;
	
	void Init(ScatterCtrl& scatter);
	virtual ~DataDlg() {};
	
	void OnTab();
	void OnArrayBar(Bar &menu);
	void ArrayCopy();
	void ArraySelect();
	void ArraySaveToFile(String fileName);
	
	class DataSourceX : public Convert {
	public:
		Value Format(const Value& q) const;
		ScatterDraw *pscatter;
		int index;
	} dataSourceX;
	
	class DataSourceY : public Convert {
	public:
		Value Format(const Value& q) const;
		ScatterDraw *pscatter;
		int index;
	};
	Upp::Array<DataSourceY> dataSourceYArr;
	
private:
	ScatterCtrl *pscatter;
	
	Upp::Array <WithDataSeries <StaticRect> > series;
};

class PropertiesDlg : public WithProperties<TopWindow> {
public:
	typedef PropertiesDlg CLASSNAME;
	
	void Init(ScatterCtrl& scatter);
	virtual ~PropertiesDlg() {};
	
	void Set(int tab);
	void OnTab();
	void Perform();
	virtual void Close() {	
		if (pscatter->GetCount() == 0)
			RejectBreak(IDOK);
		TopWindow::Close();
	}
		
private:
	ScatterCtrl* pscatter;
	MeasuresTab measures;
	TextsTab texts;
	LegendTab legend;
	SeriesTab series;
	GeneralTab general;
};


class ProcessingTab : public WithProcessingTab<StaticRect> {
public:
	typedef ProcessingTab CLASSNAME;

	ProcessingTab();
	virtual ~ProcessingTab() {};
	
	void Init(ScatterCtrl& scatter) {pscatter = &scatter;}
	void UpdateField(const String name, int id);
	void OnFFT();
	void OnOp();
	void OnAutoSensSector();
	void OnShowEquation();
	void UpdateEquations();
	void OnFromTo();
	void OnOperation();
	void OnSet();
	void OnUpdateSensitivity();
	void OnFit();
		
private:
	ScatterCtrl* pscatter;
	int id;
	
	WithProcessingTabFitLeft<StaticRect> tabFitLeft;
	WithProcessingTabFitRight<StaticRect> tabFitRight;
	Splitter splitterTabFit;
	WithProcessingTabFrequencyLeft<StaticRect> tabFreqLeft;
	WithProcessingTabFrequencyRight<StaticRect> tabFreqRight;
	Splitter splitterTabFreq;
	WithProcessingTabOpLeft<StaticRect> tabOpLeft;
	WithProcessingTabOpRight<StaticRect> tabOpRight;
	Splitter splitterTabOp;
	WithProcessingTabBestFitLeft<StaticRect> tabBestFitLeft;
	WithProcessingTabBestFitRight<StaticRect> tabBestFitRight;
	Splitter splitterTabBestFit;
	
	void ArrayCopy();
	void ArraySelect();
	void OnArrayBar(Bar &menu);
	Upp::Array<ExplicitEquation> equationTypes;
	UserEquation *userEquation;
	//GridCtrlSource ds;
	
	bool avgFirst, linearFirst, cuadraticFirst, cubicFirst, sinusFirst, sinusTendFirst, splineFirst;
	double r2Linear, r2Cuadratic, r2Cubic, r2Sinus, r2SinusTend;
	bool tabFreqFirst, tabOpFirst, tabBestFitFirst;
	
	Vector<Pointf> fft;
	AvgEquation average;
	LinearEquation linear;
	PolynomialEquation2 cuadratic;
	PolynomialEquation3 cubic;
	SinEquation sinus, sinusTend;
	SplineEquation spline;
	Vector<Pointf> upperEnvelope, lowerEnvelope;
	Vector<Pointf> movAvg, secAvg, cumAvg;
	DataXRange dataXRange;
	bool exclamationOpened;
	double newWidthMax, newWidthMin, newWidthMovAvg;
};

class ProcessingDlg : public TopWindow {
public:
	typedef ProcessingDlg CLASSNAME;

	void Init(ScatterCtrl& scatter);
	virtual ~ProcessingDlg() {};

private:
	ScatterCtrl* pscatter;
	Upp::Array<ProcessingTab> tabs;
	WithProcessingLeft<StaticRect> list; 
	WithProcessingRight<StaticRect> right;
	Splitter splitter;
	
	void UpdateFields();
};

#endif
