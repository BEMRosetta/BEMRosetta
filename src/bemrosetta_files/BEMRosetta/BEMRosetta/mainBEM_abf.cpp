#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <CtrlScroll/CtrlScroll.h>

#include <BEMRosetta/BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;

#include "main.h"


void MainABForce::Init(DataToShow _dataToShow) {
	CtrlLayout(*this);
	
	dataToShow = _dataToShow;
	
	selTab = 0;
	isFilling = false;
	tab.WhenSet = [&] {
		LOGTAB(tab);
		if (!isFilling)
			selTab = tab.Get();
	};
}

void MainABForce::Clear() {
	tab.Reset();
	selTab = 0;
}

bool MainABForce::Load(BEMData &bem) {
	TempAssign<bool> _isFilling(isFilling, true);
	try {
		tab.Reset();
		Upp::Array<HydroClass> &hydro = bem.hydros; 
		if (hydro.IsEmpty()) 
			return false;
		String format;
		switch (dataToShow) {
		case DATA_A:		format = t_("A%s");			break;		
		case DATA_B:		format = t_("B%s");			break;
		case DATA_FORCE_SC:	format = t_("Fsc%s%.1fº");	break;
		case DATA_FORCE_FK:	format = t_("Ffk%s%.1fº");	break;
		case DATA_FORCE_EX:	format = t_("Fex%s%.1fº");	break;
		case DATA_RAO:		format = t_("RAO%s%.1fº");	break;
		case DATA_STS_MA:	NEVER();
		case DATA_STS_PH:	NEVER();
		}
		int sdof = 6*bem.Nb;
		if (dataToShow == DATA_A || dataToShow == DATA_B) {
			plots.SetCount(sdof);
			for (int i = 0; i < sdof; ++i) {
				plots[i].SetCount(sdof);
				for (int j = 0; j < sdof; ++j) {
					if (!bem.onlyDiagonal || i == j) {
						plots[i][j].Init(i, j, dataToShow);
						if (plots[i][j].Load(hydro)) {
							if (i != j)
								tab.Add(plots[i][j].SizePos(), Format(format, Hydro::StrBDOF(i, j)));
							else
								tab.Add(plots[i][j].SizePos(), Format(format, Hydro::StrBDOF(i)));
						}
					}
				}
			}
		} else {
			int Nh = bem.head.GetCount();
			if (Nh < 0) 
				return false;
			
			plots.SetCount(Nh);
			for (int ih = 0; ih < Nh; ++ih) 
				plots[ih].SetCount(sdof);
			for (int i = 0; i < sdof; ++i) {
				for (int ih = 0; ih < Nh; ++ih) {
					plots[ih][i].Init(i, bem.head[ih], dataToShow);
					if (plots[ih][i].Load(hydro))
						tab.Add(plots[ih][i].SizePos(), Format(format, Hydro::StrBDOFAbrev(i), bem.head[ih]));
				}
			}
		}
		
		if (tab.GetCount() == 0)
			return false;
		else if (tab.GetCount() > selTab)	
			tab.Set(selTab);
		return true;
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
		return false;
	}
}

void MainPlot::Init() {
	if (!isInit) {
		CtrlLayout(*this);
		isInit = true;
	}
	Clear();
}

void MainPlot::Init(int _idof, double jdof_ih, DataToShow _dataToShow) {
	MainPlot::Init();
	
	idof = _idof;
	jdof = int(jdof_ih);
	heading = jdof_ih;
	dataToShow = _dataToShow;
	
	scatter.ShowAllMenus();
	String title, labelY, labelY2;
	switch (dataToShow) {
	case DATA_A:		title = Format(t_("Added mass %s"), Hydro::StrBDOF(idof, jdof));		
						labelY = t_("Added mass");				
						break;		
	case DATA_B:		title = Format(t_("Radiation damping %s"), Hydro::StrBDOF(idof, jdof));
						labelY = t_("Radiation damping");		
						break;
	case DATA_FORCE_SC:	title = Format(t_("Diffraction scattering force %s heading %.1fº"), Hydro::StrBDOF(idof), heading);
						labelY = t_("Diffraction scattering force");
						labelY2 = t_("Diffraction scattering force phase [º]");	
						break;
	case DATA_FORCE_FK:	title = Format(t_("Froude-Krylov force %s heading %.1fº"), Hydro::StrBDOF(idof), heading);
						labelY = t_("Froude-Krylov force");		
						labelY2 = t_("Froude-Krylov force phase [º]");			
						break;
	case DATA_FORCE_EX:	title = Format(t_("Excitation Force %s heading %.1fº"), Hydro::StrBDOF(idof), heading);
						labelY = t_("Excitation force");		
						labelY2 = t_("Excitation force phase [º]");				
						break;
	case DATA_RAO:		title = Format(t_("Response Amplitude Operator %s heading %.1fº"), Hydro::StrBDOF(idof), heading);
						labelY = t_("RAO []");		
						labelY2 = t_("RAO phase [º]");							
						break;
	case DATA_STS_MA:	title = Format(t_("Z magnitude B(w)+jw(A(w)-Ainf) %s"), Hydro::StrBDOFFull(idof, jdof));		
						labelY = t_("Magnitude");				
						break;
	case DATA_STS_PH:	title = Format(t_("Z phase B(w)+jw(A(w)-Ainf) %s"), Hydro::StrBDOFFull(idof, jdof));		
						labelY = t_("Phase");				
						break;
	}
	scatter.SetTitle(title).SetTitleFont(Arial());
	scatter.SetLabelY(labelY);
	if (!labelY2.IsEmpty()) {
		scatter.SetDrawY2Reticle(true);
		scatter.SetPlotAreaRightMargin(80);
		scatter.SetLabelY2(labelY2);
	}
}

bool MainPlot::Load(const Upp::Array<HydroClass> &hydro) {
	ABF_source.SetCount(hydro.GetCount());
	ABF_source2.SetCount(hydro.GetCount());
	Ainf_source.SetCount(hydro.GetCount());
	
	dim = !mbm().menuPlot.showNdim;
	markW = mbm().menuPlot.showPoints ? 10 : 0;
	show_w = mbm().menuPlot.opwT == 0;
	if (show_w)
		scatter.SetLabelX(t_("w [rad/s]"));
	else
		scatter.SetLabelX(t_("T [s]"));
	
	bool loaded = false;
	for (int id = 0; id < hydro.GetCount(); ++id) {
		const Hydro &hy = hydro[id].hd();
		LoadEach(hy, id, loaded);
	}
	if (mbm().menuPlot.autoFit)
		scatter.ZoomToFit(true, true);
	return loaded;
}

bool MainPlot::Load(const Hydro &hy) {
	scatter.RemoveAllSeries();
	ABF_source.SetCount(1);
	ABF_source2.SetCount(1);
	Ainf_source.SetCount(1);
	
	dim = !mbm().menuPlot.showNdim;
	markW = mbm().menuPlot.showPoints ? 10 : 0;
	show_w = mbm().menuPlot.opwT == 0;
	if (show_w)
		scatter.SetLabelX(t_("w [rad/s]"));
	else
		scatter.SetLabelX(t_("T [s]"));
	
	bool loaded = false;
	LoadEach(hy, 0, loaded);
	
	if (mbm().menuPlot.autoFit)
		scatter.ZoomToFit(true, true);
	return loaded;
}

void MainPlot::LoadEach(const Hydro &hy, int id, bool &loaded) {
	int ih = -1;
	if (dataToShow != DATA_A && dataToShow != DATA_B) 
		ih = hy.GetHeadId(heading);
	String nameType = Format("%s(%s)", hy.name, hy.GetCodeStrAbr());
	if (dataToShow == DATA_A) {
		Upp::Color acolor = Null;
		if (hy.IsLoadedA()) {
			if (ABF_source[id].Init(hy, idof, jdof, PLOT_A, show_w, !dim)) {
				loaded = true;
				scatter.AddSeries(ABF_source[id]).Legend(Format(t_("A_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>();
				if (dim)
					scatter.Units("Ns2/m");
				double dummy;
				scatter.GetStroke(scatter.GetCount()-1, dummy, acolor);
			}
		}
		if (hy.IsLoadedAwinf()) {
			if (Ainf_source[id].Init(hy, idof, jdof, PLOT_AINF, show_w, !dim)) {
				loaded = true;
				scatter.AddSeries(Ainf_source[id]).Legend(Format(t_("Ainf_%s"), nameType)).Dash(LINE_DOTTED).Stroke(2, acolor).NoMark();
				if (dim)
					scatter.Units("Ns2/m");
			}
		}
	} else if (dataToShow == DATA_B && hy.IsLoadedB()) {
		if (ABF_source[id].Init(hy, idof, jdof, PLOT_B, show_w, !dim)) {
			loaded = true;
			scatter.AddSeries(ABF_source[id]).Legend(Format(t_("B_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>();
			if (dim)
				scatter.Units("Ns/m");
		}
	} else if (dataToShow == DATA_FORCE_SC && hy.IsLoadedFsc() && ih >= 0) {
		if (ABF_source[id].Init(hy, idof, ih, PLOT_FORCE_SC_MA, show_w, !dim)) {
			loaded = true;
			scatter.AddSeries(ABF_source[id]).Legend(Format(t_("Fsc_ma_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>();
			if (dim)
				scatter.Units("N");
			if (ABF_source2[id].Init(hy, idof, ih, PLOT_FORCE_SC_PH, show_w, !dim)) {
				loaded = true;
				if (mbm().menuPlot.showPhase)
					scatter.AddSeries(ABF_source2[id]).Legend(Format(t_("Fsc_ph_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetDataSecondaryY().Units("rad");
			}
		}
	} else if (dataToShow == DATA_FORCE_FK && hy.IsLoadedFfk() && ih >= 0) {
		if (ABF_source[id].Init(hy, idof, ih, PLOT_FORCE_FK_MA, show_w, !dim)) {
			loaded = true;
			scatter.AddSeries(ABF_source[id]).Legend(Format(t_("Ffk_ma_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>();
			if (dim)
				scatter.Units("N");
			if (ABF_source2[id].Init(hy, idof, ih, PLOT_FORCE_FK_PH, show_w, !dim)) {
				loaded = true;
				if (mbm().menuPlot.showPhase)
					scatter.AddSeries(ABF_source2[id]).Legend(Format(t_("Ffk_ph_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetDataSecondaryY().Units("rad");
			}
		}
	} else if (dataToShow == DATA_FORCE_EX && hy.IsLoadedFex() && ih >= 0) {
		if (ABF_source[id].Init(hy, idof, ih, PLOT_FORCE_EX_MA, show_w, !dim)) {
			loaded = true;
			scatter.AddSeries(ABF_source[id]).Legend(Format(t_("Fex_ma_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>();
			if (dim)
				scatter.Units("N");
			if (ABF_source2[id].Init(hy, idof, ih, PLOT_FORCE_EX_PH, show_w, !dim)) {
				loaded = true;
				if (mbm().menuPlot.showPhase)
					scatter.AddSeries(ABF_source2[id]).Legend(Format(t_("Fex_ph_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetDataSecondaryY().Units("rad");
			}
		}
	} else if (dataToShow == DATA_RAO && hy.IsLoadedRAO() && ih >= 0) {
		if (ABF_source[id].Init(hy, idof, ih, PLOT_RAO_MA, show_w, !dim)) {
			loaded = true;
			scatter.AddSeries(ABF_source[id]).Legend(Format(t_("RAO_ma_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>();
			if (ABF_source2[id].Init(hy, idof, ih, PLOT_RAO_PH, show_w, !dim)) {
				loaded = true;
				if (mbm().menuPlot.showPhase)
					scatter.AddSeries(ABF_source2[id]).Legend(Format(t_("RAO_ph_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetDataSecondaryY().Units("rad");
			}
		}
	} else if (dataToShow == DATA_STS_MA && hy.IsLoadedA() && hy.IsLoadedB()) {
		if (ABF_source[id].Init(hy, idof, jdof, PLOT_STS_MA, show_w, !dim)) {
			loaded = true;
			scatter.AddSeries(ABF_source[id]).Legend(Format(t_("Zmag_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>();
			if (dim)
				scatter.Units("Ns/m");
		}
	} else if (dataToShow == DATA_STS_PH && hy.IsLoadedA() && hy.IsLoadedB()) {
		if (ABF_source[id].Init(hy, idof, jdof, PLOT_STS_PH, show_w, !dim)) {
			loaded = true;
			scatter.AddSeries(ABF_source[id]).Legend(Format(t_("Zph_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>();
			scatter.Units("rad");
		}
	}
}

void MainPlot::Clear() {
	ABF_source.Clear();
	ABF_source2.Clear();
	Ainf_source.Clear();
	scatter.RemoveAllSeries();
}