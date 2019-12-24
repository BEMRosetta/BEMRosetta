#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <CtrlScroll/CtrlScroll.h>

#include <BEMRosetta/BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;

#include "main.h"


void MainPlot::Init(bool vert) {
	if (!isInit) {
		scatt.SetPlotAreaLeftMargin(90).SetPlotAreaRightMargin(70).SetPlotAreaBottomMargin(50)
			   .SetTitleFont(SansSerifZ(12)).ShowAllMenus();
		scatP.SetPlotAreaLeftMargin(90).SetPlotAreaRightMargin(70).SetPlotAreaBottomMargin(50)
			   .SetTitleFont(SansSerifZ(12)).ShowAllMenus();
		scatt.LinkedWith(scatP);
		if (vert)
			splitter.Vert(scatt.SizePos(), scatP.SizePos());
		else
			splitter.Horz(scatt.SizePos(), scatP.SizePos());
		Add(splitter.SizePos());
		scatt.ShowAllMenus();
		scatP.ShowAllMenus();
		isInit = true;
	}
	Clear();
}

void MainPlot::Init(int _idf, double jdf_ih, DataToShow _dataToShow) {
	MainPlot::Init(true);
	
	idf = _idf;
	jdf = int(jdf_ih);
	heading = jdf_ih;
	dataToShow = _dataToShow;
	
	String title, title2, labelY, labelY2;
	switch (dataToShow) {
	case DATA_A:		title = Format(t_("Added mass %s"), Hydro::StrBDOF(idf, jdf));		
						labelY = t_("Added mass");
						splitter.SetPos(10000, 0);				
						break;		
	case DATA_B:		title = Format(t_("Radiation damping %s"), Hydro::StrBDOF(idf, jdf));
						labelY = t_("Radiation damping");
						splitter.SetPos(10000, 0);	
						break;
	case DATA_FORCE_SC:	title = Format(t_("Diffraction scattering force %s heading %.1fº"), Hydro::StrBDOF(idf), heading);
						labelY = t_("Diffraction scattering force");
						labelY2 = t_("Diffraction scattering force phase [rad]");
						splitter.SetPos(5000, 0);	
						break;
	case DATA_FORCE_FK:	title = Format(t_("Froude-Krylov force %s heading %.1fº"), Hydro::StrBDOF(idf), heading);
						labelY = t_("Froude-Krylov force");		
						labelY2 = t_("Froude-Krylov force phase [rad]");
						splitter.SetPos(5000, 0);			
						break;
	case DATA_FORCE_EX:	title = Format(t_("Excitation Force %s heading %.1fº"), Hydro::StrBDOF(idf), heading);
						labelY = t_("Excitation force");		
						labelY2 = t_("Excitation force phase [rad]");
						splitter.SetPos(5000, 0);				
						break;
	case DATA_RAO:		title = Format(t_("Response Amplitude Operator %s heading %.1fº"), Hydro::StrBDOF(idf), heading);
						labelY = t_("RAO []");		
						labelY2 = t_("RAO phase [rad]");
						splitter.SetPos(5000, 0);							
						break;
	case DATA_STS:		title = Format(t_("Magnitude Z = B(w)+jw(A(w)-Ainf) %s"), Hydro::StrBDOFFull(idf, jdf));		
						labelY = t_("Magnitude");				
						title2 = Format(t_("Phase Z = B(w)+jw(A(w)-Ainf) %s"), Hydro::StrBDOFFull(idf, jdf));		
						labelY2 = t_("Phase");
						splitter.SetPos(5000, 0);				
						break;
	case DATA_STS2:		title = Format(t_("Magnitude response %s"), Hydro::StrBDOF(idf, jdf));
						title2 = Format(t_("Frequency response %s"), Hydro::StrBDOF(idf, jdf));
						splitter.SetPos(5000, 0);
	}
	scatt.SetTitle(title);
	scatt.SetLabelY(labelY);
	scatP.SetTitle(title2);
	scatP.SetLabelY(labelY2);
}

bool MainPlot::Load(const Upp::Array<HydroClass> &hydro, const MainBEM &mbm, const Vector<int> &ids) {
	scatt.RemoveAllSeries();
	scatP.RemoveAllSeries();
	ABFZ_source.SetCount(ids.GetCount());
	ABFZ_source2.SetCount(ids.GetCount());
	Ainf_source.SetCount(ids.GetCount());
	A0_source.SetCount(ids.GetCount());
	TFS_source.SetCount(ids.GetCount());
	TFS_source2.SetCount(ids.GetCount());
		
	dim = !mbm.menuPlot.showNdim;
	markW = mbm.menuPlot.showPoints ? 10 : 0;
	show_w = mbm.menuPlot.opwT == 0;
	if (show_w) {
		scatt.SetLabelX(t_("w [rad/s]"));
		scatP.SetLabelX(t_("w [rad/s]"));
	} else {
		scatt.SetLabelX(t_("T [s]"));
		scatP.SetLabelX(t_("T [s]"));
	}
	bool loaded = false;
	for (int id = 0; id < ids.GetCount(); ++id) {
		const Hydro &hy = hydro[ids[id]].hd();
		LoadEach(hy, id, loaded);
	}
	if (mbm.menuPlot.autoFit) {
		scatt.ZoomToFit(true, true);
		if (mbm.menuPlot.fromY0) {
			double yRange = max<double>(0, scatt.GetYMin()) + scatt.GetYRange();
			scatt.SetXYMin(Null, 0).SetRange(Null, yRange);
		}
		scatP.ZoomToFit(true, false);
		scatP.SetXYMin(Null, -M_PI).SetRange(Null, 2*M_PI).SetMajorUnits(Null, 1);
		scatP.SetMinUnits(Null, M_PI-3);
	}
	return loaded;
}

bool MainPlot::Load(const Hydro &hy, const MainBEM &mbm) {
	scatt.RemoveAllSeries();
	scatP.RemoveAllSeries();
	ABFZ_source.SetCount(1);
	ABFZ_source2.SetCount(1);
	Ainf_source.SetCount(1);
	A0_source.SetCount(1);
	
	dim = !mbm.menuPlot.showNdim;
	markW = mbm.menuPlot.showPoints ? 10 : 0;
	show_w = mbm.menuPlot.opwT == 0;
	if (show_w) {
		scatt.SetLabelX(t_("w [rad/s]"));
		scatP.SetLabelX(t_("w [rad/s]"));
	}else {
		scatt.SetLabelX(t_("T [s]"));
		scatP.SetLabelX(t_("T [s]"));
	}
	bool loaded = false;
	int idc = ArrayModel_IdHydro(mbm.menuFOAMM.arrayModel);
	LoadEach(hy, 0, loaded, idc);
	
	if (mbm.menuPlot.autoFit) {
		scatt.ZoomToFit(true, true);
		if (mbm.menuPlot.fromY0) {
			double ymax = scatt.GetYMin() + scatt.GetYRange();
			scatt.SetXYMin(Null, 0).SetRange(Null, ymax);
		}
		scatP.ZoomToFit(true, false);
		scatP.SetXYMin(Null, -M_PI).SetRange(Null, 2*M_PI).SetMajorUnits(Null, 1);
		scatP.SetMinUnits(Null, M_PI-3);
	}
	return loaded;
}

void MainPlot::LoadEach(const Hydro &hy, int id, bool &loaded, int idc) {
	int ih = -1;
	if (dataToShow != DATA_A && dataToShow != DATA_B) 
		ih = hy.GetHeadId(heading);
	String nameType = Format("%s(%s)", hy.name, hy.GetCodeStrAbr());
	if (idc < 0)
		idc = id;
	const Upp::Color &color = GetColorId(idc);
	if (dataToShow == DATA_A) {
		Upp::Color acolor = Null;
		if (hy.IsLoadedA()) {
			if (ABFZ_source[id].Init(hy, idf, jdf, PLOT_A, show_w, !dim)) {
				loaded = true;
				scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("A_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().Stroke(2, color);
				if (dim)
					scatt.Units(t_("Ns2/m"));
				double dummy;
				scatt.GetStroke(scatt.GetCount()-1, dummy, acolor);
			}
		}
		if (hy.IsLoadedAwinf()) {
			if (Ainf_source[id].Init(hy, idf, jdf, PLOT_AINF, show_w, !dim)) {
				loaded = true;
				scatt.AddSeries(Ainf_source[id]).Legend(Format(t_("Ainf_%s"), nameType)).Dash(LINE_DOTTED).SetMarkColor(acolor).MarkStyle<SquareMarkPlot>().Stroke(2, color);
				scatt.Stroke(show_w ? 2 : 0, acolor).SetMarkWidth(show_w ? 0 : 8);
				if (dim)
					scatt.Units(t_("Ns2/m"));
			}
		}
		if (hy.IsLoadedAw0()) {
			if (A0_source[id].Init(hy, idf, jdf, PLOT_A0, show_w, !dim)) {
				loaded = true;
				scatt.AddSeries(A0_source[id]).Legend(Format(t_("A0_%s"), nameType)).Dash(LINE_DOTTED).SetMarkColor(acolor).MarkStyle<SquareMarkPlot>().Stroke(2, color);
				scatt.Stroke(!show_w ? 2 : 0, acolor).SetMarkWidth(!show_w ? 0 : 8);
				if (dim)
					scatt.Units(t_("Ns2/m"));
			}
		}
	} else if (dataToShow == DATA_B && hy.IsLoadedB()) {
		if (ABFZ_source[id].Init(hy, idf, jdf, PLOT_B, show_w, !dim)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("B_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().Stroke(2, color);
			if (dim)
				scatt.Units(t_("Ns/m"));
		}
	} else if (dataToShow == DATA_FORCE_SC && hy.IsLoadedFsc() && ih >= 0) {
		if (ABFZ_source[id].Init(hy, idf, ih, PLOT_FORCE_SC_MA, show_w, !dim)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("Fsc_ma_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().Stroke(2, color);
			if (dim)
				scatt.Units(t_("N"));
			if (ABFZ_source2[id].Init(hy, idf, ih, PLOT_FORCE_SC_PH, show_w, !dim)) {
				loaded = true;
				scatP.AddSeries(ABFZ_source2[id]).Legend(Format(t_("Fsc_ph_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().Units("rad").Stroke(2, color);
			}
		}
	} else if (dataToShow == DATA_FORCE_FK && hy.IsLoadedFfk() && ih >= 0) {
		if (ABFZ_source[id].Init(hy, idf, ih, PLOT_FORCE_FK_MA, show_w, !dim)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("Ffk_ma_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().Stroke(2, color);
			if (dim)
				scatt.Units(t_("N"));
			if (ABFZ_source2[id].Init(hy, idf, ih, PLOT_FORCE_FK_PH, show_w, !dim)) {
				loaded = true;
				scatP.AddSeries(ABFZ_source2[id]).Legend(Format(t_("Ffk_ph_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().Units("rad").Stroke(2, color);
			}
		}
	} else if (dataToShow == DATA_FORCE_EX && hy.IsLoadedFex() && ih >= 0) {
		if (ABFZ_source[id].Init(hy, idf, ih, PLOT_FORCE_EX_MA, show_w, !dim)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("Fex_ma_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().Stroke(2, color);
			if (dim)
				scatt.Units(t_("N"));
			if (ABFZ_source2[id].Init(hy, idf, ih, PLOT_FORCE_EX_PH, show_w, !dim)) {
				loaded = true;
				scatP.AddSeries(ABFZ_source2[id]).Legend(Format(t_("Fex_ph_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().Units("rad").Stroke(2, color);
			}
		}
	} else if (dataToShow == DATA_RAO && hy.IsLoadedRAO() && ih >= 0) {
		if (ABFZ_source[id].Init(hy, idf, ih, PLOT_RAO_MA, show_w, !dim)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("RAO_ma_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().Stroke(2, color);
			if (ABFZ_source2[id].Init(hy, idf, ih, PLOT_RAO_PH, show_w, !dim)) {
				loaded = true;
				scatP.AddSeries(ABFZ_source2[id]).Legend(Format(t_("RAO_ph_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().Units("rad").Stroke(2, color);
			}
		}
	} else if (dataToShow == DATA_STS && hy.IsLoadedA() && hy.IsLoadedB()) {
		if (ABFZ_source[id].Init(hy, idf, jdf, PLOT_Z_MA, show_w, !dim)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("Zmag_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().Stroke(2, color);
			if (dim)
				scatt.Units(t_("Ns/m"));
			if (ABFZ_source2[id].Init(hy, idf, jdf, PLOT_Z_PH, show_w, !dim)) {
				loaded = true;
				scatP.AddSeries(ABFZ_source2[id]).Legend(Format(t_("Zph_%s"), nameType)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().Units("rad").Stroke(2, color);
			}
		}
	} else if (dataToShow == DATA_STS2 && hy.IsLoadedStateSpace()) {
		if (ABFZ_source[id].Init(hy, idf, jdf, PLOT_Z_MA, show_w, !dim)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("Zmag %s"), hy.name)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().Stroke(2, color);//.Units("dB");
			if (ABFZ_source2[id].Init(hy, idf, jdf, PLOT_Z_PH, show_w, !dim)) {
				loaded = true;
				scatP.AddSeries(ABFZ_source2[id]).Legend(Format(t_("Zph %s"), hy.name)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().Units("rad").Stroke(2, color);
			}
		}
		if (TFS_source[id].Init(hy, idf, jdf, PLOT_TFS_MA, show_w, !dim)) {
			loaded = true;
			scatt.AddSeries(TFS_source[id]).Legend(Format(t_("TFSmag %s"), hy.name)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().Stroke(4, color).Dash(LINE_DASH_DOT);//.Units("dB");
			if (TFS_source2[id].Init(hy, idf, jdf, PLOT_TFS_PH, show_w, !dim)) {
				loaded = true;
				scatP.AddSeries(TFS_source2[id]).Legend(Format(t_("TFSph %s"), hy.name)).SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().Units("rad").Stroke(4, color).Dash(LINE_DASH_DOT);
			}
		}
	}
}

void MainPlot::Clear() {
	ABFZ_source.Clear();
	ABFZ_source2.Clear();
	Ainf_source.Clear();
	A0_source.Clear();
	scatt.RemoveAllSeries();
	scatP.RemoveAllSeries();
}