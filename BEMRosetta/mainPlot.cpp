#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <TabBar/TabBar.h>

#include <BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;

#include "main.h"


void MainPlot::Init(bool vert) {
	if (!isInit) {
		int len = StdFont().GetHeight();
		
		scatt.SetPlotAreaLeftMargin(8*len).SetPlotAreaRightMargin(len).SetPlotAreaBottomMargin(4*len)
			   .SetTitleFont(SansSerifZ(12)).ShowAllMenus();
		scatP.SetPlotAreaLeftMargin(8*len).SetPlotAreaRightMargin(len).SetPlotAreaBottomMargin(4*len)
			   .SetTitleFont(SansSerifZ(12)).ShowAllMenus();
		scatt.LinkedWith(scatP);
		
		compare.Init(scatt);
		splitCompare.Horz(scatt.SizePos(), compare.SizePos());
		splitCompare.SetPositions(7000, 10000).SetInitialPositionId(1).SetButtonNumber(1).SetButtonWidth(len);
		if (vert)
			splitter.Vert(splitCompare.SizePos(), scatP.SizePos());
		else
			splitter.Horz(splitCompare.SizePos(), scatP.SizePos());
		Add(splitter.SizePos());
		scatt.ShowAllMenus();
		scatP.ShowAllMenus();
		isInit = true;
	}
	Clear();
}

void MainPlot::Init(int _idf, double jdf_ih, Hydro::DataToShow _dataToShow) {
	MainPlot::Init(true);
	
	plot_idf = _idf;
	plot_jdf = int(jdf_ih);
	heading = jdf_ih;
	dataToShow = _dataToShow;
	
	compare.Init(dataToShow);
	
	String title, title2, labelY, labelY2;
	switch (dataToShow) {
	case Hydro::DATA_A:	title = Format(t_("Added mass %s"), Hydro::StrBDOF(plot_idf, plot_jdf));		
						labelY = t_("Added mass");
						splitter.SetPos(10000, 0);				
						break;		
	case Hydro::DATA_B:	title = Format(t_("Radiation damping %s"), Hydro::StrBDOF(plot_idf, plot_jdf));
						labelY = t_("Radiation damping");
						splitter.SetPos(10000, 0);	
						break;
	case Hydro::DATA_AINFW:	title = Format(t_("Added mass at infinity (ω) %s"), Hydro::StrBDOF(plot_idf, plot_jdf));		
						labelY = t_("Added mass at infinity (ω)");
						splitter.SetPos(10000, 0);				
						break;
	case Hydro::DATA_K:	title = Format(t_("Kirf Impulse Response Function %s"), Hydro::StrBDOF(plot_idf, plot_jdf));
						labelY = t_("Kirf Impulse Response Function");
						splitter.SetPos(10000, 0);	
						break;
	case Hydro::DATA_FORCE_SC:	title = Format(t_("Diffraction scattering force %s heading %.1fº"), Hydro::StrBDOF(plot_idf), heading);
						labelY = t_("Diffraction scattering force");
						labelY2 = t_("Diffraction scattering force phase [rad]");
						splitter.SetPos(5000, 0);	
						break;
	case Hydro::DATA_FORCE_FK:	title = Format(t_("Froude-Krylov force %s heading %.1fº"), Hydro::StrBDOF(plot_idf), heading);
						labelY = t_("Froude-Krylov force");		
						labelY2 = t_("Froude-Krylov force phase [rad]");
						splitter.SetPos(5000, 0);			
						break;
	case Hydro::DATA_FORCE_EX:	title = Format(t_("Excitation Force %s heading %.1fº"), Hydro::StrBDOF(plot_idf), heading);
						labelY = t_("Excitation force");		
						labelY2 = t_("Excitation force phase [rad]");
						splitter.SetPos(5000, 0);				
						break;
	case Hydro::DATA_RAO:	title = Format(t_("Response Amplitude Operator %s heading %.1fº"), Hydro::StrBDOF(plot_idf), heading);
						labelY = t_("RAO []");		
						labelY2 = t_("RAO phase [rad]");
						splitter.SetPos(5000, 0);							
						break;
	case Hydro::DATA_STS:	title = Format(t_("Magnitude Kr(ω) = B(ω)+jω{A(ω)-A∞} %s"), Hydro::StrBDOFFull(plot_idf, plot_jdf));		
						labelY = t_("Magnitude");				
						title2 = Format(t_("Phase Kr(ω) = B(ω)+jω{A(ω)-A∞} %s"), Hydro::StrBDOFFull(plot_idf, plot_jdf));		
						labelY2 = t_("Phase");
						splitter.SetPos(5000, 0);				
						break;
	case Hydro::DATA_STS2:	title = Format(t_("Magnitude response %s"), Hydro::StrBDOF(plot_idf, plot_jdf));
						title2 = Format(t_("Frequency response %s"), Hydro::StrBDOF(plot_idf, plot_jdf));
						splitter.SetPos(5000, 0);
	}
	scatt.SetTitle(title);
	scatt.SetLabelY(labelY);
	scatP.SetTitle(title2);
	scatP.SetLabelY(labelY2);
}

bool MainPlot::Load(const Upp::Array<HydroClass> &hydro, const MainBEM &mbm, const Upp::Vector<int> &ids) {
	scatt.RemoveAllSeries();
	scatP.RemoveAllSeries();
	ABFZ_source.SetCount(ids.size());
	ABFZ_source2.SetCount(ids.size());
	Ainf_source.SetCount(ids.size());
	A0_source.SetCount(ids.size());
	TFS_source.SetCount(ids.size());
	TFS_source2.SetCount(ids.size());
		
	dim = !mbm.menuPlot.showNdim;
	markW = mbm.menuPlot.showPoints ? 10 : 0;
	show_w = mbm.menuPlot.opwT == 0;
	if (show_w && dataToShow != Hydro::DATA_K) {
		scatt.SetLabelX(t_("ω [rad/s]"));
		scatP.SetLabelX(t_("ω [rad/s]"));
	} else {
		scatt.SetLabelX(t_("T [s]"));
		scatP.SetLabelX(t_("T [s]"));
	}
	bool loaded = false;
	for (int id = 0; id < ids.size(); ++id) {
		const Hydro &hy = hydro[ids[id]].hd();
		LoadEach(hy, id, loaded, hy.GetId());
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
	compare.Load();
	
	return loaded;
}

bool MainPlot::Load(const Hydro &hy, const MainBEM &mainBem) {
	scatt.RemoveAllSeries();
	scatP.RemoveAllSeries();
	ABFZ_source.SetCount(1);
	ABFZ_source2.SetCount(1);
	Ainf_source.SetCount(1);
	A0_source.SetCount(1);
	
	dim = !mainBem.menuPlot.showNdim;
	markW = mainBem.menuPlot.showPoints ? 10 : 0;
	show_w = mainBem.menuPlot.opwT == 0;
	if (show_w && dataToShow != Hydro::DATA_K) {
		scatt.SetLabelX(t_("ω [rad/s]"));
		scatP.SetLabelX(t_("ω [rad/s]"));
	} else {
		scatt.SetLabelX(t_("T [s]"));
		scatP.SetLabelX(t_("T [s]"));
	}
	bool loaded = false;
	int idc = ArrayModel_Id(mainBem.listLoaded);
	LoadEach(hy, 0, loaded, idc);
	
	if (mainBem.menuPlot.autoFit) {
		scatt.ZoomToFit(true, true);
		if (mainBem.menuPlot.fromY0) {
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
	if (dataToShow != Hydro::DATA_A && dataToShow != Hydro::DATA_B && dataToShow != Hydro::DATA_AINFW) 
		ih = hy.GetHeadId(heading);
	String nameType = Format("%s(%s)", hy.name, hy.GetCodeStrAbr());
	if (idc < 0)
		idc = id;
	const Upp::Color &color = GetColorId(idc);
	if (dataToShow == Hydro::DATA_A) {
		if (hy.IsLoadedA()) {
			if (ABFZ_source[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_A, show_w, !dim)) {
				loaded = true;
				scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("A_%s"), nameType)).
						SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
				if (dim)
					scatt.Units(t_("Ns2/m"));
			}
		}
		if (hy.IsLoadedAwinf()) {
			if (Ainf_source[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_AINF, show_w, !dim)) {
				loaded = true;
				scatt.AddSeries(Ainf_source[id]).Legend(Format(t_("A∞_%s"), nameType)).
						MarkStyle<SquareMarkPlot>().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_DOTTED);
				scatt.Stroke(show_w ? 2 : 0, color).SetMarkWidth(show_w ? 0 : 8);
				if (dim)
					scatt.Units(t_("Ns2/m"));
			}
		}
		if (hy.IsLoadedAw0()) {
			if (A0_source[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_A0, show_w, !dim)) {
				loaded = true;
				scatt.AddSeries(A0_source[id]).Legend(Format(t_("A0_%s"), nameType)).
						SetMarkWidth(!show_w ? 0 : 8).SetMarkColor(color).MarkStyle<SquareMarkPlot>().
						Stroke(2, color).Dash(LINE_DOTTED);
				if (dim)
					scatt.Units(t_("Ns2/m"));
			}
		}
	} else if (dataToShow == Hydro::DATA_AINFW) {
		if (hy.IsLoadedAINFW()) {
			if (ABFZ_source[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_AINFW, show_w, !dim)) {
				loaded = true;
				scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("A∞(ω)_%s"), nameType)).
						SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
				if (dim)
					scatt.Units(t_("Ns2/m"));
			}
		}
		if (hy.IsLoadedAwinf()) {
			if (Ainf_source[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_AINF, show_w, !dim)) {
				loaded = true;
				scatt.AddSeries(Ainf_source[id]).Legend(Format(t_("A∞_%s"), nameType)).
						MarkStyle<SquareMarkPlot>().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_DOTTED);
				scatt.Stroke(show_w ? 2 : 0, color).SetMarkWidth(show_w ? 0 : 8);
				if (dim)
					scatt.Units(t_("Ns2/m"));
			}
		}
	} else if (dataToShow == Hydro::DATA_B && hy.IsLoadedB()) {
		if (ABFZ_source[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_B, show_w, !dim)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("B_%s"), nameType)).
						SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
			if (dim)
				scatt.Units(t_("Ns/m"));
		}
	} else if (dataToShow == Hydro::DATA_K && hy.IsLoadedB()) {
		if (ABFZ_source[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_K, show_w, !dim)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("K_%s"), nameType)).
						SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
			if (dim)
				scatt.Units(t_("N/m"));
		}
	} else if (dataToShow == Hydro::DATA_FORCE_SC && hy.IsLoadedFsc() && ih >= 0) {
		if (ABFZ_source[id].Init(hy, plot_idf, ih, Hydro::PLOT_FORCE_SC_MA, show_w, !dim)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("Fsc_ma_%s"), nameType)).
						SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
			if (dim)
				scatt.Units(t_("N/m"));
			if (ABFZ_source2[id].Init(hy, plot_idf, ih, Hydro::PLOT_FORCE_SC_PH, show_w, !dim)) {
				loaded = true;
				scatP.AddSeries(ABFZ_source2[id]).Legend(Format(t_("Fsc_ph_%s"), nameType)).Units(t_("rad")).
						SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
			}
		}
	} else if (dataToShow == Hydro::DATA_FORCE_FK && hy.IsLoadedFfk() && ih >= 0) {
		if (ABFZ_source[id].Init(hy, plot_idf, ih, Hydro::PLOT_FORCE_FK_MA, show_w, !dim)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("Ffk_ma_%s"), nameType)).
						SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
			if (dim)
				scatt.Units(t_("N/m"));
			if (ABFZ_source2[id].Init(hy, plot_idf, ih, Hydro::PLOT_FORCE_FK_PH, show_w, !dim)) {
				loaded = true;
				scatP.AddSeries(ABFZ_source2[id]).Legend(Format(t_("Ffk_ph_%s"), nameType)).Units(t_("rad")).
						SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
			}
		}
	} else if (dataToShow == Hydro::DATA_FORCE_EX && hy.IsLoadedFex() && ih >= 0) {
		if (ABFZ_source[id].Init(hy, plot_idf, ih, Hydro::PLOT_FORCE_EX_MA, show_w, !dim)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("Fex_ma_%s"), nameType)).
						SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
			if (dim)
				scatt.Units(t_("N/m"));
			if (ABFZ_source2[id].Init(hy, plot_idf, ih, Hydro::PLOT_FORCE_EX_PH, show_w, !dim)) {
				loaded = true;
				scatP.AddSeries(ABFZ_source2[id]).Legend(Format(t_("Fex_ph_%s"), nameType)).Units(t_("rad")).
						SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
			}
		}
	} else if (dataToShow == Hydro::DATA_RAO && hy.IsLoadedRAO() && ih >= 0) {
		if (ABFZ_source[id].Init(hy, plot_idf, ih, Hydro::PLOT_RAO_MA, show_w, !dim)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("RAO_ma_%s"), nameType)).
						SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
			if (ABFZ_source2[id].Init(hy, plot_idf, ih, Hydro::PLOT_RAO_PH, show_w, !dim)) {
				loaded = true;
				scatP.AddSeries(ABFZ_source2[id]).Legend(Format(t_("RAO_ph_%s"), nameType)).Units(t_("rad")).
						SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
			}
		}
	} else if (dataToShow == Hydro::DATA_STS && hy.IsLoadedA() && hy.IsLoadedB()) {
		if (ABFZ_source[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_Z_MA, show_w, !dim)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("Zmag_%s"), nameType)).
						SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
			if (dim)
				scatt.Units(t_("Ns/m"));
			if (ABFZ_source2[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_Z_PH, show_w, !dim)) {
				loaded = true;
				scatP.AddSeries(ABFZ_source2[id]).Legend(Format(t_("Zph_%s"), nameType)).Units(t_("rad")).
						SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
			}
		}
	} else if (dataToShow == Hydro::DATA_STS2 && hy.IsLoadedStateSpace()) {
		if (ABFZ_source[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_Z_MA, show_w, !dim)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("Kr_mag %s"), hy.name)).
						SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetMarkColor(color).
						Stroke(2, color);//.Units("dB");
			if (ABFZ_source2[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_Z_PH, show_w, !dim)) {
				loaded = true;
				scatP.AddSeries(ABFZ_source2[id]).Legend(Format(t_("Kr_ph %s"), hy.name)).Units(t_("rad")).
						SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
			}
		}
		if (TFS_source[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_TFS_MA, show_w, !dim)) {
			loaded = true;
			const Upp::Color &bcolor = LtRed();
			scatt.AddSeries(TFS_source[id]).Legend(Format(t_("TFS_mag %s"), hy.name)).
						SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetMarkColor(bcolor).
						Stroke(4, bcolor).Dash(LINE_DASH_DOT);//.Units("dB");
			if (TFS_source2[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_TFS_PH, show_w, !dim)) {
				loaded = true;
				scatP.AddSeries(TFS_source2[id]).Legend(Format(t_("TFS_ph %s"), hy.name)).Units(t_("rad")).
						SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetMarkColor(bcolor).
						Stroke(2, bcolor).Dash(LINE_SOLID);
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