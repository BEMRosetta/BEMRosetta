// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <TabBar/TabBar.h>
#include <DropGrid/DropGrid.h>

#include <BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;

#include "main.h"


void MainPlot::Init(bool vert) {
	if (!isInit) {
		int len = StdFont().GetHeight();
		
		scatt.SetMargin(6*len, len, len, 4*len).SetTitleFont(SansSerifZ(12));
		scatP.SetMargin(6*len, len, len, 4*len).SetTitleFont(SansSerifZ(12));
		scatt.LinkedWith(scatP);
		
		splitCompare.Horz(scatt.SizePos(), compare.SizePos());
		splitCompare.SetPositions(7000, 10000).SetInitialPositionId(1).SetButtonNumber(1).SetButtonWidth(len);
		if (vert)
			splitter.Vert(splitCompare.SizePos(), scatP.SizePos());
		else
			splitter.Horz(splitCompare.SizePos(), scatP.SizePos());
		Add(splitter.SizePos());
		scatt.ShowAllMenus().SetSciExpTop();
		scatP.ShowAllMenus().SetSciExpTop();
		
		compare.Init(scatt, splitCompare);
		
		isInit = true;
	}
	Clear();
}

void MainPlot::Init(int _idf, double jdf_ih, Hydro::DataToShow _dataToShow, double _heading1) {
	MainPlot::Init(true);
	
	plot_idf = _idf;
	plot_jdf = int(jdf_ih);
	heading0 = jdf_ih;
	heading1 = _heading1;
	dataToShow = _dataToShow;
	
	compare.Init(dataToShow);
	
	String smp_up, smp_down;
	if (show_ma_ph) {
		smp_up = t_("mag");
		smp_down = t_("phase");
	} else {
		smp_up = t_("real");
		smp_down = t_("imag");
	}
						
	String title, title2, labelY, labelY2;
	switch (dataToShow) {
	case Hydro::DATA_A:	title = Format(t_("Added mass %s"), BEM::StrBDOF2(plot_idf, plot_jdf, false));		
						labelY = t_("Added mass");
						splitter.SetPos(10000, 0);				
						break;		
	case Hydro::DATA_B:	title = Format(t_("Radiation damping %s"), BEM::StrBDOF2(plot_idf, plot_jdf, false));
						labelY = t_("Radiation damping");
						splitter.SetPos(10000, 0);	
						break;
	case Hydro::DATA_MD:title = Format(t_("Mean drift %s heading %.1f:%.1fº"), BEM::StrBDOF(plot_idf, false), heading0, heading1);
						labelY = t_("Mean drift");
						splitter.SetPos(10000, 0);	
						break;						
	case Hydro::DATA_AINFW:	title = Format(t_("Added mass at infinity (ω) %s"), BEM::StrBDOF2(plot_idf, plot_jdf, false));		
						labelY = t_("Added mass at infinity (ω)");
						splitter.SetPos(10000, 0);				
						break;
	case Hydro::DATA_KIRF:	title = Format(t_("Kirf Impulse Response Function %s"), BEM::StrBDOF2(plot_idf, plot_jdf, false));
						labelY = t_("Kirf Impulse Response Function");
						splitter.SetPos(10000, 0);	
						break;
	case Hydro::DATA_FORCE_SC:	title = Format(t_("Diffraction scattering force %s heading %.1fº"), BEM::StrBDOF(plot_idf, false), heading0);
						title2 = title + " " + smp_down;
						title  = title + " " + smp_up;
						labelY  = t_("Diffraction scattering force") + S(" ") + smp_up;
						labelY2 = t_("Diffraction scattering force") + S(" ") + smp_down;
						splitter.SetPos(5000, 0);	
						break;
	case Hydro::DATA_FORCE_FK:	title = Format(t_("Froude-Krylov force %s heading %.1fº"), BEM::StrBDOF(plot_idf, false), heading0);
						title2 = title + " " + smp_down;
						title  = title + " " + smp_up;
						labelY  = t_("Froude-Krylov force") + S(" ") + smp_up;		
						labelY2 = t_("Froude-Krylov force") + S(" ") + smp_down;
						splitter.SetPos(5000, 0);			
						break;
	case Hydro::DATA_FORCE_EX:	title = Format(t_("Excitation Force %s heading %.1fº"), BEM::StrBDOF(plot_idf, false), heading0);
						title2 = title + " " + smp_down;
						title  = title + " " + smp_up;
						labelY  = t_("Excitation force") + S(" ") + smp_up;		
						labelY2 = t_("Excitation force") + S(" ") + smp_down;
						splitter.SetPos(5000, 0);				
						break;
	case Hydro::DATA_RAO:	title = Format(t_("Response Amplitude Operator %s heading %.1fº"), BEM::StrBDOF(plot_idf, false), heading0);
						title2 = title + " " + smp_down;
						title  = title + " " + smp_up;
						labelY  = t_("RAO force") + S(" ") + smp_up;		
						labelY2 = t_("RAO force") + S(" ") + smp_down;
						splitter.SetPos(5000, 0);							
						break;
	case Hydro::DATA_STS:	title = Format(t_("Magnitude Kr(ω) = B(ω)+jω{A(ω)-A∞} %s"), BEM::StrBDOFFull(plot_idf, plot_jdf));		
						title2 = Format(t_("Phase Kr(ω) = B(ω)+jω{A(ω)-A∞} %s"), BEM::StrBDOFFull(plot_idf, plot_jdf));		
						if (show_ma_ph) {
							labelY = t_("Magnitude");				
							labelY2 = t_("Phase");
						} else {
							
						}
						splitter.SetPos(5000, 0);				
						break;
	case Hydro::DATA_STS2:	title = Format(t_("Magnitude response %s"), BEM::StrBDOF2(plot_idf, plot_jdf, false));
						title2 = Format(t_("Frequency response %s"), BEM::StrBDOF2(plot_idf, plot_jdf, false));
						splitter.SetPos(5000, 0);
	}
	scatt.SetTitle(title);
	scatt.SetLabelY(labelY);
	scatP.SetTitle(title2);
	scatP.SetLabelY(labelY2);
	
	scatt.SetLegendWithUnits(Bem().legend_w_units);
	scatP.SetLegendWithUnits(Bem().legend_w_units);
}

bool MainPlot::Load(const UArray<HydroClass> &hydro, const MainBEM &mbm, const UVector<int> &ids) {
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
	show_ma_ph = mbm.menuPlot.opMP == 0;
	
	if (show_w && dataToShow != Hydro::DATA_KIRF) {
		scatt.SetLabelX(t_("ω [rad/s]"));
		scatP.SetLabelX(t_("ω [rad/s]"));
	} else {
		scatt.SetLabelX(t_("T [s]"));
		scatP.SetLabelX(t_("T [s]"));
	}
	bool loaded = false;
	for (int id = 0; id < ids.size(); ++id) {
		const Hydro &hy = hydro[ids[id]].hd();
		bool ld = false;
		LoadEach(hy, id, ld, hy.GetId());
		loaded |= ld;		// If only one is loaded, then loaded = true
	}
	if (mbm.menuPlot.autoFit) {
		scatt.ZoomToFit(true, true);
		if (mbm.menuPlot.fromY0) {
			double yRange = max<double>(0, scatt.GetYMin()) + scatt.GetYRange();
			scatt.SetXYMin(Null, 0).SetRange(Null, yRange);
		}
		if (show_ma_ph) {
			scatP.ZoomToFit(true, false);
			scatP.SetXYMin(Null, -M_PI).SetRange(Null, 2*M_PI).SetMajorUnits(Null, 1);
			scatP.SetMinUnits(Null, M_PI-3);
		} else if (mbm.menuPlot.fromY0) {
			double yRange = max<double>(0, scatP.GetYMin()) + scatP.GetYRange();
			scatP.SetXYMin(Null, 0).SetRange(Null, yRange);
		} else
			scatP.ZoomToFit(true, true);
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
	show_ma_ph = mainBem.menuPlot.opMP == 0;
	
	if (show_w && dataToShow != Hydro::DATA_KIRF) {
		scatt.SetLabelX(t_("ω [rad/s]"));
		scatP.SetLabelX(t_("ω [rad/s]"));
	} else {
		scatt.SetLabelX(t_("T [s]"));
		scatP.SetLabelX(t_("T [s]"));
	}
	show_ma_ph = mainBem.menuPlot.opMP;
	
	bool loaded = false;
	int idc = ArrayModel_Id(mainBem.listLoaded);
	LoadEach(hy, 0, loaded, idc);
	
	if (mainBem.menuPlot.autoFit) {
		scatt.ZoomToFit(true, true);
		if (mainBem.menuPlot.fromY0) {
			double ymax = scatt.GetYMin() + scatt.GetYRange();
			scatt.SetXYMin(Null, 0).SetRange(Null, ymax);
		}
		if (show_ma_ph) {
			scatP.ZoomToFit(true, false);
			scatP.SetXYMin(Null, -M_PI).SetRange(Null, 2*M_PI).SetMajorUnits(Null, 1);
			scatP.SetMinUnits(Null, M_PI-3);
		} else
			scatP.ZoomToFit(true, true);
	}
	return loaded;
}

void MainPlot::LoadEach(const Hydro &hy, int id, bool &loaded, int idc) {
	int ih = -1;
	if (dataToShow != Hydro::DATA_A && dataToShow != Hydro::DATA_B && dataToShow != Hydro::DATA_AINFW) {
		if (dataToShow == Hydro::DATA_MD) {
			std::complex<double> h(heading0, heading1);
			ih = hy.GetHeadIdMD(h);
		} else
			ih = hy.GetHeadId(heading0);
	}
	String nameType = hy.name;
	if (Bem().legend_w_solver)
		nameType << Format("(%s)", hy.GetCodeStrAbr());
	
	String sids = BEM::strDOFnum_sub[plot_idf];
	if (dataToShow == Hydro::DATA_A || dataToShow == Hydro::DATA_AINFW || dataToShow == Hydro::DATA_B || dataToShow == Hydro::DATA_KIRF)
		sids += BEM::strDOFnum_sub[plot_jdf];
	
	if (idc < 0)
		idc = id;
	const Upp::Color &color = GetColorId(idc);
	String st = show_ma_ph ? t_("ₘₐ") : t_("ᵣ");
	String sp = show_ma_ph ? t_("ₚₕ") : t_("ᵢₘ");
	
	if (dataToShow == Hydro::DATA_A) {
		if (hy.IsLoadedA(plot_idf, plot_jdf)) {
			if (ABFZ_source[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_A, show_w, !dim, show_ma_ph)) {
				loaded = true;
				scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("A%s %s"), sids, nameType)).
						SetMarkWidth(markW).MarkStyle<CircleMarkPlot>().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
				if (dim)
					scatt.Units(Hydro::A_units(!dim, plot_idf, plot_jdf));
			}
		}
		if (hy.IsLoadedAinf()) {
			if (Ainf_source[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_AINF, show_w, !dim, show_ma_ph)) {
				loaded = true;
				scatt.AddSeries(Ainf_source[id]).Legend(Format(t_("A∞%s %s"), sids, nameType)).
						MarkStyle<SquareMarkPlot>().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_DOTTED);
				scatt.Stroke(show_w ? 2 : 0, color).SetMarkWidth(show_w ? 0 : 8);
				if (dim)
					scatt.Units(Hydro::A_units(!dim, plot_idf, plot_jdf));
			}
		}
		if (hy.IsLoadedA0()) {
			if (A0_source[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_A0, show_w, !dim, show_ma_ph)) {
				loaded = true;
				scatt.AddSeries(A0_source[id]).Legend(Format(t_("A₀%s %s"), sids, nameType)).
						SetMarkWidth(!show_w ? 0 : 8).SetMarkColor(color).MarkStyle<SquareMarkPlot>().
						Stroke(2, color).Dash(LINE_DOTTED);
				if (dim)
					scatt.Units(Hydro::A_units(!dim, plot_idf, plot_jdf));
			}
		}
	} else if (dataToShow == Hydro::DATA_AINFW) {
		if (hy.IsLoadedAinf_w(plot_idf, plot_jdf)) {
			if (ABFZ_source[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_AINFW, show_w, !dim, show_ma_ph)) {
				loaded = true;
				scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("A∞%s(ω)%s"), sids, nameType)).
						SetMarkWidth(markW).SetMarkStyleType().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
				if (dim)
					scatt.Units(Hydro::A_units(!dim, plot_idf, plot_jdf));
			}
		}
		if (hy.IsLoadedAinf()) {
			if (Ainf_source[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_AINF, show_w, !dim, show_ma_ph)) {
				loaded = true;
				scatt.AddSeries(Ainf_source[id]).Legend(Format(t_("A∞%s %s"), sids, nameType)).
						MarkStyle<SquareMarkPlot>().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_DOTTED);
				scatt.Stroke(show_w ? 2 : 0, color).SetMarkWidth(show_w ? 0 : 8);
				if (dim)
					scatt.Units(Hydro::A_units(!dim, plot_idf, plot_jdf));
			}
		}
	} else if (dataToShow == Hydro::DATA_B && hy.IsLoadedB(plot_idf, plot_jdf)) {
		if (ABFZ_source[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_B, show_w, !dim, show_ma_ph)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("B%s %s"), sids, nameType)).
						SetMarkWidth(markW).SetMarkStyleType().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
			if (dim)
				scatt.Units(Hydro::B_units(!dim, plot_idf, plot_jdf));
		}
	} else if (dataToShow == Hydro::DATA_MD && ih >= 0 && hy.IsLoadedMD()) {
		if (ABFZ_source[id].Init(hy, plot_idf, ih, Hydro::PLOT_MD, show_w, !dim, true)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("MD%s %s %d"), sids, nameType, hy.mdtype)).
						SetMarkWidth(markW).SetMarkStyleType().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
			if (dim)
				scatt.Units(Hydro::MD_units(!dim, plot_idf));
		}
	} else if (dataToShow == Hydro::DATA_KIRF && hy.IsLoadedKirf()) {
		if (ABFZ_source[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_KIRF, show_w, !dim, show_ma_ph)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("K%s %s"), sids, nameType)).
						SetMarkWidth(markW).SetMarkStyleType().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
			if (dim)
				scatt.Units(Hydro::Kirf_units(!dim, plot_idf, plot_jdf));
		}
	} else if (dataToShow == Hydro::DATA_FORCE_SC && ih >= 0 && hy.IsLoadedFsc(plot_idf, ih)) {
		if (ABFZ_source[id].Init(hy, plot_idf, ih, Hydro::PLOT_FORCE_SC_1, show_w, !dim, show_ma_ph)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("Fsc%s%s %s"), st, sids, nameType)).
						SetMarkWidth(markW).SetMarkStyleType().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
			if (dim)
				scatt.Units(Hydro::F_units(!dim, plot_idf));
			if (ABFZ_source2[id].Init(hy, plot_idf, ih, Hydro::PLOT_FORCE_SC_2, show_w, !dim, show_ma_ph)) {
				loaded = true;
				scatP.AddSeries(ABFZ_source2[id]).Legend(Format(t_("Fsc%s%s %s"), sp, sids, nameType)).Units(t_("rad")).
						SetMarkWidth(markW).SetMarkStyleType().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
				if (dim && !show_ma_ph)
					scatP.Units(Hydro::F_units(!dim, plot_idf));
			}
		}
	} else if (dataToShow == Hydro::DATA_FORCE_FK && ih >= 0 && hy.IsLoadedFfk(plot_idf, ih)) {
		if (ABFZ_source[id].Init(hy, plot_idf, ih, Hydro::PLOT_FORCE_FK_1, show_w, !dim, show_ma_ph)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("Ffk%s%s %s"), st, sids, nameType)).
						SetMarkWidth(markW).SetMarkStyleType().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
			if (dim)
				scatt.Units(Hydro::F_units(!dim, plot_idf));
			if (ABFZ_source2[id].Init(hy, plot_idf, ih, Hydro::PLOT_FORCE_FK_2, show_w, !dim, show_ma_ph)) {
				loaded = true;
				scatP.AddSeries(ABFZ_source2[id]).Legend(Format(t_("Ffk%s%s %s"), sp, sids, nameType)).Units(t_("rad")).
						SetMarkWidth(markW).SetMarkStyleType().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
				if (dim && !show_ma_ph)
					scatP.Units(Hydro::F_units(!dim, plot_idf));
			}
		}
	} else if (dataToShow == Hydro::DATA_FORCE_EX && ih >= 0 && hy.IsLoadedFex(plot_idf, ih)) {
		if (ABFZ_source[id].Init(hy, plot_idf, ih, Hydro::PLOT_FORCE_EX_1, show_w, !dim, show_ma_ph)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("Fex%s%s %s"), st, sids, nameType)).
						SetMarkWidth(markW).SetMarkStyleType().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
			if (dim)
				scatt.Units(Hydro::F_units(!dim, plot_idf));
			if (ABFZ_source2[id].Init(hy, plot_idf, ih, Hydro::PLOT_FORCE_EX_2, show_w, !dim, show_ma_ph)) {
				loaded = true;
				scatP.AddSeries(ABFZ_source2[id]).Legend(Format(t_("Fex%s%s %s"), sp, sids, nameType)).Units(t_("rad")).
						SetMarkWidth(markW).SetMarkStyleType().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
				if (dim && !show_ma_ph)
					scatP.Units(Hydro::F_units(!dim, plot_idf));
			}
		}
	} else if (dataToShow == Hydro::DATA_RAO  && ih >= 0 && hy.IsLoadedRAO(plot_idf, ih)) {
		if (ABFZ_source[id].Init(hy, plot_idf, ih, Hydro::PLOT_RAO_1, show_w, !dim, show_ma_ph)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("RAO%s%s %s"), st, sids, nameType)).Units(Hydro::RAO_units(plot_idf)).
						SetMarkWidth(markW).SetMarkStyleType().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
			if (ABFZ_source2[id].Init(hy, plot_idf, ih, Hydro::PLOT_RAO_2, show_w, !dim, show_ma_ph)) {
				loaded = true;
				scatP.AddSeries(ABFZ_source2[id]).Legend(Format(t_("RAO%s%s %s"), sp, sids, nameType)).Units(t_("rad")).
						SetMarkWidth(markW).SetMarkStyleType().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
				if (!show_ma_ph)
					scatP.Units(Hydro::RAO_units(plot_idf));		
			}
		}
	} else if (dataToShow == Hydro::DATA_STS && hy.IsLoadedA() && hy.IsLoadedB()) {
		if (ABFZ_source[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_Z_1, show_w, !dim, show_ma_ph)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("Z%s%s_%s"), st, sids, nameType)).
						SetMarkWidth(markW).SetMarkStyleType().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
			if (dim)
				scatt.Units(t_("Ns/m"));
			if (ABFZ_source2[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_Z_2, show_w, !dim, show_ma_ph)) {
				loaded = true;
				scatP.AddSeries(ABFZ_source2[id]).Legend(Format(t_("Z%s%s %s"), sp, sids, nameType)).Units(t_("rad")).
						SetMarkWidth(markW).SetMarkStyleType().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
				if (dim && !show_ma_ph)
					scatt.Units(t_("Ns/m"));
			}
		}
	} else if (dataToShow == Hydro::DATA_STS2 && hy.IsLoadedStateSpace()) {
		if (ABFZ_source[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_Z_1, show_w, !dim, show_ma_ph)) {
			loaded = true;
			scatt.AddSeries(ABFZ_source[id]).Legend(Format(t_("Kr_%s %s"), st, hy.name)).
						SetMarkWidth(markW).SetMarkStyleType().SetMarkColor(color).
						Stroke(2, color);//.Units("dB");
			if (ABFZ_source2[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_Z_2, show_w, !dim, show_ma_ph)) {
				loaded = true;
				scatP.AddSeries(ABFZ_source2[id]).Legend(Format(t_("Kr_%s %s"), sp, hy.name)).Units(t_("rad")).
						SetMarkWidth(markW).SetMarkStyleType().SetMarkColor(color).
						Stroke(2, color).Dash(LINE_SOLID);
				if (dim && !show_ma_ph)
					scatt.Units(t_("Ns/m"));

			}
		}
		if (TFS_source[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_TFS_1, show_w, !dim, show_ma_ph)) {
			loaded = true;
			const Upp::Color &bcolor = LtRed();
			scatt.AddSeries(TFS_source[id]).Legend(Format(t_("TFS_%s %s"), st, hy.name)).
						SetMarkWidth(markW).SetMarkStyleType().SetMarkColor(bcolor).
						Stroke(4, bcolor).Dash(LINE_DASH_DOT);//.Units("dB");
			if (TFS_source2[id].Init(hy, plot_idf, plot_jdf, Hydro::PLOT_TFS_2, show_w, !dim, show_ma_ph)) {
				loaded = true;
				scatP.AddSeries(TFS_source2[id]).Legend(Format(t_("TFS_%s %s"), sp, hy.name)).Units(t_("rad")).
						SetMarkWidth(markW).SetMarkStyleType().SetMarkColor(bcolor).
						Stroke(2, bcolor).Dash(LINE_SOLID);
				if (dim && !show_ma_ph)
					scatt.Units(t_("Ns/m"));

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