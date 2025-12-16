// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <SurfaceCanvas/SurfaceCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <TabBar/TabBar.h>
#include <DropGrid/DropGrid.h>

#include <BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;

#define IMAGECLASS Img2
#define IMAGEFILE <BEMRosetta/main.iml>
#include <Draw/iml.h>

#include "main.h"
#include "clip.brc"


void MainSetupFOAMM::Init() {
	CtrlLayout(*this);
	
	arrayFoamm.SetLineCy(EditField::GetStdHeight());
	arrayFoamm.AddColumn(t_("     Sel"), 30);
	arrayFoamm.AddColumn(t_("Body"), 20);
	arrayFoamm.AddColumn(t_("Dof"), 20);	
	arrayFoamm.AddColumn(t_("Body"), 20);	
	arrayFoamm.AddColumn(t_("Dof"), 20);
	arrayFoamm.AddColumn(t_("From (rad/s)"), 40);
	arrayFoamm.AddColumn(t_("To (rad/s)"), 40);
	arrayFoamm.AddColumn(t_("Frequencies (rad/s)"), 60);
	arrayFoamm.WhenSel = [&] {WhenSelArrayCases();};
	
	selectAll << [&] {
			for (int i = 0; i < arrayFoamm.GetCount(); ++i)
				arrayFoamm.Set(i, 0, ~selectAll);
		};
	setToAll  << [&] {
			int id = arrayFoamm.GetCursor();
			if (id < 0)
				return;
			for (int i = 0; i < arrayFoamm.GetCount(); ++i) {
				if (id != i) {
					arrayFoamm.Set(i, 5, arrayFoamm.Get(id, 5));
					arrayFoamm.Set(i, 6, arrayFoamm.Get(id, 6));
					arrayFoamm.Set(i, 7, arrayFoamm.Get(id, 7));
				}
			}
		};
	
	fromFreq.WhenAction = THISBACK(WhenArrayCases);
	frameSet.Add(fromFreq.GetRectEnter());
	toFreq.WhenAction   = THISBACK(WhenArrayCases);
	frameSet.Add(toFreq.GetRectEnter());
	selector.Init(THISBACK(WhenArrayCases), frameSet);
	frameSet.WhenEnter = THISBACK(WhenFocus);
	frameSet.Set(fromFreq.GetRectEnter());
	
	rectPlots.Add(plots.SizePos());	
	
	plots.Init(true);
	
	plots.scatt.WhenPainter = THISBACK1(OnPainter, &plots.scatt);
	plots.scatt.WhenMouseClick = THISBACK1(OnMouse, &plots.scatt);
	plots.scatP.WhenPainter = THISBACK1(OnPainter, &plots.scatP);
	plots.scatP.WhenMouseClick = THISBACK1(OnMouse, &plots.scatP);
}

void MenuFOAMM::OnCursor() {
	MainBEM &mainBEM = GetParentCtrl<MainBEM>(this);
	int id = ArrayModel_IndexHydro(mainBEM.listLoaded);
	if (id < 0)
		return;
	if (ArrayCtrlSelectedGetCount(mainBEM.listLoaded) > 1)
		return;
	setup->WhenSelArrayModel(id, Bem());	
}

void MenuFOAMM::Init(MainBEM &mainBEM, MainSetupFOAMM &_setup) {
	CtrlLayout(*this);
	setup = &_setup;
	
	butLoad.WhenAction 	= [&] {
		if (OnFOAMM()) {
			int id = mainBEM.mainTab.Find(mainBEM.mainStateSpace);
			mainBEM.mainTab.GetItem(id).Enable(true);
			int idPlot = mainBEM.menuTab.Find(mainBEM.menuPlot);
			mainBEM.menuTab.GetItem(idPlot).Enable(true);
			mainBEM.mainTab.Set(0);
			mainBEM.mainTab.Set(id);
		}
	};
	
	foammLogo.Set(Img2::FOAMM());
	foammLogo.SetHyperlink("https://coer.maynoothuniversity.ie/downloads/");
	foammWorking.LoadBuffer(String(animatedStar, animatedStar_length));
	foammWorking.Hide();
	status.Hide();
	progress.Hide();
	butCancel.Hide();
	butCancel << [&] {isCancelled = true;};
}

void MainSetupFOAMM::WhenFocus() {
	plots.RefreshScatter();
}

void MainSetupFOAMM::OnPainter(Painter &w, ScatterCtrl *pscat) {
	ScatterCtrl &scat = *pscat;
	int plotW = scat.GetPlotWidth(), plotH = scat.GetPlotHeight();
	
	for (int i = 0; i < selector.size(); ++i) {
		if (IsNum(selector.Get(i))) {
			double xFreq = scat.GetPosX(selector.Get(i));
			if (selector.IsSelected(i))
				DrawLineOpa(w, xFreq, 0, xFreq, plotH, 1, 1, 2, LtCyan(), "2 2");
			else
				DrawLineOpa(w, xFreq, 0, xFreq, plotH, 1, 1, 2, LtBlue(), "2 2");
		}
	}

	if (!IsNull(fromFreq)) {
		double xFrom = scat.GetPosX(~fromFreq);
		FillRectangleOpa(w, 0, 0, xFrom, plotH, 0.5, Null, LtBlue());
	}
	if (!IsNull(toFreq)) {
		double xTo = scat.GetPosX(~toFreq);
		FillRectangleOpa(w, xTo, 0, plotW, plotH, 0.5, Null, LtBlue());
	}
}

void MainSetupFOAMM::OnMouse(Point p, dword, ScatterCtrl::MouseAction action, ScatterCtrl *pscat) {
	ScatterCtrl &scat = *pscat;
	if (action != ScatterCtrl::LEFT_DOWN && action != ScatterCtrl::LEFT_MOVE)
		return; 
		
	double freq = scat.GetRealPosX(p.x);

	if (!scat.IsEmpty()) {
		DataSource &data = scat.GetDataSource(0);
		freq = data.x(data.ClosestX(freq));
	}
	
	if (fromFreq.IsShownFrame()) {	
		fromFreq <<= freq;
		fromFreq.WhenAction();
	} else if (toFreq.IsShownFrame()) {	
		toFreq <<= freq;
		toFreq.WhenAction();
	} else {
		int id = selector.GetSelected();
		if (id >= 0) 
			selector.Set(id, freq);
	}
}

void MainSetupFOAMM::WhenSelArrayModel(int _id, BEM &bem) {
	arrayFoamm.Clear();
	options.Clear();
	
	idd = _id;
	
	ASSERT(_id < Bem().hydros.size());
	
	const Hydro &hy = Bem().hydros[_id];
	
	for (int ib0 = 0; ib0 < hy.dt.Nb; ++ib0) {
		for (int idf = 0; idf < 6; ++idf) {
			for (int ib1 = 0; ib1 < hy.dt.Nb; ++ib1) {
				for (int jdf = 0; jdf < 6; ++jdf) {
					if (!bem.onlyDiagonal || idf == jdf) {
						int _idf = ib0*6 + idf;
						int _jdf = ib1*6 + jdf;
		
						if (hy.IsLoadedA(_idf, _jdf) && hy.IsLoadedB(_idf, jdf)) {
							arrayFoamm.Add(false, ib0+1, BEM::StrDOF(idf), ib1+1, BEM::StrDOF(jdf));
							int row = arrayFoamm.GetCount()-1;
							arrayFoamm.SetCtrl(row, 0, options.Add());
							options.Top() << [=] {options[row].SetFocus();};
						}
					}
				}
			}
		}
	}
	if (arrayFoamm.GetCount() > 0)
		arrayFoamm.SetCursor(0);
}

void MainSetupFOAMM::WhenSelArrayCases() {
	try {
		if (idd < 0)
			return;
		int row = arrayFoamm.GetCursor();
		if (row < 0)
			return;
	
		bool opChoose = arrayFoamm.Get(row, 0);
		int ib0 = int(arrayFoamm.Get(row, 1)) - 1;
		int idf = BEM::DOFStr(arrayFoamm.Get(row, 2));
		int ib1 = int(arrayFoamm.Get(row, 3)) - 1;
		int jdf = BEM::DOFStr(arrayFoamm.Get(row, 4));
		fromFreq <<= arrayFoamm.Get(row, 5);
		toFreq   <<= arrayFoamm.Get(row, 6);

		selector.Clear();
			
		String freqs = arrayFoamm.Get(row, 7);
		UVector<String> afreqs = Split(freqs, ';');
		for (int i = 0; i < afreqs.size(); ++i)
			selector.AddField(ScanDouble(afreqs[i]));
		
		if (opChoose)
			Status(Check(~fromFreq, ~toFreq, ~freqs));
		
		const Hydro &hy = Bem().hydros[idd];
		
		plots.Init(idf + 6*ib0, jdf + 6*ib1, Hydro::DATA_STS);
		MainBEM &mbm = GetParentCtrl<MainBEM>(this);
		plots.Load(hy, mbm);
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
		return;
	}
}

void MainSetupFOAMM::WhenArrayCases() {
	int row = arrayFoamm.GetCursor();
	if (row < 0)
		return;
	
	arrayFoamm.Set(row, 5, ~fromFreq);
	arrayFoamm.Set(row, 6, ~toFreq);

	UVector<double> freqs;
	for (int i = 0; i < selector.size(); ++i) {
		double freq = selector.Get(i);
		if (IsNum(freq)) {
			if (!plots.scatt.IsEmpty()) { 
				DataSource &data = plots.scatt.GetDataSource(0);
				freq = data.x(data.ClosestX(freq));
			}
		}
		freqs << freq;
	}
	
	Sort(freqs);
	
	String sfreqs;
	for (int i = 0; i < freqs.size(); ++i) {
		if (!sfreqs.IsEmpty())
			sfreqs << ";";
		sfreqs << freqs[i];
	}
	
	arrayFoamm.Set(row, 7, sfreqs);
	
	Status(Check(~fromFreq, ~toFreq, sfreqs));
	
	plots.RefreshScatter();
}

String MainSetupFOAMM::Check(double fromFreq, double toFreq, String freqs) {
	if (!IsNum(fromFreq))
		return t_("From: frequency is empty");

	if (!IsNum(toFreq))
		return t_("To: frequency is empty");

	if (toFreq <= fromFreq) 
		return t_("From: frequency has to be lower than To: frequency");
			
	UVector<String> afreqs = Split(freqs, ';');
	if (afreqs.IsEmpty()) 
		return t_("No frequency has been selected");
	
	UVector<double> unique;
	for (int i = 0; i < afreqs.size(); ++i) {
		double freq = ScanDouble(afreqs[i]);
		if (freq < fromFreq || freq > toFreq) 
			return t_("Selected frequencies have to be between lower and higher limits");
		FindAddRatio(unique, freq, 0.001);
	}
	if (unique.size() != afreqs.size()) 
		return t_("Some selected frequencies are repeated");
	
	return String("");
}

bool MainSetupFOAMM::Get(UVector<int> &ib0s, UVector<int> &idfs, UVector<int> &ib1s, UVector<int> &jdfs,
		UVector<double> &froms, UVector<double> &tos, UVector<UVector<double>> &freqs) {
	for (int row = 0; row < arrayFoamm.GetCount(); ++row) {
		bool proc = arrayFoamm.Get(row, 0);
		if (proc) {
			int ib0 = int(arrayFoamm.Get(row, 1))-1;
			ib0s << ib0;
			String sidf = arrayFoamm.Get(row, 2);
			idfs << BEM::DOFStr(sidf);
			int ib1 = int(arrayFoamm.Get(row, 3))-1;
			ib1s << ib1;
			String sjdf = arrayFoamm.Get(row, 4);
			jdfs << BEM::DOFStr(sjdf);
			double from = arrayFoamm.Get(row, 5);
			double to   = arrayFoamm.Get(row, 6);
			String strfreqs = arrayFoamm.Get(row, 7);
			String err = Check(from, to, strfreqs);
			if (!err.IsEmpty()) {
				BEM::PrintError(Format(t_("Problem in body %d (%s), body %s (%s): %s"), ib0+1, sidf, ib1+1, sjdf, err));
				return false;		
			}
			froms << from;
			tos << to;
			UVector<double> &f = freqs.Add();
			UVector<String> fs = Split(strfreqs, ';');
			for (int i = 0; i < fs.size(); ++i)
				f << ScanDouble(fs[i]);
		}
	}
	if (ib0s.IsEmpty()) {
		BEM::PrintError(t_("No case has been selected"));
		return false;			
	}
	return true;
}

void MainSetupFOAMM::Clear() {
	plots.Clear();
}

void MenuFOAMM::Clear() {
//	arrayModel.Clear();
}
		
bool MenuFOAMM::OnFOAMM() {
	UVector<int> ib0s, ib1s, idfs, jdfs;
	UVector<double> froms, tos;
	UVector<UVector<double>> freqs;
	String ret;
	
	try {
		MainBEM &mainBEM = GetParentCtrl<MainBEM>(this);
		
		int idx = ArrayModel_IndexHydro(mainBEM.listLoaded);
		if (idx < 0)
			return false;
		if (mainBEM.listLoaded.GetCount() != 1 && ArrayCtrlSelectedGetCount(mainBEM.listLoaded) != 1)
			return false;
		
		if (!setup->Get(ib0s, idfs, ib1s, jdfs, froms, tos, freqs))
			return false;
		
		foammWorking.Show();
		foammWorking.Play();
		status.Show();
		progress.Show();
		butCancel.Show();
		butLoad.Disable();
		WaitCursor wait;
		isCancelled = false;
		status.SetText(t_("Starts processing"));
		Foamm &foamm = static_cast<Foamm&>(Bem().hydros[idx]);
		foamm.Get(ib0s, idfs, ib1s, jdfs, froms, tos, freqs,
			[&](String str, int pos)->bool {
				if (!str.IsEmpty())
					status.SetText(str);	
				if (IsNull(pos))
					; 
				else
					progress.Set(pos, 100);
				ProcessEvents(); 
				return isCancelled;
			}, 
			[&](String str) {
				if (!str.IsEmpty()) {
					str.Replace("\r", "");
					str.Replace("\n\n", "\n");
					BEM::PrintError(t_("FOAMM message:&") + DeQtfLf(str));
				}
				ProcessEvents(); 
			});
	} catch (Exc e) {
		ret = e;
	}
	foammWorking.Hide();
	foammWorking.Stop();
	status.Hide();
	progress.Hide();
	butCancel.Hide();
	butLoad.Enable();
	if (!ret.IsEmpty()) {
		BEM::PrintError(ret);
		return false;
	}
	return true;
}
