// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2021, the BEMRosetta author and contributors
#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <TabBar/TabBar.h>

#include <BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;

#define IMAGECLASS Img2
#define IMAGEFILE <BEMRosetta/main.iml>
#include <Draw/iml.h>

#include "main.h"
#include "clip.brc"


void MainSetupFOAMM::Init() {
	CtrlLayout(*this);
	
	arrayCases.SetLineCy(EditField::GetStdHeight());
	arrayCases.AddColumn(t_("     Sel"), 30);
	arrayCases.AddColumn(t_("Body"), 20);
	arrayCases.AddColumn(t_("Row"), 20);	
	arrayCases.AddColumn(t_("Column"), 20);
	arrayCases.AddColumn(t_("From (rad/s)"), 40);
	arrayCases.AddColumn(t_("To (rad/s)"), 40);
	arrayCases.AddColumn(t_("Frequencies (rad/s)"), 60);
	arrayCases.WhenSel = [&] {WhenSelArrayCases();};
	
	selectAll << [&] {
			for (int i = 0; i < arrayCases.GetCount(); ++i)
				arrayCases.Set(i, 0, ~selectAll);
		};
	setToAll  << [&] {
			int id = arrayCases.GetCursor();
			if (id < 0)
				return;
			for (int i = 0; i < arrayCases.GetCount(); ++i) {
				if (id != i) {
					arrayCases.Set(i, 4, arrayCases.Get(id, 4));
					arrayCases.Set(i, 5, arrayCases.Get(id, 5));
					arrayCases.Set(i, 6, arrayCases.Get(id, 6));
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
	MainBEM &mainBEM = GetDefinedParent<MainBEM>(this);
	int id = ArrayModel_IdHydro(mainBEM.listLoaded);
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
	foammLogo.SetHyperlink("http://www.eeng.nuim.ie/coer/downloads/");
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
		if (!IsNull(selector.Get(i))) {
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
	arrayCases.Clear();
	options.Clear();
	
	id = _id;
	
	ASSERT(id < Bem().hydros.size());
	
	const Hydro &hydro = Bem().hydros[id].hd();
	
	for (int ib = 0; ib < hydro.Nb; ++ib) {
		for (int idf = 0; idf < 6; ++idf) {
			for (int jdf = 0; jdf < 6; ++jdf) {
				if (!bem.onlyDiagonal || idf == jdf) {
					int _idf = hydro.GetOrder()[ib*6 + idf];
					int _jdf = hydro.GetOrder()[ib*6 + jdf];
	
					if (hydro.IsLoadedA() && hydro.IsLoadedB() && 
						!IsNull(hydro.A[_idf][_jdf][0]) && !IsNull(hydro.B[_idf][_jdf][0])) {
						arrayCases.Add(false, ib+1, BEM::StrDOF(idf), BEM::StrDOF(jdf));
						int row = arrayCases.GetCount()-1;
						arrayCases.SetCtrl(row, 0, options.Add());
						options.Top() << [=] {options[row].SetFocus();};
					}
				}
			}
		}
	}
	if (arrayCases.GetCount() > 0)
		arrayCases.SetCursor(0);
}

void MainSetupFOAMM::WhenSelArrayCases() {
	try {
		int row = arrayCases.GetCursor();
		if (row < 0)
			return;
	
		bool opChoose = arrayCases.Get(row, 0);
		fromFreq <<= arrayCases.Get(row, 4);
		toFreq   <<= arrayCases.Get(row, 5);

		selector.Clear();
			
		String freqs = arrayCases.Get(row, 6);
		Upp::Vector<String> afreqs = Split(freqs, ';');
		for (int i = 0; i < afreqs.size(); ++i)
			selector.AddField(ScanDouble(afreqs[i]));
		
		if (opChoose)
			Status(Check(~fromFreq, ~toFreq, ~freqs));
		
		const Hydro &hydro = Bem().hydros[id].hd();
		
		int ib = int(arrayCases.Get(row, 1)) - 1;
		int idf = BEM::DOFStr(arrayCases.Get(row, 2));
		int jdf = BEM::DOFStr(arrayCases.Get(row, 3));
		
		plots.Init(idf + 6*ib, jdf + 6*ib, Hydro::DATA_STS);
		MainBEM &mbm = GetDefinedParent<MainBEM>(this);
		plots.Load(hydro, mbm);
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
		return;
	}
}

void MainSetupFOAMM::WhenArrayCases() {
	int row = arrayCases.GetCursor();
	if (row < 0)
		return;
	
	arrayCases.Set(row, 4, ~fromFreq);
	arrayCases.Set(row, 5, ~toFreq);

	Upp::Vector<double> freqs;
	for (int i = 0; i < selector.size(); ++i) {
		double freq = selector.Get(i);
		if (!IsNull(freq)) {
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
	
	arrayCases.Set(row, 6, sfreqs);
	
	Status(Check(~fromFreq, ~toFreq, sfreqs));
	
	plots.RefreshScatter();
}

String MainSetupFOAMM::Check(double fromFreq, double toFreq, String freqs) {
	if (IsNull(fromFreq))
		return t_("From: frequency is empty");

	if (IsNull(toFreq))
		return t_("To: frequency is empty");

	if (toFreq <= fromFreq) 
		return t_("From: frequency has to be lower than To: frequency");
			
	Upp::Vector<String> afreqs = Split(freqs, ';');
	if (afreqs.IsEmpty()) 
		return t_("No frequency has been selected");
	
	Upp::Vector<double> unique;
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

bool MainSetupFOAMM::Get(Upp::Vector<int> &ibs, Upp::Vector<int> &idfs, Upp::Vector<int> &jdfs,
		Upp::Vector<double> &froms, Upp::Vector<double> &tos, Upp::Vector<Upp::Vector<double>> &freqs) {
	for (int row = 0; row < arrayCases.GetCount(); ++row) {
		bool proc = arrayCases.Get(row, 0);
		if (proc) {
			int ib = int(arrayCases.Get(row, 1))-1;
			ibs << ib;
			String sidf = arrayCases.Get(row, 2);
			idfs << BEM::DOFStr(sidf);
			String sjdf = arrayCases.Get(row, 3);
			jdfs << BEM::DOFStr(sjdf);
			double from = arrayCases.Get(row, 4);
			double to = arrayCases.Get(row, 5);
			String strfreqs = arrayCases.Get(row, 6);
			String err = Check(from, to, strfreqs);
			if (!err.IsEmpty()) {
				Exclamation(Format(t_("Problem in body %d (%s, %s): %s"), ib+1, sidf, sjdf, err));
				return false;		
			}
			froms << from;
			tos << to;
			Upp::Vector<double> &f = freqs.Add();
			Upp::Vector<String> fs = Split(strfreqs, ';');
			for (int i = 0; i < fs.size(); ++i)
				f << ScanDouble(fs[i]);
		}
	}
	if (ibs.IsEmpty()) {
		Exclamation(t_("No case has been selected"));
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
	Upp::Vector<int> ibs, idfs, jdfs;
	Upp::Vector<double> froms, tos;
	Upp::Vector<Upp::Vector<double>> freqs;
	String ret;
	
	try {
		MainBEM &mainBEM = GetDefinedParent<MainBEM>(this);
		
		int id = ArrayModel_IdHydro(mainBEM.listLoaded);
		if (id < 0)
			return false;
		if (mainBEM.listLoaded.GetCount() != 1 && ArrayCtrlSelectedGetCount(mainBEM.listLoaded) != 1)
			return false;
		
		if (!setup->Get(ibs, idfs, jdfs, froms, tos, freqs))
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
		Foamm &foamm = static_cast<Foamm&>(Bem().hydros[id]);
		foamm.Get(ibs, idfs, jdfs, froms, tos, freqs,
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
					Exclamation(t_("FOAMM message:&") + DeQtfLf(str));
				}
				ProcessEvents(); 
			});
	} catch (Exc e) {
		ret = DeQtfLf(e);
	}
	foammWorking.Hide();
	foammWorking.Stop();
	status.Hide();
	progress.Hide();
	butCancel.Hide();
	butLoad.Enable();
	if (!ret.IsEmpty()) {
		Exclamation(ret);
		return false;
	}
	return true;
}
