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


void MainABForce::Init(Hydro::DataToShow _dataToShow) {
	Add(tab.SizePos());
	
	dataToShow = _dataToShow;
	
	selTab = 0;
	isFilling = false;
	tab.WhenSet = [&] {
		//LOGTAB(tab);
		if (!isFilling)
			selTab = tab.Get();
	};
}

void MainABForce::Clear() {
	tab.Reset();
	selTab = 0;
}

bool MainABForce::Load(BEM &bem, const UVector<int> &ids, int ih) {
	TempAssign<bool> _isFilling(isFilling, true);
	try {
		MainBEM &mbm = GetDefinedParent<MainBEM>(this);
	
		tab.Reset();
		UArray<HydroClass> &hydros = bem.hydros; 
		if (hydros.IsEmpty() || ids.IsEmpty()) 
			return false;
		
		int sdof = 6*bem.Nb;
		int nloaded = 0;
		if (dataToShow == Hydro::DATA_A || dataToShow == Hydro::DATA_B || dataToShow == Hydro::DATA_AINFW || 
			dataToShow == Hydro::DATA_K) {
			plots.SetCount(sdof);
			for (int idf = 0; idf < sdof; ++idf) {
				plots[idf].SetCount(sdof);
				for (int jdf = 0; jdf < sdof; ++jdf) {
					if (!bem.onlyDiagonal || idf == jdf) {
						plots[idf][jdf].Init(idf, jdf, dataToShow);
						if (plots[idf][jdf].Load(hydros, mbm, ids)) {
							nloaded++;
							if (idf != jdf)
								tab.Add(plots[idf][jdf].SizePos(), Format(t_("%s"), BEM::StrBDOF2(idf, jdf, false)));
							else
								tab.Add(plots[idf][jdf].SizePos(), Format(t_("%s"), BEM::StrBDOF(idf, false)));
						}
					}
				}
			}
		} else {
			ih = max(ih, 0);
			
			if (dataToShow == Hydro::DATA_MD) {
				int Nh = bem.headAllMD.size();
				if (Nh == 0) 
					return false;
				
				plots.SetCount(1);
				plots[0].SetCount(sdof);
				for (int idf = 0; idf < sdof; ++idf) {
					const std::complex<double> &head = bem.headAllMD[ih];
					plots[0][idf].Init(idf, head.real(), dataToShow, head.imag());
					if (plots[0][idf].Load(hydros, mbm, ids))
						nloaded++;
				}
				if (nloaded > 0)
					UpdateHeadMD(bem, ih);
			} else {
				int Nh = bem.headAll.size();
				if (Nh == 0) 
					return false;
				
				plots.SetCount(1);
				plots[0].SetCount(sdof);
				for (int idf = 0; idf < sdof; ++idf) {
					plots[0][idf].Init(idf, bem.headAll[ih/*bem.orderHeadAll[ih]*/], dataToShow);
					if (plots[0][idf].Load(hydros, mbm, ids))
						nloaded++;
				}
				if (nloaded > 0)
					UpdateHead(bem, ih);
			}
		}
		
		if (nloaded == 0)
			return false;
		
		if (nloaded <= selTab)	
			selTab = 0;
		
		tab.Set(selTab);
		return true;
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
		return false;
	}
}
		
void MainABForce::UpdateHead(BEM &bem, int ih) {
	int it = max(0, tab.Get());
	{
		TempAssign<bool> _isFilling(isFilling, true);
		tab.Reset();
		for (int idf = 0; idf < 6*bem.Nb; ++idf) 
			tab.Add(plots[0][idf].SizePos(), Format("%s", BEM::StrBDOF(idf, false)));//, bem.headAll[ih/*bem.orderHeadAll[ih]*/]));
	}
	tab.Set(it);
}

void MainABForce::UpdateHeadMD(BEM &bem, int ih) {
	int it = max(0, tab.Get());
	{
		TempAssign<bool> _isFilling(isFilling, true);
		tab.Reset();
		//const std::complex<double> &h = bem.headAllMD[ih];
		for (int idf = 0; idf < 6*bem.Nb; ++idf) 
			tab.Add(plots[0][idf].SizePos(), Format("%s", BEM::StrBDOF(idf, false)));///*, h.real(), h.imag())*/);
	}
	tab.Set(it);
}