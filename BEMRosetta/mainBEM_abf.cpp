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

#include "main.h"


void MainABForce::Init(Hydro::DataToShow _dataToShow) {
	Add(tab.SizePos());
	
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

bool MainABForce::Load(BEM &bem, const Upp::Vector<int> &ids) {
	TempAssign<bool> _isFilling(isFilling, true);
	try {
		MainBEM &mbm = GetDefinedParent<MainBEM>(this);
	
		tab.Reset();
		Upp::Array<HydroClass> &hydros = bem.hydros; 
		if (hydros.IsEmpty() || ids.IsEmpty()) 
			return false;
		String format;
		switch (dataToShow) {
		case Hydro::DATA_A:			
		case Hydro::DATA_B:
		case Hydro::DATA_AINFW:
		case Hydro::DATA_K:			format = t_("%s");		break;
		case Hydro::DATA_FORCE_SC:	
		case Hydro::DATA_FORCE_FK:	
		case Hydro::DATA_FORCE_EX:	
		case Hydro::DATA_RAO:		format = t_("%s%.1fº");	break;
		case Hydro::DATA_STS:		NEVER();
		case Hydro::DATA_STS2:		NEVER();
		}
		int sdof = 6*bem.Nb;
		if (dataToShow == Hydro::DATA_A || dataToShow == Hydro::DATA_B || dataToShow == Hydro::DATA_AINFW || 
			dataToShow == Hydro::DATA_K) {
			plots.SetCount(sdof);
			for (int idf = 0; idf < sdof; ++idf) {
				plots[idf].SetCount(sdof);
				for (int jdf = 0; jdf < sdof; ++jdf) {
					if (!bem.onlyDiagonal || idf == jdf) {
						plots[idf][jdf].Init(idf, jdf, dataToShow);
						if (plots[idf][jdf].Load(hydros, mbm, ids)) {
							if (idf != jdf)
								tab.Add(plots[idf][jdf].SizePos(), Format(format, BEM::StrBDOF2(idf, jdf, false)));
							else
								tab.Add(plots[idf][jdf].SizePos(), Format(format, BEM::StrBDOF(idf, false)));
						}
					}
				}
			}
		} else {
			int Nh = bem.headAll.size();
			if (Nh < 0) 
				return false;
			
			plots.SetCount(Nh);
			for (int ih = 0; ih < Nh; ++ih) 
				plots[ih].SetCount(sdof);
			for (int idf = 0; idf < sdof; ++idf) {
				for (int ih = 0; ih < Nh; ++ih) {
					plots[ih][idf].Init(idf, bem.headAll[ih], dataToShow);
					if (plots[ih][idf].Load(hydros, mbm, ids))
						tab.Add(plots[ih][idf].SizePos(), Format(format, BEM::StrBDOF(idf, false), bem.headAll[ih]));
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
