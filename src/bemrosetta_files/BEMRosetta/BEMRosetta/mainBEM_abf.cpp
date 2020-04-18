#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <TabBar/TabBar.h>

#include <BEMRosetta/BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;

#include "main.h"


void MainABForce::Init(DataToShow _dataToShow) {
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

bool MainABForce::Load(BEMData &bem, const Vector<int> &ids) {
	TempAssign<bool> _isFilling(isFilling, true);
	try {
		MainBEM &mbm = GetDefinedParent<MainBEM>(this);
	
		tab.Reset();
		Upp::Array<HydroClass> &hydros = bem.hydros; 
		if (hydros.IsEmpty() || ids.IsEmpty()) 
			return false;
		String format;
		switch (dataToShow) {
		case DATA_A:		format = t_("%s");		break;		
		case DATA_B:		format = t_("%s");		break;
		case DATA_FORCE_SC:	format = t_("%s%.1fº");	break;
		case DATA_FORCE_FK:	format = t_("%s%.1fº");	break;
		case DATA_FORCE_EX:	format = t_("%s%.1fº");	break;
		case DATA_RAO:		format = t_("%s%.1fº");	break;
		case DATA_STS:		NEVER();
		case DATA_STS2:		NEVER();
		}
		int sdof = 6*bem.Nb;
		if (dataToShow == DATA_A || dataToShow == DATA_B) {
			plots.SetCount(sdof);
			for (int idf = 0; idf < sdof; ++idf) {
				plots[idf].SetCount(sdof);
				for (int jdf = 0; jdf < sdof; ++jdf) {
					if (!bem.onlyDiagonal || idf == jdf) {
						plots[idf][jdf].Init(idf, jdf, dataToShow);
						if (plots[idf][jdf].Load(hydros, mbm, ids)) {
							if (idf != jdf)
								tab.Add(plots[idf][jdf].SizePos(), Format(format, Hydro::StrBDOF(idf, jdf)));
							else
								tab.Add(plots[idf][jdf].SizePos(), Format(format, Hydro::StrBDOF(idf)));
						}
					}
				}
			}
		} else {
			int Nh = bem.headAll.GetCount();
			if (Nh < 0) 
				return false;
			
			plots.SetCount(Nh);
			for (int ih = 0; ih < Nh; ++ih) 
				plots[ih].SetCount(sdof);
			for (int idf = 0; idf < sdof; ++idf) {
				for (int ih = 0; ih < Nh; ++ih) {
					plots[ih][idf].Init(idf, bem.headAll[ih], dataToShow);
					if (plots[ih][idf].Load(hydros, mbm, ids))
						tab.Add(plots[ih][idf].SizePos(), Format(format, Hydro::StrBDOFAbrev(idf), bem.headAll[ih]));
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
