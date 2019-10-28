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
		case DATA_STS:		NEVER();
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
