// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <TabBar/TabBar.h>
#include <DropGrid/DropGrid.h>
#include <STEM4U/Utility.h>
#include <STEM4U/CrossCorrelation.h>

#include <BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;

#include "main.h"

UVector<CompareParameters *> CompareParameters::plist;
bool CompareParameters::comparing = false;

void CompareParameters::Init(ScatterDraw &scatter, SplitterButton &splitter) {
	CtrlLayout(*this);
	
	plist << this;
	
	pscatter = &scatter;
	psplitter = &splitter;
	
	splitter.WhenAction = [&] {
		NON_REENTRANT_V;
		if (comparing)			// To avoid recursive update of everything with everything
			return;
		comparing = true;
		int pos = psplitter->GetPos();
		for (int i = 0; i < plist.size(); ++i) {
			CompareParameters *c = plist[i];
			c->psplitter->SetPos(pos);	
		}
		comparing = false;
	};
	
	swRelative <<= 0;
	
	auto Up = [&] {
		NON_REENTRANT_V;
		for (CompareParameters *c : plist) {
			c->swRelative <<= ~swRelative;
			c->opMinRange <<= ~opMinRange;
			c->Load();
		}
	};
	
	swRelative.WhenAction = [=]() {Up();};
	opMinRange.WhenAction = [=]() {Up();};
}

CompareParameters::~CompareParameters() {
	int id = Find(plist, this);
	if (id >= 0)
		plist.Remove(id);
}

void CompareParameters::Init(Hydro::DataToShow _dataToShow) {
	this->dataToShow = _dataToShow;
	
	swRelative.Enable();
	opMinRange.Enable();
}

void CompareParameters::Load() {
	ScatterDraw &scatter = *pscatter;
	
	double wmin = std::numeric_limits<double>::min(), 
		   wmax = std::numeric_limits<double>::max();
	int idFirst = Null, idLast = Null;		// Ids first and last of the data with series
	
	for(int i = 0; i < scatter.GetCount(); i++) {
		DataSource &data = scatter.GetDataSource(i);
		if (data.GetCount() > 0) {
			String leg = scatter.GetLegend(i);
			if (leg.StartsWith("A∞(ω)") || !(leg.StartsWith("A∞") || leg.StartsWith("A₀"))) {
				wmin = max(data.MinX(), wmin);
				wmax = min(data.MaxX(), wmax);
				if (IsNull(idFirst))
					idFirst = i;
				idLast = i;
			}
		}
	}
	if (IsNull(idFirst) || IsNull(idLast))
		return;
	
	if (wmin > wmax)
		wmin = wmax;
	
	list.Reset();
	list.SetLineCy(EditField::GetStdHeight()).MultiSelect();
	list.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, list);};
	
	list.AddColumn(t_("Name"), 20);
	list.AddColumn(t_("RMS"), 10);
	if (swRelative > 0) {
		list.AddColumn(t_("RMSE"), 10);
		list.AddColumn(t_("Xcorr"), 10);
	}
	
	UVector<double> rms, rmse, xcorr, ainf, a0;
	UVector<String> leg;
	VectorXd x0, y0;		// Reference series

	if (swRelative > 0) {		// To fill x0, y0 if relative. Has to be done before next for
		VectorXd x, y;
		if (swRelative == 1) 
			scatter.GetDataSource(idFirst).CopyXY(x0, y0);
		else
			scatter.GetDataSource(idLast).CopyXY(x0, y0);
		if (opMinRange) 
			Segment(x0, y0, wmin, wmax, x0, y0);		
	}
	
	for(int i = 0; i < scatter.GetCount(); i++) {
		DataSource &data = scatter.GetDataSource(i);
		if (data.GetCount() == 0) 
			continue;
		
		String str = scatter.GetLegend(i);
		VectorXd x, y;
		data.CopyXY(x, y);
		
		if (str.StartsWith("A∞(ω)") || !(str.StartsWith("A∞") || str.StartsWith("A₀"))) {
			leg << str;
			if (opMinRange) 
				Segment(x, y, wmin, wmax, x, y);
			
			if ((swRelative == 1 && i == 0) || (swRelative == 2 && i == scatter.GetCount() - 1)) {
				rms << 1;
				rmse << 0;
				xcorr << 1;	
			} else {
				if (y.size() > 0)
					rms << RMS(y);
				
				if (swRelative > 0) {
					VectorXd yyy;
					ResampleY(x, y, x0, yyy);
					
					double mean = abs(y0.mean());
					if (mean > 0)
						rmse << RMSE(yyy, y0)/mean;
					else
						rmse << Null;
					
					VectorXd rr, lags;
					XCorr(yyy, y0, rr, lags, 'c');
					if (rr.size() == 0)
						xcorr << Null;
					else	
						xcorr << rr.maxCoeff();
				}
			}
		} else if (str.StartsWith("A∞")) {
			if (ainf.size() == 0)
				list.AddColumn(t_("A∞"), 10);
			ainf << y[0]; 
		} else if (str.StartsWith("A0")) {
			if (a0.size() == 0)
				list.AddColumn(t_("A0"), 10);
			a0 << y[0];
		}
	}	
	if (rms.size() == 0)
		return;
	
	for(int row = 0; row < rms.size(); row++) {
		int col = 0;
		list.Set(row, col++, leg[row]);	
		
		double _r = Null, _rm = Null, _x = Null, _ai = Null, _a0 = Null;
		
		if (swRelative == 0) {
			_r = rms[row];	
			if (ainf.size() > row) 
				_ai = ainf[row];		
			if (a0.size() > row) 
				_a0 = a0[row];
		} else if (swRelative == 1) {
			_r = First(rms) != 0 ? rms[row]/First(rms) : Null;	
			if (rmse.size() > row) 
				_rm = rmse[row];
			if (xcorr.size() > row) 
				_x = xcorr[row];
			if (ainf.size() > row) 
				_ai = First(ainf) != 0 ? ainf[row]/First(ainf) : Null;		
			if (a0.size() > row) 
				_a0 = First(a0) != 0 ? a0[row]/First(a0) : Null;
		} else {// if (swRelative == 2)
			_r = Last(rms) != 0 ? rms[row]/Last(rms) : Null;	
			if (rmse.size() > row) 
				_rm = rmse[row];
			if (xcorr.size() > row) 
				_x = xcorr[row];
			if (ainf.size() > row) 
				_ai = Last(ainf) != 0 ? ainf[row]/Last(ainf) : Null;		
			if (a0.size() > row) 
				_a0 = Last(a0) != 0 ? a0[row]/Last(a0) : Null;
		}
		
		if (!IsNull(_r))
			list.Set(row, col++, FDS(_r, 10, false));	
		if (!IsNull(_rm))
			list.Set(row, col++, FDS(_rm, 10, false));
		if (!IsNull(_x))
			list.Set(row, col++, FDS(_x, 10, false));
		if (!IsNull(_ai))
			list.Set(row, col++, FDS(_ai, 10, false));		
		if (!IsNull(_a0))
			list.Set(row, col++, FDS(_a0, 10, false));
	}
}

