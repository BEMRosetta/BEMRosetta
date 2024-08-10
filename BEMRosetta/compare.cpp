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
	if (opMinRange) {
		bool found = false;
		for(int i = 0; i < scatter.GetCount(); i++) {
			String leg = scatter.GetLegend(i);
			if (leg.StartsWith("A∞(ω)") || !(leg.StartsWith("A∞") || leg.StartsWith("A₀"))) {
				DataSource &data = scatter.GetDataSource(i);
				wmin = max(data.MinX(), wmin);
				wmax = min(data.MaxX(), wmax);
				found = true;
			}
		}
		if (!found)
			return;
	}
	if (wmin > wmax)
		wmin = wmax;
	
	list.Reset();
	list.SetLineCy(EditField::GetStdHeight()).MultiSelect();
	list.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, list);};
	
	list.AddColumn(t_("Name"), 20);
	list.AddColumn(t_("RMS"), 10);
	if (swRelative > 0) 
		list.AddColumn(t_("Xcorr"), 10);
	
	UVector<double> rms, xcorr, ainf, a0;
	UVector<String> leg;
	VectorXd x0, y0;

	if (swRelative > 0) {		// To fill x0, y0 if relative. Has to be done before next for
		for(int i = 0; i < scatter.GetCount(); i++) {
			String str = scatter.GetLegend(i);
			VectorXd x, y;
			scatter.GetDataSource(i).CopyXY(x, y);
				
			if (str.StartsWith("A∞(ω)") || !(str.StartsWith("A∞") || str.StartsWith("A₀"))) {
				VectorXd xx, yy;
				if (!opMinRange) {
					xx = x;
					yy = y;
				} else
					Segment(x, y, wmin, wmax, xx, yy);
				
				x0 = xx;
				y0 = yy;
				if (i == 0 && swRelative == 1) 
					break;
			}
		}
	}
	
	for(int i = 0; i < scatter.GetCount(); i++) {
		String str = scatter.GetLegend(i);
		VectorXd x, y;
		scatter.GetDataSource(i).CopyXY(x, y);
		
		if (x.size() == 0 || y.size() == 0)
			return;
	
		if (str.StartsWith("A∞(ω)") || !(str.StartsWith("A∞") || str.StartsWith("A₀"))) {
			leg << str;
			VectorXd xx, yy;
			if (!opMinRange) {
				xx = x;
				yy = y;
			} else
				Segment(x, y, wmin, wmax, xx, yy);
			
			if (yy.size() > 0)
				rms << RMS(yy);
			
			if (swRelative > 0) {
				VectorXd yyy;
				ResampleY(x, y, x0, yyy);
				VectorXd rr, lags;
				XCorr(yyy, y0, rr, lags, 'c');
				if (rr.size() == 0)
					xcorr << Null;
				else	
					xcorr << rr.maxCoeff();
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
		
		double _r, _x, _ai, _a0;
		
		if (swRelative == 0) {
			_r = rms[row];	
			if (ainf.size() > row+1) 
				_ai = ainf[row];		
			if (a0.size() > row+1) 
				_a0 = a0[row];
		} else if (swRelative == 1) {
			_r = First(rms) != 0 ? rms[row]/First(rms) : Null;	
			if (xcorr.size() > row+1) 
				_x = xcorr[row];
			if (ainf.size() > row+1) 
				_ai = First(ainf) != 0 ? ainf[row]/First(ainf) : Null;		
			if (a0.size() > row+1) 
				_a0 = First(a0) != 0 ? a0[row]/First(a0) : Null;
		} else {// if (swRelative == 2)
			_r = Last(rms) != 0 ? rms[row]/Last(rms) : Null;	
			if (xcorr.size() > row+1) 
				_x = xcorr[row];
			if (ainf.size() > row+1) 
				_ai = Last(ainf) != 0 ? ainf[row]/Last(ainf) : Null;		
			if (a0.size() > row+1) 
				_a0 = Last(a0) != 0 ? a0[row]/Last(a0) : Null;
		}
		
		list.Set(row, col++, FDS(_r, 10, false));	
		if (xcorr.size() > 0) 
			list.Set(row, col++, FDS(_x, 10, false));
		if (ainf.size() > 0) 
			list.Set(row, col++, FDS(_ai, 10, false));		
		if (a0.size() > 0) 
			list.Set(row, col++, FDS(_a0, 10, false));
	}
}

