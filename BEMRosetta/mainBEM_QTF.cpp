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
#include <BEMRosetta_cl/functions.h>

using namespace Upp;

#include "main.h"


void QTFTabDof::Init(MainQTF &par, int posSplitter, int ib, int idof) {
	this->ib = ib;
	this->idof = idof;
	Add(splitter);
	splitter.Horz(leftsplit.SizePos(), rightsplit.SizePos());
	splitter.SetPos(posSplitter, 0);
	leftsplit.Add(up.array, 0, 0).Add(down.array, 1, 0);
	rightsplit.Add(up.sc, 0, 0).Add(down.sc, 1, 0);
	
	up.sc.Add(up.surf, 0, 0).Add(up.scatter, 0, 1);
	down.sc.Add(down.surf, 0, 0).Add(down.scatter, 0, 1);
	up.sc.WhenWidths = down.sc.WhenWidths = [&](int width, int height, UVector<int> &widths) {
		widths[0] = height;
		widths[1] = max(0, width - height);
	};
	
	up.isUp = true;
	down.isUp = false;
		
	up	.surf.ShowInfo().ShowContextMenu().ShowPropertiesDlg().ShowProcessDlg().SetLeftMargin(50).SetTopMargin(25).SetBottomMargin(50);
	down.surf.ShowInfo().ShowContextMenu().ShowPropertiesDlg().ShowProcessDlg().SetLeftMargin(50).SetTopMargin(25).SetBottomMargin(50);
	up	.surf.LinkedWith(down.surf);
	
	up  .surf.WhenPainter = THISBACK(OnPainter);
	down.surf.WhenPainter = THISBACK(OnPainter);
	up  .surf.WhenDraw    = THISBACK(OnDraw);
	down.surf.WhenDraw    = THISBACK(OnDraw);
	
	up  .surf.WhenMouseClick = [&](Point p, dword keyflags, ScatterCtrl::MouseAction action) {OnClick(p, this->idof, action);};
	down.surf.WhenMouseClick = [&](Point p, dword keyflags, ScatterCtrl::MouseAction action) {OnClick(p, this->idof, action);};
	
	int len = StdFont().GetHeight();
	
	up  .surf.SetMargin(4*len, len, int(2.5*len), 4*len);
	down.surf.SetMargin(4*len, len, int(2.5*len), 4*len);
	
	up.  scatter.SetMargin(6*len, len, len, 4*len).SetTitleFont(SansSerifZ(12)).ShowAllMenus();
	down.scatter.SetMargin(6*len, len, len, 4*len).SetTitleFont(SansSerifZ(12)).ShowAllMenus();		   
	up.scatter.LinkedWith(down.scatter);
	
	parent = &par;
}

double QTFTabDof::qwT(const Hydro &hd, int id) const {
	return show_w ? hd.qw[id] : 2*M_PI/hd.qw[id];
}

VectorXd QTFTabDof::qwT(const Hydro &hd) const {
	VectorXd ret(hd.qw.size());
	for (int i = 0; i < hd.qw.size(); ++i)
		ret[i] = show_w ? hd.qw[i] : 2*M_PI/hd.qw[i];
	return ret;
}

Pointf &QTFTabDof::Pf() {
	return parent->pf;
}

void QTFTabDof::DoClick(Data &up, int idof) {
	up.dataPlot.Clear();
	up.scatter.RemoveAllSeries();
	
	String strmag;
	
	if (up.show_ma_ph) {
		if (up.isUp) {
			up.labelY = t_("Magnitude");
			up.ma_ph = t_("ma");
			strmag = t_("mag");
			up.units = idof < 3 ? t_("N/m²") : t_("N m/m²");
		} else {
			up.labelY = t_("Phase");
			up.ma_ph = t_("ph");
			strmag = t_("phase");
			up.units = t_("rad");
		}
	} else {
		if (up.isUp) {
			up.labelY = t_("Real");
			up.ma_ph = t_("re");
			strmag = t_("real");
		} else {
			up.labelY = t_("Imaginary");
			up.ma_ph = t_("im");
			strmag = t_("imag");
		}
		up.units = idof < 3 ? t_("N/m²") : t_("N m/m²");
	}
	up.scatter.SetLabelY(up.labelY);
	String saxis;
	if (typec == 'v')
		saxis = t_("Y axis");
	else
		saxis = t_("X axis");
	up.scatter.SetLabelX(Format("%s %s", saxis, show_w ? t_("ω [rad/s]") : t_("T [s]")));
	
	double avgT = 0;
	for (int i = 0; i < Bem().hydros.size(); ++i) {
		const Hydro &hd = Bem().hydros[i].hd();
		if (!hd.IsLoadedQTF(isSum))
			continue;
		 
		int idh = FindDelta(hd.qh, head, 2.);
		if (idh < 0) 
			continue;
		
		VectorXd Qw = qwT(hd);
		if (!show_w)
			Upp::Reverse(Qw);
		
		MatrixXd m = GetMat(hd, up, idh, show_w);
		
		UArray<Pointf> &d = up.dataPlot.Add();	
		
		int id = Null;
		if (IsNull(Pf())) {
			double freq = Avg(Last(Qw), First(Qw));
			id = FindClosest(Qw, freq);
			Pf().x = Pf().y = Qw(id);
		} else {
			if (snap) {
				if (typec == 'h') {
					id = FindClosest(Qw, Pf().y);
					Pf().y = Qw(id);
				} else if (typec == 'v') {			
					id = FindClosest(Qw, Pf().x);
					Pf().x = Qw(id);
				}
			}
		}
		Pointf from, to;
		
		if (typec == 'h' || typec == 'v') {
			if (snap)
				avgT += Qw(id);
			else
				avgT += Pf().x;
		} else if (typec == 'd') 
			Diagonal (Pf(), First(Qw), Last(Qw), from, to);
		else
			Conjugate(Pf(), First(Qw), Last(Qw), from, to);
		
		if (typec == 'h' || typec == 'v') {
			if (snap) {
				for (int iw = 0; iw < hd.qw.size(); ++iw) {
					double val;
					if (typec == 'h') 
						val = m(id, iw);
					else 
						val = m(iw, id);
					d << Pointf(Qw(iw), val);
				}	
			} else {
				if (typec == 'h') {
					from = Pointf(First(Qw), Pf().y);
					to = Pointf(Last(Qw), Pf().y);
				} else {
					from = Pointf(First(Qw), Pf().x);
					to = Pointf(Last(Qw), Pf().x);
				}
				if (!IsNull(from) && !IsNull(to)) {
					double dx = (to.x - from.x)/300;
					double dy = (to.y - from.y)/300;
					for (int i = 0; i <= 300; ++i)
						d << Pointf(from.x + i*dx, BilinearInterpolate(from.x + i*dx, from.y + i*dy, Qw, Qw, m));
				}	
			}
		} else {
			if (!IsNull(from) && !IsNull(to)) {
				double dx = (to.x - from.x)/300;
				double dy = (to.y - from.y)/300;
				for (int i = 0; i <= 300; ++i)
					d << Pointf(from.x + i*dx, BilinearInterpolate(from.x + i*dx, from.y + i*dy, Qw, Qw, m));
			}
		}
		int idc = hd.GetId();
		const Upp::Color &color = GetColorId(idc);
		String nameType = Format(t_("QTF %s %s(%s)"), up.ma_ph, hd.name, hd.GetCodeStrAbr());
		up.scatter.AddSeries(d).Legend(nameType).Units(up.units).SetMarkColor(color).Stroke(2, color);
		if (!showPoints)
			up.scatter.NoMark();
	}
	if (typec == 'h' || typec == 'v')
		avgT /= Bem().hydros.size();		// Average value
	
	String strw;
	if (typec == 'd')
		strw = t_("Diagonal");
	else if (typec == 'c')
		strw = t_("Conjugate");
	else 
		strw = Format("%.2f %s", avgT, show_w ? "rad/s" : "s");
	up.scatter.SetTitle(Format(t_("QTF %s %d.%s %s heading %.1f:%.1fº %s"), isSum ? "sum" : "dif", ib+1, BEM::StrDOF(idof), strw, real(head), imag(head), strmag));
	
	if (autoFit) {
		up.scatter.ZoomToFit(true, true);
		if (up.isUp || !up.show_ma_ph) {
			if (fromY0) {
				double yRange = max<double>(0, up.scatter.GetYMin()) + up.scatter.GetYRange();
				up.scatter.SetXYMin(Null, 0).SetRange(Null, yRange);
			}
		} else {
			if (up.show_ma_ph && !up.isUp) {
				up.scatter.ZoomToFit(true, false);
				up.scatter.SetXYMin(Null, -M_PI).SetRange(Null, 2*M_PI).SetMajorUnits(Null, 1);
				up.scatter.SetMinUnits(Null, M_PI-3);
			} else if (fromY0) {
				double yRange = max<double>(0, up.scatter.GetYMin()) + up.scatter.GetYRange();
				up.scatter.SetXYMin(Null, 0).SetRange(Null, yRange);
			} 
		}
	}
	up.scatter.Refresh();
}

void QTFTabDof::OnClick(Point p, int idof, ScatterCtrl::MouseAction action) {
	if (action != ScatterCtrl::LEFT_DOWN && action != ScatterCtrl::LEFT_MOVE)
		return;
	
	Pf().x = up.surf.GetRealPosX(p.x);
	Pf().y = up.surf.GetRealPosY(p.y);
	
	up.surf.Refresh();
	down.surf.Refresh();
	
	DoClick(up, idof);
	DoClick(down, idof);
}

char QTFTabDof::GetWhat(const Data &data) {
	if (data.show_ma_ph) {
		if (data.isUp)
			return 'm';
		else
			return 'p';
	} else {
		if (data.isUp)
			return 'r';
		else 
			return 'i';
	}
}
	
double QTFTabDof::GetData(const Hydro &hd, const Data &data, int idh, int ifr1, int ifr2) const {
	return hd.GetQTFVal(ib, idof, idh, ifr1, ifr2, isSum, GetWhat(data));
}

MatrixXd QTFTabDof::GetMat(const Hydro &hd, const Data &data, int idh, bool show_w) const {
	MatrixXd m = hd.GetQTFMat(ib, idof, idh, isSum, GetWhat(data));
	if (!show_w)
		ReverseX(m);
	return m;
}
			
void QTFTabDof::UpdateArray(const Hydro &hd, bool show_ma_ph, Data &data, bool opBilinear, bool opSnap) {
	data.show_ma_ph = show_ma_ph;
	
	int qtfNf = int(hd.qw.size());

	ArrayCtrl &up = data.array;
	
	up.Reset();
	up.SetLineCy(EditField::GetStdHeight()).HeaderObject().Absolute();
	up.MultiSelect();
	up.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, up);};

	up.AddColumn(show_w ? t_("ω [rad/s]") : t_("T [s]"), 60);
	for (int c = 0; c < qtfNf; ++c)
		up.AddColumn(FDS(show_w ? hd.qw[c] : 2*M_PI/hd.qw[qtfNf-1-c], 8), 90);
	for (int r = qtfNf-1; r >= 0; --r)
		up.Add		(FDS(show_w ? hd.qw[r] : 2*M_PI/hd.qw[qtfNf-1-r], 8));
	
	data.xAxis.SetCount(qtfNf);
	for (int i = 0; i < qtfNf; ++i) 
		data.xAxis[i] = (show_w ? hd.qw[i] : 2*M_PI/hd.qw[qtfNf-1-i]);
	
	
	MatrixXd mat = GetMat(hd, data, ih, show_w);
	double mn = mat.minCoeff(), mx = mat.maxCoeff();
	if (mx == mn)
		mx = Null;
			
	data.zData.Clear();
	data.zData.Reserve(qtfNf*qtfNf);

	for (int if1 = 0; if1 < qtfNf; ++if1) {
		for (int if2 = 0; if2 < qtfNf; ++if2) {
			double val = mat(if1, if2);
			if (IsNull(val)) {
				up.Set(if2, 1+if1, "-");
				data.zData << Null;
			} else {				
				data.zData << val;
				
				if (show_ma_ph && !data.isUp) 
					up.Set(if2, 1+if1, FDS(val, 10, false));
				else {
					if (IsNull(mx)) 
						up.Set(qtfNf-1-if2, 1+if1, FDS(val, 10, false));
					else {
						double rat = (val - mn)/(mx - mn);
						
						::Color backColor = GetRainbowColor(rat, White(), LtBlue(), 0);
						::Color color = Black();
						if (Grayscale(backColor) < 150)
							color = White();
						
						String str = FDS(val, 10, false);
						
						up.Set(qtfNf-1-if2, 1+if1, AttrText(str).Center().Ink(color).Paper(backColor));
					}
				}
			}
		}
	}
	if (data.xAxis[0] > data.xAxis[1]) {
		ReverseX(data.xAxis);
		
		for (int r = 0; r < qtfNf/2; ++r)
			for (int c = 0; c < qtfNf; ++c)
				Swap(data.zData[r*qtfNf + c], data.zData[(qtfNf-1-r)*qtfNf + (qtfNf-1-c)]);
		if (Odd(qtfNf)) {
			int r = qtfNf/2;
			for (int c = 0; c < qtfNf/2; ++c)
				Swap(data.zData[r*qtfNf + c], data.zData[r*qtfNf + (qtfNf-1-c)]);
		}
	}
	
	data.dataSurf.Init(data.zData, data.xAxis, data.xAxis, opBilinear ? TableInterpolate::BILINEAR : TableInterpolate::NO, false);

	data.surf.AddSurf(data.dataSurf);
	data.surf.SetRainbowPaletteTextColor(White);
	data.surf.ZoomToFitZ().ZoomToFit(true, true);
}
	
void QTFTabDof::Load(const Hydro &hd, int ib, int ih, int idof, bool ndim, bool show_w, bool show_ma_ph, bool isSum, bool opBilinear, bool opSnap, bool showPoints, bool fromY0, bool autoFit, int posSplitter, bool resetPf) {
	try {
		splitter.SetPos(posSplitter, 0);
		
		const UArray<UArray<UArray<MatrixXcd>>> &qtf = isSum ? hd.qtfsum : hd.qtfdif;
		if (qtf.size() <= ib || qtf[ib].size() <= ih || qtf[ib][ih].size() <= idof)
			return;
		
		this->isSum = isSum;
		this->ib = ib;
		this->ih = ih;
		this->idof = idof;
		this->ndim = ndim;
		this->show_w = show_w;
		this->head = hd.qh[ih];
		this->showPoints = showPoints;
		this->fromY0 = fromY0;
		this->autoFit = autoFit;
		this->snap = opSnap;
		
		if (resetPf)
			Pf() = Null;
		
		UpdateArray(hd, show_ma_ph, up, opBilinear, opSnap);
		UpdateArray(hd, show_ma_ph, down, opBilinear, opSnap);
		DoClick(up, idof);
		DoClick(down, idof);
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}	
	
void MainQTF::Init(MainBEM &parent) {
	CtrlLayout(*this);
	
	try {
		_mbm = &parent;
		ArrayCtrl &headQTF = parent.menuPlot.headQTF;
		
		headQTF.Reset();
		headQTF.NoHeader();
		headQTF.AddColumn("", 20);
		headQTF.AddColumn("", 20);
		
		opLine <<= 0;
		
		opQTF  		<< [&] {
			isSumm = opQTF.GetData() == FSUM;
			OnHeadingsSel(&headQTF, false);
		};
		opBilinear  	<< THISBACK(OnSurf);
		opSnap 			<< THISBACK2(OnHeadingsSel, &headQTF, true);
		opLine 			<< THISBACK2(OnHeadingsSel, &headQTF, true);
		
		headQTF.WhenSel << THISBACK2(OnHeadingsSel, &headQTF, false);
		tab.WhenSet 	<< THISBACK2(OnHeadingsSel, &headQTF, false);
		
		opSnap = true;
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}		
}

void MainQTF::OnHeadingsSel(ArrayCtrl *headQTF, bool resetPf) {
	if (isLoading)
		return;
	
	int ih = headQTF->GetCursor();
	if (ih < 0)
		return;

	Unload(idof);
	
	try {
		WaitCursor wait;
		
		opSnap.Enable(opLine > 1);
			
		MainBEM &mbm = *_mbm;

		int idHydro = mbm.GetIdOneSelected(false);
		if (idHydro < 0)
			return;
	
		const Hydro &hd = Bem().hydros[idHydro].hd();
	
		head = hd.qh[ih];
		bool ndim = mbm.menuPlot.showNdim;
		bool show_w = mbm.menuPlot.opwT == 0;
		bool show_ma_ph = mbm.menuPlot.opMP == 0;
		bool showPoints = mbm.menuPlot.showPoints;
		bool fromY0 = mbm.menuPlot.fromY0;
		bool autoFit = mbm.menuPlot.autoFit;
		bool isSum = opQTF.GetData() == FSUM;
		
		idof = tab.Get();
		if (idof < 0)
			return;
		
		ib = idof/6;
		idof = idof - 6*ib;

		OnSurf();
		dof[idof+6*ib].Load(hd, ib, ih, idof, ndim, show_w, show_ma_ph, isSum, ~opBilinear, ~opSnap, showPoints, fromY0, autoFit, posSplitter, resetPf);
		
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
}

void MainQTF::OnSurf() {
	switch (int(~opLine)) {
	case 0:	dof[idof+6*ib].typec = 'd';	break;
	case 1:	dof[idof+6*ib].typec = 'c';	break;
	case 2:	dof[idof+6*ib].typec = 'h';	break;
	default:dof[idof+6*ib].typec = 'v';
	}
	dof[idof+6*ib].up.dataSurf.SetInterpolate(opBilinear ? TableInterpolate::BILINEAR : TableInterpolate::NO);
	dof[idof+6*ib].up.surf.Refresh();
	dof[idof+6*ib].down.dataSurf.SetInterpolate(opBilinear ? TableInterpolate::BILINEAR : TableInterpolate::NO);
	dof[idof+6*ib].down.surf.Refresh();
}
	
bool MainQTF::Load() {
	try {	
		MainBEM &mbm = *_mbm;
		
		int idHydro;
		
		{
			TempAssign<bool> _isLoading(isLoading, true);

			tab.Reset();

			//idHydro = mbm.GetIdOneSelected(false);
			//if (idHydro < 0) 
			//	return false;
						
			dof.SetCount(6*Bem().Nb);
			for (int ib = 0; ib < Bem().Nb; ++ib) {
				for (int idf = 0; idf < 6; ++idf) {
					dof[idf+6*ib].Init(*this, posSplitter, ib, idf);
					tab.Add(dof[idf+6*ib].SizePos(), Format("%d.%s", ib+1, BEM::StrDOF(idf)));
				}
			}
			if (tab.GetCount() >= idof + 6*ib && idof >= 0)
				tab.Set(idof + 6*ib);
		}
			
		ArrayCtrl &headQTF = mbm.menuPlot.headQTF;
		
		headQTF.Clear();
	
		idHydro = -1;
		for (int row = 0; row < mbm.listLoaded.GetCount(); ++row) {
			if (mbm.listLoaded.IsSelected(row)) {
				idHydro = ArrayModel_IdHydro(mbm.listLoaded, row);
				break;
			}
		}	// Only one available => directly selected
		if (idHydro < 0 && mbm.listLoaded.GetCount() == 1)
			idHydro = ArrayModel_IdHydro(mbm.listLoaded, 0);
		//if (idHydro < 0) 
		//	return false;
		
		//if (ArrayCtrlSelectedGetCount(mbm.listLoaded) > 1) 
		//	return false;
		
		// Show the tab if any model has QTFs
		bool show = false;
		for (const HydroClass &h : Bem().hydros) {
			if (h.hd().IsLoadedQTF(true) || h.hd().IsLoadedQTF(false)) {
				show = true;
				break;
			}
		}
		if (!show)
			return false;
		////
		
		if (idHydro >= 0) {
			const Hydro &hd = Bem().hydros[idHydro].hd();
			
			opQTF.Clear();
			if (hd.IsLoadedQTF(true))
				opQTF.Add(FSUM, t_("Summation"));
			else 
				isSumm = false;
			if (hd.IsLoadedQTF(false))
				opQTF.Add(FDIFFERENCE, t_("Difference"));
			else
				isSumm = true;
			if (opQTF.GetCount() > 1)
				opQTF.SetIndex(isSumm ? FSUM : FDIFFERENCE);
			else if (opQTF.GetCount() > 0)
				opQTF.SetIndex(0);
				
			/*UArray<std::complex<double>> qh;					// Prepare qtf headings to be shown ordered
			for (const auto &c : hd.qh)
				qh << FixHeading(c, Bem().headingType);
			
			Sort(qh, [&](auto& a, auto& b)->bool const { 
				if (a.real() < b.real())
					return true; 
				else if (a.real() > b.real())
					return false;
				else
					return a.imag() < b.imag();	
			});*/
		
			for (int ih = 0; ih < hd.qh.size(); ++ih)
				headQTF.Add(hd.qh[ih].real(), hd.qh[ih].imag());
					
			if (headQTF.GetCount() > 0) {
				int id = FindClosest(hd.qh, head);
				if (id < 0)
					id = 0;
				headQTF.SetCursor(id);
			}
			mbm.menuPlotList.SetQTF();
		}
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}				
	return true;		
}

void MainQTF::Unload(int idf) {
	if (ib < 0)
		return;
	
	if (idf < 0) {
		idf = tab.Get();
		if (idf < 0)
			return;
	}
	
	posSplitter = dof[idf+6*ib].splitter.GetPos(0);
}
