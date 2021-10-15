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

using namespace Eigen;



void MainStiffness::Init() {
	CtrlLayout(*this);
	
	numDecimals <<= THISBACK(PrintData);
	opEmptyZero <<= THISBACK(PrintData);
}

void MainStiffness::Clear() {
	array.Reset();
	array.NoHeader().SetLineCy(EditField::GetStdHeight()).HeaderObject().Absolute();
	array.MultiSelect().SpanWideCells();
	array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, array);};
	
	data.Clear();
	row0s.Clear();
	col0s.Clear();
}

void MainStiffness::AddPrepare(int &row0, int &col0, String name, int icase, String bodyName, int ibody, int idc) {
	row0s << (row0 = ibody*9 + 1);
	col0s << (col0 = icase*8);
	
	array.Set(0, col0, AttrText(name).Bold());
	array.Set(2, col0, AttrText(" ").Paper(GetColorId(idc)));
	
	while (array.GetColumnCount() < col0 + 7) {
		if (icase > 0) {
			array.AddColumn("", 10);
			int icol = array.GetColumnCount() - 1;
			array.HeaderObject().Tab(icol).SetMargin(0);
		}
		array.AddColumn("", 20);
		for (int i = 0; i < 2; ++i)
			array.AddColumn("", 20);	
		for (int i = 0; i < 4; ++i)
			array.AddColumn("", 70);
	}
	array.Set(row0, col0, AttrText(Format(t_("#%d body. %s"), ibody + 1, bodyName)).Bold());
	for (int i = 0; i < 6; ++i) {
		array.Set(row0 + 1, 	col0 + i + 1, AttrText(FormatInt(i + 1)).Bold().Align(ALIGN_CENTER));
		array.Set(row0 + i + 2, col0, 	  	  AttrText(FormatInt(i + 1)).Bold().Align(ALIGN_CENTER));
	}
	for (int i = 1; i <= ibody; ++i) 
		array.SetLineCy(i*9, 10);	
}

void MainStiffness::PrintData() {
	if (IsNull(numDecimals) || numDecimals < 0 || numDecimals > 14)
		numDecimals <<= 0;
	for (int i = 0; i < data.size(); ++i) {
		int row0 = row0s[i];
		int col0 = col0s[i];
		if (data[i].size() != 0) {
			for (int r = 0; r < 6; ++r) {
				for (int c = 0; c < 6; ++c) {
					double val = data[i](r, c);
					String data = FormatF(val, numDecimals);
					if (ScanDouble(data) == 0)
						data = opEmptyZero ? "" : "0";
					array.Set(row0 + r + 2, col0 + c + 1, AttrText(data).Align(ALIGN_RIGHT));
				}
			}
		}
	}
}

void MainStiffness::Add(const Mesh &mesh, int icase, bool button) {
	String name = mesh.fileName;
	
	const MatrixXd &K = mesh.C;
	data << clone(K);
	
	int idc = mesh.GetId();
	
	int row0, col0;
	AddPrepare(row0, col0, name, icase, "", 0, idc);
	
	if (K.size() == 0)
		return;

	if (button) {
		array.CreateCtrl<Button>(row0, col0+5, false).SetLabel(t_("Save")).Tip(t_("Saves to Wamit .hst stiffness matrix format"))
			<< [&] {
				FileSel fs;
				fs.Type(t_("Wamit stiffness matrix format"), "*.hst");
				if (fs.ExecuteSaveAs(t_("Save to Wamit .hst stiffness matrix format"))) 
					static_cast<const WamitMesh &>(mesh).SaveHST(~fs, Bem().rho, Bem().g);
			};
	}
	if (button && Bem().hydros.size() > 0) {
		array.CreateCtrl<Button>(row0, col0+6, false).SetLabel(t_("Copy")).Tip(t_("Copies matrix and paste it in selected BEM Coefficients file and body"))
			<< [=] {
				WithBEMList<TopWindow> w;
				CtrlLayout(w);
				w.Title(t_("Copies matrix and paste it in selected BEM Coefficients file and body"));
				w.array.SetLineCy(EditField::GetStdHeight());
				w.array.AddColumn(t_("File"), 20);
				w.array.AddColumn(t_("Body"), 10);
				w.array.HeaderObject().HideTab(w.array.AddColumn().HeaderTab().GetIndex());
				w.array.HeaderObject().HideTab(w.array.AddColumn().HeaderTab().GetIndex());
				for (int f = 0; f < Bem().hydros.size(); ++f) {
					const Hydro &hy = Bem().hydros[f].hd();
					for (int ib = 0; ib < hy.Nb; ++ib)
						w.array.Add(hy.name, hy.names[ib].IsEmpty() ? AsString(ib+1) : hy.names[ib], f, ib);
				}
				bool cancel = true;
				w.butSelect << [&] {cancel = false;	w.Close();};
				w.butCancel << [&] {w.Close();}; 
				w.Execute();
				if (!cancel) {
					int id = w.array.GetCursor();
					if (id < 0)
						return; 
					int f = w.array.Get(id, 2);
					int ib = w.array.Get(id, 3);
					Bem().hydros[f].hd().SetC(ib, K);
				}
			 };
	}
}
	
void MainStiffness::Add(String name, int icase, String bodyName, int ibody, const Hydro &hydro, int idc) {
	data << hydro.C_dim(ibody);
	
	int row0, col0;
	AddPrepare(row0, col0, name, icase, bodyName, ibody, idc);

	if (hydro.C.IsEmpty() || hydro.C[ibody].size() == 0)
		return;
}

bool MainStiffness::Load(Upp::Array<HydroClass> &hydros, const Upp::Vector<int> &ids) {
	Clear();
	
	for (int i = 0; i < ids.size(); ++i) {
		int isurf = ids[i];
		Hydro &hydro = hydros[isurf].hd();
		for (int ibody = 0; ibody < hydro.Nb; ++ibody) 
			Add(hydro.name, i, hydro.names[ibody], ibody, hydro, hydro.GetId());
	}
	PrintData();
	return true;
}	

void MainStiffness::Load(Upp::Array<Mesh> &surfs, const Upp::Vector<int> &ids) {
	Clear();

	for (int i = 0; i < ids.size(); ++i) {
		int isurf = ids[i];
		if (isurf >= 0)	
			Add(surfs[isurf], i, true);
	}
	PrintData();
}