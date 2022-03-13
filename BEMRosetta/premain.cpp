// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include <BEMRosetta_cl/BEMRosetta.h>
#include <GridCtrl/GridCtrl.h>
#include "premain.h"
#include "auxiliar.h"

void ArrangeDOF::Init(Hydro &hydro) {
	CtrlLayout(*this);
	
	selecting = false;
	
	listOrig.SetLineCy(EditField::GetStdHeight()).HideSb();
	
	listOrig.WhenSel = [this] {
		if (!selecting) {
			selecting = true;
			dofList.SetCursor(listOrig.GetCursor());
			selecting = false;
		}
	};
	listOrig.WhenScroll = [this] {dofList.ScrollTo(listOrig.GetScroll());};
	
	listOrig.AddColumn(t_("DOF"), 40);
	listOrig.AddColumn(t_("Available"), 40).With([](One<Ctrl>& x) {
			x.Create<Option>().NoWantFocus().SetReadOnly();
		}
	);
	
	dofList.SetLineCy(EditField::GetStdHeight()).AutoHideSb();
	
	dofList.WhenDropInsert = [&](int line, PasteClip& d) { 
		if (DnDInsert(line, d)) {
			Upp::Vector<int> order;
			for (int i = 0; i < dofList.GetCount(); ++i) {
				int ib, idf;
				BEM::DOFFromStr(dofList.Get(i, 0).ToString(), ib, idf); 
				order << (idf + 6*ib);
			}
			/*Upp::Vector<int> neworder;
			for (int i = 0; i < order.size(); ++i) {
				int id = FindIndex(order, i);
				neworder << id;
			}*/
			hydro.SetOrder(order);
		}
	};
	dofList.WhenDrag = [=] { 
		if(dofList.DoDragAndDrop(InternalClip(dofList, "array")) == DND_MOVE)
			dofList.RemoveSelection();
	};
	dofList.WhenSel = [this] {
		if (!selecting) {
			selecting = true;
			listOrig.SetCursor(dofList.GetCursor());
			selecting = false;
		}
	};
	dofList.WhenScroll = [this] {listOrig.ScrollTo(dofList.GetScroll());};

	dofList.AddColumn(t_("Arrange DOF"));
	
	for (int ib = 0; ib < hydro.Nb; ++ib) {
		for (int idf = 0; idf < 6; ++idf) 
			listOrig.Add(InitCaps(BEM::StrBDOF(ib*6+idf, false)), hydro.IsAvailableDOF(ib, idf));
	}
	for (int i = 0; i < 6*hydro.Nb; ++i) {
		int id = Find(hydro.GetOrder(), i);
		dofList.Add(InitCaps(BEM::StrBDOF(id, false)));
	}
	colorMark.Color(GetColorId(hydro.GetId()));
}

bool ArrangeDOF::DnDInsert(int line, PasteClip& d) {
	if(AcceptInternal<ArrayCtrl>(d, "array")) {
		dofList.InsertDrop(line, d);
		dofList.SetFocus();
		return true;
	}
	return false;
}

