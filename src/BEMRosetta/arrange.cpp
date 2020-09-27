#include <BEMRosetta_cl/BEMRosetta.h>
#include "arrange.h"


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
	listOrig.AddColumn(t_("Available"), 60).With([](One<Ctrl>& x) {
			x.Create<Option>().NoWantFocus().SetReadOnly();
		}
	);
	
	dofList.SetLineCy(EditField::GetStdHeight()).AutoHideSb();
	
	dofList.WhenDropInsert = [&](int line, PasteClip& d) { 
		if (DnDInsert(line, d)) {
			Upp::Vector<int> order;
			for (int i = 0; i < dofList.GetCount(); ++i) {
				int ib, idf;
				Hydro::DOFFromStr(dofList.Get(i, 0).ToString(), ib, idf); 
				order << idf + 6*ib;
			}
			Upp::Vector<int> neworder;
			for (int i = 0; i < order.GetCount(); ++i) {
				int id = FindIndex(order, i);
				neworder << id;
			}
			hydro.SetOrder(neworder);
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
	
	for (int i = 0; i < 6*hydro.Nb; ++i) 
		listOrig.Add(i+1, hydro.IsAvailableDOF(0, i));
		
	for (int i = 0; i < 6*hydro.Nb; ++i) {
		int id = Find(hydro.GetOrder(), i);
		dofList.Add(InitCaps(hydro.StrBDOF(id)));
	}
}

bool ArrangeDOF::DnDInsert(int line, PasteClip& d) {
	if(AcceptInternal<ArrayCtrl>(d, "array")) {
		dofList.InsertDrop(line, d);
		dofList.SetFocus();
		return true;
	}
	return false;
}

