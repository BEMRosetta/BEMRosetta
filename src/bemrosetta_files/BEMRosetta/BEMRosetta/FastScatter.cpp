#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <RasterPlayer/RasterPlayer.h>
#include <TabBar/TabBar.h>

using namespace Upp;

#include <BEMRosetta/BEMRosetta_cl/FastOut.h>
#include "FastScatter.h"

#include "clip.brc"

void FastScatter::Init() {
	CtrlLayout(*this);
	
	CtrlLayout(left);
	CtrlLayout(right);
	splitter.Horz(left.SizePos(), right.SizePos());

	CtrlLayout(leftSearch);
	CtrlLayout(rightSearch);
	right.splitterSearch.Horz(leftSearch.SizePos(), rightSearch.SizePos());
	
	file.WhenChange = THISBACK(OnLoad);
	file.BrowseRightWidth(40).UseOpenFolder(true).BrowseOpenFolderWidth(10);
	butLoad.Tip(t_("Loads FAST out/outb file")).WhenAction = [&] {file.DoGo();};
	file.Type(t_("FAST output file"), "*.out, *.outb"); 
	
	left.scatter.ShowAllMenus().SetMode(ScatterDraw::MD_DRAW);
	
	leftSearch.array.AddColumn("");
	rightSearch.array.AddColumn("");
	rightSearch.label.SetLabel(t_("Right axis"));

	right.filter.WhenAction = THISBACK(OnFilter);
	right.paramList.AddColumn("");
	
	player.LoadBuffer(String(animatedStar, animatedStar_length));
	player.WhenAction = [&] {
		if (!player.IsRunning()) {
		 	player.Play();
			updateTime.WhenEnter();		
		} else {
			player.Stop();
		}
	};
	updateTime <<= "1:00";
	updateTime.WhenEnter = [&] {
		int seconds = int(StringToSeconds(~updateTime));
		updateTime <<= SecondsToString(seconds, 0, false, false, true, false, true);
		updateTime.CancelSelection();
	};
}

void FastScatter::OnFilter() {
	Vector<String> list = datafast.GetParameterList(~right.filter);
	right.paramList.Clear();
	for (int rw = 0; rw < list.GetCount(); ++rw)
		right.paramList.Add(list[rw]);
}

bool FastScatter::OnLoad() {
	try {
		if (!datafast.Load(~file)) {
			statusBar.Temporary(Format("File '%s' not found", ~file));
			return false;
		}
		
		WaitCursor waitcursor;
		statusBar.Temporary(t_("Loading file"));
		
		left.scatter.RemoveAllSeries();

		//dataSource.SetCount(datafast.parameters.GetCount());
		//for (int c = 0; c < dataSource.GetCount(); ++c) 
		//	dataSource[c].Init(datafast, c);
		
		//ShowSelected();
		OnFilter();
	} catch (Exc e) {
		Exclamation(Format("Error: %s", DeQtf(e)));	
	}
	return true;
}

void FastScatter::ShowSelected() {
	left.scatter.SetLabelX(t_("Time"));
	left.scatter.RemoveAllSeries();
	Vector<int> idsx, idsy, idsFixed;
	for (int rw = 0; rw < leftSearch.array.GetCount(); ++rw) {
		String param = Trim(leftSearch.array.Get(rw, 0));
		if (!param.IsEmpty()) {
			int col = datafast.FindCol(param);
			if (IsNull(col))
				statusBar.Temporary(Format("Parameter %s does not exist", param));
			else
				left.scatter.AddSeries(datafast.dataOut, 0, col, idsx, idsy, idsFixed, false).NoMark().Legend(param).Units(datafast.units[col], t_("sec")).Stroke(1);	
		}
	}
	for (int rw = 0; rw < rightSearch.array.GetCount(); ++rw) {
		String param = Trim(rightSearch.array.Get(rw, 0));
		if (!param.IsEmpty()) {
			int col = datafast.FindCol(param);
			if (IsNull(col))
				statusBar.Temporary(Format("Parameter %s does not exist", param));
			else
				left.scatter.AddSeries(datafast.dataOut, 0, col, idsx, idsy, idsFixed, false).NoMark().Legend(param).Units(datafast.units[col], t_("sec")).SetDataSecondaryY().Stroke(1);	
		}
	}
	bool rightEmpty = rightSearch.array.GetCount() == 0;
	left.scatter.SetDrawY2Reticle(!rightEmpty).SetDrawY2ReticleNumbers(!rightEmpty);
	left.scatter.SetSequentialXAll().SetFastViewX();
	left.scatter.ZoomToFit(true, true);	
}

/*Value FastScatter::DataSource::Format(const Value& q) const {
	ASSERT(datafast);
	return datafast->dataOut[col][q];
}*/

void FastScatterTabs::Init() {
	tabBar.Crosses(true).Crosses(true).ContextMenu(true);
	tabBar.SetTop();
	tabBar.WhenAction = THISBACK(OnTab);
	tabBar.WhenClose = THISBACK(OnCloseTab);

	AddFrame(tabBar);

	AddTab(t_("Empty"));
}

void FastScatterTabs::AddTab(String filename) {
	int id = tabBar.GetCursor();
	FastScatter *sct = nullptr;
	if (id >= 0 && (t_("Empty") == tabBar.GetValue(id))) { 
		int key = tabBar.GetKey(id);
		int idKey = tabKeys.Find(key);
		sct = &tabScatters[idKey];
		tabBar.SetValue(key, GetFileTitle(filename));
	} else {
		int key = tabCounter++;
		tabKeys << key;
		sct = &tabScatters.Add();
		sct->Init();
		tabBar.AddKey(key, GetFileTitle(filename));
		Add(sct->SizePos());
	}
	if (filename != t_("Empty")) {
		sct->file <<= filename;
		sct->file.WhenChange();
	}
}

void FastScatterTabs::OnTab() {	
	if (tabBar.GetCount() == 0)
		return; 
	
	int id = tabBar.GetCursor();
	if (id < 0)
		return;

	int key = tabBar.GetKey(id);
	int idKey = tabKeys.Find(key);
	for (int i = 0; i < tabKeys.GetCount(); ++i) 
		tabScatters[i].Show(i == idKey); 
}

void FastScatterTabs::OnCloseTab(Value key) {
	for (int i = 0; i < tabKeys.GetCount(); ++i) { 
		if (tabKeys[i] == key) {
			tabKeys.Remove(i);	
			tabScatters.Remove(i);
		}
	}
}

void FastScatterTabs::DragAndDrop(Point, PasteClip& d) {
	if (IsDragAndDropSource())
		return;
	if (AcceptFiles(d)) 
		LoadDragDrop(GetFiles(d));
}

bool FastScatterTabs::Key(dword key, int) {
	if (key == K_CTRL_V) 
		return LoadDragDrop(GetFiles(Ctrl::Clipboard()));
	return false;
}

bool FastScatterTabs::LoadDragDrop(const Vector<String> &files) {
	if (files.GetCount() == 0)
		return false;
		
	for (int i = 0; i < files.GetCount(); ++i) 
		AddTab(files[i]);

	return true;
}
