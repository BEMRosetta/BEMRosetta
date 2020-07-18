#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <RasterPlayer/RasterPlayer.h>
#include <TabBar/TabBar.h>

using namespace Upp;

#define IMAGECLASS ImgF
#define IMAGEFILE <BEMRosetta/BEMRosetta/FastScatter.iml>
#include <Draw/iml.h>

#include <BEMRosetta/BEMRosetta_cl/FastOut.h>
#include "auxiliar.h"
#include "FastScatter.h"

#include "clip.brc"

void FastScatter::Init(Function <void(String)> OnFile, StatusBar &_statusBar) {
	WhenFile = OnFile;
	statusBar = &_statusBar;
	CtrlLayout(*this);
	
	CtrlLayout(left);
	CtrlLayout(right);
	splitter.Horz(left.SizePos(), right.SizePos());
	splitter.SetPos(7000, 0);

	CtrlLayout(leftSearch);
	CtrlLayout(rightSearch);
	right.splitterSearch.Horz(leftSearch.SizePos(), rightSearch.SizePos());
	
	file.WhenChange = THISBACK(OnLoad);
	file.BrowseRightWidth(40).UseOpenFolder().BrowseOpenFolderWidth(10)
		.Tip(t_("Enter file path to show, or drop it from file explorer"));
	butLoad.Tip(t_("Loads FAST out/outb file")).WhenAction = [&] {file.DoGo();};
	file.Type(t_("FAST output file"), "*.out, *.outb"); 
		
	left.scatter.ShowAllMenus().SetMode(ScatterDraw::MD_DRAW);
	
	leftSearch.array.AutoHideSb().NoHeader().SetLineCy(EditField::GetStdHeight()).MultiSelect();
	rightSearch.array.AutoHideSb().NoHeader().SetLineCy(EditField::GetStdHeight()).MultiSelect();
	
	leftSearch.array.Tip(t_("Hover mouse to add new parameters to left axis"));
	rightSearch.array.Tip(t_("Hover mouse to add new parameters to right axis"));
	
	leftSearch.array.WhenDropInsert = [&] (int line, PasteClip& d) {OnDropInsert(line, d, leftSearch.array); };
	leftSearch.array.WhenDrop 		= [&] (PasteClip& d) {OnDrop(d, leftSearch.array);};
	leftSearch.array.WhenDrag 		= [&] {OnDrag(leftSearch.array, true);};
	
	rightSearch.array.WhenDropInsert= [&] (int line, PasteClip& d) {OnDropInsert(line, d, rightSearch.array); };
	rightSearch.array.WhenDrop 		= [&] (PasteClip& d) {OnDrop(d, rightSearch.array);};
	rightSearch.array.WhenDrag 		= [&] {OnDrag(rightSearch.array, true);};
	
	leftSearch.array.WhenLeftDouble = THISBACK1(WhenArrayLeftDouble, &leftSearch.array);
	leftSearch.array.WhenEnterKey   = THISBACK1(WhenArrayLeftDouble, &leftSearch.array);
	rightSearch.array.WhenLeftDouble = THISBACK1(WhenArrayLeftDouble, &rightSearch.array);	
	rightSearch.array.WhenEnterKey   = THISBACK1(WhenArrayLeftDouble, &rightSearch.array);	
	
	leftSearch.array.Removing().NoAskRemove().WhenArrayAction = THISBACK(ShowSelected);
	rightSearch.array.Removing().NoAskRemove().WhenArrayAction = THISBACK(ShowSelected);
	
	leftSearch.array.AddColumn("");
	rightSearch.array.AddColumn("");
	rightSearch.label.SetLabel(t_("Right axis"));
	
	frameSet.Add(leftSearch.array.GetRectEnter(), 1);
	frameSet.Add(rightSearch.array.GetRectEnter(), 1);
	//frameSet.WhenEnter = THISBACK(WhenFilter);
	frameSet.Set(leftSearch.array.GetRectEnter());
	
	right.filterParam.WhenAction = THISBACK1(OnFilter, true);
	right.filterParam.Tip(t_("Filter parameters to display. * are allowed"));
	right.filterUnits.WhenAction = THISBACK1(OnFilter, true);
	right.filterUnits.Tip(t_("Filter parameters to be displayed according to their units. * are allowed"));
	right.arrayParam.AddColumn(t_("Filtered fields"), 80);
	right.arrayParam.AddColumn(t_("Units"), 20);
	right.arrayParam.SetLineCy(EditField::GetStdHeight()).MultiSelect()
		 .Tip(t_("Double click to choose parameters to display"));
	right.arrayParam.WhenLeftDouble = THISBACK1(WhenArrayLeftDouble, &right.arrayParam);
	right.arrayParam.WhenEnterKey   = THISBACK1(WhenArrayLeftDouble, &right.arrayParam);
	
	right.arrayParam.WhenDropInsert= [&] (int line, PasteClip& d) {OnDropInsert(line, d, right.arrayParam); };
	right.arrayParam.WhenDrop 		= [&] (PasteClip& d) {OnDrop(d, right.arrayParam);};
	right.arrayParam.WhenDrag 		= [&] {OnDrag(right.arrayParam, false);};
	
	player.LoadBuffer(String(animatedStar, animatedStar_length));
	player.SetSpeed(0.5).Tip(t_("Click to update periodically plot from file"));
	player.WhenAction = [&] {
		if (!player.IsRunning()) {
			timeStop.Reset();
		 	player.Play();
			updateTime.WhenEnter();		
			OnTimer();
		} else 
			player.Stop();
	};
	timer.Set(int(-10*1000), THISBACK(OnTimer));
	
	updateTime.Tip(t_("Interval to update data from file"));
	updateTime <<= "1:00";
	updateTime.WhenEnter = [&] {
		int seconds = int(StringToSeconds(~updateTime));
		updateTime <<= SecondsToString(seconds, 0, false, false, true, false, true);
		updateTime.CancelSelection();
	};
}

void FastScatter::OnTimer() {
	if (!player.IsRunning())
		return;
	int time = int(StringToSeconds(~updateTime));
	
	if (timeStop.Seconds() < time)
		return;
	timeStop.Reset();
	
	OnLoad();
}

bool ArrayExists(const ArrayCtrl &a, String val) {
	for (int rw = 0; rw < a.GetCount(); ++rw)
		if (AsString(a.Get(rw, 0)) == val)
			return true;
	return false;
}

void ArrayRemoveDuplicates(const ArrayCtrl &a, int id, Vector<bool> &index) {
	if (index[id])
		return;
	for (int rw = id+1; rw < a.GetCount(); ++rw)
		if (a.Get(rw, 0) == a.Get(id, 0))
			index[rw] = true;
}

void ArrayRemoveDuplicates(ArrayCtrl &a) {
	Vector<bool> index;
	index.SetCount(a.GetCount(), false);
	for (int rw = 0; rw < a.GetCount(); ++rw)
		ArrayRemoveDuplicates(a, rw, index);
	for (int rw = a.GetCount()-1; rw >= 0; --rw)
		if (index[rw])
			a.Remove(rw);
}

void FastScatter::OnDropInsert(int line, PasteClip& d, ArrayCtrl &array) {
	if(AcceptInternal<ArrayCtrl>(d, "array")) {
		array.InsertDrop(line, d);
		array.SetFocus();
		ArrayRemoveDuplicates(array);
		ShowSelected();
	} else if(AcceptText(d) && !ArrayExists(array, GetString(d))) {
		array.Insert(line);
		array.Set(line, 0, GetString(d));
		array.SetCursor(line);
		array.SetFocus();
		ShowSelected();
	}	
}

void FastScatter::OnDrop(PasteClip& d, ArrayCtrl &array) {
	if(AcceptText(d)) {
		array.Add(GetString(d), GetString(d));
		array.SetFocus();
		ShowSelected();
	}
}

void FastScatter::OnDrag(ArrayCtrl &array, bool remove) {
	if(array.DoDragAndDrop(InternalClip(array, "array"), array.GetDragSample()) == DND_MOVE)
		if (remove)
			array.RemoveSelection();
}
	
bool FastScatter::AddParameter(String param,  ArrayCtrl *parray) {
	ArrayCtrl &arr = parray != NULL ? *parray : (leftSearch.array.IsShownFrame() ? leftSearch.array : rightSearch.array);
	for (int rw = 0; rw < arr.GetCount(); ++rw) {
		if (arr.Get(rw, 0) == param) 
			return false;
	}
	arr.Add(param);
	return true;
}

void FastScatter::WhenArrayLeftDouble(ArrayCtrl *parray) {
	ArrayCtrl &array = *parray;
	
	int id = array.GetCursor();
	if (id < 0)
		return;
	String param = array.Get(id, 0);
	if (parray == &leftSearch.array) {
		array.Remove(id);
		AddParameter(param, &rightSearch.array);	
	} else if (parray == &rightSearch.array) {
		array.Remove(id);
		AddParameter(param, &leftSearch.array);	
	} else
		AddParameter(param, nullptr);
	
	ShowSelected();
}
 
void FastScatter::OnFilter(bool show) {
	SortedVectorMap<String, String> list = datafast.GetList(~right.filterParam, ~right.filterUnits);
	if (list.GetCount() == 1) {
		if (AddParameter(list.GetKey(0), nullptr)) {
			right.filterParam.Clear();
			right.filterUnits.Clear();
			list = datafast.GetList();
			if (show)
				ShowSelected();
		}
	} 
	right.arrayParam.Clear();
	right.filterUnits.ClearList();
	for (int rw = 0; rw < list.GetCount(); ++rw) {
		right.arrayParam.Add(list.GetKey(rw), list[rw]);
		right.filterUnits.FindAddList(list[rw]);
	}
}

void FastScatter::Clear() {
	file.Clear();
	left.scatter.RemoveAllSeries();
	datafast.Clear();
	dataSource.Clear();
	OnFilter(true);
}

bool FastScatter::OnLoad() {
	try {
		WaitCursor waitcursor;
		
		if (!datafast.Load(~file)) {
			statusBar->Temporary(Format("File '%s' not found", ~file));
			return false;
		}
		
		statusBar->Temporary(t_("Loading file"));
		
		left.scatter.RemoveAllSeries();

		dataSource.SetCount(datafast.parameters.GetCount());
		for (int c = 0; c < dataSource.GetCount(); ++c) 
			dataSource[c].Init(datafast, c);
		
		LoadParams();
		OnFilter(false);
		ShowSelected();
		
		WhenFile(~file);
	} catch (Exc e) {
		Exclamation(Format("Error: %s", DeQtf(e)));	
		return false;
	}
	return true;
}

void FastScatter::ShowSelected() {
	left.scatter.SetLabelX(t_("Time"));
	left.scatter.RemoveAllSeries();
	Upp::Vector<int> idsx, idsy, idsFixed;
	for (int rw = 0; rw < leftSearch.array.GetCount(); ++rw) {
		String param = Trim(leftSearch.array.Get(rw, 0));
		if (!param.IsEmpty()) {
			int col = datafast.FindCol(param);
			if (IsNull(col))
				statusBar->Temporary(Format("Parameter %s does not exist", param));
			else
				left.scatter.AddSeries(datafast.dataOut, 0, col, idsx, idsy, idsFixed, false).NoMark().Legend(param).Units(datafast.units[col], t_("sec")).Stroke(1);	
		}
	}
	for (int rw = 0; rw < rightSearch.array.GetCount(); ++rw) {
		String param = Trim(rightSearch.array.Get(rw, 0));
		if (!param.IsEmpty()) {
			int col = datafast.FindCol(param);
			if (IsNull(col))
				statusBar->Temporary(Format("Parameter %s does not exist", param));
			else
				left.scatter.AddSeries(datafast.dataOut, 0, col, idsx, idsy, idsFixed, false).NoMark().Legend(param).Units(datafast.units[col], t_("sec")).SetDataSecondaryY().Stroke(1);	
		}
	}
	left.scatter.SetPlotAreaLeftMargin(70);
	bool rightEmpty = rightSearch.array.GetCount() == 0;
	if (!rightEmpty)
		left.scatter.SetPlotAreaRightMargin(70);
	left.scatter.SetDrawY2Reticle(!rightEmpty).SetDrawY2ReticleNumbers(!rightEmpty);
	left.scatter.SetSequentialXAll().SetFastViewX();
	left.scatter.ZoomToFit(true, true);	
}

Value FastScatter::DataSource::Format(const Value& q) const {
	ASSERT(datafast);
	return datafast->dataOut[col][q];
}

// GetAppDataFolder(), "BEMRosetta"
void FastScatterTabs::Init(String appDataFolder, StatusBar &_statusBar) {
	statusBar = &_statusBar;
	String folder = AppendFileNameX(appDataFolder, "FASTScatter");
	DirectoryCreate(folder);
	if (!DirectoryExists(folder))
		return;
		
	LoadFromJsonFile(*this, AppendFileName(folder, "FASTScatter.json"));
	
	tabBar.Crosses(true).ContextMenu(true);
	tabBar.SetTop();
	tabBar.Tip(t_("Click '+' to add new file in additional tab"));
	tabBar.WhenAction = THISBACK(OnTab);
	tabBar.WhenClose = THISBACK(OnCloseTab);
	tabBar.CancelClose = [=] (Value key) -> bool {return key == -1 || tabBar.Get(key) == t_("New Tab");};
	tabBar.CancelDragAndDrop = [=](int from, int to) {
		if (to < 0 || to > tabBar.GetCount()-1)
			return true;
		if (tabBar.CancelClose(tabBar.GetKey(from)))
			return true;
		if (tabBar.GetKey(to) == -1 || tabBar.GetValue(to) == t_("New Tab"))
			return true;
		return false;
	};
	
	styleTab = TabBar::StyleDefault();
	styleTab.crosses[0] = ImgF::VARIANT3_CR1();
	styleTab.crosses[1] = ImgF::VARIANT3_CR2();
	styleTab.crosses[2] = ImgF::VARIANT3_CR3();
	tabBar.SetStyle(styleTab);
	
	AddFrame(tabBar);

	tabBar.AddKey(-1, "", ImgF::Plus());
}

void FastScatterTabs::AddTab(String filename) {
	FastScatter *sct = nullptr;
	String title;
	if (filename == t_("New Tab"))
		title = t_("New Tab");
	else
		title = GetFileTitle(GetUpperFolder(filename)) + "/" + GetFileTitle(filename);
	
	int id = tabBar.GetCursor();
	Value key;
	if (id >= 0 && (t_("New Tab") == tabBar.GetValue(id))) { 
		key = tabBar.GetKey(id);
		int idKey = tabKeys.Find(key);
		fileNames[idKey] = filename;
		sct = &tabScatters[idKey];
		tabBar.SetValue(key, title);
	} else {
		key = tabCounter++;
		tabKeys << key;
		sct = &tabScatters.Add();
		int pos = max(tabBar.GetCount()-1, 0);
		sct->Init([=] (String filename) {
			String title = GetFileTitle(GetUpperFolder(filename)) + "/" + GetFileTitle(filename);
			tabBar.SetValue(key, title);
		}, *statusBar);
		tabBar.InsertKey(pos, key, title);
		Add(sct->SizePos());
		fileNames << filename;
	}
	AddHistory(filename);
	if (int(key) > -1 && filename != t_("New Tab")) {
		sct->file <<= filename;
		sct->file.WhenChange();
	}
}

void FastScatterTabs::AddHistory(String filename) {
	if (filename != t_("New Tab"))
		history.FindAdd(filename);
	for (int it = 0; it < tabBar.GetCount(); ++it) {
		const Value &key = tabBar.GetKey(it);
		if (int(key) > -1) {
			tabScatters[it].file.ClearHistory();
			for (int i = 0; i < history.GetCount(); ++i)
				tabScatters[it].file.AddHistory(history[i]);
		}
	}
}

void FastScatterTabs::OnTab() {	
	if (tabBar.GetCount() == 0)
		return; 
	int id = tabBar.GetCursor();
	if (id < 0)
		return;
	
	Value key = tabBar.GetKey(id);
	
	if (int(key) < 0) {		// Pressed '+'
		if (id > 0 && tabBar.FindValue(t_("New Tab")) >= 0)
			tabBar.SetCursor(tabBar.GetCursor()-1);
		else
			AddTab(t_("New Tab"));
		return;
	}
	
	int idKey = tabKeys.Find(key);
	for (int i = 0; i < tabKeys.GetCount(); ++i) 
		tabScatters[i].Show(i == idKey); 
}

void FastScatterTabs::OnCloseTab(Value key) {
	int idKey = tabKeys.Find(key);
	FastScatter &sct = tabScatters[idKey];
	sct.SaveParams();
	tabKeys.Remove(idKey);	
	tabScatters.Remove(idKey);
	fileNames.Remove(idKey);
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

bool FastScatterTabs::LoadDragDrop(const Upp::Vector<String> &files) {
	if (files.GetCount() == 0)
		return false;
		
	for (int i = 0; i < files.GetCount(); ++i) 
		AddTab(files[i]);

	return true;
}

FastScatterTabs::~FastScatterTabs() {
	for (int i = 0; i < tabScatters.GetCount(); ++i)
		tabScatters[i].SaveParams();	
	
	String file = AppendFileNameX(GetAppDataFolder(), "BEMRosetta", "FASTScatter", "FASTScatter.json");	
	StoreAsJsonFile(*this, file, true);
}

void FastScatter::Params::Get(const ArrayCtrl &aleft, const ArrayCtrl &aright) {
	left.Clear();
	for (int rw = 0; rw < aleft.GetCount(); ++rw)
		left << aleft.Get(rw, 0);
	right.Clear();
	for (int rw = 0; rw < aright.GetCount(); ++rw)
		right << aright.Get(rw, 0);
}

void FastScatter::Params::Set(ArrayCtrl &aleft, ArrayCtrl &aright) const {
	aleft.Clear();
	for (int rw = 0; rw < left.GetCount(); ++rw)
		aleft.Set(rw, 0, left[rw]);
	aright.Clear();
	for (int rw = 0; rw < right.GetCount(); ++rw)
		aright.Set(rw, 0, right[rw]);
}
	
void FastScatter::LoadParams() {
	String strpath = SHA1StringS(~file).Left(12);
	String strname = SHA1StringS(GetFileName(~file)).Left(12);
	
	String folder = AppendFileNameX(GetAppDataFolder(), "BEMRosetta", "FASTScatter");
	DirectoryCreate(folder);
	if (!DirectoryExists(folder))
		return;
	
	String fileName;
	FindFile ffpath(AppendFileName(folder, strpath + "_*.json"));
	if (ffpath) 
		fileName = ffpath.GetPath();
	else {
		FindFile ffname(AppendFileName(folder,  "*_" + strname + ".json"));
		if (ffname) 
			fileName = ffname.GetPath();
	}
	if (fileName.IsEmpty())
		return;

	Params params;
	if (!LoadFromJsonFile(params, fileName))
		return;
	
	params.Set(leftSearch.array, rightSearch.array);
}

void FastScatter::SaveParams() {
	String strpath = SHA1StringS(~file).Left(12);
	String strname = SHA1StringS(GetFileName(~file)).Left(12);
	
	String folder = AppendFileNameX(GetAppDataFolder(), "BEMRosetta", "FASTScatter");
	DirectoryCreate(folder);
	if (!DirectoryExists(folder))
		return;
	
	String fileName = AppendFileName(folder, strpath + "_" + strname + ".json");

	Params params;
	params.Get(leftSearch.array, rightSearch.array);
	StoreAsJsonFile(params, fileName, true);
}

void MainFASTW::Init(String appDataFolder, const Image &icon, const Image &largeIcon, StatusBar &statusBar) {
	fast.Init(appDataFolder, statusBar);
	Add(fast.SizePos());
	Title(t_("BEMRosetta FAST .out Reader")).Sizeable().Zoomable().Icon(icon, largeIcon);
}