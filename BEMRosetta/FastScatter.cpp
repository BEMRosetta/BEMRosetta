// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <RasterPlayer/RasterPlayer.h>
#include <TabBar/TabBar.h>
#include <DropGrid/DropGrid.h>

using namespace Upp;

#define IMAGECLASS ImgF
#define IMAGEFILE <BEMRosetta/FastScatter.iml>
#include <Draw/iml.h>

#include <BEMRosetta_cl/BEMRosetta.h>
#include <BEMRosetta_cl/FastOut.h>
#include "auxiliar.h"
#include "FastScatter.h"

#include "clip.brc"


void FastScatter::Init(Function <void(String)> OnFile, Function <void(String)> OnCopyTabs, StatusBar &statusBar) {
	int len = StdFont().GetHeight();
	
	CtrlLayout(compare);
	compare.butCalc <<= THISBACK(OnCalc);
	compare.editStart <<= 0;
	compare.array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, compare.array);};
	
	fscatter.Init(this, OnFile, OnCopyTabs, statusBar);
#ifdef flagDEBUG
	splitCompare.Horz(fscatter.SizePos(), compare.SizePos());
	splitCompare.SetPositions(3000, 10000).SetInitialPositionId(1).SetButtonNumber(1).SetButtonWidth(len);
	Add(splitCompare.SizePos());
#else
	Add(fscatter.SizePos());
#endif
}
	
void FastScatter::OnCalc() {
	try {
		compare.array.Reset();
		compare.array.SetLineCy(EditField::GetStdHeight()).MultiSelect();
		
		int idTime = fscatter.datafast.GetIdTime(~compare.editStart);
		int num = fscatter.datafast.size();
		if (idTime < 0 || idTime >= num) {
			Exclamation("Bad start time");
			return;
		}
		const UVector<UVector<String>> params = {{"+FAIRTEN1_t", ".3E", "max"},
												 {"+FAIRTEN2_t", ".3E", "max"},
												 {"+FAIRTEN3_t", ".3E", "max"},
												 {"+FAIRTEN4_t", ".3E", "max"},
												 {"+FAIRTEN5_t", ".3E", "max"},
												 {"+FAIRTEN6_t", ".3E", "max"},
												 {"+FAIRTEN7_t", ".3E", "max"},
												 {"+PtfmShift",  ".0f", "mean", "maxval"}, 
												 {"+PtfmTilt",   ".1f", "mean", "maxval"}, 
												 {"PtfmYaw",     ".1f", "mean", "maxval"}, 
												 {"+NcIMUTA",    ".1f", "max"}, 
												 {"+TwrBsShear", ".0f", "max"},
												 {"+TwrBsBend",  ".0f", "max"}
		};
		
		int col = 0;
		for (const UVector<String> &param : params) {				
			if (param.size() < 3)
				throw Exc(t_("Wrong number of parameters"));
			String strpar = param[0];
			int id = fscatter.datafast.GetParameter_throw(strpar);
			String units = fscatter.datafast.GetUnit(id);
			for (int i = 2; i < param.size(); i++) {
				if (i == 2)
					compare.array.AddColumn(Format("%s [%s]", strpar, units));
				else
					compare.array.AddColumn("");
	
				String stat = param[i];
				compare.array.Set(0, col++, stat);
			}
		}
		
		int row = 1;
		col = 0;			
		for (const UVector<String> &param : params) {
			String strpar = param[0];
			String format = param[1];
			int id = fscatter.datafast.GetParameter_throw(strpar);
			const VectorXd data = fscatter.datafast.GetVector(id).tail(num - idTime);
			
			for (int i = 2; i < param.size(); i++) {
				String stat = param[i];
				double val;
				if (stat == "mean") 
					val = data.mean();
				else if (stat == "min") 
					val = data.minCoeff();
				else if (stat == "max") 
					val = data.maxCoeff();
				else if (stat == "maxval") { 
					double mx = data.maxCoeff();
					double mn = data.minCoeff();
					if (abs(mx) > abs(mn))
						val = mx;
					else
						val = mn;
				} else if (stat == "std") 
					val = sqrt((data.array() - data.mean()).square().sum() / (data.size() - 1));
				else
					throw Exc(Format("Unknown '%s' statistic in parameter '%s'", stat, strpar));
				
				compare.array.Set(row, col++, Format("%" + format, val));
			}
		}
	} catch (const Exc &e) {
		Exclamation(Format("Error: %s", DeQtf(e)));	
	}
}
	
void FastScatterBase::Init(FastScatter *parent, Function <void(String)> OnFile, Function <void(String)> OnCopyTabs, StatusBar &_statusBar) {
	fastScatter = parent;
	WhenFile = OnFile;
	WhenCopyTabs = OnCopyTabs;
	statusBar = &_statusBar;
	CtrlLayout(*this);
	
	CtrlLayout(left);
	CtrlLayout(right);
	splitter.Horz(left.SizePos(), right.SizePos());
	splitter.SetPos(7000, 0);

	CtrlLayout(leftSearch);
	CtrlLayout(rightSearch);
	right.splitterSearch.Horz(leftSearch.SizePos(), rightSearch.SizePos());
	
	file.WhenChange = THISBACK1(OnLoad, false);
	file.BrowseRightWidth(40).UseOpenFolder().BrowseOpenFolderWidth(10)
		.Tip(t_("Enter file path to show, or drop it from file explorer"));
	butLoad.Tip(t_("Loads FAST out/outb file")) << [&] {file.DoGo();};
	file.Type(t_("FAST output file"), "*.out, *.outb"); 
	butSaveAs <<= THISBACK(OnSaveAs);
	butSaveAs.Tip(t_("Saves data file"));
	dropFormat.Add(".out").Add(".csv").Add(".csv only selected");
	dropFormat.SetIndex(0);
	
	UpdateButtons(false);
	
	left.scatter.ShowAllMenus().SetMode(ScatterDraw::MD_DRAW);
	
	leftSearch.array.AutoHideSb().NoHeader().SetLineCy(EditField::GetStdHeight()).MultiSelect();
	rightSearch.array.AutoHideSb().NoHeader().SetLineCy(EditField::GetStdHeight()).MultiSelect();
	
	leftSearch.array.Tip(t_("Hover mouse to add new parameters to left axis"));
	rightSearch.array.Tip(t_("Hover mouse to add new parameters to right axis"));
	
	leftSearch.array.WhenDropInsert = [&] (int line, PasteClip& d) {OnDropInsert(line, d, leftSearch.array); };
	leftSearch.array.WhenDrag 		= [&] {OnDrag(leftSearch.array, true);};
	
	rightSearch.array.WhenDropInsert= [&] (int line, PasteClip& d) {OnDropInsert(line, d, rightSearch.array); };
	rightSearch.array.WhenDrag 		= [&] {OnDrag(rightSearch.array, true);};
	
	leftSearch.array.WhenLeftDouble  = THISBACK1(WhenArrayLeftDouble, &leftSearch.array);
	leftSearch.array.WhenEnterKey    = THISBACK1(WhenArrayLeftDouble, &leftSearch.array);
	rightSearch.array.WhenLeftDouble = THISBACK1(WhenArrayLeftDouble, &rightSearch.array);	
	rightSearch.array.WhenEnterKey   = THISBACK1(WhenArrayLeftDouble, &rightSearch.array);	
	
	leftSearch.array.Removing().NoAskRemove().WhenArrayAction = THISBACK(ShowSelected);
	rightSearch.array.Removing().NoAskRemove().WhenArrayAction = THISBACK(ShowSelected);
	
	leftSearch.array.AddColumn("");
	leftSearch.array.NoVertGrid();
	rightSearch.array.AddColumn("");
	rightSearch.array.NoVertGrid();
	rightSearch.label.SetLabel(t_("Right axis"));
	
	frameSet.Add(leftSearch.array.GetRectEnter(), 1);
	frameSet.Add(rightSearch.array.GetRectEnter(), 1);
	frameSet.Set(leftSearch.array.GetRectEnter());
	
	right.filterParam.WhenAction = THISBACK1(OnFilter, true);
	right.filterParam.Tip(t_("Filters parameters to display. * are allowed"));
	right.filterUnits.WhenAction = THISBACK1(OnFilter, true);
	right.filterUnits.Tip(t_("Filters parameters to be displayed according to their units. * are allowed"));
	right.arrayParam.AddColumn(t_("Fields"), 80);
	right.arrayParam.AddColumn(t_("Units"), 20);
	right.arrayParam.SetLineCy(EditField::GetStdHeight()).MultiSelect()
		 .Tip(t_("Double click to choose parameters to display"));
	right.arrayParam.WhenLeftDouble = THISBACK1(WhenArrayLeftDouble, &right.arrayParam);
	right.arrayParam.WhenEnterKey   = THISBACK1(WhenArrayLeftDouble, &right.arrayParam);
	
	right.butCopy 		<<= THISBACK(SelCopy);
	right.butCopy.Tip(t_("Copy selected parameters to clipboard"));
	right.butPaste 		<<= THISBACK1(SelPaste, "");
	right.butPaste.Tip(t_("Paste selected parameters from clipboard"));
	right.butSetOnAllTabs 	<<= THISBACK(SelCopyTabs);
	right.butSetOnAllTabs.Tip(t_("Copy selected parameters on the other tabs"));
	
	right.arrayParam.WhenDropInsert= [&] (int line, PasteClip& d) {OnDropInsert(line, d, right.arrayParam); };
	right.arrayParam.WhenDrag 		= [&] {OnDrag(right.arrayParam, false);};
	
	player.LoadBuffer(String(animatedStar, animatedStar_length));
	player.SetSpeed(0.5).Tip(t_("Click to update periodically plot from file"));
	player << [&] {
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
	
	saveFolder = GetDesktopFolder();
}

String FastScatterBase::SelectedStr() {
	String strLeft;
	for (int rw = 0; rw < leftSearch.array.GetCount(); ++rw) {
		String param = Trim(leftSearch.array.Get(rw, 0));
		if (!param.IsEmpty()) {
			if (!strLeft.IsEmpty())
				strLeft << ",";
			strLeft << param;
		}
	}
	String strRight;
	for (int rw = 0; rw < rightSearch.array.GetCount(); ++rw) {
		String param = Trim(rightSearch.array.Get(rw, 0));
		if (!param.IsEmpty()) {
			if (!strRight.IsEmpty())
				strRight << ",";
			strRight << param;
		}
	}
	return strLeft + ";" + strRight;
}

void FastScatterBase::SelCopy() {
	WriteClipboardText(SelectedStr());
}

void FastScatterBase::SelPaste(String str) {
	if (str.IsEmpty())
		str = ReadClipboardText();
	UVector<String> params = Split(str, ";");
	params.SetCount(2);
	UVector<String> left = Split(params[0], ",");
	leftSearch.array.Clear();
	for (int rw = 0; rw < left.size(); ++rw) {
		if (!IsNull(datafast.GetParameter(left[rw])))
			leftSearch.array.Set(rw, 0, left[rw]);
	}
	UVector<String> right = Split(params[1], ",");
	rightSearch.array.Clear();
	for (int rw = 0; rw < right.size(); ++rw) {
		if (!IsNull(datafast.GetParameter(right[rw])))
			rightSearch.array.Set(rw, 0, right[rw]);
	}
	ShowSelected();	
}

void FastScatterBase::SelCopyTabs() {
	WhenCopyTabs(SelectedStr());
}
	
void FastScatterBase::OnTimer() {
	if (!player.IsRunning())
		return;
	int time = int(StringToSeconds(~updateTime));
	
	if (timeStop.Seconds() < time)
		return;
	timeStop.Reset();
	
	OnLoad(true);
}

bool ArrayExists(const ArrayCtrl &a, String val) {
	for (int rw = 0; rw < a.GetCount(); ++rw)
		if (AsString(a.Get(rw, 0)) == val)
			return true;
	return false;
}

void ArrayRemoveDuplicates(const ArrayCtrl &a, int id, UVector<bool> &index) {
	if (index[id])
		return;
	for (int rw = id+1; rw < a.GetCount(); ++rw)
		if (a.Get(rw, 0) == a.Get(id, 0))
			index[rw] = true;
}

void ArrayRemoveDuplicates(ArrayCtrl &a) {
	UVector<bool> index;
	index.SetCount(a.GetCount(), false);
	for (int rw = 0; rw < a.GetCount(); ++rw)
		ArrayRemoveDuplicates(a, rw, index);
	for (int rw = a.GetCount()-1; rw >= 0; --rw)
		if (index[rw])
			a.Remove(rw);
}

void FastScatterBase::OnDropInsert(int line, PasteClip& d, ArrayCtrl &array) {
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

void FastScatterBase::OnDrag(ArrayCtrl &array, bool remove) {
	if(array.DoDragAndDrop(InternalClip(array, "array"), array.GetDragSample()) == DND_MOVE)
		if (remove)
			array.RemoveSelection();
}
	
bool FastScatterBase::AddParameter(String param,  ArrayCtrl *parray) {
	ArrayCtrl &arr = parray != NULL ? *parray : (leftSearch.array.IsShownFrame() ? leftSearch.array : rightSearch.array);
	for (int rw = 0; rw < arr.GetCount(); ++rw) {
		if (arr.Get(rw, 0) == param) 
			return false;
	}
	arr.Add(param);
	return true;
}

void FastScatterBase::WhenArrayLeftDouble(ArrayCtrl *parray) {
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
 
void FastScatterBase::OnFilter(bool show) {
	SortedVectorMap<String, String> list = datafast.GetList(Trim(~right.filterParam), Trim(~right.filterUnits));
	if (list.size() == 1) {
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
	for (int rw = 0; rw < list.size(); ++rw) {
		right.arrayParam.Add(list.GetKey(rw), list[rw]);
		right.filterUnits.FindAddList(list[rw]);
	}
}

void FastScatterBase::Clear() {
	file.Clear();
	left.scatter.RemoveAllSeries();
	datafast.Clear();
	dataSource.Clear();
	OnFilter(true);
	UpdateButtons(false);
}

void FastScatterBase::UpdateButtons(bool on) {
	butSaveAs.Enable(on);
	dropFormat.Enable(on);
}

bool FastScatterBase::OnLoad(bool justUpdate) {
	try {
		WaitCursor waitcursor;

		if (!justUpdate)
			statusBar->Temporary(t_("Loading file"));
		else
			statusBar->Temporary(t_("Updating data"));
		Ctrl::ProcessEvents();
		
		left.scatter.Disable();
		
		String fileName = FastOut::GetFileToLoad(~file);
		if (IsNull(fileName)) {
			statusBar->Temporary(Format(t_("File '%s' not found"), ~file));
			left.scatter.Enable();
			UpdateButtons(false);
			return false;
		}
		if (!datafast.Load(fileName)) {
			statusBar->Temporary(Format(t_("File '%s' temporarily blocked by OpenFAST"), ~file));
			left.scatter.Enable();
			UpdateButtons(false);
			return false;
		}
		
		dataSource.SetCount(datafast.parameters.size());
		if (!justUpdate) {
			for (int c = 0; c < dataSource.size(); ++c) 
				dataSource[c].Init(datafast, c);
		
			LoadParams();
			OnFilter(false);
			ShowSelected();
		
			WhenFile(fileName);
		} else 
			ShowSelected();
		
	} catch (const Exc &e) {
		Exclamation(Format("Error: %s", DeQtf(e)));	
		left.scatter.Enable();
		UpdateButtons(false);
		return false;
	}
	left.scatter.Enable();
	
	UpdateButtons(true);
	return true;
}

void FastScatterBase::OnSaveAs() {
	try {
		statusBar->Temporary(t_("Saving data"));
		String fileType = ~dropFormat;
		
		String ext;
		FileSel fs;
		if (fileType == ".out") {
			fs.Type(t_("OpenFAST .out format"), "*.out");
			ext = ".out";
		} else {
			fs.Type(t_("CSV format"), ".csv");
			ext = ".csv";
		}
		fs.ActiveDir(saveFolder);
		fs.ActiveType(0);
		fs.Set(ForceExt(~file, ext));
		
		if (fs.ExecuteSaveAs(t_("Save register data"))) {
			WaitCursor waitcursor;
			
			if (fileType == ".out") {	
				if (!datafast.Save(~fs, ".out")) {
					Exclamation(Format("Problem saving '%s'", ~fs));
					return;
				}
			} else if (fileType == ".csv") {	
				if (!datafast.Save(~fs, ".csv", ScatterDraw::GetDefaultCSVSeparator())) {
					Exclamation(Format("Problem saving '%s'", ~fs));
					return;
				}
			} else {
				if (!left.scatter.SaveToFileData(~fs)) {
					Exclamation(Format("Problem saving '%s'", ~fs));
					return;
				}
			}
		}
		saveFolder = GetFileFolder(~fs);
	} catch (const Exc &e) {
		Exclamation(Format("Error: %s", DeQtf(e)));	
	}		
}
	
void FastScatterBase::ShowSelected() {
	WaitCursor wait;
	
	left.scatter.SetLabelX(t_("Time"));
	left.scatter.RemoveAllSeries().SetSequentialXAll().SetFastViewX();
	UVector<int> idsx, idsy, idsFixed;
	for (int rw = 0; rw < leftSearch.array.GetCount(); ++rw) {
		String param = Trim(leftSearch.array.Get(rw, 0));
		if (!param.IsEmpty()) {
			int col = datafast.GetParameter(param);
			if (IsNull(col))
				statusBar->Temporary(Format("Parameter %s does not exist", param));
			else
				left.scatter.AddSeries(datafast.dataOut, 0, col, idsx, idsy, idsFixed, false).NoMark().Legend(param).Units(datafast.units[col], t_("s")).Stroke(1);	
		}
	}
	for (int rw = 0; rw < rightSearch.array.GetCount(); ++rw) {
		String param = Trim(rightSearch.array.Get(rw, 0));
		if (!param.IsEmpty()) {
			int col = datafast.GetParameter(param);
			if (IsNull(col))
				statusBar->Temporary(Format("Parameter %s does not exist", param));
			else
				left.scatter.AddSeries(datafast.dataOut, 0, col, idsx, idsy, idsFixed, false).NoMark().Legend(param).Units(datafast.units[col], t_("s")).SetDataSecondaryY().Stroke(1);	
		}
	}
	left.scatter.SetPlotAreaLeftMargin(8*StdFont().GetHeight());
	left.scatter.SetPlotAreaBottomMargin(4*StdFont().GetHeight());
	bool rightEmpty = rightSearch.array.GetCount() == 0;
	if (!rightEmpty)
		left.scatter.SetPlotAreaRightMargin(8*StdFont().GetHeight());
	else
		left.scatter.SetPlotAreaRightMargin(2*StdFont().GetHeight());
	left.scatter.SetDrawY2Reticle(!rightEmpty).SetDrawY2ReticleNumbers(!rightEmpty);
	left.scatter.ZoomToFit(true, true);	
	
	SaveParams();
}

Value FastScatterBase::DataSource::Format(const Value& q) const {
	ASSERT(datafast);
	return datafast->dataOut[col][q];
}

void FastScatterTabs::Init(String appDataFolder, StatusBar &_statusBar) {
	statusBar = &_statusBar;
	String folder = AppendFileNameX(appDataFolder, "FASTScatter");
	if (!DirectoryCreateX(folder))
		return;
		
	LoadFromJsonFile(*this, AppendFileNameX(folder, "FASTScatter.json"));
	
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
				String title = GetFileName(GetUpperFolder(filename)) + "/" + GetFileTitle(filename);
				tabBar.SetValue(key, title);
				AddHistory(filename);
			}, [=] (String clipboard) {
				if (tabBar.GetCount() == 0)
					return; 
				int id = tabBar.GetCursor();
				if (id < 0)
					return;
				Value key = tabBar.GetKey(id);
				for (int i = 0; i < tabKeys.size(); ++i) {
					if (tabKeys[i] != key) {
						tabScatters[i].fscatter.SelPaste(clipboard);
						tabScatters[i].compare.editStart <<= ~tabScatters[id].compare.editStart;
					}
				}
			}, *statusBar);
		tabBar.InsertKey(pos, key, title);
		Add(sct->SizePos());
		fileNames << filename;
	}
	AddHistory(filename);
	if (int(key) > -1 && filename != t_("New Tab")) {
		sct->fscatter.file <<= filename;
		if (!sct->fscatter.file.WhenChange())
			OnCloseTab(key); 
	}
}

void FastScatterTabs::AddHistory(String filename) {
	if (filename != t_("New Tab")) {
		int id = history.Find(filename);
		if (id >= 0)
			history.Remove(id);
		history.Insert(0, filename);
	}
	for (int it = 0; it < tabBar.GetCount(); ++it) {
		const Value &key = tabBar.GetKey(it);
		if (int(key) > -1 && it < tabScatters.size()) {
			tabScatters[it].fscatter.file.ClearHistory();
			for (int i = history.size()-1; i >= 0; --i)
				tabScatters[it].fscatter.file.AddHistory(history[i]);
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
	for (int i = 0; i < tabKeys.size(); ++i) 
		tabScatters[i].Show(i == idKey); 
}

void FastScatterTabs::OnCloseTab(Value key) {
	int idKey = tabKeys.Find(key);
	if (idKey < 0)
		return;
	FastScatterBase &sct = tabScatters[idKey].fscatter;
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

bool FastScatterTabs::LoadDragDrop(const UVector<String> &files) {
	if (files.size() == 0)
		return false;
	
	UVector<String> fis;
	for (int i = 0; i < files.size(); ++i)
		fis.Append(FastOut::GetFilesToLoad(files[i]));
	
	if (fis.size() > 6) {
		if (!PromptYesNo(Format(t_("You are about to open %d files"), fis.size())))
			return false;
	}
	
	for (int i = 0; i < fis.size(); ++i) 
		AddTab(fis[i]);

	return true;
}

FastScatterTabs::~FastScatterTabs() {
	for (int i = 0; i < tabScatters.size(); ++i)
		tabScatters[i].fscatter.SaveParams();	
	
	String file = AppendFileNameX(GetAppDataFolder(), "BEMRosetta", "FASTScatter", "FASTScatter.json");	
	StoreAsJsonFile(*this, file, true);
}

void FastScatterBase::Params::Get(const ArrayCtrl &aleft, const ArrayCtrl &aright) {
	left.Clear();
	for (int rw = 0; rw < aleft.GetCount(); ++rw)
		left << aleft.Get(rw, 0);
	right.Clear();
	for (int rw = 0; rw < aright.GetCount(); ++rw)
		right << aright.Get(rw, 0);
}

void FastScatterBase::Params::Set(ArrayCtrl &aleft, ArrayCtrl &aright) const {
	aleft.Clear();
	for (int rw = 0; rw < left.size(); ++rw)
		aleft.Set(rw, 0, left[rw]);
	aright.Clear();
	for (int rw = 0; rw < right.size(); ++rw)
		aright.Set(rw, 0, right[rw]);
}
	
void FastScatterBase::LoadParams() {
	String strpath = SHA1StringS(~file).Left(12);
	String strname = SHA1StringS(GetFileName(~file)).Left(12);
	
	String folder = AppendFileNameX(GetAppDataFolder(), "BEMRosetta", "FASTScatter");
	if (!DirectoryCreateX(folder))
		return;
	
	String fileName;
	FindFile ffpath(AppendFileNameX(folder, strpath + "_*.json"));
	if (ffpath) 
		fileName = ffpath.GetPath();
	else {
		FindFile ffname(AppendFileNameX(folder,  "*_" + strname + ".json"));
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

void FastScatterBase::SaveParams() {
	String strpath = SHA1StringS(~file).Left(12);
	String strname = SHA1StringS(GetFileName(~file)).Left(12);
	
	String folder = AppendFileNameX(GetAppDataFolder(), "BEMRosetta", "FASTScatter");
	if (!DirectoryCreateX(folder))
		return;
	
	String fileName = AppendFileNameX(folder, strpath + "_" + strname + ".json");

	Params params;
	params.Get(leftSearch.array, rightSearch.array);
	StoreAsJsonFile(params, fileName, true);
}

void MainFASTW::Init(String appDataFolder, const Image &icon, const Image &largeIcon, StatusBar &statusBar, Function <void()> _WhenClose) {
	WhenClose = _WhenClose;
	fast.Init(appDataFolder, statusBar);
	Add(fast.SizePos());
	Title(t_("BEMRosetta FAST .out+b Reader")).Sizeable().Zoomable().Icon(icon, largeIcon);
}
