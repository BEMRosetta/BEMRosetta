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


void FastScatter::Init(Function <bool(String)> OnFile, Function <void(String)> OnCopyTabs, StatusBar &statusBar) {
	int len = StdFont().GetHeight();
	
	CtrlLayout(compare);
	compare.butCalc <<= THISBACK(OnCalc);
	compare.array.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, compare.array);};
	
	compare.fileMetrics.WhenChange  = [&] {
		bool ret = OnLoadCompare();
		if (ret) {
			FastScatterTabs *pf = GetParentCtrlP<FastScatterTabs>(this);
			pf->AddHistoryMetrics(~compare.fileMetrics);
		}
		return ret;
	};
	compare.fileMetrics.activeDir = GetDesktopFolder();
	compare.fileMetrics.BrowseRightWidth(40).UseOpenFolder().BrowseOpenFolderWidth(10)
		.Tip(t_("Enter file with set of parameters and metrics"));
	compare.fileMetrics.Type(t_("json file"), "*.json").Type(t_("All files"), "*.*"); 
	compare.butLoad.Tip(t_("Loads parameters and metrics file")) << [&] {compare.fileMetrics.DoGo();};
	compare.butSave.Tip(t_("Saves parameters and metrics")) << [&] {OnSaveCompare();};
	
	compare.arrayStats.AddColumn(t_("Parameter"), 80);
	compare.arrayStats.GetColumn(0).Edit(edits.Add()); 
	compare.arrayStats.AddColumn(t_("# decimals"), 80);
	compare.arrayStats.GetColumn(1).Edit(edits.Add()); 
	compare.arrayStats.AddColumn(t_("Statistics"), 200);
	compare.arrayStats.GetColumn(2).Edit(edits.Add()); 
	compare.arrayStats.Proportional().Clipboard().Editing().Removing().Appending().Duplicating().SetToolBar().MultiSelect().ExtraPaste().SelectRow(false);

	compare.dropOperations.AddColumn("Metric", 10);
	compare.dropOperations.AddColumn("Comment", 20);
	compare.dropOperations.Add("mean", "")
						  .Add("min", "")
						  .Add("max", "")
						  .Add("maxval", "max(abs(max), abs(min))")
						  .Add("minmean", "min - mean")
						  .Add("maxmean", "max - mean")
						  .Add("stddev", "stddev")
						  .Add("amplitude", "RAO*height")
						  .Add("rao", "RAO from FFT")
						  .Add("rao_mean", "RAO from average")
						  .Add("percentile", "Percentile (ratio)")
						  .Add("weibull", "Weibull fitting percentile (ratio)")
						  .Add("td", "Design tension (time, γmean, γdyn)");
	compare.dropOperations.Resizeable(false).Width(200);					  
	compare.dropOperations.WhenAction = [&] {
		String operation;
		int r = compare.dropOperations.GetIndex();
		if (r < 0)
			return;
		operation = compare.dropOperations.Get(r, 0);
		compare.edOperations <<= operation;
		r = compare.arrayStats.GetCursorId();
		if (r < 0)
			return;
		String str = compare.arrayStats.Get(r, 2);
		if (!str.IsEmpty())
			str << ", ";
		str << operation << "()";
		compare.arrayStats.Set(r, 2, str);
	};
	compare.edOperations <<= "mean";
	
	fscbase.Init(this, OnFile, OnCopyTabs, statusBar);
//#ifdef flagDEBUG
	splitCompare.Horz(fscbase.SizePos(), compare.SizePos());
	splitCompare.SetPositions(3000, 10000).SetInitialPositionId(1).SetButtonNumber(1).SetButtonWidth(len);
	Add(splitCompare.SizePos());
//#else
//	Add(fscbase.SizePos());
//#endif
}

bool FastScatter::OnLoadCompare() {
	String filename = Trim(~compare.fileMetrics);
	
	if (filename.IsEmpty()) {
		compare.fileMetrics.DoBrowse();
		return true;
	}
	if (!LoadFromJsonFile(params, filename)) {
		BEM::PrintError(t_("Impossible to load file"));
		return false;
	}
	ParamsToGrid();
		
	return true;
}

void FastScatter::ParamsToGrid() {
	compare.arrayStats.Clear(false);
	for (int i = 0; i < params.params.size(); ++i) {
		String stats;
		for (int is = 0; is < params.params[i].metrics.size(); ++is) {
			if (is > 0)
				stats << ", ";
			stats << params.params[i].metrics[is];
		}
		compare.arrayStats.Add(params.params[i].name, params.params[i].decimals, stats);
	}
}

void FastScatter::GridToParams() {
	params.params.Clear();
	params.params.SetCount(compare.arrayStats.GetRowCount());
	for (int r = 0; r < compare.arrayStats.GetRowCount(); ++r) {
		params.params[r].name = AsString(compare.arrayStats.Get(r, 0));
		params.params[r].decimals = ScanInt(compare.arrayStats.Get(r, 1).ToString());
		String stats = compare.arrayStats.Get(r, 2).ToString();
		int inparen = 0;
		String met;
		for (char c : stats) {		// converts "rao_mean, rao(), factor(2, 5)" into ["rao_mean", "rao", "factor(2, 5)"]
			if (c == '(') {
				if (inparen != 0)
					throw Exc("Parenthesis opened after parenthesis is not allowed");
				inparen = 1;
				met << c;
			} else if (c == ')') {
				if (inparen != 1)
					throw Exc("Parenthesis closed after parenthesis is not allowed");
				inparen = 2;
				met << c;
			} else if (c == ',' && (inparen == 0 || inparen == 2)) {
				String fun = Trim(met);
				if (!fun.IsEmpty())
					params.params[r].metrics << fun;
				met.Clear();
				inparen = 0;
			} else if (inparen < 2)
				met << c;
			// Discards everything after ')' and before ','
		}	
		if (!met.IsEmpty()) {
			String fun = Trim(met);
			if (!fun.IsEmpty())
				params.params[r].metrics << fun;
		}
	}
}

void FastScatter::OnSaveCompare() {
	GridToParams();
	
	GuiLock __;
	
	try {
		FileSel fs;
		
		fs.Type("json file", "*.json");
		
		fs.ActiveType(0);
		fs.ActiveDir(GetFileFolder(~compare.fileMetrics));
		
		if (!fs.ExecuteSaveAs(t_("Save parameters")))
			throw Exc(t_("Cancelled by the user"));
		
		String fileName = ~fs;
				
		if (!StoreAsJsonFile(params, fileName, true))
			throw Exc(t_("Impossible to save file"));
		
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void FastScatter::OnCalc() {
	auto GetEquationName = [&](String eq) {
		int id = eq.Find(":");
		if (id >= 0)
			return Trim(eq.Left(id));
		return eq;
	};
	
	try {
		compare.array.Reset();
		compare.array.SetLineCy(EditField::GetStdHeight()).MultiSelect().NoHeader();
		
		if (fscbase.IsEmpty()) {
			BEM::PrintError("No data loaded");
			return;	
		}

		GridToParams();
				
		UVector<UVector<Value>> table;
		ParameterMetrics realparams; 
		
		Progress progress(t_("Calculating"), 100);
		auto Status = [&](String str, int _pos) {
			if (!IsNull(str))
				progress.SetText(str); 
			progress.SetPos(_pos); 
			return !progress.Canceled();
		};
			
		Calc(fscbase.left.dataFast, params, realparams, fscbase.timeStart, fscbase.timeEnd, table, Status);
		
		int col = 0;
		compare.array.AddColumn("");
		compare.array.Set(0, col++, t_("File"));
		compare.array.AddColumn("");
		compare.array.Set(0, col++, t_("Begin [s]"));
		compare.array.AddColumn("");
		compare.array.Set(0, col++, t_("End [s]"));
		FastOut &fast = fscbase.left.dataFast[0];
		UVector<String> fmt;
		fmt << "s" << "s" << "s";
				
		for (const ParameterMetric &param : realparams.params) {				
			String units;
			int id;
			if ((id = fast.GetParameterX(param.name)) >= 0)
				units = fast.GetUnit(id);
			for (int i = 0; i < param.metrics.size(); i++) {
				compare.array.AddColumn("");
				if (i == 0) {
					String name = compare.opHideEquations ? GetEquationName(param.name) : param.name;
					if (name.Find("[") < 0)
						name = Format("%s [%s]", name, units);
					compare.array.Set(0, col, name);
				}
				compare.array.Set(1, col++, param.metrics[i]);
				fmt << ("." << FormatInt(param.decimals) << "f");
			}
		}
		for (int row = 0; row < table.size(); ++row) {
			for (col = 0; col < table[row].size(); ++col) {
				String val;
				if (col < 3)
					val = Format("%" + fmt[col], S(table[row][col]));
				else {
					double d = double(table[row][col]);
					if (IsNull(d))
						val = "-";
					else
						val = Format("%" + fmt[col], d);
				}
				compare.array.Set(row+2, col, val);
			}
		}
			
	} catch (const Exc &e) {
		BEM::PrintError(e);	
	}
}
	
void FastScatterBase::Init(FastScatter *parent, Function <bool(String)> OnFile, Function <void(String)> OnCopyTabs, StatusBar &_statusBar) {
	fastScatter = parent;
	
	WhenFile = OnFile;
	WhenCopyTabs = OnCopyTabs;
	statusBar = &_statusBar;
	CtrlLayout(*this);
	
	CtrlLayout(rightT);
	CtrlLayout(rightB);
	right.Add(rightT, 0, 0).Add(rightB, 1, 0);
	splitter.Horz(left.SizePos(), right.SizePos());
	splitter.SetPos(7000, 0);

	CtrlLayout(leftSearch);
	CtrlLayout(rightSearch);
	rightB.splitterSearch.Horz(leftSearch.SizePos(), rightSearch.SizePos());
	CtrlLayout(rightUp);
	CtrlLayout(rightDown);
	rightB.splitterDown.Vert(rightUp.SizePos(), rightDown.SizePos());
	rightB.splitterDown.SetPos(7000, 0);
	
	rectangleButtons.SetBackground(SColorFace());
	
	file.WhenChange = [&] {
		return WhenFile(~file);
	};
	file.BrowseRightWidth(40).UseOpenFolder().BrowseOpenFolderWidth(10)
		.Tip(t_("Enter file path to show, or drop it from file explorer"));
	file.Type(t_("All files"), "*.*").Type(t_("FAST output file"), "*.out, *.outb")
									 .Type(t_("CSV file"), "*.csv").Type(t_("DeepLines Wind file"), "*.db")
									 .Type(t_("AQWA Naut file"), "*.lis"); 
	butLoad.Tip(t_("Loads FAST out/outb file")) << [&] {file.DoGo();};
	butSaveAs <<= THISBACK(OnSaveAs);
	butSaveAs.Tip(t_("Saves data file"));
	dropFormat.Add(".out").Add(".csv").Add(".csv only selected");
	dropFormat.SetIndex(0);
	
	opTime <<= 0;
	opTime.WhenAction = [&] {
		ShowSelected(false);
	};
	
	opLoad3 <<= 0;
	opLoad3.WhenAction = [&] {
		if (opLoad3 == 0)	
			right.SetHeight(0, 0);	
		else
			right.SetHeight(0, 200);	
		ShowSelected(true);
	};
	opLoad3.WhenAction();
	
	editStart.WhenEnter = THISBACK1(ShowSelected, true);
	opFrom.WhenAction = THISBACK1(ShowSelected, true);
	editEnd.WhenEnter = THISBACK1(ShowSelected, true);
	
	UpdateButtons(false);
	
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
	
	leftSearch.array.Removing().NoAskRemove().WhenArrayAction = [&] {
		for(int i = 0; i < leftSearch.array.GetCount(); i++)
			if (leftSearch.array.IsSelected(i))
				return;
		ShowSelected(false);
	};
	rightSearch.array.Removing().NoAskRemove().WhenArrayAction = THISBACK1(ShowSelected, false);
	
	leftSearch.array.AddColumn("");
	leftSearch.array.NoVertGrid();
	rightSearch.array.AddColumn("");
	rightSearch.array.NoVertGrid();
	rightSearch.label.SetLabel(t_("Right axis"));
	
	frameSet.Add(leftSearch.array.GetRectEnter(), 1);
	frameSet.Add(rightSearch.array.GetRectEnter(), 1);
	frameSet.Set(leftSearch.array.GetRectEnter());
	
	rightB.filterParam.WhenAction = THISBACK1(OnFilter, true);
	rightB.filterParam.Tip(t_("Filters parameters to display. * are allowed"));
	rightB.filterUnits.WhenAction = THISBACK1(OnFilter, true);
	rightB.filterUnits.Tip(t_("Filters parameters to be displayed according to their units. * are allowed"));
	rightUp.arrayParam.AddColumn(t_("Fields"), 80);
	rightUp.arrayParam.AddColumn(t_("Units"), 20);
	rightUp.arrayParam.SetLineCy(EditField::GetStdHeight()).MultiSelect()
		 .Tip(t_("Double click to choose parameters to display"));
	rightUp.arrayParam.WhenLeftDouble = THISBACK1(WhenArrayLeftDouble, &rightUp.arrayParam);
	rightUp.arrayParam.WhenEnterKey   = THISBACK1(WhenArrayLeftDouble, &rightUp.arrayParam);
	
	rightB.butCopy 		<<= THISBACK(SelCopy);
	rightB.butCopy.Tip(t_("Copy selected parameters to clipboard"));
	rightB.butPaste 		<<= THISBACK1(SelPaste, "");
	rightB.butPaste.Tip(t_("Paste selected parameters from clipboard"));
	rightB.butSetOnAllTabs 	<<= THISBACK(SelCopyTabs);
	rightB.butSetOnAllTabs.Tip(t_("Copy selected parameters on the other tabs"));
	
	rightUp.arrayParam.WhenDropInsert= [&] (int line, PasteClip& d) {OnDropInsert(line, d, rightUp.arrayParam); };
	rightUp.arrayParam.WhenDrag 		= [&] {OnDrag(rightUp.arrayParam, false);};
	
	rightB.opZoomToFit <<= true;
	
	rightDown.arrayPoints.AddColumn(t_("Name"), 80);
	rightDown.arrayPoints.GetColumn(0).Edit(edits.Add()); 
	rightDown.arrayPoints.AddColumn(t_("X"), 80);
	rightDown.arrayPoints.GetColumn(1).Edit(editNums.Add()); 
	rightDown.arrayPoints.AddColumn(t_("Y"), 80);
	rightDown.arrayPoints.GetColumn(2).Edit(editNums.Add()); 
	rightDown.arrayPoints.AddColumn(t_("Z"), 80);
	rightDown.arrayPoints.GetColumn(3).Edit(editNums.Add()); 
	rightDown.arrayPoints.Proportional().Clipboard().Editing().Removing().Appending().Duplicating().SetToolBar();
	rightDown.arrayPoints.WhenInsertRow = [&] {
		for (int r = 0; r < rightDown.arrayPoints.GetRowCount(); ++r) {
			if (IsNull(rightDown.arrayPoints.Get(r, 0)))
				rightDown.arrayPoints.Set(r, 0, FormatInt(r+1));
		}
	};
	
	rightT.arrayFiles.NoHeader().SetLineCy(EditField::GetStdHeight());
	rightT.arrayFiles.AddColumn(t_("Id"), 2);
	rightT.arrayFiles.AddColumn(t_("File"), 20);
	rightT.arrayFiles.WhenSel = [&] {
		int row = rightT.arrayFiles.GetCursor();
		if (row < 0)
			return;
		file <<= rightT.arrayFiles.Get(row, 1);
	};
	
	editStart <<= "-";
	editEnd <<= "-";
	
	opFrom.MinCaseHeight(int(1.2*StdFont().GetCy()));
	opFrom = 0;
	
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
		int seconds = int(StringToSeconds(String(~updateTime)));
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
	String strPoints;
	for (int rw = 0; rw < rightDown.arrayPoints.GetRowCount(); ++rw) {
		for (int c = 0; c < 4; ++c) {
			if (!strPoints.IsEmpty())
				strPoints << ",";
			strPoints << rightDown.arrayPoints.Get(rw, c);
		}
	}
	return strLeft + ";" + strRight + ";" + strPoints;
}

void FastScatterBase::SelCopy() {
	WriteClipboardText(SelectedStr());
}

void FastScatterBase::SelPaste(String str) {
	if (str.IsEmpty())
		str = ReadClipboardText();
	UVector<String> params = Split(str, ";");
	params.SetCount(3);
	UVector<String> leftp = Split(params[0], ",");
	leftSearch.array.Clear();
	for (int rw = 0; rw < leftp.size(); ++rw) {
		if (left.dataFast[0].GetParameterX(leftp[rw]) >= 0)
			leftSearch.array.Set(rw, 0, leftp[rw]);
	}
	UVector<String> rightp = Split(params[1], ",");
	rightSearch.array.Clear();
	for (int rw = 0; rw < rightp.size(); ++rw) {
		if (left.dataFast[0].GetParameterX(rightp[rw]) >= 0)
			rightSearch.array.Set(rw, 0, rightp[rw]);
	}
	UVector<String> points = Split(params[2], ",");
	rightDown.arrayPoints.Clear();
	for (int rw = 0; rw < points.size()/4; ++rw) {
		for (int c = 0; c < 4; ++c)
			rightDown.arrayPoints.Set(rw, c, points[rw*4 + c]);
	}
	ShowSelected(false);	
}

void FastScatterBase::SelCopyTabs() {
	WhenCopyTabs(SelectedStr());
}
	
void FastScatterBase::OnTimer() {
	if (!player.IsRunning())
		return;
	int time = int(StringToSeconds(String(~updateTime)));
	
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
	disableShowSelected = true;
	if(AcceptInternal<ArrayCtrl>(d, "array")) {
		array.InsertDrop(line, d);
		array.SetFocus();
		ArrayRemoveDuplicates(array);
	} else if(AcceptText(d) && !ArrayExists(array, GetString(d))) {
		array.Insert(line);
		array.Set(line, 0, GetString(d));
		array.SetCursor(line);
		array.SetFocus();
	}	
	disableShowSelected = false;
	ShowSelected(true);
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
	
	ShowSelected(false);
}
 
void FastScatterBase::OnFilter(bool show) {
	SortedVectorMap<String, String> list = left.dataFast[0].GetList(Trim(~rightB.filterParam), Trim(~rightB.filterUnits));
	if (list.size() == 1) {
		if (AddParameter(list.GetKey(0), nullptr)) {
			rightB.filterParam.Clear();
			rightB.filterUnits.Clear();
			rightB.filterUnits.Clear();
			list = left.dataFast[0].GetList();
			if (show)
				ShowSelected(false);
		}
	} 
	rightUp.arrayParam.Clear();
	rightB.filterUnits.ClearList();
	for (int rw = 0; rw < list.size(); ++rw) {
		if (ToLower(list.GetKey(rw)) != "time") {
			rightUp.arrayParam.Add(list.GetKey(rw), list[rw]);
			rightB.filterUnits.FindAddList(list[rw]);
		}
	}
}

void FastScatterBase::Clear() {
	file.Clear();
	left.scatter.Clear();
	left.dataFast.Clear();
	left.dataSource.Clear();
	OnFilter(true);
	UpdateButtons(false);
}

void FastScatterBase::UpdateButtons(bool on) {
	butSaveAs.Enable(on);
	dropFormat.Enable(on);
}

void FastScatterBase::OnLoad() {
	for (int r = 0; r < rightT.arrayFiles.GetCount(); ++r)
		OnLoad0(rightT.arrayFiles.Get(r, 1));
}

bool FastScatterBase::OnLoad(String fileName) {
	FastScatterTabs *pf = GetParentCtrlP<FastScatterTabs>(this);
	
	if (pf && !pf->loadingDragDrop) {
		OnLoad();
		for (int r = 0; r < rightT.arrayFiles.GetCount(); ++r) {		
			String file = rightT.arrayFiles.Get(r, 1);
			if (fileName == file || 
				(ForceExt(fileName, "") == ForceExt(file, "") && 
				 PatternMatch(".out*", GetFileExt(fileName)) && 
				 PatternMatch(".out*", GetFileExt(file)))) 
				return true;
		}
	}
	return OnLoad0(fileName);
}
	
bool FastScatterBase::OnLoad0(String fileName0) {
	NON_REENTRANT(false);
	
	try {
		String fileName = FastOut::GetFileToLoad(fileName0);
		if (IsNull(fileName) || IsVoid(fileName)) {
			String error;
			if (IsVoid(fileName))
			 	error = Format(t_("File '%s' not supported"), fileName0);
			else
			 	error = Format(t_("File '%s' not found"), fileName0);
			//Exclamation(DeQtf(error));
			//statusBar->Temporary(error);
			BEM::PrintError(error);
			left.EnableX();
			UpdateButtons(false);
			return false;
		}
		
		int iff = -1;
		if (opLoad3 == 0) {
			if (left.scatterSize() > 0)
				iff = 0;
		} else
			iff = left.Find(fileName);
		
		bool justUpdate = false;
		if (iff < 0) {
			statusBar->Temporary(t_("Loading file"));
			iff = left.AddAll(opLoad3 == 2);
		} else {
			statusBar->Temporary(t_("Updating data"));
			justUpdate = true;
		}
		
		FastOut &fout = left.dataFast[iff];
		
		fout.pointParams.Clear();
		for (int r = 0; r < rightDown.arrayPoints.GetRowCount(); ++r) {
			String name;
			if (IsNull(rightDown.arrayPoints.Get(r, 0)))
				name = FormatInt(r+1);
			else
				name = rightDown.arrayPoints.Get(r, 0);
			Point3D p = Point3D(rightDown.arrayPoints.Get(r, 1), rightDown.arrayPoints.Get(r, 2), rightDown.arrayPoints.Get(r, 3));
			if (!IsNull(p) && fout.FindParam(name) < 0) 
				fout.pointParams << FastOut::PointParam(name, p, &fout);
		}
		
		String ret;
		{
			WaitCursor waitcursor;
			
			left.EnableX(false);
			
			Progress progress(t_("Loading file"), 100);
			
			auto Status = [&](String str, int _pos) {
				if (!IsNull(str))
					progress.SetText(str); 
				progress.SetPos(_pos); 
				return !progress.Canceled();
			};
			
			ret = fout.Load(fileName, Status);
		}
		if (!ret.IsEmpty()) {
			BEM::PrintError(Format(t_("Problem reading file '%s': %s"), ~file, ret));
			left.EnableX();
			UpdateButtons(false);
			return false;
		}
		if (fout.GetNumData() == 0) {
			BEM::PrintError(Format(t_("File '%s' is empty"), ~file));
			left.EnableX();
			UpdateButtons(false);
			return false;
		}
		
		for (int r = 0; r < fout.pointParams.size(); ++r) {	// Update with points from elastodyn.dat
			rightDown.arrayPoints.Set(r, 0, fout.pointParams[r].name);
			rightDown.arrayPoints.Set(r, 1, fout.pointParams[r].pos.x);
			rightDown.arrayPoints.Set(r, 2, fout.pointParams[r].pos.y);
			rightDown.arrayPoints.Set(r, 3, fout.pointParams[r].pos.z);
		}		
	
		/*if (!justUpdate) {
			UArray<ScatterLeft::DataSource> &src = left.dataSource[iff];
			src.SetCount(fout.parameters.size());
		
			for (int c = 0; c < src.size(); ++c) 
				src[c].Init(fout, c);
		
			LoadParams();
			OnFilter(false);
			ShowSelected(true);
			
			rightT.arrayFiles.Add(rightT.arrayFiles.GetCount()+1, fileName);
		} else 
			ShowSelected(false);*/
		
		UArray<ScatterLeft::DataSource> &src = left.dataSource[iff];
		src.SetCount(fout.parameters.size());
	
		for (int c = 0; c < src.size(); ++c) 
			src[c].Init(fout, c);
	
		LoadParams();
		OnFilter(false);
		ShowSelected(true);
		
		if (!justUpdate)
			rightT.arrayFiles.Add(rightT.arrayFiles.GetCount()+1, fileName);
		
	} catch (const Exc &e) {
		BEM::PrintError(Format("Error: %s", e));	
		left.EnableX();
		UpdateButtons(false);
		return false;
	}
	left.EnableX();
	
	UpdateButtons(true);
	return true;
}

void FastScatterBase::OnSaveAs() {
	try {
		if (rightT.arrayFiles.GetCount() == 0) {
			BEM::PrintError(t_("No file to save"));
			return;
		}
			
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
		fs.Set(ForceExtSafer(~file, ext));
		
		int row = rightT.arrayFiles.GetCursor();
		if (row < 0)
			row = 0;
		String fileName = rightT.arrayFiles.Get(row, 1);
		
		if (fs.ExecuteSaveAs(t_("Save register data"))) {
			Progress progress(t_("Saving file"), 100);
			
			auto Status = [&](String str, int _pos) {
				progress.SetText(str); 
				progress.SetPos(_pos); 
				return !progress.Canceled();
			};
			
			if (fileType == ".out") {	
				if (!left.dataFast[row].Save(fs.Get(), Status, ".out")) {
					BEM::PrintError(Format("Problem saving '%s'", fileName));
					return;
				}
			} else if (fileType == ".csv") {	
				if (!left.dataFast[row].Save(fs.Get(), Status, ".csv", ScatterDraw::GetDefaultCSVSeparator())) {
					BEM::PrintError(Format("Problem saving '%s'", fileName));
					return;
				}
			} else {
				int id = opLoad3 == 2 ? row : 0;
				if (!left.scatter[id].SaveToFileData(fs.Get())) {
					BEM::PrintError(Format("Problem saving '%s'", fileName));
					return;
				}
			}
		}
		saveFolder = GetFileFolder(fs.Get());
	} catch (const Exc &e) {
		BEM::PrintError(DeQtf(e));	
	}		
}
	
void FastScatterBase::ShowSelected(bool zoomtofit) {
	NON_REENTRANT_V;
	
	if (disableShowSelected)
		return;
	
	try {
		WaitCursor wait;
	
		if (opLoad3 == 0 && left.dataFast.size() > 0)
			file <<= left.dataFast[0].GetFileName();
		
		int datasize    = opLoad3 == 0 ? min(1, left.dataFast.size()) : left.dataFast.size();
		int scattersize = opLoad3 == 2 ? datasize : min(1, left.scatterSize());
		
		FastScatterTabs *pf = GetParentCtrlP<FastScatterTabs>(this);
		if (pf) {
			FastScatterTabs &f = *pf;
			
			Value key = -1;
			int id = f.tabBar.GetCursor();
			if (id >= 0) {
				key = f.tabBar.GetKey(id);
				
				if (scattersize == 1) {
					String title = left.dataFast[0].GetFileName();
					if (GetUpperFolder(title).IsEmpty())
						title = GetFileTitle(title);
					else
						title = GetFileName(GetUpperFolder(title)) + "/" + GetFileTitle(title);
					f.tabBar.SetValue(key, title);
				} else
					f.tabBar.SetValue(key, t_("Multiple files"));				
			}
		}
		
		if (zoomtofit || ~rightB.opZoomToFit) {
			left.ClearScatter();
			for (int i = 0; i < scattersize; ++i)
				left.AddScatter();
		}
		
		for (int iff = 0; iff < scattersize; ++iff) {
			auto &scat = left.scatter[iff];
			auto &fast = left.dataFast[iff];
			
			if (fast.GetNumData() > 0) {
				scat.SetLabelX(t_("time"));
				scat.RemoveAllSeries().SetSequentialXAll().SetFastViewX();
				scat.SetLeftMargin(8*StdFont().GetHeight());
				scat.SetTopMargin(StdFont().GetHeight());
				scat.SetBottomMargin(4*StdFont().GetHeight());
				bool rightEmpty = rightSearch.array.GetCount() == 0;
				if (!rightEmpty)
					scat.SetRightMargin(8*StdFont().GetHeight());
				else
					scat.SetRightMargin(2*StdFont().GetHeight());
				scat.SetDrawY2Reticle(!rightEmpty).SetDrawY2ReticleNumbers(!rightEmpty);
				if (scattersize > 1) {
					String title = fast.GetFileName();
					if (GetUpperFolder(title).IsEmpty())
						title = GetFileTitle(title);
					else
						title = GetFileName(GetUpperFolder(title)) + "/" + GetFileTitle(title);
					Font f = scat.GetTitleFont();
					f.Height(12);
					scat.SetTitleFont(f).SetTitle(title);	
				}
			}
		}

		if (!left.dataFast.IsEmpty()) {
			double maxTimeStart = First(left.dataFast).GetTimeStart();
			for (int i = 1; i < left.dataFast.size(); ++i) 
				maxTimeStart = max(maxTimeStart, left.dataFast[i].GetTimeStart());
			double numEditStart = StringToSeconds(String(~editStart));
			if (~editStart == "-" || ~editStart == "")  
				timeStart = maxTimeStart;
			else {
				if (!IsNull(numEditStart)) {
					if (numEditStart < 0)
						throw Exc("Wrong start number");
					if (numEditStart < maxTimeStart) 
						editStart <<= timeStart;	
					timeStart = max(maxTimeStart, numEditStart);
				}
			}
			
			double minTimeEnd = First(left.dataFast).GetTimeEnd();
			for (int i = 1; i < left.dataFast.size(); ++i) 
				minTimeEnd = min(minTimeEnd, left.dataFast[i].GetTimeEnd());
			double numEditEnd  = StringToSeconds(String(~editEnd));
			timeEnd = numEditEnd;
			if (~editEnd == "-" || ~editEnd == "") 
				timeEnd = minTimeEnd;
			if (!IsNull(numEditEnd)) {
				if (numEditEnd < 0)
					throw Exc("Wrong end number");
				if (opFrom == 0) {
					if (numEditEnd > minTimeEnd) 
						editEnd <<= minTimeEnd;
					timeEnd = min(minTimeEnd, numEditEnd);
				} else if (opFrom == 1) {
					timeEnd = timeStart + timeEnd;
					if (timeEnd > minTimeEnd) {
						timeEnd = minTimeEnd;
						editEnd <<= minTimeEnd;
					}
				} else {
					timeEnd = minTimeEnd - timeEnd; 
					if (timeEnd > minTimeEnd) {
						timeEnd = minTimeEnd;
						editEnd <<= 0;
					}
				}
				if (timeEnd <= timeStart)
					throw Exc("Series ends before it starts");
			}
		}
		
		double xfactor;
		String unitsx;
		if (opTime == 0) {
			xfactor = 1;
			unitsx = t_("s");
		} else if (opTime == 1) {
			xfactor = 1/60.;
			unitsx = t_("min");
		} else if (opTime == 2) {
			xfactor = 1/3600.;
			unitsx = t_("h");
		} else if (opTime == 3) {
			xfactor = 1/3600./24.;
			unitsx = t_("d");
		}
		
		for (int iff = 0; iff < datasize; ++iff) {
			auto &fast = left.dataFast[iff];
			auto &scat = opLoad3 == 2 ? left.scatter[iff] : left.scatter[0];
			
			if (fast.GetNumData() > 0) {
				int idBegin = fast.GetIdTime(timeStart);
				if (IsNull(idBegin))
					throw Exc(t_("Incorrect start time"));
				
				int id = fast.GetIdTime(timeEnd);	
				if (IsNull(id))
					throw Exc(t_("Incorrect end time"));
				
				int numData = id - idBegin + 1;
				if (numData <= 1) 
					throw Exc(t_("Incorrect start or end times"));
				
				UVector<int> idsx, idsy, idsFixed;
				for (int rw = 0; rw < leftSearch.array.GetCount(); ++rw) {
					String param = Trim(leftSearch.array.Get(rw, 0));
					if (!param.IsEmpty()) {
						int col = fast.GetParameterX(param);
						if (col < 0)
							statusBar->Temporary(Format("Parameter %s does not exist", param));
						else {
							if (left.dataFast.size() > 1 && opLoad3 != 0)
								param = Format("%d.", iff+1) + param;
							scat.AddSeries(fast.dataOut, 0, col, idsx, idsy, idsFixed, false, idBegin, numData, xfactor)
								.NoMark().Legend(param).Units(fast.units[col], unitsx).Stroke(1);	
						}
					}
				}
				for (int rw = 0; rw < rightSearch.array.GetCount(); ++rw) {
					String param = Trim(rightSearch.array.Get(rw, 0));
					if (!param.IsEmpty()) {
						int col = fast.GetParameterX(param);
						if (col < 0)
							statusBar->Temporary(Format("Parameter %s does not exist", param));
						else {
							if (left.dataFast.size() > 1 && opLoad3 != 0)
								param = Format("%d.", iff+1) + param;
							scat.AddSeries(fast.dataOut, 0, col, idsx, idsy, idsFixed, false, idBegin, numData, xfactor)
								.NoMark().Legend(param).Units(fast.units[col], unitsx).SetDataSecondaryY().Stroke(1);	
						}
					}
				}
			}
		}
		
		if (zoomtofit || ~rightB.opZoomToFit) {
			for (int iff = 0; iff < scattersize; ++iff) {
				auto &scat = left.scatter[iff];
				
				scat.ZoomToFit(true, true);	
			}
		}
		
		SaveParams();
	} catch (const Exc &e) {
		BEM::PrintError(DeQtf(e));	
	}
}

void FastScatterTabs::Init(String appDataFolder, StatusBar &_statusBar) {
	statusBar = &_statusBar;
	String folder = AFX(appDataFolder, "FASTScatter");
	if (!DirectoryCreateX(folder))
		return;
		
	LoadFromJsonFile(*this, AFX(folder, "FASTScatter.json"));
	
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

void FastScatterTabs::AddFile(String filename) {
	String title;
	if (filename == t_("New Tab"))
		title = t_("New Tab");
	
	Value key = -1;
	int idKey = -1;
	int id = tabBar.GetCursor();
	if (id >= 0) {
		key = tabBar.GetKey(id);
		idKey = tabKeys.Find(key);
	}
		
	if (id >= 0 && (t_("New Tab") == tabBar.GetValue(id))) { 
		title = GetFileName(GetUpperFolder(filename)) + "/" + GetFileTitle(filename);
		tabBar.SetValue(key, title);
	} else {
		bool addTab = true;
		
		if (idKey >= 0) {
			FastScatterBase &fscbase = tabScatters[idKey].fscbase;
		
			int idLoaded = -1;
			for (int i = 0; i < fscbase.left.dataFast.size(); ++i) {
				if (fscbase.left.dataFast[i].GetFileName() == filename) {
					idLoaded = i;
					break;
				}
			}
			if (fscbase.opLoad3 == 0 && idLoaded < 0)
				addTab = true;
			else if (idLoaded < 0)
				addTab = false;
		}
			
		if (addTab) {
			key = tabCounter++;
			tabKeys << key;
			FastScatter &sct = tabScatters.Add();
			idKey = tabScatters.size()-1;
			int pos = max(tabBar.GetCount()-1, 0);

			sct.Init([=] (String filename) {
					if (tabBar.GetCount() == 0)
						return false; 
					int id = tabBar.GetCursor();
					if (id < 0)
						return false;
					Value key = tabBar.GetKey(id);
					bool ret = false;
					if (filename == t_("New Tab"))
						ret = true;
					else {
						for (int i = 0; i < tabKeys.size(); ++i) {
							if (tabKeys[i] == key ) {
								ret = tabScatters[i].fscbase.OnLoad(filename);
								break;
							}
						}
					}
					if (ret)		
						AddHistory(filename);
					return ret;
				}, [=] (String clipboard) {
					if (tabBar.GetCount() == 0)
						return; 
					int id = tabBar.GetCursor();
					if (id < 0)
						return;
					Value key = tabBar.GetKey(id);
					for (int i = 0; i < tabKeys.size(); ++i) {
						if (tabKeys[i] != key) {
							tabScatters[i].fscbase.SelPaste(clipboard);
							tabScatters[i].fscbase.editStart <<= ~tabScatters[id].fscbase.editStart;
							tabScatters[i].fscbase.editEnd <<= ~tabScatters[id].fscbase.editEnd;
						}
					}
				}, *statusBar);
			tabBar.InsertKey(pos, key, title);
			Add(sct.SizePos());
		}
	}
	AddHistory(filename);
	AddHistoryMetrics("");
	if (int(key) > -1) {// && filename != t_("New Tab")) {
		FastScatter &sct = tabScatters[idKey];
		sct.fscbase.file <<= filename;
		if (!sct.fscbase.file.WhenChange())
			OnCloseTab(key); 
		if (filename == t_("New Tab"))
			sct.fscbase.file <<= "";
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
			tabScatters[it].fscbase.file.ClearHistory();
			for (int i = history.size()-1; i >= 0; --i)
				tabScatters[it].fscbase.file.AddHistory(history[i]);
		}
	}
}

void FastScatterTabs::AddHistoryMetrics(String filename) {
	int id = historyMetrics.Find(filename);
	if (id >= 0)
		historyMetrics.Remove(id);
	historyMetrics.Insert(0, filename);
	for (int it = 0; it < tabBar.GetCount(); ++it) {
		const Value &key = tabBar.GetKey(it);
		if (int(key) > -1 && it < tabScatters.size()) {
			tabScatters[it].fscbase.file.ClearHistory();
			for (int i = historyMetrics.size()-1; i >= 0; --i)
				tabScatters[it].compare.fileMetrics.AddHistory(historyMetrics[i]);
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
			AddFile(t_("New Tab"));
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
	FastScatterBase &sct = tabScatters[idKey].fscbase;
	sct.SaveParams();
	tabKeys.Remove(idKey);	
	tabScatters.Remove(idKey);
	//fileNames.Remove(idKey);
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
	
	UVector<String> filesDrop;
	for (int i = 0; i < files.size(); ++i)
		filesDrop.Append(FastOut::GetFilesToLoad(files[i]));
	
	Sort(filesDrop);
	
	if (filesDrop.size() > 6) {
		if (!PromptYesNo(Format(t_("You are about to open %d files"), filesDrop.size())))
			return false;
	}
	
	loadingDragDrop = true;
	
	for (auto &f : filesDrop) 
		AddFile(f);
	
	loadingDragDrop = false;

	return true;
}

FastScatterTabs::~FastScatterTabs() {
	for (int i = 0; i < tabScatters.size(); ++i)
		tabScatters[i].fscbase.SaveParams();	
	
	String file = AFX(GetAppDataFolder(), "BEMRosetta", "FASTScatter", "FASTScatter.json");	
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

void FastScatterBase::Params::Set(FastScatterBase &parent, ArrayCtrl &aleft, ArrayCtrl &aright) const {
	parent.disableShowSelected = true;
	aleft.Clear();
	for (int rw = 0; rw < left.size(); ++rw)
		aleft.Set(rw, 0, left[rw]);
	aright.Clear();
	for (int rw = 0; rw < right.size(); ++rw)
		aright.Set(rw, 0, right[rw]);
	parent.disableShowSelected = false;
	parent.ShowSelected(false);
}
	
void FastScatterBase::LoadParams() {
	String strpath = SHA1StringS(~file).Left(12);
	String strname = SHA1StringS(GetFileName(~file)).Left(12);
	
	String folder = AFX(GetAppDataFolder(), "BEMRosetta", "FASTScatter");
	if (!DirectoryCreateX(folder))
		return;
	
	String fileName;
	FindFile ffpath(AFX(folder, strpath + "_*.json"));
	if (ffpath) 
		fileName = ffpath.GetPath();
	else {
		FindFile ffname(AFX(folder,  "*_" + strname + ".json"));
		if (ffname) 
			fileName = ffname.GetPath();
	}
	if (fileName.IsEmpty())
		return;

	Params params;
	if (!LoadFromJsonFile(params, fileName))
		return;
	
	params.Set(*this, leftSearch.array, rightSearch.array);
}

void FastScatterBase::SaveParams() {
	String strpath = SHA1StringS(~file).Left(12);
	String strname = SHA1StringS(GetFileName(~file)).Left(12);
	
	String folder = AFX(GetAppDataFolder(), "BEMRosetta", "FASTScatter");
	if (!DirectoryCreateX(folder))
		return;
	
	String fileName = AFX(folder, strpath + "_" + strname + ".json");

	Params params;
	params.Get(leftSearch.array, rightSearch.array);
	StoreAsJsonFile(params, fileName, true);
}
