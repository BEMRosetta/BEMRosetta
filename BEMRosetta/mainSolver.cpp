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

using namespace Upp;

#define IMAGECLASS ImgNemoh
#define IMAGEFILE <BEMRosetta/main.iml>
#include <Draw/iml.h>

#include "main.h"


struct HeadFreqConvert : Convert {
	virtual Value Format(const Value& q) const {
		if (IsNull(q))
			return String();
		return FormatDoubleDecimals(ScanDouble(q.ToString()), 5);
	}
};

void InitGrid(GridCtrl &grid, EditDouble edit[]) {
	grid.Reset();
	grid.Absolute().Editing().Clipboard().Sorting(false);
	for (int i = 0; i < 6; ++i)
		grid.AddColumn(InitCaps(BEM::StrDOF(i))).Edit(edit[i]);
	for (int y = 0; y < 6; ++y)
		for (int x = 0; x < 6; ++x)
			grid.Set(y, x, 0.);
}

MainSolverBody::MainSolverBody() {
	const String meshFiles = ".gdf .dat .stl .pnl .msh .grd";
	String meshFilesAst = clone(meshFiles);
	meshFilesAst.Replace(".", "*.");
	
	meshFile.Type(Format("All supported mesh files (%s)", meshFiles), meshFilesAst);
	meshFile.AllFilesType();
	lidFile.Type(Format("All supported mesh files (%s)", meshFiles), meshFilesAst);
	lidFile.AllFilesType();
	
	x_g <<= 0;
	y_g <<= 0;
	z_g <<= 0;
	
	x_0 <<= 0;
	y_0 <<= 0;
	z_0 <<= 0;
	
	butC0toCg.WhenAction = [&] {
		x_g <<= ~x_0;
		y_g <<= ~y_0;
		z_g <<= ~z_0;
	};
	butCgtoC0.WhenAction = [&] {
		x_0 <<= ~x_g;
		y_0 <<= ~y_g;
		z_0 <<= ~z_g;
	};	
	InitGrid(M, editMass);
	InitGrid(Dlin, editLinear);
	InitGrid(Dquad, editQuadratic);	
	InitGrid(Cadd, editAdd);
}

void MainSolver::Init() {
	tab.Add(genScroll.AddPaneV(gen).SizePos(), "General");
	tab.Add(bodies.SizePos(), "Bodies");
	tab.Add(save.SizePos(), "Save");
	
	WithMainSolver_Body<StaticRect> &b = bodiesEach.Add();
	CtrlScroll &bscroll = bodiesEachScroll.Add();
	CtrlLayout(b);
	bodies.array.Add(bscroll.AddPane(b, true, true).SizePos(), "Body 1");
	b.name <<= "Body 1";
	b.name.WhenAction = [&]() {
		bodies.array.grid.Set(bodies.array.GetCursor(), 0, ~b.name);
	};
	bodies.array.SetWidth(60);
	
	
	const String bemFiles = ".cal .in .dat";
	String bemFilesAst = clone(bemFiles);
	bemFilesAst.Replace(".", "*.");
	loadFrom.Type(Format("All supported bem files (%s)", bemFiles), bemFilesAst);
	loadFrom.AllFilesType();
	String extView = ToLower(GetFileExt(loadFrom.GetData().ToString()));
	if (extView.IsEmpty())
		loadFrom.ActiveType(0);
	else if (bemFiles.Find(extView) >= 0)
		loadFrom.ActiveType(0);
	else
		loadFrom.ActiveType(1);
	
	loadFrom.WhenChange = [&] {return OnLoad();};
	loadFrom.BrowseRightWidth(40).UseOpenFolder(true).BrowseOpenFolderWidth(10).UseDropping();
	butLoad.WhenAction = [&] {loadFrom.DoGo();};

	save.saveTo.BrowseRightWidth(40).UseOpenFolder().BrowseOpenFolderWidth(10).UseDropping();
	save.butSave.WhenAction = [&] {OnSave();};
	
	bodies.butAdd <<= THISBACK(arrayOnAdd);
	bodies.butDuplicate <<= THISBACK(arrayOnDuplicate);
	bodies.butRemove <<= THISBACK(arrayOnRemove);
	
	
	for (int i = 0; i < Hydro::NUMBEM; ++i)
		if (Hydro::caseCanSave[i])
			save.dropSolver.Add(i, Hydro::bemStr[i]);		
			
	save.dropSolver.SetIndex(min(dropSolverVal, save.dropSolver.GetCount()-1));
	save.dropSolver.WhenAction = [&] {
		bool isNemoh = int(~save.dropSolver) >= Hydro::NEMOH && int(~save.dropSolver) <= Hydro::NEMOHv3;
		//if (isNemoh)
			//tab.Set(nemohScroll);
		//else
		//	tab.Set(hamsScroll);
		
		save.numThreads.Enable(!isNemoh);
		save.labnumThreads.Enable(!isNemoh);
		//save.nemohScroll.Enable(isNemoh);
		save.opIncludeBin.Enable(~save.dropSolver != Hydro::CAPYTAINE);
		
		//bodies.array.HeaderObject().ShowTab(1, !isNemoh);
		
		//bodies.lidFile.Enable(!isNemoh);
		
		//bodies.M.Enable(!isNemoh);
		//bodies.Dlin.Enable(!isNemoh);
		//bodies.Dquad.Enable(!isNemoh);
		//bodies.C.Enable(!isNemoh);
		//bodies.Cext.Enable(!isNemoh);
		//bodies.Cadd.Enable(false);		// Not used for now
	};
	save.dropSolver.WhenAction();
	
	gen.opInfinite.WhenAction = [&] {
		gen.height.Enable(!bool(~gen.opInfinite));
	};
	gen.opInfinite.WhenAction();
	
	//gen.gridFreq.AddColumn("w [rad/s]").SetConvert(Single<HeadFreqConvert>()).Edit(editF);
//	gen.gridFreq.Editing().MultiSelect().Removing().Clipboard().Sorting(false);
	gen.gridFreq.WhenPaste = gen.gridFreq.WhenEnter = gen.gridFreq.WhenCursor = gen.gridFreq.WhenRemoveRow = [&]() {
		UVector<double> data;
		for (int i = 0; i < gen.gridFreq.GetCount(); ++i)
			data << gen.gridFreq(i, 0);
		if (data.IsEmpty())
			return;
		Sort(data);
		gen.Nf <<= data.GetCount();
		gen.minF <<= First(data);
		gen.maxF <<= Last(data);
	};
	
	gen.minF.WhenAction = gen.maxF.WhenAction = gen.Nf.WhenAction = [&]() {
		if (IsNull(gen.Nf) || gen.Nf < 1 || IsNull(gen.minF) || gen.minF < 0 || IsNull(gen.maxF) || gen.maxF <= gen.minF)
			return;
		double delta;
		if (gen.Nf == 1)
			delta = gen.maxF - gen.minF;
		else
			delta = (gen.maxF - gen.minF)/(gen.Nf - 1);
		gen.gridFreq.Clear();
		for (int i = 0; i < gen.Nf; ++i)
			gen.gridFreq.Add(gen.minF + i*delta);	
	};
	
	//gen.gridHead.AddColumn("head [deg]").SetConvert(Single<HeadFreqConvert>()).Edit(editH);
	//gen.gridHead.Editing().MultiSelect().Removing().Clipboard().Sorting(false);
	gen.gridHead.WhenPaste = gen.gridHead.WhenEnter = gen.gridHead.WhenCursor = gen.gridHead.WhenRemoveRow = [&]() {
		UVector<double> data;
		for (int i = 0; i < gen.gridHead.GetCount(); ++i)
			data << gen.gridHead(i, 0);
		if (data.IsEmpty())
			return;
		Sort(data);
		gen.Nh <<= data.GetCount();
		gen.minH <<= First(data);
		gen.maxH <<= Last(data);
	};	
	
	gen.minH.WhenAction = gen.maxH.WhenAction = gen.Nh.WhenAction = [&]() {
		if (IsNull(gen.Nh) || gen.Nh < 1 || IsNull(gen.minH) || gen.minH < 0 || IsNull(gen.maxH) || gen.maxH <= gen.minH)
			return;
		double delta;
		if (gen.Nh == 1)
			delta = gen.maxH - gen.minH;
		else
			delta = (gen.maxH - gen.minH)/(gen.Nh - 1);
		gen.gridHead.Clear();
		for (int i = 0; i < gen.Nh; ++i)
			gen.gridHead.Add(gen.minH + i*delta);	
	};
}

void MainSolver::InitBeforeSerialize() {
	CtrlLayout(gen);
	CtrlLayout(bodies);
	CtrlLayout(save);
	CtrlLayout(*this);
	
	gen.gridFreq.AddColumn("w [rad/s]").SetConvert(Single<HeadFreqConvert>()).Edit(editF);
	gen.gridFreq.Editing().MultiSelect().Removing().Clipboard().Sorting(false);
	
	gen.gridHead.AddColumn("head [deg]").SetConvert(Single<HeadFreqConvert>()).Edit(editH);
	gen.gridHead.Editing().MultiSelect().Removing().Clipboard().Sorting(false);
}

void MainSolver::InitSerialize(bool ret) {
	if (!ret || IsNull(save.opIncludeBin)) 
		save.opIncludeBin = true;
	if (!ret || IsNull(save.numSplit)) 	
		save.numSplit = 2;	
	String manufacturer, productName, version, mbSerial;
	int numberOfProcessors;
	GetSystemInfo(manufacturer, productName, version, numberOfProcessors, mbSerial);
	if (!ret || IsNull(save.numThreads)) 	
		save.numThreads = numberOfProcessors;
	if (!ret || IsNull(~gen.xeff))
		gen.xeff <<= 0;
	if (!ret || IsNull(~gen.yeff))
		gen.yeff <<= 0;
	//if (!ret || IsNull(~bodies.cx))
	//	bodies.cx <<= 0;
	//if (!ret || IsNull(~bodies.cy))
	//	bodies.cy <<= 0;
	//if (!ret || IsNull(~bodies.cz))
	//	bodies.cz <<= 0;
	//if (!ret || IsNull(~gen.Nf))
		//gen.Nf <<= 100;
	int numF = gen.gridFreq.GetRowCount();
	if (numF == 0) {
		gen.Nf <<= 100;
		gen.minF <<= 0;
		gen.maxF <<= 4;
	} else {
		gen.Nf <<= numF;
		gen.minF <<= gen.gridFreq(0, 0);
		gen.maxF <<= gen.gridFreq(numF-1, 0);
	}
	//if (!ret || IsNull(~gen.minF))
	//	gen.minF <<= 0;
	//if (!ret || IsNull(~gen.maxF))
	//	gen.maxF <<= 4;
	/*if (!ret || IsNull(~gen.Nh))
		gen.Nh <<= 1;
	if (!ret || IsNull(~gen.minH))
		gen.minH <<= 0;
	if (!ret || IsNull(~gen.maxH))
		gen.maxH <<= 0;*/
	int numH = gen.gridHead.GetRowCount();
	if (numH == 0) {
		gen.Nh <<= 1;
		gen.minH <<= 0;
		gen.maxH <<= 0;
	} else {
		gen.Nh <<= numH;
		gen.minH <<= gen.gridHead(0, 0);
		gen.maxH <<= gen.gridHead(numH-1, 0);
	}
	if (!ret || IsNull(~gen.height)) {
		gen.height <<= Null;
		gen.opInfinite <<= true;
	}
}
	
void MainSolver::Jsonize(JsonIO &json) {
	if (json.IsLoading()) {
		save.opIncludeBin = Null;
		save.numSplit = Null;	
		gen.xeff <<= Null;
		gen.yeff <<= Null;
		//bodies.cx <<= Null;
		//bodies.cy <<= Null;
		//bodies.cz <<= Null;
		//gen.Nf <<= Null;
		//gen.minF <<= Null;
		//gen.maxF <<= Null;
		//gen.Nh <<= Null;
		//gen.minH <<= Null;
		//gen.maxH <<= Null;
	} else {
		dropSolverVal = save.dropSolver.GetIndex();
	}
	json
		("loadFrom", loadFrom)
		("saveTo", save.saveTo)
		("opIncludeBin", save.opIncludeBin)
		("numSplit", save.numSplit)
		("numThreads", save.numThreads)
		("opSplit", save.opSplit)
		("xeff", gen.xeff)
		("yeff", gen.yeff)
		//("cx", bodies.cx)
		//("cy", bodies.cy)
		//("cz", bodies.cz)
		//("Nf", gen.Nf)
		//("minF", gen.minF)
		//("maxF", gen.maxF)
		("gen_gridFreq", gen.gridFreq)
		//("Nh", gen.Nh)
		//("minH", gen.minH)
		//("maxH", gen.maxH)
		("gen_gridHead", gen.gridHead)
		("dropSolver", dropSolverVal)
		("height", gen.height)
		("opInfinite", gen.opInfinite)
		("gen_g", gen.g)
		("gen_rho", gen.rho)
	;
	if (json.IsLoading()) {
		if (IsNull(dropSolverVal) || dropSolverVal < 0)
			dropSolverVal = 0;
	}
}

bool MainSolver::OnLoad() {
	String file = ~loadFrom;
	
	try {
		Load(file);
		save.dropSolver.WhenAction();
	} catch (const Exc &e) {
		BEM::PrintError(DeQtfLf(e));
		return false;
	}
	return true;		
}		

void MainSolver::Load() {
	if (IsNull(~gen.g))
		gen.g <<= Bem().g;
	if (IsNull(~gen.rho))
		gen.rho <<= Bem().rho;	
	if (IsNull(~gen.height)) {
		gen.opInfinite <<= (Bem().depth < 0);
		gen.height.Enable(Bem().depth > 0);
		gen.height <<= (Bem().depth > 0 ? Bem().depth : Null);
	}
}
	
void MainSolver::Load(String file) {
	Hydro data;
	
	Progress progress(t_("Loading BEM files..."), 100); 
	
	data.LoadCase(file, [&](String str, int _pos) {
			progress.SetText(str); 
			progress.SetPos(_pos); 
			return !progress.Canceled();
	});
	bool isNemoh = data.hd().IsNemoh();
	
	gen.opInfinite <<= (data.hd().h < 0);
	gen.height.Enable(data.hd().h > 0);
	gen.height <<= (data.hd().h > 0 ? data.hd().h : Null);
	
	//InitArray(true);
	
	bodies.array.SetCount(data.hd().msh.size());
	bodiesEach.SetCount(data.hd().msh.size());
	for (int i = 0; i < data.hd().msh.size(); ++i) {
		const Mesh &b = data.hd().msh[i];
		WithMainSolver_Body<StaticRect> &d = bodiesEach[i];
		//WithMainSolver_Body<StaticRect> &d0 = *bodies.array.Get(i);
		d.meshFile <<= b.fileName;
		d.lidFile <<= b.lidFile;  
		d.x_0 = b.c0[0];
		d.y_0 = b.c0[1];
		d.z_0 = b.c0[2];
		d.x_g = b.cg[0];
		d.y_g = b.cg[1];
		d.z_g = b.cg[2];
		
		MatrixXdToGridCtrl(d.M, b.M);
		MatrixXdToGridCtrl(d.Dlin, b.Dlin);
		MatrixXdToGridCtrl(d.Dquad, b.Dquad);
		MatrixXdToGridCtrl(d.Cadd, b.Cadd);
	}
	bodies.array.SetCursor(0);
	save.dropSolver.WhenAction();
		
	if (!data.hd().w.IsEmpty()) {
		gen.Nf <<= data.hd().Nf;
		gen.minF <<= First(data.hd().w);
		gen.maxF <<= Last(data.hd().w);
		VectorToGridCtrl(gen.gridFreq, data.hd().w);
	}
	if (!data.hd().head.IsEmpty()) {
		gen.Nh <<= data.hd().Nh;
		gen.minH <<= First(data.hd().head);
		gen.maxH <<= Last(data.hd().head);
		VectorToGridCtrl(gen.gridHead, data.hd().head);
	}
	
	if (isNemoh) 
		save.dropSolver <<= Hydro::NEMOHv115;
	else if (data.hd().solver == Hydro::HAMS)
		save.dropSolver <<= Hydro::HAMS;
	
	gen.g <<= data.hd().g;
	gen.rho <<= data.hd().rho;
	
	gen.xeff <<= data.hd().x_w;
	gen.yeff <<= data.hd().y_w;
	
	gen.boxIrf <<= data.hd().Tirf.size() > 0;
	if (data.hd().Tirf.size() != 0) {
		gen.irfStep <<= (data.hd().Tirf[1] - data.hd().Tirf[0]);
		gen.irfDuration <<= Last(data.hd().Tirf);
	}
}
/*
void MainSolver::InitArray(bool isNemoh) {
	bodies.array.Reset();
	bodies.array.SetLineCy(EditField::GetStdHeight());
	bodies.array.AddColumn("Mesh file", 40);
	bodies.array.AddColumn("Lid file", 40);
	bodies.array.AddColumn("RotX",  10);
	bodies.array.AddColumn("RotY",  10);
	bodies.array.AddColumn("RotZ",  10);
	
	dropSolver.WhenAction();

	InitGrid(bodies.M, editMass);
	InitGrid(bodies.Dlin, editLinear);
	InitGrid(bodies.Dquad, editQuadratic);	
	InitGrid(bodies.C, editInternal);
	//InitGrid(bodies.Cext, editExternal);
	InitGrid(bodies.Cadd, editAdd);
}*/

void MainSolver::LoadMatrix(GridCtrl &grid, const Eigen::MatrixXd &mat) {
	for (int y = 0; y < 6; ++y)
		for (int x = 0; x < 6; ++x)
			grid.Set(y, x, mat(x, y));
}

bool MainSolver::Save(Hydro &hd, bool isNemoh) {
	if (!gen.opInfinite)
		hd.hd().h = ~gen.height;
	else
		hd.hd().h = -1;
	
	hd.hd().msh.SetCount(bodiesEach.size());
	for (int i = 0; i < hd.hd().msh.size(); ++i) {
		Mesh &b = hd.hd().msh[i];
		b.name = ~bodiesEach[i].name;
		b.fileName = bodiesEach[i].meshFile;
		b.lidFile  = bodiesEach[i].lidFile;
		b.c0[0]    = bodiesEach[i].x_0;
		b.c0[1]    = bodiesEach[i].y_0;
		b.c0[2]    = bodiesEach[i].z_0;
		b.cg[0] = bodiesEach[i].x_g;
		b.cg[1] = bodiesEach[i].y_g;
		b.cg[2] = bodiesEach[i].z_g;
		b.M = GridCtrlToMatrixXd(bodiesEach[i].M);
		b.Dlin = GridCtrlToMatrixXd(bodiesEach[i].Dlin);
		b.Dquad = GridCtrlToMatrixXd(bodiesEach[i].Dquad);
		b.Cadd = GridCtrlToMatrixXd(bodiesEach[i].Cadd);
	}
		
	hd.hd().Nf = ~gen.Nf;
	hd.hd().w.SetCount(~gen.Nf);
	for (int i = 0; i < ~gen.gridFreq.GetRowCount(); ++i)
		hd.hd().w[i] = ScanDouble(~gen.gridFreq.Get(i, 0));
	
	hd.hd().Nh = ~gen.Nh;
	hd.hd().head.SetCount(~gen.Nh);
	for (int i = 0; i < ~gen.gridHead.GetRowCount(); ++i)
		hd.hd().head[i] = ScanDouble(~gen.gridHead.Get(i, 0));
	
	hd.hd().g = ~gen.g;
	hd.hd().rho = ~gen.rho;
	
	hd.hd().x_w = ~gen.xeff;
	hd.hd().y_w = ~gen.yeff;
	
	//hd.irf = ~nemoh.boxIrf;
	if (~gen.boxIrf) {
		//hd.irfStep = ~nemoh.irfStep;
		//hd.irfDuration = ~nemoh.irfDuration;
		double irfStep = ~gen.irfStep;
		double irfDuration = ~gen.irfDuration;
		int n = irfDuration/irfStep; 
		LinSpaced(hd.hd().Tirf, n, 0, irfDuration);
	} else
		hd.hd().Tirf.resize(0);
		// hd.irfStep = hd.irfDuration = 0;
	
	return true;
}
/*
void MainSolver::arrayOnCursor() {
	int id = bodies.array.GetCursor();
	if (id < 0)
		return;
	
	bodies.meshFile <<= bodies.array.Get(id, 0);
	bodies.lidFile  <<= bodies.array.Get(id, 1);
	bodies.cx 		<<= bodies.array.Get(id, 2);
	bodies.cy 		<<= bodies.array.Get(id, 3);
	bodies.cz 		<<= bodies.array.Get(id, 4);
}
*/
/*
bool MainSolver::ArrayUpdateCursor() {
	try {
		bool isNemoh = ~dropSolver != Hydro::HAMS;
		
		int id = bodies.array.GetCursor();
		if (id < 0) {
			if (bodies.array.GetCount() == 0) {
				InitArray(isNemoh);
				bodies.array.Add();
				id = 0;
				arrayClear();
			} else
				id = bodies.array.GetCount()-1;
		}	
		
		bodies.array.Set(id, 0, ~bodies.meshFile);
		bodies.array.Set(id, 1, ~bodies.lidFile);
		bodies.array.Set(id, 2, ~bodies.cx);
		bodies.array.Set(id, 3, ~bodies.cy);
		bodies.array.Set(id, 4,~bodies.cz);
	
		bodies.array.Update();
		
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
		return false;
	}
	
	return true;
}
*/
void MainSolver::arrayClear() {
	bodies.array.Clear();
	bodiesEach.Clear();
}

void MainSolver::arrayOnAdd() {
	WithMainSolver_Body<StaticRect> &b = bodiesEach.Add();
	CtrlScroll &bscroll = bodiesEachScroll.Add();
	CtrlLayout(b);
	String name = Format("Body %d", bodiesEach.size());
	bodies.array.Add(bscroll.AddPane(b, true, true).SizePos(), name);
	b.name <<= name;
	b.name.WhenAction = [&]() {
		bodies.array.grid.Set(bodies.array.GetCursor(), 0, ~b.name);
	};
	
/*	if (bodies.array.GetCount() == 0) {
		bool isNemoh = ~dropSolver != Hydro::HAMS;
		InitArray(isNemoh);
	}
	bodies.array.Add();
	bodies.array.SetCursor(bodies.array.GetCount()-1);	
	arrayClear();
	ArrayUpdateCursor();
	dropSolver.WhenAction();*/
}

void MainSolver::arrayOnDuplicate() {
	if (bodies.array.size() == 0) {
		BEM::PrintError(t_("No body available to duplicate"));
		return;
	}
	int id = bodies.array.GetCursor();
	if (id < 0) {
		BEM::PrintError(t_("Please select body to duplicate"));
		return;
	}
	WithMainSolver_Body<StaticRect> &sel = bodiesEach[id];
	
	WithMainSolver_Body<StaticRect> &last = bodiesEach.Add();
	CtrlScroll &bscroll = bodiesEachScroll.Add();
	CtrlLayout(last);
	
	last.name <<= ~sel.name;
	last.x_0 <<= ~sel.x_0;
	last.y_0 <<= ~sel.y_0;
	last.z_0 <<= ~sel.z_0;
	last.x_g <<= ~sel.x_g;
	last.y_g <<= ~sel.y_g;
	last.z_g <<= ~sel.z_g;
	
	last.meshFile <<= ~sel.meshFile;
	last.lidFile <<= ~sel.lidFile;
	
	MatrixXdToGridCtrl(last.M, GridCtrlToMatrixXd(sel.M));
	MatrixXdToGridCtrl(last.Cadd, GridCtrlToMatrixXd(sel.Cadd));
	MatrixXdToGridCtrl(last.Dlin, GridCtrlToMatrixXd(sel.Dlin));
	MatrixXdToGridCtrl(last.Dquad, GridCtrlToMatrixXd(sel.Dquad));
	
	bodies.array.Add(bscroll.AddPane(last, true, true).SizePos(), ~sel.name);
	last.name.WhenAction = [&]() {
		bodies.array.grid.Set(bodies.array.GetCursor(), 0, ~sel.name);
	};	
	
	
/*	int nr = id + 1;
	bodies.array.Insert(nr);
	for (int c = 0; c < bodies.array.GetColumnCount(); ++c)
		bodies.array.Set(nr, c, bodies.array.Get(id, c));
	bodies.array.Disable();
	bodies.array.SetCursor(nr);
	bodies.array.Enable();
	arrayOnCursor();*/
}

void MainSolver::arrayOnRemove() {
	if (bodies.array.size() == 0) {
		BEM::PrintError(t_("No body available to remove"));
		return;
	}
	int id = bodies.array.GetCursor();
	if (id < 0) {
		BEM::PrintError(t_("Please select body to remove"));
		return;
	}
	bodies.array.Remove(id);
	bodiesEach.Remove(id);
	bodiesEachScroll.Remove(id);
	bodies.array.SetCursor(max(id, bodiesEach.size()-1));
	/*if (id >= bodies.array.GetCount()) {
		id = bodies.array.GetCount()-1;
		if (id < 0) {
			arrayClear();
			ArrayUpdateCursor();
			return;
		}
	} 
	bodies.array.SetCursor(id);*/
}

bool MainSolver::OnSave() {
	try {
		String folder = ~save.saveTo;
		
		bool isNemoh = ~save.dropSolver != Hydro::HAMS;
		
		Hydro data;
		
		if (!Save(data, isNemoh))
			return false;
		
		UVector<String> errors = data.Check(~save.dropSolver);
		
		if (!errors.IsEmpty()) {
			String str;
			if (errors.size() == 1)
				str << "\n " << errors[0];
			else {
				for (int i = 0; i < errors.size(); ++i)
				 	str << "\n- " << errors[i];
			}
			if (!ErrorOKCancel(Format(t_("Errors found in data:%s&Do you wish to continue?"), DeQtfLf(str))))
				return false;
		}
		if (!DirectoryExists(folder)) {
			if (!PromptYesNo(Format(t_("Folder %s does not exist.&Do you wish to create it?"), DeQtfLf(folder))))
				return false;
			RealizeDirectory(folder);
		} else {
			if (!PromptYesNo(Format(t_("Folder %s contents will be deleted.&Do you wish to continue?"), DeQtfLf(folder))))
				return false;
		}
		if (~save.opSplit) {
			if (IsNull(~save.numSplit)) {
				BEM::PrintError(t_("Please enter number of parts to split the simulation (min. is 2)"));
				return false;
			} else if (int(~save.numSplit) > data.hd().Nf) {
				if (PromptOKCancel(Format(t_("Number of split cases %d must not be higher than number of frequencies %d"), int(~save.numSplit), data.hd().Nf)
							   + S("&") + t_("Do you wish to fit the number of cases to the number of frequencies?"))) 
					save.numSplit <<= data.hd().Nf;
				else
					return false;
			}
		}
		
		WaitCursor waitcursor;
		
		if (isNemoh) 
			data.SaveFolderCase(folder, ~save.opIncludeBin, ~save.opSplit ? int(~save.numSplit) : 1, Null		 , ~save.dropSolver);
		else
			data.SaveFolderCase(folder, ~save.opIncludeBin, ~save.opSplit ? int(~save.numSplit) : 1, ~save.numThreads, ~save.dropSolver);

	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
		return false;
	}
	
	return true;
}

Eigen::MatrixXd GridCtrlToMatrixXd(const GridCtrl &grid) {
	MatrixXd mat(grid.GetColumnCount(), grid.GetRowCount());
	
	for (int x = 0; x < mat.cols(); ++x)
		for (int y = 0; y < mat.cols(); ++y)
			mat(x, y) = grid.Get(y, x);
	
	return mat;
}

void MatrixXdToGridCtrl(GridCtrl &grid, const Eigen::MatrixXd &mat) {
	for (int x = 0; x < mat.cols(); ++x)
		for (int y = 0; y < mat.cols(); ++y)
			grid.Set(y, x, mat(x, y));
}

void VectorToGridCtrl(GridCtrl &grid, const Upp::Vector<double> &mat) {
	for (int x = 0; x < mat.size(); ++x)
		grid.Set(x, 0, mat[x]);
}