// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <SurfaceCanvas/SurfaceCanvas.h>
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
	grid.Absolute().Editing().Clipboard().Sorting(false).FixedPaste().MultiSelect().SelectRow(false);
	for (int i = 0; i < 6; ++i)
		grid.AddColumn(InitCaps(BEM::StrDOF(i))).Edit(edit[i]);
	for (int y = 0; y < 6; ++y)
		for (int x = 0; x < 6; ++x)
			grid.Set(y, x, 0.);
}

MainSolverBody::MainSolverBody() {
	const String meshFiles = ".gdf .dat .stl .pnl .msh .grd .hst .mesh";
	String meshFilesAst = clone(meshFiles);
	meshFilesAst.Replace(".", "*.");
	
	fileMesh.Type(Format("All supported mesh files (%s)", meshFiles), meshFilesAst);
	fileMesh.AllFilesType();
	fileMesh.WhenChange << [&] {butMesh.WhenAction(); return true;};
	fileLid.Type(Format("All supported mesh files (%s)", meshFiles), meshFilesAst);
	fileLid.AllFilesType();
	fileLid.WhenChange  << [&] {butLid.WhenAction(); return true;};
	
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
	InitGrid(Cmoor, editMoor);
	
	UVector<int> idxs;
	
	butMesh << [&]() {
		Body::Load(mesh, ~fileMesh, Bem().rho, Bem().g, Null, Null, false, idxs);
		SetTexts(true);
	};
	butLid << [&]() {
		Body::Load(lid, ~fileLid, Bem().rho, Bem().g, Null, Null, false, idxs);
		SetTexts();
	};
	butMeshClear << [&]() {
		mesh.Clear();
		SetTexts();
	};
	butLidClear << [&]() {
		lid.Clear();
		SetTexts();
	};
}

void MainSolverBody::SetTexts(bool updateInertia) {
	if (mesh.IsEmpty())
		labMesh.SetText(t_("Not loaded")).SetFont(labMesh.GetFont().Bold(false).Italic(true));
	else
		labMesh.SetText(Format(t_("Panels: %d. Nodes: %d"), mesh.dt.mesh.panels.size(), mesh.dt.mesh.nodes.size())).SetFont(labMesh.GetFont().Bold(true).Italic(false));

	if (lid.IsEmpty())
		labLid.SetText(t_("Not loaded")).SetFont(labMesh.GetFont().Bold(false).Italic(true));
	else
		labLid.SetText(Format(t_("Panels: %d. Nodes: %d"), lid.dt.mesh.panels.size(), lid.dt.mesh.nodes.size())).SetFont(labMesh.GetFont().Bold(true).Italic(false));
	
	if (updateInertia) {
		if (mesh.dt.M.size() == 36) 
			MatrixXdToGridCtrl(M, mesh.dt.M, 6, 6, 0);
		if (!IsNull(mesh.dt.c0)) {
			x_0 <<= mesh.dt.c0.x;
			y_0 <<= mesh.dt.c0.y;
			z_0 <<= mesh.dt.c0.z;
		}
		if (!IsNull(mesh.dt.cg)) {
			x_g <<= mesh.dt.cg.x;
			y_g <<= mesh.dt.cg.y;
			z_g <<= mesh.dt.cg.z;
		}
	}
}

void MainSolver::Init() {
	tab.Add(genScroll.AddPaneV(gen).SizePos(), "1. General");
	tab.Add(bodies.SizePos(), "2. Bodies");
	tab.Add(save.SizePos(), "3. Save");
	
	MainSolverBody &b = bodiesEach.Add();
	CtrlScroll &bscroll = bodiesEachScroll.Add();
	CtrlLayout(b);
	
	bodies.array.Add(bscroll.AddPane(b, true, true).SizePos(), "Body 1");
	b.name <<= "Body 1";
	b.name.WhenAction = [&]() {
		bodies.array.grid.Set(bodies.array.GetCursor(), 0, ~b.name);
	};
	bodies.array.SetWidth(60);
	
	bodies.butAdd <<= THISBACK(arrayOnAdd);
	bodies.butDuplicate <<= THISBACK(arrayOnDuplicate);
	bodies.butRemove <<= THISBACK(arrayOnRemove);
	
	const String bemFiles = ".cal .in .dat .nc .owr .yml .out .1 .2 .3 .3sc .3fk .hst .4 .7 .8 .9 .12d .12s .cfg .frc .pot .mmx .wam";
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
	save.withPotentials.WhenAction = [&] {
		if (save.withPotentials && int(~save.dropSolver) == Hydro::CAPYTAINE_PY)
			save.withMesh <<= true;
	;};
	save.opSplit.WhenAction   = [&] {save.numSplit.Enable(save.opSplit);};
	save.opThreads.WhenAction = [&] {save.numThreads.Enable(!save.opThreads);};
	
	save.arrayDOF.Reset();
	save.arrayDOF.AddColumn(t_("DOF"));
	save.arrayDOF.AddColumn(t_("Enabled")).With( 
		[](One<Ctrl>& x) {
			x.Create<Option>().NoWantFocus();
		}
	);
	
#ifndef flagDEBUG
	save.arrayAdditional.Hide();
	save.labBoxPot.Hide();
	save.withPotentials.Hide();
	save.labAdditional.Hide();
	save.opWaveHeight.Hide();
	save.arrayArea.Hide();
#endif
	
	for (int i = 0; i < 6; ++i)
		save.arrayDOF.Add(InitCaps(BEM::StrDOF(i)), true);
	
	if (save.arrayAdditional.GetColumnCount() != 3) {
		save.arrayAdditional.Appending().Removing().Editing().Sorting(true).MultiSelect().Clipboard().ExtraPaste().SelectRow(false).Duplicating().SetToolBar();
		save.arrayAdditional.AddColumn(t_("x")).Edit(edAdditionalX);
		save.arrayAdditional.AddColumn(t_("y")).Edit(edAdditionalY);
		save.arrayAdditional.AddColumn(t_("z")).Edit(edAdditionalZ);
	}
	if (save.arrayArea.GetColumnCount() != 3) {
		save.arrayArea.Editing().MultiSelect().Clipboard().Sorting(false).ExtraPaste().SelectRow(false);
		save.arrayArea.AddColumn(t_("Data"), 30);
		save.arrayArea.AddColumn(t_("x"), 10).Edit(edAreaX);
		save.arrayArea.AddColumn(t_("y"), 10).Edit(edAreaY);
		save.arrayArea.Add(t_("Left/Top [m]"));
		save.arrayArea.Add(t_("Width/Height [m]"));
		save.arrayArea.Add(t_("Number of points"));
	}
	save.opWaveHeight.WhenAction = [&]{
		bool enabled = save.opWaveHeight.IsEnabled() && save.opWaveHeight;
		save.arrayArea.Enable(enabled);
	};
	save.opWaveHeight.WhenAction();
	
	for (int i = 0; i < Hydro::NUMBEM; ++i)
		if (Hydro::caseCanSave[i])
			save.dropSolver.Add(i, Hydro::GetBemStrCase(static_cast<Hydro::BEM_FMT>(i)));		

	save.dropSolver.SetIndex(min(dropSolverVal, save.dropSolver.GetCount()-1));
	save.dropSolver.WhenAction = [&] {
		int solver = ~save.dropSolver;
		bool isNemoh = solver >= Hydro::NEMOH && solver <= Hydro::NEMOHv3;
		
		save.symX.Enable(!isNemoh && solver != Hydro::BEMROSETTA_H5);
		save.symY.Enable(solver != Hydro::BEMROSETTA_H5);
		save.opIncludeBin.Enable(isNemoh || solver == Hydro::HAMS || solver == Hydro::ORCAWAVE_YML);
		save.opSplit.Enable(isNemoh || solver == Hydro::HAMS);
		save.numSplit.Enable(save.opSplit.IsEnabled() && save.opSplit);
		save.opThreads.Enable(solver == Hydro::ORCAWAVE_YML || solver == Hydro::AQWA_DAT || 
							  solver == Hydro::CAPYTAINE_PY || solver == Hydro::WAMIT || 
							  solver == Hydro::HAMS);
		save.numThreads.Enable(save.opThreads.IsEnabled() && !save.opThreads);
		
		save.labDOF.  Enable(isNemoh);// || solver == Hydro::WAMIT || solver == Hydro::CAPYTAINE_PY);
		save.arrayDOF.Enable(isNemoh);// || solver == Hydro::WAMIT || solver == Hydro::CAPYTAINE_PY);
		
		save.withMesh.Enable(solver == Hydro::CAPYTAINE_PY);
		save.withPotentials.Enable(solver == Hydro::ORCAWAVE_YML || solver == Hydro::AQWA_DAT || solver == Hydro::CAPYTAINE_PY || solver == Hydro::HAMS || solver == Hydro::WAMIT);
		save.arrayAdditional.Enable(solver == Hydro::ORCAWAVE_YML || solver == Hydro::AQWA_DAT || solver == Hydro::CAPYTAINE_PY || solver == Hydro::HAMS || solver == Hydro::WAMIT);
		save.opWaveHeight.Enable(solver == Hydro::WAMIT ||  solver == Hydro::HAMS);
		save.opWaveHeight.WhenAction();
		save.withQTF.Enable(solver == Hydro::ORCAWAVE_YML || solver == Hydro::AQWA_DAT || solver == Hydro::WAMIT);
	};
	save.dropSolver.WhenAction();
	
	gen.opInfinite.WhenAction = [&] {
		gen.height.Enable(!bool(~gen.opInfinite));
	};
	gen.opInfinite.WhenAction();
	
	gen.opFreq.WhenAction = [&] {
		switch (gen.opFreq) {
		case 0:		gen.listFreq.Set("Wave periods", "T", "s");			break;
		case 1:		gen.listFreq.Set("Wave frequencies", "ω", "rad/s");	break;
		default:	gen.listFreq.Set("Wave frequencies", "ω", "Hz");
		}		
	};
	gen.opFreq.WhenAction();
}

void MainSolver::InitBeforeSerialize() {
	gen.listFreq.InitBeforeSerialize("Wave frequencies", "ω", "rad/s", 0.1, 6, 100);
	gen.listHead.InitBeforeSerialize("Wave directions/headings", "Head", "deg", 0, 180, 5);
	
	CtrlLayout(gen);
	CtrlLayout(bodies);
	CtrlLayout(save);
	CtrlLayout(*this);
	
	gen.listFreq.from.Min(0.001);
	gen.listHead.from.Min(-360);
	gen.listHead.from.Max(360);
}

void MainSolver::InitAfterSerialize(bool ret) {
	gen.listFreq.InitAfterSerialize();
	gen.listHead.InitAfterSerialize();
	
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
		("gen_listFreq", gen.listFreq)
		("gen_listHead", gen.listHead)
		("dropSolver", dropSolverVal)
		("height", gen.height)
		("opInfinite", gen.opInfinite)
		("gen_g", gen.g)
		("gen_rho", gen.rho)
		("gen_opFreq", gen.opFreq)
	;
	if (json.IsLoading()) {
		if (IsNull(dropSolverVal) || dropSolverVal < 0)
			dropSolverVal = 0;
		if (IsNull(gen.opFreq))
			gen.opFreq <<= 1;
	}
}

bool MainSolver::OnLoad() {
	String file = ~loadFrom;
	
	try {
		Load(file);
		save.dropSolver.WhenAction();
	} catch (const Exc &e) {
		BEM::PrintError(e);
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
		if (gen.opInfinite == false)
			gen.opInfinite <<= (Bem().depth < 0);
		gen.height.Enable(Bem().depth > 0);
		gen.height <<= (Bem().depth > 0 ? Bem().depth : Null);
	}
}
	
void MainSolver::Load(String file) {
	Hydro tmp_hy;											// Just for loading and transferring to gen and bodies
	
	Progress progress(t_("Loading BEM files..."), 100); 

	gen.opFreq <<= 1;

	WaitCursor waitcursor;
	
	tmp_hy.LoadCase(file, [&](String str, int _pos) {
			progress.SetText(str); 
			progress.SetPos(_pos); 
			return !progress.Canceled();
	});
	
	gen.opInfinite <<= (tmp_hy.dt.h < 0);
	gen.height.Enable(tmp_hy.dt.h > 0);
	gen.height <<= (tmp_hy.dt.h > 0 ? tmp_hy.dt.h : Null);
	
	UVector<int> idxs;
	bodies.array.Clear();
	bodiesEach.Clear();
	bodiesEachScroll.Clear();
	for (int ib = 0; ib < tmp_hy.dt.msh.size(); ++ib) {
		Body &tmp_b = tmp_hy.dt.msh[ib];
		MainSolverBody &b = bodiesEach.Add();
		CtrlScroll &bscroll = bodiesEachScroll.Add();
		CtrlLayout(b);
		bodies.array.Add(bscroll.AddPane(b, true, true).SizePos(), tmp_b.dt.name);
		
		b.name <<= tmp_b.dt.name;
		b.fileMesh <<= tmp_b.dt.fileName;
		b.fileLid <<= tmp_b.dt.lidFile;  
		
		if (tmp_b.IsEmpty())
			Body::Load(b.mesh, tmp_b.dt.fileName, tmp_hy.dt.rho, tmp_hy.dt.g, Null, Null, false, idxs);
		else
			b.mesh = clone(tmp_b);

		Body::Load(b.lid, tmp_b.dt.lidFile, tmp_hy.dt.rho, tmp_hy.dt.g, Null, Null, false, idxs);
		
		b.SetTexts();
		
		b.x_0 = tmp_b.dt.c0.x;
		b.y_0 = tmp_b.dt.c0.y;
		b.z_0 = tmp_b.dt.c0.z;
		b.x_g = Nvl(tmp_b.dt.cg.x, b.mesh.dt.cg.x);
		b.y_g = Nvl(tmp_b.dt.cg.y, b.mesh.dt.cg.y);
		b.z_g = Nvl(tmp_b.dt.cg.z, b.mesh.dt.cg.z);
		
		MatrixXdToGridCtrl(b.M, tmp_b.dt.M, 6, 6, 0);
		MatrixXdToGridCtrl(b.Dlin, tmp_b.dt.Dlin, 6, 6, 0);
		MatrixXdToGridCtrl(b.Dquad, tmp_b.dt.Dquad, 6, 6, 0);
		MatrixXdToGridCtrl(b.Cadd, tmp_b.dt.Cadd, 6, 6, 0);
		MatrixXdToGridCtrl(b.Cmoor, tmp_b.dt.Cmoor, 6, 6, 0);
	}
	
	bodies.array.SetCursor(0);
	save.dropSolver.WhenAction();
		
	if (!tmp_hy.dt.w.IsEmpty()) {
		gen.listFreq.number <<= tmp_hy.dt.Nf;
		gen.listFreq.from <<= First(tmp_hy.dt.w);
		gen.listFreq.to <<= Last(tmp_hy.dt.w);
		VectorToGridCtrl(gen.listFreq.grid, tmp_hy.dt.w, 0, Null);
	}
	if (!tmp_hy.dt.head.IsEmpty()) {
		gen.listHead.number <<= tmp_hy.dt.Nh;
		gen.listHead.from <<= First(tmp_hy.dt.head);
		gen.listHead.to <<= Last(tmp_hy.dt.head);
		VectorToGridCtrl(gen.listHead.grid, tmp_hy.dt.head, 0, Null);
	}
	
	if (save.dropSolver.HasKey(tmp_hy.dt.solver))
		save.dropSolver.SetData(tmp_hy.dt.solver);
	
	gen.g <<= tmp_hy.dt.g;
	gen.rho <<= tmp_hy.dt.rho;
	
	gen.xeff <<= tmp_hy.dt.x_w;
	gen.yeff <<= tmp_hy.dt.y_w;
	
	/*gen.boxIrf <<= tmp_hy.dt.Tirf.size() > 0;
	if (tmp_hy.dt.Tirf.size() != 0) {
		gen.irfStep <<= (tmp_hy.dt.Tirf[1] - tmp_hy.dt.Tirf[0]);
		gen.irfDuration <<= Last(tmp_hy.dt5.Tirf);
	}*/
	// To be improved by extracting panel points and free surface points
	save.arrayAdditional.Clear();
	for (const auto &p : tmp_hy.listPoints)
		save.arrayAdditional.Add(p.x, p.y, p.z);
}

void MainSolver::LoadMatrix(GridCtrl &grid, const Eigen::MatrixXd &mat) {
	for (int y = 0; y < 6; ++y)
		for (int x = 0; x < 6; ++x)
			grid.Set(y, x, mat(x, y));
}

bool MainSolver::CopyHydro(Hydro &hy, UArray<Body> &lids) {
	if (!gen.opInfinite)
		hy.dt.h = ~gen.height;
	else
		hy.dt.h = -1;
	
	hy.dt.msh.SetCount(bodiesEach.size());
	lids.SetCount(bodiesEach.size());
	hy.dt.Nb = hy.dt.msh.size();
	
	for (int i = 0; i < hy.dt.Nb; ++i) {
		Body &b = hy.dt.msh[i];
		b = clone(bodiesEach[i].mesh);
		lids[i] = clone(bodiesEach[i].lid);
		b.dt.name = ~bodiesEach[i].name;
		b.dt.fileName = bodiesEach[i].fileMesh;
		b.dt.lidFile  = bodiesEach[i].fileLid;
		b.dt.c0[0]    = bodiesEach[i].x_0;
		b.dt.c0[1]    = bodiesEach[i].y_0;
		b.dt.c0[2]    = bodiesEach[i].z_0;
		b.dt.cg[0]    = bodiesEach[i].x_g;
		b.dt.cg[1]    = bodiesEach[i].y_g;
		b.dt.cg[2]    = bodiesEach[i].z_g;
		b.dt.M 		  = GridCtrlToMatrixXd(bodiesEach[i].M);
		b.dt.Dlin 	  = GridCtrlToMatrixXd(bodiesEach[i].Dlin);
		b.dt.Dquad 	  = GridCtrlToMatrixXd(bodiesEach[i].Dquad);
		b.dt.Cadd 	  = GridCtrlToMatrixXd(bodiesEach[i].Cadd);
		b.dt.Cmoor 	  = GridCtrlToMatrixXd(bodiesEach[i].Cmoor);
	}
	
	hy.dt.Nf = gen.listFreq.number;
	hy.dt.w.SetCount(gen.listFreq.number);
	for (int i = 0; i < gen.listFreq.number; ++i) {
		double dat = ScanDouble(~gen.listFreq.grid.Get(i, 0));;
		if (gen.opFreq == 0)
			dat = 2*M_PI/dat;
		else if (gen.opFreq == 2)
			dat /= 2*M_PI;
		hy.dt.w[i] = dat;	// Always in rad/s
	}
	Sort(hy.dt.w);
	
	hy.dt.Nh = gen.listHead.number;
	hy.dt.head.SetCount(gen.listHead.number);
	for (int i = 0; i < gen.listHead.number; ++i)
		hy.dt.head[i] = ScanDouble(~gen.listHead.grid.Get(i, 0));
	
	hy.dt.g = ~gen.g;
	hy.dt.rho = ~gen.rho;
	
	hy.dt.x_w = ~gen.xeff;
	hy.dt.y_w = ~gen.yeff;
	
	//hd.irf = ~nemoh.boxIrf;
	/*if (~gen.boxIrf) {
		//hd.irfStep = ~nemoh.irfStep;
		//hd.irfDuration = ~nemoh.irfDuration;
		double irfStep = ~gen.irfStep;
		double irfDuration = ~gen.irfDuration;
		int n = int(irfDuration/irfStep); 
		LinSpaced(hy.dt.Tirf, n, 0, irfDuration);
	} else*/
		hy.dt.Tirf.resize(0);
		// hd.irfStep = hd.irfDuration = 0;
	
	return true;
}

void MainSolver::LoadDragDrop() {
	GuiLock __;
	
	Sort(filesToDrop);
	for (int i = filesToDrop.size()-1; i > 0; --i)
		if (ToLower(GetFileTitle(filesToDrop[i])) == ToLower(GetFileTitle(filesToDrop[i-1])))
			filesToDrop.Remove(i);
		
	bool followWithErrors = false;
	for (int i = 0; i < filesToDrop.size(); ++i) {
		loadFrom <<= filesToDrop[i];
		Status(Format(t_("Loading '%s'"), filesToDrop[i]));
		if (!OnLoad() && !followWithErrors && filesToDrop.size() - i > 1) {
			if (!PromptYesNo(Format(t_("Do you wish to try with the pending %d files?"), filesToDrop.size() - i - 1)))
				return;
			followWithErrors = true;
		}
		ProcessEvents();
	}
}
	
void MainSolver::DragAndDrop(Point , PasteClip& d) {
	GuiLock __;
	if (IsDragAndDropSource())
		return;
	if (AcceptFiles(d)) {
		filesToDrop = GetFiles(d);
		timerDrop.Set(0, [=] {LoadDragDrop();});
		return;
	}
	timerDrop.Kill();
}

bool MainSolver::Key(dword key, int ) {
	GuiLock __;
	if (key == K_CTRL_V) {
		filesToDrop = GetFiles(Ctrl::Clipboard());
		timerDrop.Set(0, [=] {LoadDragDrop();});
		return true;
	}
	return false;
}

void MainSolver::arrayClear() {
	bodies.array.Clear();
	bodiesEach.Clear();
}

void MainSolver::arrayOnAdd() {
	MainSolverBody &b = bodiesEach.Add();
	CtrlScroll &bscroll = bodiesEachScroll.Add();
	CtrlLayout(b);
	String name = Format("Body %d", bodiesEach.size());
	bodies.array.Add(bscroll.AddPane(b, true, true).SizePos(), name);
	b.name <<= name;
	b.name.WhenAction = [&]() {
		bodies.array.grid.Set(bodies.array.GetCursor(), 0, ~b.name);
	};
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
	MainSolverBody &sel = bodiesEach[id];
	
	MainSolverBody &last = bodiesEach.Add();
	CtrlScroll &bscroll = bodiesEachScroll.Add();
	CtrlLayout(last);
	
	last.name <<= ~sel.name;
	last.x_0 <<= ~sel.x_0;
	last.y_0 <<= ~sel.y_0;
	last.z_0 <<= ~sel.z_0;
	last.x_g <<= ~sel.x_g;
	last.y_g <<= ~sel.y_g;
	last.z_g <<= ~sel.z_g;
	
	last.fileMesh <<= ~sel.fileMesh;
	last.fileLid <<= ~sel.fileLid;
	
	MatrixXdToGridCtrl(last.M, GridCtrlToMatrixXd(sel.M), 6, 6, 0);
	MatrixXdToGridCtrl(last.Cadd, GridCtrlToMatrixXd(sel.Cadd), 6, 6, 0);
	MatrixXdToGridCtrl(last.Cmoor, GridCtrlToMatrixXd(sel.Cmoor), 6, 6, 0);
	MatrixXdToGridCtrl(last.Dlin, GridCtrlToMatrixXd(sel.Dlin), 6, 6, 0);
	MatrixXdToGridCtrl(last.Dquad, GridCtrlToMatrixXd(sel.Dquad), 6, 6, 0);
	
	last.mesh = clone(sel.mesh);
	last.lid = clone(sel.lid);
	
	last.SetTexts();
	
	bodies.array.Add(bscroll.AddPane(last, true, true).SizePos(), ~sel.name);
	last.name.WhenAction = [&]() {
		bodies.array.grid.Set(bodies.array.GetCursor(), 0, ~sel.name);
	};	
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
}

bool MainSolver::OnSave() {
	try {
		String folder = ~save.saveTo;
		
		if (save.withPotentials)
			save.withMesh <<= true;
		
		Hydro hy;
		UArray<Body> lids;
		
		if (!CopyHydro(hy, lids))
			return false;
		
		Hydro::BEM_FMT solver = (Hydro::BEM_FMT)int(~save.dropSolver);
		
		UVector<String> errors = hy.Check(solver);
		
		if (!errors.IsEmpty()) {
			String str;
			if (errors.size() == 1)
				str << "\n " << errors[0];
			else {
				for (int i = 0; i < errors.size(); ++i)
				 	str << "\n- " << errors[i];
			}
			if (!PromptOKCancel(Format(t_("Problems found in data:%s&Do you wish to continue?"), DeQtfLf(str))))
				return false;
		}
		if (!DirectoryExists(folder))
			RealizeDirectory(folder);
		else {
			if (!PromptYesNo(Format(t_("Folder %s contents will be overwritten.&Do you wish to continue?"), DeQtfLf(folder))))
				return false;
		}
		if (~save.opSplit) {
			if (IsNull(~save.numSplit)) {
				BEM::PrintError(t_("Please enter number of parts to split the simulation (min. is 2)"));
				return false;
			} else if (int(~save.numSplit) > hy.dt.Nf) {
				if (PromptOKCancel(Format(t_("Number of split cases %d must not be higher than number of frequencies %d"), int(~save.numSplit), hy.dt.Nf)
							   + S("&") + t_("Do you wish to fit the number of cases to the number of frequencies?"))) 
					save.numSplit <<= hy.dt.Nf;
				else
					return false;
			}
		}
		
		UVector<bool> listDOF(6);
		for (int r = 0; r < 6; ++r) 
			listDOF[r] = save.arrayDOF.Get(r, 1);
		
		// Listing all the points to be measured
		UVector<Point3D> listPoints;
		for (int r = 0; r < save.arrayAdditional.GetRowCount(); ++r)
			listPoints << Point3D(save.arrayAdditional.Get(r, 0),
								  save.arrayAdditional.Get(r, 1),
								  save.arrayAdditional.Get(r, 2));
		double x0 = save.arrayArea.Get(0, 1),
			   y0 = save.arrayArea.Get(0, 2);
		double dx = save.arrayArea.Get(1, 1),
			   dy = save.arrayArea.Get(1, 2);
		int    nx = save.arrayArea.Get(2, 1),
			   ny = save.arrayArea.Get(2, 2);
		
		if (~save.opWaveHeight) {
			if (nx <= 1 || ny <= 1) {
				BEM::PrintError(t_("Number of points has to be higher than 1"));
				return false;
			}
			
			VectorXd xpoints, ypoints;
			LinSpaced(xpoints, nx, x0, x0+dx);
			LinSpaced(ypoints, ny, y0, y0+dy);
			for (double x : xpoints)
			 	for (double y : ypoints)
					listPoints << Point3D(x, y, 0);
		}
		WaitCursor waitcursor;
		
		int nSplit = 1;
		if (save.opSplit.IsEnabled() && save.opSplit)
			nSplit = save.numSplit;
		
		int nThreads = Null;
		if (save.opThreads.IsEnabled() && !save.opThreads)
			nThreads = save.numThreads;
		
		hy.SaveFolderCase(folder, ~save.opIncludeBin, nSplit, nThreads, solver, 
			~save.withPotentials, ~save.withMesh, ~save.withQTF, ~save.symY, ~save.symX, lids, listDOF, listPoints);

	} catch (Exc e) {
		BEM::PrintError(e);
		return false;
	}
	
	return true;
}

Eigen::MatrixXd GridCtrlToMatrixXd(const GridCtrl &grid) {
	MatrixXd mat(grid.GetColumnCount(), grid.GetRowCount());
	
	for (int x = 0; x < mat.cols(); ++x)
		for (int y = 0; y < mat.rows(); ++y)
			mat(y, x) = grid.Get(y, x);
	
	return mat;
}

void MatrixXdToGridCtrl(GridCtrl &grid, const Eigen::MatrixXd &mat, int rows, int cols, double val) {
	grid.Clear(false);
	for (int x = 0; x < cols; ++x)
		for (int y = 0; y < rows; ++y)
			grid.Set(y, x, val);
	for (int x = 0; x < mat.cols(); ++x)
		for (int y = 0; y < mat.cols(); ++y)
			grid.Set(y, x, mat(y, x));
}

void VectorToGridCtrl(GridCtrl &grid, const Upp::Vector<double> &mat, int rows, double val) {
	grid.Clear(false);
	for (int y = 0; y < rows; ++y)
		grid.Set(y, 0, val);
	for (int y = 0; y < mat.size(); ++y)
		grid.Set(y, 0, mat[y]);
}