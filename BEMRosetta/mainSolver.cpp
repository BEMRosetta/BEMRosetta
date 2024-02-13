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


void MainSolver::Init(const BEM &bem) {
	CtrlLayout(bodies);
	bodiesScroll.AddPaneV(bodies).SizePos();
	CtrlLayout(nemoh);
	CtrlLayout(hams);
	CtrlLayout(*this);
	
	tab.Add(nemohScroll.AddPaneV(nemoh).SizePos(), "Nemoh");
	//tab.Add(hamsScroll.AddPaneV(hams).SizePos(),   "HAMS");
	
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
	
	loadFrom.WhenChange = [&] {return OnLoad(bem);};
	loadFrom.BrowseRightWidth(40).UseOpenFolder(true).BrowseOpenFolderWidth(10).UseDropping();
	butLoad.WhenAction = [&] {loadFrom.DoGo();};

	saveTo.BrowseRightWidth(40).UseOpenFolder().BrowseOpenFolderWidth(10).UseDropping();
	butSave.WhenAction = [&] {OnSave(bem);};
	
	bodies.array.WhenSel = THISBACK(arrayOnCursor);
	
	bodies.butAdd <<= THISBACK(arrayOnAdd);
	bodies.butDuplicate <<= THISBACK(arrayOnDuplicate);
	bodies.butRemove <<= THISBACK(arrayOnRemove);
	bodies.butAllDOF.WhenAction = [&] {
		bool set = bodies.butAllDOF.GetLabel() == t_("All DOF");
		if (set)
			bodies.butAllDOF.SetLabel(t_("No DOF"));
		else
			bodies.butAllDOF.SetLabel(t_("All DOF"));
		bodies.surge <<= set; 	
		bodies.sway  <<= set; 	
		bodies.heave <<= set; 
		bodies.roll  <<= set; 	
		bodies.pitch <<= set; 
		bodies.yaw   <<= set;
		ArrayUpdateCursor();
	};
	
	const String meshFiles = ".gdf .dat .stl .pnl .msh";
	String meshFilesAst = clone(meshFiles);
	meshFilesAst.Replace(".", "*.");
	
	bodies.meshFile.WhenAction 	= [&] {ArrayUpdateCursor();};
	bodies.meshFile.WhenChange 	= [&] {return ArrayUpdateCursor();};
	bodies.meshFile.Type(Format("All supported mesh files (%s)", meshFiles), meshFilesAst);
	bodies.meshFile.AllFilesType();
	bodies.lidFile.WhenAction  	= [&] {ArrayUpdateCursor();};
	bodies.lidFile.WhenChange  	= [&] {return ArrayUpdateCursor();};
	bodies.lidFile.Type(Format("All supported mesh files (%s)", meshFiles), meshFilesAst);
	bodies.lidFile.AllFilesType();
	
	bodies.xcm <<= 0;
	bodies.ycm <<= 0;
	bodies.zcm <<= 0;
	
	bodies.surge.WhenAction 	= [&] {ArrayUpdateCursor();};
	bodies.sway.WhenAction 		= [&] {ArrayUpdateCursor();};
	bodies.heave.WhenAction 	= [&] {ArrayUpdateCursor();};
	bodies.roll.WhenAction 		= [&] {ArrayUpdateCursor();};
	bodies.pitch.WhenAction 	= [&] {ArrayUpdateCursor();};
	bodies.yaw.WhenAction 		= [&] {ArrayUpdateCursor();};
	bodies.cx.WhenAction 		= [&] {ArrayUpdateCursor();};
	bodies.cy.WhenAction 		= [&] {ArrayUpdateCursor();};
	bodies.cz.WhenAction 		= [&] {ArrayUpdateCursor();};
	
	bodies.butC0toCg.WhenAction = [&] {
		bodies.xcm <<= ~bodies.cx;
		bodies.ycm <<= ~bodies.cy;
		bodies.zcm <<= ~bodies.cz;
	};
	bodies.butCgtoC0.WhenAction = [&] {
		bodies.cx <<= ~bodies.xcm;
		bodies.cy <<= ~bodies.ycm;
		bodies.cz <<= ~bodies.zcm;
	};
	
	nemoh.freeSurface.Transparent(false);
	nemoh.freeSurface.WhenAction = [&] {
		nemoh.freeX.Enable(~nemoh.freeSurface);	
		nemoh.freeY.Enable(~nemoh.freeSurface);
	};
	nemoh.freeSurface.WhenAction();
	
	for (int i = 0; i < BEMCase::NUMSOLVERS; ++i)
		if (BEMCase::solverCanSave[i])
			dropSolver.Add(i, BEMCase::solverStr[i]);		
			
	dropSolver.SetIndex(dropSolverVal);
	dropSolver.WhenAction = [&] {
		bool isNemoh = ~dropSolver != BEMCase::HAMS;
		if (isNemoh)
			tab.Set(nemohScroll);
		//else
		//	tab.Set(hamsScroll);
		
		numThreads.Enable(!isNemoh);
		labnumThreads.Enable(!isNemoh);
		nemohScroll.Enable(isNemoh);
		opIncludeBin.Enable(~dropSolver != BEMCase::CAPYTAINE);
		
		bodies.array.HeaderObject().ShowTab(1, !isNemoh);
		bodies.array.HeaderObject().ShowTab(2, isNemoh);
		bodies.array.HeaderObject().ShowTab(3, isNemoh);
		bodies.array.HeaderObject().ShowTab(4, isNemoh);
		bodies.array.HeaderObject().ShowTab(5, isNemoh);
		bodies.array.HeaderObject().ShowTab(6, isNemoh);
		bodies.array.HeaderObject().ShowTab(7, isNemoh);
		
		bodies.lidFile.Enable(!isNemoh);
		bodies.butAllDOF.Enable(isNemoh);
		bodies.surge.Enable(isNemoh);
		bodies.sway.Enable(isNemoh);
		bodies.heave.Enable(isNemoh);
		bodies.roll.Enable(isNemoh);
		bodies.pitch.Enable(isNemoh);
		bodies.yaw.Enable(isNemoh);
		
		//bodies.xcm.Enable(!isNemoh);
		//bodies.ycm.Enable(!isNemoh);
		//bodies.zcm.Enable(!isNemoh);
		bodies.M.Enable(!isNemoh);
		bodies.Dlin.Enable(!isNemoh);
		bodies.Dquad.Enable(!isNemoh);
		bodies.C.Enable(!isNemoh);
		bodies.Cext.Enable(!isNemoh);
		bodies.Cadd.Enable(false);		// Not used for now
	};
	dropSolver.WhenAction();
	
	opInfinite.WhenAction = [&] {
		height.Enable(!bool(~opInfinite));
	};
	opInfinite.WhenAction();
}

void MainSolver::InitSerialize(bool ret) {
	if (!ret || IsNull(opIncludeBin)) 
		opIncludeBin = true;
	if (!ret || IsNull(numSplit)) 	
		numSplit = 2;	
	String manufacturer, productName, version, mbSerial;
	int numberOfProcessors;
	GetSystemInfo(manufacturer, productName, version, numberOfProcessors, mbSerial);
	if (!ret || IsNull(numThreads)) 	
		numThreads = numberOfProcessors;
	if (!ret || IsNull(~nemoh.xeff))
		nemoh.xeff <<= 0;
	if (!ret || IsNull(~nemoh.yeff))
		nemoh.yeff <<= 0;
	if (!ret || IsNull(~bodies.cx))
		bodies.cx <<= 0;
	if (!ret || IsNull(~bodies.cy))
		bodies.cy <<= 0;
	if (!ret || IsNull(~bodies.cz))
		bodies.cz <<= 0;
	if (!ret || IsNull(~Nf))
		Nf <<= 100;
	if (!ret || IsNull(~minF))
		minF <<= 0;
	if (!ret || IsNull(~maxF))
		maxF <<= 4;
	if (!ret || IsNull(~Nh))
		Nh <<= 1;
	if (!ret || IsNull(~minH))
		minH <<= 0;
	if (!ret || IsNull(~maxH))
		maxH <<= 0;
	if (!ret || IsNull(~height)) {
		height <<= Null;
		opInfinite <<= true;
	}
}
	
void MainSolver::Jsonize(JsonIO &json) {
	if (json.IsLoading()) {
		opIncludeBin = Null;
		numSplit = Null;	
		nemoh.xeff <<= Null;
		nemoh.yeff <<= Null;
		bodies.cx <<= Null;
		bodies.cy <<= Null;
		bodies.cz <<= Null;
		Nf <<= Null;
		minF <<= Null;
		maxF <<= Null;
		Nh <<= Null;
		minH <<= Null;
		maxH <<= Null;
	} else {
		dropSolverVal = dropSolver.GetData();
	}
	json
		("loadFrom", loadFrom)
		("saveTo", saveTo)
		("opIncludeBin", opIncludeBin)
		("numSplit", numSplit)
		("numThreads", numThreads)
		("opSplit", opSplit)
		("xeff", nemoh.xeff)
		("yeff", nemoh.yeff)
		("cx", bodies.cx)
		("cy", bodies.cy)
		("cz", bodies.cz)
		("Nf", Nf)
		("minF", minF)
		("maxF", maxF)
		("Nh", Nh)
		("minH", minH)
		("maxH", maxH)
		("dropSolver", dropSolverVal)
		("height", height)
		("opInfinite", opInfinite)
	;
	if (json.IsLoading()) {
		if (IsNull(dropSolverVal) || dropSolverVal < 0)
			dropSolverVal = 0;
	}
}

bool MainSolver::OnLoad(const BEM &bem) {
	String file = ~loadFrom;
	
	try {
		Load(file, bem);
		dropSolver.WhenAction();
	} catch (const Exc &e) {
		BEM::PrintError(DeQtfLf(e));
		return false;
	}
	return true;		
}		

void MainSolver::Load(const BEM &bem) {
	if (IsNull(~nemoh.g))
		nemoh.g <<= bem.g;
	if (IsNull(~nemoh.rho))
		nemoh.rho <<= bem.rho;	
	if (IsNull(~height)) {
		opInfinite <<= (bem.depth < 0);
		height.Enable(bem.depth > 0);
		height <<= (bem.depth > 0 ? bem.depth : Null);
	}
}
	
void MainSolver::Load(String file, const BEM &bem) {
	BEMCase data;
	
	data.Load(file, bem);
	bool isNemoh = data.IsNemoh();
	
	opInfinite <<= (data.h < 0);
	height.Enable(data.h > 0);
	height <<= (data.h > 0 ? data.h : Null);
	
	InitArray(true);
	
	for (int i = 0; i < data.bodies.size(); ++i) {
		const BEMBody &b = data.bodies[i];
		bodies.array.Add(b.meshFile, b.lidFile,  
					b.dof[BEM::SURGE], b.dof[BEM::SWAY], b.dof[BEM::HEAVE], 
					b.dof[BEM::ROLL], b.dof[BEM::PITCH], b.dof[BEM::YAW], 
					b.c0[0], b.c0[1], b.c0[2]);
	}
	bodies.array.SetCursor(0);
	dropSolver.WhenAction();
		
	Nf <<= data.Nf;
	
	minF <<= data.minF;
	maxF <<= data.maxF;	
	
	Nh <<= data.Nh;
	minH <<= data.minH;
	maxH <<= data.maxH;
	
	if (isNemoh) 
		dropSolver <<= BEMCase::NEMOHv115;
	else if (data.solver == BEMCase::HAMS)
		dropSolver <<= BEMCase::HAMS;
	
	const BEMBody &b = data.bodies[0];
	
	bodies.xcm = b.cg[0];
	bodies.ycm = b.cg[1];
	bodies.zcm = b.cg[2];
	MatrixXdToGridCtrl(bodies.M, b.M);
	MatrixXdToGridCtrl(bodies.Dlin, b.Dlin);
	MatrixXdToGridCtrl(bodies.Dquad, b.Dquad);
	MatrixXdToGridCtrl(bodies.C, b.C);
	MatrixXdToGridCtrl(bodies.Cext, b.Cext);
	MatrixXdToGridCtrl(bodies.Cadd, b.Cadd);
	
	nemoh.g <<= data.g;
	nemoh.rho <<= data.rho;
	
	nemoh.xeff <<= data.xeff;
	nemoh.yeff <<= data.yeff;
	
	nemoh.boxIrf <<= data.irf;
	nemoh.irfStep <<= data.irfStep;
	nemoh.irfDuration <<= data.irfDuration;
	
	nemoh.boxKochin <<= (data.nKochin > 0);
	nemoh.nKochin <<= data.nKochin;
	nemoh.minK <<= data.minK;
	nemoh.maxK <<= data.maxK;
	
	nemoh.showPressure <<= data.showPressure;
	
	nemoh.freeSurface <<= (data.nFreeX > 0);
	
	nemoh.freeX <<= (data.nFreeX > 0);
	nemoh.nFreeX <<= data.nFreeX;
	nemoh.domainX <<= data.domainX;
		
	nemoh.freeY <<= (data.nFreeY > 0);
	nemoh.nFreeY <<= data.nFreeY;
	nemoh.domainY <<= data.domainY;
}

void MainSolver::InitArray(bool isNemoh) {
	bodies.array.Reset();
	bodies.array.SetLineCy(EditField::GetStdHeight());
	bodies.array.AddColumn("Mesh file", 40);
	bodies.array.AddColumn("Lid file", 40);
	bodies.array.AddColumn("Surge", 10).With([](One<Ctrl>& x) {x.Create<Option>().SetReadOnly();});
	bodies.array.AddColumn("Sway",  10).With([](One<Ctrl>& x) {x.Create<Option>().SetReadOnly();});
	bodies.array.AddColumn("Heave", 10).With([](One<Ctrl>& x) {x.Create<Option>().SetReadOnly();});
	bodies.array.AddColumn("Roll",  10).With([](One<Ctrl>& x) {x.Create<Option>().SetReadOnly();});
	bodies.array.AddColumn("Pitch", 10).With([](One<Ctrl>& x) {x.Create<Option>().SetReadOnly();});
	bodies.array.AddColumn("Yaw",   10).With([](One<Ctrl>& x) {x.Create<Option>().SetReadOnly();});
	bodies.array.AddColumn("RotX",  10);
	bodies.array.AddColumn("RotY",  10);
	bodies.array.AddColumn("RotZ",  10);
	
	dropSolver.WhenAction();

	InitGrid(bodies.M, editMass);
	InitGrid(bodies.Dlin, editLinear);
	InitGrid(bodies.Dquad, editQuadratic);	
	InitGrid(bodies.C, editInternal);
	InitGrid(bodies.Cext, editExternal);
	InitGrid(bodies.Cadd, editAdd);
}

void MainSolver::InitGrid(GridCtrl &grid, EditDouble edit[]) {
	grid.Reset();
	grid.Absolute().Editing().Clipboard();
	for (int i = 0; i < 6; ++i)
		grid.AddColumn(InitCaps(BEM::StrDOF(i))).Edit(edit[i]);
	for (int y = 0; y < 6; ++y)
		for (int x = 0; x < 6; ++x)
			grid.Set(y, x, 0.);
}

void MainSolver::LoadMatrix(GridCtrl &grid, const Eigen::MatrixXd &mat) {
	for (int y = 0; y < 6; ++y)
		for (int x = 0; x < 6; ++x)
			grid.Set(y, x, mat(x, y));
}

bool MainSolver::Save(BEMCase &data, bool isNemoh) {
	if (!opInfinite)
		data.h = ~height;
	else
		data.h = -1;
	
	data.bodies.SetCount(bodies.array.GetCount());
	for (int i = 0; i < data.bodies.size(); ++i) {
		BEMBody &b = data.bodies[i];
		b.meshFile = bodies.array.Get(i, 0);
		b.lidFile = bodies.array.Get(i, 1);

		b.dof[BEM::SURGE] = bodies.array.Get(i, 2);
		b.dof[BEM::SWAY]  = bodies.array.Get(i, 3);
		b.dof[BEM::HEAVE] = bodies.array.Get(i, 4);
		b.dof[BEM::ROLL]  = bodies.array.Get(i, 5);
		b.dof[BEM::PITCH] = bodies.array.Get(i, 6);
		b.dof[BEM::YAW]   = bodies.array.Get(i, 7);
		b.c0[0] 	 		= bodies.array.Get(i, 8);
		b.c0[1] 	 		= bodies.array.Get(i, 9);
		b.c0[2] 	 		= bodies.array.Get(i, 10);
		
		if (i == 0) {
			b.cg[0] = bodies.xcm;
			b.cg[1] = bodies.ycm;
			b.cg[2] = bodies.zcm;
			GridCtrlToMatrixXd(b.M, bodies.M);
			GridCtrlToMatrixXd(b.Dlin, bodies.Dlin);
			GridCtrlToMatrixXd(b.Dquad, bodies.Dquad);
			GridCtrlToMatrixXd(b.C, bodies.C);
			GridCtrlToMatrixXd(b.Cext, bodies.Cext);
			GridCtrlToMatrixXd(b.Cadd, bodies.Cadd);
		}
	}
		
	data.Nf = ~Nf;
	data.minF = ~minF;
	data.maxF = ~maxF;	
	
	data.Nh = ~Nh;
	data.minH = ~minH;
	data.maxH = ~maxH;
	
	data.g = ~nemoh.g;
	data.rho = ~nemoh.rho;
	
	data.xeff = ~nemoh.xeff;
	data.yeff = ~nemoh.yeff;
	
	data.irf = ~nemoh.boxIrf;
	if (~nemoh.boxIrf) {
		data.irfStep = ~nemoh.irfStep;
		data.irfDuration = ~nemoh.irfDuration;
	} else
		data.irfStep = data.irfDuration = 0;
	
	if (~nemoh.boxKochin) {
		data.nKochin = ~nemoh.nKochin;
		data.minK = ~nemoh.minK;
		data.maxK = ~nemoh.maxK;
	} else {
		data.nKochin = 0;
		data.minK = data.maxK = 0;
	}
	
	data.showPressure = ~nemoh.showPressure;
	
	if (~nemoh.freeSurface) {
		data.nFreeX = ~nemoh.nFreeX;
		data.domainX = ~nemoh.domainX;
		
		data.nFreeY = ~nemoh.nFreeY;
		data.domainY = ~nemoh.domainY;	
	} else {
		data.nFreeX = data.nFreeY = 0;
		data.domainX = data.domainY = 0;
	}
	
	return true;
}

void MainSolver::arrayOnCursor() {
	int id = bodies.array.GetCursor();
	if (id < 0)
		return;
	
	bodies.meshFile <<= bodies.array.Get(id, 0);
	bodies.lidFile  <<= bodies.array.Get(id, 1);
	bodies.surge 	<<= bodies.array.Get(id, 2);
	bodies.sway 	<<= bodies.array.Get(id, 3);
	bodies.heave 	<<= bodies.array.Get(id, 4);
	bodies.roll 	<<= bodies.array.Get(id, 5);
	bodies.pitch 	<<= bodies.array.Get(id, 6);
	bodies.yaw 	 	<<= bodies.array.Get(id, 7);
	bodies.cx 		<<= bodies.array.Get(id, 8);
	bodies.cy 		<<= bodies.array.Get(id, 9);
	bodies.cz 		<<= bodies.array.Get(id, 10);
}

bool MainSolver::ArrayUpdateCursor() {
	try {
		bool isNemoh = ~dropSolver != BEMCase::HAMS;
		
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
		
		bodies.array.Set(id, 2, ~bodies.surge);
		bodies.array.Set(id, 3, ~bodies.sway);
		bodies.array.Set(id, 4, ~bodies.heave);
		bodies.array.Set(id, 5, ~bodies.roll);
		bodies.array.Set(id, 6, ~bodies.pitch);
		bodies.array.Set(id, 7, ~bodies.yaw);
		bodies.array.Set(id, 8, ~bodies.cx);
		bodies.array.Set(id, 9, ~bodies.cy);
		bodies.array.Set(id, 10,~bodies.cz);
	
		bodies.array.Update();
		
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
		return false;
	}
	
	return true;
}

void MainSolver::arrayClear() {
	bodies.meshFile <<= "";
	bodies.lidFile  <<= "";
	bodies.surge 	<<= true;
	bodies.sway 	<<= true;
	bodies.heave 	<<= true;
	bodies.roll 	<<= true;
	bodies.pitch 	<<= true;
	bodies.yaw 	 	<<= true;
	bodies.cx 		<<= 0;
	bodies.cy 		<<= 0;
	bodies.cz 		<<= 0;
}

void MainSolver::arrayOnAdd() {
	if (bodies.array.GetCount() == 0) {
		bool isNemoh = ~dropSolver != BEMCase::HAMS;
		InitArray(isNemoh);
	}
	bodies.array.Add();
	bodies.array.SetCursor(bodies.array.GetCount()-1);	
	arrayClear();
	ArrayUpdateCursor();
	dropSolver.WhenAction();
}

void MainSolver::arrayOnDuplicate() {
	if (bodies.array.GetCount() == 0) {
		BEM::PrintError(t_("No body available to duplicate"));
		return;
	}
	int id = bodies.array.GetCursor();
	if (id < 0) {
		BEM::PrintError(t_("Please select body to duplicate"));
		return;
	}
	int nr = id + 1;
	bodies.array.Insert(nr);
	for (int c = 0; c < bodies.array.GetColumnCount(); ++c)
		bodies.array.Set(nr, c, bodies.array.Get(id, c));
	bodies.array.Disable();
	bodies.array.SetCursor(nr);
	bodies.array.Enable();
	arrayOnCursor();
}

void MainSolver::arrayOnRemove() {
	if (bodies.array.GetCount() == 0) {
		BEM::PrintError(t_("No body available to remove"));
		return;
	}
	int id = bodies.array.GetCursor();
	if (id < 0) {
		BEM::PrintError(t_("Please select body to remove"));
		return;
	}
	bodies.array.Remove(id);
	if (id >= bodies.array.GetCount()) {
		id = bodies.array.GetCount()-1;
		if (id < 0) {
			arrayClear();
			ArrayUpdateCursor();
			return;
		}
	} 
	bodies.array.SetCursor(id);
}

bool MainSolver::OnSave(const BEM &bem) {
	try {
		String folder = ~saveTo;
		
		bool isNemoh = ~dropSolver != BEMCase::HAMS;
		
		BEMCase data;
		
		if (!Save(data, isNemoh))
			return false;
		
		UVector<String> errors = data.Check(~dropSolver);
		
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
		if (~opSplit) {
			if (IsNull(~numSplit)) {
				BEM::PrintError(t_("Please enter number of parts to split the simulation (min. is 2)"));
				return false;
			} else if (int(~numSplit) > data.Nf) {
				if (PromptOKCancel(Format(t_("Number of split cases %d must not be higher than number of frequencies %d"), int(~numSplit), data.Nf)
							   + S("&") + t_("Do you wish to fit the number of cases to the number of frequencies?"))) 
					numSplit <<= data.Nf;
				else
					return false;
			}
		}
		
		WaitCursor waitcursor;
		
		if (isNemoh) 
			data.SaveFolder(folder, ~opIncludeBin, ~opSplit ? int(~numSplit) : 1, Null		 , bem, ~dropSolver);
		else
			data.SaveFolder(folder, ~opIncludeBin, ~opSplit ? int(~numSplit) : 1, ~numThreads, bem, ~dropSolver);

	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
		return false;
	}
	
	return true;
}

void GridCtrlToMatrixXd(Eigen::MatrixXd &mat, const GridCtrl &grid) {
	mat.resize(grid.GetColumnCount(), grid.GetRowCount());
	
	for (int x = 0; x < mat.cols(); ++x)
		for (int y = 0; y < mat.cols(); ++y)
			mat(x, y) = grid.Get(y, x);
}

void MatrixXdToGridCtrl(GridCtrl &grid, const Eigen::MatrixXd &mat) {
	for (int x = 0; x < mat.cols(); ++x)
		for (int y = 0; y < mat.cols(); ++y)
			grid.Set(y, x, mat(x, y));
}
