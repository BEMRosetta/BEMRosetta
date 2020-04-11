#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas.h>
#include <RasterPlayer/RasterPlayer.h>

#include <BEMRosetta/BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;

#define IMAGECLASS ImgNemoh
#define IMAGEFILE <BEMRosetta/BEMRosetta/main.iml>
#include <Draw/iml.h>

#include "main.h"


void MainNemoh::Init(const BEMData &bem) {
	CtrlLayout(*this);
	
	loadFrom.WhenChange = THISBACK(OnLoad);
	loadFrom.BrowseRightWidth(40).UseOpenFolder().BrowseOpenFolderWidth(10).UseDropping();
	butLoad.WhenAction = [&] {loadFrom.DoGo();};

	//saveTo.WhenChange  = [&] {return OnSave(bem);};
	saveTo.BrowseRightWidth(40).UseOpenFolder().BrowseOpenFolderWidth(10).UseDropping();
	butSave.WhenAction = [&] {OnSave(bem);};
	
	array.WhenSel = THISBACK(arrayOnCursor);
	
	butAdd <<= THISBACK(arrayOnAdd);
	butDuplicate <<= THISBACK(arrayOnDuplicate);
	butRemove <<= THISBACK(arrayOnRemove);
	
	meshFile.WhenAction = THISBACK(arrayUpdateCursor);
	meshFile.WhenChange = [&] {arrayUpdateCursor(); return true;};
	surge.WhenAction = THISBACK(arrayUpdateCursor);
	sway.WhenAction = THISBACK(arrayUpdateCursor);
	heave.WhenAction = THISBACK(arrayUpdateCursor);
	roll.WhenAction = THISBACK(arrayUpdateCursor);
	pitch.WhenAction = THISBACK(arrayUpdateCursor);
	yaw.WhenAction = THISBACK(arrayUpdateCursor);
	cx.WhenAction = THISBACK(arrayUpdateCursor);
	cy.WhenAction = THISBACK(arrayUpdateCursor);
	cz.WhenAction = THISBACK(arrayUpdateCursor);
	
	freeSurface.Transparent(false);
	freeSurface.WhenAction = [&] {freeX.Enable(~freeSurface);	freeY.Enable(~freeSurface);};
	freeSurface.WhenAction();
	
	opSplit.WhenAction = [&] {numSplit.Enable(~opSplit);		labelSplit.Enable(~opSplit);};
	opSplit.WhenAction();
	
	dropSolver.Add(0, t_("Nemoh 115+"));
	dropSolver.Add(1, t_("Nemoh"));
	dropSolver.SetIndex(0);
	
	ecnLogo.Set(ImgNemoh::ECN());
	ecnLogo.SetHyperlink("https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-presentation-192863.kjsp");
}

void MainNemoh::InitSerialize(bool ret) {
	if (!ret || IsNull(opIncludeBin)) 
		opIncludeBin = true;
	if (!ret || IsNull(numSplit)) 	
		numSplit = 1;	
	if (!ret || IsNull(~xeff))
		xeff <<= 0;
	if (!ret || IsNull(~yeff))
		yeff <<= 0;
	if (!ret || IsNull(~cx))
		cx <<= 0;
	if (!ret || IsNull(~cy))
		cy <<= 0;
	if (!ret || IsNull(~cz))
		cz <<= 0;
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
}

void MainNemoh::Load(const BEMData &bem) {
	if (IsNull(~g))
		g <<= bem.g;
	if (IsNull(~rho))
		rho <<= bem.rho;	
	if (IsNull(~h))
		h <<= bem.depth;
}

void MainNemoh::Jsonize(JsonIO &json) {
	json
		("loadFrom", loadFrom)
		("saveTo", saveTo)
		("opIncludeBin", opIncludeBin)
		("numSplit", numSplit)
		("opSplit", opSplit)
		("xeff", xeff)
		("yeff", yeff)
		("cx", cx)
		("cy", cy)
		("cz", cz)
		("Nf", Nf)
		("minF", minF)
		("maxF", maxF)
		("Nh", Nf)
		("minH", minH)
		("maxH", maxH)
	;
}

bool MainNemoh::OnLoad() {
	String file = ~loadFrom;
	
	try {
		NemohCal data;
		
		if (!data.Load(file)) {
			Exclamation(Format("Problem loading %s file", DeQtf(file)));
			return false;
		}
		Load(data);	
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
		return false;
	}
	
	return true;
}

void MainNemoh::InitArray() {
	array.Reset();
	array.SetLineCy(EditField::GetStdHeight());
	array.AddColumn("File", 40);
	array.AddColumn("#points", 10);
	array.AddColumn("#panels", 10);
	array.AddColumn("Surge", 10).With([](One<Ctrl>& x) {x.Create<Option>().SetReadOnly();});
	array.AddColumn("Sway",  10).With([](One<Ctrl>& x) {x.Create<Option>().SetReadOnly();});
	array.AddColumn("Heave", 10).With([](One<Ctrl>& x) {x.Create<Option>().SetReadOnly();});
	array.AddColumn("Roll",  10).With([](One<Ctrl>& x) {x.Create<Option>().SetReadOnly();});
	array.AddColumn("Pitch", 10).With([](One<Ctrl>& x) {x.Create<Option>().SetReadOnly();});
	array.AddColumn("Yaw",   10).With([](One<Ctrl>& x) {x.Create<Option>().SetReadOnly();});
	array.AddColumn("RotX",  10);
	array.AddColumn("RotY",  10);
	array.AddColumn("RotZ",  10);
}

void MainNemoh::Load(const NemohCal &data) {
	g <<= data.g;
	rho <<= data.rho;
	h <<= data.h;
	xeff <<= data.xeff;
	yeff <<= data.yeff;
	
	InitArray();
	
	for (int i = 0; i < data.bodies.GetCount(); ++i) {
		const NemohBody &b = data.bodies[i];
		array.Add(b.meshFile, b.npoints, b.npanels, b.surge, b.sway, b.heave, b.roll, b.pitch, b.yaw, b.cx, b.cy, b.cz);
	}
	array.SetCursor(0);
		
	Nf <<= data.Nf;
	
	minF <<= data.minF;
	maxF <<= data.maxF;	
	
	Nh <<= data.Nh;
	minH <<= data.minH;
	maxH <<= data.maxH;
	
	irf <<= data.irf;
	irfStep <<= data.irfStep;
	irfDuration <<= data.irfDuration;
	
	kochin <<= data.nKochin > 0;
	nKochin <<= data.nKochin;
	minK <<= data.minK;
	maxK <<= data.maxK;
	
	showPressure <<= data.showPressure;
	
	freeSurface <<= data.nFreeX > 0;
	freeX <<= data.nFreeX > 0;
	nFreeX <<= data.nFreeX;
	domainX <<= data.domainX;
		
	freeY <<= data.nFreeY > 0;
	nFreeY <<= data.nFreeY;
	domainY <<= data.domainY;
}

void MainNemoh::Save(NemohCal &data) {
	data.g = ~g;
	data.rho = ~rho;
	data.h = ~h;
	data.xeff = ~xeff;
	data.yeff = ~yeff;
	
	data.bodies.SetCount(array.GetCount());
	for (int i = 0; i < data.bodies.GetCount(); ++i) {
		NemohBody &b = data.bodies[i];
		b.meshFile = array.Get(i, 0);
		b.npoints = array.Get(i, 1);
		b.npanels = array.Get(i, 2);
		if (FileExists(b.meshFile)) {
			MeshData dat;
			bool x0z;
			if (dat.LoadDatNemoh(b.meshFile, x0z).IsEmpty()) {
				b.npoints = dat.mesh.GetNumNodes();
				b.npanels = dat.mesh.GetNumPanels();
			}
		}		
		b.surge = array.Get(i, 3);
		b.sway = array.Get(i, 4);
		b.heave = array.Get(i, 5);
		b.roll = array.Get(i, 6);
		b.pitch = array.Get(i, 7);
		b.yaw = array.Get(i, 8);
		b.cx = array.Get(i, 9);
		b.cy = array.Get(i, 10);
		b.cz = array.Get(i, 11);
	}
		
	data.Nf = ~Nf;
	
	data.minF = ~minF;
	data.maxF = ~maxF;	
	
	data.Nh = ~Nh;
	data.minH = ~minH;
	data.maxH = ~maxH;
	
	data.irf = ~irf;
	if (~irf) {
		data.irfStep = ~irfStep;
		data.irfDuration = ~irfDuration;
	} else
		data.irfStep = data.irfDuration = 0;
	
	if (~kochin) {
		data.nKochin = ~nKochin;
		data.minK = ~minK;
		data.maxK = ~maxK;
	} else {
		data.nKochin = 0;
		data.minK = data.maxK = 0;
	}
	
	data.showPressure = ~showPressure;
	
	if (~freeSurface) {
		data.nFreeX = ~nFreeX;
		data.domainX = ~domainX;
		
		data.nFreeY = ~nFreeY;
		data.domainY = ~domainY;	
	} else {
		data.nFreeX = data.nFreeY = 0;
		data.domainX = data.domainY = 0;
	}
}

void MainNemoh::arrayOnCursor() {
	int id = array.GetCursor();
	if (id < 0)
		return;
	
	meshFile <<= array.Get(id, 0);
	surge <<= array.Get(id, 3);
	sway <<= array.Get(id, 4);
	heave <<= array.Get(id, 5);
	roll <<= array.Get(id, 6);
	pitch <<= array.Get(id, 7);
	yaw <<= array.Get(id, 8);
	cx <<= array.Get(id, 9);
	cy <<= array.Get(id, 10);
	cz <<= array.Get(id, 11);
}

void MainNemoh::arrayUpdateCursor() {
	int id = array.GetCursor();
	if (id < 0) {
		if (array.GetCount() == 0) {
			InitArray();
			array.Add();
			id = 0;
			cx <<= 0;
			cy <<= 0;
			cz <<= 0;
		} else
			id = array.GetCount()-1;
	}	
	array.Set(id, 0, ~meshFile);
	array.Set(id, 3, ~surge);
	array.Set(id, 4, ~sway);
	array.Set(id, 5, ~heave);
	array.Set(id, 6, ~roll);
	array.Set(id, 7, ~pitch);
	array.Set(id, 8, ~yaw);
	array.Set(id, 9, ~cx);
	array.Set(id, 10, ~cy);
	array.Set(id, 11, ~cz);
}

void MainNemoh::arrayClear() {
	meshFile <<= "";
	surge <<= false;
	sway <<= false;
	heave <<= false;
	roll <<= false;
	pitch <<= false;
	yaw <<= false;
	cx <<= 0;
	cy <<= 0;
	cz <<= 0;
}

void MainNemoh::arrayOnAdd() {
	if (array.GetCount() == 0)
		InitArray();
	array.Add();
	array.SetCursor(array.GetCount()-1);	
	arrayClear();
	arrayUpdateCursor();
}

void MainNemoh::arrayOnDuplicate() {
	if (array.GetCount() == 0) {
		Exclamation(t_("No body available to duplicate"));
		return;
	}
	int id = array.GetCursor();
	if (id < 0) {
		Exclamation(t_("Please select body to duplicate"));
		return;
	}
	int nr = id + 1;
	array.Insert(nr);
	for (int c = 0; c < array.GetColumnCount(); ++c)
		array.Set(nr, c, array.Get(id, c));
	array.Disable();
	array.SetCursor(nr);
	array.Enable();
	arrayOnCursor();
}

void MainNemoh::arrayOnRemove() {
	if (array.GetCount() == 0) {
		Exclamation(t_("No body available to remove"));
		return;
	}
	int id = array.GetCursor();
	if (id < 0) {
		Exclamation(t_("Please select body to remove"));
		return;
	}
	array.Remove(id);
	if (id >= array.GetCount()) {
		id = array.GetCount()-1;
		if (id < 0) {
			arrayClear();
			arrayUpdateCursor();
			return;
		}
	} 
	array.SetCursor(id);
}

bool MainNemoh::OnSave(const BEMData &bem) {
	try {
		String nemohFolder = ~saveTo;
		
		NemohCal data;
		Save(data);	
		Vector<String> res = data.Check();
		if (!res.IsEmpty()) {
			String str;
			for (int i = 0; i < res.GetCount(); ++i)
			 	str << "\n- " << res[i];
			if (!ErrorOKCancel(Format(t_("Errors found in Nemoh data:%s&Do you wish to continue?"), DeQtfLf(str))))
				return false;
		}
		if (!DirectoryExists(nemohFolder)) {
			if (!PromptYesNo(Format(t_("Folder %s does not exist.&Do you wish to create it?"), DeQtfLf(nemohFolder))))
				return false;
			RealizeDirectory(nemohFolder);
		} else {
			if (!PromptYesNo(Format(t_("Folder %s contents will be deleted.&Do you wish to continue?"), DeQtfLf(nemohFolder))))
				return false;
		}
		if (~opSplit) {
			if (IsNull(~numSplit)) {
				Exclamation(t_("Please enter number of parts to split the simulation (min. is 1)"));
				return false;
			} else if (int(~numSplit) > data.Nf) {
				if (PromptOKCancel(Format(t_("Number of Nemoh cases %d must not be higher than number of frequencies %d"), int(~numSplit), data.Nf)
							   + S("&") + t_("Do you wish to fit the number of cases to the number of frequencies?"))) 
					numSplit <<= data.Nf;
				else
					return false;
			}
		}
		data.SaveFolder(nemohFolder, ~opIncludeBin, ~opSplit ? int(~numSplit) : 1, bem, dropSolver.GetData());
	} catch (Exc e) {
		Exclamation(DeQtfLf(e));
		return false;
	}
	
	return true;
}

void MainNemoh::DragAndDrop(Point , PasteClip& d) {
	if (IsDragAndDropSource())
		return;
	if (AcceptFiles(d)) {
		Vector<String> files = GetFiles(d);
		for (int i = 0; i < files.GetCount(); ++i) {
			loadFrom <<= files[i];
			OnLoad();
			break;
		}
	}
}

bool MainNemoh::Key(dword key, int ) {
	if (key == K_CTRL_V) {
		Vector<String> files = GetFiles(Ctrl::Clipboard());
		for (int i = 0; i < files.GetCount(); ++i) {
			loadFrom <<= files[i];
			OnLoad();
			break;
		}
		return true;
	}
	return false;
}