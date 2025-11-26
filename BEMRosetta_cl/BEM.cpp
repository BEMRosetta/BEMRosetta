// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2024, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#include <ScatterDraw/DataSource.h>
#include <ScatterDraw/Equation.h>
#include "functions.h"
#include <SysInfo/Crash.h>
#include <MatIO/matio.h>

using namespace Upp;
using namespace Eigen;


Function <void(String)> BEM::Print 		  = [](String s) {Cout() << s;};
Function <void(String)> BEM::PrintWarning = [](String s) {Cout() << t_("Warning: ") << s;};
Function <void(String)> BEM::PrintError   = [](String s) {Cout() << t_("ERROR: ") << s;};

const char *BEM::strDOFtext[] 	 = {t_("surge"), t_("sway"), t_("heave"), t_("roll"), t_("pitch"), t_("yaw")};
const char *BEM::strDOFtextAbrev[] = {t_("s"), t_("w"), t_("h"), t_("r"), t_("p"), t_("y")};
const char *BEM::strDOFnum[] 	 = {t_("1"), t_("2"), t_("3"), t_("4"), t_("5"), t_("6")};
const char *BEM::strDOFnum_sub[] = {t_("₁"), t_("₂"), t_("₃"), t_("₄"), t_("₅"), t_("₆")};
const char *BEM::strDOFxyz[] 	 = {t_("x"), t_("y"), t_("z"), t_("rx"), t_("ry"), t_("rz")};
	
const char *BasicBEM::strDOFType[] = {t_("1,2,3,4,5,6"), t_("surge,sway,"), t_("x,y,z,rx,ry,rz"), ""};
BasicBEM::DOFType BEM::dofType = BasicBEM::DOFSurgeSway;

const char *BasicBEM::strHeadingType[] = {t_("-180->180º"), t_("0->360º"), ""};
BasicBEM::HeadingType BEM::headingType = BasicBEM::HEAD_180_180;
	

bool PrintStatus(String s, int v) {
	if (!IsNull(v) && v >= 0) {
		int done;
		if (v == 0)
			done = 0;
		else
			done = 20*v/100;
		int pending = 20 - done;
		printf("|%s%s|", ~String('*', done), ~String('-', pending));
	}
	const int totalChar = 80;
	if (!IsNull(s)) {
		String str = RemoveAccents(s);
		int num = totalChar - 1 - str.GetCount();
		printf(" %s%s", ~str, ~String(' ', num > 0 ? num : 0));
	} else
		printf("%s", ~String(' ', totalChar));
	printf("\r");
	return true;
};

	
BEM::BEM() {
	bemFilesAst = clone(bstFilesExt);	// clone(bemFilesExt);
	bemFilesAst.Replace(".", "*.");
	experimental = ToLower(GetExeTitle()).Find("experimental") >= 0;
}

int BEM::LoadBEM(String fileName, Function <bool(String, int)> Status, bool checkDuplicated) {
	Status(t_("Loading files"), 0);
	
	if (checkDuplicated) {
		for (int i = 0; i < hydros.size(); ++i) {
			if (hydros[i].dt.file == fileName) {
				BEM::Print(S("\n") + t_("Model is already loaded"));
				throw Exc(Format(t_("Model '%s' is already loaded"), fileName));
			}
		}
	}
	
	int num = Hydro::LoadHydro(hydros, fileName, Status);
	
	Bem().Nb = 0;
	for (int i = 0; i < hydros.size(); ++i) 
		Bem().Nb = max(Bem().Nb, hydros[i].dt.Nb);	
		
	UpdateHeadAll();
	UpdateHeadAllMD();
	
	return num;
}

Hydro &BEM::Join(UVector<int> &ids, Function <bool(String, int)> Status) {
	UVector<Hydro *>hydrosp;
	
	hydrosp.SetCount(ids.size());
	for (int i = 0; i < ids.size(); ++i) 
		hydrosp[i] = &hydros[ids[i]]; 
	
	Hydro hy;
	hy.Join(hydrosp);
	String error = hy.AfterLoad(Status);
	if (!error.IsEmpty()) 
		throw Exc(Format(t_("Problem joining models: '%s'\n%s"), error));	
	
	Sort(ids, StdLess<int>());
	for (int i = ids.size()-1; i >= 0; --i)
		hydros.Remove(ids[i]);
	hy.IncrementIdCount();
	hydros << hy;
	return Last(hydros);
}

void BEM::RemoveHydro(int id) {
	if (id < 0)
		return;
	hydros.Remove(id);
	if (hydros.IsEmpty())
		Hydro::ResetIdCount();
}

Hydro &BEM::Duplicate(int id) {
	Hydro &hy = hydros.Add();
	hy.Copy(hydros[id]);
	hy.IncrementIdCount();
	return hy;
}

Hydro &BEM::Average(UVector<int> &ids) {
	Hydro &hy = hydros.Add();
	hy.Average(hydros, ids);
	
	return hy;
}

void BEM::SymmetrizeForces(int id, bool xAxis) {
	hydros[id].Symmetrize_Forces(xAxis);
	hydros[id].Symmetrize_QTF(xAxis);
	hydros[id].Symmetrize_MD(xAxis);
	
	UpdateHeadAll();
	UpdateHeadAllMD();
}

void BEM::UpdateHeadAll() {
	headAll.Clear();
				
	for (int id = 0; id < hydros.size(); ++id) {
		for (int ih = 0; ih < hydros[id].dt.head.size(); ++ih) 
			FindAddDelta(headAll, FixHeading(hydros[id].dt.head[ih], headingType), 0.1);
	}
	Sort(headAll);
}

void BEM::UpdateHeadAllMD() {
	headAllMD.Clear();
				
	for (int id = 0; id < hydros.size(); ++id) {
		for (int ih = 0; ih < hydros[id].dt.mdhead.size(); ++ih) 
			FindAddDelta(headAllMD, FixHeading(hydros[id].dt.mdhead[ih], headingType), 0.1);
	}
	Sort(headAllMD, SortComplex);
}

void BEM::A0(int id) {
	hydros[id].GetA0();
}

void BEM::Kirf(int id, double maxT) {
	hydros[id].GetK_IRF(maxT, numValsA);
}

void BEM::Ainf(int id) {
	hydros[id].GetAinf();
} 

void BEM::Ainf_w(int id) {
	hydros[id].GetAinf_w();
}

void BEM::RAO(int id, double critDamp) {
	hydros[id].GetRAO(critDamp);
}

void BEM::BH(int id, int &num) {
	hydros[id].GetB_H(num);
} 

void BEM::Symmetrize(int id) {
	hydros[id].Symmetrize();
}

void BEM::OgilvieCompliance(int id, bool zremoval, bool thinremoval, bool decayingTail, UVector<int> &vidof, UVector<int> &vjdof) {
	hydros[id].GetOgilvieCompliance(zremoval, thinremoval, decayingTail, vidof, vjdof);
}

void BEM::ResetForces(int id, Hydro::FORCE force, bool forceMD, Hydro::FORCE forceQtf) {
	hydros[id].ResetForces(force, forceMD, forceQtf);
}

void BEM::MultiplyDOF(int id, double factor, const UVector<int> &idDOF, bool a, bool b, bool diag, bool f, bool md, bool qtf, bool C) {
	hydros[id].MultiplyDOF(factor, idDOF, a, b, diag, f, md, qtf, C);
}

void BEM::SwapDOF(int id, int ib1, int idof1, int ib2, int idof2) {
	if (idof1 < 0)
		hydros[id].SwapDOF(ib1, ib2);
	else	
		hydros[id].SwapDOF(ib1, idof1, ib2, idof2);
}

void BEM::FillFrequencyGapsABForces(int id, bool zero, int maxFreq) {
	hydros[id].FillFrequencyGapsABForces(zero, maxFreq);
}

void BEM::FillFrequencyGapsQTF(int id, bool zero, int maxFreq) {
	hydros[id].FillFrequencyGapsQTF(zero, maxFreq);
}

void BEM::FillFrequencyGapsABForcesZero(int id) {
	hydros[id].FillFrequencyGapsABForcesZero();
}

void BEM::FillFrequencyGapsQTFZero(int id) {
	hydros[id].FillFrequencyGapsQTFZero();
}

void BEM::CopyQTF_MD(int id) {
	hydros[id].CopyQTF_MD();
}

void BEM::DeleteBodies(int id, const UVector<int> &idBody) {
	hydros[id].DeleteBodies(idBody);
}

void BEM::DeleteHeadingsFrequencies(int id, const UVector<int> &idFreq, const UVector<int> &idFreqQTF, const UVector<int> &idHead, const UVector<int> &idHeadMD, const UVector<int> &idHeadQTF) {
	hydros[id].DeleteFrequencies(idFreq);
	hydros[id].DeleteFrequenciesQTF(idFreqQTF);
	hydros[id].DeleteHeadings(idHead);
	hydros[id].DeleteHeadingsMD(idHeadMD);
	hydros[id].DeleteHeadingsQTF(idHeadQTF);
}

void BEM::TranslationTo(int id, const MatrixXd &to, Function <bool(String, int pos)> Status) {
	hydros[id].GetTranslationTo(to, zeroIfEmpty, Status);
}

void BEM::WaveTo(int id, double xto, double yto) {
	hydros[id].GetWaveTo(xto, yto, g);
}

String BEM::SpreadNegative(int id, Function <bool(String, int)> Status) {
	return hydros[id].SpreadNegative(Status);
}

void BEM::MapMeshes(int idh, int ib, const UVector<int> &idms, bool oneCase) {
	hydros[idh].MapMeshes(hydros, ib, idms, oneCase);	
	
	Bem().Nb = 0;
	for (int i = 0; i < hydros.size(); ++i) 
		Bem().Nb = max(Bem().Nb, hydros[i].dt.Nb);	
}

int BEM::LoadBody(String fileName, Function <bool(String, int pos)> Status, bool cleanPanels, bool checkDuplicated, const UVector<int> &idxs) {
	Status(Format(t_("Loading mesh '%s'"), fileName), 10);
	
	if (checkDuplicated) {
		for (int i = 0; i < surfs.size(); ++i) {
			if (surfs[i].dt.fileName == fileName) {
				BEM::Print(S("\n") + t_("Model is already loaded"));
				throw Exc(t_("Model is already loaded"));
			}
		}
	}
	UArray<Body> meshes;
	String error = Body::Load(meshes, fileName, rho, g, cleanPanels, roundVal, roundEps, idxs);
	if (!error.IsEmpty()) {
		BEM::Print("\n" + Format(t_("Problem loading '%s'") + S("\n%s"), fileName, error));
		throw Exc(Format(t_("Problem loading '%s'") + S("\n%s"), fileName, error));
	}
	int num = meshes.size();
	for (int i = 0; i < num; ++i)
		surfs.Add(pick(meshes[i]));
	
	Status(Format(t_("Mesh '%s' loaded"), fileName), 100);
	return num;
}

void BEM::SaveBody(String fileName, const UVector<int> &ids, Body::MESH_FMT type, Body::MESH_TYPE meshType, bool symX, bool symY) {
	if (type == Body::UNKNOWN) {
		String ext = ToLower(GetFileExt(fileName));
		
		if (ext == ".gdf")
			type = Body::WAMIT_GDF;
		else if (ext == ".dat")
			type = Body::NEMOH_DAT;
		else if (ext == ".")
			type = Body::NEMOH_PRE;
		else if (ext == ".pnl")
			type = Body::HAMS_PNL;
		else if (ext == ".stl")
			type = Body::STL_TXT;
		else if (ext == ".mesh")
			type = Body::BEM_MESH;
		else
			throw Exc(Format(t_("Conversion to file type '%s' not supported"), fileName));
	}
	
	UArray<Body> bds;
	for (int i = 0; i < ids.size(); ++i)
		bds.Add(clone(surfs[ids[i]]));
	
	UVector<String> fileNames;
	if (ids.size() == 1)
		fileNames << fileName;
	else {
		for (int ib = 0; ib < ids.size(); ++ib)
			fileNames << Format("%s_%d", fileName, ib+1);
	}
	Body::SaveAs(bds, fileNames, type, meshType, rho, g, symX, symY);
}

void BEM::HealingBody(int id, bool basic, Function <bool(String, int)> Status) {
	Status(Format(t_("Healing mesh '%s'"), surfs[id].dt.fileName), 10);
	Print(S("\n\n") + Format(t_("Healing mesh '%s'"), surfs[id].dt.fileName));
	
	String ret;
	try {
		ret = surfs[id].Heal(basic, rho, g, roundVal, roundEps, Status);
	} catch (Exc e) {
		surfs.SetCount(surfs.size()-1);
		Print("\n" + Format(t_("Problem healing '%s': %s") + S("\n%s"), e));
		throw std::move(e);
	}
	if (!ret.IsEmpty()) {
		ret.Replace("\n", "\n- ");
		Print(ret);
	} else
		Print(S(". ") + t_("The mesh is in good condition"));
}

void BEM::OrientSurface(int id, Function <bool(String, int)> Status) {
	Status(Format(t_("Orienting surface mesh '%s'"), surfs[id].dt.fileName), 10);
	Print(S("\n\n") + Format(t_("Orienting surface mesh '%s'"), surfs[id].dt.fileName));
	
	try {
		surfs[id].Orient();
	} catch (Exc e) {
		surfs.SetCount(surfs.size()-1);
		Print("\n" + Format(t_("Problem orienting surface '%s': %s") + S("\n%s"), e));
		throw std::move(e);
	}
}

void BEM::UnderwaterBody(int id, Function <bool(String, int pos)> Status) {
	Status(Format(t_("Getting underwater mesh '%s'"), surfs[id].dt.fileName), 10);
	
	Body &mesh = surfs.Add();
	Body &orig = surfs[id];
	mesh.dt.fileName = orig.dt.fileName;
	
	try {
		mesh.dt.mesh.CutZ(orig.dt.mesh, -1);
	} catch (Exc e) {
		surfs.SetCount(surfs.size()-1);
		Print("\n" + Format(t_("Problem loading '%s': %s") + S("\n%s"), e));
		throw std::move(e);
	}
}

void BEM::RemoveBody(int id) {
	Body *b = &(surfs[id]);
	for (int i = 0; i < surfs.size(); ++i) {
		if (i != id) {
			for (int j = 0; j < surfs[i].cdt.damagedBodies.size(); ++j) {
				if (surfs[i].cdt.damagedBodies[j] == b) {
					surfs[i].cdt.damagedBodies.Remove(j);
					break;
				}
			}
		}
	}
	if (surfs[id].dt.GetCode() == Body::MOORING_MESH)
		fast.Clear();
	surfs.Remove(id);
	if (surfs.IsEmpty()) {
		Body::ResetIdCount();
		fast.Clear();
	}
}

void BEM::DuplicateBody(int id) {
	Body &surf = surfs.Add();
	surf.dt.SetCode(Body::EDIT);
	
	try {
		surf = clone(surfs[id]);
	} catch (Exc e) {
		surfs.SetCount(surfs.size()-1);
		Print("\n" + Format(t_("Problem duplicating '%s': %s") + S("\n%s"), e));
		throw std::move(e);
	}
}

void BEM::JoinBody(const UVector<int> &idsJoin) {
	Body &surf = surfs.Add();
	surf.dt.SetCode(Body::EDIT);
	
	try {
		surf = clone(surfs[First(idsJoin)]);
		for (int i = 1; i < idsJoin.size(); ++i) {
			surf.Append(surfs[idsJoin[i]].dt.mesh, rho, g);	
			surf.dt.fileName << "/" << surfs[idsJoin[i]].dt.fileName;
		}
	} catch (Exc e) {
		surfs.SetCount(surfs.size()-1);
		Print("\n" + Format(t_("Problem joining '%s': %s") + S("\n%s"), e));
		throw std::move(e);
	}
}

UVector<int> BEM::SplitBody(int id, Function <bool(String, int pos)> Status) {
	Status(Format(t_("Splitting mesh '%s'"), surfs[id].dt.fileName), 0);
	Body &orig = surfs[id];
	
	UVector<int> ret;
	try {
		UVector<UVector<int>> sets = orig.dt.mesh.GetPanelSets(Status);
		if (sets.size() == 1)
			return ret;
		for (int i = 0; i < sets.size(); ++i) {		
			Body &surf = surfs.Add();
			ret << surfs.size()-1-1;		// One more as id is later removed
			for (int ii = 0; ii < sets[i].size(); ++ii) 
				surf.dt.mesh.panels << clone(orig.dt.mesh.panels[sets[i][ii]]);	
			
			surf.dt.mesh.nodes = clone(orig.dt.mesh.nodes);
			surf.dt.SetCode(orig.dt.GetCode());
			surf.AfterLoad(rho, g, false, true);
		}
		RemoveBody(id);
	} catch (Exc e) {
		surfs.SetCount(surfs.size() - ret.size());
		Print("\n" + Format(t_("Problem loading '%s': %s") + S("\n%s"), e));
		throw std::move(e);
	}
	return ret;
}

void BEM::AddFlatRectangle(double x, double y, double z, double size, double panWidth, double panHeight) {
	try {
		Body &surf = surfs.Add();

		surf.dt.SetCode(Body::EDIT);
		surf.dt.mesh.AddFlatRectangle(panWidth, panHeight, size, size); 
		surf.dt.mesh.Translate(x, y, z);
		surf.dt.c0 = Point3D(0, 0, 0);
	} catch (Exc e) {
		surfs.SetCount(surfs.size() - 1);
		Print("\n" + Format(t_("Problem adding flat panel: %s"), e));
		throw std::move(e);
	}	
}

void BEM::AddRevolution(double x, double y, double z, double size, UVector<Pointf> &vals, double angle, bool close, Function <bool(String)> Prompt) {
	try {
		Body &surf = surfs.Add();

		surf.dt.SetCode(Body::EDIT);
		surf.dt.mesh.AddRevolution(vals, size, angle, close, Prompt); 
		surf.dt.mesh.Translate(x, y, z);
		surf.dt.c0 = Point3D(0, 0, 0);
	} catch (Exc e) {
		surfs.SetCount(surfs.size() - 1);
		Print("\n" + Format(t_("Problem adding revolution surface: %s"), e));
		throw std::move(e);
	}	
}

void BEM::Extract(int id, int cutXYZ, int cutPosNeg) {
	try {
		Body &surf = surfs.Add();
		Body &orig = surfs[id];
		
		surf.dt.SetCode(Body::EDIT);
		if (cutXYZ == 0)
			surf.dt.mesh.CutX(orig.dt.mesh, cutPosNeg == 0 ? 1 : -1);
		else if (cutXYZ == 1)
			surf.dt.mesh.CutY(orig.dt.mesh, cutPosNeg == 0 ? 1 : -1);
		else
			surf.dt.mesh.CutZ(orig.dt.mesh, cutPosNeg == 0 ? 1 : -1);
		surf.dt.c0 = clone(orig.dt.c0);
	} catch (Exc e) {
		surfs.SetCount(surfs.size() - 1);
		Print("\n" + Format(t_("Problem adding flat panel: %s"), e));
		throw std::move(e);
	}	
}	
	
void BEM::AddPanels(const Body &surfFrom, UVector<int> &panelList) {
	try {
		Body &surf = surfs.Add();

		surf.dt.SetCode(Body::EDIT);
		surf.dt.mesh.AddPanels(surfFrom.dt.mesh, panelList); 
		surf.dt.c0 = clone(surfFrom.dt.c0);
	} catch (Exc e) {
		surfs.SetCount(surfs.size() - 1);
		Print("\n" + Format(t_("Problem adding surface from panels: %s"), e));
		throw std::move(e);
	}	
}

void BEM::AddPolygonalPanel(double x, double y, double z, double size, UVector<Pointf> &vals, bool quads) {
	try {
		Body &surf = surfs.Add();

		surf.dt.SetCode(Body::EDIT);
		surf.dt.mesh.AddPolygonalPanel(vals, size, true, quads); 
		surf.dt.mesh.Translate(x, y, z);
		surf.dt.c0 = Point3D(0, 0, 0);
	} catch (Exc e) {
		surfs.SetCount(surfs.size() - 1);
		Print("\n" + Format(t_("Problem adding polygonal panel: %s"), e));
		throw std::move(e);
	}	
}

void BEM::Extrude(int id, double dx, double dy, double dz, bool close) {
	try {
		Body &surf = surfs[id];
	
		if (surf.dt.mesh.volumex > 0.001 && surf.dt.mesh.volumey > 0.001 && surf.dt.mesh.volumez > 0.001)
			throw Exc(t_("It is only possible to extrude a flat surface"));

		surf.dt.mesh.Extrude(dx, dy, dz, close);
		
	} catch (Exc e) {
		Print("\n" + Format(t_("Problem extruding surface: %s"), e));
		throw std::move(e);
	}	
}

void BEM::AddWaterSurface(int id, char c, double meshRatio, bool quads) {
	try {
		Body &surf = surfs.Add();

		surf.dt.SetCode(Body::EDIT);
		surf.dt.mesh.AddWaterSurface(surfs[id].dt.mesh, surfs[id].dt.under, c, roundVal, roundEps, meshRatio, quads); 
		surf.dt.c0 = surfs[id].dt.c0;
		
		if (c == 'r')
			surf.dt.name = t_("Water surface removed");
		else if (c == 'f')
			surf.dt.name = t_("Water surface");
		else if (c == 'e')
			surf.dt.name = t_("Water surface extracted");
		surf.dt.name = surfs[id].dt.name + " " + surf.dt.name;
		surf.dt.fileName =  "";
		
		if (true) {		// Clean panels
			Surface::RemoveDuplicatedPanels(surf.dt.mesh.panels);
			Surface::RemoveDuplicatedPointsAndRenumber(surf.dt.mesh.panels, surf.dt.mesh.nodes);
			Surface::RemoveDuplicatedPanels(surf.dt.mesh.panels);
		}
		surf.AfterLoad(rho, g, false, true);
	} catch (Exc e) {
		surfs.SetCount(surfs.size() - 1);
		Print("\n" + Format(t_("Problem adding water surface: %s"), e));
		throw std::move(e);
	}	
}
			
String BEM::LoadSerializeJson() {
	csvSeparator = Null;
	volWarning = volError = Null;
	roundVal = roundEps = Null;
	
	String ret;
	String folder = GetBEMRosettaDataFolder();
	if (!DirectoryCreateX(folder))
		ret = Format(t_("Impossible to create folder '%s' to store configuration file"), folder);
	else {
		String fileName = AFX(folder, "configdata.cf");
		if (!FileExists(fileName)) 
			ret = t_("First time");
		else {
			String jsonText = LoadFile(fileName);
			if (jsonText.IsEmpty())
				ret = Format(t_("Configuration file '%s' is empty"), fileName);
			else {
				ret = LoadFromJsonError(*this, jsonText);
				if (!ret.IsEmpty())
					ret = Format(t_("Problem loading configuration file '%s': %s"), fileName, ret);
			}
			if (!ret.IsEmpty()) {
				DirectoryCreateX(AFX(folder, "Errors"));
				FileCopy(fileName, AFX(folder, "Errors", "configdata.cf"));
			}
		}
	}
	
	bool ok = ret.IsEmpty();
	if (!ok || IsNull(csvSeparator))
		csvSeparator = ";";
	if (!ok || IsNull(volWarning))
		volWarning = 1;
	if (!ok || IsNull(volError))
		volError = 10;
	if (!ok || IsNull(roundVal))
		roundVal = 1;
	if (!ok || IsNull(roundEps))
		roundEps = 1E-8;
	if (!ok || IsNull(g)) 
		g = 9.80665;
	if (!ok || IsNull(depth)) 
		depth = 100;
	if (!ok || IsNull(rho)) 
		rho = 1025;
	if (!ok || IsNull(len)) 
		len = 1;
	//if (!ok || IsNull(discardNegDOF))
	//	discardNegDOF = false;
	//if (!ok || IsNull(thres)) 
	//	thres = 0.01;
	if (!ok || IsNull(calcAinf))
		calcAinf = true;
	if (!ok || IsNull(calcAinf_w))
		calcAinf_w = true;
	if (!ok || IsNull(maxTimeA))
		maxTimeA = 120;
	if (!ok || IsNull(numValsA))
		numValsA = 1000;
	if (!ok || IsNull(onlyDiagonal))
		onlyDiagonal = false;
	if (!ok || IsNull(volWarning))	
		volWarning = 1;
	if (!ok || IsNull(volError))
		volError = 10;
	if (!ok || IsNull(roundVal))
		roundVal = 1;
	if (!ok || IsNull(roundEps))
		roundEps = 1E-8;
	if (!ok || IsNull(legend_w_solver))
		legend_w_solver = true;
	if (!ok || IsNull(legend_w_units))
		legend_w_units = true;
	if (!ok || IsNull(zeroIfEmpty))
		zeroIfEmpty = true;
	if (!ok || IsNull(opT))
		opT = 0;
						
	return ret;
}

bool BEM::ClearTempFiles() {
	String folder = GetTempFilesFolder();
	DeleteFolderDeepWildcardsX(folder, "*.*");	Sleep(100);
	return DirectoryCreateX(folder);
}
	
bool BEM::StoreSerializeJson() {
	String folder = AFX(GetAppDataFolder(), "BEMRosetta");
	if (!DirectoryCreateX(folder))
		return 0;
	String fileName = AFX(folder, "configdata.cf");
	return StoreAsJsonFile(*this, fileName, false);
}

String FormatDoubleEmpty(double val) {
	if (IsNull(val))
		return t_("'empty'");
	else
		return FDS(val, 10, false);
}

String FormatIntEmpty(int val) {
	if (IsNull(val))
		return t_("'empty'");
	else
		return FormatInt(val);
}
