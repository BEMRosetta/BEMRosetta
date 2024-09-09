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
Function <void(String)> BEM::PrintWarning = [](String s) {Cout() << s;};
Function <void(String)> BEM::PrintError   = [](String s) {Cout() << s;};

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
	if (v == 0)
		printf("\n%s", ~RemoveAccents(s));
	else {
		printf("\r%s", ~RemoveAccents(s));
		int done = v/5;
		int pending = 20 - done;
		printf(" |%s%s|                    ", ~String('*', done), ~String('-', pending));
	}
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

void BEM::RAO(int id) {
	hydros[id].GetRAO();
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

void BEM::MultiplyDOF(int id, double factor, const UVector<int> &idDOF, bool a, bool b, bool diag, bool f, bool md, bool qtf) {
	hydros[id].MultiplyDOF(factor, idDOF, a, b, diag, f, md, qtf);
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

int BEM::LoadBody(String fileName, Function <bool(String, int pos)> Status, bool cleanPanels, bool checkDuplicated) {
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
	String error = Body::Load(meshes, fileName, rho, g, cleanPanels, roundVal, roundEps);
	if (!error.IsEmpty()) {
		BEM::Print("\n" + Format(t_("Problem loading '%s'") + S("\n%s"), fileName, error));
		throw Exc(Format(t_("Problem loading '%s'") + S("\n%s"), fileName, error));
	}
	int num = meshes.size();
	for (int i = 0; i < num; ++i)
		surfs.Add(pick(meshes[i]));
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
	
	Body::SaveAs(bds, fileName, type, meshType, rho, g, symX, symY);
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
	surfs.Remove(id);
	if (surfs.IsEmpty())
		Body::ResetIdCount();
}

void BEM::JoinBody(int idDest, int idOrig) {
	const Body &orig = surfs[idOrig];
	Body &dest = surfs[idDest];
	dest.dt.fileName << "/" << orig.dt.fileName;
	
	try {
		dest.Append(orig.dt.mesh, rho, g);
		RemoveBody(idOrig);
	} catch (Exc e) {
		surfs.SetCount(surfs.size()-1);
		Print("\n" + Format(t_("Problem loading '%s': %s") + S("\n%s"), e));
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

void BEM::AddRevolution(double x, double y, double z, double size, UVector<Pointf> &vals) {
	try {
		Body &surf = surfs.Add();

		surf.dt.SetCode(Body::EDIT);
		surf.dt.mesh.AddRevolution(vals, size); 
		surf.dt.mesh.Translate(x, y, z);
		surf.dt.c0 = Point3D(0, 0, 0);
	} catch (Exc e) {
		surfs.SetCount(surfs.size() - 1);
		Print("\n" + Format(t_("Problem adding revolution surface: %s"), e));
		throw std::move(e);
	}	
}

void BEM::AddPolygonalPanel(double x, double y, double z, double size, UVector<Pointf> &vals) {
	try {
		Body &surf = surfs.Add();

		surf.dt.SetCode(Body::EDIT);
		surf.dt.mesh.AddPolygonalPanel(vals, size, true); 
		surf.dt.mesh.Translate(x, y, z);
		surf.dt.c0 = Point3D(0, 0, 0);
	} catch (Exc e) {
		surfs.SetCount(surfs.size() - 1);
		Print("\n" + Format(t_("Problem adding revolution surface: %s"), e));
		throw std::move(e);
	}	
}

void BEM::Extrude(int id, double dx, double dy, double dz, bool close) {
	try {
		Body &surf = surfs[id];
	
		if (surf.dt.mesh.volume != 0)
			throw Exc(t_("It is only possible to extrude a flat surface"));

		surf.dt.mesh.Extrude(dx, dy, dz, close);
		
	} catch (Exc e) {
		Print("\n" + Format(t_("Problem extruding surface: %s"), e));
		throw std::move(e);
	}	
}

void BEM::AddWaterSurface(int id, char c) {
	try {
		Body &surf = surfs.Add();

		surf.dt.SetCode(Body::EDIT);
		surf.dt.mesh.AddWaterSurface(surfs[id].dt.mesh, surfs[id].dt.under, c, roundVal, roundEps); 
		
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
		Print("\n" + Format(t_("Problem adding revolution surface: %s"), e));
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
	return StoreAsJsonFile(*this, fileName, true);
}


String Hydro::LoadSerialization(String fileName) {
	BEM::Print("\n\n" + Format(t_("Loading '%s'"), dt.file));
	
	if (!LoadFromJsonFile(*this, dt.file)) 
		return Format(t_("Error loading '%s'"), dt.file);
	
	dt.file = fileName;
	return String();
}
	
void Hydro::SaveSerialization(String fileName) const {
	BEM::Print("\n\n" + Format(t_("Saving '%s'"), fileName));
	if (!StoreAsJsonFile(*this, fileName, true)) {
		BEM::PrintError("\n" + Format(t_("Error saving '%s'"), fileName));
		throw Exc(Format(t_("Error saving '%s'"), fileName));
	}
}

void Hydro::SaveForce(FileOut &out, const Hydro::Forces &f) const {
	const String &sep = Bem().csvSeparator;
	
	out << sep;
	for (int ib = 0; ib < dt.Nb; ++ib) {			
		for (int idf = 0; idf < 6; ++idf) 			
			out << sep << ((dt.Nb > 1 ? (FormatInt(ib+1) + "-") : S("")) + BEM::StrDOF(idf)) << sep;
	}
	out << "\n";
	
	out  << "Head [deg]" << sep << "Frec [rad/s]";
	for (int ib = 0; ib < dt.Nb; ++ib) {			
		for (int idf = 0; idf < 6; ++idf) 			
			out << sep << "mag" << sep << "phase";
	}
	out << "\n";
	
	UVector<int> ow = GetSortOrderX(dt.w);
	UVector<int> oh = GetSortOrderX(dt.head);
	
	for (int ih = 0; ih < dt.Nh; ++ih) {
		out << dt.head[oh[ih]];
		for (int ifr = 0; ifr < dt.Nf; ++ifr) {
			for (int ib = 0; ib < dt.Nb; ++ib) {		
				out << sep;
				out << dt.w[ow[ifr]];				
				for (int idf = 0; idf < 6; ++idf) { 	
					out << sep;	
					if (IsNum(f[ib][oh[ih]](ow[ifr], idf)))	{
						const std::complex<double> &c = F_dim(f, oh[ih], ow[ifr], idf, ib);
						out << FormatDouble(abs(c)) << sep << FormatDouble(ToDeg(arg(c)));
					} else
						out << sep;
				}
				out << "\n";
			}
		}
	}
}	

void Hydro::SaveMD(FileOut &out) const {
	const String &sep = Bem().csvSeparator;
	
	out  << "Head [deg]" << sep << "Frec [rad/s]";
	for (int ib = 0; ib < dt.Nb; ++ib) {			
		for (int idf = 0; idf < 6; ++idf) 			
			out << sep << ((dt.Nb > 1 ? (FormatInt(ib+1) + "-") : S("")) + BEM::StrDOF(idf));
	}
	out << "\n";
	
	UVector<int> ow = GetSortOrderX(dt.w);
	
	for (int ih = 0; ih < dt.mdhead.size(); ++ih) {
		const std::complex<double> &hh = dt.mdhead[ih];
		out << Format("%s %s", FormatDouble(hh.real()), FormatDouble(hh.imag()));
		for (int ifr = 0; ifr < dt.Nf; ++ifr) {
			for (int ib = 0; ib < dt.Nb; ++ib) {		
				out << sep;
				out << dt.w[ow[ifr]];				
				for (int idf = 0; idf < 6; ++idf) { 	
					out << sep;	
					if (IsNum(dt.md[ib][ih][idf](ifr))) 
						out << FormatDouble(Md_dim(idf, ih, ifr));
				}
				out << "\n";
			}
		}
	}
}

void Hydro::SaveC(FileOut &out) const {
	const String &sep = Bem().csvSeparator;
	
	out  << "DoF";
	for (int idf = 0; idf < 6; ++idf) 			
		out << sep << BEM::StrDOF(idf);
	out << "\n";
		
	for (int ib = 0; ib < dt.Nb; ++ib) {		
		for (int idf1 = 0; idf1 < 6; ++idf1) { 	
			out << ((dt.Nb > 1 ? (FormatInt(ib+1) + "-") : S("")) + BEM::StrDOF(idf1));
			out << sep;	
			for (int idf2 = 0; idf2 < 6; ++idf2) {
				if (IsNum(dt.msh[ib].dt.C(idf1, idf2)))	
					out << FormatDouble(C_dim(ib, idf1, idf2));
				out << sep;
			}
			out << "\n";
		}
	}
}

void Hydro::SaveM(FileOut &out) const {
	const String &sep = Bem().csvSeparator;
	
	out  << "DoF";
	for (int idf = 0; idf < 6; ++idf) 			
		out << sep << BEM::StrDOF(idf);
	out << "\n";
		
	for (int ib = 0; ib < dt.Nb; ++ib) {		
		for (int idf1 = 0; idf1 < 6; ++idf1) { 	
			out << ((dt.Nb > 1 ? (FormatInt(ib+1) + "-") : S("")) + BEM::StrDOF(idf1));
			out << sep;	
			for (int idf2 = 0; idf2 < 6; ++idf2) {
				if (IsNum(dt.msh[ib].dt.M(idf1, idf2)))	
					out << FormatDouble(dt.msh[ib].dt.M(idf1, idf2));
				out << sep;
			}
			out << "\n";
		}
	}
}
	
void Hydro::SaveCSVMat(String fileName) const {
	BEM::Print("\n\n" + Format(t_("Saving '%s'"), fileName));
	
	String folder = GetFileFolder(fileName);
	String nname = GetFileTitle(fileName);
	String ext = GetFileExt(fileName);
	
	if (IsLoadedA())  {
		String files = AFX(folder, nname + "_A" + ext);
		FileOut out(files);
		if (!out.IsOpen())
			throw Exc(Format(t_("Impossible to save '%s'. File already used."), files));


		const String &sep = Bem().csvSeparator;
		
		out  << "Frec [rad/s]" << sep << "DoF";
		for (int ib = 0; ib < dt.Nb; ++ib) {			
			for (int idf = 0; idf < 6; ++idf) 			
				out << sep << ((dt.Nb > 1 ? (FormatInt(ib+1) + "-") : S("")) + BEM::StrDOF(idf));
		}
		out << "\n";
	
		
		if (IsLoadedA0()) {
			for (int ib = 0; ib < dt.Nb; ++ib) {		
				if (ib == 0)
					out << "0";				
				for (int idf1 = 0; idf1 < 6; ++idf1) { 	
					out << sep;	
					out << ((dt.Nb > 1 ? (FormatInt(ib+1) + "-") : S("")) + BEM::StrDOF(idf1));
					out << sep;	
					for (int idf2 = 0; idf2 < 6; ++idf2) {
						if (IsNum(dt.A0(idf1 + 6*ib, idf2 + 6*ib)))	
							out << FormatDouble(A0_dim(idf1, idf2));
						out << sep;
					}
					out << "\n";
				}
			}
		}
		UVector<int> ow = GetSortOrderX(dt.w);
		
		for (int ifr = 0; ifr < dt.Nf; ++ifr) {
			for (int ib = 0; ib < dt.Nb; ++ib) {		
				if (ib == 0)
					out << dt.w[ow[ifr]];				
				for (int idf1 = 0; idf1 < 6; ++idf1) { 	
					out << sep;	
					out << ((dt.Nb > 1 ? (FormatInt(ib+1) + "-") : S("")) + BEM::StrDOF(idf1));
					out << sep;	
					for (int idf2 = 0; idf2 < 6; ++idf2) {
						if (IsNum(dt.A[idf1 + 6*ib][idf2 + 6*ib][ow[ifr]]))	
							out << FormatDouble(A_dim(ow[ifr], idf1, idf2));
						out << sep;
					}
					out << "\n";
				}
			}
		}
		if (IsLoadedAinf()) {
			for (int ib = 0; ib < dt.Nb; ++ib) {		
				if (ib == 0)
					out << "inf";				
				for (int idf1 = 0; idf1 < 6; ++idf1) { 	
					out << sep;	
					out << ((dt.Nb > 1 ? (FormatInt(ib+1) + "-") : S("")) + BEM::StrDOF(idf1));
					out << sep;	
					for (int idf2 = 0; idf2 < 6; ++idf2) {
						if (IsNum(dt.Ainf(idf1 + 6*ib, idf2 + 6*ib)))	
							out << FormatDouble(Ainf_dim(idf1, idf2));
						out << sep;
					}
					out << "\n";
				}
			}
		}	
	}
	
	if (IsLoadedB())  {
		String files = AFX(folder, nname + "_B" + ext);
		FileOut out(files);
		if (!out.IsOpen())
			throw Exc(Format(t_("Impossible to save '%s'. File already used."), files));


		const String &sep = Bem().csvSeparator;
		
		out  << "Frec [rad/s]" << sep << "DoF";
		for (int ib = 0; ib < dt.Nb; ++ib) {			
			for (int idf = 0; idf < 6; ++idf) 			
				out << sep << ((dt.Nb > 1 ? (FormatInt(ib+1) + "-") : S("")) + BEM::StrDOF(idf));
		}
		out << "\n";
	
		
		UVector<int> ow = GetSortOrderX(dt.w);
		
		for (int ifr = 0; ifr < dt.Nf; ++ifr) {
			for (int ib = 0; ib < dt.Nb; ++ib) {		
				if (ib == 0)
					out << dt.w[ow[ifr]];				
				for (int idf1 = 0; idf1 < 6; ++idf1) { 	
					out << sep;	
					out << ((dt.Nb > 1 ? (FormatInt(ib+1) + "-") : S("")) + BEM::StrDOF(idf1));
					out << sep;	
					for (int idf2 = 0; idf2 < 6; ++idf2) {
						if (IsNum(dt.B[idf1 + 6*ib][idf2 + 6*ib][ow[ifr]]))	
							out << FormatDouble(B_dim(ow[ifr], idf1, idf2));
						out << sep;
					}
					out << "\n";
				}
			}
		}
	}
	
	if (IsLoadedC())  {
		String files = AFX(folder, nname + "_C" + ext);
		FileOut out(files);
		if (!out.IsOpen())
			throw Exc(Format(t_("Impossible to save '%s'. File already used."), files));


		SaveC(out);
	}
	
	if (IsLoadedM())  {
		String files = AFX(folder, nname + "_M" + ext);
		FileOut out(files);
		if (!out.IsOpen())
			throw Exc(Format(t_("Impossible to save '%s'. File already used."), files));


		SaveM(out);
	}
	
	if (IsLoadedFex())  {
		String files = AFX(folder, nname + "_Fex" + ext);
		FileOut out(files);
		if (!out.IsOpen())
			throw Exc(Format(t_("Impossible to save '%s'. File already used."), files));

		SaveForce(out, dt.ex);
	}
	
	if (IsLoadedMD())  {
		String files = AFX(folder, nname + "_MD" + ext);
		FileOut out(files);
		if (!out.IsOpen())
			throw Exc(Format(t_("Impossible to save '%s'. File already used."), files));

		SaveMD(out);
	}		
}

void Hydro::SaveCSVTable(String fileName) const {
	BEM::Print("\n\n" + Format(t_("Saving '%s'"), fileName));
	
	String folder = GetFileFolder(fileName);
	String nname = GetFileTitle(fileName);
	String ext = GetFileExt(fileName);
	
	if (IsLoadedA())  {
		String files = AFX(folder, nname + "_A" + ext);
		FileOut out(files);
		if (!out.IsOpen())
			throw Exc(Format(t_("Impossible to save '%s'. File already used."), files));


		const String &sep = Bem().csvSeparator;
		
		out  << "Frec [rad/s]";
		for (int ib = 0; ib < dt.Nb; ++ib) {			
			for (int idf1 = 0; idf1 < 6; ++idf1) 
				for (int idf2 = 0; idf2 < 6; ++idf2) 
					out << sep << ((dt.Nb > 1 ? (FormatInt(ib+1) + "-") : S("")) + BEM::StrDOF(idf1) + "-" + BEM::StrDOF(idf2));
		}
		out << "\n";
	
		if (IsLoadedA0()) {
			for (int ib = 0; ib < dt.Nb; ++ib) {		
				if (ib == 0)
					out << "0" << sep;				
				for (int idf1 = 0; idf1 < 6; ++idf1) { 	
					for (int idf2 = 0; idf2 < 6; ++idf2) {
						if (IsNum(dt.A0(idf1 + 6*ib, idf2 + 6*ib)))	
							out << FormatDouble(A0_dim(idf1, idf2));
						out << sep;
					}
				}
			}
		}
		out << "\n";
		
		UVector<int> ow = GetSortOrderX(dt.w);
		
		for (int ifr = 0; ifr < dt.Nf; ++ifr) {
			for (int ib = 0; ib < dt.Nb; ++ib) {		
				if (ib == 0)
					out << dt.w[ow[ifr]] << sep;				
				for (int idf1 = 0; idf1 < 6; ++idf1) { 	
					for (int idf2 = 0; idf2 < 6; ++idf2) {
						if (IsNum(dt.A[idf1 + 6*ib][idf2 + 6*ib][ow[ifr]]))	
							out << FormatDouble(A_dim(ow[ifr], idf1, idf2));
						out << sep;
					}
				}
				out << "\n";
			}
		}
		if (IsLoadedAinf()) {
			for (int ib = 0; ib < dt.Nb; ++ib) {		
				if (ib == 0)
					out << "inf" << sep;				
				for (int idf1 = 0; idf1 < 6; ++idf1) { 	
					for (int idf2 = 0; idf2 < 6; ++idf2) {
						if (IsNum(dt.Ainf(idf1 + 6*ib, idf2 + 6*ib)))	
							out << FormatDouble(Ainf_dim(idf1, idf2));
						out << sep;
					}
				}
			}
		}	
	}
	
	if (IsLoadedB())  {
		String files = AFX(folder, nname + "_B" + ext);
		FileOut out(files);
		if (!out.IsOpen())
			throw Exc(Format(t_("Impossible to save '%s'. File already used."), files));


		const String &sep = Bem().csvSeparator;
		
		out  << "Frec [rad/s]";
		for (int ib = 0; ib < dt.Nb; ++ib) {			
			for (int idf1 = 0; idf1 < 6; ++idf1) 
				for (int idf2 = 0; idf2 < 6; ++idf2) 
					out << sep << ((dt.Nb > 1 ? (FormatInt(ib+1) + "-") : S("")) + BEM::StrDOF(idf1) + "-" + BEM::StrDOF(idf2));
		}
		out << "\n";
		
		UVector<int> ow = GetSortOrderX(dt.w);
		
		for (int ifr = 0; ifr < dt.Nf; ++ifr) {
			for (int ib = 0; ib < dt.Nb; ++ib) {		
				if (ib == 0)
					out << dt.w[ow[ifr]] << sep;			
				for (int idf1 = 0; idf1 < 6; ++idf1) { 	
					for (int idf2 = 0; idf2 < 6; ++idf2) {
						if (IsNum(dt.B[idf1 + 6*ib][idf2 + 6*ib][ow[ifr]]))	
							out << FormatDouble(B_dim(ow[ifr], idf1, idf2));
						out << sep;
					}
				}
				out << "\n";
			}
		}
	}
	
	if (IsLoadedC())  {
		String files = AFX(folder, nname + "_C" + ext);
		FileOut out(files);
		if (!out.IsOpen())
			throw Exc(Format(t_("Impossible to save '%s'. File already used."), files));


		SaveC(out);
	}
	
	if (IsLoadedM())  {
		String files = AFX(folder, nname + "_M" + ext);
		FileOut out(files);
		if (!out.IsOpen())
			throw Exc(Format(t_("Impossible to save '%s'. File already used."), files));


		SaveM(out);
	}
	
	if (IsLoadedFex())  {
		String files = AFX(folder, nname + "_Fex" + ext);
		FileOut out(files);
		if (!out.IsOpen())
			throw Exc(Format(t_("Impossible to save '%s'. File already used."), files));

		SaveForce(out, dt.ex);
	}	
	
	if (IsLoadedMD())  {
		String files = AFX(folder, nname + "_MD" + ext);
		FileOut out(files);
		if (!out.IsOpen())
			throw Exc(Format(t_("Impossible to save '%s'. File already used."), files));

		SaveMD(out);
	}	
}

int IsTabSpace(int c) {
	if (c == '\t' || c == ' ' || c == '!')
		return true;
	return false;
}

UVector<int> NumSets(int num, int numsets) {
	ASSERT(numsets > 0);
	UVector<int> ret;
	ret.SetCount(numsets);
	
	for (int i = 0; numsets > 0; ++i) {
		int delta = int(num/numsets);
		ret[i] = delta;
		num -= delta;
		numsets--;
	}
	return ret;
}

String FormatWam(double d) {
	if (!IsNum(d))
		return "0.0";
	return (d >= 0 ? " " : "-") + Format("%12E", abs(d));
}

void LineParserWamit::LoadWamitJoinedFields(String _line) {	
	line = _line;
	fields.Clear();
	UVector<String> prefields = Split(line, IsTabSpace, true);
	for (int id = 0; id < prefields.size(); ++id) {
		String s = prefields[id];
		String ns;
		for (int i = 0; i < s.GetCount(); ++i) {	
			int c = s[i];
			if (c == '-') {
				if (i == 0)
					ns.Cat(c);
				else if (s[i-1] == 'E')
					ns.Cat(c);
				else {
					fields << ns;
					ns.Clear();
					ns.Cat(c);
				}
			} else
				ns.Cat(c);
		}
		fields << ns;
	}
}

void Hydro::LoadCase(String fileName, Function <bool(String, int)> Status) {
	dt.file = fileName;
	
	String ret;
	if (ToLower(GetFileName(fileName)) == "nemoh.cal")
		ret = static_cast<Nemoh&>(*this).Load(fileName);
	else if (ToLower(GetFileExt(fileName)) == ".in")
		ret = static_cast<Hams&>(*this).Load(fileName, Status); 
	else if (ToLower(GetFileExt(fileName)) == ".dat") 
		ret = static_cast<Aqwa&>(*this).Load(fileName, Status);
	else if (ToLower(GetFileExt(fileName)) == ".lis") 
		ret = static_cast<Aqwa&>(*this).Load(fileName, Status);
	else if (ToLower(GetFileExt(fileName)) == ".ah1") 
		ret = static_cast<Aqwa&>(*this).Load(fileName, Status);
	else if (ToLower(GetFileExt(fileName)) == ".nc") {
		UArray<Hydro> hydros;
		int num;
		ret = CapyNC_Load(fileName, hydros, num);
		if (ret.IsEmpty() && num > 0)
			*this = pick(First(hydros));
	}
#ifdef PLATFORM_WIN32	 
	else if (ToLower(GetFileExt(fileName)) == ".owr")
		ret = static_cast<OrcaWave&>(*this).Load(fileName, Status);
#endif	
	else if (ToLower(GetFileExt(fileName)) == ".yml")
		ret = static_cast<OrcaWave&>(*this).Load(fileName, Status);
	else
		ret = t_("Unknown BEM input format");
	
	if (!ret.IsEmpty())
		throw Exc(ret);
	
	if (IsNull(dt.rho))
		dt.rho = Bem().rho;
	if (IsNull(dt.g))
		dt.g = Bem().g;	
	
	AfterLoad();
}

void Hydro::SaveFolderCase(String folder, bool bin, int numCases, int numThreads, BEM_FMT solver, 
			bool withPotentials, bool withMesh, bool withQTF, bool x0z, bool y0z, UArray<Body> &lids) {
	if (solver == Hydro::CAPYTAINE || solver == Hydro::NEMOH || solver == Hydro::NEMOHv115 || solver == Hydro::NEMOHv3 || solver == Hydro::SEAFEM_NEMOH)
		static_cast<const Nemoh &>(*this).SaveFolder(folder, bin, numCases, solver, x0z, y0z);
	else if (solver == Hydro::CAPYTAINE_PY)
		static_cast<const Nemoh &>(*this).SaveFolder_Capy(folder, withPotentials, withMesh, x0z, y0z, lids);
	else if (solver == Hydro::HAMS)
		static_cast<const Hams &>(*this).SaveFolder(folder, bin, numCases, numThreads, x0z, y0z, lids);
	else if (solver == Hydro::ORCAWAVE_YML)
		static_cast<const OrcaWave &>(*this).SaveFolder_OW_YML(folder, bin, numThreads, withPotentials, withMesh, withQTF, x0z, y0z);
	else if (solver == Hydro::AQWA_DAT)
		static_cast<const Aqwa &>(*this).SaveCaseDat(folder, numThreads, withPotentials, withQTF, x0z, y0z);
	else if (solver == Hydro::BEMROSETTA_H5) {
		dt.solver = BEMROSETTA_H5;
		if (!dt.msh.IsEmpty() && !IsLoadedPotsIncB()) 
			GetPotentialsIncident();
		if (IsLoadedPotsIncB()) 
			GetForcesFromPotentials(dt.pots_inc_bmr, dt.fk_pot_bmr);
		
		dt.fk = clone(dt.fk_pot_bmr);
		static_cast<BemioH5&>(*this).Save(AFX(folder, GetFileTitle(folder) + ".h5"));
	}
	else
		throw Exc(t_("Format is not supported"));
}

void Hydro::BeforeSaveCase(String folderBase, int numCases, bool deleteFolder) const {
	if (numCases < 1)
		throw Exc(Format(t_("Number cases must be higher than 1 (%d)"), numCases));
	
	if (numCases > dt.Nf)
		throw Exc(Format(t_("Number of cases %d must not be higher than number of frequencies %d"), numCases, dt.Nf));
	
	if (deleteFolder) {		// If called from GUI, user has been warned
		if (!DeleteFileDeepWildcardsX(folderBase))
			throw Exc(Format(t_("Impossible to clean folder '%s'. Maybe it is in use"), folderBase));
		Sleep(100);
	}
	if (!DirectoryCreateX(folderBase))
		throw Exc(Format(t_("Problem creating '%s' folder"), folderBase));
}

UVector<String> Hydro::Check(BEM_FMT type) const {
	UVector<String> ret;
	
	if (IsNull(dt.rho) || dt.rho < 0 || dt.rho > 10000)
		 ret << Format(t_("Incorrect rho %s"), FormatDoubleEmpty(dt.rho));
	if (IsNull(dt.g) || dt.g < 0 || dt.g > 100)
		ret << Format(t_("Incorrect g %s"), FormatDoubleEmpty(dt.g));
	
	if (IsNull(dt.h) || dt.h < -1 || dt.h > 100000)
		ret << Format(t_("Incorrect depth %s"), FormatDoubleEmpty(dt.h));

	if (IsNull(dt.Nf) || dt.Nf < 1 || dt.Nf > 1000)
		ret << Format(t_("Incorrect number of frequencies %s"), FormatIntEmpty(dt.Nf));
	
	if (IsNull(dt.Nh) || dt.Nh < 1 || dt.Nh > 1000)
		ret << Format(t_("Incorrect number of headings %s"), FormatIntEmpty(dt.Nh));
	
	if (type == BEM_FMT::HAMS)
		ret = static_cast<const Hams&>(*this).Check();
	
	if (First(dt.w) <= 0.01)
		ret << Format(t_("First frequency %f < 0.01 is too low"), First(dt.w));
	
	return ret;
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
