#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#include <STEM4U/Utility.h>

bool Nemoh::Load(String file, double) {
	try {
		String ext = GetFileExt(file); 

		if (ext == ".tec") {
			String folder = GetFileFolder(file);
			String folderTitle = GetFileName(folder);
			if (ToLower(folderTitle) != "results") 
				throw Exc(Format(t_(".tec file '%s' should have to be in 'results' folder"), file));
			bool found = false;
			String upperFolder = GetUpperFolder(folder);
			for (FindFile ff(AppendFileName(upperFolder, "*.*")); ff; ++ff) {
				if (ff.IsFile()) {
					if (ToLower(ff.GetName()) == "nemoh.cal") {
						file = ff.GetPath();
						found = true;
						break;
					}
				}
			}
			if (!found)
				throw Exc(Format(t_("nemoh.cal file not found in '%s' folder"), upperFolder));
		}
	
		if (ext == ".cal" || ext == ".tec")
			hd().code = Hydro::NEMOH;
		else
			hd().code = Hydro::SEAFEM_NEMOH;
	
		hd().file = file;
		hd().name = GetFileTitle(GetFileFolder(file));
		folder = GetFileFolder(file);
		hd().len = 1;
		hd().dimen = true;
		hd().Nb = Null;
	
		String fileCal;
		BEMData::Print("\n\n" + Format(t_("Loading '%s'"), file));
		if (hd().code == Hydro::NEMOH) 
			fileCal = file;
		else 
			fileCal = AppendFileName(folder, "Nemoh_output/Nemoh.cal");
		if (!Load_Cal(fileCal)) 
			throw Exc(Format(t_("File '%s' not found"), fileCal));
		
		String fileRad, folderForces;
		if (hd().code == Hydro::NEMOH) {
			BEMData::Print(S("\n- ") + t_("Hydrostatics file(s) 'Mesh/Hydrostatics*.dat'"));
			if (!Load_Hydrostatics())
				BEMData::PrintWarning(S(": **") + t_("Not found") + "**");
			BEMData::Print(S("\n- ") + t_("KH file(s) 'Mesh/KH*.dat'"));
			if (!Load_KH())
				BEMData::PrintWarning(S(": **") + t_("Not found") + "**");
			fileRad = AppendFileName(folder, AppendFileName("Results", "RadiationCoefficients.tec"));
			folderForces = folder;
		} else {
			if (!Load_Inf(file)) 
				throw Exc(Format(t_("File '%s' not found"), file));

			fileRad = AppendFileName(folder, AppendFileName("Nemoh_output/Results", "RadiationCoefficients.tec"));
			folderForces = AppendFileName(folder, "Nemoh_output");
		} 
		BEMData::Print(S("\n- ") + t_("Radiation file 'RadiationCoefficients.tec'"));
		if (!Load_Radiation(fileRad))
			BEMData::PrintWarning(S(": **") + t_("Not found") + "**");

		BEMData::Print(S("\n- ") + t_("Excitation force file 'ExcitationForce.tec'"));
		if (!Load_Excitation(folderForces))
			BEMData::PrintWarning(S(": **") + t_("Not found") + "**");
		
		BEMData::Print(S("\n- ") + t_("Diffraction force file 'DiffractionForce.tec'"));
		if (!Load_Diffraction(folderForces))
			BEMData::PrintWarning(S(": **") + t_("Not found") + "**");
		BEMData::Print(S("\n- ") + t_("Froude Krylov file 'FKForce.tec'"));
		if (!Load_FroudeKrylov(folderForces))
			BEMData::PrintWarning(S(": **") + t_("Not found") + "**");
		
		if (hd().code == Hydro::NEMOH) {
			if (!hd().dof.IsEmpty()) {
				BEMData::Print(S("\n- ") + t_("IRF file(s) 'IRF.tec'"));
				if (!Load_IRF(AppendFileName(folder, AppendFileName("Results", "IRF.tec"))))
					BEMData::PrintWarning(S(": **") + t_("Not found") + "**");
			}
		}
		if (IsNull(hd().Nb))
			return false;
	} catch (Exc e) {
		BEMData::PrintError(Format("\n%s: %s", t_("Error"), e));
		hd().lastError = e;
		return false;
	}
	return true;
}

bool Nemoh::Load_Cal(String fileName) {	
	if (!datacal.Load(fileName))
		return false;

	hd().rho = datacal.rho;
	hd().g = datacal.g;
	hd().h = datacal.h;
	if (hd().h == 0)
		hd().h = -1;
	hd().Nb = datacal.bodies.size();
	for (int i = 0; i < hd().Nb; ++i) 	
		hd().names << GetFileTitle(datacal.bodies[i].meshFile);
	hd().dof.SetCount(hd().Nb, 6);
	//for (int i = 0; i < hd().Nb; ++i)
	//	hd().dof[i] = data.bodies[i].ndof;
	hd().Nf = datacal.Nf;
	LinSpaced(hd().w, hd().Nf, datacal.minF, datacal.maxF); 
 	hd().T.SetCount(hd().Nf);
    for (int i = 0; i < hd().Nf; ++i) 
		hd().T[i] = 2*M_PI/hd().w[i];  
   
	hd().Nh = datacal.Nh;  						
    LinSpaced(hd().head, hd().Nh, datacal.minH, datacal.maxH); 		

	hd().dataFromW = true;
	
	if (datacal.irf) {
		hd().Tirf.resize(int(datacal.irfDuration/datacal.irfStep));
		for (int i = 0; i < hd().Tirf.size(); ++i) 
			hd().Tirf[i] = i*datacal.irfStep;
	}
	return true;
}

int NemohCal::GetNumArgs(const FieldSplit &f) {
	for (int i = 0; i < f.size(); ++i) {
		double num = ScanDouble(f.GetText(i));
		if (IsNull(num))
			return i;
	}
	return f.size();
}
	
void NemohCal::LoadFreeSurface(const FileInLine &in, const FieldSplit &f) {
	nFreeX = f.GetInt(0);	nFreeY = f.GetInt(1);	domainX = f.GetDouble(2);	domainY = f.GetDouble(3);
	if (nFreeX < 0)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect number of points in x direction %s"), f.GetText(0)));
	if (nFreeX > 0 && (nFreeY <= 0 || domainX < 0 || domainY < 0))
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect free surface elevation %s"), f.GetText()));	
}


void NemohCal::LoadKochin(const FileInLine &in, const FieldSplit &f) {
	nKochin = f.GetInt(0);	minK = f.GetDouble(1);	maxK = f.GetDouble(2);
	if (nKochin < 0)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect number of Kochin function directions %s"), f.GetText(0)));
	if (nKochin > 0) {
		if (minK < -360)
			throw Exc(in.Str() + "\n"  + Format(t_("Incorrect Kochin direction %s"), f.GetText(1)));
		if (maxK > 360)
			throw Exc(in.Str() + "\n"  + Format(t_("Incorrect Kochin direction %s"), f.GetText(2)));
		if (maxK <= minK)
			throw Exc(in.Str() + "\n"  + Format(t_("Minimum Kochin direction %s has to be lower than maximum direction %s"), f.GetText(1), f.GetText(2)));	
	}
}

bool NemohCal::Load(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	String line;
	FieldSplit f(in);
	f.IsSeparator = IsTabSpace;
	
	in.GetLine();
	f.Load(in.GetLine());	
	rho = f.GetDouble(0);
	if (rho < 0 || rho > 10000)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect rho %s"), f.GetText(0)));
	
	f.Load(in.GetLine());	
	g = f.GetDouble(0);
	if (g < 0 || g > 100)
		throw Exc(in.Str() + "\n" + Format(t_("Incorrect g %s"), f.GetText(0)));
	
	f.Load(in.GetLine());	
	String sh = ToLower(f.GetText(0));
	if (sh == "inf")
		h = -1;
	else {
		h = ScanDouble(sh);
		if (h < 0 || h > 100000)
			throw Exc(in.Str() + "\n"  + Format(t_("Incorrect depth %s"), f.GetText(0)));
		else if (h == 0)
			h = -1;
	}
	f.Load(in.GetLine());	xeff = f.GetDouble(0);	yeff = f.GetDouble(1);
	in.GetLine();
	f.Load(in.GetLine());	int Nb = f.GetInt(0);
	if (Nb < 1 || Nb > 100)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect number of bodies %s"), f.GetText(0)));
	bodies.SetCount(Nb);
	for (int ib = 0; ib < Nb; ++ib) {
		BemBody &body = bodies[ib];
		in.GetLine();
		f.Load(in.GetLine());	body.meshFile = f.GetText(0);
		f.Load(in.GetLine());	body.npoints = f.GetInt(0);		body.npanels = f.GetInt(1);
		String file = AppendFileName(GetFileFolder(fileName), body.meshFile);
		if (FileExists(file)) {
			body.meshFile = file;
			MeshData dat;
			bool x0z;
			if (dat.LoadDatNemoh(file, x0z).IsEmpty()) {
				body.npoints = dat.mesh.GetNumNodes();
				body.npanels = dat.mesh.GetNumPanels();
			}
		} else
			throw Exc(in.Str() + "\n"  + Format(t_("Mesh file '%s ' not found"), file));
		
		if (body.npoints < 1 || body.npoints > 100000000)
			throw Exc(in.Str() + "\n"  + Format(t_("Incorrect number of points %s"), f.GetText(0)));
		if (body.npanels < 1 || body.npanels > 100000000)
			throw Exc(in.Str() + "\n"  + Format(t_("Incorrect number of panels %s"), f.GetText(1)));	
		f.Load(in.GetLine());	
		body.ndof = f.GetInt(0);
		if (body.ndof < 0 || body.ndof > 6)
			throw Exc(in.Str() + "\n"  + Format(t_("Incorrect DOF number %s in body %d"), f.GetText(0), ib+1));
		for (int idf = 0; idf < body.ndof; ++idf) {
			f.Load(in.GetLine());
			int type = f.GetInt(0);
			bool x = f.GetDouble(1) > 0;		
			bool y = f.GetDouble(2) > 0;
			bool z = f.GetDouble(3) > 0;
			double cx = f.GetDouble(4);
			double cy = f.GetDouble(5);
			double cz = f.GetDouble(6);
			if (type == 1) {
				if (x) 
					body.dof[Hydro::SURGE] = true;
				else if (y)
					body.dof[Hydro::SWAY] = true;
				else if (z)
					body.dof[Hydro::HEAVE] = true;
			} else if (type == 2) {
				body.c0[0] = cx;
				body.c0[1] = cy;
				body.c0[2] = cz;
				if (x) 
					body.dof[Hydro::ROLL] = true;
				else if (y)
					body.dof[Hydro::PITCH] = true;
				else if (z)
					body.dof[Hydro::YAW] = true;
			} else
				throw Exc(in.Str() + "\n"  + Format(t_("Incorrect DOF type %d set in body %d"), f.GetText(0), ib+1));
		}
		f.Load(in.GetLine());	int nforces = f.GetInt(0);
		in.GetLine(nforces);	// Discarded
		f.Load(in.GetLine());	int nadditional = f.GetInt(0);	
		in.GetLine(nadditional);// Discarded	
	}
	in.GetLine();
	f.Load(in.GetLine());	Nf = f.GetInt(0);	minF = f.GetDouble(1);	maxF = f.GetDouble(2);
	if (Nf < 1 || Nf > 1000)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect number of frequencies %s"), f.GetText(0)));
	if (minF < 0)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect frequency %s"), f.GetText(1)));
	if (maxF < minF)
		throw Exc(in.Str() + "\n"  + Format(t_("Minimum frequency %s has to be lower than maximum frequency %s"), f.GetText(1), f.GetText(2)));	
	
	f.Load(in.GetLine());	Nh = f.GetInt(0);	minH = f.GetDouble(1);	maxH = f.GetDouble(2);
	if (Nh < 1 || Nh > 1000)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect number of headings %s"), f.GetText(0)));
	if (minH < -360)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect direction %s"), f.GetText(1)));
	if (maxH > 360)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect direction %s"), f.GetText(2)));
	if (maxH < minH)
		throw Exc(in.Str() + "\n"  + Format(t_("Minimum direction %s has to be lower than maximum direction %s"), f.GetText(1), f.GetText(2)));	
	
	in.GetLine();
	f.Load(in.GetLine());	irf = f.GetInt(0) > 0;	irfStep = f.GetDouble(1);	irfDuration = f.GetDouble(2);
	if (irf && irfStep <= 0)
		throw Exc(in.Str() + "\n"  + Format(t_("Incorrect IRF step %s"), f.GetText(1)));
	if (irf && irfDuration <= irfStep)
		throw Exc(in.Str() + "\n"  + Format(t_("IRF step %s has to be lower than duration %s"), f.GetText(1), f.GetText(2)));	
	f.Load(in.GetLine());	showPressure = f.GetInt(0) > 0;
	
	bool loadedFree = false, loadedKochin = false;
	f.Load(in.GetLine());	
	if (GetNumArgs(f) == 4) {
		LoadFreeSurface(in, f);
		loadedFree = true;
	} else if (GetNumArgs(f) == 3) {
		LoadKochin(in, f);
		loadedKochin = true;
	} else
		throw Exc(in.Str() + "\n"  + Format(t_("Unexpected data %s"), f.GetText()));			
	
	f.Load(in.GetLine());	
	if (GetNumArgs(f) == 4) {
		if (!loadedFree)
			LoadFreeSurface(in, f);
		else
			throw Exc(in.Str() + "\n"  + Format(t_("Free surface data is already loaded %s"), f.GetText()));				
	} else if (GetNumArgs(f) == 3) {
		if (!loadedKochin)
			LoadKochin(in, f);
		else
			throw Exc(in.Str() + "\n"  + Format(t_("Kochin data is already loaded %s"), f.GetText()));				
	} else if (GetNumArgs(f) > 0)
		throw Exc(in.Str() + "\n"  + Format(t_("Unexpected data %s"), f.GetText()));			
	else {
		if (!loadedFree) {
			nFreeX = nFreeY = 0;	
			domainX = domainY = 0;
		}
		if (!loadedKochin) {	
			nKochin = 0;	
			minK = maxK = 0;
		}
	}
	return true;
}

String FormatDoubleEmpty(double val) {
	if (IsNull(val))
		return t_("'empty'");
	else
		return FormatDouble(val);
}

String FormatIntEmpty(int val) {
	if (IsNull(val))
		return t_("'empty'");
	else
		return FormatInt(val);
}

Vector<String> BemCal::Check(bool oneBody) {
	Vector<String> ret;
	
	if (IsNull(h) || h < -1 || h > 100000)
		ret << Format(t_("Incorrect depth %s"), FormatDoubleEmpty(h));

	if (IsNull(Nf) || Nf < 1 || Nf > 1000)
		ret << Format(t_("Incorrect number of frequencies %s"), FormatIntEmpty(Nf));
	if (IsNull(minF) || minF < 0)
		ret << Format(t_("Incorrect min frequency %s"), FormatDoubleEmpty(minF));
	if (IsNull(maxF) || maxF < minF)
		ret << Format(t_("Minimum frequency %s has to be lower than maximum frequency %s"), FormatDoubleEmpty(minF), FormatDoubleEmpty(maxF));	
	
	if (IsNull(Nh) || Nh < 1 || Nh > 1000)
		ret << Format(t_("Incorrect number of headings %s"), FormatIntEmpty(Nh));
	if (IsNull(minH) || minH < -180)
		ret << Format(t_("Incorrect min heading %s"), FormatDoubleEmpty(minH));
	if (IsNull(maxH) || maxH > 180)
		ret << Format(t_("Incorrect max heading %s"), FormatDoubleEmpty(maxH));
	if (maxH < minH)
		ret << Format(t_("Minimum heading %s has to be lower than maximum heading %s"), FormatDoubleEmpty(minH), FormatDoubleEmpty(maxH));	

	if(bodies.size() < 1)
		ret << t_("The case has to include at least one body");
	if (oneBody && bodies.size() > 1)
		ret << t_("This solver just processes one body");
	
	return ret;
}
	
Vector<String> NemohCal::Check() {
	Vector<String> ret;
	
	ret.Append(BemCal::Check(false));
	
	if (IsNull(rho) || rho < 0 || rho > 10000)
		 ret << Format(t_("Incorrect rho %s"), FormatDoubleEmpty(rho));
	if (IsNull(g) || g < 0 || g > 100)
		ret << Format(t_("Incorrect g %s"), FormatDoubleEmpty(g));

	if (irf) {
		if (IsNull(irfStep) || irfStep <= 0)
			ret << Format(t_("Incorrect IRF step %s"), FormatDoubleEmpty(irfStep));
		if (IsNull(irfDuration) || irfDuration <= irfStep)
			ret << Format(t_("IRF step %s has to be lower than duration %s"), FormatDoubleEmpty(irfStep), FormatDoubleEmpty(irfDuration));	
	}
	
	if (IsNull(nFreeX) || nFreeX < 0)
		ret << Format(t_("Incorrect number of points in x direction %s (0 for no free surface calculation)"), FormatIntEmpty(nFreeX));
	if (nFreeX > 0) {
		if (IsNull(nFreeY) || nFreeY <= 0)
			ret << Format(t_("Incorrect number of points in x direction %s"), FormatIntEmpty(nFreeY));
		if (IsNull(domainX) || domainX <= 0)
			ret << Format(t_("Incorrect free surface domain X %s"), FormatDoubleEmpty(domainX));
		if (IsNull(domainY) || domainY <= 0)
			ret << Format(t_("Incorrect free surface domain Y %s"), FormatDoubleEmpty(domainY));
	}
	
	if (IsNull(nKochin) || nKochin < 0)
		ret << Format(t_("Incorrect number of Kochin function directions %s"), FormatIntEmpty(nKochin));
	if (nKochin > 0) {
		if (IsNull(minK) || minK < -180)
			ret << Format(t_("Incorrect Kochin direction %s"), FormatDoubleEmpty(minK));
		if (IsNull(maxK) || maxK < 180)
			ret << Format(t_("Incorrect Kochin direction %s"),FormatDoubleEmpty(minK));
		if (maxK <= minK)
			Format(t_("Minimum Kochin direction %s has to be lower than maximum direction %s"), FormatDoubleEmpty(minK), FormatDoubleEmpty(maxK));	
	}	
	
	return ret;
}

void NemohCal::Save_Id(String folder) const {
	String fileName = AppendFileName(folder, "ID.dat");
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to create '%s'"), fileName));
	out << "1\n.";
}

void NemohCal::Save_Bat(String folder, String batname, String caseFolder, bool bin, String preName, String solvName, String postName) const {
	String fileName = AppendFileName(folder, batname);
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to create '%s'"), fileName));
	if (!IsNull(caseFolder))
		out << "cd \"" << caseFolder << "\"\n";
	String strBin;
	if (bin)
		strBin = AppendFileName(caseFolder.IsEmpty() ? "." : "..", "bin");
	if (preName.IsEmpty()) {
		if (solvName == "capytaine")
			out << "\"" << solvName << "\"\n";
		else
			out << "\"" << AppendFileName(strBin, solvName) << "\" -all\n";
	} else
		out << "\"" << AppendFileName(strBin, preName) << "\"\n"
			<< "\"" << AppendFileName(strBin, solvName) << "\"\n"
			<< "\"" << AppendFileName(strBin, postName) << "\"";
}

void NemohCal::Save_Input(String folder) const {
	String fileName = AppendFileName(folder, "Input.txt");
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to create '%s'"), fileName));
	
	out << "--- Calculation parameters ------------------------------------------------------------------------------------" << "\n";
	out << "0          ! Indiq_solver   ! - ! Solver (0) Direct Gauss (1) GMRES (2) GMRES with FMM acceleration (2 not implemented yet)" << "\n";
	out << "20         ! IRES           ! - ! Restart parameter for GMRES" << "\n";
	out << "5.E-07     ! TOL_GMRES      ! - ! Stopping criterion for GMRES" << "\n";
	out << "100        ! MAXIT          ! - ! Maximum iterations for GMRES" << "\n";
	out << "1          ! Sav_potential  ! -	! Save potential for visualization";
}
	

String NemohHeader(String str) {
	String ret = "--- " << str << " ";
	ret << String('-', 130 - ret.GetCount());
	return ret;
}

String NemohField(String str, int length) {
	String ret = str;
	if (length > ret.GetCount())
		ret << String(' ', length - ret.GetCount());
	return ret + " ";
}

BemBody::BemBody() {
	dof.SetCount(6, false);	
	cg = Eigen::Vector3d::Zero();
	mass.setConstant(6, 6, 0);
	linearDamping.setConstant(6, 6, 0);
	quadraticDamping.setConstant(6, 6, 0);
	hydrostaticRestoring.setConstant(6, 6, 0);
	externalRestoring.setConstant(6, 6, 0);
}
	
int BemBody::GetNDOF() const {
	int ret = 0;
	for (auto &d : dof)
		ret += d;
	return ret;
}

void BemCal::BeforeSave(String folderBase, int numCases, bool deleteFolder) const {
	if (numCases < 1)
		throw Exc(Format(t_("Number cases must be higher than 1 (%d)"), numCases));
	
	if (numCases > Nf)
		throw Exc(Format(t_("Number of cases %d must not be higher than number of frequencies %d"), numCases, Nf));
	
	if (deleteFolder) {		// If called from GUI, user has been warned
		if (!DeleteFileDeepWildcardsX(folderBase))
			throw Exc(Format(t_("Impossible to clean folder '%s'. Maybe it is in use"), folderBase));
		Sleep(100);
	}
	if (!DirectoryCreateX(folderBase))
		throw Exc(Format(t_("Problem creating '%s' folder"), folderBase));
}

void NemohCal::SaveFolder(String folderBase, bool bin, int numCases, const BEMData &bem, int solver) const {
	SaveFolder0(folderBase, bin, 1, bem, true, solver);
	if (numCases > 1)
		SaveFolder0(folderBase, bin, numCases, bem, false, solver);
}

void NemohCal::SaveFolder0(String folderBase, bool bin, int numCases, const BEMData &bem, 
							bool deleteFolder, int solver) const {
	BeforeSave(folderBase, numCases, deleteFolder);
	
	#define MIN_F_NEMOH 0.01
	
	double fixminF = minF;
	if (fixminF < MIN_F_NEMOH)
		fixminF = MIN_F_NEMOH;
	
	Vector<int> valsf;
	int _nf;
	double _minf, _maxf;
	int ifr = 0;
	Vector<double> freqs;
	if (numCases > 1) { 
		LinSpaced(freqs, Nf, fixminF, maxF);
		valsf = NumSets(Nf, numCases);
	}
	
	String preName, solvName, postName, batName = "Nemoh";
	if (solver == BemCal::CAPYTAINE) {
		solvName = "capytaine";
		batName = "Capytaine_bat";
	} else if (bin) {
		String binResults = AppendFileName(folderBase, "bin");
		if (!DirectoryCreateX(binResults))
			throw Exc(Format(t_("Problem creating '%s' folder"), binResults));
		
		if (solver == BemCal::NEMOHv115) {
			solvName = GetFileName(bem.nemohPathNew);
			String destNew = AppendFileName(binResults, solvName);
			if (!FileCopy(bem.nemohPathNew, destNew)) 
				throw Exc(Format(t_("Problem copying preprocessor file from '%s'"), bem.nemohPathNew));					
		} else if (solver == BemCal::NEMOH) {
			preName = GetFileName(bem.nemohPathPreprocessor);
			String destProprocessor = AppendFileName(binResults, preName);
			if (!FileCopy(bem.nemohPathPreprocessor, destProprocessor)) 
				throw Exc(Format(t_("Problem copying preprocessor file from '%s'"), bem.nemohPathPreprocessor));		
			solvName = GetFileName(bem.nemohPathSolver);
			String destSolver = AppendFileName(binResults, solvName);
			if (!FileCopy(bem.nemohPathSolver, destSolver)) 
				throw Exc(Format(t_("Problem copying solver file from '%s'"), bem.nemohPathSolver));		
			postName = GetFileName(bem.nemohPathPostprocessor);
			String destPostprocessor = AppendFileName(binResults, postName);
			if (!FileCopy(bem.nemohPathPostprocessor, destPostprocessor)) 
				throw Exc(Format(t_("Problem copying postprocessor file from '%s'"), bem.nemohPathPostprocessor));		
		} else if (solver == BemCal::CAPYTAINE) 
			solvName = "capytaine";

	} else {
		if (solver == BemCal::NEMOHv115) 
			solvName = "nemoh";
		else {
			preName = "preprocessor";
			solvName = "solver";
			postName = "postprocessor";
		}
	}
		
	//String sumcases;
	for (int i = 0; i < numCases; ++i) {
		String folder;
		if (numCases > 1) {
			folder = AppendFileName(folderBase, Format("%s_Part_%d", batName, i+1));
			if (!DirectoryCreateX(folder))
				throw Exc(Format(t_("Problem creating '%s' folder"), folder));
			//sumcases << " " << AppendFileName(folder, "Nemoh.cal");
			_minf = freqs[ifr];
			int deltaf = valsf[i];
			_maxf = freqs[ifr + deltaf - 1];
			_nf = deltaf;
			ifr += deltaf;
		} else {
			folder = folderBase;
			_nf = Nf;
			_minf = fixminF;
			_maxf = maxF;
		}
		Save_Id(folder);
		Save_Input(folder);
		String folderMesh = AppendFileName(folder, "mesh");
		if (!DirectoryCreateX(folderMesh))
			throw Exc(Format(t_("Problem creating '%s' folder"), folderMesh));
	
		for (int ib = 0; ib < bodies.size(); ++ib) {
			String name = GetFileName(bodies[ib].meshFile);
			name = RemoveAccents(name);
			name.Replace(" ", "_");
			String dest = AppendFileName(folderMesh, name);
			if (GetFileExt(name) == ".dat") {
				if (!FileCopy(bodies[ib].meshFile, dest)) 
					throw Exc(Format(t_("Problem copying mesh file from '%s'"), bodies[ib].meshFile));
			} else {
				MeshData mesh;
				String err = mesh.Load(bodies[ib].meshFile);
				if (!err.IsEmpty())
					throw Exc(err);
				mesh.SaveDatNemoh(dest, mesh.mesh, false);
			}
		}
		Save_Cal(folder, _nf, _minf, _maxf);
		
		String folderResults = AppendFileName(folder, "results");
		if (!DirectoryCreateX(folderResults))
			throw Exc(Format(t_("Problem creating '%s' folder"), folderResults));
		
		if (bin && !GetFileName(bem.nemohPathGREN).IsEmpty()) {
			String destGREN = AppendFileName(folder, GetFileName(bem.nemohPathGREN));
			if (!FileCopy(bem.nemohPathGREN, destGREN)) 
				throw Exc(Format(t_("Problem copying gren file '%s'"), bem.nemohPathGREN));
		}
		
		if (numCases > 1) 
			Save_Bat(folderBase, Format("%s_Part_%d.bat", batName, i+1), Format("%s_Part_%d", batName, i+1), bin, preName, solvName, postName);
		else
			Save_Bat(folder, Format("%s.bat", batName), Null, bin, preName, solvName, postName);
	}
}
	
void NemohCal::Save_Cal(String folder, int _nf, double _minf, double _maxf) const {
	String fileName = AppendFileName(folder, "Nemoh.cal");
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'"), fileName));
	
	int cp = 28;	
	out << NemohHeader("Environment - Created with BEMRosetta") << "\n";
	out << NemohField(Format("%f", rho), cp) 		   << "! RHO             ! KG/M**3   ! Fluid specific volume" << "\n";
	out << NemohField(Format("%f", g), cp)   		   << "! G               ! M/S**2    ! Gravity " << "\n";
	double seth = h;
	if (h < 0)
		seth = 0;
	out << NemohField(Format("%f", seth), cp)   	   << "! DEPTH           ! M         ! Water depth" << "\n";
	out << NemohField(Format("%f %f", xeff, yeff), cp) << "! XEFF YEFF       ! M         ! Wave measurement point" << "\n";
	
	out << NemohHeader("Description of floating bodies") << "\n";
	out << NemohField(Format("%d", bodies.size()), cp) << "! Number of bodies" << "\n";
	
	for (int i = 0; i < bodies.size(); ++i) {
		const BemBody &b = bodies[i];
		out << NemohHeader(Format("Body %d", i+1)) << "\n";	
		String name = GetFileName(b.meshFile);
		name = RemoveAccents(name);
		name.Replace(" ", "_");
		String file = AppendFileName("mesh", name);
		
		out << NemohField(Format("%s", file), cp) << "! Name of mesh file" << "\n";
		out << NemohField(Format("%d %d", b.npoints, b.npanels), cp) << "! Number of points and number of panels" << "\n";	
		out << NemohField(Format("%d", b.GetNDOF()), cp) << "! Number of degrees of freedom" << "\n";	
		if (b.dof[Hydro::SURGE])
			out << NemohField("1 1. 0. 0. 0. 0. 0.", cp) << "! Surge" << "\n";	
		if (b.dof[Hydro::SWAY])
			out << NemohField("1 0. 1. 0. 0. 0. 0.", cp) << "! Sway" << "\n";	
		if (b.dof[Hydro::HEAVE])
			out << NemohField("1 0. 0. 1. 0. 0. 0.", cp) << "! Heave" << "\n";	
		if (b.dof[Hydro::ROLL])
			out << NemohField(Format("2 1. 0. 0. %.2f %.2f %.2f", b.c0[0], b.c0[1], b.c0[2]), cp) << "! Roll about a point" << "\n";	
		if (b.dof[Hydro::PITCH])
			out << NemohField(Format("2 0. 1. 0. %.2f %.2f %.2f", b.c0[0], b.c0[1], b.c0[2]), cp) << "! Pitch about a point" << "\n";	
		if (b.dof[Hydro::YAW])		
			out << NemohField(Format("2 0. 0. 1. %.2f %.2f %.2f", b.c0[0], b.c0[1], b.c0[2]), cp) << "! Yaw about a point" << "\n";	
		out << NemohField(Format("%d", b.GetNDOF()), cp) << "! Number of resulting generalised forces" << "\n";	
		if (b.dof[Hydro::SURGE])
			out << NemohField("1 1. 0. 0. 0. 0. 0.", cp) << "! Force in x direction" << "\n";	
		if (b.dof[Hydro::SWAY])
			out << NemohField("1 0. 1. 0. 0. 0. 0.", cp) << "! Force in y direction" << "\n";	
		if (b.dof[Hydro::HEAVE])
			out << NemohField("1 0. 0. 1. 0. 0. 0.", cp) << "! Force in z direction" << "\n";	
		if (b.dof[Hydro::ROLL])
			out << NemohField(Format("2 1. 0. 0. %.2f %.2f %.2f", b.c0[0], b.c0[1], b.c0[2]), cp) << "! Moment force in x direction about a point" << "\n";	
		if (b.dof[Hydro::PITCH])
			out << NemohField(Format("2 0. 1. 0. %.2f %.2f %.2f", b.c0[0], b.c0[1], b.c0[2]), cp) << "! Moment force in y direction about a point" << "\n";	
		if (b.dof[Hydro::YAW])		
			out << NemohField(Format("2 0. 0. 1. %.2f %.2f %.2f", b.c0[0], b.c0[1], b.c0[2]), cp) << "! Moment force in z direction about a point" << "\n";	
		out << NemohField("0", cp) << "! Number of lines of additional information" << "\n";
	}
	out << NemohHeader("Load cases to be solved") << "\n";
	out << NemohField(Format("%d %f %f", _nf, _minf, _maxf), cp) << "! Number of wave frequencies, Min, and Max (rad/s)" << "\n";
	out << NemohField(Format("%d %f %f", Nh, minH, maxH), cp) << "! Number of wave directions, Min and Max (degrees)" << "\n";
	
	out << NemohHeader("Post processing") << "\n";
	out << NemohField(Format("%4<d %.2f %.2f", irf ? 1 : 0, irfStep, irfDuration), cp) << "! IRF                    ! IRF calculation (0 for no calculation), time step and duration" << "\n";
	out << NemohField(Format("%d", showPressure ? 1 : 0), cp) << "! Show pressure" << "\n";	
	out << NemohField(Format("%4<d %.2f %.2f", nKochin, minK, maxK), cp) << "! Kochin function        ! Number of directions of calculation (0 for no calculations), Min and Max (degrees)" << "\n";
	out << NemohField(Format("%4<d %4<d %.2f %.2f", nFreeX, nFreeY, domainX, domainY), cp) << "! Free surface elevation ! Number of points in x direction (0 for no calculations) and y direction and dimensions of domain in x and y direction" << "\n";
	
	out << "---";
}

bool Nemoh::Load_Inf(String fileName) {
	if (hd().Nb != 1)
		throw Exc(Format(t_("SeaFEM_Nemoh only allows one body, found %d"), hd().Nb));

	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
				
	hd().cg.setConstant(3, 1, Null);
	hd().cb.setConstant(3, 1, Null);
	hd().Vo.SetCount(1, Null);
	hd().C.SetCount(1);
	hd().C[0].setConstant(6, 6, Null);   
	
	double minimumDirectionAngle = 0;
	
	String line;
	while(!in.IsEof()) {
		line = in.GetLine();
		int pos;
		if ((pos = line.FindAfter("XG [m]=")) >= 0) 
			hd().cg(0, 0) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("YG [m]=")) >= 0) 
			hd().cg(1, 0) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("ZG [m]=")) >= 0) 
			hd().cg(2, 0) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("XC [m]=")) >= 0) 
			hd().cb(0, 0) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("YC [m]=")) >= 0) 
			hd().cb(1, 0) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("ZC [m]=")) >= 0) 
			hd().cb(2, 0) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("Displacement [m3]=")) >= 0) 
			hd().Vo[0] = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [3][3] [N/m]=")) >= 0) 
			hd().C[0](2, 2) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [3][4] [N/rad]=")) >= 0) 
			hd().C[0](2, 3) = hd().C[0](3, 2) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [3][5] [N/rad]=")) >= 0) 
			hd().C[0](2, 4) = hd().C[0](4, 2) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [4][4] [Nm/rad]=")) >= 0) 
			hd().C[0](3, 3) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [4][5] [Nm/rad]=")) >= 0) 
			hd().C[0](3, 4) = hd().C[0](4, 3) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [4][6] [Nm/rad]=")) >= 0) 
			hd().C[0](3, 5) = hd().C[0](5, 3) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [5][5] [Nm/rad]=")) >= 0) 
			hd().C[0](4, 4) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("K [5][6] [Nm/rad]=")) >= 0) 
			hd().C[0](4, 5) = hd().C[0](5, 4) = ScanDouble(line.Mid(pos));
		else if ((pos = line.FindAfter("Minimum direction angle =")) >= 0)
			minimumDirectionAngle = ScanDouble(line.Mid(pos))*180/M_PI;
	}
	
	for (int ih = 0; ih < hd().head.size(); ++ih)
		hd().head[ih] -= minimumDirectionAngle;

	return true;	
}
	
bool Nemoh::Load_Hydrostatics() {
	hd().cg.setConstant(3, hd().Nb, Null);
	hd().cb.setConstant(3, hd().Nb, Null);
	hd().Vo.SetCount(hd().Nb, Null);
	String line;
	
	for (int b = 0; b < hd().Nb; ++b) {
	    String fileHydro;
	    if (hd().Nb == 1)
	        fileHydro = AppendFileName(folder, AppendFileName("Mesh", "Hydrostatics.dat"));
	    else
	        fileHydro = AppendFileName(folder, AppendFileName("Mesh", Format("Hydrostatics_%d.dat", b)));
	    
	    FileInLine in(fileHydro);
	    if (!in.IsOpen())
	        return false;
	    
	    FieldSplit f(in);
	    f.IsSeparator = IsTabSpace;
	    for (int i = 0; i < 3 && !in.IsEof(); ++i) {
			f.Load(in.GetLine());
			hd().cg(i, b) = f.GetDouble(6);
			hd().cb(i, b) = f.GetDouble(2);
	    }
		f.Load(in.GetLine());
	    hd().Vo[b] = f.GetDouble(2); 		
	}
	return true;
}

bool Nemoh::Load_KH() {
	hd().C.SetCount(hd().Nb);
	for (int ib = 0; ib < hd().Nb; ++ib) {
	    String fileKH;
		if (hd().Nb == 1) 
			fileKH = AppendFileName(folder, AppendFileName("Mesh", "KH.dat"));
		else 
			fileKH = AppendFileName(folder, AppendFileName("Mesh", Format("KH_%d.dat", ib)));
	    
	    FileInLine in(fileKH);
		if (!in.IsOpen()) 
	        return false;

		hd().C[ib].setConstant(6, 6, 0);
	    
	    FieldSplit f(in);
	    f.IsSeparator = IsTabSpace;
		for (int i = 0; i < 6 && !in.IsEof(); ++i) {
			f.Load(in.GetLine());
			for (int ifr = 0; ifr < 6; ++ifr)
				hd().C[ib](i, ifr) = f.GetDouble(ifr);
		}
	}
	return true;
}

bool Nemoh::Load_Radiation(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	String line;
	FieldSplit f(in);
	f.IsSeparator = IsTabSpace;
	in.GetLine();
	while(!in.IsEof()) {
		line = in.GetLine();
	    if (line.Find("Motion of body") >= 0 || line.Find("dof_") >= 0)
	        break;
	}
	hd().A.SetCount(6*hd().Nb);
	hd().B.SetCount(6*hd().Nb);
	for (int i = 0; i < 6*hd().Nb; ++i) {
		hd().A[i].SetCount(6*hd().Nb);
		hd().B[i].SetCount(6*hd().Nb);
		for (int j = 0; j < 6*hd().Nb; ++j) {
			hd().A[i][j].setConstant(hd().Nf, Null);	
			hd().B[i][j].setConstant(hd().Nf, Null);	
		}
	}
	for (int ibody = 0; ibody < hd().Nb; ++ibody) {
		for (int idof = 0; idof < 6; ++idof) {
			if (datacal.IsDof(ibody, idof)) {
				for (int ifr = 0; ifr < hd().Nf; ++ifr) {	
					f.Load(in.GetLine());
					int col = 1;
					for (int idof2 = 0; idof2 < 6; ++idof2) {			
						if (datacal.IsDof(ibody, idof2)) {
							hd().A[ibody*6 + idof][ibody*6 + idof2][ifr] = f.GetDouble(col++);
		        			hd().B[ibody*6 + idof][ibody*6 + idof2][ifr] = f.GetDouble(col++);
			        	}
					}
				}
		    	if (idof < 6)
		    		in.GetLine();
	    	}
		}
	}
	return true;
}

bool Nemoh::Load_Excitation(String folder) {	
	return Load_Forces(hd().ex, folder, "ExcitationForce.tec");
}

bool Nemoh::Load_Diffraction(String folder) {
	return Load_Forces(hd().sc, folder, "DiffractionForce.tec");
}

bool Nemoh::Load_FroudeKrylov(String folder) {
	return Load_Forces(hd().fk, folder, "FKForce.tec");
}

bool Nemoh::Load_Forces(Hydro::Forces &fc, String nfolder, String fileName) {
	FileInLine in(AppendFileName(nfolder, AppendFileName("Results", fileName)));
	if (!in.IsOpen())
		return false;
	String line;
	FieldSplit f(in);
	f.IsSeparator = IsTabSpace;
	in.GetLine();
	Vector<Vector<int>> dof;
	dof.SetCount(hd().Nb);
	while(!in.IsEof()) {
		line = in.GetLine();
		if (line.StartsWith("Zone") || line.StartsWith("angle"))
	        break;
	}
	hd().Initialize_Forces(fc);
	for (int ih = 0; ih < hd().Nh; ++ih) {
		int ifr = 0;
		while(!in.IsEof()) {
			line = in.GetLine();
			if (line.StartsWith("Zone") || line.StartsWith("angle"))
				break;
			f.Load(line);
			int il = 0;
			for (int ib = 0; ib < hd().Nb; ++ib) {
				for (int ibdof = 0; ibdof < 6; ++ibdof) {
					if (datacal.IsDof(ib, ibdof)) {
						if (ifr >= hd().Nf)
							throw Exc(in.Str() + "\n"  + t_("Number of frequencies higher than the defined in Nemoh.cal file"));		
						//if (ib >= hd().Nb)
						//	throw Exc(in.Str() + "\n"  + t_("Number of bodies higher than the defined in Nemoh.cal file"));		
						double ma = fc.ma[ih](ifr, ib*6+ibdof) = f.GetDouble(1 + 2*il);	
						double ph = fc.ph[ih](ifr, ib*6+ibdof) = -f.GetDouble(1 + 2*il + 1); //-Phase to follow Wamit
						fc.re[ih](ifr, ibdof) = ma*cos(ph); 
						fc.im[ih](ifr, ibdof) = ma*sin(ph);
						il++; 
					}
				}
			}
			ifr++;
		}
	}
	return true;
}

bool Nemoh::Load_IRF(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	String line;
	FieldSplit f(in);	
	f.IsSeparator = IsTabSpace;
	hd().Awinf.setConstant(hd().Nb*6, hd().Nb*6, Null);
	//int ibodydof = 0;
	hd().Kirf.SetCount(hd().Nb*6); 	// Initialize Kirf		
    for (int i = 0; i < hd().Nb*6; ++i) {
    	hd().Kirf[i].SetCount(hd().Nb*6); 			 
   		for (int j = 0; j < hd().Nb*6; ++j)
			hd().Kirf[i][j].setConstant(hd().Tirf.size(), Null);
    }
    while(!in.IsEof()) {
		line = in.GetLine();	
		if (line.Find("Zone t=") >= 0) 
			break;
	}
	for (int ib = 0; ib < hd().Nb; ++ib) {
		for (int ibdof = 0; ibdof < 6; ++ibdof) {
			if (datacal.IsDof(ib, ibdof)) {
				for (int iNt = 0; iNt < hd().Tirf.size(); ++iNt) {
					f.Load(in.GetLine());
					int col = 1;
					for (int idf = 0; idf < 6; ++idf) {
						if (datacal.IsDof(ib, idf)) {
							hd().Awinf(ibdof, idf) = f.GetDouble(col++);
							hd().Kirf[ib*6+ibdof][ib*6+idf][iNt] = f.GetDouble(col++);
						}
					}
				}
				in.GetLine();
			}
		}
	}
	return true;
}

void Nemoh::Save(String ) {
	throw Exc("Option not implemented");
}		

