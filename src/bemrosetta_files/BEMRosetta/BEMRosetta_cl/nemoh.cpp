#include "BEMRosetta.h"

bool Nemoh::Load(String file, double) {
	hd().file = file;
	hd().name = GetFileTitle(GetFileFolder(file));
	folder = GetFileFolder(file);
	hd().len = 1;
	hd().dimen = true;
	hd().Nb = Null;
	
	String ext = GetFileExt(file); 
	if (ext == ".cal")
		hd().code = Hydro::NEMOH;
	else
		hd().code = Hydro::SEAFEM_NEMOH;
	
	try {
		String fileCal;
		BEMData::Print("\n\n" + Format(t_("Loading '%s'"), file));
		if (hd().code == Hydro::NEMOH) 
			fileCal = file;
		else 
			fileCal = AppendFileName(folder, "Nemoh_output/Nemoh.cal");
		if (!Load_Cal(fileCal)) 
			throw Exc("\n" + Format(t_("File '%s' not found"), fileCal));
		
		String fileRad, folderForces;
		if (hd().code == Hydro::NEMOH) {
			BEMData::Print(x_("\n- ") + t_("Hydrostatics file(s) 'Mesh/Hydrostatics*.dat'"));
			if (!Load_Hydrostatics())
				BEMData::PrintWarning(x_(": **") + t_("Not found") + "**");
			BEMData::Print(x_("\n- ") + t_("KH file(s) 'Mesh/KH*.dat'"));
			if (!Load_KH())
				BEMData::PrintWarning(x_(": **") + t_("Not found") + "**");
			fileRad = AppendFileName(folder, AppendFileName("Results", "RadiationCoefficients.tec"));
			folderForces = folder;
		} else {
			if (!Load_Inf(file)) 
				throw Exc("\n" + Format(t_("File '%s' not found"), file));

			fileRad = AppendFileName(folder, AppendFileName("Nemoh_output/Results", "RadiationCoefficients.tec"));
			folderForces = AppendFileName(folder, "Nemoh_output");
		} 
		
		BEMData::Print(x_("\n- ") + t_("Radiation file 'RadiationCoefficients.tec'"));
		if (!Load_Radiation(fileRad))
			BEMData::PrintWarning(x_(": **") + t_("Not found") + "**");
		
		BEMData::Print(x_("\n- ") + t_("Excitation force file 'ExcitationForce.tec'"));
		if (!Load_Excitation(folderForces))
			BEMData::PrintWarning(x_(": **") + t_("Not found") + "**");
		
		BEMData::Print(x_("\n- ") + t_("Diffraction force file 'DiffractionForce.tec'"));
		if (!Load_Diffraction(folderForces))
			BEMData::PrintWarning(x_(": **") + t_("Not found") + "**");
		BEMData::Print(x_("\n- ") + t_("Froude Krylov file 'FKForce.tec'"));
		if (!Load_FroudeKrylov(folderForces))
			BEMData::PrintWarning(x_(": **") + t_("Not found") + "**");
		
		if (hd().code == Hydro::NEMOH) {
			if (!hd().dof.IsEmpty()) {
				BEMData::Print(x_("\n- ") + t_("IRF file(s) 'IRF.tec'"));
				if (!Load_IRF(AppendFileName(folder, AppendFileName("Results", "IRF.tec"))))
					BEMData::PrintWarning(x_(": **") + t_("Not found") + "**");
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
	NemohCal data;
	
	if (!data.Load(fileName))
		return false;

	hd().rho = data.rho;
	hd().g = data.g;
	hd().h = data.h;
	if (hd().h == 0)
		hd().h = -1;
	hd().Nb = data.bodies.GetCount();
	for (int i = 0; i < hd().Nb; ++i) 	
		hd().names << GetFileTitle(data.bodies[i].meshFile);
	hd().Nf = data.Nf;
	LinSpaced(hd().w, hd().Nf, data.minF, data.maxF); 
 	hd().T.SetCount(hd().Nf);
    for (int i = 0; i < hd().Nf; ++i) {
		hd().T[i] = 2*M_PI/hd().w[i];  
    }
	hd().Nh = data.Nh;  						
    LinSpaced(hd().head, hd().Nh, data.minD, data.maxD); 		

	hd().dataFromW = true;
	
	return true;
}

int NemohCal::GetNumArgs(const FieldSplit &f) {
	for (int i = 0; i < f.GetCount(); ++i) {
		double num = ScanDouble(f.GetText(i));
		if (IsNull(num))
			return i;
	}
	return f.GetCount();
}
	
void NemohCal::LoadFreeSurface(const FileInLine &in, const FieldSplit &f) {
	nFreeX = f.GetInt(0);	nFreeY = f.GetInt(1);	domainX = f.GetDouble(2);	domainY = f.GetDouble(3);
	if (nFreeX < 0)
		throw Exc(Format(t_("[%d] Incorrect number of points in x direction %s"), in.GetLineNumber(), f.GetText(0)));
	if (nFreeX > 0 && (nFreeY <= 0 || domainX < 0 || domainY < 0))
		throw Exc(Format(t_("[%d] Incorrect free surface elevation %s"), in.GetLineNumber(), f.GetText()));	
}


void NemohCal::LoadKochin(const FileInLine &in, const FieldSplit &f) {
	nKochin = f.GetInt(0);	minK = f.GetDouble(1);	maxK = f.GetDouble(2);
	if (nKochin < 0)
		throw Exc(Format(t_("[%d] Incorrect number of Kochin function directions %s"), in.GetLineNumber(), f.GetText(0)));
	if (nKochin > 0) {
		if (minK < -360)
			throw Exc(Format(t_("[%d] Incorrect Kochin direction %s"), in.GetLineNumber(), f.GetText(1)));
		if (maxK > 360)
			throw Exc(Format(t_("[%d] Incorrect Kochin direction %s"), in.GetLineNumber(), f.GetText(2)));
		if (maxK <= minK)
			throw Exc(Format(t_("[%d] Minimum Kochin direction %s has to be lower than maximum direction %s"), in.GetLineNumber(), f.GetText(1), f.GetText(2)));	
	}
}

bool NemohCal::Load(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	String line;
	FieldSplit f(in);
	
	in.GetLine();
	f.Load(in.GetLine());	rho  = f.GetDouble(0);
	if (rho < 0 || rho > 10000)
		throw Exc(Format(t_("[%d] Incorrect rho %s"), in.GetLineNumber(), f.GetText(0)));
	f.Load(in.GetLine());	g    = f.GetDouble(0);
	if (g < 0 || g > 100)
		throw Exc(Format(t_("[%d] Incorrect g %s"), in.GetLineNumber(), f.GetText(0)));
	f.Load(in.GetLine());	h    = f.GetDouble(0);
	if (h < 0 || h > 100000)
		throw Exc(Format(t_("[%d] Incorrect depth %s"), in.GetLineNumber(), f.GetText(0)));
	f.Load(in.GetLine());	xeff = f.GetDouble(0);	yeff = f.GetDouble(1);
	in.GetLine();
	f.Load(in.GetLine());	int Nb = f.GetInt(0);
	if (Nb < 1 || Nb > 100)
		throw Exc(Format(t_("[%d] Incorrect number of bodies %s"), in.GetLineNumber(), f.GetText(0)));
	bodies.SetCount(Nb);
	for (int ib = 0; ib < Nb; ++ib) {
		NemohBody &body = bodies[ib];
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
		}
		
		if (body.npoints < 1 || body.npoints > 100000000)
			throw Exc(Format(t_("[%d] Incorrect number of points %s"), in.GetLineNumber(), f.GetText(0)));
		if (body.npanels < 1 || body.npanels > 100000000)
			throw Exc(Format(t_("[%d] Incorrect number of panels %s"), in.GetLineNumber(), f.GetText(1)));	
		f.Load(in.GetLine());	int ndof = f.GetInt(0);
		if (ndof < 0 || ndof > 6)
			throw Exc(Format(t_("[%d] Incorrect DOF %s in body %d"), in.GetLineNumber(), f.GetText(0), ib+1));
		for (int idof = 0; idof < ndof; ++idof) {
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
					body.surge = true;
				else if (y)
					body.sway = true;
				else if (z)
					body.heave = true;
			} else if (type == 2) {
				body.cx = cx;
				body.cy = cy;
				body.cz = cz;
				if (x) 
					body.roll = true;
				else if (y)
					body.pitch = true;
				else if (z)
					body.yaw = true;
			} else
				throw Exc(Format(t_("[%d] Incorrect DOF type %d set in body %d"), in.GetLineNumber(), f.GetText(0), ib+1));
		}
		f.Load(in.GetLine());	int nforces = f.GetInt(0);
		in.GetLine(nforces);	// Discarded
		f.Load(in.GetLine());	int nadditional = f.GetInt(0);	
		in.GetLine(nadditional);// Discarded	
	}
	in.GetLine();
	f.Load(in.GetLine());	Nf = f.GetInt(0);	minF = f.GetDouble(1);	maxF = f.GetDouble(2);
	if (Nf < 1 || Nf > 1000)
		throw Exc(Format(t_("[%d] Incorrect number of frequencies %s"), in.GetLineNumber(), f.GetText(0)));
	if (minF <= 0)
		throw Exc(Format(t_("[%d] Incorrect frequency %s"), in.GetLineNumber(), f.GetText(1)));
	if (maxF <= minF)
		throw Exc(Format(t_("[%d] Minimum frequency %s has to be lower than maximum frequency %s"), in.GetLineNumber(), f.GetText(1), f.GetText(2)));	
	
	f.Load(in.GetLine());	Nh = f.GetInt(0);	minD = f.GetDouble(1);	maxD = f.GetDouble(2);
	if (Nh < 1 || Nh > 1000)
		throw Exc(Format(t_("[%d] Incorrect number of directions %s"), in.GetLineNumber(), f.GetText(0)));
	if (minD < -180)
		throw Exc(Format(t_("[%d] Incorrect direction %s"), in.GetLineNumber(), f.GetText(1)));
	if (maxD > 180)
		throw Exc(Format(t_("[%d] Incorrect direction %s"), in.GetLineNumber(), f.GetText(2)));
	if (maxD < minD)
		throw Exc(Format(t_("[%d] Minimum direction %s has to be lower than maximum direction %s"), in.GetLineNumber(), f.GetText(1), f.GetText(2)));	
	
	in.GetLine();
	f.Load(in.GetLine());	irf = f.GetInt(0) > 0;	irfStep = f.GetDouble(1);	irfDuration = f.GetDouble(2);
	if (irfStep <= 0)
		throw Exc(Format(t_("[%d] Incorrect IRF step %s"), in.GetLineNumber(), f.GetText(1)));
	if (irfDuration <= irfStep)
		throw Exc(Format(t_("[%d] IRF step %s has to be lower than duration %s"), in.GetLineNumber(), f.GetText(1), f.GetText(2)));	
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
		throw Exc(Format(t_("[%d] Unexpected data %s"), in.GetLineNumber(), f.GetText()));			
	
	f.Load(in.GetLine());	
	if (GetNumArgs(f) == 4) {
		if (!loadedFree)
			LoadFreeSurface(in, f);
		else
			throw Exc(Format(t_("[%d] Free surface data already loaded %s"), in.GetLineNumber(), f.GetText()));				
	} else if (GetNumArgs(f) == 3) {
		if (!loadedKochin)
			LoadKochin(in, f);
		else
			throw Exc(Format(t_("[%d] Kochin data already loaded %s"), in.GetLineNumber(), f.GetText()));				
	} else if (GetNumArgs(f) > 0)
		throw Exc(Format(t_("[%d] Unexpected data %s"), in.GetLineNumber(), f.GetText()));			
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

Vector<String> NemohCal::Check() {
	Vector<String> ret;
	
	if (rho < 0 || rho > 10000)
		 ret << Format(t_("Incorrect rho %s"), FormatDoubleEmpty(rho));
	if (g < 0 || g > 100)
		ret << Format(t_("Incorrect g %s"), FormatDoubleEmpty(g));
	if (h < 0 || h > 100000)
		ret << Format(t_("Incorrect depth %s"), FormatDoubleEmpty(h));

	if (Nf < 1 || Nf > 1000)
		ret << Format(t_("Incorrect number of frequencies %s"), FormatIntEmpty(Nf));
	if (minF <= 0)
		ret << Format(t_("Incorrect frequency %s"), FormatDoubleEmpty(minF));
	if (maxF <= minF)
		ret << Format(t_("Minimum frequency %s has to be lower than maximum frequency %s"), FormatDoubleEmpty(minF), FormatDoubleEmpty(maxF));	
	
	if (Nh < 1 || Nh > 1000)
		ret << Format(t_("Incorrect number of directions %s"), FormatIntEmpty(Nh));
	if (minD < -180)
		ret << Format(t_("Incorrect direction %s"), FormatDoubleEmpty(minD));
	if (maxD > 180)
		ret << Format(t_("Incorrect direction %s"), FormatDoubleEmpty(maxD));
	if (maxD < minD)
		ret << Format(t_("Minimum direction %s has to be lower than maximum direction %s"), FormatDoubleEmpty(minD), FormatDoubleEmpty(maxD));	
	
	if (irf) {
		if (irfStep <= 0)
			ret << Format(t_("Incorrect IRF step %s"), FormatDoubleEmpty(irfStep));
		if (irfDuration <= irfStep)
			ret << Format(t_("IRF step %s has to be lower than duration %s"), FormatDoubleEmpty(irfStep), FormatDoubleEmpty(irfDuration));	
	}
	
	if (nFreeX < 0)
		ret << Format(t_("Incorrect number of points in x direction %s (0 for no free surface calculation)"), FormatIntEmpty(nFreeX));
	if (nFreeX > 0) {
		if (nFreeY <= 0)
			ret << Format(t_("Incorrect number of points in x direction %s"), FormatIntEmpty(nFreeY));
		if (domainX <= 0)
			ret << Format(t_("Incorrect free surface domain X %s"), FormatDoubleEmpty(domainX));
		if (domainY <= 0)
			ret << Format(t_("Incorrect free surface domain Y %s"), FormatDoubleEmpty(domainY));
	}
	
	if (nKochin < 0)
		ret << Format(t_("Incorrect number of Kochin function directions %s"), FormatIntEmpty(nKochin));
	if (nKochin > 0) {
		if (minK < -180)
			ret << Format(t_("Incorrect Kochin direction %s"), FormatDoubleEmpty(minK));
		if (maxK < 180)
			ret << Format(t_("Incorrect Kochin direction %s"),FormatDoubleEmpty(minK));
		if (maxK <= minK)
			Format(t_("Minimum Kochin direction %s has to be lower than maximum direction %s"), FormatDoubleEmpty(minK), FormatDoubleEmpty(maxK));	
	}	
	
	return ret;
}

void NemohCal::CreateId(String folder) const {
	String fileName = AppendFileName(folder, "ID.dat");
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to create '%s'"), fileName));
	out << "1\n.";
}

void NemohCal::CreateBat(String folder, bool bin, String preName, String solvName, String postName) const {
	String fileName = AppendFileName(folder, "Nemoh.bat");
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to create '%s'"), fileName));
	String strBin;
	if (bin)
		strBin = "\"bin";
	out << AppendFileName(strBin, preName) << "\"\n"
		<< AppendFileName(strBin, solvName) << "\"\n"
		<< AppendFileName(strBin, postName) << "\"";
}

void NemohCal::CreateInput(String folder) const {
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
	return ret;
}

int NemohBody::GetNDOF() const {
	int ret = 0;
	if (surge)
		ret++;
	if (sway)
		ret++;	
	if (heave)
		ret++;
	if (roll)
		ret++;
	if (pitch)
		ret++;
	if (yaw)
		ret++;
	return ret;
}

void NemohCal::SaveFolder(String folder, bool bin, const BEMData &bem) const {
	CreateId(folder);
	CreateInput(folder);
	String folderMesh = AppendFileName(folder, "mesh");
	DeleteFolderDeep(folderMesh);
	if (!DirectoryCreate(folderMesh)) 
		throw Exc(Format(t_("Problem creating %s folder"), folderMesh));

	for (int i = 0; i < bodies.GetCount(); ++i) {
		String name = GetFileName(bodies[i].meshFile);
		String dest = AppendFileName(folderMesh, name);
		if (!FileCopy(bodies[i].meshFile, dest)) 
			throw Exc(Format(t_("Problem copying mesh file '%s'"), bodies[i].meshFile));
	}
	Save_Cal(folder);
	
	String folderResults = AppendFileName(folder, "results");
	DeleteFolderDeep(folderResults);
	if (!DirectoryCreate(folderResults)) 
		throw Exc(Format(t_("Problem creating '%s' folder"), folderResults));
	
	String preName = "preprocessor.exe";
	String solvName = "solver.exe";
	String postName = "postprocessor.exe";
	if (bin) {
		String binResults = AppendFileName(folder, "bin");
		DeleteFolderDeep(binResults);
		if (!DirectoryCreate(binResults)) 
			throw Exc(Format(t_("Problem creating '%s' folder"), binResults));
		
		preName = GetFileName(bem.nemohPathPreprocessor);
		String destProprocessor = AppendFileName(binResults, preName);
		if (!FileCopy(bem.nemohPathPreprocessor, destProprocessor)) 
			throw Exc(Format(t_("Problem copying preprocessor file '%s'"), bem.nemohPathPreprocessor));		
		solvName = GetFileName(bem.nemohPathSolver);
		String destSolver = AppendFileName(binResults, solvName);
		if (!FileCopy(bem.nemohPathSolver, destSolver)) 
			throw Exc(Format(t_("Problem copying solver file '%s'"), bem.nemohPathSolver));		
		postName = GetFileName(bem.nemohPathPostprocessor);
		String destPostprocessor = AppendFileName(binResults, postName);
		if (!FileCopy(bem.nemohPathPostprocessor, destPostprocessor)) 
			throw Exc(Format(t_("Problem copying postprocessor file '%s'"), bem.nemohPathPostprocessor));		
		if (!GetFileName(bem.nemohPathGREN).IsEmpty()) {
			String destGREN = AppendFileName(folder, GetFileName(bem.nemohPathGREN));
			if (!FileCopy(bem.nemohPathGREN, destGREN)) 
				throw Exc(Format(t_("Problem copying gren file '%s'"), bem.nemohPathGREN));
		}
	}
	CreateBat(folder, bin, preName, solvName, postName);
}

void NemohCal::Save_Cal(String folder) const {
	String fileName = AppendFileName(folder, "Nemoh.cal");
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'"), fileName));
		
	out << NemohHeader("Environment - Created with BEMRosetta") << "\n";
	out << NemohField(Format("%f", rho), 30) 		   << "! RHO             ! KG/M**3   ! Fluid specific volume" << "\n";
	out << NemohField(Format("%f", g), 30)   		   << "! G               ! M/S**2    ! Gravity " << "\n";
	out << NemohField(Format("%f", h), 30)   		   << "! DEPTH           ! M         ! Water depth" << "\n";
	out << NemohField(Format("%f %f", xeff, yeff), 30) << "! XEFF YEFF       ! M         ! Wave measurement point" << "\n";
	
	out << NemohHeader("Description of floating bodies") << "\n";
	out << NemohField(Format("%d", bodies.GetCount()), 30) << "! Number of bodies" << "\n";
	
	for (int i = 0; i < bodies.GetCount(); ++i) {
		const NemohBody &b = bodies[i];
		out << NemohHeader(Format("Body %d", i+1)) << "\n";	
		String file = AppendFileName("mesh", GetFileName(b.meshFile));
		
		out << NemohField(Format("%s", file), 30) << "! Name of mesh file" << "\n";
		out << NemohField(Format("%d %d", b.npoints, b.npanels), 30) << "! Number of points and number of panels" << "\n";	
		out << NemohField(Format("%d", b.GetNDOF()), 30) << "! Number of degrees of freedom" << "\n";	
		if (b.surge)
			out << NemohField("1 1. 0. 0. 0. 0. 0.", 30) << "! Surge" << "\n";	
		if (b.sway)
			out << NemohField("1 0. 1. 0. 0. 0. 0.", 30) << "! Sway" << "\n";	
		if (b.heave)
			out << NemohField("1 0. 0. 1. 0. 0. 0.", 30) << "! Heave" << "\n";	
		if (b.roll)
			out << NemohField(Format("2 1. 0. 0. %.2f %.2f %.2f", b.cx, b.cy, b.cz), 30) << "! Roll about a point" << "\n";	
		if (b.pitch)
			out << NemohField(Format("2 0. 1. 0. %.2f %.2f %.2f", b.cx, b.cy, b.cz), 30) << "! Pitch about a point" << "\n";	
		if (b.yaw)		
			out << NemohField(Format("2 0. 0. 1. %.2f %.2f %.2f", b.cx, b.cy, b.cz), 30) << "! Yaw about a point" << "\n";	
		out << NemohField(Format("%d", b.GetNDOF()), 30) << "! Number of resulting generalised forces" << "\n";	
		if (b.surge)
			out << NemohField("1 1. 0. 0. 0. 0. 0.", 30) << "! Force in x direction" << "\n";	
		if (b.sway)
			out << NemohField("1 0. 1. 0. 0. 0. 0.", 30) << "! Force in y direction" << "\n";	
		if (b.heave)
			out << NemohField("1 0. 0. 1. 0. 0. 0.", 30) << "! Force in z direction" << "\n";	
		if (b.roll)
			out << NemohField(Format("2 1. 0. 0. %.2f %.2f %.2f", b.cx, b.cy, b.cz), 30) << "! Moment force in x direction about a point" << "\n";	
		if (b.pitch)
			out << NemohField(Format("2 0. 1. 0. %.2f %.2f %.2f", b.cx, b.cy, b.cz), 30) << "! Moment force in y direction about a point" << "\n";	
		if (b.yaw)		
			out << NemohField(Format("2 0. 0. 1. %.2f %.2f %.2f", b.cx, b.cy, b.cz), 30) << "! Moment force in z direction about a point" << "\n";	
		out << NemohField("0", 30) << "! Number of lines of additional information" << "\n";
	}
	out << NemohHeader("Load cases to be solved") << "\n";
	out << NemohField(Format("%d %f %f", Nf, minF, maxF), 30) << "! Number of wave frequencies, Min, and Max (rad/s)" << "\n";
	out << NemohField(Format("%d %f %f", Nh, minD, maxD), 30) << "! Number of wave directions, Min and Max (degrees)" << "\n";
	
	out << NemohHeader("Post processing") << "\n";
	out << NemohField(Format("%4<d %.2f %.2f", irf ? 1 : 0, irfStep, irfDuration), 30) << "! IRF                    ! IRF calculation (0 for no calculation), time step and duration" << "\n";
	out << NemohField(Format("%d", showPressure ? 1 : 0), 30) << "! Show pressure" << "\n";	
	out << NemohField(Format("%4<d %.2f %.2f", nKochin, minK, maxK), 30) << "! Kochin function        ! Number of directions of calculation (0 for no calculations), Min and Max (degrees)" << "\n";
	out << NemohField(Format("%4<d %4<d %.2f %.2f", nFreeX, nFreeY, domainX, domainY), 30) << "! Free surface elevation ! Number of points in x direction (0 for no calcutions) and y direction and dimensions of domain in x and y direction" << "\n";
	
	out << "---";
}

bool Nemoh::Load_Inf(String fileName) {
	if (hd().Nb != 1)
		throw Exc(Format(t_("SeaFEM_Nemoh only allows one body, found %d"), hd().Nb));
		
	hd().cg.setConstant(3, 1, Null);
	hd().cb.setConstant(3, 1, Null);
	hd().Vo.SetCount(1, Null);
	hd().C.SetCount(1);
	hd().C[0].setConstant(6, 6, 0);   
	
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
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
	
	for (int ih = 0; ih < hd().head.GetCount(); ++ih)
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
	    
		hd().C[ib].setConstant(6, 6, 0);    
	    FileInLine in(fileKH);
	    if (!in.IsOpen()) 
	        return false;
	    
	    FieldSplit f(in);
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
	in.GetLine();
	Vector<int> dof;
	dof.SetCount(hd().Nb, 0);
	while(!in.IsEof()) {
		line = in.GetLine();
	    if (line.Find("Motion of body") >= 0)
	        break;
	    f.Load(line);
	    int ibody = f.GetInt(1) - 1;
	    int ndof = f.GetInt(2);
		dof[ibody] = max(dof[ibody], ndof);    
	}
	if (hd().dof.IsEmpty())
		hd().dof = pick(dof);
	else if (!IsEqualRange(dof, hd().dof))
		throw Format(t_("DOF does not match in %s"), fileName);
	hd().A.SetCount(hd().Nf);
	hd().B.SetCount(hd().Nf);
	for (int k = 0; k < hd().Nf; ++k) {
		hd().A[k].setConstant(hd().Nb*6, hd().Nb*6, Null);
		hd().B[k].setConstant(hd().Nb*6, hd().Nb*6, Null);
	}
	int ibodydof = 0;
	for (int ibody = 0; ibody < hd().Nb; ++ibody) {
		for (int idof = 0; idof < hd().dof[ibody]; ++idof) {
			for (int k = 0; k < hd().Nf; ++k) {	
				f.Load(in.GetLine());
				for (int df = 0; df < hd().dof[ibody]; ++df) {		
					hd().A[k](ibodydof, df) = f.GetDouble(1 + 2*df);
	        		hd().B[k](ibodydof, df) = f.GetDouble(2 + 2*df);
				}
	        }
	        ++ibodydof;
	        in.GetLine();
	    }
	}
	return true;
}

bool Nemoh::Load_Excitation(String folder) {	
	return Load_Forces(hd().ex, folder, "ExcitationForce.tec", "Diffraction force");
}

bool Nemoh::Load_Diffraction(String folder) {
	return Load_Forces(hd().sc, folder, "DiffractionForce.tec", "Diffraction force");
}

bool Nemoh::Load_FroudeKrylov(String folder) {
	return Load_Forces(hd().fk, folder, "FKForce.tec", "FKforce");
}

bool Nemoh::Load_Forces(Hydro::Forces &fc, String nfolder, String fileName, String textDelim) {
	FileInLine in(AppendFileName(nfolder, AppendFileName("Results", fileName)));
	if (!in.IsOpen())
		return false;
	String line;
	FieldSplit f(in);
	in.GetLine();
	Vector<int> dof;
	dof.SetCount(hd().Nb, 0);
	while(!in.IsEof()) {
		line = in.GetLine();
	    if (line.Find(textDelim) >= 0)
	        break;
	    f.Load(line);
	    int ibody = f.GetInt(1) - 1;
	    int ndof = f.GetInt(2);
		dof[ibody] = max(dof[ibody], ndof);    
	}
	if (hd().dof.IsEmpty())
		hd().dof = pick(dof);
	else if (!IsEqualRange(dof, hd().dof))
		throw Format(t_("[%d] DOF does not match in %s"), in.GetLineNumber(), fileName);
	hd().Initialize_Forces(fc);
	for (int h = 0; h < hd().Nh; ++h) {
		int ifr = 0;
		while(!in.IsEof()) {
			line = in.GetLine();
			if (line.Find(textDelim) >= 0)
				break;
			f.Load(line);
			int ib = 0, idof = 0, ibdof = 0;
			for (int i = 0; i < hd().Nb*6; ++i) {
				if (ifr >= hd().Nf)
					throw Exc(Format(t_("[%d] Number of frequencies higher than the defined in Nemoh.cal file"), in.GetLineNumber()));		
				if (ib >= hd().Nb)
					throw Exc(Format(t_("[%d] Number of bodies higher than the defined in Nemoh.cal file"), in.GetLineNumber()));		
				double ma = fc.ma[h](ifr, ibdof) = f.GetDouble(1 + 2*i);	
				double ph = fc.ph[h](ifr, ibdof) = -f.GetDouble(1 + 2*i + 1); //-Phase to follow Wamit
				fc.re[h](ifr, ibdof) = ma*cos(ph); 
				fc.im[h](ifr, ibdof) = ma*sin(ph); 
				idof++;
				ibdof++;
				if (idof >= hd().dof[ib]) {
					idof = 0;
					ib++;
					ibdof = 6*ib;
					if (ib >= hd().Nb)
						break;
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
	hd().Awinf.setConstant(hd().Nb*6, hd().Nb*6, Null);
	int ibodydof = 0;
	for (int ibody = 0; ibody < hd().Nb; ++ibody) {
		for (int idof = 0; idof < hd().dof[ibody]; ++idof) {
			while(!in.IsEof()) {
				line = in.GetLine();	
				if (line.Find("Zone t=") >= 0) 
					break;
			}
			line = in.GetLine();	
			f.Load(line);
			for (int df = 0; df < hd().dof[ibody]; ++df) 
				hd().Awinf(ibodydof, df) = f.GetDouble(1 + 2*df);
			
			++ibodydof;
		}
	}
	return true;
}

void Nemoh::Save(String ) {
	throw Exc("Option not implemented");
}		

