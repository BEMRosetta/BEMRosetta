// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2023, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#include <Hdf5/hdf5.h>


void HydroStar::SaveCase(String folder, bool withPotentials, bool x0z, bool y0z, 
				const UVector<bool> &listDOF, int qtftype) const {
	if (!DirectoryCreateX(folder))
		throw Exc(Format(t_("Problem creating '%s' folder"), folder));
	
	Save_HSG(AFX(folder, "input.hsg"));
	Save_MCN(AFX(folder, "input.mcn"));
	Save_QTF(AFX(folder, "input.qtf"), qtftype);
	Save_RAO(AFX(folder, "input.rao"), listDOF, qtftype > 0);
	Save_RDF(AFX(folder, "input.rdf"));
	Save_DFT(AFX(folder, "input.dft"));
	
	for (int ib = 0; ib < dt.Nb; ++ib) {
		String dest = AFX(folder, Format("Body_%d.hst", ib+1));
		Body::SaveAs(dt.msh[ib], dest, Body::HYDROSTAR_HST, Body::UNDERWATER, dt.rho, dt.g, y0z, x0z);
	}
}

void HydroStar::Save_HSG(String fileName) const {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to save '%s'. File already used."), fileName));
	
	out << "<?xml version=\"1.0\" ?>\n"
		<< "<hydrostar_project>\n"
		<< "	<InputFiles Name=\"input.mcn\"/>\n"
		<< "	<InputFiles Name=\"input.qtf\"/>\n"
		<< "	<InputFiles Name=\"input.rao\"/>\n"
		<< "	<InputFiles Name=\"input.rdf\"/>\n";
	for (int ib = 0; ib < dt.Nb; ++ib)
		out << Format("	<InputFiles Name=\"Body_%d.hst\"/>\n", ib+1);
	out	<< "	<InputFiles Name=\"input.dft\"/>\n"
	   	<< "</hydrostar_project>";
}

void HydroStar::Save_MCN(String fileName) const {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to save '%s'. File already used."), fileName));
		
	out << 	"#Diffraction results to use\n"
			"FILENAME rd1\n";
			
	out << 	"\n";
	
	out <<	"#Mass of the body (in kg)\n";
	for (int ib = 0; ib < dt.Nb; ++ib)
		out << Format("MASS_BODY  %d     %.3f\n", ib+1, dt.msh[ib].GetMass());

	out << 	"\n";
	
	out <<	"#Center of gravity (in mesh reference)\n";
	for (int ib = 0; ib < dt.Nb; ++ib)
		out << Format("COGPOINT_BODY  %d      %.3E    %.3E    %.3E\n", ib+1, dt.msh[ib].dt.cg.x, dt.msh[ib].dt.cg.y, dt.msh[ib].dt.cg.z);
		
	out << 	"\n";
	
	out <<	"#Rotational inertia\n";

	for (int ib = 0; ib < dt.Nb; ++ib) {
		double m = dt.msh[ib].GetMass();
		if (m > 0)
			out << Format("GYRADIUS_BODY  %d      %.3E    %.3E    %.3E     %.3E     %.3E  %.3E\n", ib+1, 
							sqrt(dt.msh[ib].dt.M(3, 3)/m), sqrt(dt.msh[ib].dt.M(3, 4)/m), sqrt(dt.msh[ib].dt.M(3, 5)/m),
														   sqrt(dt.msh[ib].dt.M(4, 4)/m), sqrt(dt.msh[ib].dt.M(4, 5)/m),
																						  sqrt(dt.msh[ib].dt.M(5, 5)/m));
	}
	out << 	"\n";
	
	out << "#RAIDEUR D'ANCRAGE + FREE SURFACE\n";
	for (int ib = 0; ib < dt.Nb; ++ib) {
		if (IsLoadedCMoor(ib) || IsLoadedCAdd(ib)) {
			out << Format("STIFFNESS_MATRIX TYPE 0 BODY %d %d\n", ib+1, ib+1);
			for (int idof = 0; idof < 6; ++idof)
				for (int jdof = 0; jdof < 6; ++jdof) {
					double d = 0;
					if (IsLoadedCMoor(ib, idof, jdof))
						d += CMoor_dim(ib, idof, jdof);
					if (IsLoadedCAdd(ib, idof, jdof))
						d += CAdd_dim(ib, idof, jdof);
					out << Format("%d %d %.3E\n", idof+1, jdof+1, d);
				}
			out << "ENDSTIFFNESS_MATRIX\n";
		}
	}
	
	out << 	"\n\n";
	
	for (int ib = 0; ib < dt.Nb; ++ib) {
		if (IsLoadedDlin(ib)) {
			out << "DAMPING_MATRIX  TYPE 0\n";
			for (int idof = 0; idof < 6; ++idof)
				for (int jdof = 0; jdof < 6; ++jdof) {
					double d = 0;
					if (IsLoadedDlin(ib, idof, jdof))
						d = dt.msh[ib].dt.Dlin(idof, jdof);
					out << Format("%d %d %.3E\n", 6*ib + idof+1, 6*ib + jdof+1, d);
				}
			out << "ENDDAMPING_MATRIX\n";
		}
	}
	
	out << "\n\n";	
	
	out << Format("RHO %.1f\n", rho_ndim());
	out << "\n";
	out << Format("GRAVITY %.5f\n", g_ndim());
	out << "\n";
	out << "ENDFILE";
}

void HydroStar::Save_QTF(String fileName, int qtftype) const {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to save '%s'. File already used."), fileName));
		
	out << "\n";
	
	double dw = (dt.w.Top() - dt.w[0])/(dt.w.size()-1);
	out << Format("DIFFREQUENCE   0.0   %.2f   %.3f\n", dt.w.Top(), dw);
	out << Format("WAVFREQUENCE   %.1f   %.2f    %.3f\n", dt.w[0], dt.w.Top(), dw);
	
	out << "\n";
	
	if (qtftype == 9)
		out << "TYPEFORMULE   NEAR-FIELD\n";
	else if (qtftype == 7)	
		out << "TYPEFORMULE   MIDDLE-FIELD\n"; 
	else
		throw Exc(Format("Type of QTF %d is not supported in HydroStar", qtftype));
	
	out << "\n";
	
	if (qtftype == 7) {
		out << Format("NBBOITE %d\n", dt.Nb);
		for (int ib = 0; ib < dt.Nb; ++ib) {
			Surface s = clone(dt.msh[ib].dt.under);
			s.GetSegments();
			s.GetEnvelope();
			double dz = s.GetAvgLenSegment();
			double zmin = 1.5*s.env.minZ;
			out << Format("CSFILE AUTO BODY %d %.3f  %3f  %3f\n", ib+1, zmin, 0, dz);
		}
	}
		
	out << "\n";
	
	out << "MULTIDIRECTIONELLE\n";
	
	out << "\n";
	
	out << "ENDFILE";	
}

void HydroStar::Save_RAO(String fileName, const UVector<bool> &listDOF, bool qtf) const {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to save '%s'. File already used."), fileName));
		
	out << "\n";
	
	out << "#Motion\n";
	for (int ib = 0; ib < dt.Nb; ++ib) {
		for (int idof = 0; idof < 6; ++idof)
			if (listDOF[idof]) {
				String sdof = BEM::StrDOF(idof, false);
				out << Format("%s  %d  FILE  %s.rao\n", Format("G%8<s", ToUpper(sdof)) + S("BODY"), ib+1, sdof);
			}
	}
	
	out << "\n";
	
	out << "#Added mass and damping\n";
	for (int ib = 0; ib < dt.Nb; ++ib) {
		out << Format("CA BODY  %d  FILE  ca.rao TERM", ib+1);
		for (int idof = 0; idof < 6; ++idof)
			if (listDOF[idof]) 
				out << " " << idof+1 << idof+1;
		out << "\n";
		out << Format("CM BODY  %d  FILE  cm.rao TERM", ib+1);
		for (int idof = 0; idof < 6; ++idof)
			if (listDOF[idof]) 
				out << " " << idof+1 << idof+1;
		out << "\n";
	}

	out << "\n";

	UVector<String> strdof = {"FX", "FY", "FZ", "MX", "MY", "MZ"};
		
	out << "#First order forces\n";
	for (int ib = 0; ib < dt.Nb; ++ib) {
		for (int idof = 0; idof < 6; ++idof)
			if (listDOF[idof]) 
				out << Format("%s BODY %d FILE %s.rao\n", ToUpper(strdof[idof]) + "F1ST", ib+1, ToLower(strdof[idof]) + "f1st");
	}
	
	if (qtf) {
		out << "\n";
			
		out << "#Full QTF (hsqtf)\n";
		for (int ib = 0; ib < dt.Nb; ++ib) {
			for (int idof = 0; idof < 6; ++idof)
				if (listDOF[idof]) 
					out << Format("QTF%s BODY %d FILE qtf_%s.rao\n", ToUpper(strdof[idof]), ib+1, ToLower(strdof[idof]));
		}
	}
	
	out << "\n";

	out << "ORCAFLEX";	
}

void HydroStar::Save_RDF(String fileName) const {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to save '%s'. File already used."), fileName));
		

	out << 	"#Name of the output file\n"
			"FILENAME rd1\n";
			
	out << 	"\n";
	
	out << 	"#Range of frequency\n"
			"FREQUENCY   TYPE   1\n";
	for (double ww : dt.w)
		out << Format(" %.3f", ww);
	out << 	"\n";
	out << 	"ENDFREQUENCY\n";
	
	out << 	"\n";
	
	out << 	"#Range of heading\n"
			"HEADING   TYPE   1\n";
	for (double hh : dt.head)
		out << Format(" %.3f", hh);
	out << 	"\n";
	out << 	"ENDHEADING\n";
	
	out << 	"\n";
	
	out << 	"#Waterdepth\n"
			"WATERDEPTH   " << (dt.h > 0 ? Format("%.3f", dt.h) : "inf") << "\n";
	
	out << "\n";
	
	out << "REFWAVE  0.0   0.0\n";
	for (int ib = 0; ib < dt.Nb; ++ib)
		out << Format("REFPOINT %d   %.3f   %.3f   %.3f\n", ib+1, dt.msh[ib].dt.c0.x, dt.msh[ib].dt.c0.y, dt.msh[ib].dt.c0.z);
	
	out << "\n";
	
	out << 	"SPEEDS TYPE 0\n"
			"1 0.0\n"
			"ENDSPEEDS\n";

	out << "\n";

	out << "ENDFILE";	
}

void HydroStar::Save_DFT(String fileName) const {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to save '%s'. File already used."), fileName));
	
	out <<  "NFORMULE YES\n"
			"FFORMULE YES\n"
			"MFORMULE NO\n"
			"ENDFILE";
}
	
bool HydroStar::Load_HDF(Wamit &wam, Function <bool(String, int)> Status) {
	Status("Loading HDF file", Null);
	
	FindFile ff(AFX(GetFileFolder(wam.dt.file), "*.hdf"));
	if (!ff)
		return false;
	
	Hdf5File hfile;
	hfile.Open(ff.GetPath(), H5F_ACC_RDONLY);
	
	hfile.ChangeGroup("VTKHDF");
	
	int ncell;
	ncell = hfile.GetInt("NumberOfCells");
	
	int ndata;
	UVector<double> w, head;
	{
		hfile.ChangeGroup("FieldData");
		
		hfile.GetDouble("frequency", w);
		hfile.GetDouble("heading", head);
	
		if (w.size() != head.size())
			throw Exc("'frequency' and 'heading' sizes mismatch");
				
		Upp::Index<double> w0, head0;	
		for (int i = 0; i < w.size(); ++i) {
			w0.FindAdd(w[i]);
			head0.FindAdd(head[i]);
		}
		
		if (w0.size() != wam.dt.w.size())
			throw Exc("'frequency' number of frequencies mismatch with .out");
		
		if (head0.size() != wam.dt.head.size())
			throw Exc("'heading' number of headings mismatch with .out");
	
		int nhead = head0.size();
		int nw = w0.size();
		
		//		(nh*nw (diffraction) + nh*nw (incident) + nh*nw*6dof (radiation))*2 (re/im)
		ndata = (nhead*nw            + nhead*nw         + nhead*nw*6)            *2;
		
		if (w.size() != ndata)
			throw Exc("Wrong size in 'frequency'");
		
		if (head.size() != ndata)
			throw Exc("Wrong size in 'heading'");
		
		
		hfile.UpGroup();
	}

	FindFile ff2(AFX(GetFileFolder(wam.dt.file), "*.hst"));
	if (!ff2)
		throw Exc (t_("Mesh data (.hst file) not found"));
	
	UArray<Body> mesh;
	bool y0z, x0z;
	String ret = Body::Load(mesh, ff2.GetPath(), wam.dt.rho, wam.dt.g, false, Null, Null, y0z, x0z);
	if (!ret.IsEmpty())
		throw Exc (Format(t_("Mesh data (.hst file) problem: %s"), ret));
	
	First(wam.dt.msh).dt.mesh = pick(First(mesh).dt.mesh);
	
	if (ncell != First(wam.dt.msh).dt.mesh.GetNumPanels())
		throw Exc(Format(t_("The number of panels in hdf (%d) doesn't match with the number in the hst mesh (%d)"), ncell, First(wam.dt.msh).dt.mesh.GetNumPanels()));
	
	wam.Initialize_PotsRad();
	wam.Initialize_PotsIncDiff(wam.dt.pots_inc);
	wam.Initialize_PotsIncDiff(wam.dt.pots_dif);
	
	{
		hfile.ChangeGroup("CellData");

		Status("Listing HDF data", Null);
		
		UVector<String> names = hfile.ListGroup();	// Faster than ListGroupDatasets(), and valid as there are only datasets
		UVector<double> data;
		int ib = 0;
		for (int id = 0; id < names.size(); ++id) {
			if (!(id%100))
				Status("Processing HDF panel data", (id*100)/names.size());
			
			const String &name = names[id];
			hfile.GetDouble(name, data);
			
			if (data.size()	!= ncell)
				throw Exc(Format(t_("Wrong number of items in '%s'"), name));

			bool real = name.EndsWith("_RE");
			char type;
			String stype = name.Mid(15, 3);
			if (stype == "RAD")
				type = 'r';
			else if (stype == "INC")
				type = 'i';
			else if (stype == "DIF")
				type = 'd';
			else
				throw Exc(Format(t_("Unknown type in '%s'"), name));
			
			int idf;
			if (type == 'r') {
				String sdof = name.Mid(19, 2);	
				if (sdof == "su")
					idf = 0;
				else if (sdof == "sw")
					idf = 1;
				else if (sdof == "he")
					idf = 2;
				else if (sdof == "ro")
					idf = 3;
				else if (sdof == "pi")
					idf = 4;
				else if (sdof == "ya")
					idf = 5;
				else
					throw Exc(Format(t_("Unknown dof in '%s'"), name));
			}
			int ifr = FindRoundDecimals(wam.dt.w, w[id], 3);
			if (ifr < 0)
				throw Exc(Format(t_("Frequency %f found in the .hdf but not found in the .out file"), w[id]));
			
			int ihead = FindRoundDecimals(wam.dt.head, head[id], 1);
			if (ihead < 0)
				throw Exc(Format(t_("Heading %f found in the .hdf but not found in the .out file"), head[id]));
			
			double factor = wam.dt.g/w[id];
			
			for (int ip = 0; ip < ncell; ++ip) {
				std::complex<double> *d;
				if (type == 'r') {
					d = &wam.dt.pots_rad[ib][ip][idf][ifr];
					if (real)							// p = -iρωΦ ; Φ = [Im(p) - iRe(p)]/ρω
						d->imag(-data[ip]*factor);
					else
						d->real(data[ip]*factor);
				} else if (type == 'i')	{
					d = &wam.dt.pots_inc[ib][ip][ihead][ifr];
					if (real)							// p = -iρωΦ ; Φ = [Im(p) - iRe(p)]/ρω
						d->imag(data[ip]*factor);
					else
						d->real(-data[ip]*factor);
				} else {
					d = &wam.dt.pots_dif[ib][ip][ihead][ifr];
					if (real)							// p = -iρωΦ ; Φ = [Im(p) - iRe(p)]/ρω
						d->imag(data[ip]*factor);
					else
						d->real(-data[ip]*factor);
				}
			}
		}
		hfile.UpGroup();		
	}
		
	return true;
}

// This is only for HydroStar
bool HydroStar::Load_MCN(String fileName, int nb, UVector<Point3D> &refPoint, UVector<Pointf> &refWave) {
	if (nb == 0)
		return false;
	
	String folder = GetFileFolder(fileName);
	FindFile mcn(AFX(folder, "*.mcn"));
	if (!mcn)
		return false;
	
	String mcnFile = mcn.GetPath();
	
	FileInLine in(mcnFile);
	if (!in.IsOpen())
		return false;

	LineParserWamit f(in);
	f.IsSeparator = IsTabSpace;
	
	UVector<Point3D> cog(nb, Null);
	refPoint.SetCount(nb, Null);
	refWave.SetCount(nb, Null);
	
	while(!in.IsEof()) {
		f.GetLine_discard_empty();
		if (f.IsEof())
			break;

		if (f.GetText(0) == "COGPOINT_BODY") {
			int ib = f.GetInt(1);
			if (ib < 1 || ib > nb)
				throw Exc(in.Str() + "\n"  + t_("Wrong body in COGPOINT_BODY"));
			cog[ib-1] = Point3D(f.GetDouble(2), f.GetDouble(3), f.GetDouble(4));
		} else if (f.GetText(0).StartsWith("REFPOINT")) {
			int ib = f.GetInt(1);
			if (ib < 1 || ib > nb)
				throw Exc(in.Str() + "\n"  + t_("Wrong body in REFPOINT"));
			refPoint[ib-1] = Point3D(f.GetDouble(2), f.GetDouble(3), f.GetDouble(4));
		} else if (f.GetText(0).StartsWith("REFWAVE")) 
			refWave.Set(0, Pointf(f.GetDouble(1), f.GetDouble(2)), nb);
	}
	if (IsNull(cog[0]))
		return false;
	
	if (IsNull(refPoint[0]))
		refPoint = clone(cog);
	
	return true;
}

