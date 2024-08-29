// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2023, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#include <Hdf5/hdf5.h>


bool Wamit::Load_HDF(Function <bool(String, int)> Status) {
	Status("Loading HDF file", Null);
	
	FindFile ff(AFX(GetFileFolder(dt.file), "*.hdf"));
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
		
		if (w0.size() != dt.w.size())
			throw Exc("'frequency' number of frequencies mismatch with .out");
		
		if (head0.size() != dt.head.size())
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

	FindFile ff2(AFX(GetFileFolder(dt.file), "*.hst"));
	if (!ff2)
		throw Exc (t_("Mesh data (.hst file) not found"));
	
	UArray<Body> mesh;
	bool y0z, x0z;
	String ret = Body::Load(mesh, ff2.GetPath(), dt.rho, dt.g, false, Null, Null, y0z, x0z);
	if (!ret.IsEmpty())
		throw Exc (Format(t_("Mesh data (.hst file) problem: %s"), ret));
	
	First(dt.msh).dt.mesh = pick(First(mesh).dt.mesh);
	
	if (ncell != First(dt.msh).dt.mesh.GetNumPanels())
		throw Exc(Format(t_("The number of panels in hdf (%d) doesn't match with the number in the hst mesh (%d)"), ncell, First(dt.msh).dt.mesh.GetNumPanels()));
	
	Initialize_PotsRad();
	Initialize_PotsIncDiff(dt.pots_inc);
	Initialize_PotsIncDiff(dt.pots_dif);
	
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
			int ifr = FindRoundDecimals(dt.w, w[id], 3);
			if (ifr < 0)
				throw Exc(Format(t_("Frequency %f found in the .hdf but not found in the .out file"), w[id]));
			
			int ihead = FindRoundDecimals(dt.head, head[id], 1);
			if (ihead < 0)
				throw Exc(Format(t_("Heading %f found in the .hdf but not found in the .out file"), head[id]));
			
			double factor = dt.g/w[id];
			
			for (int ip = 0; ip < ncell; ++ip) {
				std::complex<double> *d;
				if (type == 'r')
					d = &dt.pots_rad[ib][ip][idf][ifr];
				else if (type == 'i')	
					d = &dt.pots_inc[ib][ip][ihead][ifr];
				else
					d = &dt.pots_dif[ib][ip][ihead][ifr];
				
				if (real)							// p = iρωΦ		Φ = Im(p/ρω) - iRe(p/ρω)
					d->imag(-data[ip]*factor);
				else
					d->real(data[ip]*factor);
			}
		}
		hfile.UpGroup();		
	}
		
	return true;
}



