// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2025, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include <MatIO/matio.h>


String Matlab::Load(String file, double) {
	dt.file = file;
	dt.name = GetFileTitle(file);
	dt.dimen = true;
	dt.len = 1;
	dt.solver = Hydro::MATLAB;
	dt.Nb = Null;
	dt.x_w = dt.y_w = 0;
	
	try {
		BEM::Print("\n\n" + Format(t_("Loading '%s'"), file));

		BEM::Print("\n- " + S(t_("Mat file")));
		
		Load_Mat();
		
		if (IsNull(dt.Nb))
			return t_("No data found");
	
	} catch (Exc e) {
		return e;
	}
	
	return String();
}

void Matlab::Load_Mat() {
	String fileName = ForceExtSafer(dt.file, ".mat");
	
	MatFile mfile;
	
	if (!mfile.OpenRead(fileName))
		throw Exc(t_("Impossible to load .mat file"));
	
	{
		UVector<String> freq = {"w", "frequency", "frequencies"};
		int id = mfile.Exist(freq, true);
		if (id < 0)
			throw Exc(t_("Frequencies not found"));
		mfile.Get(freq[id], dt.w, true);
		dt.Nf = dt.w.size();
	}
	{
		UVector<String> head = {"head", "headings", "h"};
		int id = mfile.Exist(head, true);
		if (id < 0)
			throw Exc(t_("Headings not found"));
		mfile.Get(head[id], dt.head, true);
		dt.Nh = dt.head.size();
	}
	dt.h = mfile.Get<double>("depth");
	{
		UVector<String> g = {"g", "gravity"};
		int id = mfile.Exist(g, true);
		if (id < 0)
			throw Exc(t_("Gravity not found"));
		dt.g = mfile.Get<double>(g[id], true);
	}
	{
		UVector<String> rho = {"rho", "density"};
		int id = mfile.Exist(rho, true);
		if (id < 0)
			throw Exc(t_("Density not found"));
		dt.rho = mfile.Get<double>(rho[id], true);
	}
	dt.Nb = Null;
	
	auto LoadAB = [&](UArray<UArray<VectorXd>> &a, const char *name, const UVector<String> &sA) {
		int id = mfile.Exist(sA, true);	
		if (id >= 0) {
			MultiDimMatrix<double> A;
			mfile.Get(sA[id], A, true);
			if (A.GetAxisDim(2) != dt.Nf)
				throw Exc(Format(t_("%s dimension is not correct (Nf) (6*nBodies, 6*nBodies, nFrequencies)"), name));
			int dim = A.GetAxisDim(0);
			if (dim != A.GetAxisDim(1))
				throw Exc(Format(t_("%s dimension is not correct (Nf) (6*nBodies, 6*nBodies, nFrequencies)"), name));
			if (!(dim%6))
				throw Exc(Format(t_("%s dimension is not correct (Nf) (6*nBodies, 6*nBodies, nFrequencies)"), name));
			int Nb = dim/6;
			if (!IsNull(dt.Nb) && dt.Nb != Nb)
				throw Exc(Format(t_("Body number mismatch in %s"), name));
			dt.Nb = Nb;
			Initialize_AB(a, 0);
			for (int r = 0; r < 6*dt.Nb; ++r) 
				for (int c = 0; c < 6*dt.Nb; ++c) 
					for (int ifr = 0; ifr < dt.Nf; ++ifr) 
						a[r][c][ifr] = A(r, c, ifr);
		}
	};
	LoadAB(dt.A, "A", {"A", "added_mass", "addedmass"});
	LoadAB(dt.B, "B", {"B", "radiation_damping", "radiationdamping"});

	{
		UVector<String> s = {"Ainf"};
		int id = mfile.Exist(s, true);	
		if (id >= 0) {
			mfile.Get(s[id], dt.Ainf, true);
			if ((dt.Ainf.rows() != dt.Ainf.cols()) || (!IsNull(dt.Nb) && (6*dt.Ainf.rows() != 6*dt.Nb)))
				throw Exc(t_("Ainf dimension is not correct (6*nBodies, 6*nBodies)"));
			dt.Nb = dt.Ainf.cols()/6;
		}
	}
	{
		if (mfile.Exist("A0", true)) {
			mfile.Get("A0", dt.A0, true);
			if ((dt.A0.rows() != dt.A0.cols()) || (!IsNull(dt.Nb) && (6*dt.A0.rows() != 6*dt.Nb)))
				throw Exc(t_("A0 dimension is not correct (6*nBodies, 6*nBodies)"));
			dt.Nb = dt.A0.cols()/6;
		}
	}
	{
		if (mfile.Exist("vol", true)) {
			UVector<double> vol;
			mfile.Get("vol", vol, true);
			if (!IsNull(dt.Nb) && vol.size() != dt.Nb)
				throw Exc(t_("Vol dimension is not correct (1, nBodies)"));
			dt.Nb = vol.size();
			if (dt.msh.IsEmpty())
				dt.msh.SetCount(dt.Nb);
			for (int ib = 0; ib < dt.Nb; ++ib)
				dt.msh[ib].dt.Vo = vol[ib];
		}
	}
	{
		if (mfile.Exist("Cg", true)) {
			MatrixXd cg;
			mfile.Get("Cg", cg, true);
			if (!IsNull(dt.Nb) && cg.cols() != dt.Nb)
				throw Exc(t_("Cg dimension is not correct (cols != Nb) (3, nBodies)"));
			if (cg.rows() != 3)
				throw Exc(t_("Cg dimension is not correct (rows != 3) (3, nBodies)"));
			dt.Nb = cg.cols();
			if (dt.msh.IsEmpty())
				dt.msh.SetCount(dt.Nb);
			for (int ib = 0; ib < dt.Nb; ++ib)
				dt.msh[ib].dt.cg = Point3D(cg(0, ib), cg(1, ib), cg(2, ib));
		}
	}
	{
		if (mfile.Exist("Cb", true)) {
			MatrixXd cb;
			mfile.Get("Cb", cb, true);
			if (!IsNull(dt.Nb) && cb.cols() != dt.Nb)
				throw Exc(t_("Cb dimension is not correct (cols != Nb) (3, nBodies)"));
			if (cb.rows() != 3)
				throw Exc(t_("Cb dimension is not correct (rows != 3) (3, nBodies)"));
			dt.Nb = cb.cols();
			if (dt.msh.IsEmpty())
				dt.msh.SetCount(dt.Nb);
			for (int ib = 0; ib < dt.Nb; ++ib)
				dt.msh[ib].dt.cb = Point3D(cb(0, ib), cb(1, ib), cb(2, ib));
		}
	}
	{
		if (mfile.Exist("C0", true)) {
			MatrixXd c0;
			mfile.Get("C0", c0, true);
			if (!IsNull(dt.Nb) && c0.cols() != dt.Nb)
				throw Exc(t_("C0 dimension is not correct (cols != Nb) (3, nBodies)"));
			if (c0.rows() != 3)
				throw Exc(t_("C0 dimension is not correct (rows != 3) (3, nBodies)"));
			dt.Nb = c0.cols();
			if (dt.msh.IsEmpty())
				dt.msh.SetCount(dt.Nb);
			for (int ib = 0; ib < dt.Nb; ++ib)
				dt.msh[ib].dt.cg = Point3D(c0(0, ib), c0(1, ib), c0(2, ib));
		}
	}	
	{
		UVector<String> sB = {"Kh", "K", "stiffness_matrix", "stiffnessmatrix"};
		int id = mfile.Exist(sB, true);	
		if (id >= 0) {
			MatrixXd C;
			mfile.Get(sB[id], C, true);
			if ((C.rows() != C.cols()) || (!IsNull(dt.Nb) && (6*C.rows() != 6*dt.Nb)))
				throw Exc(t_("Kh dimension is not correct (6*nBodies, 6*nBodies)"));
			dt.Nb = C.cols()/6;
			for (int ib = 0; ib < dt.Nb; ++ib)
				dt.msh[ib].dt.C = C.block<6, 6>(6*ib, 6*ib);
		}
	}
	{
		UVector<String> sB = {"mass", "inertia", "m"};
		int id = mfile.Exist(sB, true);	
		if (id >= 0) {
			MatrixXd M;
			mfile.Get(sB[id], M, true);
			if ((M.rows() != M.cols()) || (!IsNull(dt.Nb) && (6*M.rows() != 6*dt.Nb)))
				throw Exc(t_("Kh dimension is not correct (6*nBodies, 6*nBodies)"));
			dt.Nb = M.cols()/6;
			for (int ib = 0; ib < dt.Nb; ++ib)
				dt.msh[ib].dt.M = M.block<6, 6>(6*ib, 6*ib);
		}
	}
	auto LoadF = [&](Forces &f, const char *name, const UVector<String> &sB) {
		int id = mfile.Exist(sB, true);	
		if (id >= 0) {
			MultiDimMatrix<std::complex<double>> F;
			mfile.Get(sB[id], F, true);
			if (F.GetNumAxis() != 4)
				throw Exc(Format(t_("%s dimension is not correct (6, nHeadings, nFrequencies, nBodies)"), name));
			if (F.GetAxisDim(0) != 6)
				throw Exc(Format(t_("%s dimension is not correct (6, nHeadings, nFrequencies, nBodies)"), name));
			if (!IsNull(dt.Nh) &&  dt.Nh != F.GetAxisDim(1)) 
				throw Exc(Format(t_("%s dimension is not correct (6, nHeadings, nFrequencies, nBodies)"), name));
			else
				dt.Nh = F.GetAxisDim(1);
			if (!IsNull(dt.Nf) &&  dt.Nf != F.GetAxisDim(2)) 
				throw Exc(Format(t_("%s dimension is not correct (6, nHeadings, nFrequencies, nBodies)"), name));
			else
				dt.Nf = F.GetAxisDim(2);
			if (!IsNull(dt.Nb) &&  dt.Nb != F.GetAxisDim(3)) 
				throw Exc(Format(t_("%s dimension is not correct (6, nHeadings, nFrequencies, nBodies)"), name));
			else
				dt.Nb = F.GetAxisDim(3);
			Initialize_Forces(f);
			for (int idf = 0; idf < 6; ++idf) 
				for (int ih = 0; ih < dt.Nh; ++ih) 
					for (int ib = 0; ib < dt.Nb; ++ib) 		
						for (int ifr = 0; ifr < dt.Nf; ++ifr) 
							f[ib][ih](ifr, idf) = F(idf, ih, ifr, ib);			
		}
	};
	LoadF(dt.ex, "Fex", {"excitation_force", "excitationforce", "fex"});
	LoadF(dt.fk, "Fex", {"froude_krylov", "froudekrylov", "ffk"});
	LoadF(dt.sc, "Fex", {"scattering_force", "scatteringforce", "diffraction_force", "diffractionforce", "fsc"});	
}
	
	
void Matlab::Save(String file) const {
	String fileName = ForceExtSafer(file, ".mat");
	
	MatFile mfile;
	
	if (!mfile.OpenCreate(fileName, MAT_FT_MAT73))
		throw Exc(t_("Impossible to create .mat file"));
	
	mfile.Set("w", dt.w);
	mfile.Set("head", dt.head);
	mfile.Set("depth", dt.h);
	mfile.Set("rho", dt.rho);
	mfile.Set("g", dt.g);
	
	if (IsLoadedA()) {
		MultiDimMatrix<double> A(6*dt.Nb, 6*dt.Nb, dt.Nf);
		for (int r = 0; r < 6*dt.Nb; ++r) {
			for (int c = 0; c < 6*dt.Nb; ++c) {
				for (int ifr = 0; ifr < dt.Nf; ++ifr) {
					if (dt.A[r][c].size() == dt.Nf) 
						A(r, c, ifr) = dt.A[r][c][ifr];
					else
						A(r, c, ifr) = 0;
				}
			}
		}
		mfile.Set("A", A);
	}
	if (IsLoadedB()) {
		MultiDimMatrix<double> B(6*dt.Nb, 6*dt.Nb, dt.Nf);
		for (int r = 0; r < 6*dt.Nb; ++r) {
			for (int c = 0; c < 6*dt.Nb; ++c) {
				for (int ifr = 0; ifr < dt.Nf; ++ifr) {
					if (dt.B[r][c].size() == dt.Nf) 
						B(r, c, ifr) = dt.B[r][c][ifr];
					else
						B(r, c, ifr) = 0;
				}
			}
		}
		mfile.Set("B", B);
	}
	
	if (IsLoadedAinf()) 
		mfile.Set("Ainf", dt.Ainf);
		
	if (IsLoadedA0()) 
		mfile.Set("A0", dt.A0);
	
	if ((!IsLoadedFsc() || !IsLoadedFfk()) && IsLoadedFex()) {
		MultiDimMatrix<std::complex<double>> Fex(6, dt.Nh, dt.Nf, dt.Nb);
		for (int idf = 0; idf < 6; ++idf) 
			for (int ih = 0; ih < dt.Nh; ++ih) 
				for (int ib = 0; ib < dt.Nb; ++ib) 		
					for (int ifr = 0; ifr < dt.Nf; ++ifr) 
						Fex(idf, ih, ifr, ib) = dt.ex[ib][ih](ifr, idf);
		mfile.Set("Fex", Fex);				
	}
	if (IsLoadedFfk()) {
		MultiDimMatrix<std::complex<double>> Ffk(6, dt.Nh, dt.Nf, dt.Nb);
		for (int idf = 0; idf < 6; ++idf) 
			for (int ih = 0; ih < dt.Nh; ++ih) 
				for (int ib = 0; ib < dt.Nb; ++ib) 		
					for (int ifr = 0; ifr < dt.Nf; ++ifr) 
						Ffk(idf, ih, ifr, ib) = dt.fk[ib][ih](ifr, idf);
		mfile.Set("Ffk", Ffk);	
	}
	if (IsLoadedFsc()) {
		MultiDimMatrix<std::complex<double>> Fsc(6, dt.Nh, dt.Nf, dt.Nb);
		for (int idf = 0; idf < 6; ++idf) 
			for (int ih = 0; ih < dt.Nh; ++ih) 
				for (int ib = 0; ib < dt.Nb; ++ib) 		
					for (int ifr = 0; ifr < dt.Nf; ++ifr) 
						Fsc(idf, ih, ifr, ib) = dt.sc[ib][ih](ifr, idf);
		mfile.Set("Fsc", Fsc);	
	}
	
	if (IsLoadedAnyC()) {
		MatrixXd C = MatrixXd::Zero(6*dt.Nb, 6*dt.Nb);
		for (int ib = 0; ib < dt.Nb; ++ib)
			if (IsLoadedC(ib, 5, 5))
				C.block<6, 6>(6*ib, 6*ib) = dt.msh[ib].dt.C;
		mfile.Set("Kh", C);	
	}
	
	if (IsLoadedAnyM()) {
		MatrixXd M = MatrixXd::Zero(6*dt.Nb, 6*dt.Nb);
		for (int ib = 0; ib < dt.Nb; ++ib)
			M.block<6, 6>(6*ib, 6*ib) = dt.msh[ib].dt.M;
		mfile.Set("Mass", M);	
	}

	VectorXd vol(dt.Nb);
	for (int ib = 0; ib < dt.Nb; ++ib) {
		if (!IsNull(dt.msh[ib].dt.cb)) 
			vol(ib) = dt.msh[ib].dt.Vo;
		else
			vol(ib) = NaNDouble;
	}
	mfile.Set("Vol", vol);
	
	MatrixXd Cb(3, dt.Nb);
	for (int ib = 0; ib < dt.Nb; ++ib) {
		for (int i = 0; i < 3; ++i) {
			if (!IsNull(dt.msh[ib].dt.cb)) 
				Cb(i, ib) = dt.msh[ib].dt.cb[i];
			else
				Cb(i, ib) = NaNDouble;
		}
	}
	mfile.Set("Cb", Cb);
		
	MatrixXd Cg(3, dt.Nb);
	for (int ib = 0; ib < dt.Nb; ++ib) {
		for (int i = 0; i < 3; ++i) {
			if (!IsNull(dt.msh[ib].dt.cb)) 
				Cg(i, ib) = dt.msh[ib].dt.cg[i];
			else
				Cg(i, ib) = NaNDouble;
		}
	}
	mfile.Set("Cg", Cg);
	
	MatrixXd C0(3, dt.Nb);
	for (int ib = 0; ib < dt.Nb; ++ib) {
		for (int i = 0; i < 3; ++i) {
			if (!IsNull(dt.msh[ib].dt.cb)) 
				C0(i, ib) = dt.msh[ib].dt.c0[i];
			else
				C0(i, ib) = NaNDouble;
		}
	}		
	mfile.Set("C0", C0);
	
	UVector<String> legend = {"DoF: surge, sway, heave, roll, pitch, yaw",
					"Units: m, Kg, N, s, rad",
					"Parameters:",
					"- w: Frequencies (1, nFrequencies)",
					"- head: Headings (1, nHeadings)",
					"- depth: Calculation depth",
					"- rho: Water density",
					"- g: Gravity",
					"- A: Added mass (6*nBodies, 6*nBodies, nFrequencies)",
					"- B: Radiation damping (6*nBodies, 6*nBodies, nFrequencies)",
					"- Ainf: Added mass at infinite frequency (6*nBodies, 6*nBodies)",
					"- A0: Added mass at frequency = 0 (6*nBodies, 6*nBodies)",
					"- Fex = Ffk + Fsc: Excitation force (imaginary numbers) (6, nHeadings, nFrequencies, nBodies)",
					"- Ffk: Froude-Krylov force (imaginary numbers) (6, nHeadings, nFrequencies, nBodies)",
					"- Fsc: Diffracion scattering force (imaginary numbers) (6, nHeadings, nFrequencies, nBodies)",
					"- Kh: Hydrostatic stiffness matrix (6*nBodies, 6*nBodies)",
					"- Mass: Inertia matrix (6*nBodies, 6*nBodies)",
					"- Vol: Submerged volume (1, nBodies)",
					"- Cg: Centres of gravity (3, nBodies)",
					"- Cb: Centres of buoyancy (3, nBodies)",
					"- C0: Centres of body axis (3, nBodies)",
					"Comments:",
					"- Global axis (wave origin) is at 0, 0",
					"- The incident wave system follows Wamit criteria",
					"File generated by BEMRosetta: Hydrodynamic solvers viewer and converter"};
	mfile.SetCell("Legend", legend);
}
