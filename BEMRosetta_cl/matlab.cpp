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
		MatVar var = mfile.GetVar({"w", "frequency", "frequencies"}, true);
		if (!var.IsLoaded())
			throw Exc(t_("Frequencies not found"));
		dt.w = mfile.ReadRowMajor<double>(var);
		dt.Nf = dt.w.size();
	}
	{
		MatVar head = mfile.GetVar({"head", "headings", "h"}, true);
		if (!head.IsLoaded())
			throw Exc(t_("Headings not found"));
		dt.head = mfile.ReadRowMajor<double>(head);
		dt.Nh = dt.head.size();
	}
	dt.h = mfile.ReadDouble(mfile.GetVar("depth", true));
	{
		MatVar g = mfile.GetVar({"g", "gravity"});
		if (!g.IsLoaded())
			throw Exc(t_("Gravity not found"));
		dt.g = mfile.ReadDouble(g);
	}
	{
		MatVar rho = mfile.GetVar({"rho", "density"}, true);
		if (!rho.IsLoaded())
			throw Exc(t_("Density not found"));
		dt.rho = mfile.ReadDouble(rho);
	}
	dt.Nb = Null;
	
	auto LoadAB = [&](UArray<UArray<VectorXd>> &a, const char *name, const UVector<String> &sA) {
		MatVar id = mfile.GetVar(sA, true);	
		if (id.IsLoaded()) {
			MultiDimMatrix<double> A = MatFile::ReadMultiDim<double>(id);
			if (A.GetAxisDim(2) != dt.Nf)
				throw Exc(Format(t_("%s dimension is not correct (Nf) (6*nBodies, 6*nBodies, nFrequencies)"), name));
			int dim = A.GetAxisDim(0);
			if (dim != A.GetAxisDim(1))
				throw Exc(Format(t_("%s dimension is not correct (Nf) (6*nBodies, 6*nBodies, nFrequencies)"), name));
			if (dim%6)
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
	LoadAB(dt.A, "A", {"A", "addedmass"});
	LoadAB(dt.B, "B", {"B", "radiationdamping"});

	{
		MatVar id = mfile.GetVar("Ainf", true);
		if (id.IsLoaded()) {
			dt.Ainf = MatFile::ReadEigenDouble(id);
			if ((dt.Ainf.rows() != dt.Ainf.cols()) || (!IsNull(dt.Nb) && (dt.Ainf.rows() != 6*dt.Nb)))
				throw Exc(t_("Ainf dimension is not correct (6*nBodies, 6*nBodies)"));
			dt.Nb = (int)dt.Ainf.cols()/6;
		}
	}
	{
		MatVar id = mfile.GetVar("A0", true);
		if (id.IsLoaded()) {
			dt.A0 = MatFile::ReadEigenDouble(id);
			if ((dt.A0.rows() != dt.A0.cols()) || (!IsNull(dt.Nb) && (dt.A0.rows() != 6*dt.Nb)))
				throw Exc(t_("A0 dimension is not correct (6*nBodies, 6*nBodies)"));
			dt.Nb = (int)dt.A0.cols()/6;
		}
	}
	{
		MatVar id = mfile.GetVar("vol", true);
		if (id.IsLoaded()) {
			UVector<double> vol = MatFile::ReadColMajor<double>(id);
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
		MatVar id = mfile.GetVar("Cg", true);
		if (id.IsLoaded()) {
			MatrixXd cg = MatFile::ReadEigenDouble(id);
			if (!IsNull(dt.Nb) && cg.cols() != dt.Nb)
				throw Exc(t_("Cg dimension is not correct (cols != Nb) (3, nBodies)"));
			if (cg.rows() != 3)
				throw Exc(t_("Cg dimension is not correct (rows != 3) (3, nBodies)"));
			dt.Nb = (int)cg.cols();
			if (dt.msh.IsEmpty())
				dt.msh.SetCount(dt.Nb);
			for (int ib = 0; ib < dt.Nb; ++ib)
				dt.msh[ib].dt.cg = Point3D(cg(0, ib), cg(1, ib), cg(2, ib));
		}
	}
	{
		MatVar id = mfile.GetVar("Cb", true);
		if (id.IsLoaded()) {
			MatrixXd cb = MatFile::ReadEigenDouble(id);
			if (!IsNull(dt.Nb) && cb.cols() != dt.Nb)
				throw Exc(t_("Cb dimension is not correct (cols != Nb) (3, nBodies)"));
			if (cb.rows() != 3)
				throw Exc(t_("Cb dimension is not correct (rows != 3) (3, nBodies)"));
			dt.Nb = (int)cb.cols();
			if (dt.msh.IsEmpty())
				dt.msh.SetCount(dt.Nb);
			for (int ib = 0; ib < dt.Nb; ++ib)
				dt.msh[ib].dt.cb = Point3D(cb(0, ib), cb(1, ib), cb(2, ib));
		}
	}
	{
		MatVar id = mfile.GetVar("C0", true);
		if (id.IsLoaded()) {
			MatrixXd c0 = MatFile::ReadEigenDouble(id);
			if (!IsNull(dt.Nb) && c0.cols() != dt.Nb)
				throw Exc(t_("C0 dimension is not correct (cols != Nb) (3, nBodies)"));
			if (c0.rows() != 3)
				throw Exc(t_("C0 dimension is not correct (rows != 3) (3, nBodies)"));
			dt.Nb = (int)c0.cols();
			if (dt.msh.IsEmpty())
				dt.msh.SetCount(dt.Nb);
			for (int ib = 0; ib < dt.Nb; ++ib)
				dt.msh[ib].dt.c0 = Point3D(c0(0, ib), c0(1, ib), c0(2, ib));
		}
	}	
	{
		MatVar id = mfile.GetVar({"Kh", "K", "stiffness_matrix", "stiffnessmatrix"}, true);
		if (id.IsLoaded()) {
			MatrixXd C = MatFile::ReadEigenDouble(id);
			if ((C.rows() != C.cols()) || (!IsNull(dt.Nb) && (C.rows() != 6*dt.Nb)))
				throw Exc(t_("Kh dimension is not correct (6*nBodies, 6*nBodies)"));
			dt.Nb = (int)C.cols()/6;
			for (int ib = 0; ib < dt.Nb; ++ib)
				dt.msh[ib].dt.C = C.block<6, 6>(6*ib, 6*ib);
		}
	}
	{
		MatVar id = mfile.GetVar({"mass", "inertia", "m"}, true);
		if (id.IsLoaded()) {
			MatrixXd M = MatFile::ReadEigenDouble(id);
			if ((M.rows() != M.cols()) || (!IsNull(dt.Nb) && (M.rows() != 6*dt.Nb)))
				throw Exc(t_("Kh dimension is not correct (6*nBodies, 6*nBodies)"));
			dt.Nb = (int)M.cols()/6;
			for (int ib = 0; ib < dt.Nb; ++ib)
				dt.msh[ib].dt.M = M.block<6, 6>(6*ib, 6*ib);
		}
	}
	auto LoadF = [&](Forces &f, const char *name, const UVector<String> &sB) {
		MatVar id = mfile.GetVar(sB, true);	
		if (id.IsLoaded()) {
			MultiDimMatrix<std::complex<double>> F = MatFile::ReadMultiDim<std::complex<double>>(id);
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
	LoadF(dt.ex, "Fex", {"excitationforce", "fex"});
	LoadF(dt.fk, "Ffk", {"froudekrylov", "ffk"});
	LoadF(dt.sc, "Fsc", {"scatteringforce", "diffractionforce", "fsc"});	
}
	
	
void Matlab::Save(String file) const {
	String fileName = ForceExtSafer(file, ".mat");
	
	MatFile mfile;
	
	if (!mfile.OpenCreate(fileName, MAT_FT_MAT73))
		throw Exc(t_("Impossible to create .mat file"));
	
	mfile.Write("w", dt.w);
	mfile.Write("head", dt.head);
	mfile.Write("depth", dt.h);
	mfile.Write("rho", dt.rho);
	mfile.Write("g", dt.g);
	
	if (IsLoadedA()) {
		MultiDimMatrix<double> A(6*dt.Nb, 6*dt.Nb, dt.Nf);
		for (int r = 0; r < 6*dt.Nb; ++r) {
			for (int c = 0; c < 6*dt.Nb; ++c) {
				for (int ifr = 0; ifr < dt.Nf; ++ifr) {
					if (dt.A[r][c].size() == dt.Nf) 
						A(r, c, ifr) = A_dim(ifr, r, c);
					else
						A(r, c, ifr) = 0;
				}
			}
		}
		mfile.Write("A", A);
	}
	if (IsLoadedB()) {
		MultiDimMatrix<double> B(6*dt.Nb, 6*dt.Nb, dt.Nf);
		for (int r = 0; r < 6*dt.Nb; ++r) {
			for (int c = 0; c < 6*dt.Nb; ++c) {
				for (int ifr = 0; ifr < dt.Nf; ++ifr) {
					if (dt.B[r][c].size() == dt.Nf) 
						B(r, c, ifr) = B_dim(ifr, r, c);
					else
						B(r, c, ifr) = 0;
				}
			}
		}
		mfile.Write("B", B);
	}
	
	if (IsLoadedAinf()) 
		mfile.Write("Ainf", Ainf_mat(false));
		
	if (IsLoadedA0()) 
		mfile.Write("A0", A0_mat(false));
	
	if ((!IsLoadedFsc() || !IsLoadedFfk()) && IsLoadedFex()) {
		MultiDimMatrix<std::complex<double>> Fex(6, dt.Nh, dt.Nf, dt.Nb);
		for (int idf = 0; idf < 6; ++idf) 
			for (int ih = 0; ih < dt.Nh; ++ih) 
				for (int ib = 0; ib < dt.Nb; ++ib) 		
					for (int ifr = 0; ifr < dt.Nf; ++ifr) 
						Fex(idf, ih, ifr, ib) = F_dim(dt.ex, ih, ifr, idf, ib);
		mfile.Write("Fex", Fex);				
	}
	if (IsLoadedFfk()) {
		MultiDimMatrix<std::complex<double>> Ffk(6, dt.Nh, dt.Nf, dt.Nb);
		for (int idf = 0; idf < 6; ++idf) 
			for (int ih = 0; ih < dt.Nh; ++ih) 
				for (int ib = 0; ib < dt.Nb; ++ib) 		
					for (int ifr = 0; ifr < dt.Nf; ++ifr) 
						Ffk(idf, ih, ifr, ib) = F_dim(dt.fk, ih, ifr, idf, ib);
		mfile.Write("Ffk", Ffk);	
	}
	if (IsLoadedFsc()) {
		MultiDimMatrix<std::complex<double>> Fsc(6, dt.Nh, dt.Nf, dt.Nb);
		for (int idf = 0; idf < 6; ++idf) 
			for (int ih = 0; ih < dt.Nh; ++ih) 
				for (int ib = 0; ib < dt.Nb; ++ib) 		
					for (int ifr = 0; ifr < dt.Nf; ++ifr) 
						Fsc(idf, ih, ifr, ib) = F_dim(dt.sc, ih, ifr, idf, ib);
		mfile.Write("Fsc", Fsc);	
	}
	
	if (IsLoadedAnyC()) {
		MatrixXd C = MatrixXd::Zero(6*dt.Nb, 6*dt.Nb);
		for (int ib = 0; ib < dt.Nb; ++ib)
			if (IsLoadedC(ib, 5, 5))
				C.block<6, 6>(6*ib, 6*ib) = C_mat(false, ib);
		mfile.Write("Kh", C);	
	}
	
	if (IsLoadedAnyM()) {
		MatrixXd M = MatrixXd::Zero(6*dt.Nb, 6*dt.Nb);
		for (int ib = 0; ib < dt.Nb; ++ib)
			M.block<6, 6>(6*ib, 6*ib) = dt.msh[ib].dt.M;
		mfile.Write("Mass", M);	
	}

	VectorXd vol(dt.Nb);
	bool vnan = true;
	for (int ib = 0; ib < dt.Nb; ++ib) {
		if (!IsNull(dt.msh[ib].dt.Vo)) {
			vol(ib) = dt.msh[ib].dt.Vo;
			vnan = false;
		} else
			vol(ib) = NaNDouble;
	}
	if (!vnan)
		mfile.Write("Vol", vol);
	
	MatrixXd Cb(3, dt.Nb);
	bool cbnan = true;
	for (int ib = 0; ib < dt.Nb; ++ib) {
		for (int i = 0; i < 3; ++i) {
			if (!IsNull(dt.msh[ib].dt.cb)) {
				Cb(i, ib) = dt.msh[ib].dt.cb[i];
				cbnan = false;
			} else
				Cb(i, ib) = NaNDouble;
		}
	}
	if (!cbnan)
		mfile.Write("Cb", Cb);
		
	MatrixXd Cg(3, dt.Nb);
	bool cgnan = true;
	for (int ib = 0; ib < dt.Nb; ++ib) {
		for (int i = 0; i < 3; ++i) {
			if (!IsNull(dt.msh[ib].dt.cg)) {
				Cg(i, ib) = dt.msh[ib].dt.cg[i];
				cgnan = false;
			} else
				Cg(i, ib) = NaNDouble;
		}
	}
	if (!cgnan)
		mfile.Write("Cg", Cg);
	
	MatrixXd C0(3, dt.Nb);
	bool c0nan = true;
	for (int ib = 0; ib < dt.Nb; ++ib) {
		for (int i = 0; i < 3; ++i) {
			if (!IsNull(dt.msh[ib].dt.c0)) {
				C0(i, ib) = dt.msh[ib].dt.c0[i];
				c0nan = false;
			} else
				C0(i, ib) = NaNDouble;
		}
	}	
	if (!c0nan)	
		mfile.Write("C0", C0);
	
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
	mfile.WriteCell("Legend", legend);
}
