// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2025, the BEMRosetta author and contributors
#include "BEMRosetta.h"
//#include "BEMRosetta_int.h"
#include <MatIO/matio.h>

void Matlab::Save(String file) const {
	String fileName = ForceExt(file, ".mat");
	
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
		MultiDimMatrix<std::complex<double>> Fex(6*dt.Nh, dt.Nf, dt.Nb);
		for (int idf = 0; idf < 6; ++idf) 
			for (int ih = 0; ih < dt.Nh; ++ih) 
				for (int ib = 0; ib < dt.Nb; ++ib) 		
					for (int ifr = 0; ifr < dt.Nf; ++ifr) 
						Fex(idf + 6*ih, ifr, ib) = dt.ex[ib][ih](ifr, idf);
		mfile.Set("Fex", Fex);				
	}
	if (IsLoadedFfk()) {
		MultiDimMatrix<std::complex<double>> Ffk(6*dt.Nh, dt.Nf, dt.Nb);
		for (int idf = 0; idf < 6; ++idf) 
			for (int ih = 0; ih < dt.Nh; ++ih) 
				for (int ib = 0; ib < dt.Nb; ++ib) 		
					for (int ifr = 0; ifr < dt.Nf; ++ifr) 
						Ffk(idf + 6*ih, ifr, ib) = dt.fk[ib][ih](ifr, idf);
		mfile.Set("Ffk", Ffk);	
	}
	if (IsLoadedFsc()) {
		MultiDimMatrix<std::complex<double>> Fsc(6*dt.Nh, dt.Nf, dt.Nb);
		for (int idf = 0; idf < 6; ++idf) 
			for (int ih = 0; ih < dt.Nh; ++ih) 
				for (int ib = 0; ib < dt.Nb; ++ib) 		
					for (int ifr = 0; ifr < dt.Nf; ++ifr) 
						Fsc(idf + 6*ih, ifr, ib) = dt.sc[ib][ih](ifr, idf);
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
					"- Fex = Ffk + Fsc: Excitation force (imaginary numbers) (6*nHeadings, nFrequencies, nBodies)",
					"- Ffk: Froude-Krylov force (imaginary numbers) (6*nHeadings, nFrequencies, nBodies)",
					"- Fsc: Diffracion scattering force (imaginary numbers) (6*nHeadings, nFrequencies, nBodies)",
					"- Kh: Hydrostatic stiffness matrix (6*nBodies, 6*nBodies)",
					"- Mass: Inertia matrix (6*nBodies, 6*nBodies)",
					"- Vol: Submerged volume (1, nBodies)",
					"- Cg: Centres of gravity (3, nBodies)",
					"- Cb: Centres of buoyancy (3, nBodies)",
					"- C0: Centres of body axis (3, nBodies)",
					"Comments:",
					"- Global axis (wave origin) is at 0, 0",
					"- The incident wave system follows Wamit criteria"};
	mfile.SetCell("Legend", legend);
}
