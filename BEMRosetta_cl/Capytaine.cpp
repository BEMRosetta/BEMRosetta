// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2024, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#include <NetCDF/NetCDF.h>

using namespace Upp;

String CapyNC_Load(const String &file, UArray<Hydro> &hydros, int &num) {
	try {
		BEM::Print("\n\n" + Format(t_("Loading '%s'"), file));
		BEM::Print("\n- " + S(t_("NC file")));
		
		String name = GetFileTitle(file);
	
		NetCDFFile cdf(file);

		double _g = cdf.GetDouble("g");
		
		UVector<double> _rho;
		cdf.GetDouble("rho", _rho);

		UVector<double> _h;
		cdf.GetDouble("water_depth", _h);
		for (auto &hh : _h) {
			if (hh == std::numeric_limits<double>::infinity())
				hh = -1;
		}
		
		UVector<double> _head;
		cdf.GetDouble("wave_direction", _head);
		for (auto &hh : _head)
			hh = ToDeg(hh);
		int Nh = _head.size();
		
		UVector<double> _w;
		cdf.GetDouble("omega", _w);
		int Nf = _w.size();

		//UVector<double> _T;
		//cdf.GetDouble("period", _T);

		int numaxisAB = 3;
		if (_rho.size() > 1)
			numaxisAB++;
		if (_h.size() > 1)
			numaxisAB++;
		
		MultiDimMatrixRowMajor<double> a;
		cdf.GetDouble("added_mass", a);
		ASSERT(numaxisAB == a.GetNumAxis());

		MultiDimMatrixRowMajor<double> b;
		cdf.GetDouble("radiation_damping", b);
		ASSERT(numaxisAB == b.GetNumAxis());
		
		int numaxisF = 4;
		if (_rho.size() > 1)
			numaxisF++;
		if (_h.size() > 1)
			numaxisF++;
		
		MultiDimMatrixRowMajor<double> sc;
		cdf.GetDouble("diffraction_force", sc);
		ASSERT(numaxisF == sc.GetNumAxis());
		
		MultiDimMatrixRowMajor<double> fk;
		cdf.GetDouble("Froude_Krylov_force", fk);
		ASSERT(numaxisF == fk.GetNumAxis());
					
		if (_rho.size() == 1) {		// Added to simplify handling later
			a.InsertAxis(0, 0);
			b.InsertAxis(0, 0);
			sc.InsertAxis(0, 0);
			fk.InsertAxis(0, 0);
		}
		if (_h.size() == 1) {
			a.InsertAxis(1, 0);
			b.InsertAxis(1, 0);
			sc.InsertAxis(1, 0);
			fk.InsertAxis(1, 0);
		}
		
		int Nb = a.size(3)/6;
		
		ASSERT(a.size(2) == Nf && a.size(3) == 6*Nb && a.size(4) == 6*Nb);
		ASSERT(b.size(2) == Nf && b.size(3) == 6*Nb && b.size(4) == 6*Nb);
				
		ASSERT(sc.size(2) == 2 && sc.size(3) == Nf && sc.size(4) == Nh && sc.size(5) == 6*Nb);
		ASSERT(fk.size(2) == 2 && fk.size(3) == Nf && fk.size(4) == Nh && fk.size(5) == 6*Nb);
		
		UArray<MatrixXd> M;
		if (cdf.ExistVar("inertia_matrix")) {
			MatrixXd _M;
			cdf.GetDouble("inertia_matrix", _M);
			M.SetCount(Nb);
			for (int ib = 0; ib < Nb; ++ib)
				M[ib] = _M.block(ib*6, ib*6, 6, 6);
		}
		
		UArray<MatrixXd> C;
		if (cdf.ExistVar("hydrostatic_stiffness")) {
			MatrixXd _C;
			cdf.GetDouble("hydrostatic_stiffness", _C);
			C.SetCount(Nb);
			for (int ib = 0; ib < Nb; ++ib)
				C[ib] = _C.block(ib*6, ib*6, 6, 6);
		}
		
		MatrixXd c0;
		if (cdf.ExistVar("rotation_center")) {
			cdf.GetDouble("rotation_center", c0);
			c0.transposeInPlace();
			ASSERT(c0.rows() == 3 && c0.cols() == Nb);
		} else 
			c0 = MatrixXd::Zero(3, Nb);
		
		String bodies = cdf.GetString("body_name");
		UVector<String> bds = Split(bodies, "+");			
		
		auto LoadAB = [&](const MultiDimMatrixRowMajor<double> &_a, UArray<UArray<VectorXd>> &a, int irho, int ih, int ib) {
			for (int r = 0; r < 6; ++r) 
				for (int c = 0; c < 6*Nb; ++c) 
					for (int iw = 0; iw < Nf; ++iw) 
						a[r + 6*ib][c](iw) = _a(irho, ih, iw, r, c);
		};
		auto LoadForce = [&](const MultiDimMatrixRowMajor<double> &_f, Hydro::Forces &f, int irho, int _ih, int ib) {
			for (int idf = 0; idf < 6; ++idf) 
				for (int ih = 0; ih < Nh; ++ih) 
					for (int iw = 0; iw < Nf; ++iw) 
						f.force[ih](iw, 6*ib + idf) = std::complex<double>(_f(irho, _ih, 0, iw, ih, idf + 6*ib), 
																		   -_f(irho, _ih, 1, iw, ih, idf + 6*ib));//-Imaginary to follow Wamit
		};
	
		num = _rho.size()*_h.size();
		
		for (int irho = 0; irho < _rho.size(); ++irho) {
			for (int ih = 0; ih < _h.size(); ++ih) {
				Hydro &hy = hydros.Add();
				
				hy.file = file;
				hy.name = name;
				if (!_rho.IsEmpty())
					hy.name + Format("_rho%.0f", _rho[irho]);
				if (!_h.IsEmpty())
					hy.name + Format("_h%.0f", _h[ih]);
				hy.dimen = true;
				hy.len = 1;
				hy.solver = Hydro::CAPYNC;
		
				hy.g = _g;
				hy.rho = _rho[irho];
				hy.h = _h[ih];
				
				hy.Nb = Nb;
				
				hy.dataFromW = true;
				
				hy.w = clone(_w);
				//hy.T = clone(_T);
				hy.Nf = Nf;
				hy.head = clone(_head);
				hy.Nh = Nh;
				//hy.M = clone(M);
				//hy.C = clone(C);
				//hy.c0 = clone(c0);
				hy.msh.SetCount(Nb);
				for (int ib = 0; ib < Nb; ++ib) {
					for (int i = 0; i < 3; ++i)
						hy.msh[ib].c0[i] = c0(i, ib);
					if (bds.size() > ib)
						hy.msh[ib].name = bds[ib];
					hy.msh[ib].M = clone(M[ib]);
					hy.msh[ib].C = clone(C[ib]);
				}
				
				hy.Initialize_AB(hy.A);
				hy.Initialize_AB(hy.B);
				hy.Initialize_Forces();
				
				for (int ib = 0; ib < Nb; ++ib) {
					LoadAB(a, hy.A, irho, ih, ib);
					LoadAB(b, hy.B, irho, ih, ib);
					LoadForce(sc, hy.sc, irho, ih, ib);
					LoadForce(fk, hy.fk, irho, ih, ib);
				}
				//hy.GetFexFromFscFfk();
				
				//hy.names.SetCount(Nb);
				//for (int i = 0; i < min(bds.size(), Nb); ++i)
				//	hy.names[i] = bds[i];
				
				/*hy.dof.Clear();	hy.dof.SetCount(hy.Nb, 0);
				for (int i = 0; i < hy.Nb; ++i)
					hy.dof[i] = 6;*/
			}
		}
	} catch (Exc e) {
		return e;
	}
	return String();
}


// diffraction_force		rexim x Nf x Nh x 6xNb


// omega					Nf
// period
// wave_direction			Nh
// inertia_matrix			6xNb x 6xNb
// hydrostatic_stiffness
// added_mass				2 x 2 x Nf x 6xNb x 6xNb
// radiation_damping
// diffraction_force		2 x 2 x re/im x Nf x Nh x 6xNb
// Froude_Krylov_force
// excitation_force (no esta siempre)

// A						[6*Nb][6*Nb][Nf]
// C						[Nb](6, 6)
// M
// Force					[Nh](Nf, 6*Nb) 
  


