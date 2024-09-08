// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2024, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#include <NetCDF/NetCDF.h>

using namespace Upp;

String CapyNC_Load(const String &file, UArray<Hydro> &hydros, int &num) {
	num = 0;
	
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

		int numaxisAB = 3;
		if (_rho.size() > 1)
			numaxisAB++;
		if (_h.size() > 1)
			numaxisAB++;
		
		MultiDimMatrixRowMajor<double> a;
		cdf.GetDouble("added_mass", a);
		if (!(numaxisAB == a.GetNumAxis()))
			throw Exc(t_("Wrong dimension in added_mass"));

		MultiDimMatrixRowMajor<double> b;
		cdf.GetDouble("radiation_damping", b);
		if (!(numaxisAB == b.GetNumAxis()))
			throw Exc(t_("Wrong dimension in radiation_damping"));
		
		int numaxisF = 4;
		if (_rho.size() > 1)
			numaxisF++;
		if (_h.size() > 1)
			numaxisF++;
		
		MultiDimMatrixRowMajor<double> sc;
		cdf.GetDouble("diffraction_force", sc);
		if (!(numaxisF == sc.GetNumAxis()))
			throw Exc(t_("Wrong dimension in diffraction_force"));
		
		MultiDimMatrixRowMajor<double> fk;
		cdf.GetDouble("Froude_Krylov_force", fk);
		if (!(numaxisF == fk.GetNumAxis()))
			throw Exc(t_("Wrong dimension in Froude_Krylov_force"));
			

		MultiDimMatrixRowMajor<double> rao;
		if (cdf.ExistVar("RAO")) {
			cdf.GetDouble("RAO", rao);
			if (!(numaxisF == rao.GetNumAxis()))
				throw Exc(t_("Wrong dimension in RAO"));
		}
		if (_rho.size() == 1) {		// Added to simplify handling later
			a.InsertAxis(0, 1);
			b.InsertAxis(0, 1);
			sc.InsertAxis(0, 1);
			fk.InsertAxis(0, 1);
			if (!rao.IsEmpty())
				rao.InsertAxis(0, 1);
		}
		if (_h.size() == 1) {
			a.InsertAxis(1, 1);
			b.InsertAxis(1, 1);
			sc.InsertAxis(1, 1);
			fk.InsertAxis(1, 1);
			if (!rao.IsEmpty())
				rao.InsertAxis(1, 1);
		}
		
		int Nb = a.size(3)/6;
		
		if (!(a.size(2) == Nf && a.size(3) == 6*Nb && a.size(4) == 6*Nb))
			throw Exc(t_("Wrong dimension in a"));
		if (!(b.size(2) == Nf && b.size(3) == 6*Nb && b.size(4) == 6*Nb))
			throw Exc(t_("Wrong dimension in b"));
				
		if (!(sc.size(2) == 2 && sc.size(3) == Nf && sc.size(4) == Nh && sc.size(5) == 6*Nb))
			throw Exc(t_("Wrong dimension in sc"));
		if (!(fk.size(2) == 2 && fk.size(3) == Nf && fk.size(4) == Nh && fk.size(5) == 6*Nb))
			throw Exc(t_("Wrong dimension in fk"));
		if (!rao.IsEmpty() && !(rao.size(2) == 2 && rao.size(3) == Nf && rao.size(4) == Nh && rao.size(5) == 6*Nb))
			throw Exc(t_("Wrong dimension in rao"));
		
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
			if (!(c0.rows() == 3 && c0.cols() == Nb))
				throw Exc(t_("Wrong dimension in c0"));
		} else 
			c0 = MatrixXd::Zero(3, Nb);

		MatrixXd cg;
		if (cdf.ExistVar("center_of_mass")) {
			cdf.GetDouble("center_of_mass", cg);
			cg.transposeInPlace();
			if (!(cg.rows() == 3 && cg.cols() == Nb))
				throw Exc(t_("Wrong dimension in cg"));
		} else 
			cg = MatrixXd::Zero(3, Nb);
					
		String bodies = cdf.GetString("body_name");
		UVector<String> bds = Split(bodies, "+");	
		
		
		int numPan = 0;
		UVector<int> bodyPan;
		UVector<int> bodyIdEachPan;
		
		MultiDimMatrixRowMajor<double> dofDefinition;
		if (cdf.ExistVar("dof_definition")) {
			cdf.GetDouble("dof_definition", dofDefinition);
			if (!(3 == dofDefinition.GetNumAxis()))
				throw Exc(t_("Wrong dimension in dofDefinition"));
			if (!(dofDefinition.size(0) == 6*Nb && dofDefinition.size(2) == 3))
				throw Exc(t_("Wrong dimension in dofDefinition 2"));
			
			numPan = dofDefinition.size(1);
			bodyPan.SetCount(numPan);
			bodyIdEachPan.SetCount(numPan);
			UVector<int> bodyIdPan;
			bodyIdPan.SetCount(Nb, 0);
			
			auto GetBody = [=](int ip, const MultiDimMatrixRowMajor<double> &d)->int {
				for (int ib = 0; ib < Nb; ++ib) {
					for (int idof = 0; idof < 6; ++idof)
						for (int i = 0; i < 3; ++i)	
							if (d(idof + 6*ib, ip, i) != 0)
								return ib;
				}
				NEVER();	return -1;
			};
			for (int ip = 0; ip < numPan; ++ip) {
				int ib = bodyPan[ip] = GetBody(ip, dofDefinition);	// Assigns each panel to a body
				bodyIdEachPan[ip] = bodyIdPan[ib];					// Sets the id of each panel in its body
				bodyIdPan[ib]++;									// Increments the actual id of next panel in each body
			}
		}
				
		MultiDimMatrixRowMajor<double> pan;
		MultiDimMatrixRowMajor<double> rad_press, dif_press, inc_press;

		if (cdf.ExistVar("mesh_vertices")) {
			cdf.GetDouble("mesh_vertices", pan);
			if (!(3 == pan.GetNumAxis()))
				throw Exc(t_("Wrong dimension in mesh_vertices"));		
		
			numPan = pan.size(0);
			if (!(pan.size(1) == 4 && pan.size(2) == 3))
				throw Exc(t_("Wrong dimension in mesh_vertices 2"));
		}
		if (cdf.ExistVar("radiation_pressure")) {
			cdf.GetDouble("radiation_pressure", rad_press);
			if (!(numaxisF == rad_press.GetNumAxis()))
				throw Exc(t_("Wrong dimension in radiation_pressure"));		

			if (_rho.size() == 1)  		// Added to simplify handling later
				rad_press.InsertAxis(0, 1);
			if (_h.size() == 1) 
				rad_press.InsertAxis(1, 1);
		
			if (!(rad_press.size(2) == 2 && rad_press.size(3) == Nf && rad_press.size(4) == 6*Nb && rad_press.size(5) == numPan))
				throw Exc(t_("Wrong dimension in radiation_pressure 2"));
		}
		if (cdf.ExistVar("diffraction_pressure")) {
			cdf.GetDouble("diffraction_pressure", dif_press);
			if (!(numaxisF == dif_press.GetNumAxis()))
				throw Exc(t_("Wrong dimension in diffraction_pressure"));		
		
			if (_rho.size() == 1) 		// Added to simplify handling later
				dif_press.InsertAxis(0, 1);
			if (_h.size() == 1) 
				dif_press.InsertAxis(1, 1);
			
			if (!(dif_press.size(2) == 2 && dif_press.size(3) == Nf && dif_press.size(4) == Nh && dif_press.size(5) == numPan))
				throw Exc(t_("Wrong dimension in diffraction_pressure 2"));
		}
		if (cdf.ExistVar("incident_pressure")) {
			cdf.GetDouble("incident_pressure", inc_press);
			if (!(numaxisF == inc_press.GetNumAxis()))
				throw Exc(t_("Wrong dimension in incident_pressure"));		
		
			if (_rho.size() == 1) 		// Added to simplify handling later
				inc_press.InsertAxis(0, 1);				
			if (_h.size() == 1) 
				inc_press.InsertAxis(1, 1);
			if (!(inc_press.size(2) == 2 && inc_press.size(3) == Nf && inc_press.size(4) == Nh && inc_press.size(5) == numPan))
				throw Exc(t_("Wrong dimension in incident_pressure 2"));
		}

		auto LoadAB = [&](const MultiDimMatrixRowMajor<double> &_a, UArray<UArray<VectorXd>> &a, int irho, int ih) {
			for (int r = 0; r < 6*Nb; ++r) 
				for (int c = 0; c < 6*Nb; ++c) 
					for (int iw = 0; iw < Nf; ++iw) 
						a[r][c](iw) = _a(irho, ih, iw, r, c);
		};
		auto LoadForce = [&](const MultiDimMatrixRowMajor<double> &_f, Hydro::Forces &f, int irho, int _ih, int ib) {
			for (int idf = 0; idf < 6; ++idf) 
				for (int ih = 0; ih < Nh; ++ih) 
					for (int iw = 0; iw < Nf; ++iw) 
						f[ib][ih](iw, idf) = std::complex<double>(_f(irho, _ih, 0, iw, ih, idf + 6*ib), 
																 -_f(irho, _ih, 1, iw, ih, idf + 6*ib));	//-Imaginary to follow Wamit
		};
		
		for (int irho = 0; irho < _rho.size(); ++irho) {
			for (int ih = 0; ih < _h.size(); ++ih) {
				Hydro &hy = hydros.Add();
				num++;
				
				hy.dt.file = file;
				hy.dt.name = name;
				if (_rho.size() > 1)
					hy.dt.name += Format("_rho%.0f", _rho[irho]);
				if (_h.size() > 1)
					hy.dt.name += Format("_h%.0f", _h[ih]);
				hy.dt.dimen = true;
				hy.dt.len = 1;
				hy.dt.solver = Hydro::CAPY_NC;
		
				hy.dt.x_w = hy.dt.y_w = 0;
				
				hy.dt.g = _g;
				hy.dt.rho = _rho[irho];
				hy.dt.h = _h[ih];
				
				hy.dt.Nb = Nb;
				
				//hy.dt.dataFromW = true;
				
				hy.dt.w = clone(_w);

				hy.dt.Nf = Nf;
				hy.dt.head = clone(_head);
				hy.dt.Nh = Nh;

				hy.dt.msh.SetCount(Nb);
				for (int ib = 0; ib < Nb; ++ib) {
					for (int i = 0; i < 3; ++i)
						hy.dt.msh[ib].dt.c0[i] = c0(i, ib);
					for (int i = 0; i < 3; ++i)
						hy.dt.msh[ib].dt.cg[i] = cg(i, ib);
					if (bds.size() > ib)
						hy.dt.msh[ib].dt.name = bds[ib];
					if (M.size() > ib)
						hy.dt.msh[ib].dt.M = clone(M[ib]);
					if (C.size() > ib)
						hy.dt.msh[ib].dt.C = clone(C[ib]);
				}
				
				hy.Initialize_AB(hy.dt.A);
				hy.Initialize_AB(hy.dt.B);
								
				LoadAB(a, hy.dt.A, irho, ih);
				LoadAB(b, hy.dt.B, irho, ih);
				
				hy.Initialize_Forces();
				if (!rao.IsEmpty())
					hy.Initialize_Forces(hy.dt.rao);
	
				for (int ib = 0; ib < Nb; ++ib) {
					LoadForce(sc, hy.dt.sc, irho, ih, ib);
					LoadForce(fk, hy.dt.fk, irho, ih, ib);
					if (!rao.IsEmpty())
						LoadForce(rao, hy.dt.rao, irho, ih, ib);
				}
				
				if (pan.size() > 0) {
					for (int ipall = 0; ipall < numPan; ++ipall) {
						int ib = bodyPan[ipall];
						Surface &msh = hy.dt.msh[ib].dt.mesh;
						Panel &p = msh.panels.Add();
				
						for (int i = 0; i < 4; ++i) {
							Point3D pnt(pan(ipall, i, 0), pan(ipall, i, 1), pan(ipall, i, 2));
							p.id[i] = FindAdd(msh.nodes, pnt);
						}
					}
				}
				
				if (rad_press.size() > 0) {
					hy.Initialize_PotsRad();
					
					for (int ipall = 0; ipall < numPan; ++ipall) {
						int ib = bodyPan[ipall];
						int ip = bodyIdEachPan[ipall];
						for (int ifr = 0; ifr < Nf; ++ifr) {
							double rw = hy.dt.rho*sqr(hy.dt.w[ifr]);							
							for (int idf = 0; idf < 6; ++idf) {
								double re = rad_press(irho, ih, 0, ifr, idf + 6*ib, ipall);
								double im = rad_press(irho, ih, 1, ifr, idf + 6*ib, ipall);
								hy.dt.pots_rad[ib][ip][idf][ifr] += std::complex<double>(-re, im)/rw; // p = iρωΦ
							}
						}
					}
				}
				if (inc_press.size() > 0) {
					hy.Initialize_PotsIncDiff(hy.dt.pots_inc);
					
					for (int ipall = 0; ipall < numPan; ++ipall) {
						int ib = bodyPan[ipall];
						int ip = bodyIdEachPan[ipall];
						for (int ifr = 0; ifr < Nf; ++ifr) {
							double rw = hy.dt.rho*hy.dt.w[ifr];
							for (int ihead = 0; ihead < Nh; ++ihead) {
								double re = inc_press(irho, ih, 0, ifr, ihead, ipall);
								double im = inc_press(irho, ih, 1, ifr, ihead, ipall);
								hy.dt.pots_inc[ib][ip][ihead][ifr] += std::complex<double>(im, re)/rw; // p = iρωΦ
							}
						}
					}
				}
				if (dif_press.size() > 0) {
					hy.Initialize_PotsIncDiff(hy.dt.pots_dif);
					
					for (int ipall = 0; ipall < numPan; ++ipall) {
						int ib = bodyPan[ipall];
						int ip = bodyIdEachPan[ipall];
						for (int ifr = 0; ifr < Nf; ++ifr) {
							double rw = hy.dt.rho*hy.dt.w[ifr];
							for (int ihead = 0; ihead < Nh; ++ihead) {
								double re = dif_press(irho, ih, 0, ifr, ihead, ipall);
								double im = dif_press(irho, ih, 1, ifr, ihead, ipall);
								hy.dt.pots_dif[ib][ip][ihead][ifr] += std::complex<double>(im, re)/rw; // p = iρωΦ
							}
						}
					}
				}
			}
		}
	} catch (Exc e) {
		return e;
	}
	return String();
}


//	void Initialize_PotsRad();
//	void Initialize_PotsIncDiff(UArray<UArray<UArray<UArray<std::complex<double>>>>> &pots);
	
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
// excitation_force 

// A						[6*Nb][6*Nb][Nf]
// C						[Nb](6, 6)
// M
// Force					[Nh](Nf, 6*Nb) 
  


