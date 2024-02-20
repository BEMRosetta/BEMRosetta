// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2023, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#include <Hdf5/hdf5.h>

bool BemioH5::Load(String file, double) {
	hd().file = file;
	hd().name = GetFileTitle(file);
	hd().dimen = true;
	hd().len = 1;
	hd().code = Hydro::BEMIOH5;
	hd().Nb = Null;
	hd().x_w = hd().y_w = 0;
	
	try {
		BEM::Print("\n\n" + Format(t_("Loading '%s'"), file));

		BEM::Print("\n- " + S(t_("H5 file")));
		
		Load_H5();
		
		if (IsNull(hd().Nb))
			throw Exc(t_("No data found"));
	
		hd().dof.Clear();	hd().dof.SetCount(hd().Nb, 0);
		for (int i = 0; i < hd().Nb; ++i)
			hd().dof[i] = 6;
	} catch (Exc e) {
		//BEM::PrintError(Format("\n%s: %s", t_("Error"), e));
		hd().lastError = e;
		return false;
	}
	
	return true;
}

void BemioH5::Load_H5() {
	String fileName = ForceExt(hd().file, ".h5");
	
	Hdf5File hfile;
	
	if (!hfile.Open(fileName, H5F_ACC_RDONLY))
		throw Exc(Format(t_("file %s"), fileName) + "\n" + t_("File not found or blocked"));
	
	hd().dataFromW = true;
	
	for (hd().Nb = 0; hfile.ExistGroup(Format("body%d", hd().Nb+1)); hd().Nb++) 
		;	
	
	{
		hfile.ChangeGroup("bem_data");
		String str = ToLower(hfile.GetString("code"));
		if (str.Find("wamit") >= 0)
			hd().code = Hydro::WAMIT;
		else if (str.Find("hams") >= 0)
			hd().code = Hydro::HAMS;
		else if (str.Find("nemoh") >= 0)
			hd().code = Hydro::NEMOH;
		else if (str.Find("aqwa") >= 0)
			hd().code = Hydro::AQWA;
		else if (str.Find("capytaine") >= 0)
			hd().code = Hydro::CAPYTAINE;
		
		hfile.UpGroup();
	}	
	{
		hfile.ChangeGroup("simulation_parameters");
		
		hfile.GetDouble("T", hd().T);
		hfile.GetDouble("w", hd().w);
		hd().Nf = hd().w.size();
		hfile.GetDouble("wave_dir", hd().head);
		hd().Nh = hd().head.size();
		hd().rho = hfile.GetDouble("rho");
		hd().g = hfile.GetDouble("g");
		hd().dimen = hfile.GetDouble("scaled") != 0.;
		
		if (hfile.GetType("water_depth") == H5T_STRING) {
			String str = hfile.GetString("water_depth");
			if (str == "infinite")
				hd().h = -1;
			else
				throw Exc(Format("Unknown depth '%s'", str));
		} else 
			hd().h = hfile.GetDouble("water_depth");

		hfile.UpGroup();
	}
	
	hd().Ainf.setConstant(hd().Nb*6, hd().Nb*6, NaNDouble);
	
	hd().Initialize_AB(hd().A);
	UArray<UArray<VectorXd>> a;		
	hd().Initialize_AB(a);
	
	hd().Initialize_AB(hd().B);
	UArray<UArray<VectorXd>> b;		
	hd().Initialize_AB(b);
	
	hd().Initialize_Forces();
	Hydro::Forces ex, sc, fk;	
	hd().Initialize_Forces(ex);	
	hd().Initialize_Forces(sc); 
	hd().Initialize_Forces(fk);
	
	hd().names.SetCount(hd().Nb);
	hd().Vo.SetCount(hd().Nb);
	hd().dof.SetCount(hd().Nb);
	hd().cb.resize(3, hd().Nb);
	hd().cg.resize(3, hd().Nb);
	hd().c0.setConstant(3, hd().Nb, 0);
	hd().C.SetCount(hd().Nb);
	for (int ib = 0; ib < hd().Nb; ++ib) 
		hd().C[ib].setConstant(6, 6, 0);
	hd().M.SetCount(hd().Nb);
	for (int ib = 0; ib < hd().Nb; ++ib) 
		hd().M[ib].setConstant(6, 6, 0);

	auto LoadComponentsAB = [&](UArray<UArray<VectorXd>> &a, int ib) {
		MatrixXd data;
		if (hfile.ChangeGroup("components")) {
			for (int r = 0; r < 6; ++r) {
				for (int c = 0; c < 6*hd().Nb; ++c) {	
					hfile.GetDouble(Format("%d_%d", r+1, c+1), data);
					if (data.rows() != hd().Nf && data.cols() != 2)
						throw Exc("Wrong data dimension reading added_mass components");
					for (int iw = 0; iw < hd().Nf; ++iw)
						a[ib*6 + r][c](iw) = data(iw, 1);
				}
			}
			hfile.UpGroup();	
		}		
	};
	auto LoadAllAB = [&](UArray<UArray<VectorXd>> &a, int ib) {
		MultiDimMatrixRowMajor<double> d;
		hfile.GetDouble("all", d);
		for (int r = 0; r < 6; ++r) 
			for (int c = 0; c < 6*hd().Nb; ++c) 
				for (int iw = 0; iw < hd().Nf; ++iw) 
					a[r + 6*ib][c](iw) = d(r, c, iw);
	};
	
	auto LoadForce = [&](Hydro::Forces &f, int ib) {
		MatrixXd data;
		if (hfile.ChangeGroup("components")) {
			if (hfile.ChangeGroup("re")) {
				for (int idf = 0; idf < 6; ++idf) {
					for (int ih = 0; ih < hd().Nh; ++ih) {
						hfile.GetDouble(Format("%d_%d", idf+1, ih+1), data);
						if (data.rows() != hd().Nf && data.cols() != 2)
							throw Exc("Wrong data dimension reading force real component");
						for (int iw = 0; iw < hd().Nf; ++iw)
							f.force[ih](iw, 6*ib + idf).real(data(iw, 1));
					}
				}
				hfile.UpGroup();	
			}
			if (hfile.ChangeGroup("im")) {
				for (int idf = 0; idf < 6; ++idf) {
					for (int ih = 0; ih < hd().Nh; ++ih) {
						hfile.GetDouble(Format("%d_%d", idf+1, ih+1), data);
						if (data.rows() != hd().Nf && data.cols() != 2)
							throw Exc("Wrong data dimension reading force imag component");
						for (int iw = 0; iw < hd().Nf; ++iw)
							f.force[ih](iw, 6*ib + idf).imag(data(iw, 1));
					}
				}
				hfile.UpGroup();	
			}
			hfile.UpGroup();	
		}
	};
	auto LoadForceAll = [&](Hydro::Forces &f, int ib) {
		MultiDimMatrixRowMajor<double> d;
		if (hfile.ExistDataset("re")) {
			hfile.GetDouble("re", d);
			for (int idf = 0; idf < 6; ++idf) 
				for (int ih = 0; ih < hd().Nh; ++ih) 
					for (int iw = 0; iw < hd().Nf; ++iw) 
						f.force[ih](iw, 6*ib + idf).real(d(idf, ih, iw));		
		}
		if (hfile.ExistDataset("im")) {
			hfile.GetDouble("im", d);
			for (int idf = 0; idf < 6; ++idf) 
				for (int ih = 0; ih < hd().Nh; ++ih) 
					for (int iw = 0; iw < hd().Nf; ++iw) 
						f.force[ih](iw, 6*ib + idf).imag(d(idf, ih, iw));			
		}
	};
		
	MatrixXd data;
	
	for (int ib = 0; ib < hd().Nb; ++ib) {
		if (hfile.ChangeGroup(Format("body%d", ib+1))) {
			if (hfile.ChangeGroup("properties")) {
				hd().names[ib] = hfile.GetString("name");
				hd().Vo[ib] = hfile.GetDouble("disp_vol");
				hd().dof[ib] = int(hfile.GetDouble("dof"));
				hfile.GetDouble("cb", data);
				if (data.rows() != 3 && data.cols() != 1)
					throw Exc("Wrong data dimension in cb");
				hd().cb.col(ib) = data;
				hfile.GetDouble("cg", data);
				if (data.rows() != 3 && data.cols() != 1)
					throw Exc("Wrong data dimension in cg");
				hd().cg.col(ib) = data;
			
				hfile.UpGroup();	
			}
			if (hfile.ChangeGroup("hydro_coeffs")) {
				if (hfile.ChangeGroup("added_mass")) {
					LoadAllAB(hd().A, ib);
					LoadComponentsAB(a, ib);	
						
					hfile.GetDouble("inf_freq", data);
					if (data.rows() != 6 && data.cols() != 6*hd().Nb)
						throw Exc("Wrong data dimension in inf_freq");
					hd().Ainf.block(ib*6, 0, 6, 6*hd().Nb) = data;	
						
					hfile.UpGroup();	
				}
				if (hfile.ChangeGroup("radiation_damping")) {
					LoadAllAB(hd().B, ib);
					LoadComponentsAB(b, ib);
						
					hfile.UpGroup();	
				}
				if (hfile.ChangeGroup("excitation")) {
					LoadForce(hd().ex, ib);
					LoadForceAll(ex, ib);
					if (hfile.ChangeGroup("froude-krylov")) {
						LoadForce(hd().fk, ib);
						LoadForceAll(fk, ib);
						hfile.UpGroup();		
					}
					if (hfile.ChangeGroup("scattering")) {
						LoadForce(hd().sc, ib);
						LoadForceAll(sc, ib);
						hfile.UpGroup();	
					}
					hfile.UpGroup();	
				}
				hfile.GetDouble("linear_restoring_stiffness", hd().C[ib]);
				
				hfile.UpGroup();	
			}
			hfile.UpGroup();	
		}
	}
	
	if (hd().IsLoadedA())
		hd().Compare_A(a);
	if (hd().IsLoadedB())
		hd().Compare_B(b);
	if (hd().IsLoadedFex())
		hd().Compare_F(hd().ex, ex, "Excitation");
	if (hd().IsLoadedFsc())
		hd().Compare_F(hd().sc, sc, "Scattering");
	if (hd().IsLoadedFfk())
		hd().Compare_F(hd().fk, fk, "Froude-Krylov");
}

bool BemioH5::Save(String file) {
	String fileName = ForceExt(file, ".h5");
	
	Hdf5File hfile;
	
	if (!hfile.Create(fileName))
		return false;
	
	{
		if (hfile.CreateGroup("bem_data", true)) {
			hfile.Set("code", hd().GetCodeStr());
			hfile.UpGroup();
		}
	}
	{
		if (hfile.CreateGroup("simulation_parameters", true)) {
			
			hfile.Set("T", hd().T).SetDescription("Wave frequencies").SetUnits("s");
			hfile.Set("w", hd().w).SetDescription("Wave periods").SetUnits("rad/s");
			hfile.Set("wave_dir", hd().head).SetDescription("Wave direction").SetUnits("deg");
			hfile.Set("rho", hd().rho).SetDescription("Water density").SetUnits("kg/m^3");
			hfile.Set("g", hd().g).SetDescription("Gravitational acceleration").SetUnits("m/s^2");
			hfile.Set("scaled", hd().dimen ? 1. : 0.).SetDescription("");
			if (hd().h < 0)
				hfile.Set("water_depth", "infinite");
			else
				hfile.Set("water_depth", hd().h);
			hfile.SetDescription("Water depth").SetUnits("m or infinite");
			hfile.UpGroup();
		}
	}

	auto SaveComponentsAB = [&](const UArray<UArray<VectorXd>> &a, int ib, String caption) {
		String str = Format("%s components as a function of frequency", caption);
		MatrixXd data(hd().Nf, 2);
		data.col(0) = hd().Get_T();
		if (hfile.CreateGroup("components", true)) {
			for (int r = 0; r < 6; ++r) {
				for (int c = 0; c < 6*hd().Nb; ++c) {	
					for (int iw = 0; iw < hd().Nf; ++iw)
						data(iw, 1) = a[ib*6 + r][c](iw);
					hfile.Set(Format("%d_%d", r+1, c+1), data).SetDescription(str);
				}
			}
			hfile.UpGroup();	
		}		
	};
	auto SaveAllAB = [&](const UArray<UArray<VectorXd>> &a, int ib, String caption) {
		MultiDimMatrixRowMajor<double> d(6, hd().Nb*6, hd().Nf);
		for (int r = 0; r < 6; ++r) 
			for (int c = 0; c < 6*hd().Nb; ++c) 
				for (int iw = 0; iw < hd().Nf; ++iw) 
					d(r, c, iw) = a[r + 6*ib][c](iw);
		hfile.Set("all", d).SetDescription(caption);
	};
	
	auto SaveForce = [&](Hydro::Forces &f, int ib, String caption) {
		MatrixXd data(hd().Nf, 2);
		data.col(0) = hd().Get_T();
		if (hfile.CreateGroup("components", true)) {
			if (hfile.CreateGroup("re", true)) {
				String str = Format("Real component of %s force as a function of frequency", caption);
				for (int idf = 0; idf < 6; ++idf) {
					for (int ih = 0; ih < hd().Nh; ++ih) {
						for (int iw = 0; iw < hd().Nf; ++iw)
							data(iw, 1) = f.force[ih](iw, 6*ib + idf).real();
						hfile.Set(Format("%d_%d", idf+1, ih+1), data).SetDescription(str);
					}
				}
				hfile.UpGroup();	
			}
			if (hfile.CreateGroup("im", true)) {
				String str = Format("Imaginary component of %s force as a function of frequency", caption);
				for (int idf = 0; idf < 6; ++idf) {
					for (int ih = 0; ih < hd().Nh; ++ih) {
						for (int iw = 0; iw < hd().Nf; ++iw)
							data(iw, 1) = f.force[ih](iw, 6*ib + idf).imag();
						hfile.Set(Format("%d_%d", idf+1, ih+1), data).SetDescription(str);
					}
				}
				hfile.UpGroup();	
			}
			if (hfile.CreateGroup("mag", true)) {
				String str = Format("Magnitude of %s force as a function of frequency", caption);
				for (int idf = 0; idf < 6; ++idf) {
					for (int ih = 0; ih < hd().Nh; ++ih) {
						for (int iw = 0; iw < hd().Nf; ++iw)
							data(iw, 1) = abs(f.force[ih](iw, 6*ib + idf));
						hfile.Set(Format("%d_%d", idf+1, ih+1), data).SetDescription(str);
					}
				}
				hfile.UpGroup();	
			}
			if (hfile.CreateGroup("phase", true)) {
				String str = Format("Phase of %s force as a function of frequency", caption);
				for (int idf = 0; idf < 6; ++idf) {
					for (int ih = 0; ih < hd().Nh; ++ih) {
						for (int iw = 0; iw < hd().Nf; ++iw)
							data(iw, 1) = arg(f.force[ih](iw, 6*ib + idf));
						hfile.Set(Format("%d_%d", idf+1, ih+1), data).SetDescription(str).SetUnits("rad");
					}
				}
				hfile.UpGroup();	
			}
			hfile.UpGroup();	
		}
	};
	auto SaveForceAll = [&](Hydro::Forces &f, int ib, String caption) {
		MultiDimMatrixRowMajor<double> d_m(6, hd().Nh, hd().Nf), d_p(6, hd().Nh, hd().Nf), d_r(6, hd().Nh, hd().Nf), d_i(6, hd().Nh, hd().Nf);		
		for (int idf = 0; idf < 6; ++idf) 
			for (int ih = 0; ih < hd().Nh; ++ih) 
				for (int iw = 0; iw < hd().Nf; ++iw) {
					const std::complex<double> &n = f.force[ih](iw, 6*ib + idf);		
					d_r(idf, ih, iw) = n.real();		
					d_i(idf, ih, iw) = n.imag();		
					d_m(idf, ih, iw) = abs(n);		
					d_p(idf, ih, iw) = arg(n);		
				}
				
		hfile.Set("re",    d_r).SetDescription(Format("Real component of %s force", caption));		
		hfile.Set("im",    d_i).SetDescription(Format("Imaginary component of %s force", caption));
		hfile.Set("mag",   d_m).SetDescription(Format("Magnitude of %s force", caption));
		hfile.Set("phase", d_p).SetDescription(Format("Phase angle of %s force", caption));
	};
	
	for (int ib = 0; ib < hd().Nb; ++ib) {
		if (hfile.CreateGroup(Format("body%d", ib+1), true)) {
			if (hfile.CreateGroup("properties", true)) {
				MatrixXd mat;
				hfile.Set("body_number", ib+1.).SetDescription("Number of rigid body from the BEM simulation");
				hfile.Set("name", ~hd().names[ib]);
				hfile.Set("disp_vol", hd().Vo[ib]).SetDescription("Displaced volume").SetUnits("m^3");
				hfile.Set("dof", (double)hd().dof[ib]).SetDescription("Degrees of freedom");
				hfile.Set("dof_start", 1.);
				hfile.Set("dof_end", 6.);
				mat = hd().cb.col(ib);
				hfile.Set("cb", mat).SetDescription("Centre of bouyancy").SetUnits("m");
				mat = hd().cg.col(ib);
				hfile.Set("cg", mat).SetDescription("Centre of gravity").SetUnits("m");
				
				hfile.UpGroup();	
			}
			if (hfile.CreateGroup("hydro_coeffs", true)) {
				if (hfile.CreateGroup("added_mass", true)) {
					if (hd().IsLoadedAinf()) {
						MatrixXd mat = hd().Ainf.block(ib*6, 0, 6, 6*hd().Nb);
						hfile.Set("inf_freq", mat).SetDescription("Infinite frequency added mass").SetUnits("kg");
					}
					if (hd().IsLoadedA()) {
						SaveAllAB(hd().A, ib, "Added mass");
						SaveComponentsAB(hd().A, ib, "Added mass");
					}
					hfile.UpGroup();	
				}
				if (hfile.CreateGroup("radiation_damping", true)) {
					if (hd().IsLoadedB()) {
						SaveAllAB(hd().B, ib, "Radiation damping");
						SaveComponentsAB(hd().B, ib, "Radiation damping");
					}						
					hfile.UpGroup();	
				}
				if (hfile.CreateGroup("excitation", true)) {
					if (hd().IsLoadedFex()) {
						SaveForceAll(hd().ex, ib, "excitation");
						SaveForce(hd().ex, ib, "excitation");
					}
					if (hfile.CreateGroup("froude-krylov", true)) {
						if (hd().IsLoadedFfk()) {
							SaveForceAll(hd().fk, ib, "froude-krylov");
							SaveForce(hd().fk, ib, "froude-krylov");
						}
						hfile.UpGroup();		
					}
					if (hfile.CreateGroup("scattering", true)) {
						if (hd().IsLoadedFsc()) {
							SaveForceAll(hd().sc, ib, "scattering");
							SaveForce(hd().sc, ib, "scattering");
						}
						hfile.UpGroup();	
					}
					hfile.UpGroup();	
				}
				hfile.Set("linear_restoring_stiffness", hd().C[ib]).SetDescription("Hydrostatic stiffness matrix");
				
				hfile.UpGroup();	
			}
			hfile.UpGroup();
		}
	}
	
	return true;
}
