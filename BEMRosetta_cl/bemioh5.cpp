// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2023, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#include <Hdf5/hdf5.h>

String BemioH5::Load(String file, double) {
	dt.file = file;
	dt.name = GetFileTitle(file);
	dt.dimen = true;
	dt.len = 1;
	dt.solver = Hydro::BEMIO_H5;
	dt.Nb = Null;
	dt.x_w = dt.y_w = 0;
	
	try {
		BEM::Print("\n\n" + Format(t_("Loading '%s'"), file));

		BEM::Print("\n- " + S(t_("H5 file")));
		
		Load_H5();
		
		if (IsNull(dt.Nb))
			return t_("No data found");
	
	} catch (Exc e) {
		return e;
	}
	
	return String();
}

void BemioH5::Load_H5() {
	String fileName = ForceExtSafer(dt.file, ".h5");
	
	Hdf5File hfile;
	hfile.Open(fileName, H5F_ACC_RDONLY);
	
	for (dt.Nb = 0; hfile.ExistGroup(Format("body%d", dt.Nb+1)); dt.Nb++) 
		;	
	
	{
		hfile.ChangeGroup("bem_data");
		String str = ToLower(hfile.GetString("code"));
		if (str.Find("wamit") >= 0)
			dt.solver = Hydro::WAMIT;
		else if (str.Find("hams") >= 0)
			dt.solver = Hydro::HAMS;
		else if (str.Find("nemoh") >= 0)
			dt.solver = Hydro::NEMOH;
		else if (str.Find("aqwa") >= 0)
			dt.solver = Hydro::AQWA;
		else if (str.Find("capytaine") >= 0)
			dt.solver = Hydro::CAPYTAINE;
		
		hfile.UpGroup();
	}	
	{
		hfile.ChangeGroup("simulation_parameters");
		
		hfile.GetDouble("w", dt.w);
		dt.Nf = dt.w.size();
		hfile.GetDouble("wave_dir", dt.head);
		dt.Nh = dt.head.size();
		dt.rho = hfile.GetDouble("rho");
		dt.g = hfile.GetDouble("g");
		dt.dimen = hfile.GetDouble("scaled") != 0.;
		
		if (hfile.GetType("water_depth") == H5T_STRING) {
			String str = hfile.GetString("water_depth");
			if (str == "infinite")
				dt.h = -1;
			else
				throw Exc(Format("Unknown depth '%s'", str));
		} else 
			dt.h = hfile.GetDouble("water_depth");

		hfile.UpGroup();
	}
	
	dt.Ainf.setConstant(dt.Nb*6, dt.Nb*6, NaNDouble);
	
	Initialize_AB(dt.A);
	UArray<UArray<VectorXd>> a;		
	Initialize_AB(a);
	
	Initialize_AB(dt.B);
	UArray<UArray<VectorXd>> b;		
	Initialize_AB(b);
	
	Initialize_Forces();
	Hydro::Forces ex, sc, fk;	
	Initialize_Forces(ex);	
	Initialize_Forces(sc); 
	Initialize_Forces(fk);
	
	dt.msh.SetCount(dt.Nb);
	for (int ib = 0; ib < dt.Nb; ++ib) 
		dt.msh[ib].dt.C.setConstant(6, 6, 0);
	for (int ib = 0; ib < dt.Nb; ++ib) 
		dt.msh[ib].dt.M.setConstant(6, 6, 0);

	auto LoadComponentsAB = [&](UArray<UArray<VectorXd>> &a, int ib) {
		MatrixXd data;
		if (hfile.ChangeGroup("components")) {
			for (int r = 0; r < 6; ++r) {
				for (int c = 0; c < 6*dt.Nb; ++c) {	
					hfile.GetDouble(Format("%d_%d", r+1, c+1), data);
					if (data.rows() != dt.Nf && data.cols() != 2)
						throw Exc("Wrong data dimension reading added_mass components");
					for (int iw = 0; iw < dt.Nf; ++iw)
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
			for (int c = 0; c < 6*dt.Nb; ++c) 
				for (int iw = 0; iw < dt.Nf; ++iw) 
					a[r + 6*ib][c](iw) = d(r, c, iw);
	};
	
	auto LoadForce = [&](Hydro::Forces &f, int ib) {
		MatrixXd data;
		if (hfile.ChangeGroup("components")) {
			if (hfile.ChangeGroup("re")) {
				for (int idf = 0; idf < 6; ++idf) {
					for (int ih = 0; ih < dt.Nh; ++ih) {
						hfile.GetDouble(Format("%d_%d", idf+1, ih+1), data);
						if (data.rows() != dt.Nf && data.cols() != 2)
							throw Exc("Wrong data dimension reading force real component");
						for (int iw = 0; iw < dt.Nf; ++iw)
							f[ib][ih](iw, idf).real(data(iw, 1));
					}
				}
				hfile.UpGroup();	
			}
			if (hfile.ChangeGroup("im")) {
				for (int idf = 0; idf < 6; ++idf) {
					for (int ih = 0; ih < dt.Nh; ++ih) {
						hfile.GetDouble(Format("%d_%d", idf+1, ih+1), data);
						if (data.rows() != dt.Nf && data.cols() != 2)
							throw Exc("Wrong data dimension reading force imag component");
						for (int iw = 0; iw < dt.Nf; ++iw)
							f[ib][ih](iw, idf).imag(data(iw, 1));
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
				for (int ih = 0; ih < dt.Nh; ++ih) 
					for (int iw = 0; iw < dt.Nf; ++iw) 
						f[ib][ih](iw, idf).real(d(idf, ih, iw));		
		}
		if (hfile.ExistDataset("im")) {
			hfile.GetDouble("im", d);
			for (int idf = 0; idf < 6; ++idf) 
				for (int ih = 0; ih < dt.Nh; ++ih) 
					for (int iw = 0; iw < dt.Nf; ++iw) 
						f[ib][ih](iw, idf).imag(d(idf, ih, iw));			
		}
	};
		
	MatrixXd data;
	
	for (int ib = 0; ib < dt.Nb; ++ib) {
		if (hfile.ChangeGroup(Format("body%d", ib+1))) {
			if (hfile.ChangeGroup("properties")) {
				dt.msh[ib].dt.name = hfile.GetString("name");
				dt.msh[ib].dt.Vo = hfile.GetDouble("disp_vol");
				if (hfile.ExistDataset("cb")) {
					hfile.GetDouble("cb", data);
					if (data.rows() != 3 && data.cols() != 1)
						throw Exc("Wrong data dimension in cb");
					dt.msh[ib].dt.cb = data;
				}
				if (hfile.ExistDataset("cg")) {
					hfile.GetDouble("cg", data);
					if (data.rows() != 3 && data.cols() != 1)
						throw Exc("Wrong data dimension in cg");
					dt.msh[ib].dt.cg = data;
				}
				if (hfile.ExistDataset("c0")) {
					hfile.GetDouble("c0", data);
					if (data.rows() != 3 && data.cols() != 1)
						throw Exc("Wrong data dimension in c0");
					dt.msh[ib].dt.c0 = data;
				}
				hfile.UpGroup();	
			}
			if (hfile.ChangeGroup("hydro_coeffs")) {
				if (hfile.ChangeGroup("added_mass")) {
					if (hfile.ExistDataset("all")) 
						LoadAllAB(dt.A, ib);
					if (hfile.ExistGroup("components")) 
						LoadComponentsAB(a, ib);	
					
					if (hfile.ExistDataset("inf_freq")) {
						hfile.GetDouble("inf_freq", data);
						if (data.rows() != 6 && data.cols() != 6*dt.Nb)
							throw Exc("Wrong data dimension in inf_freq");
						dt.Ainf.block(ib*6, 0, 6, 6*dt.Nb) = data;	
					}
					hfile.UpGroup();	
				}
				if (hfile.ChangeGroup("radiation_damping")) {
					if (hfile.ExistDataset("all")) 
						LoadAllAB(dt.B, ib);
					if (hfile.ExistGroup("components")) 
						LoadComponentsAB(b, ib);
						
					hfile.UpGroup();	
				}
				if (hfile.ChangeGroup("excitation")) {
					LoadForce(dt.ex, ib);
					LoadForceAll(ex, ib);
					if (hfile.ChangeGroup("froude-krylov")) {
						LoadForce(dt.fk, ib);
						LoadForceAll(fk, ib);
						hfile.UpGroup();		
					}
					if (hfile.ChangeGroup("scattering")) {
						LoadForce(dt.sc, ib);
						LoadForceAll(sc, ib);
						hfile.UpGroup();	
					}
					hfile.UpGroup();	
				}
				hfile.GetDouble("linear_restoring_stiffness", dt.msh[ib].dt.C);
				
				hfile.UpGroup();	
			}
			hfile.UpGroup();	
		}
	}
	
	if (IsLoadedA())
		Compare_A(a);
	if (IsLoadedB())
		Compare_B(b);
	if (IsLoadedFex())
		Compare_F(dt.ex, ex, "Excitation forces");
	if (IsLoadedFsc())
		Compare_F(dt.sc, sc, "Scattering forces");
	if (IsLoadedFfk())
		Compare_F(dt.fk, fk, "Froude-Krylov forces");
	
	for (int ib = 0; ib < dt.msh.size(); ++ib)
		if (IsNull(dt.msh[ib].dt.c0))
			dt.msh[ib].dt.c0 = dt.msh[ib].dt.cg;
}

void BemioH5::Save(String file) const {
	String fileName = ForceExtSafer(file, ".h5");
	
	Hdf5File hfile;
	
	hfile.Create(fileName);
	
	{
		if (hfile.CreateGroup("bem_data", true)) {
			hfile.Set("code", GetCodeStr());
			hfile.UpGroup();
		}
	}
	{
		if (hfile.CreateGroup("simulation_parameters", true)) {
			
			VectorXd T = Get_T();
			hfile.Set("T", T).SetDescription("Wave frequencies").SetUnits("s");
			hfile.Set("w", dt.w).SetDescription("Wave periods").SetUnits("rad/s");
			hfile.Set("wave_dir", dt.head).SetDescription("Wave direction").SetUnits("deg");
			hfile.Set("rho", dt.rho).SetDescription("Water density").SetUnits("kg/m^3");
			hfile.Set("g", dt.g).SetDescription("Gravitational acceleration").SetUnits("m/s^2");
			hfile.Set("scaled", dt.dimen ? 1. : 0.).SetDescription("");
			if (dt.h < 0)
				hfile.Set("water_depth", "infinite");
			else
				hfile.Set("water_depth", dt.h);
			hfile.SetDescription("Water depth").SetUnits("m or infinite");
			hfile.UpGroup();
		}
	}

	auto SaveComponentsAB = [&](const UArray<UArray<VectorXd>> &a, int ib, String caption) {
		String str = Format("%s components as a function of frequency", caption);
		MatrixXd data(dt.Nf, 2);
		data.col(0) = Get_T();
		if (hfile.CreateGroup("components", true)) {
			for (int r = 0; r < 6; ++r) {
				for (int c = 0; c < 6*dt.Nb; ++c) {	
					for (int iw = 0; iw < dt.Nf; ++iw) {
						if (a[ib*6 + r][c].size() == dt.Nf)
							data(iw, 1) = a[ib*6 + r][c](iw);
						else
							data(iw, 1) = 0;
					}
					hfile.Set(Format("%d_%d", r+1, c+1), data).SetDescription(str);
				}
			}
			hfile.UpGroup();	
		}		
	};
	auto SaveAllAB = [&](const UArray<UArray<VectorXd>> &a, int ib, String caption) {
		MultiDimMatrixRowMajor<double> d(6, dt.Nb*6, dt.Nf);
		for (int r = 0; r < 6; ++r) {
			for (int c = 0; c < 6*dt.Nb; ++c) 
				for (int iw = 0; iw < dt.Nf; ++iw) {
					if (a[r + 6*ib][c].size() == dt.Nf) 
						d(r, c, iw) = a[r + 6*ib][c](iw);
					else
						d(r, c, iw) = 0;
				}
		}
		hfile.Set("all", d).SetDescription(caption);
	};
	
	auto SaveForce = [&](const Hydro::Forces &f, int ib, String caption) {
		MatrixXd data(dt.Nf, 2);
		data.col(0) = Get_T();
		if (hfile.CreateGroup("components", true)) {
			if (hfile.CreateGroup("re", true)) {
				String str = Format("Real component of %s force as a function of frequency", caption);
				for (int idf = 0; idf < 6; ++idf) {
					for (int ih = 0; ih < dt.Nh; ++ih) {
						for (int iw = 0; iw < dt.Nf; ++iw)
							data(iw, 1) = f[ib][ih](iw, idf).real();
						hfile.Set(Format("%d_%d", idf+1, ih+1), data).SetDescription(str);
					}
				}
				hfile.UpGroup();	
			}
			if (hfile.CreateGroup("im", true)) {
				String str = Format("Imaginary component of %s force as a function of frequency", caption);
				for (int idf = 0; idf < 6; ++idf) {
					for (int ih = 0; ih < dt.Nh; ++ih) {
						for (int iw = 0; iw < dt.Nf; ++iw)
							data(iw, 1) = f[ib][ih](iw, idf).imag();
						hfile.Set(Format("%d_%d", idf+1, ih+1), data).SetDescription(str);
					}
				}
				hfile.UpGroup();	
			}
			if (hfile.CreateGroup("mag", true)) {
				String str = Format("Magnitude of %s force as a function of frequency", caption);
				for (int idf = 0; idf < 6; ++idf) {
					for (int ih = 0; ih < dt.Nh; ++ih) {
						for (int iw = 0; iw < dt.Nf; ++iw)
							data(iw, 1) = abs(f[ib][ih](iw, idf));
						hfile.Set(Format("%d_%d", idf+1, ih+1), data).SetDescription(str);
					}
				}
				hfile.UpGroup();	
			}
			if (hfile.CreateGroup("phase", true)) {
				String str = Format("Phase of %s force as a function of frequency", caption);
				for (int idf = 0; idf < 6; ++idf) {
					for (int ih = 0; ih < dt.Nh; ++ih) {
						for (int iw = 0; iw < dt.Nf; ++iw)
							data(iw, 1) = arg(f[ib][ih](iw, idf));
						hfile.Set(Format("%d_%d", idf+1, ih+1), data).SetDescription(str).SetUnits("rad");
					}
				}
				hfile.UpGroup();	
			}
			hfile.UpGroup();	
		}
	};
	auto SaveForceAll = [&](const Hydro::Forces &f, int ib, String caption) {
		MultiDimMatrixRowMajor<double> d_m(6, dt.Nh, dt.Nf), d_p(6, dt.Nh, dt.Nf), d_r(6, dt.Nh, dt.Nf), d_i(6, dt.Nh, dt.Nf);		
		for (int idf = 0; idf < 6; ++idf) 
			for (int ih = 0; ih < dt.Nh; ++ih) 
				for (int iw = 0; iw < dt.Nf; ++iw) {
					const std::complex<double> &n = f[ib][ih](iw, idf);		
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
	
	for (int ib = 0; ib < dt.Nb; ++ib) {
		if (hfile.CreateGroup(Format("body%d", ib+1), true)) {
			if (hfile.CreateGroup("properties", true)) {
				MatrixXd mat;
				hfile.Set("body_number", ib+1.).SetDescription("Number of rigid body from the BEM simulation");
				hfile.Set("name", ~dt.msh[ib].dt.name);
				//if (dt.Vo.size() > ib)
					hfile.Set("disp_vol", dt.msh[ib].dt.Vo).SetDescription("Displaced volume").SetUnits("m^3");
				hfile.Set("dof", 6/*(double)dt.dof[ib]*/).SetDescription("Degrees of freedom");
				hfile.Set("dof_start", 1.);
				hfile.Set("dof_end", 6.);
				if (!IsNull(dt.msh[ib].dt.cb)) {
					mat = dt.msh[ib].dt.cb;
					mat = mat.transpose();
					hfile.Set("cb", mat).SetDescription("Centre of buoyancy").SetUnits("m");
				}
				if (!IsNull(dt.msh[ib].dt.cg)) {
					mat = dt.msh[ib].dt.cg;
					mat = mat.transpose();
					hfile.Set("cg", mat).SetDescription("Centre of gravity").SetUnits("m");
				}
				if (!IsNull(dt.msh[ib].dt.c0)) {
					mat = dt.msh[ib].dt.c0;
					mat = mat.transpose();
					hfile.Set("c0", mat).SetDescription("Body axis").SetUnits("m");
				}
				hfile.UpGroup();	
			}
			if (hfile.CreateGroup("hydro_coeffs", true)) {
				if (hfile.CreateGroup("added_mass", true)) {
					if (IsLoadedAinf()) {
						MatrixXd mat = dt.Ainf.block(ib*6, 0, 6, 6*dt.Nb);
						hfile.Set("inf_freq", mat).SetDescription("Infinite frequency added mass").SetUnits("kg");
					}
					if (IsLoadedA()) {
						SaveAllAB(dt.A, ib, "Added mass");
						SaveComponentsAB(dt.A, ib, "Added mass");
					}
					hfile.UpGroup();	
				}
				if (hfile.CreateGroup("radiation_damping", true)) {
					if (IsLoadedB()) {
						SaveAllAB(dt.B, ib, "Radiation damping");
						SaveComponentsAB(dt.B, ib, "Radiation damping");
					}						
					hfile.UpGroup();	
				}
				if (hfile.CreateGroup("excitation", true)) {
					if (IsLoadedFex()) {
						SaveForceAll(dt.ex, ib, "excitation");
						SaveForce(dt.ex, ib, "excitation");
					}
					if (hfile.CreateGroup("froude-krylov", true)) {
						if (IsLoadedFfk()) {
							SaveForceAll(dt.fk, ib, "froude-krylov");
							SaveForce(dt.fk, ib, "froude-krylov");
						}
						hfile.UpGroup();		
					}
					if (hfile.CreateGroup("scattering", true)) {
						if (IsLoadedFsc()) {
							SaveForceAll(dt.sc, ib, "scattering");
							SaveForce(dt.sc, ib, "scattering");
						}
						hfile.UpGroup();	
					}
					hfile.UpGroup();	
				}
				hfile.Set("linear_restoring_stiffness", dt.msh[ib].dt.C).SetDescription("Hydrostatic stiffness matrix");
				
				hfile.UpGroup();	
			}
			hfile.UpGroup();
		}
	}
}
