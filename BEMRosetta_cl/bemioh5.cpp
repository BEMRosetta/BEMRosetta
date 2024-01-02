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
	
	try {
		BEM::Print("\n\n" + Format(t_("Loading '%s'"), file));

		BEM::Print("\n- " + S(t_("H5 file")));
		if (!Load_H5()) 
			BEM::PrintWarning(S(": ** H5 file ") + t_("Not found") + "**");
		
		if (IsNull(hd().Nb))
			return false;
	
		hd().dof.Clear();	hd().dof.SetCount(hd().Nb, 0);
		for (int i = 0; i < hd().Nb; ++i)
			hd().dof[i] = 6;
	} catch (Exc e) {
		BEM::PrintError(Format("\n%s: %s", t_("Error"), e));
		//hd().lastError = e;
		return false;
	}
	
	return true;
}

bool BemioH5::Load_H5() {
	String fileName = ForceExt(hd().file, ".h5");
	
	Hdf5File hfile;
	
	if (!hfile.Open(fileName, H5F_ACC_RDONLY))
		return false;
	
	hd().dataFromW = true;
	
	for (hd().Nb = 0; hfile.ExistGroup(Format("body%d", hd().Nb+1)); hd().Nb++) 
		;	
	
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
	hd().Initialize_AB(hd().B);
	hd().Initialize_Forces();
	
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
	
	MatrixXd data;
	
	for (int ib = 0; ib < hd().Nb; ++ib) {
		if (hfile.ChangeGroup(Format("body%d", ib+1))) {
			if (hfile.ChangeGroup("properties")) {
				hd().names[ib] = hfile.GetString("name");
				hd().Vo[ib] = hfile.GetDouble("disp_vol");
				hd().dof[ib] = hfile.GetDouble("dof");
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
					hfile.GetDouble("inf_freq", data);
					if (data.rows() != 6 && data.cols() != 6*hd().Nb)
						throw Exc("Wrong data dimension in inf_freq");
					hd().Ainf.block(ib*6, 0, 6, 6*hd().Nb) = data;	
		
					LoadComponentsAB(hd().A, ib);
						
					hfile.UpGroup();	
				}
				if (hfile.ChangeGroup("radiation_damping")) {
					LoadComponentsAB(hd().B, ib);
						
					hfile.UpGroup();	
				}
				if (hfile.ChangeGroup("excitation")) {
					LoadForce(hd().ex, ib);
					if (hfile.ChangeGroup("froude-krylov")) {
						LoadForce(hd().fk, ib);
						hfile.UpGroup();		
					}
					if (hfile.ChangeGroup("scattering")) {
						LoadForce(hd().sc, ib);
						hfile.UpGroup();	
					}
					hfile.UpGroup();	
				}
				hfile.UpGroup();	
			}
			hfile.UpGroup();	
		}
	}
	
	return true;
}

bool BemioH5::Save(String file) {
	String fileName = ForceExt(file, ".h5");
	
	Hdf5File hfile;
	
	if (!hfile.Create(fileName))
		return false;
	
	hfile.CreateGroup("bem_data");
	
	{
		hfile.CreateGroup("simulation_parameters", true);
		
		hfile.Set("T", hd().T).SetDescription("Wave frequencies").SetUnits("s");
		hfile.Set("w", hd().w).SetDescription("Wave periods").SetUnits("rad/s");
		hfile.Set("wave_dir", hd().head);
		hfile.Set("rho", hd().rho);
		hfile.Set("g", hd().g);
		hfile.Set("scaled", hd().dimen ? 1. : 0.);
		if (hd().h < 0)
			hfile.Set("water_depth", "infinite");
		else
			hfile.Set("water_depth", hd().h);
		hfile.UpGroup();
	}

	auto SaveComponentsAB = [&](const UArray<UArray<VectorXd>> &a, int ib) {
		MatrixXd data(hd().Nf, 2);
		data.col(0) = hd().Get_T();
		if (hfile.CreateGroup("components", true)) {
			for (int r = 0; r < 6; ++r) {
				for (int c = 0; c < 6*hd().Nb; ++c) {	
					for (int iw = 0; iw < hd().Nf; ++iw)
						data(iw, 1) = a[ib*6 + r][c](iw);
					hfile.Set(Format("%d_%d", r+1, c+1), data);
				}
			}
			hfile.UpGroup();	
		}		
	};
	
	auto SaveForce = [&](Hydro::Forces &f, int ib) {
		MatrixXd data(hd().Nf, 2);
		data.col(0) = hd().Get_T();
		if (hfile.CreateGroup("components", true)) {
			if (hfile.CreateGroup("re", true)) {
				for (int idf = 0; idf < 6; ++idf) {
					for (int ih = 0; ih < hd().Nh; ++ih) {
						for (int iw = 0; iw < hd().Nf; ++iw)
							data(iw, 1) = f.force[ih](iw, 6*ib + idf).real();
						hfile.Set(Format("%d_%d", idf+1, ih+1), data);
					}
				}
				hfile.UpGroup();	
			}
			if (hfile.CreateGroup("im", true)) {
				for (int idf = 0; idf < 6; ++idf) {
					for (int ih = 0; ih < hd().Nh; ++ih) {
						for (int iw = 0; iw < hd().Nf; ++iw)
							data(iw, 1) = f.force[ih](iw, 6*ib + idf).imag();
						hfile.Set(Format("%d_%d", idf+1, ih+1), data);
					}
				}
				hfile.UpGroup();	
			}
			if (hfile.CreateGroup("mag", true)) {
				for (int idf = 0; idf < 6; ++idf) {
					for (int ih = 0; ih < hd().Nh; ++ih) {
						for (int iw = 0; iw < hd().Nf; ++iw)
							data(iw, 1) = abs(f.force[ih](iw, 6*ib + idf));
						hfile.Set(Format("%d_%d", idf+1, ih+1), data);
					}
				}
				hfile.UpGroup();	
			}
			if (hfile.CreateGroup("phase", true)) {
				for (int idf = 0; idf < 6; ++idf) {
					for (int ih = 0; ih < hd().Nh; ++ih) {
						for (int iw = 0; iw < hd().Nf; ++iw)
							data(iw, 1) = arg(f.force[ih](iw, 6*ib + idf));
						hfile.Set(Format("%d_%d", idf+1, ih+1), data);
					}
				}
				hfile.UpGroup();	
			}
			hfile.UpGroup();	
		}
	};
	
	
	for (int ib = 0; ib < hd().Nb; ++ib) {
		if (hfile.CreateGroup(Format("body%d", ib+1), true)) {
			if (hfile.CreateGroup("properties", true)) {
				MatrixXd mat;
				hfile.Set("body_number", ib+1.);
				hfile.Set("name", hd().names[ib]);
				hfile.Set("disp_vol", hd().Vo[ib]);
				hfile.Set("dof", (double)hd().dof[ib]);
				hfile.Set("dof_start", 1.);
				hfile.Set("dof_end", 6.);
				mat = hd().cb.col(ib);
				hfile.Set("cb", mat);
				mat = hd().cg.col(ib);
				hfile.Set("cg", mat);
				
				hfile.UpGroup();	
			}
			if (hfile.CreateGroup("hydro_coeffs", true)) {
				if (hfile.CreateGroup("added_mass", true)) {
					if (hd().IsLoadedAinf()) {
						MatrixXd mat = hd().Ainf.block(ib*6, 0, 6, 6*hd().Nb);
						hfile.Set("inf_freq", mat);
					}
					if (hd().IsLoadedA())
						SaveComponentsAB(hd().A, ib);
						
					hfile.UpGroup();	
				}
				if (hfile.CreateGroup("radiation_damping", true)) {
					if (hd().IsLoadedB())
						SaveComponentsAB(hd().B, ib);
						
					hfile.UpGroup();	
				}
				if (hfile.CreateGroup("excitation", true)) {
					if (hd().IsLoadedFex())
						SaveForce(hd().ex, ib);
					if (hfile.CreateGroup("froude-krylov", true)) {
						if (hd().IsLoadedFfk())
							SaveForce(hd().fk, ib);
						hfile.UpGroup();		
					}
					if (hfile.CreateGroup("scattering", true)) {
						if (hd().IsLoadedFsc())
							SaveForce(hd().sc, ib);
						hfile.UpGroup();	
					}
					hfile.UpGroup();	
				}
				hfile.UpGroup();	
			}
			hfile.UpGroup();
		}
	}
	
	return true;
}
