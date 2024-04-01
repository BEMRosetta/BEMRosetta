// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2023, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#include <Hdf5/hdf5.h>

bool CapyNC::Load(String file, double) {
	hd().file = file;
	hd().name = GetFileTitle(file);
	hd().dimen = true;
	hd().len = 1;
	hd().code = Hydro::CAPYNC;
	hd().Nb = Null;
	hd().x_w = hd().y_w = 0;
	
	try {
		BEM::Print("\n\n" + Format(t_("Loading '%s'"), file));

		BEM::Print("\n- " + S(t_("NC file")));
		
		Load_NC();
		
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

void CapyNC::Load_NC() {
	String fileName = ForceExtSafer(hd().file, ".nc");
	
	Hdf5File hfile;
	
	hfile.Open(fileName, H5F_ACC_RDONLY);
	
	hd().dataFromW = true;
	
	hfile.GetDouble("omega", hd().w);
	hfile.GetDouble("period", hd().T);
	
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

