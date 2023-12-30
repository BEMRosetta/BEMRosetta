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
	
	if (!hfile.Load(fileName))
		return false;
	
	hfile.ChangeGroup("simulation_parameters");
	
	Eigen::VectorXd vect;

	hfile.GetDouble("T", vect);
	Copy(vect, hd().T);
	
	hfile.GetDouble("w", vect);
	Copy(vect, hd().w);	
	
	hd().rho = hfile.GetDouble("rho");
	
	if (hfile.GetType("water_depth") == H5T_STRING) {
		String str = hfile.GetString("water_depth");
		if (str == "infinite")
			hd().h = -1;
		else
			throw Exc(Format("Unknown depth '%s'", str));
	} else 
		hd().h = hfile.GetDouble("water_depth");
				
	

	
	return true;
}
