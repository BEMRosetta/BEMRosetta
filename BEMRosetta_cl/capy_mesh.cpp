// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"


String CapyBody::Load_NC(UArray<Body> &msh, String fileName, double g) {
	BEM bem;
	bem.g = g;	// This value is necessary, but discarded
	
	try {
		bem.LoadBEM(fileName, Null, false);
		
		Hydro &hy = bem.hydros[0];

		msh = pick(hy.dt.msh);
	} catch (Exc e) {
		return t_("Parsing error: ") + e;
	}
		
	return String();
}

