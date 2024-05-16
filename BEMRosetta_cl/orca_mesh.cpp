// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2024, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"

#ifdef PLATFORM_WIN32

String ORCAMesh::Load_OWR(UArray<Mesh> &msh, String fileName, double g, bool &y0z, bool &x0z) {
	y0z = x0z = false;
	
	BEM bem;
	bem.g = 9.8;	// This value is necessary, but discarded
	
	try {
		bem.LoadBEM(fileName, Null, false);
		
		Hydro &hy = bem.hydros[0];

		msh = pick(hy.dt.msh);
		y0z = hy.dt.symX;
		x0z = hy.dt.symY;
	} catch (Exc e) {
		return t_("Error loading owr: ") + e;
	}
		
	return String();
}

#endif