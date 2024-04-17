// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2024, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"

#ifdef PLATFORM_WIN32

String ORCAMesh::Load_OWR(UArray<Mesh> &mesh, String fileName, double g, bool &y0z, bool &x0z) {
	y0z = x0z = false;
	
	BEM bem;
	bem.g = 9.8;	// This value is necessary, but discarded
	
	try {
		bem.LoadBEM(fileName, Null, false);
		
		Hydro &hyd = bem.hydros[0].hd();

		mesh = pick(hyd.meshes);
		y0z = hyd.symX;
		x0z = hyd.symY;
	} catch (Exc e) {
		return t_("Parsing error: ") + e;
	}
		
	return String();
}

#endif