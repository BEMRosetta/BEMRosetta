// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"


String WamitBody::LoadDat(UArray<Body> &mesh, String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		return t_(F("Impossible to open '%s'", fileName));
	
	try {
		String line;
		LineParser f(in);
		f.IsSeparator = IsTabSpace;
		
		line = ToUpper(TrimBoth(in.GetLine()));
		if (!line.StartsWith("ZONE"))
			return in.Str() + "\n"  + t_("'ZONE' field not found");	// To detect Wamit format
	
		line.Replace("\"", "");
		line.Replace(" ", "");
		
		int pos;
		int T = Null;
		pos = line.FindAfter("T=");
		if (pos > 0) 
			T = ScanInt(line.Mid(pos));
		int I = Null;
		pos = line.FindAfter("I=");
		if (pos > 0) 
			I = ScanInt(line.Mid(pos));
		int J = Null;
		pos = line.FindAfter("J=");
		if (pos > 0) 
			J = ScanInt(line.Mid(pos));
		String F;
		pos = line.FindAfter("F=");
		if (pos > 0) 
			F = line.Mid(pos);
		
		Body &msh = mesh.Add();
		msh.dt.fileName = fileName;
		msh.dt.SetCode(Body::WAMIT_DAT);
	
		if (IsNull(T)) {
			while(!in.IsEof()) {
				int id0 = msh.dt.mesh.nodes.size();
				for (int i = 0; i < I*J; ++i) {
					line = in.GetLine();	
					f.Load(line);
					
					double x = f.GetDouble(0);	
					double y = f.GetDouble(1);	
					double z = f.GetDouble(2);	
						
					Point3D &node = msh.dt.mesh.nodes.Add();
					node.x = x;
					node.y = y;
					node.z = z;
				}
				for (int i = 0; i < I-1; ++i) {
					for (int j = 0; j < J-1; ++j) {
						Panel &panel = msh.dt.mesh.panels.Add();
						panel.id[0] = id0 + I*j     + i;
						panel.id[1] = id0 + I*j     + i+1;
						panel.id[2] = id0 + I*(j+1) + i;
						panel.id[3] = id0 + I*(j+1) + i+1;
					}
				}
				in.GetLine();
			}
		} else {
			for (int i = 0; i < I; ++i) {
				line = in.GetLine();	
				f.Load(line);
				
				double x = f.GetDouble(0);	
				double y = f.GetDouble(1);	
				double z = f.GetDouble(2);	
					
				Point3D &node = msh.dt.mesh.nodes.Add();
				node.x = x;
				node.y = y;
				node.z = z;
			}
			for (int i = 0; i < I/4; ++i) {
				line = in.GetLine();	
				f.Load(line);
				
				Panel &panel = msh.dt.mesh.panels.Add();
				for (int ii = 0; ii < 4; ++ii)
					panel.id[ii] = f.GetInt(ii) - 1;
			}
		}
	} catch (Exc e) {
		return t_("Parsing error: ") + e;
	}
	
	return String();
}


String WamitBody::LoadGGdf(LineParser &f, String folder, UArray<Body> &mesh, bool &y0z, bool &x0z) {
	try {
		f.GetLine();
		int Nb = f.GetInt(0);
		
		for (int ib = 0;ib <Nb; ++ib) {
			f.GetLine();
			String fileName = f.GetText(0);
			if (!FileExists(fileName)) {
				fileName = AFX(folder, fileName);
				if (!FileExists(fileName))
					return F(t_("File '%s' does not exist"), f.GetText(0));
			}
			double dummy_g;
			String ret = WamitBody::LoadGdf(mesh, fileName, y0z, x0z, dummy_g);
			if (!ret.IsEmpty())
				return ret;
			f.GetLine();
			double xbody = f.GetDouble(0);
			double ybody = f.GetDouble(1);
			double zbody = f.GetDouble(2);
			double rbody = f.GetDouble(3);
			
			Body &b = Last(mesh);
			if (y0z) {
				b.dt.mesh.DeployXSymmetry();
				b.dt.spline.DeployXSymmetry();
			}
			if (x0z) {
				b.dt.mesh.DeployYSymmetry();	
				b.dt.spline.DeployYSymmetry();	
			}
			b.Rotate(0, 0, -ToRad(rbody), 0, 0, 0);
			b.Translate(-xbody, -ybody, -zbody);
			
			f.GetLine();
		}
	} catch (Exc e) {
		return t_("Parsing error: ") + e;
	}	
	x0z = y0z = false;
	return String();
}

String WamitBody::LoadPot(UArray<Body> &mesh, String fileName, bool &y0z, bool &x0z, double &g) {
	Wamit wam;
	
	wam.Load_pot(fileName);
	for (int ib = 0; ib < wam.dt.msh.size(); ++ib)
		mesh << pick(wam.dt.msh[ib]);
	
	return String();	
}

	
String WamitBody::LoadGdf(UArray<Body> &_mesh, String fileName, bool &y0z, bool &x0z, double &g) {
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		return t_(F("Impossible to open '%s'", fileName));
	
	try {
		String line;
		LineParser f(in);	
		f.IsSeparator = [](int c)->int {
			if (c == '\t' || c == ' ' || c == '!' || c == '=')
				return true;
			return false;
		};
		
		in.GetLine();
		line = in.GetLine();	
		f.Load(line);
		
		double len = 1;
		g = Null;
		if (line.Find("ILOWHICSF") < 0) {
			g = f.GetDouble(1);
			if (g < 0)
				return t_("Wrong gravity in .gdf file");
			len = f.GetDouble(0);
			if (len == -1)
				return LoadGGdf(f, GetFileFolder(fileName), _mesh, y0z, x0z);
			if (len < 1)
				return t_("Wrong length scale in .gdf file");
		} else {
			if (f.GetInt(0) > 1)
				return t_("Only ILOWHICSF=0/1 is supported in .csf file");
		}
		line = in.GetLine();	
		f.Load(line);
		y0z = f.GetInt(0) != 0;
		x0z = f.GetInt(1) != 0;
		
		line = in.GetLine();	
		f.Load(line);
		int nPatches = f.GetInt(0);
		
		int igdef = 0;
		
		if (f.size() >= 2) {
			igdef = f.GetInt_nothrow(1);
			if (!IsNull(igdef)) {
				if (igdef == 0 || igdef == 1)
					;
				else if (igdef == 2)
					return t_(".gdf files represented by MultiSurf .ms2 files (IGDEF = 2) are not supported");
				else
					return t_(".gdf files represented by a special subrutine (IGDEF < 0  or > 2) are not supported");
			} else
				igdef = 0;		// There was text
		}
		if (nPatches < 1)
			return t_("Number of patches not found or zero in .gdf file");
		
		Body &body = _mesh.Add();
		body.dt.name = GetFileName(fileName);
				
		body.dt.SetCode(igdef == 0 ? Body::WAMIT_GDF : Body::WAMIT_GDF2);
		
		if (igdef == 0) {
			Surface &mesh = body.dt.mesh;
			int ids[4];
			bool npand = false;
			while(!in.IsEof()) {
				int ip = 0;
				f.GetLine();
				for (int i = 0; i < 4; ++i) {
					if (f.GetText(1) == "NPAND") { // Dipoles loaded as normal panels
						npand = true;
						break;
					}
					double x = f.GetDouble(0+ip);///len;	// Test05 from v6 manuals seems to indicate this
					double y = f.GetDouble(1+ip);///len;	
					double z = f.GetDouble(2+ip);///len;	
					
					bool found = false;
					for (int iin = 0; iin < mesh.nodes.size(); ++iin) {
						Point3D &node = mesh.nodes[iin];
						if (x == node.x && y == node.y && z == node.z) {
							ids[i] = iin;
							found = true;
							break;
						}
					}
					if (!found) {
						Point3D &node = mesh.nodes.Add();
						node.x = x;
						node.y = y;
						node.z = z;
						ids[i] = mesh.nodes.size() - 1;
					} 
					// To read multi point rows
					if (f.size() >= ip + 6 && f.IsDouble(ip + 3) && f.IsDouble(ip + 4) && f.IsDouble(ip + 5))
						ip += 3;
					else if (i < 3)
						f.GetLine();
				}
				if (!npand) {
					Panel &panel = mesh.panels.Add();
					for (int i = 0; i < 4; ++i)
						panel.id[i] = ids[i];
				}
				if (mesh.panels.size() == nPatches)
					break;
			} 
		} else if (igdef == 1) {
			SurfaceBSpline &spline = body.dt.spline;
			
			while(!in.IsEof()) {
				BSplinePatch patch;
				
				f.GetLine_discard_empty();
				patch.nug = f.GetInt(0);
				patch.nvg = f.GetInt(1);
				f.GetLine_discard_empty();
				patch.kug = f.GetInt(0);
				patch.kvg = f.GetInt(1);
				
				int nua = patch.nug + 2*patch.kug - 1;
				int nva = patch.nvg + 2*patch.kvg - 1;
				int nb = patch.GetNUBasis() * patch.GetNVBasis();
				
				patch.uKnots.SetCount(nua);
				for (int read = 0; read < nua;) {
					f.GetLine_discard_empty();
					for (int iline = 0; iline < f.size() && read < nua; ++iline)
						patch.uKnots[read++] = f.GetDouble(iline);
				}
				patch.vKnots.SetCount(nva);
				for (int read = 0; read < nva;) {
					f.GetLine_discard_empty();
					for (int iline = 0; iline < f.size() && read < nva; ++iline)
						patch.vKnots[read++] = f.GetDouble(iline);
				}
				patch.controlPoints.SetCount(nb);
				for (int j = 0; j < nb; j++) {
					f.GetLine_discard_empty();
					double x = f.GetDouble(0);
					double y = f.GetDouble(1);
					double z = f.GetDouble(2);
					patch.controlPoints[j] = Point3D(x, y, z);
				}
				spline.Add(pick(patch));
				
				if (spline.patches.size() == nPatches)
					break;
			}
		}
		//if (mesh.panels.size() != nPatches)
		//	return t_("Wrong number of patches in .gdf file");
	} catch (Exc e) {
		return t_("Parsing error: ") + e;
	}
		
	return String();
}

void WamitBody::SaveGdf(String fileName, const Surface &surf, double g, bool y0z, bool x0z, bool iscsf) {
	if (iscsf)
		ForceExt(fileName, ".csf");
	
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(F("Impossible to open '%s'", fileName));	
	
	const UVector<Panel> &panels = surf.panels;
	const UVector<Point3D> &nodes = surf.nodes;
	
	out << "BEMRosetta GDF mesh file export\n";
	if (!iscsf)
		out << F("  %12d   %12f 	ULEN GRAV\n", 1, g);
	else
		out << "  0                           	ILOWHICSF\n";
	out << F("  %12d   %12d 	ISX  ISY\n", y0z ? 1 : 0, x0z ? 1 : 0);
	out << F("  %12d\n", panels.size());
	for (int ip = 0; ip < panels.size(); ++ip) {
		for (int i = 0; i < 4; ++i) {
			int id = panels[ip].id[i];
			const Point3D &p = nodes[id]; 
			out << F("  % 014.7E   %0 14.7E   % 014.7E\n", p.x, p.y, p.z);
		}
	}	 
}

void WamitBody::SaveHST(String fileName, double rho, double g) const {
	Wamit::Save_hst_static(dt.C, fileName, rho, g);
}
	