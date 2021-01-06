#include "BEMRosetta.h"

String MeshData::LoadDatWamit(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		return Format(t_("Impossible to open file '%s'"), fileName);
	
	this->fileName = fileName;
	SetCode(MeshData::WAMIT_DAT);
	
	try {
		String line;
		FieldSplit f(in);
		f.IsSeparator = IsTabSpace;
		
		line = ToUpper(TrimBoth(in.GetLine()));
		if (!line.StartsWith("ZONE"))
			return in.Str() + "\n"  + t_("'ZONE' field not found");
	
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
		
		if (IsNull(T)) {
			while(!in.IsEof()) {
				int id0 = mesh.nodes.size();
				for (int i = 0; i < I*J; ++i) {
					line = in.GetLine();	
					f.Load(line);
					
					double x = f.GetDouble(0);	
					double y = f.GetDouble(1);	
					double z = f.GetDouble(2);	
						
					Point3D &node = mesh.nodes.Add();
					node.x = x;
					node.y = y;
					node.z = z;
				}
				for (int i = 0; i < I-1; ++i) {
					for (int j = 0; j < J-1; ++j) {
						Panel &panel = mesh.panels.Add();
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
					
				Point3D &node = mesh.nodes.Add();
				node.x = x;
				node.y = y;
				node.z = z;
			}
			for (int i = 0; i < I/4; ++i) {
				line = in.GetLine();	
				f.Load(line);
				
				Panel &panel = mesh.panels.Add();
				for (int ii = 0; ii < 4; ++ii)
					panel.id[ii] = f.GetInt(ii) - 1;
			}
		}
	} catch (Exc e) {
		return t_("Parsing error: ") + e;
	}
	
	return String();
}
	
String MeshData::LoadGdfWamit(String fileName, bool &y0z, bool &x0z) {
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		return Format(t_("Impossible to open file '%s'"), fileName);
	
	SetCode(MeshData::WAMIT_GDF);
	
	try {
		String line;
		FieldSplit f(in);	
		f.IsSeparator = IsTabSpace;
		
		in.GetLine();
		line = in.GetLine();	
		f.Load(line);
		double len = f.GetDouble(0);
		if (len < 1)
			return t_("Wrong length scale in .gdf file");
					
		line = in.GetLine();	
		f.Load(line);
		y0z = f.GetInt(0) != 0;
		x0z = f.GetInt(1) != 0;
		
		line = in.GetLine();	
		f.Load(line);
		int nPatches = f.GetInt(0);
		if (nPatches < 1)
			return t_("Number of patches not found in .gdf file");
				
		mesh.Clear();
		
		while(!in.IsEof()) {
			int ids[4];
			bool npand = false;
			for (int i = 0; i < 4; ++i) {
				line = in.GetLine();	
				f.Load(line);
				
				if (f.GetText(1) == "NPAND") { // Dipoles loaded as normal panels
					npand = true;
					break;
				}
				double x = f.GetDouble(0)*len;	
				double y = f.GetDouble(1)*len;	
				double z = f.GetDouble(2)*len;	
				
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
			}
			if (!npand) {
				Panel &panel = mesh.panels.Add();
				for (int i = 0; i < 4; ++i)
					panel.id[i] = ids[i];
			}
			if (mesh.panels.size() == nPatches)
				break;
		}
		//if (mesh.panels.size() != nPatches)
		//	return t_("Wrong number of patches in .gdf file");
	} catch (Exc e) {
		return t_("Parsing error: ") + e;
	}
		
	return String();
}

void MeshData::SaveGdfWamit(String fileName, const Surface &surf, double g, bool y0z, bool x0z) {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'"), fileName));	
	
	const Vector<Panel> &panels = surf.panels;
	const Vector<Point3D> &nodes = surf.nodes;
	
	out << "BEMRosetta GDF mesh file export\n";
	out << Format("  %12d   %12f 	ULEN GRAV\n", 1, g);
	out << Format("  %12d   %12d 	ISX  ISY\n", y0z ? 1 : 0, x0z ? 1 : 0);
	out << Format("  %12d\n", panels.size());
	for (int ip = 0; ip < panels.size(); ++ip) {
		for (int i = 0; i < 4; ++i) {
			int id = panels[ip].id[i];
			const Point3D &p = nodes[id]; 
			out << Format("  % 014.7E   %0 14.7E   % 014.7E\n", p.x, p.y, p.z);
		}
	}	 
}
	