#include "BEMRosetta.h"


String MeshData::LoadStlTxt(String fileName, bool &isText) {
	isText = false;
	
	FileInLine in(fileName);
	if (!in.IsOpen()) 
		return Format(t_("Impossible to open file '%s'"), fileName);
	
	this->file = fileName;
	SetCode(MeshData::STL_TXT);
	
	try {
		String line;
		FieldSplit f(in);	
		
		line = ToLower(TrimBoth(in.GetLine()));
		if (!line.StartsWith("solid"))
			return Format(t_("[%d] 'solid' field not found"), in.GetLineNumber());
		
		isText = true;
		
	    while (!in.IsEof()) {	
			line = ToLower(TrimBoth(in.GetLine()));
			if (line.IsEmpty())
				continue;
			
			if (line.StartsWith("endsolid"))
				break;
			
			if (!line.StartsWith("facet normal"))
				return Format(t_("[%d] 'facet normal' field not found"), in.GetLineNumber());

			line = ToLower(TrimBoth(in.GetLine()));
			if (!line.StartsWith("outer loop"))
				return Format(t_("[%d] 'outer loop' field not found"), in.GetLineNumber());			

			int ids[5];
			for (int i = 0; i < 5; ++i) {
				f.Load(in.GetLine());
				String label = f.GetText(0);
				
				if (label == "vertex") {		
					if (i == 4)
						return Format(t_("[%d] Too much vertex in facet"), in.GetLineNumber());			
					Point3D node(f.GetDouble(1), f.GetDouble(2), f.GetDouble(3));
					mesh.nodes0 << node;
					ids[i] = mesh.nodes0.GetCount() - 1;
				} else if (label == "endloop") {
					if (i < 3)
						return Format(t_("[%d] Too few vertex in facet"), in.GetLineNumber());
					Panel &panel = mesh.panels.Add();
					panel.id[0] = ids[0];
					panel.id[1] = ids[1];
					panel.id[2] = ids[2];
					if (i == 3) 
						panel.id[3] = ids[0];
					else  
						panel.id[3] = ids[3];
				} else if (label == "endfacet") 
					break;
				else
					return Format(t_("[%d] Label '%s' not handled in facet"), in.GetLineNumber(), label);	
			}
	    }
	} catch (Exc e) {
		return e;
	}
	
	return String();
}
			
static void STLFacetTxtOut(FileOut &out, const Point3D &p0, const Point3D &p1, const Point3D &p2, 
						const Point3D &p3, const Vector3D &normal) {
	out << "facet normal " << normal.x << " " << normal.y << " " << normal.z << "\n";
	out << "   outer loop" << "\n";
	out << "      vertex " << p0.x << " " << p0.y << " " << p0.z << "\n";
	out << "      vertex " << p1.x << " " << p1.y << " " << p1.z << "\n";
	out << "      vertex " << p2.x << " " << p2.y << " " << p2.z << "\n";
	if (!IsNull(p3))
		out << "      vertex " << p3.x << " " << p3.y << " " << p3.z << "\n";	
	out << "   endloop" << "\n";
	out << "endfacet" << "\n";
}

void MeshData::SaveStlTxt(String fileName, const Vector<Panel> &panels, const Vector<Point3D> &nodes) {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'\n"), fileName));	
	
	bool forceTriangles = true;
	
	out << "solid BEMRosetta STL mesh export" << "\n";

	for (int i = 0; i < panels.GetCount(); ++i) {
		const Panel &panel = panels[i];
		const Point3D &p0 = nodes[panel.id[0]];
		const Point3D &p1 = nodes[panel.id[1]];
		const Point3D &p2 = nodes[panel.id[2]];
		if (forceTriangles) {
			STLFacetTxtOut(out, p0, p1, p2, Null, panel.normal0);
			if (!panel.IsTriangle()) {
				const Point3D &p3 = nodes[panel.id[3]];
				STLFacetTxtOut(out, p2, p3, p0, Null, panel.normal1);	
			}	
		} else {
			if (panel.IsTriangle()) 
				STLFacetTxtOut(out, p0, p1, p2, Null, panel.normalPaint);
			else {
				const Point3D &p3 = nodes[panel.id[3]];
				STLFacetTxtOut(out, p0, p1, p2, p3, panel.normalPaint);	
			}
		}
	}
	out << "endsolid BEMRosetta STL mesh export" << "\n";
}

String MeshData::LoadStlBin(String fileName) {
	FileInData in(fileName);
	if (!in.IsOpen())
		return Format(t_("Impossible to open file '%s'"), fileName);
	
	this->file = fileName;
	SetCode(MeshData::STL_BIN);
	
	try {	
		StringBuffer headerB(80);
    	in.Read(headerB, 80);
    	this->header = headerB;
    	
    	if (this->header.StartsWith("solid"))
    		return t_("Binary stl must not begin with 'solid' text");
    	
    	int32 numFacets = in.Read<int32>();
    	
    	while (!in.IsEof()) {
	    	Vector3D normal;
	    	normal.x = static_cast<double>(in.Read<float>());
	    	normal.y = static_cast<double>(in.Read<float>());
	    	normal.z = static_cast<double>(in.Read<float>());
	
			Panel &panel = mesh.panels.Add();

			Point3D &node0 = mesh.nodes0.Add();
			node0.x = static_cast<double>(in.Read<float>());
			node0.y = static_cast<double>(in.Read<float>());
			node0.z = static_cast<double>(in.Read<float>());
			panel.id[0] = mesh.nodes0.GetCount()-1;
			
			Point3D &node1 = mesh.nodes0.Add();
			node1.x = static_cast<double>(in.Read<float>());
			node1.y = static_cast<double>(in.Read<float>());
			node1.z = static_cast<double>(in.Read<float>());
			panel.id[1] = mesh.nodes0.GetCount()-1;
			
			Point3D &node2 = mesh.nodes0.Add();			
			node2.x = static_cast<double>(in.Read<float>());
			node2.y = static_cast<double>(in.Read<float>());
			node2.z = static_cast<double>(in.Read<float>());
			panel.id[2] = mesh.nodes0.GetCount()-1;
			panel.id[3] = panel.id[0];		
					
			int16 attributeByteCount = in.Read<int16>();
    	}
	} catch (Exc e) {
		return e;
	}
	return String();
}

static void STLFacetBinNodeOut(FileOutData &out, const Point3D &node) {
	out.Write(static_cast<float>(node.x));
	out.Write(static_cast<float>(node.y));
	out.Write(static_cast<float>(node.z));	
}

void MeshData::SaveStlBin(String fileName, const Vector<Panel> &panels, const Vector<Point3D> &nodes) {
	FileOutData out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'\n"), fileName));	

	String header = "BEMRosetta STL mesh export";
	header << String(' ', 80 - header.GetCount());
	out.Put64(header, 80);
	
	out.Write(static_cast<int32>(panels.GetCount()));
	
	for (int i = 0; i < panels.GetCount(); ++i) {
		const Panel &panel = panels[i];
		
		out.Write(static_cast<float>(panel.normal0.x));
		out.Write(static_cast<float>(panel.normal0.y));
		out.Write(static_cast<float>(panel.normal0.z));
		STLFacetBinNodeOut(out, nodes[panel.id[0]]);
		STLFacetBinNodeOut(out, nodes[panel.id[1]]);
		STLFacetBinNodeOut(out, nodes[panel.id[2]]);
		out.Write(static_cast<int16>(0));
		
		if (!panel.IsTriangle()) {
			out.Write(static_cast<float>(panel.normal1.x));
			out.Write(static_cast<float>(panel.normal1.y));
			out.Write(static_cast<float>(panel.normal1.z));
			STLFacetBinNodeOut(out, nodes[panel.id[2]]);
			STLFacetBinNodeOut(out, nodes[panel.id[3]]);
			STLFacetBinNodeOut(out, nodes[panel.id[0]]);	
			out.Write(static_cast<int16>(0));
		}
	}
}
