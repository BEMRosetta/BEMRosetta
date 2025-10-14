// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2025, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"

#include <Core/Core.h>
#include <SysInfo/SysInfo.h>
#include <Surface/Surface.h>
#include <STEM4U/Mooring.h>


using namespace Upp;


void Mooring::Jsonize(JsonIO &json) {
	json
		("lineTypes", lineTypes)
		("lineProperties", lineProperties)
		("connections", connections)
		("vessels", vessels)
		("depth", depth)
	;	
}

void Mooring::LineType::Jsonize(JsonIO &json) {
	if (json.IsLoading()) {
		mass = diameter = bl = 0;
		ea = "0";
		bazeta = "-1";
	}
	json
		("name", name)
		("mass", mass)
		("diameter", diameter)
		("BL", bl)
		("EA", ea)
	;
}

void Mooring::LineProperty::Jsonize(JsonIO &json) {
	if (json.IsLoading()) {
		length = 0;
	}
	json
		("name", name)
		("nameType", nameType)
		("length", length)
		("from", from)
		("to", to)
		("numseg", numseg)
	;
}

void Mooring::Connection::Jsonize(JsonIO &json) {
	if (json.IsLoading())
		x = y = z = 0;
	json
		("name", name)
		("where", where)
		("x", x)
		("y", y)
		("z", z)
	;
}

void Mooring::Vessel::Jsonize(JsonIO &json) {
	json
		("name", name)
	;
}

bool Mooring::Load(String file) {
	try {
		if (!LoadFromJsonFile(*this, ~file))
			return false;
	} catch (const Exc &e) {
		BEM::PrintError(e);
		return false;
	}	
	return true;
}

bool Mooring::Save(String file) {
	try {
		if (!StoreAsJsonFile(*this, file, true)) {
			BEM::PrintError(Format("Problem loading %s file", file));
			return false;
		}
	} catch (const Exc &e) {
		BEM::PrintError(e);
		return false;
	}
	
	return true;
}

bool Mooring::LoadMoordyn(String file) {
	UVector<LineType> lT;
	UVector<LineProperty> lP;	
	UVector<Connection> connect;
	
	FileInLine in(file);
	if (!in.IsOpen()) 
		return false;
	
	depth = Null;
	
	LineParser f(in);
	f.IsSeparator = IsTabSpace;
	
	double zmin = std::numeric_limits<double>::max();
	char status = '\0';
	while (!f.IsEof()) {
		if (f.GetLineNumber() > 1000 && (lT.IsEmpty() && lP.IsEmpty() && connect.IsEmpty()))
			return false;
		String line = Trim(f.GetLine());
		String lowline = ToLower(line);
		if (lowline.StartsWith("---")) {
			if (lowline.Find("type") > 0) {
				status = 't';
				f.GetLine(2);
			} else if (lowline.Find("point") > 0) {
				status = 'p';
				f.GetLine(2);
			} else if (lowline.Find("line") > 0) {
				status = 'l';
				f.GetLine(2);
			} else if (lowline.Find("option") > 0) {
				status = 'o';
			} else
				status = '\0';
		} else {
			if (status == 't') {
				LineType &line = lT.Add();
				line.name = f.GetText(0);
				line.diameter = f.GetDouble(1);
				line.mass = f.GetDouble(2);
				line.ea = f.GetText(3);
				line.bazeta = f.GetText(4);
				line.bl = Null;
			} else if (status == 'p') {
				Connection &conn = connect.Add();
				conn.name = f.GetText(0);
				conn.where = ToLower(f.GetText(1));
				const Upp::Index<String> attachments = {"fixed", "anchor", "coupled", "vessel", "free", "point"};
				if (attachments.Find(conn.where) < 0) {
					if (!conn.where.StartsWith("body"))
	 					throw Exc(Format(t_("Unknown type of attachment '%s' in point '%s'"), f.GetText(1), conn.name)); 
				}
				conn.x = f.GetDouble(2);
				conn.y = f.GetDouble(3);
				conn.z = f.GetDouble(4);
				zmin = min(zmin, conn.z);
			} else if (status == 'l') {
				LineProperty &prop = lP.Add();
				prop.name = f.GetText(0);
				prop.nameType = f.GetText(1);
				prop.from = f.GetText(2);
				prop.to = f.GetText(3);
				prop.length = f.GetDouble(4);
				prop.numseg = f.GetInt(5);
			} else if (status == 'o') {
				String param = ToLower(f.GetText(1));
				if (f.size() >= 2) {
					if (param == "dtm")
						dtM = f.GetDouble(0);
					else if (param == "kbot")
						kbot = f.GetDouble(0);
					else if (param == "cbot")
						cbot = f.GetDouble(0);
					else if (param == "dtic")
						dtIC = f.GetDouble(0);
					else if (param == "tmaxic")
						TmaxIC = f.GetDouble(0);
					else if (param == "cdscaleic")
						CdScaleIC = f.GetDouble(0);
					else if (param == "threshic")
						threshIC = f.GetDouble(0);
					else if (param == "wtrdpth")
						depth = f.GetDouble(0);
				}
			}
		}
	}
	if (IsNull(depth))
		depth = abs(zmin);
	
	if (lT.IsEmpty() && lP.IsEmpty() && connect.IsEmpty())
		return false;
	
	lineTypes = pick(lT);
	lineProperties = pick(lP);	
	connections = pick(connect);
	
	return true;
}

bool Mooring::SaveMoordyn(String file, int mooring, bool fairleads, bool anchors) {
	double ratioSegsToPlot = 0.1;
	int minSegsToPlot = 10;
	
	FileOut out(file);
	if (!out.IsOpen()) 
		return false;	
	
	out <<  "--------------------- MoorDyn Input File ------------------------------------\n"
			"BEMRosetta generated basic mooring file\n"
			"FALSE    Echo      - echo the input file data (flag)\n";
			
	out <<  "----------------------- LINE TYPES ------------------------------------------\n"
		<<  Format(" %12=s %7=s %8=s %13=s %8=s %6=s %6=s %6=s %6=s %6=s\n", "Name", "Diam", "MassDen", "EA ", "BA/-zeta", "EI ", "Cd ", "Ca ", "CdAx", "CaAx")
		<<  Format(" %12=s %7=s %8=s %13=s %8=s %6=s %6=s %6=s %6=s %6=s\n", "(-) ", "(m) ", "(kg/m) ", "(N)", "(N-s/-) ", "(-)", "(-)", "(-)", "(-) ", "(-) ");
	for (int i = 0; i < lineTypes.size(); ++i) {
		const LineType &line = lineTypes[i];
		out << Format(" %12=s %7.4f %8.1f %13s %8s %6.1f %6.4f %6.1f %6.4f %6.1f\n", line.name, line.diameter, line.mass, line.ea, line.bazeta, 0, 1.3325, 1, 0.2032, 0.5);
	}
			
	out <<  "----------------------- POINTS ----------------------------------------------\n"
		<<  Format(" %3=s %8=s %12=s %12=s %12=s %6=s %6=s %6=s %6=s\n", "ID ", "Type", " X ", " Y ", " Z ", " M  ", "  V  ", " CdA ", "CA ")
		<<  Format(" %3=s %8=s %12=s %12=s %12=s %6=s %6=s %6=s %6=s\n", "(-)", "(-) ", "(m)", "(m)", "(m)", "(kg)", "(m^3)", "(m^2)", "(-)");
	for (int i = 0; i < connections.size(); ++i) {		
		const Connection &conn = connections[i];
		out << Format(" %3d %8=s %12.4f %12.4f %12.4f %6.1f %6.1f %6.1f %6.1f\n", i+1, InitCaps(conn.where), conn.x, conn.y, conn.z, 0, 0, 0, 0);
	}
	out <<  "----------------------- LINES -----------------------------------------------\n"
		<<	Format(" %3=s %12=s %8=s %8=s %12=s %8=s %8=s\n", "ID ", "LineType", "AttachA", "AttachB", "UnstrLen", "NumSegs", "Outputs")
		<<	Format(" %3=s %12=s %8=s %8=s %12=s %8=s %8=s\n", "(-)", " (name) ", "  (#)  ", "  (#)  ", "   (m)  ", "  (-)  ", "  (-)  ");
	
	String sfair, sanch;
	for (int i = 0; i < connections.size(); ++i) {
		const Connection &conn = connections[i];	
		if (fairleads && conn.where == "vessel") {
			sfair << "\"";
			sfair << Format("Point%d", i+1) + "FX, ";
			sfair << Format("Point%d", i+1) + "FY, ";
			sfair << Format("Point%d", i+1) + "FZ";
			sfair << "\"\n";
		} else if (anchors && (conn.where == "fixed" || conn.where == "anchor")) {
			sanch << "\"";
			sanch << Format("Point%d", i+1) + "FX, ";
			sanch << Format("Point%d", i+1) + "FY, ";
			sanch << Format("Point%d", i+1) + "FZ";
			sanch << "\"\n";
		}
	}
	for (int i = 0; i < lineProperties.size(); ++i) {		
		const LineProperty &line = lineProperties[i];
		int from = -1, to = -1;
		double zfrom, zto;
		for (int i = 0; i < connections.size(); ++i) {
			const Connection &conn = connections[i];
			if (line.from == conn.name) {
				from = i;
				zfrom = conn.z;
			} else if (line.to == conn.name) {
				to = i;
				zto = conn.z;
			}
		}
		if (from < 0) {
			BEM::PrintError(Format(t_("Point from %s of line %s not found"), line.from, line.name));
			return false;
		}
		if (to < 0) {
			BEM::PrintError(Format(t_("Point to %s of line %s not found"), line.to, line.name));
			return false;
		}
		if (zfrom > zto)
			Swap(from, to);		// MoorDyn expects first anchor and second fairlead
		out << Format(" %3d %12=s %8d %8d %12.3f %8d %8s\n", i+1, line.nameType, from+1, to+1, line.length, line.numseg, mooring == 2 ? "p" : "-");
	}
	
	double _dtM 		= IsNull(dtM) 		? 0.001 : dtM;
	double _kbot 		= IsNull(kbot) 		? 3E6 : kbot;
	double _cbot 		= IsNull(cbot) 		? 3E5 : cbot;
	double _dtIC 		= IsNull(dtIC) 		? 1 : dtIC;
	double _TmaxIC 		= IsNull(TmaxIC) 	? 60 : TmaxIC;
	double _CdScaleIC 	= IsNull(CdScaleIC) ? 4 : CdScaleIC;
	double _threshIC 	= IsNull(threshIC) 	? 0.001 : threshIC;
	
	out <<	"---------------------- SOLVER OPTIONS ---------------------------------------\n";
	if (!IsNull(_dtM))
		out << Format("%-22.4f %11<s- time step to use in mooring integration (s). In case of NaNs reduce < 0.001, although slows down speed\n", _dtM, "dtM");
	if (!IsNull(_kbot))
		out << Format("%-22.1e %11<s- bottom stiffness (Pa/m)\n", _kbot, "kbot");
	if (!IsNull(_cbot))
		out << Format("%-22.1e %11<s- bottom damping (Pa-s/m)\n", _cbot, "cbot");
	if (!IsNull(_dtIC))
		out << Format("%-22.1f %11<s- time interval for analyzing convergence during IC gen (s)\n", _dtIC, "dtIC");
	if (!IsNull(_TmaxIC))
		out << Format("%-22.1f %11<s- max time for ic gen (s)\n", _TmaxIC, "TmaxIC");
	if (!IsNull(_CdScaleIC))
		out << Format("%-22.1f %11<s- factor by which to scale drag coefficients during dynamic relaxation (-)\n", _CdScaleIC, "CdScaleIC");
	if (!IsNull(_threshIC))
		out << Format("%-22.3f %11<s- threshold for IC convergence (-)\n", _threshIC, "threshIC");
	
	out << Format("%-22.1f %11<s- water depth (m)\n", depth, "WtrDpth");
	out << "0                      OutSwitch  - enable .MD output (0/1)\n";
	
	if (fairleads || anchors || mooring == 1) {
		out <<  "------------------------ OUTPUTS --------------------------------------------\n";
		out << sfair << sanch;

		if (mooring == 1) {
			for (int i = 0; i < lineProperties.size(); ++i) {
				const LineProperty &line = lineProperties[i];
				int numseg = max(int(line.numseg*ratioSegsToPlot), min(line.numseg, minSegsToPlot));
				UVector<int> points(numseg+1);
				double step = 1;
				if (line.numseg > numseg)
					step = line.numseg/numseg;
				for (int i = 0; i < numseg; ++i) 
		        	points[i] = (int)(i * step);
				points[numseg] = line.numseg;
				
				out << "\"";
				for (int ip = 0; ip < points.size(); ip++) {
					if (ip > 0)
						out << ", ";
					String str = "L" + FormatInt(i+1) + "N" + FormatInt(points[ip]) + "P";
					out << str << "X, " << str << "Y, " << str << "Z";
				}
				out << "\"\n";	
			}
		}
		out << "END\n";
	}
	out <<  "------------------------- need this line ------------------------------------";
	
	return true;
}

String Mooring::Test() {
	for (const auto &line : lineProperties)
		if (line.from != "" && line.from == line.to)
			return Format(t_("Line '%s' ends '%s' have to be different"), line.name, line.from);
		
	return String();
}

const Mooring::Connection *Mooring::GetConnectionP(String name) const {
	for (const Connection &con : connections) {
		if (con.name == name)
			return &con;
	}
	return nullptr;
}

const Mooring::Connection &Mooring::GetConnection(String name) const {
	const Mooring::Connection *p = GetConnectionP(name);	
	if (p) 
		return *p;
	throw Exc(t_(Format("Connection '%s' is not found", name)));
}

Mooring::LineType &Mooring::GetLineType(String name) {
	for (LineType &line : lineTypes) {
		if (line.name == name)
			return line;
	}
	throw Exc(t_(Format("Line type '%s' is not found", name)));
}

bool Mooring::FindClosest(Mooring::ClosestInfo &info) {
	auto Distance3D = [](double x1, double y1, double z1, double x2, double y2, double z2)->double {
	    return ::sqrt(sqr(x2 - x1) + sqr(y2 - y1) + sqr(z2 - z1));
	};
	
	auto DistanceClosestConnector = [&](double x, double y, double z)->double {
		double mind = std::numeric_limits<double>::max();
		for (Connection &c : connections) 
			mind = min(mind, Distance3D(x, y, z, c.x, c.y, c.z));
		return mind;
	};
	
	if (lineProperties.IsEmpty()) {
		info.distance = Null;
		return false;
	}
	
	info.distance = std::numeric_limits<double>::max();

    for(int i = 0; i < lineProperties.size(); i++) {
        const LineProperty& lineA = lineProperties[i];

        for(int j = i + 1; j < lineProperties.size(); j++) {
            const LineProperty& lineB = lineProperties[j];

            for(int p1 = 0; p1 < lineA.x.size(); p1++) {
                if (DistanceClosestConnector(lineA.x[p1], lineA.y[p1], lineA.z[p1]) > 10) {				// Discards closest connector
	                for(int p2 = 0; p2 < lineB.x.size(); p2++) {
	                    if (DistanceClosestConnector(lineB.x[p2], lineB.y[p2], lineB.z[p2]) > 10) {		// Discards closest connector
		                    double dist = Distance3D(lineA.x[p1], lineA.y[p1], lineA.z[p1],
		                                             lineB.x[p2], lineB.y[p2], lineB.z[p2]);
		                    if(dist < info.distance) {
		                        info.distance = dist;
		                        info.line1 = i;
		                        info.point1 = p1;
		                        info.line2 = j;
		                        info.point2 = p2;
		                    }
	                    }
	                }
                }
            }
        }	
    }
    if (info.distance == std::numeric_limits<double>::max()) {
        info.distance = Null;
        return false;
    }
    return true;
}

bool Mooring::Calc(double rho_water, int num, double rho_m3) {
	if (IsNull(num) || num <= 0)
		return false;
	
	for (LineProperty &line : lineProperties) {
		UVector<double> vpos;
		
		const LineType &linetype = GetLineType(line.nameType);
		const Connection &from = GetConnection(line.from);
		const Connection &to = GetConnection(line.to);
		double zanchor = Null, zvessel = Null;
		double fromx = from.x;
		double tox = to.x;
		double fromy = from.y;
		double toy = to.y;
		
		bool reverse = false;
		int id;
		id = FindVessel(from.where);
		if (id >= 0) {
			fromx += vessels[id].dx;  
			fromy += vessels[id].dy;
			zvessel = from.z;
			zanchor = to.z;
		}
		id = FindVessel(to.where);
		if (id >= 0) {
			tox += vessels[id].dx;
			toy += vessels[id].dy;
			zvessel = to.z;
			zanchor = from.z;
			reverse = true;
		}
		if (IsNull(zanchor)) { // Both are fixed
			if (to.z < from.z) {
				zvessel = to.z;
				zanchor = from.z;
				reverse = true;
			} else {
				zvessel = from.z;
				zanchor = to.z;
			}
		}
		
		double xanchorvessel = sqrt(sqr(fromx - tox) + sqr(fromy - toy));
		
		line.status = Catenary(linetype.mass, rho_m3, rho_water, line.length, linetype.bl, 
					xanchorvessel, zanchor + depth, zvessel + depth, line.fanchorvessel, line.fVanchor, line.fVvessel, 
					line.lenonfloor, vpos, line.z, num);	
		
		if (reverse) 
			::Reverse(vpos);

		line.fVanchor = -line.fVanchor;
		line.theta = atan2(toy - fromy, tox - fromx);
		line.x.SetCount(vpos.size());		
		line.y.SetCount(vpos.size());		
		for (int i = 0; i < vpos.size(); ++i) {
			line.x[i] = fromx + vpos[vpos.size()-i-1]*cos(line.theta);
			line.y[i] = fromy + vpos[vpos.size()-i-1]*sin(line.theta);
		}
		for (double &z : line.z)
			z -= depth;
		if (reverse)
			line.theta += M_PI;
	}
	return true;
}
	
	
	

