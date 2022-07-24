// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"

#include <Core/Core.h>
#include <SysInfo/SysInfo/SysInfo.h>
#include <Surface/Surface/Surface.h>
#include <STEM4U/STEM4U/Mooring.h>


using namespace Upp;


void Mooring::Jsonize(JsonIO &json) {
	json
		("lineTypes", lineTypes)
		("lineProperties", lineProperties)
		("connections", connections)
		("depth", depth)
	;	
}

void Mooring::LineType::Jsonize(JsonIO &json) {
	if (json.IsLoading()) {
		mass = diameter = bl = 0;
	}
	json
		("name", name)
		("mass", mass)
		("diameter", diameter)
		("BL", bl)
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
	;
}

void Mooring::Connection::Jsonize(JsonIO &json) {
	if (json.IsLoading()) {
		vessel = true;
		x = y = z = 0;
	}
	json
		("name", name)
		("vessel", vessel)
		("x", x)
		("y", y)
		("z", z)
	;
}

bool Mooring::Load(String file) {
	try {
		if (!LoadFromJsonFile(*this, ~file)) {
			BEM::PrintError(Format("Problem loading %s file", file));
			return false;
		}
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

String Mooring::Test() {
	for (const auto &line : lineProperties)
		if (line.from != "" && line.from == line.to)
			return Format(t_("Line '%s' ends '%s' have to be different"), line.name, line.from);
		
	return String();
}

Mooring::Connection &Mooring::GetConnection(String name) {
	for (Connection &con : connections) {
		if (con.name == name)
			return con;
	}
	throw Exc(t_(Format("Connection '%s' is not found", name)));
}

Mooring::LineType &Mooring::GetLineType(String name) {
	for (LineType &line : lineTypes) {
		if (line.name == name)
			return line;
	}
	throw Exc(t_(Format("Line type '%s' is not found", name)));
}

bool Mooring::Calc(double x, double y, double rho_water) {
	for (LineProperty &line : lineProperties) {
		UVector<double> vpos;
		int num = 50;
		double rho_m3 = 7850;	// Steel
		
		const LineType &linetype = GetLineType(line.nameType);
		const Connection &from = GetConnection(line.from);
		const Connection &to = GetConnection(line.to);
		double zanchor, zvessel;
		double fromx = from.x;
		double tox = to.x;
		double fromy = from.y;
		double toy = to.y;
		if (!from.vessel) {
			zanchor = to.z;
			zvessel = from.z;
			fromx += x;  
			fromy += y;
		} else {
			zanchor = from.z;
			zvessel = to.z;
			tox += x;  
			toy += y;
		}
		double xanchorvessel = sqrt(sqr(fromx - tox) + sqr(fromy - toy));
		
		line.status = Catenary(linetype.mass, rho_m3, rho_water, line.length, linetype.bl, 
					xanchorvessel, zanchor + depth, zvessel + depth, line.fanchorvessel, line.fVvessel, line.fVanchor, 
					line.lenonfloor, vpos, line.z, num);	
		
		line.fVanchor = -line.fVanchor;
		line.theta = atan2(toy - fromy, tox - fromx);
		line.x.SetCount(vpos.size());		
		line.y.SetCount(vpos.size());		
		for (int i = 0; i < vpos.size(); ++i) {
			line.x[i] = fromx + vpos[i]*cos(line.theta);
			line.y[i] = fromy + vpos[i]*sin(line.theta);
		}
	}
	return true;
}
	
	
	

