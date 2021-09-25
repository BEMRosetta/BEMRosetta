// SPDX-License-Identifier: GPL-3.0-or-later
#include "BEMRosetta.h"

#include <Core/Core.h>
#include <SysInfo/SysInfo.h>
#include <Surface/Surface.h>

using namespace Upp;


void Mooring::Jsonize(JsonIO &json) {
	json
		("lineTypes", lineTypes)
		("lineProperties", lineProperties)
		("connections", connections)
	;	
}

void Mooring::LineType::Jsonize(JsonIO &json) {
	if (json.IsLoading()) {
		mass = diameter = 0;
	}
	json
		("name", name)
		("mass", mass)
		("diameter", diameter)
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
			BEMData::PrintError(Format("Problem loading %s file", file));
			return false;
		}
	} catch (const Exc &e) {
		BEMData::PrintError(e);
		return false;
	}	
	return true;
}

bool Mooring::Save(String file) {
	try {
		if (!StoreAsJsonFile(*this, file, true)) {
			BEMData::PrintError(Format("Problem loading %s file", file));
			return false;
		}
	} catch (const Exc &e) {
		BEMData::PrintError(e);
		return false;
	}
	
	return true;
}

String Mooring::Test() {
	for (const auto &line : lineProperties)
		if (line.from == line.to)
			return Format(t_("Line '%s' ends '%s' have to be different"), line.name, line.from);
		
	return String();
}

