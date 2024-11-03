// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2023, the BEMRosetta author and contributors
#include "BEMRosetta.h"

#include <Core/Core.h>
#include <plugin/sqlite3/Sqlite3.h>

using namespace Upp;

#include "FastOut.h"

String FastOut::LoadDb(String file, Function <bool(String, int)> Status) {
	Clear();
	
	parameters << "TIME_SCALE";
	units << "seg";
	descriptions << "Time";
	
	Sqlite3Session sqlite3;

	if(!sqlite3.Open(file)) 
		return t_(Format("Impossible to open '%s'", file));

	fileName = file;
	
	Status(t_("Loading structure"), 10);
	{
		Sql sqlvariable(sqlite3);
		if (!sqlvariable.Execute("SELECT long_name, unit, description FROM 'variable'"))
			throw Exc(t_("'variable' table not found in .db file"));
	
		while(sqlvariable.Fetch()) {
			parameters << sqlvariable[0];
			units << sqlvariable[1];
			descriptions << sqlvariable[2];
		}
	}
	
	Status(t_("Loading wave series"), 20);
	int idwave = -1;
	UVector <double> timewave, wave;
	{
		Sql sqlwave(sqlite3);
		if (!sqlwave.Execute("SELECT * FROM 'evolution_wave_serie'"))
			throw Exc(t_("'evolution_wave_serie' table not found in .db file"));
		
		while(sqlwave.Fetch()) {
			timewave << sqlwave[0];
			wave << sqlwave[1];
		}
		String wavename = ToLower(sqlwave.GetColumnInfo(1).name);
		idwave = FindFunction(parameters, [&](const String &param)->bool {return ToLower(param) == wavename;});
	}
	
	int numCol = parameters.size();
	
	dataOut.SetCount(numCol+calcParams.size());
	
	Sql sqldata(sqlite3);
	if (!sqldata.Execute("SELECT * FROM 'time_serie'"))
		throw Exc(t_("'time_serie' table not found in .db file"));
	
	UVector<int> dataCols(sqldata.GetColumns());
	for (int i = 0; i < sqldata.GetColumns(); ++i) {
		String par = ToLower(sqldata.GetColumnInfo(i).name); 
		dataCols[i] = FindFunction(parameters, [&](const String &param)->bool {return ToLower(param) == par;});
	}

	Status(t_("Loading data"), 50);
	while(sqldata.Fetch()) {
		for (int i = 0; i < dataCols.size(); ++i) {
			if (dataCols[i] >= 0)
				dataOut[dataCols[i]] << sqldata[i];
		}
	}
	
	if (idwave >= 0) {
		dataOut[idwave].SetCount(GetNumData(), Null);
		for (int i = 0; i < GetNumData(); ++i) {
			int id = FindFunction(timewave, [&](const double &v)->bool {return abs(v - dataOut[0][i]) < 0.001;});
			if (id >= 0)
				dataOut[idwave][i] = wave[id];
		}
	}
	
	parameters[0] = "time";
	for (int i = parameters.size()-1; i >= 0; --i) {
		if (dataOut[i].IsEmpty()) {
			dataOut.Remove(i);
			parameters.Remove(i);
			units.Remove(i);
			descriptions.Remove(i);
		}
	}
	
	return "";
}