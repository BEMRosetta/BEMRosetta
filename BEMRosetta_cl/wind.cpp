// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2024, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"

#include <iostream>
#include <fstream>
#include <unsupported/Eigen/CXX11/Tensor>

using namespace Eigen;


String Wind::Load(String fileName, String ext) {
	const UVector<String> extStr = {".bts", ".wnd"};
	UVector<String> extS;
	
	fileName = Trim(fileName);
	if (*fileName.Last() == '.')
		fileName = fileName.Left(fileName.GetCount()-1);		// Removes a last point, if it exists
	
	String fileExt = GetFileExt(fileName);
	bool hasExt = !IsEmpty(fileExt) && fileExt.GetCount() <= 4;	// An extension larger than ... is not an extension
	
	if (IsEmpty(ext)) {
		ext = fileExt;
		if (IsEmpty(ext))	
			extS = clone(extStr);
		else
			extS << ext;
	} else
		extS << ext;

	String ret = Format(t_("File '%s' not found"), fileName);
	for (const String &ext : extS) {
		String file;
		if (hasExt)		
 			file = ForceExt(fileName, ext);
		else
			file = fileName + ext;
		if (FileExists(file)) {
			if (ext == ".bts") {
				if(IsEmpty(ret = static_cast<BTSWind&>(*this).LoadBTS(file)))
					break;	
			} else if (ext == ".wnd") {
				if(IsEmpty(ret = static_cast<WNDWind&>(*this).LoadWND(file)))
					break;	
			} else
				return Format(t_("'%s' format not supported for reading"), fileName);
		}
	}
	if (IsEmpty(ret))
		SetYZ();
	return ret;
}

void Wind::SetYZ() {
	Arange(yPos, 0, ny-1);
	yPos = yPos.array()*dy - dy*(ny-1)/2;
	
	Arange(zPos, 0, nz-1);
	zPos = zPos.array()*dz + zGrid;
}

String Wind::Save(String fileName, String ext) {
	const UVector<String> extStr = {".bts"};
	UVector<String> extS;
	
	fileName = Trim(fileName);
	if (*fileName.Last() == '.')
		fileName = fileName.Left(fileName.GetCount()-1);	// Removes a last point, if it exists
	
	String fileExt = GetFileExt(fileName);
	bool hasExt = !IsEmpty(fileExt) && fileExt.GetCount() <= 4;		// An extension larger than ... is not an extension
	
	if (IsEmpty(ext)) {
		ext = fileExt;
		if (IsEmpty(ext))	
			extS = clone(extStr);
		else
			extS << ext;
	} else
		extS << ext;

	String ret = Format(t_("Impossible to save '%s'"), fileName);
	for (const String &ext : extS) {
		String file;
		if (hasExt)		
 			file = ForceExt(fileName, ext);
		else
			file = fileName + ext;
		
		if (ext == ".bts") {
			if(IsEmpty(ret = static_cast<BTSWind&>(*this).SaveBTS(file)))
				return "";	
		} else
			return Format(t_("'%s' format not supported for writing"), fileName);
	}
	return ret;
}

const char *Wind::GetWindTypeStr() const {
	switch(fc) {
	case 1:	return "1-component von Karman";
	case 2:	return "1-component Kaimal";
	case 3:	return "3-component von Karman";
	case 4:	return "Improved von Karman";
	case 5:	return "IEC-2 Kaimal";
	case 7:	return "General Kaimal";
	case 8:	return "Mann model";
	default:return "Unknown";
	}
}

void Wind::Report(Grid &grid) const {
	grid.HeaderRows(1);
	grid.HeaderCols({"Parameter", "Unit"}, {20, 6});
	grid.AddCol({GetFileName(fileName)}, 20);
	
	grid.Add({"File name", ""}, 		fileName);	
	grid.Add({"Description", ""}, 		description);
	grid.Add({"Duration", "s"},			Grid::Nvl(nt, Format("%.1f", nt*dt)));
	grid.Add({"Time step", "s"},		Grid::Nvl(dt, Format("%.6f", dt)));
	grid.Add({"Number of Z points", ""},Grid::Nvl(nz, FormatInt(nz)));
	grid.Add({"Number of Y points", ""},Grid::Nvl(ny, FormatInt(nz)));
	grid.Add({"Number of tower points", ""},Grid::Nvl(ntwr, FormatInt(nz)));
	grid.Add({"Grid height", "m"}, 		Grid::Nvl(nz, Format("%.1f", (nz-1)*dz)));
	grid.Add({"Grid width", "m"}, 		Grid::Nvl(ny, Format("%.1f", (ny-1)*dy)));
	grid.Add({"Hub height", "m"}, 	 	Grid::Nvl(zHub, Format("%.1f", zHub)));
	grid.Add({"Z bottom grid", "m"}, 	Grid::Nvl(zGrid, Format("%.2f", zGrid)));
	grid.Add({"Wind type", ""}, 	 	fc > 0 ? String(GetWindTypeStr()) : String());
}


void Wind::GetPos(double z, double y, int &idz, int &idy) {
	idz = FindClosest(zPos, z);
	idy = FindClosest(yPos, y);
}

VectorXd Wind::Get(int idz, int idy) {
	VectorXd ret(nt);
	
	for (int it = 0; it < nt; ++it)
		ret(it) = Norm(velocity(it, 0, idy, idz), velocity(it, 1, idy, idz), velocity(it, 2, idy, idz));
	
	return ret;
}

VectorXd Wind::GetTime() {
	VectorXd ret(nt);
	
	LinSpaced(ret, nt, 0, nt*dt);
	
	return ret;
}
	
void ArrayWind::Report(Grid &grid) {
	for (const Wind	&w : *this)
		w.Report(grid);
}
		