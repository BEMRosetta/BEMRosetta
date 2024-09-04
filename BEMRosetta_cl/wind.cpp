// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2024, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"

#include <iostream>
#include <fstream>
#include <unsupported/Eigen/CXX11/Tensor>

using namespace Eigen;


String Wind::Load(String _fileName, String ext) {
	const UVector<String> extStr = {".bts", ".wnd"};
	UVector<String> extS;
	
	fileName = Trim(_fileName);
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
	for (const String &eext : extS) {
		String file;
		if (hasExt)		
 			file = ForceExtSafer(fileName, eext);
		else
			file = fileName + eext;
		if (FileExists(file)) {
			if (eext == ".bts") {
				if(IsEmpty(ret = static_cast<BTSWind&>(*this).LoadBTS(file)))
					break;	
			} else if (eext == ".wnd") {
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

String Wind::Save(String fileSave, String ext) const {
	const UVector<String> extStr = {".bts"};
	UVector<String> extS;
	
	fileSave = Trim(fileSave);
	if (*fileSave.Last() == '.')
		fileSave = fileSave.Left(fileSave.GetCount()-1);	// Removes a last point, if it exists
	
	String fileExt = GetFileExt(fileSave);
	bool hasExt = !IsEmpty(fileExt) && fileExt.GetCount() <= 4;		// An extension larger than ... is not an extension
	
	if (IsEmpty(ext)) {
		ext = fileExt;
		if (IsEmpty(ext))	
			extS = clone(extStr);
		else
			extS << ext;
	} else
		extS << ext;

	String ret = Format(t_("Impossible to save '%s'"), fileSave);
	for (const String &eext : extS) {
		String file;
		if (hasExt)		
 			file = ForceExt(fileSave, eext);
		else
			file = fileSave + eext;
		
		if (eext == ".bts") {
			if(IsEmpty(ret = static_cast<const BTSWind&>(*this).SaveBTS(file)))
				return "";	
		} else
			return Format(t_("'%s' format not supported for writing"), fileSave);
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
	grid.SetCol(0).SetRow(0);
	grid.ColWidths({20, 6, 50});
	grid.SetRow({"Parameter", 			"Unit", "Data"});	
	grid.SetRow({"File name", 			"", 	GetFileName(fileName)});	
	grid.SetRow({"Description", 		"", 	description});
	grid.SetRow({"Duration", 			"s", 	Grid::Nvl(nt, Format("%.1f", nt*dt))});
	grid.SetRow({"Time step", 			"s", 	Grid::Nvl(dt, Format("%.6f", dt))});
	grid.SetRow({"Number of Z points", 	"", 	Grid::Nvl(nz, FormatInt(nz))});
	grid.SetRow({"Number of Y points", 	"", 	Grid::Nvl(ny, FormatInt(nz))});
	grid.SetRow({"Number of tower points", "", 	Grid::Nvl(ntwr, FormatInt(nz))});
	grid.SetRow({"Grid height", 		"m", 	Grid::Nvl(nz, Format("%.1f", (nz-1)*dz))});
	grid.SetRow({"Grid width", 			"m", 	Grid::Nvl(ny, Format("%.1f", (ny-1)*dy))});
	grid.SetRow({"Hub height", 			"m", 	Grid::Nvl(zHub, Format("%.1f", zHub))});
	grid.SetRow({"Z bottom grid", 		"m", 	Grid::Nvl(zGrid, Format("%.2f", zGrid))});
	grid.SetRow({"Wind type", 			"", 	fc > 0 ? String(GetWindTypeStr()) : String()});
}


void Wind::GetPos(double z, double y, int &idz, int &idy) {
	idz = FindClosest(zPos, z);
	idy = FindClosest(yPos, y);
}

VectorXd Wind::GetNorm(int idz, int idy) {
	VectorXd ret(nt);
	
	for (int it = 0; it < nt; ++it)
		ret(it) = Norm(velocity(it, 0, idy, idz), velocity(it, 1, idy, idz), velocity(it, 2, idy, idz));
	
	return ret;
}

VectorXd Wind::Get(int ic, int idz, int idy) {
	if (ic >= velocity.dimension(1))
		throw Exc(Format(t_("Component %d is not available"), ic+1));
	
	VectorXd ret(nt);
	
	for (int it = 0; it < nt; ++it)
		ret(it) = velocity(it, ic, idy, idz);
	
	return ret;
}

VectorXd Wind::GetTime() {
	VectorXd ret(nt);
	
	LinSpaced(ret, nt, 0, nt*dt);
	
	return ret;
}

int Wind::GetTimeId(double time) {
	return int(time/dt);
}

void Wind::SetTI(int uvw, float ti, float &tiOld, float offset) {
	if (uvw >= nffc)
		throw Exc(Format(t_("Component %d is higher than available components %d"), uvw, nffc));
	
	ASSERT(nt   == velocity.dimension(0));
	ASSERT(nffc == velocity.dimension(1));
	ASSERT(ny   == velocity.dimension(2));
	ASSERT(nz   == velocity.dimension(3));
	
	for (int it = 0; it < nt; ++it) {
	    for (int iz = 0; iz < nz; ++iz)
	        for (int iy = 0; iy < ny; ++iy) 
	            velocity(it,uvw,iy,iz) = (velocity(it,uvw,iy,iz) - offset)*(ti/tiOld) + offset;
	}
	tiOld = ti;       
}

void Wind::SetPowerLaw(float pl, float zh) {
	ASSERT(nt   == velocity.dimension(0));
	ASSERT(nffc == velocity.dimension(1));
	ASSERT(ny   == velocity.dimension(2));
	ASSERT(nz   == velocity.dimension(3));
	ASSERT(nz	== zPos.size());
	
	ASSERT(nffc > 0);	// At least 1 component, u
	int ic = 0;			// Only for u component
	
	for (int it = 0; it < nt; ++it) 
	    for (int iz = 0; iz < nz; ++iz)
	        for (int iy = 0; iy < ny; ++iy) 
	            velocity(it,ic,iy,iz) += mffws * (pow(zPos[iz]/zh, pl) - 1);
}

void Wind::SetFactor(float fu, float fv, float fw) {
	SetFactor(0, fu);
	SetFactor(1, fv);
	SetFactor(2, fw);
}

void Wind::SetFactor(int ic, float f) {
	if (ic >= nffc)
		return;
	
	ASSERT(nt   == velocity.dimension(0));
	ASSERT(nffc == velocity.dimension(1));
	ASSERT(ny   == velocity.dimension(2));
	ASSERT(nz   == velocity.dimension(3));
	ASSERT(nz	== zPos.size());	
	
	for (int it = 0; it < nt; ++it) 
	    for (int iz = 0; iz < nz; ++iz)
	        for (int iy = 0; iy < ny; ++iy) 
	            velocity(it,ic,iy,iz) *= f;
}

void ArrayWind::Report(Grid &grid) {
	for (const Wind	&w : *this)
		w.Report(grid);
}
		