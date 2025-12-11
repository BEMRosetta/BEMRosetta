// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2025, the BEMRosetta author and contributors
#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <SurfaceCanvas/SurfaceCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <TabBar/TabBar.h>
#include <DropGrid/DropGrid.h>

#include <BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;

#include "main.h"


MenuProcessInertia::MenuProcessInertia() {
	CtrlLayout(*this);

	opInertia = 0;
	opMass = 0;
	
	mass.WhenAction    = [&] {
		if (Bem().surfs[_idx].dt.mesh.volume > 0)
			density <<= mass/Bem().surfs[_idx].dt.mesh.volume;
		Action();
	};
	density.WhenAction = [&] {
		mass <<= density*Bem().surfs[_idx].dt.mesh.volume; 
		Action();
	};
	
	opMass.WhenAction    = [&]{OpMass_WhenAction(true);};
	opInertia <<= THISBACK(Action);
	
	edDivideZ <<= THISBACK(Action);
	opDivide  <<= THISBACK(Action);
	
	opInertia.MinCaseHeight(int(1.5*StdFont().GetHeight())).SetVertical();
	opMass.MinCaseHeight(int(1.5*StdFont().GetHeight())).SetVertical();
	
	ok << [&] {Close();};
	
	butSetC0Vol << [&] {
		x_0 = ~x_g;
		y_0 = ~y_g;
		z_0 = ~z_g;
		Action();
	};
	
	butCopy  <<= THISBACK(CopyToBody);
	
	x_0 <<= THISBACK(Action);
	y_0 <<= THISBACK(Action);
	z_0 <<= THISBACK(Action);
	x_g <<= THISBACK(Action);
	y_g <<= THISBACK(Action);
	z_g <<= THISBACK(Action);
	
	grid.MultiSelect().Removing(false).Clipboard().Sorting(false).Absolute().Editing().SelectRow(false);
	for (int i = 0; i < 6; ++i) {
		edit[i].NotNull();
		grid.AddColumn(BEM::StrDOF(i), 50).Edit(edit[i]);
	}
}

void MenuProcessInertia::CopyToBody() {
	Body &msh = Bem().surfs[_idx];
	
	msh.dt.c0.x = _mb->menuProcess.x_0 = ~x_0;
	msh.dt.c0.y = _mb->menuProcess.y_0 = ~y_0;
	msh.dt.c0.z = _mb->menuProcess.z_0 = ~z_0;
	
	msh.dt.cg.x = _mb->menuProcess.x_g = ~x_g;
	msh.dt.cg.y = _mb->menuProcess.y_g = ~y_g;
	msh.dt.cg.z = _mb->menuProcess.z_g = ~z_g;
	
	//msh.dt.cg0 = msh.dt.cg;
	
	if (opMass < 3) {
		msh.dt.M.resize(6, 6);
		for (int r = 0; r < 6; ++r)		
			for (int c = 0; c < 6; ++c)
				msh.dt.M(r, c) = ScanDouble(grid.Get(r, c).ToString());
		_mb->menuProcess.mass <<= msh.GetMass();
	}
	
	msh.AfterLoad(Bem().rho, Bem().g, false, false);
	_mb->UpdateLast(_idx);
	_mb->mainView.FullRefresh(*_mb);
	
	Close();
}

void MenuProcessInertia::Init(MainBody &b, int idx) {
	_mb = &b;
	_idx = idx;
	
	opInertia = 0;
	opMass.DisableCase(3);
	opMass.DisableCase(4);
	opMass = 0;
	
	Body &mesh = Bem().surfs[idx];
	edVolume <<= mesh.dt.mesh.volume;
	
	x_0 = mesh.dt.c0.x;
	y_0 = mesh.dt.c0.y;
	z_0 = mesh.dt.c0.z;
	
	x_g = mesh.dt.cg.x;
	y_g = mesh.dt.cg.y;
	z_g = mesh.dt.cg.z;
	
	mass <<= mesh.GetMass();
	if (!IsNull(mesh.dt.mesh.volume) && mesh.dt.mesh.volume > 0)
		density <<= mesh.GetMass()/mesh.dt.mesh.volume;
	
	grid.Clear();
	if (mesh.dt.M.size() == 36) {
		for (int r = 0; r < 6; ++r)		
			for (int c = 0; c < 6; ++c)
				grid.Set(r, c, mesh.dt.M(r, c));
	}
	Action();
}

void MenuProcessInertia::Action() {
	int opmass = ~opMass;
	
	grid.Ready(false);
	if (opmass <= 3) {
		for (int c = 0; c < 6; ++c) 
			grid.GetColumn(c).Width(50);
	} else {
		for (int c = 0; c < 3; ++c)
			grid.GetColumn(c).Hidden();
	}
	grid.Ready(true);
		
	auto EmptyGrid = [&] {
		grid.Clear();
		if (opmass <= 3) {
			for (int r = 0; r < 6; ++r)		
				for (int c = 0; c < 6; ++c)
					grid.Set(r, c, "-");
		} else {
			for (int r = 0; r < 3; ++r)	
				for (int c = 0; c < 3; ++c) 
					grid.Set(r, c+3, "-");
		}	
	};
	try {
		Point3D c0(~x_0, ~y_0, ~z_0);
		if (IsNull(c0)) {
			EmptyGrid();
			return;	
		}
		Point3D cg(~x_g, ~y_g, ~z_g);
		
		Body &mesh = Bem().surfs[_idx];
		
		grid.Editing(opInertia == 0);
		x_g.SetEditable(~opInertia == 0);
		y_g.SetEditable(~opInertia == 0);
		z_g.SetEditable(~opInertia == 0 || ~opDivide == true);

		for (int i = 0; i < 6; ++i) 
			edit[i].SetEditable(opInertia == 0);
		
		opDivide.Show(opInertia > 0);
		edDivideZ.Show(opInertia > 0);
		edDivideZ.Enable(~opDivide == true);
				
		if (opInertia == 0) {
			opMass.DisableCase(3);
			opMass.DisableCase(4);
		} else {
			opMass.EnableCase(3);
			opMass.EnableCase(4);
		}
		if (opMass == 3 || opMass == 4) 
			opInertia.DisableCase(0);
		else
			opInertia.EnableCase(0);
	
		bool isvol;
		if (opInertia == 0) {
			EmptyGrid();
			return;	
		}
		
		if (opInertia == 1) {
			isvol = true;
			Point3D ccg = mesh.dt.mesh.GetCentreOfBuoyancy();
			cg.x = ccg.x;	cg.y = ccg.y;
			if (IsNull(cg.z)) {
				cg.z = ccg.z;
				z_g <<= cg.z;
			}
		} else {
			isvol = false;
			Point3D ccg = mesh.dt.mesh.GetCentreOfGravity_Surface();
			cg.x = ccg.x;	cg.y = ccg.y;
			if (IsNull(cg.z)) {
				cg.z = ccg.z;
				z_g <<= cg.z;
			}
		}
		x_g <<= cg.x;
		y_g <<= cg.y;

		
		if (isvol && mesh.dt.mesh.VolumeMatch(Bem().volError, Bem().volError) < 0) {
			opInertia = 0;
			throw Exc(t_("Incomplete mesh or wrongly oriented panels"));
		}
		
		if (IsNull(cg)) {
			EmptyGrid();
			return;	
		}
		
		Matrix3d inertia3;
		MatrixXd inertia6;
		if (!opDivide) {
			if (!mesh.dt.mesh.GetInertia33(inertia3, c0, isvol, true)) {
				EmptyGrid();
				return;	
			}
		} else {	// Divide body in 2
			if (IsNull(edDivideZ)) {
				EmptyGrid();
				return;
			}
			
			Surface msh = clone(mesh.dt.mesh);
			msh.Translate(0, 0, -edDivideZ);
			
			Surface msha, mshb;
			msha.CutZ(msh,  1).GetPanelParams().GetVolume();
			mshb.CutZ(msh, -1).GetPanelParams().GetVolume();
			
			Point3D ca = msha.GetCentreOfBuoyancy();
			Point3D cb = mshb.GetCentreOfBuoyancy();
			
			if (IsNull(ca) || IsNull(cb)) {
				EmptyGrid();
				return;	
			}
			
			Surface mshlid = clone(msha);
			mshlid.Flatten('z');
			
			mshb.Append(mshlid);
			mshlid.Orient();
			msha.Append(mshlid);
			
			Point3D c00 = c0;
			c00.z -= edDivideZ;
			
			Matrix3d ia, ib;
			if (!msha.GetInertia33(ia, c00, isvol, true) || !mshb.GetInertia33(ib, c00, isvol, true)) {
				EmptyGrid();
				return;	
			}
			
			ca.z += edDivideZ;
			cb.z += edDivideZ;
			
			double m = 1;
			double ma = m*(cg.z - cb.z)/(ca.z - cb.z);
			double mb = m - ma;
			
			inertia3 = ia*ma + ib*mb;
		}
		Surface::GetInertia66(inertia6, inertia3, cg, c0, true);
		
		if (opmass == 3)
			inertia6 *= mesh.dt.mesh.volume;
		else if (opmass < 3) {
			double m;
			if (!IsNull(mass) && (opmass == 0 || opmass == 1))
				m = mass;
			else {
				if (!IsNull(density))
					m = density*mesh.dt.mesh.volume; 
			}
			inertia6 *= m;
		}
		grid.Clear();
		if (opmass <= 3) {
			for (int r = 0; r < 6; ++r)		
				for (int c = 0; c < 6; ++c)
					grid.Set(r, c, inertia6(r, c));
		} else {
			for (int r = 0; r < 3; ++r)	{	
				for (int c = 0; c < 3; ++c) {
					double val = inertia3(r, c);
					int sign = Sign(val);
					grid.Set(r, c+3, sign*sqrt(abs(val)));
				}
			}
		}
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}	
}

void MenuProcessInertia::OpMass_WhenAction(bool action) {
	Body &mesh = Bem().surfs[_idx];
	
	mass.Enable(opMass == 1);
	density.Enable(opMass == 2);
	
	if (opMass == 0) {
		mass = mesh.GetMass();
		if (!IsNull(mass) && mesh.dt.mesh.volume > 0) 
			density = mass/mesh.dt.mesh.volume;
	} else if (opMass == 1) {
		if (IsNull(mass)) {
			if (IsNull(density))
				mass = mesh.GetMass(); 
			else
				mass = density*mesh.dt.mesh.volume;
		}
	} else if (opMass == 2) {
		if (IsNull(density)) {
			if (IsNull(mass))
				density = mesh.GetMass()/mesh.dt.mesh.volume; 
			else if (mesh.dt.mesh.volume > 0)
				density = mass/mesh.dt.mesh.volume; 
		}
	} 
	
	String str;
	if (opMass == 3)
		str = t_("Volume moments of inertia");
	else if (opMass == 4) 
		str = t_("Radii of gyration");
	else				
		str = t_("Moments of inertia");
	
	butCopy.Enable(opMass <= 2);
	
	labInertia.SetLabel(str);
		
	Action();
}	