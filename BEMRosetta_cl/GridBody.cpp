// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2024, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
//#include <ScatterDraw/DataSource.h>
//#include <ScatterDraw/Equation.h>
//#include "functions.h"
//#include <SysInfo/Crash.h>
//#include <STEM4U/SeaWaves.h>
//#include <MatIO/matio.h>

using namespace Upp;
using namespace Eigen;


void GridBody::TableHeaders(const Hydro &hy, double head, bool pot, UVector<String> &str, UVector<char> &group, UVector<int> &col) {
	const char *xyz[] = {"x", "y", "z"};
	
	str.Clear();			group.Clear();		col.Clear();
	
	str << t_("#");			group << 'n';		col << 0;
	
	if (hy.IsLoadedMesh()) {
		int n = 0;
		for (int c = 0; c < 4; ++c) {
			str << Format(t_("Node %d"), c+1);			group << 'm';	col << n++;
		}		
		str << t_("Area");								group << 'm';	col << n++;
		for (int c = 0; c < 3; ++c) {
			str << Format(t_("Center %s"), xyz[c]);		group << 'm';	col << n++;
		}
		for (int c = 0; c < 3; ++c) {
			str << Format(t_("Normal %s"), xyz[c]);		group << 'm';	col << n++;
		}
	}
	if (hy.IsLoadedPotsRad()) {
		int n = 0;
		for (int i = 0; i < 6; ++i) {
			if (pot) {
				str << Format(t_("|Φrad| %s"), 		BEM::StrDOF(i));	group << 'r'; col << n++;
				str << Format(t_("arg(Φrad) %s"), 	BEM::StrDOF(i));	group << 'r'; col << n++;
			} else {
				str << Format(t_("|Prad| %s"), 		BEM::StrDOF(i));	group << 'r'; col << n++;
				str << Format(t_("arg(Prad) %s"), 	BEM::StrDOF(i));	group << 'r'; col << n++;	
			}
		}
		if (Bem().onlyDiagonal) {
			for (int r = 0; r < 6; ++r) {
				str << Format(t_("A_%s"), BEM::StrDOF(r));		group << 'r';		col << n++;
			}
		} else {
			for (int r = 0; r < 6; ++r) 
				for (int c = 0; c < 6; ++c) {
					str << Format(t_("A_%s_%s"),BEM::StrDOF(r), BEM::StrDOF(c));	group << 'r';	col << n++;
				}
		}
		if (Bem().onlyDiagonal) {
			for (int r = 0; r < 6; ++r) {
				str << Format(t_("B_%s"), BEM::StrDOF(r)); 		group << 'r'; col << n++;
			}
		} else {
			for (int r = 0; r < 6; ++r) 
				for (int c = 0; c < 6; ++c) {
					str << Format(t_("B_%s_%s"),BEM::StrDOF(r), BEM::StrDOF(c)); 	group << 'r'; 	col << n++;
				}
		}
	} 
	if (hy.IsLoadedPotsInc()) {
		int n = 0;
		if (pot) {
			str << Format(t_("|Φinc| %.0f"), head); 	group << 'i'; 	col << n++;
			str << Format(t_("arg(Φinc) %.0f"), head);	group << 'i'; 	col << n++;
		} else {
			str << Format(t_("|Pinc| %.0f"), head); 	group << 'i'; 	col << n++;
			str << Format(t_("arg(Pinc) %.0f"), head);	group << 'i'; 	col << n++;	
		}
		for (int idf = 0; idf < 6; ++idf) {
			str << Format(t_("|Ffk| %s %.0f"), BEM::StrDOF(idf), head); 	group << 'i'; col << n++;
			str << Format(t_("arg(Ffk) %s %.0f"), BEM::StrDOF(idf), head); 	group << 'i'; col << n++;
		}
	}
	if (hy.IsLoadedPotsIncB()) {
		int n = 0;
		if (pot) {
			str << Format(t_("|Φinc_bmr| %.0f"), head), 	group << 'b'; 	col << n++;
			str << Format(t_("arg(Φinc_bmr) %.0f"), head),	group << 'b'; 	col << n++;
		} else {
			str << Format(t_("|Pinc_bmr| %.0f"), head), 	group << 'b'; 	col << n++;
			str << Format(t_("arg(Pinc_bmr) %.0f"), head),	group << 'b'; 	col << n++;			
		}
		for (int idf = 0; idf < 6; ++idf) {
			str << Format(t_("|Ffk_bmr| %s %.0f"), BEM::StrDOF(idf), head); 	group << 'b'; col << n++;
			str << Format(t_("arg(Ffk_bmr) %s %.0f"), BEM::StrDOF(idf), head); 	group << 'b'; col << n++;
		}
	}
	if (hy.IsLoadedPotsDif()) {
		int n = 0;
		if (pot) {
			str << Format(t_("|Φdif|"), head); 	group << 'd'; col << n++;
			str << Format(t_("arg(Φdif) %.0f"), head); 	group << 'd'; col << n++;
		} else {
			str << Format(t_("|Pdif|"), head); 	group << 'd'; col << n++;
			str << Format(t_("arg(Pdif) %.0f"), head); 	group << 'd'; col << n++;	
		}
	}
}
		

void GridBody::Load(int idx, int ib, int &numNodes, int &numPanels) {
	const char *xyz[] = {"x", "y", "z"};
	
	const Hydro &hy = Bem().hydros[idx];
	{
		numNodes = hy.dt.msh[ib].dt.mesh.nodes.size();
		
		grdNodes.Clear();
		grdNodes.SetVirtualCount(numNodes);
		
		dataSourceNodes.Clear();
		if (hy.IsLoadedMesh()) {
			grdNodes.AddVirtualCol(t_("#"), dataSourceNodes.Add().Init(hy.dt.msh[ib].dt.mesh, -2), 60);
			for (int c = 0; c < 3; ++c) 
				grdNodes.AddVirtualCol(Format("%s", xyz[c]), dataSourceNodes.Add().Init(hy.dt.msh[ib].dt.mesh, -2), 80);
		}
	}{
		UVector<int> nnum;
		UVector<String> snum;
		if (hy.IsLoadedMesh()) {
			nnum << hy.dt.msh[ib].dt.mesh.panels.size();
			snum << "mesh";
		}
		if (hy.IsLoadedPotsRad(ib)) {
			nnum << hy.dt.pots_rad[ib].size();
			snum << "radiation";
		}
		if (hy.IsLoadedPotsDif(ib)) {
			nnum << hy.dt.pots_dif[ib].size();
			snum << "diffraction";
		}
		if (hy.IsLoadedPotsInc(ib)) {
			nnum << hy.dt.pots_inc[ib].size();
			snum << "incident";
		}
		if (hy.IsLoadedPotsIncB(ib)) {
			nnum << hy.dt.pots_inc_bmr[ib].size();
			snum << "incident_bemr";
		}
		for (int i = 0; i < nnum.size()-1; ++i)
			if (nnum[i] != nnum[i+1]) 
				throw Exc(Format("Number of panels doesn't match between %s and %s", snum[i], snum[i+1]));
			
		numPanels = First(nnum);
			
		grdPanels.Clear();
		grdPanels.SetVirtualCount(numPanels);
		
		dataSourcePanels.Clear();
		
		double head = 0;
		
		UVector<String> str;
		UVector<char> group;
		UVector<int> col;
		TableHeaders(hy, head, true, str, group, col);
		
		for (int i = 0; i < str.size(); ++i)
			grdPanels.AddVirtualCol(str[i], dataSourcePanels.Add().Init(idx, ib, group[i], col[i]), 60);		
	}
}


void GridBody::UpdatePanelHeaders() {
	if (dataSourcePanels.IsEmpty())
		return;
	
	int idx = First(dataSourcePanels).idx;
	int ih = First(dataSourcePanels).ih;
	bool pot = First(dataSourcePanels).pot;
	
	const Hydro &hy = Bem().hydros[idx];
	
	double head = hy.dt.head[ih];
	
	UVector<String> str;
	UVector<char> group;
	UVector<int> col;
	TableHeaders(hy, head, pot, str, group, col);
	
	for (int i = 0; i < str.size(); ++i)
		grdPanels.SetVirtualHeader(i, str[i]);	
}

Value GridBody::DataSourceNodes::Format(const Value& q) const {
	ASSERT(pmesh);
	int iq = q;
	if (pmesh->nodes.size() <= iq)
		return Null;
	
	const Point3D &p = pmesh->nodes[iq];
	switch (xyz) {
	case -2:	return iq + 1;		// id
	case  0:	return p.x;			// x
	case  1:	return p.y;			// y
	case  2:	return p.z;			// z
	default: NEVER();return Null;
	}
}

Value GridBody::DataSourcePanels::Format(const Value& q) const {
	ASSERT(idx >= 0);
	const Hydro &hy = Bem().hydros[idx];
	int ip = q;
	const Surface &s = hy.dt.msh[ib].dt.mesh;
	if (s.panels.size() <= ip)
		return Null;
	
	if (group == 'n')
		return ip + 1;		// id
	else if (group == 'm') {
		switch (col) {
		case 0:	return s.panels[ip].id[0];
		case 1:	return s.panels[ip].id[1];
		case 2:	return s.panels[ip].id[2];
		case 3:	return s.panels[ip].id[3];
		case 4:	return s.panels[ip].surface0 + s.panels[ip].surface1;
		case 5:	return s.panels[ip].centroidPaint.x;
		case 6:	return s.panels[ip].centroidPaint.y;
		case 7:	return s.panels[ip].centroidPaint.z;
		case 8:	return s.panels[ip].normalPaint.x;
		case 9:	return s.panels[ip].normalPaint.y;
		case 10:return s.panels[ip].normalPaint.z;
		}
	} else if (group == 'r') {
		if (col < 2*6) {
			int idf = col/2;
			if (pot) {
				if (col%2 == 0)
					return abs		(hy.dt.pots_rad[ib][ip][idf][ifr]);
				else
					return ToDeg(arg(hy.dt.pots_rad[ib][ip][idf][ifr]));
			} else {
				if (col%2 == 0)
					return abs		(hy.P_rad(ib, ip, idf, ifr));
				else
					return ToDeg(arg(hy.P_rad(ib, ip, idf, ifr)));
			}
		}
		int cl = col - 12;
		if (Bem().onlyDiagonal) {
			if (cl < 6) {
				int col2 = cl;	
				return hy.dt.Apan(ib, ip, col2, col2, ifr);
			} else {
				int col2 = cl - 6;	
				return hy.B_pan(ib, ip, col2, col2, ifr);
			}
		} else {
			if (cl < 36) {
				int row = cl/6;
				int col2 = cl - 6*row;
				return hy.dt.Apan(ib, ip, row, col2, ifr);
			} else {
				cl -= 36;	
				int row = cl/6;
				int col2 = cl - 6*row;
				return hy.B_pan(ib, ip, row, col2, ifr);
			}
		}
	} else if (group == 'i') {
		if (col < 2) {
			if (pot) {
				if (col == 0)
					return abs		(hy.dt.pots_inc[ib][ip][ih][ifr]);
				else
					return ToDeg(arg(hy.dt.pots_inc[ib][ip][ih][ifr]));
			} else {
				if (col == 0)
					return abs		(hy.P_inc(ib, ip, ih, ifr));
				else
					return ToDeg(arg(hy.P_inc(ib, ip, ih, ifr)));
			}
		}
		int cl = col - 2;
		int idf = cl/2;
		if (cl%2 == 0)
			return 		 abs(hy.Ffk_pan(ib, ip, ih, idf, ifr));
		else
			return ToDeg(arg(hy.Ffk_pan(ib, ip, ih, idf, ifr)));
	} else if (group == 'b') {
		if (col < 2) {
			if (pot) {
				if (col == 0)
					return abs		(hy.dt.pots_inc_bmr[ib][ip][ih][ifr]);
				else
					return ToDeg(arg(hy.dt.pots_inc_bmr[ib][ip][ih][ifr]));
			} else {
				if (col == 0)
					return abs		(hy.P_inc_bmr(ib, ip, ih, ifr));
				else
					return ToDeg(arg(hy.P_inc_bmr(ib, ip, ih, ifr)));	
			}
		}
		int cl = col - 2;
		int idf = cl/2;
		if (cl%2 == 0)
			return 		 abs(hy.Ffk_pan_bmr(ib, ip, ih, idf, ifr));
		else
			return ToDeg(arg(hy.Ffk_pan_bmr(ib, ip, ih, idf, ifr)));
	} else if (group == 'd') {
		if (col < 2) {
			if (pot) {
				if (col == 0)
					return abs		(hy.dt.pots_dif[ib][ip][ih][ifr]);
				else
					return ToDeg(arg(hy.dt.pots_dif[ib][ip][ih][ifr]));
			} else {
				if (col == 0)
					return abs		(hy.P_dif(ib, ip, ih, ifr));
				else
					return ToDeg(arg(hy.P_dif(ib, ip, ih, ifr)));	
			}
		}
	}
	return Null;
}

