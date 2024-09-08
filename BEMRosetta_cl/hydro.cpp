// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include <STEM4U/Integral.h>
#include <STEM4U/Utility.h>
#include <STEM4U/SeaWaves.h>
#include "functions.h"
#include "heal.h"
#include <Xlnt/Xlnt.h>
#include <Npy/Npy.h>

using namespace Eigen;

double Hydro::GetK_IRF_MaxT(const UVector<double> &w) {
	if (w.size() < 2)
		return -1;
	double delta = 0;
	int num = 0;
	for (int iw = 1; iw < w.size(); ++iw)
		if (w[iw] != w[iw-1]) {
			delta += w[iw] - w[iw-1];
			num++;
		}
	delta = delta/num;
		
	return M_PI/delta;		// (2*M_PI/delta)/2;
}

double Hydro::GetK_IRF_MaxT() const {
	return GetK_IRF_MaxT(dt.w);
}

void Hydro::GetK_IRF(double maxT, int numT) {
	if (dt.Nf == 0 || dt.B.IsEmpty())
		return;
	
    dt.Kirf.SetCount(dt.Nb*6); 			
    for (int i = 0; i < dt.Nb*6; ++i) {
    	dt.Kirf[i].SetCount(dt.Nb*6); 			 
   		for (int j = 0; j < dt.Nb*6; ++j)
			dt.Kirf[i][j].setConstant(numT, NaNDouble);
    }
		
	GetTirf(dt.Tirf, numT, maxT);
	
	UVector<double> y(dt.Nf);
  	for (int idf = 0; idf < dt.Nb*6; ++idf) {
    	for (int jdf = 0; jdf < dt.Nb*6; ++jdf) { 
			if (dt.B[idf][jdf].size() == 0 || !IsNum(dt.B[idf][jdf][0])) 
				continue;
			if (dt.dimen)
				GetKirf(dt.Kirf[idf][jdf], dt.Tirf, Get_w(), dt.B[idf][jdf]);
			else {
				GetKirf(dt.Kirf[idf][jdf], dt.Tirf, Get_w(), B_dim(idf, jdf));
				dt.Kirf[idf][jdf] /= g_rho_dim();
			}
    	}
  	}
}  

void Hydro::GetAinf() {
	if (dt.Nf == 0 || dt.A.size() < dt.Nb*6 || !IsLoadedKirf())
		return;	
	
	dt.Ainf.setConstant(dt.Nb*6, dt.Nb*6, NaNDouble);
	
    for (int i = 0; i < dt.Nb*6; ++i) 
        for (int j = 0; j < dt.Nb*6; ++j) 
    		if (IsNum(dt.Kirf[i][j][0]))
		    	dt.Ainf(i, j) = ::GetAinf(dt.Kirf[i][j], dt.Tirf, Get_w(), dt.A[i][j]);
}

void Hydro::GetRAO() {
	if (dt.Nf == 0 || dt.A.size() < dt.Nb*6 || dt.B.size() < dt.Nb*6)
		throw Exc(t_("Insufficient data to get RAO: Added mass and Radiation damping are required"));	

	for (int ib = 0; ib < dt.Nb; ++ib) 
		if (dt.msh[ib].dt.C.rows() < 6 || dt.msh[ib].dt.C.cols() < 6 || 
			dt.msh[ib].dt.M.rows() < 6 || dt.msh[ib].dt.M.cols() < 6) 
			throw Exc(t_("Insufficient data to get RAO: Hydrostatic stiffness and Inertia matrix are required"));   
			      
	Initialize_Forces(dt.rao);

	MatrixXd D = MatrixXd::Zero(6, 6);		// Unused for now
	MatrixXd D2 = MatrixXd::Zero(6, 6);
	
	for (int ib = 0; ib < dt.Nb; ++ib) {
		MatrixXd C = C_(false, ib);
		const MatrixXd &M_ = dt.msh[ib].dt.M;
		for (int ih = 0; ih < dt.Nh; ++ih) {	
			for (int ifr = 0; ifr < dt.Nf; ++ifr) {
				VectorXcd RAO = GetRAO(dt.w[ifr], A_mat(false, ifr, ib, ib), B_mat(false, ifr, ib, ib), 
								F_(false, dt.ex, ih, ifr, ib), C, M_, D, D2);
				for (int idf = 0; idf < 6; ++idf)
					dt.rao[ib][ih](ifr, idf) = RAO[idf];
			}
		}
	}
}

VectorXcd Hydro::GetRAO(double w, const MatrixXd &Aw, const MatrixXd &Bw, const VectorXcd &Fwh, 
		const MatrixXd &C, const MatrixXd &M, const MatrixXd &D, const MatrixXd &D2) {
	const std::complex<double> j = std::complex<double>(0, 1);

	MatrixXd Aw0 = Aw.unaryExpr([](double x){return IsNum(x) ? x : 0;}),		// Replaces Null with 0
			 Bw0 = Bw.unaryExpr([](double x){return IsNum(x) ? x : 0;});
	
	MatrixXcd m = C - sqr(w)*(M + Aw0) + j*w*(Bw0 + D);
	if (!FullPivLU<MatrixXcd>(m).isInvertible())
	   throw Exc(t_("Problem solving RAO"));
	
	return Fwh.transpose()*m.inverse();
}
	
void Hydro::InitAinf_w() {
	dt.Ainf_w.SetCount(dt.Nb*6); 			
    for (int i = 0; i < dt.Nb*6; ++i) {
    	dt.Ainf_w[i].SetCount(dt.Nb*6); 			 
   		for (int j = 0; j < dt.Nb*6; ++j)
			dt.Ainf_w[i][j].setConstant(dt.Nf, NaNDouble);
    }
}

void Hydro::GetAinf_w() {
	if (dt.Nf == 0 || dt.A.size() < dt.Nb*6 || !IsLoadedKirf())
		return;	
	
	InitAinf_w();
    
    for (int idf = 0; idf < dt.Nb*6; ++idf)
        for (int jdf = 0; jdf < dt.Nb*6; ++jdf) {
    		if (!IsLoadedB(idf, jdf)) 
                continue;
    		if (dt.dimen)
		    	::GetAinf_w(dt.Ainf_w[idf][jdf], dt.Kirf[idf][jdf], dt.Tirf, Get_w(), dt.A[idf][jdf]);
            else {
                ::GetAinf_w(dt.Ainf_w[idf][jdf], dt.Kirf[idf][jdf]*g_rho_dim(), dt.Tirf, Get_w(), A_dim(idf, jdf));
                dt.Ainf_w[idf][jdf] *= (1/(rho_dim()*pow(dt.len, GetK_AB(idf, jdf))));
            }
        }
}

void Hydro::GetB_H(int &num) {
	if (!IsLoadedFex())
		throw Exc(t_("The excitation force is not loaded"));
	if (dt.Nf == 0)
		throw Exc(t_("No frecuencies loaded"));
	if (num > dt.Nh)
		throw Exc(t_("Number of headings is higher than available"));
	if (num <= 0)
		throw Exc(t_("Not enough headings for Haskind calculation"));
	if (IsNull(dt.h))
		throw Exc(t_("Unknown depth"));

	UVector<double> head360 = clone(dt.head);

	double angle = 360/num;
	double nextAngle = head360[0] + angle;
	UVector<int> idRemove;
	for (int i = 1; i < head360.size(); ++i) {
		if (head360[i] < nextAngle) 
			idRemove << i;
		else
			nextAngle += angle;
	}
	num = dt.head.size() - idRemove.size();
	for (int i = idRemove.size()-1; i >= 0; --i)
		head360.Remove(idRemove[i]);
	
	String shead = FormatF(head360[0], 1);
	for (int i = 1; i < head360.size(); ++i) 
		shead << ", " + FormatF(head360[i], 1);
	BEM::Print("\n" + Format(t_("Haskind got for %d headings %s"), num, shead));

	enum RangeType {R_0_360, R_x_360, R_0_x, R_x_x};	// There must be data from 0 to 360 deg
	RangeType rangeType;
	double hd0 = head360[0],
		   hdl = Last(head360);
	if (hd0 == 0 && hdl < 360) {
		head360 << 360;
		rangeType = R_0_x;
	} else if (hd0 > 0 && hdl == 360) {
		head360.Insert(0, 0);
		rangeType = R_x_360;
	} else if (hd0 > 0 && hdl < 360) {
		head360.Insert(0, 0);
		head360 << 360;
		rangeType = R_x_x;
	} else	
		rangeType = R_0_360;	

	VectorXd val(dt.Nf);
	for (int ifr = 0; ifr < dt.Nf; ++ifr)
	 	val[ifr] = dt.w[ifr]*SeaWaves::WaveNumber_w(dt.w[ifr], -1, g_dim(), true)/(4*M_PI*rho_dim()*g_dim()*g_dim());
	
	Initialize_AB(dt.B_H);
	
	for (int ib = 0; ib < dt.Nb; ++ib) {
	    for (int idf = 0; idf < 6; ++idf) {
			if (!IsLoadedFex(idf, 0, ib)) 		
	            continue;
			
			VectorXd b(dt.Nf);
			for (int ifr = 0; ifr < dt.Nf; ++ifr) {
				UVector<double> F2(dt.Nh);
				for (int ih = 0; ih < dt.Nh; ++ih) 
					F2[ih] = sqr(F_dim(abs(dt.ex[ib][ih](ifr, idf)), idf));
				
				for (int i = idRemove.size()-1; i >= 0; --i)
					F2.Remove(idRemove[i]);
				
				if (rangeType == R_0_x) 					// There must be data from 0 to 360 deg
					F2 << F2[0];
				else if (rangeType == R_x_360) 
					F2.Insert(0, Last(F2));
				else if (rangeType == R_x_x) {
					double f = LinearInterpolate(360., hdl, 360 + hd0, Last(F2), F2[0]);
					F2.Insert(0, f);
					F2 << f;
				}
				b(ifr) = Integral(head360, F2, SIMPSON_1_3)*val[ifr]*M_PI/180;
			}
			if (dt.dimen)
				dt.B_H[idf][idf] = b*rho_ndim()/rho_dim();
			else {
				for (int ifr = 0; ifr < dt.Nf; ++ifr)
					b(ifr) /= (rho_dim()*pow(dt.len, GetK_AB(idf, idf))*dt.w[ifr]);
				dt.B_H[idf][idf] = b;
			}
	    }
    }
}

void Hydro::GetOgilvieCompliance(bool zremoval, bool thinremoval, bool decayingTail, UVector<int> &vidof, UVector<int> &vjdof) {
	vidof.Clear();
	vjdof.Clear();
	if (dt.Nf == 0 || dt.A.size() < dt.Nb*6)
		return;	
	
	HealBEM data;
	
	if (dt.Ainf.size() == 0) 
		dt.Ainf.setConstant(6*dt.Nb, 6*dt.Nb, NaNDouble); 
	
	if (dt.Ainf_w.size() == 0) {
		dt.Ainf_w.SetCount(dt.Nb*6); 			
	    for (int idf = 0; idf < dt.Nb*6; ++idf) {
    		dt.Ainf_w[idf].SetCount(dt.Nb*6); 			 
   			for (int jdf = 0; jdf < dt.Nb*6; ++jdf)
				dt.Ainf_w[idf][jdf].setConstant(dt.Nf, NaNDouble);
	    }
    }
    double maxT = min(Bem().maxTimeA, Hydro::GetK_IRF_MaxT(dt.w));
    int numT = Bem().numValsA;
    
    if (dt.Kirf.size() == 0) {
        dt.Kirf.SetCount(dt.Nb*6); 			
	    for (int idf = 0; idf < dt.Nb*6; ++idf) {
	    	dt.Kirf[idf].SetCount(dt.Nb*6); 			 
	   		for (int jdf = 0; jdf < dt.Nb*6; ++jdf)
				dt.Kirf[idf][jdf].setConstant(numT, NaNDouble);
	    }
	}
		
    for (int idf = 0; idf < dt.Nb*6; ++idf) {
        MatrixXd ex_hf(dt.Nh, dt.Nf);
        
        for (int jdf = 0; jdf < dt.Nb*6; ++jdf) {
    		if (dt.B[idf][jdf].size() == 0 || !IsNum(dt.B[idf][jdf][0])) 
                ;
            else {
                bool done;
	    		if (data.Load(Get_w(), A_dim(idf, jdf), Ainf_dim(idf, jdf), B_dim(idf, jdf), numT, maxT, ex_hf) &&
					data.Heal(zremoval, thinremoval, decayingTail, done)) {
	            	data.Save(dt.A[idf][jdf], dt.Ainf_w[idf][jdf], dt.Ainf(idf, jdf), dt.B[idf][jdf], dt.Tirf, dt.Kirf[idf][jdf]); 
	            	if (done) {
		            	vidof << idf;
		            	vjdof << jdf;
	            	}
	    		} else
	    			data.Reset(dt.A[idf][jdf], dt.Ainf_w[idf][jdf], dt.Ainf(idf, jdf), dt.B[idf][jdf], dt.Kirf[idf][jdf]);
	    		if (dt.dimen) {
	    			dt.dimen = false;
	    			dt.A[idf][jdf] 	= pick(A_ndim(idf, jdf));
	    			dt.Ainf_w[idf][jdf] *= (rho_ndim()/rho_dim());
	    			dt.Ainf(idf, jdf)   *= (rho_ndim()/rho_dim());
	    			dt.B[idf][jdf] 	= pick(B_ndim(idf, jdf));
	    			dt.Kirf[idf][jdf]  = pick(Kirf_ndim(idf, jdf));
	    			dt.dimen = true;
	    		} else {
	    			dt.dimen = true;
	    			dt.A[idf][jdf] 	= pick(A_ndim(idf, jdf));
	    			dt.Ainf_w[idf][jdf] *= (1/(rho_ndim()*pow(dt.len, GetK_AB(idf, jdf))));
	    			dt.Ainf(idf, jdf)   *= (1/(rho_ndim()*pow(dt.len, GetK_AB(idf, jdf))));
	    			dt.B[idf][jdf] 	= pick(B_ndim(idf, jdf));
	    			dt.Kirf[idf][jdf] 	= pick(Kirf_ndim(idf, jdf));
	    			dt.dimen = false;
	    		}
            }
        }
    }
    dt.rao.Clear();	// Previous RAO is now invalid
}

void Hydro::GetWaveTo(double xto, double yto, double g) {
	double dx = xto - dt.x_w,
		   dy = yto - dt.y_w;
	
	for (int ib = 0; ib < dt.Nb; ++ib) 
		AddWave(ib, dx, dy, g);		// This g is imposed by loaded file when loading, and by BEM when saving
	
	dt.x_w = xto;
	dt.y_w = yto;
}

String Hydro::SpreadNegative(Function <bool(String, int)> Status) {
	String ret;
	UVector<String> errors;

	int numT = 0, num = 0;
	for (int ib = 0; ib < dt.Nb; ++ib) 
		numT += dt.pots_rad[ib].size();

	for (int ib = 0; ib < dt.Nb; ++ib) {
		for (int idp = 0; idp < dt.pots_rad[ib].size(); ++idp) {
			int adv = 100*num/numT;
			if (Status && !(adv%2)) {
				if (!Status(t_("Spreading negative values in diagonal"), adv))
					throw Exc(t_("Stop by user"));
			}
			num++;	
				
			UVector<int> panIDs;
			for (int ifr = 0; ifr < dt.Nf; ++ifr) {
				for (int idf = 0; idf < 6; ++idf) {		// Only diagonal
					double &apan = dt.Apan(ib, idp, idf, idf, ifr);
					if (apan < 0) {
						if (panIDs.IsEmpty())			// It is only get if the added mass in any dof is negative
							dt.msh[ib].dt.mesh.GetClosestPanels(idp, panIDs);
						for (int i = 0; i < panIDs.size(); ++i) {
							double &apan_i = dt.Apan(ib, panIDs[i], idf, idf, ifr);
							if (apan_i > 0) {
								if (apan_i + apan >= 0) {	// apan_i has enough mass
									apan_i += apan;
									apan = 0;
									break;
								} else {					// apan_i has not enough mass
									apan += apan_i;
									apan_i = 0;
								}
							}
						}
						if (apan < 0) {			// Impossible to spread
							FindAdd(errors, Format(t_("%d.%s Freq %.3f rad/s"), ib+1, BEM::strDOFtext[idf], dt.w[ifr]));
							// Resets this dof and frequency for all panels
							for (int i = 0; i < dt.pots_rad[ib].size(); ++i) 
								dt.Apan(ib, i, idf, idf, ifr) = 0;
						}
					}
				}
			}
		}
	}
	Sort(errors);
	for (const String &s : errors) {
		if (!ret.IsEmpty())
			ret << "\n";
		ret << s;
	}
	return ret;
}

void Hydro::MapNodes(int ib, UVector<Point3D> &points, Tensor<double, 4> &Apan_nodes, Tensor<double, 4> &Bpan_nodes) const {
	if (points.IsEmpty())
		throw Exc(t_("No points to map"));
	
	Apan_nodes.resize(points.size(), 6, 6, dt.Nf);
	Apan_nodes.setZero();
	Bpan_nodes.resize(points.size(), 6, 6, dt.Nf);
	Bpan_nodes.setZero();
	
	auto GetClosest = [](const Point3D &p, const UVector<Point3D> &points)->int {
		double dmin = std::numeric_limits<double>::max();
		int ipmin = -1;
		for (int ip = 0; ip < points.size(); ++ip) {
			double d = Distance(p, points[ip]);
			if (d < dmin) {
				dmin = d;
				ipmin = ip;
			}
		}
		return ipmin;
	};
	
	for (int idp = 0; idp < dt.pots_rad[ib].size(); ++idp) {
		const Point3D &p = dt.msh[ib].dt.mesh.panels[idp].centroidPaint;
		int ip = GetClosest(p, points);
		for (int ifr = 0; ifr < dt.Nf; ++ifr) {
			for (int idf1 = 0; idf1 < 6; ++idf1) {		
				for (int idf2 = 0; idf2 < 6; ++idf2) {
					Apan_nodes(ip, idf1, idf2, ifr) += dt.Apan(ib, idp, idf1, idf2, ifr);
					Bpan_nodes(ip, idf1, idf2, ifr) += B_pan(ib, idp, idf1, idf2, ifr);
				}
			}
		}
	}
}

void Hydro::SaveMap(Grid &g, int ifr, bool onlyDiagonal, const UVector<int> &ids, const UVector<Point3D> &points, 
		const Tensor<double, 4> &Apan, const Tensor<double, 4> &Bpan) const {
	g.SetNumHeaderRows(1);
	g.SetRow(0);
	int col = 0;
	g.Set(Null, col++, t_("Id"));	g.AddCol(30);
	g.Set(Null, col++, t_("x"));	g.AddCol(60);
	g.Set(Null, col++, t_("y"));	g.AddCol(60);
	g.Set(Null, col++, t_("z"));	g.AddCol(60);

	if (Apan.size() > 0) {
		if (onlyDiagonal) {
			for (int c = 0; c < 6; ++c) {
				g.Set(Null, col++, Format(t_("A_%s_%s"), BEM::StrDOF(c), BEM::StrDOF(c)));		g.AddCol(60);
			}
			for (int c = 0; c < 6; ++c) {
				g.Set(Null, col++, Format(t_("B_%s_%s"), BEM::StrDOF(c), BEM::StrDOF(c)));		g.AddCol(60);
			}
		} else {
			for (int r = 0; r < 6; ++r) 
				for (int c = 0; c < 6; ++c) {
					g.Set(Null, col++, Format(t_("A_%s_%s"), BEM::StrDOF(r), BEM::StrDOF(c)));	g.AddCol(60);
				}
			for (int r = 0; r < 6; ++r) 
				for (int c = 0; c < 6; ++c) {
					g.Set(Null, col++, Format(t_("B_%s_%s"), BEM::StrDOF(r), BEM::StrDOF(c)));	g.AddCol(60);
				}
		}
	}
			
	for (int row = 0; row < ids.size(); ++row) {
		col = 0;
		g.SetRow(row+1);
		g.Set(Null, col++, ids[row]);
		g.Set(Null, col++, points[row].x);
		g.Set(Null, col++, points[row].y);
		g.Set(Null, col++, points[row].z);
		
		if (Apan.size() > 0) {
			if (onlyDiagonal) {
				for (int c = 0; c < 6; ++c) 	
					g.Set(Null, col++, Apan(row, c, c, ifr));
				for (int c = 0; c < 6; ++c) 	
					g.Set(Null, col++, Bpan(row, c, c, ifr));
			} else {
				//for (int row = 0; row < ids.size(); ++row) {
					for (int r = 0; r < 6; ++r) 
						for (int c = 0; c < 6; ++c)
							g.Set(Null, col++, Apan(row, r, c, ifr));	
					for (int r = 0; r < 6; ++r) 
						for (int c = 0; c < 6; ++c)
							g.Set(Null, col++, Bpan(row, r, c, ifr));		
				//}
			}
		}
	}
}
		
void Hydro::SaveMap(String fileName, String type, int ifr, bool onlyDiagonal, const UVector<int> &ids, const UVector<Point3D> &points, 
		const Tensor<double, 4> &Apan, const Tensor<double, 4> &Bpan) const {
	if (IsNull(ifr)) {
		UArray<Grid> grid(dt.Nf);
		for (ifr = 0; ifr < dt.Nf; ++ifr) 
			SaveMap(grid[ifr], ifr, onlyDiagonal, ids, points, Apan, Bpan);
		
		if (type == ".csv") {
			for (ifr = 0; ifr < dt.Nf; ++ifr) { 
				String folder = GetFileFolder(fileName);
				String name = GetFileTitle(fileName);
				String ext = GetFileExt(fileName);
				String fname = AFX(folder, Format("%s_%.3f%s", name, dt.w[ifr], ext));
				SaveFile(fname, grid[ifr].AsString(false, false, ScatterDraw::GetDefaultCSVSeparator()));
			}
		} else if (type == ".xlsx") {
			xlnt::workbook wb;
			for (ifr = 0; ifr < dt.Nf; ++ifr) { 
				xlnt::worksheet ws;
				String title = Format("%.3f", dt.w[ifr]); 
				if (ifr == 0)
					ws = wb.active_sheet();	
				else
					ws = wb.create_sheet();
				ws.title(~title);
				XlsxFill(ws, grid[ifr], false);
			}
			wb.save(~fileName); 
		}
	} else {
		Grid grid;
		SaveMap(grid, ifr, onlyDiagonal, ids, points, Apan, Bpan);
		
		if (type == ".csv") 
			SaveFile(fileName, grid.AsString(false, false, ScatterDraw::GetDefaultCSVSeparator()));
		else if (type == ".xlsx") {
			xlnt::workbook wb;
			xlnt::worksheet ws = wb.active_sheet();	
			ws.title("Data");
			XlsxFill(ws, grid, false);
			wb.save(~fileName); 
		}
	}
}

void Hydro::SaveAkselos(int ib, String file) const {
	if (dt.msh.IsEmpty() || dt.msh[0].dt.mesh.panels.IsEmpty())
		throw Exc(t_("No mesh is available"));
	if (!IsLoadedPotsRad(ib)) 
		throw Exc(Format(t_("No radiation potentials/pressures are available for body %d"), ib+1));

	String folder = GetFileFolder(file);
	file = GetFileTitle(file);
	
	{
		Grid g;
		g.SetRow({"Id", "heading"});
		for (int ih = 0; ih < dt.Nh; ++ih) 
			g.SetRow({ih, dt.head[ih]});
		SaveFile(AFX(folder, Format("%s_Headings.csv", file)), g.AsString(false, false, ","));
	}{
		Grid g;
		g.SetRow({"Id", "period"});
		VectorXd T = Get_T();
		Sort(T);
		for (int iT = 0; iT < dt.Nf; ++iT) 
			g.SetRow({iT, T[iT]});
		SaveFile(AFX(folder, Format("%s_Periods.csv", file)), g.AsString(false, false, ","));
	}{
		Grid g;
		g.SetRow({"Id", "period"});
		for (int idf = 0; idf < dt.Nf; ++idf) 
			g.SetRow({idf, dt.w[idf]/2/M_PI});
		SaveFile(AFX(folder, Format("%s_Freq.csv", file)), g.AsString(false, false, ","));
	}{
		Grid g;
		for (int r = 0; r < 6; ++r)
			for (int c = 0; c < 6; ++c)
				g.Set(r, c, dt.msh[ib].dt.C(r, c));
		SaveFile(AFX(folder, "Hydrostatic_Results.csv"), g.AsString(false, false, ","));
	}{
		Grid g;
		const UVector<String> sdof = {"x", "y", "z"};
		UVector<Value> header = {"Panel index", "Id", "Area", "centroidX", "centroidY", "centroidZ", "normalX", "normalY", "normalZ"};
		for (int iv = 0; iv < 4; ++iv)
			for (int idof = 0; idof < 3; ++idof)	
				header << Format("Vertex_%d_%s", iv+1, sdof[idof]);
		g.SetRow(header);
		for (int ip = 0; ip < dt.msh[ib].dt.mesh.panels.size(); ++ip) {
			const Panel &p = dt.msh[ib].dt.mesh.panels[ip];
			const UVector<Point3D> &nodes = dt.msh[ib].dt.mesh.nodes;
			g.Set(ip).Set(ib).Set(p.surface0 + p.surface1);
			g.Set(p.centroidPaint.x).Set(p.centroidPaint.y).Set(p.centroidPaint.z);
			g.Set(p.normalPaint.x).Set(p.normalPaint.y).Set(p.normalPaint.z);
			for (int iv = 0; iv < 4; ++iv) {
				const Point3D &n = nodes[p.id[iv]];
				g.Set(n.x).Set(n.y).Set(n.z);
			}
			g.NextRowLF();
		}
		SaveFile(AFX(folder, Format("%s_Panel.csv", file)), g.AsString(false, false, ","));
	}
	
	MultiDimMatrixRowMajor<std::complex<double>> rad(6, dt.Nf, dt.msh[ib].dt.mesh.panels.size());
	//rad.SetZero();
	for (int idof = 0; idof < 6; ++idof)
		for (int idf = 0; idf < dt.Nf; ++idf)	
			for (int ip = 0; ip < dt.msh[ib].dt.mesh.panels.size(); ++ip)	
				rad(idof, idf, ip) = P_rad(ib, ip, idof, dt.Nf - idf - 1);	// Inverse order to save in period order

	MultiDimMatrixRowMajor<std::complex<double>> dif(dt.Nh, dt.Nf, dt.msh[ib].dt.mesh.panels.size());
	//dif.SetZero();
	for (int ih = 0; ih < dt.Nh; ++ih)
		for (int idf = 0; idf < dt.Nf; ++idf)	
			for (int ip = 0; ip < dt.msh[ib].dt.mesh.panels.size(); ++ip)	
				dif(ih, idf, ip) = P_dif(ib, ip, ih, dt.Nf - idf - 1);		// Inverse order to save in period order

	Npz npz;
	Npy &nrad = npz.Add("radiation");
	nrad.Set(rad);
	Npy &ndif = npz.Add("diffraction");
	ndif.Set(dif);
	npz.Save(AFX(folder, Format("%s_Pressure.npz", file)));
}

void Hydro::AddWave(int ib, double dx, double dy, double g) {
	if (dx == 0 && dy == 0)
		return;
  	auto CalcF = [&](Forces &ex, const UVector<double> &k) {
    	Forces exforce = clone(ex);
    	
	    for (int ih = 0; ih < dt.Nh; ++ih) {
	        double angle = ToRad(dt.head[ih]);
			double dist = dx*cos(angle) + dy*sin(angle);
		
			for (int ifr = 0; ifr < dt.Nf; ++ifr) {
				double ph = k[ifr]*dist;
				for (int idf = 0; idf < 6; ++idf) 
					AddPhase(exforce[ib][ih](ifr, idf), ph);		// Add the phase
			}
	    }
		ex = pick(exforce);
    };
    
    UVector<double> k(dt.Nf);
	for (int ifr = 0; ifr < dt.Nf; ++ifr) 
		k[ifr] = SeaWaves::WaveNumber_w(dt.w[ifr], dt.h, g);
    	
	if (IsLoadedFex())
		CalcF(dt.ex, k);
	if (IsLoadedFsc())
		CalcF(dt.sc, k);
	if (IsLoadedFfk())
		CalcF(dt.fk, k);
	if (IsLoadedFsc_pot())
		CalcF(dt.sc_pot, k);
	if (IsLoadedFfk_pot())
		CalcF(dt.fk_pot, k);
	if (IsLoadedFfk_pot_bmr())
		CalcF(dt.fk_pot_bmr, k);
	
	if (IsLoadedRAO())
		CalcF(dt.rao, k);
	
    auto CalcQTF = [&](UArray<UArray<UArray<MatrixXcd>>> &qtf, const UVector<double> &qk, bool isSum) {
        int sign = isSum ? 1 : -1;
		{
	        for (int ih = 0; ih < dt.qhead.size(); ++ih) {
	            double angle = ToRad(dt.qhead[ih].imag());
	            double dist = dx*cos(angle) + dy*sin(angle);
				for (int ifr1 = 0; ifr1 < dt.qw.size(); ++ifr1) {
					for (int ifr2 = 0; ifr2 < dt.qw.size(); ++ifr2) {
						double ph = (qk[ifr2] + sign*qk[ifr1])*dist;
						for (int idf = 0; idf < 6; ++idf) 
							AddPhase(qtf[ib][ih][idf](ifr1, ifr2), -ph);
					}
				}
	        }
		}
    };

    UVector<double> qk(dt.Nf);
	for (int ifr = 0; ifr < dt.qw.size(); ++ifr) 
		qk[ifr] = SeaWaves::WaveNumber_w(dt.qw[ifr], dt.h, g_dim());
			
	if (IsLoadedQTF(true)) 
		CalcQTF(dt.qtfsum, qk, true);		
	if (IsLoadedQTF(false))	
		CalcQTF(dt.qtfdif, qk, false);	
}


void Hydro::GetTranslationTo(const MatrixXd &to, bool force, Function <bool(String, int pos)> Status) {
	MatrixXd delta(3, dt.Nb);
	for (int ib = 0; ib < dt.Nb; ++ib) 
		for (int idf = 0; idf < 3; ++idf) 	
			delta(idf, ib) = to(idf, ib) - dt.msh[ib].dt.c0[idf];
	
	auto Nvl0 = [](Matrix3d &mat) {
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				mat(i, j) = ::Nvl(mat(i, j), 0.);		
	};
	auto CopyFrom = [](const UArray<UArray<VectorXd>> &a, int i0, int j0, int iif)->Matrix3d {
		Matrix3d ret;
		
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				ret(i, j) = a[i0+i][j0+j][iif];
		return ret;
	};
	auto CopyTo = [](const Matrix3d &from, UArray<UArray<VectorXd>> &a, int i0, int j0, int iif){
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				a[i0+i][j0+j][iif] = from(i, j);
	};
		
	auto TransAB = [&](UArray<UArray<VectorXd>> &A) {
        UArray<UArray<VectorXd>> An;
        Initialize_AB(An);

		for (int ib = 0; ib < dt.Nb; ++ib) {
			int ib6 = ib*6;
			double dx = delta(0, ib);
			double dy = delta(1, ib);
			double dz = delta(2, ib);
			
			Matrix3d Rg {{  0, -dz,  dy},
						 { dz,   0, -dx},
						 {-dy,  dx,   0}};
			
			for (int jb = 0; jb < dt.Nb; ++jb) {
				int jb6 = jb*6;
				
				if (!force) {
					for (int idof = 0; idof < 6; ++idof) {		// All dof are available?
						for (int jdof = 0; jdof < 6; ++jdof)
							if (!IsNum(A[ib6 + idof][jb6 + jdof][0])) 
								throw Exc("Coefficient translations require all DOFs to be available");
					}
				}
				
				for (int iif = 0; iif < dt.Nf; ++iif) {
				  	Matrix3d Q11 = CopyFrom(A, ib6 + 0, jb6 + 0, iif);
	  				Matrix3d Q12 = CopyFrom(A, ib6 + 0, jb6 + 3, iif);
					Matrix3d Q21 = CopyFrom(A, ib6 + 3, jb6 + 0, iif);
					Matrix3d Q22 = CopyFrom(A, ib6 + 3, jb6 + 3, iif);
					
					Nvl0(Q11);	Nvl0(Q12);	Nvl0(Q21);	Nvl0(Q22);
					
					CopyTo(Q11, 							  An, ib6 + 0, jb6 + 0, iif);
					CopyTo(Q12 + Q11*Rg, 					  An, ib6 + 0, jb6 + 3, iif);
					CopyTo(Q21 - Rg*Q11, 					  An, ib6 + 3, jb6 + 0, iif);
					CopyTo(Q22 - Rg*Q12 + Q21*Rg - Rg*Q11*Rg, An, ib6 + 3, jb6 + 3, iif);
				}
			}
		}
		A = pick(An);
    };
	
	Status(t_("Translating A"), 10);
	if (IsLoadedA())
		TransAB(dt.A);
	if (IsLoadedAinf_w())
		TransAB(dt.Ainf_w);

	Status(t_("Translating B"), 20);
	if (IsLoadedB())
		TransAB(dt.B);
	if (IsLoadedB_H())
		TransAB(dt.B_H);
	
    auto TransA = [&](MatrixXd &A) {
        MatrixXd An(6,6);
	
		for (int ib = 0; ib < dt.Nb; ++ib) {
			int ib6 = ib*6;
			double dx = delta(0, ib);
			double dy = delta(1, ib);
			double dz = delta(2, ib);
			
			Matrix3d Rg {{  0, -dz,  dy},
						 { dz,   0, -dx},
						 {-dy,  dx,   0}};
			
			for (int jb = 0; jb < dt.Nb; ++jb) {
				int jb6 = jb*6;
				
				if (!force) {
					for (int idof = 0; idof < 6; ++idof) {
						for (int jdof = 0; jdof < 6; ++jdof)
							if (!IsNum(A(ib6 + idof, jb6 + jdof))) 
								throw Exc("Coefficient translations require all DOFs to be available");
					}
				}
				
  				Matrix3d Q11 = A.block<3,3>(ib6 + 0, jb6 + 0); 
  				Matrix3d Q12 = A.block<3,3>(ib6 + 0, jb6 + 3);
				Matrix3d Q21 = A.block<3,3>(ib6 + 3, jb6 + 0);
				Matrix3d Q22 = A.block<3,3>(ib6 + 3, jb6 + 3);
				
				Nvl0(Q11);	Nvl0(Q12);	Nvl0(Q21);	Nvl0(Q22);
				
				An.block<3,3>(ib6 + 0, jb6 + 0) = Q11;
				An.block<3,3>(ib6 + 0, jb6 + 3) = Q12 + Q11*Rg;
				An.block<3,3>(ib6 + 3, jb6 + 0) = Q21 - Rg*Q11;
				An.block<3,3>(ib6 + 3, jb6 + 3) = Q22 - Rg*Q12 + Q21*Rg - Rg*Q11*Rg;
			} 
		}
		A = pick(An);
    };
    
    if (IsLoadedA0())
		TransA(dt.A0);
    if (IsLoadedAinf())
		TransA(dt.Ainf);
		    
    auto TransF = [&](Forces &ex) {
    	Forces exforce = clone(ex);
    	
	    for (int ih = 0; ih < dt.Nh; ++ih) {
	    	for (int ib = 0; ib < dt.Nb; ++ib) {
	    		double dx = delta(0, ib);
				double dy = delta(1, ib);
				double dz = delta(2, ib);
			
				for (int ifr = 0; ifr < dt.Nf; ++ifr) {
					exforce[ib][ih](ifr, 3) += -dy*exforce[ib][ih](ifr, 2) + dz*exforce[ib][ih](ifr, 1);
	    			exforce[ib][ih](ifr, 4) += -dz*exforce[ib][ih](ifr, 0) + dx*exforce[ib][ih](ifr, 2);
	    			exforce[ib][ih](ifr, 5) += -dx*exforce[ib][ih](ifr, 1) + dy*exforce[ib][ih](ifr, 0);
				}
	    	}
	    }
		ex = pick(exforce);
    };
    
	Status(t_("Translating Forces"), 30);	
	if (IsLoadedFex())
		TransF(dt.ex);
	if (IsLoadedFsc())
		TransF(dt.sc);
	if (IsLoadedFfk())
		TransF(dt.fk);

    auto TransMD = [&]() {
    	UArray<UArray<UArray<VectorXd>>> mdn = clone(dt.md);
    	
    	for (int ib = 0; ib < dt.Nb; ++ib) {
    		double dx = delta(0, ib);
			double dy = delta(1, ib);
			double dz = delta(2, ib);    		
    		for (int ih = 0; ih < dt.mdhead.size(); ++ih) {
    			for (int ifr = 0; ifr < dt.Nf; ++ifr) {
    				mdn[ib][ih][3](ifr) += -dy*mdn[ib][ih][2](ifr) + dz*mdn[ib][ih][1](ifr);
    				mdn[ib][ih][4](ifr) += -dz*mdn[ib][ih][0](ifr) + dx*mdn[ib][ih][2](ifr);
    				mdn[ib][ih][5](ifr) += -dx*mdn[ib][ih][1](ifr) + dy*mdn[ib][ih][0](ifr);
				}
	    	}
	    }
		dt.md = pick(mdn);
    };
    
	if (IsLoadedMD()) {
		Status(t_("Translating MD"), 40);
		TransMD();
	}
	
    auto TransQTF = [&](UArray<UArray<UArray<MatrixXcd>>> &qtf, bool isSum) {
		for (int ib = 0; ib < dt.Nb; ++ib) {
			double dx = delta(0, ib);
			double dy = delta(1, ib);
			double dz = delta(2, ib);
	        for (int ih = 0; ih < dt.qhead.size(); ++ih) {
				for (int ifr1 = 0; ifr1 < dt.qw.size(); ++ifr1) {
					for (int ifr2 = 0; ifr2 < dt.qw.size(); ++ifr2) {
						std::complex<double> &v0 = qtf[ib][ih][0](ifr1, ifr2),
							 				 &v1 = qtf[ib][ih][1](ifr1, ifr2),
							 				 &v2 = qtf[ib][ih][2](ifr1, ifr2),
							 				 &v3 = qtf[ib][ih][3](ifr1, ifr2),
							 				 &v4 = qtf[ib][ih][4](ifr1, ifr2),
							 				 &v5 = qtf[ib][ih][5](ifr1, ifr2);
						
						v3 += -dy*v2 + dz*v1;
						v4 += -dz*v0 + dx*v2;
						v5 += -dx*v1 + dy*v0;
					}
				}
	        }
		}
    };

	// QTF translation only valid for same headings. Crossed headings are deleted
	for (int ih = int(dt.qhead.size())-1; ih >= 0; --ih) { 
		if (dt.qhead[ih].real() != dt.qhead[ih].imag()) {
			Remove(dt.qhead, ih);
			for (int ib = 0; ib < dt.Nb; ++ib) {
				if (IsLoadedQTF(true)) 
					dt.qtfsum[ib].Remove(ih);
				if (IsLoadedQTF(false)) 
					dt.qtfdif[ib].Remove(ih);
			}
		}
	}
		
	if (IsLoadedQTF(true)) {
		Status(t_("Translating QTF"), 50);
		TransQTF(dt.qtfsum, true);		
	}
	if (IsLoadedQTF(false))	{
		Status(t_("Translating QTF"), 60);
		TransQTF(dt.qtfdif, false);
	}
	
	if (IsLoadedPotsRad()) {
		for (int ib = 0; ib < dt.Nb; ++ib) {
			for	(int ip = 0; ip < dt.pots_rad[ib].size(); ++ip) {
				UArray<UArray<std::complex<double>>> &pot = dt.pots_rad[ib][ip];
				for	(int ifr = 0; ifr < dt.Nf; ++ifr) {	
					pot[3][ifr] -= (pot[2][ifr]*delta(1, ib) - pot[1][ifr]*delta(2, ib));
					pot[4][ifr] -= (pot[0][ifr]*delta(2, ib) - pot[2][ifr]*delta(0, ib));
					pot[5][ifr] -= (pot[1][ifr]*delta(0, ib) - pot[0][ifr]*delta(1, ib));
				}
			}
		}
	}
	
	if (IsLoadedM()) {
		for (int ib = 0; ib < dt.Nb; ++ib) 
			Surface::TranslateInertia66(dt.msh[ib].dt.M, dt.msh[ib].dt.cg, dt.msh[ib].dt.c0, Point3D(to.col(ib)));
	}
	
	if (IsLoadedKirf()) {
		double maxT = GetK_IRF_MaxT();
		if (maxT < 0)
			maxT = Bem().maxTimeA;
		else if (Bem().maxTimeA < maxT) 
			maxT = Bem().maxTimeA;

		GetK_IRF(maxT, Bem().numValsA);
	}
	
	// Some previous data are now invalid
	dt.rao.Clear();	
	
	for (int ib = 0; ib < dt.Nb; ++ib) {
		dt.msh[ib].dt.C = MatrixXd();
		dt.msh[ib].dt.Dlin = MatrixXd();
	}
	
	for (int ib = 0; ib < dt.Nb; ++ib)
		dt.msh[ib].dt.c0 = to.col(ib);
	
	String error = AfterLoad();
	if (!error.IsEmpty())
		throw Exc(Format(t_("Problem translating model: '%s'\n%s"), error));	
}

void Hydro::CompleteForces1st() {
	if (!IsLoadedFex() && IsLoadedFsc() && IsLoadedFfk()) 
		GetFexFromFscFfk();
	if (!IsLoadedFsc() && IsLoadedFex() && IsLoadedFfk()) 
		GetFscFromFexFfk();
	if (!IsLoadedFfk() && IsLoadedFex() && IsLoadedFsc()) 
		GetFfkFromFexFsc();	
}

void Hydro::ResetForces1st(Hydro::FORCE force) {
	if (force == Hydro::FK) {
		if (!IsLoadedFfk())
			return;
		if (!IsLoadedFsc() && !IsLoadedFex())
			return;
		
		if (IsLoadedFsc()) 
			dt.ex = clone(dt.sc);
		else {
			for (int ib = 0; ib < dt.Nb; ++ib) 
				for (int ih = 0; ih < dt.Nh; ++ih) 
					for (int ifr = 0; ifr < dt.Nf; ++ifr) 
						for (int i = 0; i < 6; ++i) 
							if (IsNum(dt.fk[ib][ih](ifr, i))) 
								dt.ex[ib][ih](ifr, i) = dt.ex[ib][ih](ifr, i) - dt.fk[ib][ih](ifr, i);
		}
		dt.fk.Clear();
	} else if (force == Hydro::SCATTERING) {
		if (!IsLoadedFsc())
			return;
		if (!IsLoadedFfk() && !IsLoadedFex())
			return;
		
		if (IsLoadedFfk()) 
			dt.ex = clone(dt.fk);
		else {
			for (int ib = 0; ib < dt.Nb; ++ib) 
				for (int ih = 0; ih < dt.Nh; ++ih) 
					for (int ifr = 0; ifr < dt.Nf; ++ifr) 
						for (int i = 0; i < 6; ++i) 
							if (IsNum(dt.sc[ib][ih](ifr, i))) 
								dt.ex[ib][ih](ifr, i) = dt.ex[ib][ih](ifr, i) - dt.sc[ib][ih](ifr, i);
		}
		dt.sc.Clear();		
	} else {
		dt.ex.Clear();		
		dt.sc.Clear();		
		dt.fk.Clear();		
	}
	dt.md.Clear();
}

void Hydro::ResetForces(Hydro::FORCE force, bool forceMD, Hydro::FORCE forceQtf) {
	if (force != Hydro::NONE)
		Hydro::ResetForces1st(force);

	if (forceMD)
		dt.md.Clear();
	
	if (forceQtf == Hydro::ALL || forceQtf == Hydro::QTFSUM) 
		dt.qtfsum.Clear();
	if (forceQtf == Hydro::ALL || forceQtf == Hydro::QTFDIF) 
		dt.qtfdif.Clear();
}

void Hydro::MultiplyDOF(double factor, const UVector<int> &_idDOF, bool a, bool b, bool diag, bool f, bool ismd, bool qtf) {
	if (_idDOF.size() == 0) 
		return;
	
	UVector<int> idDOF;
	for (int idof = 0; idof < _idDOF.size(); ++idof)
		for (int ib = 0; ib < dt.Nb; ++ib)
			idDOF << _idDOF[idof] + ib*6;
	
	auto MultiplyAB = [&](UArray<UArray<VectorXd>> &A) {
		for (int idf = 0; idf < 6*dt.Nb; ++idf) {
			for (int jdf = 0; jdf < 6*dt.Nb; ++jdf) {
				for (int idof = 0; idof < idDOF.size(); ++idof) {
					if (( diag &&  idf == idDOF[idof] && jdf == idDOF[idof]) ||
					    (!diag && (idf == idDOF[idof] || jdf == idDOF[idof]))) {
						A[idf][jdf] *= factor;		
						break;
					}
				}
			}
		}
    };		
	if (a && IsLoadedA())
		MultiplyAB(dt.A);
	if (a && IsLoadedAinf_w())
		MultiplyAB(dt.Ainf_w);
	if (b && IsLoadedB())
		MultiplyAB(dt.B);
	
	if (a && IsLoadedA_P())
		MultiplyAB(dt.A_P);
	if (b && IsLoadedB_H())
		MultiplyAB(dt.B_H);
	if (b && IsLoadedB_P())
		MultiplyAB(dt.B_P);

	auto MultiplyAinfA0 = [&](MatrixXd &A) {
		for (int idf = 0; idf < 6*dt.Nb; ++idf) {
			for (int jdf = 0; jdf < 6*dt.Nb; ++jdf) {
				for (int idof = 0; idof < idDOF.size(); ++idof) {
					if (( diag &&  idf == idDOF[idof] && jdf == idDOF[idof]) ||
					    (!diag && (idf == idDOF[idof] || jdf == idDOF[idof]))) {
						A(idf, jdf) *= factor;		
						break;
					}
				}
			}
		}	
    };	
	if (a && IsLoadedAinf()) 
		MultiplyAinfA0(dt.Ainf);
	if (a && IsLoadedA0()) 
		MultiplyAinfA0(dt.A0);
		
	auto MultiplyF = [&](Forces &ex) {
		for (int ib = 0; ib < dt.Nb; ++ib) 
			for (int ih = 0; ih < dt.Nh; ++ih) 
				for (int ifr = 0; ifr < dt.Nf; ++ifr) 
					for (int idof = 0; idof < _idDOF.size(); ++idof) 
						ex[ib][ih](ifr, _idDOF[idof]) *= factor;
	};
	if (f && IsLoadedFex())
		MultiplyF(dt.ex);
	if (f && IsLoadedFsc())
		MultiplyF(dt.sc);
	if (f && IsLoadedFfk())
		MultiplyF(dt.fk);	
	if (f && IsLoadedRAO())
		MultiplyF(dt.rao);

	if (f && IsLoadedFsc_pot())
		MultiplyF(dt.sc_pot);
	if (f && IsLoadedFfk_pot())
		MultiplyF(dt.fk_pot);
	if (f && IsLoadedFfk_pot_bmr())
		MultiplyF(dt.fk_pot_bmr);

	auto MultiplyMD = [&]() {
		for (int ib = 0; ib < dt.Nb; ++ib) 
    		for (int ih = 0; ih < dt.mdhead.size(); ++ih) 
    			for (int idof = 0; idof < idDOF.size(); ++idof) 
	    			for (int ifr = 0; ifr < dt.Nf; ++ifr) 
	    				dt.md[ib][ih][idDOF[idof]](ifr) *= factor;
	};
	if (ismd && IsLoadedMD()) 
		MultiplyMD();
			
	auto MultiplySumDif = [&](UArray<UArray<UArray<MatrixXcd>>> &qtf) {
		for (int ib = 0; ib < dt.Nb; ++ib)
	        for (int ih = 0; ih < dt.qhead.size(); ++ih) 
				for (int idof = 0; idof < _idDOF.size(); ++idof) 
					qtf[ib][ih][_idDOF[idof]] *= factor;													
	};
	if (qtf && IsLoadedQTF(true)) 
		MultiplySumDif(dt.qtfsum);
	if (qtf && IsLoadedQTF(false))
		MultiplySumDif(dt.qtfdif);
	
	// Some previous data is now invalid
	dt.Kirf.Clear();
	
	String error = AfterLoad();
	if (!error.IsEmpty())
		throw Exc(Format(t_("Problem reseting DOF: '%s'\n%s"), error));
}

void Hydro::SwapDOF(int ib1, int ib2) {
	for (int idof = 0; idof < 6; ++idof)	
		SwapDOF(ib1, idof, ib2, idof);	
	
	Swap(dt.msh[ib1].dt.c0, dt.msh[ib2].dt.c0);
	Swap(dt.msh[ib1].dt.cg, dt.msh[ib2].dt.cg);
	Swap(dt.msh[ib1].dt.cb, dt.msh[ib2].dt.cb);
	Swap(dt.msh[ib1].dt.Vo, dt.msh[ib2].dt.Vo);
	Swap(dt.msh[ib1].dt.name, dt.msh[ib2].dt.name);
}

void Hydro::SwapDOF(int ib1, int idof1, int ib2, int idof2) {
	auto SwapAB = [&](UArray<UArray<VectorXd>> &A) {
		UArray<UArray<VectorXd>> An(6*dt.Nb);
		for (int i = 0; i < 6*dt.Nb; ++i) 
			An[i].SetCount(6*dt.Nb);
		
		for (int idof = 0; idof < 6*dt.Nb; ++idof) {
			for (int jdof = 0; jdof < 6*dt.Nb; ++jdof) {
				int idofn = idof, jdofn = jdof;
				if (idofn == idof1+6*ib1)
					idofn = idof2+6*ib2;
				else if (idofn == idof2+6*ib2)
					idofn = idof1+6*ib1;
				if (jdofn == idof1+6*ib1)
					jdofn = idof2+6*ib2;
				else if (jdofn == idof2+6*ib2)
					jdofn = idof1+6*ib1;	 
				An[idofn][jdofn] = pick(A[idof][jdof]);
			}
		}
    	A = pick(An);
    };		
	if (IsLoadedA())
		SwapAB(dt.A);
	if (IsLoadedAinf_w())
		SwapAB(dt.Ainf_w);
	if (IsLoadedB())
		SwapAB(dt.B);
	
	if (IsLoadedA_P())
		SwapAB(dt.A_P);
	if (IsLoadedB_H())
		SwapAB(dt.B_H);
	if (IsLoadedB_P())
		SwapAB(dt.B_P);
				
	auto SwapAinfA0 = [&](MatrixXd &A) {
		Swap(A, idof1+6*ib1, idof2+6*ib2);
    };	
	if (IsLoadedAinf()) 
		SwapAinfA0(dt.Ainf);
	if (IsLoadedA0()) 
		SwapAinfA0(dt.A0);
			  
	auto SwapF = [&](Forces &ex) {
		for (int ih = 0; ih < dt.Nh; ++ih) 
			for (int ifr = 0; ifr < dt.Nf; ++ifr)
				Swap(ex[ib2][ih](ifr, idof2), ex[ib1][ih](ifr, idof1));
	};
	if (IsLoadedFex())
		SwapF(dt.ex);
	if (IsLoadedFsc())
		SwapF(dt.sc);
	if (IsLoadedFfk())
		SwapF(dt.fk);	
	if (IsLoadedRAO())
		SwapF(dt.rao);

	if (IsLoadedFsc_pot())
		SwapF(dt.sc_pot);
	if (IsLoadedFfk_pot())
		SwapF(dt.fk_pot);
	if (IsLoadedFfk_pot_bmr())
		SwapF(dt.fk_pot_bmr);

	auto SwapMD = [&]() {
		for (int ih = 0; ih < dt.mdhead.size(); ++ih) 
			Swap(dt.md[ib1][ih][idof1], dt.md[ib2][ih][idof2]);
	};
	if (IsLoadedMD(true)) 
		SwapMD();
		
	auto SwapSumDif = [&](UArray<UArray<UArray<MatrixXcd>>> &qtf) {
        for (int ih = 0; ih < dt.qhead.size(); ++ih) 
			Swap(qtf[ib1][ih][idof1], qtf[ib2][ih][idof2]); 		
	};
	if (IsLoadedQTF(true)) 
		SwapSumDif(dt.qtfsum);
	if (IsLoadedQTF(false))
		SwapSumDif(dt.qtfdif);

	if (IsLoadedC()) 
		Swap(dt.msh[ib1].dt.C, dt.msh[ib2].dt.C, idof1, idof2);
	
	if (IsLoadedM()) 
		Swap(dt.msh[ib1].dt.M, dt.msh[ib2].dt.M, idof1, idof2);

	if (IsLoadedKirf()) // Kirf structure is like A and B
		SwapAB(dt.Kirf);
		
	String error = AfterLoad();
	if (!error.IsEmpty())
		throw Exc(Format(t_("Problem swaping DOF: '%s'\n%s"), error));	
}

void Hydro::DeleteFrequencies(const UVector<int> &idFreq) {
	if (idFreq.IsEmpty()) 
		return;
	
	auto DeleteAB = [&](UArray<UArray<VectorXd>> &A) {
        UArray<UArray<VectorXd>> An;
	
		An.SetCount(6*dt.Nb);
		for (int idof = 0; idof < 6*dt.Nb; ++idof) {
			An[idof].SetCount(6*dt.Nb);
			for (int jdof = 0; jdof < 6*dt.Nb; ++jdof) {
				An[idof][jdof].resize(dt.Nf - idFreq.size());	
				int i = 0, j = 0;
				for (int iif = 0; iif < dt.Nf; ++iif) {
					if (j >= idFreq.size() || iif != idFreq[j])
						An[idof][jdof][i++] = A[idof][jdof][iif];		
					else 
						j++;
				}
			}
		}
		A = pick(An);
    };
	
	if (IsLoadedA())
		DeleteAB(dt.A);
	if (IsLoadedAinf_w())
		DeleteAB(dt.Ainf_w);
	if (IsLoadedB())
		DeleteAB(dt.B);
	
	if (IsLoadedA_P())
		DeleteAB(dt.A_P);
	if (IsLoadedB_H())
		DeleteAB(dt.B_H);
	if (IsLoadedB_P())
		DeleteAB(dt.B_P);

	auto DeleteF = [&](Forces &ex) {
        Forces _ex;
	
		_ex.SetCount(dt.Nb);
		for (int ib = 0; ib < dt.Nb; ++ib) {
			_ex[ib].SetCount(dt.Nh);
		    for (int ih = 0; ih < dt.Nh; ++ih) {
		        _ex[ib][ih].resize(dt.Nf - idFreq.size(), 6);
		    	for (int idof = 0; idof < 6; ++idof) {
					int i = 0, j = 0;
					for (int iif = 0; iif < dt.Nf; ++iif) {
						if (j >= idFreq.size() || iif != idFreq[j]) 
							_ex[ib][ih](i++, idof) = ex[ib][ih](iif, idof);
						else 
							j++;
					}
		    	}
		    }
	    }
	    ex = pick(_ex);
    };	

	if (IsLoadedFex())
		DeleteF(dt.ex);
	if (IsLoadedFsc())
		DeleteF(dt.sc);
	if (IsLoadedFfk())
		DeleteF(dt.fk);	
	if (IsLoadedRAO())
		DeleteF(dt.rao);
	
	if (IsLoadedFsc_pot())
		DeleteF(dt.sc_pot);
	if (IsLoadedFfk_pot())
		DeleteF(dt.sc_pot);
	if (IsLoadedFfk_pot_bmr())
		DeleteF(dt.sc_pot);

	auto DeleteMD = [&]() {
		UArray<UArray<UArray<VectorXd>>> mdn;

		mdn.SetCount(dt.Nb);
		for (int ib = 0; ib < dt.Nb; ++ib) {
			mdn[ib].SetCount(int(dt.mdhead.size()));
    		for (int ih = 0; ih < dt.mdhead.size(); ++ih) {
    			mdn[ib][ih].SetCount(6);
    			for (int idf = 0; idf < 6; ++idf) {
    				mdn[ib][ih][idf].setConstant(dt.Nf - idFreq.size());
    				int i = 0, j = 0;
					for (int iif = 0; iif < dt.Nf; ++iif) {
						if (j >= idFreq.size() || iif != idFreq[j])
							dt.md[ib][ih][idf](i++) = dt.md[ib][ih][idf](iif);		
						else 
							j++;
					}
    			}
    		}
		}
	    dt.md = pick(mdn);			
    };

	if (IsLoadedMD())
		DeleteMD();
			
	int j = idFreq.size()-1;	
	for (int i = dt.w.size()-1; i >= 0 && j >= 0; --i) {
		if (i == idFreq[j]) {	
			dt.w.Remove(i);
			j--;
		}
	}
	dt.Nf = dt.w.size();
}

void Hydro::DeleteFrequenciesQTF(const UVector<int> &idFreqQTF) {
	if (idFreqQTF.size() > 0) {
		UVector<int> vids;
		LinSpaced(vids, int(dt.qw.size()), 0, int(dt.qw.size())-1);
		for (int i = idFreqQTF.size()-1; i >= 0; --i) 
			vids.Remove(idFreqQTF[i]);
		VectorXi ids;
		::Copy(vids, ids);
		dt.qw = VectorXd(dt.qw(ids));
		
		auto DeleteSumDif = [&](UArray<UArray<UArray<MatrixXcd>>> &qtf) {
			for (int ib = 0; ib < dt.Nb; ++ib)
		        for (int ih = 0; ih < dt.qhead.size(); ++ih) 
					for (int idf = 0; idf < 6; ++idf) {
						MatrixXcd &m = qtf[ib][ih][idf];
						m = MatrixXcd(m(indexing::all, ids));
						m = MatrixXcd(m(ids, indexing::all));
					}
		};
		if (IsLoadedQTF(true)) 
			DeleteSumDif(dt.qtfsum);
		if (IsLoadedQTF(false))
			DeleteSumDif(dt.qtfdif);
	}
}

void Hydro::DeleteHeadings(const UVector<int> &idHead) {
	if (idHead.size() > 0) {
		auto DeleteF = [&](Forces &ex) {
			int j = idHead.size()-1;	
			for (int i = dt.head.size()-1; i >= 0 && j >= 0; --i) {
				if (i == idHead[j]) {	
					ex.Remove(i);
					j--;
				}
			}
	    };	
	
		if (IsLoadedFex())
			DeleteF(dt.ex);
		if (IsLoadedFsc())
			DeleteF(dt.sc);
		if (IsLoadedFfk())
			DeleteF(dt.fk);	
		if (IsLoadedRAO())
			DeleteF(dt.rao);
	
		if (IsLoadedFsc_pot())
			DeleteF(dt.sc_pot);	
		if (IsLoadedFfk_pot())
			DeleteF(dt.fk_pot);	
		if (IsLoadedFfk_pot_bmr())
			DeleteF(dt.fk_pot_bmr);	
		
		int j = idHead.size()-1;	
		for (int i = dt.head.size()-1; i >= 0 && j >= 0; --i) {
			if (i == idHead[j]) {	
				dt.head.Remove(i);
				j--;
			}
		}
		dt.Nh = dt.head.size();
	}
}

void Hydro::DeleteHeadingsMD(const UVector<int> &idHead) {
	if (idHead.size() > 0) {
		auto DeleteMD = [&]() {
			for (int ib = 0; ib < dt.Nb; ++ib) {
	    		int j = idHead.size()-1;	
				for (int i = int(dt.mdhead.size())-1; i >= 0 && j >= 0; --i) {
	    			if (i == idHead[j]) {	
						dt.md[ib].Remove(i);
						j--;
					}
				}
			}
	    };

		if (IsLoadedMD())
			DeleteMD();
		
		UArray<std::complex<double>> mdh;
		::Copy(dt.mdhead, mdh);
		int j = idHead.size()-1;	
		for (int i = mdh.size()-1; i >= 0 && j >= 0; --i) {
			if (i == idHead[j]) {	
				mdh.Remove(i);
				j--;
			}
		}
		::Copy(mdh, dt.mdhead);
	}
}

void Hydro::DeleteHeadingsQTF(const UVector<int> &idHeadQTF) {
	if (idHeadQTF.size() > 0) {
		UVector<int> vids;
		LinSpaced(vids, int(dt.qhead.size()), 0, int(dt.qhead.size())-1);
		for (int i = idHeadQTF.size()-1; i >= 0; --i) 
			vids.Remove(idHeadQTF[i]);
		VectorXi ids;
		::Copy(vids, ids);
		dt.qhead = VectorXcd(dt.qhead(ids));
			
		auto DeleteSumDif = [&](UArray<UArray<UArray<MatrixXcd>>> &qtf) {
			for (int ib = 0; ib < dt.Nb; ++ib)
		        for (int ih = idHeadQTF.size()-1; ih >= 0; --ih)
		        	qtf[ib].Remove(idHeadQTF[ih]);
		};
		if (IsLoadedQTF(true)) 
			DeleteSumDif(dt.qtfsum);
		if (IsLoadedQTF(false))
			DeleteSumDif(dt.qtfdif);
	}
}

void Hydro::FillFrequencyGapsABForces(bool zero, int maxFreq) {
	if (dt.w.size() == 0)
		return;

	VectorXd w_, nw;
	::Copy(dt.w, w_);
	
	UVector<int> idsx, w0x;
	GapFillingAxisParams(w_, maxFreq, idsx, w0x, nw);
	
	auto FillAB = [&](UArray<UArray<VectorXd>> &A) {
		for (int idof = 0; idof < 6*dt.Nb; ++idof) {
			for (int jdof = 0; jdof < 6*dt.Nb; ++jdof) {
				VectorXd nm;
				const VectorXd &m = A[idof][jdof];
				GapFilling(w_, m, idsx, w0x, nw, nm, zero, maxFreq);					
				A[idof][jdof] = pick(nm);
			}
		}
    };
		
	if (IsLoadedA())
		FillAB(dt.A);
	if (IsLoadedAinf_w())
		FillAB(dt.Ainf_w);
	if (IsLoadedB())
		FillAB(dt.B);
	
	if (IsLoadedA_P())
		FillAB(dt.A_P);
	if (IsLoadedB_H())
		FillAB(dt.B_H);
	if (IsLoadedB_P())
		FillAB(dt.B_P);
	
	auto FillF = [&](Forces &ex) {
		for (int ib = 0; ib < dt.Nb; ++ib) {
		    for (int ih = 0; ih < dt.Nh; ++ih) {
		        MatrixXcd nmn(nw.size(), 6);
		    	for (int idof = 0; idof < 6; ++idof) {
		    		VectorXcd nm;
		    		const VectorXcd &m = ex[ib][ih].col(idof);
		    		GapFilling(w_, m, idsx, w0x, nw, nm, zero, maxFreq);					
					nmn.col(idof) = nm;
		    	}
		    	ex[ib][ih] = pick(nmn);
		    }
	    }
    };	

	if (IsLoadedFex())
		FillF(dt.ex);
	if (IsLoadedFsc())
		FillF(dt.sc);
	if (IsLoadedFfk())
		FillF(dt.fk);	
	if (IsLoadedRAO())
		FillF(dt.rao);	
	
	if (IsLoadedFsc_pot())
		FillF(dt.sc_pot);
	if (IsLoadedFfk_pot())
		FillF(dt.fk_pot);	
	if (IsLoadedFfk_pot_bmr())
		FillF(dt.fk_pot_bmr);	
	
	auto FillMD = [&]() {
		for (int ib = 0; ib < dt.Nb; ++ib) {
    		for (int ih = 0; ih < dt.mdhead.size(); ++ih) {
    			for (int idf = 0; idf < 6; ++idf) {
    				VectorXd nm(nw.size());
	    			VectorXd &m = dt.md[ib][ih][idf];
					GapFilling(w_, m, idsx, w0x, nw, nm, zero, maxFreq);
					dt.md[ib][ih][idf] = pick(nm);
    			}
    		}
		}
    };		
	
	if (IsLoadedMD())
		FillMD();		
	
	dt.Nf = int(nw.size());
	::Copy(nw, dt.w);
}

void Hydro::FillFrequencyGapsQTF(bool zero, int maxFreq) {
	if (dt.qw.size() == 0)
		return;
	
	VectorXd nw;

	UVector<int> idsx, w0x;
	GapFillingAxisParams(dt.qw, maxFreq, idsx, w0x, nw);
	
	auto FillSumDif = [&](UArray<UArray<UArray<MatrixXcd>>> &qtf) {
		for (int ib = 0; ib < dt.Nb; ++ib) {
	        for (int ih = 0; ih < dt.qhead.size(); ++ih) {
				for (int idof = 0; idof < 6; ++idof) {
					MatrixXcd nm, &m = qtf[ib][ih][idof];
					GapFilling(dt.qw, dt.qw, m, idsx, w0x, idsx, w0x, nw, nw, nm, zero, maxFreq);					
					m = pick(nm);
				}
	        }
		}
	};

	if (IsLoadedQTF(true)) 
		FillSumDif(dt.qtfsum);
	if (IsLoadedQTF(false)) 
		FillSumDif(dt.qtfdif);
	
	dt.qw = pick(nw);
}

void Hydro::FillFrequencyGapsABForcesZero() {
	if (dt.w.size() == 0)
		return;

	auto FillAB = [&](UArray<UArray<VectorXd>> &A) {
		for (int idof = 0; idof < 6*dt.Nb; ++idof) {
			for (int jdof = 0; jdof < 6*dt.Nb; ++jdof) {
				VectorXd &a = A[idof][jdof];
				if (a.size() == 0 || !IsNum(a(0)))
					a = VectorXd::Zero(dt.Nf);
			}
		}
    };
		
	if (IsLoadedA())
		FillAB(dt.A);
	if (IsLoadedAinf_w())
		FillAB(dt.Ainf_w);
	if (IsLoadedB())
		FillAB(dt.B);

	if (IsLoadedA_P())
		FillAB(dt.A_P);
	if (IsLoadedB_H())
		FillAB(dt.B_H);
	if (IsLoadedB_P())
		FillAB(dt.B_P);
	
	auto FillA = [&](MatrixXd &A) {
		if (A.size() == 0)
			A = MatrixXd::Zero(6*dt.Nb, 6*dt.Nb);
		else 
			A = A.unaryExpr([](double x){return IsNum(x) ? x : 0;});		// Replace NaN with 0
    };
		
	if (IsLoadedAinf())
		FillA(dt.Ainf);
	if (IsLoadedA0())
		FillA(dt.A0);
	//if (IsLoadedDlin())
		//FillA(Dlin);
	
	auto FillF = [&](Forces &ex) {
		for (int ib = 0; ib < dt.Nb; ++ib) {
		    for (int ih = 0; ih < dt.Nh; ++ih) {
		        MatrixXcd nmn(dt.Nf, 6);
		    	for (int idof = 0; idof < 6; ++idof) {
		    		const VectorXcd &m = ex[ib][ih].col(idof);
		    		if (!IsNum(m(0))) 
		    			nmn.col(idof) = VectorXcd::Zero(dt.Nf);
		    		else 
		    			nmn.col(idof) = m;
		    	}
		    	ex[ib][ih] = pick(nmn);
		    }
	    }
    };	

	if (IsLoadedFex())
		FillF(dt.ex);
	if (IsLoadedFsc())
		FillF(dt.sc);
	if (IsLoadedFfk())
		FillF(dt.fk);	
	if (IsLoadedRAO())
		FillF(dt.rao);	
	
	if (IsLoadedFsc_pot())
		FillF(dt.sc_pot);
	if (IsLoadedFfk_pot())
		FillF(dt.fk_pot);
	if (IsLoadedFfk_pot_bmr())
		FillF(dt.fk_pot_bmr);
	
	auto FillMD = [&]() {
		for (int ib = 0; ib < dt.Nb; ++ib) {
    		for (int ih = 0; ih < dt.mdhead.size(); ++ih) {
    			for (int idf = 0; idf < 6; ++idf) {
    				if (dt.md[ib][ih][idf].size() == 0 || !IsNum(dt.md[ib][ih][idf][0]))
						dt.md[ib][ih][idf] = VectorXd::Zero(dt.Nf);
    			}
    		}
		}
    };		
	
	if (IsLoadedMD())
		FillMD();			
}

void Hydro::FillFrequencyGapsQTFZero() {
	if (dt.qw.size() == 0)
		return;
	
	Eigen::Index nf = dt.qw.size();
	
	auto FillSumDif = [&](UArray<UArray<UArray<MatrixXcd>>> &qtf) {
		for (int ib = 0; ib < dt.Nb; ++ib) {
	        for (int ih = 0; ih < dt.qhead.size(); ++ih) {
				for (int idof = 0; idof < 6; ++idof) {
					MatrixXcd &m = qtf[ib][ih][idof];
					if (m.size() == 0 || !IsNum(m(0,0)))
						m = MatrixXcd::Zero(nf, nf);
				}
	        }
		}
	};

	if (IsLoadedQTF(true)) 
		FillSumDif(dt.qtfsum);
	if (IsLoadedQTF(false)) 
		FillSumDif(dt.qtfdif);
}

void Hydro::CopyQTF_MD() {
	dt.mdtype = 9;
	::Copy(dt.qhead, dt.mdhead);
	
	Initialize_MD(dt.md, dt.Nb, int(dt.mdhead.size()), dt.Nf);

	VectorXd ww;
	::Copy(dt.w, ww);
	
	for (int ib = 0; ib < dt.Nb; ++ib) {
        for (int ih = 0; ih < dt.qhead.size(); ++ih) {
			for (int idof = 0; idof < 6; ++idof) {
				const MatrixXcd &m = dt.qtfdif[ib][ih][idof];
				VectorXd diag = Eigen::abs(m.diagonal().array());
				ResampleY(dt.qw, diag, ww, dt.md[ib][ih][idof]);
            }
        }
	}
}

VectorXd AvgSafe(const VectorXd &a, const VectorXd &b) {
	ASSERT(a.size() == b.size());
	VectorXd r(a.size());
	for (int i = 0; i < a.size(); ++i) 
		r[i] = AvgSafe(a[i], b[i]);
	return r;
}

// Forces the symmetry in values that have to be symmetric
void Hydro::Symmetrize() {
	auto SymmetrizeAB = [&](UArray<UArray<VectorXd>> &A) {
		for (int idf = 0; idf < 6*dt.Nb; ++idf) 
			for (int jdf = idf+1; jdf < 6*dt.Nb; ++jdf) 
				A[idf][jdf] = A[jdf][idf] = AvgSafe(A[idf][jdf], A[jdf][idf]);
    };		
	if (IsLoadedA())
		SymmetrizeAB(dt.A);
	if (IsLoadedAinf_w())
		SymmetrizeAB(dt.Ainf_w);
	if (IsLoadedB())
		SymmetrizeAB(dt.B);

	if (IsLoadedA_P())
		SymmetrizeAB(dt.A_P);
	if (IsLoadedB_H())
		SymmetrizeAB(dt.B_H);
	if (IsLoadedB_P())
		SymmetrizeAB(dt.B_P);


	auto SymmetrizeAinfA0 = [&](MatrixXd &A) {
		for (int idf = 0; idf < 6*dt.Nb; ++idf) 
			for (int jdf = idf+1; jdf < 6*dt.Nb; ++jdf) 
				A(idf, jdf) = A(jdf, idf) = AvgSafe(A(idf, jdf), A(jdf, idf));
    };	
	if (IsLoadedAinf()) 
		SymmetrizeAinfA0(dt.Ainf);
	if (IsLoadedA0()) 
		SymmetrizeAinfA0(dt.A0);
		
	auto SymmetrizeSumDif = [&](UArray<UArray<UArray<MatrixXcd>>> &qtf, bool isSum) {
		for (int ib = 0; ib < dt.Nb; ++ib) {
	        for (int ih = 0; ih < dt.qhead.size(); ++ih) {
				for (int idf = 0; idf < 6; ++idf) { 
					MatrixXcd &c = qtf[ib][ih][idf];
					Eigen::Index rows = c.rows();
					for (int iw = 0; iw < rows; ++iw) {
						for (int jw = iw+1; jw < rows; ++jw) {
							std::complex<double> &cij = c(iw, jw), &cji = c(jw, iw);
							if (isSum)
								cij = cji = AvgSafe(cij, cji);
							else {
								std::complex<double> cji_ = std::complex<double>(cji.real(), -cji.imag());
								cij = cji_ = AvgSafe(cij, cji_);
								cji = std::complex<double>(cji_.real(), -cji_.imag());
							}
						}
					}
				}
	        }
		}
	};
	if (IsLoadedQTF(true)) 
		SymmetrizeSumDif(dt.qtfsum, true);
	if (IsLoadedQTF(false))
		SymmetrizeSumDif(dt.qtfdif, false);
	
	String error = AfterLoad();
	if (!error.IsEmpty())
		throw Exc(Format(t_("Problem symmetrizing data: '%s'\n%s"), error));	
	
}

double Hydro::GetQTFVal(int ib, int idof, int idh, int ifr1, int ifr2, bool isSum, char what, bool getDim) const {
	const UArray<UArray<UArray<MatrixXcd>>> &qtf = isSum ? dt.qtfsum : dt.qtfdif;
	if (qtf.IsEmpty())
		return Null;
	
	const MatrixXcd &m = qtf[ib][idh][idof];
	
	if (IsNull(m(ifr1, ifr2)))
		return Null;
	
	switch (what) {
	case 'm':	return F_(!getDim, abs(m(ifr1, ifr2)), idof);
	case 'p':	return arg(m(ifr1, ifr2));	 
	case 'r':	return F_(!getDim, m(ifr1, ifr2).real(), idof);
	case 'i':	return F_(!getDim, m(ifr1, ifr2).imag(), idof);
	}
	NEVER();	
	return Null;
}

MatrixXd Hydro::GetQTFMat(int ib, int idof, int idh, bool isSum, char what, bool getDim) const {
	MatrixXd ret;
	
	const UArray<UArray<UArray<MatrixXcd>>> &qtf = isSum ? dt.qtfsum : dt.qtfdif;
	if (qtf.IsEmpty())
		return ret;
	
	const MatrixXcd &m = qtf[ib][idh][idof];
	
	if (m.size() == 0)
		return ret;
	
	ret.resize(m.rows(), m.cols());

	switch (what) {
	case 'm':	for (int i = 0; i < m.size(); ++i)	ret(i) = F_(!getDim, abs(m(i)), idof);	break;
	case 'p':	for (int i = 0; i < m.size(); ++i)	ret(i) = arg(m(i));						break;
	case 'r':	for (int i = 0; i < m.size(); ++i)	ret(i) = F_(!getDim, m(i).real(), idof);break;
	case 'i':	for (int i = 0; i < m.size(); ++i)	ret(i) = F_(!getDim, m(i).imag(), idof);break;
	default: NEVER();
	}
	return ret;
}	

int Hydro::LoadHydro(UArray<Hydro> &hydros, String file, Function <bool(String, int)> Status) {
	String ext = ToLower(GetFileExt(file));
	String ret;
	
	int num = 1;
	
	Status(t_("Loading BEM file"), -1);
	
	if (ext == ".nc")
		ret = CapyNC_Load(file, hydros, num);
	else {
		Hydro &hyd = hydros.Add();	
	
		if (ext == ".cal" || ext == ".tec" || ext == ".inf") 
			ret = static_cast<Nemoh&>(hyd).Load(file, Status);
		else if (ext == ".out" || ext == ".hdf" || ext == ".mcn") 
			ret = static_cast<Wamit&>(hyd).Load(file, Status);
		else if (ext == ".in") 
			ret = static_cast<Hams&>(hyd).Load(file, Status);
		else if (ext == ".dat" || ext == ".fst") 
			ret = static_cast<Fast&>(hyd).Load(file, Status);
		else if (ext == ".1" || ext == ".2" || ext == ".3" || ext == ".3sc" || ext == ".3fk" || ext == ".7" || ext == ".8" || ext == ".9" ||
				   ext == ".hst" || ext == ".4" || ext == ".12s" || ext == ".12d" || ext == ".frc" || ext == ".pot" || ext == ".mmx") 
			ret = static_cast<Wamit&>(hyd).Load(file, Status);
		else if (ext == ".ah1" || ext == ".lis" || ext == ".qtf") 
			ret = static_cast<Aqwa&>(hyd).Load(file, Status);
		else if (ext == ".hdb") 
			ret = static_cast<Diodore&>(hyd).Load(file, Status);
		else if (ext == ".yml")
			ret = static_cast<OrcaWave&>(hyd).Load(file, Status);
	#ifdef PLATFORM_WIN32	
		else if (ext == ".owr") 
			ret = static_cast<OrcaWave&>(hyd).Load(file, Status);
	#endif
		else if (ext == ".mat") 
			ret = static_cast<Foamm&>(hyd).Load(file);
		else if (ext == ".bem") 
			ret = static_cast<Hydro&>(hyd).LoadSerialization(file);
		else if (ext == ".h5") 
			ret = static_cast<BemioH5&>(hyd).Load(file, Status);
		else if (ext == ".owd") 
			ret = t_("OrcaWAVE .owd binary format is not supported.\nHowever OrcaFLEX .yml is supported.\nTo get it, load the .owd file in OrcaFlex and save it as .yml");
		else 
			ret = Format(t_("Unknown BEM file extension in '%s'"), file);	
	}
	if (!ret.IsEmpty()) {
		hydros.SetCount(hydros.size() - num);
		throw Exc(ret);//Format(t_("Problem loading '%s'\n%s"), file, error));	
	}
	
	for (int i = hydros.size() - num; i < hydros.size(); ++i) {
		Hydro &hyd = hydros[i];

		ret = hyd.AfterLoad(Status);
		if (!ret.IsEmpty()) {
			//String error = RemoveAccents(ret);
			hydros.SetCount(hydros.size() - num);
			throw Exc(Format(t_("Problem processing '%s'\n%s"), file, ret));	
		}
		
		hyd.IncrementIdCount();
	}
	
	return num;
}

void Heal();
void Load(const VectorXd &w, const VectorXd &A, const VectorXd &B, double maxT, int num);
void Save(const VectorXd &w, VectorXd &A, VectorXd &Ainfw, double &ainf, VectorXd &B, 
			VectorXd &Tirf, VectorXd &Kinf);
			   				




	