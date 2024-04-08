// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include <STEM4U/Integral.h>
#include <STEM4U/Utility.h>
#include <STEM4U/SeaWaves.h>
#include "functions.h"
#include "heal.h"
#include <Xlnt/Xlnt.h>

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
	return GetK_IRF_MaxT(w);
}

void Hydro::GetK_IRF(double maxT, int numT) {
	if (Nf == 0 || B.IsEmpty())
		return;
	
    Kirf.SetCount(Nb*6); 			
    for (int i = 0; i < Nb*6; ++i) {
    	Kirf[i].SetCount(Nb*6); 			 
   		for (int j = 0; j < Nb*6; ++j)
			Kirf[i][j].setConstant(numT, NaNDouble);
    }
		
	GetTirf(Tirf, numT, maxT);
	
	UVector<double> y(Nf);
  	for (int idf = 0; idf < Nb*6; ++idf) {
    	for (int jdf = 0; jdf < Nb*6; ++jdf) { 
			if (B[idf][jdf].size() == 0 || !IsNum(B[idf][jdf][0])) 
				continue;
			if (dimen)
				GetKirf(Kirf[idf][jdf], Tirf, Get_w(), B[idf][jdf]);
			else {
				GetKirf(Kirf[idf][jdf], Tirf, Get_w(), B_dim(idf, jdf));
				Kirf[idf][jdf] /= g_rho_dim();
			}
    	}
  	}
}  

void Hydro::GetAinf() {
	if (Nf == 0 || A.size() < Nb*6 || !IsLoadedKirf())
		return;	
	
	Ainf.setConstant(Nb*6, Nb*6, NaNDouble);
	
    for (int i = 0; i < Nb*6; ++i) 
        for (int j = 0; j < Nb*6; ++j) 
            if (IsNum(Kirf[i][j][0]))
		    	Ainf(i, j) = ::GetAinf(Kirf[i][j], Tirf, Get_w(), A[i][j]);
}

void Hydro::GetRAO() {
	if (Nf == 0 || A.size() < Nb*6 || B.size() < Nb*6)
		throw Exc(t_("Insufficient data to get RAO: A and B are required"));	

	for (int ib = 0; ib < Nb; ++ib) 
		if (C.size() < ib+1 || C[ib].rows() < 6 || C[ib].cols() < 6 || 
			  M.size() < ib+1 || M[ib].rows() < 6 || M[ib].cols() < 6) 
			throw Exc(t_("Insufficient data to get RAO: C and M are required"));   
			      
	Initialize_Forces(rao);

	MatrixXd D = MatrixXd::Zero(6, 6);		// Unused for now
	MatrixXd D2 = MatrixXd::Zero(6, 6);
	
	for (int ib = 0; ib < Nb; ++ib) {
		MatrixXd C = C_(false, ib);
		const MatrixXd &M_ = M[ib];
		for (int ih = 0; ih < Nh; ++ih) {	
			for (int ifr = 0; ifr < Nf; ++ifr) {
				VectorXcd RAO = GetRAO(w[ifr], A_mat(false, ifr, ib, ib), B_mat(false, ifr, ib, ib), 
								F_(false, ex, ih, ifr), C, M_, D, D2);
				for (int idf = 0; idf < 6; ++idf)
					rao.force[ih](ifr, idf+6*ib) = RAO[idf];
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
	
	return Fwh.transpose()*M.inverse();
}
	
void Hydro::InitAinf_w() {
	Ainf_w.SetCount(Nb*6); 			
    for (int i = 0; i < Nb*6; ++i) {
    	Ainf_w[i].SetCount(Nb*6); 			 
   		for (int j = 0; j < Nb*6; ++j)
			Ainf_w[i][j].setConstant(Nf, NaNDouble);
    }
}

void Hydro::GetAinf_w() {
	if (Nf == 0 || A.size() < Nb*6 || !IsLoadedKirf())
		return;	
	
	InitAinf_w();
    
    for (int idf = 0; idf < Nb*6; ++idf)
        for (int jdf = 0; jdf < Nb*6; ++jdf) {
            if (!IsLoadedB(idf, jdf)) 
                continue;
            if (dimen)
		    	::GetAinf_w(Ainf_w[idf][jdf], Kirf[idf][jdf], Tirf, Get_w(), A[idf][jdf]);
            else {
                ::GetAinf_w(Ainf_w[idf][jdf], Kirf[idf][jdf]*g_rho_dim(), Tirf, Get_w(), A_dim(idf, jdf));
                Ainf_w[idf][jdf] *= (1/(rho_dim()*pow(len, GetK_AB(idf, jdf))));
            }
        }
}

void Hydro::GetB_H(int &num) {
	if (!IsLoadedFex())
		throw Exc(t_("The excitation force is not loaded"));
	if (Nf == 0)
		throw Exc(t_("No frecuencies loaded"));
	if (num > Nh)
		throw Exc(t_("Number of headings is higher than available"));
	if (num <= 0)
		throw Exc(t_("Not enough headings for Haskind calculation"));
	if (IsNull(h))
		throw Exc(t_("Unknown depth"));

	UVector<double> head360 = clone(head);

	double angle = 360/num;
	double nextAngle = head360[0] + angle;
	UVector<int> idRemove;
	for (int i = 1; i < head360.size(); ++i) {
		if (head360[i] < nextAngle) 
			idRemove << i;
		else
			nextAngle += angle;
	}
	num = head.size() - idRemove.size();
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

	VectorXd val(Nf);
	for (int ifr = 0; ifr < Nf; ++ifr)
	 	val[ifr] = w[ifr]*SeaWaves::WaveNumber_w(w[ifr], -1, g_dim(), true)/(4*M_PI*rho_dim()*g_dim()*g_dim());
	
	Initialize_AB(B_H);
	
    for (int idof = 0; idof < Nb*6; ++idof) {
        if (!IsLoadedFex(idof)) 		
            continue;
		
		VectorXd b(Nf);
		for (int ifr = 0; ifr < Nf; ++ifr) {
			UVector<double> F2(Nh);
			for (int ih = 0; ih < Nh; ++ih) 
				F2[ih] = sqr(F_dim(abs(ex.force[ih](ifr, idof)), idof));
			
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
		if (dimen)
			B_H[idof][idof] = b*rho_ndim()/rho_dim();
		else {
			for (int ifr = 0; ifr < Nf; ++ifr)
				b(ifr) /= (rho_dim()*pow(len, GetK_AB(idof, idof))*w[ifr]);
			B_H[idof][idof] = b;
		}
    }
}

void Hydro::GetOgilvieCompliance(bool zremoval, bool thinremoval, bool decayingTail, UVector<int> &vidof, UVector<int> &vjdof) {
	vidof.Clear();
	vjdof.Clear();
	if (Nf == 0 || A.size() < Nb*6)
		return;	
	
	HealBEM data;
	
	if (Ainf.size() == 0) 
		Ainf.setConstant(6*Nb, 6*Nb, NaNDouble); 
	
	if (Ainf_w.size() == 0) {
		Ainf_w.SetCount(Nb*6); 			
	    for (int idf = 0; idf < Nb*6; ++idf) {
    		Ainf_w[idf].SetCount(Nb*6); 			 
   			for (int jdf = 0; jdf < Nb*6; ++jdf)
				Ainf_w[idf][jdf].setConstant(Nf, NaNDouble);
	    }
    }
    double maxT = min(bem->maxTimeA, Hydro::GetK_IRF_MaxT(w));
    int numT = bem->numValsA;
    
    if (Kirf.size() == 0) {
        Kirf.SetCount(Nb*6); 			
	    for (int idf = 0; idf < Nb*6; ++idf) {
	    	Kirf[idf].SetCount(Nb*6); 			 
	   		for (int jdf = 0; jdf < Nb*6; ++jdf)
				Kirf[idf][jdf].setConstant(numT, NaNDouble);
	    }
	}
		
    for (int idf = 0; idf < Nb*6; ++idf) {
        MatrixXd ex_hf(Nh, Nf);
        
        for (int jdf = 0; jdf < Nb*6; ++jdf) {
            if (B[idf][jdf].size() == 0 || !IsNum(B[idf][jdf][0])) 
                ;
            else {
                bool done;
	    		if (data.Load(Get_w(), A_dim(idf, jdf), Ainf_dim(idf, jdf), B_dim(idf, jdf), numT, maxT, ex_hf) &&
					data.Heal(zremoval, thinremoval, decayingTail, done)) {
	            	data.Save(Get_w(), A[idf][jdf], Ainf_w[idf][jdf], Ainf(idf, jdf), B[idf][jdf], Tirf, Kirf[idf][jdf]); 
	            	if (done) {
		            	vidof << idf;
		            	vjdof << jdf;
	            	}
	    		} else
	    			data.Reset(Get_w(), A[idf][jdf], Ainf_w[idf][jdf], Ainf(idf, jdf), B[idf][jdf], Tirf, Kirf[idf][jdf]);
	    		if (dimen) {
	    			dimen = false;
	    			A[idf][jdf] 	= pick(A_ndim(idf, jdf));
	    			Ainf_w[idf][jdf] *= (rho_ndim()/rho_dim());
	    			Ainf(idf, jdf)   *= (rho_ndim()/rho_dim());
	    			B[idf][jdf] 	= pick(B_ndim(idf, jdf));
	    			Kirf[idf][jdf]  = pick(Kirf_ndim(idf, jdf));
	    			dimen = true;
	    		} else {
	    			dimen = true;
	    			A[idf][jdf] 	= pick(A_ndim(idf, jdf));
	    			Ainf_w[idf][jdf] *= (1/(rho_ndim()*pow(len, GetK_AB(idf, jdf))));
	    			Ainf(idf, jdf)   *= (1/(rho_ndim()*pow(len, GetK_AB(idf, jdf))));
	    			B[idf][jdf] 	= pick(B_ndim(idf, jdf));
	    			Kirf[idf][jdf] 	= pick(Kirf_ndim(idf, jdf));
	    			dimen = false;
	    		}
            }
        }
    }
    rao.Clear();	// Previous RAO is now invalid
}

void Hydro::GetWaveTo(double xto, double yto) {
	double dx = xto - x_w,
		   dy = yto - y_w;
	
	for (int ib = 0; ib < Nb; ++ib) 
		AddWave(ib, dx, dy);
	
	x_w = xto;
	y_w = yto;
}

String Hydro::SpreadNegative() {
	String ret;
	UVector<String> errors;
	for (int ib = 0; ib < Nb; ++ib) {
		for (int idp = 0; idp < pots[ib].size(); ++idp) {
			UVector<int> panIDs;
			for (int ifr = 0; ifr < Nf; ++ifr) {
				for (int idf = 0; idf < 6; ++idf) {		// Only diagonal
					if (Apan(ib, idp, idf, idf, ifr) < 0) {
						if (panIDs.IsEmpty())
							meshes[ib].mesh.GetClosestPanels(idp, panIDs);
						for (int i = 0; i < panIDs.size(); ++i) {
							if (Apan(ib, i, idf, idf, ifr) > 0) {
								double delta = Apan(ib, i, idf, idf, ifr) + Apan(ib, idp, idf, idf, ifr);
								if (delta >= 0) {
									Apan(ib, i, idf, idf, ifr) += Apan(ib, idp, idf, idf, ifr);
									Apan(ib, idp, idf, idf, ifr) = 0;
									break;
								} else {
									Apan(ib, idp, idf, idf, ifr) += Apan(ib, i, idf, idf, ifr);
									Apan(ib, i, idf, idf, ifr) = 0;
								}
							}
						}
						if (Apan(ib, idp, idf, idf, ifr) < 0) {
							FindAdd(errors, Format(t_("%d.%s Freq %d (%f rad/s)"), idf+1, BEM::strDOFtext[idf], ifr, w[ifr]));
							// Resets this dof and frequency for all panels
							for (int i = 0; i < pots[ib].size(); ++i) 
								Apan(ib, i, idf, idf, ifr) = 0;
						}
					}
				}
			}
		}
	}
	Sort(errors);
	for (const String s : errors) {
		if (!ret.IsEmpty())
			ret << "\n";
		ret << s;
	}
	return ret;
}

void Hydro::MapNodes(int ib, UVector<Point3D> &points, Tensor<double, 4> &Apan_nodes, Tensor<double, 4> &Bpan_nodes) const {
	Apan_nodes.resize(points.size(), 6, 6, Nf);
	Apan_nodes.setZero();
	Bpan_nodes.resize(points.size(), 6, 6, Nf);
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
	
	for (int idp = 0; idp < pots[ib].size(); ++idp) {
		const Point3D &p = meshes[ib].mesh.panels[idp].centroidPaint;
		int ip = GetClosest(p, points);
		for (int ifr = 0; ifr < Nf; ++ifr) {
			for (int idf1 = 0; idf1 < 6; ++idf1) {		
				for (int idf2 = 0; idf2 < 6; ++idf2) {
					Apan_nodes(ip, idf1, idf2, ifr) += Apan(ib, idp, idf1, idf2, ifr);
					Bpan_nodes(ip, idf1, idf2, ifr) += Bpan(ib, idp, idf1, idf2, ifr);
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
			for (int c = 0; c < 6; ++c) 
				g.Set(Null, col++, Format(t_("A_%s_%s"), BEM::StrDOF(c), BEM::StrDOF(c)));		g.AddCol(60);
			for (int c = 0; c < 6; ++c) 
				g.Set(Null, col++, Format(t_("B_%s_%s"), BEM::StrDOF(c), BEM::StrDOF(c)));		g.AddCol(60);
		} else {
			for (int r = 0; r < 6; ++r) 
				for (int c = 0; c < 6; ++c) 
					g.Set(Null, col++, Format(t_("A_%s_%s"), BEM::StrDOF(r), BEM::StrDOF(c)));	g.AddCol(60);
			for (int r = 0; r < 6; ++r) 
				for (int c = 0; c < 6; ++c) 
					g.Set(Null, col++, Format(t_("B_%s_%s"), BEM::StrDOF(r), BEM::StrDOF(c)));	g.AddCol(60);
		}
	}
			
	for (int row = 0; row < ids.size(); ++row) {
		int col = 0;
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
				for (int row = 0; row < ids.size(); ++row) {
					for (int r = 0; r < 6; ++r) 
						for (int c = 0; c < 6; ++c)
							g.Set(Null, col++, Apan(row, r, c, ifr));	
					for (int r = 0; r < 6; ++r) 
						for (int c = 0; c < 6; ++c)
							g.Set(Null, col++, Bpan(row, r, c, ifr));		
				}
			}
		}
	}
}
		
void Hydro::SaveMap(String fileName, String type, int ifr, bool onlyDiagonal, const UVector<int> &ids, const UVector<Point3D> &points, 
		const Tensor<double, 4> &Apan, const Tensor<double, 4> &Bpan) const {
	if (IsNull(ifr)) {
		UArray<Grid> grid(Nf);
		for (int ifr = 0; ifr < Nf; ++ifr) 
			SaveMap(grid[ifr], ifr, onlyDiagonal, ids, points, Apan, Bpan);
		
		if (type == ".csv") {
			for (int ifr = 0; ifr < Nf; ++ifr) { 
				String folder = GetFileFolder(fileName);
				String name = GetFileTitle(fileName);
				String ext = GetFileExt(fileName);
				String fname = AFX(folder, Format("%s_%.3f%s", name, w[ifr], ext));
				SaveFile(fname, grid[ifr].GetString(false, false, ScatterDraw::GetDefaultCSVSeparator()));
			}
		} else if (type == ".xlsx") {
			xlnt::workbook wb;
			for (int ifr = 0; ifr < Nf; ++ifr) { 
				xlnt::worksheet ws;
				String title = Format("%.3f", w[ifr]); 
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
			SaveFile(fileName, grid.GetString(false, false, ScatterDraw::GetDefaultCSVSeparator()));
		else if (type == ".xlsx") {
			xlnt::workbook wb;
			xlnt::worksheet ws = wb.active_sheet();	
			ws.title("Data");
			XlsxFill(ws, grid, false);
			wb.save(~fileName); 
		}
	}
}

void Hydro::AddWave(int ib, double dx, double dy) {
	if (dx == 0 && dy == 0)
		return;
  	auto CalcF = [&](Forces &ex, const UVector<double> &k) {
    	UArray<MatrixXcd> exforce = clone(ex.force);
    	
	    for (int ih = 0; ih < Nh; ++ih) {
	        double angle = ToRad(head[ih]);
			double dist = dx*cos(angle) + dy*sin(angle);
		
    		int ib6 = ib*6;
			for (int ifr = 0; ifr < Nf; ++ifr) {
				double ph = k[ifr]*dist;
				for (int idf = 0; idf < 6; ++idf) 
					AddPhase(exforce[ih](ifr, idf + ib6), ph);		// Add the phase
			}
	    }
		ex.force = pick(exforce);
    };
    
    UVector<double> k(Nf);
	for (int ifr = 0; ifr < Nf; ++ifr) 
		k[ifr] = SeaWaves::WaveNumber(T[ifr], h, g_dim());
    	
	if (IsLoadedFex())
		CalcF(ex, k);
	if (IsLoadedFsc())
		CalcF(sc, k);
	if (IsLoadedFfk())
		CalcF(fk, k);

    auto CalcQTF = [&](UArray<UArray<UArray<MatrixXcd>>> &qtf, const UVector<double> &qk, bool isSum) {
        int sign = isSum ? 1 : -1;
		{
	        for (int ih = 0; ih < qh.size(); ++ih) {
	            double angle = ToRad(qh[ih].imag());
	            double dist = dx*cos(angle) + dy*sin(angle);
				for (int ifr1 = 0; ifr1 < qw.size(); ++ifr1) {
					for (int ifr2 = 0; ifr2 < qw.size(); ++ifr2) {
						std::complex<double> &v0 = qtf[ib][ih][0](ifr1, ifr2),
							 				 &v1 = qtf[ib][ih][1](ifr1, ifr2),
							 				 &v2 = qtf[ib][ih][2](ifr1, ifr2),
							 				 &v3 = qtf[ib][ih][3](ifr1, ifr2),
							 				 &v4 = qtf[ib][ih][4](ifr1, ifr2),
							 				 &v5 = qtf[ib][ih][5](ifr1, ifr2);
						
						double ph = (qk[ifr2] + sign*qk[ifr1])*dist;
						for (int idf = 0; idf < 6; ++idf) 
							AddPhase(qtf[ib][ih][idf](ifr1, ifr2), -ph);
					}
				}
	        }
		}
    };

    UVector<double> qk(Nf);
	for (int ifr = 0; ifr < qw.size(); ++ifr) 
		qk[ifr] = SeaWaves::WaveNumber_w(qw[ifr], h, g_dim());
			
	if (IsLoadedQTF(true)) 
		CalcQTF(qtfsum, qk, true);		
	if (IsLoadedQTF(false))	
		CalcQTF(qtfdif, qk, false);	
}


void Hydro::GetTranslationTo(const MatrixXd &to) {
	MatrixXd delta(3, Nb);
	for (int ib = 0; ib < Nb; ++ib) 
		for (int idf = 0; idf < 3; ++idf) 	
			delta(idf, ib) = to(idf, ib) - c0(idf, ib);
		
	
	auto CalcAB = [&](auto &A) {
        auto An = clone(A);
	
		for (int ib = 0; ib < Nb; ++ib) {
			int ib6 = ib*6;
			double dx = delta(0, ib);
			double dy = delta(1, ib);
			double dz = delta(2, ib);
			
			for (int jb = 0; jb < Nb; ++jb) {
				int jb6 = jb*6;
				
				for (int idof = 0; idof < 6; ++idof) {		// All dof are available?
					for (int jdof = 0; jdof < 6; ++jdof)
						if (!IsNum(A[ib6 + idof][jb6 + jdof][0])) 
							throw Exc("Coefficient translations require all DOFs to be available");
				}
				
				for (int iif = 0; iif < Nf; ++iif) {
					An[ib6 + 0][jb6 + 3][iif] += - dy*A[ib6 + 0][jb6 + 2][iif] + dz*A[ib6 + 0][jb6 + 1][iif];
					An[ib6 + 1][jb6 + 3][iif] += - dy*A[ib6 + 1][jb6 + 2][iif] + dz*A[ib6 + 1][jb6 + 1][iif];
					An[ib6 + 2][jb6 + 3][iif] += - dy*A[ib6 + 2][jb6 + 2][iif] + dz*A[ib6 + 2][jb6 + 1][iif];

					An[ib6 + 0][jb6 + 4][iif] += - dz*A[ib6 + 0][jb6 + 0][iif] + dx*A[ib6 + 0][jb6 + 2][iif];
					An[ib6 + 1][jb6 + 4][iif] += - dz*A[ib6 + 1][jb6 + 0][iif] + dx*A[ib6 + 1][jb6 + 2][iif];
					An[ib6 + 2][jb6 + 4][iif] += - dz*A[ib6 + 2][jb6 + 0][iif] + dx*A[ib6 + 2][jb6 + 2][iif];

					An[ib6 + 0][jb6 + 5][iif] += - dx*A[ib6 + 0][jb6 + 1][iif] + dy*A[ib6 + 0][jb6 + 0][iif];
					An[ib6 + 1][jb6 + 5][iif] += - dx*A[ib6 + 1][jb6 + 1][iif] + dy*A[ib6 + 1][jb6 + 0][iif];
					An[ib6 + 2][jb6 + 5][iif] += - dx*A[ib6 + 2][jb6 + 1][iif] + dy*A[ib6 + 2][jb6 + 0][iif];

					An[ib6 + 3][jb6 + 0][iif] += - dy*A[ib6 + 2][jb6 + 0][iif] + dz*A[ib6 + 1][jb6 + 0][iif];	
					An[ib6 + 3][jb6 + 1][iif] += - dy*A[ib6 + 2][jb6 + 1][iif] + dz*A[ib6 + 1][jb6 + 1][iif];
					An[ib6 + 3][jb6 + 2][iif] += - dy*A[ib6 + 2][jb6 + 2][iif] + dz*A[ib6 + 1][jb6 + 2][iif];

					An[ib6 + 4][jb6 + 0][iif] += - dz*A[ib6 + 0][jb6 + 0][iif] + dx*A[ib6 + 2][jb6 + 0][iif];	
					An[ib6 + 4][jb6 + 1][iif] += - dz*A[ib6 + 0][jb6 + 1][iif] + dx*A[ib6 + 2][jb6 + 1][iif];
					An[ib6 + 4][jb6 + 2][iif] += - dz*A[ib6 + 0][jb6 + 2][iif] + dx*A[ib6 + 2][jb6 + 2][iif];

					An[ib6 + 5][jb6 + 0][iif] += - dx*A[ib6 + 1][jb6 + 0][iif] + dy*A[ib6 + 0][jb6 + 0][iif];	
					An[ib6 + 5][jb6 + 1][iif] += - dx*A[ib6 + 1][jb6 + 1][iif] + dy*A[ib6 + 0][jb6 + 1][iif];
					An[ib6 + 5][jb6 + 2][iif] += - dx*A[ib6 + 1][jb6 + 2][iif] + dy*A[ib6 + 0][jb6 + 2][iif];

					An[ib6 + 3][jb6 + 3][iif] += -    2*dy*A[ib6 + 2][jb6 + 3][iif] +  2*dz*A[ib6 + 1][jb6 + 3][iif]
												 +   dy*dy*A[ib6 + 2][jb6 + 2][iif] + dz*dz*A[ib6 + 1][jb6 + 1][iif]
												 - 2*dy*dz*A[ib6 + 1][jb6 + 2][iif];   

		    		An[ib6 + 4][jb6 + 4][iif] += -    2*dz*A[ib6 + 0][jb6 + 4][iif] +  2*dx*A[ib6 + 2][jb6 + 4][iif]
												 +   dz*dz*A[ib6 + 0][jb6 + 0][iif] + dx*dx*A[ib6 + 2][jb6 + 2][iif]
												 - 2*dz*dx*A[ib6 + 0][jb6 + 2][iif];

					An[ib6 + 5][jb6 + 5][iif] += -    2*dx*A[ib6 + 1][jb6 + 5][iif] +  2*dy*A[ib6 + 0][jb6 + 5][iif]
												 +   dx*dx*A[ib6 + 1][jb6 + 1][iif] + dy*dy*A[ib6 + 0][jb6 + 0][iif]
												 - 2*dx*dy*A[ib6 + 0][jb6 + 1][iif];
				}
			}
		}
		A = pick(An);
    };
	
	if (IsLoadedA())
		CalcAB(A);
	if (IsLoadedAinf_w())
		CalcAB(Ainf_w);
	if (IsLoadedB())
		CalcAB(B);
    
    auto CalcA = [&](auto &A) {
        auto An = clone(A);
	
		for (int ib = 0; ib < Nb; ++ib) {
			int ib6 = ib*6;
			double dx = delta(0, ib);
			double dy = delta(1, ib);
			double dz = delta(2, ib);
			
			for (int jb = 0; jb < Nb; ++jb) {
				int jb6 = jb*6;
				
				for (int idof = 0; idof < 6; ++idof) {
					for (int jdof = 0; jdof < 6; ++jdof)
						if (!IsNum(A(ib6 + idof, jb6 + jdof))) 
							throw Exc("Coefficient translations require all DOFs to be available");
				}
				
				An(ib6 + 0, jb6 + 3) += - dy*A(ib6 + 0, jb6 + 2) + dz*A(ib6 + 0, jb6 + 1);
				An(ib6 + 1, jb6 + 3) += - dy*A(ib6 + 1, jb6 + 2) + dz*A(ib6 + 1, jb6 + 1);
				An(ib6 + 2, jb6 + 3) += - dy*A(ib6 + 2, jb6 + 2) + dz*A(ib6 + 2, jb6 + 1);
	
				An(ib6 + 0, jb6 + 4) += - dz*A(ib6 + 0, jb6 + 0) + dx*A(ib6 + 0, jb6 + 2);
				An(ib6 + 1, jb6 + 4) += - dz*A(ib6 + 1, jb6 + 0) + dx*A(ib6 + 1, jb6 + 2);
				An(ib6 + 2, jb6 + 4) += - dz*A(ib6 + 2, jb6 + 0) + dx*A(ib6 + 2, jb6 + 2);
	
				An(ib6 + 0, jb6 + 5) += - dx*A(ib6 + 0, jb6 + 1) + dy*A(ib6 + 0, jb6 + 0);
				An(ib6 + 1, jb6 + 5) += - dx*A(ib6 + 1, jb6 + 1) + dy*A(ib6 + 1, jb6 + 0);
				An(ib6 + 2, jb6 + 5) += - dx*A(ib6 + 2, jb6 + 1) + dy*A(ib6 + 2, jb6 + 0);
	
				An(ib6 + 3, jb6 + 0) += - dy*A(ib6 + 2, jb6 + 0) + dz*A(ib6 + 1, jb6 + 0);	
				An(ib6 + 3, jb6 + 1) += - dy*A(ib6 + 2, jb6 + 1) + dz*A(ib6 + 1, jb6 + 1);
				An(ib6 + 3, jb6 + 2) += - dy*A(ib6 + 2, jb6 + 2) + dz*A(ib6 + 1, jb6 + 2);
	
				An(ib6 + 4, jb6 + 0) += - dz*A(ib6 + 0, jb6 + 0) + dx*A(ib6 + 2, jb6 + 0);	
				An(ib6 + 4, jb6 + 1) += - dz*A(ib6 + 0, jb6 + 1) + dx*A(ib6 + 2, jb6 + 1);
				An(ib6 + 4, jb6 + 2) += - dz*A(ib6 + 0, jb6 + 2) + dx*A(ib6 + 2, jb6 + 2);
	
				An(ib6 + 5, jb6 + 0) += - dx*A(ib6 + 1, jb6 + 0) + dy*A(ib6 + 0, jb6 + 0);	
				An(ib6 + 5, jb6 + 1) += - dx*A(ib6 + 1, jb6 + 1) + dy*A(ib6 + 0, jb6 + 1);
				An(ib6 + 5, jb6 + 2) += - dx*A(ib6 + 1, jb6 + 2) + dy*A(ib6 + 0, jb6 + 2);
				
				An(ib6 + 3, jb6 + 3) += -    2*dy*A(ib6 + 2, jb6 + 3) +  2*dz*A(ib6 + 1, jb6 + 3)
									 	+   dy*dy*A(ib6 + 2, jb6 + 2) + dz*dz*A(ib6 + 1, jb6 + 1)
										- 2*dy*dz*A(ib6 + 1, jb6 + 2);   
	
	    		An(ib6 + 4, jb6 + 4) += -    2*dz*A(ib6 + 0, jb6 + 4) +  2*dx*A(ib6 + 2, jb6 + 4)
										+   dz*dz*A(ib6 + 0, jb6 + 0) + dx*dx*A(ib6 + 2, jb6 + 2)
										- 2*dz*dx*A(ib6 + 0, jb6 + 2);
	
				An(ib6 + 5, jb6 + 5) += -    2*dx*A(ib6 + 1, jb6 + 5) +  2*dy*A(ib6 + 0, jb6 + 5)
										+   dx*dx*A(ib6 + 1, jb6 + 1) + dy*dy*A(ib6 + 0, jb6 + 0)
										- 2*dx*dy*A(ib6 + 0, jb6 + 1);
			} 
		}
		A = pick(An);
    };
    
    if (IsLoadedA0())
		CalcA(A0);
    if (IsLoadedAinf())
		CalcA(Ainf);
	if (IsLoadedDlin())
		CalcA(Dlin);
		    
    auto CalcF = [&](Forces &ex) {
    	UArray<MatrixXcd> exforce = clone(ex.force);
    	
	    for (int ih = 0; ih < Nh; ++ih) {
	        double angle = ToRad(head[ih]);
	    	for (int ib = 0; ib < Nb; ++ib) {
	    		double dx = delta(0, ib);
				double dy = delta(1, ib);
				double dz = delta(2, ib);
				double dist = dx*cos(angle) + dy*sin(angle);
			
	    		int ib6 = ib*6;
				for (int ifr = 0; ifr < Nf; ++ifr) {
					exforce[ih](ifr, 3 + ib6) += -dy*exforce[ih](ifr, 2 + ib6) + dz*exforce[ih](ifr, 1 + ib6);
	    			exforce[ih](ifr, 4 + ib6) += -dz*exforce[ih](ifr, 0 + ib6) + dx*exforce[ih](ifr, 2 + ib6);
	    			exforce[ih](ifr, 5 + ib6) += -dx*exforce[ih](ifr, 1 + ib6) + dy*exforce[ih](ifr, 0 + ib6);
				}
	    	}
	    }
		ex.force = pick(exforce);
    };
    	
	if (IsLoadedFex())
		CalcF(ex);
	if (IsLoadedFsc())
		CalcF(sc);
	if (IsLoadedFfk())
		CalcF(fk);

    auto CalcMD = [&]() {
    	UArray<UArray<UArray<VectorXd>>> mdn = clone(md);
    	
    	for (int ib = 0; ib < Nb; ++ib) {
    		double dx = delta(0, ib);
			double dy = delta(1, ib);
			double dz = delta(2, ib);    		
    		for (int ih = 0; ih < mdhead.size(); ++ih) {
    			for (int ifr = 0; ifr < Nf; ++ifr) {
    				mdn[ib][ih][3](ifr) += -dy*mdn[ib][ih][2](ifr) + dz*mdn[ib][ih][1](ifr);
    				mdn[ib][ih][4](ifr) += -dz*mdn[ib][ih][0](ifr) + dx*mdn[ib][ih][2](ifr);
    				mdn[ib][ih][5](ifr) += -dx*mdn[ib][ih][1](ifr) + dy*mdn[ib][ih][0](ifr);
				}
	    	}
	    }
		md = pick(mdn);
    };
    
	if (IsLoadedMD())
		CalcMD();
	
    auto CalcQTF = [&](UArray<UArray<UArray<MatrixXcd>>> &qtf, bool isSum) {
        int sign = isSum ? 1 : -1;
		for (int ib = 0; ib < Nb; ++ib) {
			double dx = delta(0, ib);
			double dy = delta(1, ib);
			double dz = delta(2, ib);
	        for (int ih = 0; ih < qh.size(); ++ih) {
	            double angle = ToRad(qh[ih].imag());
	            double dist = dx*cos(angle) + dy*sin(angle);
				for (int ifr1 = 0; ifr1 < qw.size(); ++ifr1) {
					for (int ifr2 = 0; ifr2 < qw.size(); ++ifr2) {
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
	for (int ih = int(qh.size())-1; ih >= 0; --ih) { 
		if (qh[ih].real() != qh[ih].imag()) {
			Remove(qh, ih);
			for (int ib = 0; ib < Nb; ++ib) {
				if (IsLoadedQTF(true)) 
					qtfsum[ib].Remove(ih);
				if (IsLoadedQTF(false)) 
					qtfdif[ib].Remove(ih);
			}
		}
	}
			
	if (IsLoadedQTF(true)) 
		CalcQTF(qtfsum, true);		
	if (IsLoadedQTF(false))	
		CalcQTF(qtfdif, false);
	
	if (IsLoadedM()) {
		for (int ib = 0; ib < Nb; ++ib) 
			Surface::TranslateInertia66(M[ib], Point3D(cg.col(ib)), Point3D(c0.col(ib)), Point3D(to.col(ib)));
	}
	
	if (IsLoadedKirf()) {
		double maxT = GetK_IRF_MaxT();
		if (maxT < 0)
			maxT = bem->maxTimeA;
		else if (bem->maxTimeA < maxT) 
			maxT = bem->maxTimeA;

		GetK_IRF(maxT, bem->numValsA);
	}
	
	// Some previous data are now invalid
	rao.Clear();	
	C.Clear();
	
/*	for (int ib = 0; ib < Nb; ++ib) {
		double dx = delta(0, ib);
		double dy = delta(1, ib);			// Previous translation has moved implicitly the wave origin
		AddWave(ib, -dx, -dy);				// Translate the wave back to the original position
	}*/
	
	c0 = clone(to);
		
	if (!AfterLoad()) {
		String error = GetLastError();
		throw Exc(Format(t_("Problem translating model: '%s'\n%s"), error));	
	}
}

void Hydro::ResetForces1st(Hydro::FORCE force) {
	if (force == Hydro::FK) {
		if (!IsLoadedFfk())
			return;
		if (!IsLoadedFsc() && !IsLoadedFex())
			return;
		
		if (IsLoadedFsc()) 
			ex = clone(sc);
		else {
			for (int ih = 0; ih < Nh; ++ih) {
				for (int ifr = 0; ifr < Nf; ++ifr) 
					for (int i = 0; i < Nb*6; ++i) 
						if (IsNum(sc.force[ih](ifr, i))) 
							ex.force[ih](ifr, i) = ex.force[ih](ifr, i) - fk.force[ih](ifr, i);
			}		
		}
		fk.Clear();
	} else if (force == Hydro::SCATTERING) {
		if (!IsLoadedFsc())
			return;
		if (!IsLoadedFfk() && !IsLoadedFex())
			return;
		
		if (IsLoadedFfk()) 
			ex = clone(fk);
		else {
			for (int ih = 0; ih < Nh; ++ih) 
				for (int ifr = 0; ifr < Nf; ++ifr) 
					for (int i = 0; i < Nb*6; ++i) 
						if (IsNum(sc.force[ih](ifr, i))) 
							ex.force[ih](ifr, i) = ex.force[ih](ifr, i) - sc.force[ih](ifr, i);
		}
		sc.Clear();		
	} else {
		ex.Clear();		
		sc.Clear();		
		fk.Clear();		
	}
	md.Clear();
}

void Hydro::ResetForces(Hydro::FORCE force, bool forceMD, Hydro::FORCE forceQtf) {
	if (force != Hydro::NONE)
		Hydro::ResetForces1st(force);

	if (forceMD)
		md.Clear();
	
	if (forceQtf == Hydro::ALL || forceQtf == Hydro::QTFSUM) 
		qtfsum.Clear();
	if (forceQtf == Hydro::ALL || forceQtf == Hydro::QTFDIF) 
		qtfdif.Clear();
}

void Hydro::MultiplyDOF(double factor, const UVector<int> &_idDOF, bool a, bool b, bool diag, bool f, bool ismd, bool qtf) {
	if (_idDOF.size() == 0) 
		return;
	
	UVector<int> idDOF;
	for (int idof = 0; idof < _idDOF.size(); ++idof)
		for (int ib = 0; ib < Nb; ++ib)
			idDOF << _idDOF[idof] + ib*6;
	
	auto MultiplyAB = [&](UArray<UArray<VectorXd>> &A) {
		for (int idf = 0; idf < 6*Nb; ++idf) {
			for (int jdf = 0; jdf < 6*Nb; ++jdf) {
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
		MultiplyAB(A);
	if (a && IsLoadedAinf_w())
		MultiplyAB(Ainf_w);
	if (b && IsLoadedB())
		MultiplyAB(B);

	auto MultiplyAinfA0 = [&](MatrixXd &A) {
		for (int idf = 0; idf < 6*Nb; ++idf) {
			for (int jdf = 0; jdf < 6*Nb; ++jdf) {
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
		MultiplyAinfA0(Ainf);
	if (a && IsLoadedA0()) 
		MultiplyAinfA0(A0);
		
	auto MultiplyF = [&](Forces &ex) {
		for (int ih = 0; ih < Nh; ++ih) 
			for (int ifr = 0; ifr < Nf; ++ifr) 
				for (int idof = 0; idof < idDOF.size(); ++idof) 
					ex.force[ih](ifr, idDOF[idof]) *= factor;
	};
	if (f && IsLoadedFex())
		MultiplyF(ex);
	if (f && IsLoadedFsc())
		MultiplyF(sc);
	if (f && IsLoadedFfk())
		MultiplyF(fk);	
	if (f && IsLoadedRAO())
		MultiplyF(rao);

	auto MultiplyMD = [&]() {
		for (int ib = 0; ib < Nb; ++ib) 
    		for (int ih = 0; ih < mdhead.size(); ++ih) 
    			for (int idof = 0; idof < idDOF.size(); ++idof) 
	    			for (int ifr = 0; ifr < Nf; ++ifr) 
	    				md[ib][ih][idDOF[idof]](ifr) *= factor;
	};
	if (ismd && IsLoadedMD()) 
		MultiplyMD();
			
	auto MultiplySumDif = [&](UArray<UArray<UArray<MatrixXcd>>> &qtf) {
		for (int ib = 0; ib < Nb; ++ib)
	        for (int ih = 0; ih < qh.size(); ++ih) 
				for (int idof = 0; idof < _idDOF.size(); ++idof) 
					qtf[ib][ih][_idDOF[idof]] *= factor;													
	};
	if (qtf && IsLoadedQTF(true)) 
		MultiplySumDif(qtfsum);
	if (qtf && IsLoadedQTF(false))
		MultiplySumDif(qtfdif);
	
	// Some previous data is now invalid
	Kirf.Clear();
	
	if (!AfterLoad()) {
		String error = GetLastError();
		throw Exc(Format(t_("Problem reseting DOF: '%s'\n%s"), error));	
	}
}

void Hydro::SwapDOF(int ib1, int ib2) {
	for (int idof = 0; idof < 6; ++idof)	
		SwapDOF(ib1, idof, ib2, idof);	
	
	c0.col(ib1).swap(c0.col(ib2));
	cg.col(ib1).swap(cg.col(ib2));
	cb.col(ib1).swap(cb.col(ib2));
	Swap(Vo[ib1], Vo[ib2]);
	Swap(names[ib1], names[ib2]);
}

void Hydro::SwapDOF(int ib1, int idof1, int ib2, int idof2) {
	auto SwapAB = [&](UArray<UArray<VectorXd>> &A) {
		UArray<UArray<VectorXd>> An(6*Nb);
		for (int i = 0; i < 6*Nb; ++i) 
			An[i].SetCount(6*Nb);
		
		for (int idof = 0; idof < 6*Nb; ++idof) {
			for (int jdof = 0; jdof < 6*Nb; ++jdof) {
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
		SwapAB(A);
	if (IsLoadedAinf_w())
		SwapAB(Ainf_w);
	if (IsLoadedB())
		SwapAB(B);

		
	auto SwapAinfA0 = [&](MatrixXd &A) {
		Swap(A, idof1+6*ib1, idof2+6*ib2);
    };	
	if (IsLoadedAinf()) 
		SwapAinfA0(Ainf);
	if (IsLoadedA0()) 
		SwapAinfA0(A0);
			  
	auto SwapF = [&](Forces &ex) {
		for (int ih = 0; ih < Nh; ++ih) {
			MatrixXcd n(Nf, 6*Nb);
			for (int idof = 0; idof < 6*Nb; ++idof) {
				int idofn = idof;
				if (idofn == idof1+6*ib1)
					idofn = idof2+6*ib2;
				else if (idofn == idof2+6*ib2)
					idofn = idof1+6*ib1;
				
	    		const VectorXcd &m = ex.force[ih].col(idof);
				n.col(idofn) = m;
	    	}
	    	ex.force[ih] = pick(n);
		}
	};
	if (IsLoadedFex())
		SwapF(ex);
	if (IsLoadedFsc())
		SwapF(sc);
	if (IsLoadedFfk())
		SwapF(fk);	
	if (IsLoadedRAO())
		SwapF(rao);

	auto SwapMD = [&]() {
		for (int ih = 0; ih < mdhead.size(); ++ih) 
			Swap(md[ib1][ih][idof1], md[ib2][ih][idof2]);
	};
	if (IsLoadedMD(true)) 
		SwapMD();
		
	auto SwapSumDif = [&](UArray<UArray<UArray<MatrixXcd>>> &qtf) {
        for (int ih = 0; ih < qh.size(); ++ih) 
			Swap(qtf[ib1][ih][idof1], qtf[ib2][ih][idof2]); 		
	};
	if (IsLoadedQTF(true)) 
		SwapSumDif(qtfsum);
	if (IsLoadedQTF(false))
		SwapSumDif(qtfdif);

	if (IsLoadedC()) 
		Swap(C[ib1], C[ib2], idof1, idof2);
	
	if (IsLoadedM()) 
		Swap(M[ib1], M[ib2], idof1, idof2);

	if (IsLoadedKirf()) // Kirf structure is like A and B
		SwapAB(Kirf);
		
	if (!AfterLoad()) {
		String error = GetLastError();
		throw Exc(Format(t_("Problem swaping DOF: '%s'\n%s"), error));	
	}
}

void Hydro::DeleteFrequencies(const UVector<int> &idFreq) {
	if (idFreq.size() > 0) {
		auto DeleteAB = [&](UArray<UArray<VectorXd>> &A) {
	        UArray<UArray<VectorXd>> An;
		
			An.SetCount(6*Nb);
			for (int idof = 0; idof < 6*Nb; ++idof) {
				An[idof].SetCount(6*Nb);
				for (int jdof = 0; jdof < 6*Nb; ++jdof) {
					An[idof][jdof].resize(Nf - idFreq.size());	
					int i = 0, j = 0;
					for (int iif = 0; iif < Nf; ++iif) {
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
			DeleteAB(A);
		if (IsLoadedAinf_w())
			DeleteAB(Ainf_w);
		if (IsLoadedB())
			DeleteAB(B);
	
		auto DeleteF = [&](Forces &ex) {
	        Forces _ex;
		
			_ex.force.SetCount(Nh);
		    for (int ih = 0; ih < Nh; ++ih) {
		        _ex.force[ih].resize(Nf-idFreq.size(), 6*Nb);
		    	for (int idof = 0; idof < 6*Nb; ++idof) {
					int i = 0, j = 0;
					for (int iif = 0; iif < Nf; ++iif) {
						if (j >= idFreq.size() || iif != idFreq[j]) 
							_ex.force[ih](i++, idof) = ex.force[ih](iif, idof);
						else 
							j++;
					}
		    	}
		    }
		    ex = pick(_ex);
	    };	
	
		if (IsLoadedFex())
			DeleteF(ex);
		if (IsLoadedFsc())
			DeleteF(sc);
		if (IsLoadedFfk())
			DeleteF(fk);	
		if (IsLoadedRAO())
			DeleteF(rao);

		auto DeleteMD = [&]() {
			UArray<UArray<UArray<VectorXd>>> mdn;

			mdn.SetCount(Nb);
			for (int ib = 0; ib < Nb; ++ib) {
				mdn[ib].SetCount(int(mdhead.size()));
	    		for (int ih = 0; ih < mdhead.size(); ++ih) {
	    			mdn[ib][ih].SetCount(6);
	    			for (int idf = 0; idf < 6; ++idf) {
	    				mdn[ib][ih][idf].setConstant(Nf-idFreq.size());
	    				int i = 0, j = 0;
						for (int iif = 0; iif < Nf; ++iif) {
							if (j >= idFreq.size() || iif != idFreq[j])
								md[ib][ih][idf](i++) = md[ib][ih][idf](iif);		
							else 
								j++;
						}
	    			}
	    		}
			}
		    md = pick(mdn);			
	    };

		if (IsLoadedMD())
			DeleteMD();
				
		int j = idFreq.size()-1;	
		for (int i = w.size()-1; i >= 0 && j >= 0; --i) {
			if (i == idFreq[j]) {	
				w.Remove(i);
				T.Remove(i);
				j--;
			}
		}
		Nf = w.size();
	}
}

void Hydro::DeleteFrequenciesQTF(const UVector<int> &idFreqQTF) {
	if (idFreqQTF.size() > 0) {
		UVector<int> vids;
		LinSpaced(vids, int(qw.size()), 0, int(qw.size())-1);
		for (int i = idFreqQTF.size()-1; i >= 0; --i) 
			vids.Remove(idFreqQTF[i]);
		VectorXi ids;
		::Copy(vids, ids);
		qw = VectorXd(qw(ids));
		
		auto DeleteSumDif = [&](UArray<UArray<UArray<MatrixXcd>>> &qtf) {
			for (int ib = 0; ib < Nb; ++ib)
		        for (int ih = 0; ih < qh.size(); ++ih) 
					for (int idf = 0; idf < 6; ++idf) {
						MatrixXcd &m = qtf[ib][ih][idf];
						m = MatrixXcd(m(indexing::all, ids));
						m = MatrixXcd(m(ids, indexing::all));
					}
		};
		if (IsLoadedQTF(true)) 
			DeleteSumDif(qtfsum);
		if (IsLoadedQTF(false))
			DeleteSumDif(qtfdif);
	}
}

void Hydro::DeleteHeadings(const UVector<int> &idHead) {
	if (idHead.size() > 0) {
		auto DeleteF = [&](Forces &ex) {
			int j = idHead.size()-1;	
			for (int i = head.size()-1; i >= 0 && j >= 0; --i) {
				if (i == idHead[j]) {	
					ex.force.Remove(i);
					j--;
				}
			}
	    };	
	
		if (IsLoadedFex())
			DeleteF(ex);
		if (IsLoadedFsc())
			DeleteF(sc);
		if (IsLoadedFfk())
			DeleteF(fk);	
		if (IsLoadedRAO())
			DeleteF(rao);
	
		int j = idHead.size()-1;	
		for (int i = head.size()-1; i >= 0 && j >= 0; --i) {
			if (i == idHead[j]) {	
				head.Remove(i);
				j--;
			}
		}
		Nh = head.size();
	}
}

void Hydro::DeleteHeadingsMD(const UVector<int> &idHead) {
	if (idHead.size() > 0) {
		auto DeleteMD = [&]() {
			for (int ib = 0; ib < Nb; ++ib) {
	    		int j = idHead.size()-1;	
				for (int i = int(mdhead.size())-1; i >= 0 && j >= 0; --i) {
	    			if (i == idHead[j]) {	
						md[ib].Remove(i);
						j--;
					}
				}
			}
	    };

		if (IsLoadedMD())
			DeleteMD();
		
		UArray<std::complex<double>> mdh;
		::Copy(mdhead, mdh);
		int j = idHead.size()-1;	
		for (int i = mdh.size()-1; i >= 0 && j >= 0; --i) {
			if (i == idHead[j]) {	
				mdh.Remove(i);
				j--;
			}
		}
		::Copy(mdh, mdhead);
	}
}

void Hydro::DeleteHeadingsQTF(const UVector<int> &idHeadQTF) {
	if (idHeadQTF.size() > 0) {
		UVector<int> vids;
		LinSpaced(vids, int(qh.size()), 0, int(qh.size())-1);
		for (int i = idHeadQTF.size()-1; i >= 0; --i) 
			vids.Remove(idHeadQTF[i]);
		VectorXi ids;
		::Copy(vids, ids);
		qh = VectorXcd(qh(ids));
			
		auto DeleteSumDif = [&](UArray<UArray<UArray<MatrixXcd>>> &qtf) {
			for (int ib = 0; ib < Nb; ++ib)
		        for (int ih = idHeadQTF.size()-1; ih >= 0; --ih)
		        	qtf[ib].Remove(idHeadQTF[ih]);
		};
		if (IsLoadedQTF(true)) 
			DeleteSumDif(qtfsum);
		if (IsLoadedQTF(false))
			DeleteSumDif(qtfdif);
	}
}

void Hydro::FillFrequencyGapsABForces(bool zero, int maxFreq) {
	if (w.size() == 0)
		return;

	VectorXd w_, nw;
	::Copy(w, w_);
	
	UVector<int> idsx, w0x;
	GapFillingAxisParams(w_, maxFreq, idsx, w0x, nw);
	
	auto FillAB = [&](UArray<UArray<VectorXd>> &A) {
		for (int idof = 0; idof < 6*Nb; ++idof) {
			for (int jdof = 0; jdof < 6*Nb; ++jdof) {
				VectorXd nm;
				const VectorXd &m = A[idof][jdof];
				GapFilling(w_, m, idsx, w0x, nw, nm, zero, maxFreq);					
				A[idof][jdof] = pick(nm);
			}
		}
    };
		
	if (IsLoadedA())
		FillAB(A);
	if (IsLoadedAinf_w())
		FillAB(Ainf_w);
	if (IsLoadedB())
		FillAB(B);
	
	auto FillF = [&](Forces &ex) {
	    for (int ih = 0; ih < Nh; ++ih) {
	        MatrixXcd nmn(nw.size(), 6*Nb);
	    	for (int idof = 0; idof < 6*Nb; ++idof) {
	    		VectorXcd nm;
	    		const VectorXcd &m = ex.force[ih].col(idof);
	    		GapFilling(w_, m, idsx, w0x, nw, nm, zero, maxFreq);					
				nmn.col(idof) = nm;
	    	}
	    	ex.force[ih] = pick(nmn);
	    }
    };	

	if (IsLoadedFex())
		FillF(ex);
	if (IsLoadedFsc())
		FillF(sc);
	if (IsLoadedFfk())
		FillF(fk);	
	if (IsLoadedRAO())
		FillF(rao);	
	
	auto FillMD = [&]() {
		for (int ib = 0; ib < Nb; ++ib) {
    		for (int ih = 0; ih < mdhead.size(); ++ih) {
    			for (int idf = 0; idf < 6; ++idf) {
    				VectorXd nm(nw.size());
	    			VectorXd &m = md[ib][ih][idf];
					GapFilling(w_, m, idsx, w0x, nw, nm, zero, maxFreq);
					md[ib][ih][idf] = pick(nm);
    			}
    		}
		}
    };		
	
	if (IsLoadedMD())
		FillMD();		
	
	Nf = int(nw.size());
	::Copy(nw, w);
	T.SetCount(Nf);
	for (int i = 0; i < Nf; ++i) 
		T[i] = 2*M_PI/w[i];
}

void Hydro::FillFrequencyGapsQTF(bool zero, int maxFreq) {
	if (qw.size() == 0)
		return;
	
	VectorXd nw;

	UVector<int> idsx, w0x;
	GapFillingAxisParams(qw, maxFreq, idsx, w0x, nw);
	
	auto FillSumDif = [&](UArray<UArray<UArray<MatrixXcd>>> &qtf) {
		for (int ib = 0; ib < Nb; ++ib) {
	        for (int ih = 0; ih < qh.size(); ++ih) {
				for (int idof = 0; idof < 6; ++idof) {
					MatrixXcd nm, &m = qtf[ib][ih][idof];
					GapFilling(qw, qw, m, idsx, w0x, idsx, w0x, nw, nw, nm, zero, maxFreq);					
					m = pick(nm);
				}
	        }
		}
	};

	if (IsLoadedQTF(true)) 
		FillSumDif(qtfsum);
	if (IsLoadedQTF(false)) 
		FillSumDif(qtfdif);
	
	qw = pick(nw);
}

void Hydro::FillFrequencyGapsABForcesZero() {
	if (w.size() == 0)
		return;

	auto FillAB = [&](UArray<UArray<VectorXd>> &A) {
		for (int idof = 0; idof < 6*Nb; ++idof) {
			for (int jdof = 0; jdof < 6*Nb; ++jdof) {
				VectorXd &a = A[idof][jdof];
				if (a.size() == 0 || !IsNum(a(0)))
					a = VectorXd::Zero(Nf);
			}
		}
    };
		
	if (IsLoadedA())
		FillAB(A);
	if (IsLoadedAinf_w())
		FillAB(Ainf_w);
	if (IsLoadedB())
		FillAB(B);

	auto FillA = [&](MatrixXd &A) {
		if (A.size() == 0)
			A = MatrixXd::Zero(6*Nb, 6*Nb);
		else 
			A = A.unaryExpr([](double x){return IsNum(x) ? x : 0;});		// Replace NaN with 0
    };
		
	if (IsLoadedAinf())
		FillA(Ainf);
	if (IsLoadedA0())
		FillA(A0);
	if (IsLoadedDlin())
		FillA(Dlin);
	
	auto FillF = [&](Forces &ex) {
	    for (int ih = 0; ih < Nh; ++ih) {
	        MatrixXcd nmn(Nf, 6*Nb);
	    	for (int idof = 0; idof < 6*Nb; ++idof) {
	    		const VectorXcd &m = ex.force[ih].col(idof);
	    		if (!IsNum(m(0))) 
	    			nmn.col(idof) = VectorXcd::Zero(Nf);
	    		else 
	    			nmn.col(idof) = m;
	    	}
	    	ex.force[ih] = pick(nmn);
	    }
    };	

	if (IsLoadedFex())
		FillF(ex);
	if (IsLoadedFsc())
		FillF(sc);
	if (IsLoadedFfk())
		FillF(fk);	
	if (IsLoadedRAO())
		FillF(rao);	
	
	auto FillMD = [&]() {
		for (int ib = 0; ib < Nb; ++ib) {
    		for (int ih = 0; ih < mdhead.size(); ++ih) {
    			for (int idf = 0; idf < 6; ++idf) {
    				if (md[ib][ih][idf].size() == 0 || !IsNum(md[ib][ih][idf][0]))
						md[ib][ih][idf] = VectorXd::Zero(Nf);
    			}
    		}
		}
    };		
	
	if (IsLoadedMD())
		FillMD();			
}

void Hydro::FillFrequencyGapsQTFZero() {
	if (qw.size() == 0)
		return;
	
	Eigen::Index nf = qw.size();
	
	auto FillSumDif = [&](UArray<UArray<UArray<MatrixXcd>>> &qtf) {
		for (int ib = 0; ib < Nb; ++ib) {
	        for (int ih = 0; ih < qh.size(); ++ih) {
				for (int idof = 0; idof < 6; ++idof) {
					MatrixXcd &m = qtf[ib][ih][idof];
					if (m.size() == 0 || !IsNum(m(0,0)))
						m = MatrixXcd::Zero(nf, nf);
				}
	        }
		}
	};

	if (IsLoadedQTF(true)) 
		FillSumDif(qtfsum);
	if (IsLoadedQTF(false)) 
		FillSumDif(qtfdif);
}

void Hydro::CopyQTF_MD() {
	mdtype = 9;
	::Copy(qh, mdhead);
	
	Initialize_MD(md, Nb, int(mdhead.size()), Nf);

	VectorXd ww;
	::Copy(w, ww);
	
	for (int ib = 0; ib < Nb; ++ib) {
        for (int ih = 0; ih < qh.size(); ++ih) {
			for (int idof = 0; idof < 6; ++idof) {
				const MatrixXcd &m = qtfdif[ib][ih][idof];
				VectorXd diag = Eigen::abs(m.diagonal().array());
				ResampleY(qw, diag, ww, md[ib][ih][idof]);
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
		for (int idf = 0; idf < 6*Nb; ++idf) 
			for (int jdf = idf+1; jdf < 6*Nb; ++jdf) 
				A[idf][jdf] = A[jdf][idf] = AvgSafe(A[idf][jdf], A[jdf][idf]);
    };		
	if (IsLoadedA())
		SymmetrizeAB(A);
	if (IsLoadedAinf_w())
		SymmetrizeAB(Ainf_w);
	if (IsLoadedB())
		SymmetrizeAB(B);

	auto SymmetrizeAinfA0 = [&](MatrixXd &A) {
		for (int idf = 0; idf < 6*Nb; ++idf) 
			for (int jdf = idf+1; jdf < 6*Nb; ++jdf) 
				A(idf, jdf) = A(jdf, idf) = AvgSafe(A(idf, jdf), A(jdf, idf));
    };	
	if (IsLoadedAinf()) 
		SymmetrizeAinfA0(Ainf);
	if (IsLoadedA0()) 
		SymmetrizeAinfA0(A0);
		
	auto SymmetrizeSumDif = [&](UArray<UArray<UArray<MatrixXcd>>> &qtf, bool isSum) {
		for (int ib = 0; ib < Nb; ++ib) {
	        for (int ih = 0; ih < qh.size(); ++ih) {
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
		SymmetrizeSumDif(qtfsum, true);
	if (IsLoadedQTF(false))
		SymmetrizeSumDif(qtfdif, false);
	
	if (!AfterLoad()) {
		String error = GetLastError();
		throw Exc(Format(t_("Problem symmetrizing data: '%s'\n%s"), error));	
	}
}

double Hydro::GetQTFVal(int ib, int idof, int idh, int ifr1, int ifr2, bool isSum, char what, bool getDim) const {
	const UArray<UArray<UArray<MatrixXcd>>> &qtf = isSum ? qtfsum : qtfdif;
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
	
	const UArray<UArray<UArray<MatrixXcd>>> &qtf = isSum ? qtfsum : qtfdif;
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

void Heal();
void Load(const VectorXd &w, const VectorXd &A, const VectorXd &B, double maxT, int num);
void Save(const VectorXd &w, VectorXd &A, VectorXd &Ainfw, double &ainf, VectorXd &B, 
			VectorXd &Tirf, VectorXd &Kinf);
			   				

	
	