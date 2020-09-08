#include "BEMRosetta.h"


bool Wamit::Load(String file) {
	hd().code = Hydro::WAMIT;
	hd().file = file;	
	hd().name = GetFileTitle(file);
		
	try {
		if (GetFileExt(file) == ".out") {
			BEMData::Print("\n\n" + Format(t_("Loading out file '%s'"), file));
			if (!Load_out()) {
				BEMData::PrintWarning("\n" + Format(t_("File '%s' not found"), file));
				return false;
			}
			String fileSC = ForceExt(file, ".3sc");
			BEMData::Print("\n- " + Format(t_("Scattering file '%s'"), GetFileName(fileSC)));
			if (!Load_Scattering(fileSC))
				BEMData::PrintWarning(S(": **") + t_("Not found") + "**");
			String fileFK = ForceExt(file, ".3fk");
			BEMData::Print("\n- " + Format(t_("Froude-Krylov file '%s'"), GetFileName(fileFK)));
			if (!Load_FK(fileFK))
				BEMData::PrintWarning(S(": **") + t_("Not found") + "**");
		} else if (S(".1.3.hst.4.12d.12s").Find(GetFileExt(file)) >= 0) {
			String file1 = ForceExt(file, ".1");
			BEMData::Print("\n- " + Format(t_("Hydrodynamic coefficients A and B .1 file '%s'"), GetFileName(file1)));
			if (!Load_1(file1))
				BEMData::PrintWarning(S(": **") + t_("Not found") + "**");
			
			String file3 = ForceExt(file, ".3");
			BEMData::Print("\n- " + Format(t_("Diffraction exciting .3 file '%s'"), GetFileName(file3)));
			if (!Load_3(file3))
				BEMData::PrintWarning(S(": **") + t_("Not found") + "**");
			
			String fileHST = ForceExt(file, ".hst");
			BEMData::Print("\n- " + Format(t_("Hydrostatic restoring file '%s'"), GetFileName(fileHST)));
			if (!Load_hst(fileHST))
				BEMData::PrintWarning(S(": **") + t_("Not found") + "**");
		
			String fileRAO = ForceExt(file, ".4");
			BEMData::Print("\n- " + Format(t_("RAO file '%s'"), GetFileName(fileRAO)));
			if (!Load_4(fileRAO))
				BEMData::PrintWarning(S(": **") + t_("Not found") + "**");

			if (IsNull(hd().Nh))
				hd().Nh = 0;
			if (IsNull(hd().Nf))
				hd().Nf = 0;

			String file12s = ForceExt(file, ".12s");
			BEMData::Print("\n- " + Format(t_("Second order sum coefficients .12s file '%s'"), GetFileName(file12s)));
			if (!Load_12(file12s, true))
				BEMData::PrintWarning(S(": **") + t_("Not found") + "**");
			
			String file12d = ForceExt(file, ".12d");
			BEMData::Print("\n- " + Format(t_("Second order mean drift coefficients .12d file '%s'"), GetFileName(file12d)));
			if (!Load_12(file12d, false))
				BEMData::PrintWarning(S(": **") + t_("Not found") + "**");
		}
		
		if (IsNull(hd().Nb))
			return false;
		
		hd().dof.Clear();	hd().dof.SetCount(hd().Nb, 0);
		for (int i = 0; i < hd().Nb; ++i)
			hd().dof[i] = 6;
	} catch (Exc e) {
		BEMData::PrintError(Format("\n%s: %s", t_("Error"), e));
		hd().lastError = e;
		return false;
	}
	
	return true;
}

void Wamit::Save(String file, bool force_T, int qtfHeading) {
	try {
		if (hd().IsLoadedA() && hd().IsLoadedB()) {
			String file1 = ForceExt(file, ".1");
			BEMData::Print("\n- " + Format(t_("Hydrodynamic coefficients A and B file '%s'"), GetFileName(file1)));
			Save_1(file1, force_T);
		}
		if (hd().IsLoadedFex()) {
			String file3 = ForceExt(file, ".3");
			BEMData::Print("\n- " + Format(t_("Diffraction exciting file '%s'"), GetFileName(file3)));
			Save_3(file3, force_T);
		}
		if (hd().IsLoadedC()) {
			String fileHST = ForceExt(file, ".hst");
			BEMData::Print("\n- " + Format(t_("Hydrostatic restoring file '%s'"), GetFileName(fileHST)));
			Save_hst(fileHST);
		}
		if (hd().IsLoadedRAO()) {
			String fileRAO = ForceExt(file, ".4");
			BEMData::Print("\n- " + Format(t_("RAO file '%s'"), GetFileName(fileRAO)));
			Save_4(fileRAO, force_T);
		}
		if (hd().IsLoadedQTF()) {
			String fileQTFs = ForceExt(file, ".12s");
			BEMData::Print("\n- " + Format(t_("QTF file '%s'"), GetFileName(fileQTFs)));
			Save_12(fileQTFs, true, force_T, true, qtfHeading);
			String fileQTFd = ForceExt(file, ".12d");
			BEMData::Print("\n- " + Format(t_("QTF file '%s'"), GetFileName(fileQTFd)));
			Save_12(fileQTFd, false, force_T, true, qtfHeading);
		}
	} catch (Exc e) {
		BEMData::PrintError(Format("\n%s: %s", t_("Error"), e));
		hd().lastError = e;
	}
}


bool Wamit::Load_out() {
	hd().Nb = 0;
	hd().Nf = 0;
	hd().Nh = 0;
	int pos;
	int ibody = -1;
	hd().dimen = false;
	
	FileInLine in(hd().file);
	if (!in.IsOpen())
		return false;
	String line;
	FieldSplitWamit f(in);
	f.IsSeparator = IsTabSpace;
	
	hd().names.Clear();
	while(!in.IsEof()) {
		line = in.GetLine();
		f.Load(line);
		if (line.Find("N=") >= 0 && line.Find("Body number:") < 0) {
			hd().Nb++;
			hd().names << GetFileTitle(f.GetText(2));
		} else if ((pos = line.FindAfter("Input from Geometric Data File:")) >= 0) {
			hd().Nb = 1;
			hd().names << GetFileTitle(TrimBoth(line.Mid(pos)));
		} else if (line.Find("POTEN run date and starting time:") >= 0) {
			hd().cg.setConstant(3, hd().Nb, Null);
			hd().cb.setConstant(3, hd().Nb, Null);
			hd().Vo.SetCount(hd().Nb, Null);
			hd().C.SetCount(hd().Nb);
		} else if (line.Find("Gravity:") >= 0) {
			hd().g = f.GetDouble(1);
			hd().len = f.GetDouble(4);
		} else if (line.Find("Water depth:") >= 0) {
			if (ToLower(f.GetText(2)) == "infinite")
				hd().h = -1;
			else {
				hd().h = f.GetDouble(2);
				if (hd().h < 0)
					throw Exc(in.Str() + "\n" +  t_("Water depth has to be positive"));
			}
			if (line.Find("Water density:") >= 0) 
				hd().rho = f.GetDouble(5);			
		} else if (line.Find("XBODY =") >= 0) {
			ibody++;
			if (ibody >= hd().Nb)
				throw Exc(in.Str() + "\n"  + Format(t_("Found additional bodies over %d"), hd().Nb));
			if (hd().cg.rows() < 3 || hd().cg.cols() < hd().Nb)
			 	throw Exc(in.Str() + "\n"  + t_("cg matrix is not dimensioned"));
			hd().cg(0, ibody) = f.GetDouble(2);
			hd().cg(1, ibody) = f.GetDouble(5);
			hd().cg(2, ibody) = f.GetDouble(8);
		} else if ((pos = line.FindAfter("Volumes (VOLX,VOLY,VOLZ):")) >= 0) {
			if (hd().Vo.GetCount() < hd().Nb)
			 	throw Exc(in.Str() + "\n"  + t_("Vo matrix is not dimensioned"));		
			hd().Vo[ibody] = ScanDouble(line.Mid(pos));
		} else if (line.Find("Center of Buoyancy (Xb,Yb,Zb):") >= 0) {
			if (hd().cb.rows() < 3 || hd().cg.cols() < hd().Nb)
			 	throw Exc(in.Str() + "\n"  + t_("cb matrix is not dimensioned"));
			hd().cb(0, ibody) = f.GetDouble(4) + hd().cg(0, ibody);
			hd().cb(1, ibody) = f.GetDouble(5) + hd().cg(1, ibody);
			hd().cb(2, ibody) = f.GetDouble(6) + hd().cg(2, ibody);
		} else if (line.Find("Hydrostatic and gravitational") >= 0) {
			if (hd().C.GetCount() < hd().Nb)
			 	throw Exc(in.Str() + "\n"  + t_("C matrix is not dimensioned"));
			hd().C[ibody].setConstant(6, 6, 0);
			f.LoadWamitJoinedFields(in.GetLine());
			hd().C[ibody](2, 2) = f.GetDouble(1);
			hd().C[ibody](2, 3) = hd().C[ibody](3, 2) = f.GetDouble(2);
			hd().C[ibody](2, 4) = hd().C[ibody](4, 2) = f.GetDouble(3);
			f.LoadWamitJoinedFields(in.GetLine());
			hd().C[ibody](3, 3) = f.GetDouble(1);
			hd().C[ibody](3, 4) = hd().C[ibody](4, 3) = f.GetDouble(2);
			hd().C[ibody](3, 5) = hd().C[ibody](5, 3) = f.GetDouble(3);
			f.LoadWamitJoinedFields(in.GetLine());
			hd().C[ibody](4, 4) = f.GetDouble(1);
			hd().C[ibody](4, 5) = hd().C[ibody](5, 4) = f.GetDouble(2);
		} else if (line.Find("Output from  WAMIT") >= 0) {
			hd().head.Clear();
			FileInLine::Pos fpos = in.GetPos();
			
			bool foundNh = false;
			while (!in.IsEof()) {
				line = in.GetLine();
				if (line.Find("Wave period (sec)") >= 0) {
					++hd().Nf;
					if (hd().Nh > 0 && !foundNh)
						foundNh = true;
				} else if (!foundNh) {
					if (hd().Nh > 0 && (line.Find("*********************") >= 0 ||
								   line.Find("FORCES AND MOMENTS") >= 0)) 
						foundNh = true;
					else if (line.Find("Wave Heading (deg) :") >= 0) {
						f.Load(line);
						hd().head << f.GetDouble(4);
						++hd().Nh;
					}
				}
			}
			if (hd().Nb == 0 || hd().Nh == 0 || hd().Nf == 0)
				throw Exc(in.Str() + "\n"  + Format(t_("Wrong format in Wamit file '%s'"), hd().file));
		
			hd().T.SetCount(hd().Nf);
			hd().w.SetCount(hd().Nf);
			
			in.SeekPos(fpos);
			while (in.GetLine().Find("Wave period = infinite") < 0 && !in.IsEof())
				; 
			if (!in.IsEof()) {
				hd().Aw0.setConstant(hd().Nb*6, hd().Nb*6, Null);
				Load_A(in, hd().Aw0);
			}
			in.SeekPos(fpos);
			while (in.GetLine().Find("Wave period = zero") < 0 && !in.IsEof())
				; 
			if (!in.IsEof()) {
				hd().Awinf.setConstant(hd().Nb*6, hd().Nb*6, Null);
				Load_A(in, hd().Awinf);
			}
			
			in.SeekPos(fpos);
			
			int ifr = -1;
			while (!in.IsEof()) {
				line = in.GetLine();
				while ((line = in.GetLine()).Find("Wave period (sec)") < 0 && !in.IsEof())
					; 
				if (in.IsEof())
					return true;
				
				f.Load(line);
				
				ifr++;
				if (OUTB(ifr, hd().Nf))
					throw Exc(in.Str() + "\n" + Format(t_("Found additional frequencies over %d"), hd().Nf));
				
	            hd().T[ifr] = f.GetDouble(4);  			
	            hd().w[ifr] = fround(2*M_PI/hd().T[ifr], 8);
	            hd().dataFromW = false;
	            
	            bool nextFreq = false;
	            while (!in.IsEof() && !nextFreq) {
	            	line = in.GetLine();
	            	if (line.Find("ADDED-MASS AND DAMPING COEFFICIENTS") >= 0) {
	            		if (hd().A.IsEmpty()) {
							hd().A.SetCount(hd().Nf);
							hd().B.SetCount(hd().Nf);
						}
						in.GetLine(2);
						if (hd().A.GetCount() < hd().Nf)
			 				throw Exc(in.Str() + "\n"  + t_("A matrix is not dimensioned"));
						if (hd().B.GetCount() < hd().Nf)
			 				throw Exc(in.Str() + "\n"  + t_("B matrix is not dimensioned"));
						hd().A[ifr].setConstant(hd().Nb*6, hd().Nb*6, Null);
		            	hd().B[ifr].setConstant(hd().Nb*6, hd().Nb*6, Null);
		            
			            while (!in.IsEof()) {
							line = TrimBoth(in.GetLine());
							if (line.IsEmpty())
			                	break;
							f.Load(line);
							int i = f.GetInt(0) - 1;
							int j = f.GetInt(1) - 1;
							double Aij = f.GetDouble(2);
							double Bij = f.GetDouble(3);
							if (OUTB(i, hd().Nb*6) || OUTB(j, hd().Nb*6))
								throw Exc(in.Str() + "\n"  + Format(t_("Index (%d, %d) out of bounds"), i, j));
							hd().A[ifr](i, j) = Aij;
							hd().B[ifr](i, j) = Bij;
						}
						hd().GetBodyDOF();
	            	} else if (line.Find("DIFFRACTION EXCITING FORCES AND MOMENTS") >= 0) {
						if (hd().ex.ma.IsEmpty()) 
							hd().Initialize_Forces(hd().ex);
						
						int ih = 0;
						while (!in.IsEof()) {		
							line = in.GetLine();
							if (line.Find("Wave Heading (deg) :") >= 0) {
								in.GetLine(3); 
								while (!TrimBoth(line = in.GetLine()).IsEmpty()) {
									f.Load(line);
									double ma = f.GetDouble(1);
									double ph = f.GetDouble(2)*M_PI/180;
									double re = ma*cos(ph);
									double im = ma*sin(ph);
									int i = abs(f.GetInt(0)) - 1;
									if (OUTB(ih, hd().Nh) || OUTB(ifr, hd().Nf) || OUTB(i, hd().Nb*6))
										throw Exc(in.Str() + "\n"  + Format(t_("Index [%d](%d, %d) out of bounds"), ih, ifr, i));
									hd().ex.ma[ih](ifr, i) = ma;	
									hd().ex.ph[ih](ifr, i) = ph;	
									hd().ex.re[ih](ifr, i) = re;	
									hd().ex.im[ih](ifr, i) = im;	
								}
								ih++;
							} else if (line.Find("RESPONSE AMPLITUDE OPERATORS") >= 0 ||
									   line.Find("SURGE, SWAY & YAW DRIFT FORCES") >= 0 ||
									   line.Find("SURGE, SWAY, HEAVE, ROLL, PITCH & YAW DRIFT FORCES") >= 0 ||
									   line.Find("HYDRODYNAMIC PRESSURE IN FLUID DOMAIN") >= 0 ||
									   line.Find("*************************************") >= 0) {
								nextFreq = true;
								break;
							}
						}
					}
				}
			}
		}
	}

	return true;
}


void Wamit::Load_A(FileInLine &in, Eigen::MatrixXd &A) {
	in.GetLine(6);
	while (!in.IsEof()) {
		String line = TrimBoth(in.GetLine());
		if (line.IsEmpty())
           	break;
		FieldSplit f(in);
		f.IsSeparator = IsTabSpace;
		f.Load(line);
		int i = f.GetInt(0) - 1;
		int j = f.GetInt(1) - 1;
		double Aij = f.GetDouble(2);
		if (OUTB(i, A.rows()) || OUTB(j, A.cols()))
			throw Exc(in.Str() + "\n"  + Format(t_("Index (%d, %d) out of bounds"), i, j));
		A(i, j) = Aij;
	}
}

bool Wamit::Load_Scattering(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	FieldSplit f(in);
	f.IsSeparator = IsTabSpace;
	
	hd().Initialize_Forces(hd().sc);
	
	in.GetLine();
    for (int ifr = 0; ifr < hd().Nf; ++ifr) {
        for (int ih = 0; ih < hd().Nh; ++ih) {
            for (int k = 0; k < hd().Nb*6; ++k) {		// Number of non-zero dof
        		f.Load(in.GetLine());
        		
        		int i = f.GetInt(2) - 1;
        		if (OUTB(i, hd().Nb*6))
        			throw Exc(in.Str() + "\n"  + Format(t_("Index (%d) out of bounds"), i));
                hd().sc.ma[ih](ifr, i) = f.GetDouble(3);
                hd().sc.ph[ih](ifr, i) = f.GetDouble(4)*M_PI/180;
                hd().sc.re[ih](ifr, i) = f.GetDouble(5);
                hd().sc.im[ih](ifr, i) = f.GetDouble(6);
            }
        }
    }
    return true;
}
		
bool Wamit::Load_FK(String fileName) {
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	FieldSplit f(in);
	f.IsSeparator = IsTabSpace;
	
	hd().Initialize_Forces(hd().fk);
	
	in.GetLine();	
    for (int ifr = 0; ifr < hd().Nf; ++ifr) {
        for (int ih = 0; ih < hd().Nh; ++ih) {
            for (int k = 0; k < hd().Nb*6; ++k) {		// Number of non-zero dof
                f.Load(in.GetLine());
        		
        		int i = f.GetInt(2) - 1;
        		if (OUTB(i, hd().Nb*6))
        			throw Exc(in.Str() + "\n"  + Format(t_("Index (%d) out of bounds"), i));
                hd().fk.ma[ih](ifr, i) = f.GetDouble(3);
                hd().fk.ph[ih](ifr, i) = f.GetDouble(4)*M_PI/180;
                hd().fk.re[ih](ifr, i) = f.GetDouble(5);
                hd().fk.im[ih](ifr, i) = f.GetDouble(6);
            }
        }
    }	
    return true;
}

bool Wamit::Load_1(String fileName) {
	hd().dimen = false;
	hd().len = 1;
	
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	FieldSplit f(in);
 	f.IsSeparator = IsTabSpace;
 	
 	FileInLine::Pos fpos;
 	while (IsNull(ScanDouble(in.GetLine())) && !in.IsEof())
 		fpos = in.GetPos();
	if (in.IsEof())
		throw Exc("Error in file format");
	
	Vector<double> T; 	
    Vector<double> w; 
    
	in.SeekPos(fpos);
	
	int maxDof = 0;
	bool thereIsAw0 = false, thereIsAwinf = false; 
	while (!in.IsEof()) {
		f.Load(in.GetLine());
		double freq = f.GetDouble(0);
		if (freq < 0)
			thereIsAw0 = true;
		else if (freq == 0)
			thereIsAwinf = true;
		else
			FindAdd(w, freq);
		
		int dof = f.GetInt(1);
		if (dof > maxDof)
			maxDof = dof-1;
	}
	
	int Nb = 1 + int(maxDof/6);
	if (!IsNull(hd().Nb) && hd().Nb < Nb)
		throw Exc(in.Str() + "\n"  + Format(t_("Number of bodies loaded is lower than previous (%d != %d)"), hd().Nb, Nb));
	hd().Nb = Nb;	
	
	int Nf = w.GetCount();
	if (!IsNull(hd().Nf) && hd().Nf != Nf)
		throw Exc(in.Str() + "\n"  + Format(t_("Number of frequencies loaded is different than previous (%d != %d)"), hd().Nf, Nf));
	hd().Nf = Nf;
	
	if (hd().Nb == 0 || hd().Nf < 2)
		throw Exc(in.Str() + "\n"  + Format(t_("Wrong format in Wamit file '%s'"), hd().file));
	
	if (w[0] > w[1]) {
		hd().dataFromW = false;
		T = pick(w);
		w.SetCount(hd().Nf);	
	} else {
		hd().dataFromW = true;
		T.SetCount(hd().Nf);
	}
	
	hd().A.SetCount(hd().Nf);
	hd().B.SetCount(hd().Nf);	
	if (thereIsAw0)
		hd().Aw0.setConstant(hd().Nb*6, hd().Nb*6, Null);
	if (thereIsAwinf)
		hd().Awinf.setConstant(hd().Nb*6, hd().Nb*6, Null);

	for (int ifr = 0; ifr < hd().Nf; ++ifr) {
		if (hd().dataFromW)
			T[ifr] = fround(2*M_PI/w[ifr], 8);
		else
			w[ifr] = fround(2*M_PI/T[ifr], 8);
		hd().A[ifr].setConstant(hd().Nb*6, hd().Nb*6, Null);
	  	hd().B[ifr].setConstant(hd().Nb*6, hd().Nb*6, Null);
	}
	hd().names.SetCount(Nb);
	
	if (hd().w.IsEmpty()) {
		hd().w = pick(w);
		hd().T = pick(T);
	} else if (!CompareRatio(hd().w, w, 0.001))
		throw Exc(in.Str() + "\n"  + Format(t_("Frequencies loaded are different than previous\nPrevious: %s\nSeries:   %s"), ToString(hd().w), ToString(w)));
	else if (!CompareRatio(hd().T, T, 0.001))
		throw Exc(in.Str() + "\n"  + Format(t_("Periods loaded are different than previous\nPrevious: %s\nSeries:   %s"), ToString(hd().T), ToString(T)));
				
	in.SeekPos(fpos);
	
	while (!in.IsEof()) {
		f.Load(in.GetLine());
		double freq = f.GetDouble(0);
 		int i = f.GetInt(1) - 1;
 		int j = f.GetInt(2) - 1;
 		if (i >= Nb*6 || i < 0 || j >= Nb*6 || j < 0)
			throw Exc(in.Str() + "\n"  + Format(t_("DOF # does not match (%d, %d)"), i+1, j+1));
 		
 		double Aij = f.GetDouble(3);
 		
 		if ((freq < 0)) {
 			if (!thereIsAw0)
				throw Exc(in.Str() + "\n"  + t_("A[w=inf] is not expected"));
			hd().Aw0(i, j) = Aij;
		} else if (freq == 0) {
			if (!thereIsAwinf)
				throw Exc(in.Str() + "\n"  + t_("A[w=0] is not expected"));				
			hd().Awinf(i, j) = Aij;
		} else {
			int ifr;
			if (hd().dataFromW)
				ifr = FindRatio(hd().w, freq, 0.001);
			else
				ifr = FindRatio(hd().T, freq, 0.001);
			if (ifr < 0)
				throw Exc(in.Str() + "\n"  + Format(t_("Frequency %f is unknown"), freq));
		
		  	hd().A[ifr](i, j) = Aij;    
		  	hd().B[ifr](i, j) = f.GetDouble(4);   	
		}
	}
	return true;	
}
 
bool Wamit::Load_3(String fileName) {
	hd().dimen = false;
	hd().len = 1;
	
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	FieldSplit f(in);
 	f.IsSeparator = IsTabSpace;
 
 	FileInLine::Pos fpos;
 	while (IsNull(ScanDouble(in.GetLine())) && !in.IsEof())
 		fpos = in.GetPos();
	if (in.IsEof())
		throw Exc(t_("Error in file format"));
	
	Vector<double> T; 	
    Vector<double> w;
    
	in.SeekPos(fpos);
	
	int maxDof = 0;	
	hd().head.Clear();
	while (!in.IsEof()) {
		f.Load(in.GetLine());
		double freq = f.GetDouble(0);
		double head = f.GetDouble(1);
		FindAdd(w, freq);
		FindAdd(hd().head, head);
		
		int dof = f.GetInt(2);
		if (dof > maxDof)
			maxDof = dof-1;
	}
	
	int Nb = 1 + int(maxDof/6);
	if (!IsNull(hd().Nb) && hd().Nb < Nb)
		throw Exc(in.Str() + "\n"  + Format(t_("Number of bodies loaded is lower than previous (%d != %d)"), hd().Nb, Nb));
	hd().Nb = Nb;
	if (hd().names.IsEmpty())
		hd().names.SetCount(hd().Nb);
		
	if (hd().head.GetCount() == 0)
		throw Exc(in.Str() + "\n" + Format(t_("Wrong format in Wamit file '%s'"), hd().file));
	
	if (!IsNull(hd().Nh) && hd().Nh != hd().head.GetCount())
		throw Exc(in.Str() + "\n" + "\n"  + Format(t_("Number of headings is different than previous (%d != %d)"), hd().Nh, hd().head.GetCount()));
	hd().Nh = hd().head.GetCount();
	
	int Nf = w.GetCount();
	if (!IsNull(hd().Nf) && hd().Nf != Nf)
		throw Exc(in.Str() + "\n"  + Format(t_("Number of frequencies loaded is different than previous (%d != %d)"), hd().Nf, Nf));
	hd().Nf = Nf;
			
	if (w[0] > w[1]) {
		hd().dataFromW = false;
		T = pick(w);
		w.SetCount(hd().Nf);	
	} else {
		hd().dataFromW = true;
		T.SetCount(hd().Nf);
	}
	
	for (int ifr = 0; ifr < hd().Nf; ++ifr) {
		if (hd().dataFromW)
			T[ifr] = fround(2*M_PI/w[ifr], 8);
		else
			w[ifr] = fround(2*M_PI/T[ifr], 8);
	}

	if (hd().w.IsEmpty()) {
		hd().w = pick(w);
		hd().T = pick(T);
	} else if (!CompareRatio(hd().w, w, 0.001))
		throw Exc(in.Str() + "\n"  + Format(t_("Frequencies loaded are different than previous\nPrevious: %s\nSeries:   %s"), ToString(hd().w), ToString(w)));
	else if (!CompareRatio(hd().T, T, 0.001))
		throw Exc(in.Str() + "\n"  + Format(t_("Periods loaded are different than previous\nPrevious: %s\nSeries:   %s"), ToString(hd().T), ToString(T)));
	
	in.SeekPos(fpos);
	
	hd().Initialize_Forces(hd().ex);
	
	while (!in.IsEof()) {
		f.Load(in.GetLine());
		double freq = f.GetDouble(0);
		int ifr;
		if (hd().dataFromW)
		 	ifr = FindRatio(hd().w, freq, 0.001);
		else
			ifr = FindRatio(hd().T, freq, 0.001);
		if (ifr < 0)
			throw Exc(in.Str() + "\n"  + Format(t_("Frequency %f is unknown"), freq));
		double head = f.GetDouble(1);
		int ih = FindRatio(hd().head, head, 0.001);
		if (ih < 0)
			throw Exc(in.Str() + "\n"  + Format(t_("Heading %f is unknown"), head));
			
		int i = f.GetInt(2) - 1;		
		
       	hd().ex.ma[ih](ifr, i) = f.GetDouble(3);
     	hd().ex.ph[ih](ifr, i) = f.GetDouble(4)*M_PI/180;
        hd().ex.re[ih](ifr, i) = f.GetDouble(5);
        hd().ex.im[ih](ifr, i) = f.GetDouble(6);
	}
		
	return true;
}

bool Wamit::Load_hst(String fileName) {
	hd().dimen = false;
	if (IsNull(hd().len))
		hd().len = 1;
	
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	FieldSplit f(in);
	f.IsSeparator = IsTabSpace;
 
 	FileInLine::Pos fpos;
 	while (IsNull(ScanDouble(in.GetLine())) && !in.IsEof())
 		fpos = in.GetPos();
	if (in.IsEof())
		throw Exc(t_("Error in file format"));
	
	in.SeekPos(fpos);
	
	hd().C.SetCount(hd().Nb);
	for(int ib = 0; ib < hd().Nb; ++ib)
		hd().C[ib].setConstant(6, 6, 0);

	while (!in.IsEof()) {
		f.Load(in.GetLine());	
		int i = f.GetInt(0) - 1;
		int ib_i = i/6;
		i = i - ib_i*6;
		int j = f.GetInt(1) - 1;
		int ib_j = j/6;
		j = j - ib_j*6;
		if (ib_i == ib_j) 
			hd().C[ib_i](i, j) = f.GetDouble(2);
	}
		
	return true;
}

bool Wamit::Load_4(String fileName) {
	hd().dimen = false;
	if (IsNull(hd().len))
		hd().len = 1;
	
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	FieldSplit f(in);
	f.IsSeparator = IsTabSpace;
 
 	FileInLine::Pos fpos;
 	while (IsNull(ScanDouble(in.GetLine())) && !in.IsEof())
 		fpos = in.GetPos();
	if (in.IsEof())
		throw Exc("Error in file format");
	
	Vector<double> T; 	
    Vector<double> w;
    
	in.SeekPos(fpos);
	
	int maxDof = 0;
	hd().head.Clear();
	while (!in.IsEof()) {
		f.Load(in.GetLine());
		double freq = f.GetDouble(0);
		FindAdd(w, freq);
		
		int dof = f.GetInt(2);
		if (dof > maxDof)
			maxDof = dof-1;
		
		double head = f.GetDouble(1);
		
		FindAdd(hd().head, head);
	}

	if (hd().head.GetCount() == 0)
		throw Exc(in.Str() + "\n"  + Format(t_("Wrong format in Wamit file '%s'"), hd().file));
	
	if (!IsNull(hd().Nh) && hd().Nh != hd().head.GetCount())
		throw Exc(in.Str() + "\n"  + Format(t_("Number of headings loaded is different than previous (%d != %d)"), hd().Nh, hd().head.GetCount()));
	hd().Nh = hd().head.GetCount();
	
	int Nb = 1 + int(maxDof/6);
	if (!IsNull(hd().Nb) && hd().Nb < Nb)
		throw Exc(in.Str() + "\n"  + Format(t_("Number of bodies loaded is lower than previous (%d != %d)"), hd().Nb, Nb));
	hd().Nb = Nb;
	if (hd().names.IsEmpty())
		hd().names.SetCount(hd().Nb);
		
	int Nf = w.GetCount();
	if (!IsNull(hd().Nf) && hd().Nf != Nf)
		throw Exc(in.Str() + "\n"  + Format(t_("Number of frequencies loaded is different than previous (%d != %d)"), hd().Nf, Nf));
	hd().Nf = Nf;
	
	if (hd().Nb == 0 || hd().Nf < 2)
		throw Exc(in.Str() + "\n"  + Format(t_("Wrong format in Wamit file '%s'"), hd().file));
	
	if (w[0] > w[1]) {
		hd().dataFromW = false;
		T = pick(w);
		w.SetCount(hd().Nf);	
	} else {
		hd().dataFromW = true;
		T.SetCount(hd().Nf);
	}
	
	for (int ifr = 0; ifr < hd().Nf; ++ifr) {
		if (hd().dataFromW)
			T[ifr] = fround(2*M_PI/w[ifr], 8);
		else
			w[ifr] = fround(2*M_PI/T[ifr], 8);
	}
	
	if (hd().w.IsEmpty()) {
		hd().w = pick(w);
		hd().T = pick(T);
	} else if (!CompareRatio(hd().w, w, 0.01))
		throw Exc(in.Str() + "\n"  + Format(t_("Frequencies loaded are different than previous\nPrevious: %s\nSeries:   %s"), ToString(hd().w), ToString(w)));
	else if (!CompareRatio(hd().T, T, 0.01))
		throw Exc(Format(t_("[%s] Periods loaded are different than previous\nPrevious: %s\nSeries:   %s"), ToString(hd().T), ToString(T)));
	
	if (hd().names.IsEmpty())
		hd().names.SetCount(hd().Nb);
	
	in.SeekPos(fpos);
	
	hd().Initialize_RAO();
		
	while (!in.IsEof()) {
		f.Load(in.GetLine());
		double freq = f.GetDouble(0);
		int ifr;
		if (hd().dataFromW)
		 	ifr = FindRatio(hd().w, freq, 0.01);
		else
			ifr = FindRatio(hd().T, freq, 0.01);
		if (ifr < 0)
			throw Exc(in.Str() + "\n"  + Format(t_("Frequency %f is unknown"), freq));
		double head = f.GetDouble(1);
		int ih = FindRatio(hd().head, head, 0.001);
		if (ih < 0)
			throw Exc(in.Str() + "\n"  + Format(t_("Heading %f is unknown"), head));
			
		int i = f.GetInt(2) - 1;		
		
       	hd().rao.ma[ih](ifr, i) = f.GetDouble(3);
     	hd().rao.ph[ih](ifr, i) = f.GetDouble(4)*M_PI/180;
        hd().rao.re[ih](ifr, i) = f.GetDouble(5);
        hd().rao.im[ih](ifr, i) = f.GetDouble(6);
	}
		
	return true;
}

bool Wamit::Load_12(String fileName, bool isSum) {
	hd().dimen = false;
	if (IsNull(hd().len))
		hd().len = 1;
	
	FileInLine in(fileName);
	if (!in.IsOpen())
		return false;
	
	int nrows = hd().Nh*hd().Nf*hd().Nf;
	Upp::Array<Hydro::QTF> &qtfList = isSum ? hd().qtfsum : hd().qtfdif;
	
	qtfList.Clear();
	qtfList.Reserve(hd().Nb*nrows);
		
	FieldSplit f(in);
	f.IsSeparator = IsTabSpace;
	
	in.GetLine();		// Avoid header
	
	FileInLine::Pos fpos = in.GetPos();
	
    Vector<double> w, T;
    Vector<double> head;
    
    bool isRad = true;
	int maxDof = 0;	
	while (!in.IsEof()) {
		f.Load(in.GetLine());
		double freq = f.GetDouble(0);
		double hd = f.GetDouble(2);
		FindAdd(w, freq);
		FindAdd(head, hd);
		
		int dof = f.GetInt(4);
		if (dof > maxDof)
			maxDof = dof-1;
		
		if (abs(f.GetDouble(6)) > M_PI + 0.01)
			isRad = false;
	}
	
	int Nb = 1 + int(maxDof/6);
	if (!IsNull(hd().Nb) && hd().Nb < Nb)
		throw Exc(Format(t_("Number of bodies loaded is lower than previous (%d != %d)"), hd().Nb, Nb));
	hd().Nb = Nb;	
	if (hd().names.IsEmpty())
		hd().names.SetCount(hd().Nb);
		
	int Nh = head.GetCount();
	if (Nh == 0)
		throw Exc(Format(t_("Wrong format in Wamit file '%s'. No headings found"), hd().file));
	if (hd().qtfhead.IsEmpty()) 
		hd().qtfhead = pick(head);
	else {
		for(int i = 0; i < head.GetCount(); i++) {
			bool found = false;
			for(int j = 0; j < hd().qtfhead.GetCount(); j++) {
				if (EqualRatio(head[i], hd().qtfhead[j], 0.001, 0.0001)) {
					found = true;
					break;
				}
			}
			if (!found)
				throw Exc(Format(t_("Heading %f not found in previous %s"), head[i], ToString(hd().head)));
		}
	}

	int Nf = w.GetCount();
	
	if (w[0] > w[1]) {
		hd().qtfdataFromW = false;
		T = pick(w);
		w.SetCount(Nf);	
	} else {
		hd().qtfdataFromW = true;
		T.SetCount(Nf);
	}
	
	for (int ifr = 0; ifr < Nf; ++ifr) {
		if (hd().qtfdataFromW)
			T[ifr] = fround(2*M_PI/w[ifr], 8);
		else
			w[ifr] = fround(2*M_PI/T[ifr], 8);
	}
	
	if (hd().qtfw.IsEmpty()) {
		hd().qtfw = pick(w);
		hd().qtfT = pick(T);
	} else if (!CompareRatio(hd().qtfw, w, 0.01))
		throw Exc(in.Str() + "\n"  + Format(t_("Frequencies loaded are different than previous\nPrevious: %s\nSeries:   %s"), ToString(hd().w), ToString(w)));
	else if (!CompareRatio(hd().qtfT, T, 0.01))
		throw Exc(Format(t_("[%s] Periods loaded are different than previous\nPrevious: %s\nSeries:   %s"), ToString(hd().T), ToString(T)));
	
	in.SeekPos(fpos);
	
	while (!in.IsEof()) {
		f.Load(in.GetLine());
	
		double wT1 = f.GetDouble(0);
		int iwT1;
		if (hd().qtfdataFromW) {
			iwT1 = FindClosest(hd().qtfw, wT1);
			if (!EqualRatio(hd().qtfw[iwT1], wT1, 0.01, 0.001))
				throw Exc(in.Str() + "\n"  + Format(t_("Frequency 1 %f not found"), wT1));
		} else {
			iwT1 = FindClosest(hd().qtfT, wT1);
			if (!EqualRatio(hd().qtfT[iwT1], wT1, 0.01, 0.001))
				throw Exc(in.Str() + "\n"  + Format(t_("Period 1 %f not found"), wT1));
		}
		double wT2 = f.GetDouble(1);
		int iwT2;
		if (hd().qtfdataFromW) {
			iwT2 = FindClosest(hd().qtfw, wT2);
			if (!EqualRatio(hd().qtfw[iwT2], wT2, 0.01, 0.001))
				throw Exc(in.Str() + "\n"  + Format(t_("Frequency 2 %f not found"), wT2));
		} else {
			iwT2 = FindClosest(hd().qtfT, wT2);
			if (!EqualRatio(hd().qtfT[iwT2], wT2, 0.01, 0.001))
				throw Exc(in.Str() + "\n"  + Format(t_("Period 2 %f not found"), wT2));
		}
		double h1 = f.GetDouble(2);
		int ih1 = FindClosest(hd().qtfhead, h1);
		if (!EqualRatio(hd().qtfhead[ih1], h1, 0.01, 0.001))
			throw Exc(in.Str() + "\n"  + Format(t_("Heading 1 %f not found"), h1));
		double h2 = f.GetDouble(3);
		int ih2 = FindClosest(hd().qtfhead, h2);
		if (!EqualRatio(hd().qtfhead[ih2], h2, 0.01, 0.001))
			throw Exc(in.Str() + "\n"  + Format(t_("Heading 2 %f not found"), h2));

	    int idof = f.GetInt(4)-1;
	    if (idof < 0 || idof > hd().Nb*6 -1)
	        throw Exc(in.Str() + "\n"  + Format(t_("Wrong DOF id %d"), idof+1));
	        		
		int ib = int(idof/6);
		idof = idof - ib*6;
		
		Hydro::QTF *pqtf;
		int id = hd().GetQTFId(qtfList, ib, ih1, ih2, iwT1, iwT2);
		if (id >= 0)
			pqtf = &qtfList[id];
		else {
			pqtf = &qtfList.Add();
			pqtf->Set(ib, ih1, ih2, iwT1, iwT2);
		}
	    
	    Hydro::QTF &qtf = *pqtf;
	    
	    qtf.fma[idof] = f.GetDouble(5);
	    double ph = f.GetDouble(6);
	    qtf.fph[idof] = isRad ? ph : ToRad(ph);
	    qtf.fre[idof] = f.GetDouble(7);
	    qtf.fim[idof] = f.GetDouble(8);    
	        
		double fre = qtf.fma[idof]*cos(qtf.fph[idof]);
		double fim = qtf.fma[idof]*sin(qtf.fph[idof]);
		
		if (abs(fre - qtf.fre[idof]) > 0.1)  
			throw Exc(in.Str() + "\n"  + Format(t_("Real force %f does not match with magnitude %f and phase %f (%f)"), 
										qtf.fre[idof], qtf.fma[idof], ph, fre));
		if (abs(fim - qtf.fim[idof]) > 0.1)  
			throw Exc(in.Str() + "\n"  + Format(t_("Imaginary force %f does not match with magnitude %f and phase %f (%f)"), 
										qtf.fim[idof], qtf.fma[idof], ph, fim));
	}
	return true;
}

void Wamit::Save_1(String fileName, bool force_T) {
	if (!(hd().IsLoadedA() && hd().IsLoadedB())) 
		return;
		
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'"), fileName));
	
	if (hd().IsLoadedAw0()) {
		for (int i = 0; i < hd().Nb*6; ++i)  
			for (int j = 0; j < hd().Nb*6; ++j)
				if (!IsNull(hd().Aw0(i, j))) 
					out << Format(" %s %5d %5d %s\n", FormatWam(-1), i+1, j+1,
													  FormatWam(hd().Aw0_ndim(i, j)));
	}
	if (hd().IsLoadedAwinf()) {
		for (int i = 0; i < hd().Nb*6; ++i)  
			for (int j = 0; j < hd().Nb*6; ++j)
				if (!IsNull(hd().Awinf(i, j))) 
					out << Format(" %s %5d %5d %s\n", FormatWam(0), i+1, j+1,
													  FormatWam(hd().Awinf_ndim(i, j)));
	}
	
	if (hd().Nf < 2)
		throw Exc(t_("No enough data to save (at least 2 frequencies)"));
		
	Vector<double> *pdata;
	if (force_T)
		pdata = &hd().T;
	else if (hd().dataFromW) 
		pdata = &hd().w;
	else
		pdata = &hd().T;
	Vector<double> &data = *pdata;
	
	int ifr0, ifrEnd, ifrDelta;
	bool growing = data[1] > data[0];
	if ((growing && (hd().dataFromW && !force_T)) || (!growing && !(hd().dataFromW && !force_T))) {
		ifr0 = 0;
		ifrEnd = hd().Nf;
		ifrDelta = 1;
	} else {
		ifr0 = hd().Nf - 1;
		ifrEnd = -1;
		ifrDelta = -1;
	}
	
	for (int ifr = ifr0; ifr != ifrEnd; ifr += ifrDelta)
		for (int i = 0; i < hd().Nb*6; ++i)  
			for (int j = 0; j < hd().Nb*6; ++j)
				if (!IsNull(hd().A[ifr](i, j)) && !IsNull(hd().B[ifr](i, j))) 
					out << Format(" %s %5d %5d %s %s\n", FormatWam(data[ifr]), i+1, j+1,
														 FormatWam(hd().A_ndim(ifr, i, j)), 
														 FormatWam(hd().B_ndim(ifr, i, j)));
}

void Wamit::Save_3(String fileName, bool force_T) {
	if (!hd().IsLoadedFex()) 
		return;
	
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'"), fileName));

	if (hd().Nf < 2)
		throw Exc(t_("No enough data to save (at least 2 frequencies)"));

	Vector<double> *pdata;
	int ifr0, ifrEnd, ifrDelta;
	if (hd().dataFromW && !force_T) 
		pdata = &hd().w;
	else
		pdata = &hd().T;
	Vector<double> &data = *pdata;
	
	if (((data[1] > data[0]) && (hd().dataFromW && !force_T)) || ((data[1] < data[0]) && !(hd().dataFromW && !force_T))) {
		ifr0 = 0;
		ifrEnd = hd().Nf;
		ifrDelta = 1;
	} else {
		ifr0 = hd().Nf - 1;
		ifrEnd = -1;
		ifrDelta = -1;
	}
	
	for (int ifr = ifr0; ifr != ifrEnd; ifr += ifrDelta)
		for (int ih = 0; ih < hd().Nh; ++ih)
			for (int i = 0; i < hd().Nb*6; ++i)
				if (!IsNull(hd().ex.ma[ih](ifr, i))) 
					out << Format(" %s %s %5d %s %s %s %s\n", 
									FormatWam(data[ifr]), FormatWam(hd().head[ih]), i+1,
									FormatWam(hd().F_ma_ndim(hd().ex, ih, ifr, i)), 
									FormatWam(hd().ex.ph[ih](ifr, i)*180/M_PI),
									FormatWam(hd().F_re_ndim(hd().ex, ih, ifr, i)), 
									FormatWam(hd().F_im_ndim(hd().ex, ih, ifr, i)));
}

void Wamit::Save_hst(String fileName) {
	if (!hd().IsLoadedC()) 
		return;
		
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'"), fileName));

	for (int i = 0; i < 6*hd().Nb; ++i)  
		for (int j = 0; j < 6*hd().Nb; ++j) {
			int ib_i = i/6;
			int ii = i - ib_i*6;
			int ib_j = j/6;
			int jj = j - ib_j*6;
			if (!IsNull(hd().C[ib_i](ii, jj))) 
				out << Format(" %5d %5d  %s\n", i+1, j+1, 
								FormatWam(hd().C_ndim(ib_i, ii, jj)));
		}
}

// K is supposed to be dimensionalized
void Wamit::Save_hst_static(Eigen::MatrixXd C, String fileName, double rho, double g) {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'"), fileName));

	for (int i = 0; i < 6; ++i)  
		for (int j = 0; j < 6; ++j) {
			if (!IsNull(C(i, j))) 
				out << Format(" %5d %5d  %s\n", i+1, j+1, 
								FormatWam(C(i, j)/rho/g));
		}
}

void Wamit::Save_4(String fileName, bool force_T) {
	if (!hd().IsLoadedRAO()) 
		return;
		
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'"), fileName));
	
	if (hd().Nf < 2)
		throw Exc(t_("No enough data to save (at least 2 frequencies)"));
		
	Vector<double> *pdata;
	int ifr0, ifrEnd, ifrDelta;
	if (hd().dataFromW && !force_T) 
		pdata = &hd().w;
	else
		pdata = &hd().T;
	Vector<double> &data = *pdata;
	
	if (((data[1] > data[0]) && (hd().dataFromW && !force_T)) || ((data[1] < data[0]) && !(hd().dataFromW && !force_T))) {
		ifr0 = 0;
		ifrEnd = hd().Nf;
		ifrDelta = 1;
	} else {
		ifr0 = hd().Nf - 1;
		ifrEnd = -1;
		ifrDelta = -1;
	}
	
	for (int ifr = ifr0; ifr != ifrEnd; ifr += ifrDelta)
		for (int ih = 0; ih < hd().Nh; ++ih)
			for (int i = 0; i < hd().Nb*6; ++i)
				if (!IsNull(hd().rao.ma[ih](ifr, i))) 
					out << Format(" %s %s %5d %s %s %s %s\n", 
									FormatWam(data[ifr]), FormatWam(hd().head[ih]), i+1,
									FormatWam(hd().R_ma_ndim(hd().rao, ih, ifr, i)), 
									FormatWam(hd().rao.ph[ih](ifr, i)*180/M_PI),
									FormatWam(hd().R_re_ndim(hd().rao, ih, ifr, i)), 
									FormatWam(hd().R_im_ndim(hd().rao, ih, ifr, i)));
}
	
void Wamit::Save_12(String fileName, bool isSum, bool force_T, bool force_Deg, int qtfHeading) {
	if (!hd().IsLoadedQTF()) 
		return;
	
	if (isSum)
		fileName = ForceExt(fileName, ".12s");
	else
		fileName = ForceExt(fileName, ".12d");
	
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format(t_("Impossible to open '%s'"), fileName));
	
	int Nf = hd().qtfw.GetCount();
	int Nh = hd().qtfhead.GetCount(); 
	
	if (Nf < 2)
		throw Exc(t_("No enough data to save (at least 2 frequencies)"));
			
	Upp::Array<Hydro::QTF> &qtfList = isSum ? hd().qtfsum : hd().qtfdif;
		
	out << " WAMIT Numeric Output -- Filename  " << Format("%20<s", GetFileName(fileName)) << "  " << Format("%", GetSysTime()) << "\n";
	
	Vector<double> *pdata;
	int ifr0, ifrEnd, ifrDelta;
	if (hd().qtfdataFromW && !force_T) 
		pdata = &hd().qtfw;
	else
		pdata = &hd().qtfT;
	Vector<double> &data = *pdata;
	
	if (((data[1] > data[0]) && (hd().dataFromW && !force_T)) || ((data[1] < data[0]) && !(hd().dataFromW && !force_T))) {
		ifr0 = 0;
		ifrEnd = hd().qtfw.GetCount();
		ifrDelta = 1;
	} else {
		ifr0 = hd().qtfw.GetCount() - 1;
		ifrEnd = -1;
		ifrDelta = -1;
	}
	
	static int idf12[] = {1, 3, 5, 2, 4, 6};
	
	for (int ifr1 = ifr0; ifr1 != ifrEnd; ifr1 += ifrDelta) 
		for (int ifr2 = ifr0; ifr2 != ifrEnd; ifr2 += ifrDelta) 
			for (int ih1 = 0; ih1 < Nh; ++ih1) {
				if (!IsNull(qtfHeading) && ih1 != qtfHeading)
					continue;
				for (int ih2 = 0; ih2 < Nh; ++ih2) {
					if (!IsNull(qtfHeading) && ih1 != qtfHeading)
						continue;	 
					for (int ib = 0; ib < hd().Nb; ++ib) {
						int id = hd().GetQTFId(qtfList, ib, ih1, ih2, ifr1, ifr2);
						if (id >= 0) {
							double qtfhead1, qtfhead2;
							if (!IsNull(qtfHeading)) 
								qtfhead1 = qtfhead2 = 0;		// Just 0º
							else {
								qtfhead1 = hd().qtfhead[ih1];
								qtfhead2 = hd().qtfhead[ih2];
							}							
							for (int _idf = 0; _idf < 6; ++_idf) {
								int idf = idf12[_idf]-1;
								Hydro::QTF &qtf = qtfList[id];
								out << Format("   %s   %s   %s   %s   %2d   %s   %s   %s   %s\n", 
										FormatWam(data[ifr1]),
										FormatWam(data[ifr2]), 
										FormatWam(qtfhead1),
										FormatWam(qtfhead2),
										ib*6 + idf + 1,
										FormatWam(hd().F_ndim(qtf.fma[idf], idf)), 
										FormatWam(!force_Deg ? qtf.fph[idf] : ToDeg(qtf.fph[idf])),
										FormatWam(hd().F_ndim(qtf.fre[idf], idf)), 
										FormatWam(hd().F_ndim(qtf.fim[idf], idf)));
							}
						}
					}
				}
			}
}
