#include "BEMRosetta.h"


bool Wamit::Load(String file, double rho) {
	code = WAMIT;
	this->file = file;	
	this->name = GetFileTitle(file);
	
	this->rho = rho;
	
	try {
		Print("\n\n" + Format("Loading '%s'", file));
		if (!Load_out()) {
			Print("\n" + Format("File '%s' not found", file));
			return false;
		}
		String fileSC = ForceExt(file, ".3sc");
		Print("\n" + Format("- Scattering file '%s'", GetFileName(fileSC)));
		if (!Load_Scattering(fileSC))
			Print(": **Not found**");
		String fileFK = ForceExt(file, ".3fk");
		Print("\n" + Format("- Froude-Krylov file '%s'", GetFileName(fileFK)));
		if (!Load_FK(fileFK))
			Print(": **Not found**");
		
		AfterLoad();
	} catch (Exc e) {
		PrintError("\nError: " + e);
		lastError = e;
		return false;
	}
	
	return true;
}

bool Wamit::Load_out() {
	Nb = 0;
	Nf = 0;
	Nh = 0;
	int pos;
	int ibody = -1;
	
	FileIn in(file);
	if (!in.IsOpen())
		return false;
	String line;
	FieldSplit f;
	names.Clear();
	while(!in.IsEof()) {
		line = in.GetLine();
		f.Load(line);
		if (line.Find("N=") >= 0 && line.Find("Body number:") < 0) {
			Nb++;
			names << GetFileTitle(f.GetText(2));
		} else if ((pos = line.FindAfter("Input from Geometric Data File:")) >= 0) {
			Nb = 1;
			names << GetFileTitle(TrimBoth(line.Mid(pos)));
		} else if (line.Find("POTEN run date and starting time:") >= 0) {
			cg.setConstant(3, Nb, nan(""));
			cb.setConstant(3, Nb, nan(""));
			Vo.SetCount(Nb, nan(""));
			C.SetCount(Nb);
		} else if (line.Find("Gravity:") >= 0) {
			g = f.GetDouble(1);
			len = f.GetDouble(4);
		} else if (line.Find("Water depth:") >= 0) {
			if (ToLower(f.GetText(2)) == "infinite")
				h = INFINITY;
			else
				h = f.GetDouble(2);
			if (line.Find("Water density:") >= 0) 
				rho = f.GetDouble(5);			
		} else if (line.Find("XBODY =") >= 0) {
			ibody++;
			cg(0, ibody) = f.GetDouble(2);
			cg(1, ibody) = f.GetDouble(5);
			cg(2, ibody) = f.GetDouble(8);
		} else if ((pos = line.FindAfter("Volumes (VOLX,VOLY,VOLZ):")) >= 0) 
			Vo[ibody] = ScanDouble(line.Mid(pos));
		else if (line.Find("Center of Buoyancy (Xb,Yb,Zb):") >= 0) {
			cb(0, ibody) = f.GetDouble(4) + cg(0, ibody);
			cb(1, ibody) = f.GetDouble(5) + cg(1, ibody);
			cb(2, ibody) = f.GetDouble(6) + cg(2, ibody);
		} else if (line.Find("Hydrostatic and gravitational") >= 0) {
			C[ibody].setConstant(6, 6, 0);
			f.LoadWamitJoinedFields(in.GetLine());
			C[ibody](2, 2) = f.GetDouble(1);
			C[ibody](2, 3) = C[ibody](3, 2) = f.GetDouble(2);
			C[ibody](2, 4) = C[ibody](4, 2) = f.GetDouble(3);
			f.LoadWamitJoinedFields(in.GetLine());
			C[ibody](3, 3) = f.GetDouble(1);
			C[ibody](3, 4) = C[ibody](4, 3) = f.GetDouble(2);
			C[ibody](3, 5) = C[ibody](5, 3) = f.GetDouble(3);
			f.LoadWamitJoinedFields(in.GetLine());
			C[ibody](4, 4) = f.GetDouble(1);
			C[ibody](4, 5) = C[ibody](5, 4) = f.GetDouble(2);
		} else if (line.Find("Output from  WAMIT") >= 0) {
			head.Clear();
			int64 fpos = in.GetPos();
			
			bool foundNh = false;
			while (!in.IsEof()) {
				line = in.GetLine();
				if (line.Find("Wave period (sec)") >= 0) {
					++Nf;
					if (Nh > 0 && !foundNh)
						foundNh = true;
				} else if (!foundNh) {
					if (Nh > 0 && (line.Find("*********************") >= 0 ||
								   line.Find("FORCES AND MOMENTS") >= 0)) 
						foundNh = true;
					else if (line.Find("Wave Heading (deg) :") >= 0) {
						f.Load(line);
						head << f.GetDouble(4);
						++Nh;
					}
				}
			}
			if (Nb == 0 || Nh == 0 || Nf == 0)
				throw Exc(Format("Wrong format in Wamit file '%s'", file));
		
			T.SetCount(Nf);
			w.SetCount(Nf);
			
			in.Seek(fpos);
			while (in.GetLine().Find("Wave period = infinite") < 0 && !in.IsEof())
				; 
			if (!in.IsEof()) {
				Aw0.setConstant(Nb*6, Nb*6, nan(""));
				Load_A(in, Aw0);
			}
			in.Seek(fpos);
			while (in.GetLine().Find("Wave period = zero") < 0 && !in.IsEof())
				; 
			if (!in.IsEof()) {
				Awinf.setConstant(Nb*6, Nb*6, nan(""));
				Load_A(in, Awinf);
			}
			
			in.Seek(fpos);
			
			int ifr = -1;
			while (!in.IsEof()) {
				line = in.GetLine();
				while ((line = in.GetLine()).Find("Wave period (sec)") < 0 && !in.IsEof())
					; 
				if (in.IsEof())
					return true;
				
				f.Load(line);
				
				ifr++;
	            T[ifr] = f.GetDouble(4);  			// Wave periods
	            w[ifr] = 2*M_PI/T[ifr];  			// Wave frequencies
	            
	            bool nextFreq = false;
	            while (!in.IsEof() && !nextFreq) {
	            	line = in.GetLine();
	            	if (line.Find("ADDED-MASS AND DAMPING COEFFICIENTS") >= 0) {
	            		if (A.IsEmpty()) {
							A.SetCount(Nf);
							B.SetCount(Nf);
						}
						in.GetLine();	in.GetLine();
						A[ifr].setConstant(Nb*6, Nb*6, nan(""));
		            	B[ifr].setConstant(Nb*6, Nb*6, nan(""));
		            
			            while (!in.IsEof()) {
							line = TrimBoth(in.GetLine());
							if (line.IsEmpty())
			                	break;
							f.Load(line);
							int i = f.GetInt(0) - 1;
							int j = f.GetInt(1) - 1;
							double Aij = f.GetDouble(2);
							double Bij = f.GetDouble(3);
							A[ifr](i, j) = Aij;
							B[ifr](i, j) = Bij;
						}
						GetBodyDOF();
	            	} else if (line.Find("DIFFRACTION EXCITING FORCES AND MOMENTS") >= 0) {
						if (ex.ma.IsEmpty()) 
							Initialize_Forces(ex);
						
						int ih = 0;
						while (!in.IsEof()) {		
							line = in.GetLine();
							if (line.Find("Wave Heading (deg) :") >= 0) {
								in.GetLine(); in.GetLine(); in.GetLine();
								int idof = 0;
								while (!TrimBoth(line = in.GetLine()).IsEmpty()) {
									f.Load(line);
									double ma = f.GetDouble(1);
									double ph = f.GetDouble(2)*M_PI/180;
									double re = ma*cos(ph);
									double im = ma*sin(ph);
									int i = abs(f.GetInt(0)) - 1;
									ex.ma[ih](ifr, i) = ma;	
									ex.ph[ih](ifr, i) = ph;	
									ex.re[ih](ifr, i) = re;	
									ex.im[ih](ifr, i) = im;	
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


void Wamit::Load_A(FileIn &in, MatrixXd &A) {
	in.GetLine(); in.GetLine(); in.GetLine(); in.GetLine(); in.GetLine(); in.GetLine();
	while (!in.IsEof()) {
		String line = TrimBoth(in.GetLine());
		if (line.IsEmpty())
           	break;
		FieldSplit f;
		f.Load(line);
		int i = f.GetInt(0) - 1;
		int j = f.GetInt(1) - 1;
		double Aij = f.GetDouble(2);
		A(i, j) = Aij;
	}
}

bool Wamit::Load_Scattering(String fileName) {
	FileIn in(fileName);
	if (!in.IsOpen())
		return false;
	FieldSplit f;
	
	Initialize_Forces(sc);
	
	in.GetLine();
    for (int ifr = 0; ifr < Nf; ++ifr) {
        for (int ih = 0; ih < Nh; ++ih) {
            for (int k = 0; k < Nb*6; ++k) {		// Number of non-zero dof
        		f.Load(in.GetLine());
        		
        		int i = f.GetInt(2) - 1;
                sc.ma[ih](ifr, i) = f.GetDouble(3);
                sc.ph[ih](ifr, i) = f.GetDouble(4)*M_PI/180;
                sc.re[ih](ifr, i) = f.GetDouble(5);
                sc.im[ih](ifr, i) = f.GetDouble(6);
            }
        }
    }
    return true;
}
		
bool Wamit::Load_FK(String fileName) {
	FileIn in(fileName);
	if (!in.IsOpen())
		return false;
	FieldSplit f;
	
	Initialize_Forces(fk);
	
	in.GetLine();	
    for (int ifr = 0; ifr < Nf; ++ifr) {
        for (int ih = 0; ih < Nh; ++ih) {
            for (int k = 0; k < Nb*6; ++k) {		// Number of non-zero dof
                f.Load(in.GetLine());
        		
        		int i = f.GetInt(2) - 1;
                fk.ma[ih](ifr, i) = f.GetDouble(3);
                fk.ph[ih](ifr, i) = f.GetDouble(4)*M_PI/180;
                fk.re[ih](ifr, i) = f.GetDouble(5);
                fk.im[ih](ifr, i) = f.GetDouble(6);
            }
        }
    }	
    return true;
}
		
void Wamit::Save(String file) {

}