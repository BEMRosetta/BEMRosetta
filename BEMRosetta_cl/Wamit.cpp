#include "BEMRosetta.h"


bool Wamit::Load(String file, double rho) {
	hd().code = Hydro::WAMIT;
	hd().file = file;	
	hd().name = GetFileTitle(file);
	
	hd().rho = rho;
	
	try {
		hd().Print("\n\n" + Format("Loading '%s'", file));
		if (!Load_out()) {
			hd().Print("\n" + Format("File '%s' not found", file));
			return false;
		}
		String fileSC = ForceExt(file, ".3sc");
		hd().Print("\n" + Format("- Scattering file '%s'", GetFileName(fileSC)));
		if (!Load_Scattering(fileSC))
			hd().Print(": **Not found**");
		String fileFK = ForceExt(file, ".3fk");
		hd().Print("\n" + Format("- Froude-Krylov file '%s'", GetFileName(fileFK)));
		if (!Load_FK(fileFK))
			hd().Print(": **Not found**");
		
		hd().AfterLoad();
	} catch (Exc e) {
		hd().PrintError("\nError: " + e);
		hd().lastError = e;
		return false;
	}
	
	return true;
}

bool Wamit::Load_out() {
	hd().Nb = 0;
	hd().Nf = 0;
	hd().Nh = 0;
	int pos;
	int ibody = -1;
	
	FileIn in(hd().file);
	if (!in.IsOpen())
		return false;
	String line;
	FieldSplit f;
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
			hd().cg.setConstant(3, hd().Nb, nan(""));
			hd().cb.setConstant(3, hd().Nb, nan(""));
			hd().Vo.SetCount(hd().Nb, nan(""));
			hd().C.SetCount(hd().Nb);
		} else if (line.Find("Gravity:") >= 0) {
			hd().g = f.GetDouble(1);
			hd().len = f.GetDouble(4);
		} else if (line.Find("Water depth:") >= 0) {
			if (ToLower(f.GetText(2)) == "infinite")
				hd().h = INFINITY;
			else
				hd().h = f.GetDouble(2);
			if (line.Find("Water density:") >= 0) 
				hd().rho = f.GetDouble(5);			
		} else if (line.Find("XBODY =") >= 0) {
			ibody++;
			hd().cg(0, ibody) = f.GetDouble(2);
			hd().cg(1, ibody) = f.GetDouble(5);
			hd().cg(2, ibody) = f.GetDouble(8);
		} else if ((pos = line.FindAfter("Volumes (VOLX,VOLY,VOLZ):")) >= 0) 
			hd().Vo[ibody] = ScanDouble(line.Mid(pos));
		else if (line.Find("Center of Buoyancy (Xb,Yb,Zb):") >= 0) {
			hd().cb(0, ibody) = f.GetDouble(4) + hd().cg(0, ibody);
			hd().cb(1, ibody) = f.GetDouble(5) + hd().cg(1, ibody);
			hd().cb(2, ibody) = f.GetDouble(6) + hd().cg(2, ibody);
		} else if (line.Find("Hydrostatic and gravitational") >= 0) {
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
			int64 fpos = in.GetPos();
			
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
				throw Exc(Format("Wrong format in Wamit file '%s'", hd().file));
		
			hd().T.SetCount(hd().Nf);
			hd().w.SetCount(hd().Nf);
			
			in.Seek(fpos);
			while (in.GetLine().Find("Wave period = infinite") < 0 && !in.IsEof())
				; 
			if (!in.IsEof()) {
				hd().Aw0.setConstant(hd().Nb*6, hd().Nb*6, nan(""));
				Load_A(in, hd().Aw0);
			}
			in.Seek(fpos);
			while (in.GetLine().Find("Wave period = zero") < 0 && !in.IsEof())
				; 
			if (!in.IsEof()) {
				hd().Awinf.setConstant(hd().Nb*6, hd().Nb*6, nan(""));
				Load_A(in, hd().Awinf);
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
	            hd().T[ifr] = f.GetDouble(4);  			// Wave periods
	            hd().w[ifr] = 2*M_PI/hd().T[ifr];  		// Wave frequencies
	            
	            bool nextFreq = false;
	            while (!in.IsEof() && !nextFreq) {
	            	line = in.GetLine();
	            	if (line.Find("ADDED-MASS AND DAMPING COEFFICIENTS") >= 0) {
	            		if (hd().A.IsEmpty()) {
							hd().A.SetCount(hd().Nf);
							hd().B.SetCount(hd().Nf);
						}
						in.GetLine();	in.GetLine();
						hd().A[ifr].setConstant(hd().Nb*6, hd().Nb*6, nan(""));
		            	hd().B[ifr].setConstant(hd().Nb*6, hd().Nb*6, nan(""));
		            
			            while (!in.IsEof()) {
							line = TrimBoth(in.GetLine());
							if (line.IsEmpty())
			                	break;
							f.Load(line);
							int i = f.GetInt(0) - 1;
							int j = f.GetInt(1) - 1;
							double Aij = f.GetDouble(2);
							double Bij = f.GetDouble(3);
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
								in.GetLine(); in.GetLine(); in.GetLine();
								int idof = 0;
								while (!TrimBoth(line = in.GetLine()).IsEmpty()) {
									f.Load(line);
									double ma = f.GetDouble(1);
									double ph = f.GetDouble(2)*M_PI/180;
									double re = ma*cos(ph);
									double im = ma*sin(ph);
									int i = abs(f.GetInt(0)) - 1;
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
	
	hd().Initialize_Forces(hd().sc);
	
	in.GetLine();
    for (int ifr = 0; ifr < hd().Nf; ++ifr) {
        for (int ih = 0; ih < hd().Nh; ++ih) {
            for (int k = 0; k < hd().Nb*6; ++k) {		// Number of non-zero dof
        		f.Load(in.GetLine());
        		
        		int i = f.GetInt(2) - 1;
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
	FileIn in(fileName);
	if (!in.IsOpen())
		return false;
	FieldSplit f;
	
	hd().Initialize_Forces(hd().fk);
	
	in.GetLine();	
    for (int ifr = 0; ifr < hd().Nf; ++ifr) {
        for (int ih = 0; ih < hd().Nh; ++ih) {
            for (int k = 0; k < hd().Nb*6; ++k) {		// Number of non-zero dof
                f.Load(in.GetLine());
        		
        		int i = f.GetInt(2) - 1;
                hd().fk.ma[ih](ifr, i) = f.GetDouble(3);
                hd().fk.ph[ih](ifr, i) = f.GetDouble(4)*M_PI/180;
                hd().fk.re[ih](ifr, i) = f.GetDouble(5);
                hd().fk.im[ih](ifr, i) = f.GetDouble(6);
            }
        }
    }	
    return true;
}
		
void Wamit::Save(String file) {
	throw Exc("Option not implemented");
}

bool Wamit::LoadMesh(String fileName) {
	FileIn in(fileName);
	if (!in.IsOpen()) {
		hd().PrintError("\n" + Format("Impossible to open file '%s'", fileName));
		mh().lastError = Format("Impossible to open file '%s'", fileName);
		return false;
	}
	mh().file = fileName;
	
	String line;
	FieldSplit f;	
	
	try {
		in.GetLine();
		line = in.GetLine();	
		f.Load(line);
		double scale = f.GetDouble(0);
		
		line = in.GetLine();	
		f.Load(line);
		mh().y0z = f.GetInt(0) != 0;
		mh().x0z = f.GetInt(1) != 0;
		
		line = in.GetLine();	
		f.Load(line);
		int nPatches = f.GetInt(0);
		
		mh().nodes.Clear();
		mh().panels.Clear();
		
		while(!in.IsEof()) {
			int ids[4];
			for (int i = 0; i < 4; ++i) {
				line = in.GetLine();	
				f.Load(line);
				
				double x = f.GetDouble(0)*scale;	
				double y = f.GetDouble(1)*scale;	
				double z = f.GetDouble(2)*scale;	
				
				bool found = false;
				for (int in = 0; in < mh().nodes.GetCount(); ++in) {
					Point3D &node = mh().nodes[in];
					if (x == node.x && y == node.y && z == node.z) {
						ids[i] = in;
						found = true;
						break;
					}
				}
				if (!found) {
					Point3D &node = mh().nodes.Add();
					node.x = x;
					node.y = y;
					node.z = z;
					ids[i] = mh().nodes.GetCount() - 1;
				}
			}
			Panel &panel = mh().panels.Add();
			panel.id0 = ids[0];
			panel.id1 = ids[1];
			panel.id2 = ids[2];
			panel.id3 = ids[3];
		}
		if (mh().panels.GetCount() != nPatches)
			throw Exc("Wrong number of patches in .gdf file");
		if (mh().Check())
			throw Exc("Wrong nodes found in Wamit .gdf mesh file");
		mh().GetLimits();
	} catch (Exc e) {
		hd().PrintError("\nError: " + e);
		mh().lastError = e;
		return false;
	}
	
	return true;
}