#include "BEMRosetta.h"

#include <plugin/matio/matio.h>

bool Foamm::Load(String file) {
	hd().code = Hydro::FOAMM;
	hd().file = file;	
	hd().name = GetFileTitle(file);
	
	try {
		if (GetFileExt(file) == ".mat") {
			hd().Print("\n\n" + Format(t_("Loading mat file '%s'"), file));
			if (!Load_mat(file)) {
				hd().PrintWarning("\n" + Format(t_("File '%s' not found"), file));
				return false;
			}
		}
		hd().dof.Clear();	hd().dof.SetCount(hd().Nb, 0);
		for (int i = 0; i < hd().Nb; ++i)
			hd().dof[i] = 1;
	
		hd().AfterLoad();
	} catch (Exc e) {
		hd().PrintError(Format("\n%s: %s", t_("Error"), e));
		hd().lastError = e;
		return false;
	}
	
	return true;
}

void Foamm::Save(String file) {
	try {
		String file1 = ForceExt(file, ".mat");
		hd().Print("\n- " + Format(t_("FOAMM file '%s'"), GetFileName(file1)));
		Save_mat(file1);
	} catch (Exc e) {
		hd().PrintError(Format("\n%s: %s", t_("Error"), e));
		hd().lastError = e;
	}
}

bool Foamm::Load_mat(String file) {
	hd().Nb = 1;
	hd().dof.SetCount(1);
	hd().dof[0] = 1;
	hd().Nh = 1;
	hd().head.SetCount(1);
	hd().head[0] = 0;
	
	MatFile mat;
	
	if (!mat.OpenRead(file)) 
		return false;

	MatMatrix<double> w = mat.VarReadMat<double>("w");	
	if (w.GetCount() == 0)
		throw Exc("Vector w not found");
	
	hd().Nf = w.GetCount();
	
	hd().w.SetCount(hd().Nf);
	for (int ifr = 0; ifr < hd().Nf; ++ifr) 
		hd().w[ifr] = w[ifr];
	
	MatMatrix<double> A = mat.VarReadMat<double>("A");	
	if (A.GetCount() == 0)
		throw Exc("Vector A not found");
	
	MatMatrix<double> B = mat.VarReadMat<double>("B");	
	if (B.GetCount() == 0)
		throw Exc("Vector B not found");
	
	if (w.GetCount() != A.GetCount())
		throw Exc("Vectors w and A size does not match");
	if (w.GetCount() != B.GetCount())
		throw Exc("Vectors w and B size does not match");
	
	hd().A.SetCount(hd().Nf);
	hd().B.SetCount(hd().Nf);
	for (int ifr = 0; ifr < hd().Nf; ++ifr) {
		hd().A[ifr].setConstant(hd().Nb*6, hd().Nb*6, Null);
		hd().B[ifr].setConstant(hd().Nb*6, hd().Nb*6, Null);	
	}
	for (int ifr = 0; ifr < hd().Nf; ++ifr) {
		hd().A[ifr](0, 0) = A[ifr];
		hd().B[ifr](0, 0) = B[ifr];
	}
	
	hd().names << "Body";
	
	double Mu = mat.VarRead<double>("Mu");
	if (!IsNull(Mu)) {
		hd().Awinf.setConstant(hd().Nb*6, hd().Nb*6, Null);
		hd().Awinf(0, 0) = Mu;
	}

	return true;
}

void Foamm::Save_mat(String file) {
	
}
