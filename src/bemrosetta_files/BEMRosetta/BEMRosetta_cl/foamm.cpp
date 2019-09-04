#include "BEMRosetta.h"

#include <plugin/matio/matio.h>

bool Foamm::Load(String file) {
	hd().code = Hydro::FOAMM;
	hd().file = file;	
	hd().name = GetFileTitle(file);
	hd().len = 1;
	hd().dimen = true;
	hd().Nb = Null;
	
	try {
		if (GetFileExt(file) == ".mat") {
			BEMData::Print("\n\n" + Format(t_("Loading mat file '%s'"), file));
			if (!Load_mat(file)) {
				BEMData::PrintWarning("\n" + Format(t_("File '%s' not found"), file));
				return false;
			}
		}
		if (IsNull(hd().Nb))
			return false;
		
		hd().dof.Clear();	hd().dof.SetCount(hd().Nb, 0);
		for (int i = 0; i < hd().Nb; ++i)
			hd().dof[i] = 1;
	} catch (Exc e) {
		BEMData::PrintError(Format("\n%s: %s", t_("Error"), e));
		hd().lastError = e;
		return false;
	}
	
	return true;
}

void Foamm::Save(String file) {
	try {
		String file1 = ForceExt(file, ".mat");
		BEMData::Print("\n- " + Format(t_("FOAMM file '%s'"), GetFileName(file1)));
		Save_mat(file1);
	} catch (Exc e) {
		BEMData::PrintError(Format("\n%s: %s", t_("Error"), e));
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
		throw Exc(x_("\n") + t_("Vector w not found"));
		
	hd().Nf = w.GetCount();
	
	hd().w.SetCount(hd().Nf);
	hd().T.SetCount(hd().Nf);
	hd().dataFromW = true;
	for (int ifr = 0; ifr < hd().Nf; ++ifr) {
		hd().w[ifr] = w[ifr];
		hd().T[ifr] = 2*M_PI/w[ifr];
	}
	
	MatMatrix<double> A = mat.VarReadMat<double>("A");	
	if (A.GetCount() == 0)
		BEMData::Print(x_("\n") + t_("Vector A not found"));
	else {
		hd().A.SetCount(hd().Nf);
		if (hd().Nf != A.GetCount())
			throw Exc(x_("\n") + t_("Vectors w and A size does not match"));
		for (int ifr = 0; ifr < hd().Nf; ++ifr) 
			hd().A[ifr].setConstant(hd().Nb*6, hd().Nb*6, Null);
		for (int ifr = 0; ifr < hd().Nf; ++ifr) 
			hd().A[ifr](0, 0) = A[ifr];
	}
	
	MatMatrix<double> B = mat.VarReadMat<double>("B");	
	if (B.GetCount() == 0)
		BEMData::Print(x_("\n") + t_("Vector B not found"));
	else {
		hd().B.SetCount(hd().Nf);
		if (hd().Nf != B.GetCount())
			throw Exc(x_("\n") + t_("Vectors w and B size does not match"));
		for (int ifr = 0; ifr < hd().Nf; ++ifr) 
			hd().B[ifr].setConstant(hd().Nb*6, hd().Nb*6, Null);	
		for (int ifr = 0; ifr < hd().Nf; ++ifr) 
			hd().B[ifr](0, 0) = B[ifr];
	}
	
	hd().names << "Body";
	
	double Mu = mat.VarRead<double>("Mu");
	if (!IsNull(Mu)) {
		hd().Awinf.setConstant(hd().Nb*6, hd().Nb*6, Null);
		hd().Awinf(0, 0) = Mu;
	}

	MatMatrix<std::complex<double>> Z = mat.VarReadMat<std::complex<double>>("Z");	
	if (Z.GetCount() == 0)
		BEMData::Print(x_("\n") + t_("Vector Z not found"));
	else {
		hd().Z.SetCount(hd().Nf);
		if (hd().Nf != Z.GetCount())
			throw Exc(x_("\n") + t_("Vectors w and Z size does not match"));
		for (int ifr = 0; ifr < hd().Nf; ++ifr) 
			hd().Z[ifr] = Z[ifr];
	}

	MatMatrix<std::complex<double>> TFSResponse = mat.VarReadMat<std::complex<double>>("TFSResponse");	
	if (TFSResponse.GetCount() == 0)
		BEMData::Print(x_("\n") + t_("Vector TFSResponse not found"));
	else {
		hd().TFSResponse.SetCount(hd().Nf);
		if (hd().Nf != TFSResponse.GetCount())
			throw Exc(x_("\n") + t_("Vectors w and TFSResponse size does not match"));
		for (int ifr = 0; ifr < hd().Nf; ++ifr) 
			hd().TFSResponse[ifr] = TFSResponse[ifr];
	}
	
	return true;
}

void Foamm::Save_mat(String ) {
	throw Exc("Option not implemented");
}
