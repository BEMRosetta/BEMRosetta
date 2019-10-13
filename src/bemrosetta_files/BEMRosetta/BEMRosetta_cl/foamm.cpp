#include "BEMRosetta.h"

#include <plugin/matio/matio.h>

bool Foamm::Load(String file) {
	hd().code = Hydro::FOAMM;
	hd().file = file;	
	hd().name = GetFileTitle(file);
	hd().len = 1;
	hd().dimen = true;
	hd().Nb = 1;
	hd().dof.SetCount(1);
	hd().dof[0] = 1;
	hd().Nh = 1;
	hd().head.SetCount(1);
	hd().head[0] = 0;
	
	try {
		if (GetFileExt(file) == ".mat") {
			BEMData::Print("\n\n" + Format(t_("Loading mat file '%s'"), file));
			if (!Load_mat(file, 0, 0, true)) {
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

bool Foamm::Load_mat(String file, int idof, int jdof, bool loadCoeff) {
	MatFile mat;
	
	if (!mat.OpenRead(file)) 
		return false;
	
	if (loadCoeff) {
		MatMatrix<double> w = mat.VarReadMat<double>("w");	
		if (w.GetCount() == 0)
			throw Exc(S("\n") + t_("Vector w not found"));
			
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
			BEMData::Print(S("\n") + t_("Vector A not found"));
		else {
			hd().A.SetCount(hd().Nf);
			if (hd().Nf != A.GetCount())
				throw Exc(S("\n") + t_("Vectors w and A size does not match"));
			for (int ifr = 0; ifr < hd().Nf; ++ifr) 
				hd().A[ifr].setConstant(hd().Nb*6, hd().Nb*6, Null);
			for (int ifr = 0; ifr < hd().Nf; ++ifr) 
				hd().A[ifr](idof, jdof) = A[ifr];
		}
		
		MatMatrix<double> B = mat.VarReadMat<double>("B");	
		if (B.GetCount() == 0)
			BEMData::Print(S("\n") + t_("Vector B not found"));
		else {
			hd().B.SetCount(hd().Nf);
			if (hd().Nf != B.GetCount())
				throw Exc(S("\n") + t_("Vectors w and B size does not match"));
			for (int ifr = 0; ifr < hd().Nf; ++ifr) 
				hd().B[ifr].setConstant(hd().Nb*6, hd().Nb*6, Null);	
			for (int ifr = 0; ifr < hd().Nf; ++ifr) 
				hd().B[ifr](idof, jdof) = B[ifr];
		}
		
		hd().names << "Body";
		
		double Mu = mat.VarRead<double>("Mu");
		if (!IsNull(Mu)) {
			hd().Awinf.setConstant(hd().Nb*6, hd().Nb*6, Null);
			hd().Awinf(idof, jdof) = Mu;
		}
	}
	
	hd().InitializeSts();
	Hydro::StateSpace &sts = hd().sts[idof][jdof];
	
	MatMatrix<std::complex<double>> Z = mat.VarReadMat<std::complex<double>>("Z");	
	if (Z.GetCount() == 0)
		BEMData::Print(S("\n") + t_("Vector Z not found"));
	else {
		sts.Z.SetCount(hd().Nf);
		if (hd().Nf != Z.GetCount())
			throw Exc(S("\n") + t_("Vectors w and Z size does not match"));
		for (int ifr = 0; ifr < hd().Nf; ++ifr) 
			sts.Z[ifr] = Z[ifr];
	}

	MatMatrix<std::complex<double>> TFSResponse = mat.VarReadMat<std::complex<double>>("TFSResponse");	
	if (TFSResponse.GetCount() == 0)
		BEMData::Print(S("\n") + t_("Vector TFSResponse not found"));
	else {
		sts.TFSResponse.SetCount(hd().Nf);
		if (hd().Nf != TFSResponse.GetCount())
			throw Exc(S("\n") + t_("Vectors w and TFSResponse size does not match"));
		for (int ifr = 0; ifr < hd().Nf; ++ifr) 
			sts.TFSResponse[ifr] = TFSResponse[ifr];
	}
	
	MatMatrix<double> A_ss = mat.VarReadMat<double>("A_ss");	
	if (A_ss.GetCount() == 0)
		BEMData::Print(S("\n") + t_("Matrix A_ss not found"));
	else {
		sts.A_ss.setConstant(A_ss.GetRows(), A_ss.GetCols(), Null);
		for (int r = 0; r < A_ss.GetRows(); ++r)
			for (int c = 0; c < A_ss.GetCols(); ++c)
				sts.A_ss(r, c) = A_ss(r, c);
	}
	
	MatMatrix<double> B_ss = mat.VarReadMat<double>("B_ss");	
	if (B_ss.GetCount() == 0)
		BEMData::Print(S("\n") + t_("Matrix B_ss not found"));
	else {
		sts.B_ss.setConstant(B_ss.GetRows(), Null);
		for (int r = 0; r < B_ss.GetRows(); ++r)
			sts.B_ss(r) = B_ss(r, 0);
	}

	MatMatrix<double> C_ss = mat.VarReadMat<double>("C_ss");	
	if (C_ss.GetCount() == 0)
		BEMData::Print(S("\n") + t_("Matrix C_ss not found"));
	else {
		sts.C_ss.setConstant(C_ss.GetCols(), Null);
		for (int c = 0; c < C_ss.GetCols(); ++c)
			sts.C_ss(c) = C_ss(0, c);
	}
	
	MatMatrix<double> ssFrequencies = mat.VarReadMat<double>("Frequencies");	
	if (ssFrequencies.GetCols() == 0)
		BEMData::Print(S("\n") + t_("Matrix Frequencies not found"));
	else {
		sts.ssFrequencies.setConstant(ssFrequencies.GetCols(), Null);
		for (int c = 0; c < ssFrequencies.GetCols(); ++c)
			sts.ssFrequencies[c] = ssFrequencies(0, c);
	}

	MatMatrix<double> ssFreqRange = mat.VarReadMat<double>("FreqRange");	
	if (ssFreqRange.GetCols() == 0)
		BEMData::Print(S("\n") + t_("Matrix FreqRange not found"));
	else {
		sts.ssFreqRange.setConstant(ssFreqRange.GetCols(), Null);
		for (int c = 0; c < ssFreqRange.GetCols(); ++c)
			sts.ssFreqRange[c] = ssFreqRange(0, c);
	}

	sts.ssMAE = mat.VarRead<double>("MAE");		
			
	return true;
}

void Foamm::Get(const Vector<int> &ibs, const Vector<int> &idofs, const Vector<int> &jdofs,
		const Vector<double> &froms, const Vector<double> &tos, const Vector<Vector<double>> &freqs, 
		Function <bool(String, int)> Status, Function <void(String)> FOAMMMessage) {
	for (int i = 0; i < ibs.GetCount(); ++i) {
		Status(Format(t_("Processing case %d"), i+1), int((100*i)/ibs.GetCount()));
		Get_Each(ibs[i], idofs[i], jdofs[i], froms[i], tos[i], freqs[i], Status, FOAMMMessage);
	}
}

void Foamm::Get_Each(int ibody, int idof, int jdof, double from, double to, const Vector<double> &freqs, 
		Function <bool(String, int)> Status, Function <void(String)> FOAMMMessage) {
	Uuid id = Uuid::Create();
	String folder = AppendFileName(BEMData::GetTempFilesFolder(), Format(id));
	if (!DirectoryCreate(folder))
		throw Exc(Format(t_("Problem creating temporary FOAMM folder '%s'"), folder));			
	String file = AppendFileName(folder, "temp_file.mat");
	
	MatFile mat;
	
	if (!mat.OpenCreate(file, MAT_FT_MAT5)) 
		throw Exc(Format(t_("Problem creating FOAMM file '%s'"), file));

	int idf = ibody*6 + idof;
	int jdf = ibody*6 + jdof;

	MatMatrix<double> matA(hd().Nf, 1);
	for (int ifr = 0; ifr < hd().Nf; ++ifr)
		matA(ifr, 0) = hd().A_dim(ifr, idf, jdf);
 	if (!mat.VarWrite("A", matA))
 		throw Exc(Format(t_("Problem writing %s to file '%s'"), "A", file));

 	if (!mat.VarWrite<double>("Mu", hd().Awinf_dim(idf, jdf)))
 		throw Exc(Format(t_("Problem writing %s to file '%s'"), "Mu", file));
 		 	
	MatMatrix<double> matB(hd().Nf, 1);
	for (int ifr = 0; ifr < hd().Nf; ++ifr)
		matB(ifr, 0) = hd().B_dim(ifr, idf, jdf);
	if (!mat.VarWrite("B", matB))
 		throw Exc(Format(t_("Problem writing %s to file '%s'"), "B", file));
	
	MatMatrix<double> matw(1, hd().Nf);
	for (int ifr = 0; ifr < hd().Nf; ++ifr)
		matw(0, ifr) = hd().w[ifr];
	if (!mat.VarWrite("w", matw))
 		throw Exc(Format(t_("Problem writing %s to file '%s'"), "w", file));
	
	/*MatMatrix<double> matDof(1, 6);
	for (int i = 0; i < 6; ++i) {
		if (i == idof)
			matDof(0, i) = 1;
		else
			matDof(0, i) = 0;
	}
	if (!mat.VarWrite("Dof", matDof))
 		throw Exc(Format(t_("Problem writing %s to file '%s'"), "Dof", file));*/
	
	Vector<String> optionsVars;
	optionsVars << "Mode" << "Method" << "FreqRangeChoice" << "FreqChoice";
	MatVar options("Options", 1, 1, optionsVars);
	
	options.VarWriteStruct<double>("Mode", 0);
	options.VarWriteStruct<double>("Method", 0);
	//options.VarWriteStruct("FreqRangeChoice", "G");
	//options.VarWriteStruct("FreqChoice", "G");
	
	MatMatrix<double> freqRangeChoice(1, 2);
	freqRangeChoice(0, 0) = from;
	freqRangeChoice(0, 1) = to;
	options.VarWriteStruct<double>("FreqRangeChoice", freqRangeChoice);

	MatMatrix<double> freqChoice(1, freqs.GetCount());
	for (int i = 0; i < freqs.GetCount(); ++i) 
		freqChoice(0, i) = freqs[i];
	options.VarWriteStruct<double>("FreqChoice", freqChoice);
	
	Vector<String> optimVars;
	optimVars << "InitCond" << "Tol" << "maxEval" << "maxIter" << "StepTol" << "ThresRel" << "ThresAbs";
	MatVar optim("Optim", 1, 1, optimVars);
	optim.VarWriteStruct<double>("InitCond", 50);
	optim.VarWriteStruct<double>("Tol", 1E-5);
	optim.VarWriteStruct<double>("maxEval", 1E3);
	optim.VarWriteStruct<double>("maxIter", 200);
	optim.VarWriteStruct<double>("StepTol", 1E-6);
	optim.VarWriteStruct<double>("ThresRel", 0.03);
	optim.VarWriteStruct<double>("ThresAbs", 0.1);
	
	options.VarWriteStruct("Optim", optim);
	mat.VarWrite(options);
	
	mat.Close();
	
	LocalProcess process;
	if (!process.Start(hd().GetBEMData().foammPath, NULL, folder))
		throw Exc(Format(t_("Problem launching FOAMM from '%s'"), file));

	String msg, reso, rese;
	bool endProcess = false;
	while (process.IsRunning()) {
		if (!endProcess) {
			if (process.Read2(reso, rese)) {
				msg.Clear();
				if (!reso.IsEmpty())
					msg << reso;
				if (!rese.IsEmpty()) {
					if (!msg.IsEmpty())
						msg << ";";
					msg << rese;
				}
				if (!msg.IsEmpty() && !msg.StartsWith("MAE:"))
					FOAMMMessage(msg);
			}
			endProcess = Status("", Null); 
			if (endProcess)
				process.Kill();
		}
		Sleep(200);
	}
	if (endProcess) {
		DeleteFolderDeep(folder);
		throw Exc(t_("Process ended by user"));
	}
	Load_mat(file, idof, jdof, false);
	DeleteFolderDeep(folder);
}

