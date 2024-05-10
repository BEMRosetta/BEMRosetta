// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"

#include <MatIO/matio.h>

String Foamm::Load(String file) {
	hd().solver = Hydro::FOAMM;
	hd().file = file;	
	hd().name = GetFileTitle(file);
	hd().len = 1;
	hd().dimen = true;
	hd().Nb = 1;
	hd().msh.SetCount(1);
	//hd().dof.SetCount(1);
	//hd().dof[0] = 1;
	hd().Nh = 1;
	hd().head.SetCount(1);
	hd().head[0] = 0;
	
	try {
		if (GetFileExt(file) == ".mat") {
			BEM::Print("\n\n" + Format(t_("Loading mat file '%s'"), file));
			
			Load_mat(file, 0, 0, true);
		}
		if (IsNull(hd().Nb))
			return t_("No data found");
		
		/*hd().dof.Clear();	hd().dof.SetCount(hd().Nb, 0);
		for (int i = 0; i < hd().Nb; ++i)
			hd().dof[i] = 1;*/
	} catch (Exc e) {
		//BEM::PrintError(Format("\n%s: %s", t_("Error"), e));
		//hd().lastError = Format(t_("file %s "), file) + e;
		return e;
	}
	
	return String();
}

void Foamm::Load_mat(String file, int idf, int jdf, bool loadCoeff) {
	MatFile mat;
	
	if (!mat.OpenRead(file)) 
		throw Exc(S("\n") + t_("File not found or blocked"));
	
	hd().stsProcessor = "FOAMM by COER (http://www.eeng.nuim.ie/coer/)";
	
	if (loadCoeff) {
		MatMatrix<double> w = mat.VarReadMat<double>("w");	
		if (w.size() == 0)
			throw Exc(S("\n") + t_("Vector w not found"));
			
		hd().Nf = w.size();
		
		hd().w.SetCount(hd().Nf);
		//hd().T.SetCount(hd().Nf);
		hd().dataFromW = true;
		for (int ifr = 0; ifr < hd().Nf; ++ifr) {
			hd().w[ifr] = w[ifr];
			//hd().T[ifr] = 2*M_PI/w[ifr];
		}
		
		MatMatrix<double> A = mat.VarReadMat<double>("A");	
		if (A.size() == 0)
			BEM::Print(S("\n") + t_("Vector A not found"));
		else {
			if (hd().Nf != A.size())
				throw Exc(S("\n") + t_("Vectors w and A size does not match"));
			hd().A[idf][jdf].setConstant(hd().Nf);	
			for (int ifr = 0; ifr < hd().Nf; ++ifr) 
				hd().A[idf][jdf][ifr] = A[ifr];
		}
	
		MatMatrix<double> B = mat.VarReadMat<double>("B");	
		if (B.size() == 0)
			BEM::Print(S("\n") + t_("Vector B not found"));
		else {
			if (hd().Nf != B.size())
				throw Exc(S("\n") + t_("Vectors w and A size does not match"));
			hd().B[idf][jdf].setConstant(hd().Nf);	
			for (int ifr = 0; ifr < hd().Nf; ++ifr) 
				hd().B[idf][jdf][ifr] = B[ifr];
		}
		
		hd().msh[0].name = "Body";
		
		double Mu = mat.VarRead<double>("Mu");
		if (!IsNull(Mu)) {
			hd().Ainf.setConstant(hd().Nb*6, hd().Nb*6, 0);
			hd().Ainf(idf, jdf) = Mu;
		}
	}
	
	hd().Initialize_Sts();
	Hydro::StateSpace &sts = hd().sts[idf][jdf];

	MatMatrix<std::complex<double>> TFS = mat.VarReadMat<std::complex<double>>("TFSResponse");	
	if (TFS.size() == 0)
		BEM::Print(S("\n") + t_("Vector TFSResponse not found"));
	else {
		sts.TFS.SetCount(hd().Nf);
		if (hd().Nf != TFS.size())
			throw Exc(S("\n") + t_("Vectors w and TFSResponse size does not match"));
		for (int ifr = 0; ifr < hd().Nf; ++ifr) 
			sts.TFS[ifr] = TFS[ifr];
	}
	
	MatMatrix<double> A_ss = mat.VarReadMat<double>("A_ss");	
	if (A_ss.size() == 0)
		BEM::Print(S("\n") + t_("Matrix A_ss not found"));
	else {
		sts.A_ss.setConstant(A_ss.GetRows(), A_ss.GetCols(), Null);
		for (int r = 0; r < A_ss.GetRows(); ++r)
			for (int c = 0; c < A_ss.GetCols(); ++c)
				sts.A_ss(r, c) = A_ss(r, c);
	}
	
	MatMatrix<double> B_ss = mat.VarReadMat<double>("B_ss");	
	if (B_ss.size() == 0)
		BEM::Print(S("\n") + t_("Matrix B_ss not found"));
	else {
		sts.B_ss.setConstant(B_ss.GetRows(), Null);
		for (int r = 0; r < B_ss.GetRows(); ++r)
			sts.B_ss(r) = B_ss(r, 0);
	}

	MatMatrix<double> C_ss = mat.VarReadMat<double>("C_ss");	
	if (C_ss.size() == 0)
		BEM::Print(S("\n") + t_("Matrix C_ss not found"));
	else {
		sts.C_ss.setConstant(C_ss.GetCols(), Null);
		for (int c = 0; c < C_ss.GetCols(); ++c)
			sts.C_ss(c) = C_ss(0, c);
	}
	
	MatMatrix<double> ssFrequencies = mat.VarReadMat<double>("Frequencies");	
	if (ssFrequencies.GetCols() == 0)
		BEM::Print(S("\n") + t_("Matrix Frequencies not found"));
	else {
		sts.ssFrequencies.setConstant(ssFrequencies.GetCols(), Null);
		for (int c = 0; c < ssFrequencies.GetCols(); ++c)
			sts.ssFrequencies[c] = ssFrequencies(0, c);
	}

	MatMatrix<double> ssFreqRange = mat.VarReadMat<double>("FreqRange");	
	if (ssFreqRange.GetCols() == 0)
		BEM::Print(S("\n") + t_("Matrix FreqRange not found"));
	else {
		sts.ssFreqRange.setConstant(ssFreqRange.GetCols(), Null);
		for (int c = 0; c < ssFreqRange.GetCols(); ++c)
			sts.ssFreqRange[c] = ssFreqRange(0, c);
	}

	sts.ssMAE = mat.VarRead<double>("MAE");	
	
	hd().dimenSTS = true;	
}

void Foamm::Get(const UVector<int> &ibs, const UVector<int> &idfs, const UVector<int> &jdfs,
		const UVector<double> &froms, const UVector<double> &tos, const UVector<UVector<double>> &freqs, 
		Function <bool(String, int)> Status, Function <void(String)> FOAMMMessage) {
	if (!FileExists(Bem().foammPath))
		throw Exc(t_("FOAMM not found. Please set FOAMM path in Options"));
	for (int i = 0; i < ibs.size(); ++i) {
		Status(Format(t_("Processing case %d"), i+1), int((100*i)/ibs.size()));
		Get_Each(ibs[i], idfs[i], jdfs[i], froms[i], tos[i], freqs[i], Status, FOAMMMessage);
	}
}

void Foamm::Get_Each(int ibody, int _idf, int _jdf, double from, double to, const UVector<double> &freqs, 
		Function <bool(String, int)> Status, Function <void(String)> FOAMMMessage) {
	Uuid id = Uuid::Create();
	String folder = AFX(BEM::GetTempFilesFolder(), Format(id));
	if (!DirectoryCreateX(folder))
		throw Exc(Format(t_("Problem creating temporary FOAMM folder '%s'"), folder));			
	String file = AFX(folder, "temp_file.mat");
	
	MatFile mat;
	
	if (!mat.OpenCreate(file, MAT_FT_MAT5)) 
		throw Exc(Format(t_("Problem creating FOAMM file '%s'"), file));

	int idf = ibody*6 + _idf;
	int jdf = ibody*6 + _jdf;

	MatMatrix<double> matA(hd().Nf, 1);
	for (int ifr = 0; ifr < hd().Nf; ++ifr)
		matA(ifr, 0) = hd().A_dim(ifr, idf, jdf);
 	if (!mat.VarWrite("A", matA))
 		throw Exc(Format(t_("Problem writing %s to file '%s'"), "A", file));

 	if (!mat.VarWrite<double>("Mu", hd().Ainf_dim(idf, jdf)))
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
		if (i == idf)
			matDof(0, i) = 1;
		else
			matDof(0, i) = 0;
	}
	if (!mat.VarWrite("Dof", matDof))
 		throw Exc(Format(t_("Problem writing %s to file '%s'"), "Dof", file));*/
	
	UVector<String> optionsVars;
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

	MatMatrix<double> freqChoice(1, freqs.size());
	for (int i = 0; i < freqs.size(); ++i) 
		freqChoice(0, i) = freqs[i];
	options.VarWriteStruct<double>("FreqChoice", freqChoice);
	
	UVector<String> optimVars;
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
	if (!process.Start(Bem().foammPath, NULL, folder))
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
		DeleteFolderDeep(folder);	Sleep(100);
		throw Exc(t_("Process ended by user"));
	}
	if (process.GetExitCode() != 0) {
		DeleteFolderDeep(folder);	Sleep(100);
		throw Exc(t_("FOAMM ended with error"));
	}
	Load_mat(file, idf, jdf, false);
	DeleteFolderDeep(folder);		Sleep(100);
}

