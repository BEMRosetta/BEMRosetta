// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"

#include <MatIO/matio.h>

String Foamm::Load(String file) {
	dt.solver = Hydro::FOAMM;
	dt.file = file;	
	dt.name = GetFileTitle(file);
	dt.len = 1;
	dt.dimen = true;
	dt.Nb = 1;
	dt.msh.SetCount(1);
	//dt.dof.SetCount(1);
	//dt.dof[0] = 1;
	dt.Nh = 1;
	dt.head.SetCount(1);
	dt.head[0] = 0;
	
	try {
		if (GetFileExt(file) == ".mat") {
			BEM::Print("\n\n" + Format(t_("Loading mat file '%s'"), file));
			
			Load_mat(file, 0, 0, true);
		}
		if (IsNull(dt.Nb))
			return t_("No data found");
		
		/*dt.dof.Clear();	dt.dof.SetCount(dt.Nb, 0);
		for (int i = 0; i < dt.Nb; ++i)
			dt.dof[i] = 1;*/
	} catch (Exc e) {
		//BEM::PrintError(Format("\n%s: %s", t_("Error"), e));
		//dt.lastError = Format(t_("file %s "), file) + e;
		return e;
	}
	
	return String();
}

void Foamm::Load_mat(String file, int idf, int jdf, bool loadCoeff) {
	MatFile mat;
	
	if (!mat.OpenRead(file)) 
		throw Exc(S("\n") + t_("File not found or blocked"));
	
	dt.stsProcessor = "FOAMM by COER (http://www.eeng.nuim.ie/coer/)";
	
	if (loadCoeff) {
		if (!mat.Exist("w"))
			throw Exc(S("\n") + t_("Vector w not found"));
		mat.Get<double>("w", dt.w);	
			
		dt.Nf = dt.w.size();
		
		if (!mat.Exist("A"))
			throw Exc(S("\n") + t_("Vector A not found"));		
		UVector<double> A;
		mat.Get<double>("A", A);	
		if (dt.Nf != A.size())
			throw Exc(S("\n") + t_("Vectors w and A size does not match"));
		dt.A[idf][jdf].setConstant(dt.Nf);	
		for (int ifr = 0; ifr < dt.Nf; ++ifr) 
			dt.A[idf][jdf][ifr] = A[ifr];
	
		if (!mat.Exist("B"))
			throw Exc(S("\n") + t_("Vector B not found"));		
		UVector<double> B;
		mat.Get<double>("B", B);	
		if (dt.Nf != B.size())
			throw Exc(S("\n") + t_("Vectors w and B size does not match"));
		dt.B[idf][jdf].setConstant(dt.Nf);	
		for (int ifr = 0; ifr < dt.Nf; ++ifr) 
			dt.B[idf][jdf][ifr] = B[ifr];
		
		dt.msh[0].dt.name = "Body";
		
		double Mu = mat.Get<double>("Mu");
		if (!IsNull(Mu)) {
			dt.Ainf.setConstant(dt.Nb*6, dt.Nb*6, 0);
			dt.Ainf(idf, jdf) = Mu;
		}
	}
	
	Initialize_Sts();
	Hydro::StateSpace &sts = dt.sts[idf][jdf];

	if (!mat.Exist("TFSResponse"))
		throw Exc(S("\n") + t_("Vector TFSResponse not found"));
	UVector<std::complex<double>> TFS;
	mat.Get<std::complex<double>>("TFSResponse", TFS);	
	sts.TFS.SetCount(dt.Nf);
	if (dt.Nf != TFS.size())
		throw Exc(S("\n") + t_("Vectors w and TFSResponse size does not match"));
	for (int ifr = 0; ifr < dt.Nf; ++ifr) 
		sts.TFS[ifr] = TFS[ifr];
	
	if (!mat.Exist("A_ss"))
		throw Exc(S("\n") + t_("Matrix A_ss not found"));
	MatrixXd A_ss;
	mat.Get("A_ss", A_ss);	
	sts.A_ss.setConstant(A_ss.rows(), A_ss.cols(), Null);
	for (int r = 0; r < A_ss.rows(); ++r)
		for (int c = 0; c < A_ss.cols(); ++c)
			sts.A_ss(r, c) = A_ss(r, c);
	
	if (!mat.Exist("B_ss"))
		throw Exc(S("\n") + t_("Matrix B_ss not found"));
	MatrixXd B_ss;
	mat.Get("B_ss", B_ss);	
	sts.B_ss.setConstant(B_ss.rows(), Null);
	for (int r = 0; r < B_ss.rows(); ++r)
		sts.B_ss(r) = B_ss(r, 0);

	if (!mat.Exist("C_ss"))
		throw Exc(S("\n") + t_("Matrix C_ss not found"));
	MatrixXd C_ss;
	mat.Get("C_ss", C_ss);	
	sts.C_ss.setConstant(C_ss.cols(), Null);
	for (int c = 0; c < C_ss.cols(); ++c)
		sts.C_ss(c) = C_ss(0, c);
	
	if (!mat.Exist("Frequencies"))
		throw Exc(S("\n") + t_("Matrix Frequencies not found"));
	MatrixXd ssFrequencies;
	mat.Get("Frequencies", ssFrequencies);	
	sts.ssFrequencies.setConstant(ssFrequencies.cols(), Null);
	for (int c = 0; c < ssFrequencies.cols(); ++c)
		sts.ssFrequencies[c] = ssFrequencies(0, c);

	if (!mat.Exist("FreqRange"))
		throw Exc(S("\n") + t_("Matrix FreqRange not found"));
	MatrixXd ssFreqRange;
	mat.Get("FreqRange", ssFreqRange);	
	sts.ssFreqRange.setConstant(ssFreqRange.cols(), Null);
	for (int c = 0; c < ssFreqRange.cols(); ++c)
		sts.ssFreqRange[c] = ssFreqRange(0, c);

	sts.ssMAE = mat.Get<double>("MAE");	
	
	dt.dimenSTS = true;	
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

	MatrixXd matA(dt.Nf, 1);
	for (int ifr = 0; ifr < dt.Nf; ++ifr)
		matA(ifr, 0) = A_dim(ifr, idf, jdf);
 	if (!mat.Set("A", matA))
 		throw Exc(Format(t_("Problem writing %s to file '%s'"), "A", file));

 	if (!mat.Set<double>("Mu", Ainf_dim(idf, jdf)))
 		throw Exc(Format(t_("Problem writing %s to file '%s'"), "Mu", file));
 		 	
	MatrixXd matB(dt.Nf, 1);
	for (int ifr = 0; ifr < dt.Nf; ++ifr)
		matB(ifr, 0) = B_dim(ifr, idf, jdf);
	if (!mat.Set("B", matB))
 		throw Exc(Format(t_("Problem writing %s to file '%s'"), "B", file));
	
	MatrixXd matw(1, dt.Nf);
	for (int ifr = 0; ifr < dt.Nf; ++ifr)
		matw(0, ifr) = dt.w[ifr];
	if (!mat.Set("w", matw))
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
	
	options.SetStruct<double>("Mode", 0);
	options.SetStruct<double>("Method", 0);
	//options.VarWriteStruct("FreqRangeChoice", "G");
	//options.VarWriteStruct("FreqChoice", "G");
	
	MatrixXd freqRangeChoice(1, 2);
	freqRangeChoice(0, 0) = from;
	freqRangeChoice(0, 1) = to;
	options.SetStruct("FreqRangeChoice", freqRangeChoice);

	MatrixXd freqChoice(1, freqs.size());
	for (int i = 0; i < freqs.size(); ++i) 
		freqChoice(0, i) = freqs[i];
	options.SetStruct("FreqChoice", freqChoice);
	
	UVector<String> optimVars;
	optimVars << "InitCond" << "Tol" << "maxEval" << "maxIter" << "StepTol" << "ThresRel" << "ThresAbs";
	MatVar optim("Optim", 1, 1, optimVars);
	optim.SetStruct<double>("InitCond", 50);
	optim.SetStruct<double>("Tol", 1E-5);
	optim.SetStruct<double>("maxEval", 1E3);
	optim.SetStruct<double>("maxIter", 200);
	optim.SetStruct<double>("StepTol", 1E-6);
	optim.SetStruct<double>("ThresRel", 0.03);
	optim.SetStruct<double>("ThresAbs", 0.1);
	
	options.SetStruct("Optim", optim);
	mat.Set(options);
	
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

