#include "BEMRosetta.h"
#include "data.brc"
#include <plugin/bz2/bz2.h>
#include <plugin/lzma/lzma.h>
#include <plugin/lz4/lz4.h>
#include <plugin/zstd/zstd.h>

bool Fast::Load(String file, double g) {
	this->file = file;	
	this->name = GetFileTitle(file);
	
	this->g = g;
	
	try {
		String hydroFile;
		bool isFast = false;
		if (GetFileExt(file) == ".dat") {
			isFast = true;
			Print("\n\n" + Format("Loading '%s'", file));
			if (!Load_dat()) {
				Print("\n" + Format("File '%s' not found", file));
				return false;
			}
			hydroFile = AppendFileName(GetFileFolder(file), AppendFileName(hydroFolder, name));
			code = FAST_WAMIT;
		} else {
			Print("\n\n" + Format("Loading '%s'", file));
			hydroFile = file;
			code = WAMIT_1_3;
		}
		
		String file1 = ForceExt(hydroFile, ".1");
		Print("\n" + Format("- Hydrodynamic coefficients A and B file '%s'", GetFileName(file1)));
		if (!Load_1(file1))
			Print(": **Not found**");
		else if (isFast && Nb > 1)
			throw Exc(Format("FAST does not support more than one body in file '%s'", file));
		
		String file3 = ForceExt(hydroFile, ".3");
		Print("\n" + Format("- Diffraction exciting file '%s'", GetFileName(file3)));
		if (!Load_3(file3))
			Print(": **Not found**");
		String fileHST = ForceExt(hydroFile, ".hst");
		Print("\n" + Format("- Hydrostatic restoring file '%s'", GetFileName(fileHST)));
		if (!Load_hst(fileHST))
			Print(": **Not found**");
		
		AfterLoad();
	} catch (Exc e) {
		PrintError("\nError: " + e);
		lastError = e;
		return false;
	}
	
	return true;
}

bool Fast::Load_dat() {
	FileIn in(file);
	if (!in.IsOpen())
		return false;

	Nb = WaveNDir = 1;
	Vo.SetCount(Nb, 0);
	rho = h = len = WaveDirRange = Null;
	
	FieldSplit f;
	while (!in.IsEof()) {
		f.Load(in.GetLine());
		if (f.GetText(1) == "WtrDens") 
			rho = f.GetDouble(0);
		else if (f.GetText(1) == "WtrDpth") 
			h = f.GetDouble(0);
		else if (f.GetText(1) == "WAMITULEN") 
			len = f.GetDouble(0);
		else if (f.GetText(1) == "PtfmVol0") 
			Vo[0] = f.GetDouble(0);
		else if (f.GetText(1) == "WaveNDir") 
			WaveNDir = f.GetInt(0);
		else if (f.GetText(1) == "WaveDirRange") 
			WaveDirRange = f.GetDouble(0);
		else if (f.GetText(1) == "PotFile") { 
			String path = f.GetText(0);
			path.Replace("\"", "");
			hydroFolder = GetFileFolder(path);
			name = GetFileName(path);
		}
	}
	if (IsNull(rho) || IsNull(h) || IsNull(len))
		throw Exc(Format("Wrong format in Fast file '%s'", file));
	
	return true;
}

bool Fast::Load_1(String fileName) {
	w.Clear();
	Nf = Null;
		
	FileIn in(fileName);
	if (!in.IsOpen())
		return false;
	FieldSplit f;
 
 	int64 fpos = 0;
 	while (IsNull(ScanDouble(in.GetLine())) && !in.IsEof())
 		fpos = in.GetPos();
	if (in.IsEof())
		throw Exc("Error in file format");
	
	in.Seek(fpos);
	
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
	
	Nb = 1 + int(maxDof/6);
	Nf = w.GetCount();
	if (Nb == 0 || Nf < 2)
		throw Exc(Format("Wrong format in Wamit file '%s'", file));
	
	if (w[0] > w[1]) {
		readW = false;
		T = pick(w);
		w.SetCount(Nf);	
	} else {
		readW = true;
		T.SetCount(Nf);
	}
	
	A.SetCount(Nf);
	B.SetCount(Nf);	
	if (thereIsAw0)
		Aw0.setConstant(Nb*6, Nb*6, nan(""));
	if (thereIsAwinf)
		Awinf.setConstant(Nb*6, Nb*6, nan(""));

	for (int ifr = 0; ifr < Nf; ++ifr) {
		if (readW)
			T[ifr] = 2*M_PI/w[ifr];
		else
			w[ifr] = 2*M_PI/T[ifr];
		A[ifr].setConstant(Nb*6, Nb*6, nan(""));
	  	B[ifr].setConstant(Nb*6, Nb*6, nan(""));
	}
			
	in.Seek(fpos);
	
	while (!in.IsEof()) {
		f.Load(in.GetLine());
		double freq = f.GetDouble(0);
 		int i = f.GetInt(1) - 1;
 		int j = f.GetInt(2) - 1;
 		double Aij = f.GetDouble(3);
 		
 		if ((freq < 0 && readW) || (freq == 0 && !readW))
			Awinf(i, j) = Aij;
		else if (freq <= 0)
			Aw0(i, j) = Aij;
		else {
			int ifr;
			if (readW)
				ifr = FindIndex(w, freq);
			else
				ifr = FindIndex(T, freq);
		  	A[ifr](i, j) = Aij;    
		  	B[ifr](i, j) = f.GetDouble(4);   	
		}
	}		
	return true;	
}
 
bool Fast::Load_3(String fileName) {
	Nh = Null;
	
	FileIn in(fileName);
	if (!in.IsOpen())
		return false;
	FieldSplit f;
 
 	int64 fpos = 0;
 	while (IsNull(ScanDouble(in.GetLine())) && !in.IsEof())
 		fpos = in.GetPos();
	if (in.IsEof())
		throw Exc("Error in file format");
	
	in.Seek(fpos);
	
	head.Clear();
	while (!in.IsEof()) {
		f.Load(in.GetLine());
		double h = f.GetDouble(1);
		
		FindAdd(head, h);
	}
	
	if (head.GetCount() == 0)
		throw Exc(Format("Wrong format in Wamit file '%s'", file));
	
	if (!IsNull(Nh) && Nh != head.GetCount())
		Print(Format("Number of headings in .dat does not match with .3 (%d != %d)", Nh, head.GetCount()));
	Nh = head.GetCount();
	
	if (abs(head[0]) != abs(head[head.GetCount()-1]))
		Print(Format("Fast requires simetric wave headings. .3 file headings found from %f to %f", head[0], head[head.GetCount()-1])); 
		
	Initialize_Forces();
	
	in.Seek(fpos);
	
	while (!in.IsEof()) {
		f.Load(in.GetLine());
		double freq = f.GetDouble(0);
		int ifr;
		if (readW)
		 	ifr = FindIndex(w, freq);
		else
			ifr = FindIndex(T, freq);
		double h = f.GetDouble(1);
		int ih = FindIndex(head, h);
		int i = f.GetInt(2) - 1;		
		
       	ex.ma[ih](ifr, i) = f.GetDouble(3);
     	ex.ph[ih](ifr, i) = f.GetDouble(4);
        ex.re[ih](ifr, i) = f.GetDouble(5);
        ex.im[ih](ifr, i) = f.GetDouble(6);
	}
		
	return true;
}

bool Fast::Load_hst(String fileName) {
	FileIn in(fileName);
	if (!in.IsOpen())
		return false;
	FieldSplit f;
 
 	int64 fpos = 0;
 	while (IsNull(ScanDouble(in.GetLine())) && !in.IsEof())
 		fpos = in.GetPos();
	if (in.IsEof())
		throw Exc("Error in file format");
	
	in.Seek(fpos);
	
	C.SetCount(Nb);
	for(int ibody = 0; ibody < Nb; ++ibody)
		C[ibody].setConstant(6, 6, 0);

	while (!in.IsEof()) {
		f.Load(in.GetLine());	
		int i = f.GetInt(0) - 1;
		int ib_i = i/6;
		i = i - ib_i*6;
		int j = f.GetInt(1) - 1;
		int ib_j = j/6;
		j = j - ib_j*6;
		if (ib_i == ib_j) 
			C[ib_i](i, j) = f.GetDouble(2);
	}
		
	return true;
}

void Fast::Save(String file, bool isFast) {
	try {
		String hydroFile;
		if (isFast) {
			file = ForceExt(file, ".dat");
			Print("\n\n" + Format("Saving '%s'", file));
			Save_dat(file, true);
			hydroFile = AppendFileName(AppendFileName(GetFileFolder(file), hydroFolder), name);
			DirectoryCreate(AppendFileName(GetFileFolder(file), hydroFolder));
		} else {
			hydroFile = ForceExt(file, ".1");
			Print("\n\n" + Format("Saving '%s'", file));
		}
		
		String file1 = ForceExt(hydroFile, ".1");
		Print("\n" + Format("- Hydrodynamic coefficients A and B file '%s'", GetFileName(file1)));
		Save_1(file1);
		
		String file3 = ForceExt(hydroFile, ".3");
		Print("\n" + Format("- Diffraction exciting file '%s'", GetFileName(file3)));
		Save_3(file3);

		String fileHST = ForceExt(hydroFile, ".hst");
		Print("\n" + Format("- Hydrostatic restoring file '%s'", GetFileName(fileHST)));
		Save_hst(fileHST);
	} catch (Exc e) {
		PrintError("\nError: " + e);
		lastError = e;
	}
}

void Fast::Save_dat(String fileName, bool force) {
	String strFile;
	if (FileExists(fileName)) {
		double lVo = Null, lrho = Null, lh = Null, llen = Null, lWaveDirRange = Null;
		int lWaveNDir = Null;
		
		FileIn in(fileName);
		if (!in.IsOpen())
			return;
	
		FieldSplit f;
		while (!in.IsEof()) {
			f.Load(in.GetLine());
			if (f.GetText(1) == "WtrDens") 
				lrho = f.GetDouble(0);
			else if (f.GetText(1) == "WtrDpth") 
				lh = f.GetDouble(0);
			else if (f.GetText(1) == "WAMITULEN") 
				llen = f.GetDouble(0);
			else if (f.GetText(1) == "PtfmVol0") 
				lVo = f.GetDouble(0);
			else if (f.GetText(1) == "WaveNDir") 
				lWaveNDir = f.GetInt(0);
			else if (f.GetText(1) == "WaveDirRange") 
				lWaveDirRange = f.GetDouble(0);
		}
		in.Close();
		
		if (Nb != 1)
			throw Exc("Number of bodies different to 1 incompatible with Fast");
		if (IsNull(lVo))
			throw Exc(Format("Volume (PtfmVol0) not found in Fast file '%s'", fileName));
		if (IsNull(lrho))
			throw Exc(Format("Density (WtrDens) not found in Fast file '%s'", fileName));
		if (IsNull(lh))
			throw Exc(Format("Water depth (WtrDpth) not found in Fast file '%s'", fileName));
		if (IsNull(llen))
			throw Exc(Format("Length scale (WAMITULEN) not found in Fast file '%s'", fileName));
		if (IsNull(lWaveNDir))
			throw Exc(Format("Number of wave directions (WaveNDir) not found in Fast file '%s'", fileName));
		if (IsNull(lWaveDirRange))
			throw Exc(Format("Range of wave directions (WaveDirRange) not found in Fast file '%s'", fileName));
		
		strFile = LoadFile(fileName);
		int poslf, pos;
			
		if (!force) {
			if (lVo != Vo[0])
				throw Exc(Format("Different volume (%f != %f) in Fast file '%s'", Vo[0], lVo, file));
			if (lrho != rho)
				throw Exc(Format("Different density (%f != %f) in Fast file '%s'", rho, lrho, file));
			if (lh != h)
				throw Exc(Format("Different water depth (%f != %f) in Fast file '%s'", h, lh, file));
			if (llen != len)
				throw Exc(Format("Different length scale (%f != %f) in Fast file '%s'", len, llen, file));
			if (!IsNull(WaveNDir) && lWaveNDir != WaveNDir)
				throw Exc(Format("Different wave headings (%d != %d) in Fast file '%s'", WaveNDir, lWaveNDir, file));
			if (!IsNull(WaveDirRange) && lWaveDirRange != WaveDirRange)
				throw Exc(Format("Headings range do not match (%f != %f) in Fast file '%s'", WaveDirRange, lWaveDirRange, file));
		} else {		
			pos   = strFile.Find("WtrDens");
			poslf = strFile.ReverseFind("\n", pos);
			if (pos < 0 || poslf < 0)
				throw Exc(Format("Bad format parsing Fast file '%s' for WtrDens", file));
			strFile = strFile.Left(poslf+1) + Format("%14>f   ", rho) + strFile.Mid(pos);
			pos   = strFile.Find("WtrDpth");
			poslf = strFile.ReverseFind("\n", pos);
			if (pos < 0 || poslf < 0)
				throw Exc(Format("Bad format parsing Fast file '%s' for WtrDpth", file));
			strFile = strFile.Left(poslf+1) + Format("%14>f   ", h) + strFile.Mid(pos);
			pos   = strFile.Find("WAMITULEN");
			poslf = strFile.ReverseFind("\n", pos);
			if (pos < 0 || poslf < 0)
				throw Exc(Format("Bad format parsing Fast file '%s' for WAMITULEN", file));
			strFile = strFile.Left(poslf+1) + Format("%14>f   ", len) + strFile.Mid(pos);
			pos   = strFile.Find("PtfmVol0");
			poslf = strFile.ReverseFind("\n", pos);
			if (pos < 0 || poslf < 0)
				throw Exc(Format("Bad format parsing Fast file '%s' for PtfmVol0", file));
			strFile = strFile.Left(poslf+1) + Format("%14>f   ", Vo[0]) + strFile.Mid(pos);
			pos   = strFile.Find("WaveNDir");
			poslf = strFile.ReverseFind("\n", pos);
			if (pos < 0 || poslf < 0)
				throw Exc(Format("Bad format parsing Fast file '%s' for WaveNDir", file));
			if (IsNull(WaveNDir))
				strFile = strFile.Left(poslf+1) + Format("%14>d   ", Nh) + strFile.Mid(pos);
			else
				strFile = strFile.Left(poslf+1) + Format("%14>d   ", WaveNDir) + strFile.Mid(pos);
			pos   = strFile.Find("WaveDirRange");
			poslf = strFile.ReverseFind("\n", pos);
			if (pos < 0 || poslf < 0)
				throw Exc(Format("Bad format parsing Fast file '%s' for WaveDirRange", file));
			if (IsNull(WaveDirRange))
				strFile = strFile.Left(poslf+1) + Format("%14>f   ", (head[Nh-1] - head[0])/2) + strFile.Mid(pos);
			else
				strFile = strFile.Left(poslf+1) + Format("%14>f   ", WaveDirRange) + strFile.Mid(pos);
		}
		pos   = strFile.Find("PotFile");
		poslf = strFile.ReverseFind("\n", pos);
		if (pos < 0 || poslf < 0)
			throw Exc(Format("Bad format parsing Fast file '%s' for PotFile", file));
		
		if (hydroFolder.IsEmpty())
			hydroFolder = "HydroData";
		String folder = AppendFileName(hydroFolder, name);
		strFile = strFile.Left(poslf+1) + Format("\"%s\" ", folder) + strFile.Mid(pos);
	} else {
		strFile = ZstdDecompress(hydroDyn, hydroDyn_length);
		strFile.Replace("[WtrDens]", FormatDouble(rho));
		strFile.Replace("[WtrDpth]", FormatDouble(h));
		strFile.Replace("[WAMITULEN]", FormatDouble(len));
		strFile.Replace("[PtfmVol0]", FormatDouble(Vo[0]));
		if (IsNull(WaveNDir))
			strFile.Replace("[WaveNDir]", FormatInt(Nh));
		else
			strFile.Replace("[WaveNDir]", FormatInt(WaveNDir));
		if (IsNull(WaveDirRange))
			strFile.Replace("[WaveDirRange]", FormatDouble((head[Nh-1] - head[0])/2));
		else
			strFile.Replace("[WaveDirRange]", FormatDouble(WaveDirRange));
		strFile.Replace("[PotFile]", Format("\"%s\"", AppendFileName(hydroFolder, name)));
	}
	if (!SaveFile(fileName, strFile))
		throw Exc(Format("Imposible to save file '%s'", file));
}

void Fast::Save_1(String fileName) {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format("Impossible to open '%s'", fileName));
	
	if (IsLoadedAw0()) {
		for (int i = 0; i < Nb*6; ++i)  
			for (int j = 0; j < Nb*6; ++j)
				if (!IsNaN(Aw0(i, j))) 
					out << Format(" %s %5d %5d %s\n", FormatWam(-1), i+1, j+1,
													  FormatWam(Aw0(i, j)/(rho*pow(len, Hydro::GetK_AB(i, j)))));
	}
	if (IsLoadedAwinf()) {
		for (int i = 0; i < Nb*6; ++i)  
			for (int j = 0; j < Nb*6; ++j)
				if (!IsNaN(Awinf(i, j))) 
					out << Format(" %s %5d %5d %s\n", FormatWam(0), i+1, j+1,
													  FormatWam(Awinf(i, j)/(rho*pow(len, Hydro::GetK_AB(i, j)))));
	}
	if (IsLoadedA() && IsLoadedB()) {
		for (int ifr = 0; ifr < Nf; ++ifr)
			for (int i = 0; i < Nb*6; ++i)  
				for (int j = 0; j < Nb*6; ++j)
					if (!IsNaN(A[ifr](i, j)) && !IsNaN(B[ifr](i, j))) 
						out << Format(" %s %5d %5d %s %s\n", FormatWam(T[ifr]), i+1, j+1,
										FormatWam(A[ifr](i, j)/(rho*pow(len, Hydro::GetK_AB(i, j)))), 
										FormatWam(B[ifr](i, j)/(rho*pow(len, Hydro::GetK_AB(i, j))*w[ifr])));
	}
}

void Fast::Save_3(String fileName) {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format("Impossible to open '%s'", fileName));

	if (IsLoadedFex()) {
		for (int ifr = 0; ifr < Nf; ++ifr)
			for (int ih = 0; ih < Nh; ++ih)
				for (int i = 0; i < Nb*6; ++i)
					if (!IsNaN(ex.ma[ih](ifr, i))) {
						double k = g*rho*pow(len, GetK_F(i));
						out << Format(" %s %s %5d %s %s %s %s\n", FormatWam(T[ifr]), FormatWam(head[ih]), i+1,
										FormatWam(ex.ma[ih](ifr, i)/k), FormatWam(ex.ph[ih](ifr, i)),
										FormatWam(ex.re[ih](ifr, i)/k), FormatWam(ex.im[ih](ifr, i)/k));
					}
	}
}

void Fast::Save_hst(String fileName) {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format("Impossible to open '%s'", fileName));

	if (IsLoadedC()) {
		for (int i = 0; i < 6*Nb; ++i)  
			for (int j = 0; j < 6*Nb; ++j) {
				int ib_i = i/6;
				int ii = i - ib_i*6;
				int ib_j = j/6;
				int jj = j - ib_j*6;
				out << Format(" %5d %5d  %s\n", i+1, j+1, FormatWam(C[ib_i](ii, jj)/(g*rho*pow(len, GetK_C(i, j)))));
			}
	}
}
