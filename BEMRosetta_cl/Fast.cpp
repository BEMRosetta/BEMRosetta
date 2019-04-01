#include "BEMRosetta.h"
#include "data.brc"
#include <plugin/bz2/bz2.h>
#include <plugin/lzma/lzma.h>
#include <plugin/lz4/lz4.h>
#include <plugin/zstd/zstd.h>

bool Fast::Load(String file, double g) {
	hd().file = file;	
	hd().name = GetFileTitle(file);
	
	hd().g = g;
	
	try {
		String hydroFile;
		bool isFast = false;
		if (GetFileExt(file) == ".dat") {
			isFast = true;
			hd().Print("\n\n" + Format("Loading '%s'", file));
			if (!Load_dat()) {
				hd().Print("\n" + Format("File '%s' not found", file));
				return false;
			}
			hydroFile = AppendFileName(GetFileFolder(file), AppendFileName(hydroFolder, hd().name));
			hd().code = Hydro::FAST_WAMIT;
		} else {
			hd().Print("\n\n" + Format("Loading '%s'", file));
			hydroFile = file;
			hd().code = Hydro::WAMIT_1_3;
		}
		
		String file1 = ForceExt(hydroFile, ".1");
		hd().Print("\n" + Format("- Hydrodynamic coefficients A and B file '%s'", GetFileName(file1)));
		
		hd().w.Clear();
		hd().Nf = Null;
		hd().Nh = Null;
		
		if (!Load_1(file1))
			hd().Print(": **Not found**");
		else if (isFast && hd().Nb > 1)
			throw Exc(Format("FAST does not support more than one body in file '%s'", file));
		
		String file3 = ForceExt(hydroFile, ".3");
		hd().Print("\n" + Format("- Diffraction exciting file '%s'", GetFileName(file3)));
		if (!Load_3(file3))
			hd().Print(": **Not found**");
		
		String fileHST = ForceExt(hydroFile, ".hst");
		hd().Print("\n" + Format("- Hydrostatic restoring file '%s'", GetFileName(fileHST)));
		if (!Load_hst(fileHST))
			hd().Print(": **Not found**");
		
		String fileRAO = ForceExt(file, ".4");
		hd().Print("\n" + Format("- RAO file '%s'", GetFileName(fileRAO)));
		if (!Load_4(fileRAO))
			hd().Print(": **Not found**");
		
		hd().AfterLoad();
	} catch (Exc e) {
		hd().PrintError("\nError: " + e);
		hd().lastError = e;
		return false;
	}
	
	return true;
}

bool Fast::Load_dat() {
	FileIn in(hd().file);
	if (!in.IsOpen())
		return false;

	hd().Nb = WaveNDir = 1;
	hd().Vo.SetCount(hd().Nb, 0);
	hd().rho = hd().h = hd().len = WaveDirRange = Null;
	
	FieldSplit f;
	while (!in.IsEof()) {
		f.Load(in.GetLine());
		if (f.GetText(1) == "WtrDens") 
			hd().rho = f.GetDouble(0);
		else if (f.GetText(1) == "WtrDpth") 
			hd().h = f.GetDouble(0);
		else if (f.GetText(1) == "WAMITULEN") 
			hd().len = f.GetDouble(0);
		else if (f.GetText(1) == "PtfmVol0") 
			hd().Vo[0] = f.GetDouble(0);
		else if (f.GetText(1) == "WaveNDir") 
			WaveNDir = f.GetInt(0);
		else if (f.GetText(1) == "WaveDirRange") 
			WaveDirRange = f.GetDouble(0);
		else if (f.GetText(1) == "PotFile") { 
			String path = f.GetText(0);
			path.Replace("\"", "");
			hydroFolder = GetFileFolder(path);
			hd().name = GetFileName(path);
		}
	}
	if (IsNull(hd().rho) || IsNull(hd().h) || IsNull(hd().len))
		throw Exc(Format("Wrong format in Fast file '%s'", hd().file));
	
	return true;
}

bool Fast::Load_1(String fileName) {
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
			FindAdd(hd().w, freq);
		
		int dof = f.GetInt(1);
		if (dof > maxDof)
			maxDof = dof-1;
	}
	
	int Nb = 1 + int(maxDof/6);
	if (!IsNull(hd().Nb) && hd().Nb < Nb)
		throw Exc(Format("Number of bodies loaded is lower than found in .4 (%d != %d)", hd().Nb, Nb));
	hd().Nb = Nb;	
	
	int Nf = hd().w.GetCount();
	if (!IsNull(hd().Nf) && hd().Nf != Nf)
		throw Exc(Format("Number of frequencies loaded is different than found in .4 (%d != %d)", hd().Nf, Nf));
	hd().Nf = Nf;
	
	if (hd().Nb == 0 || hd().Nf < 2)
		throw Exc(Format("Wrong format in Wamit file '%s'", hd().file));
	
	if (hd().w[0] > hd().w[1]) {
		readW = false;
		hd().T = pick(hd().w);
		hd().w.SetCount(hd().Nf);	
	} else {
		readW = true;
		hd().T.SetCount(hd().Nf);
	}
	
	hd().A.SetCount(hd().Nf);
	hd().B.SetCount(hd().Nf);	
	if (thereIsAw0)
		hd().Aw0.setConstant(hd().Nb*6, hd().Nb*6, nan(""));
	if (thereIsAwinf)
		hd().Awinf.setConstant(hd().Nb*6, hd().Nb*6, nan(""));

	for (int ifr = 0; ifr < hd().Nf; ++ifr) {
		if (readW)
			hd().T[ifr] = 2*M_PI/hd().w[ifr];
		else
			hd().w[ifr] = 2*M_PI/hd().T[ifr];
		hd().A[ifr].setConstant(hd().Nb*6, hd().Nb*6, nan(""));
	  	hd().B[ifr].setConstant(hd().Nb*6, hd().Nb*6, nan(""));
	}
			
	in.Seek(fpos);
	
	while (!in.IsEof()) {
		f.Load(in.GetLine());
		double freq = f.GetDouble(0);
 		int i = f.GetInt(1) - 1;
 		int j = f.GetInt(2) - 1;
 		double Aij = f.GetDouble(3);
 		
 		if ((freq < 0 && readW) || (freq == 0 && !readW))
			hd().Awinf(i, j) = Aij;
		else if (freq <= 0)
			hd().Aw0(i, j) = Aij;
		else {
			int ifr;
			if (readW)
				ifr = FindIndex(hd().w, freq);
			else
				ifr = FindIndex(hd().T, freq);
			if (ifr < 0)
				throw Exc(Format("Unknown frequency %f", freq));
		
		  	hd().A[ifr](i, j) = Aij;    
		  	hd().B[ifr](i, j) = f.GetDouble(4);   	
		}
	}		
	return true;	
}
 
bool Fast::Load_3(String fileName) {
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
	
	hd().head.Clear();
	while (!in.IsEof()) {
		f.Load(in.GetLine());
		double head = f.GetDouble(1);
		
		FindAdd(hd().head, head);
	}
	
	if (hd().head.GetCount() == 0)
		throw Exc(Format("Wrong format in Wamit file '%s'", hd().file));
	
	if (!IsNull(hd().Nh) && hd().Nh != hd().head.GetCount())
		throw Exc(Format("Number of headings loaded do not match with .3 (%d != %d)", hd().Nh, hd().head.GetCount()));
	hd().Nh = hd().head.GetCount();
	
	if (abs(hd().head[0]) != abs(hd().head[hd().head.GetCount()-1]))
		hd().Print(Format("Fast requires simetric wave headings. .3 file headings found from %f to %f", hd().head[0], hd().head[hd().head.GetCount()-1])); 
		
	hd().Initialize_Forces();
	
	for (int ifr = 0; ifr < hd().Nf; ++ifr) {
		if (readW)
			hd().T[ifr] = 2*M_PI/hd().w[ifr];
		else
			hd().w[ifr] = 2*M_PI/hd().T[ifr];
	}
	
	in.Seek(fpos);
	
	while (!in.IsEof()) {
		f.Load(in.GetLine());
		double freq = f.GetDouble(0);
		int ifr;
		if (readW)
		 	ifr = FindIndex(hd().w, freq);
		else
			ifr = FindIndex(hd().T, freq);
		if (ifr < 0)
			throw Exc(Format("Unknown frequency %f", freq));
		double head = f.GetDouble(1);
		int ih = FindIndex(hd().head, head);
		if (ih < 0)
			throw Exc(Format("Unknown heading %f", head));
			
		int i = f.GetInt(2) - 1;		
		
       	hd().ex.ma[ih](ifr, i) = f.GetDouble(3);
     	hd().ex.ph[ih](ifr, i) = f.GetDouble(4);
        hd().ex.re[ih](ifr, i) = f.GetDouble(5);
        hd().ex.im[ih](ifr, i) = f.GetDouble(6);
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
	
	hd().C.SetCount(hd().Nb);
	for(int ibody = 0; ibody < hd().Nb; ++ibody)
		hd().C[ibody].setConstant(6, 6, 0);

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

bool Fast::Load_4(String fileName) {
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
	hd().head.Clear();
	while (!in.IsEof()) {
		f.Load(in.GetLine());
		double freq = f.GetDouble(0);
		FindAdd(hd().w, freq);
		
		int dof = f.GetInt(2);
		if (dof > maxDof)
			maxDof = dof-1;
		
		double head = f.GetDouble(1);
		
		FindAdd(hd().head, head);
	}

	if (hd().head.GetCount() == 0)
		throw Exc(Format("Wrong format in Wamit file '%s'", hd().file));
	
	if (!IsNull(hd().Nh) && hd().Nh != hd().head.GetCount())
		throw Exc(Format("Number of headings loaded do not match with .4 (%d != %d)", hd().Nh, hd().head.GetCount()));
	hd().Nh = hd().head.GetCount();
	
	int Nb = 1 + int(maxDof/6);
	if (!IsNull(hd().Nb) && hd().Nb < Nb)
		throw Exc(Format("Number of bodies loaded is lower than found in .4 (%d != %d)", hd().Nb, Nb));
	hd().Nb = Nb;
	
	int Nf = hd().w.GetCount();
	if (!IsNull(hd().Nf) && hd().Nf != Nf)
		throw Exc(Format("Number of frequencies loaded is different than found in .4 (%d != %d)", hd().Nf, Nf));
	hd().Nf = Nf;
	
	if (hd().Nb == 0 || hd().Nf < 2)
		throw Exc(Format("Wrong format in Wamit file '%s'", hd().file));
	
	bool readW;
	if (hd().w[0] > hd().w[1]) {
		readW = false;
		hd().T = pick(hd().w);
		hd().w.SetCount(hd().Nf);	
	} else {
		readW = true;
		hd().T.SetCount(hd().Nf);
	}
	
	hd().Initialize_RAO();
	
	for (int ifr = 0; ifr < hd().Nf; ++ifr) {
		if (readW)
			hd().T[ifr] = 2*M_PI/hd().w[ifr];
		else
			hd().w[ifr] = 2*M_PI/hd().T[ifr];
	}
	
	in.Seek(fpos);
	
	while (!in.IsEof()) {
		f.Load(in.GetLine());
		double freq = f.GetDouble(0);
		int ifr;
		if (readW)
		 	ifr = FindIndex(hd().w, freq);
		else
			ifr = FindIndex(hd().T, freq);
		if (ifr < 0)
			throw Exc(Format("Unknown frequency %f", freq));
		double head = f.GetDouble(1);
		int ih = FindIndex(hd().head, head);
		if (ih < 0)
			throw Exc(Format("Unknown heading %f", head));
			
		int i = f.GetInt(2) - 1;		
		
       	hd().rao.ma[ih](ifr, i) = f.GetDouble(3);
     	hd().rao.ph[ih](ifr, i) = f.GetDouble(4);
        hd().rao.re[ih](ifr, i) = f.GetDouble(5);
        hd().rao.im[ih](ifr, i) = f.GetDouble(6);
	}
		
	return true;
}

void Fast::Save(String file, bool isFast) {
	try {
		String hydroFile;
		if (isFast) {
			file = ForceExt(file, ".dat");
			hd().Print("\n\n" + Format("Saving '%s'", file));
			Save_dat(file, true);
			hydroFile = AppendFileName(AppendFileName(GetFileFolder(file), hydroFolder), hd().name);
			DirectoryCreate(AppendFileName(GetFileFolder(file), hydroFolder));
		} else {
			hydroFile = ForceExt(file, ".1");
			hd().Print("\n\n" + Format("Saving '%s'", file));
		}
		
		String file1 = ForceExt(hydroFile, ".1");
		hd().Print("\n" + Format("- Hydrodynamic coefficients A and B file '%s'", GetFileName(file1)));
		Save_1(file1);
		
		String file3 = ForceExt(hydroFile, ".3");
		hd().Print("\n" + Format("- Diffraction exciting file '%s'", GetFileName(file3)));
		Save_3(file3);

		String fileHST = ForceExt(hydroFile, ".hst");
		hd().Print("\n" + Format("- Hydrostatic restoring file '%s'", GetFileName(fileHST)));
		Save_hst(fileHST);
	} catch (Exc e) {
		hd().PrintError("\nError: " + e);
		hd().lastError = e;
	}
}

void Fast::Save_dat(String fileName, bool force) {
	String strFile;
	
	if (hydroFolder.IsEmpty())
		hydroFolder = "HydroData";
	
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
		
		if (hd().Nb != 1)
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
			if (lVo != hd().Vo[0])
				throw Exc(Format("Different volume (%f != %f) in Fast file '%s'", hd().Vo[0], lVo, hd().file));
			if (lrho != hd().rho)
				throw Exc(Format("Different density (%f != %f) in Fast file '%s'", hd().rho, lrho, hd().file));
			if (lh != hd().h)
				throw Exc(Format("Different water depth (%f != %f) in Fast file '%s'", hd().h, lh, hd().file));
			if (llen != hd().len)
				throw Exc(Format("Different length scale (%f != %f) in Fast file '%s'", hd().len, llen, hd().file));
			if (!IsNull(WaveNDir) && lWaveNDir != WaveNDir)
				throw Exc(Format("Different wave headings (%d != %d) in Fast file '%s'", WaveNDir, lWaveNDir, hd().file));
			if (!IsNull(WaveDirRange) && lWaveDirRange != WaveDirRange)
				throw Exc(Format("Headings range do not match (%f != %f) in Fast file '%s'", WaveDirRange, lWaveDirRange, hd().file));
		} else {		
			pos   = strFile.Find("WtrDens");
			poslf = strFile.ReverseFind("\n", pos);
			if (pos < 0 || poslf < 0)
				throw Exc(Format("Bad format parsing Fast file '%s' for WtrDens", hd().file));
			strFile = strFile.Left(poslf+1) + Format("%14>f   ", hd().rho) + strFile.Mid(pos);
			pos   = strFile.Find("WtrDpth");
			poslf = strFile.ReverseFind("\n", pos);
			if (pos < 0 || poslf < 0)
				throw Exc(Format("Bad format parsing Fast file '%s' for WtrDpth", hd().file));
			strFile = strFile.Left(poslf+1) + Format("%14>f   ", hd().h) + strFile.Mid(pos);
			pos   = strFile.Find("WAMITULEN");
			poslf = strFile.ReverseFind("\n", pos);
			if (pos < 0 || poslf < 0)
				throw Exc(Format("Bad format parsing Fast file '%s' for WAMITULEN", hd().file));
			strFile = strFile.Left(poslf+1) + Format("%14>f   ", hd().len) + strFile.Mid(pos);
			pos   = strFile.Find("PtfmVol0");
			poslf = strFile.ReverseFind("\n", pos);
			if (pos < 0 || poslf < 0)
				throw Exc(Format("Bad format parsing Fast file '%s' for PtfmVol0", hd().file));
			double hdVo0 = 0;
			if (hd().Vo.GetCount() > 0) 				
				hdVo0 = hd().Vo[0];
			strFile = strFile.Left(poslf+1) + Format("%14>f   ", hdVo0) + strFile.Mid(pos);
			pos   = strFile.Find("WaveNDir");
			poslf = strFile.ReverseFind("\n", pos);
			if (pos < 0 || poslf < 0)
				throw Exc(Format("Bad format parsing Fast file '%s' for WaveNDir", hd().file));
			if (IsNull(WaveNDir))
				strFile = strFile.Left(poslf+1) + Format("%14>d   ", hd().Nh) + strFile.Mid(pos);
			else
				strFile = strFile.Left(poslf+1) + Format("%14>d   ", WaveNDir) + strFile.Mid(pos);
			pos   = strFile.Find("WaveDirRange");
			poslf = strFile.ReverseFind("\n", pos);
			if (pos < 0 || poslf < 0)
				throw Exc(Format("Bad format parsing Fast file '%s' for WaveDirRange", hd().file));
			if (IsNull(WaveDirRange))
				strFile = strFile.Left(poslf+1) + Format("%14>f   ", (hd().head[hd().Nh-1] - hd().head[0])/2) + strFile.Mid(pos);
			else
				strFile = strFile.Left(poslf+1) + Format("%14>f   ", WaveDirRange) + strFile.Mid(pos);
		}
		pos   = strFile.Find("PotFile");
		poslf = strFile.ReverseFind("\n", pos);
		if (pos < 0 || poslf < 0)
			throw Exc(Format("Bad format parsing Fast file '%s' for PotFile", hd().file));
		
		String folder = AppendFileName(hydroFolder, hd().name);
		strFile = strFile.Left(poslf+1) + Format("\"%s\" ", folder) + strFile.Mid(pos);
	} else {
		strFile = ZstdDecompress(hydroDyn, hydroDyn_length);
		strFile.Replace("[WtrDens]", FormatDouble(hd().rho));
		strFile.Replace("[WtrDpth]", FormatDouble(hd().h));
		strFile.Replace("[WAMITULEN]", FormatDouble(hd().len));
		double hdVo0 = 0;
		if (hd().Vo.GetCount() > 0) 				
			hdVo0 = hd().Vo[0];
		strFile.Replace("[PtfmVol0]", FormatDouble(hdVo0));
		if (IsNull(WaveNDir))
			strFile.Replace("[WaveNDir]", FormatInt(hd().Nh));
		else
			strFile.Replace("[WaveNDir]", FormatInt(WaveNDir));
		if (IsNull(WaveDirRange))
			strFile.Replace("[WaveDirRange]", FormatDouble((hd().head[hd().Nh-1] - hd().head[0])/2));
		else
			strFile.Replace("[WaveDirRange]", FormatDouble(WaveDirRange));
		strFile.Replace("[PotFile]", Format("\"%s\"", AppendFileName(hydroFolder, hd().name)));
	}
	if (!SaveFile(fileName, strFile))
		throw Exc(Format("Imposible to save file '%s'", hd().file));
}

void Fast::Save_1(String fileName) {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format("Impossible to open '%s'", fileName));
	
	if (hd().IsLoadedAw0()) {
		for (int i = 0; i < hd().Nb*6; ++i)  
			for (int j = 0; j < hd().Nb*6; ++j)
				if (!IsNaN(hd().Aw0(i, j))) 
					out << Format(" %s %5d %5d %s\n", FormatWam(-1), i+1, j+1,
													  FormatWam(hd().Aw0(i, j)/(hd().rho*pow(hd().len, Hydro::GetK_AB(i, j)))));
	}
	if (hd().IsLoadedAwinf()) {
		for (int i = 0; i < hd().Nb*6; ++i)  
			for (int j = 0; j < hd().Nb*6; ++j)
				if (!IsNaN(hd().Awinf(i, j))) 
					out << Format(" %s %5d %5d %s\n", FormatWam(0), i+1, j+1,
													  FormatWam(hd().Awinf(i, j)/(hd().rho*pow(hd().len, Hydro::GetK_AB(i, j)))));
	}
	if (hd().IsLoadedA() && hd().IsLoadedB()) {
		for (int ifr = 0; ifr < hd().Nf; ++ifr)
			for (int i = 0; i < hd().Nb*6; ++i)  
				for (int j = 0; j < hd().Nb*6; ++j)
					if (!IsNaN(hd().A[ifr](i, j)) && !IsNaN(hd().B[ifr](i, j))) 
						out << Format(" %s %5d %5d %s %s\n", FormatWam(hd().T[ifr]), i+1, j+1,
										FormatWam(hd().A[ifr](i, j)/(hd().rho*pow(hd().len, Hydro::GetK_AB(i, j)))), 
										FormatWam(hd().B[ifr](i, j)/(hd().rho*pow(hd().len, Hydro::GetK_AB(i, j))*hd().w[ifr])));
	}
}

void Fast::Save_3(String fileName) {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format("Impossible to open '%s'", fileName));

	if (hd().IsLoadedFex()) {
		for (int ifr = 0; ifr < hd().Nf; ++ifr)
			for (int ih = 0; ih < hd().Nh; ++ih)
				for (int i = 0; i < hd().Nb*6; ++i)
					if (!IsNaN(hd().ex.ma[ih](ifr, i))) {
						double k = hd().g*hd().rho*pow(hd().len, Hydro::GetK_F(i));
						out << Format(" %s %s %5d %s %s %s %s\n", FormatWam(hd().T[ifr]), FormatWam(hd().head[ih]), i+1,
										FormatWam(hd().ex.ma[ih](ifr, i)/k), FormatWam(hd().ex.ph[ih](ifr, i)),
										FormatWam(hd().ex.re[ih](ifr, i)/k), FormatWam(hd().ex.im[ih](ifr, i)/k));
					}
	}
}

void Fast::Save_hst(String fileName) {
	FileOut out(fileName);
	if (!out.IsOpen())
		throw Exc(Format("Impossible to open '%s'", fileName));

	if (hd().IsLoadedC()) {
		for (int i = 0; i < 6*hd().Nb; ++i)  
			for (int j = 0; j < 6*hd().Nb; ++j) {
				int ib_i = i/6;
				int ii = i - ib_i*6;
				int ib_j = j/6;
				int jj = j - ib_j*6;
				out << Format(" %5d %5d  %s\n", i+1, j+1, FormatWam(hd().C[ib_i](ii, jj)/(hd().g*hd().rho*pow(hd().len, Hydro::GetK_C(i, j)))));
			}
	}
}
