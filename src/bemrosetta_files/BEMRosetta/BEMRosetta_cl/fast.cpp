#include "BEMRosetta.h"
#include "data.brc"
#include <ScatterDraw/Unpedantic.h>
#include <plugin/bz2/bz2.h>
#include <plugin/lzma/lzma.h>
#include <plugin/lz4/lz4.h>
#include <plugin/zstd/zstd.h>
#include <ScatterDraw/Pedantic.h>

bool Fast::Load(String file, double g) {
	hd().file = file;	
	hd().name = GetFileTitle(file);
	
	hd().g = g;
	
	try {
		String hydroFile;
		if (GetFileExt(file) != ".dat") 
			throw Exc("\n" + Format(t_("File '%s' is not of FAST type"), file));
			
		BEMData::Print("\n\n" + Format(t_("Loading '%s'"), file));
		if (!Load_dat()) 
			throw Exc("\n" + Format(t_("File '%s' not found"), file));

		hydroFile = AppendFileName(GetFileFolder(file), AppendFileName(hydroFolder, hd().name));
		hd().code = Hydro::FAST_WAMIT;
		
		if (!Wamit::Load(ForceExt(hydroFile, ".hst"))) {
			if (!Wamit::Load(ForceExt(hydroFile, ".1"))) 
				return false;
		}
		if (IsNull(hd().Nb))
			return false;
		
		if (hd().Nb > 1)
			throw Exc(Format(t_("FAST does not support more than one body in file '%s'"), file));	
		if (hd().head.IsEmpty())
			throw Exc(t_("No wave headings found in Wamit file"));
		if (abs(hd().head[0]) != abs(hd().head[hd().head.GetCount()-1]))
			throw Exc(Format(t_("FAST requires simetric wave headings. .3 file headings found from %f to %f"), hd().head[0], hd().head[hd().head.GetCount()-1])); 
	} catch (Exc e) {
		BEMData::PrintError("\nError: " + e);
		hd().lastError = e;
		return false;
	}
	
	return true;
}

bool Fast::Load_dat() {
	FileInLine in(hd().file);
	if (!in.IsOpen())
		return false;

	hd().Nb = WaveNDir = 1;
	hd().Vo.SetCount(hd().Nb, 0);
	hd().rho = hd().h = hd().len = WaveDirRange = Null;
	
	FieldSplit f(in);
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
		throw Exc(Format(t_("Wrong format in FAST file '%s'"), hd().file));
	
	return true;
}


void Fast::Save(String file) {
	try {
		String hydroFile;
		file = ForceExt(file, ".dat");
		Save_dat(file, true);
		hydroFile = AppendFileName(AppendFileName(GetFileFolder(file), hydroFolder), hd().name);
		DirectoryCreate(AppendFileName(GetFileFolder(file), hydroFolder));
	
		Wamit::Save(hydroFile);
	} catch (Exc e) {
		BEMData::PrintError(Format("\n%s: %s", t_("Error"), e));
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
		
		FileInLine in(fileName);
		if (!in.IsOpen())
			return;
	
		FieldSplit f(in);
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
			throw Exc("Number of bodies different to 1 incompatible with FAST");
		if (IsNull(lVo))
			throw Exc(Format("Volume (PtfmVol0) not found in FAST file '%s'", fileName));
		if (IsNull(lrho))
			throw Exc(Format("Density (WtrDens) not found in FAST file '%s'", fileName));
		if (IsNull(lh))
			throw Exc(Format("Water depth (WtrDpth) not found in FAST file '%s'", fileName));
		if (IsNull(llen))
			throw Exc(Format("Length scale (WAMITULEN) not found in FAST file '%s'", fileName));
		if (IsNull(lWaveNDir))
			throw Exc(Format("Number of wave directions (WaveNDir) not found in FAST file '%s'", fileName));
		if (IsNull(lWaveDirRange))
			throw Exc(Format("Range of wave directions (WaveDirRange) not found in FAST file '%s'", fileName));
		
		strFile = LoadFile(fileName);
		int poslf, pos;
			
		if (!force) {
			if (lVo != hd().Vo[0])
				throw Exc(Format(t_("Different %s (%f != %f) in FAST file '%s'"), t_("volume"), hd().Vo[0], lVo, hd().file));
			if (lrho != hd().rho)
				throw Exc(Format(t_("Different %s (%f != %f) in FAST file '%s'"), t_("density"), hd().rho, lrho, hd().file));
			if (lh != hd().h)
				throw Exc(Format(t_("Different %s (%f != %f) in FAST file '%s'"), t_("water depth"), hd().h, lh, hd().file));
			if (llen != hd().len)
				throw Exc(Format(t_("Different %s (%f != %f) in FAST file '%s'"), t_("length scale"), hd().len, llen, hd().file));
			if (!IsNull(WaveNDir) && lWaveNDir != WaveNDir)
				throw Exc(Format(t_("Different %s (%d != %d) in FAST file '%s'"), t_("number of wave headings"), WaveNDir, lWaveNDir, hd().file));
			if (!IsNull(WaveDirRange) && lWaveDirRange != WaveDirRange)
				throw Exc(Format(t_("Different %s (%f != %f) in FAST file '%s'"), t_("headings range"), WaveDirRange, lWaveDirRange, hd().file));
		} else {		
			pos   = strFile.Find("WtrDens");
			poslf = strFile.ReverseFind("\n", pos);
			if (pos < 0 || poslf < 0)
				throw Exc(Format(t_("Bad format parsing FAST file '%s' for %s"), hd().file, "WtrDens"));
			strFile = strFile.Left(poslf+1) + Format("%14>f   ", hd().rho) + strFile.Mid(pos);
			pos   = strFile.Find("WtrDpth");
			poslf = strFile.ReverseFind("\n", pos);
			if (pos < 0 || poslf < 0)
				throw Exc(Format(t_("Bad format parsing FAST file '%s' for %s"), hd().file, "WtrDpth"));
			strFile = strFile.Left(poslf+1) + Format("%14>f   ", hd().h) + strFile.Mid(pos);
			pos   = strFile.Find("WAMITULEN");
			poslf = strFile.ReverseFind("\n", pos);
			if (pos < 0 || poslf < 0)
				throw Exc(Format(t_("Bad format parsing FAST file '%s' for %s"), hd().file, "WAMITULEN"));
			strFile = strFile.Left(poslf+1) + Format("%14>f   ", hd().len) + strFile.Mid(pos);
			pos   = strFile.Find("PtfmVol0");
			poslf = strFile.ReverseFind("\n", pos);
			if (pos < 0 || poslf < 0)
				throw Exc(Format(t_("Bad format parsing FAST file '%s' for %s"), hd().file, "PtfmVol0"));
			double hdVo0 = 0;
			if (hd().Vo.GetCount() > 0) 				
				hdVo0 = hd().Vo[0];
			strFile = strFile.Left(poslf+1) + Format("%14>f   ", hdVo0) + strFile.Mid(pos);
			pos   = strFile.Find("WaveNDir");
			poslf = strFile.ReverseFind("\n", pos);
			if (pos < 0 || poslf < 0)
				throw Exc(Format(t_("Bad format parsing FAST file '%s' for %s"), hd().file, "WaveNDir"));
			if (IsNull(WaveNDir))
				strFile = strFile.Left(poslf+1) + Format("%14>d   ", hd().Nh) + strFile.Mid(pos);
			else
				strFile = strFile.Left(poslf+1) + Format("%14>d   ", WaveNDir) + strFile.Mid(pos);
			pos   = strFile.Find("WaveDirRange");
			poslf = strFile.ReverseFind("\n", pos);
			if (pos < 0 || poslf < 0)
				throw Exc(Format(t_("Bad format parsing FAST file '%s' for %s"), hd().file, "WaveDirRange"));
			if (IsNull(WaveDirRange))
				strFile = strFile.Left(poslf+1) + Format("%14>f   ", (hd().head[hd().Nh-1] - hd().head[0])/2) + strFile.Mid(pos);
			else
				strFile = strFile.Left(poslf+1) + Format("%14>f   ", WaveDirRange) + strFile.Mid(pos);
		}
		pos   = strFile.Find("PotFile");
		poslf = strFile.ReverseFind("\n", pos);
		if (pos < 0 || poslf < 0)
			throw Exc(Format(t_("Bad format parsing FAST file '%s' for %s"), hd().file, "PotFile"));
		
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
		throw Exc(Format(t_("Imposible to save file '%s'"), hd().file));
}


