// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"


typedef void (__stdcall *TProgressHandlerProc)(HINSTANCE, LPCWSTR, int *);

class Orca {
public:
	Orca() {}
	~Orca() {
		if (dll.GetHandle() != 0) {
			int lpStatus;
			
			if (wave) {
				DestroyDiffraction(wave, &lpStatus);
				if (lpStatus != 0)
					throwError("DestroyDiffraction");
			}
			
			if (flex) {
				DestroyModel(flex, &lpStatus);
				if (lpStatus != 0)
					throwError("DestroyModel");
			}
			FinaliseLibrary(&lpStatus);
			if (lpStatus != 0)
				throwError("FinaliseLibrary");
		}
	}
	
	bool Init(String _dllFile) {
		dllFile = _dllFile;
		if (!dll.Load(dllFile))
			return false;

		CreateModel  = DLLGetFunction(dll, void, C_CreateModel, (HINSTANCE *handle, HWND hCaller, int *lpStatus));
		DestroyModel = DLLGetFunction(dll, void, C_DestroyModel,(HINSTANCE handle, int *lpStatus));
		LoadData     = DLLGetFunction(dll, void, C_LoadDataW,   (HINSTANCE handle, LPCWSTR wcs, int *lpStatus));
		SaveData     = DLLGetFunction(dll, void, C_SaveDataW,   (HINSTANCE handle, LPCWSTR wcs, int *lpStatus));
		
		CreateDiffraction      = DLLGetFunction(dll, void, C_CreateDiffraction, 	  (HINSTANCE *handle, int *lpStatus));
		DestroyDiffraction     = DLLGetFunction(dll, void, C_DestroyDiffraction,	  (HINSTANCE handle, int *lpStatus));
		LoadDiffractionData    = DLLGetFunction(dll, void, C_LoadDiffractionDataW,    (HINSTANCE handle, LPCWSTR wcs, int *lpStatus));
		SaveDiffractionData    = DLLGetFunction(dll, void, C_SaveDiffractionDataW,    (HINSTANCE handle, LPCWSTR wcs, int *lpStatus));
		CalculateDiffraction   = DLLGetFunction(dll, void, C_CalculateDiffractionW,	  (HINSTANCE handle, TProgressHandlerProc proc, int *lpStatus));
		SaveDiffractionResults = DLLGetFunction(dll, void, C_SaveDiffractionResultsW, (HINSTANCE handle, LPCWSTR wcs, int *lpStatus));
		
		SetModelThreadCount = DLLGetFunction(dll, void, C_SetModelThreadCount, (HINSTANCE handle, int threadCount, int *lpStatus));
		GetModelThreadCount = DLLGetFunction(dll, int, C_GetModelThreadCount,  (HINSTANCE handle, int *lpStatus));
		
		GetLastErrorString = DLLGetFunction(dll, int,  C_GetLastErrorStringW, (LPCWSTR wcs));
		FinaliseLibrary    = DLLGetFunction(dll, void, C_FinaliseLibrary,     (int *lpStatus));
		
		return true;
	}
	
	bool FindInit() {
		UArray<SoftwareDetails> orcadata = GetSoftwareDetails("*OrcaFlex*");	// Get installed versions
		if (orcadata.IsEmpty())
			return false;
		
		int iversion = 0;
		UVector<int> version = orcadata[0].GetVersion();
		for (int i = 1; i < orcadata.size(); ++i) {							// Get the newest version from the installed
			UVector<int> each = orcadata[i].GetVersion();
			if (SoftwareDetails::IsHigherVersion(each, version)) {			// Version are numbers separated by .
				iversion = i;
				version = pick(each);
			}
		}
		String arch;
	#ifdef CPU_64
		arch = "Win64";
	#else
		arch = "Win32";
	#endif
		String path = AppendFileNameX(orcadata[iversion].path, "OrcFxAPI", arch, "OrcFxAPI.dll");	// Assembles the path
		
		return Init(path);
	}
	
	bool IsLoaded()	{return dll;}

	void LoadFlex(String owryml) {
		if (!dll && !FindInit())
			throw Exc("Orca DLL not loaded");
		
		int lpStatus;
		LPCWSTR wcs;
		
		if (flex) {
			DestroyModel(flex, &lpStatus);
			if (lpStatus != 0)
				throwError("DestroyModel");
		}		
		
		CreateModel(&flex, 0, &lpStatus);
		if (lpStatus != 0)
			throwError("CreateModel");
		
		if (!StringToWide(owryml, wcs))
			throwError("StringToWide LoadData");
		LoadData(flex, wcs, &lpStatus);
		if (lpStatus != 0)
			throwError("LoadData");		
	}
	
	void SaveFlex(String owryml) {
		if (!dll && !FindInit())
			throw Exc("Orca DLL not loaded");
		
		int lpStatus;
		LPCWSTR wcs;
					
		if (!StringToWide(owryml, wcs))
			throwError("StringToWide SaveData");
		SaveData(flex, wcs, &lpStatus);
		if (lpStatus != 0)
			throwError("SaveData");
	}
	
	void LoadWave(String owryml) {
		if (!dll && !FindInit())
			throw Exc("Orca DLL not loaded");
			
		if (!wave) {
			int lpStatus;	
			CreateDiffraction(&wave, &lpStatus);
			if (lpStatus != 0)
				throwError("CreateDiffraction");
		}
		
		int lpStatus;
		LPCWSTR wcs;
		
		if (!StringToWide(owryml, wcs))
			throwError("StringToWide LoadDiffractionData");
		LoadDiffractionData(wave, wcs, &lpStatus);
		if (lpStatus != 0)
			throwError("LoadDiffractionData");	
			
		String owd = ForceExt(owryml, ".owd");
		
		if (!StringToWide(owd, wcs))
			throwError("StringToWide LoadDiffractionData SaveDiffractionData");
		SaveDiffractionData(wave, wcs, &lpStatus);
		if (lpStatus != 0)
			throwError("LoadDiffractionData SaveDiffractionData");
	}
	
	static void __stdcall ProgressHandlerProc(HINSTANCE handle, LPCWSTR lpProgress, BOOL *lpCancel) {
		static int lastPerc = -1;
		String msg = WideToString(lpProgress);
		int perc = -1;
		int idp = msg.Find("%");
		if (idp >= 0) {
			int idbp = msg.ReverseFind(" ", idp);
			if (idbp >= 0) 
				perc = ScanInt(msg.Mid(idbp+1, idp-idbp));
		}
		if (lastPerc == perc) 
			return;
		
		lastPerc = perc;
		
		Time et;
		if (perc == 0)
			et = Null;
		else {
			int64 sec = GetSysTime() - startCalc;
			int64 estsec = int64(100*sec/double(perc));
			et = GetSysTime() + estsec;
		}
		*lpCancel = WhenWave(msg, perc, et);
	}

	void RunWave() {
		if (!dll && !FindInit())
			throw Exc("Orca DLL not loaded");
		if (!wave) {
			int lpStatus;
			CreateDiffraction(&wave, &lpStatus);
			if (lpStatus != 0)
				throwError("CreateDiffraction");
		}
		
		int lpStatus;
		startCalc = GetSysTime();
		
		CalculateDiffraction(wave, ProgressHandlerProc, &lpStatus);
		if (lpStatus != 0)
			throwError("CalculateDiffraction");	
	}
	
	void SaveWaveResults(String owryml) {
		if (!dll && !FindInit())
			throw Exc("Orca DLL not loaded");
		if (!wave) {
			int lpStatus;
			CreateDiffraction(&wave, &lpStatus);
			if (lpStatus != 0)
				throwError("CreateDiffraction");
		}
		
		int lpStatus;
		LPCWSTR wcs;
		
		String owr = ForceExt(owryml, ".owr");
		
		if (!StringToWide(owr, wcs))
			throwError("StringToWide SaveData");
		SaveDiffractionResults(wave, wcs, &lpStatus);
		if (lpStatus != 0)
			throwError("SaveDiffractionResults");		
		
		if (GetFileExt(owryml) == ".yml") {
			HINSTANCE temp;
			
			CreateModel(&temp, 0, &lpStatus);
			if (lpStatus != 0)
				throwError("CreateModel");
			
			if (!StringToWide(owr, wcs))
				throwError("StringToWide LoadData");
			LoadData(temp, wcs, &lpStatus);
			if (lpStatus != 0)
				throwError("LoadData");	
			
			if (!StringToWide(owryml, wcs))
				throwError("StringToWide SaveData");
			SaveData(temp, wcs, &lpStatus);
			if (lpStatus != 0)
				throwError("SaveData");
			
			DestroyModel(temp, &lpStatus);
			if (lpStatus != 0)
				throwError("DestroyModel");
		}
	}
	
	void SetThreadCount(int nth) {
		if (!dll && !FindInit())
			throw Exc("Orca DLL not loaded");
		if (!wave) {
			int lpStatus;
			CreateDiffraction(&wave, &lpStatus);
			if (lpStatus != 0)
				throwError("CreateDiffraction");
		}
		
		int lpStatus;

		SetModelThreadCount(wave, nth, &lpStatus);
		if (lpStatus != 0)
			throwError("SetThreadCount");		
	}

	int GetThreadCount() {
		if (!dll && !FindInit())
			throw Exc("Orca DLL not loaded");
		if (!wave) {
			int lpStatus;
			CreateDiffraction(&wave, &lpStatus);
			if (lpStatus != 0)
				throwError("CreateDiffraction");
		}
		
		int lpStatus;

		int nth = GetModelThreadCount(wave, &lpStatus);
		if (lpStatus != 0)
			throwError("GetThreadCount");	
			
		return nth;	
	}
	
	String GetDLLPath() const		{return dllFile;}	
	static Function<bool(String, int, const Time &)> WhenWave;
	static Time startCalc;
	
private:
	Dl dll;
	String dllFile;
	
	HINSTANCE wave = 0, flex = 0;
	
	void (*CreateModel)(HINSTANCE *handle, HWND hCaller, int *lpStatus);
	void (*DestroyModel)(HINSTANCE handle, int *lpStatus);
	void (*LoadData)(HINSTANCE handle, LPCWSTR wcs, int *lpStatus);
	void (*SaveData)(HINSTANCE handle, LPCWSTR wcs, int *lpStatus);
	
	void (*CreateDiffraction)(HINSTANCE *handle, int *lpStatus);
	void (*DestroyDiffraction)(HINSTANCE handle, int *lpStatus);
	void (*LoadDiffractionData)(HINSTANCE handle, LPCWSTR wcs, int *lpStatus);
	void (*SaveDiffractionData)(HINSTANCE handle, LPCWSTR wcs, int *lpStatus);
	void (*CalculateDiffraction)(HINSTANCE handle, TProgressHandlerProc proc, int *lpStatus);
	void (*SaveDiffractionResults)(HINSTANCE handle, LPCWSTR wcs, int *lpStatus);

	void (*SetModelThreadCount)(HINSTANCE handle, int threadCount, int *lpStatus);
	int (*GetModelThreadCount)(HINSTANCE handle, int *lpStatus);
	
	int (*GetLastErrorString)(LPCWSTR wcs);
	
	void (*FinaliseLibrary)(int *lpStatus);
	
	String GetErrorString() {
		LPWSTR wcs;
		int len = GetLastErrorString(NULL);
		Buffer<wchar_t> rw(len);
		wcs = (LPWSTR)rw.begin();
		GetLastErrorString(wcs);
		return WideToString(wcs, len);
	}
	
	void throwError(String where) {
		throw Exc(Format("'%s': %s", where, GetErrorString()));
	}
};


