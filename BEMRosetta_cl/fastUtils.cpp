#include "BEMRosetta.h"
#include "BEMRosetta_int.h"


String GetSpaces(int num) {return String(' ', num);}

class FASTFiles {
public:
	void Load(String file) {
		String path = GetFileFolder(file);
		
		fast.fileName = file;
		elastodyn.fileName = AppendFileName(path, fast.GetString("EDFile"));
		hydrodyn.fileName = AppendFileName(path, fast.GetString("HydroFile"));
	}
	void Save() {
		fast.Save();
		elastodyn.Save();
		hydrodyn.Save();
	}
	
private:
	class File {
	public:
		String fileName;
		String fileText;
		bool isChanged;
		
		void Save() const {
			if (!isChanged)
				return;
			if (!SaveFile(fileName, fileText))
				throw Exc(Format(t_("Impossible to save file '%s'"), fileName));
		}
		
		String GetString(String var) {
			if (fileText.IsEmpty()) {
				fileText = LoadFile(fileName);
				if (fileText.IsEmpty())
					throw Exc(Format(t_("Impossible to read file '%s'"), fileName));
			}
			String res;
			Vector<String> vars = Split(var, "/");
			if (vars.size() == 2)
				res = GetFASTVar(fileText, vars[1], vars[0]);
			else if (vars.size() == 1)		
				res = GetFASTVar(fileText, vars[0]);
			else
				throw Exc(Format(t_("Wrong variable '%s' in GetString"), var));
			
			if (res[0] == '\"')		// Remove quotes
				res = res.Mid(1);
			if (res[res.GetCount()-1] == '\"')
				res = res.Left(res.GetCount()-1);
			return res;
		}
		double GetDouble(String var) {
			double ddata = ScanDouble(GetString(var));
			if (IsNull(ddata))
				throw Exc(Format(t_("Wrong variable '%s' in GetDouble"), var));
			return ddata;
		}
		bool GetBool(String var) {
			String data = ToLower(GetString(var));
			if (data == "true")
				return true;
			if (data == "false")
				return true;
			int idata = ScanInt(data);
			if (idata == 1)
				return true;
			if (idata == 0)
				return false;
			throw Exc(Format(t_("Wrong variable '%s' in GetBool"), var));
		}
		double GetMatrix(String var, int row, int col) {
			int posIni, posEnd;
			GetMatrixIds(var, row, col, posIni, posEnd);
			
			String data = fileText.Mid(posIni, posEnd-posIni);
			double ddata = ScanDouble(data);
			if (IsNull(ddata))
				throw Exc(Format(t_("Problem reading variable '%s' in GetMatrix %d, %d"), var, row, col));
			return ddata;
		}
		void SetMatrix(String var, int row, int col, double val) {
			int posIni, posEnd;
			GetMatrixIds(var, row, col, posIni, posEnd);
			
			int delta = posEnd-posIni-1;
			
			fileText = fileText.Left(posIni) + S(" ") + FormatDoubleSize(val, delta) + fileText.Mid(posEnd);
		}
		
		void SetString(String var, String val) {
			val = S("\"") + val + S("\"");
			SetString0(var, val);
		}
		void SetDouble(String var, double val, int numDec = 4) {
			SetString0(var, FormatDouble(val, 4, FD_CAP_E));
		}
		void SetBool(String var, bool val) {
			SetString0(var, val ? "True" : "False");
		}
		
	private:
		void SetString0(String var, String val) {
			if (fileText.IsEmpty()) {
				fileText = LoadFile(fileName);
				if (fileText.IsEmpty())
					throw Exc(Format(t_("Impossible to read file '%s'"), fileName));
			}
			Vector<String> vars = Split(var, "/");
			if (vars.size() == 2) {
				isChanged = true;
				return SetFASTVar(fileText, vars[1], val, vars[0]);
			} else if (vars.size() == 1) {
				isChanged = true;	
				return SetFASTVar(fileText, vars[0], val);
			} else
				throw Exc(Format(t_("Wrong variable '%s' in SetString"), var));
		}
		void GetMatrixIds(String var, int row, int col, int &posIni, int &posEnd) {
			if (fileText.IsEmpty()) {
				fileText = LoadFile(fileName);
				if (fileText.IsEmpty())
					throw Exc(Format(t_("Impossible to read file '%s'"), fileName));
			}
			int id;
			if ((id = fileText.Find(var)) < 0)
				throw Exc(Format(t_("Wrong variable '%s' in GetMatrixIds"), var));
			if ((id = fileText.ReverseFindAfter("\n", id)) < 0)
				throw Exc(Format(t_("Problem reading variable '%s' in GetMatrixIds"), var));
			
			for (int i = 0; i < row; ++i) {
				if ((id = fileText.FindAfter("\n", id)) < 0)
					throw Exc(Format(t_("Problem reading variable '%s' row %d in GetMatrixIds"), var, row));
			}
			for (int ic = 0; ic <= col; ++ic) {
				posIni = id;
				while (id < fileText.GetCount() && IsSpace(fileText[id]))
					id++;
				while (id < fileText.GetCount() && !IsSpace(fileText[id]))
					id++;
				posEnd = id;
			}
			if (id == fileText.GetCount())
				throw Exc(Format(t_("Problem reading variable '%s' col %d in GetMatrixIds"), var, col));				
		}
	};

public:
	File fast, elastodyn, hydrodyn;
};

void TestFast() {
	FASTFiles data;
	
	data.Load("C:\\Users\\0203853\\Desktop\\Master 1350\\_fst.fst");
	
	double HubMass = data.elastodyn.GetDouble("HubMass");
	double HubIner = data.elastodyn.GetDouble("HubIner");
	
	double dat = data.hydrodyn.GetMatrix("AddBLin", 3, 1);
	data.hydrodyn.SetMatrix("AddBLin", 4, 0, 2.5423423324324345);
	data.hydrodyn.SetMatrix("AddBLin", 4, 1, 23445.5345);
	data.hydrodyn.SetMatrix("AddBLin", 4, 2, 23445454.5345);
	data.hydrodyn.SetMatrix("AddBLin", 4, 3, 2344545445.5345);
	data.hydrodyn.SetMatrix("AddBLin", 4, 4, 2344545445.53434543545);
	data.hydrodyn.SetMatrix("AddBLin", 4, 5, 2344545445332);
	
	data.Save();
}


bool GetFASTVarLine(const String &strFile, String varName, String paragraph, int &posIni, int &pos, int pos0 = 0) {
	pos = pos0;
	if (!paragraph.IsEmpty()) {
		pos = strFile.FindAfter(paragraph, pos);
		if (pos < 0) 
			return false;
	}
	while (true) {
		pos = strFile.Find(varName, pos);
		if (pos < 0) 
			return false;
		
		// Var name is surrounded by spaces
		if (!IsSpace(strFile[pos-1]) || !IsSpace(strFile[pos + varName.GetCount()])) {
			pos += varName.GetCount();
			continue;
		}
		posIni = strFile.ReverseFind("\n", pos);
		if (posIni < 0) 
			throw Exc(Format("Bad format parsing FAST file for variable %s", varName));
		posIni++;	

		String str = Trim(strFile.Mid(posIni, pos - posIni));
		if (str.IsEmpty())
			return true;
		
		// Is a text
		if (str[0] == '"' && str[str.GetCount()-1] == '"') 
			return true;
		
		// Spaces between... varname found, but in a comment
		bool found = true;
		for (int i = 0; i < str.GetCount(); ++i) {
			if (IsSpace(str[i])) {
				found = false;
				break;
			}
		}
		if (found)
			return true;
		pos += varName.GetCount();
	}
	
	return true;
}

void SetFASTVar(String &strFile, String varName, String value, String paragraph) {
	int posIni, pos, pos0 = 0;
	while (GetFASTVarLine(strFile, varName, paragraph, posIni, pos, pos0)) {
		String sactual = strFile.Mid(posIni, pos - posIni);
		sactual.Replace("\t", "    ");
		int posEndValue;
		for (posEndValue = sactual.GetCount()-1; posEndValue >= 0; --posEndValue)
			if (!IsSpace(sactual[posEndValue]))
				break;
		
		int nright, nleft;
		if (posEndValue+1 >= value.GetCount()) {		// Value added maintaining right margin
			nright = sactual.GetCount() - posEndValue - 1;
			nleft = sactual.GetCount() - nright - value.GetCount();
		} else if (sactual.GetCount()-1 > value.GetCount()) {// Value added maintaining varname position
			int navail = sactual.GetCount() - value.GetCount();	
			nright = navail/2;
			nleft = navail - nright;
		} else 											// Value added moving all
			nright = nleft = 1;

		strFile = strFile.Left(posIni) + GetSpaces(nleft) + value + GetSpaces(nright) + strFile.Mid(pos);

		pos0 = pos + varName.GetCount();
	}
}

String GetFASTVar(const String &strFile, String varName, String paragraph) {
	int posIni, pos;
	if (!GetFASTVarLine(strFile, varName, paragraph, posIni, pos, 0))
		return Null;
	return Trim(strFile.Mid(posIni, pos - posIni));	
}
