// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"
#include <random>

String GetSpaces(int num) {
	ASSERT(num >= 0);
	return String(' ', num);
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
			posIni = 0;		// It is the first line
			//throw Exc(Format("Bad format parsing FAST file for variable %s", varName));
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
			//nright = navail/2;
			//nleft = navail - nright;
			nleft = 0;									// Values tend to be at the left
			nright = navail - nleft;
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

void GetFASTMatrixIds(const String &strFile, String var, int row, int col, int &posIni, int &posEnd) {
	int id;
	if ((id = strFile.Find(var)) < 0)
		throw Exc(Format(t_("Wrong variable '%s' in GetMatrixIds"), var));
	if ((id = strFile.ReverseFindAfter("\n", id)) < 0)
		throw Exc(Format(t_("Problem reading variable '%s' in GetMatrixIds"), var));
	
	for (int i = 0; i < row; ++i) {
		if ((id = strFile.FindAfter("\n", id)) < 0)
			throw Exc(Format(t_("Problem reading variable '%s' row %d in GetMatrixIds"), var, row));
	}
	for (int ic = 0; ic <= col; ++ic) {
		posIni = id;
		while (id < strFile.GetCount() && IsSpace(strFile[id]))
			id++;
		while (id < strFile.GetCount() && !IsSpace(strFile[id]))
			id++;
		posEnd = id;
	}
	if (id == strFile.GetCount())
		posEnd = id;
		//throw Exc(Format(t_("Problem reading variable '%s' col %d in GetFASTMatrixIds"), var, col));				
}

double GetFASTMatrixVal(const String &strFile, String var, int row, int col) {
	int posIni, posEnd;
	GetFASTMatrixIds(strFile, var, row, col, posIni, posEnd);
	
	String data = strFile.Mid(posIni, posEnd-posIni);
	double ddata = ScanDouble(data);
	if (!IsNum(ddata))
		throw Exc(Format(t_("Problem reading variable '%s' in GetFASTMatrixVal %d, %d"), var, row, col));
	return ddata;
}

Eigen::MatrixXd GetFASTMatrix(const String &strFile, String var, int rows, int cols) {
	Eigen::MatrixXd ret(rows, cols);
	
	for (int r = 0; r < rows; ++r)
		for (int c = 0; c < cols; ++c)
			ret(r, c) = GetFASTMatrixVal(strFile, var, r, c);
			
	return ret;
}
		
UVector<UVector<String>> GetFASTArray(const String &strFile, String var, String paragraph) {
	UVector<UVector<String>> ret;
	int posIni, pos;
	if (!GetFASTVarLine(strFile, var, paragraph, posIni, pos, 0))
		return ret;
	
	int num = ScanInt(strFile.Mid(posIni, pos - posIni));	
	
	for (int i = 0; i < 3; ++i) {
		pos = strFile.FindAfter("\n", pos);
		if (pos < 0)
			throw Exc(Format(t_("Problem reading variable '%s.%s' in GetFASTArray"), paragraph, var));
	}
	for (int i = 0; i < num; ++i) {
		UVector<String> &str = ret.Add();
		int npos = strFile.FindAfter("\n", pos);
		if (npos < 0)
			throw Exc(Format(t_("Problem reading variable '%s.%s' in GetFASTArray"), paragraph, var));
		String line = strFile.Mid(pos, npos-pos);
		line.Replace("\t", " ");
		line = Trim(line);
		str = Split(line, ' ', true);
		pos = npos;
	}
	return ret;
}
