// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"


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
