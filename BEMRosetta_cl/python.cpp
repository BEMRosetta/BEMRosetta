// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2026, the BEMRosetta author and contributors
#include "BEMRosetta.h"


static void ListArgsCFunction(const String &strargs, const UVector <String> &ctypes, 
						UVector<int> &argTypeId, UVector<String> &argVars) {
	UVector<String> args = Split(strargs, ",");
	argTypeId.Clear();
	argVars.Clear();
	for (auto &arg : args) {
		arg = Trim(arg);
		String argtype = arg;
		// Search [] to convert int dim[3] into int *
		if (arg.Find("[") > 0) {
			if (argtype.StartsWith("int"))
				argtype = "int *";	
			else if (argtype.StartsWith("const int"))
				argtype = "const int *";	
		}
		int j;
		for (j = 0; j < ctypes.size(); ++j) {
			if (argtype.StartsWith(ctypes[j])) 
				break;
		}
		if (j >= ctypes.size())
			throw Exc(F("Type in argument '%s' not found", arg));
		argTypeId << j;
		if (argtype == "int *")
			argVars << Trim(arg.Mid(4));
		else if (argtype == "const int *")
			argVars << Trim(arg.Mid(10));
		else
			argVars << Trim(arg.Mid(ctypes[j].GetCount()));
	}
}

static String ToArgs(const UVector<String> &args) {
	String ret;
	for (int i = 0; i < args.size(); ++i) {
		if (i > 0)
			ret << ", ";
		ret << args[i];
	}
	return ret;	
}

String GetPythonDeclaration(const String &name, const String &prefix, const String &include) {
	const UVector<String> ctypes = {"void", "double **", 									   "int *", 					   "double *", 					      "int",    	  "double",   		  "const char *", 	"bool",   		 "const double *",					"const int *"}; 
	const UVector<String> ptypes = {"None", "ctypes.POINTER(ctypes.POINTER(ctypes.c_double))", "ctypes.POINTER(ctypes.c_int)", "ctypes.POINTER(ctypes.c_double)", "ctypes.c_int", "ctypes.c_double", "ctypes.c_char_p", "ctypes.c_bool", "ctypes.POINTER(ctypes.c_double)", "ctypes.POINTER(ctypes.c_int)"}; 
	const UVector<bool> isPy_C   = {true,   false, 											   false, 						   false, 						      false,    	  false, 	   		  false,          	false, 	         true, 								false};
	const UVector<bool> isC_Py   = {true,   true, 											   true, 						   true, 						      false,    	  false, 	   		  false,          	false, 	 	     false, 							true};
	
	String str;
	
	str << "# " << name << " python functions list\n"
		   "import os\n"
		   "import ctypes\n"
		   "import numpy as np\n\n";
	 
	UVector<String> strIn;
	UVector<String> strOut;
	UVector<String> strSubNames;
	strSubNames.Add();
	UVector<int> strSubIds;
	strSubIds << 0;
	Upp::Index<String> subnamespaces, subnamespaces_name;
	subnamespaces << name;
	subnamespaces_name << name;
	
	String cleaned = CleanCFromDeclaration(include);
					
	UVector<String> lines = Split(cleaned, "\n");
	for (const auto &line : lines) {
		int pospar = line.Find("(");
		String function, outputType;
		
		for (int i = 0; i < ctypes.size(); ++i) {
			outputType = ctypes[i];
			if (line.StartsWith(outputType)) {
				function = Trim(line.Mid(outputType.GetCount(), pospar - outputType.GetCount()));
				strOut << F("self.libc.%s.restype = %s", function, ptypes[i]);
				break;
			} 
		}
		
		if (function.IsEmpty())
			continue;
			
		int posparout = line.Find(")");
		String strargs = line.Mid(pospar+1, posparout - pospar-1);
		
		UVector<int> argTypeId;
		UVector<String> argVars;
		ListArgsCFunction(strargs, ctypes, argTypeId, argVars);
							
		UVector<String> pargs, cargs, pargTypes;
		String pre, post, returns;
		int idata = 0;
		bool nextIsIntp = false;
		String prevct;
		for (int i = 0; i < argTypeId.size(); ++i) {
			String ctp = ctypes[argTypeId[i]];
			String ct = ctp;
			ct.Replace("*", "");
			ct = Trim(ct);
			String ptp = ptypes[argTypeId[i]];
			pargTypes << ptp;
			String var = argVars[i];
			if (nextIsIntp) {
				String strdim = var;
				int dim = 1;
				int pos = strdim.FindAfter("[");
				if (pos >= 0)
					dim = ScanInt(strdim.Mid(pos));
				
				cargs<< F("ctypes.byref(_data%d), _size%d", idata, idata);
				pre  << F("        _data%d = ctypes.POINTER(ctypes.c_%s)()\n", idata, prevct)
        			 << F("        _size%d = (ctypes.c_int * %d)()\n", idata, dim);
        		
        		String nptype;
        		if (prevct == "double")
        			nptype = "np.float64";
        		else if (prevct == "float")
        			nptype = "np.float32";
        		else if (prevct == "int")
        			nptype = "np.int64";
        		else if (prevct == "int32")
        			nptype = "np.int32";
        		
        		String dims, mults;
        		for (int idim = 0; idim < dim; ++idim) {
        			if (!dims.IsEmpty()) {
        				mults << "*";
        				dims << ", ";
        			}
        			mults << F("_size%d[%d]", idata, idim);
        			dims  << F("_size%d[%d]", idata, idim);
        		}
        		if (dim == 1)
        			dims << ",";	// This forces to be a tuple of 1 element...
        		post << F("        if %s == 0:\n", mults)
        			 << F("            return np.empty((%s), dtype=%s)\n", dims, nptype);
				post << F("        %s = np.ctypeslib.as_array(_data%d, shape=(%s))\n", argVars[i-1], idata, dims);
				nextIsIntp = false;
			} else if (ctp.Find("**") > 0) {
        		if (!returns.IsEmpty())
        			 returns << ", ";
        		returns << var;
        		nextIsIntp = true;
        		prevct = ct;
			} else if (ctp.Find("*") > 0 && ct.Find("const") < 0) {
				cargs << F("ctypes.byref(%s)", var);
				pre  << F("        %s = ctypes.c_%s()\n", var, ct);
        		if (!returns.IsEmpty())
        			 returns << ", ";
        		returns << var << ".value";
			} else if (i > 0 && isC_Py[argTypeId[i-1]]) {
				cargs << F("ctypes.byref(_size%d)", idata);
				idata++;
			} else if (i > 0 && isPy_C[argTypeId[i-1]]) {
				if (argVars[i].Find("[2]") > 0)
					cargs << F("(ctypes.c_int * 2)(*%s.shape)", argVars[i-1]);
				else
					cargs << F("%s.size", argVars[i-1]);
				pargs << argVars[i-1];
			} else if (ctypes[argTypeId[i]] == "const double *") {
				cargs << F("%s.ctypes.data_as(ctypes.POINTER(ctypes.c_double))", var);
				if (argVars[i+1].Find("[2]") > 0) {
					pre << F("        %s = np.asarray(%s, dtype=np.float64)\n"
						   		  "        if %s.ndim != 2:\n"
	    				   		  "            raise ValueError('Function expects a 2D array')\n"
	    				   		  "        %s = np.ascontiguousarray(%s)\n", var, var, var, var, var);
				} else {
					pre << F("        %s = np.asarray(%s, dtype=np.float64)\n"
						   		  "        if %s.ndim != 1:\n"
	    				   		  "            raise ValueError('Function expects a 1D array')\n", var, var, var);
				}
			} else if (ctypes[argTypeId[i]] == "const char *") {
				cargs << F("os.fspath(%s).encode('utf-8')", var);
				pargs << var;
			} else {
				cargs << var;
				pargs << var;
			}
		}
		strIn << F("self.libc.%s.argtypes = [%s]", function, ToArgs(pargTypes));
				
		String fname = function;
		fname.Replace(prefix + "_", "");
		
		UVector<String> fnames = Split(fname, '_');
		if (fnames.IsEmpty())
			throw Exc(F("Wrong function '%s'", fname));
		
		if (fnames[0] == "")
			fnames.Remove(0);
		
		fname = Last(fnames);
		String subname;
		int idsubname = 0;
		for (int iname = 0; iname < fnames.size()-1; ++iname) {
			String parent = subname;
			subname << fnames[iname];
			idsubname = subnamespaces.Find(subname);
			if (idsubname < 0) {
				idsubname = subnamespaces.size();
				subnamespaces << subname;
				subnamespaces_name << fnames[iname];
				strSubNames.Add();
				int idParent = subnamespaces.Find(parent);
				strSubIds << (idParent < 0 ? 0 : idParent);
			}
		}
		/*
		int pos1 = fname.Find("_");
		int idsubname = 0;
		if (pos1 > 0) {
			String subname = fname.Left(pos1);
			idsubname = subnamespaces.Find(subname);
			if (idsubname < 0) {
				idsubname = subnamespaces.size();
				subnamespaces << subname;
				strSubnames.Add();
			}
			fname = fname.Mid(pos1+1);
		}
		*/
		strSubNames[idsubname] << "    def " << fname << "(self";
		if (!pargs.IsEmpty())
			strSubNames[idsubname] << ", " << ToArgs(pargs);
		strSubNames[idsubname] << "):\n";
		
		strSubNames[idsubname] << pre ;
		String sret = outputType != "void" ? "_ret = " : "";
		strSubNames[idsubname] << F("        %s%s%s(%s)\n", sret, "self.libc.", function, ToArgs(cargs));
		strSubNames[idsubname] << "        self._raise_if_error()\n";
		if (!post.IsEmpty()) 
			strSubNames[idsubname] << post;
		else if (outputType != "void") {
			String ret;
			if (outputType == "const char *")
				ret = "_ret.decode('UTF-8', errors=\"replace\")";
			else	
				ret = "_ret";
			if (!returns.IsEmpty())
				returns.Insert(0, ", ");
			returns.Insert(0, ret);
		}
		if (!returns.IsEmpty())
			strSubNames[idsubname] << "        return " << returns << "\n";
		
		strSubNames[idsubname] << "\n";
	}
	
	strSubNames[0] <<
		F("    def _raise_if_error(self):\n"
			   "        err_ptr = self.libc.%s_GetLastError()\n"
			   "        if err_ptr:\n"
			   "            msg = err_ptr.decode('UTF-8', errors=\"replace\")\n"
			   "            raise RuntimeError(msg)\n\n", prefix);

	strSubNames[0].Insert(0, "\n");
	strSubNames[0].Insert(0, "        self.Init()\n");
	//strSubNames[0].Insert(0, "\n");
	for (int i = strSubNames.size()-1; i > 0; --i)
		strSubNames[strSubIds[i]].Insert(0, F("        self.%s = _%s(self.libc, self._raise_if_error)\n", subnamespaces_name[i], subnamespaces[i]));
	
	for (int i = strIn.size()-1; i >= 0; --i) {
		strSubNames[0].Insert(0, "        " << strOut[i] << "\n\n");
		strSubNames[0].Insert(0, "        " << strIn[i] << "\n");
	}
		
	for (int i = 0; i < strSubNames.size(); ++i) {
		String sinit = "class " << F(i == 0 ? "" : "_") << subnamespaces[i] << ":\n";
		if (i == 0)
			sinit <<
		   		"    def __init__(self, path_dll):\n"
		   		"        self.libc = ctypes.CDLL(path_dll)\n\n";
		else
			sinit <<
			    "    def __init__(self, lib, raise_if_error):\n"
        		"        self.libc = lib\n"
        		"        self._raise_if_error = raise_if_error\n\n";
        
		strSubNames[i].Insert(0, sinit);
	}
	
	for (int i = 0; i < strSubNames.size(); ++i) 
		str << /*"\n" << */strSubNames[i];
	
	return str = Trim(str);	
}