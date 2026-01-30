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
		}
		int j;
		for (j = 0; j < ctypes.size(); ++j) {
			if (argtype.StartsWith(ctypes[j])) 
				break;
		}
		if (j >= ctypes.size())
			throw Exc(Format("Type in argument '%s' not found", arg));
		argTypeId << j;
		argVars   << Trim(arg.Mid(ctypes[j].GetCount()));
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

String GetPythonDeclaration(const String &name, const String &include) {
	static String str;
	const UVector<String> ctypes = {"void", "double **", 									   "int *", 					   "double *", 					      "int",    	   "double",   		  "const char *", 	 "bool",   		  "const double *"}; 
	const UVector<String> ptypes = {"None", "ctypes.POINTER(ctypes.POINTER(ctypes.c_double))", "ctypes.POINTER(ctypes.c_int)", "ctypes.POINTER(ctypes.c_double)", "ctypes.c_int", "ctypes.c_double", "ctypes.c_char_p", "ctypes.c_bool", "np.ctypeslib.ndpointer(dtype=np.float64)"}; 
	const UVector<bool> isPy_C   = {true,   false, 											   false, 						   false, 						      false,    	   false, 	   		  false,          	 false, 	      true};
	const UVector<bool> isC_Py   = {true,   true, 											   true, 						   true, 						      false,    	   false, 	   		  false,          	 false, 	 	  false};
	
	str << "# " << name << " python functions list\n"
		   "import ctypes\n"
		   "import numpy as np\n";
	   
	UVector<String> strIn;
	UVector<String> strOut;
	UVector<String> strSubnames;
	Upp::Index<String> subnamespaces;
	strSubnames.Add();
	subnamespaces << name;
	
	String cleaned = CleanCFromDeclaration(include);
					
	UVector<String> lines = Split(cleaned, "\n");
	for (const auto &line : lines) {
		int pospar = line.Find("(");
		String function, outputType;
		
		for (int i = 0; i < ctypes.size(); ++i) {
			outputType = ctypes[i];
			if (line.StartsWith(outputType)) {
				function = Trim(line.Mid(outputType.GetCount(), pospar - outputType.GetCount()));
				strOut << Format("self.libc.%s.restype = %s", function, ptypes[i]);
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
			if (nextIsIntp) {
				String strdim = argVars[i];
				int dim = 1;
				int pos = strdim.FindAfter("[");
				if (pos >= 0)
					dim = ScanInt(strdim.Mid(pos));
				
				cargs<< Format("ctypes.byref(_data%d), _size%d", idata, idata);
				pre  << Format("        _data%d = ctypes.POINTER(ctypes.c_%s)()\n", idata, prevct)
        			 << Format("        _size%d = (ctypes.c_int * %d)()\n", idata, dim);
        		
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
        			mults << Format("_size%d[%d]", idata, idim);
        			dims  << Format("_size%d[%d]", idata, idim);
        		}
        		if (dim == 1)
        			dims << ",";	// This forces to be a tuple of 1 element...
        		post << Format("        if %s == 0:\n", mults)
        			 << Format("            return np.empty((%s), dtype=%s)\n", dims, nptype);
				post << Format("        %s = np.ctypeslib.as_array(_data%d, shape=(%s))\n", argVars[i-1], idata, dims);
				nextIsIntp = false;
			} else if (ctp.Find("**") > 0) {
        		if (!returns.IsEmpty())
        			 returns << ", ";
        		returns << argVars[i];
        		nextIsIntp = true;
        		prevct = ct;
			} else if (ctp.Find("*") > 0 && ct != "const char") {
				cargs << Format("ctypes.byref(%s)", argVars[i]);
				pre  << Format("        %s = ctypes.c_%s()\n", argVars[i], ct);
        		if (!returns.IsEmpty())
        			 returns << ", ";
        		returns << argVars[i] << ".value";
			} else {
				if (i > 0 && isC_Py[argTypeId[i-1]]) {
					cargs << Format("ctypes.byref(_size%d)", idata);
					idata++;
				} else if (i > 0 && isPy_C[argTypeId[i-1]])
					cargs << Format("len(%s)", argVars[i-1]);
				else {
					if (ctypes[argTypeId[i]] == "const char *")
						cargs << Format("%s.encode('UTF-8')", argVars[i]);
					else
						cargs << argVars[i];
					pargs << argVars[i];
				}
			}
		}
		strIn << Format("self.libc.%s.argtypes = [%s]", function, ToArgs(pargTypes));
				
		String fname = function;
		fname.Replace("DLL_", "");
		
		int pos = fname.Find("_");
		int idsubname = 0;
		if (pos > 0) {
			String subname = fname.Left(pos);
			idsubname = subnamespaces.Find(subname);
			if (idsubname < 0) {
				idsubname = subnamespaces.size();
				subnamespaces << subname;
				strSubnames.Add();
			}
			fname = fname.Mid(pos+1);
		}
		
		strSubnames[idsubname] << "    def " << fname << "(self";
		if (!pargs.IsEmpty())
			strSubnames[idsubname] << ", " << ToArgs(pargs);
		strSubnames[idsubname] << "):\n";
		
		strSubnames[idsubname] << pre ;
		String sret = outputType != "void" ? "_ret = " : "";
		strSubnames[idsubname] << Format("        %s%s%s(%s)\n", sret, "self.libc.", function, ToArgs(cargs));
		strSubnames[idsubname] << "        self._raise_if_error()\n";
		if (!post.IsEmpty()) 
			strSubnames[idsubname] << post;
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
			strSubnames[idsubname] << "        return " << returns << "\n";
		
		strSubnames[idsubname] << "\n";
	}
	
	strSubnames[0] <<
		"    def _raise_if_error(self):\n"
		"        err_ptr = self.libc.DLL_GetLastError()\n"
		"        if err_ptr:\n"
		"            msg = err_ptr.decode('UTF-8', errors=\"replace\")\n"
		"            raise RuntimeError(msg)\n\n";

	strSubnames[0].Insert(0, "\n\n");
	strSubnames[0].Insert(0, "        self.Init()\n");
	strSubnames[0].Insert(0, "\n");
	for (int i = strSubnames.size()-1; i > 0; --i)
		strSubnames[0].Insert(0, Format("        self.%s = _%s(self.libc, self._raise_if_error)\n", subnamespaces[i], subnamespaces[i]));
	
	for (int i = strIn.size()-1; i >= 0; --i) {
		strSubnames[0].Insert(0, "        " << strOut[i] << "\n\n");
		strSubnames[0].Insert(0, "        " << strIn[i] << "\n");
	}
		
	for (int i = 0; i < strSubnames.size(); ++i) {
		String sinit = "class " << S(i == 0 ? "" : "_") << subnamespaces[i] << ":\n";
		if (i == 0)
			sinit <<
		   		"    def __init__(self, path_dll):\n"
		   		"        self.libc = ctypes.CDLL(path_dll)\n\n";
		else
			sinit <<
			    "    def __init__(self, lib, raise_if_error):\n"
        		"        self.libc = lib\n"
        		"        self._raise_if_error = raise_if_error\n\n";
        
		strSubnames[i].Insert(0, sinit);
	}
	
	for (int i = 0; i < strSubnames.size(); ++i) 
		str << "\n" << strSubnames[i];
	
	return str = Trim(str);	
}