// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2026, the BEMRosetta author and contributors
#include "BEMRosetta.h"

String CleanCFromDeclaration(const String &include, bool removeSemicolon) {
	String str = include;
	
	str.Replace("	__declspec(dllexport) ", "");
	str.Replace("	L_EXPORT ", "");
	str.Replace("extern \"C\" {", "");
	str.Replace("};", "");
	str.Replace("\r\n\r\n", "\r\n");
	str.Replace(" noexcept", "");
	str.Replace("  ", " ");
	str.Replace("\t", " ");
	str.Replace(" ;", ";");
	str.Replace(" (", "(");
	str.Replace(" )", ")");
	
	if (removeSemicolon) 
		str.Replace(");", ")");
	
	return str;
}

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
	
	String cleaned = CleanCFromDeclaration(include, true);
					
	UVector<String> lines = Split(cleaned, "\n");
	
	bool infunctions = false;
	
	for (String line : lines) {
		Replace(line, "\r", "");
		line = Trim(line);
		if (!infunctions && line.Find("// DLL functions") >= 0)
			infunctions = true;
		else if (infunctions) {
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
	        		returns << var << ".copy()";
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
	}
	
	strSubNames[0] <<
		F("    def _raise_if_error(self):\n"
			   "        err_ptr = self.libc.%s_GetLastError()\n"
			   "        if err_ptr:\n"
			   "            msg = err_ptr.decode('UTF-8', errors=\"replace\")\n"
			   "            raise RuntimeError(msg)\n\n", prefix);

	strSubNames[0].Insert(0, "\n");
	strSubNames[0].Insert(0, "        self.Init()\n");
	
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
		str << strSubNames[i];
	
	return str = Trim(str);	
}

String BMR_CFunctions_List(const String &include, bool isC) {
	String ret;
	
	String sbegin;
	if (isC)
		sbegin = "// C DLL functions";
	else
		sbegin = "// DLL functions";
		
	UVector<String> lines = Split(include, '\n');
	
	bool infunctions = false;
	
	for (String line : lines) {
		Replace(line, "\r", "");
		line = Trim(line);
		if (!infunctions && line.Find(sbegin) >= 0)
			infunctions = true;
		else if (infunctions) {
			if (line.StartsWith("//")) {
				ret << line << "\n";
			} else if (line.StartsWith("L_EXPORT")) {
				line.Replace("L_EXPORT", "");	
				line.Replace("NOEXCEPT", "");	
				line.Replace(";", "");	
				line = Trim(line);
				bool isvoid = line.StartsWith("void");
				int _pos = line.Find("_");
				String retType = Trim(line.Left(_pos));
				line = line.Mid(_pos);
				retType.Replace("void", "int");	
				if (!retType.EndsWith("*"))
					retType << " ";
				ret << retType << line.Mid(1) << "\n\n";
			}
		}
	}
	return ret;
}

struct FuncNode : Moveable<FuncNode> {
    VectorMap<String, FuncNode> children;
    bool isFunc = false;
    String retType;
    String fullName;   			// e.g. "BMR_Bem_Data_Get"
    UVector<String> args;       // Argument list
};

String GetArgs(const UVector<String> &args) {
    String sargs;
    for (int i = 0; i < args.size(); i += 2) {
        if (i > 0)
            sargs << ", ";
        String s = Replace(args[i], "#", args[i+1]);
        sargs << s;
    }
    return sargs;
}

String GetArgNames(const UVector<String> &args) {
    String sargs;
    for (int i = 0; i < args.size(); i += 2) {
        if (i > 0)
            sargs << ", ";
        String sarg = args[i+1];
        sarg.Replace("[1]", "");	sarg.Replace("[2]", "");
        sargs << sarg;
    }
    return sargs;	
}


static String EmitNode(const String& name, const FuncNode& node, int depth, const String& parentClass) {
    String ind(' ', depth*4);
    String s;

    if (node.isFunc) {
        String retType = node.retType;
        if (!retType.EndsWith("*"))
            retType << " ";
        s << ind << "private:\n";
        s << ind << "typedef " << retType << " (*" << name << "_t)(" << GetArgs(node.args) << ");\n";
        s << ind << "#ifdef BEMROSETTA_DYNAMIC\n";
        s << ind << name << "_t _" << name << " = nullptr;\n";
        //s << ind << "#else\n";
        //s << ind << name << "_t _" << name << " = &::" << node.fullName << ";\n";
        s << ind << "#endif\n";
        s << ind << "public:\n";
        s << ind << retType << name << "(" << GetArgs(node.args) << ") {\n";
        s << ind << "#ifdef BEMROSETTA_DYNAMIC\n";
	        s << ind << "    ";
	        if (node.retType != "void")
	        	s << retType << "ret = ";
        	s << "_" << name << "(" << GetArgNames(node.args) << ");\n";
        s << ind << "#else\n";
            s << ind << "    ";
	        if (node.retType != "void")
	        	s << retType << "ret = ";
        	s << node.fullName << "(" << GetArgNames(node.args) << ");\n";
        s << ind << "#endif\n";
        s << ind << "    if (_BMR_GetLastError())\n";
        s << ind << "        throw std::runtime_error(_BMR_GetLastError());\n";
        if (node.retType != "void")
        	s << ind << "    return ret;\n";
        s << ind << "}\n";
    } else {
        String className = name + "_t";
        s << ind << "class " << className << " {\n";
        String fri;
        if (parentClass == "BMR")
            fri = "BEMRosetta";
        else
            fri = parentClass;
        s << ind << "    friend class " << fri << ";\n";
        s << ind << "public:\n";
        for (int i = 0; i < node.children.size(); i++)
            s << EmitNode(node.children.GetKey(i), node.children[i], depth + 1, className);
        s << ind << "#ifdef BEMROSETTA_DYNAMIC\n";
        s << ind << "private:\n";
        s << ind << "    void LoadDllFunction(DLL_HANDLE dll) {\n";
        for (int i = 0; i < node.children.size(); i++) {
            const String& cname = node.children.GetKey(i);
            const FuncNode& child = node.children[i];
            if (child.isFunc)
                s << ind << "        _" << cname << " = (" << cname << "_t)DLL_SYM(dll, \"" << child.fullName << "\");\n";
            else
                s << ind << "        " << cname << ".LoadDllFunction(dll);\n";
        }
        s << ind << "    }\n";
        s << ind << "#endif\n";
        s << ind << "} " << name << ";\n";
    }
    return s;
}

static FuncNode BuildTree(const UVector<String>& retTypes, const UVector<String>& functions, const UVector<UVector<String>>& arguments) {
    FuncNode root;
    for (int i = 0; i < functions.size(); i++) {
        UVector<String> parts = Split(functions[i], '_');
        FuncNode *node = &root;
        for (int p = 0; p < parts.size(); p++) {
            node = &node->children.GetAdd(parts[p]);
            if (p == parts.size() - 1) {
                node->isFunc   = true;
                node->retType  = retTypes[i];
                node->fullName = functions[i];
                node->args     = clone(arguments[i]);
            } else
                node->isFunc   = false;
        }
    }
    return root;
}

static String EmitClassBody(const FuncNode& node, int depth, const String& className) {
    String s;
    for (int i = 0; i < node.children.size(); i++)
        s << EmitNode(node.children.GetKey(i), node.children[i], depth + 1, className);
    return s;
}

static String GenerateHeaderCpp(FuncNode& root) {
    ASSERT(root.children.size() == 1);   // single common prefix, e.g. BMR_
    String className = root.children.GetKey(0);
    const FuncNode& top = root.children[0];

    String s;
    s << "#ifdef _WIN32\n"
      << "    #include <windows.h>\n"
      << "    #define DLL_HANDLE HMODULE\n"
      << "    #define DLL_LOAD(file) LoadLibraryA(file)\n"
      << "    #define DLL_SYM(handle, name) GetProcAddress(handle, name)\n"
      << "    #define DLL_FREE(handle) FreeLibrary(handle)\n"
      << "#else\n"
      << "    #include <dlfcn.h>\n"
      << "    #define DLL_HANDLE void*\n"
      << "    #define DLL_LOAD(file) dlopen(file, RTLD_NOW)\n"
      << "    #define DLL_SYM(handle, name) dlsym(handle, name)\n"
      << "    #define DLL_FREE(handle) dlclose(handle)\n"
      << "#endif\n\n";

    s << "class " << "BEMRosetta" << " {\n"
      << "public:\n"
      << EmitClassBody(top, 1, className);

    s << "#ifdef BEMROSETTA_DYNAMIC\n"
      << "    DLL_HANDLE dll = nullptr;\n\n"
      << "    void LoadDll(const char* file_dll) {\n"
      << "        dll = DLL_LOAD(file_dll);\n"
      << "        if (!dll)\n"
      << "            throw std::runtime_error(\"DLL '\" + std::string(file_dll) + \"' not found\");\n"
      << "        LoadDllFunction(dll);\n";
    for (int i = 0; i < top.children.size(); i++) {
        if (!top.children[i].children.IsEmpty())
        	s << "        " << top.children.GetKey(i) << ".LoadDllFunction(dll);\n";
    }
    s << "    }\n\n"
      << "    BEMRosetta(const char* file_dll) {LoadDll(file_dll);	Init();}\n"
      << "    ~BEMRosetta() {if(dll) DLL_FREE(dll);}\n";
    
	String ind(' ', 4); 
    s << ind << "private:\n";
    s << ind << "    void LoadDllFunction(DLL_HANDLE dll) {\n";
    for (int i = 0; i < top.children.size(); i++) {
        const String& cname = top.children.GetKey(i);
        const FuncNode& child = top.children[i];
        if (child.isFunc)
            s << ind << "        _" << cname << " = (" << cname << "_t)DLL_SYM(dll, \"" << child.fullName << "\");\n";
    }
    s << ind << "        _BMR_GetLastError = _GetLastError;\n";
    s << ind << "    }\n";
 
 	s << "#else\n";
 	s << "    BEMRosetta() {Init();}\n";
    s << "#endif\n";
    s << "};\n";
    return s;
}

UVector<String> SplitCArguments(const String& input) {
    UVector<String> args;
    String current;
    int parens = 0;
    
    for (int i = 0; i < input.GetCount(); i++) {
        char c = input[i];
        
        if (c == ',' && parens == 0) {
            String trimmed = Trim(current);
            if(!trimmed.IsEmpty())
                args << trimmed;
            
            current.Clear();
        } else {
	        if(c == '(')
	            parens++;
	        else if(c == ')')
	            parens--;
	       
	        current << c;
        }
    }
    String trimmed = Trim(current);
    if(!trimmed.IsEmpty())
        args << trimmed;
    
    return args;
}

void GetFunctionsList(const String &include, bool isC, UVector<String> &retTypes, UVector<String> &functions, UVector<UVector<String>> &arguments, String &before, String &declaration, String &after) {
	UVector<String> lines = Split(include, '\n');
	
	String sbegin;
	if (isC)
		sbegin = "// C DLL functions";
	else
		sbegin = "// DLL functions";
	
	int status = 0;
	
	for (String line : lines) {
		line.Replace("\r", "");
		String oline = line;
		line = Trim(line);
		if (status == 1) {
			if (line.StartsWith("L_EXPORT")) {
				line.Replace("L_EXPORT", "");	
				line.Replace("NOEXCEPT", "");	
				line.Replace(";", "");	
				line = Trim(line);
				bool isvoid = line.StartsWith("void");
				int _pos = line.Find("_");
				retTypes << Trim(line.Left(_pos));
				line = line.Mid(_pos);
				int posPar = line.Find("(");
				functions << Trim(line.Left(posPar));
				line = line.Mid(posPar+1);
				int lasPar = line.ReverseFind(")");
				line = line.Left(lasPar);			
				line = Trim(line);
				UVector<String> &aa = arguments.Add();
				UVector<String> args = SplitCArguments(line);
				
				
				for (String a : args) {
					a.Replace("\t", " ");
					a.Replace("  ", " ");
					String variable = a;
					variable.Replace("(", "");	variable.Replace(")", "");	variable.Replace("*", "");
					variable.Replace("[1]", "");variable.Replace("[2]", "");
					variable.Replace("const ", " ");	variable.Replace("int ", " ");	variable.Replace("char ", " ");
					variable.Replace("double ", " ");	variable.Replace("void ", " ");	variable.Replace(",", "");
					variable.Replace("bool ", " ");	
					variable = Trim(variable);
					a.Replace(variable, "#");
					aa << a;
					aa << variable;
				}
			} else if (line.Find("// End DLL functions") >= 0)
				status = 2;
		}
		if (status == 0)
			before << oline << "\n";
		else if (status == 1)
			declaration << oline << "\n";
		else if (status == 2)
			after << oline << "\n";
		
		if (status == 0 && line.Find(sbegin) >= 0)
			status = 1;
	}		
}

String GenerateFunctionPointers(const UVector<String> &retTypes, const UVector<String> &functions, const UVector<UVector<String>> &arguments) {
	String s;
	for (int i = 0; i < functions.size(); i++) 
		s << "\tstatic " << retTypes[i] << " (*" << functions[i] << ")(" << GetArgs(arguments[i]) << ") = 0;\n";
	return s;
}

String BMR_strCppDeclaration(const String &include) {
	UVector<String> retTypes, functions;
	UVector<UVector<String>> arguments;
	String before, declaration, after;
	
	GetFunctionsList(include, false, retTypes, functions, arguments, before, declaration, after);	
	FuncNode root = BuildTree(retTypes, functions, arguments);
	
	String ret;
	ret << before
	    << "#ifdef BEMROSETTA_DYNAMIC\n"
	    << "    static const char * (*_BMR_GetLastError)() = 0;\n"
	    //<< GenerateFunctionPointers(retTypes, functions, arguments) << "\n"
	    << "#else\n"
	    << declaration << "\n"
	    << "#endif\n"
	    << after
	    << "\n\n// C++ functions\n\n" 
	    << GenerateHeaderCpp(root);

	return ret;
}

void CollectLeaves(const String& fieldPath, const FuncNode& node, UVector<String>& fieldPaths, UVector<String>& fullNames) {
	for (int i = 0; i < node.children.size(); i++) {
		const String& cname = node.children.GetKey(i);
		const FuncNode& child = node.children[i];
		String fp = fieldPath.IsEmpty() ? cname : fieldPath + "." + cname;
		if (child.isFunc) {
			fieldPaths.Add(fp);
			fullNames.Add(child.fullName);
		} else
			CollectLeaves(fp, child, fieldPaths, fullNames);
	}
}

String EmitLeaf(const String& path, const FuncNode& node, const UVector<String>& skipErrorCheck, const String& symMacro, const String& errorCheckFn) {
	bool isVoid = TrimBoth(node.retType) == "void";
	bool skip = FindIndex(skipErrorCheck, path) >= 0; // e.g. BMR_GetLastError itself
	String sym = symMacro + "(" + path + ")";         // -> p_<path> (dynamic) or _<path> (static)

	String s;
	s << "typedef " << node.retType << " (*" << path << "_t)(" << GetArgs(node.args) << ");\n";

	s << "#ifndef BEMROSETTA_DYNAMIC\n"
	  << "extern " << node.retType << " " << sym << "(" << GetArgs(node.args) << ");\n"
	  << "#endif\n";

	String errCall = symMacro + "(" + errorCheckFn + ")()";
	s << "static inline " << node.retType << " " << path << "(" << GetArgs(node.args) << ") {\n";
	if (isVoid) {
		s << "    " << sym << "(" << GetArgNames(node.args) << ");\n";
		if (!skip)
			s << "    if (" << errCall << ")\n        BMR_TriggerError();\n";
	} else {
		s << "    " << node.retType << " ret = " << sym << "(" << GetArgNames(node.args) << ");\n";
		if (!skip)
			s << "    if (" << errCall << ")\n        BMR_TriggerError();\n";
		s << "    return ret;\n";
	}
	s << "}\n\n";
	return s;
}

const FuncNode* FindNode(const FuncNode& node, const UVector<String>& segments, int idx) {
	if (idx == segments.size())
		return &node;
	int f = node.children.Find(segments[idx]);
	if(f < 0)
		return NULL;
	return FindNode(node.children[f], segments, idx + 1);
}

String EmitTypedefs(const String& path, const FuncNode& node, const UVector<String>& skipErrorCheck, const String& symMacro, const String& errorCheckFn) {
	String s;
	for (int i = 0; i < node.children.size(); i++) {
		const String& cname = node.children.GetKey(i);
		String cpath = path.IsEmpty() ? cname : path + "_" + cname;
		s << EmitTypedefs(cpath, node.children[i], skipErrorCheck, symMacro, errorCheckFn);
	}
	if (node.isFunc)
		s << EmitLeaf(path, node, skipErrorCheck, symMacro, errorCheckFn);
	else {
		s << "typedef struct {\n";
		for (int i = 0; i < node.children.size(); i++) {
			const String& cname = node.children.GetKey(i);
			String cpath = path.IsEmpty() ? cname : path + "_" + cname;
			s << "    " << cpath << "_t " << cname << ";\n";
		}
		if (path == "BMR")
			s << "} " << "BEMRosetta" << ";\n\n";
		else
			s << "} " << path << "_t;\n\n";
	}
	return s;
}

String EmitInit(const FuncNode& node) {
	if (node.isFunc)
		return node.fullName;
	String s = "{ ";
	for (int i = 0; i < node.children.size(); i++) {
		if (i) 
			s << ", ";
		s << "." << node.children.GetKey(i) << " = " << EmitInit(node.children[i]);
	}
	return s << " }";
}

String GenerateHeaderC(FuncNode& root, const UVector<String>& skipErrorCheck = UVector<String>{ "BMR_GetLastError" }) {
	ASSERT(root.children.size() == 1);   // all functions share one prefix, e.g. BMR_
	String name = root.children.GetKey(0);
	const FuncNode& top = root.children[0];

	UVector<String> fieldPaths, fullNames;
	CollectLeaves("", top, fieldPaths, fullNames);

	String symMacro = name + "_SYM";
	String errorCheckFn = name + "_GetLastError";

	String s;
	s << "/* Auto-generated by bmr_gen — do not edit by hand. */\n\n";
	s << "#ifdef __cplusplus\nextern \"C\" {\n#endif\n\n";

	s << "#ifdef _WIN32\n"
	  << "    #include <windows.h>\n"
	  << "    #define DLL_HANDLE HMODULE\n"
	  << "    #define DLL_LOAD(file) LoadLibraryA(file)\n"
	  << "    #define DLL_SYM(handle, name) GetProcAddress(handle, name)\n"
	  << "    #define DLL_FREE(handle) FreeLibrary(handle)\n"
	  << "#else\n"
	  << "    #include <dlfcn.h>\n"
	  << "    #define DLL_HANDLE void*\n"
	  << "    #define DLL_LOAD(file) dlopen(file, RTLD_NOW)\n"
	  << "    #define DLL_SYM(handle, name) dlsym(handle, name)\n"
	  << "    #define DLL_FREE(handle) dlclose(handle)\n"
	  << "#endif\n\n";

	s  << "#define " << symMacro << "(fn) _##fn\n";

	s << "void " << name << "_TriggerError(void);\n\n";

	s << EmitTypedefs(name, top, skipErrorCheck, symMacro, errorCheckFn);

	s << "#ifdef BEMROSETTA_DYNAMIC\n"
	  << "#define BEMRosetta_Init(file) \\\n"
      << "     ( BMR_Load(file), \\\n"
	  << "     (BEMRosetta)" << EmitInit(top) << ");BMR_Init()\n"
	  << "#else\n"
	  << "#define BEMRosetta_Init() \\\n"
	  << "    (BEMRosetta)" << EmitInit(top) << ";BMR_Init()\n"
	  << "#endif\n";

	s << "#ifdef BEMROSETTA_DYNAMIC\n"
	  << "static DLL_HANDLE " << name << "_dll = 0;\n\n"
	  << "static inline int " << name << "_Load(const char* file_dll)\n{\n"
	  << "    if(" << name << "_dll)\n        return 0;   /* already loaded */\n"
	  << "    " << name << "_dll = DLL_LOAD(file_dll);\n"
	  << "    if(!" << name << "_dll)\n        return 0;\n";
	for (int i = 0; i < fullNames.size(); i++) {
		String cleanPath = fieldPaths[i];
		cleanPath.Replace(".", "_");
		cleanPath = name + "_" + cleanPath;   // fieldPaths start below `top`, EmitTypedefs' path always includes `name`
		s << "    " << symMacro << "(" << cleanPath << ") = (" << cleanPath
		  << "_t)DLL_SYM(" << name << "_dll, \"" << fullNames[i] << "\");\n";  // fullNames[i] already has its "_"
	}
	s << "    return 1;\n}\n\n"
	  << "static inline void BEMRosetta_Free() {if(BMR_dll) {DLL_FREE(BMR_dll); BMR_dll = 0;}}\n"
	  << "#else\n"
	  << "static inline void BEMRosetta_Free() {;}\n"
	  << "#endif\n\n";

	s << "#ifdef __cplusplus\n}\n#endif\n";
	return s;
}

String BMR_strCDeclaration(const String &include) {
	UVector<String> retTypes, functions;
	UVector<UVector<String>> arguments;
	String before, declaration, after;
	
	GetFunctionsList(include, true, retTypes, functions, arguments, before, declaration, after);	
	FuncNode root = BuildTree(retTypes, functions, arguments);
	
	String ret;
	ret << before
	    << "#ifdef BEMROSETTA_DYNAMIC\n"
	    << GenerateFunctionPointers(retTypes, functions, arguments) << "\n"
	    << "#else\n"
	    << declaration << "\n"
	    << "#endif\n"
	    << after
	    << "\n\n// C++ functions\n\n" 
	    << "static inline void BMR_TriggerError() {\n"
    	<< "	error_callback_t cb = _BMR_GetErrorHandlerCallback();\n"
    	<< "	if (cb)\n"
        << "		cb(_BMR_GetLastError(), _BMR_GetErrorHandlerData());\n"
        << "}\n\n"
	    << GenerateHeaderC(root);

	return ret;
}