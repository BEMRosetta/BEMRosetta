# BEMRosetta python functions list

class BEMRosetta:
    def __init__(self, path_dll):
        self.libc = ctypes.CDLL(dll)

        # INPUT TYPES
        self.libc.DemoVectorPy_C.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64), c_int]
        self.libc.DemoVectorC_Py.argtypes = [ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.c_int)]
        self.libc.DLL_FAST_Load.argtypes = [c_char_p]
        self.libc.DLL_FAST_GetParameterName.argtypes = [c_int]
        self.libc.DLL_FAST_GetUnitName.argtypes = [c_int]
        self.libc.DLL_FAST_GetParameterId.argtypes = [c_char_p]
        self.libc.DLL_FAST_GetTime.argtypes = [c_int]
        self.libc.DLL_FAST_GetData.argtypes = [c_int, c_int]
        self.libc.DLL_FAST_GetAvg.argtypes = [c_char_p]
        self.libc.DLL_FAST_LoadFile.argtypes = [c_char_p]
        self.libc.DLL_FAST_SaveFile.argtypes = [c_char_p]
        self.libc.DLL_FAST_SetVar.argtypes = [c_char_p, c_char_p, c_char_p]
        self.libc.DLL_FAST_GetVar.argtypes = [c_char_p, c_char_p]

        # OUTPUT TYPES
        self.libc.DemoVectorPy_C.restype = ctypes.c_double
        self.libc.DemoVectorC_Py.restype = ctypes.c_int
        self.libc.DLL_Version.restype = ctypes.c_char_p
        self.libc.DLL_strListFunctions.restype = ctypes.c_char_p
        self.libc.DLL_strPythonDeclaration.restype = ctypes.c_char_p
        self.libc.DLL_FAST_Load.restype = ctypes.c_int
        self.libc.DLL_FAST_GetParameterName.restype = ctypes.c_char_p
        self.libc.DLL_FAST_GetUnitName.restype = ctypes.c_char_p
        self.libc.DLL_FAST_GetParameterId.restype = ctypes.c_int
        self.libc.DLL_FAST_GetParameterCount.restype = ctypes.c_int
        self.libc.DLL_FAST_GetLen.restype = ctypes.c_int
        self.libc.DLL_FAST_GetTimeInit.restype = ctypes.c_double
        self.libc.DLL_FAST_GetTimeEnd.restype = ctypes.c_double
        self.libc.DLL_FAST_GetTime.restype = ctypes.c_double
        self.libc.DLL_FAST_GetData.restype = ctypes.c_double
        self.libc.DLL_FAST_GetAvg.restype = ctypes.c_double
        self.libc.DLL_FAST_LoadFile.restype = ctypes.c_int
        self.libc.DLL_FAST_SaveFile.restype = ctypes.c_int
        self.libc.DLL_FAST_SetVar.restype = ctypes.c_int
        self.libc.DLL_FAST_GetVar.restype = ctypes.c_char_p

    def DemoVectorPy_C(self, v):
        self.libc.DemoVectorPy_C(v, len(v))

    def DemoVectorC_Py(self, v):
        # Argument preparation
        _data0 = ctypes.POINTER(ctypes.c_double)()
        _size0 = ctypes.c_int()
        # DLL function call
        ret = self.libc.DemoVectorC_Py(ctypes.byref(_data0), ctypes.byref(_size0))
        # Vector processing
        _arraySize0 = ctypes.c_double * _size0.value
        _data0_pointer = ctypes.cast(_data0, ctypes.POINTER(_arraySize0))
		v = np.frombuffer(_data0_pointer.contents)
        return ret

    def FAST_Load(self, filename):
        self.libc.DLL_FAST_Load(str.encode(filename, 'UTF-8'))

    def FAST_GetParameterName(self, id):
        self.libc.DLL_FAST_GetParameterName(id)

    def FAST_GetUnitName(self, id):
        self.libc.DLL_FAST_GetUnitName(id)

    def FAST_GetParameterId(self, name):
        self.libc.DLL_FAST_GetParameterId(str.encode(name, 'UTF-8'))

    def FAST_GetTime(self, idtime):
        self.libc.DLL_FAST_GetTime(idtime)

    def FAST_GetData(self, idtime, idparam):
        self.libc.DLL_FAST_GetData(idtime, idparam)

    def FAST_GetAvg(self, param):
        self.libc.DLL_FAST_GetAvg(str.encode(param, 'UTF-8'))

    def FAST_LoadFile(self, file):
        self.libc.DLL_FAST_LoadFile(str.encode(file, 'UTF-8'))

    def FAST_SaveFile(self, file):
        self.libc.DLL_FAST_SaveFile(str.encode(file, 'UTF-8'))

    def FAST_SetVar(self, name, paragraph, value):
        self.libc.DLL_FAST_SetVar(str.encode(name, 'UTF-8'), str.encode(paragraph, 'UTF-8'), str.encode(value, 'UTF-8'))

    def FAST_GetVar(self, name, paragraph):
        self.libc.DLL_FAST_GetVar(str.encode(name, 'UTF-8'), str.encode(paragraph, 'UTF-8'))