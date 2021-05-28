# LIBRARY DECLARATION
libc = ctypes.CDLL(dll)

# INPUT TYPES
libc.DLL_FAST_Load.argtypes = [ctypes.c_char_p]
libc.DLL_FAST_GetParameterName.argtypes = [ctypes.c_int]
libc.DLL_FAST_GetUnitName.argtypes = [ctypes.c_int]
libc.DLL_FAST_GetParameterId.argtypes = [ctypes.c_char_p]
libc.DLL_FAST_GetTime.argtypes = [ctypes.c_int]
libc.DLL_FAST_GetData.argtypes = [ctypes.c_int, ctypes.c_int]
libc.DLL_FAST_GetAvg.argtypes = [ctypes.c_char_p]
libc.DLL_FAST_LoadFile.argtypes = [ctypes.c_char_p]
libc.DLL_FAST_SaveFile.argtypes = [ctypes.c_char_p]
libc.DLL_FAST_SetVar.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]
libc.DLL_FAST_GetVar.argtypes = [ctypes.c_char_p, ctypes.c_char_p]

# OUTPUT TYPES
libc.DLL_Version.restype = ctypes.c_char_p
libc.DLL_strListFunctions.restype = ctypes.c_char_p
libc.DLL_strPythonDeclaration.restype = ctypes.c_char_p
libc.DLL_FAST_Load.restype = ctypes.c_int
libc.DLL_FAST_GetParameterName.restype = ctypes.c_char_p
libc.DLL_FAST_GetUnitName.restype = ctypes.c_char_p
libc.DLL_FAST_GetParameterId.restype = ctypes.c_int
libc.DLL_FAST_GetParameterCount.restype = ctypes.c_int
libc.DLL_FAST_GetLen.restype = ctypes.c_int
libc.DLL_FAST_GetTimeInit.restype = ctypes.c_double
libc.DLL_FAST_GetTimeEnd.restype = ctypes.c_double
libc.DLL_FAST_GetTime.restype = ctypes.c_double
libc.DLL_FAST_GetData.restype = ctypes.c_double
libc.DLL_FAST_GetAvg.restype = ctypes.c_double
libc.DLL_FAST_LoadFile.restype = ctypes.c_int
libc.DLL_FAST_SaveFile.restype = ctypes.c_int
libc.DLL_FAST_SetVar.restype = ctypes.c_int
libc.DLL_FAST_GetVar.restype = ctypes.c_char_p