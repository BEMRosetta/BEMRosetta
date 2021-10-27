"C:\Program Files (x86)\WiX Toolset v3.11\bin\candle" BEMRosetta.wxs  -arch x64
"C:\Program Files (x86)\WiX Toolset v3.11\bin\light" BEMRosetta.wixobj
move BEMRosetta.msi ../bin
@pause