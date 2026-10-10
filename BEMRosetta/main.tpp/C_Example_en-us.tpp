topic "C Example";
[ $$0,0#00000000000000000000000000000000:Default]
[{_} 
[s0; [*4 C Example]&]
[s0;*4 &]
[s0; [C@(160.80.0)2 #include <stdio.h>]&]
[s0; [C@(160.80.0)2 #include <stdlib.h> ]&]
[s0; [C@(160.80.0)2 #include `".test`\`\libbemrosetta.h`"]&]
[s0;C2 &]
[s0; [C@(0.0.220)2 void][C2  my`_error`_handler(][C@(0.0.220)2 const][C2  
][C@(0.0.220)2 char][C2  `*message, ][C@(0.0.220)2 void][C2  `*dummy) 
`{]&]
[s0; [C2 -|printf(][C@(163.21.21)2 `"`\nError %s`"][C2 , message);]&]
[s0; [C2 -|printf(][C@(163.21.21)2 `"`\nClick Enter to end`"][C2 );]&]
[s0; [C2 -|getchar();-|]&]
[s0; [C2 -|exit(`-][C@(128.0.128)2 1][C2 );]&]
[s0; [C2 `}]&]
[s0;C2 &]
[s0; [C@(0.0.220)2 int][C2  main() `{]&]
[s0; [C2 -|]&]
[s0; [C2 -|-|printf(][C@(163.21.21)2 `"BEMRosetta C demo`\n`"][C2 );]&]
[s0;C2 &]
[s0; [C2 -|][C@(160.80.0)2 #ifdef BEMROSETTA`_DYNAMIC-|-|// DLL is loaded 
in runtime]&]
[s0; [C2 -|-|BEMRosetta bmr `= BEMRosetta`_Init(][C@(163.21.21)2 `"libbemrosetta.dll`"][C2 );
]&]
[s0; [C2 -|][C@(160.80.0)2 #else-|-|-|-|-|-|-|// DLL has to be in the .exe folder 
or in the PATH]&]
[s0; [C2 -|-|BEMRosetta bmr `= BEMRosetta`_Init();]&]
[s0; [C2 -|][C@(160.80.0)2 #endif]&]
[s0; [C2 -|]&]
[s0; [C2 -|-|bmr.SetErrorHandler(my`_error`_handler, ][C@(128.0.128)2 0][C2 );]&]
[s0; [C2 -|-|printf(][C@(163.21.21)2 `"`\nBEMRosetta version is %s`\n`"][C2 , 
bmr.Version());]&]
[s0;C2 &]
[s0; [C2 -|-|printf(][C@(163.21.21)2 `"`\n`- Mesh handling`"][C2 );]&]
[s0; [C2 -|-|][C@(0.0.220)2 const][C2  ][C@(0.0.220)2 char][C2  `*meshFile `= 
][C@(163.21.21)2 `"../examples/capytaine/Orca/Body`_1.gdf`"][C2 ;]&]
[s0; [C2 -|-|][C@(0.0.220)2 int][C2  idMesh `= bmr.Mesh.Load(meshFile);]&]
[s0; [C2 -|-|printf(][C@(163.21.21)2 `"`\nLoaded mesh `'%s`'`"][C2 , meshFile);]&]
[s0;C2 &]
[s0; [C2 -|-|][C@(0.0.220)2 double][C2  volx, voly, volz;]&]
[s0; [C2 -|-|bmr.Mesh.UnderwaterVolume.Get(`&volx, `&voly, `&volz);]&]
[s0; [C2 -|-|printf(][C@(163.21.21)2 `"`\nUnderwater volume %f`"][C2 , (volx 
`+ voly `+ volz)/][C@(128.0.128)2 3][C2 );]&]
[s0; [C2 -|-|][C@(0.0.220)2 double][C2  x, y, z;]&]
[s0; [C2 -|-|bmr.Mesh.Centre.Volume.Get(`&x, `&y, `&z);]&]
[s0; [C2 -|-|printf(][C@(163.21.21)2 `"`\nCentre of volume is %f, %f, %f`"][C2 , 
x, y, z);]&]
[s0;C2 &]
[s0; [C2 -|-|printf(][C@(163.21.21)2 `"`\n`- BEM case generation`"][C2 );]&]
[s0; [C2 -|-|bmr.Mesh.Cg.Set(][C@(128.0.128)2 0][C2 , ][C@(128.0.128)2 0][C2 , 
`-][C@(128.0.128)2 1][C2 );]&]
[s0; [C2 -|-|bmr.Mesh.C0.Set(][C@(128.0.128)2 0][C2 , ][C@(128.0.128)2 0][C2 , 
`-][C@(128.0.128)2 1][C2 );]&]
[s0; [C2 -|-|][C@(0.0.220)2 double][C2  M`[`] `= `{][C@(128.0.128)2 9.8E5][C2 , 
    ][C@(128.0.128)2 0][C2 ,     ][C@(128.0.128)2 0][C2 ,   ][C@(128.0.128)2 0][C2 , 
  ][C@(128.0.128)2 0][C2 ,   ][C@(128.0.128)2 0][C2 ,]&]
[s0; [C2 -|-|-|-|-|-|  ][C@(128.0.128)2 0][C2 , ][C@(128.0.128)2 9.8E5][C2 ,     
][C@(128.0.128)2 0][C2 ,   ][C@(128.0.128)2 0][C2 ,   ][C@(128.0.128)2 0][C2 , 
  ][C@(128.0.128)2 0][C2 ,]&]
[s0; [C2 -|-|-|-|-|-|  ][C@(128.0.128)2 0][C2 ,     ][C@(128.0.128)2 0][C2 , ][C@(128.0.128)2 9.8E5][C2 ,
   ][C@(128.0.128)2 0][C2 ,   ][C@(128.0.128)2 0][C2 ,   ][C@(128.0.128)2 0][C2 ,]&]
[s0; [C2 -|-|-|-|-|-|  ][C@(128.0.128)2 0][C2 ,     ][C@(128.0.128)2 0][C2 ,     
][C@(128.0.128)2 0][C2 , ][C@(128.0.128)2 1E7][C2 ,   ][C@(128.0.128)2 0][C2 , 
  ][C@(128.0.128)2 0][C2 ,]&]
[s0; [C2 -|-|-|-|-|-|  ][C@(128.0.128)2 0][C2 ,     ][C@(128.0.128)2 0][C2 ,     
][C@(128.0.128)2 0][C2 ,   ][C@(128.0.128)2 0][C2 , ][C@(128.0.128)2 1E7][C2 , 
  ][C@(128.0.128)2 0][C2 ,]&]
[s0; [C2 -|-|-|-|-|-|  ][C@(128.0.128)2 0][C2 ,     ][C@(128.0.128)2 0][C2 ,     
][C@(128.0.128)2 0][C2 ,   ][C@(128.0.128)2 0][C2 ,   ][C@(128.0.128)2 0][C2 , 
][C@(128.0.128)2 1E7][C2 `};]&]
[s0; [C2 -|-|][C@(0.0.220)2 int][C2  dim`[`] `= `{][C@(128.0.128)2 6][C2 , ][C@(128.0.128)2 6][C2 `};
]&]
[s0; [C2 -|-|bmr.Mesh.Inertia.Set(M, dim);]&]
[s0;C2 &]
[s0; [C2 -|-|bmr.Bem.New();]&]
[s0; [C2 -|-|bmr.Bem.depth.Set(][C@(128.0.128)2 50][C2 );]&]
[s0; [C2 -|-|bmr.Bem.g.Set(][C@(128.0.128)2 9.81][C2 );]&]
[s0; [C2 -|-|bmr.Bem.rho.Set(][C@(128.0.128)2 1025][C2 );]&]
[s0; [C2 -|-|][C@(0.0.220)2 double][C2  w`[`] `= `{][C@(128.0.128)2 0.1][C2 , 
][C@(128.0.128)2 0.5][C2 , ][C@(128.0.128)2 1][C2 , ][C@(128.0.128)2 1.5][C2 , 
][C@(128.0.128)2 2][C2 `};]&]
[s0; [C2 -|-|bmr.Bem.w.Set(w, ][C@(0.0.220)2 sizeof][C2 (w)/][C@(0.0.220)2 sizeof][C2 (][C@(0.0.220)2 d
ouble][C2 ));]&]
[s0; [C2 -|-|][C@(0.0.220)2 double][C2  head`[`] `= `{][C@(128.0.128)2 0][C2 , 
][C@(128.0.128)2 45][C2 , ][C@(128.0.128)2 90][C2 `};]&]
[s0; [C2 -|-|bmr.Bem.headings.Set(head, ][C@(0.0.220)2 sizeof][C2 (head)/][C@(0.0.220)2 sizeof][C2 (
][C@(0.0.220)2 double][C2 ));]&]
[s0;C2 &]
[s0; [C2 -|-|bmr.Bem.Mesh.Load(][C@(128.0.128)2 0][C2 , idMesh);]&]
[s0; [C2 -|-|]&]
[s0; [C2 -|-|printf(][C@(163.21.21)2 `"`\nSaving it in Capytaine format`"][C2 );]&]
[s0; [C2 -|-|bmr.Bem.SaveCase(][C@(163.21.21)2 `"../unittest/.test/Capy`"][C2 , 
][C@(163.21.21)2 `"Capytaine .py`"][C2 , ][C@(0.0.220)2 false][C2 , ][C@(0.0.220)2 false][C2 ,
 ][C@(0.0.220)2 true][C2 , ][C@(0.0.220)2 true][C2 , ][C@(163.21.21)2 `"No`"][C2 , 
][C@(0.0.220)2 false][C2 , ][C@(0.0.220)2 false][C2 , ][C@(128.0.128)2 1][C2 , 
][C@(128.0.128)2 4][C2 , ][C@(0.0.220)2 false][C2 , ][C@(0.0.220)2 false][C2 );]&]
[s0;C2 &]
[s0; [C2 -|-|printf(][C@(163.21.21)2 `"`\nRunning it in Capytaine`"][C2 );]&]
[s0; [C2 -|-|system(][C@(163.21.21)2 `"cd /d ..`\`\unittest`\`\.test`\`\Capy 
`&`& capytaine.bat`"][C2 );]&]
[s0; [C2 -|-|]&]
[s0; [C2 -|-|printf(][C@(163.21.21)2 `"`\nLoading the results in .nc format`"][C2 );]&]
[s0; [C2 -|-|bmr.Bem.Load(][C@(163.21.21)2 `"..`\`\unittest`\`\.test`\`\Capy`\`\capytaine.n
c`"][C2 );]&]
[s0; [C2 -|-|printf(][C@(163.21.21)2 `"`\nSaving the results in .h5 format`"][C2 );]&]
[s0; [C2 -|-|bmr.Bem.Save(][C@(163.21.21)2 `"..`\`\unittest`\`\.test`\`\Capy`\`\capytaine.h
5`"][C2 );]&]
[s0; [C2 -|]&]
[s0; [C2 -|printf(][C@(163.21.21)2 `"`\nProgram ended`\n`"][C2 );]&]
[s0;C2 &]
[s0; [C@(160.80.0)2 #ifdef BEMROSETTA`_DYNAMIC]&]
[s0; [C2     BEMRosetta`_Free();]&]
[s0; [C@(160.80.0)2 #endif]&]
[s0; [C2 -|-|]&]
[s0; [C2 -|][C@(0.0.220)2 return][C2  ][C@(128.0.128)2 0][C2 ;]&]
[s0; [C2 `}]]]