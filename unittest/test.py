import os
import sys
sys.path.append(".test")
from libbemrosetta import BEMRosetta

print("libbemrosetta.py test")

try:
    bemr = BEMRosetta("./.test/libbemrosetta.dll")

    print(bemr.Version())
    
    bemr.Mesh.Load("../examples/hydrostar/Mesh/Ship.hst")
    bemr.Mesh.Save("./.test/kk.gdf", ".gdf", 0, 0)
    bemr.Mesh.Load("./.test/kk.gdf")
    volx, voly, volz = bemr.Mesh.Volume.Get()
    print(f"Volume            x: {volx}, y: {voly}, z: {volz}")
    volx, voly, volz  = bemr.Mesh.UnderwaterVolume.Get()
    print(f"Underwater volume x: {volx}, y: {voly}, z: {volz}")
    print(f"Surface            : {bemr.Mesh.Surface.Get()}")
    print(f"Underwater surface : {bemr.Mesh.UnderwaterSurface.Get()}")
    print(f"Stiffness matrix   : {bemr.Mesh.HydrostaticStiffness.Get()}")
    
    
    
    os.remove("./.test/kk.gdf")
    
except Exception as e:
    print(f"\nAn error occurred: {e}")
    sys.exit(1)