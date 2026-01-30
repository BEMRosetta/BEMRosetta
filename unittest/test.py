import os
import sys
sys.path.append(".test")
from libbemrosetta import BEMRosetta

print("libbemrosetta.py test")

try:
    bemr = BEMRosetta("./.test/libbemrosetta.dll")

    print(bemr.Version())
    
    bemr.Mesh.Input("../examples/hydrostar/Mesh/Ship.hst")
    bemr.Mesh.Convert("./.test/kk.gdf", ".gdf", 0, 0)
    bemr.Mesh.Input("./.test/kk.gdf")
    volx, voly, volz = bemr.Mesh.GetVolume()
    print(f"Volume            x: {volx}, y: {voly}, z: {volz}")
    volx, voly, volz  = bemr.Mesh.GetUnderwaterVolume()
    print(f"Underwater volume x: {volx}, y: {voly}, z: {volz}")
    print(f"Surface            : {bemr.Mesh.GetSurface()}")
    print(f"Underwater surface : {bemr.Mesh.GetUnderwaterSurface()}")
    print(f"Stiffness matrix   : {bemr.Mesh.GetHydrostaticStiffness()}")
    
    
    
    os.remove("./.test/kk.gdf")
    
except Exception as e:
    print(f"\nAn error occurred: {e}")
    sys.exit(1)