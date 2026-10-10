import os
import sys
sys.path.append(".test")
from libbemrosetta import BEMRosetta

print("libbemrosetta.py test")

try:
    bmr = BEMRosetta("./.test/libbemrosetta.dll")

    print(bmr.Version())
    
    bmr.Mesh.Load("../examples/hydrostar/Mesh/Ship.hst")
    bmr.Mesh.Save("./.test/kk.gdf", ".gdf", 0, 0)
    bmr.Mesh.Load("./.test/kk.gdf")
    _, volx, voly, volz = bmr.Mesh.Volume.Get()
    print(f"Volume            x: {volx}, y: {voly}, z: {volz}")
    _, volx, voly, volz  = bmr.Mesh.UnderwaterVolume.Get()
    print(f"Underwater volume x: {volx}, y: {voly}, z: {volz}")
    print(f"Surface            : {bmr.Mesh.Surface.Get()}")
    print(f"Underwater surface : {bmr.Mesh.UnderwaterSurface.Get()}")
    bmr.Mesh.Cg.Set(0, 0, -1)
    bmr.Mesh.C0.Set(0, 0, -1)
    print(f"Stiffness matrix   : {bmr.Mesh.HydrostaticStiffness.Get()}")
    
    os.remove("./.test/kk.gdf")
    
except Exception as e:
    print(f"\nAn error occurred: {e}")
    sys.exit(1)