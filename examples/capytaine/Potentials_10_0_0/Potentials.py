import numpy as np
import capytaine as cpt
from capytaine.io.xarray import problems_from_dataset
from capytaine.bem.airy_waves import airy_waves_pressure
import xarray as xr
import os

mesh_1 = cpt.load_mesh("./mesh/Body_1.dat", file_format="nemoh")

body_1 = cpt.FloatingBody(mesh=mesh_1, dofs=cpt.rigid_body_dofs(rotation_center=(10.000000, 0.000000, 0.000000)), center_of_mass=(0.000000, 0.000000, 0.000000), name="Body 1")

body_1.inertia_matrix = body_1.compute_rigid_body_inertia()
body_1.hydrostatic_stiffness = body_1.compute_hydrostatic_stiffness()

list_of_bodies = [body_1]
all_bodies = cpt.FloatingBody.join_bodies(*list_of_bodies)
test_matrix = xr.Dataset(coords={
    "omega": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5],
    "wave_direction": [0, 1.5707963267949, 4.71238898038469],
    "radiating_dof": list(all_bodies.dofs),
    "water_depth": 100,
    "rho": 1025
})

solver = cpt.BEMSolver()
pbs = problems_from_dataset(test_matrix, all_bodies)
results = solver.solve_all(pbs, keep_details=True)
ds = cpt.assemble_dataset(results)

mesh = all_bodies.mesh

ds.coords["space_coordinate"] = ["x", "y", "z"]
ds["mesh_vertices"] = (["face", "vertices_of_face", "space_coordinate"], mesh.vertices[mesh.faces])
ds["mesh_faces_center"] = (["face", "space_coordinate"], mesh.faces_centers)
ds["dof_definition"] = (["radiating_dof", "face", "space_coordinate"], np.array([all_bodies.dofs[dof] for dof in all_bodies.dofs]))

ds["incident_pressure"] = (
    ["omega", "wave_direction", "face"],
    np.zeros((ds.sizes["omega"], ds.sizes["wave_direction"], all_bodies.mesh.nb_faces,), dtype=np.complex128),
)
ds["diffraction_pressure"] = (
    ["omega", "wave_direction", "face"],
    np.zeros((ds.sizes["omega"], ds.sizes["wave_direction"], all_bodies.mesh.nb_faces), dtype=np.complex128),
)
ds["radiation_pressure"] = (
    ["omega", "radiating_dof", "face"],
    np.zeros((ds.sizes["omega"], ds.sizes["radiating_dof"], all_bodies.mesh.nb_faces), dtype=np.complex128),
)

for res in results:
    if isinstance(res.problem, cpt.DiffractionProblem):
        ds["diffraction_pressure"].loc[dict(omega=res.omega, wave_direction=res.wave_direction)] = res.pressure
        ds["incident_pressure"].loc[dict(omega=res.omega, wave_direction=res.wave_direction)] = airy_waves_pressure(mesh, res)
    elif isinstance(res.problem, cpt.RadiationProblem):
        ds["radiation_pressure"].loc[dict(omega=res.omega, radiating_dof=res.radiating_dof)] = res.pressure

ds.coords["rigid_body_component"] = [body.name for body in list_of_bodies]
ds["rotation_center"] = (["rigid_body_component", "point_coordinates"], [body.rotation_center for body in list_of_bodies])
ds["center_of_mass"] = (["rigid_body_component", "point_coordinates"], [body.center_of_mass for body in list_of_bodies])

# Export to NetCDF file
from capytaine.io.xarray import separate_complex_values
separate_complex_values(ds).to_netcdf("Potentials.nc",
                                          encoding={"radiating_dof": {"dtype": "U"},
                                                    "influenced_dof": {"dtype": "U"}})
