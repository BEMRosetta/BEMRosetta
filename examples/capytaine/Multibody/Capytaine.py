import numpy as np
import xarray as xr
import capytaine as cpt

body_1 = cpt.FloatingBody(
            mesh=cpt.mesh_sphere(center=(0, 0, 0)),
            dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
            center_of_mass=(0, 0, 0)
        )
body_1.inertia_matrix = body_1.compute_rigid_body_inertia()
body_1.hydrostatic_stiffness = body_1.immersed_part().compute_hydrostatic_stiffness()
# If you have several rigid bodies, copy the code above to define "body_2", "body_3", etc.

body_2 = cpt.FloatingBody(
            mesh=cpt.mesh_sphere(center=(0, 10, 0)),
            dofs=cpt.rigid_body_dofs(rotation_center=(0, 10, 0)),
            center_of_mass=(0, 10, 0)
        )
body_2.inertia_matrix = body_2.compute_rigid_body_inertia()
body_2.hydrostatic_stiffness = body_2.immersed_part().compute_hydrostatic_stiffness()

list_of_bodies = [body_1, body_2]  # Replace "[body_1]" by "[body_1, body_2, body_3]" for multibody problem.

all_bodies = cpt.FloatingBody.join_bodies(*list_of_bodies).immersed_part()

# Set up paramaters
test_matrix = xr.Dataset({
        "omega": np.linspace(0.1, 6.0, 50),  # Can also specify "period", "wavelength" or "wavenumber"
        "wave_direction": np.linspace(0, np.pi, 3),
        "radiating_dof": list(all_bodies.dofs),
        "water_depth": [10.0],  # np.inf],
        "rho": [1025],          # 1025],
        })

# Do the resolution
solver = cpt.BEMSolver()
dataset = solver.fill_dataset(test_matrix, all_bodies)

# from capytaine.post_pro.rao import rao
# dataset["rao"] = rao(dataset)

dataset.coords["rigid_body_component"] = [body.name for body in list_of_bodies]
dataset["rotation_center"] = (["rigid_body_component", "point_coordinates"], [body.rotation_center for body in list_of_bodies])

# Export to NetCDF file
from capytaine.io.xarray import separate_complex_values
separate_complex_values(dataset).to_netcdf("dataset.nc",
                                           encoding={'radiating_dof': {'dtype': 'U'},
                                                     'influenced_dof': {'dtype': 'U'}})
                                                     
# Export to Nemoh file
from capytaine.io.legacy import write_dataset_as_tecplot_files
write_dataset_as_tecplot_files("./results", dataset)
