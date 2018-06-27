"""Common routines for using p2nfft with Espresso"""

from espressomd.magnetostatics import Scafacos
import numpy as np


def get_p2nfft_for_config(system, config):
    """Creates the instance of the Scafacos interaction handle and sets
       it up for p2nfft with the parameters taken from the dictionary
       containing the tuning output.
    """

    # The 3d mesh size is passed as string of comma-separated numbres
    mesh = config["mesh"]
    # The scalar mesh is applied to the longest box direction.
    # The mesh in the other direction is calculated to give a comparable
    # density of mesh points
    mesh_str = get_mesh_str_for_box(system, config)

    # Dictionary containing the method parameters as expected by the
    # p2nfft method in the Scafacos library
    method_params = {
        "p2nfft_verbose_tuning": 0,
      "pnfft_N": mesh_str,
      "pnfft_n": mesh_str,
      "pnfft_direct": 0,
      "p2nfft_r_cut": config["r_cut"],
      "p2nfft_alpha": config["alpha"],
      "p2nfft_intpol_order": "3",
      "p2nfft_reg_kernel_name": "ewald",
      "p2nfft_p": 12,
      "p2nfft_ignore_tolerance": "1",
      "pnfft_window_name": "bspline",
      "pnfft_m": config["assignment_order"],
      "pnfft_diff_ik": "0",
      "p2nfft_epsB": config["eps"]}
    scafacos = Scafacos(
        prefactor=1, method_name="p2nfft", method_params=method_params)

    return scafacos


def get_mesh_str_for_box(system, config):
    """Calculate the number of mesh points in each Cartesian direction.
       The scalar mesh is applied to the longest box direction. The toher
       mesh sizes are calculated to give a comparable mesh density.
       If the p2nfft configuration contains a box length for which the tuning
       was done, then this is used to calculate the mesh density. Otherwise,
       the current simulation box from Espresso is used.
       Returns a string of the form "a,b,c", where a,b,c are the mesh sizes
       for the Cartesian directions.

       Parameters:
         system: Espresso system handle
         config: p2nfft method configuration

    """

    # Do we have a box_l in the p2nfft configuration?
    if "box_l" in config:
        max_box_l = np.amax(config["box_l"])
    else:
        max_box_l = np.amax(system.box_l)

    # Mesh for the longest direction
    max_mesh = int(config["mesh"])

    # Mesh density in the longest direction
    mesh_dens = float(max_mesh / max_box_l)

    # Target mesh
    mesh = [0., 0., 0.]

    # Mesh for different Cartesian directions, based on the mesh density.
    for i in range(3):
        mesh[i] = int(system.box_l[i] * mesh_dens)
        # Force an even number of mesh sites per processor
        while mesh[i] % (2 * system.cell_system.get_state()["node_grid"][i]) != 0:
            mesh[i] += 1

    # Return the mesh as string as required by the Scafacos library
    mesh_str = "%d,%d,%d" % (mesh[0], mesh[1], mesh[2])
    return mesh_str


def adjust_p2nfft_mesh(system, config):
    """Expands orh shrinks the p2nfft mesh after a box size change maintaining
       the mesh density.

    """

    # Get the active scafacos interaction
    scafacos = None
    for a in system.actors:
        if isinstance(a, Scafacos):
            scafacos = a
    if scafacos is None:
        raise Exception("Scafacos is not active.")

    # Get the new mesh size
    mesh_str = get_mesh_str_for_box(system, config)
    # Get method parameters from active scafacos handle
    params = scafacos.get_params()["method_params"]

    # Set the new mesh size
    params["pnfft_N"] = mesh_str
    params["pnfft_n"] = mesh_str

    # Re-add the short range cutoff to work around an inconsistency in the
    # scafacos p2nfft interface
    params["p2nfft_r_cut"] = config["r_cut"]

    # Set the new method parameters
    scafacos.set_params({
                        "prefactor": 1, "method_name": "p2nfft", "method_params": params
                        })


def set_box_and_particles(system, particle_config):
    """Sets up the box size and particles from the provided particle configuration
       Note that existing particles are deleted.

       Parameters:
         system: Handle of the Espresso system
         particle_config: particle configuration

    """
    # Box size
    system.box_l = particle_config["box_l"]

    # PArticles
    system.part.clear()
    system.part.add(
        id=particle_config["id"], pos=particle_config["pos"], dip=particle_config["dip"])
    # Enable rotational degrees of freedom
    system.part[:].rotation = 1, 1, 1
