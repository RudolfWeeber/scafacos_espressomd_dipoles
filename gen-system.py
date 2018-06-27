"""This creates the initial system and thermlizes it without
   dipolar interactions. The resulting configuration is written to the
   file given on the command line.
"""

import model

import pickle
import sys
import argparse
from time import time

import numpy as np

# Command linke handling
parser = argparse.ArgumentParser(
    description="Tenerate gel model and thermalize it.")
parser.add_argument("--output_file", dest="output_file", required=True,
                    help="Name of the file to write the thermalized configuration to.")
parser.add_argument(
    "--steps", dest="steps", default=100, type=int, required=False,
    help="Number of thermlaization steps (100k time steps, each)")
args = parser.parse_args()


# Instance and setup the gel model.
g = model.GelModel()
g.setup_system()
g.setup_particles()
# Adjust the simulation box and move the gel to its center
g.adjust_box_and_shift_system()

# Thermalization loop
warmup_loops = args.steps
warmup_steps = 10000
for i in range(warmup_loops):
    # Tune the skin parameter of the cell system
    print(g.s.cell_system.tune_skin(
        min_skin=0.4, max_skin=4, tol=0.2, int_steps=50))

    # Time the integration
    start = time()
    g.s.integrator.run(100)
    print(i, (time() - start) / 100, "seconds per time step")

    # Run the integration
    g.s.integrator.run(warmup_steps)
    # Adjust the simulation box and move the gel to its center
    g.adjust_box_and_shift_system()

# Save the resulting configuration
f = open(args.output_file, "w")
pickle.dump({
            "box_l": g.s.box_l,
            "id": g.s.part[:].id,
            "pos": g.s.part[:].pos,
            "dip": g.s.part[:].dip}, f)
f.close()

print("Volume fraction:", len(g.s.part)
      / 6 * np.pi * g.lj_sigma**3 / g.s.volume())
