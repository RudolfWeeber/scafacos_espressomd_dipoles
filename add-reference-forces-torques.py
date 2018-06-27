"""This loads a particle configuration from a file specified on the command
   line, computes dipolar forces and torques via direct summation and stores
   them in the loaded file."

"""

import pickle
import argparse
import sys
import numpy as np
from time import time
from p2nfft_common import set_box_and_particles

import espressomd
from espressomd.magnetostatics import DipolarDirectSumCpu, Scafacos


# Command linke handling
# Command linke handling
parser = argparse.ArgumentParser(
    description="Add reference forces and torques obtained from direct summation to a particle configuration")
parser.add_argument("particle_config_file", nargs=1,
                    help="Name of the file containing the particle configuration")
args = parser.parse_args()

# Setup Espresso system
s = espressomd.System(box_l=[1,1,1]) # Resizing occurs later
s.periodicity = 0, 0, 0
s.time_step = 0.01
s.cell_system.skin = 0

# Load particle configuration
config = pickle.load(open(args.particle_config_file[0]))

# Set particle positions and dipole moments
set_box_and_particles(s, config)

# Setup dipolar direct summation
# dds=DipolarDirectSumCpu(prefactor=1)

# Alternatively use scafacos for direct summation
dds = Scafacos(method_name="direct",
               prefactor=1, method_params={"direct_cutoff": 0})

s.actors.add(dds)

# Switch of thermal noise
s.thermostat.turn_off()

# Calculate
start = time()
s.integrator.run(0)
print(time() - start, "seconds for direct summation")


# Save reference data
config["reference_force"] = s.part[:].f
config["reference_torque"] = s.part[:].torque_lab

f = open(args.particle_config_file[0], "w")
pickle.dump(config, f)
f.close()
