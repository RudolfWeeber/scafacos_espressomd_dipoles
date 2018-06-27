"""This tunes the p2nfft method parameters based on a particle configuration which includes reference fores and torques."""

from __future__ import print_function

import espressomd
import espressomd.magnetostatics as magnetostatics

import sys
import argparse
import pickle
from time import time

import numpy as np
from scipy.optimize import fmin_l_bfgs_b


from p2nfft_common import get_p2nfft_for_config, set_box_and_particles, adjust_p2nfft_mesh


def setup_es():
        """Setup the Espresso sytem and return the handle"""
        s=espressomd.System(box_l=[1,1,1]) # Will be resized later

        s.time_step= 0.01
        s.cell_system.skin =.4
        s.thermostat.turn_off()
        s.box_l =1,1,1
        s.periodicity = 0,0,0
        return s

def load_configs(fnames):
    """Loads the particle configurations to base the tuning on and
       returns them as a list of dictionaries.
       Also ensures that these configurations contain reference forces and 
       torques.
       
       Parameters:
         fnames: list of file names
      """
    configs=[]

    for fname in fnames:
        config=pickle.load(open(fname))
        required="box_l","id","pos","dip","reference_force", "reference_torque"
        for k in required:
            if k not in config:
                print("The configuration in ",fname,"lacks the key",k)
                sys.exit(1)
        configs.append(config)

    return configs



def run(system, p2nfft_config,reference_configs):
    """Calculates errors for the given p2nfft configuration and particle configurations.
       Returns (rms_force_error,rms_torque_error)
       
       Parameters:
         system: the Espresso system handle
         p2nfft_config: The dictionary containing the p2nfft method parameters
         reference_configs: Reference particle configurations as provided by load_configs()
      """
    
    # Force and torque errors
    err_f=0.
    err_t=0.
    
    # Setup a valid particle configuration so p2nfft is tuned for the correct
    # box size
    set_box_and_particles(system,reference_configs[0])
    
    # Get the scafacos handle and add it to the Espresso system
    scafacos =get_p2nfft_for_config(system,p2nfft_config)
    s.actors.add(scafacos)
    
    # Iterte over reference configurations
    first=True
    for config in reference_configs:
        # Set particles according to configuration
        set_box_and_particles(system,config)
        # Number of particles
        n=len(system.part)
        
        # Calculate forces/torques
        s.integrator.run(0)

        
        # Calculate rms force and torque errors
        err_f+=np.sqrt(np.sum(np.sum((system.part[:].f-config["reference_force"])**2,1),0)/n)
        err_t+=np.sqrt(np.sum(np.sum((system.part[:].torque_lab-config["reference_torque"])**2,1),0)/n)


    
    # Normalize by number of configurations considered
    n_configs=len(reference_configs)
    err_f/=n_configs
    err_t/=n_configs

    # De-activate p2nfft
    system.actors.remove(scafacos)

    print( p2nfft_config["r_cut"],p2nfft_config["alpha"],p2nfft_config["eps"],":", err_f,err_t)
    return err_f,err_t

    
def time_calculation(system, p2nfft_config,reference_config):
    """Times the p2nfft calculation for a single reference particle configuration.
    Returns the calculation time in seconds

    Parameters:
      system: the Espresso system handle
      p2nfft_config: dictionary containing the p2nfft method parameters
      reference_config: a single reference particle configuration. One element from the list provided by load_configs()
    """
    
    # Load the reference configuration into the system
    set_box_and_particles(system,reference_config)
    
    # Get scafacos handle and add it to the Espresso sytem
    scafacos =get_p2nfft_for_config(system,p2nfft_config)
    s.actors.add(scafacos)

    # Run once to get initial tuning done
    system.integrator.run(0)

    # Run again and measure the time
    start=time()
    system.integrator.run(0,recalc_forces=True)
    system.actors.remove(scafacos)

    return time()-start

    

def cost_function(x,system,configs, p):
  """Cost function passed to the scipy optimizer. 
     Returns the average of rms_force_error and rms_torque_error provided by run()

     Parameters:
       x: list containing the p2nfft method parameters [alpha,epsilon] tried in this step
       system: Espresso system handle
       configs: reference particle configurations as provided by load_configs()
       p: dictionary containing fixed p2nfft method parameters
    """
     
  # Add current values tried by the optimizer to the p2nfft config
  p["alpha"]=x[0]
  p["eps"]=x[1]

  res= run(system,p,configs)
  return 0.5*(res[0]+res[1])


def minimize_via_alpha_and_eps(system,configs,p):
  """Numerically minimizes the rms force and torque errors by varying alpha and epsilon
     Returns the result object of fmin_l_bfgs_b() method from scipy

     Parameters:
       system: the Espresso system handle
       configs: the reference particle configurations ad provide dby load_configs()
       p: dictionary containing the p2nfft method parameters
  """
  min_alpha=0.5/p["r_cut"]
  max_alpha=10/p["r_cut"]
  res= fmin_l_bfgs_b(cost_function,(p["alpha"],p["eps"]), args=(system,configs,p),
      bounds=((min_alpha,max_alpha),(0.001,0.15)),approx_grad=1,epsilon=0.01,maxfun=50,factr=1E9)
  print( res)
  return res

def tune_r_cut(system,mesh=None,min_r_cut=None,max_r_cut=None,r_cut_tolerance=None, target_accuracy=None,assignment_order=None, particle_configurations=None):
  """Tune the real space cutoff via a bisection.
     Return the p2nfft parameters for the smallest r_cut that is consistent with the target accuracy. The calculation time and box size for which the tuning was done are also stored.


     Parameters:
       system: the Espresso system handle
       mesh: p2nfft mesh for the longest box axis
       min_r_cut: Minimum real space cutoff for the bisection
       max_r_cut: maximum cutoff for the bisection
       r_cut_tolerance: tolerance up to which r_cut is determined
       assignment_order: p2nfft assignment order
       target_accuracy: average of rms force and torque errors to tune to
   """

  # Load particle configurations
  configs=load_configs(particle_configurations)
  
  # Bisect cutoff to find the lowest one meeting the accuracy requirement
  # Current best result 
  best_err=1E20
  best_r_cut=max_r_cut
  best_p=None

  # Parameter dictionary for p2nfft
  p={}
  # Fixed parameters
  p["mesh"] =mesh
  p["assignment_order"]=assignment_order


  
  # Initial guesses for alpha and epsilon
  p["eps"]=0.1
  p["alpha"]=1


  print("Mesh",mesh,"Bisecting r_cut in the range",min_r_cut,max_r_cut)

  
  
  iter =0
  
  # This has to run at least once to get an accuracy even if max_r_cut - min_r_cut <r_cut_tolerance
  while (max_r_cut-min_r_cut)>=r_cut_tolerance or iter==0:
    # Current r_cut
    r_cut=0.5*(max_r_cut+min_r_cut)
    
    # Add cutoff to p2nfft parameters
    p["r_cut"] =r_cut
    
    # Minimize the error for this r_cut by varying alpha end epsilon
    res=minimize_via_alpha_and_eps(system,configs,p)
    
    iter+=1
    
    # Process result
    err=res[1]
    res_p=res[0]
    
    # Is the result better than the previous best result?
    if err< target_accuracy and r_cut<best_r_cut:
      best_err=err
      best_p=p
      best_p["alpha"]=res_p[0]
      best_p["eps"]=res_p[1]
      best_r_cut=r_cut


    # Adjust bounds for the bisection
    if err>target_accuracy:
      # r_cut too small. Adjust minimum 
      min_r_cut=r_cut
    else:
      # r_cut large enough. Adjust maximum
      if err<target_accuracy:
        max_r_cut=r_cut
      else:
        # Exactly at target accuracy. SStop
        break
  
  # Was a viable r_cut found?
  if not best_p is None:
      # Time it
      best_p["time"]=time_calculation(system,best_p,configs[0])
      # Save the box size to be able to calculate the mesh density for which the tuning was done
      best_p["box_l"]=system.box_l
      print("Calculation time",best_p["time"])
  
  return best_p

        



def parse_command_line(): 
    """Parse the command line.
    Return the result provided by argparse.ArgumentParser.parse_args()
    """
    
    parser = argparse.ArgumentParser(description="Tune P2NFFT parameters for a set of particle configurations with reference forces and torques")
    
    parser.add_argument("--target_accuracy",dest="target_accuracy",
        type=float, required=True,
        help="Accuracy to which the method is tuned. The mean of rms forces and torqeus error is used.")
    parser.add_argument("--mesh",dest="meshes",
        action="append",type=int, required=True,
        help="Number of fft mesh sites for the longest box dimension")
    parser.add_argument("--assignment_order",dest="assignment_order",
        type=int, required=True,
        help="P2NFFT assignment order")
    parser.add_argument("--min_r_cut",dest="min_r_cut",
        type=float, required=True,
        help="Minimum real space cutoff to consider")
    parser.add_argument("--max_r_cut",dest="max_r_cut",
        type=float, required=True,
        help="Maximum real space cutoff to consider")
    parser.add_argument("--r_cut_tolerance",dest="r_cut_tolerance",
        type=float, required=True,
        help="Tolerance for determining the lowest possible real space cutoff")
    parser.add_argument("--output_file",dest="output_file",
        required=True,
        help="Name of the file to write tuning result to ")
    parser.add_argument("particle_configurations",
        nargs="+",
        help="Particle configurations with reference data")
    
    
    args=parser.parse_args()

    # Verify numbers >0
    for arg in "target_accuracy","assignment_order","min_r_cut","max_r_cut","r_cut_tolerance":
        if not getattr(args,arg) >0:
            print(arg,"ust be >0")
            sys.exit(1)
    
    for mesh in args.meshes:
        if not mesh >0:
            print("Mesh must be >0")
            sys.exit(1)

    # Verify cutoff parameters
    if not args.min_r_cut <args.max_r_cut:
        print("min_r_cut must be < max_r_cut")
        sys.exit(1)
    if not args.r_cut_tolerance <args.max_r_cut-args.min_r_cut:
        print("r_cut_tolerance must be < max_r_cut - min_r_cut")
        sys.exit(1)
    return args


# Command line
args=parse_command_line()    



# prepare Espresso sytem handle
s=setup_es()


# Sotrage for tuning results
all_results=[]

# Minimum and maximum cutoffs for the bisection
min_r_cut=args.min_r_cut
max_r_cut=args.max_r_cut



# Tune for highest mesh first to find lower bound for the real space cutoff
meshes=sorted(args.meshes)
mesh=meshes[-1]

# Tune
result=tune_r_cut(s,mesh=mesh,min_r_cut=min_r_cut,max_r_cut=max_r_cut, r_cut_tolerance=args.r_cut_tolerance, target_accuracy=args.target_accuracy,assignment_order=args.assignment_order,particle_configurations=args.particle_configurations)

# Can we tune for the largest mesh?
if result is None:
    print("Accuracy cannot be met. Increase mesh or max_r_cut")
    sys.exit(1)

# Store result
all_results.append(result)
# The r_cut for the highest mesh is the lower bound for smaller meshes
min_r_cut=result["r_cut"]
meshes.remove(mesh)

# Tune remaining meshes (sorted ascending)
while len(meshes)>0:
    # Consume lowest mesh
    mesh=meshes[0]
    meshes.remove(mesh)
    
    # Tune
    result=tune_r_cut(s,mesh=mesh,min_r_cut=min_r_cut,max_r_cut=max_r_cut, r_cut_tolerance=args.r_cut_tolerance, target_accuracy=args.target_accuracy,assignment_order=args.assignment_order,particle_configurations=args.particle_configurations)
    
    # Could we tune for this mesh?
    if not result is None:
        # Store result
        all_results.append(result)
        
        # The r_cut for this mesh is the upper bound for all higher meshes
        max_r_cut=result["r_cut"]

# Select fastest p2nfft configuration
all_results=sorted(all_results, key=lambda x: x["time"])

# Show summary
for r in all_results:
  print("mesh:",r["mesh"],"r_cut:",r["r_cut"],"time:",r["time"])

# Store fastest configuration
fastest=all_results[0]
f=open(args.output_file,"w")
pickle.dump(fastest,f)
f.close()

