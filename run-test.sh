mpirun=mpirun
pypresso=../es/pypresso

set -e -x -v
# Generate thermalized system
$mpirun -np 8 $pypresso gen-system.py  --output_file test_system.pcl --steps 3 

# Add reference forces and torques via direct summation
$pypresso add-reference-forces-torques.py  test_system.pcl 

# Tune p2nfft parameters
$mpirun -np 8 $pypresso tune.py  --target_accuracy 1E-3 --mesh 96  --mesh 128 --min_r_cut 8 --max_r_cut 12 --r_cut_tolerance 1 --assignment_order 5  --output test_p2nfft_config.pcl test_system.pcl

# Run a simulation with p2nfft
$mpirun -np 8 $pypresso run.py --particle_config_file test_system.pcl --p2nfft_config_file test_p2nfft_config.pcl --field 0 --output_base test_out --start_step 0 --end_step 11

# Verify the accuracy for a particle configuration produced by the simulation
$pypresso add-reference-forces-torques.py  test_out_0_0010.pcl
$mpirun -np 8 $pypresso get-accuracy.py --required_accuracy 2E-3 --particle_config_file test_out_0_0010.pcl --p2nfft_config_file test_p2nfft_config.pcl



