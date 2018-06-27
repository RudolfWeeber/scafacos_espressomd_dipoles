mpirun=mpirun
pypresso=../es/pypresso

set -e -x -v
# Generate thermalized system
$mpirun -np 8 $pypresso gen-system.py  --output_file full_system.pcl --steps 150 

# Add reference forces and torques via direct summation
$pypresso add-reference-forces-torques.py  full_system.pcl 

# Tune p2nfft parameters
$mpirun -np 8 $pypresso tune.py  --target_accuracy 1E-4 --mesh 96  --mesh 128 --mesh 160 --mesh 192 --min_r_cut 8 --max_r_cut 18 --r_cut_tolerance 0.25 --assignment_order 5  --output full_p2nfft_config.pcl full_system.pcl

# Run a simulation with p2nfft. External field 0
$mpirun -np 8 $pypresso run.py --particle_config_file full_system.pcl --p2nfft_config_file full_p2nfft_config.pcl --field 0 --output_base full_out --start_step 0 --end_step 10001 
# Verify accuracy of final result
$pypresso add-reference-forces-torques.py  full_out_0_10000.pcl
$mpirun -np 8 $pypresso get-accuracy.py --required_accuracy 2E-3 --particle_config_file full_out_0_10000.pcl --p2nfft_config_file full_p2nfft_config.pcl


# Simulation for field 10
$mpirun -np 8 $pypresso run.py --particle_config_file full_system.pcl --p2nfft_config_file full_p2nfft_config.pcl --field 10 --output_base full_out --start_step 0 --end_step 5001
# Retune after 5000 steps
$pypresso add-reference-forces-torques.py  full_out_10_5000.pcl
$mpirun -np 8 $pypresso tune.py  --target_accuracy 1E-4 --mesh 160 --mesh 192 --mesh 256 --mesh 320 --mesh 384 --min_r_cut 10 --max_r_cut 25 --r_cut_tolerance 0.25 --assignment_order 5  --output full_p2nfft_config_field_10.pcl full_out_10_5000.pcl
# Simulate remaining 5000 steps
$mpirun -np 8 $pypresso run.py --particle_config_file full_out_10_5000.pcl --p2nfft_config_file full_p2nfft_config_field_10.pcl --field 10 --output_base full_out --start_step 5000 --end_step 10001

# Verify accuracy of final result
$pypresso add-reference-forces-torques.py  full_out_10_10000.pcl
$mpirun -np 8 $pypresso get-accuracy.py --required_accuracy 2E-3 --particle_config_file full_out_10_10000.pcl --p2nfft_config_file full_p2nfft_config_field_10.pcl



