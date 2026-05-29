#!/bin/bash
#SBATCH --job-name=PO
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=60G
#SBATCH --time=160:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=xuanhong@princeton.edu
module load openmpi/gcc/4.1.0
module load fftw/gcc/3.3.9

srun ./PO_NG_Adjoint_FS_L_BFGS_B_memory
