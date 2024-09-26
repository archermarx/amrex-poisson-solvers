#!/bin/bash
#SBATCH --job-name=amrex-poisson-install
#SBATCH --account=ramanvr
#SBATCH --partition=ramanvr
#SBATCH --gpus=1
#SBATCH --cpus-per-gpu=20
#SBATCH --mem-per-cpu=4g
#SBATCH --time=2:00:00
#SBATCH --output=output.log
#SBATCH --error=output.log
#SBATCH --mail-type=END,FAIL

# Load required modules
module purge
module load gcc/10.3.0
module load cuda/12.1.1
module load cmake/3.26.3
module load openblas/0.3.23
module load openmpi/4.1.6-cuda

cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=RelWithDbgInfo \
  -DAMReX_MPI=OFF \
  -DAMReX_SPACEDIM=3 \
  -DAMReX_TINY_PROFILE=ON \
  -DAMReX_BASE_PROFILE=OFF

cmake --build build -j 20

