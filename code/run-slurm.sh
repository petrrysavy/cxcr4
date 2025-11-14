#!/bin/bash
#SBATCH --mem 96G
#SBATCH -p gpulong
#SBATCH --gres gpu:1

ml cuDNN/9.5.0.50-CUDA-12.4.0
#ml jax/0.4.35-gfbf-2023b-CUDA-12.4.0
#ml jax/0.4.26-gfbf-2023b-CUDA-12.4.0

source ~/.conda/envs/it/bin/activate

python3 --version
~/.conda/envs/it/bin/python3 --version

~/.conda/envs/it/bin/python3 main.py
