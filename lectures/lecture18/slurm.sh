#!/bin/bash
#SBATCH --mail-user=fred@fredhutch.org
#SBATCH --mail-type=END
#SBATCH --array=[1-15]
#SBATCH --nodes=1
#SBATCH --output=Rout/par-%J.out
#SBATCH --error=Rout/par-%J.err
#SBATCH --cpus-per-task=1

# Title/description here
# Usage instructions
