#!/bin/bash

#SBATCH -o sim1_uniform_150.log
#SBATCH -e sim1_uniform_150.err
#SBATCH -p mmmi -q mmmi
#SBATCH --mem=8000
#SBATCH --mail-user=maylis.layan@pasteur.fr --mail-type=END

module unload beast/v1.10.4
module load beast/v2.6.1

srun beast -beagle -beagle_CPU -statefile sim1_uniform_150.states sim1_uniform_150.xml

exit 0