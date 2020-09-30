#!/bin/bash
#SBATCH -o .log
#SBATCH -e .err
#SBATCH -p mmmi -q mmmi
#SBATCH --mem=5000
#SBATCH --mail-user=maylis.layan@pasteur.fr --mail-type=END

logcombiner -trees -burnin 2000000 sim1_uniform_150.trees.txt sim1_uniform_150_noburnin.trees.txt
treeannotator -heights median sim1_uniform_150_noburnin.trees.txt sim1_uniform_150.mcc.tree
rm sim1_uniform_150_noburnin.trees.txt

exit 0
