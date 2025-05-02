#!/bin/bash -l
#SBATCH --nodes=1 
#SBATCH --time=150:00:00   
#SBATCH --partition=smp-rh7
#SBATCH --cpus-per-task=8 
#SBATCH --mem=46gb


module load stacks/2.65 

denovo_map.pl --samples ../sh_samples/ --popmap ../popmap_sh.txt --out-path ../denovom3_M1_n1/ --paired -m 3 -M 1 -n 1 -T 16 --catalog-popmap ./catalog_sh.txt -X "populations: --vcf"

