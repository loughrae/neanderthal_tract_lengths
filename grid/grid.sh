#!/bin/sh
#SBATCH -J batchTLnov2021
#SBATCH -t 12:00:00
#SBATCH --array=1-208
#SBATCH --account=ccmb-condo
#SBATCH --mem-per-cpu=8000
module load anaconda/3-5.2.0
ID=$SLURM_ARRAY_TASK_ID

echo "hello ${ID}"
echo "sel t1 t2" > logit_results_nov_${ID}.txt
echo "sel t1 t2" > numeric_results_nov_${ID}.txt #these are mainly to replace existing files with this name in case I run the program multiple times (needed because I append later on). Would just remove
#them but that's messy in the case where they don't already exist

while read -r sel t1 t2
do

logit=$(python3 ../../../eloughra/vladjan21_2/DAIM/tract_length.py logit -p 0 None $sel -p $t1 0.03 0 -p $t2 0.03 0 2>&1) #think this should put errors in the output table as well
echo "${logit/$'\n'} ${sel} ${t1} ${t2}" >> logit_results_nov_${ID}.txt

delay="$(($t2-$t1))"
python3 ../../../eloughra/vladjan21_2/DAIM/precomp_traj.py 0.03 $sel 10000 $t2 -Tn $delay > temp_AFs_nov_${ID}.txt #to avoid overlaps from the parallel
numerics=$(python3 ../../../eloughra/vladjan21_2/DAIM/tract_length.py  precomp temp_AFs_nov_${ID}.txt 2>&1)
echo "${numerics/$'\n'} ${sel} ${t1} ${t2}" >> numeric_results_nov_${ID}.txt

done < allcombos_nov2021_${ID}.txt
