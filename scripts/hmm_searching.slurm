#!/bin/bash 

#SBATCH --job-name=hmm_searching_array
#SBATCH --mail-type=END
#SBATCH --mail-user=jgh8@kent.ac.uk
#SBATCH --partition=biosoc2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=14000
#SBATCH --error=./hmm_searching_array/%a.out
#SBATCH --out=./hmm_searching_array/%a.out
#SBATCH --time=02:00:00
#SBATCH --array=1-64%64

# conda needs special stuff to start via slurm
source /shared/home/jgh8/miniconda3/etc/profile.d/conda.sh

conda activate hmmer

# get all hmms in an array
## should be ~64
dir="../data/hmm_building/"
all_hmms=($(find ${dir} -type f -name "train*hmm"))
filenumber=$(($SLURM_ARRAY_TASK_ID - 1))

filepath=${all_hmms[$filenumber]}
gene=$(basename "$(dirname "$filepath")")

filename=${filepath##*/} #remove path up to filename start
file=${filename%.*}      # remove extension
repeat=${file##*_}       # get number

# make test dataset just for this 
test_loci_filepath=$(find ${dir}/${gene} -type f -name "test*${repeat}.txt")
test_loci_filename=${test_loci_filepath##*/}
test_loci_file=${test_loci_filename%.*}

# extract test loci into a file
echo "Extracting loci for ${test_loci_filepath}..."
seqkit grep -f ${test_loci_filepath} ${dir}/proteome.faa -o ${dir}/${gene}/${test_loci_file}.faa

# run hmm against it
echo "Searching ${filepath} against ${dir}/${gene}/${test_loci_file}.faa..."
hmmsearch -o /dev/null \
	      --tblout ${dir}/${gene}/${test_loci_file}_results.out \
	      --cpu 6 \
  	      ${filepath} \
	      ${dir}/${gene}/${test_loci_file}.faa

# delete old files
echo "Tidying up..."
rm ${dir}/${gene}/${test_loci_file}.faa
echo -e "Complete!\n"

