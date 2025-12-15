#!/bin/bash 

#SBATCH --job-name=final_hmm_array
#SBATCH --mail-type=END
#SBATCH --mail-user=jgh8@kent.ac.uk
#SBATCH --partition=biosoc2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=14000
#SBATCH --error=./final_hmm_array/%a.out
#SBATCH --out=./final_hmm_array/%a.out
#SBATCH --time=02:00:00
#SBATCH --array=1-16

# conda needs special stuff to start via slurm
source /shared/home/jgh8/miniconda3/etc/profile.d/conda.sh

conda activate hmmer

# setup
dir="../data/hmm_building/"
all_loci=($(find ${dir} -type f -name "all_loci.txt"))

# get loci ID for this array number
filenumber=$(($SLURM_ARRAY_TASK_ID - 1))
hmm_loci=${all_loci[$filenumber]}

# Make hmm	
gene=$(basename "$(dirname "$hmm_loci")")
full_filename="${hmm_loci##*/}"
filename="${full_filename%.*}"

echo "Working on ${gene}"

## pull out sequences
seqkit grep -f ${hmm_loci} ${dir}/proteome.faa -o ${dir}/${gene}/${gene}.faa

## make MSA
mafft --thread -4 --auto ${dir}/${gene}/${gene}.faa > ${dir}/${gene}/${gene}.msa

## convert to stockholm
msaconverter -i ${dir}/${gene}/${gene}.msa -o ${dir}/${gene}/${gene}.sto -p fasta -q stockholm -t protein

## make hmms
hmmbuild --amino --cpu 4 ${dir}/${gene}/${gene}.hmm ${dir}/${gene}/${gene}.sto
echo -e "${gene} complete!\n"

# search against proteome 
hmmsearch -o /dev/null \
	      --tblout ${dir}/${gene}/results.out \
	      --cpu 8 \
  	      ${dir}/${gene}/${gene}.hmm \
	      ${dir}/proteome.faa	

