#!/bin/bash
# Script to build lots of hmms from lists of loci IDs

loci_names_dir="../data/jackhmmer/loci_names/*"

for file in $loci_names_dir
do
	filename="${file##*/}"
	gene="${filename%.*}"

	echo "Working on ${gene}"

	# pull out sequences
	seqkit grep -f ${file} ../data/jackhmmer/proteome.faa -o ../data/jackhmmer/loci_sequences/${gene}.faa
	
	# make MSA
	mafft --thread -1 --auto ../data/jackhmmer/loci_sequences/${gene}.faa > ../data/jackhmmer/sequence_alignments/${gene}.msa

	# convert to stockholm
	msaconverter.py -i ../data/jackhmmer/sequence_alignments/${gene}.msa -o ../data/jackhmmer/sequence_alignments/${gene}.sto -p fasta -q stockholm -t protein

	# make hmms
	hmmbuild --amino --cpu 1 .../data/jackhmmer/hmms/${gene}.hmm .../data/jackhmmer/sequence_alignments/${gene}.sto

	echo -e "${gene} complete!\n"

done

echo "making hmm database"
cat ../data/jackhmmer/hmms/*.hmm > ../data/jackhmmer/all_hmm_database.hmm

hmmpress ../data/jackhmmer/all_hmm_database.hmm

echo "hmm building complete, searching proteome against database"
hmmscan -o ../data/jackhmmer/hmmsearch.out --tblout ../data/jackhmmer/hmmsearch.txt --cpu 16 ../data/jackhmmer/all_hmm_database.hmm ../data/jackhmmer/proteome.faa 

echo "Script complete, see ../data/jackhmmer/hmmsearch.txt for results"