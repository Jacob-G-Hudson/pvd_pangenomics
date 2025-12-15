
# conda needs special stuff to start via slurm
source /shared/home/jgh8/miniconda3/etc/profile.d/conda.sh

conda activate NRPS

# setup
## move into NRPS_transformer dir as lots of its sub-scripts need to to run
cd "../programs/NRPSTransformer"

in_dir="../../data/NRPS_transformer/run_3"
out_dir="./run_3"
all_loci=($(find ${in_dir} -type f -name "*nrps_loci.txt"))

# get loci ID for this array number
for nrps_loci in "${all_loci[@]}"
do
	
	# get bateria name
	full_filename="${nrps_loci##*/}"
	bacteria="${full_filename%_nrps_loci.txt}"

	echo "Working on ${bacteria}"

	# Check if already completed
	if [ -f ${out_dir}/results/${bacteria}.csv ]; then
		echo "$bacteria already completed! Skipping..."
		continue	
	fi
	

	# pull out sequences
	seqkit grep -f ${nrps_loci} ../../data/jackhmmer/proteome.faa -o ${out_dir}/loci_sequences/${bacteria}.faa

	# run NRPS transformer (note this is entirely copied from the NRPS transformer run.sh script)
	input_file=${out_dir}/loci_sequences/${bacteria}.faa

	## Check if the specified input file exists
	if [ ! -f "$input_file" ]; then
		echo "Error: The input file '$input_file' does not exist."
		exit 1
	fi
	echo setting input file to "$input_file"

	## Ask user to input the result file path
	# read -p "Please enter the result file path (Default: $DEFAULT_RESULT_PATH): "
	result_path="${out_dir}/results/${bacteria}.csv"
	mkdir -p "${out_dir}/results/"

	echo "setting result path to $result_path"

	## First Run HMM Scan
	if ! ./hmm/scan.sh "$input_file"; then
		echo "Error: HMM scan failed. Please check the input file and parameters."
		exit 1
	fi

	## Then Parse dbtl
	if ! python ./hmm/parse_dbtl.py; then
		echo "Error: Parsing dbtl failed. Please check the output of HMM scan."
		exit 1
	fi

	## Then generate the domains
	if ! python ./hmm/gen_domains.py --fasta_path "$input_file"; then
		echo "Error: Generating domains failed. Please check the parsed dbtl file."
		exit 1
	fi

	## Finally, run the final AI prediction
	if ! python inference.py --result_path "$result_path"; then
		echo "Warning: Final AI prediction failed. Please check the input data and parameters."
	fi

	## remove input file to save space
	rm ${out_dir}/loci_sequences/${bacteria}.faa

	echo "Successfully finished all steps!"
	echo "The result is saved to $result_path"

done