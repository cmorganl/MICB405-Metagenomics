#!/bin/bash

if [ $# -lt 6 ]; then
        echo -e "\n\tUSAGE: assembly_list.txt reads_list.txt interleaved[y|n] output_dir prefix_name n_threads\n"
	echo "'assembly_list.txt' is a file listing the FASTA files containing contigs to be binned."
	echo "                    They will be concatenated, length-filtered and deduplicated prior to binning."
	echo "'reads_list.txt'    is a file listing the FASTQ files to be used for generating the abundance profiles."
	echo "                    Interleaved FASTQ will need a single file on one line; split FASTQ needs R1 and R2 tab-delimited."
	echo "                    A split-FASTQ file list can be created with '$ ls /path/to/*/trim/*pe.?.fq.gz | paste - - >reads_list.txt'"
	echo "'interleaved'       specifies whether the FASTQ files of each sample are split (n) or interleaved (y)."
	echo "'output_dir'        is the directory for writing output files."
	echo "'prefix_name'       is the prefix for naming output files. Must be less than 20 characters."
	echo "'n_threads'         is the number of threads to be used by bwa, samtools, and MetaBAT 2.11.2."
        exit
fi

input_contig=$1
reads_list=$2
fq_form=$3
prefix_dir=$4
sid=$5
# f_ext=$( echo $input_contig | awk -F. '{ print $NF }')
# sid=$( basename $input_contig .$f_ext)
n_thread=$6

min_length=1500
dedup_id=99
log_f=$prefix_dir/${sid}_gt${min_length}_BINME_log.txt
abunds=$prefix_dir/${sid}_gt${min_length}_abund_list.txt
bin_input_contigs=$prefix_dir/${sid}_gt${min_length}_dedup.fasta
output_dir=$prefix_dir/MetaBAT2_$sid/
bad_dir=$output_dir/bad_name_bins
rpkm_f=$prefix_dir/${sid}_binned.rpkm.csv
header_map=$output_dir/binned_header_map.tsv
checkm_stats=$output_dir/MetaBAT2_${sid}_min${min_length}_checkM_stdout.tsv

# Ensure the prefix is not too long - essential for compatibility with Prokka
if [ $(echo $sid| wc -c ) -gt 21 ]; then
	echo "ERROR: 'prefix_name' ($sid) must be shorter than 20 characters! Exiting now."
	exit
fi

# Ensure formatting of FASTQ files is correct
if [ $fq_form != 'y' ] && [ $fq_form != 'n' ]; then
	echo "ERROR: Expected a 'y' or 'n' for the interleaved argument. Cannot recognize '$fq_form'! "
	exit
fi
while read line
do
	nfq=$( echo "$line" | awk --field-separator="\t" '{ print NF }' )
	if [ $nfq -eq 2 ] && [ $fq_form == 'y' ]; then
		echo "ERROR: 2 FASTQs found when only one is expected in $reads_list."
		exit
	elif [ $nfq -eq 1 ] && [ $fq_form == 'n' ]; then
		echo "ERROR: 1 FASTQ found when 2 are expected in $reads_list."
		exit
	fi
done<$reads_list

printf "\n################################################## BIN_MULTI ##################################################\n"

printf "\nPARAMETERS:\n"
printf "\tPREFIX = $sid \n"
printf "\tOUTPUT DIRECTORY = $prefix_dir \n"
printf "\tMIN_CONTIG_LENGTH = $min_length \n"
printf "\tLOG = $log_f \n"
printf "\tTHREADS = $n_thread \n"
printf "\tMIN_DEDUP_ID = $dedup_id \n"
echo

if [ ! -d $prefix_dir ]; then
	mkdir $prefix_dir
fi

if [ -f $log_f ]; then
	rm $log_f
	touch $log_f
fi 

###################################################################################################
# Prepare the reference metagenome for binning by concatenating, trimming and deduplicating
###################################################################################################
if [ ! -f $bin_input_contigs ]; then
	printf "[STATUS] Concatenating FASTA files... "
	# concatenate_FAs.py ensures there are no duplicate contig names
	concatenate_FAs.py -l $input_contig -o $prefix_dir/$sid.fasta
	printf "done.\n"
	
	printf "[STATUS] Removing all sequences shorter than ${min_length}... "
	splitFASTA -i $prefix_dir/$sid.fasta -m $min_length -o $prefix_dir/$sid\_gt$min_length.fasta
	printf "done.\n"

	printf "[STATUS] Deduplicating sequences at ${dedup_id} percent similarity... "
	dedupe.sh -Xmx20g in=$prefix_dir/$sid\_gt$min_length.fasta out=$bin_input_contigs threads=$n_thread minidentity=99 >>$log_f 2>&1
	printf "done.\n"
	printf "[STATUS] Removing intermediate FASTA files... "
	rm $prefix_dir/$sid.fasta $prefix_dir/$sid\_gt$min_length.fasta
	printf "done.\n"
else
	echo "[INFO] Using $bin_input_contigs for binning"
fi

##################################################

if [ ! -f $bin_input_contigs ]; then
	echo "ERROR: $bin_input_contigs was not created!"
	exit
fi

###################################################################################################
# Calculate each contig's abundance across samples for MetaBAT
###################################################################################################
if [ ! -f $abunds ]; then	
	printf "[STATUS] Indexing $bin_input_contigs ... "
	bwa index $bin_input_contigs 1>$prefix_dir/bwa_stdout.txt 2>$prefix_dir/bwa_stderr.txt
	printf "done.\n"

	while read line; do
		fq_sid=$( basename $line | sed 's/.fastq.*\|.fq.*//g' )
		printf "[STATUS] Aligning reads in $fq_sid ... "
		if [ $fq_form == 'y' ]; then
			bwa mem -t $n_thread -p $bin_input_contigs $line 2>$prefix_dir/bwa_stderr.txt | \
				samtools view -@ 4 -bS - | samtools sort -@ 4 - 1>$bin_input_contigs.$fq_sid.bam 2>>$log_f
		else
			fwfq=$( echo "$line" | awk -F"\t" '{ print $1 }' )
			refq=$( echo "$line" | awk -F"\t" '{ print $2 }' )
			bwa mem -t $n_thread $bin_input_contigs $fwfq $refq 2>$prefix_dir/bwa_stderr.txt | \
				samtools view -@ 4 -bS - | samtools sort -@ 4 - 1>$bin_input_contigs.$fq_sid.bam 2>>$log_f
		fi
		printf "done.\n"
	done<$reads_list
	printf "[STATUS] Building the abundance profiles... "
	jgi_summarize_bam_contig_depths --outputDepth $abunds $bin_input_contigs.*.bam >>$log_f 2>&1
	# Remove the bwa index files and BAMs
	rm $bin_input_contigs.*
	printf "done.\n"
fi

##################################################

if [ ! -f $abunds ]; then	
	echo "ERROR: Sample-wise abundance file $abunds was not created!"
	exit
else
	echo "[INFO] Using the abundance profiles found in $abunds"
fi

###################################################################################################
# Bin the reference metagenome using MetaBAT2
###################################################################################################

if [ ! -d $output_dir ]; then
	printf "[STATUS] Binning with MetaBAT 2.11.2... "
	mkdir $output_dir
	echo "metabat2 --inFile $bin_input_contigs --outFile $output_dir/$sid\_raw --abdFile $abunds --numThreads $n_thread --seed 888" >>$log_f
	metabat2 \
--inFile $bin_input_contigs \
--outFile $output_dir/$sid\_raw \
--abdFile $abunds \
--numThreads $n_thread \
--seed 888 >>$log_f 2>&1

	printf "done.\n"

##################################################

###################################################################################################
# Re-assemble each of the MAGs using minimus2 to further reduce sequence redundancy
###################################################################################################
	printf "[STATUS] Re-assembling MAGs with minimus2... "
	for asm in $output_dir/$sid\_raw*fa;
	do
		prefix=$( echo $asm | sed 's/_raw/_minimus2/g' | sed 's/.fa//g' )
		final_ff=$( echo $asm | sed 's/_raw//g' )
		echo "~/bin/amos-3.1.0/bin/toAmos -s $asm -o $prefix.afg" >>$log_f
		~/bin/amos-3.1.0/bin/toAmos -s $asm -o $prefix.afg >>$log_f
		echo "~/bin/amos-3.1.0/bin/minimus2 $prefix -D REFCOUNT=0 -D OVERLAP=200 -D MINID=95" >>$log_f
		~/bin/amos-3.1.0/bin/minimus2 $prefix -D REFCOUNT=0 -D OVERLAP=200 -D MINID=95 >>$log_f
		echo "cat $prefix.fasta $prefix.singletons.seq >$final_ff" >>$log_f
		cat $prefix.fasta $prefix.singletons.seq >$final_ff
		rm -r $prefix.runAmos.log $prefix.afg $prefix.OVL $prefix.singletons $prefix.contig $prefix.ovl $prefix.singletons.seq $prefix.coords $prefix.qry.seq $prefix.delta $prefix.bnk $prefix.ref.seq $prefix.fasta
	done
	printf "done.\n"

##################################################


###################################################################################################
# Archive the raw bins generated by MetaBAT2
###################################################################################################
	printf "[STATUS] Tar-balling raw MetaBAT MAGs... "
	mkdir $output_dir/Raw_MAGs/
	mv $output_dir/$sid\_raw.*fa $output_dir/Raw_MAGs/
	cd $output_dir/
	tar -czf $output_dir/Raw_MAGs.tar.gz Raw_MAGs/
	rm -r $output_dir/Raw_MAGs/
	printf "done.\n"
	cd -
fi
##################################################

###################################################################################################
# Run checkM to evaluate the completeness and contamination of the MAGs, summarise
###################################################################################################
if [ ! -f $checkm_stats ]; then
	echo "checkm lineage_wf --tab_table -x .fa --threads $n_thread --pplacer_threads $n_thread $output_dir $output_dir/checkM_output/ >$checkm_stats" >>$log_f
	checkm lineage_wf \
--tab_table \
-x .fa \
--threads $n_thread \
--pplacer_threads $n_thread \
$output_dir $output_dir/checkM_output/ 1>$checkm_stats
fi

num_hq=$( gawk -F"\t" '{ if ($12>90 && $13<5) print $0 }' $checkm_stats | wc -l)
echo "[INFO] Number of bins with >90%% completion and <5%% contamination: $num_hq"
num_mpq=$( gawk -F"\t" '{ if ($12>50 && $13<10) print $0 }' $checkm_stats | wc -l)
echo "[INFO] Number of bins with >50%% completion and <10%% contamination: $num_mpq"

##################################################

###################################################################################################
# Fix the header names in the bin .fa files (added by Ryan 10.31.2018)
###################################################################################################
printf "[STATUS] Renaming MAG headers to common format... "
echo "Renaming MAG headers to common format... " >>$log_f

mkdir $bad_dir
mv $output_dir*.fa $bad_dir/

for f in $bad_dir/*.fa;
do
	fname=$( basename $f )
	bin_hmap=$output_dir/$( basename $fname .fa )\_map.tsv
	python3 /usr/local/bin/rename_contig_headers_after_binning.py $f $output_dir/$fname >>$log_f
	if [ ! -f $bin_hmap ]; then
		echo "ERROR: Header map $bin_hmap does not exist!"
		exit
	else
		cat $bin_hmap >>$header_map
		rm $bin_hmap
	fi
done
rm -r $bad_dir/
printf "done.\n"

##################################################


###################################################################################################
# Calculate RPKM for all binned contigs
###################################################################################################
printf "[STATUS] Creating BWA index for $output_dir/binned_sequences.fasta... "
cat $output_dir/*fa >$output_dir/binned_sequences.fasta
if [ ! -f $output_dir/binned_sequences.fasta ]; then
	echo "ERROR: Unable to concatenate bin FASTA files in $output_dir/"
	exit
fi
bwa index $output_dir/binned_sequences.fasta 1>$prefix_dir/bwa_stdout.txt 2>$prefix_dir/bwa_stderr.txt
printf "done.\n"

touch $rpkm_f
echo "Sample,Sequence,RPKM" >>$rpkm_f
while read line; do
	fq_sid=$( basename $line | sed 's/.fastq.*\|.fq.*//g' )
	printf "[STATUS] Aligning reads in $fq_sid ... "
	if [ $fq_form == 'y' ]; then
	bwa mem -t $n_thread -p $output_dir/binned_sequences.fasta $line 1>$prefix_dir/binned_sequences.$fq_sid.sam 2>$prefix_dir/bwa_stderr.txt
	else
		fwfq=$( echo "$line" | awk -F"\t" '{ print $1 }' )
		refq=$( echo "$line" | awk -F"\t" '{ print $2 }' )
	bwa mem -t $n_thread $output_dir/binned_sequences.fasta $fwfq $refq 1>$prefix_dir/binned_sequences.$fq_sid.sam 2>$prefix_dir/bwa_stderr.txt
	fi
	printf "done.\n"
	printf "[STATUS] Calculating RPKM... "
	rpkm -c $output_dir/binned_sequences.fasta -a $prefix_dir/binned_sequences.$fq_sid.sam -o $prefix_dir/$sid\_binned.$fq_sid.rpkm.csv 1>>$log_f 2>>$log_f
	cat $prefix_dir/$sid\_binned.$fq_sid.rpkm.csv | sed "s/^/$fq_sid,/g" >>$rpkm_f
	rm $prefix_dir/$sid\_binned.$fq_sid.rpkm.csv
	printf "done.\n"
done<$reads_list
rm $prefix_dir/binned_sequences.*.sam $output_dir/binned_sequences.fasta*

###################################################################################################
# Split the wheat from the chaff, archiving the low-quality MAGs
###################################################################################################
printf "[STATUS] Archiving low quality MAGs... "
cd $output_dir
mkdir MedQPlus_MAGs/
mkdir LowQ_MAGs/
for f in $(gawk -F"\t" '{ if ($12>50 && $13<10) print $1 }' $checkm_stats)
do
	mv $f.fa MedQPlus_MAGs/
done
mv *.fa LowQ_MAGs/
tar -czf LowQ_MAGs.tar.gz LowQ_MAGs/
rm -r LowQ_MAGs/
printf "done.\n"
cd -

##################################################

###################################################################################################
# Use the Genome-taxonomy database toolkit (gtbdtk) to phylogenetically classify the MAGs
###################################################################################################
gtdbtk classify_wf --genome_dir $output_dir/MedQPlus_MAGs/ --out_dir $output_dir/gtdbtk_output/ -x .fa --cpus $n_thread

if [ ! -f $output_dir/gtdbtk_output/gtdbtk.ar*.classification_pplacer.tsv ] && [ ! -f $output_dir/gtdbtk_output/gtdbtk.bac*.classification_pplacer.tsv ]; then
	echo "ERROR: gtdbtk did not run successfully. No classification files were found in $output_dir/gtdbtk_output/"
	exit
fi
##################################################


file_list=$(gawk -F"\t" '{ if ($12>50 && $13<10) print $1 }' $checkm_stats | sed -r 's/$/.fa /g')
###################################################################################################
# Run prokka on each HQB (added by Ryan 10.31.2018) using the GTDB-tk Kingdom annotation
###################################################################################################
prokka_output_dir=$output_dir/prokka_annos
IFS=' ' read -r -a array <<< "$file_list"
for f in "${array[@]}";
do
	bin_num=$(cut -d'.' -f2 <<< $f)
	kingy=$(grep -w "$bin_num" $output_dir/gtdbtk_output/gtdbtk.*.classification_pplacer.tsv | gawk '{ print $2 }' | gawk -F';' '{ print $1 }' | sed 's/d__//g')
	printf "[STATUS] Annotating $f with ${kingy}l info using Prokka... "
	prokka \
	--kingdom $kingy \
	--outdir $prokka_output_dir/$bin_num/ \
	--force \
	$output_dir/$f
	printf "done.\n"
done

##################################################


