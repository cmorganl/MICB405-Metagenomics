#!/bin/bash

if [ $# -eq 0 ]; then
	echo "USAGE: $0 project_path sample threads interleaved[y|n] gzipped[y|n]"
	echo "To use the pipeline, create a \$project_path dir (e.g., /home/cmorganlang/path_to_project/Coal_beds) and within it a \$sample dir (e.g., CO182)."
	echo "Then deposit all forward sense reads into /path/to/\$project_path/\$sample/forward_dir and reverse strand reads in /path/to/\$project_path/\$sample/reverse_dir"
	echo "If the reads are interleaved, leave them in the '\$project_path' directory as \$sample.fastq"
	echo "Finally, run the script as described by USAGE"
	exit
fi

project_path=$1
sample=$2
threads=$3
interleaved=$4
gz=$5

project_dir=$project_path/$sample/
if [ ! -d $project_dir ]; then
	mkdir $project_dir
fi

forward_dir=$project_dir/forward_dir
reverse_dir=$project_dir/reverse_dir
trim_dir=$project_dir/trim
assembly_dir=$project_dir/assembly

Forward_fq=$forward_dir/$sample\_R1.fastq
Reverse_fq=$reverse_dir/$sample\_R2.fastq

if [ $interleaved == "y" ]; then
       	printf "Deinterleaving files... "
        interleaved_input=$project_dir/$sample.fastq
        if [ ! -f $interleaved_input ] && [ $gz == 'n' ]; then
                echo "ERROR: Unable to find interleaved fastq $interleaved_input"
                exit
        elif [ ! -f $interleaved_input.gz ] && [ $gz == 'y' ]; then
                echo "ERROR: Unable to find interleaved fastq $interleaved_input.gz"
                exit
        fi
        if [ ! -d $forward_dir ]; then
                mkdir $forward_dir
                mkdir $reverse_dir
        fi
        if [ $gz == "y" ]; then
                zcat $interleaved_input.gz |paste - - - - - - - - | tee >(cut -f 1-4 | tr "\t" "\n" > $Forward_fq) | cut -f 5-8 | tr "\t" "\n" > $Reverse_fq
        fi
        if [ $gz == "n" ]; then
                cat $interleaved_input |paste - - - - - - - - | tee >(cut -f 1-4 | tr "\t" "\n" > $Forward_fq) | cut -f 5-8 | tr "\t" "\n" > $Reverse_fq
        fi
        printf "done.\n"
        gz="n"
fi

num_forward=$( ls $forward_dir | wc -l)
num_reverse=$( ls $reverse_dir | wc -l)

if [ $num_forward == 0 ] || [ $num_reverse == 0 ]; then
	echo "ERROR: FASTQ files are not in their proper place!"
	echo "This pipeline expects gzipped fastq files to be in $forward_dir and $reverse_dir"
	echo "If the fastq is interleaved, it should be placed in $project_dir/interleaved"
	exit
fi

if [ -d $assembly_dir ]; then
	rm -r $assembly_dir
fi
mkdir $assembly_dir

cd $project_dir

if [ ! -f $Forward_fq.gz ] && [ $gz == 'y' ]; then
	#There is not a single gzipped file and the inputs are gzipped
	zcat $forward_dir/* >$Forward_fq
elif [ ! -f $Forward_fq ] && [ $gz == 'n' ]; then
	#There is not a single fastq file and the inputs are not gzipped
	cat $forward_dir/* >$Forward_fq
elif [ -f $Forward_fq ]; then
	echo "Feeding $Forward_fq to trimmomatic"
elif [ -f $Forward_fq.gz ] && [ $gz == 'y' ]; then
	#Case where the input is a single gzipped file and name is formatted correctly
	echo "Feeding $Forward_fq.gz to trimmomatic"
	gunzip -c $Forward_fq.gz | head -n 100000 >$Forward_fq
	offset=$( phred_offset_detector.py -f $Forward_fq )
	rm $Forward_fq
else
	echo "Input is unsupported"
	exit
fi

if [ ! -f $Reverse_fq.gz ] && [ $gz == 'y' ]; then
	#There is not a single gzipped file and the inputs are gzipped
	zcat $reverse_dir/* >$Reverse_fq
elif [ ! -f $Reverse_fq ] && [ $gz == 'n' ]; then
	#There is not a single fastq file and the inputs are not gzipped
	cat $reverse_dir/* >$Reverse_fq
elif [ -f $Reverse_fq ]; then
	echo "Feeding $Reverse_fq to trimmomatic"
elif [ -f $Reverse_fq.gz ]; then
	echo "Feeding $Reverse_fq.gz to trimmomatic"
else
	echo "Input is unsupported"
	exit
fi

############################################ Trimmomatic ############################################

if [ -z $offset ]; then
	offset=$( phred_offset_detector.py -f $Forward_fq )
	if [ "$offset" == '-' ]; then
		echo "ERROR: phred_offset_detector.py failed. It probably needs more reads."
        	exit
	fi
	echo "Phred offset detected to be $offset"
fi

STARTs=$(date +%s)

if [ ! -d $trim_dir ]; then
    echo "No trimmed files detected in $trim_dir. Proceeding with trimmomatic."
    mkdir $trim_dir
    mkdir $trim_dir/log
    java -jar /usr/local/bin/trimmomatic-0.35.jar PE \
-threads $threads \
$offset \
$Forward_fq* $Reverse_fq* \
$trim_dir/$sample\_pe.1.fq \
$trim_dir/$sample\_se.1.fq \
$trim_dir/$sample\_pe.2.fq \
$trim_dir/$sample\_se.2.fq \
ILLUMINACLIP:/usr/local/share/Trimmomatic-0.35/adapters/TruSeq3-PE.fa:2:3:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
1>$trim_dir/log/trim.stdout \
2>$trim_dir/log/trim.stderr

else
	echo "Using the files in $trim_dir:"
	ls -thor $trim_dir
fi

# Check to make sure the trimming process has not already completed
if [ ! -f $trim_dir/singletons.fq.gz ]; then

	printf "Counting the number of reads... "
	numreads=$(cat $trim_dir/*fq | wc -l | gawk '{ print $1/4 }')
	printf "done.\n"

	# printf "Compressing raw fastq files... "
	# if [ ! -f $Forward_fq.gz ]; then
	# 	gzip $Forward_fq
	# fi
	# if [ ! -f $Reverse_fq.gz ]; then
	# 	gzip $Reverse_fq
	# fi
	# printf "done.\n"

	########################################### fq2fa ##################################################
	cd $trim_dir
	if [ ! -f $trim_dir/singletons.fq ]; then
		cat $sample\_se.?.fq >singletons.fq
		rm $trim_dir/$sample\_se.*.fq
	fi
	
	##########################################################################################
	#fq2fa - Convert Fastq sequences to Fasta sequences.
	#Usage: fq2fa tmp.fq tmp.fa [...] 
	#       fq2fa --paired tmp.fq tmp.fa
	#       fq2fa --merge tmp_1.fq tmp_2.fq tmp.fa
	#Allowed Options: 
	#      --paired                           if the reads are paired-end in one file
	#      --merge                            if the reads are paired-end in two files
	#      --filter                           filter out reads containing 'N'
	###################################### run Megahit assembly ########################################
else
	gunzip $trim_dir/*.fq.gz
	numreads=$(cat $trim_dir/*fq | wc -l | gawk '{ print $1/4 }' )
fi

head -1000 $trim_dir/singletons.fq | fastq-to-fasta >tmp.fa
max_read=$( splitFASTA -i tmp.fa -o tmp.fa2 -v -m 50 | grep "Longest sequence" | gawk -F" = " '{ print $2 }' )
if [ "$max_read" -lt 81 ]; then
        k_max=77
elif [ "$max_read" -lt 117 ]; then
        k_max=97
elif [ "$max_read" -lt 150 ]; then
	k_max=117
else
	k_max=147
fi
rm tmp.fa tmp.fa2

cd $assembly_dir/

STARTs=$(date +%s)
STARTdate=$(date)
echo "Beginning Megahit assembly of $sample at: $STARTdate"

megahit \
-1 $trim_dir/$sample\_pe.1.fq \
-2 $trim_dir/$sample\_pe.2.fq \
--read $trim_dir/singletons.fq \
--k-min 27 \
--k-max $k_max \
--k-step 10 \
--memory 0.30 \
--min-contig-len 500 \
--out-dir $sample \
--num-cpu-threads $threads
# --presets meta-sensitive \

####################################################################################################
#Usage:
#  megahit [options] {-1 <pe1> -2 <pe2> | --12 <pe12> | -r <se>} [-o <out_dir>]
#
#  Input options that can be specified for multiple times (supporting plain text and gz/bz2 extensions)
#    -1                       <pe1>          comma-separated list of fasta/q paired-end #1 files, paired with files in <pe2>
#    -2                       <pe2>          comma-separated list of fasta/q paired-end #2 files, paired with files in <pe1>
#    --12                     <pe12>         comma-separated list of interleaved fasta/q paired-end files
#    -r/--read                <se>           comma-separated list of fasta/q single-end files
#
#  Input options that can be specified for at most ONE time (not recommended):
#    --input-cmd              <cmd>          command that outputs fasta/q reads to stdout; taken by MEGAHIT as SE reads 
#
#Optional Arguments:
#  Basic assembly options:
#    --min-count              <int>          minimum multiplicity for filtering (k_min+1)-mers, default 2
#    --k-min                  <int>          minimum kmer size (<= 127), must be odd number, default 21
#    --k-max                  <int>          maximum kmer size (<= 127), must be odd number, default 99
#    --k-step                 <int>          increment of kmer size of each iteration (<= 28), must be even number, default 10
#    --k-list                 <int,int,..>   comma-separated list of kmer size (all must be odd, in the range 15-127, increment <= 28);
#                                            override `--k-min', `--k-max' and `--k-step'
#
#  Advanced assembly options:
#    --no-mercy                              do not add mercy kmers
#    --no-bubble                             do not merge bubbles
#    --merge-level            <l,s>          merge complex bubbles of length <= l*kmer_size and similarity >= s, default 20,0.98
#    --prune-level            <int>          strength of local low depth pruning (0-2), default 2
#    --low-local-ratio        <float>        ratio threshold to define low local coverage contigs, default 0.2
#    --max-tip-len            <int>          remove tips less than this value; default 2*k for iteration of kmer_size=k
#    --no-local                              disable local assembly
#    --kmin-1pass                            use 1pass mode to build SdBG of k_min
#
#  Presets parameters:
#    --presets                <str>          override a group of parameters; possible values:
#                                            meta            '--min-count 2 --k-list 21,41,61,81,99'             (generic metagenomes, default)
#                                            meta-sensitive  '--min-count 2 --k-list 21,31,41,51,61,71,81,91,99' (more sensitive but slower)
#                                            meta-large      '--min-count 2 --k-list 27,37,47,57,67,77,87'       (large & complex metagenomes, like soil)
#                                            bulk            '--min-count 3 --k-list 31,51,71,91,99 --no-mercy'  (experimental, standard bulk sequencing with >= 30x depth)
#                                            single-cell     '--min-count 3 --k-list 21,33,55,77,99,121 --merge_level 20,0.96' (experimental, single cell data)
#
#  Hardware options:
#    -m/--memory              <float>        max memory in byte to be used in SdBG construction; default 0.9
#                                            (if set between 0-1, fraction of the machine's total memory)
#    --mem-flag               <int>          SdBG builder memory mode, default 1
#                                            0: minimum; 1: moderate; others: use all memory specified by '-m/--memory'.
#    --use-gpu                               use GPU
#    --gpu-mem                <float>        GPU memory in byte to be used. Default: auto detect to use up all free GPU memory. 
#    -t/--num-cpu-threads     <int>          number of CPU threads, at least 2. Default: auto detect to use all CPU threads.
#
#  Output options:
#    -o/--out-dir             <string>       output directory, default ./megahit_out
#    --out-prefix             <string>       output prefix (the contig file will be OUT_DIR/OUT_PREFIX.contigs.fa)
#    --min-contig-len         <int>          minimum length of contigs to output, default 200
#    --keep-tmp-files                        keep all temporary files
#
#Other Arguments:
#    --continue                              continue a MEGAHIT run from its last available check point.
#                                            please set the output directory correctly when using this option.
#    -h/--help                               print the usage message
#    -v/--version                            print version
#    --verbose                               verbose mode
################################ Clean up tmp files ################################################
ENDs=$(date +%s)
ENDdate=$(date)
echo "Assembly completed at: $ENDdate"
DIFF=$(($ENDs - $STARTs))
echo "Assembly of $sample required $DIFF seconds"
cd $assembly_dir/$sample
if [ -f final.contigs.fa ]; then
        mv final.contigs.fa $sample\_contig.fa
	mv log Megahit.log
	rm -r intermediate_contigs/ opts.txt
	if [ -f $sample\_contig.fa ]; then
		echo "Assembly completed successfully"
	fi
else
        echo "Something went wrong during Megahit assembly!"
fi
# cd $trim_dir
# ls singletons.fq $sample\_pe.*.fq >tmp
# parallel gzip <tmp
# rm tmp
echo "$sample contains $numreads reads"
