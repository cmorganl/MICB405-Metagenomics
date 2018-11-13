#!/bin/bash

# data_dir="/home/cmorganlang/Hallam_projects/Genome_assembly/Hydrocarbon/RawData/"

if [ $# -eq 0 ]; then
	echo "$0 sample_name working_dir forward.fastq reverse.fastq"
	exit
fi

sample=$1
working_dir=$2
forward=$3
reverse=$4

if [ ! -d $working_dir ]; then
	echo "$working_dir does not exist."
	exit
fi

# Example of a $line: IX0866/D1CP6ACXX_8/fastq/IX0866_D1CP6ACXX_8_AACCCC 
echo $working_dir/$sample
if [ ! -d $working_dir/$sample ]; then
	mkdir $working_dir/$sample
fi

#If reverse is not provided, assume reads are interleaved and split the interleaved FASTQ:
if [ -z $reverse ]; then
	cat $forward |paste - - - - - - - - | tee >(cut -f 1-4 | tr "\t" "\n" > $working_dir/$sample/$sample\_R1.fastq) | cut -f 5-8 | tr "\t" "\n" > $working_dir/$sample/$sample\_R2.fastq
	forward=$working_dir/$sample/$sample\_R1.fastq
	reverse=$working_dir/$sample/$sample\_R2.fastq
fi

phred_offset=$( python ~/bin/MeatLoaf/phred_offset_detector.py -f $reverse )
echo $phred_offset
po=""
if [ $phred_offset == "-phred33" ]; then
	po="--phred33"
fi

##################################### Run EMIRGE ###################################################

emirge.py $working_dir/$sample -1 $forward -2 $reverse \
-f /home/cmorganlang/bin/EMIRGE/SSURef_111_candidate_db.fasta \
-b /home/cmorganlang/bin/EMIRGE/SSU_candidate_db_btindex \
--max_read_length=151 \
--insert_mean=500 \
--insert_stddev=250 \
--processors=4 \
--iterations=20 \
--min_depth=2 \
$po

############## Rename the FASTA output and return high probability sequences  #####################
emirge_rename_fasta.py --record_prefix=$sample \
--prob_min=0.00001 \
$working_dir/$sample/iter.20 > $working_dir/$sample/$sample.emirge.fasta

############################################################ EMIRGE help ##########################
#Usage: emirge.py DIR <required_parameters> [options]
#
#This version of EMIRGE (emirge.py) attempts to reconstruct rRNA SSU genes from
#Illumina metagenomic data.
#DIR is the working directory to process data in.
#Use --help to see a list of required and optional arguments
#
#Additional information:
#https://groups.google.com/group/emirge-users
#https://github.com/csmiller/EMIRGE/wiki
#
#If you use EMIRGE in your work, please cite these manuscripts, as appropriate.
#
#Miller CS, Baker BJ, Thomas BC, Singer SW, Banfield JF (2011)
#EMIRGE: reconstruction of full-length ribosomal genes from microbial community short read sequencing data.
#Genome biology 12: R44. doi:10.1186/gb-2011-12-5-r44.
#
#Miller CS, Handley KM, Wrighton KC, Frischkorn KR, Thomas BC, Banfield JF (2013)
#Short-Read Assembly of Full-Length 16S Amplicons Reveals Bacterial Diversity in Subsurface Sediments.
#PloS one 8: e56018. doi:10.1371/journal.pone.0056018.
#
#
#Options:
#  -h, --help            show this help message and exit
#
#  Required flags:
#    These flags are all required to run EMIRGE, and may be supplied in any
#    order.
#
#    -1 reads_1.fastq[.gz]
#                        path to fastq file with \1 (forward) reads from
#                        paired-end sequencing run, or all reads from single-
#                        end sequencing run.  File may optionally be gzipped.
#                        EMIRGE expects ASCII-offset of 64 for quality scores.
#                        (Note that running EMIRGE with single-end reads is
#                        largely untested.  Please let me know how it works for
#                        you.)
#    -f FASTA_DB, --fasta_db=FASTA_DB
#                        path to fasta file of candidate SSU sequences
#    -b BOWTIE_DB, --bowtie_db=BOWTIE_DB
#                        precomputed bowtie index of candidate SSU sequences
#                        (path to appropriate prefix; see --fasta_db)
#    -l MAX_READ_LENGTH, --max_read_length=MAX_READ_LENGTH
#                        length of longest read in input data.
#
#  Required flags for paired-end reads:
#    These flags are required to run EMIRGE when you have paired-end reads
#    (the standard way of running EMIRGE), and may be supplied in any
#    order.
#
#    -2 reads_2.fastq    path to fastq file with \2 (reverse) reads from
#                        paired-end run.  File must be unzipped for mapper.
#                        EMIRGE expects ASCII-offset of 64 for quality scores.
#    -i INSERT_MEAN, --insert_mean=INSERT_MEAN
#                        insert size distribution mean.
#    -s INSERT_STDDEV, --insert_stddev=INSERT_STDDEV
#                        insert size distribution standard deviation.
#
#  Optional parameters:
#    Defaults should normally be fine for these options in order to run
#    EMIRGE
#
#    -n ITERATIONS, --iterations=ITERATIONS
#                        Number of iterations to perform.  It may be necessary
#                        to use more iterations for more complex samples
#                        (default=40)
#    -a PROCESSORS, --processors=PROCESSORS
#                        Number of processors to use in the mapping steps.  You
#                        probably want to raise this if you have the
#                        processors. (default: 1)
#    -m MAPPING, --mapping=MAPPING
#                        path to precomputed initial mapping (bam file).  If
#                        not provided, and initial mapping will be run for you.
#    -p SNP_FRACTION_THRESH, --snp_fraction_thresh=SNP_FRACTION_THRESH
#                        If fraction of variants in a candidate sequence
#                        exceeds this threhold, then split the candidate into
#                        two sequences for next iteration.  See also
#                        --variant_fraction_thresh. (default: 0.04)
#    -v VARIANT_FRACTION_THRESH, --variant_fraction_thresh=VARIANT_FRACTION_THRESH
#                        minimum probability of second most probable base at a
#                        site required in order to call site a variant.  See
#                        also --snp_fraction_thresh.  (default: 0.1)
#    -j JOIN_THRESHOLD, --join_threshold=JOIN_THRESHOLD
#                        If two candidate sequences share >= this fractional
#                        identity over their bases with mapped reads, then
#                        merge the two sequences into one for the next
#                        iteration.  (default: 0.97; valid range: [0.95, 1.0] )
#    -c MIN_DEPTH, --min_depth=MIN_DEPTH
#                        minimum average read depth below which a candidate
#                        sequence is discarded for next iteration(default: 3)
#    --nice_mapping=NICE_MAPPING
#                        If set, during mapping phase, the mapper will be
#                        "niced" by the Linux kernel with this value (default:
#                        no nice)
#    --phred33           Illumina quality values in fastq files are the (fastq
#                        standard) ascii offset of Phred+33.  This is the new
#                        default for Illumina pipeline >= 1.8. DEFAULT is still
#                        to assume that quality scores are Phred+64
#    -e SAVE_EVERY, --save_every=SAVE_EVERY
#                        every SAVE_EVERY iterations, save some information
#                        about the program's state.  This is solely for
#                        debugging information, and is NOT required to resume a
#                        run (see --resume_from below).  (default=none)
#
#  Resuming iterations:
#    These options allow you to resume iterations from a previously
#    completed EMIRGE iteration.  This requires that directories for the
#    iteration to resume from and the previous iteration both be present.
#    It is STRONGLY recommended that other options set on the command line
#    be identical to the original run.  Note that EMIRGE does not check
#    this for you!
#
#    -r RESUME_FROM, --resume_from=RESUME_FROM
#                        Resume iterations from COMPLETED iteration specified.
#                        Requires that the iteration and previous iteration
#                        fully completed, i.e. a priors file, bam file, and
#                        fasta file are all present in the iteration directory.
