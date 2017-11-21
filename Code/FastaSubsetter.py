#!/usr/bin/env python

import sys
import re
import argparse

__author__ = "Connor Morgan-Lang"


def set_arguments():
    parser = argparse.ArgumentParser()
    extract = parser.add_argument_group("extract")
    parser.add_argument("-i", "--fasta", help="The fasta file to be subsetted.", required=True)
    parser.add_argument("-l", "--list",
                        help="The list of contig headers to be subsetted from the fasta file", required=True)
    parser.add_argument("-o", "--output",
                        help="The output fasta file [DEFAULT = subset.fasta]", required=False, default="subset.fasta")
    parser.add_argument("-N", "--Nsplit", help="Split the contigs on ambiguity character 'N's [DEFAULT = False]",
                        default=False, action="store_true")
    parser.add_argument("-m", "--min_length", help="The minimum contig length [DEFAULT = 200].",
                        required=False,
                        default=200,
                        type=int)
    parser.add_argument("-c", "--min_coverage", default=0, required=False, type=int,
                        help="If the contig/scaffold coverage is included in the header name (by *_cov_X.Y_*),"
                             "then filter the contigs by their coverage [ DEFAULT = 0, no filtering ]")
    parser.add_argument("-x", "--longest",
                        help="Take the 'x' longest contigs. Writes all contigs by default.",
                        default=0, type=int)
    parser.add_argument("--inverse",
                        help="Write all contigs that are NOT in the list to output [DEFAULT = False]",
                        action="store_true",
                        default=False)
    parser.add_argument("--partition",
                        help="Specify the number of partitions to divide the input FASTA file.",
                        required=False, type=int, default=1)
    parser.add_argument("-v", "--verbose", required=False, action="store_true",
                        help="Verbosity flag to control runtime messages")
    extract.add_argument("-s", "--slice",
                         help="Flag indicating a specific sequence should be extracted from a specific header.",
                         required=False, action="store_true")
    extract.add_argument("-p", "--pos", help="Positions (start,end) of the provided contig to extract", required=False)
    args = parser.parse_args()
    return args


def subset_fasta(fasta, headers, invert):
    """
    Finds each header in `headers` within the `fasta` argument.
    NOTE: headers within the fasta file are split on spaces!
    :return: A dictionary with the contig headers of interest as keys and their respective sequences as their values.
    """
    subset = dict()
    header = ""
    add = 0
    with open(fasta) as fas:
        line = fas.readline()
        while line:
            if line[0] == '>':
                # header = line[1:].strip().split(' ')[0]
                header = line[1:].strip()
                if headers[0] == "all":
                    headers.append(header)
                if header in headers and invert is False:
                    subset[header] = ""
                    add = 1
                elif header not in headers and invert is True:
                    subset[header] = ""
                    add = 1
                else:
                    add = 0
            elif add == 1:
                subset[header] += line.strip()
            else:
                pass
            line = fas.readline()
    return subset, len(subset.keys())


def load_list(header_list):
    headers = list()
    try:
        LoH = open(header_list)
    except IOError:
        sys.exit("ERROR: Unable to find " + header_list + ". Exiting now!")
    for line in LoH:
        if line[0] == '>':
            headers.append(line[1:].strip())
        else:
            headers.append(line.strip())
    LoH.close()
    return headers


def extract_seq(contig_seq, start, end):
    """
    :return: A sub-string of the contig sequence as determined by the pos tuple.
    """
    start_pos = int(start)
    end_pos = int(end)
    return str(contig_seq[start_pos:end_pos])


def Nsplit_scaffolds(subset, min_length):
    """
    Splits all scaffolds containing 'N's on 'N' and
    adds these contigs to the dictionary to be returned.
    Removes all sub-contigs that are smaller than min_length.
    """
    unambiguous_subset = dict()
    acc = 1
    for head in subset.keys():
        if 'N' in subset[head]:
            contigs = subset[head].split('N')
            for c in contigs:
                if len(c) > min_length:
                    contig_name = head + "_" + str(acc)
                    acc += 1
                    unambiguous_subset[contig_name] = c
        else:
            contig_name = head + "_" + str(acc)
            acc += 1
            unambiguous_subset[contig_name] = subset[head]
        acc = 1
    return unambiguous_subset


def slice_contigs(subset, slice, start, end):
    contigs = dict()
    for contig in subset:
        seq = extract_seq(subset[contig], start, end)
        contigs[contig] = seq
    return contigs


def parse_contig_coverage(headers, min_coverage):
    spades_coverage_re = re.compile(r'_cov_')
    good_coverage_seqs = list()
    for header in headers:
        if spades_coverage_re.search(header):
            coverage = spades_coverage_re.split(header)[1]
            cov_int = float(coverage.split('_')[0])
            if cov_int >= min_coverage:
                good_coverage_seqs.append(header)
    return good_coverage_seqs


def write_subset_to_output(subset, output, min_length, longest, partition, min_coverage):
    total_seqs_written = 1
    fasta_length = 0
    if partition > 1:
        headers = subset.keys()
        print min_coverage
        if min_coverage > 0:
            headers = parse_contig_coverage(headers, min_coverage)
        num_subset_seqs = len(headers)
        if partition > num_subset_seqs:
            sys.exit("ERROR: Provided number of partitions is greater than number of sequences in input FASTA!")
        num_seqs_per = num_subset_seqs / partition
        seqs_parsed = 0
        while partition >= 1:
            output_file = output + "_" + str(partition)
            print output_file,
            seqs_written = 0
            with open(output_file, 'w') as fa_out:
                # TODO: Ensure the last sequence is being written
                # TODO: Handle the end of file case
                while seqs_written <= num_seqs_per and seqs_parsed < num_subset_seqs:
                    contig = headers[seqs_parsed]
                    if len(subset[contig]) > min_length:
                        fasta_length += len(subset[contig])
                        fa_out.write('>' + str(contig) + "\n")
                        fa_out.write(subset[contig] + "\n")
                        total_seqs_written += 1
                        seqs_written += 1
                    seqs_parsed += 1
            partition -= 1
    else:
        print output
        headers = sorted(subset, key=lambda contig: len(subset[contig]), reverse=True)
        if min_coverage > 0:
            headers = parse_contig_coverage(headers, min_coverage)
        with open(output, 'w') as fa_out:
            for contig in headers:
                if len(subset[contig]) > min_length and total_seqs_written <= longest:
                    fasta_length += len(subset[contig])
                    fa_out.write('>' + str(contig) + "\n")
                    fa_out.write(subset[contig] + "\n")
                    total_seqs_written += 1
    return fasta_length, total_seqs_written


def main():
    print "Beginning", sys.argv[0]
    args = set_arguments()
    headers = ["all"]
    if args.list != "all":
        headers = load_list(args.list)
    print "Finding the sequences in", args.list
    subset, num_contigs = subset_fasta(args.fasta, headers, args.inverse)
    if args.Nsplit:
        print "Splitting the sequences on 'N's..."
        subset = Nsplit_scaffolds(subset, args.min_length)
    if args.slice:
        start, end = args.pos.split(',')
        subset = slice_contigs(subset, args.slice, start, end)
    print "Writing the contig subset to",
    if args.longest == 0:
        args.longest = num_contigs
    fasta_length, n = write_subset_to_output(subset, args.output, args.min_length, args.longest,
                                             args.partition, args.min_coverage)
    if args.verbose:
        print "Genome length =", fasta_length
    return 0


main()
