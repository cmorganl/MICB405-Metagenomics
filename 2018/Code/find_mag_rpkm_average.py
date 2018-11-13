#!/usr/bin/env python

import argparse
import os
import sys

__author__ = 'Connor Morgan-Lang'


def get_options():
    """
    Function that uses argparse to parse the command-line arguments and return an object with the parameters
    :return: args - an object containing parameters parsed from the command-line
    """
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--mag_fasta", required=False,
                        help="Path to a single FASTA file for a MAG/bin or genome")
    parser.add_argument("-l", "--list", required=False,
                        help="A file listing the FASTA files to be analyzed")
    parser.add_argument("-r", "--rpkm", required=True,
                        help="The csv output file from `rpkm`")
    parser.add_argument("-o", "--output", required=False, default="MAG_RPKM_averages.csv",
                        help="Name of file to write the output.")
    parser.add_argument("-v", "--verbose", required=False, default=False, action="store_true",
                        help="Be verbose in the progress printed to stdout")
    args = parser.parse_args()

    if not args.mag_fasta and not args.list:
        sys.stderr.write("ERROR: Either the `--mag_fasta` or `--list` must be used to provide ")
    return args


def read_rpkm_csv(rpkm_file):
    rpkm_dict = dict()
    rpkm_data = open(rpkm_file, 'r')
    record = rpkm_data.readline()
    while record:
        contig, rpkm = record.split(',')
        rpkm_dict[contig] = float(rpkm)
        record = rpkm_data.readline()

    rpkm_data.close()

    return rpkm_dict


def find_mag_headers(fasta):
    headers = list()
    fasta_records = open(fasta, 'r')
    record = fasta_records.readline()
    while record:
        if record[0] == '>':
            # This is a header line
            headers.append(record.strip()[1:])
        else:
            pass
        record = fasta_records.readline()
    fasta_records.close()
    return headers


def get_averages(mag_dict, rpkm_dict):
    mag_wise_rpkm = dict()
    all_rpkm_headers = list(rpkm_dict.keys())
    # Welcome to my very inefficient for loop
    for mag in mag_dict:
        rpkm_sum = 0
        headers = mag_dict[mag]
        for header in headers:
            if header in all_rpkm_headers:
                rpkm_sum += rpkm_dict[header]
            else:
                sys.stderr.write("ERROR: Unable to find header " + header + " in rpkm file!")
                sys.exit(2)
        mag_wise_rpkm[mag] = rpkm_sum/len(headers)
    return mag_wise_rpkm


def write_averages(mag_wise_rpkm, output_csv):
    rpkm_averages_csv = open(output_csv, 'w')
    header = "MAG,RPKM\n"
    csv_text = header
    for mag in mag_wise_rpkm:
        csv_text += mag + ',' + str(mag_wise_rpkm[mag]) + '\n'
    rpkm_averages_csv.write(csv_text)

    rpkm_averages_csv.close()

    return


def main():
    args = get_options()
    rpkm_dict = read_rpkm_csv(args.rpkm)
    mag_dict = dict()
    if args.list:
        fasta_list = open(args.list, 'r')
        for fasta in fasta_list:
            mag_dict[os.path.basename(fasta.strip())] = find_mag_headers(fasta.strip())
        fasta_list.close()
    else:
        mag_dict[os.path.basename(args.mag_fasta)] = find_mag_headers(args.mag_fasta)
    mag_wise_rpkm = get_averages(mag_dict, rpkm_dict)
    write_averages(mag_wise_rpkm, args.output)


main()
