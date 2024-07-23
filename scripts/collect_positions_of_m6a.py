#!/usr/bin/env python3
#
# collect positions

# script to parse 

import argparse
import os
import csv
from collections import defaultdict

def get_args():
    parser = argparse.ArgumentParser(description="Parse exon data from file", add_help=False)
    file_directory = os.path.realpath(__file__).split("parse_exons.py")[0]
    optional = parser.add_argument_group('optional arguments')

    optional.add_argument("--file", dest='file',
                          action="store",
                          nargs='+',  # This allows multiple arguments
                          required=True,
                          default=os.path.join(file_directory, "data", "test.site_proba.csv"),
                          type=str,
                          help="List of input files to be parsed e.g. --file file1.tsv file2.tsv file3.tsv")
 
    optional.add_argument("--thread", dest='threads',
                          action="store", default="1",
                          type=str,
                          help="number of threads: currently does nothing yet")
    
    return parser.parse_args()

def parse_file(file_path):
    data = defaultdict(lambda: defaultdict(list))
    
    with open(file_path, mode='r') as file:
        csv_reader = csv.DictReader(file, delimiter='\t')
        for row in csv_reader:
            transcript_id = row['transcript_id']
            exon_number = row['exon_number']
            position = row['position']
            data[transcript_id][exon_number].append(position)
    
    return data

def main():
    args = get_args()
    file_paths = args.file
    
    all_data = defaultdict(lambda: defaultdict(list))
    
    for file_path in file_paths:
        parsed_data = parse_file(file_path)
        for transcript_id, exons in parsed_data.items():
            for exon_number, positions in exons.items():
                all_data[transcript_id][exon_number].extend(positions)
    
    for transcript_id, exons in all_data.items():
        print(f"{transcript_id}: {dict(exons)}")

if __name__ == "__main__":
    main()
