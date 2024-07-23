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
 
    optional.add_argument("--output", dest='output',
                          action="store", default="m6a_positions.tsv",
                          type=str,
                          help="Path to the output file (default: output.tsv)")
    
    optional.add_argument("--thread", dest='threads',
                          action="store", default="1",
                          type=str,
                          help="number of threads: currently does nothing yet")
    
    return parser.parse_args()


def parse_file(file_path, data):
    """
    Parse exon data from a given file and update the provided dictionary.

    This function reads a tab-separated values (TSV) file containing exon data,
    and updates the given defaultdict with transcript IDs, exon numbers, and positions.

    Args:
        file_path (str): The path to the input file to be parsed.
        data (defaultdict): A defaultdict where each key is a transcript ID,
                            and each value is another defaultdict. The inner
                            defaultdict's keys are exon numbers, and the values
                            are lists of positions.

    Returns:
        None: The function updates the provided defaultdict in place.
    """
    with open(file_path, mode='r') as file:
        csv_reader = csv.DictReader(file, delimiter='\t')
        for row in csv_reader:
            transcript_id = row['transcript_id']
            exon_number = row['exon_number']
            position = row['position']
            data[transcript_id][exon_number].append(position)


def main():
    args = get_args()
    file_paths = args.file
    output_file = args.output
    
    all_data = defaultdict(lambda: defaultdict(list))
    
    for file_path in file_paths:
        parse_file(file_path, all_data)
    
    with open(output_file, 'w', newline='') as out_file:
        writer = csv.writer(out_file, delimiter='\t')
        writer.writerow(['transcript_id', 'exon_number', 'positions', 'num_positions'])
        
        for transcript_id, exons in all_data.items():
            for exon_number, positions in exons.items():
                num_positions = len(positions)
                positions_str = ','.join(positions)
                writer.writerow([transcript_id, exon_number, positions_str, num_positions])
    
    print(f"Data successfully written to {output_file}")

if __name__ == "__main__":
    main()
