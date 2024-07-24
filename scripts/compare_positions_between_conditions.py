#!/usr/bin/env python3
#
# you need to run  collect_positions_of_m6a.py first and then compare the ouputs wit hthis script

import argparse
import os
import csv
from collections import defaultdict

def get_args():
    parser = argparse.ArgumentParser(description="Compare exon modification data from multiple files", add_help=False)
    file_directory = os.path.realpath(__file__).split("compare_exons.py")[0]
    optional = parser.add_argument_group('optional arguments')

    optional.add_argument("--files", dest='files',
                          action="store",
                          nargs='+',  # This allows multiple arguments
                          required=True,
                          default=[os.path.join(file_directory, "data", "test1.tsv"), os.path.join(file_directory, "data", "test2.tsv")],
                          type=str,
                          help="List of input files to be parsed and compared e.g. --files file1.tsv file2.tsv")
 
    optional.add_argument("--output", dest='output',
                          action="store", default=None,
                          type=str,
                          help="Path to the output file (default: derived from input filenames)")
    
    return parser.parse_args()

def parse_file(file_path):
    """
    Parse exon data from a given file and return a dictionary.

    This function reads a tab-separated values (TSV) file containing exon data,
    and returns a dictionary with transcript IDs, exon numbers, and positions.

    Args:
        file_path (str): The path to the input file to be parsed.

    Returns:
        dict: A dictionary where each key is a transcript ID, and each value
              is another dictionary. The inner dictionary's keys are exon
              numbers, and the values are lists of positions.
    """
    data = defaultdict(lambda: defaultdict(list))
    
    with open(file_path, mode='r') as file:
        csv_reader = csv.DictReader(file, delimiter='\t')
        for row in csv_reader:
            transcript_id = row['transcript_id']
            exon_number = row['exon_number']
            positions = row['positions'].split(',')
            data[transcript_id][exon_number].extend(positions)
            print(f"Parsed {transcript_id} {exon_number} with positions: {positions}")  # Debug print
    
    return data

def compare_files(file_data, filenames):
    """
    Compare exon modification data from multiple files.

    This function takes a list of dictionaries containing exon modification data
    and compares the transcripts and positions of modifications.

    Args:
        file_data (list): A list of dictionaries where each dictionary contains
                          exon modification data from a file.
        filenames (list): A list of filenames corresponding to the data.

    Returns:
        dict: A dictionary containing the comparison results.
    """
    comparison_results = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    
    all_transcripts = set(file_data[0].keys())
    for data in file_data[1:]:
        all_transcripts.intersection_update(data.keys())
    
    for transcript_id in all_transcripts:
        exon_numbers = set(file_data[0][transcript_id].keys())
        for data in file_data[1:]:
            exon_numbers.intersection_update(data[transcript_id].keys())
        
        for exon_number in exon_numbers:
            positions_sets = [set(data[transcript_id][exon_number]) for data in file_data]
            print(f"Positions sets for {transcript_id} {exon_number}: {positions_sets}")  # Debug print
            common_positions = set.intersection(*positions_sets)
            unique_positions = [set.difference(pos_set, common_positions) for pos_set in positions_sets]
            
            comparison_results[transcript_id][exon_number]['common'] = sorted(list(common_positions), key=int)
            for idx, unique in enumerate(unique_positions):
                comparison_results[transcript_id][exon_number][f'unique_{os.path.basename(filenames[idx])}'] = sorted(list(unique), key=int)
                print(f"Unique positions for {transcript_id} {exon_number} in file {filenames[idx]}: {unique}")  # Debug print
    
    return comparison_results

def main():
    args = get_args()
    file_paths = args.files
    output_file = args.output
    
    if output_file is None:
        base_names = [os.path.splitext(os.path.basename(f))[0] for f in file_paths]
        output_file = "_vs_".join(base_names) + "_comparison.tsv"
    
    file_data = [parse_file(file_path) for file_path in file_paths]
    
    comparison_results = compare_files(file_data, file_paths)
    
    with open(output_file, 'w', newline='') as out_file:
        writer = csv.writer(out_file, delimiter='\t')
        header = ['transcript_id', 'exon_number', 'common_positions', 'num_common_positions']
        for filename in file_paths:
            basename = os.path.basename(filename)
            header.extend([f'unique_positions_{basename}', f'num_unique_positions_{basename}'])
        writer.writerow(header)
        
        for transcript_id, exons in comparison_results.items():
            for exon_number, positions_data in exons.items():
                common_positions = positions_data['common']
                num_common_positions = len(common_positions)
                row = [transcript_id, exon_number, ', '.join(common_positions), num_common_positions]
                
                for filename in file_paths:
                    unique_key = f'unique_{os.path.basename(filename)}'
                    unique_positions = positions_data.get(unique_key, [])
                    num_unique_positions = len(unique_positions)
                    row.extend([', '.join(unique_positions), num_unique_positions])
                
                writer.writerow(row)
    
    print(f"Comparison data successfully written to {output_file}")

if __name__ == "__main__":
    main()
