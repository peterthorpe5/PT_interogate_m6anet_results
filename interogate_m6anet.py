#!/usr/bin/env python3
#
# interogate_m6anet.py 

import os
import sys
import shutil
import errno
import time
import argparse
from collections import defaultdict
import logging
import logging.handlers
import argparse
import matplotlib.pyplot as plt
import pandas as pd
from interogate.parse_gtf import parse_gff_gft
from interogate.return_dict import generate_transcript_coordinates
from interogate.parse_m6a_site_proba import identify_methylated_sites, query_transcript_exon
from interogate.plot import plot_methylation_distribution
from interogate.summary_stats import summarize_methylation_sites

from scipy.stats import chi2_contingency


print(" ...   libs loaded ...")

def get_args():
    parser = argparse.ArgumentParser(description="m6anet interogater:  " +
                                     "data for methylation ",
                                     add_help=False)
    file_directory = os.path.realpath(__file__).split("interogate_m6anet.py")[0]                        
    optional = parser.add_argument_group('optional arguments')

    optional.add_argument("--m6a", dest='m6a',
                          action="store",
                          nargs='+',  # This allows multiple arguments
                          required=True,
                          default=os.path.join(file_directory, "data", 
                                               "test.site_proba.csv"),
                          type=str,
                          help="List of m6anet result files to be parsed e.g. --m6a file1.csv file2.csv file3.csv")
 
    optional.add_argument("--thread", dest='threads',
                          action="store", default="1",
                          type=str,
                          help="number of threads: currently does nothing yet")
    

    optional.add_argument("--threshold", dest='threshold',
                          action="store", default=0.9,
                          type=float,
                          help="theshold for m6a dat filtering. Default is recommended 0.9")
    
    optional.add_argument("-o", "--out", dest='out',
                          action="store",
                          default="gene_exon_locations.out",
                          type=str,
                          help="outfile name")
        
    optional.add_argument("--gtf", dest='gtf',
                          action="store",
                          default=os.path.join(file_directory, "data", 
                                               "test.gtf"),
                          type=str,
                          help="input gtf file to get the transcript coordinates")
    
    optional.add_argument("-l", "--logfile", dest='logfile',
                          action="store",
                          default="pipeline.log",
                          type=str,
                          help="log file name")
    return parser.parse_args()
    

def main():
    args = get_args()

    # Set up logging
    logger = logging.getLogger('interogate_m6anet')
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    logger.addHandler(err_handler)

    try:
        logstream = open(args.logfile, 'w')
        err_handler_file = logging.StreamHandler(logstream)
        err_handler_file.setFormatter(err_formatter)
        # logfile is always verbose
        err_handler_file.setLevel(logging.INFO)
        logger.addHandler(err_handler_file)
    except:
        logger.error(f"Could not open {args.logfile} for logging")
        sys.exit(1)

    # Report input arguments
    logger.info(sys.version_info)
    logger.info("Command-line: %s", ' '.join(sys.argv))
    logger.info("Starting processing: %s", time.asctime())

    # Example usage
    logger.info("Starting processing: %s", args.gtf )
    file_path = args.gtf  # Replace with the path to your GFF or GTF file

    # Parse GTF file and generate transcript coordinates
    file_path = args.gtf
    features = parse_gff_gft(file_path)
    transcript_dict, transcript_exon_counts, gene_exon_counts, \
         last_exon_for_transcript = generate_transcript_coordinates(features)
    


    #with open(args.out, 'w') as out_file:
    #    for transcript, exons in transcript_dict.items():
    #        for exon, coordinates in exons.items():
    #            out_data = f'{transcript} exon {exon}: {coordinates}'
    #            out_file.write(out_data + '\n')
    #            print(out_data)

   # Process each m6A result file
    
    for m6a_file in args.m6a:
        logger.info("Starting processing: %s", m6a_file)
        threshold = args.threshold

        methylated_sites = identify_methylated_sites(m6a_file, threshold)
        # print(methylated_sites)

        # Determine exon/UTR location for each methylation site
        results = []
        for index, row in methylated_sites.iterrows():
            transcript_id = row['transcript_id']
            position = row['transcript_position']
            exon_number, total_exons_in_transcript = query_transcript_exon(transcript_dict,
                                                                            transcript_id, 
                                                                            position)

            if exon_number is not None:
                gene_id = transcript_id.split('.')[0]
                total_exons_in_gene = gene_exon_counts.get(gene_id, 'Unknown')
                is_last_exon = exon_number == last_exon_for_transcript.get(transcript_id)
                result = {
                    'transcript_id': transcript_id,
                    'position': position,
                    'exon_number': exon_number,
                    'total_exons_in_transcript': total_exons_in_transcript,
                    'total_exons_in_gene': total_exons_in_gene,
                    'is_last_exon': is_last_exon
                }
            else:
                result = {
                    'transcript_id': transcript_id,
                    'position': position,
                    'exon_number': 'UTR',
                    'total_exons_in_transcript': total_exons_in_transcript,
                    'total_exons_in_gene': 'Unknown',
                    'is_last_exon': False
                }
            results.append(result)

        results_df = pd.DataFrame(results)
        # print(results_df)

        # Print and save the result
        output_file = f"{os.path.splitext(m6a_file)[0]}_exon_annotated.tab"
        results_df.to_csv(output_file, index=False, sep="\t")
        print(f"Results saved to {output_file}")


        # plot out the data usage
        output_plot = f"{os.path.splitext(m6a_file)[0]}_m6a_distribution.pdf"
        plot_methylation_distribution(results_df, output_plot)

        # write out a summary per transcript usage
        output_summary = f"{os.path.splitext(m6a_file)[0]}_summary_per_transcript.tab"
        summarize_methylation_sites(results_df, output_summary)


    logger.info("Processing finished: %s", time.asctime())

if __name__ == '__main__':
    main()


