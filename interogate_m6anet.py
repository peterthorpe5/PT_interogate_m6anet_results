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
from interogate.return_dict import generate_transcript_coordinates, query_transcript_exon
from interogate.parse_m6a_site_proba import identify_methylated_sites
from interogate.plot import plot_methylation_distribution
from interogate.summary_stats import summarise_methylation_sites
from interogate.parse_trans_len import parse_transcript_lengths
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
                          help="(data.site_proba.csv) List of m6anet result files to be parsed e.g. --m6a file1.csv file2.csv file3.csv")
 
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
    
    optional.add_argument("--len", dest='trans_len',
                          action="store",
                          default=os.path.join(file_directory, "data", 
                                               "Araport11_genes.201606.cdna.len"),
                          type=str,
                          help="input gtf file to get the transcript coordinates")
    
    optional.add_argument("-l", "--logfile", dest='logfile',
                          action="store",
                          default="pipeline.log",
                          type=str,
                          help="log file name")
    optional.add_argument("--test", dest='test',
                          action="store",
                          default=False,
                          type=str,
                          help="extra printing for testing, add true if required")
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

    # get the transcript lenghts
    transcript_lengths = parse_transcript_lengths(args.trans_len)
    logger.info("processed: %s", args.trans_len)

    # Example usage
    logger.info("Starting processing: %s", args.gtf )
    file_path = args.gtf  # Replace with the path to your GFF or GTF file

    # Parse GTF file and generate transcript coordinates
    file_path = args.gtf
    features = parse_gff_gft(file_path)
    transcript_dict, transcript_exon_counts, gene_exon_counts, last_exon_for_transcript, \
            transcript_strands = generate_transcript_coordinates(features, transcript_lengths)
    

    if args.test:
        with open(args.out, 'w') as out_file:
            for transcript, exons in transcript_dict.items():
                for exon, coordinates in exons.items():
                    out_data = f'{transcript} exon {exon}: {coordinates}'
                    out_file.write(out_data + '\n')
                    print(out_data)

   # Process each m6A result file
    
    for m6a_file in args.m6a:
        # These were used in testing - useful. but not for real data. 
        # print("Transcript Dictionary:", transcript_dict)
        # print("Exon Counts per Transcript:", transcript_exon_counts)
        # print("Gene Exon Counts:", gene_exon_counts)
        ## print("Last Exon for Each Transcript:", last_exon_for_transcript)
        # print("Transcript Lengths:", transcript_lengths)

        try:
            logger.info("Starting processing: %s", m6a_file)
            threshold = args.threshold
            # Identify methylated sites
            methylated_sites = identify_methylated_sites(m6a_file, threshold)

            # Filter out invalid transcripts
            valid_transcripts = set(transcript_dict.keys())
            methylated_sites = methylated_sites[methylated_sites['transcript_id'].isin(valid_transcripts)]
            if methylated_sites.empty:
                logger.warning(f"No valid methylated sites after filtering for file: {m6a_file}")
                continue
        
            # Annotate the sites
            results = []
            for index, row in methylated_sites.iterrows():
                transcript_id = row['transcript_id']
                position = row['transcript_position']
                exon_number, total_exons = query_transcript_exon(transcript_dict, transcript_id, position)
                is_last_exon = (exon_number == last_exon_for_transcript.get(transcript_id, None))
                if exon_number is not None:
                    result = {
                        'transcript_id': transcript_id,
                        'position': position,
                        'exon_number': exon_number,
                        'total_exons_in_transcript': total_exons,
                        'total_exons_in_gene': gene_exon_counts.get(transcript_id.split('.')[0], 'Unknown'),
                        'is_last_exon': is_last_exon
                    }
                else:
                    result = {
                        'transcript_id': transcript_id,
                        'position': position,
                        'exon_number': 'UTR',
                        'total_exons_in_transcript': total_exons,
                        'total_exons_in_gene': gene_exon_counts.get(transcript_id.split('.')[0], 'Unknown'),
                        'is_last_exon': False
                    }
                results.append(result)

            results_df = pd.DataFrame(results)
            print("Results DataFrame:", results_df)
            logger.info("Results DataFrame: ")
            logger.info(results_df)

            output_file = f"{os.path.splitext(m6a_file)[0]}_exon_annotated.tab"
            results_df.to_csv(output_file, index=False, sep="\t")
            print(f"Results saved to {output_file}")
            output_plot = f"{os.path.splitext(m6a_file)[0]}_m6a_distribution.pdf"

            try:
                logger.info(f"Plot saved to {output_plot}")
                plot_methylation_distribution(results_df, output_plot, 
                                              transcript_lengths, transcript_strands,
                                              m6a_file, logger)
            except Exception as e:
                logger.error(f"An error occurred while plotting the methylation distribution: {e}")
            # Continue with the rest of the script

            # write out a summary per transcript usage
            output_summary = f"{os.path.splitext(m6a_file)[0]}_summary_per_transcript.tab"
            summarise_methylation_sites(results_df, output_summary, 
                                        logger)
        

            logger.info("Processing finished: %s", time.asctime())
            logger.info("########################\n")
        except Exception as m6a_file_e:
            logger.error(f"An error occurred while processing the file {m6a_file}: {m6a_file_e}")
            continue  # Skip to the next file

if __name__ == '__main__':
    main()


# a spiked real UTR value AT1G01010.1  position > 1290 .. 
# seem to only return UTR at AT1G01010.1,2433,2433,0.999,AAACA,0.1111111111111111
#    this should be UTR as len transcript is  1290 ???
#    6 exons in total. UTR genomic coord 
info = """
1       Araport11       mRNA    3631    5899    .       +       .       ID=AT1G01010.1;Parent=AT1G01010;Name=AT1G01010.1;Note=NAC domain containing protein 1;conf_class=2;symbol=NAC001;full_name=NAC domain containing protein 1;computational_description=NAC domain containing protein 1;conf_rating=****;gene=2200934,UniProt=Q0WV96;curator_summary=Member of the NAC domain containing family of plant specific transcriptional regulators.
1       Araport11       CDS     3760    3913    .       +       0       ID=AT1G01010:CDS:1;Parent=AT1G01010.1;Name=NAC001:CDS:1;Note=NAC domain containing protein 1;curator_summary=Member of the NAC domain containing family of plant specific transcriptional regulators.;computational_description=NAC domain containing protein 1
1       Araport11       CDS     3996    4276    .       +       2       ID=AT1G01010:CDS:2;Parent=AT1G01010.1;Name=NAC001:CDS:2;Note=NAC domain containing protein 1;curator_summary=Member of the NAC domain containing family of plant specific transcriptional regulators.;computational_description=NAC domain containing protein 1
1       Araport11       CDS     4486    4605    .       +       0       ID=AT1G01010:CDS:3;Parent=AT1G01010.1;Name=NAC001:CDS:3;Note=NAC domain containing protein 1;curator_summary=Member of the NAC domain containing family of plant specific transcriptional regulators.;computational_description=NAC domain containing protein 1
1       Araport11       CDS     4706    5095    .       +       0       ID=AT1G01010:CDS:4;Parent=AT1G01010.1;Name=NAC001:CDS:4;Note=NAC domain containing protein 1;curator_summary=Member of the NAC domain containing family of plant specific transcriptional regulators.;computational_description=NAC domain containing protein 1
1       Araport11       CDS     5174    5326    .       +       0       ID=AT1G01010:CDS:5;Parent=AT1G01010.1;Name=NAC001:CDS:5;Note=NAC domain containing protein 1;curator_summary=Member of the NAC domain containing family of plant specific transcriptional regulators.;computational_description=NAC domain containing protein 1
1       Araport11       CDS     5439    5630    .       +       0       ID=AT1G01010:CDS:6;Parent=AT1G01010.1;Name=NAC001:CDS:6;Note=NAC domain containing protein 1;curator_summary=Member of the NAC domain containing family of plant specific transcriptional regulators.;computational_description=NAC domain containing protein 1
1       Araport11       exon    3631    3913    .       +       .       ID=AT1G01010:exon:1;Parent=AT1G01010.1;Name=AT1G01010:exon:1;Note=NAC domain containing protein 1;curator_summary=Member of the NAC domain containing family of plant specific transcriptional regulators.;computational_description=NAC domain containing protein 1
1       Araport11       exon    3996    4276    .       +       .       ID=AT1G01010:exon:2;Parent=AT1G01010.1;Name=AT1G01010:exon:2;Note=NAC domain containing protein 1;curator_summary=Member of the NAC domain containing family of plant specific transcriptional regulators.;computational_description=NAC domain containing protein 1
1       Araport11       exon    4486    4605    .       +       .       ID=AT1G01010:exon:3;Parent=AT1G01010.1;Name=AT1G01010:exon:3;Note=NAC domain containing protein 1;curator_summary=Member of the NAC domain containing family of plant specific transcriptional regulators.;computational_description=NAC domain containing protein 1
1       Araport11       exon    4706    5095    .       +       .       ID=AT1G01010:exon:4;Parent=AT1G01010.1;Name=AT1G01010:exon:4;Note=NAC domain containing protein 1;curator_summary=Member of the NAC domain containing family of plant specific transcriptional regulators.;computational_description=NAC domain containing protein 1
1       Araport11       exon    5174    5326    .       +       .       ID=AT1G01010:exon:5;Parent=AT1G01010.1;Name=AT1G01010:exon:5;Note=NAC domain containing protein 1;curator_summary=Member of the NAC domain containing family of plant specific transcriptional regulators.;computational_description=NAC domain containing protein 1
1       Araport11       exon    5439    5899    .       +       .       ID=AT1G01010:exon:6;Parent=AT1G01010.1;Name=AT1G01010:exon:6;Note=NAC domain containing protein 1;curator_summary=Member of the NAC domain containing family of plant specific transcriptional regulators.;computational_description=NAC domain containing protein 1
1       Araport11       five_prime_UTR  3631    3759    .       +       .       ID=AT1G01010:five_prime_UTR:1;Parent=AT1G01010.1;Name=NAC001:five_prime_UTR:1;Note=NAC domain containing protein 1;curator_summary=Member of the NAC domain containing family of plant specific transcriptional regulators.;computational_description=NAC domain containing protein 1
1
"""
