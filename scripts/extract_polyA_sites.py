import pysam
import re
import pandas as pd
import argparse
import os


def get_args():
    parser = argparse.ArgumentParser(description="Extract poly(A) sites from nanopore direct RNAseq data",
                                     add_help=False)
    file_directory = os.path.realpath(__file__).split("extract_polyA_sites.py")[0]                        
    optional = parser.add_argument_group('optional arguments')

    optional.add_argument("--bam", dest='bam',
                          action="store",
                          nargs='+',  # This allows multiple arguments
                          required=True,
                          type=str,
                          help="List of BAM files to be parsed e.g. --bam file1.bam file2.bam file3.bam")
    
    optional.add_argument("--output", dest='output',
                          action="store",
                          type=str,
                          required=True,
                          help="Output TSV file to store the poly(A) sites")

    optional.add_argument("--fasta", dest='fasta',
                          action="store",
                          type=str,
                          required=True,
                          help="FASTA file of the reference genome")

    return parser.parse_args()


def extract_polyA_sites(bam_file, fasta_file):
    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")
    # Open the FASTA file
    fasta = pysam.FastaFile(fasta_file)
    
    polyA_sites = []

    # Iterate over each read in the BAM file
    for read in bam.fetch():
        if not read.is_unmapped:
            # Get the sequence of the read
            seq = read.query_sequence
            
            # Search for a poly(A) tail at the end of the sequence
            match = re.search(r'(A{10,})$', seq)
            
            if match:
                # Get the read name
                read_name = read.query_name
                # Get the start position of the poly(A) tail
                polyA_start = read.reference_start + match.start()
                # Calculate the length of the poly(A) tail
                polyA_length = len(match.group(0))
                
                # Map to transcriptome
                transcript_id = read.reference_name
                
                # Get the genomic coordinate
                chrom = read.reference_name
                coordinate = read.reference_start + match.start()
                
                # Get the 20bp region before the poly(A) site
                if coordinate >= 20:
                    region_start = coordinate - 20
                    region_end = coordinate
                    pre_polyA_seq = fasta.fetch(chrom, region_start, region_end)
                else:
                    pre_polyA_seq = fasta.fetch(chrom, 0, coordinate)
                
                # Store the result
                polyA_sites.append([read_name, transcript_id, coordinate, polyA_start, polyA_length, pre_polyA_seq])
    
    return polyA_sites


def main():
    args = get_args()

    all_polyA_sites = []

    for bam_file in args.bam:
        polyA_sites = extract_polyA_sites(bam_file, args.fasta)
        all_polyA_sites.extend(polyA_sites)

    # Convert results to a DataFrame
    polyA_df = pd.DataFrame(all_polyA_sites, columns=['Read_Name', 'TranscriptID', 'Genomic_Coordinate', 'PolyA_Start', 'PolyA_Length', 'Pre_PolyA_Sequence'])
    
    # Save to TSV
    polyA_df.to_csv(args.output, sep='\t', index=False)
    print(f"Poly(A) sites have been extracted and saved to {args.output}")

if __name__ == "__main__":
    main()
