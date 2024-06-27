#!/usr/bin/env python3
import os
from collections import defaultdict
import re


def generate_transcript_coordinates(features):
    """
    Generate transcript coordinates with continuous nucleotide positions for exons.
    Also, mark the last exon for each transcript.

    Parameters:
    features (list): A list of tuples, each containing the fields of a feature.

    Returns:
    tuple: A nested dictionary mapping each transcript ID to a dictionary of exons,
           where each exon maps to a list of nucleotide positions,
           a dictionary mapping each transcript ID to the number of exons in the transcript,
           and a dictionary mapping each gene ID to the total number of unique exons,
           and a dictionary marking the last exon for each transcript.
    """
    transcript_dict = defaultdict(lambda: defaultdict(list))
    transcript_exon_sets = defaultdict(set)
    gene_exon_sets = defaultdict(set)
    nucleotide_counter = defaultdict(int)
    last_exon_for_transcript = {}
    transcript_lengths = defaultdict(int)

    for feature in features:
        seqname, source, feature_type, start, end, score, strand, frame, attribute = feature
        
        if feature_type in ['exon', 'three_prime_UTR']:
            transcript_id = None
            exon_number = None
            attributes = attribute.split(';')
            for attr in attributes:
                if 'Parent' in attr:
                    transcript_id = attr.split('=')[1].strip() if '=' in attr else attr.split()[1].strip().strip('"')
                if 'ID' in attr:
                    exon_match = re.search(r'exon:(\d+)', attr)
                    if exon_match:
                        exon_number = int(exon_match.group(1))

            if transcript_id and exon_number:
                gene_id = transcript_id.split('.')[0]
                exon_positions = []
                if strand == '+':
                    for pos in range(start, end + 1):
                        nucleotide_counter[transcript_id] += 1
                        exon_positions.append(nucleotide_counter[transcript_id])
                else:
                    for pos in range(start, end + 1):
                        nucleotide_counter[transcript_id] += 1
                        exon_positions.append(nucleotide_counter[transcript_id])
                
                transcript_dict[transcript_id][exon_number] = exon_positions
                transcript_exon_sets[transcript_id].add(exon_number)
                gene_exon_sets[gene_id].add(exon_number)

                if transcript_id not in last_exon_for_transcript or exon_number > last_exon_for_transcript[transcript_id]:
                    last_exon_for_transcript[transcript_id] = exon_number

    gene_exon_counts = {gene: len(exons) for gene, exons in gene_exon_sets.items()}
    trans_exon_counts = {trans: len(exons) for trans, exons in transcript_exon_sets.items()}
    last_exon_for_transcript = {trans: max(exons) for trans, exons in transcript_exon_sets.items() if exons}

    return transcript_dict, trans_exon_counts, gene_exon_counts, last_exon_for_transcript, transcript_lengths



def query_transcript_exon(transcript_dict, transcript_id, position):
    """
    Query the exon and total number of exons for a given transcript ID and coordinate.

    Parameters:
    transcript_dict (dict): A nested dictionary mapping each transcript ID to a dictionary of exons,
                            where each exon maps to a list of nucleotide positions.
    transcript_id (str): The ID of the transcript to query.
    position (int): The nucleotide position to query.

    Returns:
    tuple: The exon number that the coordinate belongs to and the total number of exons for the transcript.
    """
    if transcript_id in transcript_dict:
        for exon_number, coordinates in transcript_dict[transcript_id].items():
            if position in coordinates:
                total_exons = len(transcript_dict[transcript_id])
                return exon_number, total_exons
    return None, None