#!/usr/bin/env python3
import os
from collections import defaultdict


def generate_transcript_coordinates(features):
    """
    Generate transcript coordinates with continuous nucleotide positions for exons.

    Parameters:
    features (list): A list of tuples, each containing the fields of a feature.

    Returns:
    dict: A nested dictionary mapping each transcript ID to a dictionary of exons,
          where each exon maps to a list of nucleotide positions.
    """
    transcript_dict = defaultdict(lambda: defaultdict(list))
    exon_nucleotide_counters = defaultdict(int)  # To count the nucleotide positions within exons
    exon_counters = defaultdict(int)  # To count the number of exons per transcript
    
    for feature in features:
        seqname, source, feature_type, start, end, score, strand, frame, attribute = feature
        
        # Process only exon features
        if feature_type == 'exon':
            # Extract transcript ID from the Parent attribute field
            attributes = attribute.split(';')
            transcript_id = None
            for attr in attributes:
                if 'Parent' in attr:
                    transcript_id = attr.split('=')[1].strip() if '=' in attr else attr.split()[1].strip().strip('"')
                    break  # Exit the loop once Parent is found
            
            # If transcript ID is found, add the coordinates to the dictionary
            if transcript_id:
                exon_counters[transcript_id] += 1  # Increment exon counter
                exon_number = exon_counters[transcript_id]  # Current exon number
                
                for pos in range(start, end + 1):
                    transcript_dict[transcript_id][exon_number].append(exon_nucleotide_counters[transcript_id] + 1)
                    exon_nucleotide_counters[transcript_id] += 1
    
    return transcript_dict


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

