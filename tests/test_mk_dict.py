#!/usr/bin/env python

"""Tests of wrapper code in pycits: clean_up"""

import os
import shutil
from collections import defaultdict
import matplotlib.pyplot as plt
import unittest
from interogate.parse_gtf import parse_gff_gft
from interogate.return_dict import generate_transcript_coordinates
from interogate.tools import NotExecutableError

# Ensure the `query_transcript_exon` function is defined
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
                total_exons = len(transcript_dict[transcript_id])  # This gives the number of exons in the transcript
                return exon_number, total_exons
    return None, None

# get the data and the dict
file_path = 'data/test.gtf'  # Replace with the path to your GFF or GTF file
features = parse_gff_gft(file_path)
transcript_dict, transcript_exon_counts, gene_exon_counts, \
     last_exon_for_transcript, transcript_lengths = generate_transcript_coordinates(features)

class TestTranscriptDict(unittest.TestCase):

    def test_dict_1(self):
        """Test 1: AT1G01010.1 at pos 3 should be exon 1."""
        query_transcript = "AT1G01010.1"
        query_position = 3
        exon_number, total_exons = query_transcript_exon(transcript_dict, query_transcript, query_position)

        if exon_number is not None:
            print(f'Transcript {query_transcript} position {query_position} is in exon {exon_number}.')
            print(f'Transcript {query_transcript} has {total_exons} exons.')
        else:
            print(f'Position {query_position} is not found in transcript {query_transcript}.')
        self.assertEqual(exon_number, 1)
        self.assertEqual(total_exons, transcript_exon_counts[query_transcript])

    def test_dict_1b(self):
        """Test 1b: AT1G01010.1 at pos 3 should be exon 1."""
        query_transcript = "AT1G01010.1"
        query_position = 3
        exon_number, total_exons = query_transcript_exon(transcript_dict, query_transcript, query_position)
        self.assertEqual(exon_number, 1)
        self.assertEqual(total_exons, transcript_exon_counts[query_transcript])

    def test_dict_2(self):
        """Test 2: AT1G01020.4 at pos 426 should be exon 4. Should have 7 exons."""
        query_transcript = "AT1G01020.4"
        query_position = 426
        exon_number, total_exons = query_transcript_exon(transcript_dict, query_transcript, query_position)
        print(f"{query_transcript} should have {total_exons} exons and position {query_position} should be exon {exon_number}")
        if exon_number is not None:
            print(f'Transcript {query_transcript} position {query_position} is in exon {exon_number}.')
            print(f'Transcript {query_transcript} has {total_exons} exons.')
        else:
            print(f'Position {query_position} is not found in transcript {query_transcript}.')
        self.assertEqual(exon_number, 4)
        self.assertEqual(total_exons, transcript_exon_counts[query_transcript])

    def test_dict_3(self):
        """Test 3: AT1G01020.4 at pos 860 should be exon 6. Should have 7 exons. NEGATIVE CODING GENE. DIFF TRANSCRIPT"""
        query_transcript = "AT1G01020.4"
        query_position = 860
        exon_number, total_exons = query_transcript_exon(transcript_dict, query_transcript, query_position)
        print(f"{query_transcript} should have {total_exons} exons and position {query_position} should be exon {exon_number}")
        if exon_number is not None:
            print(f'Transcript {query_transcript} position {query_position} is in exon {exon_number}.')
            print(f'Transcript {query_transcript} has {total_exons} exons.')
        else:
            print(f'Position {query_position} is not found in transcript {query_transcript}.')
        self.assertEqual(exon_number, 6)
        self.assertEqual(total_exons, transcript_exon_counts[query_transcript])

    def test_dict_4(self):
        """Test 4: AT1G01020.4 at pos 1550 should BREAK. Should have 7 exons."""
        query_transcript = "AT1G01020.4"
        query_position = 1550
        exon_number, total_exons = query_transcript_exon(transcript_dict, query_transcript, query_position)
        print(f"{query_transcript} should have {total_exons} exons and position {query_position} should BREAK")
        if exon_number is not None:
            print(f'Transcript {query_transcript} position {query_position} is in exon {exon_number}.')
            print(f'Transcript {query_transcript} has {total_exons} exons.')
        else:
            print(f'Position {query_position} is not found in transcript {query_transcript}.')
        self.assertIsNone(exon_number)
        self.assertEqual(total_exons, transcript_exon_counts[query_transcript])

if __name__ == '__main__':
    unittest.main()
