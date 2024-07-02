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
from interogate.parse_trans_len import parse_transcript_lengths


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
                total_exons = len(transcript_dict[transcript_id])
                return exon_number, total_exons
    return "UTR", len(transcript_dict[transcript_id])

# get the data and the dict
file_path = 'data/test.gtf'  # Replace with the path to your GFF or GTF file
features = parse_gff_gft(file_path)


 # get the transcript lenghts
transcript_lengths = parse_transcript_lengths('data/Araport11_genes.201606.cdna.len')

transcript_dict, transcript_exon_counts, gene_exon_counts, last_exon_for_transcript, \
        transcript_strands = generate_transcript_coordinates(features, transcript_lengths)


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
        """Test 2: AT1G01020.4 at pos 426 should be exon 7. Should have 7 exons."""
        query_transcript = "AT1G01020.4"
        query_position = 426
        exon_number, total_exons = query_transcript_exon(transcript_dict, query_transcript, query_position)
        print(f"{query_transcript} should have {total_exons} exons and position {query_position} should be exon {exon_number}")
        if exon_number is not None:
            print(f'Transcript {query_transcript} position {query_position} is in exon {exon_number}.')
            print(f'Transcript {query_transcript} has {total_exons} exons.')
        else:
            print(f'Position {query_position} is not found in transcript {query_transcript}.')
        self.assertEqual(exon_number, 7)
        self.assertEqual(total_exons, transcript_exon_counts[query_transcript])

# UTR
    def test_dict_3(self):
        """Test 3: AT1G01020.4 at pos 860 should be UTR. Should have 7 exons in the tans but 12 in the gene. NEGATIVE CODING GENE. DIFF TRANSCRIPT"""
        query_transcript = "AT1G01020.4"
        query_position = 860
        exon_number, total_exons = query_transcript_exon(transcript_dict, query_transcript, query_position)
        print(f"{query_transcript} should have {total_exons} exons and position {query_position} should be exon {exon_number}")
        if exon_number is not None:
            print(f'Transcript {query_transcript} position {query_position} is in exon {exon_number}.')
            print(f'Transcript {query_transcript} has {total_exons} exons.')
        else:
            print(f'Position {query_position} is not found in transcript {query_transcript}.')
        self.assertEqual(exon_number, "UTR")
        self.assertEqual(total_exons, transcript_exon_counts[query_transcript])


    def test_dict_4(self):
        """Test 4: TEST_neg_last_exon.1 at pos 5 should be exon 3(last one).
        Should have 3 exons. NEGATIVE CODING GENE. DIFF TRANSCRIPT"""
        query_transcript = "TEST_neg_last_exon.1"
        query_position = 5
        exon_number, total_exons = query_transcript_exon(transcript_dict, query_transcript, query_position)
        print(f"{query_transcript} should have {total_exons} exons and position {query_position} should be exon {exon_number}")
        if exon_number is not None:
            print(f'Transcript {query_transcript} position {query_position} is in exon {exon_number}.')
            print(f'Transcript {query_transcript} has {total_exons} exons.')
        else:
            print(f'Position {query_position} is not found in transcript {query_transcript}.')
        self.assertEqual(exon_number, 3)
        self.assertEqual(total_exons, 3)
        self.assertEqual(total_exons, transcript_exon_counts[query_transcript])


    def test_dict_5(self):
        """Test 5: TEST_1_3_UTR.1 at pos 5 should be exon 1. Should have 3 exons
        """
        query_transcript = "TEST_1_3_UTR.1"
        query_position = 5
        exon_number, total_exons = query_transcript_exon(transcript_dict, query_transcript, query_position)
        print(f"{query_transcript} should have {total_exons} exons and position {query_position} should be exon {exon_number}")
        if exon_number is not None:
            print(f'Transcript {query_transcript} position {query_position} is in exon {exon_number}.')
            print(f'Transcript {query_transcript} has {total_exons} exons.')
        else:
            print(f'Position {query_position} is not found in transcript {query_transcript}.')
        self.assertEqual(exon_number, 1)
        self.assertEqual(total_exons, 3)
        self.assertEqual(total_exons, transcript_exon_counts[query_transcript])


    def test_dict_6(self):
        """Test 6: TEST_1_3_UTR.1 at pos 25 should be exon 3(last one). 
        Should have 3 exons."""
        query_transcript = "TEST_1_3_UTR.1"
        query_position = 25
        exon_number, total_exons = query_transcript_exon(transcript_dict, query_transcript, query_position)
        print(f"{query_transcript} should have {total_exons} exons and position {query_position} should be exon {exon_number}")
        if exon_number is not None:
            print(f'Transcript {query_transcript} position {query_position} is in exon {exon_number}.')
            print(f'Transcript {query_transcript} has {total_exons} exons.')
        else:
            print(f'Position {query_position} is not found in transcript {query_transcript}.')
        self.assertEqual(exon_number, 3)
        self.assertEqual(total_exons, 3)
        self.assertEqual(total_exons, transcript_exon_counts[query_transcript])


    def test_dict_7(self):
        """Test 7: TEST_1_3_UTR.1 at pos 50 should be UTR. Should have 3 exons. 
        """
        query_transcript = "TEST_1_3_UTR.1"
        query_position = 50
        exon_number, total_exons = query_transcript_exon(transcript_dict, query_transcript, query_position)
        print(f"{query_transcript} should have {total_exons} exons and position {query_position} should be exon {exon_number}")
        if exon_number is not None:
            print(f'Transcript {query_transcript} position {query_position} is in exon {exon_number}.')
            print(f'Transcript {query_transcript} has {total_exons} exons.')
        else:
            print(f'Position {query_position} is not found in transcript {query_transcript}.')
        self.assertEqual(exon_number, "UTR")
        self.assertEqual(total_exons, 3)
        self.assertEqual(total_exons, transcript_exon_counts[query_transcript])


    def test_dict_8(self):
        """Test 8: TEST_exon1.1 at pos 100 should be exon1. Should have 3 exons."""
        query_transcript = "TEST_exon1.1"
        query_position = 100
        exon_number, total_exons = query_transcript_exon(transcript_dict, query_transcript, query_position)
        print(f"{query_transcript} should have {total_exons} exons and position {query_position} should be exon {exon_number}")
        if exon_number is not None:
            print(f'Transcript {query_transcript} position {query_position} is in exon {exon_number}.')
            print(f'Transcript {query_transcript} has {total_exons} exons.')
        else:
            print(f'Position {query_position} is not found in transcript {query_transcript}.')
        self.assertEqual(exon_number, 1)
        self.assertEqual(total_exons, 3)
        self.assertEqual(total_exons, transcript_exon_counts[query_transcript])


    def test_dict_9(self):
        """Test 9: TEST_last_exon.1 at pos 1048 should be exon2(last). Should have 2 exons. 
        """
        query_transcript = "TEST_exon1.1"
        query_position = 100
        exon_number, total_exons = query_transcript_exon(transcript_dict, query_transcript, query_position)
        print(f"{query_transcript} should have {total_exons} exons and position {query_position} should be exon {exon_number}")
        if exon_number is not None:
            print(f'Transcript {query_transcript} position {query_position} is in exon {exon_number}.')
            print(f'Transcript {query_transcript} has {total_exons} exons.')
        else:
            print(f'Position {query_position} is not found in transcript {query_transcript}.')
        self.assertEqual(exon_number, 1)
        self.assertEqual(total_exons, 3)
        self.assertEqual(total_exons, transcript_exon_counts[query_transcript])


if __name__ == '__main__':
    unittest.main()
