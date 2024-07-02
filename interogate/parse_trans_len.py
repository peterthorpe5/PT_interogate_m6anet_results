#!/usr/bin/env python3

import os
from collections import defaultdict


def test_line(line):
    """returns true lines. Not comments or blank line"""
    if not line.strip():
        return False  # if the last line is blank
    if line.startswith("#"):
        return False  # comment line
    if line.startswith("    # "):  # swarm result file
        return False  # comment line
    if line.startswith("		p"):
        return False  # comment line
    return line.rstrip()


def parse_transcript_lengths(length_file):
    """
    Parse the length file and return a defaultdict mapping transcript IDs to their lengths.

    Parameters:
    length_file (str): Path to the length file.

    Returns:
    defaultdict: defaultdict mapping transcript IDs to their lengths.
    """
    transcript_lengths = defaultdict(int)
    with open(length_file, 'r') as file:
        for line in file:
            if test_line(line):
                line = test_line(line)
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    transcript_id, length = parts
                    transcript_lengths[transcript_id] = int(length)

    return transcript_lengths