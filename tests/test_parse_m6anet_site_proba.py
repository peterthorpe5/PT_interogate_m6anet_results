#!/usr/bin/env python

"""Tests of m6anet site proba parsing  functionality"""

import os
import unittest
import pandas as pd
import tempfile
import os
from interogate.parse_m6a_site_proba import identify_methylated_sites


# 3 exon 1, last exon and UTR
test_data = """TEST.1,5,5,0.9662596747279167,TGACA,0.1395348837209302
TEST.1,35,35,0.9662596747279167,TGACA,0.1395348837209302
TEST.1,35,35,0.0662596747279167,TGACA,0.1395348837209302
TEST.1,50,50,0.0662596747279167,TGACA,0.1395348837209302""".split()



class TestIdentifyMethylatedSites(unittest.TestCase):

    def setUp(self):
        # Reduced sample data to create a temporary CSV file
        self.data = """transcript_id,transcript_position,n_reads,probability_modified,kmer,mod_ratio
AT1G01100.2,170,43,0.0662596747279167,TGACA,0.1395348837209302
AT1G01100.2,262,50,0.2530084550380707,TGACT,0.2400000000000000
AT1G01100.2,428,51,0.1141563579440117,TGACT,0.2941176470588235
AT1G01100.2,475,47,0.9452261924743652,AGACT,0.6170212765957447
AT1G01100.2,525,36,0.9956660866737366,GGACT,0.8055555555555556
AT1G01100.2,559,36,0.6470550298690796,TAACT,0.4722222222222222
AT1G01100.2,577,25,0.7123230099678040,AAACC,0.5200000000000000
AT1G01090.1,278,170,0.0189953744411469,AGACT,0.0294117647058824
AT1G01090.1,375,238,0.0949225723743439,TAACC,0.0504201680672269
AT1G01090.1,433,236,0.2292162925004959,AGACA,0.3559322033898305
AT1G01090.1,541,223,0.1027148067951202,TGACT,0.0762331838565022
AT1G01090.1,565,236,0.2193052470684052,TGACC,0.2415254237288136
AT1G01090.1,875,234,0.0939915999770164,AAACT,0.1452991452991453
AT1G01090.1,904,258,0.1621683388948441,TAACT,0.2286821705426356
AT1G01090.1,946,240,0.0135922869667411,TGACC,0.0416666666666667
AT1G01090.1,1092,281,0.0827838256955147,AGACT,0.1743772241992882
AT1G01090.1,1174,267,0.0537529326975346,AGACC,0.0262172284644195
AT1G01090.1,1551,225,0.8826547265052795,AAACT,0.7955555555555556
AT1G01090.1,1600,280,0.9126547265052795,AAACT,0.8955555555555556
AT1G01090.1,1700,290,0.9626547265052795,AAACT,0.9455555555555556"""

        # Create a temporary CSV file with the sample data
        self.temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.csv')
        self.temp_file.write(self.data.encode())
        self.temp_file.close()
        
    def tearDown(self):
        # Remove the temporary file
        os.remove(self.temp_file.name)

    def test_identify_methylated_sites(self):
        threshold = 0.9
        expected_result = pd.DataFrame({
            'transcript_id': ['AT1G01100.2', 'AT1G01100.2', 'AT1G01090.1', 'AT1G01090.1'],
            'transcript_position': [475, 525, 1600, 1700]
        })
        
        result = identify_methylated_sites(self.temp_file.name, threshold)
        
        # Print the result for debugging purposes
        print(result)

        pd.testing.assert_frame_equal(result.reset_index(drop=True), expected_result)

if __name__ == '__main__':
    unittest.main()
