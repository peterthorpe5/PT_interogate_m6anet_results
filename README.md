# interogate_m6anet_results


Interogate m6Anet results is a Python package for analysing m6A (N6-methyladenosine) modification data and overlaying this with the GTF file. The package provides tools for parsing GTF files, generating transcript coordinates, and querying specific transcript information.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Testing](#testing)
- [Contributing](#contributing)
- [License](#license)

## Installation

### Prerequisites

- Python 3.7 or higher
- Required Python packages: `numpy`, `matplotlib`, `pandas`, `nose2`

### Clone the Repository

```bash
git clone https://github.com/peterthorpe5/interogate_m6anet_results.git
cd interogate_m6anet_results


## Install Required Packages
to test

pip install numpy matplotlib pandas nose nose2

```

# what does this do and how

1) this parses a gft fie (tested) gff3 (not yet tested) and sets up a LARGE dictionary of transcript exon number to
coordinates for the nucleotide sequence.  For exmaple:

AT1G01020.4 exon 2: [283, 284, 285]

In the gtf, the cooridnates are genomic locations, these dont directly help when mapping to the transcriptome. 

```bash
python interogate_m6anet.py

```

