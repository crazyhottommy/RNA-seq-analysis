
#!/usr/bin/env python3


import json
from glob import glob
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--fastq_dir", help="Required. the full path to the fastq folder")
args = parser.parse_args()

assert args.fastq_dir is not None, "please provide the path to the fastq folder"

# the fastq_dir is the fastq folder contain many sub-folders. each subfolder's name is the sample name, and within each subfolder are the 
# fastq files for that sample
# check if folder exists:
if not os.path.exists(args.fastq_dir):
	print("fastq folder does not exist")
	exit()

# glob all the fastq files
fastqs = glob(args.fastq_dir + '/*/*fastq.gz')
FILES = {}

# Change this line to extract a sample name from each filename, the folder name is the sample name in my case
SAMPLES = [fastq.split('/')[-2] for fastq in fastqs]

for sample in SAMPLES:
    # Change 'R1' and 'R2' to match the way your mate pairs are marked.
    mate1 = lambda fastq: sample in fastq and 'R1' in fastq
    mate2 = lambda fastq: sample in fastq and 'R2' in fastq
    FILES[sample] = {}
    FILES[sample]['R1'] = sorted(filter(mate1, fastqs))
    FILES[sample]['R2'] = sorted(filter(mate2, fastqs))

js = json.dumps(FILES, indent = 4, sort_keys=True)
open('samples.json', 'w').writelines(js)
