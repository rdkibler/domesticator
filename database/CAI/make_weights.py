

import argparse
from CAI import relative_adaptiveness
import os
import json
from Bio import SeqIO


parser=argparse.ArgumentParser()
parser.add_argument("ref_fasta",help="fasta file containing reference sequences for given species")
args = parser.parse_args()

outfile = os.path.splitext(os.path.basename(args.ref_fasta))[0] + ".wts"

sequence = [str(seq.seq) for seq in SeqIO.parse(args.ref_fasta, "fasta") if len(seq) % 3 == 0] #bit of a hack
wts_dict = relative_adaptiveness(sequences = sequence)
with open(outfile, 'w+') as out:
	json.dump(wts_dict, out)