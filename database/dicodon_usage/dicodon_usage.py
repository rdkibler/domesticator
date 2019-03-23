from dnachisel import Specification, SpecEvaluation
import json
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio import SeqIO
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA
import itertools
from collections import defaultdict
import copy
from math import exp,log


def compute_dicodon_usage(sequence, table):
	"""return the sum of scores for every pair of adjacent codons"""
	seq = Seq(sequence, IUPACUnambiguousDNA)
	return(score(sequence,table))
	#return exp(score(sequence,table))

def score(sequence,table):
	seq_by_codon = ["".join(codon) for codon in grouper(3,sequence)]
	log_total = 0.0

	for first, second in zip(seq_by_codon, seq_by_codon[1:]):
		pair_freq = table[(first,second)]
		log_pair_freq = log(pair_freq)
		log_total += log_pair_freq
	return log_total
	
def grouper(n, iterable, fillvalue=None):
	args = [iter(iterable)] * n
	return itertools.zip_longest(*args, fillvalue=fillvalue)


def compute_all_possible_dicodons():
	bases = ["A","T","C","G"]
	codons = []
	for p1 in bases:
		for p2 in bases:
			for p3 in bases:
				codons.append(p1+p2+p3)
	stop_codons = ["TAA", "TGA", "TAG"]
	all_possible_dicodons = []
	for codon1 in codons:
		if codon1 not in stop_codons:
			for codon2 in codons:
				all_possible_dicodons.append((codon1,codon2))
	return all_possible_dicodons


def construct_dicodon_usage_table(ref_filename):

	standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
	standard_table.forward_table['TAG'] = '*'
	standard_table.forward_table['TAA'] = '*'
	standard_table.forward_table['TGA'] = '*'
	codon2aa = standard_table.forward_table
	aa2codon = defaultdict(lambda:[])
	for k in codon2aa.keys():
		aa2codon[codon2aa[k]].append(k)


	sequences = [str(seq.seq) for seq in SeqIO.parse(ref_filename, "fasta") if len(seq) % 3 == 0] #bit of a hack

	seqs_by_codon = [["".join(codon) for codon in grouper(3,seq)] for seq in sequences]

	#count_dicodon = initialize_dicodon_dict()
	count_dicodon = defaultdict(lambda:1)

	#count_diaa = defaultdict(lambda:defaultdict(lambda:1))
	count_diaa = {}
	for aa1 in aa2codon.keys():
		if aa1 != '*':
			for aa2 in aa2codon.keys():
				init_dicodon_count = {}
				for codon1 in aa2codon[aa1]:
					for codon2 in aa2codon[aa2]:
						init_dicodon_count[(codon1,codon2)] = 1
				count_diaa[(aa1,aa2)] = copy.deepcopy(init_dicodon_count)

	for seq_by_codon in seqs_by_codon:
		for first, second in zip(seq_by_codon, seq_by_codon[1:]):
			aa1 = codon2aa[first]
			aa2 = codon2aa[second]

			if (first,second) not in compute_all_possible_dicodons():
				message = "codon combination (%s, %s) not possible!" % (first,second)
				print(message)
				continue

			count_diaa[(aa1,aa2)][(first,second)] += 1
			#count_dicodon[(first,second)] += 1

	#print(count_diaa)


	dicodon_freq_table = {}
	for diaa, dicodons in count_diaa.items():
		#print(diaa)
		#print(dicodons)
		max_count = max(dicodons.values())
		dicodon_freqs = {k:float(v)/max_count for k,v in dicodons.items()}
		dicodon_freq_table.update(dicodon_freqs)
	#print(len(count_dicodon))
	#huh, this picks up 3 more than should be possible. Why?
	#print(max(count_dicodon.values()))
	return dicodon_freq_table

def load_dicodon_usage_table(table_filename):
	with open(table_filename, 'r') as fin:
		pretable = json.load(fin)
		table = {tuple(k.split("_")):v for k,v in pretable.items()}
		return table

def save_dicodon_usage_table(table, table_filename):
	with open(table_filename, 'w+') as fout:
		posttable = {"_".join(k):v for k,v in table.items()}
		json.dump(posttable, fout)


DICODON_TABLE_ECOLI = '/home/rdkibler/projects/domesticator/database/dicodon_usage/ecoli.heg.dicodon.wts'

#an objective for DNAChisel
class OptimizeDicodonUsage(Specification):
	"""Maximize a dicodon usage score."""

	def __init__(self, species, boost=1.0):
		self.species = species
		if self.species == "e_coli":
			self.table = load_dicodon_usage_table()
		else:
			exit("no dicodon usage table found for " + self.species)
		self.boost = boost

	def evaluate(self, problem):
		"""Return a score representing dicodon usage of a thing"""
		sequence = problem.sequence
		score = score(self.sequence,self.table)
	
		return SpecEvaluation(
			self, problem,
			score=score,
			locations=[Location(0, len(sequence))],
			message="Score: %.02f" % (score)
		)

	def __str__(self):
		"""String representation."""
		return "OptimizeDicodonUsage" % self.k


if __name__ == "__main__":
	import argparse
	import os

	parser=argparse.ArgumentParser(prog='dicodon_usage', 
		description='TODO: write description')

	parser.add_argument('--version', action='version', 
		version='%(prog)s 0.1')

	parser.add_argument("--compute_input","-i", type=str)
	parser.add_argument("--ref_fasta", type=str)
	parser.add_argument("--dicodon_usage_table", type=str)

	args = parser.parse_args()

	if not args.compute_input and not (args.ref_fasta or args.dicodon_usage_table):
		exit("bad input")

	table={}

	if args.ref_fasta:
		table = construct_dicodon_usage_table(args.ref_fasta)

		if not args.dicodon_usage_table:
			args.dicodon_usage_table = os.path.splitext(os.path.basename(args.ref_fasta))[0] + ".dicodon.wts"
		fout = args.dicodon_usage_table

		save_dicodon_usage_table(table, args.dicodon_usage_table)

	elif args.dicodon_usage_table:
		table = load_dicodon_usage_table(args.dicodon_usage_table)

	else:
		exit("need a ref_fasta or dicodon_usage_table")


	if args.compute_input:
		args.compute_input = args.compute_input.upper()
		print(compute_dicodon_usage(args.compute_input, table))