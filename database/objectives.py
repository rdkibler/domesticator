from dnachisel import *
from collections import Counter
from CAI import CAI
import json

class MinimizeKmerScore(Specification):
	"""Minimize a "no-kmers" score."""

	def __init__(self, k=9, boost=1.0):
		self.k = k
		self.boost = boost

	def evaluate(self, problem):
		"""Return a customized kmer score for the problem's sequence"""
		sequence = problem.sequence
		all_kmers = [sequence[i:i + self.k] for i in range(len(sequence) - self.k)]
		number_of_non_unique_kmers = sum([
			count
			for kmer, count in Counter(all_kmers).items()
			if count > 1
		])
		score = - (float(self.k) * number_of_non_unique_kmers) / len(sequence)
		return SpecEvaluation(
			self, problem,
			score=score,
			locations=[Location(0, len(sequence))],
			message="Score: %.02f (%d non-unique %d-mers)" % (
				score, number_of_non_unique_kmers, self.k
			)
		)

	def __str__(self):
		"""String representation."""
		return "Minimize%dmerScore" % self.k


#I think the best way to write this CAI thing is to pre-compute the weights file
#for each supported species and then just point to that in the database. I'll need
#a "compute_CAI_ref_weights.py" thing which computes weights from CAI.relative_adaptiveness

E_COLI_WTS = "/home/rdkibler/projects/domesticator/database/CAI/ecoli.heg.wts"
H_SAPIENS_WTS = "/home/rdkibler/projects/domesticator/database/CAI/hgTables.wts"
S_CEREVISIAE_WTS = "/home/rdkibler/projects/domesticator/database/CAI/scTable.wts"

#procedure for adding new species: go here:
#https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=716062121_ScFSjaHgLy4HsQVgDHKgdwv3ZJ2U&clade=other&org=0&db=0&hgta_group=genes&hgta_track=knownGene&hgta_table=knownGene&hgta_regionType=genome&position=&hgta_outputType=sequence&hgta_outFileName=
#and select your genome. Enter an output filename, hit "get output", select genomic & submit, then make sure only CDS EXONS is selected and finally "get sequence".
#Things ##should## just just work. 

#Then unzip and run that though make_weights.py in /home/rdkibler/projects/domesticator/database/CAI/make_weights.py


class MaximizeCAI(Specification):
	"""Maximizes CAI. Uses https://github.com/Benjamin-Lee/CodonAdaptationIndex/. """
	def __init__(self, species, boost=1.0):
		self.boost = boost
		if species == "e_coli":
			self.weights = json.load(open(E_COLI_WTS))
		elif species == "h_sapiens":
			self.weights = json.load(open(H_SAPIENS_WTS))
		elif species == "s_cerevisiae":
			self.weights = json.load(open(S_CEREVISIAE_WTS))

	def evaluate(self, problem):
		""" return a CAI for the problem's sequence"""
		sequence = problem.sequence

		score = CAI(sequence, weights=self.weights)

		return SpecEvaluation(
			self, problem,
			score=score,
			locations=[Location(0, len(sequence))],
			message="CAI: %.02f" % (score)
		)

	def __str__(self):
		"""String representation."""
		return "CAI"
