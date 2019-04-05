##WASTEFUL to import *
from dnachisel import * 
from collections import Counter
from CAI import CAI
import json
from dnachisel.Location import Location


class MinimizeKmerScore(Specification):
	"""Minimize a "no-kmers" score."""

	def __init__(self, k=9, location=None, boost=1.0):
		self.location = location
		self.k = k
		self.boost = boost

	def initialize_on_problem(self, problem, role=None):
		return self._copy_with_full_span_if_no_location(problem)
		# if self.location is None:
		#     location = Location(0, len(problem.sequence))
		#     return self.copy_with_changes(location=location)
		# else:
		#     return self

	def evaluate(self, problem):
		"""Return a customized kmer score for the problem's sequence"""
		sequence = self.location.extract_sequence(problem.sequence)
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
			locations=[self.location],
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
	"""Maximizes CAI. Uses https://github.com/Benjamin-Lee/CodonAdaptationIndex/. 
	
	"""
	def __init__(self, species, location=None, boost=1.0):
		self.boost = boost
		self.location = location
		if species == "e_coli":
			self.weights = json.load(open(E_COLI_WTS))
		elif species == "h_sapiens":
			self.weights = json.load(open(H_SAPIENS_WTS))
		elif species == "s_cerevisiae":
			self.weights = json.load(open(S_CEREVISIAE_WTS))


	def initialize_on_problem(self, problem, role=None):
		return self._copy_with_full_span_if_no_location(problem)
		# if self.location is None:
		#     location = Location(0, len(problem.sequence))
		#     return self.copy_with_changes(location=location)
		# else:
		#     return self

	def evaluate(self, problem):
		""" return a CAI for the problem's sequence"""

		sequence = self.location.extract_sequence(problem.sequence)

		score = CAI(sequence, weights=self.weights)

		return SpecEvaluation(
			self, problem,
			score=score,
			locations=[self.location],
			message="CAI: %.02f" % (score)
		)

	def __str__(self):
		"""String representation."""
		return "CAI"


class MaximizeDicodonAdaptiveIndex(Specification):
	def __init__(self):
		raise NotImplementedError("AvoidAlternativeStarts not implmented")


#import RNA
try:
	from RNA import fold_compound, OPTION_MFE, OPTION_WINDOW
except:
	exit("please install RNAlib (https://www.tbi.univie.ac.at/RNA/)")

class MinimizeSecondaryStructure(Specification):
	"""Uses ViannaRNA to calculate minimum free energy (MFE) over the specified location

	max_energy
		How strong can the feature be before it gets flagged

	window
	  Length of the sliding window, in nucleotides, for local MFE.
	  If not provided, the global MFE is considered

	temperature
		parameter for ViennaRNA. The temperature to predict secondary structure at

	location
	  Location objet indicating that the Specification only applies to a
	  subsegment of the sequence. Make sure it is bigger than ``window``
	  if both parameters are provided
	  """
	def __init__(self, max_energy=-5.0, location=None, boost=1.0):
		self.max_e = max_energy
		self.boost = boost
		self.window = window

		if isinstance(location, tuple):
			location = Location.from_tuple(location)
		if location is not None and (location.strand == -1):
			location = Location(location.start, location.end, 1)
		self.location = location

	def initialize_on_problem(self, problem, role=None):
		##I will need to implement this better probably?

		#will this do the whole seq? Also, will RNA tolerate DNA?
		self._copy_with_full_span_if_no_location(problem)
		return self
   
	def mfe_window_callback(self,start,end,structure,energy,data=None):
		if energy < self.max_e:
			data.append({'structure':structure,'start':start,'end':end,'energy':energy})

	def evaluate(self, problem):
		"""Return deltaG over the location."""
		fc = fold_compound(problem.sequence, None, OPTION_MFE | OPTION_WINDOW)
		
		hairpins = []

		mfe = fc.mfe_window_cb(mfe_window_callback, hairpins)

		message = "overall mfe: %0.2f kcal/mol; num violations: %d" % (mfe, len(hairpins))
		
		#locations are defined by lists of the start and stop of the segment with the problem
		hairpin_locations = [(hp['start'], hp['end']) for hp in hairpins]
		
		return SpecEvaluation(self, problem, mfe,
							  locations=hairpin_locations,
							  message=message)

	#do I need to implement localized?

	def label_parameters(self):

		self.max_e
		self.boost

		return (
					[('max_energy', str(self.max_e))]
		)