from dnachisel import *
import json
from objectives import E_COLI_WTS, H_SAPIENS_WTS, S_CEREVISIAE_WTS
from CAI import CAI

#how to write constriaints?

class AvoidHiddenStops(Specification):
	
	def __init__(self):
		raise NotImplementedError("AvoidHiddenStops not implmented")

##need constraint for alternative start sites
class AvoidAlternativeStarts(Specification):
	"""Make sure there aren't any RBS-like regions upstream of ATGs. Supports E. coli and yeast
	
	Parameters
	----------
	species
		A string defining which translation inition site consensus to use

	"""

	def __init__(self, species):
		raise NotImplementedError("AvoidAlternativeStarts not implmented")
		self.species = species

	def detect_ecoli(self, problem):
		"""
		This should just be looking for A/G rich a bit up stream of ATG
		"""
		pass

	def detect_yeast(self, problem):
		"""
		the TIS in budding yeast has a narrow consensus sequence

		Consensus0.5,−0.5 = N A\\T A\\(C|T) N\\G A\\T AUG G C\\(A|T))

		defined by the frequency weight matrix for HiConf TISs (Fig. 1), where A\\(G|U) denotes ‘A and not G and not T’
		"""
		pass

	def evaluate(self, problem):
		pass

	def __str__(self):
		return "Avoiding alternative starts for " + self.species


class ConstrainCAI(Specification):
	"""Enforce that the given sequence meets codon adaptation index requirements. Uses https://github.com/Benjamin-Lee/CodonAdaptationIndex/.

	Parameters
	----------
	species
		A string defining which wts file to use. Supported are "e_coli", "h_sapiens", and "s_cerevisiae"
	minimum
		The lower bound for CAI score. Upper bound is 1.0. Normally lower bound is 0.8. 
	location
		NOT YET IMPLEMENTED. Location of the DNA segment on which to enforce the CAI e.g.
	  ``Location(10, 45, 1)``
	"""
	best_possible_score = 1.0


	def __init__(self, species, location, minimum):
		#raise NotImplementedError("ConstrainCAI not implmented")
		self.location = location
		self.minimum = minimum
		self.species = species
		if species == "e_coli":
			self.weights = json.load(open(E_COLI_WTS))
		elif species == "h_sapiens":
			self.weights = json.load(open(H_SAPIENS_WTS))
		elif species == "s_cerevisiae":
			self.weights = json.load(open(S_CEREVISIAE_WTS))

	def evaluate(self, problem):
		""" return a CAI for the problem's sequence"""
		sequence = self.location.extract_sequence(problem.sequence)

		cai = CAI(sequence, weights=self.weights)
		
		#a constraint passes evaluation if the score is >= 0, so I need to map cai, 
		#which goes from 0 to 1 (and thus always passing) to XXX to 1 and store this as score
		#I have two points, (min, 0)_1 and (1,1)_2
		#The slope between these is (y2-y1)/(x2-x1) = (1-0)/(1-min) = 1/(1-min)
		#the line is y-y1 = m(x-x1), so y = (1/(1-min))(x-min)

		score = (1.0 / ( 1.0 - self.minimum)) * (cai - self.minimum)

		if cai >= self.minimum: 
			message = "Passed. CAI for %s meets requirements." % self.species
		else:
			message = "Failed. CAI for %s is %0.2f but must meet or exceed %0.2f" % (self.species, cai, self.minimum)

		return SpecEvaluation(
			self, problem,
			score=score,
			locations=[self.location],
			message=message
		)

	def __str__(self):
		"""String representation."""
		return "CAI for " + self.species