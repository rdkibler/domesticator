from dnachisel import *
from collections import Counter
from CAI import CAI
import json
from dnachisel.Location import Location


class MinimizeKmerScore(Specification):
	"""Minimize a "no-kmers" score."""

	def __init__(self, k=9, boost=1.0):
		self.k = k
		self.boost = boost

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
		
		sequence = self.location.extract_sequence(problem.sequence)

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

#viennaRNA/RNAlib
import RNA

class MinimizeSecondaryStructure(Specification):
	"""Uses ViannaRNA to calculate minimum free energy (MFE) over the specified location

	window
      Length of the sliding window, in nucleotides, for local SS stability.
      If not provided, the global SS stability is considered

    location
      Location objet indicating that the Specification only applies to a
      subsegment of the sequence. Make sure it is bigger than ``window``
      if both parameters are provided
      """
	def __init__(self, window=None, location=None, boost=1.0):
		self.boost = boost
		self.window = window
		if isinstance(location, tuple):
            location = Location.from_tuple(location)
        if location is not None and (location.strand == -1):
            location = Location(location.start, location.end, 1)
        self.location = location

    #def initialize_on_problem(self, problem, role=None):
    #    return self._copy_with_full_span_if_no_location(problem)
    #^ copied from EnforceGCContent. May be beneficial to calcualte an initial value



    def evaluate(self, problem):
        """Return deltaG over the location."""
        sequence = self.location.extract_sequence(problem.sequence)


        gc = gc_content(sequence, window_size=self.window)
        
        breaches = (np.maximum(0, self.mini - gc) +
                    np.maximum(0, gc - self.maxi))
        score = - breaches.sum()
        breaches_starts = wstart + (breaches > 0).nonzero()[0]

        if len(breaches_starts) == 0:
            breaches_locations = []
        elif len(breaches_starts) == 1:
            if self.window is not None:
                start = breaches_starts[0]
                breaches_locations = [[start, start + self.window]]
            else:
                breaches_locations = [[wstart, wend]]
        else:
            segments = [(bs, bs + self.window) for bs in breaches_starts]
            groups = group_nearby_segments(
                segments,
                max_start_spread=max(1,  self.locations_span))
            breaches_locations = [
                (group[0][0], group[-1][-1])
                for group in groups
            ]

        if breaches_locations == []:
            message = "Passed !"
        else:
            breaches_locations = [Location(*loc) for loc in breaches_locations]
            message = ("Out of bound on segments " +
                       ", ".join([str(l) for l in breaches_locations]))
        return SpecEvaluation(self, problem, score,
                              locations=breaches_locations,
                              message=message)






    def localized(self, location, problem=None, with_righthand=True):
        """Localize the GC content evaluation.
        For a location, the GC content evaluation will be restricted
        to [start - window, end + window]
        """
        # if self.location is not None:
        if self.window is None:
            # NOTE: this makes sense, but could be refined by computing
            # how much the local bounds should be in order not to outbound
            return self
        new_location = self.location.overlap_region(location)
        if new_location is None:
            return None 
 # VoidSpecification(parent_specification=self)
        else:
            extension = 0 if self.window is None else self.window - 1
            extended_location = location.extended(
                extension, right=with_righthand)
            new_location = self.location.overlap_region(extended_location)
        # else:
        #     if self.window is not None:
        #         new_location = location.extended(self.window + 1)
        #     else:
        #         new_location = None
        return self.copy_with_changes(location=new_location)

    def label_parameters(self):
        show_mini = self.mini is not None
        show_maxi = self.maxi is not None
        show_target = not (show_mini or show_maxi)
        show_window = self.window is not None

        return (
            show_mini * [('mini', "%.2f" % self.mini)] +
            show_maxi * [('maxi', "%.2f" % self.maxi)] +
            show_target * [('target',
                           ("%.2f" % self.target) if self.target else '')] +
            show_window * [('window',
                            ("%d" % self.window) if self.window else '')]
        )