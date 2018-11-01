from dnachisel import *
from collections import Counter

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




