#!/home/rdkibler/.conda/envs/domesticator/bin/python




import unittest

class TestExample(self):
	#using assertEquals, assertTrue, etc from unittest is nice b/c the error messages are good

	def setUp(self): 
		#runs before every test

	def tearDown(self):
		#runs after every test

	def test_basic(self):
		self.assertEquals(True,True)

##also look into pytest






#New features I want:
#	only optimize coding sequences
#	optimize against hairpins globally (even between coding and non-coding sequences)



#What I want this to do is take a protein sequence and convert it into an orderable gblock
#to do this, it'll need 
#	1. the protein sequence, 
#	2. the cloning scheme being used (MoClo-YTK, EMMA, BioBricks, etc)
#	3. specific identification required by the clonging scheme (ie MoClo-YTK type) 
#	Optional override to cloning scheme and identifiaction if you just specify the left and right flanking sequences 
#	4. the destination organism (to optimize codon usage), 
#	5. a list of sequences to avoid (like restriction sites) (default to avoid the type-II restriction sites or whatever sites used by whatever standards)
#	6. the name of the sequence (has a default value)
#	6. a flag specifiying output to be in the form of a fasta file (in which case you must also specifiy the file name), gb file, or fasta print to terminal
#		if you write to gb file, then it'll include the cloning system you specify, the part type (or whatever), and specify wether or not the overhangs are custom


#default behavior will be to codon optimize a cds
#I plan on optimizing I definitely want to be able to optimize a sequence for synthesis/GC content while preserving patterns, however, I will not be implementing that right away

print("importing libraries")
import sys

database='/home/rdkibler/projects/domesticator/database/'
sys.path.insert(0, database)

import json
from collections import Counter
from dnachisel import *
import argparse
import warnings
from Bio import SeqIO
from Bio.SeqFeature import *
from Bio.Seq import (MutableSeq, Seq)
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import constraints
import objectives


print("parsing arguments")
parser=argparse.ArgumentParser()

#Input Arguments
parser.add_argument("--dna_sequence","-d",dest="dna_sequence",help="DNA sequence you want to optimize. Required and mutaually exclusive with --protein_sequence, --dna_fasta, and --protein_fasta.")
parser.add_argument("--protein_sequence","-p",dest="protein_sequence",help="protein sequence you want to back translate and optimize. Required and mutually exclusive with --dna_sequence, dna_fasta, and protein_fasta")
parser.add_argument("--dna_fasta","-D",dest="dna_fasta",help="fasta file containing dna sequences you want to optimize")
parser.add_argument("--protein_fasta","-P",dest="protein_fasta",help="fasta file containing protein sequences you want to back translate and optimize")
parser.add_argument("--sequence_name", "-n",dest="name",help="Use with --dna_sequence or --protein_sequence. Specify the name to be associated with the inputted sequence. default to %(default)s", default="gblock")
parser.add_argument("--cds",dest="is_cds",action="store_true",default=True, help="I assume that any DNA sequence given to me is a CoDing Sequence (CDS) unless specified")

#Cloning Arguments
parser.add_argument("--cloning_scheme", "-c", dest="scheme", help="Specifies the overhang type and the specific patterns/restriction sites to avoid. Do not set flag if you don\'t want to apply a cloning scheme. Options are: \'MoClo-YTK\', [more coming soon] TODO: pet29_for_idt")
#parser.add_argument("--custom_cloning", dest="cloning_path",help="supply specifications for a custom cloning scheme")
parser.add_argument("--part_type", "-t",dest="type",help="The type to domesticate the input sequence as according to the cloning_scheme selected")
parser.add_argument("--right_flank",dest="right_flank", help="overrides the selected cloning scheme\'s 3\' flank")
parser.add_argument("--left_flank",dest="left_flank", help="overrides the selected cloning scheme\'s 5\' flank")
#parser.add_argument("--primers_only",dest="primers_only", action="store_true",default=False,help="overrides all sequence optimization, outputs primers that will convert an existing sequence into a usable part. Incompatible with --protein_sequence and --protien_fasta")
#parser.add_argument("--conversion_primers",dest="types", help="report primers that can be used to convert the designed part into other parts, specified as a comma-separated list")
#parser.add_argument("--max_primer_len",dest="primer_len", type=int,default=60,help="use with --primers_only flag. default is %(default)dC")
#parser.add_argument("--primer_min_tm",dest="primer_min_tm", type=float,default=58.0,help="TODO:implement")

#Optimization Arguments
parser.add_argument("--avoid_kmers", type=int, dest="avoid_kmers", default=9, help="avoid repeated sequences greater than k bases to help with gene synthesis. Turn off by setting 0. Default to %(default)s")
parser.add_argument("--avoid_kmers_boost", type=float, dest="kmer_boost", default=10.0, help="TODO: fill out this help. Default to %(default)s")
parser.add_argument("--avoid_patterns", dest="avoid_patterns", help="DNA sequence patterns to avoid, listed as a comma separated list (no spaces!)")
parser.add_argument("--avoid_restriction_sites", dest="avoid_restriction_sites", help="The names of the enzymes whose restriction sites you wish to avoid, such as EcoRI or BglII")
parser.add_argument("--species", "-s", dest="species", default="e_coli", help="specifies the codon bias table to use as a comma separated list. Options are: e_coli, s_cerevisiae. Defaults to %(default)s")
parser.add_argument("--codon_bias_table", dest="codon_bias_table_filename", help="overrides species, gives a custom codon bias table. See <NEED EXAMPLE> for example of formatting")
parser.add_argument("--harmonized", dest="harmonized", help="This will tell the algorithm to choose codons with the same frequency as they appear in nature, otherwise it will pick the best codon as often as possible. Defaults to %(default)s", default=False, action="store_true")
parser.add_argument("--avoid_homopolymers", dest="avoid_homopolymers", action="store_true", default=True, help="homopolymers can complicate synthesis. We minimize them by default, but you can turn this off")
parser.add_argument("--avoid_hairpints", dest="avoid_hairpins", action="store_true", default=True, help="hairpins can cause problems during synthesis so this gives the option to avoid them. Default to %(default)s)")
parser.add_argument("--optimize_terminal_GC_content", action="store_true", default=True, help="TODO: need to implement this AND write a helpful help")
#some argument for avoiding other kmer repeats?
#some argument for avoiding homopolymer?

#Output Arguments
parser.add_argument("--output_mode", dest="output_mode", default="print", help="supported output modes: print (to terminal) (default), fasta, or genbank. If no output_filename is supplied but a non-print option is selected, it'll just use the name of the input (from the fasta file) or \"gblock###\" if none specified.")
parser.add_argument("--output_filename", dest="output_filename", help="defaults to %(default)s.fasta or %(default)s.gb", default="gblocks")


args = parser.parse_args()







#Handle Argument Requirements
if args.primers_only is True and (args.protein_sequence or args.protein_fasta):
	sys.exit("Error: can't create primers with protein sequences")
num_of_required_set = sum(((args.dna_sequence is not None), (args.protein_sequence is not None), (args.protein_fasta is not None), (args.dna_fasta is not None)))
if num_of_required_set == 0:
	sys.exit("I know \"required options\" are kind of dumb, but I use them so deal with it. I suppose I should implement this as a group."
		 "One of the following must be supplied: --dna_sequence, --protein_sequence, --dna_fasta, or --protein_fasta. See --help for "
		 "more information")
if num_of_required_set >1:
	sys.exit("Only one of the following may be supplied:  --dna_sequence, --protein_sequence, --dna_fasta, or --protein_fasta. "
		 "You supplied %d. See --help for more information" % num_of_required_set)
if args.cloning_path:
	sys.exit("custom cloning schemes are not yet supported")

if args.protein_sequence or args.protein_fasta:
	args.is_cds = True

print("loading sequences")
global_constraints = []
global_objectives = []
SeqRecords = []

#Load in the sequences
if args.dna_sequence:
	SeqRecords.append(SeqRecord(Seq(args.dna_sequence, IUPAC.unambiguous_dna),name=args.name,id=args.name))
elif args.protein_sequence:
	SeqRecords.append(SeqRecord(Seq(reverse_translate(args.protein_sequence), IUPAC.unambiguous_dna),name=args.name,id=args.name))
elif args.dna_fasta:
	SeqRecords = list(SeqIO.parse(args.dna_fasta,"fasta",alphabet=IUPAC.unambiguous_dna))
elif args.protein_fasta:
	SeqRecords = list(SeqIO.parse(args.protein_fasta,"fasta",alphabet=IUPAC.protein))
	for sr in SeqRecords:
		sr.seq = Seq(reverse_translate(sr.seq),alphabet=IUPAC.unambiguous_dna)

#TODO for coding sequences, store the protein translation and compare it before and after optimization. Not sure if it's getting changed

local_constraints = [None] * len(SeqRecords)

for sr in SeqRecords:
	#TODO check if it has a name. If not, give it a default name and a number
	if args.is_cds:
		print("SHOULD I BE USING FEATURES OR ANNOTATIONS? Dude. Fix this.")
		sr.features.append(SeqFeature(location=FeatureLocation(0,len(sr.seq)),type="CDS",strand=1,id=sr.id)) #indicate cds
		sr.annotations["species"]= "optimized for %s" % (args.species)
		sr.annotations["source"] = "de novo"
		if len(sr.seq) % 3 is not 0:
			sys.exit("ERROR: cds flag is set but %s (len = %d) is not divisible by 3" % (sr.name, len(sr.seq)))

print("loaded: ")
for sr in SeqRecords:
	print(sr.id)



#load cloning scheme
cloning_file = ""
if args.cloning_path is not None:
	cloning_file = args.cloning_path
elif args.scheme is not None:
	cloning_file = database + "cloning_schema/" + args.scheme + ".json"
if cloning_file is not "":
	print("loading cloning scheme")
	with open(cloning_file) as clone_json:
		scheme = json.load(clone_json)
		print("scheme loaded from " + cloning_file)
		if scheme['type_required'] and args.type is None:
			print("parts type required for cloning scheme %s" % cloning_file)
			print("options are: %s" % scheme['parts'])
			exit()
		part_type = scheme['parts'][args.type]
		left_flank = part_type["left_flank"]
		right_flank = part_type["right_flank"]
		#override flanks if requested
		if args.right_flank: right_flank = args.right_flank
		if args.left_flank: left_flank = args.left_flank

		patterns_to_avoid = scheme['patterns_to_avoid'] + [enzyme_pattern(enzyme) for enzyme in scheme['enzymes_to_avoid']]
		#for pattern in scheme['patterns_to_avoid']:
		#	global_constraints.append(AvoidPattern(pattern))
		#for enzyme in scheme['enzymes_to_avoid']:
		#	global_constraints.append(AvoidPattern(enzyme_pattern(enzyme)))

		for sr in SeqRecords:
			#TODO
			#store part type in annotation/whatever
			#adjust cds (if cds) in record
			#NOTE: try using "SHIFT" in the SeqFeature class
			#apply flanking sequences, keeping track of CDS and shifting location as necessary
#			print(sr.seq)
			sr.seq = left_flank + sr.seq + right_flank
			for f in sr.features:
				if f.type == "CDS":
					f.location = f.location + len(left_flank)

		#need to tell the optimizer it's not allowed to touch the flanks
		#This is probably not the best way to do it but it's what I got. 
		local_constraints = [[AvoidChanges(Location(0,len(left_flank))), AvoidChanges(Location(len(sr.seq) - len(right_flank),len(sr.seq)))] + [AvoidPattern(pattern=pattern,location=Location(len(left_flank),len(sr.seq) - len(right_flank))) for pattern in patterns_to_avoid] for sr in SeqRecords]
		
		####MUST DO add EnforceTranslation(location=Location(len(left_flank),len(sr.seq) - len(right_flank))) which works with the type, if it is a type that gets translated. Shit


	print("scheme applied")
else:
	print("no cloning scheme selected")
	if args.is_cds:
		global_constraints.append(EnforceTranslation()) #Does this need a location??

print("loading global constraints")
#load constraints
if args.avoid_homopolymers:
	global_constraints += [
                AvoidPattern(homopolymer_pattern("A",6)),
                AvoidPattern(homopolymer_pattern("T",6)),
                AvoidPattern(homopolymer_pattern("G",6)),
                AvoidPattern(homopolymer_pattern("C",6))]
if args.avoid_hairpins:
	global_constraints += [AvoidHairpins()]

global_constraints += [
                EnforceGCContent(0.4,0.65), #global
                EnforceGCContent(0.25,0.8,window=50)] #local

if args.avoid_restriction_sites:
	print("avoiding:")
	for enzy in args.avoid_restriction_sites.split(","):
		print(enzy)
		global_constraints.append(AvoidPattern(enzyme_pattern(enzy)))

print("loading global objectives")
#load objectives
if args.avoid_kmers > 0:
	global_objectives.append(objectives.MinimizeKmerScore(k = args.avoid_kmers, boost = args.kmer_boost))
if args.is_cds:
	print("optimizing for:")
	for species in args.species.split(","):
		print(species)
		if args.harmonized:
			global_objectives.append(CodonOptimize(species,mode='harmonized')) #NEEDS LOCATION WITH CDS
		else:
			global_objectives.append(CodonOptimize(species,mode='best_codon')) #NEEDS LOCATION WITH CDS
	#global_objectives.append(CodonOptimize(species=args.species,mode='best_codon')) #NEEDS LOCATION WITH CDS
print("begin optimization")
#start optimization


#print(global_constraints)
#print(global_objectives)

print("NOTE: I can easily parallelize this I think")
print("NOTE: think about implementing ideas from https://www.nature.com/articles/nbt.4238")

for sr, lc in zip(SeqRecords, local_constraints):
	print("Optimizing %s" % sr.name)
	print(len(sr.seq))
	this_constraints = global_constraints
	if(lc is not None):
		this_constraints += lc
	this_objectives = global_objectives

	print("objectives: %s" % this_objectives)
	print("constraints: %s" % this_constraints)

	problem = DnaOptimizationProblem(
	    sequence=str(sr.seq),
	    constraints=this_constraints,
	    objectives=this_objectives
	)
	
	#print(problem.sequence)
	#print ("\nBefore optimization:\n")
	print (problem.constraints_text_summary(failed_only=True))
	#print (problem.objectives_text_summary())


	#This seems like a great place for GPU calculations
	problem.resolve_constraints()
	problem.optimize()
	problem.resolve_constraints(final_check=True)
	
	#print(problem.sequence)
	if  args.is_cds:
		original_prot = sr.seq.translate()
		optimize_prot = Seq(problem.sequence).translate()
		if(original_prot != optimize_prot):
			print("protein sequence changed before and after optimization")
			print(original_prot)
			print(optimize_prot)
			exit()
	print("Old seq: %s" % sr.seq)
	sr.seq = MutableSeq(problem.sequence, alphabet=IUPAC.unambiguous_dna)
	print("New seq: %s" % sr.seq)
	#print ("\nAfter optimization:\n")
	print (problem.constraints_text_summary(failed_only=True))
	
	#do something else if it's failed!
	#print (problem.objectives_text_summary())
print("finish optimization")
print("outputting")

outname = args.output_filename

args.output_mode = args.output_mode.strip().lower()

#outext = ""
#TODO this is messy. Use BioPython's built in type checker and just clean up input a bit (lower case, etc). Catch errors and then default to terminal.
if args.output_mode == "fasta":
	if outname is None:
		outname = "gblocks.fasta"
	#outext = ".fasta"
	for sr in SeqRecords:
		with open(outname, 'w') as f:
			f.write(sr.format("fasta"))
elif args.output_mode == "genbank":
	if outname is None:
		outname = "gblocks.gb"
	#outext = ".gb"
	sys.exit("not implemented")
else:
	if args.output_mode != "print":
		print("output mode not recognized")
	print("Printing to terminal")
	for sr in SeqRecords:
		print(sr.format("fasta"))
