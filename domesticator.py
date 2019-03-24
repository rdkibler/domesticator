#!/home/rdkibler/.conda/envs/domesticator/bin/python



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
from Bio.PDB import PDBParser, PPBuilder
import constraints
import objectives
import os


# NON_DNA_SYMBOLS = ['R','N','D','Q','E','H','I','L','K','M','F','P','S','W','Y','V', '*']
# NON_DNA_NON_PROTEIN_SYMBOLS = ['J','O','X','Z']

# def detect_input_mode(input_string):
# 	try:
# 		f = open(input_string, 'r')
# 		input_string = f.readline()
# 		if ">" in input_string:
# 			input_string = f.readline()
# 		is_file = True
# 		f.close()
# 	except FileNotFoundError:
# 		is_file = False
# 	except OSError:
# 		is_file = False

# 	input_type = None
# 	for char in input_string.upper():
# 		if char in NON_DNA_NON_PROTEIN_SYMBOLS:
# 			exit("input not recognized: %s" % input_string)
# 		if char in NON_DNA_SYMBOLS:
# 			if is_file:
# 				input_type = "protein_fasta"
# 			else:
# 				input_type = "protein_seq"
# 			break
# 	if not input_type:
# 		if is_file:
# 			input_type = "dna_fasta"
# 			exit("dna_fasta unimplemented")
# 		else:
# 			input_type = "dna_seq"

# 	return input_type



# def input_reader(input_string, input_mode):
# 	if input_mode == 'automatic':
# 		input_mode = detect_input_mode(input_string)
# 	elif input_mode == 'protein_seq':
# 		return Seq(reverse_translate(input_string), IUPAC.unambiguous_dna)
# 	elif input_mode == 'dna_seq':
# 		return Seq(input_string, IUPAC.unambiguous_dna)
# 	elif input_mode == 'protein_fasta':
# 		exit("protein_fasta not implemented")
# 	elif input_mode == 'dna_fasta':
# 		exit("dna_fasta not implemented")
# 	else:
# 		exit("unknown input mode")

	
		
# 		#gotta determine type



# 		exit()


# 	#before returning, make sure the seq is divisible by 3
# 	return #type should be Bio.Seq I think 



def load_template(filename):
	''' func descriptor '''

	vector = SeqIO.read(filename, "genbank")
	feats = [feat.qualifiers for feat in vector.features]
	keywords = {kw.strip("\"").split(":")[0] : kw.strip("\"").split(":")[1:] for kw in vector.annotations['keywords'][0].split(' ')}

	#prevent everything except destination from changing

	#Need to fined the dom_destination feature type which has the same label as destination_label
	#scream loudly if this wasn't found.

	#replace the destination sequence with the insert_seq

	#did the other features move around? 

	#find the constraints and objectives. How are these stored?
	#use keywords for global constraints/objectives
	#can I use features for local constraints/objectives?


	#This seq should be a SeqRecord object
	return vector



def replace_sequence_in_record(record, location, new_seq):
	#print(record, location, new_seq)
	insert = new_seq.seq
	if location.strand == 1:

		#print(dir(location))
		#print(location.extract(record.seq))
		#print(record.seq[location.start:location.end])
		adjusted_seq = record.seq[:location.start] + new_seq.seq + record.seq[location.end:]
		#exit(adjusted_seq)
		record.seq = adjusted_seq

		#print(help(location))
		#exit(dir(location))

		seq_diff = len(new_seq) - len(location)

		#adjust all features
		for feat in record.features:
			loc = feat.location
			
			earliest = min(loc.start, loc.end)
			latest = max(loc.start, loc.end)

			##gonna have to be careful about this

			if earliest > location.end and :
				feat.location += seq_diff

			if feat.location.start in location or feat.location.end in location:
				#they're overlapping
				#will have to create a whole new feature. Can't just modify this one unfortunately. 
				exit('unimplemented')

		return record

	else:
		exit("not implemented")


def insert_into_vector(vector, destination, new_seq):
	destination_annotation = None
	for feat in vector.features:
		if feat.qualifiers['label'][0] == destination:
			destination_annotation = feat

			#I will be replacing it so remove it:
			vector.features.remove(feat)
			break
	if not destination_annotation:
		exit("destination not found: " + destination)
	#print(destination_annotation)

	location = destination_annotation.location

	#print(vector)
	#print(vector.features)
	#print(dir(vector.features))
	vector = replace_sequence_in_record(vector, location, new_seq)


	#put in a new feature for this thing:
	#exit("not implemented")
	return vector



def vector_writer():
	pass


parser=argparse.ArgumentParser(prog='domesticator', description='TODO: write description')

parser.add_argument('--version', action='version', version='%(prog)s 0.2')

#Input Arguments
#parser.add_argument("--dna_sequence","-d",dest="dna_sequence",help="DNA sequence you want to optimize. Required and mutaually exclusive with --protein_sequence, --dna_fasta, and --protein_fasta.")
#parser.add_argument("--protein_sequence","-p",dest="protein_sequence",help="protein sequence you want to back translate and optimize. Required and mutually exclusive with --dna_sequence, dna_fasta, and protein_fasta")
#parser.add_argument("--dna_fasta","-D",dest="dna_fasta",help="fasta file containing dna sequences you want to optimize")
#parser.add_argument("--protein_fasta","-P",dest="protein_fasta",help="fasta file containing protein sequences you want to back translate and optimize")
#parser.add_argument("--sequence_name", "-n",dest="name",help="Use with --dna_sequence or --protein_sequence. Specify the name to be associated with the inputted sequence. default to %(default)s", default="gblock")
#parser.add_argument("--cds",dest="is_cds",action="store_true",default=True, help="I assume that any DNA sequence given to me is a CoDing Sequence (CDS) unless specified")
input_parser = parser.add_argument_group(title="Input Options", description=None)
input_parser.add_argument("input",							 			type=str, 	default=None, 			nargs="+",	help="DNA or protein sequence(s) or file(s) to be optimized. Valid inputs are full DNA or protein sequences or fasta or genbank files. Types are automatically determined. If this fails, set --input_mode to the input type.")
input_parser.add_argument("--vector", "-v", 		dest="vector", 		type=str, 	default=None, 			metavar="pEXAMPLE.gb",		help="HELP MESSAGE")
input_parser.add_argument("--destination", "-d", 	dest="destination", type=str, 	default="INSERT1", 		metavar="NAME",			help="TODO: flesh this out. Matches the dom_destination feature in the vector")
input_parser.add_argument("--input_mode", 			dest="input_mode", 	type=str, 	default="protein_fasta_file", 	help="Input mode. Protein fasta file by default.", choices=["PDB", "DNA_fasta_file", "protein_fasta_file", "DNA_sequence", "protein_sequence"])

#Cloning Arguments
#parser.add_argument("--cloning_scheme", "-c", dest="scheme", help="Specifies the overhang type and the specific patterns/restriction sites to avoid. Do not set flag if you don\'t want to apply a cloning scheme. Options are: \'MoClo-YTK\', [more coming soon] TODO: pet29_for_idt")
#parser.add_argument("--custom_cloning", dest="cloning_path",help="supply specifications for a custom cloning scheme")
#parser.add_argument("--part_type", "-t",dest="type",help="The type to domesticate the input sequence as according to the cloning_scheme selected")
#parser.add_argument("--right_flank",dest="right_flank", help="overrides the selected cloning scheme\'s 3\' flank")
#parser.add_argument("--left_flank",dest="left_flank", help="overrides the selected cloning scheme\'s 5\' flank")
#parser.add_argument("--primers_only",dest="primers_only", action="store_true",default=False,help="overrides all sequence optimization, outputs primers that will convert an existing sequence into a usable part. Incompatible with --protein_sequence and --protien_fasta")
#parser.add_argument("--conversion_primers",dest="types", help="report primers that can be used to convert the designed part into other parts, specified as a comma-separated list")
#parser.add_argument("--max_primer_len",dest="primer_len", type=int,default=60,help="use with --primers_only flag. default is %(default)dC")
#parser.add_argument("--primer_min_tm",dest="primer_min_tm", type=float,default=58.0,help="TODO:implement")

optimizer_parser = parser.add_argument_group(title="Optimizer Options", description=None)
#Optimization Arguments
optimizer_parser.add_argument("--avoid_kmers", type=int, dest="avoid_kmers", default=9, help="avoid repeated sequences greater than k bases to help with gene synthesis. Turn off by setting 0. Default to %(default)s")
optimizer_parser.add_argument("--avoid_kmers_boost", type=float, dest="kmer_boost", default=10.0, help="TODO: fill out this help. Default to %(default)s")
optimizer_parser.add_argument("--avoid_patterns", dest="avoid_patterns", help="DNA sequence patterns to avoid, listed as a comma separated list (no spaces!)")
optimizer_parser.add_argument("--avoid_restriction_sites", dest="avoid_restriction_sites", help="The names of the enzymes whose restriction sites you wish to avoid, such as EcoRI or BglII")
optimizer_parser.add_argument("--species", "-s", dest="species", default="e_coli", help="specifies the codon bias table to use as a comma separated list. Defaults to %(default)s", choices=["e_coli", "s_cerevisiae", "h_sapiens"])
#parser.add_argument("--codon_bias_table", dest="codon_bias_table_filename", help="overrides species, gives a custom codon bias table. See <NEED EXAMPLE> for example of formatting")
optimizer_parser.add_argument("--harmonized", dest="harmonized", help="This will tell the algorithm to choose codons with the same frequency as they appear in nature, otherwise it will pick the best codon as often as possible. Defaults to %(default)s", default=False, action="store_true")
optimizer_parser.add_argument("--avoid_homopolymers", dest="avoid_homopolymers", action="store_true", default=True, help="homopolymers can complicate synthesis. We minimize them by default, but you can turn this off")
optimizer_parser.add_argument("--avoid_hairpins", dest="avoid_hairpins", action="store_true", default=True, help="hairpins can cause problems during synthesis so this gives the option to avoid them. Default to %(default)s)")
optimizer_parser.add_argument("--optimize_terminal_GC_content", action="store_true", default=True, help="TODO: need to implement this AND write a helpful help")
optimizer_parser.add_argument("--constrain_CAI", dest="CAI_lower_bound", type=float, default=0.8, help="TODO: fill out help") ##UNCLEAR if this actually helps since CodonOptimize is this exact thing
#some argument for avoiding other kmer repeats?
#some argument for avoiding homopolymer?

#NOTE: using action="store_true" and default=True makes no sense. This is not an option

output_parser = parser.add_argument_group(title="Output Options", description=None)
#Output Arguments
output_parser.add_argument("--output_mode", dest="output_mode", default="print", help="supported output modes: print (to terminal) (default), fasta, or genbank. If no output_filename is supplied but a non-print option is selected, it'll just use the name of the input (from the fasta file) or \"domesticator_optimization\" if none specified.")
output_parser.add_argument("--output_filename", dest="output_filename", help="defaults to %(default)s.fasta or %(default)s.gb", default="domesticator_output")

args = parser.parse_args()

rec_counter = 1
inserts_to_optimize = []
if args.input_mode == "protein_fasta_file":
	for input_filename in args.input:
		for record in SeqIO.parse(input_filename, 'fasta'):
			record.seq = Seq(reverse_translate(record.seq), IUPAC.unambiguous_dna)
			inserts_to_optimize.append(record)
elif args.input_mode == "DNA_fasta_file":
	for input_filename in args.input:
		for record in SeqIO.parse(input_filename, 'fasta'):
			assert(len(record.seq) % 3 == 0)
			record.seq = Seq(str(record.seq), IUPAC.unambiguous_dna)
			inserts_to_optimize.append(record)
elif args.input_mode == "protein_sequnce":
	for input_sequence in args.input:
		record = SeqRecord(Seq(reverse_translate(input_sequence),IUPAC.unambiguous_dna), id="unknown_seq%d" % rec_counter, name="unknown_seq%d" % rec_counter, description="domesticator-optimized DNA sequence")
		rec_counter += 1
		inserts_to_optimize.append(record)

elif args.input_mode == "DNA_sequence":
	for input_sequence in args.input:
		record = SeqRecord(Seq(input_sequence,IUPAC.unambiguous_dna), id="unknown_seq%d" % rec_counter, name="unknown_seq%d" % rec_counter, description="domesticator-optimized DNA sequence")
		rec_counter += 1
		inserts_to_optimize.append(record)

elif args.input_mode == "PDB":
	chain="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
	parser = PDBParser()
	ppb=PPBuilder()
	for input_pdb in args.input:
		for chain_num, polypeptide in enumerate(ppb.build_peptides(parser.get_structure('name', input_pdb))):
			seq = Seq(reverse_translate(polypeptide.get_sequence()), IUPAC.unambiguous_dna)
			name = os.path.splitext(os.path.basename(input_pdb))[0] + "_" + chain[chain_num]
			record = SeqRecord(seq, id=name, name=name, description="domesticator-optimized DNA sequence")
			inserts_to_optimize.append(record)
else:
	exit("input mode not recognized: " + args.input_mode)
	
if args.vector:
	vector = load_template(args.vector)	
else:
	vector = SeqRecord()
	#empty vector

#now load all the custom global constraints and objectives?


for insert in inserts_to_optimize:
	vector = insert_into_vector(vector, args.destination, insert)
	print(vector)
	SeqIO.write([vector], "output.gb", "genbank")
	print("done")

	exit()





	constraints = []
	objectives = []
	#if a vector is provided, read it in, parse out the vector/cloning scheme-specific constraints, and insert initial sequence











exit()











#Handle Argument Requirements
#if args.primers_only is True and (args.protein_sequence or args.protein_fasta):
#	sys.exit("Error: can't create primers with protein sequences")
num_of_required_set = sum(((args.dna_sequence is not None), (args.protein_sequence is not None), (args.protein_fasta is not None), (args.dna_fasta is not None)))
if num_of_required_set == 0:
	sys.exit("I know \"required options\" are kind of dumb, but I use them so deal with it. I suppose I should implement this as a group."
		 "One of the following must be supplied: --dna_sequence, --protein_sequence, --dna_fasta, or --protein_fasta. See --help for "
		 "more information")
if num_of_required_set >1:
	sys.exit("Only one of the following may be supplied:  --dna_sequence, --protein_sequence, --dna_fasta, or --protein_fasta. "
		 "You supplied %d. See --help for more information" % num_of_required_set)


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



#note, since we're chainging from stupid json cloning schemes to genbank, make use of DnaChisel's
#built in find_specification_in_feature function under biotools. 


#load cloning scheme
cloning_file = ""
if args.scheme:
	cloning_file = database + "cloning_schema/" + args.scheme + ".json"
if cloning_file is not "":
	print("loading cloning scheme")
	with open(cloning_file) as clone_json:
		scheme = json.load(clone_json)
		print("scheme loaded from " + cloning_file)
		#if scheme['type_required'] and args.type is None:
		#	print("parts type required for cloning scheme %s" % cloning_file)
		#	print("options are: %s" % scheme['parts'])
		#	exit()
		#part_type = scheme['parts'][args.type]
		#left_flank = part_type["left_flank"]
		#right_flank = part_type["right_flank"]
		#override flanks if requested
		#if args.right_flank: right_flank = args.right_flank
		#if args.left_flank: left_flank = args.left_flank

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
			#sr.seq = left_flank + sr.seq + right_flank
			for f in sr.features:
				if f.type == "CDS":
					#f.location = f.location + len(left_flank)
					f.location = f.location

		#need to tell the optimizer it's not allowed to touch the flanks
		#This is probably not the best way to do it but it's what I got. 
		#local_constraints = [[AvoidChanges(Location(0,len(left_flank))), AvoidChanges(Location(len(sr.seq) - len(right_flank),len(sr.seq)))] + [AvoidPattern(pattern=pattern,location=Location(len(left_flank),len(sr.seq) - len(right_flank))) for pattern in patterns_to_avoid] for sr in SeqRecords]
		local_constraints = [[AvoidPattern(pattern=pattern,location=Location(len(left_flank),len(sr.seq) - len(right_flank))) for pattern in patterns_to_avoid] for sr in SeqRecords]

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
		if args.CAI_lower_bound > 0.0:
			global_constraints.append(constraints.ConstrainCAI(species=species, minimum=args.CAI_lower_bound))
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
