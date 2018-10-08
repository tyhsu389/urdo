#!/usr/bin/python

"""
The goal of this script is to create the template needed for the primer3 program.
"""

# import
import argparse
from Bio import SeqIO

# Arguments
def get_args():
	"""
	Get arguments passed to script
	"""
	parser = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
		)
	parser.add_argument(
		"-c", "--consensus",
		help="consensus sequence for file."
		)
	parser.add_argument(
		"-o", "--output",
		help="output file for primer3."
		)
	parser.add_argument(
		"-p1", "--path",
		help="need path to libraries in primer3.",
		action="store_true"
		)
	parser.add_argument(
		"-p2", "--path_specified",
		help="path to libraries in primer3",
		default="/home/user/Downloads/programs/primer3-2.4.0/src/primer3_config/", #for CLC
		#default="/scicomp/home/s0th/Downloads/primer3-2.3.4/src/primer3_config/" for CDC
		)
	args = parser.parse_args()
	return args

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------
def main():
	args = get_args()
	
	fh = open( args.output, 'w' )

	# Parse consensus (if multiple)
	for seq_record in SeqIO.parse( args.consensus, "fasta" ):
		# Variables
		fh.write( "SEQUENCE_ID=" + seq_record.id + '\n' )
		fh.write( "SEQUENCE_TEMPLATE=" + str(seq_record.seq) + '\n' )

		# Constants
		list_of_vars = ["PRIMER_TASK=generic",
				"PRIMER_PRODUCT_SIZE_RANGE=150-230",
				"PRIMER_MIN_TM=53",
				"PRIMER_OPT_TM=55",
				"PRIMER_MAX_TM=58",
				"PRIMER_MAX_NS_ACCEPTED=1",
				"PRIMER_LIBERAL_BASE=1",
				"PRIMER_EXPLAIN_FLAG=1",
				"PRIMER_NUM_RETURN=20"]
		if args.path:
			list_of_vars.append( "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" + args.path_specified )
		for item in list_of_vars:
			fh.write( item + '\n' )
		# Split between consensus
		fh.write( "=" + '\n' )
	fh.close()

if __name__ == "__main__":
    main()
