#!/usr/bin/python

"""
The goal of this script is to get primers from the primer3 output.
"""

# import
import argparse
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

# Arguments
def get_args():
	"""
	Get arguments passed to script
	"""
	parser = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
		)
	parser.add_argument(
		"-i", "--input",
		help="Primer3 output",
		nargs="+",
		)
	parser.add_argument(
		"-n", "--num",
		help="Number of 'N's to add between primers",
		default=4,
		type=int,
		)
	parser.add_argument(
		"-e", "--emboss",
		help="Output file location for EMBOSS primer_search",
		default="emboss_primers.txt",
		)
	parser.add_argument(
		"-b", "--blast",
		help="Output fasta file location for BLAST",
		default="blast_primers.txt",
		)
	args = parser.parse_args()
	return args

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------
def main():
	args = get_args()

	# Read in primer3 output
	dict_primer3 = {}
	rec_list = []
	fh_emboss = open( args.emboss, "w" )
	for primers in args.input:
		condition = primers.split('/')[-1]
		for astrline in open( primers ):
			aastrline = astrline.strip().split('=')
			dict_primer3[ aastrline[0] ] = aastrline[1]

		# Get what we need
		num_primers = int( dict_primer3["PRIMER_LEFT_NUM_RETURNED"] )

		# Loop through primers
		for i in range( num_primers ):
			primer_name = "primer_" + condition +  '_' + str(i)
			left_sequence = dict_primer3["PRIMER_LEFT_" + str(i) + "_SEQUENCE"]
			right_sequence = dict_primer3["PRIMER_RIGHT_" + str(i) + "_SEQUENCE"]
			fh_emboss.write( '\t'.join( [primer_name, left_sequence, right_sequence] ) + '\n' )
			new_sequence = left_sequence + str( args.num*'N' ) + right_sequence 
			new_rec = SeqRecord( Seq( new_sequence, IUPACAmbiguousDNA() ),
					id=primer_name, name=primer_name )
			rec_list.append( new_rec )

	# Output fasta
	SeqIO.write( rec_list, args.blast, "fasta" )
	
if __name__ == "__main__":
	main()
