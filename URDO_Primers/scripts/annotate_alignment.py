#!/usr/bin/python

"""
The goal is to read in the alignment and annotated the headers with the hrv species.
"""

import sys, re, argparse
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
                "-i", "--fasta_to_annotate",
                help="Alignment in fasta format to annotate",
                )
	parser.add_argument(
		"-t", "--taxids",
		help="Vrial taxids from NCBI for annotation",
		)
	parser.add_argument(
		"-o", "--output",
		help="Output fasta with new headers",
		)
        args = parser.parse_args()
        return args

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------
def main():
        args = get_args()

	# Read in taxids for annotation
	dict_taxa = {}
	for astrline in open( args.taxids, "r" ):
		aastrline = astrline.strip().split('\t')
		if not re.search( '^##', astrline ):
			refseq, rep, taxonomy = aastrline[0], aastrline[1], aastrline[3]
			last_taxa = taxonomy.split(',')[-1]
			last_taxa_n = "hrv_" + last_taxa.split(' ')[-1] #These lines must be changed 
			if re.search( "[Rr]hinovirus", last_taxa ): #If doing a different organism
				dict_taxa[rep] = last_taxa_n 
				dict_taxa[refseq] = last_taxa_n

	# Read in fasta file to be annotated, and store seqrecords for writing out
	list_of_records = []
	i = 0
	for seqrecord in SeqIO.parse( args.fasta_to_annotate, "fasta" ):
		name = seqrecord.id.split('|')[3].split('.')[0]
		annot = ""
		if name in dict_taxa:
			annot = dict_taxa[name]
		else:
			annot = "weird"
		seqrecord.id = '|'.join( [annot] + seqrecord.id.split('|') )
		list_of_records.append( seqrecord )

	SeqIO.write( list_of_records, args.output, "fasta" )

if __name__ == "__main__":
        main()
