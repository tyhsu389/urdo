#!/usr/bin/python

"""
The goal of this script is to generate a text file for Krona input.

Format:
Count	Genus	Species	Type(characterization,detection,pan)	kmer(if possible)

To get this, take as input:

"""

import sys
import re

# ---------------
# cli
# ---------------

# ---------------
# main
# ---------------
def main( ):
	# import kmer counts for matches and partial matches (Hamming)
	
	all_info = []
	total_kmers = 0
	for astrline in open( sys.argv[1] ): # get matches
		kmer, count, genus, species, ctype = astrline.strip().split('\t')
		newline = [count, genus, species, ctype, kmer]
		all_info.append( newline )
		total_kmers += int( count )

	for bstrline in open( sys.argv[2] ): # get hamming matches
		taxon, avgHD, count = bstrline.strip().split('\t')
		genus, species = "", ""
		if re.search( "Pan ", taxon ): #this is a pan organism
			genus = taxon.replace( "Pan ", "" )
			species = "-"
		else:
			genus, species = taxon.rsplit(' ', 1)
		newline = [count, genus, species, 'hamming', '']
		all_info.append( newline )
		total_kmers += int( count )

	dict_totals = {} # remainder will be unclassified
	for totals in open( sys.argv[3] ):
		sample, total = totals.strip().split('\t')
		if re.search( sample, sys.argv[1] ):
			mismatched = int( total ) - total_kmers
			all_info.append( [mismatched, '', '', '', ''] )

	# write out txt file
	for line in all_info:
		print( '\t'.join( str(x) for x in line ) )


if __name__ == "__main__":
	main( )
