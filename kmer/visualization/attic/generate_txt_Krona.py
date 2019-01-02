#!/usr/bin/python

"""
The goal of this script is to generate a text file for Krona input.

Format:
Count	Genus	Species	Type(characterization,detection,pan)	Char_Type(characterization value)

To get this, take as input:
1. detection_meta (counts for detection)
2. hamming (counts for hamming)
3. characterization (counts for characterization types)
4. total counts (to determine uncharacterized values)

"""
import argparse
import sys
import re

# ---------------
# cli
# ---------------
def get_args( ):
	""" Get arguments passed to script """
	parser = argparse.ArgumentParser(
		formatter_class = argparse.ArgumentDefaultsHelpFormatter
		)
	parser.add_argument( 
		"-a", "--abundances",
		required=True,
		help="Counts for exact matches, usually in /abundances/*_extended"
		)
	parser.add_argument(
		"-m", "--hamming",
		help="Counts for hamming matches, usually in /hamming/*_hamming"
		)
	parser.add_argument(
		"-c", "--character",
		help="Counts for characterization matches, usually in /characterization/*_extended"
		)
	parser.add_argument(
		"-t", "--total",
		help="Counts for total number of kmers from reads minus duplicates"
		)
	args = parser.parse_args()
	return args

def find_genus_sp( taxon ):
	genus, species = "", ""
	if re.search( "Pan ", taxon ): #this is a pan taxon
		genus = taxon.replace( "Pan ", "" )
		species = "pan"
	else:
		genus, species = taxon.rsplit( ' ', 1 )
	return genus, species

def add_to_dict( dict_count, genus, species, ctype, ptype, count ):
	if genus in dict_count:
		if species in dict_count[genus]:
			if ctype == "characterization":
				if ctype in dict_count[genus][species]:
					dict_count[genus][species][ctype].setdefault( ptype, [] ).append( float( count ) )
				else:
					dict_count[genus][species][ctype] = { ptype: [float( count )] }
			else:
				if ctype in dict_count[genus][species]:
					print( 'weird' )
				dict_count[genus][species][ctype] = float( count )
		else:
			if ctype == "characterization":
				dict_count[genus][species] = {ctype: {ptype: [float( count )]}}
			else:
				dict_count[genus][species] = {ctype: float( count )}
	else:
		if ctype == 'characterization':
			dict_count[genus] = {species: {ctype: {ptype: [float( count )]} } }
		else:
			dict_count[genus] = {species: {ctype: float( count )}}
	return dict_count

# ---------------
# main
# ---------------
def main( ):
	args = get_args()
		
	all_info = []
	total_kmers = 0
	dict_counts = {}
	for astrline in open( args.abundances ): # get detection counts /abundances/*_extended file
		perc, count, taxid, taxon, ctype = astrline.strip().split('\t')
		genus, species = find_genus_sp( taxon )
		dict_counts = add_to_dict( dict_counts, genus, species, "exact", "", count )
		total_kmers += float( count )

	if args.hamming: # if hamming exists
		for bstrline in open( args.hamming ): # get hamming matches /hamming/*_hamming file
			taxon, avgHD, count = bstrline.strip().split('\t')
			genus, species = find_genus_sp( taxon )
			dict_counts = add_to_dict( dict_counts, genus, species, "hamming", "", count )
			total_kmers += float( count )

	if args.character: #if characterization exists
		for cstrline in open( args.character ): # get characterization matches /characterization/*_extended file
			firstelement = cstrline.strip().split('\t')[0]
			if re.search( "[ATGC]+", firstelement ):
				kmer, taxon, ptype, chardetail, count = cstrline.strip().split('\t')
				genus, species = find_genus_sp( taxon )
				dict_counts = add_to_dict( dict_counts, genus, species, "characterization", ptype, count )
			
	mismatched = 0 # remainder will be unclassified
	for totals in open( args.total ): # we know 1 and 2 argv definitely exist
		sample, total = totals.strip().split('\t')
		if re.search( sample, args.abundances ):
			mismatched = float( total ) - total_kmers
			print( '\t'.join( str(x) for x in [mismatched] ) )

	# calculate final values and write out file
	for genus in dict_counts:
		for species in dict_counts[genus]:
			char, exact, ham = 0, dict_counts[genus][species]["exact"], 0
			# adjust if characterization was included
			if "characterization" in dict_counts[genus][species]:
				char = 0
				for ptype in dict_counts[genus][species]["characterization"]:
					char_count = sum( dict_counts[genus][species]["characterization"][ptype] )
					print( '\t'.join( str(x) for x in [char_count, genus, species, "characterization", ptype] ) )
					char += char_count
				exact = exact - char
				print( '\t'.join( str(x) for x in [exact, genus, species, "exact"] ) )
			else:
				print( '\t'.join( str(x) for x in [exact,genus,species,"exact"] ) )

			if "hamming" in dict_counts[genus][species]:
				ham = dict_counts[genus][species]["hamming"]
				print( '\t'.join( str(x) for x in [ham, genus, species, "hamming"] ) )
		

if __name__ == "__main__":
	main( )
