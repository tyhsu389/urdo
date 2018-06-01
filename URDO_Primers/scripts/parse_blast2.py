#!/usr/bin/python

"""
The goal of this script is to parse primer BLAST output.
"""

import argparse, re
import numpy as np
from collections import Counter
from scipy import stats
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Arguments
def get_args():
	"""
	Get arguments passed to script
	"""
	parser = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
		)
	parser.add_argument(
		"-i", "--blast",
		help="BLAST output in outfmt 6",
		)
	parser.add_argument(
		"-pf", "--primerfasta",
		help="FASTA file used to as query to BLAST",
		)
	parser.add_argument(
		"-rf", "--reffasta",
		help="FASTA file used to generate BLAST db",
		)
	parser.add_argument(
		"-t", "--taxids",
		help="Viral taxids from NCBI for annotation",
		)
	parser.add_argument(
		"-o", "--output",
		help="Output fasta with new headers",
		)
	parser.add_argument(
		"-e", "--exclude",
		help="Remove BLAST hits associated with this regex",
		type=str,
		)
	parser.add_argument(
		"-a", "--annotate",
		help="Annotate with types using taxids",
		action="store_true"
		)
	args = parser.parse_args()
	return args

def orient( sstart, send ):
	strand = "-"
	if send - sstart > 0:
		strand = "+"
	return strand

def facing_primers( p1start, p1end, p2start, p2end, p1orient, p2orient ):
	status = False
	if p1orient == "+" and p2orient == "-":
		# p1 could be 1 to 20 (start/end), and p2 is 40 and 21 (start/end)
		if p2end > p1end:
			status = True
	elif p1orient == '-' and p2orient == '+':
		# if p1 is 20 to 1 (start/end), and p2 is 21 and 40 (start/end), these would face away
		# if p1 is 40 to 21 (start/end), and p2 is 1 to 20 (start/end), these would face each other
		if p2start < p1end:
			status = True
	return status

def find_type( taxon, dict_annotation ):
	list_types = []
	subject = taxon.split('|')[-2].split('.')[0]
	if subject in dict_annotation:
		list_types.append( dict_annotation[subject] )
	else:
		list_types.append( subject )
	return list_types 

def type_single_product( p1, dict_results, dict_annotation ):
	list_types = []
	if p1 in dict_results.keys():
		for taxon in dict_results[p1]:
			list_types += find_type( taxon, dict_annotation )
	return Counter( list_types )

def primer_counts( p1, dict_results, def_filter ):
	p1_taxa, p1_filttaxa = 0, 0
	p1_hits, p1_filthits = 0, 0
	if p1 in dict_results.keys():
		p1_taxa = len( dict_results[p1].keys() ) # number of taxon hits
		for taxon in dict_results[p1]:
			p1_hits += len( dict_results[p1][taxon] )
			desc = list( dict_results[p1][taxon] )[0].split('\t')[2]

			# count only those that pass filter
			if def_filter:
				if not re.search( def_filter, desc ):
					p1_filthits += len( dict_results[p1][taxon] )
					p1_filttaxa += 1
	if not def_filter: #if no filter assign same values
		p1_filttaxa = p1_taxa
		p1_filthits = p1_hits
				
	return p1_taxa, p1_hits, p1_filttaxa, p1_filthits

def primer_status( p1, p2, dict_results, def_filter ):
	p1_taxa, p1_hits, p1_filttaxa, p1_filthits = primer_counts( p1, dict_results, def_filter )
	p2_taxa, p2_hits, p2_filttaxa, p2_filthits = primer_counts( p2, dict_results, def_filter )
	status = ""
	if( p1 in dict_results.keys() ) and( p2 not in dict_results.keys() ):
		status = "F hits, R no_hits"
	elif( p1 not in dict_results.keys() ) and( p2 in dict_results.keys() ):
		status = "F no_hits, R hits"
	elif( p1 not in dict_results.keys() ) and( p2 not in dict_results.keys() ):
		status = "F no_hits, R no_hits"
	elif( p1 in dict_results.keys() ) and( p2 in dict_results.keys() ):
		status = "F hits, R hits"
	return p1_taxa, p1_hits, p1_filttaxa, p1_filthits, p2_taxa, p2_hits, p2_filttaxa, p2_filthits, status

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------
def main():
	args = get_args()

	# Get each reference and its sequence/properties
	allrefs = {}
	if args.reffasta:
		for seqrec in SeqIO.parse( args.reffasta, "fasta" ):
			allrefs[ str(seqrec.id) ] = seqrec

	# Primer pairs (define, keep all in F/R orientation)
	dict_primers = {}
	dict_primerlen = {}
	for seqrec in SeqIO.parse( args.primerfasta, "fasta" ):
		dict_primerlen[ seqrec.id ] = seqrec.seq
		if re.search( "_F$", seqrec.id ): #forward primer
			reverse = seqrec.id.replace( "_F", "_R" )
			dict_primers[ seqrec.id ] = reverse
		elif re.search( "_R$", seqrec.id ): #reverse primer
			forward = seqrec.id.replace("_R", "_F" )
			dict_primers[ forward ] = seqrec.id
		else: #weird (should not exist)
			print( "strange" )

	# Get taxids
	dict_annotate = {}
	if args.annotate:
		for astrline in open( args.taxids ):
			aastrline = astrline.strip().split('\t')
			if not re.search( '^#', astrline ):
				NC_value, rep_value = aastrline[0], aastrline[1]
				taxonomy = aastrline[3]
				dict_annotate[NC_value] = taxonomy
				dict_annotate[rep_value] = taxonomy


	# Generate output files
	fh_primers = open( args.output + ".primerstats", "w" )
	fh_primers.write( '\t'.join( ["PRIMER_PAIR","STATUS", "FILTER", 
				"F_TAXA", "F_HITS", "F_TYPE", "F_TYPE_COUNT", "F_TAXA_FILT", "F_HITS_FILT", 
				"R_TAXA", "R_HITS", "R_TYPE", "R_TYPE_COUNT", "R_TAXA_FILT", "R_HITS_FILT",
				"NUM_PRODUCTS", "PAIR_TYPE", "PAIR_TYPE_COUNT", "NUM_PRODUCTS_FIL", "MIN", "MEAN", "MEDIAN", "MAX", "VAR"] ) + '\n' )
	fh_products = open( args.output + ".products", "w" )
	fh_products.write( '\t'.join( ["PRIMER_PAIR", "TAXAID", "DESC", "AMPLEN", "START", "END", "AMPLEN2", 
				"F_LEN", "R_LEN", "F_ID", "R_ID", "F_ALIGNED_LEN", "R_ALIGNED_LEN", 
				"F_MISMATCH", "R_MISMATCH", "F_GAP", "R_GAP", "F_STRAND" ,"R_STRAND", 
				"F_SEQ", "R_SEQ", "AMP_SEQ"] ) + '\n' )
	fh_products_filt = open( args.output + ".products_filt", "w" )
	fh_products_filt.write( '\t'.join( ["PRIMER_PAIR", "TAXAID", "DESC", "AMPLEN", "START", "END", "AMPLEN2",
				"F_LEN", "R_LEN", "F_ID", "R_ID", "F_ALIGNED_LEN", "R_ALIGNED_LEN",
				"F_MISMATCH", "R_MISMATCH", "F_GAP", "R_GAP", "F_STRAND" ,"R_STRAND", 
				"F_SEQ", "R_SEQ", "AMP_SEQ"] ) + '\n' )

	# Loop through BLAST results
	dict_results = {} # create dictionary in which we have dict[query][subject] = set( unique_hits )
	for astrline in open( args.blast ):
		aastrline = astrline.strip().split('\t')
		qaccver, saccver, stitle, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = aastrline
		# for each primer, group the hits by subject
		if qaccver in dict_results:
			if saccver in dict_results[qaccver]:
				dict_results[qaccver][saccver].add( astrline.strip() )
			else:
				dict_results[qaccver][saccver] = set( [astrline.strip()] )
		else:
			dict_results[qaccver] = {saccver: set( [astrline.strip()] )}

	# For each primer pair, get:
	#1) # of subjects hit
	#2) # of unique BLAST hits		
	#3) Status (did both F/R have hits)
	for p1 in dict_primers:
		p2 = dict_primers[p1]

		# Get information about p1 and p2
		p1_taxa, p1_hits, p1_filttaxa, p1_filthits, p2_taxa, p2_hits, p2_filttaxa, p2_filthits, status = primer_status( p1, p2, dict_results, args.exclude )
	
		# Type each primer
		p1_types, p2_types = {}, {}
		if args.annotate:
			p1_types = type_single_product( p1, dict_results, dict_annotate )
			p2_types = type_single_product( p2, dict_results, dict_annotate )
	
		# If both primers have hits, check for potential products
		# Do this by: 
		# 1) Looking for similar organisms they hit
		# 2) Checking primers are facing each other
		list_records_per_pair, list_records_per_pair_filt = [], []
		all_pairs, all_pairs_filt = [], []
		if p1 in dict_results.keys() and p2 in dict_results.keys():
			#print( p1, p2 )
			taxap1 = set( dict_results[p1].keys() )
			taxap2 = set( dict_results[p2].keys() )
			both_taxa = taxap1 & taxap2 # shared hits

			# For subjects/taxa both primers share, loop through unique hits
			for taxa in both_taxa: # each taxon they share
				for hit in dict_results[p1][taxa]: # hits for F primer 
					start, stop = int( hit.split('\t')[9] ), int( hit.split('\t')[10] )
					hitorient = orient( start, stop )
					for hit2 in dict_results[p2][taxa]: # hits for R primer
						#print( hit )
						#print( hit2 )
						sstart, sstop = int( hit2.split('\t')[9] ), int( hit2.split('\t')[10] )
						sorient = orient( sstart, sstop )
						
						# Check orientation
						face_status = facing_primers( start, stop, sstart, sstop, hitorient, sorient )
						if face_status == True: # in opposite directions
							coordlist = [ start, stop, sstart, sstop ]
							potential_product = max( coordlist ) - min( coordlist ) + 1

							# Annotate with references and try to get amplicon
							# Start with default from BLAST db
							ampseq = ""
							taxaid = taxa
							desc = hit.split('\t')[2]
							if len( allrefs ) > 0: # Or use fasta used to generate BLAST db
								ampseq = allrefs[taxa].seq[ min(coordlist)-1: max(coordlist) ]
								taxaid = allrefs[taxa].id
								desc = allrefs[taxa].description

							# Get stats for each primer
							primer_stats = [ len( dict_primerlen[p1] ), len( dict_primerlen[p2] ) ]
							for i, j in zip( hit.split('\t')[3:7], hit2.split('\t')[3:7] ):
								primer_stats += [i, j]
							primer_stats += [hitorient, sorient]

							# Combine information
							info = [p1 + ";" + p2, 
								taxaid, 
								desc,
								len( ampseq ), 
								min(coordlist), 
								max(coordlist), 
								max( coordlist ) - min( coordlist ) + 1
								]
							sequences = [dict_primerlen[p1],
								dict_primerlen[p2],
								ampseq,
								]
							results = info + primer_stats + sequences

							# Implement exclusion filter if needed
							if args.exclude:
								if not re.search( args.exclude, desc ):
									fh_products_filt.write( '\t'.join( str(x) for x in results ) + '\n' )
									list_records_per_pair_filt.append( max( coordlist ) - min( coordlist ) + 1 )
							else:
								list_records_per_pair_filt.append( max( coordlist ) - min( coordlist ) + 1 )
							# Print and join information for potential products
							fh_products.write( '\t'.join( str(x) for x in results ) + '\n' ) 
							list_records_per_pair.append( max( coordlist ) - min( coordlist ) + 1 )
							if args.annotate:
								all_pairs += find_type( taxaid, dict_annotate )
		
		# Write results for each primer pair
		nobs, minmax, mean, variance, median = 0, (0,0), 'NA', 'NA', 0
		if len( list_records_per_pair ) > 0:
			nobs, minmax, mean, variance = list( stats.describe( np.array( list_records_per_pair ) ) )[0:4]
			median = np.median( np.array( list_records_per_pair ) )
		filt_nobs = len( list_records_per_pair_filt )
		
		# Get types if appropriate
		p1_type, p2_type, pair_type = 'NA', 'NA', 'NA'
		p1_annot, p2_annot, pair_annot = "NA", "NA", "NA"
		if args.annotate:
			p1_type, p2_type, pair_type = len( p1_types.keys() ), len( p2_types.keys() ), len( Counter( all_pairs ).keys() )
			p1_annot = ",".join( [':'.join( str(y) for y in [x, p1_types[x]] ) for x in p1_types] )
			p2_annot = ",".join( [':'.join( str(y) for y in [x, p2_types[x]] ) for x in p2_types] )
			pair_annot = ",".join( [':'.join( str(y) for y in [x, Counter( all_pairs )[x]] ) for x in Counter( all_pairs )] )
		primer_pair_info = [p1 + ";"+ p2, 
			status,
			args.exclude,
			p1_taxa, 
			p1_hits,
			p1_type,
			p1_annot, 
			p1_filttaxa, 
			p1_filthits,
			p2_taxa, 
			p2_hits,
			p2_type,
			p2_annot,
			p2_filttaxa, 
			p2_filthits, 
			nobs,
			pair_type,
			pair_annot,
			filt_nobs,
			minmax[0],
			mean,
			median,
			minmax[1],
			variance
			]
		print( primer_pair_info )
		fh_primers.write( "\t".join( str(x) for x in primer_pair_info ) + '\n' )

	fh_primers.close()
	fh_products.close()
	fh_products_filt.close()

if __name__ == "__main__":
	main()
