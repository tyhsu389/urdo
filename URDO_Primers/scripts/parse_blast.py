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

def filter_pcr_product( hit, hit2, primer_stats, dict_primerlen, p1, p2 ):
	k, status = 0, True
	Fprimer_len, Rprimer_len = len( dict_primerlen[p1] ), len( dict_primerlen[p2] )
	for i, j in zip( hit.split('\t')[3:7], hit2.split('\t')[3:7] ):
		if k == 0: #%ID
			#pass
			if float( i ) < 90 or float( j ) < 90: 
				status = False
		elif k == 1: #length
			#pass
			if float( i ) < Fprimer_len*(0.90) or float( j ) < Rprimer_len*(0.90):
				status = False
		elif k== 2: #mismatch
			#pass
			if float( i ) > Fprimer_len*(0.10) or float( j ) > Rprimer_len*(0.10):
				status = False
		else: #gapopen	
			#pass
			if float( i ) > Fprimer_len*(0.10) or float( j ) > Rprimer_len*(0.10):
				status = False
		k += 1
		primer_stats += [i, j]
	return primer_stats, status

def find_pcr_product( taxa, hit, hit2, hitorient, sorient, coordlist, allrefs, dict_primerlen, p1, p2 ):
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
	# This is where you could filter pcr products
	primer_stats = [ len( dict_primerlen[p1] ), len( dict_primerlen[p2] ) ]
	primer_stats, filter_status = filter_pcr_product( hit, hit2, primer_stats, dict_primerlen, p1, p2 )
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
	return results, filter_status

def find_type( taxon, dict_annotation ):
	""" Find the taxon annotation """
	list_types = []
	print( taxon )
	subject = taxon.split('|')[-2].split('.')[0]
	if subject in dict_annotation:
		list_types.append( dict_annotation[subject] )
	else:
		list_types.append( subject )
	return list_types 

def identify_type( annotate, p1_types, p2_types, all_pairs ):
	p1_type, p2_type, pair_type = 'NA', 'NA', 'NA'
	p1_annot, p2_annot, pair_annot = "NA", "NA", "NA"
	if annotate:
		p1_count, p2_count, pair_count = Counter( p1_types ), Counter( p2_types ), Counter( all_pairs )
		p1_type = len( p1_count.keys() )
		p2_type = len( p2_count.keys() )
		pair_type = len( pair_count.keys() )
		p1_annot = ",".join( [':'.join( str(y) for y in [x, p1_count[x]] ) for x in p1_count] )
		p2_annot = ",".join( [':'.join( str(y) for y in [x, p2_count[x]] ) for x in p2_count] )
		pair_annot = ",".join( [':'.join( str(y) for y in [x, pair_count[x]] ) for x in pair_count] )
	return p1_type, p2_type, pair_type, p1_annot, p2_annot, pair_annot 

def primer_counts( p1, dict_results, dict_annotation ):
	""" Return the number of unique taxa and hits per primer """
	p1_taxa, p1_hits, p1_types = 0, 0, []
	if p1 in dict_results.keys():
		p1_taxa = len( dict_results[p1].keys() ) # number of taxon hits
		for taxon in dict_results[p1]:
			p1_hits += len( dict_results[p1][taxon] )
			p1_types += find_type( taxon, dict_annotation )
	return p1_taxa, p1_hits, p1_types

def primer_status( p1, p2, dict_results ):
	""" Determine whether both/one/none of the F/R primers have hits """
	status = ""
	if( p1 in dict_results.keys() ) and( p2 not in dict_results.keys() ):
		status = "F hits, R no_hits"
	elif( p1 not in dict_results.keys() ) and( p2 in dict_results.keys() ):
		status = "F no_hits, R hits"
	elif( p1 not in dict_results.keys() ) and( p2 not in dict_results.keys() ):
		status = "F no_hits, R no_hits"
	elif( p1 in dict_results.keys() ) and( p2 in dict_results.keys() ):
		status = "F hits, R hits"
	return status

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

	# Define the primer pairs (keep all in F/R orientation)
	dict_primers = {} #{ forward_primer: reverse_primer }
	dict_primerlen = {} # {primer: sequence}
	for seqrec in SeqIO.parse( args.primerfasta, "fasta" ):
		dict_primerlen[ seqrec.id ] = seqrec.seq
		if re.search( "_F$", seqrec.id ): #forward primer
			reverse = seqrec.id.replace( "_F", "_R" )
			dict_primers[ seqrec.id ] = reverse
		elif re.search( "_R$", seqrec.id ): #reverse primer
			forward = seqrec.id.replace("_R", "_F" )
			dict_primers[ forward ] = seqrec.id

	# Get taxids for annotation if relevant
	# Eventually change this to a two column file
	dict_annotate = {}
	if args.annotate:
		for astrline in open( args.taxids ):
			aastrline = astrline.strip().split('\t')
			if not re.search( '^#', astrline ):
				NC_value, rep_value = aastrline[0], aastrline[1]
				taxonomy = aastrline[3]
				dict_annotate[NC_value] = taxonomy
				dict_annotate[rep_value] = taxonomy

	# Loop through BLAST results and store information
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

	# Output for summary file
	fh_summary = open( args.output + ".primer_summary", "w" )

	# For each primer pair, get:
	#1) singleton hits
	#2) primer pair products
	#3) status
	# Loop through forward primers
	for p1 in dict_primers:
		p2 = dict_primers[p1] #Get reverse primer
		primername = p1 + '_' + p2

		# Generate output for each primer
		fh_products = open( args.output + primername + ".products", "w" )
		# Get status of primer
		status = primer_status( p1, p2, dict_results )
	
		# Get information about singletons
		p1_taxa, p1_hits, p1_types = primer_counts( p1, dict_results, dict_annotate )
		p2_taxa, p2_hits, p2_types = primer_counts( p2, dict_results, dict_annotate )
	
		# Get information about potential products
		# 1) Looking for similar organisms they hit
		# 2) Checking primers are facing each other
		list_records_per_pair, all_pairs = [], []
		if p1 in dict_results.keys() and p2 in dict_results.keys():
			taxap1 = set( dict_results[p1].keys() )
			taxap2 = set( dict_results[p2].keys() )
			both_taxa = taxap1 & taxap2 # shared hits which are possible products

			# For subjects/taxa both primers share, loop through unique hits
			for taxa in both_taxa: # each taxon they share
				for hit in dict_results[p1][taxa]: # hits for F primer 
					fstart, fstop = int( hit.split('\t')[9] ), int( hit.split('\t')[10] )
					forient = orient( fstart, fstop )
					for hit2 in dict_results[p2][taxa]: # hits for R primer
						rstart, rstop = int( hit2.split('\t')[9] ), int( hit2.split('\t')[10] )
						rorient = orient( rstart, rstop )
						
						# Check orientation and get product
						face_status = facing_primers( fstart, fstop, rstart, rstop, forient, rorient )
						if face_status == True: # in opposite directions
							coordlist = [ fstart, fstop, rstart, rstop ]
							potential_product = max( coordlist ) - min( coordlist ) + 1
							results, filter_status = find_pcr_product( taxa,
									hit, hit2,
									forient, rorient,
									coordlist, allrefs, 
									dict_primerlen, p1, p2 )

							# Print and join information for potential products
							if filter_status == True:
								fh_products.write( '\t'.join( str(x) for x in results ) + '\n' ) 
								list_records_per_pair.append( potential_product )
							if args.annotate:
								all_pairs += find_type( taxa, dict_annotate )
	
		# Close products file
		fh_products.close()
	
		# Write summary results for products within each primer pair
		nobs, minmax, mean, variance, median = 0, (0,0), 'NA', 'NA', 0
		if len( list_records_per_pair ) > 0:
			nobs, minmax, mean, variance = list( stats.describe( np.array( list_records_per_pair ) ) )[0:4]
			median = np.median( np.array( list_records_per_pair ) )
		
		# Get types if appropriate
		p1_type, p2_type, pair_type, p1_annot, p2_annot, pair_annot = identify_type( args.annotate, p1_types, p2_types, all_pairs )
		primer_pair_info = [p1 + ";"+ p2, 
			status,
			p1_taxa, #number of unique taxa for F primer
			p1_hits, #number of unique hits for F
			p1_type, #number of unique taxa/type
			p1_annot, #names of unique taxa/type
			p2_taxa, #same as above for R primer
			p2_hits,
			p2_type,
			p2_annot,
			nobs, #number of pcr products
			pair_type,
			pair_annot,
			minmax[0], #pcr product size stats
			mean, 
			median,
			minmax[1],
			variance
			]
		fh_summary.write( "\t".join( str(x) for x in primer_pair_info ) + '\n' )
	fh_summary.close()

if __name__ == "__main__":
	main()
