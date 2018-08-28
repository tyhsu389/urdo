#!/usr/bin/python

import argparse
import os, re
import lxml.etree

# ---------------
# cli
# ---------------
def get_args( ):
	""" Get arguments passed to script """
	parser = argparse.ArgumentParser(
		formatter_class = argparse.ArgumentDefaultsHelpFormatter
		)	
	parser.add_argument(
		"-d", "--refcoords",
		default="/scicomp/home/s0th/shatavia/URDO_Kmer/Databases/refCoords3.csv",
		help="The reference coordinates used to generate the database."
		)
	parser.add_argument(
		"-a", "--abundances_extended",
		help="Directory with extended abundances (before hamming)."
		)
	parser.add_argument(
		"-ham", "--hamming",
		help="Directory with hamming abundances."
		)
	parser.add_argument(
		"-c", "--characterization",
		help="Directory with characterization abundances."
		)
	parser.add_argument(
		"-t", "--total",
		help="Total number of kmers (not including duplicates."
		)
	args = parser.parse_args()
	return args

def parse_coords( refcoords ):
	""" Create a dictionary in which: { genus: {species: characterized} } from Devika's database """
	dict_hierarchy = {}
	pan_set, sp_set = set([]), set([])
	for astrline in open( refcoords ):
		aastrline = astrline.strip().split(',')
		refname, genus, species, pcoords, dcoords, ccoords = aastrline #defined by Devika's script
		fullname = genus.strip() + ' ' + species.strip()
		
		pname, dname, cname = genus, fullname, fullname	
		if pcoords == '': #if no pan-coordinates
			pname = "NA"
		if dcoords == '':
			dname = "NA"
		if ccoords == '':
			cname = "NA"
		pan_set.add( pname )
		sp_set.add( fullname )

		if genus in dict_hierarchy:
			dict_hierarchy[genus].setdefault( dname, set([]) ).add( cname )
		else:
			dict_hierarchy[genus] = {dname: set([cname])}
	return dict_hierarchy, pan_set, sp_set

def classify_taxon( taxon, dict_hierarchy ):
	""" Get genus and species names from database """
	genus, species = "", ""
	if re.search( "Pan ", taxon ):
		genus = taxon.replace( "Pan ", "" )
		species = "pan"
	else:
		genus, species = taxon.rsplit( ' ', 1 )
		species = taxon
		# below is to make naming same as Devika's database (db)
		# unfortunately this was inconsistent
		if genus in dict_hierarchy: #genera in db is the split taxon name
			genus = genus
			for sp in dict_hierarchy[genus]:
				if re.search( species, sp ):
					species = sp
		elif taxon in dict_hierarchy: #genera in db is the full taxon name
			genus = taxon
			for sp in dict_hierarchy[taxon]:
				if re.search( species, sp ):
					species = sp
		else: #genera in db is some section of taxon name, but species name is full taxon name
			found = False
			for genera in dict_hierarchy:
				if re.search( genera, taxon ):
					for sp in dict_hierarchy[genera]:
						if sp == taxon:
							genus = genera
							species = sp
							found = True
			if found == False:
				print( 'Never found species, check database used' )
	return genus, species

def store_counts( dict_samplecounts, sample, genus, species, counttype, ptype, counts ):
	""" Get counts for ptype {sample: {genus: {species: {counttype: {chartype: counts}}}} """
	# define second to last elemeent
	if counttype == 'characterization':
		lastelement = {ptype: [float( counts )]}
	else:
		lastelement = [float( counts )]

	if sample in dict_samplecounts:
		if genus in dict_samplecounts[sample]:
			if species in dict_samplecounts[sample][genus]:
				if counttype == 'characterization':
					if counttype in dict_samplecounts[sample][genus][species]:
						dict_samplecounts[sample][genus][species][counttype].setdefault( ptype, [] ).append( float( counts ) )
					else:
						dict_samplecounts[sample][genus][species][counttype] = {ptype: [float(counts)]}
				else:			
					dict_samplecounts[sample][genus][species].setdefault( counttype, [] ).append( float(counts) )
			else:
				dict_samplecounts[sample][genus][species] = {counttype: lastelement}
		else:
			dict_samplecounts[sample][genus] = {species: {counttype: lastelement}}
	else:
		dict_samplecounts[sample] = {genus: {species: {counttype: lastelement}}}
	return dict_samplecounts

def record_counts( dict_samplecounts, g_sample_set, sp_sample_set, sample, taxon, counttype, ptype, counts, dict_hierarchy ):
	""" Record taxa counted """
	genus, species = classify_taxon( taxon, dict_hierarchy )
	g_sample_set.add( genus )
	sp_sample_set.add( species )
	dict_samplecounts = store_counts( dict_samplecounts, sample, genus, species, counttype, ptype, counts )
	return dict_samplecounts, g_sample_set, sp_sample_set

def check_values(dict_samplecounts, sample, otherlevels ):
	kmer_count = 0
	if len( otherlevels ) == 1: #genus only
		genus = otherlevels[0]
		if genus in dict_samplecounts[sample]:
			kmer_count = dict_samplecounts[sample][genus]['value']
	elif len( otherlevels ) == 2: #genus, species
		genus, species = otherlevels
		if genus in dict_samplecounts[sample]:
			if species in dict_samplecounts[sample][genus]:
				kmer_count = dict_samplecounts[sample][genus][species]['value']
	elif len( otherlevels ) == 3: #genus, species, ctype
		genus, species, ctype = otherlevels
		if genus in dict_samplecounts[sample]:
			if species in dict_samplecounts[sample][genus]:
				if ctype in dict_samplecounts[sample][genus][species]:
					if ctype == 'characterization':
						kmer_count = dict_samplecounts[sample][genus][species][ctype]['value']
					else:
						kmer_count = sum(dict_samplecounts[sample][genus][species][ctype])
	elif len( otherlevels ) == 4:
		genus, species, ctype, ptype = otherlevels
		if genus in dict_samplecounts[sample]:
			if species in dict_samplecounts[sample][genus]:
				if ctype == 'characterization':
					if ptype in dict_samplecounts[sample][genus][species][ctype]:
						kmer_count = sum( dict_samplecounts[sample][genus][species][ctype][ptype] )
	return kmer_count

def add_to_xml( xml_prevnodes, dict_samplecounts, dict_totalcounts, samplelist, otherlevels ):
	newname = otherlevels[-1]
	xml_nodes = lxml.etree.SubElement( xml_prevnodes, 'node', name=newname )
	xml_count = lxml.etree.SubElement( xml_nodes, 'kmer_count' )
	# get value for every sample
	for sample in samplelist:
		# check value, if present, get kmer_count
		kmer_count = check_values( dict_samplecounts, sample, otherlevels )
		if newname == 'unclassified':
			kmer_count = dict_totalcounts[sample] - dict_samplecounts[sample]['value']
		kmer_val = str( int( kmer_count ) )
		xml_kval = lxml.etree.SubElement( xml_count, "val" )
		xml_kval.text = kmer_val
	return xml_nodes

# ---------------
# main
# ---------------
def main( ):
	args = get_args()

	# get all potential taxa
	dict_hierarchy, pan_set, sp_set = parse_coords( args.refcoords )

	# loop through and get sample information
	dict_samplecounts = {} #{sample1: {genus: {species: {'exact': count, 'hamming': count, 'char': count}...}}
	g_sample_set, sp_sample_set, char_sample_set = set([]), set([]), set([])
	counttype = "exact"
	for afile in os.listdir( args.abundances_extended ):
		if re.search( "_extended", afile ):
			sample_abbr = afile.replace( "_extended", "")
			for astrline in open( args.abundances_extended + '/' + afile ):
				aastrline = astrline.strip().split('\t')
				abundance, counts, taxid, taxon, chartype = aastrline
				dict_samplecounts, pan_sample_set, sp_sample_set = record_counts( dict_samplecounts, g_sample_set, sp_sample_set, sample_abbr, taxon, counttype, None, counts, dict_hierarchy )
	
	# loop through and get hamming information
	counttype = "hamming"
	if args.hamming: #may not have a folder/results
		for bfile in os.listdir( args.hamming ):
			sample_abbr = bfile.replace( "_hamming", "" )
			for bstrline in open( args.hamming + '/' + bfile ):
				bbstrline = bstrline.strip().split('\t')
				taxon, avgham, counts = bbstrline
				dict_samplecounts, pan_sample_set, sp_sample_set = record_counts( dict_samplecounts, g_sample_set, sp_sample_set, sample_abbr, taxon, counttype, None, counts, dict_hierarchy )

	# loop through and get characterization information
	counttype = "characterization"
	if args.characterization: #may not have a folder/results
		for cfile in os.listdir( args.characterization ):
			if re.search( "_extended", cfile ):
				sample_abbr = cfile.replace( "_extended", "" )
				i = 0
				for cstrline in open( args.characterization + "/" + cfile ):
					ccstrline = cstrline.strip().split('\t')
					if i == 0: #first line, get type
						pass #ignoring for now
					else:
						kmer, taxon, ptype, character, counts = ccstrline
						dict_samplecounts, pan_sample_set, sp_sample_set = record_counts( dict_samplecounts, g_sample_set, sp_sample_set, sample_abbr, taxon, counttype, ptype, counts, dict_hierarchy )
					i += 1

	# loop through and get uncounted information
	dict_totalcounts = {}
	for dstrline in open( args.total ):
		sample, counts = dstrline.strip().split('\t')
		dict_totalcounts[ sample ] = float( counts )

	# break up dictionary for easier summing
	# note - need to throw in the "unclassified here at the top level
	ptype_set = set([])
	for sample in dict_samplecounts:
		sample_sum = 0
		for genus in dict_samplecounts[sample]:
			genus_sum = 0
			for species in dict_samplecounts[sample][genus]:
				sp_sum = 0
				for ctype in dict_samplecounts[sample][genus][species]:
					if ctype == 'characterization':
						ctype_sum = 0
						for ptype in dict_samplecounts[sample][genus][species][ctype]:
							ctype_sum += sum( dict_samplecounts[sample][genus][species][ctype][ptype] )
							ptype_set.add( ptype )
						dict_samplecounts[sample][genus][species][ctype]['value'] = ctype_sum
						sp_sum += ctype_sum
					else:
						kmer_count = sum( dict_samplecounts[sample][genus][species][ctype] )
						sp_sum += kmer_count
				dict_samplecounts[sample][genus][species]['value'] = sp_sum
				genus_sum += sp_sum
			dict_samplecounts[sample][genus]['value'] = genus_sum
			sample_sum += genus_sum
		dict_samplecounts[sample]['value'] = sample_sum


	# build xml file
	root = lxml.etree.Element( 'krona' )
	# get attributes
	attributes = lxml.etree.SubElement( root, 'attributes', magnitude="kmer_count" )
	k_attr = lxml.etree.SubElement( attributes, 'attribute', display="kmer_count" )
	k_attr.text = "kmer_count"
	# get colors

	# get samples and nodes (total values )
	xml_samples = lxml.etree.SubElement( root, 'datasets' )
	xml_nodes = lxml.etree.SubElement( root, 'node', name="genus" )
	xml_skcount = lxml.etree.SubElement( xml_nodes, "kmer_count")
	samplelist = list( dict_samplecounts.keys() )
	for sample in samplelist:
		# add samples to list
		xml_ds = lxml.etree.SubElement( xml_samples, 'dataset' )
		xml_ds.text = sample
		# add total sample values to nodes
		kmer_val = str( int( dict_totalcounts[sample] ) )
		xml_skval = lxml.etree.SubElement( xml_skcount, "val" )
		xml_skval.text = kmer_val
	
	# get remaining nodes 
	## genus values
	g_sample_set.add( "unclassified" )
	for genus in g_sample_set:
		xml_gnodes = add_to_xml( xml_nodes, dict_samplecounts, dict_totalcounts, samplelist, [genus] )
		
		# species values
		if genus != "unclassified":
			# look at all species within this genus + samples, including pan species
			species = set( list( dict_hierarchy[genus].keys() ) ) & sp_sample_set
			species.add( 'pan' ) #just in case
			for p_sp in species: 
				xml_spnodes = add_to_xml( xml_gnodes, dict_samplecounts, dict_totalcounts, samplelist, [genus, p_sp] )
							
				# exact, hamming, characterization values
				# look at all types within this species
				for counttype in ["exact", "hamming", "characterization"]:
					xml_cnodes = add_to_xml( xml_spnodes, dict_samplecounts, dict_totalcounts, samplelist, [genus, p_sp, counttype] )

					# phenotypes within characterization values
					if counttype == 'characterization':
						for ptype in ptype_set:
							xml_pnodes = add_to_xml( xml_cnodes, dict_samplecounts, dict_totalcounts, samplelist, [genus, p_sp, counttype, ptype] )

	# print the xml file
	print( lxml.etree.tostring( root, pretty_print=True) )

if __name__ == "__main__":
	main( )
