# URDO k-mer pipeline visualization

## Introduction  
The URDO project currently has multiple bioinformatics pipelines in development 
for pathogen detection and characterization.  This visualization script was 
specifically written for one of these pipelines, which matches and counts 
k-mers from the multiplexed PCR assay to a reference pathogen k-mer database.

The visualization script takes the k-mer counts (as input) and outputs an xml 
file understood by [Krona](https://github.com/marbl/Krona/wiki).

## Input
We take 5 files/directories from the bioinformatic pipeline, which include:
1. `ref_coords.tsv` (file) - This file contains annotations to generate the 
reference pathogen k-mer database.
2. `abundances` (directory) - This directory contains k-mer counts for 
pathogens with exact matches to the database for each sample. 
3. `ham_out` (directory) - This directory contains k-mer counts for pathogens 
with inexact matches to the database for each sample. 
4. `characterization` (directory) - This directory contains k-mer counts for 
pathogen characteristics for each sample.
5. `total_kmers` (file) - This file contains the total k-mer count for each 
sample. 
Longer descriptions of the files are found below.

## Running the script
Assuming that the pipeline has been run and is located at `path/to/results`, 
we can run the Python script by: 
```
# Run the python script to generate the xml file
$ generate_xml_krona.py -r /path/to/results/ref_coords.tsv \
-a /path/to/results/abundances \
-ham /path/to/results/ham_out \
-c /path/to/results/characterization \
-t /path/to/results/total_kmers > /path/to/results/myxml.xml
```
To run the demo, replace `path/to/results` with `path/to/demo`. 

## Creating the Krona visualization
```
# Load the krona module and create the html file
$ module load krona/2.7
$ ktImportXML /path/to/results/myxml.xml -o xml.krona.html
```

## Detailed file descriptions  
1. `ref_coords.tsv`
  a. This file is comma-delimited, with 6 fields in the order of:
1) Name of fasta file containing genome or amplicon region, 2) Genus, 
3) Species, 4) Sequence coordinates for genus level taxonomic annotation, or 
'WS' for "whole sequence", 5) Sequence coordinates for species level taxonomic 
annotation, and 6) Sequence coordinates for trait characterization (such as 
antibiotic resistance.
  b. This file changes as other panel references are added.
2. `abundances` (directory) - This directory contains k-mer counts for
pathogens with exact matches to the database for each sample.
  a. This folder contains files with the format of `samplenumber_extended`
  b. Each file is tab-delimited, with 5 fields in the order of: 1) Percentage 
of taxon (based on kmer-count), 2) K-mer count, 3) Taxon ID, 4) Taxon 
annotation, and 5) Characterization type.
3. `ham_out` - This directory contains k-mer counts for pathogens
with inexact matches to the database for each sample.
  a. This folder contains files named `*_hamming`, and does not exist as 
an output if no inexact matches are found.
  b. Each file is tab-delimited, with 3 fields in the order of 1) Taxon 
annotation, 2) Inexact match distance, 3) K-mer count
4. `characterization`
  a. This folder contains files named `*_extended`.
  b. Each file is tab-delimited, with 5 fields in the order of 1) K-mer 
sequence (defining the trait), 2) Taxon annotation, 3) Phenotype, 
4) Detailed characterization of phenotype, 5) K-mer count
5. `total_kmers`
  a. This file contains the total number of k-mers per sample.
