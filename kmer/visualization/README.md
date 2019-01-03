# URDO k-mer pipeline visualization

## Introduction
The URDO project currently has multiple bioinformatics pipelines in development for pathogen detection and characterization. \
This visualization script was specifically written for one of these pipelines, which matches and counts k-mers from the multiplexed \
PCR assay to a reference pathogen k-mer database. 

The visualization script takes the k-mer counts (as input) and outputs an xml file understood by [Krona](https://github.com/marbl/Krona/wiki).

## Input
We take 5 files/directories from the bioinformatic pipeline, which include:
1. `ref_coords.tsv` (file) - This file contains annotations to generate the reference pathogen k-mer database.
2. `abundances` (directory) - This directory contains k-mer counts for pathogens with exact matches to the database for each sample. 
3. `ham_out` (directory) - This directory contains k-mer counts for pathogens with inexact matches to the database for each sample. 
4. `characterization` (directory) - This directory contains k-mer counts for pathogen characteristics for each sample.
5. `total_kmers` (file) - This file contains the total k-mer count for each sample.  

## Running the script
Assuming that the pipeline has been run and is located at `path/to/results`, we can run the Python script by: 
```
# Run the python script to generate the xml file
$ generate_xml_krona.py -r /path/to/ref_coords3.tsv \
-a /path/to/results/abundances \
-ham /path/to/results/ham_out \
-c /path/to/results/characterization \
-t /path/to/results/total_kmers > /path/to/results/myxml.xml
```

## Creating the Krona visualiztion
```
# Load the krona module and create the html file
$ module load krona/2.7
$ ktImportXML /path/to/results/myxml.xml -o xml.krona.html
```

## Notes
* ref_coords.tsv may change as other references are added
* ham_out may be empty, and is not necessary to include
* characterization does not always exist: it only exists if something was found
