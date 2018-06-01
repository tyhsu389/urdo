URDO viral primer pipeline
=============
Introduction
-------------
The primer pipeline should generate primers for viruses.
The pipeline currently works by:

1. Download genomes from NCBI.
2. Orgznie sequences by species/serotype.
3. Align genomes using ClustalO (all, as well as by species/serotype).
3. Convert the alignment into a consensus sequence.
4. Run Primer3 and DEGEPRIME to choose primers.
5. Align each primer set against all organism sequences (from step 1) using EMBOSS's primersearch. Determine:
  a. Primer set specificity to one or multiple species/serotype (species or pan-primer)
  b. Product size distribution
6. For primer sets that pass step 5: 
  a. BLAST each primer against "nr".
  b. BLAST each primer against kmer pipeline reference database.
  c. Perform primer-BLAST on NCBI.
7. Choose only primers that pass steps 5 and 6 for wet laboratory testing.

Related Documents
------------

Prerequisites
------------
1. Python 2.7
2. BioPython 1.7 
3. EMBOSS suite (tested with 6.6.0)
4. clustalo
5. primer3 (tested with v2.4.0)
6. DEGEPRIME
7. NCBI-blastn

Installation
------------
All of the above are currently installed as modules on the CDC servers.

Basic Usage (Procedure)
------------
**Step 1: Download genomes from NCBI.**
1. Download genomes using [Batch Entrez](https://www.ncbi.nlm.nih.gov/sites/batchentrez) at NCBI.
     For example, for hrv-c:  
     - We first visited [NCBI Viral Genomes](https://www.ncbi.nlm.nih.gov/genome/viruses/) and
        downloaded the file under "Accession list of all viral genomes".  
     - We then accessed the downloaded file named `taxid10239.nbr`. We can then 
        "grep" for the virus of interest.  
     - The "grepped" list can be used as the "list of identifiers" for the next step.
2. Input list of identifiers (a text file) under "File" and hit "Retrieve".  
3. After reaching the page with genomes, click on "Send to:" dropdown, and:  
     - Select "Complete Record".  
     - Under "Choose Destination", selection "File." for format select "FASTA"
        and check "Show GI". Then click "Create File."  
     - Save the file, which will automatically be named `sequence.fasta`.    
    
```sh
## All the below are for Step 1.1 ##
# Choose genomes named "rhinovirus C"
$ cat taxid10239.nbr | grep 'Rhinovirus C' > rhinovirusConly.txt

# Add back the header
$ head -n2 taxid10239.nbr | cat - rhinovirusConly.txt > rhinovirusConly_h.txt

# Take the second column
$ cut -f2 rhinovirusConly_h.txt > rhinovirusCrep_list.txt

# Add the reference sequence back by manually editing in vim
$ vi rhinovirusCrep_list.txt
```


**Step 2: Align viral genomes using clustalo.**
1. Check whether genome lengths are approximately the same. You cannot cluster sequences that are very 
different in length.
2. Create fasta files for the groups you would like to cluster. 
3. Cluster via clustalo.  
4. Visualize alignment either in JalView or CLC Bio (whatever you have available).

```sh
# cluster all hrv species together
$ clustalo -i YOURFASTA.fasta -o YOURFASTA.clustalo.fasta

# rename fasta files headers in order to tell between types
$ python annotate_alignment.py -i YOURFASTA.clustalo.fasta -t taxid10239.nbr -o YOURFASTA.clustalo.annotated.fasta
```


**Step 3: Create consensus sequence and run Primer3.**

Step 3.1: Convert into a consensus sequence.
Step 3.2: Create a Primer3 template file.
Step 3.3: Run Primer3.
Primer3 has a multitude of options documented [online](http://primer3.sourceforge.net/primer3_manual.htm).

```sh
## Split the alignment into each species ##

## Convert to consensus sequence ##
# Create the consensus sequence with ambiguous bases using the EMBOSS tool for each alignment.
$ module load EMBOSS/6.5.7 # do if on cluster/BioLinux
$ consambig -sequence YOURFASTA.clustalo.annotated.fasta -outseq YOURFASTA.consambig

# You also have the option to create a consensus sequence without ambiguous bases (only Ns) if you use:
$ cons -sequence YOURFASTA.clustalo.annotated.fasta -outseq YOURFASTA.cons

## Create a Primer3 template file. ##
# Primer3 works by taking in a template of the parameters you want:
$ python create_primer3_template.py -c YOURGENOMESNAME.cons -o YOURGENOMESNAME.CONS.primer3 -p1

# The input goes under the -c flag for the consensus sequence.
# The p1 flag specifies where the Primer3 config is. Use -p2 to specify a different
# path other than the default.

## Generate primers. ##
$ primer3_core < YOURFASTA.cons.primer3 > YOURFASTA.cons.primer3.out

# Combine primer output for each gene region and format for blast and primersearch.
$ get_primers.py -i YOURFASTA.cons.primer3.out -n 4 -b YOURFASTA.blast -e YOURFASTA.emboss

# The .emboss file is all the stats from primer3
# There are 2 .blast files created:
# 1) YOURFASTA.blast 
# 2) single_YOURFASTA.blast.
# The former is the forward primer stitched with the reverse complement, with 
# "Ns" in between, which can be specified with "-n"
# The latter is the two primers separate.
```

**Step 4: Trim alignment and run DEGEPRIME.**

This is currently optional. DEGEPRIME is better than Primer3 
at generating degenerate primers, but currently you need to 
manually search for primers.

Step 4.1: "Trim" the alignment.
Step 4.2: Run DEGEPRIME.

```sh
##Trim the alignment ##
$ perl ~/programs/DEGEPRIME/TrimAlignment.pl -i YOURFASTA.clustalo.annotated.fasta -min 0.5 
-o YOURFASTA.trimmed_align_file
## Create a Primer3 template file. ##

## Run DEGEPRIME
$ perl ~/programs/DEGEPRIME/DegePrime.pl -i YOURFASTA.trimmed_align_file -d 6 -l 20 -o YOURFASTA.output

## Search for primers
#Currently, we sort YOURFASTA.output by # of sequences with a matching primer, and then
#by degeneracy.
```

**Step 5: Check your primers for specificity to your organism. [Optimizing more parameters]**  
```sh
# BLAST against your taxon of interest using "human rhinovirus" as an example.
$ blastn -query YOURGENOMESNAME.CONS.BLAST -out BLASTOUTPUT_TAXSPECIFIC.FMT1.txt 
-db nr -task blastn-short -word_size 7 -max_target_seqs 50000 -reward 1 -penalty -1 
-entrez_query "human rhinovirus" -remote

# BLAST against everything EXCEPT taxon of interest.
$ blastn -query YOURGENOMESNAME.CONS.BLAST -out BLASTOUTPUT_TAXSPECIFIC.FMT1.txt 
-db nr -task blastn-short -word_size 7 -max_target_seqs 50000 -reward 1 -penalty -1 
-entrez_query "NOT human rhinovirus" -remote

# For custom, parse-able output in Step 8, use outfmt 6
$ blastn -query YOURGENOMESNAME.CONS.BLAST -out BLASTOUTPUT_TAXSPECIFIC.FMT6.txt 
-db nr -task blastn-short -word_size 7 -max_target_seqs 50000 -reward 1 -penalty -1 
-entrez_query "human rhinovirus" -remote -outfmt "6 qaccver sallseqid stitle pident 
nident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore"
```

  
**Step 8: Parse BLAST output. [Improve by adding summaries]**  
```sh
$ parse_blast.py -i BLASTOUTPUT_TAXSPECIFIC.FMT6.txt -f YOURGENOMESNAME.CONS.BLAST 
> BLASTOUTPUT_TAXSPECIFIC.PARSED.txt
```


**Step 9: BLAST amplicons sequences against nr. [Yet to be implemented]**  


Advanced Options (Procedure)
------------

Future Directions (Questions)
------------
Current Approach
- Deciding on representative genomes for whole genome alignment
  - Need a way to choose groups
- Why do consensus sequences with Ns generate more primers?
  - Is Primer3 not as good with degenerate sequences or avoid them?
- Deciding on BLAST parameters: 
  - Current parameters based on some of PrimerBLAST's parameters, but have not managed to recreate their results.
  - Need to decide on summary output that would be useful for deciding on primer pairs.
  - Perhaps output products in separate fasta file for BLAST
  - Why do mugsy primers fail in BLAST?
- Installation onto CLC Bio
  - Everything is currently implemented on the CDC Servers
  - Would like to install on CLC, but not sure if we should go for a module system

Other Approaches
- PrimerBLAST has an initial step of megaBLAST-ing your sequence of interest against all other organisms, and masking aligned sections before calling Primer3.
  - It might be possible to BLAST the genome-of-interest consensus sequence against nr and remove any section that has a hit before proceeding with the remainder of the pipeline
  - May be a more "appropriate use of BLAST" that works better
- HYDEN and DEGEPRIME probably generate better degenerate primers, but are much less helpful than Primer3
- With the BLAST script written, may be able to generate product sizes even if we use something like DEGEPRIME and Hyden
