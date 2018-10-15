This README was written by Tiffany Hsu on 10/15/2018.

In order to run this script and create Krona visualizations:
```
# Assume that the kmer pipeline has been run and is located at /path/to/results

# Run the python script to generate the xml file
$ generate_xml_krona.py -r /path/to/ref_coords3.tsv -a /path/to/results/abundances -ham /path/to/results/ham_out -c /path/to/results/characterization -t /path/to/results/total_kmers > /path/to/results/myxml.xml

# Load the krona module and create the html file
$ module load krona/2.7
$ ktImportXML /path/to/results/myxml.xml -o xml.krona.html
```

Note that:
* ref_coords.tsv may change as other references are added
* ham_out may be empty, and is not necessary to include
* characterization does not always exist: it only exists if something was found
