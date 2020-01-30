## Biases in Genome Reconstruction from Metagenomic Data
Contains the code supporting the manuscript Biases in Genome Reconstruction from Metagenomic Data

Doi: TBA


## Variation-Parsing.py for analyzing genome datasets
Will perform analysis on a set of genomes to determine the variation in G+C% and tetranucleotide frequency and the percentage of a genome in repeat regions as determined by NUCmer

### Variation-Parsing Dependencies
Python dependencies
* `seaborn`
* `matplotlib`
* `pandas`
* `scipy`
* `biopython`
* `argparse`
* `numpy`

### Variation-Parsing Install Requiremenets
Tested using the following programs:
* NUCmer (NUCleotide MUMmer) version 3.1
* show-coords

### Variation-Parsing Usage
#### Set-up
* Create a file containing a single column list of genomes to be analyzed 
  * Genomes should be in nucleotide FASTA format with a shared suffix
  * The suffix should not be included in the genome name as part of the list
* Multiple sets of genomes can be run simultaneously, each in a separate file
* If multiple genome sets are provided a boxplot figure will be generated

#### Usage
* `-i` flag designates the name(s) of the files that contain genome lists
* `-f` flag sets the suffix of the FASTA files
* Each analysis should be run separately due to an issue in boxplot generation
* Automated checkpoint files are generated to make it possible to make multiple iterations of boxplot figures quickly without performing a complete re-analysis which can be time intensive
* To re-run any part of the analysis (or set of genomes separately), be sure to remove the automatedly generated outputs that end in the suffixes:
  * `*.tetravar` for tetranucleotide frequency variance
  * `*.gcvar` for G+C% variance
  * `*.coords` for NUCmer repeat regions
```
Example
`python3 Variation-Parsing -i HotLakeMAGs,HotLakeGenomes,RefSeqGenomes,TullyTaraMAGs -f fna --runGC`
```
```
Once complete, removing file names from the -i parameter will quickly regenerate the boxplot figure with the desired genome set
`python3 Variation-Parsing -i RefSeqGenomes,TullyTaraMAGs -f fna --runGC`
```


