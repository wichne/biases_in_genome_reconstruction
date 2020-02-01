## Biases in Genome Reconstruction from Metagenomic Data
Contains the code supporting the manuscript Biases in Genome Reconstruction from Metagenomic Data

Doi: TBA

## nucmer_to_CR_and_NR.pl
Parses binned and not-binned regions from MUMmer show-coords formatted output

### nucmer_to_CR_and_NR dependencies
none

### nucmer_to_CR_and_NR Usage
nucmer_to_CR_and_NR.pl reference_vs_MAG.coords

reference_vs_MAG.coords - output from searching a genome sequence against
              a cognate MAG using nucmer (doi://) with the maxmatch option,
              and then converting the delta file to coordinates using the
              show-coords script with the options '-lrHT'
              IMPORTANT: This script uses a cutoff of 99% id and 200 bp
              for alignment regions.

OUTPUT:       Two files, one for CRs, one for NR with the format:
              genome_accession <tab> low_coord <tab> high_coord <tab> region_length <tab> MAG_accession

## gc_analysis.pl
Will perform %G+C analysis and comparison of genome with binned and not-binned regions.

### gc_analysis dependencies
Bio::SeqIO
Statistics::Basic

### gc_analysis Usage
gc_analysis.pl genome.fa CR.coords NR.coords

genome.fa - fasta file containing the genome sequence

CR.coords - a tab-delimited file describing the correctly binned regions of the
            genome. It must contain of the following fields:
            0 sequence accession (must match an accession in the fasta file)
            1 low coord
            2 high coord

NR.coords - a file with the same format as above, but describing the not-binned
            regions of the genome.

## self_cov.pl parses per-base intragenome repetitiveness

### self_cov Dependencies
none

### self_cov Usage
self_cov.pl self.coords > output.tsv

self.coords - output from searching a genome sequence against itself
              using nucmer (doi://) with the maxmatch option, and then
              converting the delta file to coordinates using the
              show-coords script with the options '-lrHT'
              IMPORTANT: This script screens the identity alignment (ie the
              whole molecule aligning to itself) and any alignment
              that is <97% identical.

output.tsv - output is a tab-delimited file in the format:
             accession <tab> position <tab> repeat depth

## genome_redundancy_analysis.pl for analyzing genome and binned and not-binned regions repetitiveness
Parses output from self_cov.pl and compares calculated repetitiveness of CR and NR to genome.

### genome_redundancy_analysis Dependencies
Bio::SeqIO
Statistics::Basic

### genome_redundancy_analysis Usage
genome_redundancy_analysis.pl self_cov.tsv genome_stats.txt CR.coords NR.coords

self_cov.tsv     - output from self_cov.pl

genome_stats.txt - descriptive statistics for each molecule in the genome fasta file
                   with the format:
                   >seq_acc description length=NNNN G+C=NN.NN A=NNNN T=NNNN C=NNNN G=NNNN
                      where the seq_acc should match the accessions in the other files,
                      length in bp, G+C as a percentage, and counts of any symbols found
                      in the sequence (including ambiguous and gap characters).

CR.coords        - a tab-delimited file describing the correctly binned regions of the
                   genome. It must contain of the following fields:
                   0 sequence accession (must match an accession in the fasta file)
                   1 low coord
                   2 high coord
                   !!! If the high coord is lower than the low coord, the program assumes
                   !!! the chromosome is circular and the region spans the origin.

NR.coords        - a file with the same format as above, but describing the not-binned
                   regions of the genome.

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


