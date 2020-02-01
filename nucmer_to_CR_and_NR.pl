#!/usr/bin/perl
use strict;

if (! @ARGV ||
    $ARGV[0] eq "-h") {
    print "IDENTIFY CR and NR REGIONS - this program takes nucmer output and identifies the correctly binned regions (CR) and not-binned regions (NR) of a genome

USAGE: $0 reference_vs_MAG.coords

reference_vs_MAG.coords - output from searching a genome sequence against
              a cognate MAG using nucmer (doi://) with the maxmatch option,
              and then converting the delta file to coordinates using the
              show-coords script with the options '-lrHT'
              IMPORTANT: This script uses a cutoff of 99% id and 200 bp 
              for alignment regions.

OUTPUT:       Two files, one for CRs, one for NR with the format:
             genome_accession <tab> low_coord <tab> high_coord <tab> region_length <tab> MAG_accession
";
    exit();
}

# open the input and output files
my $infile = $ARGV[0];
open (IN, $infile) or die "No file $infile? : $!\n";
my $CRfile = $infile . ".CR.coords";
open (CR, ">$CRfile") or die "No write $CRfile? : $!\n";
my $NRfile = $infile . ".NR.coords";
open (NR, ">$NRfile") or die "No write $NRfile? : $!\n";

# step through input
my ($this_lo, $this_hi, $this_qry, $last_hi, $last_ref, $last_reflen);
while (my $l = <IN>) {
    chomp $l;
    my @f = split/\t/, $l;
    
    # check to see if the alignment is good - spurious short repeat regions (99-150 bp) can be found in IDBA_ud assemblies
    if ($f[6] > 99 && $f[5] > 200) {
	# sometimes a MAG contig hit is split; account for this
	if ($this_qry eq $f[10]) {
	    $this_hi = $f[1] if ($f[1] > $last_hi);
	} else {
	    # print out the CR
	    printf CR "%s\t%i\t%i\t%i\t%s\n", ($last_ref, $this_lo, $this_hi, ($this_hi - ($this_lo - 1)), $this_qry) if ($this_lo > 0);
	    $this_lo = $f[0];
	    $this_hi = $f[1];
	    $this_qry = $f[10];
	}

	# Now to sort out the NR
	# if we have switched reference molecules since the last entry, print out
	if ($f[9] ne $last_ref) {
	    # need to make a NR that goes to the end of the molecule
	    # if the last hi coord on the reference was not the molecule length
	    if ($last_reflen > $last_hi) {
		printf NR "%s\t%i\t%i\t%i\n", ($last_ref, $last_hi+1, $last_reflen, $last_reflen-$last_hi);
	    }
	    # set the last hi to 0 for the new molecule
	    $last_hi=0;
	}

	# If the lo coord is higher than the last hi
	# print a NR region if there is a gap between alignments
	if ($f[0] > $last_hi + 1) { 
	    printf NR "%s\t%i\t%i\t%i\n", ($f[9], $last_hi+1, $f[0]-1,($f[0]-1)-$last_hi);
	}

	# sometimes a short alignment region will not extend past the previous alignement region
	if ($f[1] >= $last_hi) { 
	    $last_hi = $f[1] 
	}
    }
    $last_reflen = $f[7];
    $last_ref = $f[9];
}

# print out the final entry
if ($last_reflen > $last_hi) {
    printf NR "%s\t%i\t%i\t%i\n", ($last_ref, $last_hi+1, $last_reflen, $last_reflen-$last_hi);
}
