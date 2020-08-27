#!/usr/bin/perl

use strict;

if (! @ARGV ||
    $ARGV[0] eq "-h") {
    print "
This program takes nucmer output and calculates per-base coverage values for the reference sequence.

USAGE: $0 nucmer.coords > output.tsv

nucmer.coords - output from searching a genome sequence against itself
                using nucmer (https://doi.org/10.1186/gb-2004-5-2-r12) 
                with the maxmatch option, and then
                converting the delta file to coordinates using the
                show-coords script with the options '-lrHT'
                IMPORTANT: This script screens the identity alignment (ie the
                whole molecule aligning to itself) and any alignment
                that is <97% identical.

output.tsv    - output is a tab-delimited file in the format:
                accession <tab> position <tab> repeat depth
";
    exit();
}
my $file = $ARGV[0];
open (my $in, $file) or die "Can't open input file $file: $!\n";

my %data;
# go through each hit and increment at each genome position within the alignment region
while (my $line = <$in>) {
    if ($line =~ /^#/) { next }  # skip comments
    chomp $line;
    my @f = split/\t/, $line;
    # fields
    # 0 ref low
    # 1 ref high
    # 2 qry low
    # 3 qry high
    # 4 ref match region length
    # 5 qry match region length
    # 6 %id
    # 7 length of ref
    # 8 length of qry
    # 9 ref acc
    # 10 qry acc
    if ($f[6] < 97) {
	print STDERR "Skip low %id region (<97%)\n\t$line\n";
    } else {
	# tally depth across alignment region
	for (my $i = $f[0]; $i <= $f[1]; $i++) {
	    $data{$f[9]}->[$i]++;
	}
    }
}
# foreach sequence, print out the coverage at each position
print "#accession\tposition\trepeat_depth\n";
foreach my $acc (keys %data) {
    my $depth = $data{$acc};
    for (my $pos = 1; $pos<@$depth; $pos++) {
	if (!($depth->[$pos])) { $depth->[$pos] = 0 }
	print "$acc\t$pos\t$depth->[$pos]\n";
    }
}
