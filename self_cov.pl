#!/usr/bin/perl

use strict;

if (! @ARGV ||
    $ARGV[0] eq "-h") {
    print "INTRAGENOME REPETITIVENESS - this program takes nucmer output and identifies repetitive elements

USAGE: $0 self.coords > output.tsv

self.coords - output from searching a genome sequence against itself
              using nucmer (doi://) with the maxmatch option, and then
              converting the delta file to coordinates using the
              show-coords script with the options '-lrHT'
              IMPORTANT: This script screens the identity alignment (ie the
              whole molecule aligning to itself) and any alignment
              that is <97% identical.

output.tsv - output is a tab-delimited file in the format:
             accession <tab> position <tab> repeat depth
";
    exit();
}
my $file = $ARGV[0];
open (my $in, $file) or die "Can't open input file $file: $!\n";

my %data;
# first we go through the file to increment at each position matching each alignment region
# There is no screen of results here, so that must be done (if at all) outside of this program.
while (my $line = <$in>) {
    if ($line =~ /^#/) { next }  # skip comments
    chomp $line;
    my @f = split/\t/, $line;
    if (($f[0] == 1 && $f[1] == $f[7]) || $f[6] < 97) { print STDERR "Skipping $line\n"; }
    # tally depth across alignment region
    for (my $i = $f[0]; $i <= $f[1]; $i++) {
	$data{$f[9]}->[$i]++;
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
