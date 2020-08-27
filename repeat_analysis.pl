#!/usr/bin/perl

use Bio::SeqIO;
use Statistics::Basic qw(:all);
use strict;

$| = 1;

if (! @ARGV ||
    $ARGV[0] eq "-h") {
    print "REPETITIVENESS ANALYSIS

This program determines how much of a sequence is comprised of repeats by parsing search
results to calculate an average per-base coverage. For example, if a 1000bp region
contained a 200 bp segment which was repeated 10 times in the genome, with the rest 
of the region only occuring once in the genome, the repetitiveness score would be:

   ((800 x 1) + (200 x 10)) / 1000 = 2.8

The same situation on a region of 10,000bp would result in a score of 0.28. Thus the 
score reflects the degree of repetitiveness of segments within the sequence region 
and the amount of the sequence region comprised of repeated sequence.

USAGE: $0 cov.tsv genome_stats.txt CR.coords NR.coords

cov.tsv          - a tab-delimited file, such as is produced by nucmer_to_cov.pl,
                   with the following fields:
                   0 sequence accession
                   1 position
                   2 repeat depth


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

";
    exit();
}

# The process:
# 1. calculate the genome coverage for the genome
# 2. calculate the average coverage of the binned contigs (get coords from binned.coords)
# 2a. do 1000 reps of pull random seqs of same lengths
# 2b. calculate distance and p-value
# 3. calculate the average coverage of the unbinned contigs (get coords from missing.coords)
# 3a. do 1000 reps of pull random seqs of same lengths
# 3b. calculate distance and p-value

# Output is a tab-delimited text file:
# 1. average coverage of the reference genome
# 2. mean coverage for the binned contigs
# 3. std dev of coverage for the binned contigs
# 4. distance (difference) b/t bin and reference genome average coverages
# 5. bootstrap p-value for mean coverage
# 6. mean distance b/t individ CR coverage and genome coverage
# 7. std dev of above
# 8. bootstrap p-value for above
# 9-15 same, but for for the unbinned regions

# set variables
my $n = 1000; # number of bootstraps datasets to generate
my $length_cutoff = 1000; # minimum length of region (CR or NR) to evaluate

# handle command line options
my $self_cov = shift @ARGV;
my $genome_stats = shift @ARGV;
my $CR_coords_file = shift @ARGV;
my $NR_coords_file = shift @ARGV;

# read the self_cov file
open(my $in, $self_cov) || die "Can't parse self_cov file '$self_cov': $!\n";
our %COV;
while (my $line = <$in>) {
    next if ($line =~ /^#/);
    chomp $line;
    my ($id, $pos, $cov) = split/\s+/, $line;
    if (! $cov) { $cov = 0 }
    $COV{$id}[$pos] = $cov;
}

# initialize the DATA structure length to the length of the molecule.
# (parse molecule length(s) from the genome_stats file)
my %SEQ;
open (my $stats, $genome_stats) or die "Can't open stats file '$genome_stats': $!\n";
while (my $line=<$stats>) {
    chomp $line;
    my ($id, $lo, $hi, @therest);
    if ($line =~ />(\S+)\s.*length\=(\d+)/) {
	($id, $lo, $hi) = ($1, 1, $2);
    }
    $SEQ{$id}->{'length'} = $hi;
}
close $stats;

# get the coords for the CRs
# we do this first to save time in case there's a data problem
my $CR = &parse_coords($CR_coords_file);

# get the coords for the NRs
my $NR = &parse_coords($NR_coords_file);

# header for the output
print "#accession\tG %rpt\tG cov"
    . "\tCR %rpt\tCR %rpt dist\tp-value"
    . "\tCR cov\tCR cov dist\tp-value"
    . "\tNR %rpt\tNR %rpt dist\tp-value"
    . "\tNR cov\tNR cov dist\tp-value"
    . "\n";

# Now let's get computational
our $genome_rpt = 0;
foreach my $id (keys %COV) {
    print STDERR "Doing $id\n";
    # first find the average coverage for the whole genome
    my($genome_avg, $genome_sd) = &genome_avg_coverage([1, $SEQ{$id}->{'length'}], $COV{$id});
    my $cutoff = $genome_avg + (2*$genome_sd);
    
    # now determine how much of the genome has higher than average coverage (avg+sd)
    my($genome_rpt_a, $genome_rpt_dist_a) = &calc_rpt_content([[1, $SEQ{$id}->{'length'}]], $COV{$id}, $SEQ{$id}->{'length'}, $genome_rpt, $cutoff);
    $genome_rpt = $genome_rpt_a->[0];
    
    print STDERR "\tCR...\n";
    # now the CRs
    my ($CR_rpt, $CR_rpt_dist,
	$CR_cov, $CR_cov_dist) = &calc_rpt_content($CR->{$id},
						   $COV{$id},
						   $SEQ{$id}->{'length'},
						   $genome_rpt,
						   $cutoff,
						   $genome_avg);
    # calculate the average/stdev
    my $CR_rpt_m = mean(@$CR_rpt);
    my $CR_rpt_sd = stddev(@$CR_rpt);
    my $CR_rpt_dist_m = mean(@$CR_rpt_dist);
    my $CR_rpt_dist_sd = stddev(@$CR_rpt_dist);
    my $CR_cov_m = mean(@$CR_cov);
    my $CR_cov_sd = stddev(@$CR_cov);
    my $CR_cov_dist_m = mean(@$CR_cov_dist);
    my $CR_cov_dist_sd = stddev(@$CR_cov_dist);
#    print STDERR "Sneak peek:\n\tgenome: $genome_rpt\n\tCR: $CR_m +/- $CR_sd\t$CR_dist_m +/- $CR_dist_sd\n";
    # to determine distance significance, perform a bootstrap
    print STDERR "\tCR bootstrap...\n";
    my $CR_boot_rpt_count = 0;
    my $CR_boot_cov_count = 0;
    for (my $b=0; $b<$n; $b++) {
	# an equivalent set of regions is defined based on a random starting coordinate
	my $boot_coords = &random_coords($CR->{$id}, $SEQ{$id}->{'length'});
	# cov is calculated for the regions
	my ($boot_rpt, $boot_rpt_dist,
	    $boot_cov, $boot_cov_dist) = &calc_rpt_content($boot_coords,
						       $COV{$id},
						       $SEQ{$id}->{'length'},
						       $genome_rpt,
						       $cutoff,
						       $genome_avg);
	# the mean of the distances is calculated
	my $boot_rpt_dist_m = mean(@$boot_rpt_dist);
	my $boot_cov_dist_m = mean(@$boot_cov_dist);
	# score if the mean is bigger than or equal to that of the CR
	if ($boot_rpt_dist_m >= $CR_rpt_dist_m) { $CR_boot_rpt_count++ }
	if ($boot_cov_dist_m >= $CR_cov_dist_m) { $CR_boot_cov_count++ }
    }
    # calculate the estimated p-value
    my $CR_rpt_dist_p = $CR_boot_rpt_count/$n;
    my $CR_cov_dist_p = $CR_boot_cov_count/$n;

    print STDERR "\tNR...\n";
    # now the NRs
    my ($NR_rpt, $NR_rpt_dist,
	$NR_cov, $NR_cov_dist) = &calc_rpt_content($NR->{$id},
						   $COV{$id},
						   $SEQ{$id}->{'length'},
						   $genome_rpt,
						   $cutoff,
						   $genome_avg);
    # calculate the average/stdev
    my $NR_rpt_m = mean(@$NR_rpt);
    my $NR_rpt_sd = stddev(@$NR_rpt);
    my $NR_rpt_dist_m = mean(@$NR_rpt_dist);
    my $NR_rpt_dist_sd = stddev(@$NR_rpt_dist);
    my $NR_cov_m = mean(@$NR_cov);
    my $NR_cov_sd = stddev(@$NR_cov);
    my $NR_cov_dist_m = mean(@$NR_cov_dist);
    my $NR_cov_dist_sd = stddev(@$NR_cov_dist);
#    print STDERR "\tNR: $NR_m +/- $NR_sd\t$NR_dist_m +/- $NR_dist_sd\n";
    
    # to determine distance significance, perform a bootstrap
    print STDERR "\tNR bootstrap...\n";
    my $NR_boot_rpt_count = 0;
    my $NR_boot_cov_count = 0;
    for (my $b=0; $b<$n; $b++) {
	# an equivalent set of regions is defined based on a random starting coordinate
	my $boot_coords = &random_coords($NR->{$id}, $SEQ{$id}->{'length'});
	# cov is calculated for the regions
	my ($boot_rpt, $boot_rpt_dist,
	    $boot_cov, $boot_cov_dist) = &calc_rpt_content($boot_coords,
						       $COV{$id},
						       $SEQ{$id}->{'length'},
						       $genome_rpt,
						       $cutoff,
						       $genome_avg);
	# the mean of the distances is calculated
	my $boot_rpt_dist_m = mean(@$boot_rpt_dist);
	my $boot_cov_dist_m = mean(@$boot_cov_dist);
	# score if the mean is bigger than or equal to that of the NR
	if ($boot_rpt_dist_m >= $NR_rpt_dist_m) { $NR_boot_rpt_count++ }
	if ($boot_cov_dist_m >= $NR_cov_dist_m) { $NR_boot_cov_count++ }
    }
    # calculate the estimated p-value
    my $NR_rpt_dist_p = $NR_boot_rpt_count/$n;
    my $NR_cov_dist_p = $NR_boot_cov_count/$n;

    # print the output
    printf "$id\t%.4f\t%.2f\t" . 
    "%.4f +/- %.4f\t%.4f +/- %.4f\t%.3f\t" . 
    "%.2f +/- %.2f\t%.2f +/- %.2f\t%.3f\t" . 
    "%.4f +/- %.4f\t%.4f +/- %.4f\t%.3f\t" . 
    "%.2f +/- %.2f\t%.2f +/- %.2f\t%.3f\n" ,
    ($genome_rpt, $genome_avg,
     $CR_rpt_m, $CR_rpt_sd, $CR_rpt_dist_m, $CR_rpt_dist_sd, $CR_rpt_dist_p,
     $CR_cov_m, $CR_cov_sd, $CR_cov_dist_m, $CR_cov_dist_sd, $CR_cov_dist_p,
     $NR_rpt_m, $NR_rpt_sd, $NR_rpt_dist_m, $NR_rpt_dist_sd, $NR_rpt_dist_p,
     $NR_cov_m, $NR_cov_sd, $NR_cov_dist_m, $NR_cov_dist_sd, $NR_cov_dist_p);
}

exit();

sub genome_avg_coverage {
    my $coordset = shift;
    my $covref = shift;
    my @val;
    my ($lo, $hi) = @$coordset;
    my $avg = mean(@{$covref}[$lo..$hi]);
    my $sd = stddev(@{$covref}[$lo..$hi]);
    return ($avg, $sd);
}

sub calc_rpt_content {
    my $coord_ref = shift;
    my $cov_ref = shift;
    my $seqlen = shift;
    my $genome_rpt = shift;
    my $cutoff = shift;
    my $genome_avg = shift;

    my @REG_RPT;
    my @REG_RPT_DIST;
    my @REG_COV;
    my @REG_COV_DIST;
    foreach my $coordset(@$coord_ref) {
	my ($lo, $hi) = @$coordset;
	my $lgth = $hi - ($lo - 1);
	my $rpt = 0; # stores the summed length of repeat regions
	my $avg_cov;
	# if hi is larger than seqlen, assume a circular molecule
	if ($hi > $seqlen) {
	    for(my $i=$lo; $i<=$seqlen; $i++) {
		if ($cov_ref->[$i] > $cutoff) {
		    $rpt++;
		}
	    }
	    for(my $i=1; $i<=($hi - $seqlen); $i++) {
		if ($cov_ref->[$i] > $cutoff) {
		    $rpt++;
		}
	    }
	    $avg_cov = mean(@{$cov_ref}[$lo..$seqlen],@{$cov_ref}[1..($hi-$seqlen)]);
	} elsif ($hi < $lo ) {
	    for(my $i=$lo; $i<=$seqlen; $i++) {
		if ($cov_ref->[$i] > $cutoff) {
		    $rpt++;
		}
	    }
	    for(my $i=1; $i<=$hi; $i++) {
		if ($cov_ref->[$i] > $cutoff) {
		    $rpt++;
		}
	    }
	    $avg_cov = mean(@{$cov_ref}[$lo..$seqlen],@{$cov_ref}[1..$hi]);	
	} else { # otherwise, just examine the region
	    for(my $i=$lo; $i<=$hi; $i++) {
		if ($cov_ref->[$i] > $cutoff) {
		    $rpt++;
		}
	    }
	    $avg_cov = mean(@{$cov_ref}[$lo..$hi]);
	}
	# calculate the percent repeat for the region
	my $per_rpt = $rpt/$lgth;
	push @REG_RPT, $per_rpt;
	push @REG_RPT_DIST, abs($genome_rpt - $per_rpt);

	push @REG_COV, $avg_cov;
	push @REG_COV_DIST, abs($genome_avg - $avg_cov);
    }
    return(\@REG_RPT, \@REG_RPT_DIST, \@REG_COV, \@REG_COV_DIST);
}

sub parse_coords {
    my $coords_file = shift;
    my %B;
    open(my $c, $coords_file) || die "Can't open $coords_file: $!\n";
    while (my $line = <$c>) {
	chomp $line;
	my ($seq_id, $low, $high, $length, $contig) = split(/\t/, $line);
	
	# make sure the coordinates are numeric
	if ($low !~ /^\d+$/ || $high !~ /^\d+$/) { die "Bad format for file $coords_file? Low coord: $low, High coord: $high\n"; }
	
	# make sure the seq_id is present in the genome fasta
	if (! defined $COV{$seq_id}) { die "Ids don't match between file ($coords_file) and genome coverage file ($self_cov): Can't find '$seq_id' in file\n"; }

	# warn about regions spanning the origin
	if ($high < $low) { warn "!!!FOR $coords_file:: $seq_id, low ($low) is greater than high ($high)\n"; }

	# region length screen
	if ($high - $low + 1 < $length_cutoff) { warn "!!!FOR $coords_file:: $seq_id, $low/$high is < $length_cutoff nt\n"; next; }

	push @{$B{$seq_id}}, [$low, $high];
    }
    return \%B;
}

sub random_coords {
    # take a list of coordinate sets
    # generate a random start site
    # modify all coord sets to be relative to that start site
    # this guarantees non-overlapping regions of the same length and number
    my $in_coords = shift;
    my $seq_length = shift;

    # generate the random starting point
    my $base = int rand ($seq_length);

    my @rand_coords;

    # modify each coordinate set to be relative to the random starting point
    for my $coordset(@$in_coords) {
	my $length = abs($coordset->[1] - $coordset->[0]) + 1;

	# remember to screen length
	if ($length < 1000) { die "Bad length? $length\n"; }

	my $low = $coordset->[0] + $base;
	my $high = $low - 1 + $length;

	# handle alterations that span the origin
	if ($low > $seq_length) {
	    $low = $low - $seq_length;
	}
	if ($high > $seq_length) {
	    $high = $high - $seq_length;
	}

	push @rand_coords, [$low, $high];
    }
    return(\@rand_coords);
}
