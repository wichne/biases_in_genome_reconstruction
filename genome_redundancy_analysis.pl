#!/usr/bin/perl
use Bio::SeqIO;
use Statistics::Basic qw(:all);
use strict;

$| = 1;

if (! @ARGV ||
    $ARGV[0] eq "-h") {
    print "GENOME REDUNDANCY ANALYSIS

USAGE: $0 self_cov.tsv genome_stats.txt CR.coords NR.coords

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

";
    exit();
}

# The process:
# 1. calculate the genome redundancy for the genome (parse from fastaStats?)
# 2. calculate the genome redundancy of the binned contigs (get coords from binned.coords)
# 2a. do 1000 reps of pull random seqs of same lengths
# 2b. calculate distance and p-value
# 3. calculate the redundancy of the unbinned contigs (get coords from missing.coords)
# 3a. do 1000 reps of pull random seqs of same lengths
# 3b. calculate distance and p-value

# Output is a tab-delimited text file:
# 1. redundancy of the reference genome
# 2. mean redundancy for the binned sequences
# 3. std dev of redundancy for the binned sequences
# 4. distance (difference) b/t bin and reference genome redundancies
# 5. bootstrap p-value for mean redundancy
# 6. mean distance b/t individ CDR redundancy and genome redundancy
# 7. std dev of above
# 8. bootstrap p-value for above
# 9-15 same, but for for the unbinned regions

my $n = 1000; # number of bootstraps datasets to generate
my $length_cutoff = 1000; # minimum length of region (CR or NR) to evaluate

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

# get the coords for the whole genome
# set the DATA structure length to the length of the molecule.
my %SEQ;
my %REG;
open (my $stats, $genome_stats) or die "Can't open stats file '$genome_stats': $!\n";
while (my $line=<$stats>) {
    chomp $line;
    my ($id, $lo, $hi, @therest);
    if ($line =~ />(\S+)\s.*length\=(\d+)/) {
	($id, $lo, $hi) = ($1, 1, $2);
    }
    $SEQ{$id}->{'length'} = $hi;
    push @{$REG{$id}}, [$lo, $hi];
}
close $stats;

# get the coords for the binned
# we do this first to save time in case there's a data problem
my $CR = &parse_coords($CR_coords_file);

# get the coords for the unbinned
my $NR = &parse_coords($NR_coords_file);

# Now let's get computational
our $redund_genome;
print "#accession\tgenome_mean\tCR_mean\tCR_sd\tCR_dist_m\tCR_dist_sd\tCR_dist_p\tNR_mean\tNR_sd\tNR_dist_m\tNR_dist_sd\tNR_dist_p\n";
foreach my $id (keys %COV) {
    print STDERR "Doing $id\n";
    # first the whole genome
    $redund_genome = &calc_redundancy($REG{$id}, $COV{$id});

    print STDERR "\tCR...\n";
    # now the CR
    my ($CR_redund, $CR_redund_dist) = &calc_reg_redundancy($CR->{$id}, $COV{$id}, $SEQ{$id}->{'length'});
    # calculate the average/stdev
    my $CR_m = mean(@$CR_redund);
    my $CR_sd = stddev(@$CR_redund);
    my $CR_dist_m = mean(@$CR_redund_dist);
    my $CR_dist_sd = stddev(@$CR_redund_dist);
    
    # to determine distance significance, perform a bootstrap
    print STDERR "\tCR bootstrap...\n";
    my $CR_boot_count = 0;
    for (my $b=0; $b<$n; $b++) {
	# an equivalent set of regions is defined based on a random starting coordinate
	my $boot_coords = &random_coords($CR->{$id}, $SEQ{$id}->{'length'});
	# redund is calculated for the regions
	my ($boot_redund, $boot_redund_dist) = &calc_reg_redundancy($boot_coords, $COV{$id}, $SEQ{$id}->{'length'});
	# the mean of the distances is calculated
	my $boot_dist_m = mean(@$boot_redund_dist);
	# score if the mean is less than or equal to that of the CR
	if ($boot_dist_m >= $CR_dist_m) { $CR_boot_count++ }
    }
    # calculate the estimated p-value
    my $CR_dist_p = $CR_boot_count/$n;

    print STDERR "\tNR...\n";
    # now the NR
    my ($NR_redund, $NR_redund_dist) = &calc_reg_redundancy($NR->{$id}, $COV{$id}, $SEQ{$id}->{'length'});
    # calculate the average/stdev
    my $NR_m = mean(@$NR_redund);
    my $NR_sd = stddev(@$NR_redund);
    my $NR_dist_m = mean(@$NR_redund_dist);
    my $NR_dist_sd = stddev(@$NR_redund_dist);

    # to determine distance significance, perform a bootstrap
    print STDERR "\tNR bootstrap...\n";
    my $NR_boot_count = 0;
    for (my $b=0; $b<$n; $b++) {
	# an equivalent set of regions is defined based on a random starting coordinate
	my $boot_coords = &random_coords($NR->{$id}, $SEQ{$id}->{'length'});
	# redund is calculated for the regions
	my ($boot_redund, $boot_redund_dist) = &calc_reg_redundancy($boot_coords, $COV{$id}, $SEQ{$id}->{'length'});
	# the mean of the distances is calculated
	my $boot_dist_m = mean(@$boot_redund_dist);
	# score if the mean is less than or equal to that of the NR
	if ($boot_dist_m >= $NR_dist_m) { $NR_boot_count++ }
    }
    # calculate the estimated p-value
    my $NR_dist_p = $NR_boot_count/$n;

    printf "$id\t%.2f\t" . 
    "%.2f\t%.2f\t%.2f\t%.2f\t%.3f\t" . 
    "%.2f\t%.2f\t%.2f\t%.2f\t%.3f\n" ,
    ($redund_genome,
     $CR_m, $CR_sd, $CR_dist_m, $CR_dist_sd, $CR_dist_p,
     $NR_m, $NR_sd, $NR_dist_m, $NR_dist_sd, $NR_dist_p);
}

exit();

sub calc_redundancy {
    my $coordref = shift;
    my $covref = shift;
    my @val;
    foreach my $coordset (@$coordref) {
	my ($lo, $hi) = @$coordset;
	my $len = $hi - $lo + 1;
	my @reg = @{$covref}[$lo..$hi];
	push @val, @reg;
    }
    return mean(@val);
}

sub calc_reg_redundancy {
    my $coord_ref = shift;
    my $cov_ref = shift;
    my $seqlen = shift;
    my @RDND;
    my @RAND;
    my @DIST;

    # for each coordinate set describing a binned or not-binned region,
    # calculate the redundancy
    # calculate the distance (absolute difference) between the region redundancy and the complete genome value
    foreach my $coordset (@$coord_ref) {
	my($low, $high) = @$coordset;
	my $rdnd;
	
	# handle origin-spanning regions
	# assume a high coordinate higher than the molecule length indicates a spanning region
	if ($high > $seqlen) {
	    $high -= $seqlen;
	    $rdnd = mean(@{$cov_ref}[$low..$seqlen],@{$cov_ref}[1..$high]);
	} elsif ($high < $low ) {
	    $rdnd = mean(@{$cov_ref}[$low..$seqlen],@{$cov_ref}[1..$high]);
	} else {
	    $rdnd = mean(@{$cov_ref}[$low..$high]);
	}
	push @RDND, $rdnd;
	push @DIST, abs($rdnd - $redund_genome);
    }

    return(\@RDND, \@DIST);
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
