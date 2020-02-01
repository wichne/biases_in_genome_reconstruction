#!/usr/bin/perl
use Bio::SeqIO;
use Statistics::Basic qw(:all);
use strict;

if (! @ARGV ||
    $ARGV[0] eq "-h") {
    print "GC_ANALYSIS

USAGE: $0 genome.fa CR.coords NR.coords

genome.fa - fasta file containing the genome sequence

CR.coords - a tab-delimited file describing the correctly binned regions of the
            genome. It must contain of the following fields:
            0 sequence accession (must match an accession in the fasta file)
            1 low coord
            2 high coord

NR.coords - a file with the same format as above, but describing the not-binned
            regions of the genome.

";
    exit();
}

my $n = 1000; # number of bootstrap datasets to generate
my $length_cutoff = 1000; # minimum length of region (CR or NR) to evaluate


my $genome_fasta_file = shift @ARGV;
my $binned_coords_file = shift @ARGV;
my $unbinned_coords_file = shift @ARGV;

my $gfo = Bio::SeqIO->new(-file=>$genome_fasta_file);

# get the whole genome
my %SEQ;
while (my $seq_obj = $gfo->next_seq) {
    my $id = $seq_obj->display_id;
    $SEQ{$id} = $seq_obj;
}

# get the coords for the binned
# we do this first to save time in case there's a data problem
open(my $bc, $binned_coords_file) || die "Can't open $binned_coords_file: $!\n";
my %BIN;
while (my $line = <$bc>) {
    chomp $line;
    my ($seq_id, $low, $high, $length, $contig) = split(/\t/, $line);
    if ($high < $low) { warn "!!!FOR BIN $seq_id, low ($low) is greater than high ($high)\n"; next; }
    if ($high - $low + 1 < $length_cutoff) { warn "!!!FOR BIN $seq_id, $low/$high is < $length_cutoff nt\n"; next; }
    if (! defined $SEQ{$seq_id}) { die "Ids don't match between binned.coords file ($binned_coords_file) and genome fasta file ($genome_fasta_file): Can't find '$seq_id' in fasta file\n"; }
    if ($low !~ /^\d+$/ || $high !~ /^\d+$/) { die "Bad format for binned.coords file ($binned_coords_file)? Low coord: $low, High coord: $high\n"; }
    push @{$BIN{$seq_id}}, [$low, $high];
}

# get the coords for the unbinned
open(my $bc, $unbinned_coords_file) || die "Can't open $unbinned_coords_file: $!\n";
my %UNBIN;
while (my $line = <$bc>) {
    chomp $line;
    my ($seq_id, $low, $high, $length, $contig) = split(/\t/, $line);
    if ($high < $low) { warn "!!!FOR UNBIN $seq_id, low ($low) is greater than high ($high)\n"; next; }
    if ($high - $low + 1 < $length_cutoff) { warn "!!!FOR UNBIN $seq_id, $low/$high is < $length_cutoff nt\n"; next; }
    if (! defined $SEQ{$seq_id}) { die "Ids don't match between unbinned.coords file ($unbinned_coords_file) and genome fasta file ($genome_fasta_file): Can't find '$seq_id' in fasta file\n"; }
    if ($low !~ /^\d+$/ || $high !~ /^\d+$/) { die "Bad format for unbinned.coords file ($unbinned_coords_file)? Low coord: $low, High coord: $high\n"; }
    push @{$UNBIN{$seq_id}}, [$low, $high];
}

# output header line
print "Molecule id\t%G+C\tCR %G+C mean\tCR %G+C stdev\tCR dist mean\tCR dist stdev\tCR dist p\tNR %G+C mean\tNR %G+C stdev\tNR dist mean\tNR dist stdev\tNR dist p\n";
# Now let's get computational
our $GC_genome;
foreach my $id (keys %SEQ) {
    print STDERR "Doing $id\n";
    # first the whole genome
    $GC_genome = &calc_gc($SEQ{$id});

    print STDERR "\tCR...\n";
    # calculations for the correctly binned
    my ($CR_gc, $CR_gc_dist) = &get_region_gc($BIN{$id}, $SEQ{$id});

    # calculate the average/stdev %G+C and the avg/stdev distance
    my $CR_gc_m = mean(@$CR_gc);
    my $CR_gc_sd = stddev(@$CR_gc);
    my $CR_dist_m = mean(@$CR_gc_dist);
    my $CR_dist_sd = stddev(@$CR_gc_dist);
    
    # to determine distance significance, perform a bootstrap

    print STDERR "\tCR bootstrap...\n";
    my $CR_boot_count = 0;
    for (my $b=0; $b<$n; $b++) {
	my $boot_coords = &random_coords($BIN{$id}, $SEQ{$id}->length);
	my($boot_gc, $boot_gc_dist) = &get_region_gc($boot_coords, $SEQ{$id});
	my $boot_dist_m = mean(@$boot_gc_dist);
	if ($boot_dist_m >= $CR_dist_m) { $CR_boot_count++ }
    }
    my $CR_dist_p = $CR_boot_count/$n;
    
    print STDERR "\tNR...\n";
    # now the not-binned
    my ($NR_gc, $NR_gc_dist) = &get_region_gc($UNBIN{$id}, $SEQ{$id});

    # calculate the average/stdev %G+C and the avg/stdev distance
    my $NR_gc_m = mean(@$NR_gc);
    my $NR_gc_sd = stddev(@$NR_gc);
    my $NR_dist_m = mean(@$NR_gc_dist);
    my $NR_dist_sd = stddev(@$NR_gc_dist);
    
    # to determine distance significance, perform a bootstrap
    print STDERR "\tNR bootstrap...\n";
    my $NR_boot_count = 0;
    for (my $b=0; $b<$n; $b++) {
	my $boot_coords = &random_coords($UNBIN{$id}, $SEQ{$id}->length);
	my($boot_gc, $boot_gc_dist) = &get_region_gc($boot_coords, $SEQ{$id});
	my $boot_dist_m = mean(@$boot_gc_dist);
	if ($boot_dist_m >= $NR_dist_m) { $NR_boot_count++ }
    }
    my $NR_dist_p = $NR_boot_count/$n;
    
    # Print the results
    printf "$id\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.3f\t%.2f\t%.2f\t%.2f\t%.2f\t%.3f\n",
    ($GC_genome,
     $CR_gc_m, $CR_gc_sd, $CR_dist_m, $CR_dist_sd, $CR_dist_p,
     $NR_gc_m, $NR_gc_sd, $NR_dist_m, $NR_dist_sd, $NR_dist_p);
}

exit();

sub calc_gc {
    my $seq_obj = shift;
    my $seq = $seq_obj->seq;
    $seq =~ s/[^GATCgatc]+//g;
    my $len = length($seq);
    if ($len == 0) { warn "seq had no length after screening for N's\n"; return -1; }
    my %CHAR;
#    while ($seq =~ /([GCgc])/g) { $CHAR{$1}++ }
    $seq =~ s/[^GCgc]//g;
#    return sprintf "%.2f",(($CHAR{G} + $CHAR{C} + $CHAR{g} + $CHAR{c})/length($seq) * 100);
    return sprintf "%.2f", (length($seq)/$len * 100);
}

sub get_region_gc {
    my $coord_ref = shift;
    my $seq_ref = shift;
    my @GC;
    my @GCdist;
    my @RAND;

    # for each coordinate set describing a binned or not-binned region,
    # calculate the %G+C
    # calculate the distance (difference) between the region %G+C and the complete genome value
    # also make random coord sets of the same length for bootstrapping
    foreach my $coordset (@$coord_ref) {
	my ($low, $high) = @$coordset;
	my $gc;
	if ($high > $seq_ref->length) {
	    $high -= $seq_ref->length;
	    my $seq = $seq_ref->subseq($low, $seq_ref->length);
	    $seq .= $seq_ref->subseq(1, $high);
	    my $lref = Bio::Seq->new(-seq => $seq, -display_id => join("_", ($seq_ref->display_id, $low, $high)));
	    $gc = &calc_gc($lref);
	} elsif ($high < $low ) {
	    #warn("High $high less than low $low");
	    my $seq = $seq_ref->subseq($low, $seq_ref->length);
	    $seq .= $seq_ref->subseq(1, $high);
	    my $lref = Bio::Seq->new(-seq => $seq, -display_id => join("_", ($seq_ref->display_id, $low, $high)));
	    $gc = &calc_gc($lref);
	    
	} else {
	    $gc = &calc_gc($seq_ref->trunc($low, $high));
	}
	push @GC, $gc;
	push @GCdist, abs($GC_genome - $gc);
    }

    return (\@GC, \@GCdist);
}

sub random_coords {
    # take coords set
    # generate a random start site
    # modify all coords to be relative to that start site
    # this guarantees non-overlapping regions of the same length and number
    my $in_coords = shift;
    my $seq_length = shift;

    my $base = int rand ($seq_length);
    my @rand_coords;
    for my $coordset(@$in_coords) {
	my $length = abs($coordset->[1] - $coordset->[0]) + 1;
	if ($length < 1000) { die "Bad length? $length\n"; }
	my $low = $coordset->[0] + $base;
	if ($low > $seq_length) {
	    $low = $low - $seq_length;
	}
	my $high = $low - 1 + $length;
	if ($high > $seq_length) {
	    $high = $high - $seq_length;
	}
	push @rand_coords, [$low, $high];
    }

    return(\@rand_coords);
}
