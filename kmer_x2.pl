#!/usr/bin/perl
use Bio::SeqIO;
use Statistics::Basic qw(:all);
use Number::Format qw(:subs);
use strict;

$| = 1;

if (! @ARGV ||
    $ARGV[0] eq "-h") {
    print "TETRANUCLEOTIDE CHI-SQUARED ANALYSIS

Input is the genome sequence in fasta format, and two files listing the genome coordinates of the binned
and not-binned regions (described below).
The program uses the genome tetranucleotide frequencies as the expected to calculate
 1) the mean x2 value for the binned regions (and standard deviation)
 2) a boot-strapped p-value (n=1000) for the binned x2
 3) the mean x2 value for the not-binned regions (and standard deviation)
 4) a boot-strapped p-value (n=1000) for the not-binned x2
The program writes to STDOUT.
The output is tab-delimited and the format is defined by the first line of the output.

USAGE: $0 genome.fa CR.coords NR.coords > output.tsv

genome.fa - fasta file containing the genome sequence

CR.coords - a tab-delimited file describing the correctly binned regions of the
            genome. It must contain of the following fields:
            0 sequence accession (must match an accession in the fasta file)
            1 low coord
            2 high coord
            !!! If the high coord is lower than the low coord, the program assumes
            !!! the chromosome is circular and the region spans the origin.

NR.coords - a file with the same format as above, but describing the not-binned
            regions of the genome.

";
    exit();
}

my $n = 1000; # number of bootstrap datasets to generate
my $length_cutoff = 1000; # minimum length of region (CR or NR) to evaluate
our $kmer = 4; # length of kmer to evaluate

my $genome_fasta_file = shift @ARGV;
my $CR_coords_file = shift @ARGV;
my $NR_coords_file = shift @ARGV;

my $gfo = Bio::SeqIO->new(-file=>$genome_fasta_file);

# get the whole genome
my %SEQ;
while (my $seq_obj = $gfo->next_seq) {
    my $id = $seq_obj->display_id;
    $SEQ{$id} = $seq_obj;
}

# get the coords for the correctly binned
# we do this first to save time in case there's a data problem
my $CR = &parse_coords($CR_coords_file);

# get the coords for the unbinned
my $NR = &parse_coords($NR_coords_file);

# print output header line
print "Molecule id\tCR x2 mean\tCR x2 stdev\tCR x2 p\tNR x2 mean\tNR x2 stdev\tNR x2 p\n";

# Now let's get computational
foreach my $id (keys %SEQ) {
    print STDERR "Doing $id\n";
    # count kmers for the genome
    my $gen = &kmer_freq($SEQ{$id});

    print STDERR "\tCR...\n";
    # calculate the chi-squared for each region
    my $CR_x2s = &get_region_x2($CR->{$id}, $SEQ{$id}, $gen);

    # calculate the average/stdev x2
    my $CR_m = mean(@$CR_x2s);
    my $CR_sd = stddev(@$CR_x2s);
    
    # to determine distance significance, perform a bootstrap
    print STDERR "\tCR bootstrap...\n";
    my $CR_boot_count = 0;
    for (my $b=0; $b<$n; $b++) {
	# an equivalent set of regions is defined based on a random starting coordinate
	my $boot_coords = &random_coords($CR->{$id}, $SEQ{$id}->length);
	# x2 is calculated for the regions
	my $boot_x2s = &get_region_x2($boot_coords, $SEQ{$id}, $gen);
	# the mean of the distances is calculated
	my $boot_m = mean(@$boot_x2s);
	# score if the mean is less than or equal to that of the CR
	if ($boot_m >= $CR_m) { $CR_boot_count++ }
    }
    # calculate the estimated p-value
    my $CR_p = $CR_boot_count/$n;
    
    print STDERR "\tNR...\n";
    # now the not-binned
    my $NR_x2s = &get_region_x2($NR->{$id}, $SEQ{$id}, $gen);

    # calculate the average/stdev %G+C and the avg/stdev distance
    my $NR_m = mean(@$NR_x2s);
    my $NR_sd = stddev(@$NR_x2s);
    
    # to determine distance significance, perform a bootstrap
    print STDERR "\tNR bootstrap...\n";
    my $NR_boot_count = 0;
    for (my $b=0; $b<$n; $b++) {
	my $boot_coords = &random_coords($NR->{$id}, $SEQ{$id}->length);
	my $boot_x2s = &get_region_x2($boot_coords, $SEQ{$id}, $gen);
	my $boot_m = mean(@$boot_x2s);
	if ($boot_m >= $NR_m) { $NR_boot_count++ }
    }
    my $NR_p = $NR_boot_count/$n;
    
    # Print the results
    printf "$id\t" .
    "%.4f\t%.4f\t%.4f\t" . 
    "%.4f\t%.4f\t%.4f\n",
    ($CR_m, $CR_sd, $CR_p,
     $NR_m, $NR_sd, $NR_p);
}

exit();

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
	if (! defined $SEQ{$seq_id}) { die "Ids don't match between file ($coords_file) and genome fasta file ($genome_fasta_file): Can't find '$seq_id' in fasta file\n"; }

	# warn about regions spanning the origin
	if ($high < $low) { warn "!!!FOR $coords_file:: $seq_id, low ($low) is greater than high ($high)\n"; }

	# region length screen
	if ($high - $low + 1 < $length_cutoff) { warn "!!!FOR $coords_file:: $seq_id, $low/$high is < $length_cutoff nt\n"; next; }

	push @{$B{$seq_id}}, [$low, $high];
    }
    return \%B;
}

sub count_kmer {
    my $seq = shift;
    my $C = shift;
    while ($seq =~ /(.{$kmer})/g) {
	my $k = $1;
	$C->{$k}++ if ($k !~ /[^GATCgatc]/); # ignore k-mers with ambiguous symbols.
    }
}

sub kmer_freq {
    my $seqo = shift;
    # generate the reverse strand
    my $revo = $seqo->revcom();

    my $seq = $seqo->seq;
    my $revseq = $revo->seq;
    
    # save this for the frequency calculation later
    my $seqlen = length($seq);

    my %KMER;
    # forward
    count_kmer($seq, \%KMER);
 
    # reverse
    count_kmer($revseq, \%KMER);

    my %FREQ;
    foreach my $mer (keys %KMER) {
	# since we counted forward and reverse, each kmer count is divided by
	# 2 times the number of kmers in the genome (ie, genome/kmer) 
	$FREQ{$mer} = $KMER{$mer}/(2*$seqlen/$kmer);
    }
    return \%FREQ;
}

sub chi_squared {
    my $exp = shift;
    my $obs = shift;

    my $x2;
    foreach my $t (keys %$exp) {
	$x2 += (($obs->{$t} - $exp->{$t})**2)/$exp->{$t};
    }
    return $x2;
}

sub get_region_x2 {
    my $coord_ref = shift;
    my $seq_ref = shift;
    my $EXP = shift;
    
    # for each coordinate set describing a binned or not-binned region,
    # calculate the trinucleotide frequency
    # calculate the chi-squared statistic using the genome value as the expected frequency
    my @chi_data;
    
    foreach my $coordset (@$coord_ref) {
	my ($low, $high) = @$coordset;
	my $gc;

	# handle origin-spanning regions
	# assume a high coordinate higher than the molecule length indicates a spanning region
	my $lref;
	if ($high > $seq_ref->length) {
	    $high -= $seq_ref->length;
	    my $seq = $seq_ref->subseq($low, $seq_ref->length);
	    $seq .= $seq_ref->subseq(1, $high);
	    $lref = Bio::Seq->new(-seq => $seq, -display_id => join("_", ($seq_ref->display_id, $low, $high)));

	# assume a high coordinate lower than a low coordinate indicates a spanning region
	} elsif ($high < $low ) {
	    #warn("High $high less than low $low");
	    my $seq = $seq_ref->subseq($low, $seq_ref->length);
	    $seq .= $seq_ref->subseq(1, $high);
	    $lref = Bio::Seq->new(-seq => $seq, -display_id => join("_", ($seq_ref->display_id, $low, $high)));
	} else {
	    $lref = $seq_ref->trunc($low, $high);
	}

	# calculate the frequencies and calculate the chi-squared
	my $kmerfreq = &kmer_freq($lref);
	my $x2 = &chi_squared($EXP, $kmerfreq);

	# collect the data
	push @chi_data, $x2;
    }
    return (\@chi_data);
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
