#!/usr/bin/env perl
#
# modified from fasta2pwise.pl
#
# calculates pairwise snp counts for an aligned multifasta

use strict;
use warnings;
use Data::Dumper;

my $usage = "usage: $0 input_file\n";

my $in_name = $ARGV[0] or die $usage;
chomp $in_name;

unless (open (INFILE, $in_name)) {
    die "File $in_name not found.\n";
}

my %fasta_db;
my $name = '';

foreach my $in (<INFILE>) {
    chomp $in;
    if ($in =~ /^\s*$/) {
    } elsif ($in =~ /^>(\S+)/) {
        $name = $1;
        $fasta_db{$name} = '';
    } else {
	    $fasta_db{$name} .= $in;
    }
}
close INFILE;

#print Dumper (\%fasta_db);

my @seq_names = sort (keys %fasta_db);

#my %pair_sim;

print "Taxon";
foreach my $seq_name (@seq_names) {
    print "\t$seq_name";
}
print "\n";

foreach my $seq_name1 (@seq_names) {
    print "$seq_name1";
    foreach my $seq_name2 (@seq_names) {
        my @x = split ('', $fasta_db{$seq_name1});
        my @y = split ('', $fasta_db{$seq_name2});
        #my $length = 0;
        my $snps = 0;
        while (@x) {
            my $xbase = shift @x;
            my $ybase = shift @y;
            if ($xbase eq '-' or $ybase eq '-' or $xbase eq 'N' or $ybase eq 'N') {
                #++$length;
            } elsif ($xbase ne $ybase) {
                    ++$snps;
            }
        }
        #$pair_sim{$seq_name1}{$seq_name2} = $idents/$length;
        #my $sim = sprintf("%.4f", ($idents/$length));
        print "\t$snps";
    }
    print "\n";
}

#print Dumper (\%pair_sim);
            



exit;
