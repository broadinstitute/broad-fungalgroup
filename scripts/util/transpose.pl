#!/usr/bin/env perl
#
# transpose.pl
#
# converts a file of columns to rows

use strict;
use warnings;
# use Data::Dumper;

my $usage = "usage: $0 input_file > output_file \n";

my $infilename = $ARGV[0] or die $usage;
chomp $infilename;

unless (open (INFILE, $infilename)) {
    die "File $infilename not found.\n";
}

my $num_cols = 0;
my $num_rows = 0;
my %data;

foreach my $in (<INFILE>) {
        chomp $in;
        my @x = split (/\t/, $in);
        my $cols = scalar(@x);
        if ($cols > $num_cols) {
            $num_cols = $cols;
        }
        ++$num_rows;
        for (my $c = 1; $c < ($num_cols + 1); ++$c) {
            $data{$c}{$num_rows} = $x[($c-1)];
        }
}
close INFILE;

for (my $i = 1; $i < ($num_cols + 1); ++$i) {
    for (my $j = 1; $j < ($num_rows + 1); ++$j) {
        if (defined $data{$i}{$j}) {
            print "$data{$i}{$j}";
        }
        if ($j < $num_rows) {
            print "\t";
        }
    }
    print "\n";
}



exit;

