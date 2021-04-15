#!/usr/bin/perl
#
# addColumn.pl
#
# adds tab delimited info to tab delimited file
# 

use strict;
use warnings;
# use Data::Dumper;

my $usage = "usage: $0 base_file tab_file base_label_column tab_label_column tab_info_column > output_file\n";

my $pfilename = $ARGV[0] or die $usage;
my $tfilename = $ARGV[1] or die $usage;
my $bcol = $ARGV[2] or die $usage;
my $lcol = $ARGV[3] or die $usage;
my $icol = $ARGV[4] or die $usage;
chomp $icol;

if ($bcol eq 'zero') {
    $bcol = 0;
}
if ($icol eq 'zero') {
    $icol = 0;
}
if ($lcol eq 'zero') {
    $lcol = 0;
}

unless (open (TFILE, $tfilename)) {
    die "File not found.\n";
}

my %tinfo;

foreach my $t (<TFILE>) {
    chomp $t;
    my @x = split(/\t/, $t);
    $tinfo{$x[$lcol]} = $x[$icol];
}

#print Dumper(\%tinfo);

close TFILE;

unless (open (PFILE, $pfilename)) {
    die "File not found.\n";
}

foreach my $p (<PFILE>) {
    chomp $p;
    my @x = split(/\t/, $p);
    print "$p";
    if (exists $tinfo{$x[$bcol]}) {
        print "\t$tinfo{$x[$bcol]}\n";
    } else {
        print STDERR"no info for $x[$bcol]\n";
        print "\t.\n";
    }
}
 
close PFILE;





exit;
