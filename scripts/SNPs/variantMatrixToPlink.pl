#!/usr/bin/perl
#

use strict;
use warnings;
use Data::Dumper;
my $usage = "usage: $0 matrix_file\n";

my $matfile = $ARGV[0] or die $usage;

### Read in matrix ####
 
unless (open (MATFILE, $matfile)) {
   die "File $matfile not found.\n";
}
 
my %cols_to_features;
my @ordered_features;
my %profiles;
my $line = 1;
my @genome_list;
 
foreach my $in (<MATFILE>) {
   chomp $in;
   my @x = split(/\t/, $in);
   if ($line == 1) {
       shift @x;
       my $i = 1;
       foreach my $x (@x) {
           $cols_to_features{$i} = $x;
           push(@ordered_features, $x);
           ++$i;
       }
       ++$line;
   } elsif ($in =~ /^#annotation/) {
   } else {
       my $genome = shift @x;
       push (@genome_list, $genome);
       my $i = 1;
       foreach my $x (@x) {
           $profiles{$genome}{$cols_to_features{$i}} = $x;
           ++$i;
       }        
   }
}
close MATFILE;

my $pedfilename = "$matfile.ped";
open (PEDFILE, ">$pedfilename");

foreach my $genome (@genome_list) {
	print PEDFILE "$genome 1";
	foreach my $ordered_feature (@ordered_features) {
		if ($profiles{$genome}{$ordered_feature} eq '0') {	
			print PEDFILE "  1 1";
		} elsif ($profiles{$genome}{$ordered_feature} eq '1') {
			print PEDFILE "  2 2";	
		} elsif ($profiles{$genome}{$ordered_feature} eq '-') {
			print PEDFILE "  0 0";
		} else {
			print PEDFILE "  0 0";
			print "Invalid coding: $profiles{$genome}{$ordered_feature}, logged as missing data.\n";
		}
	}
	print PEDFILE "\n";
}

close PEDFILE;

my $mapfilename = "$matfile.map";
open (MAPFILE, ">$mapfilename");
foreach my $ordered_feature (@ordered_features) {
	print MAPFILE "0 $ordered_feature 0 1\n";
}
close MAPFILE;