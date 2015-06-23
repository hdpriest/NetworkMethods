#!/usr/bin/perl
use warnings;
use strict;

my $usage = "perl $0 <cyto csv file> <cutoff>\n\n";

die $usage unless $#ARGV==1;

my $file=$ARGV[0];
my $cutoff=$ARGV[1];

open(FILE,"<",$file) || die "cannot open $file!\n$!\nexiting...\n";
until(eof(FILE)){
	my $line=<FILE>;
	next if $line=~m/weight/i;
	chomp $line;
	my @line=split(/\s/,$line);
	next if $line[0] eq $line[1];
	print $line."\n" if abs($line[2])>=$cutoff;
}
close FILE;
