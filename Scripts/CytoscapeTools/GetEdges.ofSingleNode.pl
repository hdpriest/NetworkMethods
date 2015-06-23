#!/usr/bin/perl
use warnings;
use strict;

use lib '/home/hpriest/Scripts/Library';
use hdpTools;

my $usage="perl $0 <csv edge file> <node to find edges of> <threshold for value>\n";
die $usage unless $#ARGV==2;

my $node = $ARGV[1];

open(NET,"<",$ARGV[0]) || die "Cannot open $ARGV[0]!\n$!\nexiting...\n";
until(eof(NET)){
	my $line=<NET>;
	chomp $line;
	my @line=split(/\s/,$line);
	if((defined $line[0] eq $node) || ($line[1] eq $node) || ($line[4] eq $node) || ($line[5] eq $node)){
		print $line."\n" if abs($line[2]) >= $ARGV[2];
	}else{
	}
}
close NET;
