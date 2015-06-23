#!/usr/bin/perl
use warnings;
use strict;

use lib '/home/hpriest/Scripts/Library';
use hdpTools;

my $usage="perl $0 <csv edge file> <set of nodes to extract edges between> <threshold for value>\n";
die $usage unless $#ARGV==2;

my @nodes=@{hdpTools->LoadFile($ARGV[1])};
my %n;
map {$n{$_}=1} @nodes;

open(NET,"<",$ARGV[0]) || die "Cannot open $ARGV[0]!\n$!\nexiting...\n";
until(eof(NET)){
	my $line=<NET>;
	chomp $line;
	my @line=split(/\s/,$line);
	if((defined $n{$line[0]})&&(defined $n{$line[1]})){
		print $line."\n" if abs($line[2]) >= $ARGV[2];
	}else{
	}
}
close NET;
