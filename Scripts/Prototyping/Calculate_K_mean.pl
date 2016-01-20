#!/usr/bin/perl
use warnings;
use strict;

my $frame = $ARGV[0];

die "need a file of adjacency values\n" unless -e $frame;

my $sum=0;
my $N=0;
open(KD,"<",$frame) || die "Cannot open $frame!\n$!\nexiting...\n";
until(eof(KD)){
	my $line = <KD>;
	chomp $line;
	my @line = split(/\t/,$line);
	$sum += $line[2];
	$N++;
}
my $mean = $sum/$N;
warn "Total = $sum\nN = $N\nMean = $mean\n";
