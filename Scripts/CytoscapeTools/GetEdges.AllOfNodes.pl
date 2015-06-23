#!/usr/bin/perl
use warnings;
use strict;

use lib '/home/hpriest/Scripts/Library';
use hdpTools;

my $usage="perl $0 <csv edge file> <set of nodes> <threshold for value>\n\nGet Edges involving nodes of list\n\n";
die $usage unless $#ARGV==2;

my @nodes=@{hdpTools->LoadFile($ARGV[1])};
my %n1;
map {$n1{$_}=1} @nodes;

open(NET,"<",$ARGV[0]) || die "Cannot open $ARGV[0]!\n$!\nexiting...\n";
until(eof(NET)){
	my $line=<NET>;
	chomp $line;
	my @line=split(/\s/,$line);
	if((defined $n1{$line[0]})||(defined $n1{$line[1]})||(defined($n1{$line[4]}))||(defined($n1{$line[5]}))){
		print $line."\n" if abs($line[2])>= $ARGV[2];
	}
}
close NET;
