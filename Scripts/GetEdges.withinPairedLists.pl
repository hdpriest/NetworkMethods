#!/usr/bin/perl
use warnings;
use strict;

use lib '/home/hpriest/Scripts/Library';
use hdpTools;

my $usage="perl $0 <csv edge file> <set of nodes 1> <set of nodes 2>  <threshold for value>\n\nGet Edges between nodes of set A and nodes of set B";
die $usage unless $#ARGV==3;

my @nodes=@{hdpTools->LoadFile($ARGV[1])};
my %n1;
map {$n1{$_}=1} @nodes;
@nodes=@{hdpTools->LoadFile($ARGV[2])};
my %n2;
map {$n2{$_}=1} @nodes;

open(NET,"<",$ARGV[0]) || die "Cannot open $ARGV[0]!\n$!\nexiting...\n";
until(eof(NET)){
	my $line=<NET>;
	chomp $line;
	my @line=split(/\,/,$line);
	if((defined $n1{$line[0]})&&(defined $n2{$line[1]})){
		print $line."\n" if $line[2] >= $ARGV[3];
	}elsif((defined $n2{$line[0]})&&(defined $n1{$line[1]})){
		print $line."\n" if $line[2] >= $ARGV[3];
	}else{
	}
}
close NET;
