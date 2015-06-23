#!/usr/bin/perl
use warnings;
use strict;
use FindBin;
use lib "$FindBin::Bin/../../Lib";
use Tools;
use Configuration;

die "usage: perl $0 <cytoscape file> <plasticity direction> <stanza>\n\n" unless $#ARGV==2;

my $cyto = $ARGV[0];
my $plas = $ARGV[1];
my $stan = $ARGV[2];

my @values;
my $total=0;
open(CYTO,"<",$cyto) || die "cannot open $cyto!\n$!\nexiting...\n";
until(eof(CYTO)){
	my $line=<CYTO>;
	chomp $line;
	my ($n1,$n2,$rv,$av,$a1,$a2) = split(/\t/,$line);
	push @values, $rv;
	$total++;
}
close CYTO;

my %Bins=%{Tools->binValues(1,@values)};
print "Source,Plasticity,Bin,Value\n";
for(my$i=-1;$i<=1;$i+=0.1){
	my $I = sprintf("%.1f",$i);
	if(defined($Bins{$I})){
		my $V = $Bins{$I}/$total;
		print $stan.",".$plas.",".$I.",".$V."\n";
	}else{
		print $stan.",".$plas.",".$I.",0\n";
	}
}
