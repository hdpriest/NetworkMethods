#!/usr/bin/perl
use warnings;
use strict;

use lib '/home/hpriest/Scripts/Library';
use hdpTools;

my $usage="perl $0 <csv edge file> <set of nodes to extract edges between> <threshold for value>\n";
die $usage unless $#ARGV==2;

my @edges=@{hdpTools->LoadFile($ARGV[1])};

my %E;
foreach my $edge (@edges){
	my ($n1,$n2)=split(/\s/,$edge);
	if(defined($E{$n1})){
		$E{$n1}{$n2}=1;
	}else{
		$E{$n1}={};
		$E{$n1}{$n2}=1;
	}
	if(defined($E{$n2})){
		$E{$n2}{$n1}=1;
	}else{
		$E{$n2}={};
		$E{$n2}{$n1}=1;
	}
}

open(NET,"<",$ARGV[0]) || die "Cannot open $ARGV[0]!\n$!\nexiting...\n";
until(eof(NET)){
	my $line=<NET>;
	chomp $line;
	my @line=split(/\s/,$line);
	if((defined($E{$line[0]}{$line[1]}))||(defined($E{$line[1]}{$line[0]}))){
		print $line."\n" if abs($line[2]) >= $ARGV[2];
	}else{
	}
}
close NET;
