#!/usr/bin/perl
use warnings;
use strict;

my $usage = "perl $0 <cyto csv file>\n\n";

die $usage unless $#ARGV==0;

my $file	=$ARGV[0];
my %local;
my $max=0;
open(FILE,"<",$file) || die "cannot open $file!\n$!\nexiting...\n";
until(eof(FILE)){
	my $line=<FILE>;
	chomp $line;
	my @line=split(/\s/,$line);
	next if $line[0] eq $line[1];
	if(defined($local{$line[0]})){
		$local{$line[0]}+=abs($line[2]);
		$max=$local{$line[0]} if $local{$line[0]} > $max;
	}else{
		$local{$line[0]}=abs($line[2]);
	}
	if(defined($local{$line[1]})){
		$local{$line[1]}+=abs($line[2]);
		$max=$local{$line[1]} if $local{$line[1]} > $max;
	}else{
		$local{$line[1]}=abs($line[2]);
	}
}
close FILE;

my %kdist;
my $MK=0;
foreach my $key (keys %local){
	my $k = int($local{$key});
	my $K = $k/$max;
	$K = sprintf("%.2f",$K);
	$kdist{$K}+=1;
	$MK++;
}
for(my$i=0;$i<=1;$i+=0.01){
	my $I=sprintf("%.2f",$i);
	my $k=0;
	if(defined($kdist{$I})){
		$k=$kdist{$I};
	}else{
	}
	my $V=$k/$MK;
	print $I."\t".$V."\n";
}
