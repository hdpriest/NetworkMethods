#!/usr/bin/perl
use warnings;
use strict;

my $usage = "perl $0 <cyto csv file> <mask level>\n\n";

die $usage unless $#ARGV==1;

my $file	=$ARGV[0];
my $mask	=$ARGV[1];
my %local;
open(FILE,"<",$file) || die "cannot open $file!\n$!\nexiting...\n";
until(eof(FILE)){
	my $line=<FILE>;
	chomp $line;
	my @line=split(/\t/,$line);
	next if $line[0] eq $line[1];
	next if $line[3] < $mask;
	if(defined($local{$line[0]})){
		$local{$line[0]}{$line[1]}=1;
	}else{
		$local{$line[0]}={};
		$local{$line[0]}{$line[1]}=1;
	}
	if(defined($local{$line[1]})){
		$local{$line[1]}{$line[0]}=1;
	}else{
		$local{$line[1]}={};
		$local{$line[1]}{$line[0]}=1;
	}
}
close FILE;

my %kdist;

foreach my $key (keys %local){
	my $k = scalar(keys(%{$local{$key}}));
	if($k>2000){
		print $key."\n";
	}
	$kdist{$k}+=1;
}
foreach my $k (sort {$a <=> $b} keys %kdist){
#	print $k."\t".$kdist{$k}."\n";
}
