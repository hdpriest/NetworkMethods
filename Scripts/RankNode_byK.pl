#!/usr/bin/perl
use warnings;
use strict;

my $usage = "perl $0 <cyto csv file>\n\n";

die $usage unless $#ARGV==0;

my $file	=$ARGV[0];
my %local;
open(FILE,"<",$file) || die "cannot open $file!\n$!\nexiting...\n";
until(eof(FILE)){
	my $line=<FILE>;
	chomp $line;
	my @line=split(/\s/,$line);
	next if $line[0] eq $line[1];
	if(defined($local{$line[0]})){
		$local{$line[0]}+=abs($line[2]);
	}else{
		$local{$line[0]}=abs($line[2]);
	}
	if(defined($local{$line[1]})){
		$local{$line[1]}+=abs($line[2]);
	}else{
		$local{$line[1]}=abs($line[2]);
	}
	if(defined($local{$line[4]})){
		$local{$line[4]}+=abs($line[2]);
	}else{
		$local{$line[4]}=abs($line[2]);
	}
	if(defined($local{$line[5]})){
		$local{$line[5]}+=abs($line[2]);
	}else{
		$local{$line[5]}=abs($line[2]);
	}
}
close FILE;

my %kdist;
my %knode;
my $max=0;
foreach my $key (keys %local){
#	my $k = int($local{$key});
	my $k = $local{$key};
	$kdist{$k}+=1;
	if(defined($knode{$k})){
		push @{$knode{$k}}, $key;
	}else{
		$knode{$k}=[];
		push @{$knode{$k}}, $key;
	}
	$max++;
}
my $N=0;
foreach my $k (sort {$a <=> $b} keys %kdist){
	foreach my $n (@{$knode{$k}}){
		$N++;
		my $R=$max-$N+1;
		print $n."\t".$R."\t".$k."\n";
	}
}
