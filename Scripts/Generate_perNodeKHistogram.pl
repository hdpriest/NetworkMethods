#!/usr/bin/perl
use warnings;
use strict;

use FindBin;
use lib "$FindBin::Bin/../Lib";
use Tools;

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
		push @{$local{$line[0]}}, $line[2];
	}else{
		$local{$line[0]}=[];
		push @{$local{$line[0]}}, $line[2];
	}
	if(defined($local{$line[1]})){
		push @{$local{$line[1]}}, $line[2];
	}else{
		$local{$line[1]}=[];
		push @{$local{$line[1]}}, $line[2];
	}
}
close FILE;


foreach my $key (sort {$a cmp $b} keys %local){
	my @k_values = @{$local{$key}};
	my %histo = %{Tools->binValues(2,@k_values)};
	my @values;
	for(my $i=-1;$i<=1;$i+=0.01){
		my $I=sprintf("%.2f",$i);
		if(defined($histo{$I})){
			push @values, $histo{$I};
		}else{
			push @values, 0;
		}
	}
	print $key."\t".join(",",@values)."\n";
}
