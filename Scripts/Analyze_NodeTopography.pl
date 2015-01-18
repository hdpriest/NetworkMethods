#!/usr/bin/perl
use warnings;
use strict;

my $usage = "perl $0 <cyto csv file> <node> <mask level> <sig cutoff>\n\n";

die $usage unless $#ARGV==3;

my $file	=$ARGV[0];
my $node	=$ARGV[1];
my $mask	=$ARGV[2];
my $cutoff	=$ARGV[3];
my %local;
$local{all}={};
$local{sig}={};
open(FILE,"<",$file) || die "cannot open $file!\n$!\nexiting...\n";
until(eof(FILE)){
	my $line=<FILE>;
	chomp $line;
	my @line=split(/\t/,$line);
	next if $line[0] eq $line[1];
	next if $line[3] < $mask;
	next unless(($line[0] eq $node)||($line[1] eq $node));
	my @ts = ($line[0],$line[1]);
	@ts = sort {$a cmp $b} @ts;
	my $key = join("-",@ts);
	$local{all}{$key}=1;
	if($line[3] > $cutoff){
		$local{sig}{$key}=1;
	}
}
close FILE;

my $tot = scalar(keys %{$local{all}});
my $sig = scalar(keys %{$local{sig}});

print "$tot total nodes connected to $node\n";
print "$sig nodes with sig. plastic edges to $node\n";
