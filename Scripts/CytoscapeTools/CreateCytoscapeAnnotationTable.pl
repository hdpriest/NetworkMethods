#!/usr/bin/perl
use warnings;
use strict;

use FindBin;
use lib "$FindBin::Bin/../../Lib";
use Tools;

my $usage = "usage: perl $0 <cytoscape file> <gene list file> <name of column>\n\nThis script takes a list of genes, and creates a file of two fields, in which every gene in the network is assigned an interger value signifying if they are a member of that list, 1 or 0. This is useful for cytoscape graphical manipulation.\n\n";

die $usage unless $#ARGV==2;
my $cyto = $ARGV[0];
my $list = $ARGV[1];
my $name = $ARGV[2];

my %genes;
open(CYTO,"<",$cyto) || die "cannot open $cyto!\n$!\nexiting...\n";
until(eof(CYTO)){
	my $line=<CYTO>;
	chomp $line;
	my @line=split(/\t/,$line);
	$genes{$line[0]}=1;
	$genes{$line[1]}=1;
	$genes{$line[4]}=1;
	$genes{$line[5]}=1;
}
close CYTO;

my @list=@{Tools->LoadFile($list)};
my %list;
map {$list{$_}=1} @list;
print "Gene\t$name\n";
foreach my $gene (sort {$a cmp $b} keys %genes){
	if(defined($list{$gene})){
		print $gene."\t".1 ."\n";	
	}else{
		print $gene."\t".0 ."\n";	
	}
}

