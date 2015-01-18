#!/usr/bin/perl
use warnings;
use strict;

use hdpTools;


my $Dir=$ARGV[0];
my $list=$ARGV[1];
die "usage: perl $0 <input directory> <gene list>\n\n" unless $#ARGV==1;

opendir(DIR,$Dir) || die "cannot open directory $Dir\n";
my @Files=grep {m/Cluster/} readdir(DIR);
closedir DIR;

my @list=@{hdpTools->LoadFile($list)};
my %list;
map {$list{$_}=1} @list;

my @Output;
foreach my $file (@Files){
	my $Path=$Dir."/".$file;
	my $mod=$file;
	$mod=~s/\.txt//;
	$mod=~s/Cluster\.//;
	my @List=@{hdpTools->LoadFile($Path)};
	foreach my $gene (@List){	
		next unless defined $list{$gene};
		my $output=$gene."\t".$mod;
		push @Output, $output;
	}
}

print join("\n",@Output)."\n";
