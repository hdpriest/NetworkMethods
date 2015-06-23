#!/usr/bin/perl
use warnings;
use strict;

use hdpTools;


my $Dir=$ARGV[0];
die "usage: perl $0 <input directory>\n\n" unless $#ARGV==0;

opendir(DIR,$Dir) || die "cannot open directory $Dir\n";
my @Files=grep {m/Cluster/} readdir(DIR);
closedir DIR;
my @Output;
foreach my $file (@Files){
	my $Path=$Dir."/".$file;
	my $mod=$file;
	$mod=~s/\.txt//;
	$mod=~s/Cluster\.//;
	my @List=@{hdpTools->LoadFile($Path)};
	foreach my $gene (@List){	
		my $output=$gene."\t".$mod;
		push @Output, $output;
	}
}
print "Gene\tModule\n";
print join("\n",@Output)."\n";
