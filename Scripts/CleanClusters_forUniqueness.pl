#!/usr/bin/perl
use warnings;
use strict;

use lib '/home/hpriest/Scripts/Library';
use hdpTools;

my $InputDir=$ARGV[0];
my $outputDir=$ARGV[1];
die "usage: perl $0 <Directory of input file lists> <output dir>\n\nSometimes, it can be useful to repeat genes in the network input file (to allow for multiple orthology, or similar applications).\nIn these cases, some modules may contain the same gene twice, and must be \'cleaned\'.\n" unless $#ARGV==1;

opendir(DIR,$InputDir);
my @Files=grep {m/txt$/} readdir(DIR);
closedir DIR;

foreach my $file (@Files){
	my $path=$InputDir."/".$file;
	my $pathOut=$outputDir."/".$file;
	my @file=@{hdpTools->LoadFile($path)};
	my %IDs;
	for(my$i=0;$i<=$#file;$i++){
		$IDs{$file[$i]}=1;
	}
	my @ids = keys %IDs;
	hdpTools->printToFile($pathOut,\@ids);
}



