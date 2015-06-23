#!/usr/bin/perl
use warnings;
use strict;

use lib '/home/hpriest/Scripts/Library';
use hdpTools;

my $dir=$ARGV[0];

die "usage: perl $0 <input directory>\n\n" unless $#ARGV==0;

my @files = @{hdpTools->LoadDir($dir)};
my %g;
foreach my $file (@files){
	my $path=$dir."/".$file;
	my @file = @{hdpTools->LoadFile($path)};
	foreach my $line (@file){
		my ($GO,$P,$FDR,$desc)=split(/\t/,$line);
		if($FDR<0.05){
			$g{$GO}=1;
		}
	}
}

my $N=scalar(keys(%g));
print $N." unique go terms\n";
