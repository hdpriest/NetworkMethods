#!/usr/bin/perl
use warnings;
use strict;

use FindBin;
use lib "$FindBin::Bin/../../Lib";
use Tools;

my $usage = "perl $0 <directory of clusters> <output directory> <annotation file (gene<TAB>annotation)> \n\n";
die $usage unless $#ARGV==2;

my $root_in	= $ARGV[0];
my $root_out= $ARGV[1];

die "Can't find input directory : $root_in\n" unless -e $root_in;
die "Can't find output directory : $root_out\n" unless -e $root_out;
die "Can't find annotation file : $ARGV[2]\n" unless -e $ARGV[2];

my @Clusters = grep {m/txt$/} @{Tools->LoadDir($root_in)};
my @Annotation = @{Tools->LoadFile($ARGV[2])};

my %A;
foreach my $line (@Annotation){
	my ($id,$annotation)=split(/\t/,$line);
	$A{$id}=$annotation;
}

foreach my $file (@Clusters){
	my $path_in = $root_in."/".$file;
	my $path_out= $root_out."/".$file;
	my @output;
	warn "Processing $path_in\n";
	my @file=@{Tools->LoadFile($path_in)};
	my ($t,$a,$u)=(0,0,0);
	foreach my $line (@file){
		if(defined($A{$line})){
			push @output, $line."\t".$A{$line};
			$a++;
		}else{
			push @output, $line."\tnone";
			$u++;
		}
		$t++;
	}
	my $u_rat = $u/$t;
	$u_rat = sprintf("%.2f",$u_rat)*100;
	warn "$u_rat\% unannotated\n";
	Tools->printToFile($path_out,\@output);
}
