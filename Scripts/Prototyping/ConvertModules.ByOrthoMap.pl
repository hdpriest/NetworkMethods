#!/usr/bin/perl
use warnings;
use strict;

use FindBin;
use lib "$FindBin::Bin/../../Lib";
use Tools;

my $map=$ARGV[0];
my $ModuleDir=$ARGV[1];
my $moduleOut=$ARGV[2];

die "usage: perl $0 <map file> <module input directory> <module output directory>\n\n" unless $#ARGV==2;

my %map;
my @map=@{Tools->LoadFile($map)};
foreach my $line (@map){
	my ($source,$target)=split(/\t/,$line);
	$map{$source}=$target;
}

my @files = grep {m/txt$/} @{Tools->LoadDir($ModuleDir)};
foreach my $file (@files){
	my $path=$ModuleDir."/".$file;
	my @content=@{Tools->LoadFile($path)};
	my $output=$moduleOut."/".$file;
	my @out;
	foreach my $line (@content){
		chomp $line;
		if(defined($map{$line})){
			push @out, $map{$line};
		}else{
			warn $line." has no map entry!\n";
		}
	}
	Tools->printToFile($output,\@out);
}
