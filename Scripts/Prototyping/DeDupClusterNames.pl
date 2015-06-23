#!/usr/bin/perl
use warnings;
use strict;
use lib '/home/hpriest/Scripts/Library';
use hdpTools;

my $dir=$ARGV[0];

die " usage: perl $0 <directory (./*.txt)>\n\n" unless defined $dir;

opendir(DIR,$dir) || die "cannot open $dir\n$!\nexiting...\n";
my @files = grep {m/txt$/} readdir(DIR);
closedir DIR;

foreach my $file (@files){
	my $path=$dir."/".$file;
	my @file=@{hdpTools->LoadFile($path)};
	my %ids;
	map{$ids{$_}=1} @file;
	my @output=keys %ids;
	hdpTools->printToFile($path,\@output);
}
