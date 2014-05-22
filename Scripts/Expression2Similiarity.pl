#!/usr/bin/perl
use warnings;
use strict;
use FindBin;
use lib "$FindBin::Bin/../Lib";
use DataFrame;
use Tools;
use NetworkTools;
use Statistics::R;

my $R=Statistics::R->new();
$R->stop();
$R->start();
$R->send("library(dynamicTreeCut)");
$R->send("library(WGCNA");

my $usage=" perl $0 <Expression Frame For analysis> <Output file for similarity Matrix>\n\n";
die $usage unless $#ARGV==1;

my $expFrame=$ARGV[0];
my $outFile =$ARGV[1];

my $g = $NetworkTools->sendFileToRFrame($R,$expFrame,"DF","\t");

exit;
