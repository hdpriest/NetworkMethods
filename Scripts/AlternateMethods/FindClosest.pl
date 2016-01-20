#!/usr/bin/perl
use warnings;
use strict;

use FindBin;
use lib "$FindBin::Bin/../../Lib";
use Tools;
use DataFrame;

my $inputFrame=$ARGV[0];
my $GOI = $ARGV[1];

die "usage: perl $0 <input data frame> <gene of interest>\n\n" unless $#ARGV==1;

my $DF = DataFrame->new();
$DF->loadFile($inputFrame,"\t");

my @goi_data = @{$DF->getDataByID($GOI)};
foreach my $id (@{$DF->getAllIDs()}){
	next if $id eq $GOI;
	my @target_data = @{$DF->getDataByID($id)};
	my $pcc = Tools->pearsonsR(\@goi_data,\@target_data);
	print $GOI."\t".$id."\t".$pcc."\n";
}
