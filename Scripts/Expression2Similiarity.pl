#!/usr/bin/perl
use warnings;
use strict;
use FindBin;
use lib "$FindBin::Bin/../Lib";
use DataFrame;
use Tools;
use NetworkTools;
use Statistics::R;
use Cwd;
my $R=Statistics::R->new("log_dir" => cwd(),"tmp_dir" => cwd());
$R->stop();
$R->start();
$R->send("library(dynamicTreeCut)");
$R->send("library(WGCNA)");

my $usage=" perl $0 <Expression Frame For analysis> <Output file root for similarity Matrix>\n\n";
die $usage unless $#ARGV==1;

my $expFrame=$ARGV[0];
my $outFile =$ARGV[1];

my $DF=DataFrame->new();
$DF->loadFile($expFrame,"\t");

my $rf = NetworkTools->sendFileToRFrame($R,$expFrame,"DF","\t");
my $command="TDF=t(DF)";
$R->send($command);
$command="Sim=cor(TDF)";
$R->send($command);
$command="Sim";
$R->send($command);
my @header=@{NetworkTools->getMatrixHeader($R,"Sim")};
print "GeneID\t".join("\t",@header)."\n";
my $r=1;
while(my $ID=$DF->getThisID()){
	my ($id,$headRef,$valRef)=NetworkTools->getMatrixRow($R,"Sim",$r);
	my @values=@$valRef;
	print $id."\t".join("\t",@values)."\n";
	$r++;
	$DF->iterateRow();
}

#NetworkTools->saveRobjToAscii($R,"DF",$outFile.".throughReadTable.R");

exit;
