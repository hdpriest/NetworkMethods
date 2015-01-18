#!/usr/bin/perl
use warnings;
use strict;

use lib '/home/hpriest/Scripts/Library';
use hdpTools;
use WrapR;
use Statistics::Multtest qw(:all);
my $R=WrapR->new("GOrun.R.log");


my $file=$ARGV[0];
my $cut =$ARGV[1];
my $usage="perl $0 <file of permut> <cut>\n\n";

die $usage unless $#ARGV==1;

my @file=@{hdpTools->LoadFile($file)};
my @ids;
my @P;
my @obs;
my @mean;
my @stdev;
my @Z;
foreach my $line (@file){
	my @line=split(/\t/,$line);
	my $id =shift @line;
	my $num=shift @line;
	my $mean=hdpTools->mean(@line);
	my $stdev=hdpTools->stdev(@line);
	my $Z=($num-$mean)/$stdev;
	my $P = $R->ZtoP($Z);
	push @ids, $id;
	push @obs, $num;
	push @mean, $mean;
	push @stdev, $stdev;
	push @Z, $Z;
	push @P, $P;
}

my @FDR = @{BH(\@P)};
my $lt=0;
print "Obs\tmeean\tstdev\tZ\tP\tFDR\n";
for(my$i=0;$i<$#obs;$i++){
	my $line = $ids[$i]."\t".$obs[$i]."\t".$mean[$i]."\t".$stdev[$i]."\t".$Z[$i]."\t".$P[$i]."\t".$FDR[$i]."\n";
	print $line if $FDR[$i] <$cut;
	$lt++ if $FDR[$i] < $cut;
}
warn "total hits: $lt\n";
