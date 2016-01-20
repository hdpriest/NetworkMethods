#!/usr/bin/perl
use warnings;
use strict;

use lib '/home/hpriest/Scripts/Library';
use Tools;

my @list=@ARGV;
my %cuts;
print "Permutations,PositiveCut,NegativeCut\n";
foreach my $dir (@list){
	my $perm = $dir."/FDR_Calculations_Similarity.tab";
	$dir=~m/_p(\d+)_/;
	my $id=$1;
	my ($poscut,$negcut)=(0,0);
	($poscut,$negcut) = @{_parsePerms($perm)};
	my $type="perm";
	if(defined($cuts{$id})){
		push @{$cuts{$id}{pcut}}, $poscut;
		push @{$cuts{$id}{ncut}}, $negcut;
	}else{
		$cuts{$id}={};
		$cuts{$id}{pcut}=[];
		$cuts{$id}{ncut}=[];
		push @{$cuts{$id}{pcut}}, $poscut;
		push @{$cuts{$id}{ncut}}, $negcut;
	}
	#print $id.",".$poscut.",".$negcut.",".$type."\n";
}

foreach my $key (sort {$a <=> $b} keys %cuts){
	foreach my $ncut (@{$cuts{$key}{'ncut'}}){
		print $key.",negative,".$ncut."\n";
	}
	foreach my $pcut (@{$cuts{$key}{'pcut'}}){
		print $key.",positive,".$pcut."\n";
	}
}

sub _parsePerms {
	my $file=shift;
	my @file=@{Tools->LoadFile($file)};
	my $head = shift @file;
	my $pcut=2;
	my $ncut=-2;
	foreach my $line (@file){
		my @line=split(/\t/,$line);
		if($line[0]<0){
			next if $line[3] > 0.05;
			$ncut = $line[0] if $line[0] > $ncut;
		}else{
			next if $line[3] > 0.05;
			$pcut = $line[0] if $line[0] < $pcut;
		}
	}
	my @result=($pcut,$ncut);
	return \@result;
}

sub _parseSF {
	my $file=shift;
	my @file=@{Tools->LoadFile($file)};
	my $head = shift @file;
	my $cut=undef;
	foreach my $line (@file){
		my @line=split(/\t/,$line);
		if(($line[1] > 0.75) && ($line[2] < -0.8)){
			$cut=$line[0];
			last;
		}
	}
	return $cut;
}
