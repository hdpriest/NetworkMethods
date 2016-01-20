#!/usr/bin/perl
use warnings;
use strict;

my $usage = "perl $0 <cyto csv file> <file2>\n\n";

die $usage unless $#ARGV==1;

my $file=$ARGV[0];
my $file2=$ARGV[1];
my ($F,$S,$B)=(0,0,0);
my %D;
open(FILE,"<",$file) || die "cannot open $file!\n$!\nexiting...\n";
until(eof(FILE)){
	my $line=<FILE>;
	next if $line=~m/weight/i;
	chomp $line;
	my @line=split(/\s/,$line);
	next if $line[0] eq $line[1];
	my @k = ($line[0],$line[1]);
	@k = sort {$a cmp $b} @k;
	my $k = join("-",@k);
	$D{$k}=1;
}
close FILE;
open(FILE,"<",$file2) || die "cannot open $file2!\n$!\nexiting...\n";
until(eof(FILE)){
	my $line=<FILE>;
	next if $line=~m/weight/i;
	chomp $line;
	my @line=split(/\s/,$line);
	next if $line[0] eq $line[1];
	my @k = ($line[0],$line[1]);
	@k = sort {$a cmp $b} @k;
	my $k = join("-",@k);
	if(defined($D{$k})){
		$B++;
		delete $D{$k};
	}else{
		$S++;
	}
}

foreach my $k (keys %D){
	$F++;
}

warn "First: $F\nSecond: $S\nBoth: $B\n";
