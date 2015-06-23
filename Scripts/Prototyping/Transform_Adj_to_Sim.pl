#!/usr/bin/perl
use warnings;
use strict;

my $frame = $ARGV[0];
my $alpha = $ARGV[1];
my $mu    = $ARGV[2];

die "usage: perl $0 <adjacency matrix> <alpha> <mu>\n" unless $#ARGV==2;

my %transform = %{_approx_Sigmoid($alpha,$mu)};

my $N = _findN($frame);

my @K;
my @ID;
my %K;
open(KD,"<",$frame) || die "Cannot open $frame!\n$!\nexiting...\n";
until(eof(KD)){
	my $line = <KD>;
	chomp $line;
	my @line = split(/\t/,$line);
	next if $line[0] eq "";
	my $id = shift @line;
	push @ID, $id;
}
close KD;

sub _approx_Sigmoid {
	my $alpha = shift;
	my $mu = shift;
	my @Adjacencies;
	for(my $i=-1;$i<=1;$i+=0.001){
		my $I=sprintf("%.3f",$i);
		my $Adj = 1 / (1 + exp($alpha * -1 * (abs($I) - $mu)));
		$Adj = sprintf("%.3f",$Adj);
		print $I."\t".$Adj."\n";
		#(1.0f / (1 + Math.exp(A * -1 * (Math.abs(V) - M))))
	}
	exit;
}


sub _sum {
	my @array=@_;
	my $s=0;
	map {$s+=$_} @array;
	return $s;
}

sub _findN {
	my $file=shift;
	my $n=0;
	open(NH,"<",$file) || die "Cannot find $file!\n$!\nexiting...\n";
	until(eof(NH)){
		my $line=<NH>;
		$n++;
	}
	close NH;
	$n--; # has a header;
	return $n;
}
