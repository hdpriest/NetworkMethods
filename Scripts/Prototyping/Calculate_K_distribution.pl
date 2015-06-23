#!/usr/bin/perl
use warnings;
use strict;

my $frame = $ARGV[0];

die "need a file of adjacency values\n" unless -e $frame;

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
	my $sum = _sum(@line);
	my $rnd = sprintf("%d", $sum);
#	push @K, log($sum)/log(10);
	push @K, $rnd;
	$K{$rnd}++;
}
close KD;

foreach my $key (sort {$a <=> $b} keys %K){
	my $f=$K{$key}/$N;
	my $LK = log($key)/log(10);
	my $LF = log($f)/log(10);
	print $key."\t".$K{$key}."\t".$f."\t".$LK."\t".$LF."\n";
}

warn "Size of ids: ".scalar(@ID)."\n";
warn "Size of K: ".scalar(@K)."\n";
exit;
for (my$i=0;$i<=$#K;$i++){
	print $ID[$i]."\t".$K[$i]."\n";
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
