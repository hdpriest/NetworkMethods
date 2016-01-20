#!/usr/bin/perl
use warnings;
use strict;

use hdpTools;


my $Dir=$ARGV[0];
my $list=$ARGV[1];
die "usage: perl $0 <input directory> <gene list>\n\n" unless $#ARGV==1;

opendir(DIR,$Dir) || die "cannot open directory $Dir\n";
my @Files=readdir(DIR);
closedir DIR;

my @list=@{hdpTools->LoadFile($list)};
my %list;
map {$list{$_}="none"} @list;

my @Output;
foreach my $file (@Files){
	next unless (($file=~m/Cluster/i)||($file=~m/Module/i));
	my $Path=$Dir."/".$file;
	my $mod=$file;
	$mod=~s/\.txt//;
	$mod=~s/Cluster\.//;
	my @List=@{hdpTools->LoadFile($Path)};
	foreach my $gene (@List){	
		if(defined($list{$gene})){
			$list{$gene}=$mod;
		}else{
		#	my $output=$gene."\t"."none";
		#	push @Output, $output;
		}
	}
}
map {push @Output, "$_\t".$list{$_} } keys %list;
print join("\n",@Output)."\n";
