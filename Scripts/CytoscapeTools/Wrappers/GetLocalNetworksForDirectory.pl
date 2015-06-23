#!/usr/bin/perl
use warnings;
use strict;

use FindBin;
use lib "$FindBin::Bin/../../../Lib";
use Tools;

my $directory = $ARGV[0];
my $cytoscape = $ARGV[1];
my $cutoff = $ARGV[2];

my $usage = "Usage:
perl $0 <directory of gene sets (*.txt\$)> <cytoscape file to query> <cutoff for absolute edge strength>

This script wraps ../GetEdges.AllOfNodes.pl - for each gene set file in the input directory,
this script finds all edges above the provided edge value involving the provided nodes. Each
gene set, the local network is printed.

In an effort to avoid holding the whole network in perl-memory, this just iterates through the file for each gene set.

So, it is not fast, but also not a memory hog.

";

die $usage unless $#ARGV==2;

my $script = "$FindBin::Bin/../GetEdges.AllOfNodes.pl";

my @Dir = grep {m/\.txt$/} @{Tools->LoadDir($directory)};
foreach my $file (@Dir) {
	my $inPath = $directory."/".$file;
	my $out = $file;
	$out=~s/.txt/\.min$cutoff\.cyto\.tab/;
	$out = $directory."/".$out;
	my $cmd = $script." $cytoscape $inPath $cutoff > $out";
	`$cmd`;
}
