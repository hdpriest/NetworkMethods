#!/usr/bin/perl
use warnings;
use strict;

use lib '/home/hpriest/Scripts/Library';
use hdpTools;
use Statistics::Descriptive;
die "usage: perl $0 <EXP file> <min FC to include gene profile> <space-del list of module gene lists>\n\n" unless $#ARGV>=2;
my $exp=shift @ARGV;
my $cutoff=shift @ARGV;
my @EXP=@{hdpTools->LoadFile($exp)};

my %Data;
my $Head=shift @EXP;
my @Head=split(/\t/,$Head);
my $prefix = shift @Head;
foreach my $line (@EXP){
	my @Line=split(/\t/,$line);
	die "Warning: $Line[0] in file twice\n\n" if defined ($Data{$Line[0]});
	$Data{$Line[0]}=\@Line;
}

print "GeneID,Sample,Expression\n";
foreach my $file (@ARGV){
	my @D=@{hdpTools->LoadFile($file)};
	my $ID=$file;
	$ID=~s/\.sorted\.list//;
	$ID=~s/\.txt//;
	$ID=~s/Module//i;
	my $X=0;
	my @Bins;
	for(my$i=1;$i<=$#Head;$i++){
		$Bins[$i]=[];
	}
	my $N=scalar(@D);
	foreach my $gene (@D){
		if(defined($Data{$gene})){
			my @Data=@{$Data{$gene}};
			my $id=shift @Data;
			my $FC = hdpTools->max(@Data)/hdpTools->min(@Data);
			next unless $FC > $cutoff;
			my @norm=@{hdpTools->averageNormalizeArray(\@Data)};
			for(my$i=0;$i<=$#norm;$i++){
				print $id.",".$Head[$i].",".$norm[$i]."\n";
			}
		}
	}
}


