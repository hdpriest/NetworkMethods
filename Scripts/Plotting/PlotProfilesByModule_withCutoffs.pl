#!/usr/bin/perl
use warnings;
use strict;
use FindBin;
use lib "$FindBin::Bin/../../Lib";
use Tools;
use Configuration;
use Statistics::R;
use DataFrame;

die "usage: perl $0 <Merged Expression Frame> <min FC cutoff> <type: bar or line> <Directory of Modules> <output Directory>\n\n" unless $#ARGV==4;


my $Frame=$ARGV[0];
my $FC  =$ARGV[1];
my $type=$ARGV[2];
my $oDir=$ARGV[4];
my $iDir=$ARGV[3];

warn "Starting R...\n";
my $R=Statistics::R->new();
$R->startR();
$R->send("library(ggplot2)");

warn "Loading data...\n";
my $Data=DataFrame->new();
$Data->loadFile($Frame,"\t");
my @header=@{$Data->getHeader()};

my @DIR = grep {m/txt$/} @{Tools->LoadDir($iDir)};
foreach my $module (@DIR){
	my $ipath = $iDir."/".$module;
	my $opath = $oDir."/".$module;
	$opath=~s/txt//;
	die "Cannot find file: $ipath\n" unless -e $ipath;
	my @Genes = @{Tools->LoadFile($ipath)};
	my $RFrame= _genR_frame($Data,\@Genes,\@header,$opath,$FC);
	die "Could not create $RFrame.\n" unless -e $RFrame;
	_RFrameToPNG($RFrame,$opath,$type,\@header);	
}
warn "Done.\n";

exit(0);

sub _genR_frame {
	my $Frame = shift;
	my @Genes = @{$_[0]};
	my @header= @{$_[1]};
	my $pathRoot = $_[2];
	my $cut = $_[3];
	my $outFile = $pathRoot."frame.csv";
	my @output;
	push @output, "GeneID,Sample,Expression";
	foreach my $gene (@Genes){
		if($Frame->getDataByID($gene)){
			my @row=@{$Frame->getDataByID($gene)};
			my $min=Tools->min(@row);
			$min=1 if $min==0;
			my $max=Tools->max(@row);
			next unless (($max/$min)>$cut);
			my @mean = @{Tools->averageNormalizeArray(\@row)};
			my @sorted = sort {$a <=> $b} @mean;
			$sorted[0]=1 if $sorted[0]==0;
			next unless (($sorted[$#sorted]/$sorted[0])>$cut);
			for(my$i=0;$i<=$#row;$i++){
			#	push @output, $gene.",".$header[$i].",".$row[$i];
				push @output, $gene.",".$header[$i].",".$mean[$i];
			}
		}else{
		}
	}
	Tools->printToFile($outFile,\@output);
	return $outFile;
}

sub _RFrameToPNG {
	my $file =shift;
	my $path=shift;
	my $type=shift;
	my @header=@{$_[0]};
	my $xlab = "Sample";
	my $ylab = "Mean-normalized Expression";
	my $title=$path;
	$path.="png";
	my $cmd="DF=as.data.frame(read.table(\"$file\",sep=\",\",header=TRUE))";
	$R->send($cmd);
	$cmd="png(file=\"$path\")";
	$R->send($cmd);
	if($type eq "line"){
		$cmd="ggplot(data=DF,aes(x=Sample,y=Expression,group=GeneID)) + geom_line() + theme(axis.text.x = element_text(angle=90,hjust=1)) + ggtitle(\"$title\") + xlab(\"$xlab\") + ylab(\"$ylab\") + xlim(\"".join("\",\"",@header)."\")";
	}elsif($type eq "bar"){
		$cmd="ggplot(data=DF,aes(x=Sample,y=Expression,group=GeneID)) + geom_boxplot(stat = \"boxplot\", position = \"dodge\" ) + theme(axis.text.x = element_text(angle=90,hjust=1)) + ggtitle(\"$title\") + xlab(\"$xlab\") + ylab(\"$ylab\") + xlim(\"".join("\",\"",@header)."\")";
	}
	$R->send($cmd);
	my $result = $R->read();
	warn $result."\n";
	$cmd="dev.off()";
	$R->send($cmd);
	return 1;
}

