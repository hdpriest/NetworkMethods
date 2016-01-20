#!/usr/bin/perl
use warnings;
use strict;

use FindBin;
use lib "$FindBin::Bin/../../Lib";
use Cwd;
use Tools;
use DataFrame;
use Configuration;
our $cwd = getcwd;

die "usage: perl $0 <configuration file>\n\n" unless $#ARGV==0;

our $Config=Configuration->new($ARGV[0]);
_checkConfig($Config);
my $Prefix = $Config->get("OPTIONS","SQL_Prefix");
our $NetDir = $Config->get("OPTIONS","NetworkDir");
my $inputExpression = $NetDir."/InputExpression.matrix.tab";
my $phytozomeAnnotation = $Config->get("OPTIONS","PhytozomeAnnotation");
my $anno_id_col  = $Config->get("OPTIONS","Annotation_id_column");
my $anno_name_col= $Config->get("OPTIONS","Annotation_name_column");
my $anno_desc_col= $Config->get("OPTIONS","Annotation_description_column");
my $Adj		= $NetDir."/Adjacency.cytoscape.raw.tab";
my $ClusterDir	= $NetDir."/Clusters";
my $Dbase	= $Config->get("OPTIONS","DBName");
my $Mask	= $Config->get("OPTIONS","MaskLevel");
my $regex = $Config->get("OPTIONS","RegularExpression");
my $displayName = $Config->get("OPTIONS","DisplayName");
my @Exp = @{Tools->LoadFile($inputExpression)};
my @ExpFields = split(/\t/,$Exp[0]);
my $geneSpace = shift @ExpFields;
@Exp=[];
my @Genes  = sort {$a cmp $b} @{Tools->LoadFile($Config->get("OPTIONS","MasterGeneList"))};
my %Genes;
for(my$i=0;$i<=$#Genes;$i++){
	my $I=$i+1;
	$Genes{$Genes[$i]}=$I;
#### for unique universal gene id
}
warn "Making Create Statements...\n";
my $Create = _GetCreateStatement($Dbase,$Prefix);
my $Statement = _createExpressionTable($Dbase,$Prefix,@ExpFields);
$Statement .= _createMetricTable($Dbase,$Prefix,@ExpFields);
$Statement .= _prepClusterTable($Dbase,$Prefix,$ClusterDir);

print $Create."\n";
print $Statement."\n";

warn "Creating Gene Load Inserts...\n";
my $GeneLoadCommand = _loadGenes($Dbase,$Prefix,$regex,$displayName,\%Genes);
print $GeneLoadCommand."\n";
#print "DO SLEEP(5);\n";
warn "Creating Expression Load Inserts...\n";
my $ExpressionLoadCommand = _loadExpression($Dbase,$Prefix,$inputExpression,\@ExpFields,\%Genes);
print $ExpressionLoadCommand."\n";
warn "Creating Annotation Load Inserts...\n";
my $AnnotationLoadCommand = _loadAnnotation($Dbase,$Prefix,$phytozomeAnnotation,$anno_id_col,$anno_name_col,$anno_desc_col,\%Genes);
print $AnnotationLoadCommand."\n";
#print "DO SLEEP(5);\n";
warn "Creating Cluster Load Inserts...\n";
my $ClusterLoadCommand = _loadClusters($Dbase,$Prefix,$ClusterDir,\%Genes);
print $ClusterLoadCommand."\n";
#print "DO SLEEP(5);\n";
warn "Creating Adjacency Load Inserts...\n";
my $AdjacencyLoadCommand = _loadAdjacency($Dbase,$Prefix,$Adj,$Mask,\%Genes);
print $AdjacencyLoadCommand."\n";
#print "DO SLEEP(5);\n";
warn "Creating Metrics Load Inserts...\n";
my $MetricLoadCommand = _prepNetworkMetrics(\%Genes,$inputExpression,$Dbase,$Prefix,$Adj,$ClusterDir,$Mask);
print $MetricLoadCommand."\n";
warn "Done.\n";

sub _loadAdjacency {
	my $DB = shift;
	my $Prefix = shift;
	my $AdjFile = shift;
	my $mask = shift;
	my %genes=%{$_[0]};
	open(Adj,"<",$AdjFile) || die "Cannot open $AdjFile!\n$!\nexiting...\n";
	my $edgeCSV = "$cwd/$Prefix.EdgeLoad.csv";
	my $adjCSV  = "$cwd/$Prefix.AdjLoad.csv";
	my $edge_num=1;
	my @EdgeCSV;
	my @AdjCSV;
	my $adjCols = "(id,edgeId)";
	my $edgeCols = "(edgeId,adjacency,node1,node2)";
	until(eof(Adj)){
		my $line=<Adj>;
		chomp $line;
		my ($id_A,$id_B,$edge,$absEdge,$alt_A,$alt_B)=split(/\t/,$line);
		next if abs($edge) < $mask;
		$edge=sprintf("%.3f",$edge);
		die "$id_A undefined in master gene list\n" unless defined $genes{$id_A};
		die "$id_B undefined in master gene list\n" unless defined $genes{$id_B};
		my $ida = $genes{$id_A};
		my $idb = $genes{$id_B};
		push @AdjCSV, $ida.",".$edge_num;
		push @AdjCSV, $idb.",".$edge_num;
		push @EdgeCSV, "$edge_num,$edge,$ida,$idb";
		$edge_num++;
		if($edge_num==10000){
			#warn "Debugging bypass in place. Stopping at edge #: 10000\n";
			#last if $edge_num==10000;
		}
	}
	Tools->printToFile($edgeCSV,\@EdgeCSV);
	Tools->printToFile($adjCSV,\@AdjCSV);
	my $insert = "LOAD DATA INFILE '$edgeCSV' INTO TABLE `$Dbase`.`$Prefix"."_Edges` fields terminated by ',' enclosed by '' lines terminated by '\\n' $edgeCols;\n";
	$insert .= "LOAD DATA INFILE '$adjCSV' INTO TABLE `$Dbase`.`$Prefix"."_Adjacency` fields terminated by ',' enclosed by '' lines terminated by '\\n' $adjCols;";
	close Adj;
	return $insert;
}

sub _PrepGenes {
	my $file=shift;
	my $Expression=shift;
	my @Genes = keys %{$file};
	my %Genes;
	foreach my $gene (@Genes){
		$Genes{$gene}={};
		$Genes{$gene}{modular_k}=0;
		$Genes{$gene}{modular_k_rank}=0;
		$Genes{$gene}{modular_exp_rank_mean}=0;
		$Genes{$gene}{module}=0;
		$Genes{$gene}{k}=0;
		$Genes{$gene}{k_rank}=0;
		$Genes{$gene}{exp_mean}=0;
		$Genes{$gene}{exp_rank_mean}=0;
	}
	my $DataFrame = DataFrame->new();
	$DataFrame->loadFile($Expression,"\t");
	$DataFrame->initRowIterator();
	foreach my $gene (@Genes){
		my @Data = @{$DataFrame->getDataByID_withZeros($gene)};
		my $mean = Tools->mean(@Data);
		$Genes{$gene}{exp_mean}=$mean;
	}
	return \%Genes;
}

sub _PrepModules {
	my $ClusterDir = shift;
	my $ref = shift;
	my %Genes=%$ref;
	my %Modules;
	my %used; ## doubly ensures single-cluster membership
	my @Files = grep {m/txt$/} @{Tools->LoadDir($ClusterDir)};
	foreach my $file (@Files){
		my $path=$ClusterDir."/".$file;
		$file=~m/Cluster\.(\d+)\.txt/;
		my $c=$1;
		$Modules{$c}=[];
		my @file=@{Tools->LoadFile($path)};
		foreach my $id (@file){
			next if defined $used{$id};
			die "Cannot find $id in gene hash!\n" unless defined $Genes{$id};
			$Genes{$id}{module}=$c;
			push @{$Modules{$c}}, $id;
			$used{$id}=1;
		}
	}
	my @refs = (\%Genes,\%Modules);
	return \@refs;
}

sub _prepNetworkMetrics {
	my $MasterList = shift;
	my $Expression = shift;
	my $DB = shift;
	my $Prefix = shift;
	my $AdjPath= shift;
	my $ClusterDir = shift;
	my $Mask = shift;
	my $outFile ="$cwd/".$Prefix.".Metrics.csv";
	my @output;
	my %g = %{$MasterList};
	my @Genes = keys %g;
	my %Genes = %{_PrepGenes($MasterList,$Expression)};
	my ($r1,$r2) = @{_PrepModules($ClusterDir,\%Genes)};
	%Genes = %$r1;
	my %Modules = %$r2;
	%Genes = %{_PrepAdj($AdjPath,\%Genes,$Mask)};
	%Genes = %{_Rank(\%Genes,\%Modules)};
	foreach my $gene (keys %Genes){
		my @line;
		die "$gene not defined in master list\n" unless defined $g{$gene};
		push @line, $g{$gene};
		push @line, sprintf("%.3f",$Genes{$gene}{modular_k});
		push @line, $Genes{$gene}{modular_k_rank};
		push @line, sprintf("%.3f",$Genes{$gene}{modular_exp_rank_mean});
		push @line, $Genes{$gene}{module};
		push @line, sprintf("%.3f",$Genes{$gene}{k});
		push @line, $Genes{$gene}{k_rank};
		push @line, sprintf("%.3f",$Genes{$gene}{exp_mean});
		push @line, $Genes{$gene}{exp_rank_mean};
		my $line = join(",",@line);
		push @output, $line;
	}
	Tools->printToFile($outFile,\@output);
	my $columns = "(id,modular_k,modular_k_rank,modular_mean_exp_rank,module,k,k_rank,mean_exp,mean_exp_rank)";
	my $insert = "LOAD DATA INFILE '$outFile' INTO TABLE `$Dbase`.`$Prefix"."_Metrics` fields terminated by ',' enclosed by '' lines terminated by '\\n' $columns;";
	return $insert;
}

sub _Rank {
	my %Genes=%{$_[0]};
	my %Modules=%{$_[1]};
	foreach my $module (keys %Modules){
		my @genes = @{$Modules{$module}};
		my $G = scalar(@genes);
		my %this_module;
		for(my$i=0;$i<=$#genes;$i++){
			$this_module{$genes[$i]}=$Genes{$genes[$i]}{modular_k};
		}
		foreach my $gene (sort {$this_module{$a} <=> $this_module{$b}} keys %this_module){ ## low to high sort
			my $rank = $G;
			$Genes{$gene}{modular_k_rank}=$rank;
			$G--;
		}
		for(my$i=0;$i<=$#genes;$i++){
			$this_module{$genes[$i]}=$Genes{$genes[$i]}{exp_mean};
		}
		$G = scalar(@genes);
		foreach my $gene (sort {$this_module{$a} <=> $this_module{$b}} keys %this_module){ ## low to high sort
			my $rank = $G;
			$Genes{$gene}{modular_exp_rank_mean}=$rank;
			$G--;
		}
		
	}
	my $N = scalar(keys(%Genes));
	foreach my $gene (sort {$Genes{$a}{k} <=> $Genes{$b}{k}} keys %Genes){
		my $rank = $N;
		$Genes{$gene}{k_rank}=$rank;
		$N--;
	}
	$N = scalar(keys(%Genes));
	foreach my $gene (sort {$Genes{$a}{exp_mean} <=> $Genes{$b}{exp_mean}} keys %Genes){
		my $rank = $N;
		$Genes{$gene}{exp_rank_mean}=$rank;
		$N--;
	}
	return \%Genes;
}

sub _PrepAdj {
	my $AdjPath = shift;
	my %Genes = %{$_[0]};
	my $Mask = $_[1];
	open(ADJ,"<",$AdjPath) || die "Cannot open $AdjPath!\n$!\nexiting...\n";
	until(eof(ADJ)){
		my $line=<ADJ>;
		chomp $line;
		my ($IDA,$IDB,$edge,$absEdge,$alt_IDA,$alt_IDB)=split(/\s/,$line);
		next if $absEdge <= $Mask;
		die "Cannot find $IDA in Gene table!\n" unless defined $Genes{$IDA};
		die "Cannot find $IDB in Gene table!\n" unless defined $Genes{$IDB};
		if(($Genes{$IDA}{module} == $Genes{$IDB}{module}) && ($Genes{$IDA}{module}!=0)){
			$Genes{$IDA}{modular_k}+=$absEdge;
			$Genes{$IDB}{modular_k}+=$absEdge;
		}else{
		}
		$Genes{$IDA}{k}+=$absEdge;
		$Genes{$IDB}{k}+=$absEdge;
	}
	close ADJ;
	return \%Genes;
}

sub _loadClusters {
	my $Dbase=shift;
	my $prefix=shift;
	my $ClusterDir = shift;
	my %Genes = %{$_[0]};
	my @Files = grep {m/txt$/} @{Tools->LoadDir($ClusterDir)};
	my %used;
	my @output;
	my $outFile = "$cwd/".$prefix.".clusterLoad.csv";
	foreach my $file (@Files){
		my $path=$ClusterDir."/".$file;
		$file=~m/Cluster\.(\d+)\.txt/;
		my $c=$1;
		my @IDs = @{Tools->LoadFile($path)};
		foreach my $gene (@IDs){
			die "gene $gene not in master list\n"  unless defined $Genes{$gene};
			next if defined($used{$gene});
			my $id = $Genes{$gene};
			my $line = "$id,$c";
			push @output, $line;
			$used{$gene}=1;
		}
	}
	Tools->printToFile($outFile,\@output);
	my $statement = "LOAD DATA INFILE '$outFile' INTO TABLE `$Dbase`.`$prefix"."_Clusters_Genes` fields terminated by ',' enclosed by '' lines terminated by '\\n';";
	return $statement;
}

sub _loadAnnotation {
#	_loadAnnotation($Dbase,$Prefix,$phytozomeAnnotation,\%Genes);
	my $DB = shift;
	my $Prefix = shift;
	my $phyto = shift;
	my $id_col=shift;
	my $name_col=shift;
	my $desc_col=shift;
	my %Genes=%{$_[0]};
	my @input = @{Tools->LoadFile($phyto)};
	my %A;
	foreach my $line (@input){
		my @line=split(/\t/,$line);
		my $id = $line[$id_col];
		my $name = $line[$name_col];
		$name = "None" unless defined $name;
		my $desc = $line[$desc_col];
		$desc = "None" unless defined $desc; 
		if(defined($A{$id})){
			if($name eq $A{$id}{'name'}){
			}elsif($A{$id}{'desc'} eq $desc){
			}else{
	#			warn "Multiple different deflines for $id found\n";
			}
		}else{
			$A{$id}={};
			$A{$id}{'name'}=$name;
			$A{$id}{'desc'}=$desc;
		}
	}
	my $statement="";
	foreach my $gene (keys %Genes){
		my $id = $Genes{$gene};
		my $name = "None";
		my $desc = "None";
		if(defined($A{$gene})){
			$name=$A{$gene}{'name'};
			$desc=$A{$gene}{'desc'};
		}
		$name = "None" if $name eq "";
		$name=~s/\,/ /g;
		$name=~s/\'//g;
		$desc=~s/\,/ /g;
		$desc=~s/\'//g;
		my $insert = "INSERT INTO `$DB`.`$Prefix"."_Annotation` (`id`,`locus`,`name`,`description`) VALUES ($id,\'$gene\',\'$name\',\'$desc\');\n";
		$statement.=$insert;
	}
	return $statement;
}

sub _loadExpression {
	my $DB = shift;
	my $Prefix = shift;
	my $DF = shift;
	my @ExpFields=@{$_[0]};
	my %Genes=%{$_[1]};
	my $DataFrame = DataFrame->new();
	$DataFrame->loadFile($DF,"\t");
	$DataFrame->initRowIterator();
	my $statement="";
	foreach my $field (@ExpFields){
		my $biogrp = $Config->get("EXPRESSION_FIELDS",$field);
		my $display= $Config->get("DISPLAY_NAMES",$biogrp);
		my $insert = "INSERT INTO `$DB`.`$Prefix"."_Expression_Map` (`biosample_id`,`field_name`,`display_name`) VALUES ($biogrp,\'$field\',\'$display\');\n";
		$statement.=$insert;
	}
	foreach my $gene (keys %Genes){
		if($DataFrame->checkID($gene)){
			my @Data = @{$DataFrame->getDataByID_withZeros($gene)};
			my $id = $Genes{$gene};
			my $headers="`id`";
			my $values="$id";
			for(my $i=0;$i<=$#Data;$i++){
				$headers.=", `$ExpFields[$i]`";
				my $v = sprintf("%.3f",$Data[$i]);
				$values .=", $v";
			}
			my $insert = "INSERT INTO `$DB`.`$Prefix"."_Expression` ($headers) VALUES ($values);\n";
			$statement.=$insert;
		}else{
		}
	}
	return $statement;
}

sub _getMapStatement{
	my $Dbase=shift;
	my $regex=shift;
	my $Prefix=shift;
	my $dispname=shift;
	my $statement="INSERT INTO `$Dbase`.`MapReference` (`prefix`,`display_name`,`regex`) VALUES ('$Prefix','$dispname','$regex');\n";
	return $statement;
}

sub _loadGenes {
	my $Dbase = shift;
	my $Prefix= shift;
	my $regex = shift;
	my $displayName = shift;
	my %Genes = %{$_[0]};
	my $File = "$cwd/".$Prefix.".geneLoad.csv";
	my $mapStatement = _getMapStatement($Dbase,$regex,$Prefix,$displayName);
	my @output;
	foreach my $key (keys %Genes){
		my $line="$Genes{$key},$key";
		push @output, $line;
	}
	Tools->printToFile($File,\@output);
	my $insert = "LOAD DATA INFILE '$File' INTO TABLE `$Dbase`.`$Prefix"."_Genes` fields terminated by ',' enclosed by '' lines terminated by '\\n';";
	$insert .= "\n".$mapStatement;
	return $insert;
}

sub _prepClusterTable {
	my $Dbase=shift;
	my $prefix=shift;
	my $ClusterDir = shift;
	my @Files = grep {m/txt$/} @{Tools->LoadDir($ClusterDir)};
	my %used;
	my $statement="";
	foreach my $file (@Files){
		my $path=$ClusterDir."/".$file;
		$file=~m/Cluster\.(\d+)\.txt/;
		my $c=$1;
		$statement .= "INSERT INTO `$Dbase`.`$prefix"."_Clusters` (`clusterID`) VALUES ($c);\n";
	}
	return $statement;
}

sub _createMetricTable {
	my $DB = shift;
	my $PF = shift;
	my $table = $PF."_Metrics";
	my $Create ="
use $DB;
CREATE TABLE IF NOT EXISTS `$table` (
`id` INT NOT NULL,
`modular_k` DOUBLE NOT NULL,
`modular_k_rank` INT(8) NOT NULL,
`modular_mean_exp_rank` INT(8) NOT NULL,
`module` INT(8) NOT NULL,
`k` DOUBLE NOT NULL,
`k_rank` INT(8) NOT NULL,
`mean_exp` DOUBLE NOT NULL,
`mean_exp_rank` INT(8) NOT NULL,
PRIMARY KEY (`id`),
UNIQUE INDEX `id_UNIQUE` (`id` ASC),
CONSTRAINT `gene_metrics`
	FOREIGN KEY (`id`)
	REFERENCES `$DB`.`$PF"."_Genes` (`id`)
	ON DELETE NO ACTION
	ON UPDATE NO ACTION)
ENGINE = MyISAM;
";
	return $Create;
}

sub _createExpressionTable {
	my $dbase=shift;
	my $Prefix = shift;
	my $exptable = $Prefix."_Expression";
	my $gentable = $Prefix."_Gene";
	my @Fields=@_;
	my $statement="
-- -----------------------------------------------------
-- Table `$dbase`.`$Prefix"."_Expression`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `$dbase`.`$Prefix"."_Expression` (
  `id` INT NOT NULL,";
  	foreach my $field (@Fields){
  		$statement.=" `$field` DOUBLE NOT NULL,\n";
	}
	$statement.="
  PRIMARY KEY (`id`),
  UNIQUE INDEX `id_UNIQUE` (`id` ASC),
  CONSTRAINT `gene_exp`
    FOREIGN KEY (`id`)
    REFERENCES `$dbase`.`$Prefix"."_Genes` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = MyISAM;

-- -----------------------------------------------------
-- Table `$dbase`.`$Prefix"."_Expression_Map`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `$dbase`.`$Prefix"."_Expression_Map` (
  `biosample_id` INT NOT NULL,
  `field_name` VARCHAR(45) NOT NULL,
  `display_name` VARCHAR(45) NOT NULL,
  INDEX `bioid_idx` (`biosample_id` ASC),
  INDEX `fieldname_idx` (`field_name` ASC))
ENGINE = MyISAM;
";
	return $statement;
}

sub _GetCreateStatement {
	my $dbase = shift;
	my $prefix= shift;
	my $monster ="
-- MySQL Script generated by MySQL Workbench
-- 05/18/15 18:43:58
-- Model: New Model    Version: 1.0
-- MySQL Workbench Forward Engineering

SET \@OLD_UNIQUE_CHECKS=@\@UNIQUE_CHECKS, UNIQUE_CHECKS=0;
SET \@OLD_FOREIGN_KEY_CHECKS=@\@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0;
SET \@OLD_SQL_MODE=@\@SQL_MODE, SQL_MODE='TRADITIONAL,ALLOW_INVALID_DATES';

-- -----------------------------------------------------
-- Schema mydb
-- -----------------------------------------------------

-- -----------------------------------------------------
-- Schema mydb
-- -----------------------------------------------------
CREATE SCHEMA IF NOT EXISTS `$dbase` DEFAULT CHARACTER SET utf8 COLLATE utf8_general_ci ;
USE `$dbase` ;

-- -----------------------------------------------------
-- Table `$dbase`.`$Prefix"."_Genes`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `$dbase`.`$Prefix"."_Genes` (
  `id` INT NOT NULL,
  `name` VARCHAR(45) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE INDEX `gene.id_UNIQUE` (`id` ASC),
  UNIQUE INDEX `name_UNIQUE` (`name` ASC))
ENGINE = MyISAM;

-- -----------------------------------------------------
-- Table `$dbase`.`MapReference`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `$dbase`.`MapReference` (
  `prefix` VARCHAR(45) NOT NULL,
  `display_name` VARCHAR(45) NOT NULL,
  `regex` VARCHAR(45) NOT NULL,
  INDEX `prefix_idx` (`prefix` ASC),
  INDEX `name_idx` (`display_name` ASC),
  INDEX `regex_idx` (`regex` ASC))
ENGINE = MyISAM;

-- -----------------------------------------------------
-- Table `$dbase`.`$Prefix"."_Clusters`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `$dbase`.`$Prefix"."_Clusters` (
  `clusterID` INT NOT NULL,
  PRIMARY KEY (`clusterID`),
  UNIQUE INDEX `clusterID_UNIQUE` (`clusterID` ASC))
ENGINE = MyISAM;


-- -----------------------------------------------------
-- Table `$dbase`.`$Prefix"."_Edges`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `$dbase`.`$Prefix"."_Edges` (
  `edgeId` INT NOT NULL,
  `adjacency` DOUBLE NOT NULL,
  `node1` VARCHAR(45) NOT NULL,
  `node2` VARCHAR(45) NOT NULL,
  PRIMARY KEY (`edgeId`),
  UNIQUE INDEX `edgeId_UNIQUE` (`edgeId` ASC))
ENGINE = MyISAM;


-- -----------------------------------------------------
-- Table `$dbase`.`$Prefix"."_Adjacency`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `$dbase`.`$Prefix"."_Adjacency` (
  `id` INT NOT NULL,
  `edgeId` INT NOT NULL,
  INDEX `id_idx` (`id`),
  INDEX `edge_idx` (`edgeId` ASC),
  CONSTRAINT `id`
    FOREIGN KEY (`id`)
    REFERENCES `$dbase`.`$Prefix"."_Genes` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `edge`
    FOREIGN KEY (`edgeId`)
    REFERENCES `$dbase`.`$Prefix"."_Edges` (`edgeId`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = MyISAM;


-- -----------------------------------------------------
-- Table `$dbase`.`$Prefix"."_Clusters_Genes`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `$dbase`.`$Prefix"."_Clusters_Genes` (
  `id` INT NOT NULL,
  `Cluster` INT NOT NULL,
  PRIMARY KEY (`id`),
  INDEX `ind_cluster` (`Cluster` ASC),
  CONSTRAINT `gene`
    FOREIGN KEY (`id`)
    REFERENCES `$dbase`.`$Prefix"."_Genes` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `cluster`
    FOREIGN KEY (`Cluster`)
    REFERENCES `$dbase`.`$Prefix"."_Clusters` (`clusterID`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = MyISAM;




-- -----------------------------------------------------
-- Table `$dbase`.`$Prefix"."_Annotation`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `$dbase`.`$Prefix"."_Annotation` (
  `id` INT NOT NULL,
  `locus` VARCHAR(45) NOT NULL,
  `name` VARCHAR(45) NOT NULL,
  `description` VARCHAR(4096) NOT NULL, 
  PRIMARY KEY (`id`),
  UNIQUE INDEX `id_UNIQUE` (`id` ASC),
  CONSTRAINT `annotation_id`
    FOREIGN KEY (`id`)
    REFERENCES `$dbase`.`$Prefix"."_Genes` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = MyISAM;

SET SQL_MODE=\@OLD_SQL_MODE;
SET FOREIGN_KEY_CHECKS=\@OLD_FOREIGN_KEY_CHECKS;
SET UNIQUE_CHECKS=\@OLD_UNIQUE_CHECKS;
";
	return $monster;
}

sub _checkConfig { 
	my $Config = shift;
	my $Prefix = $Config->get("OPTIONS","SQL_Prefix");
	our $NetDir = $Config->get("OPTIONS","NetworkDir");
	die "Cannot find OPTIONS->NetworkDir : $NetDir!\n" unless -e $NetDir;
	my $inputExpression = $NetDir."/InputExpression.matrix.tab";
	die "Cannot find OPTIONS->inputExpression : $inputExpression!\n" unless -e $inputExpression;
	my $phytozomeAnnotation = $Config->get("OPTIONS","PhytozomeAnnotation");
	die "Cannot find OPTIONS->PhytozomeAnnotation : $phytozomeAnnotation!\n" unless -e $phytozomeAnnotation;
	my $anno_id_col  = $Config->get("OPTIONS","Annotation_id_column");
	die "Cannot find OPTIONS->Annotation_id_column : $anno_id_col!\n" if $anno_id_col eq "-1";
	my $anno_name_col= $Config->get("OPTIONS","Annotation_name_column");
	die "Cannot find OPTIONS->Annotation_name_column : $anno_name_col!\n" if $anno_name_col eq "-1"; 
	my $anno_desc_col= $Config->get("OPTIONS","Annotation_description_column");
	die "Cannot find OPTIONS->Annotation_description_column : $anno_desc_col!\n" if $anno_desc_col eq "-1";
	my $Adj		= $NetDir."/Adjacency.cytoscape.raw.tab";
	die "Cannot find ".$NetDir."/Adjacency.cytoscape.raw.tab\n" unless -e $NetDir."/Adjacency.cytoscape.raw.tab";
	my $ClusterDir	= $NetDir."/Clusters";
	die "Cannot find ". $NetDir."/Clusters\n" unless -e  $NetDir."/Clusters";
	my $Dbase	= $Config->get("OPTIONS","DBName");
	die "Cannot find option for Options->DBName!\n" if $Dbase eq "-1";
	my $Mask	= $Config->get("OPTIONS","MaskLevel");
	die "Cannot find option for Options->MaskLevel!\n" if $Mask == -1;
	die "Options->MaskLevel cannot be set to zero!\n" if $Mask == 0;
	die "Options->MaskLevel cannot be set to a negative value!\n" if $Mask<0;
	my $regex = $Config->get("OPTIONS","RegularExpression");
	die "Cannot find option for Options->RegularExpression!\n" if $regex eq "-1";
	my $displayName = $Config->get("OPTIONS","DisplayName");
	die "Cannot find option for Options->DisplayName!\n" if $displayName eq "-1";
	return 1;
}
