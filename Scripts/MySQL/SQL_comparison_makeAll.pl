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
my $Prefix = $Config->get("OPTIONS","Comp_Prefix");
our $NetDir = $Config->get("OPTIONS","NetworkDir");
my $Configuration_One = $Config->get("OPTIONS","Config_One");
my $Configuration_Two = $Config->get("OPTIONS","Config_Two");
my $C1 = Configuration->new($Configuration_One);
my $C2 = Configuration->new($Configuration_Two);
my $Prefix_One = $C1->get("OPTIONS","SQL_Prefix");
my $Prefix_Two = $C2->get("OPTIONS","SQL_Prefix");
my $Exp_One = $C1->get("OPTIONS","NetworkDir")."/InputExpression.matrix.tab";
my $Exp_Two = $C2->get("OPTIONS","NetworkDir")."/InputExpression.matrix.tab";
my $Clusters_Pos	= $NetDir."/Clusters_Pos";
my $Clusters_Neg	= $NetDir."/Clusters_Neg";
my $Dbase	= $Config->get("OPTIONS","DBName");
my $regex = $Config->get("OPTIONS","RegularExpression");
my $displayName = $Config->get("OPTIONS","DisplayName");
our $Mask = $Config->get("OPTIONS","MaskLevel");

my $Adj	= $NetDir."/Adjacency.negPlasticity.cytoscape.tab";
my @Adjacency=@{Tools->LoadFile($Adj)};
$Adj		= $NetDir."/Adjacency.posPlasticity.cytoscape.tab";
push @Adjacency, @{Tools->LoadFile($Adj)};

my @Genes  = sort {$a cmp $b} @{Tools->LoadFile($C1->get("OPTIONS","MasterGeneList"))};
warn "Loading ".$C1->get("OPTIONS","MasterGeneList") ." as master gene list\n";
my %Genes;
for(my$i=0;$i<=$#Genes;$i++){
	my $I=$i+1;
	$Genes{$Genes[$i]}=$I;
#### for unique universal gene id
}


my @AltGenes  = sort {$a cmp $b} @{Tools->LoadFile($C2->get("OPTIONS","MasterGeneList"))};
warn "Loading ".$C2->get("OPTIONS","MasterGeneList") ." as master gene list\n";
my %AltGenes;
for(my$i=0;$i<=$#AltGenes;$i++){
	my $I=$i+1;
	$AltGenes{$AltGenes[$i]}=$I;
#### for unique universal gene id
}

warn "Making Create Statements...\n";
my $Create = _GetCreateStatement($Dbase,$Prefix,$Prefix_One,$Prefix_Two);
$Create .= _createMetricTable($Dbase,$Prefix);
print $Create."\n";
my $Map = _getCompareReferenceStatement($Dbase,$regex,$Prefix,$Prefix_One,$Prefix_Two,$displayName);
print $Map."\n";
my $Statement = _prepClusterTable($Dbase,$Prefix,$Clusters_Pos,"_PositiveClusters")."\n";
$Statement .= _prepClusterTable($Dbase,$Prefix,$Clusters_Neg,"_NegativeClusters");

print $Statement."\n";

warn "Creating Cluster Load Inserts...\n";
my $ClusterLoadCommand = _loadClusters($Dbase,$Prefix,"_NegativeClusters_Genes",$Clusters_Neg,\%Genes)."\n";
$ClusterLoadCommand .= _loadClusters($Dbase,$Prefix,"_PositiveClusters_Genes",$Clusters_Pos,\%Genes);
print $ClusterLoadCommand."\n";
my $ClusterAnnotateCommand = _loadClusterAnnotation($Dbase,$Prefix,$Config);
print $ClusterAnnotateCommand."\n";

warn "Creating Adjacency Load Inserts...\n";
my $AdjacencyLoadCommand = _loadAdjacency($Dbase,$Prefix,\%Genes,\%AltGenes,\@Adjacency);
print $AdjacencyLoadCommand."\n";
#print "DO SLEEP(5);\n";
warn "Creating Metrics Load Inserts...\n";
my $MetricLoadCommand = _prepNetworkMetrics(\%Genes,\%AltGenes,\@Adjacency,$Dbase,$Prefix,$Clusters_Pos,$Clusters_Neg,$Exp_One,$Exp_Two);
print $MetricLoadCommand."\n";
warn "Done.\n";

sub _loadClusterAnnotation {
	my $statement="";
	my $Dbase = shift;
	my $Prefix = shift;
	my $Config = shift;
	return $statement unless $Config->get("CLUSTER_ANNOTATION","Load_Annotations") == 1;
	my $p_dir 	= $Config->get("OPTIONS","NetworkDir")."/".$Config->get("CLUSTER_ANNOTATION","Positive_subdirectory");
	my $n_dir 	= $Config->get("OPTIONS","NetworkDir")."/".$Config->get("CLUSTER_ANNOTATION","Negative_subdirectory");
	my $del 	= $Config->get("CLUSTER_ANNOTATION","delimiter");
	my $term_col= $Config->get("CLUSTER_ANNOTATION","term_column");
	my $desc_col= $Config->get("CLUSTER_ANNOTATION","description_column");
	my $fdrp_col= $Config->get("CLUSTER_ANNOTATION","FDR_column");
	my $fdr_cut = $Config->get("CLUSTER_ANNOTATION","FDR_maximum_cutoff");
	my $file_reg= $Config->get("CLUSTER_ANNOTATION","File_Regex");
	my $TableName = "`$Dbase`.`$Prefix"."_NegativeClusters_Annotations`";
	my @Files = grep {m/$file_reg/} @{Tools->LoadDir($n_dir)};
	$statement=_getClusterState($TableName,$n_dir,\@Files,$del,$term_col,$desc_col,$fdrp_col,$fdr_cut);
	$TableName = "`$Dbase`.`$Prefix"."_PositiveClusters_Annotations`";
	@Files = grep {m/$file_reg/} @{Tools->LoadDir($p_dir)};
	$statement.=_getClusterState($TableName,$p_dir,\@Files,$del,$term_col,$desc_col,$fdrp_col,$fdr_cut);
	return $statement;
}

sub _getClusterState {
	my $tableName=shift;
	my $dir=shift;
	my @Files=@{$_[0]};
	my $del = $_[1];
	my $term_col = $_[2];
	my $desc_col = $_[3];
	my $fdrp_col = $_[4];
	my $fdr_cut=$_[5];
	my $statement="";
	foreach my $file (@Files){
		my $path = $dir."/".$file;
		$file=~m/\.(\d+)\./;
		my $cid = $1;
		my @content = @{Tools->LoadFile($path)};
		foreach my $line (@content){
			my @line=split($del,$line);
			my $term = $line[$term_col];
			my $desc = $line[$desc_col];
			my $fdrp = $line[$fdrp_col];
			next if ($fdrp > $fdr_cut);
			$fdrp = 0 if $fdrp eq "NA";
			$fdrp = sprintf("%.5f",$fdrp);
			$statement .= "INSERT INTO $tableName (`Cluster`,`term`,`description`,`FDR_p`) VALUES ($cid,'$term','$desc','$fdrp');\n";
		}
	}
	return $statement;
}

sub _loadAdjacency {
	my $DB = shift;
	my $Prefix = shift;
	my %genes=%{$_[0]};
	my %altgenes=%{$_[1]};
	my @Adj=@{$_[2]};
	my $edgeCSV = "$cwd/$Prefix.EdgeLoad.csv";
	my $adjCSV  = "$cwd/$Prefix.AdjLoad.csv";
	my $edge_num=1;
	my @EdgeCSV;
	my @AdjCSV;
	my $adjCols = "(id,edgeId)";
	my $edgeCols = "(edgeId,adjacency,node1,node2,alt_node1,alt_node2)";
	foreach my $line (@Adj){
		chomp $line;
		my ($id_A,$id_B,$edge,$absEdge,$alt_A,$alt_B)=split(/\t/,$line);
		$edge=sprintf("%.3f",$edge);
		next if abs($edge) < $Mask;
		die "$id_A undefined in master gene list\n" unless defined $genes{$id_A};
		die "$id_B undefined in master gene list\n" unless defined $genes{$id_B};
		my $ida = $genes{$id_A};
		my $idb = $genes{$id_B};
		my $alt_ida = $altgenes{$alt_A};
		my $alt_idb = $altgenes{$alt_B};
		push @AdjCSV, $ida.",".$edge_num;
		push @AdjCSV, $idb.",".$edge_num;
		push @EdgeCSV, "$edge_num,$edge,$ida,$idb,$alt_ida,$alt_idb";
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
	return $insert;
}

sub _PrepModules {
	my $ClusterDir = shift;
	my $Gref = shift;
	my $Mref = shift;
	my $field= shift;
	my %Genes= %$Gref;
	my %Modules = %$Mref;
	my %used; ## doubly ensures single-cluster membership
	warn "Working on $ClusterDir\n";
	my @Files = grep {m/txt$/} @{Tools->LoadDir($ClusterDir)};
	foreach my $file (@Files){
		my $path=$ClusterDir."/".$file;
#		warn "working on $path\n";
		$file=~m/Cluster\.(\d+)\.txt/;
		my $c=$1;
		$c = $c * -1 if $field=~m/neg/i;
#		warn "Cluster=$c\n";
		$Modules{$c}=[];
		my @file=@{Tools->LoadFile($path)};
		foreach my $id (@file){
			next if defined $used{$id};
			die "Cannot find $id in gene hash!\n" unless defined $Genes{$id};
			$Genes{$id}{$field}=$c;
			push @{$Modules{$c}}, $id;
			$used{$id}=1;
		}
	}
	my @refs = (\%Genes,\%Modules);
	return \@refs;
}

sub _prepNetworkMetrics {
#my $MetricLoadCommand = _prepNetworkMetrics(\%Genes,\%Adjacency,$Dbase,$Prefix,$Clusters_Pos,$Clusters_Neg);
	my $MasterList = shift;
	my $AltList = shift;
	my $Adjacency = shift;
	my $DB = shift;
	my $Prefix = shift;
	my $Clusters_Pos = shift;
	my $Clusters_Neg = shift;
	my $E1 = shift;
	my $E2 = shift;
	my $outFile ="$cwd/".$Prefix.".Metrics.csv";
	my @output;
	my %g = %{$MasterList};
	my @Genes = keys %g;
	warn "Prepping Genes\n";
	my %Genes = %{_PrepGenes($MasterList,$E1,$E2)};
	my %Modules;
	warn "Prepping Modules (1)\n";
	my ($r1,$r2) = @{_PrepModules($Clusters_Pos,\%Genes,\%Modules,"pos_module")};
	%Genes = %$r1;
	%Modules = %$r2;
	warn "Prepping Modules (2)\n";
	($r1,$r2) = @{_PrepModules($Clusters_Neg,\%Genes,\%Modules,"neg_module")};
	%Genes = %$r1;
	%Modules = %$r2;
	warn "Prepping Adjacency\n";
	%Genes = %{_PrepAdj($Adjacency,\%Genes)};
	warn "Calculating Metrics\n";
	%Genes = %{_Rank(\%Genes,\%Modules)};
	foreach my $gene (keys %Genes){
		my @line;
		die "$gene not defined in master list\n" unless defined $g{$gene};
		push @line, $g{$gene};
		push @line, sprintf("%.3f",$Genes{$gene}{neg_modular_k});
		push @line, $Genes{$gene}{neg_modular_k_rank};
		push @line, sprintf("%.3f",$Genes{$gene}{pos_modular_k});
		push @line, $Genes{$gene}{pos_modular_k_rank};
		push @line, sprintf("%.3f",$Genes{$gene}{pos_mod_mean_exp_rank_1});
		push @line, sprintf("%.3f",$Genes{$gene}{pos_mod_mean_exp_rank_2});
		push @line, sprintf("%.3f",$Genes{$gene}{neg_mod_mean_exp_rank_1});
		push @line, sprintf("%.3f",$Genes{$gene}{neg_mod_mean_exp_rank_2});
		push @line, $Genes{$gene}{pos_module};
		push @line, $Genes{$gene}{neg_module};
		push @line, sprintf("%.3f",$Genes{$gene}{k});
		push @line, $Genes{$gene}{k_rank};
		push @line, sprintf("%.3f",$Genes{$gene}{exp_mean_1});
		push @line, sprintf("%.3f",$Genes{$gene}{exp_mean_2});
		push @line, $Genes{$gene}{exp_rank_1};
		push @line, $Genes{$gene}{exp_rank_2};
		my $line = join(",",@line);
		push @output, $line;
	}
	Tools->printToFile($outFile,\@output);
	my $columns = "(id,";
	$columns .= "neg_modular_k,";
	$columns .= "neg_modular_k_rank,";
	$columns .= "pos_modular_k,";
	$columns .= "pos_modular_k_rank,";
	$columns .= "pos_modular_mean_exp_rank_1,";
	$columns .= "pos_modular_mean_exp_rank_2,";
	$columns .= "neg_modular_mean_exp_rank_1,";
	$columns .= "neg_modular_mean_exp_rank_2,";
	$columns .= "pos_module,";
	$columns .= "neg_module,";
	$columns .= "k,";
	$columns .= "k_rank,";
	$columns .= "exp_mean_1,";
	$columns .= "exp_mean_2,";
	$columns .= "exp_rank_1,";
	$columns .= "exp_rank_2)";
	my $insert = "LOAD DATA INFILE '$outFile' INTO TABLE `$Dbase`.`$Prefix"."_Metrics` fields terminated by ',' enclosed by '' lines terminated by '\\n' $columns;";
	return $insert;
}

sub _Rank {
	my %Genes=%{$_[0]};
	my %Modules=%{$_[1]};
	foreach my $module (sort {$a <=> $b} keys %Modules){
		my @genes = @{$Modules{$module}};
		my $G = scalar(@genes);
		my %this_module;
		for(my$i=0;$i<=$#genes;$i++){
			if($module<0){
				$this_module{$genes[$i]}=$Genes{$genes[$i]}{neg_modular_k};
			}else{
				$this_module{$genes[$i]}=$Genes{$genes[$i]}{pos_modular_k};
			}
		}
		foreach my $gene (sort {$this_module{$a} <=> $this_module{$b}} keys %this_module){ ## low to high sort
			my $rank = $G;
			if($module<0){
				$Genes{$gene}{neg_modular_k_rank}=$rank;
			}else{
				$Genes{$gene}{pos_modular_k_rank}=$rank;
			}
			$G--;
		}
		for(my$i=0;$i<=$#genes;$i++){
			$this_module{$genes[$i]}=$Genes{$genes[$i]}{exp_mean_1};
		}
		$G = scalar(@genes);
		foreach my $gene (sort {$this_module{$a} <=> $this_module{$b}} keys %this_module){ ## low to high sort
			my $rank = $G;
			if($module<0){
				$Genes{$gene}{neg_mod_mean_exp_rank_1}=$rank;
			}else{
				$Genes{$gene}{pos_mod_mean_exp_rank_1}=$rank;
			}
			$G--;
		}
		for(my$i=0;$i<=$#genes;$i++){
			$this_module{$genes[$i]}=$Genes{$genes[$i]}{exp_mean_2};
		}
		$G = scalar(@genes);
		foreach my $gene (sort {$this_module{$a} <=> $this_module{$b}} keys %this_module){ ## low to high sort
			my $rank = $G;
			if($module<0){
				$Genes{$gene}{neg_mod_mean_exp_rank_2}=$rank;
			}else{
				$Genes{$gene}{pos_mod_mean_exp_rank_2}=$rank;
			}
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
	foreach my $gene (sort {$Genes{$a}{exp_mean_1} <=> $Genes{$b}{exp_mean_1}} keys %Genes){
		my $rank = $N;
		$Genes{$gene}{exp_rank_1}=$rank;
		$N--;
	}
	$N = scalar(keys(%Genes));
	foreach my $gene (sort {$Genes{$a}{exp_mean_2} <=> $Genes{$b}{exp_mean_2}} keys %Genes){
		my $rank = $N;
		$Genes{$gene}{exp_rank_2}=$rank;
		$N--;
	}
	return \%Genes;
}

sub _PrepAdj {
	my $Adjacency = shift;
	my %Genes = %{$_[0]};
	foreach my $line (@$Adjacency){
		chomp $line;
		my ($IDA,$IDB,$edge,$absEdge,$alt_IDA,$alt_IDB)=split(/\s/,$line);
		next if abs($edge)<$Mask;
		die "Cannot find $IDA in Gene table!\n" unless defined $Genes{$IDA};
		die "Cannot find $IDB in Gene table!\n" unless defined $Genes{$IDB};
		if(($Genes{$IDA}{pos_module} == $Genes{$IDB}{pos_module}) && ($Genes{$IDA}{pos_module}!=0)){
			$Genes{$IDA}{pos_modular_k}+=$absEdge;
			$Genes{$IDB}{pos_modular_k}+=$absEdge;
		}else{
		}
		if(($Genes{$IDA}{neg_module} == $Genes{$IDB}{neg_module}) && ($Genes{$IDA}{neg_module}!=0)){
			$Genes{$IDA}{neg_modular_k}+=$absEdge;
			$Genes{$IDB}{neg_modular_k}+=$absEdge;
		}else{
		}
		$Genes{$IDA}{k}+=$absEdge;
		$Genes{$IDB}{k}+=$absEdge;
	}
	return \%Genes;
}

sub _loadClusters {
	my $Dbase=shift;
	my $prefix=shift;
	my $table=shift;
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
	my $statement = "LOAD DATA INFILE '$outFile' INTO TABLE `$Dbase`.`$prefix"."$table` fields terminated by ',' enclosed by '' lines terminated by '\\n';";
	return $statement;
}

sub _PrepGenes {
	my $file=shift;
	my $E1=shift;
	my $E2=shift;
	my @Genes = keys %{$file};
	my %Genes;
	foreach my $gene (@Genes){
		$Genes{$gene}={};
		$Genes{$gene}{exp_mean_1}=0;
		$Genes{$gene}{exp_mean_2}=0;
		$Genes{$gene}{k}=0;
		$Genes{$gene}{k_rank}=0;
		$Genes{$gene}{exp_rank_1}=0;
		$Genes{$gene}{exp_rank_2}=0;
		$Genes{$gene}{pos_module}=0;
		$Genes{$gene}{pos_modular_k}=0;
		$Genes{$gene}{pos_modular_k_rank}=0;
		$Genes{$gene}{pos_mod_mean_exp_rank_1}=0;
		$Genes{$gene}{pos_mod_mean_exp_rank_2}=0;
		$Genes{$gene}{neg_module}=0;
		$Genes{$gene}{neg_modular_k}=0;
		$Genes{$gene}{neg_modular_k_rank}=0;
		$Genes{$gene}{neg_mod_mean_exp_rank_1}=0;
		$Genes{$gene}{neg_mod_mean_exp_rank_2}=0;
	}
	my $DF1 = DataFrame->new();
	$DF1->loadFile($E1,"\t");
	$DF1->initRowIterator();
	$DF1->createIdByIndex();
	my $DF2 = DataFrame->new();
	$DF2->loadFile($E2,"\t");
	$DF2->initRowIterator();
	foreach my $gene (@Genes){
		my @Data = @{$DF1->getDataByID_withZeros($gene)};
		my $mean = Tools->mean(@Data);
		$Genes{$gene}{exp_mean_1}=$mean;
		my $index = $DF1->getIndexOfGene($gene);
		if(defined($index)){
			my $ID2=$DF2->getIDbyIndex($index);
			$ID2="Null" unless defined $ID2;
			@Data = @{$DF2->getDataByID_withZeros($ID2)};
			$mean = Tools->mean(@Data);
			$Genes{$gene}{exp_mean_2}=$mean;
		}else{
		}
	}
	return \%Genes;
}

sub _getCompareReferenceStatement{
	my $Dbase=shift;
	my $regex=shift;
	my $pfx=shift;
	my $pfx1=shift;
	my $pfx2=shift;
	my $dispname=shift;
	my $statement="INSERT INTO `$Dbase`.`MapCompare` (`prefix`,`pfx_control`,`pfx_compare`,`display_name`,`regex`) VALUES ('$pfx','$pfx1','$pfx2','$dispname','$regex');\n";
	return $statement;
}

sub _prepClusterTable {
	my $Dbase=shift;
	my $prefix=shift;
	my $ClusterDir = shift;
	my $table = shift;
	my @Files = grep {m/txt$/} @{Tools->LoadDir($ClusterDir)};
	my %used;
	my $statement="";
	foreach my $file (@Files){
		my $path=$ClusterDir."/".$file;
		$file=~m/Cluster\.(\d+)\.txt/;
		my $c=$1;
		$statement .= "INSERT INTO `$Dbase`.`$prefix"."$table` (`clusterID`) VALUES ($c);\n";
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
`neg_modular_k` DOUBLE NOT NULL,
`neg_modular_k_rank` INT(8) NOT NULL,
`pos_modular_k` DOUBLE NOT NULL,
`pos_modular_k_rank` INT(8) NOT NULL,
`pos_modular_mean_exp_rank_1` INT(8) NOT NULL,
`pos_modular_mean_exp_rank_2` INT(8) NOT NULL,
`neg_modular_mean_exp_rank_1` INT(8) NOT NULL,
`neg_modular_mean_exp_rank_2` INT(8) NOT NULL,
`pos_module` INT(8) NOT NULL,
`neg_module` INT(8) NOT NULL,
`k` DOUBLE NOT NULL,
`k_rank` INT(8) NOT NULL,
`exp_mean_1` DOUBLE NOT NULL,
`exp_mean_2` DOUBLE NOT NULL,
`exp_rank_1` INT(8) NOT NULL,
`exp_rank_2` INT(8) NOT NULL,
PRIMARY KEY (`id`),
UNIQUE INDEX `id_UNIQUE` (`id` ASC),
CONSTRAINT `gene_comp_metrics`
	FOREIGN KEY (`id`)
	REFERENCES `$DB`.`$PF"."_Genes` (`id`)
	ON DELETE NO ACTION
	ON UPDATE NO ACTION)
ENGINE = MyISAM;
";
	return $Create;
}

sub _GetCreateStatement {
	my $dbase = shift;
	my $prefix= shift;
	my $pfx1=shift;
	my $pfx2=shift;
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
-- Table `$dbase`.`MapCompare`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `$dbase`.`MapCompare` (
  `prefix` VARCHAR(45) NOT NULL,
  `pfx_control` VARCHAR(45) NOT NULL,
  `pfx_compare` VARCHAR(45) NOT NULL,
  `display_name` VARCHAR(45) NOT NULL,
  `regex` VARCHAR(1064) NOT NULL,
  INDEX `prefix_idx` (`prefix` ASC),
  INDEX `control_idx` (`pfx_control` ASC),
  INDEX `compare_idx` (`pfx_compare` ASC),
  INDEX `name_idx` (`display_name` ASC),
  INDEX `regex_idx` (`regex` ASC))
ENGINE = MyISAM;

-- -----------------------------------------------------
-- Table `$dbase`.`$Prefix"."_PositiveClusters`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `$dbase`.`$Prefix"."_PositiveClusters` (
  `clusterID` INT NOT NULL,
  PRIMARY KEY (`clusterID`),
  UNIQUE INDEX `clusterID_UNIQUE` (`clusterID` ASC))
ENGINE = MyISAM;

-- -----------------------------------------------------
-- Table `$dbase`.`$Prefix"."_NegativeClusters`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `$dbase`.`$Prefix"."_NegativeClusters` (
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
  `alt_node1` VARCHAR(45) NOT NULL,
  `alt_node2` VARCHAR(45) NOT NULL,
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
  CONSTRAINT `id_comp`
    FOREIGN KEY (`id`)
    REFERENCES `$dbase`.`$Prefix"."_Genes` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `edge_comp`
    FOREIGN KEY (`edgeId`)
    REFERENCES `$dbase`.`$Prefix"."_Edges` (`edgeId`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = MyISAM;


-- -----------------------------------------------------
-- Table `$dbase`.`$Prefix"."_PositiveClusters_Genes`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `$dbase`.`$Prefix"."_PositiveClusters_Genes` (
  `id` INT NOT NULL,
  `Cluster` INT NOT NULL,
  PRIMARY KEY (`id`),
  INDEX `ind_cluster` (`Cluster` ASC),
  CONSTRAINT `posgene`
    FOREIGN KEY (`id`)
    REFERENCES `$dbase`.`$Prefix"."_Genes` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `poscluster`
    FOREIGN KEY (`Cluster`)
    REFERENCES `$dbase`.`$Prefix"."_PositiveClusters` (`clusterID`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = MyISAM;

-- -----------------------------------------------------
-- Table `$dbase`.`$Prefix"."_PositiveClusters_Annotations`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `$dbase`.`$Prefix"."_PositiveClusters_Annotations` (
  `Cluster` INT NOT NULL,
  `term` VARCHAR(45) NOT NULL,
  `description` VARCHAR(4096) NOT NULL,
  `FDR_p` DOUBLE NOT NULL,
  INDEX `ind_cluster` (`Cluster` ASC),
  CONSTRAINT `cluster_ann`
    FOREIGN KEY (`Cluster`)
    REFERENCES `$dbase`.`$Prefix"."_PositiveClusters` (`clusterID`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = MyISAM;

-- -----------------------------------------------------
-- Table `$dbase`.`$Prefix"."_NegativeClusters_Genes`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `$dbase`.`$Prefix"."_NegativeClusters_Genes` (
  `id` INT NOT NULL,
  `Cluster` INT NOT NULL,
  PRIMARY KEY (`id`),
  INDEX `ind_cluster` (`Cluster` ASC),
  CONSTRAINT `neggene`
    FOREIGN KEY (`id`)
    REFERENCES `$dbase`.`$Prefix"."_Genes` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `negcluster`
    FOREIGN KEY (`Cluster`)
    REFERENCES `$dbase`.`$Prefix"."_NegativeClusters` (`clusterID`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = MyISAM;

-- -----------------------------------------------------
-- Table `$dbase`.`$Prefix"."_NegativeClusters_Annotations`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `$dbase`.`$Prefix"."_NegativeClusters_Annotations` (
  `Cluster` INT NOT NULL,
  `term` VARCHAR(45) NOT NULL,
  `description` VARCHAR(4096) NOT NULL,
  `FDR_p` DOUBLE NOT NULL,
  INDEX `ind_cluster` (`Cluster` ASC),
  CONSTRAINT `cluster_ann`
    FOREIGN KEY (`Cluster`)
    REFERENCES `$dbase`.`$Prefix"."_NegativeClusters` (`clusterID`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = MyISAM;
";
	return $monster;
}

sub _checkConfig { 
	my $Config = shift;
	my $Prefix = $Config->get("OPTIONS","Comp_Prefix");
	die "Cannot find option for OPTIONS->Comp_Prefix!\n" if $Prefix eq "-1";
	my $Configuration_One = $Config->get("OPTIONS","Config_One");
	my $Configuration_Two = $Config->get("OPTIONS","Config_Two");
	die "Cannot find OPTIONS->Config_One\n" unless -e $Configuration_One;
	die "Cannot find OPTIONS->Config_Two\n" unless -e $Configuration_Two;
	my $C1 = Configuration->new($Configuration_One);
	my $C2 = Configuration->new($Configuration_Two);
	my $Exp_One = $C1->get("OPTIONS","NetworkDir")."/InputExpression.matrix.tab";
	my $Exp_Two = $C2->get("OPTIONS","NetworkDir")."/InputExpression.matrix.tab";
	die "Cannot find input expression matrix from Config_One: $Exp_One\n" unless -e $Exp_One;
	die "Cannot find input expression matrix from Config_Two: $Exp_Two\n" unless -e $Exp_Two;
	
	my $NetDir = $Config->get("OPTIONS","NetworkDir");
	die "Cannot find OPTIONS->NetworkDir : $NetDir!\n" unless -e $NetDir;

	my $Clusters_Pos	= $NetDir."/Clusters_Pos";
	my $Clusters_Neg	= $NetDir."/Clusters_Neg";
	my $Dbase	= $Config->get("OPTIONS","DBName");
	my $regex = $Config->get("OPTIONS","RegularExpression");
	my $displayName = $Config->get("OPTIONS","DisplayName");
	my $Adj	= $NetDir."/Adjacency.negPlasticity.cytoscape.tab";
	die "Cannot find negative differential adjacency file: $Adj\n" unless -e $Adj;
	$Adj		= $NetDir."/Adjacency.posPlasticity.cytoscape.tab";
	die "Cannot find positive differential adjacency file: $Adj\n" unless -e $Adj;
	die "Cannot find positive cluster directory: $Clusters_Pos\n" unless -e $Clusters_Pos;
	die "Cannot find negative cluster directory: $Clusters_Neg\n" unless -e $Clusters_Neg;
	die "Cannot find option for Options->DBName!\n" if $Dbase eq "-1";
	die "Cannot find option for Options->RegularExpression!\n" if $regex eq "-1";
	die "Cannot find option for Options->DisplayName!\n" if $displayName eq "-1";
	return 1;
}
