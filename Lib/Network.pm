#!/usr/bin/perl
use warnings;
use strict;
use FindBin;
use lib "$FindBin::Bin/../Library";
package WGCNA;
use Benchmark;
use threads;
use Statistics::R;


my $R=Statistics::R->new();
$R->stop();
$R->start();
$R->send("library(dynamicTreeCut, lib.loc=\"/home/hpriest/R/x86_64-redhat-linux-gnu-library/\")");
$R->send("library(WGCNA, lib.loc=\"/home/hpriest/R/x86_64-redhat-linux-gnu-library/\")");


sub new {
        my $class=shift;
        my $self = {
	  	nObs	=> {},
		nGenes	=> {},
		nThreads	=> shift,
        };
        bless $self, $class;
	  $self->{nThreads}=2 unless defined $self->{nThreads};
	  my $cmd="enableWGCNAThreads(nThreads = ".$self->{nThreads}.")";
	  $R->send($cmd);
        return $self;
}

sub addExpressionData {
	my $self=shift;
	my @info=@_;
	my $Eid=shift @info;
	my $id=shift @info;
	$self->{IDs}={} unless defined ($self->{IDs});
	$self->{IDs}{$Eid}=[] unless defined ($self->{IDs}{$Eid});
	$self->{Expr}={} unless defined($self->{Expr});
	$self->{Expr}{$Eid}={} unless defined($self->{Expr}{$Eid}); 
	if(defined($self->{Expr}{$Eid}{$id})){
		map {push @{$self->{Expr}{$Eid}{$id}}, $_} @info;
	}else{
		$self->{Expr}{$Eid}{$id}=\@info;
		push @{$self->{IDs}{$Eid}}, $id;
	}
	$self->{nObs}{$Eid}=scalar(@{$self->{Expr}{$Eid}{$id}});
	$self->{nGenes}{$Eid}=scalar(keys %{$self->{Expr}{$Eid}});
	return 1;
}


############# FILTERING - min/max variance, min/max presence, min/max expression, etci

sub _initMatrix {
	my $n=shift;
	my $G=shift;
	my @array;
	for(my$i=0;$i<=$G-1;$i++){
		$array[$i]=[];
		for(my$j=0;$j<=$n-1;$j++){
			$array[$i][$j]=0;
		}
	}
	return \@array;
}

sub getR {
	return $R;
}

sub startR {
	$R->startR();
	return 1;
}

sub stopR {
	$R->stopR();
	return 1;
}

sub stopRlog {
	my $self=shift;
	close RLOG;
	return 1;
}

sub _readFromR {
	return $R->read();
}

sub _sendToR {
	my $line=shift;
	$R->send($line);
	return 1;
}


sub generateDendroColorPlot {
	my $self=shift;
	my $Tree=shift;
	my $file=shift;
	my $N=shift;
	my $Rcmd="png(file=\"$file\",width=900,height=900,units=\"px\",pointsize=12,bg=\"white\",res=NA)";
	$R->send($Rcmd);
	print RLOG $Rcmd."\n";
	$Rcmd="plotDendroAndColors($Tree,rainbow($N))";
	$R->send($Rcmd);
	print RLOG $Rcmd."\n";
	$Rcmd="dev.off()";
	$R->send($Rcmd);
	print RLOG $Rcmd."\n";
	return 1;
}

sub getTreeMerge {
	my $self=shift;
	my $Tree=shift;
	my $Rcmd="$Tree\[1\]";
	print RLOG $Rcmd."\n";
	$R->send($Rcmd);
	return $R->read();
}

sub getTreeHeight {
	my $self=shift;
	my $Tree=shift;
	my $Rcmd="$Tree\[2\]";
	print RLOG $Rcmd."\n";
	$R->send($Rcmd);
	return $R->read();
}

sub _getTreeOrder {
	my $Tree=shift;
	my $Rcmd="$Tree\[3\]";
	print RLOG $Rcmd."\n";
	$R->send($Rcmd);
	return $R->read();
}

sub getTreeOrder {
	my $self=shift;
	my $Tree=shift;
	return _getTreeOrder($Tree);
}

sub exportCytoscape {
	my $self=shift;
	my $Adj=shift;
	my $Nfile=shift;
	my $Efile=shift;
	my $NameVar=shift;
	my $Rcmd="exportNetworkToCytoscape($Adj, edgeFile = \"$Efile\", nodeFile = \"$Nfile\", nodeNames=$NameVar)";
	$R->send($Rcmd);
	return 1;
}

sub getDynamicMods {
	my $self=shift;
	my $Tree=shift;
	my $DissMat=shift;
	my $Rcmd="minModuleSize = 5;
	dynamicMods = cutreeDynamic(dendro = $Tree, distM = $DissMat, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);
	table(dynamicMods)
	";
	$R->send($Rcmd);
	print RLOG $Rcmd."\n";
	my $return=$R->read();
	return $return;
}

sub generateTree {
	my $self=shift;
	my $TOM=shift;
	my $DissMat=shift;
	my $Tree=shift;
	my $Rcmd="$DissMat <- $TOM";
	$R->send($Rcmd);
	print RLOG $Rcmd."\n";
	warn "Calculating Tree Cluster Solution....\n";
#	$Rcmd="$Tree <- flashClust(as.dist($DissMat),method=\"centroid\")";
	$Rcmd="$Tree <- flashClust(as.dist($DissMat),method=\"average\")";
	my $t0=Benchmark->new();
	$R->send($Rcmd);
	my $t1=Benchmark->new();
	my $td= timediff($t1,$t0);
	warn "Elapsed time: ".timestr($td)."\n";
	print RLOG $Rcmd."\n";
	return $Tree;
}

sub _cutreeStaticColor {
	my $Tree=shift;
	my $height=shift;
	my $minSize=shift;
	my $colorVect=shift;
	my $Rcmd="$colorVect <- cutreeStaticColor($Tree, cutHeight = $height, minSize = $minSize)";
	$R->send($Rcmd);
	print RLOG $Rcmd;
	return $colorVect;
}

sub generateTOMplot {
	my $self=shift;
	my $TOMvar=shift;
	my $Tree=shift;
	my $file=shift;
	my $cut=shift;
	my $CutHeight=shift;
	my $minSize=shift;
	die "Can't call ". (caller(0))[3]." unless a Tree (hclust/flashclust) variable is passed!\n\n" unless $Tree;
	die "Can't call ". (caller(0))[3]." unless a TOM (TOMdist) variable is passed!\n\n" unless $TOMvar;
	die "Can't call ". (caller(0))[3]." unless a File (.png) variable is passed!\n\n" unless $file;
	my $Rcmd="TempDissMat <- $TOMvar";
	$R->send($Rcmd);
	print RLOG $Rcmd."\n";
	$Rcmd="DissPlot <- TempDissMat^7";
	$R->send($Rcmd);
	print RLOG $Rcmd."\n";
	if($cut==1){
		$Rcmd="ColorVect <- cutreeStaticColor($Tree, cutHeight = $CutHeight, minSize = $minSize)";
		$R->send($Rcmd);
		print RLOG $Rcmd."\n";
	}
	$Rcmd="png(file=\"$file\",width=900,height=900,units=\"px\",pointsize=12,bg=\"white\",res=NA)";
	$R->send($Rcmd);
	print RLOG $Rcmd."\n";
	warn "Generating TOM plot...\n";
	my $t0=Benchmark->new();
	$Rcmd="TOMplot(DissPlot,$Tree,colors=ColorVect,main=\"Network heatmap plot, all genes\")";
	#$Rcmd="TOMplot($TOMvar,$Tree,colors=ColorVect,main=\"Network heatmap plot, all genes\")";
	$R->send($Rcmd);
	warn "Done.\n";
	print RLOG $Rcmd."\n";
	$Rcmd="dev.off()";
	$R->send($Rcmd);
	my $t1=Benchmark->new();
	my $td= timediff($t1,$t0);
	warn "Elapsed time: ".timestr($td)."\n";
	print RLOG $Rcmd."\n";
	return 1;
}

sub generateTOMplotDynamic {
	my $self=shift;
	my $TOMvar=shift;
	my $Tree=shift;
	my $file=shift;
	my $cut=shift;
	my $CutHeight=shift;
	my $minSize=shift;
	die "Can't call ". (caller(0))[3]." unless a Tree (hclust/flashclust) variable is passed!\n\n" unless $Tree;
	die "Can't call ". (caller(0))[3]." unless a TOM (TOMdist) variable is passed!\n\n" unless $TOMvar;
	die "Can't call ". (caller(0))[3]." unless a File (.png) variable is passed!\n\n" unless $file;
	my $Rcmd="minModuleSize = $minSize;
	dynamicMods = cutreeDynamic(dendro = $Tree, distM = $TOMvar, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);
	table(dynamicMods);
	dyanmicColors = labels2colors(dynamicMods);
	";
	$R->send($Rcmd);
	$Rcmd="png(file=\"$file\",width=900,height=900,units=\"px\",pointsize=12,bg=\"white\",res=NA)";
	$R->send($Rcmd);
	print RLOG $Rcmd."\n";
	warn "Generating TOM plot...\n";
	my $t0=Benchmark->new();
	$Rcmd="TOMplot($TOMvar^4,$Tree,colors=dynamicColors,main=\"Network heatmap plot, all genes\")";
#	main = "TOM heatmap plot, module genes" )
	$R->send($Rcmd);
	warn "Done.\n";
	print RLOG $Rcmd."\n";
	$Rcmd="dev.off()";
	$R->send($Rcmd);
	my $t1=Benchmark->new();
	my $td= timediff($t1,$t0);
	warn "Elapsed time: ".timestr($td)."\n";
	print RLOG $Rcmd."\n";
	return 1;
}


sub getAdjFromPower {
	my $self=shift;
	my $Sim=shift;
	my $B=shift;
	my $adjVar=shift;
	warn "Calculating Adjacency Matrix....\n";
	my $t0=Benchmark->new();
	my $cmd="$adjVar <- adjacency.fromSimilarity($Sim,type=\"unsigned\",power=$B)";
	print RLOG $cmd."\n";
	$R->send($cmd);
	my $t1=Benchmark->new();
	my $td= timediff($t1,$t0);
	warn "Elapsed time: ".timestr($td)."\n";
	return $adjVar;
}

sub getAdjFromSigMoid {
	my $self=shift;
	my $cvar=shift;
	my $T=shift;
	my $adjVar=shift;
	warn "Calculating Adjacency Matrix....\n";
	my $t0=Benchmark->new();
	my $mu=.8;
	my $alpha=20;
	my $cmd="$adjVar <- sigmoidAdjacencyFunction($cvar,mu=$mu,alpha=$alpha)";
	print RLOG $cmd."\n";
	$R->send($cmd);
	my $t1=Benchmark->new();
	my $td= timediff($t1,$t0);
	warn "Elapsed time: ".timestr($td)."\n";
	return $adjVar;
}

sub getAdjFromSigNum {
	my $self=shift;
	my $cvar=shift;
	my $T=shift;
	warn "Calculating Adjacency Matrix....\n";
	my $t0=Benchmark->new();
	my $adjVar="AdjMat";
	my $cmd="$adjVar <- signumAdjacencyFunction($cvar,$T)";
	print RLOG $cmd."\n";
	$R->send($cmd);
	my $t1=Benchmark->new();
	my $td= timediff($t1,$t0);
	warn "Elapsed time: ".timestr($td)."\n";
	return $adjVar;
}

sub getRobjFromFile {
	my $self=shift;
	my $file=shift;
	my $Rcmd="load(\"$file\");";
	print RLOG $Rcmd."\n";
	warn "Reading from $file...\n";
	my $t0=Benchmark->new();
	$R->send($Rcmd);
	my $t1=Benchmark->new();
	my $td= timediff($t1,$t0);
	warn "Elapsed time: ".timestr($td)."\n";
	return 1;
}

sub saveRobjToAscii {
	my $self=shift;
	my $file=shift;
	my $Rvar=shift;
#	my $Rcmd="FILE <- file(description=\"$file\");";
#	$R->send($Rcmd);
	my $Rcmd="save($Rvar, file=\"$file\", ascii=TRUE);";
	warn "Writing $Rvar to ASCII $file...\n";
	my $t0=Benchmark->new();
	$R->send($Rcmd);
	my $t1=Benchmark->new();
	my $td= timediff($t1,$t0);
	warn "Elapsed time: ".timestr($td)."\n";
	return 1;
}

sub saveRobjToFile {
	my $self=shift;
	my $file=shift;
	my $Rvar=shift;
#	my $Rcmd="FILE <- file(description=\"$file\");";
#	$R->send($Rcmd);
	my $Rcmd="save($Rvar, file=\"$file\");";
	warn "Writing $Rvar to $file...\n";
	my $t0=Benchmark->new();
	$R->send($Rcmd);
	my $t1=Benchmark->new();
	my $td= timediff($t1,$t0);
	warn "Elapsed time: ".timestr($td)."\n";
	return 1;
}

sub genTOM {
	my $self=shift;
	my $adjVar=shift;
	my $TOMvar=shift;
#	my $cmd="$TOMvar <- TOMsimilarity($adjVar, TOMType=\"unsigned\",TOMDenom=\"min\",verbose=0)";
	my $cmd="$TOMvar <- TOMdist($adjVar, TOMType=\"unsigned\",TOMDenom=\"min\",verbose=0,indent=0)";
	warn "Calculating Topological Overlap Matrix....\n";
	print RLOG $cmd."\n";
	my $t0=Benchmark->new();
	$R->send($cmd);
	my $t1=Benchmark->new();
	my $td= timediff($t1,$t0);
	warn "Elapsed time: ".timestr($td)."\n";
	return $TOMvar;
}

sub _getRvalueOf {
	my $val=shift;
	my $cmd="$val";
	$R->send($cmd);
	return $R->read();
}

sub getRvalueOf {
	my $self=shift;
	my $val=shift;
	return _getRvalueOf($val);
}

sub getMatrixEntryOf {
	my $self=shift;
	my $val=shift;
	my $i=shift;
	my $j=shift;
	my $cmd="$val\[$i,$j\]";
	$R->send($cmd);
	return $R->read();
}

sub getMatrixRow {
	my $self=shift;
	my $val=shift;
	my $row=shift;
	my $cmd="$val\[$row,\]";
	$R->send($cmd);
	return $R->read();
}

sub startRlog {
	my $self=shift;
	my $logFile=shift;
	open(RLOG,">",$logFile) || die "cannot open $logFile!\n$!\nexiting...\n";
	return 1;
}

####### The cor() function (for example) in R expects to calculate correlation between columns of a matrix, not rows, thus, we need to make a correctly formatted matrix:
sub populateExprMatrix {
	my $self=shift;
	my $Eid=shift;
	my $G=$self->{nGenes}{$Eid};
	my $n=$self->{nObs}{$Eid};
	my @ids=@{$self->{IDs}{$Eid}}; #### this is important - it maps the IDs in the expression hash to the columns of the expression matrix!
	my ($t0,$t1,$td);
	warn "Populating expression matrix...\n";
	$t0=Benchmark->new();
	$self->{ExprMatrix}={} unless defined($self->{ExprMatrix});
	$self->{ExprMatrix}{$Eid}=_initMatrix($n,$G); ### function returns an array REF of the appropriate size
	for(my$g=0;$g<=$G-1;$g++){
		for(my$j=0;$j<=$n-1;$j++){
			$self->{ExprMatrix}{$Eid}[$g][$j]=${$self->{Expr}{$Eid}{$ids[$g]}}[$j];
		}
	}
	$t1=Benchmark->new();
	$td= timediff($t1,$t0);
	warn "Elapsed time: ".timestr($td)."\n";
	return 1;
}

sub _genMatrixToRline {
	my $self=shift;
	my $type=shift;
	my $nrow=shift;
	my $ncol=shift;
	my $Eid=shift;
	my $matrixString="$type <- matrix(c(";
	my $matrixFoot="), nrow=$nrow,ncol=$ncol,byrow=TRUE)";
	my $matrixBody;
	for(my$j=0;$j<=$nrow-1;$j++){
		for(my$i=0;$i<=$ncol-1;$i++){
			if(defined($matrixBody)){
				$matrixBody.=",".$self->{ExprMatrix}{$Eid}[$i][$j];
			}else{
				$matrixBody.=$self->{ExprMatrix}{$Eid}[$i][$j];
			}
		}
	}
	return $matrixString.$matrixBody.$matrixFoot;
}

sub pickSoftValues {
	my $self=shift;
	my $type=shift;
	my $ret=shift;
	my $threads=$self->{nThreads};
	my $cmd="powers = c(c(1:10), seq(from = 12, to=20, by=2))";
	$R->send($cmd);
	print RLOG $cmd."\n";
	$cmd="$ret = pickSoftThreshold($type, powerVector = powers, verbose = 5)";
	print RLOG $cmd."\n";
	warn "picking soft threasholding values (power).....\n";
	my $t0=Benchmark->new();
	$R->send($cmd);
	my $t1=Benchmark->new();
	my $td=timediff($t1,$t0);
	warn "Elapsed time: ".timestr($td)."\n";
	return $ret;
}

sub graphSoftThresh {
	my $self=shift;
	my $sftVar=shift;
	my $ID=shift;
	my $outF=$ID.".Sft.pwr.png";
	my $cmd="png(file=\"$outF\");";
	$R->send($cmd);
	my $plotCmd="sizeGrWindow(9, 5)
	par(mfrow = c(1,2));
	cex1 = 0.9;
	plot($sftVar\$fitIndices[,1], -sign($sftVar\$fitIndices[,3])*$sftVar\$fitIndices[,2],
	xlab=\"Soft Threshold (power)\",ylab=\"Scale Free Topology Model Fit,signed R^2\",type=\"n\",
	main = paste(\"Scale independence\"));
	text($sftVar\$fitIndices[,1], -sign($sftVar\$fitIndices[,3])*$sftVar\$fitIndices[,2],
	labels=powers,cex=cex1,col=\"red\");
	abline(h=0.90,col=\"red\")
	plot($sftVar\$fitIndices[,1], $sftVar\$fitIndices[,5],
	xlab=\"Soft Threshold (power)\",ylab=\"Mean Connectivity\", type=\"n\",
	main = paste(\"Mean connectivity\"))
	text($sftVar\$fitIndices[,1], $sftVar\$fitIndices[,5], labels=powers, cex=cex1,col=\"red\")";
	$R->send($plotCmd);
	$R->send("dev.off();");
	return 1;
}


sub runBiCor {
	my $self=shift;
	my $type=shift;
	my $ret=shift;
	my $threads=$self->{nThreads};
	my $cmd="$ret <- bicor($type,nThreads=$threads)";
	print RLOG $cmd."\n";
	warn "calculating bicorrelation matrix.....\n";
	my $t0=Benchmark->new();
	$R->send($cmd);
	my $t1=Benchmark->new();
	my $td=timediff($t1,$t0);
	warn "Elapsed time: ".timestr($td)."\n";
	return $ret;
}


sub runCorrelation {
	my $self=shift;
	my $type=shift;
	my $ret=shift;
	my $cmd="$ret <- cor($type,y=NULL,method=\"pearson\",use=\"all.obs\")";
	print RLOG $cmd."\n";
	warn "running correlations...\n";
	my $t0=Benchmark->new();
	$R->send($cmd);
	my $t1=Benchmark->new();
	my $td=timediff($t1,$t0);
	warn "Elapsed time: ".timestr($td)."\n";
	return $ret;
}

sub sendExprMatrixToR {
	my $self=shift;
	my $Eid=shift;
	my $ExpVar=shift;
	my $G=$self->{nGenes}{$Eid};
	my $N=$self->{nObs}{$Eid};
	my $line=_genMatrixToRline($self,$ExpVar,$N,$G,$Eid);
	print RLOG "generated a matrix from matrix type $ExpVar, sending to R...\n";
	my ($t0,$t1,$td);
	warn "Sending expression matrix...\n";
	$t0=Benchmark->new();
	$R->send($line);
	$t1=Benchmark->new();
	$td= timediff($t1,$t0);
	warn "Elapsed time: ".timestr($td)."\n";
	return $ExpVar;
}

sub getModulesFromTree {
	my $self=shift;
	my $Tree=shift;
	my $height=shift;
	my $minSize=shift;
	my $colorVect="ColorVect";
	_cutreeStaticColor($Tree,$height,$minSize,$colorVect);
	my @order=@{_parseVector(_getTreeOrder($Tree))};
	my @colors=@{_parseVector(_getRvalueOf($colorVect))};
	return [\@order,\@colors];
}

sub _parseVector { ##Victor
	my $vector=shift;
	my @vector;
	my @lines=split(/\n/,$vector);
	my $header=shift @lines if $lines[0]=~m/\$/;
	foreach my $line (@lines){
		$line=~s/\s*\[\d+\]\s+//;
		my @line=split(/\s+/,$line);
	#	@line = map{$_=~s/\"//g} @line;
		map {push @vector, $_} @line;
	}
	return \@vector;
}

sub deconvoluteTree {
	my $self=shift;
	my $Tree=shift;
	my $Height=shift;
	my $hcut=shift;
	my @matrix=@{_parseMergeMatrix($Tree)};
	my @height=@{_parseHeightMatrix($Height)};
	my %Tree=%{_convertMergeToTree(\@matrix)};
	my %used;
	my $modN=1;
	my %Modules;
	for(my$i=$#height;$i>=0;$i--){
		next unless $height[$i]<$hcut;
		my $ind=$i+1;
		next if $used{$ind};
		my @list=@{$Tree{$ind}};
		$Modules{$modN}=[];
		foreach my $ele (@list){
			push @{$Modules{$modN}}, $ele unless $used{$ele};
			$used{$ele}=1;
		}
		$modN++;
	}
	return \%Modules;
}

sub _convertMergeToTree {
	my @matrix=@{$_[0]};
	my %tree;
	foreach my $element (@matrix){
		my $ind=shift @$element;
		if(defined($tree{$ind})){
			die "tree found $ind twice - that is so wrong.\n";
		}else{
			foreach my $e (@$element){
				if($e<0){
					my $E=$e*-1;
					push @{$tree{$ind}}, $E;
				}else{
					map {push @{$tree{$ind}}, $_} @{$tree{$e}};
				}
			}
		}
	}
	return \%tree;
}

sub _parseHeightMatrix {
	my $Height=shift;
	my @Heights=split(/\n/,$Height);
	my @ret;
	my $header=shift @Heights;
	foreach my $line (@Heights){
		$line=~s/\s*\[\d+\]\s+//;
		my @line=split(/\s+/,$line);
		map {push @ret, $_} @line
	}
	return \@ret;
}

sub _parseMergeMatrix {
	my $Matrix=shift;
	my @Matrix=split(/\n/,$Matrix);
	my @ret;
	my $header=shift @Matrix;
	my $indexHead=shift @Matrix;
	for(my$i=0;$i<=$#Matrix;$i++){
		my $line=$Matrix[$i];
		if($line=~m/\[(\d+)\,\]\s+(-?\d+)\s+(-?\d+)/){
			my @this=($1,$2,$3);
			push @ret, \@this;
		}else{
			die "Could not parse line: $line\n";
		}
	}
	return \@ret;
}


sub getGeneExprProfile {
	my $self=shift;
	my $id=shift;
	return $self->{Expr}{$id};
}

sub sendArrayToVector {

######### The curious structure of this command circumnavigates some weird behavior in the
#R bridge. Suffice to say, trying to join the array on "\",\"" and send it to R caused
#truncation in the R statement, and failure. This creates the intended data structure in
#small increments, and seems to work fine.
	my $self=shift;
	my $array=shift;
	my $name=shift;
	my @array=@$array;
	for(my$i=0;$i<=$#array;$i++){
		my $pname=$array[$i];
		if($i == 0){
			my $statement="$name <- c(\"$pname\");";
			$R->send($statement);
		}else{
			my $statement="$name <- append($name, \"$pname\");";
			$R->send($statement);
		}
	}
	sleep(5);
	return 1;
}

1;

