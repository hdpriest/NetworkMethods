#!/usr/bin/perl
use warnings;
use strict;
use FindBin;
use lib "$FindBin::Bin/../Lib";
package NetworkTools;

use Benchmark;
use threads;

###### No constructor. This isn't a network object, it is a method box. 
our $R="bogus";

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
	my $R=shift;
	my $Rvar=shift;
	my $file=shift;
	my $Rcmd="save($Rvar, file=\"$file\", ascii=TRUE);";
	warn  $Rcmd."\n";
	$R->send($Rcmd);
	return 1;
}

sub saveRobjToFile {
	my $self=shift;
	my $file=shift;
	my $Rvar=shift;
	my $Rcmd="save($Rvar, file=\"$file\");";
	$R->send($Rcmd);
	return 1;
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

sub getMatrixHeader {
	my $self=shift;
	my $R=shift;
	my $val=shift;
	my $cmd="colnames($val)";
	$R->send($cmd);
	return _parseVector($R->read());
}

sub getMatrixRow {
	my $self=shift;
	my $R=shift;
	my $val=shift;
	my $row=shift;
	my $cmd="$val\[$row,\]";
	$R->send($cmd);
	return _parseMatrixRow($R->read());
}





sub _parseVector { ##Victor
	my $vector=shift;
	my @vector;
	my @lines=split(/\n/,$vector);
	my $header=shift @lines if $lines[0]=~m/\$/;
	foreach my $line (@lines){
		$line=~s/\"//g;
		$line=~s/\s*\[\d+\]\s+//;
		my @line=split(/\s+/,$line);
	#	@line = map{$_=~s/\"//g} @line;
		map {push @vector, $_} @line;
	}
	return \@vector;
}

sub _parseMatrixRow {
	my $row=shift;
	my @data=split("\n",$row);
	my @headers;
	my @values;
	my $ID;
	for(my$i=0;$i<=$#data-1;$i+=2){
		my $head=$data[$i];
		my $data=$data[$i+1];
		$head=~s/^\s+//;
		@headers=(@headers,split(/\s+/,$head));
		my @row=split(/\s+/,$data);
		$ID=shift @row;
		@values=(@values,@row);

	}
	return ($ID,\@headers,\@values);
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

