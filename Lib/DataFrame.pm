#!/usr/bin/perl
use warnings;
use strict;
use threads;
use threads::shared;
use FindBin;
use lib "$FindBin::Bin/../Lib";
use Configuration;
use Tools;
package DataFrame;




sub new {
	my $class=shift;
	my $self = {
		header	=> [],
		data		=> {},
		rowIDs	=> [], #### Row identifiers, in order. NEVER ALTER.
		IDindex 	=> {},
		currentRow	=> undef,
      };
      bless $self, $class;
      return $self;
}

sub lookaheadRow {
	my $self=shift;
	return 0 unless defined $self->{rowIDs}[$self->{currentRow}+1];
	return 1;
}

sub initRowIterator {
	my $self=shift;
	$self->{currentRow}=0;
	return 1;
}

sub createIdByIndex {
	my $self=shift;
	my %h;
	for(my $i=0;$i<=$#{$self->{rowIDs}};$i++){
		my $id = ${$self->{rowIDs}}[$i];
		$h{$id}=$i;
	}
	$self->{IDindex}=\%h;
	return 1;
}

sub getIndexOfGene {
	my $self=shift;
	my $gene=shift;
	if(defined($self->{IDindex}{$gene})){
		return $self->{IDindex}{$gene};
	}else{
		return undef;
	}
	return undef;
}

sub getIDbyIndex {
	my $self=shift;
	my $index=shift;
	if($index>$#{$self->{rowIDs}}){
		return undef;
	}else{
		my $ID=$self->{rowIDs}[$index];
		return $ID;
	}
	return undef;
}

sub getAllIDs {
	my $self=shift;
	return $self->{rowIDs};
}

sub getThisID {
	my $self=shift;
	die "cannot call getThisID without first loading a file!\n"  unless defined $self->{currentRow};
	if($self->{currentRow}>$#{$self->{rowIDs}}){
		return undef;
	}else{
		my $nextID=$self->{rowIDs}[$self->{currentRow}];
		return $nextID;
	}
	return undef;
}

sub getThisRow {
	my $self=shift;
	die "cannot call getNextRow without first loading a file!\n" unless defined $self->{currentRow};
	if($self->{currentRow}>$#{$self->{rowIDs}}){
		return undef;
	}else{
		my $ID=$self->{rowIDs}[$self->{currentRow}];
		my @row=@{$self->{data}{$ID}};
		return \@row;
	}
	return undef;
}

sub iterateRow {
	my $self=shift;
	$self->{currentRow}+=1;
	return $self->{currentRow};
}

sub getHistogramOfValuesByColumn {
	my $self= shift;
	my $col = shift;
	my %bins;
	foreach my $row (@{$self->{rowIDs}}){
		my @row=@{$self->{data}{$row}};
		my $value=sprintf("%.2f",$row[$col]);
		$bins{$value}+=1;
	}
	return \%bins;
}

sub getHistogramOfValues {
	my $self=shift;
	my %bins;
	foreach my $row (@{$self->{rowIDs}}){
		my @row=@{$self->{data}{$row}};
		foreach my $entry (@row){
			my $value=sprintf("%.2f",$entry);
			$bins{$value}+=1;
		}
	}
	return \%bins;
}

sub getIDsByMinValueRatio {
	my $self=shift;
	my $minValue=shift;
	my $ratio=shift;
	my @GRs;
	foreach my $row (@{$self->{rowIDs}}){
		die "$row undefined.\n" unless defined ($self->{data}{$row});
		if(_checkRowValueRatio($minValue,$ratio,$self->{data}{$row})){
			push @GRs, $row;
		}
	}
	return \@GRs;
}

sub getDataByColumn {
	my $self=shift;
	my @col_values;
	my $column;
	foreach my $row (@{$self->{rowIDs}}){
		my @row=@{$self->{data}{$row}};
		push @col_values, $row[$column];
	}
	return \@col_values;
}

sub getColumnTotals {
	my $self=shift;
	my @cols;
	foreach my $row (@{$self->{rowIDs}}){
		my @row=@{$self->{data}{$row}};
		for(my$i=0;$i<=$#row;$i++){
			$cols[$i]+=$row[$i];
		}
	}
	return \@cols;
}

sub printPerturbedFrame {
	my $self=shift;
	my $out =shift;
	my $del =shift;
	my @names=@{$self->{rowIDs}};
	my @output;
	push @output, "rowID".$del.join($del,@{$self->{header}});
	foreach my $row (@names){
		die "Cannot find $row in data hash. Ruh roh.\n\n" unless defined $self->{perturbation}{$row};
		my $line=$row.$del.join($del,@{$self->{perturbation}{$row}});
		push @output, $line;
	}
	Tools->printToFile($out,\@output);
	return scalar(@output);
}

sub cleanRowIDsforWGCNA {
	my $self=shift;
	$self->{perturbation}={};
	my @Rows=@{$self->{rowIDs}};
	my @Names;
	for(my$i=0;$i<=$#Rows;$i++){
		my @D=@{$self->{data}{$Rows[$i]}};
		my $n=$Rows[$i];
		$n=~s/\|/-/;
		$self->{perturbation}{$n}=\@D;
		push @Names, $n;
	}
	return \@Names;
}

sub logNTransform {
	my $self=shift;
	$self->{perturbation}={};
	my @Rows=@{$self->{rowIDs}};
	for(my$i=0;$i<=$#Rows;$i++){
		my @D=@{$self->{data}{$Rows[$i]}};
		my @P;
		foreach my $entry (@D){
			my $log=log($entry);
			push @P, $log;
		}
		$self->{perturbation}{$Rows[$i]}=\@P;
	}
	return \@Rows;
}

sub unLog2Transform {
	my $self=shift;
	$self->{perturbation}={};
	my @Rows=@{$self->{rowIDs}};
	for(my$i=0;$i<=$#Rows;$i++){
		my @D=@{$self->{data}{$Rows[$i]}};
		my @P;
		foreach my $entry (@D){
			my $log2;
			if($entry==0){
				$log2=0;
			}else{
				$log2=2**$entry;
			}
			push @P, $log2;
		}
		$self->{perturbation}{$Rows[$i]}=\@P;
	}
	return \@Rows;
}

sub log2Transform {
	my $self=shift;
	$self->{perturbation}={};
	my @Rows=@{$self->{rowIDs}};
	for(my$i=0;$i<=$#Rows;$i++){
		my @D=@{$self->{data}{$Rows[$i]}};
		my @P;
		foreach my $entry (@D){
			my $log2;
			if($entry==0){
				$log2=0;
			}else{
				$log2=Tools->log2($entry);
			}
			push @P, $log2;
		}
		$self->{perturbation}{$Rows[$i]}=\@P;
	}
	return \@Rows;
}

sub checkID {
	my $self=shift;
	my $id=shift;
	if(defined($self->{data}{$id})){
		return 1;
	}else{
		return 0;
	}
	return 0;
}

sub getDataByID_withZeros {
	my $self=shift;
	my $id=shift;
	my @output;
	if(defined($self->{data}{$id})){
		@output=@{$self->{data}{$id}};
		my $max = Tools->max(@output);
		warn "$id has all-zero entries\n" if $max == 0;
	}else{
		my @empty;
		my $N=scalar(@{$self->{data}{${$self->{rowIDs}}[0]}});
	#		warn "$N found to be the cardinality of ${$self->{rowIDs}}[0]\n";
#		warn "Could not find value for $id!\n";
		for(my$i=0;$i<$N;$i++){
			push @empty, 0;
		}
		@output=@empty;
	}
	return \@output;
}

sub getDataByID {
	my $self=shift;
	my $id=shift;
	my @output;
	if(defined($self->{data}{$id})){
		@output=@{$self->{data}{$id}};
		my $max = Tools->max(@output);
		warn "$id has all-zero entries\n" if $max == 0;
	}else{
		my @empty;
		my $N=scalar(@{$self->{data}{${$self->{rowIDs}}[0]}});
	#		warn "$N found to be the cardinality of ${$self->{rowIDs}}[0]\n";
#		warn "Could not find value for $id!\n";
		for(my$i=0;$i<$N;$i++){
			push @empty, 0;
		}
#		@output=@empty;
		return undef;
	}
	return \@output;
}

sub getDataByRowID {
	my $self=shift;
	my @Rows=@{$_[0]};
	my @output;
	foreach my $row (@Rows){
		if(defined($self->{data}{$row})){
			push @output, $self->{data}{$row};
		}else{
			my @empty;
			my $N=scalar(@{$self->{data}{${$self->{rowIDs}}[0]}});
			warn "Could not find value for $row!\n";
	#		warn "$N found to be the cardinality of ${$self->{rowIDs}}[0]\n";
			for(my$i=0;$i<$N;$i++){
				push @empty, 0;
			}
			push @output, \@empty;
		}
	}
	return \@output;
}

sub _checkRowValueRatio {
	my $minValue=shift;
	my $ratio=shift;
	my @D=@{$_[0]};
	my $N=scalar(@D);
	my $bad=0;
	for(my$i=0;$i<=$#D;$i++){
		$bad++ if $D[$i]<$minValue;
	}
	my $rat=$bad/$N;
	return 1 if $rat<=$ratio;
	return 0;
}

sub getHeader {
	my $self=shift;
	if(defined($self->{header})){
		return $self->{header};
	}else{
		die "No header found. did you load a file?\n";
	}
	die "Shouldn't have gotten here.\n";
}

sub setHeader {
	my $self=shift;
	my $ref=shift;
	$self->{header}=$ref;
	return 1;
}

sub addRow {
	my $self=shift;
	my $id  =shift;
	my $ref =shift;
	push @{$self->{rowIDs}}, $id;
	$self->{data}{$id}=$ref;
	return 1;
}

sub loadFile {
	my $self=shift;
	my $file=shift;
	my $del =shift;
	my @File=@{Tools->LoadFile($file)};
	my @Header=split($del,shift @File);
	my $junk=shift @Header;
	$self->{header}=\@Header;
	foreach my $line (@File) {
		my @line=split($del,$line);
		my $rowID=shift @line;
		if(defined($self->{data}{$rowID})){
			my $max = Tools->max(@{$self->{data}{$rowID}});
			my $t_max = Tools->max(@line);
			if($t_max > $max){
				$self->{data}{$rowID} = \@line;
#				warn "Encountered Duplicate Row IDs, Updating $rowID (stored maximum lt current maximum)\n";
			}else{
#				warn "Encountered Duplicate Row IDs. Not updating $rowID (stored maximum gt current max)\n";
			}
		}else{
			push @{$self->{rowIDs}}, $rowID;
			$self->{data}{$rowID}=\@line;
		}
	}
	$self->{currentRow}=0;
	return scalar(keys(%{$self->{data}}));
}

sub printFile {
	my $self=shift;
	my $out =shift;
	my $del =shift;
	my @output;
	push @output, join($del,@{$self->{header}});
	foreach my $row (@{$self->{rowIDs}}){
		die "Cannot find $row in data hash. Uh oh.\n\n" unless defined $self->{data}{$row};
		my $line=$row.$del.join($del,@{$self->{data}{$row}});
		push @output, $line;
	}
	Tools->printToFile($out,\@output);
	return scalar(@output);
}

sub _launchJob {
	my $command=shift;
	`$command`;
	return 1;
}

1;

