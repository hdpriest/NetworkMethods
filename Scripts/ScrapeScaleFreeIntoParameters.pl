#!/usr/bin/perl
use warnings;
use strict;

use lib '/home/hpriest/Scripts/Library';
use Tools;

my @list=@ARGV;
my %cuts;
my $min_r = 0.8; ## minimum correlation to log-log
my $max_m = -0.8; ## maximum regression slope

my $java = "/home/hpriest/Java/Install/jre1.7.0_65/bin/java -jar ";
my $jar  = "/home/hpriest/Java/NewPlas/NetworkDev.jar ";
my $cmd  = "permute ";
my $mL = 0.65;
my $mH = 0.96;
my $aH = 26;
my $aL = 16;
my $threads = 80;
my $options = "-c pcc ";
for(my $i=0;$i<=$#Dir;$i++){

print "Permutations,PositiveCut,NegativeCut\n";
foreach my $file (@list){
	my %info = %{_parseSF($file)};
	my @dir = ("Pos","Neg");
	my @types = ("Corr","Slope");
	my @m = keys %info;
	my @a = keys %{$info{$m[0]}};
	foreach my $dir (@dir){
		$cuts{$file}={};
		$cuts{$file}{$dir}={};
		my ($final_A,$final_M);
		for(my $m=0;$m<=$#m;$m++){
			my $M=$m[$m]; ## confused yet?
			for(my $a=0;$a<=$#a;$a++){
				next if((defined $final_A) && (defined $final_M));
				my $A=$a[$a];
#				$I{$mu}{$A}{$direction}{$type}=$V;
				my $r = $info{$M}{$A}{$dir}{$types[0]};
				my $s = $info{$M}{$A}{$dir}{$types[1]};
				if(($r>=$min_r) && ($s<=$max_m)){
					$final_A=$A;
					$final_M=$M;
				}
			}
		}
		if((defined $final_A) && (defined $final_M)){
			$cuts{$file}{$dir}{'A'}=$final_A;
			$cuts{$file}{$dir}{'M'}=$final_M;
		}else{
			warn "Could not find cutoff for file: $file\n";
		}
	}
}

foreach my $file (keys %cuts){
	my $p_string = " -pa ".$cuts{$file}{"Pos"}{'A'}." -pm ".$cuts{$file}{"Pos"}{'M'};
	my $n_string = " -na ".$cuts{$file}{"Neg"}{'A'}." -nm ".$cuts{$file}{"Neg"}{'M'};
}

sub _parseSF {
	my $file=shift;
	my @file=@{Tools->LoadFile($file)};
	my $type;
	my $direction;
	my @alphas;
	my %I;
	foreach my $line (@file){
		if($line=~m/^Negative/){
			$direction="Neg";
		}elsif($line=~m/^Positive/){
			$direction="Pos";
		}elsif($line=~m/^Correlation/){
			$type="Corr";
		}elsif($line=~m/^Regression/){
			$type="Slope";
		}elsif($line=~m/^Mean/){
			$type="K";
		}elsif($line=~m/^,/){
			@alphas = split(/\,/,$line);
			shift @alphas;
		}elsif($line=~m/^[01]/){
			my @values = split(/\,/,$line);
			my $mu = shift @values;
			for(my $a=0;$a<=$#alphas;$a++){
				my $V = $values[$a];
				my $A = $alphas[$a];
				if(defined($I{$mu})){
					if(defined($I{$mu}{$A})){
						if(defined($I{$mu}{$A}{$direction})){
							if(defined($I{$mu}{$A}{$direction}{$type})){
								die "you parsed wrong\n";
							}else{
								$I{$mu}{$A}{$direction}{$type}=$V;
							}
						}else{
							$I{$mu}{$A}{$direction}={};
							$I{$mu}{$A}{$direction}{$type}=$V;
						}
					}else{
						$I{$mu}{$A}={};
						$I{$mu}{$A}{$direction}={};
						$I{$mu}{$A}{$direction}{$type}=$V;
					}
				}else{
					$I{$mu}={};
					$I{$mu}{$A}={};
					$I{$mu}{$A}{$direction}={};
					$I{$mu}{$A}{$direction}{$type}=$V;
				}
			}
		}else{
		}
	}
	return \%I;
}
