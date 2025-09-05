#!/usr/bin/perl
use warnings; use strict;

open FAM, "/home/chuan/hdd1/TE_exp/Data/families.tsv" or die $!;
my %fams;
while(<FAM>){
	chomp;
	next if (/^#/);
	my @ar = split /\t/;
	my ($name) =  $ar[1] =~ /(.+)-int/;
	foreach my $ltr ( split /,/, $ar[2] ) {
		$fams{$name}{$ltr} = 1;
		$fams{$ltr}{$name} = 1;
	}
}
close FAM;

srand(12345);
my $max_extension = 300;  


my %rcd;
while(<>){
	chomp;

	my @ar = split /\s+/;
	if($ar[7] =~ /ERVK/){
		$ar[3] =~ s/#.+//g;
		$ar[3] =~ s/$ar[6]\|//;
		my $type;
		my $name ;
		if($ar[3] =~ /(.+)-int/){
			$type = "INT";
			$name = $1;
		}elsif($ar[3] =~ /LTR/){
			$name = $ar[3];
			$type = "LTR";
		}else{
			$name = $ar[3];
			$type = "Other";
		}
		#print "type $type\n";
		
		if($type ne "INT" and defined $rcd{$ar[9]}{N} ){
			1;
		}else{
			$rcd{$ar[9]}{N} = $name if $name;
		}

		$rcd{$ar[9]}{$type} ++;
		$rcd{$ar[9]}{CHR} = $ar[0];
		$rcd{$ar[9]}{L} += $ar[2] - $ar[1];
		$rcd{$ar[9]}{D} = $ar[5] ;

		if(defined $rcd{$ar[9]}{S}){
			$rcd{$ar[9]}{E} = $ar[2] ;
		}else{
			$rcd{$ar[9]}{S} = $ar[1];
			$rcd{$ar[9]}{E} = $ar[2];
		}
		#print "$ar[4]\t$rcd{$ar[14]}{S}\t$rcd{$ar[14]}{E}\t$ar[14]\n";


	}
}

my @ids = sort {$a <=> $b} keys %rcd;

for (my $i  = 0; $i < @ids; $i ++){
	my $id = $ids[$i];
	next unless ($rcd{$id});	
	my $lable;
	if( defined $rcd{$id}{INT} and defined $rcd{$id}{LTR} ){
		$lable = "STR";
	}elsif($rcd{$id}{L} > 1000 ){
		$lable = "LEN";
	}else{
		$lable = "PART";
	}

	if($lable ne "PART" ){
		my $name = $rcd{$id}{N} ? $rcd{$id}{N} : "HERVK";
		my $chr = $rcd{$id}{CHR};
		my $start = $rcd{$id}{S};
		my $end = $rcd{$id}{E};
		
		#print "GROW:$lable,$id,$start,$end,$rcd{$id}{N}\n" ;
		
		#print "look behind\n";
		for (my $j = 1; ; $j ++){
			my $ij = $i - $j;
			if($ij <  0){
				#print "index exhausted\n";
				last;
			}
			my $lid = $ids[$ij];

			next unless $rcd{$lid};
			#print "$i $j $ij $id $lid $rcd{$lid}{N}\n" ;

			my $lend = $rcd{$lid}{E} ;
			if($start - $lend > 500 or  $start - $lend < 0){
				#print "long..\n";
				last;
			}else{
				#print "close..\n";
			}

			my $lnm = $rcd{$lid}{N};
			if($rcd{$lid}{CHR} eq $chr ){
				foreach my $pair (keys %{$fams{$name}} ){
					if($lnm =~ /$pair/){
						#print "join $lid\n";
						$start = $rcd{$lid}{S};
						delete $rcd{$lid};
						last;
					}
				}
			}
		}
		#print "look foreard\n";
		for (my $j = 1; ; $j ++){
			my $ij = $i + $j;
			if($ij > @ids - 1 ){
				last;
			}
			
			my $lid = $ids[$ij];
			
			next unless ($rcd{$lid});
			#print "$i $j $ij $id $lid $rcd{$lid}{N}\n" ;

			my $lstart = $rcd{$lid}{S};
			if( $lstart - $end > 500 or $lstart - $end < 0 ){
				#print "long...\n";
				last;
			}

			my $lnm = $rcd{$lid}{N};
			if($rcd{$lid}{CHR} eq $chr ){
				foreach my $pair (keys %{$fams{$name}} ){
					if($lnm =~ /$pair/){
						#print "join $lid\n";
						$end = $rcd{$lid}{E};
						delete $rcd{$lid};
						last;
					}
				}
			}
		}
		my $extend_left  = int(rand($max_extension));
		my $extend_right = int(rand($max_extension));
		$start = $start - $extend_left;
		$start = 0 if $start < 0;

		$end = $end + $extend_right;

		$chr =~ s/GRCh38#0#//;
		print "$chr\t$start\t$end\t$name:$id:$lable\t.\t$rcd{$id}{D}\n";
		print STDERR  "$chr\tsimulate\tgene\t$start\t$end\t.\t$rcd{$id}{D}\t.\tgene_id \"$name:$id:$lable\"; locus \"$name:$id:$lable\"\n";
		delete $rcd{$id};
	}
}

#exp_rand{
#	my $lambda = 2;
#	my $exp_random = -log(rand()) / $lambda * 100;
#}
