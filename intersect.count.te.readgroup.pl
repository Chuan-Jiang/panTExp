#!/usr/bin/perl
use warnings; use strict;
use List::Util qw(min max);
use Storable;
#use Bit::Vector;

my $debug = 0;

my $strand_type = "n";

my @strand;
if($strand_type eq "f"){
	@strand = (0, 1,-1);
}elsif($strand_type eq "r") {
	@strand = (0,-1,1);
}else{
	@strand = (0,0,0);
}


##### annotations   ###
my %anno;

# three situation

#my %inside;
#my %suspend;
#my %outside;
my %n2r;
my %r2g;
my %f2r;

my %bad;
my %rgraph;

while(<>){
	chomp;
	print "\n$_\n" if $debug ;
	my @ar = split /\t/;

	my @map = @ar[0..15];
	my @bed = @ar[16..$#ar -1];

	#my @bed = @ar[0..7];
	#my @map = @ar[8..$#ar-1];

	print ">>BED: @bed\n" if $debug ;
	print ">>MAP: @map\n" if $debug ;

	########################################
	# get informaton from reads and anno ###
	# ######################################
	
	#VVVVVVVVVVVVVVVVVVVVVVV
	## from map	
	my @pair = split /,/, $map[12];
	my $which = $pair[0] == 2? 1 : 2;		
	my $read = $map[3];
	my $rlen = $map[4];
	my $dv = $map[5];
	my $as = $map[6];
	my $cs = $map[7];
	my $mpath = $map[8];
	my $plen = $map[9];
	my $phead = $map[10];
	my $ptail = $map[11];

	my $mapn = $map[13];
	my $chkn = $map[14];

	############################################
	### remove reads that have bad alignment ###
	if( $dv eq "NA" or $dv > 0.2 ){
		print "large diversity, next ..\n" if $debug;
		next;
	}
	if( $cs =~ /^\+([ATCGatcg]+)/ ){
		if(length($1) > 20){
			print "long overhang, to read form here.. next\n" if $debug;
			next;
		}
	}
	if($cs =~ /\+([ATCGatcg]+)$/){
		if(length($1) > 20){
			print "long overhang, to read form here.. next\n" if $debug;
			next;
		}
	}

	### node hash
	# can make a univeral graph, reduce memory . 
	my ($nhr) = nodehash($mpath,$phead,$ptail,$plen,$rlen,$read,$mapn,$which);

	#VVVVVVVVVVVVVVVVVVVVVVVV
	# from anno	
	#######################
	my $gene = $bed[3];
	
	# no overlap, then put in into %ouside <<<<<<< type 3
	if($gene eq "."){
		print "$read no overlap with genes, put it in outside..\n" if $debug;	
		$r2g{"$read##$mapn"}{$which}{NS} = $nhr;
		$r2g{"$read##$mapn"}{$which}{AS} = $as;
		foreach my $n (keys %{$nhr}){
			push @{$n2r{$n}}, "$read##$mapn";
		}
		next;
	}

	my @end_ns = $bed[-1] =~ /^.(\d+).+?(\d+):\d+,\d+,\d+$/;
	$gene = join "__", $gene,@end_ns;
	#print "new gene name $gene\n" if $debug;
	## extract annos , save involved reads
	unless($anno{$gene} ){
		my @nodes = split /;/, $bed[-1];
		#print "N:@nodes\n" if $debug;
		foreach my $n (@nodes){
			my($d,$id,$l,$s,$e) = $n =~ /(.)(\d+):(\d+),(\d+),(\d+)/;
			if($d eq ">"){
				$d = 1;
			}elsif($d eq "<"){
				$d = -1;
			}else{
				$d = 0;
			}
			#print "$gene $id $d $s $e\n" if $debug;
			$anno{$gene}{$id}{D} = $d;
			$anno{$gene}{$id}{L} = $l;
			$anno{$gene}{$id}{S} = $s;
			$anno{$gene}{$id}{E} = $e;
		}
	}	
	print "overlap test $read...\n" if $debug;
	#next if (defined $anchored{$gene}{$read}{$mpatch} or defined $anchored{$gene}{$read}{$pair[1]});
	# to check whether reads are fully overlapped with annotation and also its maping direction
	
	my ($olen, $nc,$susnr)  = check_overlap($nhr,$gene,$which);

	if($olen < 5){
		print "outside of annotations.\n" if $debug;
		$r2g{"$read##$mapn"}{"$which"}{NS} = $nhr;         
		$r2g{"$read##$mapn"}{"$which"}{AS} = $as;
		foreach my $n (keys %{$nhr}){
			push @{$n2r{$n}}, "$read##$mapn";
		}
	}else{
		print "hit in gene body: $gene $read $mapn $which ..\n" if $debug; 
	
		if (%{$susnr} ){
			print "overhang, have suspend reads\n" if $debug;
			$r2g{"$read##$mapn"}{"$which"}{NS} = $susnr;         
			foreach my $n (keys %{$susnr}){
				push @{$n2r{$n}}, "$read##$mapn";
			}
		}
		$r2g{"$read##$mapn"}{"$which"}{AS} = $as;
		$r2g{"$read##$mapn"}{"$which"}{AN} = $gene;   

		$f2r{$gene}{"$read##$mapn"} = $which ;
	}
	my @anchors = keys %{$susnr};
	print "ANCHOR:$olen @anchors\n" if $debug;

}

store([\%n2r, \%r2g, \%f2r], 'anchors.dat');  


########################
# subfunctuinos  #

sub nodehash{
	my($path,$ps,$pe,$plen,$rlen,$read,$mapn,$which) = @_;
	my %nh;
	my @nodes = $path =~ /(<|>)(\d+)/g ;
	for(my $i = 0; $i < @nodes; $i = $i + 2){
		my $d = $nodes[$i] eq ">"? 1 : -1;
		my $n = $nodes[$i+1];

		#sub keys: D, I, P, N, HE, TE
		# direciotn, index, Previour, Next, Headleft, Tailleft
		my @ni;
		$ni[0] = $d;
		$ni[1] = $i/2;
		$ni[2] = @nodes/2 - 1;
		my ($ss,$se) = (0,0);
		if($i == 0){
			if($d == 1){
				$ss = $ps;
			}else{
				$se = $ps;
			}
			$ni[3] = 0;
		}else{
			$ni[3] = $nodes[$i-1];
		}

		if($i == $#nodes - 1){
			if($d == 1){
				$se = $plen - $pe;
			}else{
				$ss = $plen - $pe;
			}
			$ni[4]= 1;
		}else{
			$ni[4] = $nodes[$i+3];
		}
		@ni[5,6] = ($ss,$se);

		$nh{$n} = \@ni;
	}

	return (\%nh);
}


sub check_overlap{
	# even sharing one node          , will be count as overlap. 
	my ($nhr,$gene,$which ) = @_;
	print "checking overlap with gene: $gene, read: $which\n" if $debug;	
	my $olen = 0; # total overlapped length with node annotation
	my %susp ;

	my $nodesc = keys %{$nhr};  # count of nodes 
	# iterate each node
	foreach my $n ( sort { $$nhr{$a}[1] <=> $$nhr{$b}[1] } keys %{$nhr} ){
		my $d = $$nhr{$n}[0]; 
		print "\tN:$n D:$d\n" if $debug;    
		#print "Extract:$gene node:$n S:$anno{$gene}{$n}{S} E:$anno{$gene}{$n}{E} \n" if $debug;
		if( defined  $anno{$gene}{$n} ){ 
			if( $d * $strand[$which] * $anno{$gene}{$n}{D} >= 0 ){
				my $as = $anno{$gene}{$n}{S};
				my $ae = $anno{$gene}{$n}{E};
				my $nl = $anno{$gene}{$n}{L};
				
				my ($rs,$re) =  ($$nhr{$n}[5], $$nhr{$n}[6]);
			
				my $clen = 0;   
				if($rs == 0 and $re ==0  and $as == 0 and $ae == 0){
					$clen =  $nl;
					print "\t\tall full cover, len is node lenght $nl\n" if $debug;
				}else{
					if($rs + 1 > $nl - $ae or $as + 1 > $nl -$re){
						print "\t\toverlap is zero $nl $rs $re $as $ae\n" if $debug;
					}else{
						my $overlap = min($nl-$ae, $nl - $re) - max($rs, $as);
						my $or = $overlap / ($nl - $rs - $re);
						print "\t\tpartial overlap $nl $rs $re $as $ae $overlap $or\n" if $debug;
						$clen =  $overlap;
					}
				}
				if( $clen == 0 or $clen <  ($nl - $rs -$re)  ){
					print "\t\tparitial have overhang, $n  to susp\n" if $debug;
					$susp{$n} = $$nhr{$n}; # susp means read pair  that extend outside of annotation 
				}else{
					print "\t\tfull overlap. not suspends\n" if $debug;
				}
				$olen += $clen;
			}else{
				print "\tdirectino mismatch $n to suspends  $d * $strand[$which] * $anno{$gene}{$n}{D} \n" if $debug;
				$susp{$n} = $$nhr{$n};  
			}
		}else{
			print "\tsuspent end: $n...\n" if $debug;
			$susp{$n} = $$nhr{$n};
		}
	}
	return ($olen, $nodesc ,\%susp);
}



