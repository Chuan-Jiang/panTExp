#!/usr/bin/perl
use warnings; use strict;
use List::Util qw(min max);
use Storable;
#use Bit::Vector;

my $debug = 1;

my $strand_type = "r";

my $SD ;
if($strand_type eq "f"){
	$SD = 1;
}elsif($strand_type eq "r") {
	$SD = -1;
}else{
	$SD = 0;
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

	my @map = @ar[0..7];
	my @bed = @ar[8..$#ar -1];

	#my @bed = @ar[0..7];
	#my @map = @ar[8..$#ar-1];

	print ">>BED: @bed\n" if $debug ;
	print ">>MAP: @map\nvvvvvvvvvvvv\n" if $debug ;
	########################################
	# get informaton from reads and anno ###
	# ######################################
	
	#VVVVVVVVVVVVVVVVVVVVVVV
	## from map	
	my $read = $map[3];
	my $as = $map[4];
	my $path = $map[5];
	my $mapn = $map[6];
	my $chkn = $map[7];

	### node hash
	# can make a univeral graph, reduce memory . 
	
	my ($nhr) = nodehash($path);

	#VVVVVVVVVVVVVVVVVVVVVVVV
	# from anno	
	#######################
	my $gene = $bed[3];
	
	# no overlap, then put in into %ouside , using for extension  <<<<<<< type 3
	if($gene eq "."){
		print "$read no overlap with genes, put it in ouside..\n" if $debug;	
		$r2g{"$read##$mapn"}{AS} = $as;
		$r2g{"$read##$mapn"}{AN}{N} = 1;
		#foreach my $n (keys %{$nhr}){
		#	push @{$n2r{$n}}, "$read##$mapn";
		#}
		next;
	}

	my ($fn) = $bed[-1] =~ /^.(\d+)/;
	my ($sn) = $bed[-1] =~ /(\d+):\d+,\d+,\d+$/;

	## extract annos , save involved reads
	unless($anno{$gene} ){
		my @nodes = split /;/, $bed[-1];
		#print "N:@nodes\n" if $debug;
		foreach my $n (@nodes){
			my($d,$id,$l,$s,$e) = $n =~ /(.)(\d+):(\d+),(\d+),(\d+)/;
			if($d eq "+"){
				$d = 1;
			}elsif($d eq "-"){
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
	
	my ($olen, $nc,$susnr)  = check_overlap($nhr,$gene);
	print "overlap results:$olen, $nc\n" if $debug;

	# if the same reads mapped on the gene two times, they will be assign two distince $mapn number
	$r2g{"$read##$mapn"}{AS} = $as;
	$r2g{"$read##$mapn"}{AN}{$gene} = $olen if $olen;
	
	#foreach my $n (keys %{$susnr}){
	#	push @{$n2r{$n}}, "$read##$mapn";
	#}

	#my @anchors = keys %{$susnr};
	#print "ANCHOR:$olen @anchors\n" if $debug;

}

#store(\%n2r, 'n2r.dat');
store(\%r2g, "r2g.dat");
store(\%f2r, "f2r.dat");

#store([\%n2r, \%r2g, \%f2r], 'anchors.dat');  


########################
# subfunctuinos  #

sub nodehash{
	#my($path,$ps,$pe,$plen,$rlen,$read,$mapn,$which) = @_;
	my ($path) = @_;
	my %nh;
	print ">>NodeHash: $path\n" if $debug;

	my @gs = split /\|/, $path;
	my @nm = ("R1","M","R2");
	my $ni = 0;
	for(my $j = 0; $j < @nm; $j ++){
		my @nodes = split /;/, $gs[$j];
		my $gn = $nm[$j];
		
		#my @nodes = $path =~ /(<|>)(\d+)/g ;
		for( my $i = 0; $i < @nodes; $i ++ ){
			my($d,$id,$l,$s,$e) = $nodes[$i]  =~ /(.)(\d+):(\d+),(\d+),(\d+)/;
			if($d eq "+"){
				$d = 1;
			}elsif($d eq "-"){
				$d = -1;
			}else{
				$d = 0;
			}
			my @ninf  = ($d,$ni,$gn,$l,$s,$e);
			print "SPLIT: $id, $d $ni @ninf\n" if $debug;

			$nh{$id} = \@ninf;
			$ni ++;
			#sub keys: D, I, P, N, HE, TE
			# direciotn, index, Previour, Next, Headleft, Tailleft
		}
	}
	return (\%nh);
}


sub check_overlap{
	# even sharing one node          , will be count as overlap. 
	my ($nhr,$gene ) = @_;
	print "checking overlap with gene: $gene\n" if $debug;	
	my $olen = 0; # total overlapped length with node annotation
	my %susp;

	my %nodeh = %{$nhr};	
	my $nodesc = keys %nodeh;   
	# iterate each node
	foreach my $n ( sort { $nodeh{$a}[1] <=> $nodeh{$b}[1] } keys %nodeh ){
		my $d = $$nhr{$n}[0]; 
		my $ni = $$nhr{$n}[1];
		print "\tN:$n D:$d $ni\n" if $debug;    
		#print "Extract:$gene node:$n S:$anno{$gene}{$n}{S} E:$anno{$gene}{$n}{E} \n" if $debug;
		if( defined  $anno{$gene}{$n} ){ 
			if( $d * $SD  * $anno{$gene}{$n}{D} >= 0 ){
				my $as = $anno{$gene}{$n}{S};
				my $ae = $anno{$gene}{$n}{E};
				my $nl = $anno{$gene}{$n}{L};

				my ($rs,$re) =  ($$nhr{$n}[4], $$nhr{$n}[5]);
				
				my $clen = 0;   
				if($rs == 0 and $re ==0  and $as == 0 and $ae == 0){
					$clen =  $nl;
					print "\t\tall full cover, len is node lenght $nl\n" if $debug;
				}else{
					if($rs + 1 > $nl - $ae or $as + 1 > $nl -$re){
						print "\t\toverlap is zero $nl $rs $re $as $ae\n" if $debug;
					}else{
						my $overlap = min($nl-$ae, $nl - $re) - max($rs, $as);
						print "\t\tpartial overlap $nl $rs $re $as $ae $overlap\n" if $debug;
						my $or = $overlap / ($nl - $rs - $re);
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
				print "\tdirectino mismatch $n to suspends  $d * $SD * $anno{$gene}{$n}{D} \n" if $debug;
				$susp{$n} = $$nhr{$n};  
			}
		}else{
			print "\tsuspent end: $n...\n" if $debug;
			$susp{$n} = $$nhr{$n};
		}
	}
	return ($olen, $nodesc ,\%susp);
}



