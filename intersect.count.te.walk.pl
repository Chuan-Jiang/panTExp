#!/usr/bin/perl
use warnings; use strict;
use List::Util qw(min max);
use Storable;
#use Bit::Vector;

my $prefix = shift @ARGV;
open DIS, ">$prefix.discarded.tsv" or die $!;

my $debug = 0;

my $strand_type = "n";

my @strand;
if($strand_type eq "f"){
	@strand = (0,1,-1);
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
		print  DIS "$read\tdiversity\n";
		next;
	}
	if( $cs =~ /^\+([ATCGatcg]+)/ ){
		if(length($1) > 20){
			print "long overhang, to read form here.. next\n" if $debug;
			print DIS "$read\t-$1\n";
			next;
		}
	}
	if($cs =~ /\+([ATCGatcg]+)$/){
		if(length($1) > 20){
			print "long overhang, to read form here.. next\n" if $debug;
			print DIS "$read\t+$1\n";
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
		print "$read no overlap with genes, put it in ouside..\n" if $debug;	
		#$r2g{"$read##$mapn"}{$which}{NS} = $nhr;
		$r2g{"$read##$mapn"}{$which}{AS} = $as;
		$r2g{"$read##$mapn"}{$which}{AN}{N} = 1;
		#foreach my $n (keys %{$nhr}){
		#	push @{$n2r{$n}}, "$read##$mapn";
		#}
		next;
	}

	my ($fn) = $bed[-1] =~ /^.(\d+)/;
	my ($sn) = $bed[-1] =~ /(\d+):\d+,\d+,\d+$/;

	$gene = join "__", $gene, ($fn,$sn);
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

	if(0){
	#if(%$susnr){
		if($r2g{"$read##$mapn"}{"$which"}{NS}){
			print "$read have two interwsect\n" if $debug;
			foreach my $n ( keys %{ $r2g{"$read##$mapn"}{"$which"}{NS} } ){
				unless(exists $$susnr{$n} ){
					delete $r2g{"$read##$mapn"}{"$which"}{NS}{$n};
				}
			}
		}else{
			$r2g{"$read##$mapn"}{"$which"}{NS} = $susnr;
		}
	}
	# if the same reads mapped on the gene two times, they will be assign two distince $mapn number
	$r2g{"$read##$mapn"}{"$which"}{AS} = $as;
	$r2g{"$read##$mapn"}{"$which"}{AN}{$gene} = $olen if $olen;
	
	#foreach my $n (keys %{$susnr}){
	#	push @{$n2r{$n}}, "$read##$mapn";
	#}
	$f2r{$gene}{"$read##$mapn"} = $which  if $olen; # have bedtool inserseciton, might have zero overlap

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
	my($path,$ps,$pe,$plen,$rlen,$read,$mapn,$which) = @_;
	my %nh;
	print "NodeHash: $path\n" if $debug;

	my @nodes = split /;/, $path;
	#my @nodes = $path =~ /(<|>)(\d+)/g ;
	for(my $i = 0; $i < @nodes; $i ++){
		my($n,$d,$l) = $nodes[$i] =~ /(\d+):(-?\d),(\d+)/;


		#sub keys: D, I, P, N, HE, TE
		# direciotn, index, Previour, Next, Headleft, Tailleft
		my @ni;
		$ni[0] = $d;
		$ni[1] = $i;
		my ($ss,$se) = (0,0);
		if($i == 0){
			if($d == 1){
				$ss = $ps;
			}else{
				$se = $ps;
			}
		}

		if($i == $#nodes ){
			if($d == 1){
				$se = $plen - $pe;
			}else{
				$ss = $plen - $pe;
			}
		}
		push @ni, ($ss,$se,$l);
		
		print "\t@ni\n" if $debug;
		$nh{$n} = \@ni;
	}

	return (\%nh);
}


sub check_overlap{
	# even sharing one node          , will be count as overlap. 
	my ($nhr,$gene,$which ) = @_;
	print "checking overlap with gene: $gene, read: $which\n" if $debug;	
	my $olen = 0; # total overlapped length with node annotation
	my %susp;

	my $nodesc = keys %{$nhr};   
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
				push @{$$nhr{$n}}, $nl;

				my ($rs,$re) =  ($$nhr{$n}[2], $$nhr{$n}[3]);
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



