#!/usr/bin/perl
use warnings; use strict;
use List::Util qw(min max);
use Storable;
use Bio::DB::HTS::Tabix;
use IO::Handle;
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

##### Unions ######

my %sets;
my %union_count;

##### annotations   ###
my %anno;

# three situation

#my %inside;
#my %suspend;
#my %outside;
my %r2g;
my %f2r;

my %bad;
my %rgraph;

my %nodeh;
my %links;

my %seen;

my $count = 0;

while(<>){
	chomp;
	if($count % 100000 == 0){
		print "LINE $count\n";
	}
	$count ++;

	print "\n$_\n" if $debug ;
	my @ar = split /\t/;

	my @map = @ar[0..8];
	my @bed = @ar[9..$#ar -1];

	#my @bed = @ar[0..7];
	#my @map = @ar[8..$#ar-1];

	print ">>BED: @bed\n" if $debug ;
	print ">>MAP: @map\nvvvvvvvvvvvv\n" if $debug ;
	########################################
	# get informaton from reads and anno ###
	# ######################################
	
	#VVVVVVVVVVVVVVVVVVVVVVV
	## from map	
	my $as = $map[3];
	my $path = $map[4];
	my ($read,$mapn)  = $map[5] =~ /(.+):(\d+)$/;
	my $sample = $map[6];
	my ($chkp,$chkn) = @map[7,8];

	### node hash
	# can make a univeral graph, reduce memory . 

	my ($nhr) = nodehash($path,$map[1],$map[2], $read);

	#VVVVVVVVVVVVVVVVVVVVVVVV
	# from anno	
	#######################
	my $gene = $bed[3];
	
	# no overlap, then put in into %ouside , using for extension  <<<<<<< type 3
	if($gene eq "."){
		print "$read no overlap with genes, put it in ouside..\n" if $debug;	
		$r2g{"$read##$mapn"}{AS} = $as;
		$r2g{"$read##$mapn"}{AN}{N} = 1;
		#$r2g{"$read##$mapn"}{NS} = $nhr;
		next;
	}
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
	# to check whether reads are fully overlapped with annotation and also its maping direction
	
	my ($olen, $nc,$susnr)  = check_overlap($nhr,$gene);
	print "overlap results:$olen, $nc\n" if $debug;
	
	# if the same reads mapped on the gene two times, they will be assign two distince $mapn number
	$r2g{"$read##$mapn"}{AS} = $as;
	$r2g{"$read##$mapn"}{AN}{$gene} = $olen if $olen;
	#$r2g{"$read##$mapn"}{NS} = $susnr;

	$f2r{$gene}{"$read##$mapn"} = $as  if $olen;
	#foreach my $n (keys %{$susnr}){
}

store(\%sets,"sets.dat");
store(\%union_count,"union_count.dat");

store(\%r2g, "r2g.dat");
store(\%f2r, "f2r.dat");

########################
# subfunctuinos  #

sub nodehash{
	#my($path,$ps,$pe,$plen,$rlen,$read,$mapn,$which) = @_;
	my ($path, $sn, $en,$read) = @_;
	my %nh;
	print ">>NodeHash: $path\n" if $debug;

	my @nodes = split /;/, $path;
	#my @nodes = $path =~ /(<|>)(\d+)/g ;
	
	my @units;
	for( my $i = 0; $i < @nodes; $i ++ ){
		my($id,$l,$s,$e,$t) = $nodes[$i]  =~ /(.\d+):(\d+),(\d+),(\d+),(\d+)/;

		## union
		my @us;
		if($l > 20){
			my $fras = int($l / 20)  + 1;
			my $start = int($s/20);
			my $distance = int(($l - ($start - 1) * 20 -$e)/20)  ;
			for(my $i = 0;$i < $distance; $i ++){
				my $indx = $start + $i;
				my $pos = $indx * 20;
				print "U:$id $start, $indx  $pos  $distance\n" if $debug;
				push @us, "$id:$indx"
			}
			if($id < 0){
				@us = reverse (@us);
			}
		}elsif($l > 10) {
			@us = ($id);
		}
		push @units, @us;

		if( abs($id) > $en or abs($id) < $sn){
			next;
		}

		my @ninf  = ($i,$l,$s,$e,$t);
		print "Nodes: $id @ninf\n" if $debug;
		

		$nh{$id} = \@ninf;
		#sub keys: D, I, P, N, HE, TE
		# direciotn, index, Previour, Next, Headleft, Tailleft
	}
	if($seen{$read}){
		next;
	}else{
		for(my $i = 0 ; $i < @units - 1 ; $i ++){
			union($units[$i], $units[$i + 1] );
		}
	}
	###
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
	foreach my $n ( sort { $nodeh{$a}[0] <=> $nodeh{$b}[0] } keys %nodeh ){
		#print "Extract:$gene node:$n S:$anno{$gene}{$n}{S} E:$anno{$gene}{$n}{E} \n" if $debug;
		my $n_t = abs($n);
		if( defined  $anno{$gene}{$n_t} ){
			print "$n $SD $anno{$gene}{$n_t}{D} ..\n" if $debug;
			if( $n * $SD  * $anno{$gene}{$n_t}{D} >= 0 ){
				my $as = $anno{$gene}{$n_t}{S};
				my $ae = $anno{$gene}{$n_t}{E};
				my $nl = $anno{$gene}{$n_t}{L};

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
				print "\tdirectino mismatch $n to suspends  $n * $SD * $anno{$gene}{$n_t}{D} \n" if $debug;
				$susp{$n} = $$nhr{$n};  
			}
		}else{
			print "\tsuspent end: $n...\n" if $debug;
			$susp{$n} = $$nhr{$n};
		}
	}
	return ($olen, $nodesc ,\%susp);
}




#### union-finder

# Find the root of a node with path compression

sub find {
    my ($x) = @_;
    $sets{$x} = $x unless exists $sets{$x};
    return $x if $sets{$x} eq $x;
    $sets{$x} = find($sets{$x});
    return $sets{$x};
}

# Perform union and count number of union attempts
# Only merge when the count reaches a threshold (default: 3)
sub union {
	my ($a, $b, $threshold) = @_;
    $threshold //= 3;  # default to 3 if not given
	print "UNI:$a, $b, $threshold\n" if $debug ;

    my $ra = find($a);
    my $rb = find($b);
	if($ra eq $rb){
		print "$a and $b in same group\n" if $debug;
		return;
	}

    my $key = join('|', sort ($ra, $rb));
    $union_count{$key}++;

    if ($union_count{$key} >= $threshold) {
        print "$ra and $rb to same sets\n" if $debug;
		$sets{$rb} = $ra;
    }else{
		print "not excess $threshold\n" if $debug;
	}
}

# Get the number of times a pair has been unioned
sub union_weight {
    my ($a, $b) = @_;
    my $key = join('|', sort (find($a), find($b)));
    return $union_count{$key} || 0;
}

# Print groupings (for debug)
sub print_groups {
    my %groups;
    for my $x (keys %sets) {
        push @{ $groups{find($x)} }, $x;
    }
    for my $root (keys %groups) {
        print "Group [$root]: @{$groups{$root}}\n";
    }
}

