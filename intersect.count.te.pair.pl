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

##### Unions ######

my %sets;
my %union_count;

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
	
	my ($nhr) = nodehash($path,$map[1],$map[2], $read,\%n2r);

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

store(\%n2r, 'n2r.dat');
store(\%r2g, "r2g.dat");
store(\%f2r, "f2r.dat");

########################
# subfunctuinos  #

sub nodehash{
	#my($path,$ps,$pe,$plen,$rlen,$read,$mapn,$which) = @_;
	my ($path, $sn, $en,$read, $n2r) = @_;
	my %nh;
	print ">>NodeHash: $path\n" if $debug;

	my @nodes = split /;/, $path;
	#my @nodes = $path =~ /(<|>)(\d+)/g ;
	for( my $i = 0; $i < @nodes; $i ++ ){
		my($id,$l,$s,$e,$t) = $nodes[$i]  =~ /(.\d+):(\d+),(\d+),(\d+),(\d+)/;

		## union
		


		if( abs($id) > $en or abs($id) < $sn){
			next;
		}

		my @ninf  = ($i,$l,$s,$e,$t);
		print "Nodes: $id @ninf\n" if $debug;
		
		print "N2R: $id  = $read $l $s $e \n" if $debug ;	
		push @{$$n2r{$id}}, [ $read, $l, $s, $e ];

		$nh{$id} = \@ninf;
		#sub keys: D, I, P, N, HE, TE
		# direciotn, index, Previour, Next, Headleft, Tailleft
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

    my $ra = find($a);
    my $rb = find($b);

    return if $ra eq $rb;  # already in the same group

    my $key = join('|', sort ($ra, $rb));
    $union_count{$key}++;

    if ($union_count{$key} >= $threshold) {
        $sets{$rb} = $ra;
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

