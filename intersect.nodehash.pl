#!/usr/bin/perl
use warnings; use strict;
use List::Util qw(min max);
use Storable;
use Bit::Vector;

my $debug = 1 ;
$| = 1;

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

# three situation

#my %inside;
#my %suspend;
#my %outside;
my %f2r;


my %links;
my %anno;
my %r2g;
my $closed = 0;
#################
my $count = 0;

while(<>){
	chomp;
	if($count % 100000 == 0){
		warn  "LINE $count\n";
	}
	$count ++;

	my @ar = split /\t/;

	my @map = @ar[0..8];
	my @bed = @ar[9..$#ar -1];

	my $gene = $bed[3];	
	########################################
	# get informaton from reads and anno ###
	# ######################################
	my $pmax =  max(map { abs($_) } keys %links);
	if( $pmax and $map[1] > $pmax and !$anno{$gene} ){
		print "\nfound a disconnection... $map[1] $pmax\n" if $debug ;
		
		## check_overlap
		check_overlap(\%links,\%anno);

		undef %anno;
		undef %links;
	}else{
		if($debug){
			print "continue collect reads:\n";
			print "$pmax $map[1] $anno{$gene}\n";
		}
	}
	print ">>BED: @bed\n" if $debug ;
	print ">>MAP: @map\nvvvvvvvvvvvv\n" if $debug ;
	#VVVVVVVVVVVVVVVVVVVVVVV
	## from map	
	my $as = $map[3];
	my $path = $map[4];
	my($read,$mapn) = $map[5] =~ /(.+):(\d+)$/;
	my $sample = $map[6];
	my ($chkn,$chkt) = @map[7,8];

	### node hash
	# can make a univeral graph, reduce memory . 

	my ($nhr) = nodehash($path, $map[1], $map[2], "$read##$mapn##$sample", \%links);

	#VVVVVVVVVVVVVVVVVVVVVVVV
	# from anno	
	#######################
	$r2g{"$read##$mapn"}{AS} = $as;
	# no overlap, then put in into %ouside , using for extension  <<<<<<< type 3	
	if($gene eq "."){
		next;
	}
	## extract annos , save involved reads
	my @nodes = split /;/, $bed[-1];
	my $gl = $bed[-3];

	#print "N:@nodes\n" if $debug;
	foreach my $n (@nodes){
		my($id,$l,$s,$e,$t) =  split /,/, $n ; #=~ /(.)(\d+):(\d+),(\d+),(\d+)/;
		
		my $d;	
		if($id <  0 ){
			$d = -1;
		}elsif($id > 0 ){
			$d = 1;
		}else{
			$d = 0;
		}
		$id = abs($id);
		print "save: $gene $id $d $s $e\n" if $debug;
		$anno{$gene}{$id}{D} = $d;
		$anno{$gene}{$id}{L} = $gl;
		$anno{$gene}{$id}{S} = $s;
		$anno{$gene}{$id}{E} = $e;
	}
}

# analys the last records
check_overlap(\%links,\%anno);

########################
# subfunctuinos  #
sub check_overlap{
	# even sharing one node          , will be count as overlap. 
	my ($links, $anno) = @_;
	print "checking overlap \n" if $debug;	
	unless ( keys %{ $anno }){
		print "no anno... skip\n\n" if $debug;
		return;
	}
	
	my %olen;  # reads overlapped??? with annotation	
	my %gc; # T = total genes overlapped. M = directionn mathed genes overlap


	foreach my $id ( sort { $a  <=> $b } keys %{$links} ){
		my @rs = @{$$links{$id}{R}};
		my @pns = keys %{$$links{$id}{P}};
		my @nns = keys %{$$links{$id}{N}};
		my $idc = scalar @rs;
		if($debug){
			print "ID $id C  $idc = @rs\n";
			print "ID $id P @pns\n";
			print "ID $id N @nns\n";
		}
		my $ida =  $id ;
		foreach my $g ( keys %$anno ){
			if( exists $$anno{$g}{$ida} ) {
				$gc{T}{$g} = 1;	
				if( $$anno{$g}{$ida}{D} * $SD * $id >= 0){ ## test anno node direction
					$gc{M}{$g} = 1;
					print "\t$g overlap $ida $$anno{$g}{$id}{D} * $SD * $id\n" if $debug;
					my $count = 0; ## count of overlapped positions
					my $total = 0; ## total coverage, that is depth
				
					my $nl = @{$$links{$id}{A}};
					if($debug){
						print "range $nl $ida = $$anno{$g}{$id}{S}.. $nl - $$anno{$g}{$id}{E} - 1\n";
						print "nodearray @{$$links{$id}{A}}\n";
					}	
					my @rix; # have all reads overlapped with gene
					foreach my $ai ( $anno{$g}{$id}{S}.. $nl - $$anno{$g}{$id}{E} - 1 ){ # iterate all overlapped position with gene
						my $v = $$links{$id}{A}[$ai];
						my $vl = @{$v}; ## all reads indexes 
						print "INDEX $ai lenth $vl  @$v\n" if $debug ;
						#        ------
						#        --     count will be 6. total will be 11. 
						#        ---
						#    ----------------------------- 
						if($vl > 0 ){
							$count ++; # how many position overlapped
						} 
						$total += $vl;  # total coverage 
						
						push @rix, @{$v};
					}
					@rix = uniq(@rix);

					$olen{$g}{$id} = [$count,$total,\@rix];	
					print "\t\tMATCH: count $count total $total\n" if $debug ;

				}else{
					print "\t$g overlap, but in reverse strand...\n" if $debug ;
				}
			}else{
				print "\t $g no overlaop with node $g $ida $id\n" if $debug;
			}
		}
	}

	if($debug){
		print "Total gene : ", keys %{$gc{T}}, "\n";
		print "Match gene : ", keys %{$gc{M}}, "\n";
	}

	print "findoverlaps:\n" if $debug ;
	foreach my $g (keys %olen){
		foreach my $id (keys %{$olen{$g}}){
			print "$g $id @{$olen{$g}{$id}} @{$links{$id}{R}}[@{$olen{$g}{$id}[2]}]\n";
		}
	}
	print "checking overlap finished...\n\n" if $debug ;
}

###############################
###############################
sub nodehash{
	#my($path,$ps,$pe,$plen,$rlen,$read,$mapn,$which) = @_;
	my ($path, $sn, $en,$read,$links) = @_;
	print "\t>>NodeHash: $path\n" if $debug;

	my @nodes = split /;/, $path;
	#my @nodes = $path =~ /(<|>)(\d+)/g ;
	@nodes = map{ [split /,/]  } @nodes;
	my @units;
	for( my $i = 0; $i < @nodes; $i ++ ){

		my($id,$l,$s,$e,$t) =  @{$nodes[$i]}  ; 
		#######################
		###### find links #####
		my($pn,$nn);
		if($i > 0){
			($pn) = $nodes[$i-1] -> [0]; 
		}else{
			$pn = 0;
		}
		if($i < $#nodes){
			($nn) = $nodes[$i + 1] -> [0] ;
		}else{
			$nn = 0;
		}
		##########
		$$links{$id}{P}{$pn} ++;
		$$links{$id}{N}{$nn} ++;
		push @{$$links{$id}{R}}, $read;
		my $ri = $#{$$links{$id}{R}} ; 	
		#print "read index $ri for node $id\n" if $debug ;
		if(exists $$links{$id}{A}){
			my $aref = $$links{$id}{A};
			push @{$_}, $ri for @{$aref}[$s .. $l - $e - 1];
			#$_++ for @{$$links{$id}{A}}[$s..$l-$e-1];
		}else{
			#my @na = (0) x $l;
			#$_ ++ for @na[$s..$l-$e-1];
			my @na = map { [] } 1..$l;
			push @{$_}, $ri for @na[$s..$l-$e-1];
			$$links{$id}{A} = \@na;
		}
		#print "\tUPDATE:$id @{$$links{$id}{A}} \n";
		###### current nodes to hash:
		#if( abs($id) > $en or abs($id) < $sn){
		#	next;
		#}
	}
}

sub uniq {
	my %seen;
    return grep { !$seen{$_}++ } @_;
}

sub array_range {
	my ($min, $max) = (undef, undef);
	for my $val (@_) {
		$min = $val if !defined($min) || $val < $min;
		$max = $val if !defined($max) || $val > $max;
	}
	return ($min, $max);
}
