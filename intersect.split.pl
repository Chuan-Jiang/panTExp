#!/usr/bin/perl
use warnings; use strict;
use List::Util qw(min max);
use Storable;
use Bio::DB::HTS::Tabix;

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



my %f2r;

my %links;
my %anno;
my %r2g;
my $closed = 0;

my $noover = "../02_quant/CF1.NoOverlap.tsv.gz";
#my $tabix = Bio::DB::HTS::Tabix->new(filename => $nodenei_f);            
#my $iter = $tabix -> query($win); 

#################

my $count = 0;

while(<STDIN>){
	chomp;
	if($count % 100000 == 0){
		print "LINE $count\n";
	}
	$count ++;

	my @ar = split /\t/;

	my @map = @ar[0..7];
	my @bed = @ar[8..$#ar -1];

	########################################
	# get informaton from reads and anno ###
	# ######################################
	#my $pmax =  max(map { abs($_) } keys %links);
	#if($map[1] > $ran_aln[1] ){
	#	print "\nfound a disconnection... $map[1] $pmax\n";
		
	#	my @ran_link = array_range(keys %links);
		
		## check_overlap
		#	check_overlap(\%links,\%anno);

		#	undef %anno;
		#undef %links;
		#}

	print ">>BED: @bed\n" if $debug ;
	print ">>MAP: @map\nvvvvvvvvvvvv\n" if $debug ;
	#VVVVVVVVVVVVVVVVVVVVVVV
	## from map	
	my $read = $map[3];
	my $as = $map[4];
	my $path = $map[5];
	my $mapn = $map[6];
	my $chkn = $map[7];

	### node hash
	# can make a univeral graph, reduce memory . 


	my ($nhr) = nodehash($path, $map[1], $map[2], "$read##$mapn", \%links);

	#VVVVVVVVVVVVVVVVVVVVVVVV
	# from anno	
	#######################
	my $gene = $bed[3];	
	$r2g{"$read##$mapn"}{AS} = $as;
	# no overlap, then put in into %ouside , using for extension  <<<<<<< type 3	
	#if($gene eq "."){
	#		next;
	#}
	## extract annos , save involved reads
	my @nodes = split /;/, $bed[-1];
	my $gl = $bed[-3];

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
		$anno{$id}{$gene}{D} = $d;
		$anno{$id}{$gene}{L} = $gl;
		$anno{$id}{$gene}{S} = $s;
		$anno{$id}{$gene}{E} = $e;
	}
}

########################
# subfunctuinos  #
sub check_overlap{
	# even sharing one node          , will be count as overlap. 
	my ($links, $anno) = @_;
	print "checking overlap \n" if $debug;	
	unless ( keys %{ $anno }){
		print "no anno... skip\n\n";
		return;
	}
	
	my %olen;	
	my %gc;


	foreach my $id ( sort { abs($a) <=> abs($b) } keys %{$links} ){
		my @rs = @{$$links{$id}{R}};
		my @pns = keys %{$$links{$id}{P}};
		my @nns = keys %{$$links{$id}{N}};
		my $idc = scalar @rs;
		print "ID $id C  $idc = @rs\n";
		print "ID $id P @pns\n";
		print "ID $id N @nns\n";
		
		my $ida = abs($id);
		if( exists $$anno{$ida} ){
			foreach my $g ( keys %{$$anno{$ida}} ){
				$gc{T}{$g} = 1;	
				if( $$anno{$ida}{$g}{D} * $SD * $id >= 0){ ## test anno node direction
					$gc{M}{$g} = 1;
					print "\t$g overlap $ida $$anno{$ida}{$g}{D} * $SD * $id\n";
					my $count = 0; ## count of overlapped positions
					my $total = 0; ## total coverage, that is depth
				
					my $nl = @{$$links{$id}{A}};
					print "range $nl $ida = $$anno{$ida}{$g}{S}.. $nl - $$anno{$ida}{$g}{E} - 1\n";
					print "nodearray @{$$links{$id}{A}}\n";
					
					my @rix;
					foreach my $ai ( $anno{$ida}{$g}{S}.. $nl - $$anno{$ida}{$g}{E} - 1 ){
						my $v = $$links{$id}{A}[$ai];
						my $vl = @{$v};
						print "INDEX $ai lenth $vl  @$v\n";
						if($vl > 0 ){
							$count ++;
						}
						$total += $vl;
						push @rix, @{$v};
					}
					@rix = uniq(@rix);

					$olen{$g}{$id} = [$count,$total,\@rix];	
					print "\t\tMATCH: count $count total $total\n";

				}else{
					print "\t$g overlap, but in reverse strand...\n" if $debug ;
				}
			}
		}else{
			print "\tno annotation\n" if $debug;
		}
	}

	if($debug){
		print "Total gene : ", keys %{$gc{T}}, "\n";
		print "Match gene : ", keys %{$gc{M}}, "\n";
	}

	print "findoverlaps:\n";
	foreach my $g (keys %olen){
		foreach my $id (keys %{$olen{$g}}){
			print "$g $id @{$olen{$g}{$id}} @{$links{$id}{R}}[@{$olen{$g}{$id}[2]}]\n";
		}
	}
	print "checking overlap finished...\n\n";
}

###############################
###############################
sub nodehash{
	#my($path,$ps,$pe,$plen,$rlen,$read,$mapn,$which) = @_;
	my ($path, $sn, $en,$read,$links) = @_;
	print "\t>>NodeHash: $path\n" if $debug;

	my @nodes = split /;/, $path;
	#my @nodes = $path =~ /(<|>)(\d+)/g ;
	
	my @units;
	for( my $i = 0; $i < @nodes; $i ++ ){

		my($id,$l,$s,$e,$t) = $nodes[$i]  =~ /(.\d+):(\d+),(\d+),(\d+),(\d+)/;
		#######################
		###### find links #####
		my($pn,$nn);
		if($i > 0){
			($pn) = $nodes[$i - 1] =~ /^(.\d+)/;
		}else{
			$pn = 0;
		}
		if($i < $#nodes){
			($nn) = $nodes[$i + 1] =~ /^(.\d+)/;
		}else{
			$nn = 0;
		}
		$$links{$id}{P}{$pn} ++;
		$$links{$id}{N}{$nn} ++;
		push @{$$links{$id}{R}}, $read;
		my $ri = $#{$$links{$id}{R}} ; 	
		print "read idnex $ri for node $id\n";
		if(exists $$links{$id}{A}){
			#$_++ for @{$$links{$id}{A}}[$s..$l-$e-1];
			push @{$_}, $ri  for @{$$links{$id}{A}}[$s..$l-$e-1];
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
		$min = abs($val) if !defined($min) || abs($val) < $min;
		$max = abs($val) if !defined($max) || abs($val) > $max;
	}
	return ($min, $max);
}
