#!/usr/bin/perl
use warnings; use strict;

my ($file_in, $pre, $maxd) = @ARGV;

my %ndeps;
my @lines;
my $lidx = 0;
my $node_last ;
my $filter_rd_iden = 0.8;
my %deadnode;  # this is used to filter nodes that have very high coverage
my %deadnode_p;

(my $depth_filter,$maxd)  = $maxd  ? (1, $maxd) : (0,0) ;

warn "depthfilter $depth_filter maxd $maxd\n";

my $debug = 0;
my $all_dead = 0;

my $fh;
if ($file_in =~ /\.gz$/) {
    open($fh, "-|", "gzip -dc $file_in") or die "Cannot open $file_in: $!";
} else {
    open($fh, "<", $file_in) or die "Cannot open $file_in: $!";
}

open DISC, ">$pre.discarded.tsv" or die $!;
open FILT, ">$pre.filtered.tsv" or die $!;

my $next_line  = <$fh>;
while($next_line){
	my $line = $next_line;
	$next_line  = <$fh>;
	chomp($line);
	
	print "NEW: $line\n" if $debug;
	my @ar  = split /\t/, $line ;
	my @inf1 = split /\|\|/, $ar[8];
	my @inf2 = split /\|\|/, $ar[9]; 
	
	my $csr1 = cs($inf1[2])/$inf1[0] ;
	my $csr2 = cs($inf2[2])/$inf2[0] ;

	# quality contriol
	my @dtag ;
	if($ar[6] eq "*"){
		push @dtag  , "unmap1";
	}
	if($ar[7] eq "*"){ 
		push @dtag  , "unmap2";
	}
	if($csr1  < $filter_rd_iden){
		push @dtag , "cs1=$csr1"
	}
	if($csr2  < $filter_rd_iden ){
		push @dtag , "cs2=$csr2"
	}
	if(@dtag){
		print DISC "$line\tdisc:s:", join(",", @dtag) , "\n";
		print "discarded for @dtag skip \n" if $debug;
		next;
	}
    
	# V collect head tail and length to %ninfs from a specific read. if all are dead nodes, then will get empty hash 
	# which indicated discard this line.

	my %ninfs;
	my $encounter = 0;
	foreach my $p ( $ar[6], $ar[7]){
		foreach my $ni  ( split /;/, $p){
			my @na = split /,/, $ni;
			my $n = abs($na[0]);
			if( exists $deadnode{$n} or exists $deadnode_p{$n} ){
				$encounter ++;
				next;
			}
			$ninfs{$n}{L} = $na[1];
			if($ninfs{$n}{H}){
				$ninfs{$n}{H} = $ninfs{$n}{H} < $na[2] ? $ninfs{$n}{H} : $na[2];
			}else{
				$ninfs{$n}{H} = $na[2];
			}

			if( $ninfs{ $n }{T} ){
				$ninfs{ $n }{T} = $ninfs{ $n }{T} < $na[3] ? $ninfs{ $n }{T} : $na[3];
			}else{
				$ninfs{ $n }{T} = $na[3];
			}
		}
	}

	unless(%ninfs){
		print DISC "$line\tdisc:s:HiDup\n";
		print "all deadnode found skiippp\n" if $debug;
		$all_dead = 1;
		next;
	}else{
		my @rns = keys %ninfs;
		print "remaine nodes: @rns\n" if $debug;
	}

	
	### decide if it is a break and output previous	
	if( ($node_last and $node_last <= $ar[1] ) or ( ! $encounter and $all_dead) or !$next_line)  { # the second condition is for all previous lines are all dead nodes
		
		if($node_last > $ar[1]){	
			my @dn = keys %deadnode;
			print " exit of deadnode $#dn: @dn\n" if $debug ;
		}
		
		my %disc;
		print "have a break, output lines\n" if $debug ;
		my @dns ;
		if($depth_filter){
			print "depth filter on\n" if $debug ;
			foreach my $n ( keys %ndeps ){
				if($deadnode_p{$n} or $deadnode{$n}){
					print "collect index on deadnode $n \n" if $debug; 		
					push @dns, $n; 
					$disc{$_} = 1 for  map{ @$_  } @{$ndeps{$n}{A}};
					next;
				}

				if($ndeps{$n}{C} < $maxd ){
					next;
				}
				for(my $i = 0; $i <  @{$ndeps{$n}{A}}; $i ++){
					my @p =  @{$ndeps{$n}{A}[$i]};
					my $rc = scalar @p;
					if($rc  > $maxd ){
						#print "$rc reads removed at $n $i >>> @p  \n";
						foreach my $li (@p){
							$disc{$li} = 1;
						}
					}
				}
			}
		}
		for(my $i = 0; $i < @lines; $i ++){
			unless($disc{$i}){
				print FILT "$lines[$i]\n" ;
			}else{
				print DISC "$lines[$i]\tdisc:s:HiDup\n";
			}
		}

		undef %ndeps;
		undef @lines;

		%deadnode_p = %deadnode;
		$deadnode_p{$_} = 1 for @dns ;
		undef %deadnode;
		$all_dead = 0;
		$lidx = 0;
	}else{
		print "continue to collect\n" if $debug ;
	}

	
	##################################
	### colect read coverage 	######
	##################################
	
	# V assign coverage  # pao hui 
	push @lines, $line;
	$node_last =  $ar[2];
	
	if($depth_filter){
		foreach my $n ( keys %ninfs ){
			print "pile node $n\n" if $debug ;
			my $h =  $ninfs{$n}{H};
			my $t =  $ninfs{$n}{T};
			my $l =  $ninfs{$n}{L};
			
			unless($ndeps{$n} ){
				my @par = map { [] } 1..$l;
				$ndeps{$n}{A} = \@par;
				$ndeps{$n}{c} = 0;
			}

			$ndeps{$n}{C} ++;
			
			foreach my $p ($h .. $l-$t-1){
				#print "save $n $p <- $lidx to  @{$ndeps{$n}{A}[$p]} \n";
				
				if( @{$ndeps{$n}{A}[$p]} <= $maxd){
					push @{$ndeps{$n}{A}[$p]}, $lidx;
				}else{
					print "found a positons $n $p full ..\n" if $debug ;
					unless($ndeps{$n}{A}[$p][-1] eq "TOP"){
						$ndeps{$n}{c} ++;
						print "full rate $ndeps{$n}{c} $l\n" if $debug ;
						push  @{$ndeps{$n}{A}[$p]}, "TOP";
					}
				}
			}
			if ($ndeps{$n}{c} == $l){
				print "create a dead node $n\n" if $debug ;
				$deadnode{$n} = 1;
			}else{
				print "node $n not dead yet $ndeps{$n}{c} < $l \n" if $debug ;
			}
		}
		$lidx ++;
	}
}

sub cs {
    my ($str) = shift @_;
    my $m = 0;
    while($str =~ /:(\d+)/g ){
        $m += $1;
    }
    return $m;
}
