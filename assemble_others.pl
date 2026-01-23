#!/usr/bin/perl
use warnings; use strict;
use Bit::Vector;
use Bio::DB::HTS::Tabix;


my($other_file , $g2l_file ) = @ARGV;

my $debug = 1;
my $step = 10;
#########
my $tabix_g = Bio::DB::HTS::Tabix -> new(filename => $g2l_file);
my $tabix_o = Bio::DB::HTS::Tabix -> new(filename => $other_file );

open OR, "zcat $other_file | " or die $!;

my $start_node;
my $end_node;

my %bnodes;
my @lines;
my %seen;
my @regions;

while(<OR>){
	chomp;
	my @ar = split /\t/;
	
	next if($seen{"$ar[1]+$ar[2]"} );
	
	if($end_node and $ar[1] >= $end_node) {
		# a break, then move look around for branches

		### Block : branche looks up  <<<<<
		push @regions, [$start_node, $end_node];
		print "\nstart region: $start_node = $end_node\n";
		look_branches(\@regions, \%bnodes, \@lines);
		process_lines(\@lines);
		
		undef @lines;
		undef %bnodes;
		undef $start_node;
		undef $end_node;
	}
	
	## collect branch/nodes seen in reads
	if( $ar[-2]  > 1){
		my $path = $ar[5];
		for my $item (split /;/, $path) {
			my $k = abs( (split /,/, $item)[0] );
			if($k < $ar[1] + 1 or $k > $ar[2]){
				$bnodes{$k} = 1;
			}
		}
	}
	
	$end_node = $ar[2];
	$start_node = $ar[1] + 1 unless defined $start_node;

}


sub process_lines {
	print "--> Enter process lines module...\n";
	my ($lines) = @_;
	foreach my $l (@$lines){
		print "\tPP: $l\n";
	}
}
=head
		unless($Vnode{$nid}){
			my $vec = Bit::Vector -> new( $nid );
			$Vnode{$nid}{V} = $vec;
		}
		$Vnode{$nid}{V} -> Interval_Fill($head, $nlen - $tail - 1);
		
		#####
		if($lnid){
			$Vnode{$lnid}{N}{$nid}  ++;
		}
		my $rinf = [$ar[6],$ar[3],$ar[4], $ar[7]);
		$Vnode{$nid}{R}{$head}{ = $rinf;

		push @{$Reads{$read}} ,  $rinf;

		$lnid = $nid;
	}

	if($last_node and $ar[1] >= $last_node){
		# a possible break points
		### output....

		my @pnodes = grep { $_ <= $ar[1] } keys %Vnode;
		
		foreach my $pn (@pnodes){

	}

}
=cut 

sub look_branches {
	print "--> Enter look around model...\n" if $debug ;
	my($regions, $bnodes,$lines) = @_;
	
	my($start_node, $end_node) ; 
	my $loopc = 0;
	while(1){
		#### real branch nodes
		my @bns = sort {$a <=> $b } keys %$bnodes ;
		if(@bns){
			print "\tstart finding branch @bns\n";
		}else{
			print "no branches..\n";
			return;
		}
		
		### test if node have seen or not
		my $head = shift @bns ;
		delete $$bnodes{$head};

		my $skip = 0;
		foreach my $r (@$regions){
			if($head  >= $$r[0] and $head  <= $$r[1]){
				$skip = 1;
				last;
			}
		}
		if($skip){
			print "\t--$head skip..\n";
			next;
		}

		###
		$loopc ++;

		print "checking branch $head. $#bns + 1 left\n";	
		if(! defined $start_node){
			($start_node, $end_node) = ($head,$head);
		}elsif ($head >= $start_node and $head <= $end_node and $end_node > $start_node ){
			print "\tskip encouted nodes.. $head \n";
			next;
		}

		#####################
		##### backwards #####
		#####################
		my $uloop = 1;
		
		print "backwards...$head\n" if $debug ;
		while(1){
			my $ls = $head - $uloop * $step + 1;
			my $le = $head - ($uloop - 1) * $step;
			if($le <= 0){
				last;
			}
			$ls  = $ls < 0 ? 0 : $ls;

			print "\tB $ls -> $le\n";
			
			my $ui  = $tabix_o -> query("G:$ls-$le");
			my @ul; 
			while( my $l = $ui -> next ){
				push @ul, $l;
			}
			print "\t===backward $uloop lines $#ul : $ls - $le\n" if $debug ;

			my $break = 1;
			foreach my $l ( reverse @ul ){
				print "\tback branch:$l\n";
				my  @ar = split /\t/, $l;
				if( $ar[2] >= $start_node ){
					push @$lines, $l;
					$start_node = $ar[1] + 1;
					$break = 0;
				}else{
					$break = 1;
					last;
				}

				if( $ar[-2]  > 1){
					my $path = $ar[5];
					for my $item (split /;/, $path) {
						my $k = abs( (split /,/, $item)[0] );
						if($k < $ar[1] + 1 or $k > $ar[2]){
							$$bnodes{$k} = 1;
						}
					}
				}
			}
			if($break ){
				last;
			}
			$uloop ++;
		}
		#####################
		#####  Forward  #####
		#####################

		my $floop = 1;
		print "forward...$head\n" if $debug;
		while(1){
			my $le = $head + $floop * $step  - 1;
			my $ls = $head + ($floop - 1) * $step ;
			print "\tF $ls -> $le\n";
			my $fi  = $tabix_o  -> query("G:$ls-$le");
			my @fl ;
			while( my $l = $fi -> next ){
				push @fl, $l;
			}
			print "\t===forward $floop lines $#fl : $ls - $le\n";	
			my $break = 1;
			foreach my $l ( @fl ){
				print "\tforw branch:$l\n";
				my @ar = split /\t/, $l;
				if($ar[1]+1 <= $end_node ){
					push @$lines, $l;
					$end_node = $ar[2];
					$break = 0;
				}else{
					$break = 1;
					last;
				}
				if( $ar[-2]  > 1){
					my $path = $ar[5];
					for my $item (split /;/, $path) {
						my $k = abs( (split /,/, $item)[0] );
						if($k < $ar[1] + 1 or $k > $ar[2]){
							$$bnodes{$k} = 1;
						}
					}
				}
			}
			if( $break ){
				last;
			}
			$floop ++;
		}
		push @$regions, [$start_node,$end_node];
		undef $start_node;
		undef $end_node;
	}
	print "\tread end\n";
	#look_branches($regions, $bnodes,$lines) ; 
	
}


sub nodes_ranges {

    my $debug = 0;
    my @nums = @_;

    @nums = sort { $a <=> $b } @nums;

    my %ranges;
    my %includs;
    $includs{$nums[0]} = 1;

    my $start = $nums[0];
    my $end   = $nums[0];

    for my $n (@nums[1 .. $#nums]) {
        if ($n - $end < 10) {
            # extend the current range
            $includs{$n} = 1;
            $end = $n;
        } else {
            # save previous range
            my $rg = "G:".($start-1)."-$end";
            my %inc = %includs;
            $ranges{$rg} = \%inc;
            my @ins = keys %includs;
            # start a new range
            $start = $n;
            $end   = $n;
            %includs = ($n , 1);
        }
    }

    # push last range
    my $rg = "G:". ($start-1) . "-$end";
    $ranges{$rg} = \%includs;

    return (%ranges);

}








