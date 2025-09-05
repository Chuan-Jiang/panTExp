#!/usr/bin/perl
use warnings; use strict;
use List::Util qw(min max sum);

### 
my %rds_map; ####
my %as;  ####
my %feas;

while(<>){
	chomp;
	my @ar = split /\t/;
	# reads  feature  S
	push @{$rds_map{$ar[1]}{F}{$ar[0]}{Q}}, $ar[2];
	push  @{$as{$ar[1]}}, $ar[2];
	$feas{$ar[0]}{P} = 1;    
	$feas{$ar[0]}{T} = 1;
	if(0){
		unless ($rds_map{$ar[1]}{F}{$ar[0]}{Q} ){
			$rds_map{$ar[1]}{F}{$ar[0]}{Q} = $ar[2];
			push @{$as{$ar[1]}}, $ar[2];
		}elsif( $ar[2] > $rds_map{$ar[1]}{F}{$ar[0]}{Q} ){
			$rds_map{$ar[1]}{F}{$ar[0]}{Q} = $ar[2];
			push @{$as{$ar[1]}}, $ar[2];
		}
	}
}

###############
# initiation
# #############
my $pi_prior = 0;
my $theta_prior = 2000;
my $k = scalar(keys %feas);
my $rk = scalar (keys %rds_map );
my $pi = 1 / $k;
my $theta = $pi;
my $pi_t = 1; 

my $total_weight;    # similar of the count of reads, 
my $total_weight_a;

###
# need to do:  lenth should be considered, curretly set 300 as default values.
foreach my $r (keys %rds_map){
	my ($softr,$max_q,$exp_total) = softmax( $as{$r} );
	my @features = keys %{$rds_map{$r}{F}}; 
	$rds_map{$r}{M} = $max_q;
	$total_weight += $max_q;
	$total_weight_a += $max_q if (@features > 1);
	
	#print "Read $r\n";
	foreach my $f ( @features ){
		#print "\tFeature $f $rds_map{$r}{F}{$f}{Q}\n";
		my $multi = @features == 1 ? 0 : 1; 
		$rds_map{$r}{F}{$f}{Y} =  $multi; 
		if($multi){
			@{ $rds_map{$r}{F}{$f}{Q} } = @$softr{ @{$rds_map{$r}{F}{$f}{Q}} } ;
		}else{
			@{ $rds_map{$r}{F}{$f}{Q} } = max( @$softr{ @{$rds_map{$r}{F}{$f}{Q}} } );  
		}
		#print "\t\t $r $f = $rds_map{$r}{F}{$f}{Q} $rds_map{$r}{F}{$f}{Y}\n";
	}
}
#print "prepare reading finished\n";

my $loop = 0;

### EM loops
while($loop < 100){
	##
	$loop ++;
	my $rd_tr  = estep();
	my $terminate = mstep($rd_tr);
	if($terminate){
		last;
	}
}
foreach my $f (keys %feas){
	my $ct = $feas{$f}{P}/$k * $rk;
	$ct = sprintf("%.0f", $ct);
	print "$f\t$ct\n";
}
#############################
#############################
sub mstep {
	my ($rd_tr) = shift @_;
	
	#print "MSTEP\n";
	my %feas_old = %feas;
	undef %feas; 
	foreach my $r (keys %rds_map){	
		foreach my $f ( keys %{$rds_map{$r}{F}} ){
			$feas{$f}{P} += ( $_ / $$rd_tr{$r} * $rds_map{$r}{M} + $pi_prior)/($total_weight + $pi_prior * $k) * $k  for @{$rds_map{$r}{F}{$f}{Z}} ;
			if($rds_map{$r}{F}{$f}{Y} ){
				$feas{$f}{T} += ( $_ / $$rd_tr{$r} * $rds_map{$r}{M}  + $theta_prior )/ ($total_weight_a + $theta_prior * $k) * $k for @{$rds_map{$r}{F}{$f}{Z}} ;
			}
		}
	}
	my $dif;
	my $total;
	foreach my $f ( keys %feas ){
		#warn "M:P\t$feas{$f}{P}\t$feas_old{$f}{P}\t$f\n";
		if(exists $feas{$f}{T} ){
			#warn "M:T\t$feas{$f}{T}\t$feas_old{$f}{T}\t$f\n";
		}
		$total +=  $feas{$f}{P};
		$dif += abs($feas{$f}{P} - $feas_old{$f}{P});
	}
	warn  "DIFF $dif $total $k \n";
	if($dif < 1e-6) {
		return 1;
	}else{
		return 0;
	}
}

## based on mapping quality and feature ratios for each fragment/reads, estimate total amount of features(%fea_t) and amount of a reads(from 0 to 1) to a specific freature : Z in mapr . 
sub estep {
	my %rd_t;
	#print "ESTEP start..\n";
	foreach my $r ( keys %rds_map ){
		foreach my $f (keys %{$rds_map{$r}{F}} ){
			my @z; 
			if($rds_map{$r}{F}{$f}{Y}){
				@z = map { $feas{$f}{P} * $feas{$f}{T} * $_ } @{$rds_map{$r}{F}{$f}{Q}} ;
			}else{
				@z = map { $feas{$f}{P} * $_ } @{$rds_map{$r}{F}{$f}{Q}} ;
			}
			#print "E:$r\t$f\t$z\n";
			$rds_map{$r}{F}{$f}{Z} = \@z;
			$rd_t{$r} += $_ for @z;
		}	
		#print "E:$r\ttotal\t$rd_t{$r}\n\n";
	}
	#print "ESTEP finish..\n";
	return \%rd_t;
}

######softmax function

sub softmax {
	my ($nums) = @_;
	my %nums_out;
	my $max_q;

	my $max_num = max(@$nums);
	#print "max $max_num\n";	
	my $exp_total = 0;
	for (my $i = 0; $i <@$nums; $i ++) {
		$exp_total += exp( $nums->[$i] /  $max_num  * 100 ) - 1;
	}	
 
	for (my $i = 0; $i <@$nums; $i ++) {
		
		my $r = ( exp( $nums->[$i] / $max_num * 100  ) - 1 ) ;
		my $z= ( exp( $nums->[$i] / $max_num * 100  ) - 1 ) / $exp_total;
		$nums_out{$nums->[$i]} = $r;
		if($nums -> [$i] == $max_num){
			$max_q = $r;
		}
	}
	return (\%nums_out,$max_q,$exp_total);
}

