#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Storable;
use Parallel::ForkManager;
use List::MoreUtils qw(upper_bound lower_bound);
use List::Util qw(min max uniq sum);
use feature 'say';
use Bio::DB::HTS::Tabix;
use Graph::Directed;
use File::Temp qw(tempfile);


use Time::HiRes qw(gettimeofday tv_interval);

my $dis_f;
my $gaf_input;
my $pre;
my $debug;
my $num_processes = 1; # Number of parallel processes
my $step = 20;
my $nodelen_f;
my $nodenei_f;
my $CST = 0.8;
my $lib = 1000;
my $line_bunch = 100000;

GetOptions(
    "i|input=s"  => \$gaf_input,
    "p|prefix=s" => \$pre,
	"b|neib=s" => \$nodenei_f,
    "s|step=i"  => \$step,
	"l|nodelen=s"    => \$nodelen_f,
    "g|bug"     => \$debug,
	"n|num_processes=i" => \$num_processes # Allow setting number of processes
) or die "Invalid options!\n";

#######
my $pm = Parallel::ForkManager->new($num_processes);
warn "paralell $num_processes\n";

open BED, " > $pre.gaf2.bed" or die $!;
open DIS, " > $pre.discard.gaf" or die $!;
#### gaf 
my $gaf_fh;
if ($gaf_input) {
    if ($gaf_input eq "-") {
        $gaf_fh = *STDIN;
    } else {
        open my $fh, "zcat $gaf_input |" or die $!;
        $gaf_fh = $fh;
    }
} else {
    $gaf_fh = *STDIN;
}

################
## create temp files

my @tempfiles;
my %fh_assign;
for my $i (0 .. $num_processes ) {
    my ($fhb, $filenameb) = tempfile("childbed_$i-XXXX", UNLINK => 1, DIR => ".");
    my ($fhg, $filenameg) = tempfile("childgaf_$i-XXXX", UNLINK => 1, DIR => ".");
    push @tempfiles, { fhb => $fhb, fileb => $filenameb,fhg => $fhg, fileg => $filenameg };
}

##########
print "geting node len...\n" if $debug;
my $nodelen = retrieve($nodelen_f); 

print "got node len...\n" if $debug;

### 
$pm->run_on_finish(sub {
    my ($pid, $exit_code, $ident, $exit_signal, $core_dump) = @_;
	 if ($exit_signal) {
        print "Child $pid died from signal $exit_signal\n";
    }
	print "$ident released by $pid..\n";
	delete $fh_assign{$ident};
});
##############
my @lines;
my $lastr = 0;
my $lc = 0;
while (<$gaf_fh>) {
	$lc ++;
    chomp;
    my @ar = split /\t/; 

	if ( $lc >= $line_bunch and $ar[0] ne $lastr ) { # Process in batches of 1000 lines
        print "processing line\n" if $debug;
		
		process_lines([@lines]);
		$lc = 0;
        @lines = (); # Clear the buffer
    }
	push @lines, $_;
	$lastr = $ar[0];
}

process_lines([@lines]) if @lines; # Process remaining lines

$pm->wait_all_children; # Ensure all processes complete

foreach my $tf (@tempfiles) {
    open my $inb, '<', $tf->{fileb} or die "Can't open $tf->{fileb}: $!";
    print "--- Contents of $tf->{fileb} ---\n";
    while (<$inb>) {
        print BED $_;
    }
    close $inb;
    unlink $tf->{fileb};  # Clean up

	open my $ing, '<', $tf->{fileg} or die "Can't open $tf->{fileg}: $!";
	print "--- Contents of $tf->{fileg} ---\n";
	while (<$ing>) {
		print DIS  $_;
	}
	close $ing;
	unlink $tf->{fileg};  # Clean up
}



sub process_lines {
    my ($lr) = @_;	

	my $fhi;
	#### choose file handle
	for(my $i = 0 ;$i < @tempfiles; $i ++){
		unless (  defined $fh_assign{$i} ){
			$fhi = $i;
			last;
		}
	}

	my $pid = $pm->start($fhi);
	if($pid){
		print "child $pid ... FHI:$fhi\n";
		$fh_assign{$fhi} = 1;
		return;
	}

	my $tabix = Bio::DB::HTS::Tabix->new(filename => $nodenei_f);  
	my $fhb = $tempfiles[$fhi] -> {fhb};
	my $fhg = $tempfiles[$fhi] -> {fhg};
	my @lines = @{$lr};

	my %tout; # out put for the total lines
	my $batch = scalar @lines;
	print "BATCH $batch\n";
	### get node leng
	my %multimap;
	for (my $i = 0; $i < @lines; $i += 2 ){	
		my @ar1 =  split /\t/, $lines[$i]; # @{ $lines[$i] };
		my @ar2 =  split /\t/, $lines[$i + 1]; # @{ $lines[$i + 1] };
		print "ARRAY1:@ar1\n" if $debug ;
		print "ARRAY2:@ar2\n" if $debug ;
	
		my %tag1;
		my %tag2;
		foreach my $t (@ar1[12..$#ar1]){
			my ($k,$v) = $t =~ /([^:]+):\w:(.+)/;
			$tag1{$k} = $v;
		}
		foreach my $t (@ar2[12..$#ar2]){
			my ($k,$v) = $t =~ /([^:]+):\w:(.+)/;
			$tag2{$k} = $v;
		}


		#  which mapping?
		if($multimap{$ar1[0]}){
			$multimap{$ar1[0]} ++;
		}else{
			$multimap{$ar1[0]} = 1;
		}
		# test mapped or not 
		my $discard  = "";
		if ( $ar1[5] eq "*" ){
			$discard .= "(1:U)";
		}
		if( $ar2[5] eq "*" ){
			$discard .= "(2:U)";
		}
		if($discard){	
			#push @{$tout{DIS}}, join "\t", @ar1,"rm:z:$discard";
			#push @{$tout{DIS}}, join "\t", @ar2,"rm:z:$discard";
			say $fhg join "\t", @ar1, "rm:z:$discard";
			say $fhg join "\t", @ar2, "rm:z:$discard";
			next;
		}
		# test short target
		if($ar1[6] / $ar1[1] < 0.9){
			$discard .= "(1:T)";
		}
		if($ar2[6] / $ar2[1] < 0.9){
			$discard .= "(2:T)";
		}
		if($discard ){
			say $fhg join "\t", @ar1, "rm:z:$discard";
			say $fhg join "\t", @ar2, "rm:z:$discard";
			next;
		}

		### tesing diversity
		my $m1 = 0;
		while($tag1{cs} =~ /:(\d+)/g ){
			$m1 += $1;
		}
		my $m2 = 0;
		while($tag2{cs} =~ /:(\d+)/g ){
			$m2 += $1;
		}
		if( $m1 / $ar1[1] < $CST ){
			$discard .= "(1:D)";
		}
		if( $m2 / $ar2[1] < $CST ){
			$discard .= "(2:D)";
		}
		if($discard){
			#push @{$tout{DIS}}, join "\t", @ar1,"rm:z:$discard";
			#push @{$tout{DIS}}, join "\t", @ar2,"rm:z:$discard";
			#if($debug){
				say $fhg join "\t", @ar1, "rm:z:$discard";
				say $fhg join "\t", @ar2, "rm:z:$discard";
			#}
			next;
		}

		#################################
		###### test if properly mapped ###
		my @allnodes;

		## nodes on reads
		my @nodes1;
		while($ar1[5] =~ /(<|>)(\d+)/g){
			push @nodes1, $1 eq ">" ? $2 : -$2;
		}
		my @nodes2;
		while($ar2[5] =~ /(<|>)(\d+)/g){
			push @nodes2, $1 eq ">" ? -$2: $2;
		}
		@nodes2 = reverse @nodes2;
		
		print "N1 @nodes1\n" if $debug;
		print "N2 @nodes2\n" if $debug ;
		my $at1 = $ar1[8] eq "*" ? 0 : $ar1[8];
		my $at2 = $ar2[8] eq "*" ? 0 : $ar2[8];
		my $r1_lft = $ar1[6] - $at1 + 1;
		my $r2_lft = $ar2[6] - $at2 + 1;
	
		###  COMPARE nodes   ####	
		my $pos = 0;
		for (my $i = 1; $i <= @nodes1 && $i <= @nodes2; $i++) {
			my @tail =  @nodes1[-$i..-1];
			my @head =  @nodes2[0..$i-1];
			print "@tail \\\ @head\n" if $debug ;
			if ("@tail" eq "@head") {
				$pos = $i;
			}elsif($pos){
				last;
			}
		}

		if($debug){
			print "OVERLAP  POS $pos \n"  ;
		}
		## have node overlapesssss
		if($pos ){
			if( $pos  == 1){
				my $nl = $$nodelen{ abs($nodes2[0]) } ? $$nodelen{ abs($nodes2[0]) } : 1;
				my $nodeleft = $r1_lft + $r2_lft - $nl ;
				if($nodeleft +  $ar1[1] + $ar2[1] > $lib ){
					$discard .= "(Fl)";
				}elsif($nodeleft + $ar1[1] < 0 and $nodeleft + $ar2[1] < 0){
					$discard .= "(Fs)";
				}
			}else{
				### need finish when using real dataset
			}
			if($discard){
				#push @{$tout{DIS}}, join "\t", @ar1,"rm:z:$discard";
				#push @{$tout{DIS}}, join "\t", @ar2,"rm:z:$discard";
				#if($debug){
					say $fhg join "\t", @ar1, "rm:z:$discard";
					say $fhg join "\t", @ar2, "rm:z:$discard";
				#}
				next;
			}
			my @mn = $pos == 1? $nodes2[0] : @nodes2[0..$pos-1];
			print "comb:@nodes1,,,OVER,,, @nodes2\n" if $debug;
			@allnodes = (\@nodes1, [], \@nodes2);
		##  don's have overlap. need to estimate distance  
		}else{
			## do we need to go linking nodes
			my $left = $lib - $ar1[1] - $ar2[1] - $r1_lft - $r2_lft;
			if($left < 0 ){
				print "left $left. negtive ..\n" if $debug;
				$discard .= "(Fl)";
			}else{
				my $start_time = [gettimeofday]; 
				my ($pathl,$pathr) = connect_node ($nodes1[-1] , $nodes2[0], $tabix, $left );	
				my $elapsed = tv_interval($start_time); 
				print "Execution time: $elapsed seconds\n" if $debug ;
				
				if($pathl < 0 or $pathl > $left){
					print "left $left short $pathl. negtive ..\n" if $debug;
					$discard .= "(Fl)";
				}else{
					my @mn = @{$pathr};
					print "comb:@nodes1,,,MN @mn,,, @nodes2\n" if $debug;
					@allnodes = (\@nodes1, $pathr, \@nodes2);
				}
			}
			if($discard){
				#if($debug){
					say $fhg join "\t", @ar1, "rm:z:$discard";
					say $fhg join "\t", @ar2, "rm:z:$discard";
				#}
				#push @{$tout{DIS}}, join "\t", @ar1,"rm:z:$discard";
				#push @{$tout{DIS}}, join "\t", @ar2,"rm:z:$discard";
				next;
			}
		}
		
		my ($r1_h,$r1_t)  = ($ar1[7], $ar1[6] - $ar1[8]);
		my ($r2_h,$r2_t)  = ($ar2[6] - $ar2[8], $ar2[7]);
		
		my $as = $tag1{AS} + $tag2{AS}; 
		print "ALL:$as,$multimap{$ar1[0]}, @allnodes\n" if $debug ;
		out ( \@allnodes, $nodelen, $ar1[0], $r1_h,$r1_t, $r2_h,$r2_t,  $as, $multimap{$ar1[0]} , $fhb );
	} 
	print "subprocess, finished\n" if $debug;
	$pm->finish(0, \%tout); 
}


sub connect_node {
	my ($n1,$n2, $tabix, $left) = @_;
	
	print "connect: $n1 $n2\n" if $debug ;
	my $g1 = Graph::Directed->new;     # D

	my ($start,$end)  = sort {$a <=> $b} (abs($n1), abs($n2));	
	my $distance = $end - $start  + 1;
	my @wins;
	my $ext = 20;
	if($distance < 2 *  $ext ){
		my $ws = $start  - $ext < 0 ? 0 : $start  - $ext ;
		my $we = $end + $ext;
		@wins = "G:$ws-$we";
	}else{
		my $ws1 = $start - $ext < 0 ? 0 : $start - $ext;
		my $we1 = $start + $ext;
		my $ws2 = $end - $ext <0 ? 0 : $end - $ext;
		my $we2 = $end + $ext;
		@wins = ("G:$ws1-$we1", "G:$ws2-$we2") ;
	}
	print "WIN:@wins\n" if $debug;
	## put each line into a hash
	foreach my $win(@wins){
		my $iter = $tabix -> query($win);
		while (my $line = $iter->next) {
			#print "$line\n" if $debug ;  # VCF or other tabix-formatted line
			my @ar = split /\t/, $line;
			my $cn = $ar[2];
			while($ar[3] =~ /(-?\d+)/g){
				my $n = $1;
				$g1 -> add_weighted_edge(-$cn,$n,  $$nodelen{abs($n)} ? $$nodelen{abs($n)} : 1 );
			}
			while($ar[4] =~ /(-?\d+)/g){
				my $n = $1;
				$g1 -> add_weighted_edge($cn,$n,  $$nodelen{abs($n)} ? $$nodelen{abs($n)} : 1 );
			}		
		}
	}
	my $pathl = 0;
	my @path = $g1->SP_Dijkstra($n1, $n2);
	return  -1  unless @path;
	shift @path;
	pop @path;
	for(my $i = 0; $i < @path ; $i ++){
		my $p = $path[$i];
		$pathl += $$nodelen{abs($p)} ? $$nodelen{abs($p)} : 1
	}
	print "PP1:$n1 $n2 $pathl = @path\n" if $debug;
	return ($pathl, \@path);
}

sub form_node{
	my ($ref,$nodelen, $hr,$tr) = @_;
	my $path = "";
	my @ns = @{$ref};
	for (my $i = 0; $i < @ns; $i ++){	
		my($h,$t) = (0,0) ;
		if($i == 0){
			$h  = $hr;
		}
		if($i == scalar @ns  - 1){
			$t = $tr;
		}
		if($ns[$i] < 0){
			($h,$t) = ($t,$h);
		}else{
			$ns[$i] = "+$ns[$i]";
		}
		my $l = $$nodelen{abs($ns[$i])} ? $$nodelen{abs($ns[$i])} : 1;
		$path = $path? "$path;$ns[$i]:$l,$h,$t" : "$ns[$i]:$l,$h,$t";
	}
	return $path;
}

sub out {
    my ($nodes,$nodelen,$read,$r1_h,$r1_t, $r2_h,$r2_t,$as, $multi, $fhb )  = @_;
	
	my $pathr1 = form_node($$nodes[0],$nodelen,$r1_h,$r1_t);
	my $pathr2 = form_node($$nodes[2],$nodelen,$r2_h,$r2_t);
	my $pathm =  form_node($$nodes[1],$nodelen,0,0) ;

	my $anno = "$read\t$as\t$pathr1|$pathm|$pathr2";	
	my @pn = sort { abs($a) <=> abs($b) } (@{$$nodes[0]}, @{$$nodes[1]}, @{$$nodes[2]}) ;
	@pn = uniq(@pn);	
	my @chunks;
	my $startn;
	for (my $k = 0; $k < @pn; $k ++){
		my $n = abs($pn[$k]);
		unless($startn){
			$startn = $n;
			next;
		}
		my $n_l = abs($pn[$k - 1]);
		if(abs($n - $n_l) < $step){
			next;
		}else{
			push @chunks, "$startn\t$n_l";
			$startn = $n;
		}
	}
	if($startn){
		my $ln = abs($pn[-1]);
		push @chunks, "$startn\t$ln";
	}
	print "\tCC:@chunks \n" if $debug ;
	
	#push @{$beds{R}}, $read;
	my $chunk_count = scalar @chunks;
	print "\t$chunk_count = @chunks\n" if $debug ;	
	for(my $j = 1; $j <= $chunk_count; $j ++){
		my $c = $chunks[$j-1];
		print "\t===>G\t$c\t$anno\t$multi\t$j/$chunk_count\n\n" if $debug;
		say $fhb "G\t$c\t$anno\t$multi\t$j/$chunk_count";
		#push @{ $$tout{BED} }, "G\t$c\t$anno\t$multi\t$j/$chunk_count";
	}
}

