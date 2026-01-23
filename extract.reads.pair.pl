#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Storable;
use Parallel::ForkManager;
use List::MoreUtils qw(upper_bound lower_bound);
use List::Util qw(min max uniq sum);
use feature 'say';
#use Graph::Directed;
use File::Temp qw(tempfile);
use File::Basename;
use File::Path qw(make_path);
use IO::Handle;

use Time::HiRes qw(gettimeofday tv_interval);

my $gaf_input;
my $pre;
my $debug;

my $num_processes = 1; # Number of parallel processes
my $CST = 0.8;
my $lib = 1000;
my $line_bunch = 500000;
my $nodelen_f;
my $score_diff = 10;

GetOptions(
    "i|input=s"  => \$gaf_input,
    "p|prefix=s" => \$pre,
    "g|bug"     => \$debug,
	"l|nodelen=s" => \$nodelen_f,
	"n|num_processes=i" => \$num_processes # Allow setting number of processes
) or die "Invalid options!\n";

#######
my $pm = Parallel::ForkManager->new($num_processes);
warn "paralell $num_processes\n";

# Extract directory part from prefix
my $dir = dirname($pre);

# If directory doesn't exist, create it
unless (-d $dir) {
    make_path($dir) or die "Failed to create directory $dir: $!";
}

my $bedout = "$pre.gaf2.bed" ;
#open BED, " > $pre.gaf2.bed" or die $!;

#### gaf 
my $gaf_fh;
if ($gaf_input) {
    if ($gaf_input eq "-") {
        $gaf_fh = *STDIN;
    } elsif( $gaf_input =~ /.gz/ ){
        open my $fh, "zcat $gaf_input |" or die $!;
        $gaf_fh = $fh;
    }else{
		open my $fh , " $gaf_input " or die $!;
		$gaf_fh = $fh;
	}
} else {
    $gaf_fh = *STDIN;
}

##########
print "geting node len...\n" if $debug;
my $nodelen = retrieve($nodelen_f);

print "got node len...\n" if $debug;


################
## create temp files

my @tempfiles;
my %fh_assign;
for my $i (0 .. $num_processes ) {
    my ($fhb, $filenameb) = tempfile("childbed_$i-XXXX", UNLINK => 1, DIR => ".");
    push @tempfiles, { fhb => $fhb, fileb => $filenameb};
}

##########

### 
$pm->run_on_finish(sub {
    my ($pid, $exit_code, $ident, $exit_signal, $core_dump) = @_;
	 if ($exit_signal) {
        die "Child $pid died from signal $exit_signal\n";
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

	my $sid = $ar[0];
	if($sid =~ /(.+)\/[1|2]$/){
		$sid = $1;
	}

	if ( $lc >= $line_bunch and $sid  ne $lastr ) { # Process in batches of 1000 lines
        print "processing line\n" if $debug;
		process_lines([@lines]);
		$lc = 0;
        @lines = (); # Clear the buffer
    }
	push @lines, $_;
	$lastr = $sid;
}

process_lines([@lines]) if @lines; # Process remaining lines

$pm->wait_all_children; # Ensure all processes complete


my @childsfile = map { $_ -> {fileb} } @tempfiles;
warn "concatename and sort @childsfile\n";
system("sort -k 2,2n -k 3,3n -k 4,4 -k 5,5nr  --parallel=$num_processes --buffer-size=2G  -T ./  @childsfile  > $bedout ");

#foreach my $tf (@tempfiles) {
#    open my $inb, '<', $tf->{fileb} or die "Can't open $tf->{fileb}: $!";
#    print "--- Contents of $tf->{fileb} ---\n";
#    while (<$inb>) {
#        print BED $_;
#    }
#    close $inb;
#    unlink $tf->{fileb};  # Clean up
#}

sub process_lines {
    my ($lr) = @_;	

	#### choose file handlei
	my ($fhi) =  getfh(\@tempfiles,\%fh_assign) ;  	
	

	my @lines_sub = @{$lr};
	
	my $pid = $pm->start($fhi);
	if($pid){
		print "child $$ parent $pid ... FHI:$fhi\n";
		$fh_assign{$fhi} = 1;
		return;
	}

	my $fhb = $tempfiles[$fhi] -> {fhb};

	my %tout; # out put for the total lines
	my $batch = scalar @lines_sub;
	print "BATCH $batch\n";
	### get node leng
	my %bestfriends;
	my $lc = scalar @lines_sub;

	my %as_rcd;
	
	for (my $h = 0; $h < @lines_sub ; $h += 2 ){	
		my @ar1 =  split /\t/, $lines_sub[$h]; # @{ $lines[$h] };
		my @ar2 =  split /\t/, $lines_sub[$h + 1]; # @{ $lines[$h + 1] };
		print "ARRAY1:@ar1\n" if $debug ;
		print "ARRAY2:@ar2\n" if $debug ;
		if($ar1[5] eq "*" or  $ar2[5] eq "*"){
			next;
		}

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
	
		my $rid = $ar1[0];
	
		my $as1 = $tag1{AS} // 0;
		my $as2 = $tag2{AS} // 0;
		
		my $as = $as1 + $as2;
		if( $as_rcd{$rid} ){
			### filter for each read, only keep the best alignmes - 10;
			if($as_rcd{$rid}{S} > $as + $score_diff  ){
				## already have a better pair
				next;
			}
		}else{
			### output  last rid;
			my ($lrid) = keys %as_rcd;
			if($lrid){
				my @ils = @{$as_rcd{$lrid}{L}};	
				my $lc = @ils;
				foreach my $l ( @ils ){
					print $fhb "$l\t$lc\n";
				}
			}

			%as_rcd = ();
			$as_rcd{$rid}{S} = $as;
			# look ahead the next line,to see it is the same read
			if( $h+2 < @lines_sub and $lines_sub[$h + 2] =~ /^$rid/){
				$as_rcd{$rid}{U} = 1;
			}else{
				$as_rcd{$rid}{U} = 0;
			}
		}

		## mapping statistics for each read
		my $r1 = join('||', 
			$ar1[1],
			$as1,
			$tag1{cs} // '*',
			$tag1{dv} // '*'
		);
		my $r2 = join('||',
			$ar2[1],
			$as2,
			$tag2{cs} // '*',
			$tag2{dv} // '*'
		);
			
		## nodes on reads	
		## nodes on reads
		my($sn1,$en1, $sn2,$en2) ;
		my @nodes1;
		while($ar1[5] =~ /(<|>)(\d+)/g){
			$sn1 = $2 if !defined($sn1) || $2 < $sn1;
			$en1 = $2 if !defined($en1) || $2 > $en1;
			push @nodes1, $1 eq ">" ? $2 : "-$2";
		}
		my @nodes2;
		while($ar2[5] =~ /(<|>)(\d+)/g){
			$sn2 = $2 if !defined($sn2) || $2 < $sn2;
			$en2 = $2 if !defined($en2) || $2 > $en2;
			push @nodes2, $1 eq ">" ? "-$2" : $2;
		}
		@nodes2 = reverse @nodes2;

		# format nodessss, foreach node , include length ,, head and tail ?
		my $p1;
		my $p2;
		if($ar1[7] eq "*" or $ar1[8] eq "*" ){
			warn "@ar1 \n @ar2\n";
			$p1 = "*";
		}else{
			$p1 = form_node(\@nodes1,$nodelen, $ar1[7]   ,$ar1[6] - $ar1[8]  );
		}
		if($ar2[7] eq "*" or $ar2[8] eq "*" ){
			$p2 = "*";
		}else{
			$p2 = form_node(\@nodes2,$nodelen, $ar2[6] - $ar2[8]  , $ar2[7]   );
		}

		print "N1 $p1  @nodes1\n" if $debug;
		print "N2 $p2 @nodes2\n" if $debug ;
		
		#my ($sn, $en) = sort {$a <=> $b} ($nodes1[-1], $nodes2[0]);
		my ($snr,$enr) = (sort { $a <=> $b } grep { defined } ($sn1, $sn2, $en1, $en2))[0, -1];
		$snr --;	

		### overlap
		my $overlap = 0;
		print "test overlap @nodes1 = @nodes2\n" if $debug;
		if(@nodes1 and @nodes2){
			my $ii =  @nodes1 < @nodes2 ? scalar @nodes1 : scalar @nodes2;
			print "itr $ii\n" if $debug ;
			for my $len (1 .. $ii ) {
				# take last $len elems from @a1 and first $len elems from @a2
				my @tail = @nodes1[-$len .. -1];
				my @head = @nodes2[0 .. $len-1];
				print "TTTT:@tail = @head\n" if $debug;
				if ("@tail" eq "@head") {
					$overlap = $len;
				 }
			}
		}
		my $rang = "$snr\t$enr"; ### out put line the node range / coordinates
		
		push @{$as_rcd{$rid}{L}}, "G\t$rang\t$rid\t$as\t$overlap\t$p1\t$p2\t$r1\t$r2";
		#print $fhb "G\t$rang\t$rid:$as_rcd{$rid}{U}\t$as\t$overlap\t$p1\t$p2\t$r1\t$r2\n";	
		#print "G\t$rang\t$r1\t$r2\n" if $debug; 	
		#out ( \@allnodes, $nodelen, $ar1[0], $ar1[7], $ar2[7], $as, $multimap{$ar1[0]} , $fhb );
	} 
	print "subprocess, finished\n" if $debug;
	$fhb -> flush;
	$pm->finish(0, \%tout); 
}

#######################
sub form_node{
    my ($ref ,$nodelen, $hr,$tr) = @_;
    my $path = "";

    my @ns = @{ $ref };

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
        }

        my $l = $$nodelen{abs($ns[$i])} ? $$nodelen{abs($ns[$i])} : 1;

        $path = $path? "$path;$ns[$i],$l,$h,$t" : "$ns[$i],$l,$h,$t";
    }
    return $path;
}


sub getfh { 
    my ($tempfiles,$fh_assign) =   @_;

    my $fhi;
    while(! defined $fhi){
        for(my $i = 0 ;$i < @$tempfiles; $i ++){
            unless (  defined $$fh_assign{$i} ){
                $fhi = $i;
                $$fh_assign{$i} = 1;
                last;
            }
        }
    }
    return ($fhi);
}





