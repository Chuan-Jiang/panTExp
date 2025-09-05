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


use Time::HiRes qw(gettimeofday tv_interval);

my $gaf_input;
my $pre;
my $debug;

my $num_processes = 1; # Number of parallel processes
my $CST = 0.8;
my $lib = 1000;
my $line_bunch = 500000;
my $nodelen_f;

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
my %fh_assign;

# Extract directory part from prefix
my $dir = dirname($pre);

# If directory doesn't exist, create it
unless (-d $dir) {
    make_path($dir) or die "Failed to create directory $dir: $!";
}

open BED, " > $pre.gaf2.bed" or die $!;

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



################
## create temp files

my @tempfiles;
for my $i (0 .. $num_processes ) {
    my ($fhb, $filenameb) = tempfile("childbed_$i-XXXX", UNLINK => 1, DIR => ".");
    push @tempfiles, { fhb => $fhb, fileb => $filenameb};
}

##########

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

foreach my $tf (@tempfiles) {
    open my $inb, '<', $tf->{fileb} or die "Can't open $tf->{fileb}: $!";
    print "--- Contents of $tf->{fileb} ---\n";
    while (<$inb>) {
        print BED $_;
    }
    close $inb;
    unlink $tf->{fileb};  # Clean up

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

	my @lines_sub = @{$lr};
	
	my $pid = $pm->start($fhi);
	if($pid){
		print "child $pid ... FHI:$fhi\n";
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
	for (my $h = 0; $h < @lines_sub ; $h += 2 ){	
		my @ar1 =  split /\t/, $lines_sub[$h]; # @{ $lines[$h] };
		my @ar2 =  split /\t/, $lines_sub[$h + 1]; # @{ $lines[$h + 1] };
		print "ARRAY1:@ar1\n" if $debug ;
		print "ARRAY2:@ar2\n" if $debug ;
		if($ar1[5] eq "*" and $ar2[5] eq "*"){
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

		my $r1 = join('||', 
			$ar1[1],
			$tag1{AS} // '*',
			$tag1{cs} // '*',
			$tag1{dv} // '*'
		);
		my $r2 = join('||',
			$ar2[1],
			$tag2{AS} // '*',
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

		my $p1;
		my $p2;
		if($ar1[7] eq "*" or $ar1[8] eq "*" ){
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
		my $max = 0 ; # best index should have the higest overlap
		for ( my $i = 1; $i <= @nodes1 && $i <= @nodes2; $i ++ ){
			my $total = 0;
			my $match = 0;
			for (my $j = 0; $j < $i; $j ++){
				my $n1 = $nodes1[-$i + $j];
				my $n2 = $nodes2[$j];
				my $n1l = $$nodelen{abs($n1)} // 1;
				my $n2l = $$nodelen{abs($n2)} // 1;
				if($n1 eq $n2){
					$total += $n1l * 2;
					$match += $n1l * 2;
				}else{
					$total += ($n1l + $n2l) ;
				}
			}
			if($match/$total > $max){
				$overlap = $i;
				$max = $match/$total;
			}
		}

		my $rang = "$snr\t$enr";

		print $fhb "G\t$rang\t$overlap\t$p1\t$p2\t$rid\t$r1\t$r2\n";	
		print "G\t$rang\t$r1\t$r2\n" if $debug; 	
		#out ( \@allnodes, $nodelen, $ar1[0], $ar1[7], $ar2[7], $as, $multimap{$ar1[0]} , $fhb );
	} 
	print "subprocess, finished\n" if $debug;
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

