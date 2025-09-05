#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Parallel::ForkManager;
use List::MoreUtils qw(upper_bound lower_bound);
use List::Util qw(min max);
use Storable ;
use File::Temp qw/tempdir tempfile/;


my $bed_f;
my $dis_f;
my $gaf_input;
my $pre;
my $debug;
my $num_processes = 1; # Number of parallel processes
my $step = 20;
my $nodelen_f;
my $nodegap_f;

GetOptions(
    "i|input=s"  => \$gaf_input,
    "p|prefix=s" => \$pre,
	"a|nodegap=s" => \$nodegap_f, ### determine wheter reads are properly mapped
    "s|step=i"  => \$step,
	"l|nodelen=s"    => \$nodelen_f,
	"b|bed"      => \$bed_f,
    "d|dis"      => \$dis_f,
    "g|bug"     => \$debug,
	"n|num_processes=i" => \$num_processes # Allow setting number of processes
) or die "Invalid options!\n";

#######
my $pm = Parallel::ForkManager->new($num_processes);
warn "paralell $num_processes\n";

if ($bed_f) {
    open BED, "| pigz  -p 4 > $pre.gaf2.bed.gz" or die $!;
}
if ($dis_f) {
    open GAF, "| pigz -p 4 > $pre.discord.gaf.gz" or die $!;
    open RDS, "| pigz  > $pre.discord.ids.gz" or die $!;
}

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

print "geting node len...\n" if $debug;
my $nodelen = retrieve($nodelen_f); 


print "got node len...\n" if $debug;

print "get node gaps...\n" if $debug;
my %nodegap;
open GAP , "$nodegap_f" or die $!;
while(<GAP>){
	chomp;
	my @ar = split /\t/;
	#my $number = int($ar[1] /1000) + 1;
	my($fn,$sn) = sort {$a <=> $b} @ar[2,3];
	$nodegap{$fn}{$sn} = $sn - $fn;
}
my @nodegap_ks = sort {$a <=> $b} keys %nodegap;
print "node gaps obtained\n" if $debug;

close GAP;

### 
$pm->run_on_finish(sub {
    my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $tout) = @_;
    if (defined $tout) {
		foreach my $id (keys %{$tout}) {
			if($$tout{$id}{ID}){
				print RDS "$id\n";
			}
			if($$tout{$id}{GAF}){
				print GAF join "\n", @{$$tout{$id}{GAF}};
				print GAF "\n";
			}
			if($$tout{$id}{BED}){
				if($debug){
					my $count = scalar @{$$tout{$id}{BED}};
					foreach my $l ( @{$$tout{$id}{BED}} ){
						print "BED$id:S$l E\n";
					}
				}
				if( @{$$tout{$id}{BED}} ){
					print BED join "" , @{$$tout{$id}{BED}};
				}
			}
		}

    }
});
##############

my @lines;
my $lastr = 0;
my $lc = 0;
while (<$gaf_fh>) {
	$lc ++;
    chomp;
    my @ar = split /\t/; 

	if ( $lc >= 100000 and $ar[0] ne $lastr ) { # Process in batches of 1000 lines
        print "processing line\n" if $debug;
		process_lines([@lines]);
		$lc = 0;
        @lines = (); # Clear the buffer
    }
	push @lines, \@ar;
	$lastr = $ar[0];
}

process_lines([@lines]) if @lines; # Process remaining lines

sub process_lines {
    my ($lr) = @_;
    $pm->start and return; # Fork a new procesfff

	my @lines = @{$lr};

	my %tout; # out put for the total lines
    my %pair;

	### get node leng
	for (my $i = 0; $i < @lines; $i ++){	
		my @ar = @{ $lines[$i] };
		print "lines: @ar\n" if $debug; 	
		my $next;
		my @ns_p;
		if($i % 2 == 0){
			$next = join ",", @{$lines[$i + 1]}[5,6,7,8];
			(@ns_p) =$lines[$i+1][5] =~ /(\d+)/g;
		}else{
			$next = join ",", @{$lines[$i - 1]}[5,6,7,8];
		}
	
		# test if two ends are close or not
		if(@ns_p){	
			my @ns = $lines[$i][5] =~ /(\d+)/g;	
			my($tip1,$tip2) = sort {$a <=> $b} ($ns_p[-1],$ns[-1]);
		
			@ns_p = sort {$a <=> $b} @ns_p;
			@ns = sort {$a <=> $b} @ns;
			my $close  = 0;
			
			my($k, $j) = (0,0);
			my $min_diff = 'inf';
			while ( $k < @ns && $j < @ns_p ) {
				my $diff = abs($ns[$k] - $ns_p[$j]);
				$min_diff = $diff if $diff < $min_diff;
				if($min_diff < 100 ){
					print "looks resonable.. skip.\n" if $debug;
					$close = 1;
					last;
				}
				# Move the pointer of the smaller value
				if ($ns[$k] < $ns_p[$j]) {
					$k++;
				} else {
					$j++;
				}
			}

			if( $close == 0  ){	
				print "need test..\n" if $debug;
				my $thr = 100;
				my $idx1 = upper_bound{ $_ <=>  $tip1 - $thr } @nodegap_ks;
				my $idx2 = lower_bound{ $_ <=>  $tip1 + $thr } @nodegap_ks;
				my $nd;
				for my $i ($idx1 .. $idx2){
					print "$tip1 $tip2, $i, $idx1 $idx2, $nodegap_ks[$idx1] $nodegap_ks[$idx2]\n" if $debug ;
					my $n = $nodegap_ks[$i];
					$nd = $n - $tip1;
					if($nd > 100){
						print "first attemp faled\n" if $debug;
					}else{
						my @linked = sort {$a<=>$b} keys %{$nodegap{$n}};
						my @dis = map { abs($_ - $tip2) } @linked;
						my $min_d = min(@dis);
						if($min_d + $nd < 100){
							print "have jump.. \n" if $debug;
							$close = 1; 
							last;
						}else{
							print "distance not bridege...\n" if $debug;
						}
					}

				}	
			}

			########################
			##################################
			###############  skip distantly mapped reads  #######
			unless ($close == 1){
				$i ++;
				print "this alignment skipp..\n\n" if $debug;
				next;
			}else{
				print "this alignment keepped ..\n\n" if $debug;
			}
		}


		my %tags;
        for (my $i = 12; $i < @ar; $i++) { # the last hashreference, shoult not ini for loop
            my ($t, $v) = $ar[$i] =~ /^([^:]*):\w:(.+)/;
            $tags{$t} = $v;
        }
        $tags{fn} ? push @ar, "2,$next"  : $tags{fp} ? push @ar, "1,$next" : undef;
		foreach my $k ("cs","AS","dv"){
			push @ar, defined $tags{$k} ? $tags{$k} : "NA";
		}

		#die "$line\n" unless $tags{cs};
        my $discor = 0;

        #die "$_\n" unless (defined $tags{pd});
        if($ar[2] eq "*"){
            $discor = 1;
        }elsif( defined $tags{pd} and  $tags{pd} == 0 ){ # Properly paired flag, by Xian
            $discor = 1;
        }elsif($tags{cs} =~ /^\+([ATCGatcg]+)/){
            $discor = 1 if (length($1) > 20);
        }elsif ($tags{cs} =~ /\+([ATCGatcg]+)$/){
            $discor = 1 if (length($1) > 20);
        }elsif(! defined $tags{dv} ){
            $discor = 1;
        }elsif($tags{dv} > 0.2){
            $discor = 1;
		}

		# decide whether new pair encountered.
        if ($pair{$ar[0]}) {
			print "collect...\n" if $debug;
			#push @{$pair{$ar[0]}{D}}, $discor;
			#push @{$pair{$ar[0]}{P}}, \@ar;
			push @{$pair{$ar[0]}}, [$discor,@ar];
        } else {
            print "\nnew line encounter. output previous...\n" if $debug;
			out(\%pair,\%tout,$nodelen);
            %pair = ();
			#push @{$pair{$ar[0]}{D}}, $discor;
			#push @{$pair{$ar[0]}{P}}, \@ar;
			push @{$pair{$ar[0]}}, [$discor,@ar];
        }
    }
    out(\%pair, \%tout,$nodelen );
	$pm->finish(0, \%tout); 
}

$pm->wait_all_children; # Ensure all processes complete

sub out {
    my ($pair,$tout,$nodelen)  = @_;
    return unless %{$pair};
	
    my ($id) = keys %{$pair};
	print "\tOUTID:$id\n" if $debug;
	
	my @lines = @{$$pair{$id}};
	for (my $i = 0; $i < @lines; $i ++){
		my $ii =  int($i/2) + 1;
		my($discor,@lar) = @{$lines[$i]};
		### to find discordant reads
		if($discor and $dis_f){
			$$tout{$id}{ID} = 1;
			push @{$$tout{$id}{GAF}}, join "\t",@lar;
		}
		### transform into bed format	
		if ($bed_f) {
			my $read =	$lar[-4]; # which read, 1 or 2?
			
			# extract node len
			my $path = "";
			while($lar[5] =~ /(.)(\d+)/g){
				my $d = $1 eq ">" ? 1 : -1;
				my $l ;
				if(%{$nodelen}){
					$l = $$nodelen{$2} ? $$nodelen{$2} : 1;
				}else{
					$l = 0;
				}

				$path = $path ? "$path;$2:$d,$l" : "$2:$d,$l";
			}
			my $anno = "$lar[0]\t$lar[1]\t$lar[-1]\t$lar[-2]\t$lar[-3]\t$path\t$lar[6]\t$lar[7]\t$lar[8]";
            #print "LINE $i @lar\n";
            my %ns ;
			while($lar[5] =~ /(<|>)(\d+)/g){
			    my ($id,$in) = ($1,$2);
				$ns{$in} = $id;
			}
			
			my @pn = sort {$a <=> $b} keys %ns;
			my @chunks;
			my $startn;
			for (my $k = 0; $k < @pn; $k ++){
				my $n = $pn[$k];
				unless($startn){
					$startn = $n;
					next;
				}
                my $n_l = $pn[$k - 1];
                if(abs($n - $n_l) < $step){
					next;
				}else{
					push @chunks, "$startn\t$n_l";
                    $startn = $n;
				}
			}
			if($startn){
				push @chunks, "$startn\t$pn[-1]";
			}
			print "\tCC:@chunks >> $read >> $anno >> $read\n" if $debug ;
			#push @{$beds{R}}, $read;
			my $bed_s;
			my $map_ct = scalar @lines;
			my $chunk_count = scalar @chunks;
			print "\t$chunk_count = @chunks\n" if $debug ;	
			for(my $j = 1; $j <= $chunk_count; $j ++){
				my $c = $chunks[$j-1];

				#G	start	end	anno(readid,dv,AS,cs,path,ln,start,end) howmanyreadsmapped  read which_chunk how_many _chunks
				print "\t===>G\t$c\t$anno\t$read\t$ii\t$j/$chunk_count\t$discor\n\n" if $debug;
				$bed_s .= "G\t$c\t$anno\t$read\t$ii\t$j/$chunk_count\t$discor\n";
			}
			push @{$$tout{$id}{BED}}, $bed_s if $bed_s;
		}
	}
}

sub tmp_f{
    make_path("./tmp") unless -d "./tmp";
        my %options = (
        TEMPLATE => 'mytempXXXXX',  # Default template
        DIR      => './tmp', # Default directory
        SUFFIX   => '.tmp',         # Default suffix
        UNLINK   => 1,              # Auto-delete file when closed
        @_,                        # Override defaults with any provided arguments
    );
    return tempfile(%options);
}



