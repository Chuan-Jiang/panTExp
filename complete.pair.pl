#!/usr/bin/perl
use warnings; use strict;
use Parallel::ForkManager;
use Getopt::Long;
use Bio::DB::HTS::Tabix;
use IO::Handle;

use List::Util qw(min max uniq sum);
use File::Temp qw(tempfile);
use File::Basename;
use File::Path qw(make_path);




my $num_processes = 1;
my $line_bunch = 100000;
my $nodenei_f ;
my $input_f;
my $pre;
my $frag_len = 1000;
my $debug;
my $samplename = "NA";
#$| = 1;

GetOptions(
    "i|input=s"  => \$input_f,
    "p|prefix=s" => \$pre,
    "b|neib=s" => \$nodenei_f,
    "d|debug"     => \$debug,
    "n|num_processes=i" => \$num_processes, # Allow setting number of processes
	"s|samplename=s" => \$samplename
) or die "Invalid options!\n";

###
# Extract directory part from prefix
my $pdir = dirname($pre);

# If directory doesn't exist, create it
unless (-d $pdir) {
    make_path($pdir) or die "Failed to create directory $pdir: $!";
}
#

open KEP, "> $pre.keep.tsv" or die $!;
open DIS, "> $pre.disc.tsv" or die $!;

srand(1);

my $pm = Parallel::ForkManager->new($num_processes);

##### determine input format  
my $in_fh;
if ($input_f) {
    if ($input_f eq "-") {
        $in_fh = *STDIN;
    } else {
        open my $fh, "zcat $input_f |" or die $!;
        $in_fh = $fh;
    }
} else {
    $in_fh = *STDIN;
}


##### create temp files

my @tempfiles;
my %fh_assign;
for my $i (0 .. $num_processes ) {
    my ($fhk, $filenamek) = tempfile("childkeep_$i-XXXX", UNLINK => 0, DIR => ".");
    my ($fhd, $filenamed) = tempfile("childdisc_$i-XXXX", UNLINK => 0, DIR => ".");
    push @tempfiles, { fhk => $fhk, filek => $filenamek ,fhd => $fhd, filed => $filenamed };
}


###### clean subprocess 
$pm->run_on_finish(sub {
    my ($pid, $exit_code, $ident, $exit_signal, $core_dump) = @_;
    if ($exit_signal) {
        print "Child $pid died from signal $exit_signal\n";
    }
    print "file handle $ident released by $pid..\n";
    delete $fh_assign{$ident};
});


################
###  MAIN   ####

my @lines;
my %lnode;

my $lastr = 0;   
my $lc = 0;
my $tc = 0;

while(<$in_fh>){
	chomp;
	my $rl = $_;

	$lc ++;
    $tc ++;
    my @ar = split /\t/;
	my @ns = map { /(?:^|;)(-?\d+)/g } @ar[4,5];
	#print "reads matched reads: @ar[4,5] => @ns\n" if $debug;

    if ( $lc >= $line_bunch and $ar[1] >  $lastr   ) { # Process in batches of 1000 lines
		my $gap = 1;
		foreach my $n (@ns ){
			if($lnode{ $n }){
				$gap = 0 ;
				last;
			}
		}
		if($gap){
			print "processing line\n" if $debug;
			
			my ($fhi, $fhk, $fhd) = getfh(\@tempfiles,\%fh_assign) ;
			
			###
			my $pid = $pm->start($fhi);
			if($pid){
				print "PARENT:child $pid with fhi $fhi... \n";
				$lc = 0;
				@lines = (); # Clear the buffer
				%lnode = ();
				#next;
			}else{
				print "\t^---CHILD: assiged with jobs, get filehandle $fhi  \n";
				print "inchild $lc $tc in line to process\n" if $debug;	
				process_lines(\@lines, $fhk,$fhd);
				$pm -> finish;
			}
		}
    }
    push @lines, $rl;
    map { $lnode{ $_ } = 1 } @ns; 
	$lastr = $ar[2];
}

if(@lines){
	my ($fhi, $fhk, $fhd) = getfh(\@tempfiles,\%fh_assign) ;
	process_lines(\@lines, $fhk , $fhd ) if @lines;
}

$pm -> wait_all_children;


### LAST combine outputs ###########
foreach my $tf (@tempfiles) {
    open my $ink, '<', $tf->{filek} or die "Can't open $tf->{filek}: $!";
    print "--- Contents of $tf->{filek} ---\n";
    while (<$ink>) {
        print KEP  $_;
    }
    close $ink;
	unlink $tf->{filek};  # Clean up

    open my $ind, '<', $tf->{filed} or die "Can't open $tf->{filed}: $!";
    print "--- Contents of $tf->{filed} ---\n";
    while (<$ind>) {
        print DIS  $_;
    }
    close $ind;
	unlink $tf->{filed};  # Clean up
}



### SUB routine ####

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
	my $fhk; my $fhd;
	$fhk = $$tempfiles[$fhi] -> {fhk};
	$fhd = $$tempfiles[$fhi] -> {fhd};
	return ($fhi, $fhk,$fhd);
}

sub process_lines {
	my ($lr,$fhk,$fhd) = @_;

	my @lines_s = @{$lr};
	## choose hanle index

	my %graph;
	my %complete;
	my %lnode;

	my $last_p;
	my $init_p;

	my $tabix = Bio::DB::HTS::Tabix->new(filename => $nodenei_f)  if $nodenei_f ; 
	
	for(my $i = 0; $i < @lines_s; $i ++){
		my $line_w = $lines_s[$i];
		print "NEW ==> $lines_s[$i]\n" if $debug;

		my @ar = split /\t/, $lines_s[$i];
		my @ns = map { /(?:^|;)(-?\d+)/g } @ar[4,5];
		if($last_p and $ar[1] > $last_p ){	
			
			my $gap = 1;
			foreach my $n (@ns){
				if($lnode{$n}){
					$gap = 0;
					last;
				}
			}
			if($gap){

				## finish completess
				# check graph file
				print "gap $last_p  $ar[1] ..\n" if $debug;
				
				inside_tabix($tabix, $init_p, $last_p ,\%graph);
				
				my $lap = 0;
				foreach my $n ($ar[1] .. $ar[2]){
					if(exists $graph{$n} ){
						$lap ++;
					}
				}
				print "over $lap\n" if $debug;
				unless ( $lap ){
					link_it (\%complete,\%graph,$tabix,$init_p,$last_p, $fhk,$fhd);
					
					undef %graph ;
					undef %complete;
					$last_p = undef;
					$init_p = undef;
				}else{
					$last_p = $ar[2];
				}
			}
		}else{
			$last_p = $ar[2];
		}
		
		$init_p = $init_p // $ar[1];
		
		map { $lnode{$_} = 1 } @ns;

		my @nodes1 = @{ build($ar[4],\%graph) };
		my @nodes2 = @{ build($ar[5],\%graph) };

		my @inf1 = split /\|\|/, $ar[7];
		my @inf2 = split /\|\|/, $ar[8];
	
		###### quality checgg
		# cs check
		my $disc = "";
		my $csr1 = cs($inf1[2])/$inf1[0] ;
		my $csr2 = cs($inf2[2])/$inf2[0] ;

		if($csr1  < 0.9 or $csr2 < 0.9 ){
			$disc .= "(cs=$csr1,$csr2)";
		}
		# map stat check, target are shorter than read
		if($ar[5] eq "*"){
			$disc .= "(t1=*)";
		}
		if($ar[6] eq "*"){
			$disc .= "(t2=*)";
		}
		
		# one node = overlap check 
		if($ar[3] ==  1){
			my @n1 = @{$nodes1[-1]};
			my @n2 = @{$nodes2[0]};

			my $nl = $n1[1];
			my $nodeleft = $n1[3] + $n2[2] - $nl ; # negtive means have overlap.
			print "\tnodeleft $nodeleft $n1[3] + $n2[2] - $nl\n" if $debug ;
			
			if($nodeleft + $inf1[0] + $inf2[0] > $frag_len  ){
				$disc = ".(M)";
				#}elsif($nodeleft + $inf1[0] < 0 and $nodeleft + $inf2[0] < 0){
				#$disc = ".(D)";
			}
		}
		
		###### firt round quality check finished : reads alignment based
		unless($disc){
			my $read_num = $tc + $i;
			my $read = "$ar[6]:$read_num";
			my $score = $inf1[1] + $inf2[1];
			my $rlen  = $inf1[0] + $inf2[0];

			my ($ns,$ne) = ($ar[1],$ar[2]);
			
			if($ar[3] > 0 ){
				# under this situation, should output to keep file : no $disc and have overlap
				# merge nodes, make node string
				my @path_parts;
				my @path_nodes;
				my $frag  = 0;

				for(my $j = 0; $j < @nodes1; $j ++){
					my $over = $ar[3] - (@nodes1 - $j );
					if($over < 0){
						my @n1 = @{ $nodes1[$j] };
						push @path_parts , join ",", @n1, 1; 
						push @path_nodes , $n1[0];
						$frag += $n1[1] - $n1[2] - $n1[3]; 
					}else{
						my @n1 = @{ $nodes1[$j] } ;
						my @n2 = @{ $nodes2[$over] } ;
						my $h = min($n1[2],$n2[2]);
						my $t = min($n1[3],$n2[3]);
						my $mn = "$n1[0],$n1[1],$h,$t,3";
						push @path_parts, $mn; 
						push @path_nodes,  $n1[0] ;
						$frag += $n1[1] - $h - $t;
					}
				}
				for (my $k = $ar[3]; $k < @nodes2; $k ++ ){
					my @n2 = @{ $nodes2[$k] };
					push @path_parts, join ",", @n2, 2 ;
					push @path_nodes,  $n2[0] ;
					$frag += $n2[1] - $n2[2] - $n2[3] ;
				}

				if($debug){
					print "\toverlaps nodes:\n";
					print "\t@path_parts\n";
					print "\t@path_nodes\n";
				}
				## merge path
				my $path_m = join ";", @path_parts;
				my @chunks = @{ chunks (\@path_nodes) };
				my $chunks_c = @chunks;
				for(my $m = 1; $m  <= $chunks_c; $m ++){
					print $fhk "G\t$chunks[$m-1]\t$score\t$frag\t$path_m\t$read\t$samplename\t$m\t$chunks_c\n";
					print " ===> G\t$chunks[$m-1]\t$score\t$frag\t$path_m\t$read\t$samplename\t$m\t$chunks_c\n" if $debug ;
				}
			}else{
				# have no overlap, and not sure about the middle situation. 
				print "\t need to be completed\n" if $debug;
				$complete{$nodes1[-1][0]}{$read} = [\@nodes1,\@nodes2 , $score, $rlen, $ns,$ne , $lines_s[$i] ];
			}
		}else{
			# have $disc defined, so this annotation must go to discard file
			print "\tDISCARD  $lines_s[$i]\t$disc\n" if $debug;
			print $fhd "$line_w\t$disc\n";
		}
	}
	if(%complete){
		link_it (\%complete,\%graph,$tabix,$init_p,$last_p, $fhk,$fhd);
	}
	##
	$fhk -> flush();
	$fhd -> flush();

	print "\tbachof lines preprocess finished\n" if $debug;
}

####
sub dest_checkin {
	my ( $complete, $bn, $back, $energy, $des , $record) = @_;
	
	foreach my $r (keys %{$$complete{$bn}} ){
		if ($$record{$r}){
			next;
		}else{
			$$record{$r} = 1;
		}

		foreach my $p ( keys %{ $$back{$bn} } ){
			my $en = $$complete{$bn}{$r}[1][0][0];
			# the bridge can not have length larger than $left
			my $left  = $frag_len - $$complete{$bn}{$r}[3];
			if($bn > 0){
				$left -= $$complete{$bn}{$r}[0][-1][3];
			}else{
				$left -= $$complete{$bn}{$r}[0][-1][2];
			}
			if($en > 0){
				$left -= $$complete{$bn}{$r}[1][0][2];
			}else{
				$left -= $$complete{$bn}{$r}[1][0][3];
			}
			####
			if($left > 0){
				print "\t$r registed at $p energy : $left aim from $bn -> $en \n" if $debug;
				$$des{$r}{A} = $en;
				$$des{$r}{S} = $bn;
				$$energy{$r}{$p} = $left ;
			}
		}
	}
}

#####
sub link_it {
	my($complete,$graph, $tabix,$init_p, $last_p, $fhk, $fhd )  = @_;
	print "####### \n Entering link model \n########\n" if $debug;	
	
	### collect numbers
	my @gns = keys %{$graph};
	my $nodes = keys %{$graph};
	my $comp = keys %{$complete};
	my $overs = 0;
	my $overe = 0;
	foreach my $c (keys %{$complete} ){
		if (exists $$graph{$c}){
			$overs ++;
			foreach my $e (keys %{$$complete{$c}}){	
				my $en = $$complete{$c}{$e}[1][0][0];
				if(exists $$graph{$en}){
					$overe ++;
				}
			}
		}
	}

	print "range $init_p = $last_p ; nodes in graph $nodes complete ends: $comp = $overs $overe @gns \n" if $debug;

	my %seen; # record checked nodes
	my %records; # record checked reads
	# find the start and end nodes
	foreach my $sn ( sort {$a <=> $b} keys %{$complete} ){	
		if ( $seen{$sn} ) {
			print "$sn checked. next" if $debug;
			next;
		}
		### checkin destination nodes
		my %back; #  node -> path with this node as the last one
		my @paths = ("$sn,$$graph{abs($sn)}{L},0,0,0");  # collectio of pathes
		$back{$sn}{0} = 1;  # the last nodes point to path index
	
		my %energy ; # the left length of one read. key are path index => reads => length
		my %des; # have key as "endnode", value as read name	
	
		print "\ninit process with node $sn  $$graph{abs($sn)}{L},0,0,0 \n" if $debug;

		print "NOLENGTH $sn\n" unless $$graph{abs($sn)}{L} and  $sn; 
		# %back have all the leading nodes as keys.
		while( %back ){
			########		
			my @iters = keys %back;
			print "nodes at end @iters\n" if $debug;

			### too many pathes, ignore....
			if(@paths > 100){
				foreach my $r ( keys %energy ){
					my $start = $des{$r}{S};
					my $line = $$complete{$start}{$r}[6] ;
					print $fhd "$line\t(CP)\n";
					delete $energy{$r};
					delete $$complete{$start}{$r};
				}
				last;
			}

			foreach my $bn (@iters){  # $b, base node
				# begain to growth path
				if($debug){
					print "\ngrowth node :  $bn\n";	
				}	
				# if a leading edge nodes are the edge of reads pair.  then check_in /
				if($$complete{$bn}){
					print "register new tadpole \n" if $debug;
					dest_checkin ( $complete, $bn , \%back , \%energy, \%des , \%records); # completehash, base_node, des_hash, paths_hash, curent_paths
				}
				$seen{$bn}  = 1;

				# extending graph according to requirement
				unless (exists $$graph{$bn} ){
					print "$bn not in reads, need check graph???\n" if $debug;
					extend_tabix($bn,$tabix,$graph);
				}

				# test if SNP bubble were identified. if it is a snp, then merge it as the same node group. 
				print "\ttest if  SNP bubble been found\n" if $debug ;
				my %ens ; # extending nodes , the representive one as the first key, all are used the second level keys. values are depthes. 
				my %forwards;
				foreach my $nn ( keys %{$$graph{$bn}{C}} ){
					#if(abs($nn) >= $last_p or abs($nn) <=  $init_p){
					#	print "outside of boundary $nn \n";
					#	next;
					#}
					if($seen{$nn}){
						# print "previous nodes encountered.. skip\n";
						#next;
					}

					unless (defined $$graph{abs($nn)}{L}){
						print  "$nn length  not exists\n" if $debug;
						next;
					}
					if( $$graph{abs($nn)}{L} == 1 ){ ### only consider SNPs, not INDELs
						my $nnn = join ",", (sort keys %{$$graph{$nn}{C}});
						print "$nn have length 1 => $nnn ,\n" if $debug;
						$forwards{$nnn}{$nn} = $$graph{$bn}{C}{$nn} ;
					}else{
						print "$nn have length > 1 $$graph{abs($nn)}{L}   \n" if $debug;
						$ens{$nn}{$nn} = $$graph{$bn}{C}{$nn} ;
						$ens{$nn}{D} += $$graph{$bn}{C}{$nn};
					}
				}
				
				foreach my $f ( keys %forwards ){
					my @nnns = sort keys %{$forwards{$f}};
					my @depth = map { $$graph{$bn}{C}{$_} } @nnns;
					my $ci = choose (@depth);

					if($debug){
						print "Depth for next nodes: @nnns = @depth\n";
						print "$f == @nnns choose $ci \n";
					}
					map {$ens{$nnns[$ci]}{D}  += $forwards{$f}{$_} } @nnns;
				}	
					
				###### get idxes for a specific leading node
				my @path_idxs = keys %{$back{$bn}}; # path index with this $bn at the end
				
				
				print "iterating path index : @path_idxs extending from $bn\n" if $debug;	
				## iterating each previous paths
				my %tadpole_alive;
				foreach my $pi ( @path_idxs ){
					my $last = $#paths ;
				
					my @ens_rep = keys %ens;
				
					# iterating next steps, each potential nodes will have a potential path index
					for(my $j = 0; $j < @ens_rep ; $j ++){
						### determine the path index
						my $ib = $j * @ens_rep;

						my @nidx = ( $ib + $pi, $ib + $last+1 .. $ib + $last + $#ens_rep ) ; # current idex and series of new idex in case of @ens have more then one element
						print "extend or create new path idex : (@nidx) for new ends: ( @ens_rep ) based on original idx $pi \n" if $debug;

						### collect infor about this extending node 
						# node id
						my $en = $ens_rep[$j];	
						unless(exists $$graph{$en} ){
							print "$en at the bobundary..\n" if $debug;
							next;
						}
						# type of nodes, "" form graph file 0 from other reads
						my $type = "";
						if($ens{$en}{D} > 0){
							$type = 0;
						}
						# groups of nodes associateed with $en
						my @en_as = keys %{$ens{$en}};
						# the string represented by this node
						my $en_s = join "|", @en_as;

						# node length
						my $en_l = $$graph{ abs($en) }{ L };
						
						# the path index the node will go into
						my $en_i = $nidx[$j];
						
						print "EXTENDING NODE: $en, $en_s, $en_l, $en_i\n" if $debug ;

						print "Test if read are still valide if extending to $en\n" if $debug;
						#delete $$graph{$bn}{C}{$en}; # delete seen path, avoid endless loop	
						## each read is a tadpole, it is trying find its mom (dest_node). tadpole have energy to a maximum distance that is determing the the library size.
						foreach my $r ( keys %energy ){
							unless($energy{$r}{$pi}){
								next;
							}
							print "energy left  $r $energy{$r}{$pi} - $en_l \n" if $debug;
							if(  $des{$r}{A} eq  $en ){
								my ($ipath) = $paths[$pi] =~ /(?:$des{$r}{S},\d+,\d+,\d+,\S;)(.*)/;
								$ipath = "" unless $ipath ;
								my $frag = 0;
								while($ipath =~ /\d+,(\d+),\d+,\d+,\d*;?/g){
									$frag += $1;
								}

								print "found mom  $r  $des{$r}{S} => $en:$ipath:$paths[$pi]\n" if $debug ;
								my @path_parts;
								my @path_nodes;
								@path_nodes = $ipath =~ /(?:;|^)(-?\d+)/g if $ipath ; 
								print "iner path nodes: @path_nodes <<< $ipath <<< $paths[$pi] \n" if $debug;

								my @nodes1 = @{$$complete{ $des{$r}{S} }{$r}[0] };
								my @nodes2 = @{$$complete{ $des{$r}{S} }{$r}[1] };
								for(my $i = 0; $i < @nodes1  ; $i ++){
									my @n1 = @{ $nodes1[$i] };
									push @path_parts , join ",", @n1, 1; 
									push @path_nodes , $n1[0];
									$frag  += $n1[1] - $n1[2] - $n1[3] ;
								}
	
								if ($ipath){
									push @path_parts, $ipath;
								}

								for(my $i = 0; $i < @nodes2  ; $i ++){
									my @n1 = @{ $nodes2[$i] };
									push @path_parts , join ",", @n1, 2;
									push @path_nodes,  $n1[0] ;
									$frag  += $n1[1] - $n1[2] - $n1[3] ;
								}

								if($debug){
									print "\toverlaps nodes:\n";
									print "\t@path_parts\n";
									print "\t@path_nodes\n";
								}
								## merge path
								my $path_m = join ";", @path_parts;
								my @chunks = @{ chunks (\@path_nodes) };
								my $chunks_c = @chunks;
								for(my $i = 1; $i  <= $chunks_c; $i ++){
									print $fhk "G\t$chunks[$i-1]\t$$complete{ $des{$r}{S} }{$r}[2]\t$frag\t$path_m\t$r\t$samplename\t$i\t$chunks_c\n";
									print " ===> G\t$chunks[$i-1]\t$$complete{ $des{$r}{S} }{$r}[2]\t$frag\t$path_m\t$r\t$samplename\t$i\t$chunks_c\n" if $debug ;
								}
								delete 	$des{$r};
								delete $energy{$r};
							}elsif( $energy{$r}{$pi} - $en_l  < 0 ){
								print "energy used up for   $r to extend $en ...discard this reads\n" if $debug;	
								delete $energy{$r}{$en_i};
							}else{
								print "keep finding...\n" if $debug;
								$energy{$r}{$en_i} -= $en_l;
							}
						}
						$paths[$en_i] =  "$paths[$pi];$en,$$graph{abs($en)}{L},0,0,$type"  ;
						$back{$en}{$en_i} = 1;
					}
				}
				#####
				delete $back{$bn};
			}
			unless (keys %des ){
				last;
			}
		}
		foreach my $r ( keys %energy ){
			my $start = $des{$r}{S};
			my $line = $$complete{$start}{$r}[6] ;
			print $fhd "$line\t(DIS)\n";
		}
	}
	print "Link Module finished\n############\n" if $debug;
}


sub choose{
	my @size = @_;

	if(@size == 1){
		return 0;
	}

	my $total = sum(@size);
	if($total <=0){
		return 0;
	}

	my $rnd = rand($total);

	my $acc = 0;
	my $chosen;
	for(my $i = 0; $i < @size; $i ++){
		$acc += $size[$i];
		if ($rnd < $acc) {
			return $i;
		}
    }
}


### insize_tabix
sub inside_tabix {
	my ($tabix,$min, $max, $graph) = @_;
	
	my $win = "G:$min-$max";
	print "inside check $win..\n" if $debug;

	my $iter = $tabix -> query($win);

	#####
	while (my $line = $iter -> next) {
		my @ar = split /\t/, $line;
		
		my $cn = $ar[2];
		$$graph{$cn}{L} = $ar[5];
		while( $ar[3] =~ /(-?\d+)/g){                                                                   
			my $n = $1 * -1; 
			unless ( defined $$graph{-$cn}{C}{$n} ){		
				$$graph{-$cn}{C}{$n} = 0;
			}
		}
		while( $ar[4] =~ /(-?\d+)/g){
			my $n = $1 ;
			unless( defined $$graph{$cn}{C}{$n} ) {
				$$graph{$cn}{C}{$n} = 0;
			}
		}
	}
}


### 
sub extend_tabix {
	my($bn,$tabix,$graph) = @_;

	print "\textending $bn...\n" if $debug;	
	my $start  = abs($bn) - 20;
	my $end = abs($bn) + 20;

	my $win = "G:$start-$end";
	my $iter = $tabix -> query($win);


	my $found = 0;
	my @lines;

	while (my $line = $iter -> next) {
		my @ar = split /\t/, $line;
		push @lines, $line;
		
		while( $ar[3] =~ /(-?\d+)/g){
			my $n = $1  * -1 ;
			if(defined $$graph{$n}){
				$found ++;
			}
		}
		while( $ar[4] =~ /(-?\d+)/g){
			my $n = $1 ;
			if(defined $$graph{$n}){
			   $found ++;
			}
		}
	}
	if($found == 0){
		return 0;
	}
	foreach my $l (@lines){
		my @ar = split /\t/, $l;
		my $cn = $ar[2];
		$$graph{$cn}{L} = $ar[5];
		while( $ar[3] =~ /(-?\d+)/g){                                                                   
			my $n = $1  * -1 ; #> 0 ?  "+$1" : $1 ;
			unless ( defined $$graph{-$cn}{C}{$n} ){
				$$graph{-$cn}{C}{$n} = 0;
			}
			#$g1 -> add_weighted_edge(-$cn,$n,  $$nodelen{abs($n)} ? $$nodelen{abs($n)} : 1 );
		}
		while( $ar[4] =~ /(-?\d+)/g){
			my $n = $1 ; #> 0 ? "+$1" : $1 ;
			unless( defined $$graph{$cn}{C}{$n} ) {
				$$graph{$cn}{C}{$n} = 0;
			}
			#$g1 -> add_weighted_edge("+$cn" ,$n,  $$nodelen{abs($n)} ? $$nodelen{abs($n)} : 1 );
		}
	}
	return $found ;
}

sub build {
	my($path,$graph) = @_;	
	if($path eq "*"){
		return ([]);
	}

	my %pinf ;
	my @nodes = split /;/, $path;
	@nodes = map{ [split /,/] } @nodes;
	my @nodes_r = reverse @nodes;
	for (my $i = 0; $i < @nodes  ; $i ++ ) {
		if($i < $#nodes){
			# forward
			$$graph{$nodes[$i][0]}{C}{$nodes[$i+1][0]} ++;
			# reverse
			$$graph{-1*$nodes_r[$i][0]}{C}{-1*$nodes_r[$i+1][0]} ++;
		}
		if($nodes[$i][0] eq "*"){
			print "********$path\n";
		}
		$$graph{abs($nodes[$i][0])}{L} = $nodes[$i][1];
	}
	return (\@nodes);
}

#####
sub chunks {
	my $step = 100;
	my ($nr)  = @_;

	my @pn = map { abs ( $_ ) } @$nr ;
	@pn = sort { $a  <=> $b  } @pn;
	@pn = uniq(@pn);	
	my @chunks;
	my $startn;
	for (my $k = 0; $k < @pn; $k ++){
		my $n = abs($pn[$k]);
		unless( defined $startn){
			$startn = $n - 1;
			next;
		}
		my $n_l = abs($pn[$k - 1]);
		if(abs($n - $n_l) < $step){
			next;
		}else{
			push @chunks, "$startn\t$n_l";
			$startn = $n - 1;
		}
	}
	if($startn){
		my $ln = abs($pn[-1]);
		push @chunks, "$startn\t$ln";
	}
	return \@chunks;
}
#####
sub cs {
	my ($str) = shift @_;
	my $m = 0;
	while($str =~ /:(\d+)/g ){
		$m += $1;
	}
	return $m;
}

