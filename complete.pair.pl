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

my $gline_min = 2000;
my $gline_max = 5000;
my $bline_max = 100000;

my $nodenei_f ;
my $input_f;
my $l2g_ref_f;
my $g2l_ref_f;

my $pre;
my $debug ;
my $samplename = "NA";


### filters
my $filter_rd_iden = 0.9;
my $frag_len = 1000;

my %stat;
GetOptions(
    "i|input=s"  => \$input_f,
	"p|prefix=s" => \$pre,
    "b|neib=s" => \$nodenei_f,
    "l|l2g=s" => \$l2g_ref_f,
    "g|g2l=s" => \$g2l_ref_f,
	"d|debug"     => \$debug,
    "n|num_processes=i" => \$num_processes, # Allow setting number of processes
	"s|samplename=s" => \$samplename,
	"iden=f" => \$filter_rd_iden # the required percentage of identity with pan-genome for read

) or die "Invalid options!\n";

###
# Extract directory part from prefix
my $pdir = dirname($pre);

# If directory doesn't exist, create it
unless (-d $pdir) {
    make_path($pdir) or die "Failed to create directory $pdir: $!";
}
#
my $keep_file = "$pre.keep.tsv";
my $disc_file = "$pre.disc.tsv";
my $stat_file = "$pre.stat.tsv";
open STAT, "> $stat_file " or die $!;

srand(1);

my $pm = Parallel::ForkManager->new($num_processes);

##### determine input format
my $in_fh;
if ($input_f) {
    if ($input_f eq "-") {
        $in_fh = *STDIN;
    } elsif($input_f =~ /.gz/) {
        open my $fh, "zcat $input_f |" or die $!;
        $in_fh = $fh;
    } else {
		open my $fh, "$input_f" or die $!;
		$in_fh = $fh;
	}
} else {
    $in_fh = *STDIN;
}

##### create temp files

my @tempfiles;
my %fh_assign;
my @keeps ;
my @discs ;
my @mrcds ; 
for my $i (0 .. $num_processes ) {
	# fhd destination fhk keep fhr
    my ($fhk, $filenamek) = tempfile("childkeep_$i-XXXX", UNLINK => 1, DIR => ".");
	my ($fhd, $filenamed) = tempfile("childdisc_$i-XXXX", UNLINK => 1, DIR => ".");
    my ($fhr, $filenamer) = tempfile("childmrcd_$i-XXXX", UNLINK => 1, DIR => ".");
	push @keeps, $filenamek;
	push @discs, $filenamed;
	push @mrcds, $filenamer;
	push @tempfiles, { fhr => $fhr, filer => $filenamer, fhk => $fhk, filek => $filenamek ,fhd => $fhd, filed => $filenamed };
}

###### clean subprocess
$pm->run_on_finish(sub {
    my ($pid, $exit_code, $ident, $exit_signal, $core_dump,$data_ref) = @_;
    if ($exit_signal) {
        die "Child $pid died from signal $exit_signal $exit_code \n";
    }
    if (defined $data_ref) {
      my %hash = %$data_ref;
		#my %stat_c; # read decision statistics. D=discard, O=overlap, W=wander rescue , G=gadd`R=reference rescue
      $stat{D} += $hash{D};
      $stat{O} += $hash{O};
      $stat{W} += $hash{W};
      $stat{G} += $hash{G};
      $stat{R} += $hash{R};
    }
    print "\t===file handle $ident released by $pid..\n";
    delete $fh_assign{$ident};
});


################
###  MAIN   ####

my @lines;
my %lnode;

my $lastr = 0;
my $lc = 0;
my $gc = 0;
my $tc = 0;

while(<$in_fh>){
	chomp;
	my $rl = $_;

	my @ar = split /\t/;
	$lc ++;
	$gc ++ if $ar[5] == 0; # means no overlap between paired reads
	$tc ++;
	my @ns = map { /(?:^|;)(-?\d+)/g } @ar[6,7];
	#print "reads matched reads: @ar[4,5] => @ns\n" if $debug;

	my $gap = 1; # determing wheter current line have gap with the previous node
	foreach my $n (@ns ){
		if($n eq $ar[1] or $n eq $ar[2]){
			next;
		}
		if($lnode{ $n }){
			$gap = 0;
			last;
		}
	}
	# determine wheter it is a chunk of linked reads
	my $ana = 0;
	if ( $gc >= $gline_min and $ar[1] >= $lastr and $gap ) { # Process in batches of 1000 lines
		$ana = 1 ;
	}elsif($gc > $gline_max  and $ar[1] >= $lastr and $gap ){
		$ana = 2;
	}elsif($gc > $gline_max * 1.5 and $ar[1] >= $lastr){
		$ana = 3;
	}elsif($lc > $bline_max and $ar[1] >= $lastr and $gap  ){
		$ana = 4;
	}elsif($lc > $bline_max * 1.5 and $ar[1] >= $lastr){
		$ana = 5;
	}elsif( $lc > $bline_max * 2){
		$ana = 6;
	}elsif ($gc > $gline_max * 2){
		$ana = 7;
	}

	if($ana ){
		%lnode = ();
		process_lines(\@lines,$tc,$lc,$gc,$ana);
		@lines = ();
		$lc = 0;
		$gc = 0;
	}
	push @lines, $rl;
	map { $lnode{ $_ } = 1 } @ns;
	$lastr = $ar[2];
}

if(@lines){
	process_lines(\@lines,$tc,$lc,$gc,2 );
}

$pm -> wait_all_children;

### LAST combine outputs ###########

###################################
## collect multiple aligned reads  source is %rdrm 
my %muls;
foreach my $mf (@mrcds){
	open my $fh, $mf or die $!;
	
	while(<$fh>){
		chomp;
		my @ar = split /\t/;
		if( defined $muls{$ar[0]} ){
			my $r = $muls{$ar[0]} - $ar[2];
			if($r == 0){
				delete $muls{$ar[0]};
			}else{
				$muls{$ar[0]} = $r;
			}
		}else{
			$muls{$ar[0]} = $ar[1] - $ar[2];
		}
	}
	close $fh;
}


## prepare file handle
warn " final sort | sort -k 2,2n -k 3,3n  -k 5,5n -k 6,6 --parallel=$num_processes --buffer-size=2G  -T ./  > $keep_file \n";
open KEP,  " | sort -k 2,2n -k 3,3n  -k 5,5n -k 6,6 --parallel=$num_processes --buffer-size=2G  -T ./  > $keep_file " or die $!;
#########

my %chks;

foreach my $kf (@keeps){
	warn "reading and count $kf\n";
	open my $fh , $kf or die $!;
	while(<$fh>){
		chomp;
		my @ar = split /\t/;
		my $mul;
		my($r,$rn) = $ar[6] =~ /(.+):(\d+)$/; # guess: $rn is the uniq read id.
		
		my $mark;
		if ( !$muls{$r} ) {
			$mark = "0";
		}elsif ($ar[9] > 1 && $ar[8] > 1) {
			$mark = $chks{$rn};
			delete $chks{$rn} if $ar[8] == $ar[9];
		}else {
			$mark = $muls{$r};
			$chks{$rn} = $mark if $ar[9] > 1;   # only assign to chks when ar[9] > 1
			$muls{$r} = "$muls{$r}-1";
		}
		$ar[6] = "$r:$mark";
		print KEP join("\t",@ar),"\n";

	}
	close $fh;
}


my $deleted_count = unlink @keeps ;
print "$deleted_count files deleted\n";

warn "summarize disc file\n";
system("cat @discs > $disc_file");

print STAT "Total:\nLowQuality:$stat{D}\nOverlapRescue:$stat{O}\nWanderRescue:$stat{W}\nGap:$stat{G}\nReferRescue:$stat{R}\n";


warn "all done. output files: $keep_file , $disc_file , $stat_file \n";

# Catch Ctrl-C (SIGINT)
unless($debug){
	$SIG{INT} = sub {
		print "\nCaught Ctrl-C, cleaning up...\n";
		unlink @keeps, @discs or warn "Couldn't delete : $!";
	  	exit 1;
	};
}

#############################
#############################
### SUB routine #############
#############################

sub process_lines {
	my ($lr,$tc,$lc,$gc,$ana ) = @_;

	my @lines  = @{$lr};

	#### choose file handlei
	my ($fhi) =  getfh(\@tempfiles,\%fh_assign) ;
	print "enter child $tc $lc $gc $ana\n" if $debug;

	my $pid = $pm->start($fhi);
	if($pid){
		print "\t^--parent $$ child  $pid ... FHI:$fhi $tc $lc $gc $ana\n";
		$fh_assign{$fhi} = 1;
		return;
	}

	my $tabix = Bio::DB::HTS::Tabix -> new(filename => $nodenei_f)  if $nodenei_f ;
	my $tabix_l = Bio::DB::HTS::Tabix -> new(filename => $l2g_ref_f) ;
	my $tabix_g = Bio::DB::HTS::Tabix -> new(filename => $g2l_ref_f) ;
	
	my $fhk = $tempfiles[$fhi] -> {fhk};
	my $fhd = $tempfiles[$fhi] -> {fhd};
	my $fhr = $tempfiles[$fhi] -> {fhr};

	my %graph;
	my %collect; # first key is read name
	my %complete; # first key if node is, it used for link_it to fill gaps
	my %nodesb;  # collect all nodes in a lines block
	my %stat_c; # read decision statistics. D=discard, O=overlap, W=wander rescue , G=gap, R=reference rescue
	my %rdrm ; # read multiple mapping record
	
	
	my $last_p;
	my $init_p;
	for(my $i = 0; $i < @lines; $i ++){
		my $line_w = $lines[$i];
		print "NEW ==> $lines[$i]\n" if $debug;

		my @ar = split /\t/, $line_w;
		my $mul_c = pop @ar;
		my $rd = $ar[3];
		
		$rdrm{$rd}{R} = $mul_c ;
		
		##############################################################
		### code block to decide whether init line_it module  ########
		if($last_p and $ar[1] > $last_p ){	
			#VVVVVVVVVVVVVVVVVVVVVVVVVVVVV
			# if a gap, then should link and output
			## finish completess
			# check graph file
			print "gap $last_p  $ar[1] ..\n" if $debug;

			# collect relations from graph and update them into %graph
			inside_tabix($tabix, $init_p, $last_p ,\%graph);
			my $init_pr = $init_p;
			$init_p = undef;

			my $lap = 0;
			foreach my $n ($ar[1] .. $ar[2]){
				if(exists $graph{$n} ){
					$lap ++;
					last unless $debug;
				}
			}
			print "over $lap\n" if $debug;
			unless ( $lap ){
				if( keys %collect ){
					link_it (\%complete, \%collect, \%graph, \%nodesb, $tabix, $tabix_l,$tabix_g, $fhk, $fhd,\%stat_c  ,\ %rdrm);
				}else{
					print "no collected reads, skip link_it\n" if $debug;;
				}

				undef %graph ;
				undef %complete;
				undef %nodesb;
				undef %collect;

				$last_p = undef;
				$init_p = undef;
			}else{
				$last_p = $ar[2];
			}
		}

		$init_p = $init_p // $ar[1];
		my @inf1 = split /\|\|/, $ar[8];
		my @inf2 = split /\|\|/, $ar[9];

		###########################
		###### quality check ######
		###########################

		# sub cs check
		my $disc = "";
		
		my $csr1 = cs($inf1[2])/$inf1[0] ;
		my $csr2 = cs($inf2[2])/$inf2[0] ;

		if($csr1  < $filter_rd_iden  or $csr2 < $filter_rd_iden ){
			$disc .= "(cs=$csr1,$csr2)";
		}
		# map stat check, target are shorter than read
		if($ar[6] eq "*"){
			$disc .= "(t1=*)";
		}
		if($ar[7] eq "*"){
			$disc .= "(t2=*)";
		}

		# one node = overlap check, fragment length check
		
		my @n1 = $ar[6] =~ /(-?\d+),(\d+),(\d+),(\d+)$/;  #@{$nodes1[-1]};
		my @n2 = $ar[7] =~ /^(-?\d+),(\d+),(\d+),(\d+)/;  # @{$nodes2[0]};
		my $n1_rest = $n1[0] > 0 ? $n1[3] : $n1[2];
		my $n2_rest = $n2[0] > 0 ? $n2[2] : $n2[3];

		if($ar[5] ==  1){
			my $nl = $n1[1];
			my $nodeleft = $n1_rest + $n2_rest  - $nl ; # negtive means have overlap.

			if($nodeleft + $inf1[0] + $inf2[0] > $frag_len  ){
				print "\tdiscard nodeleft $nodeleft $inf1[0] $inf2[0]\n" if $debug ;
				$disc .= "(M)";
				#}elsif($nodeleft + $inf1[0] < 0 and $nodeleft + $inf2[0] < 0){
				#$disc = ".(D)";
			}
		}elsif($ar[5] == 0){
			if($n1_rest + $n2_rest + $inf1[0] + $inf2[0] > $frag_len ){
				print "\tdiscard gap $n1_rest + $n2_rest \n" if $debug ;
				$disc .= "(M)";
			}
		}


		my %uniq_nodes;
		my ($as, $ov, $l1 , $l2);
	
		#######################################
		## redundant mapping detection ########
		if($ar[6] ne "*" and $ar[7] ne "*"){
			## redundant detection
			$uniq_nodes{1} =  parse_node($ar[6]);
			$uniq_nodes{2} =  parse_node($ar[7]);

			$as = $ar[4];
			$ov = $ar[5];
			$l1 = $inf1[0];
			$l2 = $inf2[0];

			if(defined $collect{$rd}){
				foreach my $read_num (keys %{ $collect{$rd}  }){
					my %collectr = %{$collect{$rd}{$read_num}};

					my ($len1) = set_overlap($uniq_nodes{1}, $collectr{R1});
					my ($len2) = set_overlap($uniq_nodes{2}, $collectr{R2});
					my ($lenx1) = set_overlap($uniq_nodes{1}, $collectr{R2});
					my ($lenx2) = set_overlap($uniq_nodes{2}, $collectr{R1});

					print "found duplicates reads $read_num = $rd $len1 $len2 $lenx1 $lenx2 \n" if $debug;
					if( ( $len1/$l1 > 0.8 and  $len2/$l2 > 0.8)  or ( $lenx1/$l1 > 0.8 and  $lenx2/$l2 > 0.8)  ){
						### same alignment, decide which one should be kept
						print " overlapped ....\n" if $debug ;
						my $substitute = 0 ;
						if($as > $collectr{AS}){
							$substitute = 1;
						}elsif($as == $collectr{AS}){
							my @uniq1 = uniq(keys %{ $uniq_nodes{1} }, keys %{ $uniq_nodes{2} });
							my @uniq2 = uniq(keys %{ $collectr{R1} }, keys %{ $collectr{R2} });
							if(@uniq1 < @uniq2){
								$substitute = 1;
							}
						}
						if($substitute){
							print "discard previous. $rd\n" if $debug;
							print $fhd join "\t" , (@{$collectr{AR}}, "(DUP)") , "\n";
							$rdrm{$rd}{D} ++; 
							$stat_c{D} ++;
							if(defined $collectr{COM}){
								my ($n,$r) = @{$collectr{COM}} ; # [$nodes1[-1][0],$read];
								delete $complete{$n}{$r};
							}
							delete $collect{$rd}{$read_num} ;
						}else{
							print "current one discard... create DISC: $disc tag\n" if $debug ;
							$disc .= "(DUP)";
						}
						last;
					}else{
						print "not significantly overlapped..\n" if $debug;
					}
				}
			}
		}
		###### firt round quality check finished : reads alignment based
		# $disc have discard information. if not discard , then  connect, transform and output
		
		unless($disc) {
			## prepare	
			my $read_num = $tc  + $i;  
			print "\t=$rd $read_num collected\n" if $debug;	
			map { $nodesb{$_} = 1 } keys %{$uniq_nodes{1}};
			map { $nodesb{$_} = 1 } keys %{$uniq_nodes{2}};

			$collect{$rd}{$read_num}{R1} = $uniq_nodes{1};
			$collect{$rd}{$read_num}{R2} = $uniq_nodes{2};
			$collect{$rd}{$read_num}{AS} = $as;
			$collect{$rd}{$read_num}{LN} = $l1 + $l2; 
			$collect{$rd}{$read_num}{AR} = \@ar ;
			print "save line @ar to read $rd:$read_num \n" if $debug;
			my @nodes1 = @{ build($ar[6],\%graph) };
			my @nodes2 = @{ build($ar[7],\%graph) }; # it return an array of node annotation sub-array
			
			####
			my $read = "$rd:$read_num";
			my $score = $inf1[1] + $inf2[1];
			my $rlen  = $inf1[0] + $inf2[0];

			my ($ns,$ne) = ($ar[1],$ar[2]);

			if($ar[5] > 0 ){
				# under this situation, should output to keep file : no $disc and have overlap
				# merge nodes, make node string
				$collect{$rd}{$read_num}{HUG} = 1; # two reads hugged with each other
			}else{
				# have no overlap, and not sure about the middle situation.
				print "\t need to be completed $read_num\n" if $debug;
				$complete{$nodes1[-1][0]}{$read} = [ $nodes1[-1], $nodes2[0],  $rd,$read_num ];

				$collect{$rd}{$read_num}{COM} = [$nodes1[-1][0],$read]; # COM means this is a read need to be filled , it can link with %complete, it is important
			}
		}else{
			# have $disc defined, so this annotation must go to discard file
			print "\tDISCARD  $lines[$i]\t$disc\n" if $debug;
			print $fhd "$line_w\t$disc\n";
			$stat_c{D} ++;
			$rdrm{$rd}{D} ++ ;
		}
		#### update iteration related information
		$last_p = $ar[2];
	}
	######
	if(%collect){
		link_it (\%complete, \%collect, \%graph, \%nodesb, $tabix, $tabix_l,$tabix_g, $fhk, $fhd,\%stat_c  ,\ %rdrm);
	}
	foreach my $r (keys %rdrm ){
		my $rt = $rdrm{$r}{R};
		my $rd = $rdrm{$r}{D} // 0;
		#print "save mul:$r\t$rt\t$rd\n";
		if($rt - $rd >  1){
			print $fhr "$r\t$rt\t$rd\n";
		}
	}


	##
	$fhk -> flush();
	$fhd -> flush();
	$fhr -> flush();
	$fhk -> sync();
	$fhd -> sync();
	$fhr -> sync();

	print "\tbachof lines preprocess finished\n" if $debug;
	$stat_c{D} //= 0;
	$stat_c{O} //= 0;
	$stat_c{W} //= 0;
	$stat_c{R} //= 0;
	$stat_c{G} //= 0;

	print "child $$:\n\tLowQuality:$stat_c{D}\n\tOverlapRescue:$stat_c{O}\n\tReferRescue:$stat_c{R}\n\tWanderRescue:$stat_c{W}\n\tGap:$stat_c{G}\n";
	$pm -> finish(0,\%stat_c);
	#$pm -> finish(0);
}

sub link_it {
	#link_it (\%complete, \%collect, \%graph, \%nodesb, $tabix, $tabix_l,$tabix_g, $fhk, $fhd,\%stat_c  ,\ %rdrm);
    my( $complete, $collect, $graph, 
		$nodesb, 
		$tabix, $tabix_l, $tabix_g, 
		$fhk, $fhd ,
		$stat_c, $rdrm )  = @_;
    print "####### \n Entering link model \n########\n" if $debug;
	
	print "\n>>>link_it step1 : prepare linear coordinates \n" if $debug;
	
	## STEP1: collect associated position in linear references default is GRCh38
	my %g2l;
	my @nodes_rs = keys %$nodesb;
	#unless (@nodes_rs){
		#my @emptys_r =  keys %$collect; 
		#print "\t\tempty nodes in this link sets\n@emptys_r????\n";
		#die;
	#}
	my %nranges = nodes_ranges( @nodes_rs );
	
	foreach my $r ( sort keys %nranges ){
		my $iter = $tabix_g -> query($r);

	    print "\tquery $r..\n" if $debug ;
		while (my $aline = $iter -> next) {
			my @ar = split /\t/, $aline;
			if($$nodesb{$ar[2]}){
				my($chr,$cs,$ce) = $ar[3] =~ /#(.+):(\d+)-(\d+)/;
				$g2l{$ar[2]} = [$chr,$cs,$ce];
			}
		}
	}

	print "\n>>>link_it step2 ...\n" if $debug;
	## STEP2:  2.1, transform to linear positioin for each reads and save them into LP and FL tag of collect	
	my %reaches; # similar as %complete	
	foreach my $r ( sort keys %$collect ){
        foreach my $rn ( sort keys %{ $$collect{$r} } ){
			print "\n\tlink_it: test $r $rn \n" if $debug ; 
			print "\t@{$$collect{$r}{$rn}{AR}}\n" if $debug ;
			### find location in linear reference
			# VVVVVVVVVVV default is GRCh38
			my $collectr = $$collect{$r}{$rn};

			my ($flen_r, $lp ) = linear_coor($collect, $r,$rn, \%g2l, $tabix_l, \%reaches  );

			$$collect{$r}{$rn}{LP} = $lp;
			$$collect{$r}{$rn}{FL} = $flen_r;
			print "\tstep2:need to link $r $rn\n" if $debug;
		}
	}
			
		
	## STEP3: fill gaps 
	#
	#                        {R1} = %nodes
	#						 {R2} = %nodes
	#						 {LN} = toltal read length
	#=>$collect{$rd}{$rd_num}{AS} = $as
	#                        {AR} = \@ar the whole line
	#                        {HUG} = 1 when has overlapped path 
	#						 {COM} = [$lastnode_read1 , $rd:$read_num]
	#						 {LP}  =  overhang length and linear positions....
	#						 {FL}  = total fragment length ? -1 means not avalable
	#
	#
	#=>$complete{$lastnode_read1}{$rd:$read_num} = [ @read1_lastnode, @read2_firstnode, $readname, $read_num]
	#
	#=>%reaches{chr}{edge1}{edge2} = [ $path1, $path2, $sn, $en ,$r, $rn,  $dir, $trunked_length ]
	# my @meta = ($path1, $path2, $sn, $en, $r, $rn,$dir);

	overlap_rescue ($collect,$fhk,$stat_c);
	refer_rescue($complete, $collect, \%reaches, \%g2l, $tabix_l, $fhk, $stat_c);
	wander_rescue($complete, $collect, $graph, $tabix, $fhk, $stat_c);

	## STEP4: collect lost reads 
    foreach my $node ( keys %$complete ){
		foreach my $read ( keys %{ $$complete{$node} } ){
			my ($n1,$n2,$r,$rn) = @{ $$complete{$node}{$read} };
			my $collectr = $$collect{$r}{$rn};
			print "\t$read lost in linking ..\n" if $debug ;
			print $fhd join("\t", @{$$collectr{AR}} , "(GAP)") , "\n";
			$$stat_c{G} ++;
			$$rdrm{$r}{D} ++ ;
		}
	}
    print "Link Module finished\n############\n\n" if $debug;

}

sub overlap_rescue {
	my( $collect, $fhk, $stat_c )  = @_;
	foreach my $r ( sort keys %$collect ){
		foreach my $rn ( sort keys %{ $$collect{$r} } ){
			my $collectr = $$collect{$r}{$rn};

			if($$collectr{HUG}){
				my $flen_r = $$collectr{FL};
				my $lp = $$collectr{LP};
				print "\tstep3:overlap rescue $r $rn \n" if $debug ;
				my ($path_m, $path_nodes,$frag)  = overlapped_pair( $collectr );    
				format_out($path_m,$path_nodes,$frag, $$collectr{AS}, $r, $rn, $fhk,"$lp;$flen_r");
				
				delete $$collect{$r}{$rn};
				$$stat_c{O} ++;
			}
		}
	}
}

sub refer_rescue {
	my($complete,$collect,$reaches,$g2l,$tabix_l ,$fhk, $stat_c) = @_;

	print "\n>>>refer_rescue ...\n" if $debug;
	foreach my $chr ( keys %$reaches ){
		my %eles; # this is a subset of %$reaches	
		my($cs,$ce);

		foreach my $edge1 (sort {$a <=> $b} keys %{$$reaches{$chr}}  ){
			foreach my $edge2 (sort {$b <=> $a} keys %{$$reaches{$chr}{$edge1}}){
				my @paths = @{$$reaches{$chr}{$edge1}{$edge2}};
				foreach my $p (@paths){
					print "\trefer rescue $chr:$edge1-$edge2 @$p\n" if $debug;					
				}
				if($cs and $ce ){
					if ($edge1 > $ce + 5000 ){
						## extract the collected reaches
						print "range $chr:$cs-$ce collected ..\n" if $debug;
						refer_rescue_worker(\%eles, "$chr:$cs-$ce", $tabix_l, $collect, $complete, $fhk, $stat_c);

						undef %eles;
						my ($cs,$ce) = ($edge1,$edge2);
					}else{
						$ce = $edge2 if $edge2 > $ce;
						$eles{$edge1}{$edge2} = $$reaches{$chr}{$edge1}{$edge2};
					}
				}else{
					($cs,$ce) = ($edge1,$edge2);
					$eles{$edge1}{$edge2} = $$reaches{$chr}{$edge1}{$edge2};
				}
			}
		}
		print "range $chr:$cs-$ce collected ..\n" if $debug;
		refer_rescue_worker(\%eles, "$chr:$cs-$ce", $tabix_l, $collect, $complete, $fhk, $stat_c);
	}
}

sub refer_rescue_worker{
	my ($eles,$win,$tabix,$collect,$complete, $fhk, $stat_c ) = @_;
	
	print "\n>>>refer_rescue_worker for $win ...\n" if $debug;
	my $iter = $tabix -> query($win);
	
	my %accumulate;
	my @lines;
	my $lct = 0;
	while ( my $aline = $iter -> next ) {
		print "\n\tlinear to graph: $aline\n" if $debug;
		my @ar = split /\t/, $aline;
		my ($id,$dir,$len) =  $ar[3] =~ /(\d+)(<|>)(\d+)/;
		$dir = $dir eq ">" ? 1 : -1;
		push @lines, [$id,$dir,$len];
		
		# initiate accumulation
		print "\t<<1>>initiate accumulation...\n" if $debug;
		foreach my $edge1 ( sort {$a <=> $b} keys %$eles ){
			if($edge1 <= $ar[2] and $edge1 > $ar[1]){
				print "\tinitiate accumulation for edge1:$edge1\n" if $debug;
				$accumulate{$edge1} = $lct;
			}else{
				#print "\tedge1:$edge1 not in range $ar[1]-$ar[2]\n" if $debug;
			}
		}

		# checking accumulation
		print "\t<<2>>checking accumulated edges ...\n" if $debug;
		foreach my $edge1 ( sort  keys %accumulate ){
			print "\t.from edge1 $edge1\n" if $debug;
			foreach my $edge2 ( keys %{$$eles{$edge1}} ){
				print "\ttest aim: $edge1-$edge2\n" if $debug;
				if($edge2 <= $ar[1] or  $edge2 > $ar[2]){
					print "\tofftarget edge2:$edge2 $ar[1] $ar[2]\n" if $debug;
					next;
				}else{	
					print "\ton target edge2:$edge2 $ar[1] $ar[2]\n" if $debug;
				}
				my @gmetas = @{$$eles{$edge1}{$edge2}};
				my $starti = $accumulate{$edge1};
				
				my @path_fill  =  @lines[$starti+1   .. $lct-1 ];
				foreach my $gm (@gmetas){
					if($debug){
						print "\t\\\--->rescue $edge1-$edge2 :  @$gm\n";
					}
					my ($path1, $path2, $sn, $en, $r, $rn,$pdir,$rlen) = @$gm;
					if($pdir  == -1){
						@path_fill = reverse @path_fill;
						$_->[1] *= -1 for @path_fill;
					}
					my $ipath = join ";",map { "$$_[0],$$_[1],$$_[2],0,0,4" } @path_fill;
					my $fraglen = abs($edge2 - $edge1) + $rlen ; 	
						
					print "\tinner path: abs($edge2-$edge1) + $rlen. $ipath\n" if $debug;
					
					my @path_nodes =  (keys %{ $$collect{$r}{$rn}{R1} }, keys %{ $$collect{$r}{$rn}{R2} }, map { $_->[0] } @path_fill ) ;

					my @path_parts;
					if($path1){
						push @path_parts, $path1;
					}
					if($ipath ){
						push @path_parts,$ipath;
					}
					if($path2){
						push @path_parts, $path2;
					}

					my $mpath =  join ";", @path_parts;
					
					my $lp = $$collect{$r}{$rn}{LP};
					my $flen_r = $$collect{$r}{$rn}{FL};
					format_out( $mpath, \@path_nodes, $fraglen, $$collect{$r}{$rn}{AS}, $r, $rn , $fhk,"$lp;$flen_r") ;
			
					$$stat_c{R} ++;	
					my ($lastnode,$read ) = @{$$collect{$r}{$rn}{COM}} ;
					delete $$collect{$r}{$rn};
					delete $$complete{$lastnode}{$read};
				}
			}
		}
		$lct ++;
	}
}
	
sub wander_rescue {
	my($complete,$collect,$graph,$tabix,$fhk,$stat_c) = @_;

    print "\n@@@@@@@@@@@@@, WANDER, @@@@@@@@@ \n \titeration from small to large nodes to link gaps\n" if $debug;
    
	#################
	# Iteration of all nodes that are last one of read1
	# find the start and end nodes
	#
	# $sn = start_node
    foreach my $sn ( sort { $a <=> $b } keys %{$complete} ){
		
        ### checkin destination nodes
		# both hashs %back and %paths are like "Condor Heroes":
		# one used to record the last/curent nodes, one use to record all pathes.	
		#
		my %back; #  node -> path with this node as the last one
        my %paths;
		#_____________________
		#
		my $snd = $sn < 0 ? -1 : 1;
		my $snid = abs($sn);
		$back{$sn}{0} = [$sn];  # the last nodes point to path index in @paths
        $paths{0} = "$snid,$snd,$$graph{$snid}{L},0,0,0";

        #push @{$paths{$sn}},  "$sn,$$graph{abs($sn)}{L},0,0,0";
        #my @paths = ("$sn,$$graph{abs($sn)}{L},0,0,0");  # collectio of pathes === path recorder

        my %energy ; # the left length of one read. key are path index => reads => length
        my %des ; # have key as "endnode", value as read name

        print "\ninit process with node $sn  $$graph{ $snid }{L},0,0,0 \n" if $debug;
        # %back have all the leading nodes as keys.
        while( %back ){
            ########
            if($debug){
                my @iters = keys %back;
                print "nodes at end @iters\n";
            }
			#################################
			## >>> condition to break <<< ##
			### too many pathes, ignore....
            if(keys %paths > 20){
                print "too many pathes encountered, last..\n" if $debug;
                last;
            }
	
            foreach my $bn  ( keys %back ){  # $bn, base node, or the current last nodes of each path.
                my ($bn_i, $bn_d) = $bn > 0 ? ($bn,1) : ($bn * -1, -1); # $bn the current base node
                # begain to growth path
                if($debug){
                    print "\ngrowth node :  $bn $bn_i $bn_d \n";
                }
                # if a leading edge nodes are the edge of reads pair.  then check_in /
                if($$complete{$bn}){
                    print "register new tadpole \n" if $debug;
                    my $rc = dest_checkin ( $complete, $collect, $bn , \%back , \%energy, \%des  ); # completehash, base_node, des_hash, paths_hash, curent_paths
				}

                # extending graph according to requirement
                if ( ! defined $$graph{ $bn_i }  ){
                    print "$bn not in reads, need check graph???\n" if $debug;
                    extend_tabix($bn,$tabix,$graph, 30  ); # get surrounding nodes within range -20 +20
                }

				######### WORKING ON GRAPH HASH ###########
                # test if SNP bubble were identified. if it is a snp, then merge it as the same node group.
				# extending nodes , the representive one as the first key, then key C and D .  all are used the second level keys. values are depthes.
                my %ens ; 
                my @bn_neis = keys %{$$graph{$bn_i}{$bn_d} }; # two keys are nodeid and direction 
                if(@bn_neis == 1){
                    $ens{$bn_neis[0]}{C} = \@bn_neis;
                    $ens{$bn_neis[0]}{D} = $$graph{$bn_i}{$bn_d}{$bn_neis[0]}; # is means a weight of edge, depth of coverage
                }else{
                    ### might have bubbles
                    my %forwards;
                    foreach my $nn ( sort  { $$graph{$bn_i}{$bn_d}{$b} <=> $$graph{$bn_i}{$bn_d}{$a} }  @bn_neis){
                        my($nn_i,$nn_d) = $nn > 0 ? ($nn, 1) : (-1 * $nn, -1);
                        unless (defined $$graph{$nn_i}){
							next;
						}

						if( $$graph{$nn_i}{L} <= 10 ){ ### only consider SNPs, not INDELs
                            my $nnn = join ",", (sort keys %{$$graph{ $nn_i }{ $nn_d }});
                            print "$nn $nn_i $nn_d have length less than 10  => $nnn \n" if $debug;
                            $forwards{$nnn}{$nn} = $$graph{$bn_i}{$bn_d}{$nn} ;
                        }else{ ## node larger than 10, incorparate into path iteration
                            print "$nn have length > 1 $$graph{abs($nn)}{L}   \n" if $debug;
                            $ens{$nn}{C} = [$nn] ; # {$nn} = $$graph{$bn_i}{$bn_d}{$nn} ;
                            $ens{$nn}{D} = $$graph{$bn_i}{$bn_d}{$nn};
                        }
                    }

                    # one node bubble, the most common one
                    foreach my $f ( keys %forwards ){
                        my @mid_ns = sort keys %{$forwards{$f}}; # this @mid_ns is realted to $nn in %forwards hash
                        my @depth = map { $$graph{$bn_i}{$bn_d}{$_} } @mid_ns ;
                        my $ci = choose (@depth);

                        if($debug and @mid_ns > 1){
                            print "found bubble from $bn -> @mid_ns\n";
                            print "Depth for next nodes: @mid_ns  = @depth\n";
                            print " @mid_ns  leat to nn $f choose $ci \n";
                        }
                        map { $ens{$mid_ns[$ci]}{D}  += $forwards{$f}{$_} } @mid_ns ;
                        $ens{$mid_ns[$ci]}{C}  = \@mid_ns;
                    }
                }
                
				###### get idxes for a specific leading node
                ## iterating each previous paths
                #
                #------ 0 \                                    |
                #           \   $bn                            |
                #--------1 ---- N12   ->  N13, N14, N15        |  iter 1 for $bn
                #   @path_idx /              @ens_rep          |
                #---------2 /                                  |
                #
                #--------3  ---- N11    -> N16,N18             |   iter 2 for another $bn
                
				my @path_idxs = keys %{$back{$bn}}; # path index with this $bn at the end
                my @ens_rep = keys %ens;###  ======= this is as nodes should be extended into. more than 1 nodes, will increase number of pathes.

                if($debug){
                    if(@ens_rep){
                        print "next nodes: @ens_rep\n";
                    }else{
                        print "no nodes extendinge..\n";
                    }
					print "Pathes indexs: @path_idxs ended with $bn . extend to @ens_rep \n" if $debug;
                }

                foreach my $pi ( @path_idxs ){
                    unless(@ens_rep){ # dead end for this path
                        delete $paths{$pi};
                        next;
                    }

					# create new path index .
                    my $last = (sort {$a <=> $b} keys %paths )[-1] ;    #paths ; # last index of all paths
                    my @nidx = ( $pi, $last+1 .. $last + $#ens_rep ) ; # current idex and series of new idex in case of @ens have more then one element
                    if($debug){
                        print "extend or create new path idex : (@nidx) for new ends: ( @ens_rep ) the current lasti index  is $last\n";
                    }

                    # iterating next steps, each potential nodes will have a potential path index
                    for(my $j = 0; $j < @ens_rep ; $j ++){
                        ### collect infor about this extending node
                        # node id
                        my $en = $ens_rep[$j];
                        my ($en_i, $en_d) = $en > 0 ? ($en, 1) :  (-1 * $en, -1);
                        my %en_g = map {$_ => 1}  @{$ens{$en}{C}}; # nodes from the same bubble 
                        
						unless(exists $$graph{$en_i} ){
                            print "$en at the bobundary..\n" if $debug;
                            next;
                        }else{
							if($debug){
								my %n = %{$$graph{$en_i}};
								foreach my $k (keys %n) {
									if($k ne "L"){
										my @ik = keys %{$n{$k}};
										print "$k =>  @ik\n";
									}else{
										print "LENGH $n{$k}\n";
									}
								}
							}
						}
                        # type of nodes, "" form graph file 0 from other reads
                        my $type = "";  # empty means the node is recovered from graph, not sequenced
                        if($ens{$en}{D} > 0){
                            $type = 0; # 0 means, the extends nodes have overlap in graph.
                        }

                        # node length
                        my $en_l = $$graph{ $en_i }{ L };
						print "node:$en_i have length $en_l\n" if $debug;
						print "node:$en_i why no length\n" unless ($en_l);

                        # the path index the node will go into
                        my $path_i = shift @nidx; #$nidx[$j];

                        print "EXTENDING NODE: $en, lenth: $en_l, path index: $path_i\n" if $debug ;

                        print "Test if read are still valide if extending to $en\n" if $debug;
                        #delete $$graph{$bn}{C}{$en}; # delete seen path, avoid endless loop
                        ## each read is a tadpole, it is trying find its mom (dest_node). tadpole have energy to a maximum distance that is determing the the library size.
                        foreach my $r (sort  keys %energy ){
                            unless($energy{$r}{$pi}){
                                print "tadpole $r not on this path $pi\n" if $debug ;
                                next;
                            }
                            my $energy_status = $energy{$r}{$pi};
                            print "energy left  $r $energy_status  - $en_l \n" if $debug;
                            if(  $en_g{ $des{$r}{A} }  ) { # if a aim/dest node exists in extending node groups . 
								my $sn = $des{$r}{S}; ;
								my ($sn_i, $sn_d) = $sn > 0 ? ($sn,1) : (-1 * $sn, -1);
								
								my ($ipath) = $paths{$pi} =~ /(?:$sn_i,$sn_d,\d+,\d+,\d+,\S;)(.*)/;
                                $ipath = "" unless $ipath ;
                                my $frag = 0;
                                while($ipath =~ /\d+,(?:1|-1),(\d+),\d+,\d+,\d*;?/g){
                                    $frag += $1;
                                }

                                print "found mom  $r  source $des{$r}{S} => end $en <> $ipath <> $paths{$pi}\n" if $debug ;
                                my @path_parts;
                                my @ends;

								my @path_nodes;
                                @path_nodes = $ipath =~ /(?:;|^)(\d+)/g if $ipath ;
                                print "iner path nodes: @path_nodes <<< $ipath <<< $paths{$pi} \n" if $debug;
								
								my($node1r,$node2r,$rd,$read_num) = @{$$complete{ $des{$r}{S} }{$r}};	
								my $collectr = $$collect{$rd}{$read_num};

								my %r1 = %{$$collectr{R1}};
								my %r2 = %{$$collectr{R2}};
								
								my @nodes1 = sort { $r1{$a}{I} <=> $r1{$b}{I} } keys %r1;
								my @nodes2 = sort { $r2{$a}{I} <=> $r2{$b}{I} } keys %r2;
								
								#### first read nodes
                                for(my $i = 0; $i < @nodes1  ; $i ++){
                                    my $n1 = $nodes1[$i];
									push @path_nodes , $n1;
									my $ht ;
									if($i == $#nodes1 ){
										if($r1{$n1}{D} > 0){
											$ht = "$r1{$n1}{H},0";
										}elsif($r1{$n1}{D} < 0){
											$ht = "0,$r1{$n1}{T}";
										}
									}else{
										$ht = "$r1{$n1}{H},$r1{$n1}{T}";
									}

									push @path_parts, join ",", $n1,$r1{$n1}{D},$r1{$n1}{L},$ht,1;
									
									my $n_contrib = 0;
									if($i == $#nodes1 ){
										if($r1{$n1}{D} > 0){
											$n_contrib = $r1{$n1}{L} - $r1{$n1}{H};
										}elsif($r1{$n1}{D} < 0){
											$n_contrib = $r1{$n1}{L} - $r1{$n1}{T};
										}
									}else{
										$n_contrib = $r1{$n1}{L} - $r1{$n1}{H} - $r1{$n1}{T} ;
									}

									$frag += $n_contrib;;
                                }
            
								### inner path nodes
                                if ($ipath){
                                    push @path_parts, $ipath;
                                }
              

								#### second read nodes
								for(my $i = 0; $i < @nodes2  ; $i ++){
									my $n2 = $nodes2[$i];
									push @path_nodes , $n2;
									my $ht;
									if($i == 0 ){
										if($r2{$n2}{D} > 0){
											$ht = "0,$r2{$n2}{T}";
										}elsif($r2{$n2}{D} < 0){
											$ht = "$r2{$n2}{H},0";
										}
									}else{
										$ht = "$r2{$n2}{H},$r2{$n2}{T}";
									}

									push @path_parts, join ",", $n2,$r2{$n2}{D},$r2{$n2}{L},$ht ,2;
									
									my $n_contrib = 0;
									if($i == 0 ){
										if($r2{$n2}{D} > 0){
											$n_contrib = $r2{$n2}{L} - $r2{$n2}{T};
										}elsif($r2{$n2}{D} < 0){
											$n_contrib = $r2{$n2}{L} - $r2{$n2}{H};
										}
									}else{
										$n_contrib = $r2{$n2}{L} - $r2{$n2}{H} - $r2{$n2}{T} ;
									}
									$frag  += $r2{$n2}{L} - $r2{$n2}{H} - $r2{$n2}{T} ;
								}


                                if($debug){
                                    print "\tpaths: @path_parts\n";
                                    print "\tnodes: @path_nodes\n";                                                                                                                        
                                }
                                ## merge path
                                my $path_m = join ";", @path_parts;
										
								my $lp = $$collectr{LP};
								my $flen_r = $$collectr{FL};
								format_out($path_m,\@path_nodes, $frag, $$collectr{AS}, $rd, $read_num, $fhk, "$lp;$flen_r" );

								$$stat_c{W} ++;	
								
								delete $$collect{$rd}{$read_num};
								delete $$complete{ $des{$r}{S} }{$r};
								#delete $$complete{ $des{$r}{S} }{$r};
                                delete $des{$r};
                                delete $energy{$r};
                            }else{
                                my $new_energy = $energy_status - $en_l ;
                                print "energy updated $new_energy\n" if $debug;
                                if($new_energy < 0){
									delete $energy{$r}{$path_i};
								}else{
									$energy{$r}{$path_i} = $new_energy;
								}
                            }
                        }
                        $paths{$path_i} =  "$paths{$pi};$en_i,$en_d,$$graph{abs($en)}{L},0,0,$type"  ;
                        $back{$en}{$path_i} = 1;
                    }
                }
                #####
                delete $back{$bn};
            }
            
			unless (keys %des ){
                print "luckily, all gaps have be filled\n" if $debug ;
                last;
            }
        }
    }
}

sub dest_checkin {
	my ($complete, $collect, $bn, # input data
		$back, # backtrack of paths based on the last node 
		$energy, $des ) = @_;
	
	my $rc = 0;
	foreach my $r (keys %{$$complete{$bn}} ){
		my ($snr,$enr,$rd,$read_num )  = @{$$complete{$bn}{$r}};

		my @sn = @$snr;
		my @en = @$enr;

		my @lar = @{$$collect{$rd}{$read_num}{AR}};
		my $lensum = 0;
		for my $i ( -1, -2){
			my $l = (split /\|\|/, $lar[$i])[0];
			unless ($l =~ /^\d+$/){
				print "why @lar\n";
			}
			$lensum += $l;
		}

		foreach my $p ( keys %{ $$back{$bn} } ){
			# the bridge can not have length larger than $left
			my $left  = $frag_len - $lensum;
			if($bn > 0){
				$left -= $sn[3];
			}else{
				$left -= $sn[2];
			}
			if($en[0] > 0){
				$left -= $en[2];
			}else{
				$left -= $en[3];
			}
			####
			if($left > 0){
				$rc ++;
				$$des{$r}{A} = $en[0];
				$$des{$r}{S} = $bn;
				$$energy{$r}{$p} = $left ;
			}
		}
		$rc ++;
	}
	return ($rc);
}

sub nodes_ranges {
    my $debug = 0;
    my @nums = @_;

    @nums = sort { $a <=> $b } @nums; 
    my %ranges;
    my %includs;
    if(! defined $nums[0] ){
		return ;
	}
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

sub linear_coor{
	#convert_linear($collect, $r, $rn, \%g2l);
	my $prefix = "\t###LC###";
	print "\n$prefix enter linear_coor\n" if $debug;
	my($collect, $r, $rn ,$g2l,$tabix_l,$reaches ) = @_;

	my $collectr = $$collect{$r}{$rn};
	
	my $rlen = $$collectr{LN}; 
	my %r1 = %{$$collectr{R1}};
	my %r2 = %{$$collectr{R2}};
	
	my @nodes1 = sort { $r1{$a}{I} <=> $r1{$b}{I} } keys %r1;
	my @nodes2 = sort { $r2{$a}{I} <=> $r2{$b}{I} } keys %r2;

	#####################
	### PART 1: 
	#####################
	## >>>>>>> infer the linear position. 
	#project graph mapping to linear coordinate
	#
	my ($hc,$tc) = (0,0); # head and tail soft clipped length
	my ($lsc,$lec) = ("","");
	my ($ls,$le) = (0,0);
	## first read
	foreach my $n (@nodes1){
		if($$g2l{$n}){
			my($chr, $lhp,$ltp) = project_node_to_linear($$g2l{$n} , $r1{$n}); 
			print "$prefix\\->read1: $n projected to $chr : $lhp - $ltp \n" if $debug;
			$lsc = $chr;
			$ls = $r1{$n}{D}> 0 ? $lhp : $ltp;
			last;
		}else{
			$hc += $r1{$n}{L} - $r1{$n}{H} - $r1{$n}{T};
			if($hc > 30){
				last;
			}
		}
	}
	# second read	
	foreach my $n (reverse @nodes2){
		if($$g2l{$n}){
			my($chr, $lhp,$ltp) = project_node_to_linear($$g2l{$n} , $r2{$n}) ;
			print "$prefix\\->read2: $n projected to $chr : $lhp - $ltp \n" if $debug;
			$lec = $chr; 
			$le = $r2{$n}{D} > 0? $ltp : $lhp;
			last;
		}else{
			$tc += $r2{$n}{L} - $r2{$n}{H} - $r2{$n}{T};
			if($tc > 30){
				last;
			}
		}
	}
	
	my $rpos = "" ;
	my $rflen  = 0;
	print "$prefix:conditions:$lsc $lec $ls $le $frag_len - $hc - $tc\n" if $debug;
	if( $lsc and ($lsc eq $lec) and (abs($ls - $le) < ( $frag_len - $hc - $tc )) ){
		$rpos = "$lsc:$ls-$le";
		$rflen = abs($ls - $le) + 1 ;
	}
	#return ($rflen, $rpos );

	###############
	### PART 2: 
	################
	## collect the reach ends of two reads
	## use linear reference to fill gap
	## 
	if($rpos and ! $$collectr{HUG} ){
		print "\n$prefix>>collect reads reaches ...\n" if $debug;
		my ($sn, $en, 
			$path1, $path2,
			$edge1, $edge2 # real first  edge positon of linear
		);
		my ($sns,$ens)  = (0,0);
		
		for (my $i = $#nodes1; $i >= 0; $i--) {
			my $n = $nodes1[$i];
			if( $$g2l{$n} ){
				## accumulate path;
				my($chr, $lhp,$ltp) = project_node_to_linear($$g2l{$n} , $r1{$n});
				if($chr eq $lsc){
					$sn = $n;
					print "\t\t-->first node $n projected to $chr : $lhp - $ltp \n" if $debug;
					my $en1 = $nodes1[$i];
					if($i == 0){
						if($r1{$en1}{D} > 0){
							$path1 = "$en1,$r1{$en1}{D},$r1{$en1}{L},$r1{$en1}{H},0,1";
						}else{
							$path1 = "$en1,$r1{$en1}{D},$r1{$en1}{L},0,$r1{$en1}{T},1";
						}
					}else{
						$path1 = join ";", map{ "$_,$r1{$_}{D},$r1{$_}{L},$r1{$_}{H},$r1{$_}{T},1" } @nodes1[0..$i-1];
						$path1 .=  ";$en1,$r1{$en1}{D},$r1{$en1}{L},0,0,1";
					}

					$edge1 = $r1{$n}{D} > 0 ? $ltp : $lhp;
					last;
				}
			}else{
				$sns +=  $r1{$n}{L} - $r1{$n}{H} - $r1{$n}{T};
				if($sns > 2 ){
					last;
				}
			}
		}
		
		for (my $i = 0; $i < @nodes2 ; $i++ ) {
			my $n = $nodes2[$i];
			if( $$g2l{$n} ){
				my($chr, $lhp,$ltp) = project_node_to_linear($$g2l{$n},$r2{$n});
				if($chr eq $lsc){
					$en = $n;
					print "\t\t-->second node $n projected to $chr : $lhp - $ltp \n" if $debug;
					my $en2 = $nodes2[$i];
					if($i == $#nodes2){
						if($r2{$en2}{D} > 0){
							$path2 = "$en2,$r2{$en2}{D},$r2{$en2}{L},0,$r2{$en2}{T},2";
						}else{
							$path2 = "$en2,$r2{$en2}{D},$r2{$en2}{L},$r2{$en2}{H},0,2";
						}
					}else{
						$path2 = "$en2,$r2{$en2}{D},$r2{$en2}{L},0,0,2;";
						$path2 .= join ";", map{ "$_,$r2{$_}{D},$r2{$_}{L},$r2{$_}{H},$r2{$_}{T},2" }@nodes2[$i+1 .. $#nodes2]; 
					}
					$edge2 = $r2{$n}{D} > 0 ? $lhp : $ltp;
					last;
				}
			}else{
				$ens +=  $r2{$n}{L} - $r2{$n}{H} - $r2{$n}{T};
				if($ens > 2 ){
					last;
				}
			}
		}
		if($path1 and $path2 ){
			if(abs($edge1 - $edge2) > $rflen + 40 ){
				print "$prefix>>reach edges too far apart: $edge1 - $edge2 > $rflen + 40\n" if $debug;
			}elsif(($edge1 - $ls)	* ($edge2 - $le)  >= 0 ){
				print "$prefix>>reach edges not spanning the gap: ($edge1 - $ls)	* ($edge2 - $le)  >= 0 \n" if $debug;
			}elsif( ($edge1 - $edge2) * ($ls - $le)  <= 0 ){
				print "$prefix>>reach edges not in correct order: ($edge1 - $edge2) * ($ls - $le)  <= 0 \n" if $debug;
			}else{
				print "$prefix>>reach edges collected: $edge1 - $edge2 \n" if $debug;
				my $dir ;
				if($edge1 < $edge2){
					$dir = 1;
				}else{
					$dir = -1;
					($edge1,$edge2) = ($edge2,$edge1);
				}
				my @meta = ($path1	, $path2, $sn, $en, $r, $rn,$dir, $rlen - $sns - $ens );
				push @{$$reaches{$lsc}{$edge1}{$edge2}},\@meta;
			}
		}
	}else{
		print "$prefix>>hugged, no need to collect reaches\n" if $debug;
	}
	print "$prefix leave linear_coor\n\n" if $debug;
	return ($rflen, $rpos );
}
sub project_node_to_linear{
	my($lp,$rp) = @_;
	my @lpos = @{$lp} ; # always are for head and tail 
	my $ld = $lpos[1] < $lpos[2] ? 1 : -1;
	
	my %gpos = %$rp;
	my $nd = $gpos{D};
	
	my ($rhp,$rtp) ; # real linear position of node head and tail

	$rhp = $lpos[1] + $ld * ($gpos{H} + 1); # $rhp is 1 based position
	$rtp = $lpos[2] - $ld * $gpos{T}; # $rtp is 1 based position

	
	print "\tproject_linear: @lpos, l d : $ld, node d : $nd, head $gpos{H} tail $gpos{T}\n" if $debug;
	return ($lpos[0],$rhp,$rtp);
}

sub overlapped_pair{
	my($collectr) = @_;

	my @path_parts;
	my @ends;	
	my @path_nodes;
	my $frag  = 0;

	my %r1 = %{$$collectr{R1}};
	my %r2 = %{$$collectr{R2}};
	my @ar = @{$$collectr{AR}}; 
	
	my @nodes1 = sort { $r1{$a}{I} <=> $r1{$b}{I} } keys %r1;
	my @nodes2 = sort { $r2{$a}{I} <=> $r2{$b}{I} } keys %r2;
	
	for(my $j = 0; $j < @nodes1; $j ++){
		my $over = $ar[5] - (@nodes1 - $j );
		my $nid1 = $nodes1[$j];
		my @n1 = ($nid1, @{$r1{$nid1}}{qw(D L H T)} )  ;
		if($over < 0){
			push @path_parts , join ",", @n1, 1;
			push @path_nodes , $n1[0];
			$frag += $n1[2] - $n1[3] - $n1[4];
		}else{
			my $nid2 = $nodes2[$over];
			my @n2 = ($nid2, @{$r2{$nid2}}{qw(D L H T)} ); 

			my $h = min($n1[3],$n2[3]);
			my $t = min($n1[4],$n2[4]);
			
			my $mn = "$n1[0],$n1[1],$n1[2],$h,$t,3";
			push @path_parts, $mn;
			push @path_nodes,  $n1[0] ;
			$frag += $n1[2] - $h - $t;
		}
	}
	for (my $k = $ar[5]; $k < @nodes2; $k ++ ){
		my $nid2 = $nodes2[$k];
		my @n2 = ($nid2, @{$r2{$nid2}}{qw(D L H T)} ); 
		
		push @path_parts, join ",", @n2, 2 ;
		push @path_nodes,  $n2[0] ;
		$frag += $n2[2] - $n2[3] - $n2[4] ;
	}

	## merge path
	my $path_m = join ";", @path_parts;
	return ($path_m,\@path_nodes,$frag);
}

sub format_out {
	my($path_m, $path_nodes, $frag,$as , $r,$rn, $fh,$rstr) = @_;
	#print $fhk "$out_str\t$lp;$flen_r\n";
	my @chunks = @{ chunks ($path_nodes) };
	my @chunks_s = map { join ",", @$_  } @chunks;
	my $chunks_c = @chunks;

	for(my $m = 1; $m  <= $chunks_c; $m ++){
		my @rests = @chunks_s;
		splice(@rests, $m - 1 , 1);         # remove 1 element at index 2
		my $joined = join(';', @rests);    # join remaining elements
		$joined = $joined? $joined : ".";		
		print $fh "G\t$chunks[$m-1][0]\t$chunks[$m-1][1]\t$as\t$frag\t$path_m\t$r:$rn\t$samplename\t$m\t$chunks_c\t$joined\t$rstr\n";
	}
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

sub choose{
	my @size = @_;

	if(@size == 1){
		return 0;
	}

	my $total = sum(@size);
	if($total <= 0 ){
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
			$$graph{$cn}{-1}{$n} += 0;
		}
		while( $ar[4] =~ /(-?\d+)/g){
			my $n = $1 ;
			$$graph{$cn}{1}{$n} += 0;
		}
	}
}

###
sub extend_tabix {
	my($bn,$tabix,$graph,$step ) = @_;

	my $start  = abs($bn) - $step ;
	my $end = abs($bn) + $step ;

	my $win = "G:$start-$end";
	print "\textending base $bn window $win ...\n" if $debug;
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
		print "\tno known nodes meeted, extending failed...\n" if $debug;
		return 0;
	}
	foreach my $l (@lines){
		my @ar = split /\t/, $l;
		my $cn = $ar[2];
		$$graph{$cn}{L} = $ar[5];
		while( $ar[3] =~ /(-?\d+)/g){
			my $n = $1  * -1 ; #> 0 ?  "+$1" : $1 ;
			$$graph{$cn}{-1}{$n} += 0;
		}
		while( $ar[4] =~ /(-?\d+)/g){
			my $n = $1 ; #> 0 ? "+$1" : $1 ;
			$$graph{$cn}{1}{$n} += 0;
		}
	}
	return $found ;
}

sub build {	
	my($path, $graph ) = @_;
	
	if($debug){
		print "build read path to graph: $path\n";
	}

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
			my $fn1 = $nodes[$i][0];
			my $fn2 = $nodes[$i+1][0];

			my $m = $fn1 > 0 ? 1 : -1;
			$$graph{$fn1 * $m }{$m }{$fn2} ++;
			print "B:$fn1 * $m = $m = $fn2\n" if $debug;

			my $rn1 = $nodes_r[$i][0] * -1;
			my $rn2 = $nodes_r[$i+1][0] * -1 ;

			my $mr = $rn1 > 0 ? 1 : -1;
			$$graph{ $rn1 * $mr }{$mr}{$rn2} ++;
			print "B:$rn1 * $mr = $mr = $rn2\n" if $debug;
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
			push @chunks, [$startn, $n_l];
			$startn = $n - 1;
		}
	}
	if($startn){
		my $ln = abs($pn[-1]);
		push @chunks, [$startn,$ln];
	}
	return \@chunks;
}
######

sub cs {
	my ($str) = shift @_;
	my $m = 0;
	while($str =~ /:(\d+)/g ){
		$m += $1;
	}
	return $m;
}

sub set_overlap {
	my ($h1,$h2) = @_;

	my $len = 0;
	foreach my $elem (keys %$h1) {
		if (exists $$h2{$elem}){
			my $ol = overlap_length($$h1{$elem}{L},$$h1{$elem}{H}, $$h1{$elem}{T}, $$h2{$elem}{H}, $$h2{$elem}{T});
			$len += $ol;
		}
	}
	return $len;
	#return ($inter / keys %$h1, $inter / keys %$h2 )
}

sub overlap_length {
	my ($L, $a_start, $a_end, $b_start, $b_end) = @_;

	my $a1 = $a_start;
	my $a2 = $L - $a_end;

	my $b1 = $b_start;
	my $b2 = $L - $b_end;

	my $left  = $a1 > $b1 ? $a1 : $b1;
	my $right = $a2 < $b2 ? $a2 : $b2;

	return $right > $left ? $right - $left : 0;
}

sub parse_node{
	my $p = shift @_;
	my %nodes;

	my @ns = split /;/, $p;
	for (my $i = 0 ; $i < @ns; $i ++){
		my $na = $ns[$i];
		my @infs = split /,/, $na;
		my $d = $infs[0] > 0 ? 1 : -1;
		my $n = abs($infs[0]);
		$nodes{$n}{L} = $infs[1];
		$nodes{$n}{H} = $infs[2];
		$nodes{$n}{T} = $infs[3];
		$nodes{$n}{I} = $i;
		$nodes{$n}{D} = $d;
	}
	return \%nodes;
}



=head
my %mreads;
my %mreads_confirmed;
foreach my $kf ( @keeps){
	warn "reading and count $kf\n";
	open my $fh , $kf or die $!;
	while(<$fh>){
		chomp;
		my @ar = split /\t/;	
		my $rd = $ar[6];
		my ($rds,$rep,$num) = $rd =~ /^(.*):(\d+):(\d+)$/;
		if($rep == 0){
			$ar[6] = "$rds:0";
			print KEP join("\t",@ar),"\n";
		}else{
			if( defined $mreads_confirmed{$rds} ){
				$mreads_confirmed{$rds} ++ ;
				$ar[6] = "$rds:$mreads_confirmed{$rds}";
				print KEP join("\t",@ar),"\n";
			}elsif( defined $mreads{$rds} ){
				my @preves = @{$mreads{$rds}};
				
				$preves[6] = "$rds:1";
				print KEP join("\t",@preves),"\n";
				$ar[6] = "$rds:2";
				print KEP join("\t",@ar),"\n";
				
				$mreads_confirmed{$rds} = 2  ;
				delete $mreads{$rds} ;
			}else{
				$mreads{$rds}  =  \@ar;
			}
		}
	}
	close $fh;
}
warn "output multiple reads\n";
foreach my $rds ( keys %mreads ){
	my @ls = @{$mreads{$rds}};
	$ls[6] = "$rds:0";
	print KEP join("\t",@ls),"\n";
}
close KEP;

=cut


=head
	for my $n1 (keys %$complete ){
		foreach my $r ( keys %{$$complete{$n1}} ){
			my ($n1r,$n2r,$rd,$rn) = @{$$complete{$n1}{$r}};

			my %r1 = $$collect{$rd}{$rn}{R1};
			my %r2 = $$collect{$rd}{$rn}{R2};

			my @nodes1 = sort { $r1{$a}{I} <=> $r1{$b}{I} } keys %r1;
			my @nodes2 = sort { $r2{$b}{I} <=> $r2{$a}{I} } keys %r2;
	
			my ($hc,$tc) = (0,0); # head and tail soft clipped length
			my ($lsc,$ls,$lec,$le) = ('','');
			## first read
			foreach my $n (@nodes1){
				if($$g2l{$n}){
					my @lpos = @{$$g2l{$n}};
					my $ld = $lpos[1] < $lpos[2] ? 1 : -1;
					
					my $nd = $r1{$n}{D};
					
					my $rp ; # real position in the middle of node
					if($nd > 0){
						$rp = $lpos[1] + $ld * $r1{$n}{H};
					}else{
						$rp = $lpos[2] - $ld * $r1{$n}{T};
					}
					($lsc,$ls) = ($lpos[0],$rp);
					last;
				}else{
					$hc += $r1{$n}{L} - $r1{$n}{H} - $r1{$n}{T};
				}
			}
			# second read
			
			foreach my $n (@nodes2){
				if($$g2l{$n}){
					my @lpos = @{$$g2l{$n}};
					my $ld = $lpos[1] < $lpos[2] ? 1 : -1;
					
					my $nd = $r2{$n}{D};
					
					my $rp ; # real position in the middle of node
					if($nd > 0){
						$rp = $lpos[2] - $ld * $r1{$n}{T};
					}else{
						$rp = $lpos[1] + $ld * $r1{$n}{H};
					}
					($lec,$le )  = ($lpos[0],$rp);
					last;
				}else{
					$tc += $r1{$n}{L} - $r1{$n}{H} - $r1{$n}{T};
				}
			}
		}
	}



		my $reads_wgap = 0;
		foreach my $r (keys %{$$complete{$sn}}){	
			my($rd,$rn) = $r =~ s/(.+):(\d+)$//;
			if($filled{$rd}{$rn}){
				next;
			}else{
				$reads_wgap = 1;
				last;
			}
		}
		unless ($reads_wgap){
			print "all associated reads checked..\n";
			next;
		}




