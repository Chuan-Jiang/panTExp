#!/usr/bin/perl
use warnings; use strict;
use List::Util qw(min max);
use Storable  qw(dclone) ;
use Bit::Vector;
use Parallel::ForkManager;
use Getopt::Long;
use File::Temp qw( tempfile );
use IO::Handle;
use File::Basename;
use File::Path qw( make_path );
use Fcntl;
use Bio::DB::HTS::Tabix;

$| = 1;


print "Running: perl $0 @ARGV\n";
Getopt::Long::Configure("no_ignore_case");
####
my $input_f;
my $input_mf;
my $input_mb;
my $infor_f;
my $nodenei_f;
my $sampleinfo ;
my $batch_size = 100 ; # one batch have at least 10 genes
my $num_processes = 4; # number of processes to run in parallel
my $strand_type = "r";
my $pre ;
my $debug = 0 ;
my $help;
my $OLEN = 20 ;
my $orat = 0.8 ;
my $asam = 2;
my $dist_gmer = 2000;
my $nread_gmer = 2;
my $nsamp_gmer = 2;
my $hit_clus   = 100 ; # the cluster distance to merge expression alignemnt hits. used only for inter-annotation region
my $hps_c = 2; # required number of haps to support a merge
my $uniflg = 0;

#my $run0 = Bit::Vector->new($hit_clus);

GetOptions(
    "o|overlap=s"  => \$input_f,
    "m|missed=s"  => \$input_mf,
	"k|break=s"   => \$input_mb,
    "i|infor=s"   => \$infor_f,
	"N|nodenei=s" => \$nodenei_f,
	"p|prefix=s" => \$pre,
	"b|batch=s" => \$batch_size,
    "d|debug"     => \$debug,
    "n|num_processes=i" => \$num_processes, # Allow setting number of processes
    "s|sampleinfo=s" => \$sampleinfo,
	"t|strand=s" => \$strand_type, # strand type: f, r, or n
	"olen" => \$OLEN,
	"orat" => \$orat,
	"asam" => \$asam, # expressed freatures should exist in at least $asam samples
	"dist_gmer" => \$dist_gmer, # negtive value meaning disable merge
	"nread_gmer" => \$nread_gmer,
	"nsamp_gmer" => \$nsamp_gmer,
	"h|help" => \$help
) or die "Invalid options!\n";

# If --help is called or no arguments were passed
if ($help ) {
    print_help();
    exit;
}

# your script logic here...
sub print_help {
    print <<"EOF";
Usage: script.pl -i <input_file> -p <prefix> [options]

Required arguments:
    -i, --input <file>         Path to input file (default: STDIN if unspecified)
    -p, --prefix <string>      Prefix for output files
	-I, --infor <file>         pangenono infor about annotation
Optional arguments:
    -b, --batch <number>       Batch size of genes for processing (default: 10)
    -n, --num_processes <num>  Number of processes to run in parallel (default: 4)
    -d, --debug                Enable debug mode (default: off)
    -s, --sampleinfo <file>    Sample information file
    -t, --strand <type>        Strand type: f (forward), r (reverse), or n (none) (default: r)
    --olen <number>            Overlap length (default: 20)
    --orat <number>            ORAT threshold (default: 0.8)
    --asam <number>            Minimum samples expressing features (default: 2)
    --dist_gmer <number>       Distance for gMer merging (negative to disable) (default: 2000)
    --nread_gmer <number>      Minimum reads for gMer (default: 2)
    --nsamp_gmer <number>      Minimum samples for gMer (default: 2)
    -h, --help                 Show this help message and exit
EOF
}


###
# Extract directory part from prefix
my $pdir = dirname($pre);

# If directory doesn't exist, create it
unless (-d $pdir) {
    make_path($pdir) or die "Failed to create directory $pdir: $!";
}
#


################
my $use_fork  = 1;
if($num_processes == 1){
	$use_fork = 0;
}elsif($num_processes == 0){
	die "use at least 1 processor\n";
}else{
	print "Using $num_processes processes for parallel execution.\n";
}

warn "use fork: $use_fork \n";
##### determine input format
my $in_fh;
if ($input_f) {
    if ($input_f eq "-") {
        $in_fh = *STDIN;
    } elsif($input_f =~ /\.gz$/) {
        open my $fh, "zcat $input_f |" or die $!;
        $in_fh = $fh;
    } else {
		open my $fh, "<", $input_f or die "Cannot open $input_f: $!";
		$in_fh = $fh;
	}
} else {
	warn "reading input from pipe\n";
    $in_fh = *STDIN;
}

###############################
## creat sample file handles ##
my %samples_fhs;
open IN, $sampleinfo  or die "Cannot open $sampleinfo: $!";
while(<IN>){
	chomp;
	next if /^#/;
	my @ar = split /\t/;
	open my $fh, "> $pre.$ar[0].assign" or die "Cannot open $pre.$ar[0].assign : $!";
	$samples_fhs{$ar[0]} = $fh;
}
close IN;

###############################
#### configurations ###########
my $SD ;
if($strand_type eq "f"){
	$SD = 1;
}elsif($strand_type eq "r") {
	$SD = -1;
}else{
	$SD = 0;
}

warn "sequencing direction coefficient: $SD ($strand_type)\n";

#########################################
##### create temp files for child #######

my @tempfiles;
my %fh_assign;
my @tmps;
for my $i (0 .. $num_processes  ) {
	my ($fh, $filename) = tempfile("childassign_$i-XXXX", UNLINK => 1, DIR => ".");
    push @tmps, $filename;
	my %fls = ( fh => $fh, file => $filename);
	if($debug){
		my ($fhd,$filenamed) = tempfile("childdestiny_$i-XXXX", UNLINK => 1, DIR => ".") ;
		$fls{fhd} = $fhd;
		$fls{filed} = $filenamed;
		push @tmps, $filenamed;
	}
    push @tempfiles, \%fls ;
}

 # Catch Ctrl-C (SIGINT)
 #unless($debug){
 #   $SIG{INT} = sub {
 #       print "\nCaught Ctrl-C, cleaning up...\n";
 #       unlink @tmps  or warn "Couldn't delete : $!";
 #       exit 1;
 #   };
 #}


##################
##################
srand(1);
my $pm = Parallel::ForkManager->new($num_processes);
$pm->run_on_finish(sub {
    my ($pid, $exit_code, $ident, $exit_signal, $core_dump) = @_;

    if ($exit_signal) {
        warn "Child $pid (ident=$ident) died from signal $exit_signal\n";
    }
    elsif ($core_dump) {
        warn "Child $pid (ident=$ident) dumped core!\n";
    }
    elsif ($exit_code != 0) {
        warn "Child $pid (ident=$ident) exited with error code $exit_code\n";
    }
	print "file handle $ident released by $pid..\n";
	delete $fh_assign{$ident};
});


my $tabix_m = Bio::DB::HTS::Tabix->new(filename => $input_mf) ;
#my $tabix_nn = Bio::DB::HTS::Tabix->new(filename => $nodenei_f) ;


### read in brekpoit in non-annotaiton reads files
my @brks1;
my @brks2;
open BRK, $input_mb or die $!;
while(<BRK>){
	chomp;
	my @ar = split /\t/;
	push @brks1, $ar[0];
	push @brks2, $ar[1];
}
close BRK;


#############################
##### MAIN PROCESS ##########
##### read input flow==######
##############################
my $count = 0;
my $count_b = 0;
my @block = ( [] ) ;

my $node_last ;
my %gene_seen;
my $path_rep = 0;
#my %vec;
my $gap_c = 0 ;

my $line = <$in_fh>;
while( defined $line ){
	chomp ( $line );
	my $next_line = <$in_fh>;

	
	$count ++;
	$count_b ++;

	my @ar = split /\t/, $line;
	#my @map = @ar[0..11];
	#my @bed = @ar[14..$#ar -1];

	my @genes  = split /,/, $ar[-1] ;

	# have seen this genes before?
	my $old_gene = grep { exists $gene_seen{$_} } @genes;
	my $gc = scalar keys %gene_seen;

	################################################################
	# # ====> first condition
	# the decision to issue a process_lines subfunction are based on these conditions:
	# 1, encounter genes never seen before.
	# 2, the current nodes is not within the previous node id range
	# 3, if it is the last line, also need to run process_line
	#####
	# ======> second condition 
	# the previous condition only observed in reads have overlaps with annotaion. 
	# if the three conditions matched, then, i should checked reads that has no overlaps with any annotation. ATTENSYION< here: only
	# no-overlapped reads that in the downstream were check. up-stream should be checked in process_lines().
	##
	# so $tabix_m -> query($win) command is executaed to backtrack the reads before the current nodes. for each island, only down stream was checked.
	#
	# 1, if there is a  node id gap  between two consective lines.
	# 2, save node coverage to a hash %vec. it is binary based . accumlation coverage of nodes. with  aim to find gaps in long read
	# 3, 
	#
	# ===> decisions
	#	 if gap encoutheed, then the node id upstream of gap will be save to %block;
	#  
	#  this is a block have four islands:
	#	  [----?????AAAAAAA++++---????AAAA???AAA++++-----????AAA?AAA+++-----??AAA??AAA+++-----]
	#           ^            ^    ^                 ^      ^           ^     ^           ^
	#           - : gap
	#           ? : need to be confirmed in process_lines. the unknown region might covered by no-annotation reads
	#           + : covered by read lines don't have overlaps with annotation
	#           A : covered by read lines have overlaps with annotation
	#           ^ start node and end note for each island
	#  
	#	   for each elements in @block is a reference for an array:
	#	    ---[start_node, Readline1-A, Readline2-+, ...., Readline3-A, linei-+, end_node]---- 
	#
	# $gap_c have the count how many times of gaps encoutered, that equals the current number of islands in block
	# if number of island is larger than a threshold, then initatiate process
	#
	#
	# ===> collect information
	# prepare several blocks, each block is sperate island separate by gaps. several island makes block for submmision
	#
	#
	
	#====> first condition
	if( ( (!$old_gene )  and $node_last and $ar[1]  >= $node_last ) or !$next_line  ){
		#my $end = $next_line ?  $ar[1] + 1 : ""; 
		my $gap = 0;
	
		my $end = $ar[1] + 1;
		if($end < $node_last){
			$end = $node_last;
		}
		my $tail = $node_last;	
		my $fi = first_index_gt($node_last, \@brks1);
	
		if($fi){
			if( $brks1[$fi]  < $end ){
				$gap = 1;
				$tail = $brks1[$fi];
			}elsif( $brks2[$fi - 1] > $node_last ) {
				$gap = 1;
				$tail = $brks2[$fi - 1] > $end ? $end : $brks2[$fi - 1];
			}
		}


		### -=====> decidions, 
		if( $gap == 1 or !$next_line ){
			$gap_c ++;
			#print "gap found, save it to block tail $tail\n" ;
			push @{ $block[-1] }, $tail;

			# test if inita a worker......
			if( $count_b > 1000 ){
				## assign worker
				print "\t--->collect and send a block at $count $count_b | $gc $gap_c \n\n";
				process_lines(\@block);
				undef @block;
				undef %gene_seen;
				#undef %vec;
				$count_b = 0;
				$gap_c = 0;
			}
			push @block, [ ];
		}else{
			print "no gap found ..\n" if $debug ;
		}
	}
	######################################
	# =====> collect information
	# add start node for block
	unless( @{ $block[-1] } ){
		push @{ $block[-1] } , $ar[1] ++;
	}
	# update path_vec
	# path_vec($ar[5],\%vec);
	# save lines
	push @{ $block[-1] }, $line;
	$node_last = $ar[2];
	map { $gene_seen{$_} ++ } @genes;
	$line = $next_line;
	
	if(!$next_line){
		push @{ $block[-1] }, $ar[2];
	}
}

#############################################
# Process any remaining lines after the loop
if(@block){
	process_lines( \@block );
}

if($use_fork){
	print "Parent is waiting childs\n";
	$pm->wait_all_children; # Wait for all child processes to finish
	print "All child processes finished.\n";
}

######## CLEANUP #######
### LAST combine outputs ###########
open my $fh_destiny , ">", "$pre.destiny.tsv" or die "Cannot open $pre.destiny.tsv: $!" if $debug ;
foreach my $tf (@tempfiles) {
	#my $fh = $tf -> {fh};
	#my $fhd = $tf -> {fhd};
	#$fh -> flush;
	#$fhd -> flush;

	if($debug){
		open my $ind, "<", $tf->{filed} or die "Can't open $tf->{file}: $!";
		print "--- Contents of $tf->{filed} ---\n";
		while (<$ind>) {
			print $fh_destiny $_; # write to destiny file;
		}
		close $ind;
		unlink $tf->{filed}; # delete temp file
	}

    open my $in, '<', $tf->{file} or die "Can't open $tf->{file}: $!";
    print "--- Contents of $tf->{file} ---\n";
    while (<$in>) {
		chomp;
		my @ar = split /\t/;
		print "=:@ar\n" if $debug;
		my $sample = $ar[-1];
		my $ofh = $samples_fhs{$sample};
		die "$sample $tf->{file} $_ \n" unless $ofh;
		#warn "$ofh" ,  join("\t", @ar[0..$#ar]), "\n";
		print $ofh join("\t", @ar[0..$#ar]), "\n"; # write to sample file
    }
	close $in;
	unlink $tf->{file}; # delete temp file
}

#################################################
################ SUB FUNCTIONS ##################
#################################################
sub process_lines {
	my($blkr) = @_;
	my $block = dclone($blkr);
	undef @$blkr ;

	#### statistics

	### initiate working horse subprocess: single or parallel:
	my $fhi;
	if($use_fork){
		($fhi ) = getfh(\@tempfiles,\%fh_assign) ;
		my $pid = $pm->start($fhi);

		if($pid){
			return;
		}
		my $ppid = getppid();
		print "\t---CHILD $$ parent $ppid : get filehandle $fhi with Count $count  count_b: $count_b block_gaps: $gap_c \n";
	}else{
		$fhi = 0;
		print "\t---SINGLE  $$ : get filehandle $fhi with Count $count  count_b: $count_b block_gaps: $gap_c \n";
	}

	my $fh = $tempfiles[$fhi]->{fh};
	my $fhd = $tempfiles[$fhi]->{fhd} if $debug ;

	my $tabix_i = Bio::DB::HTS::Tabix -> new(filename => $infor_f) ;
	my $tabix_m = Bio::DB::HTS::Tabix -> new(filename => $input_mf) ;

	# iteration each island:
	foreach my $bl  (@$block){
		my @lines = @$bl;

		my $head = shift @lines;
		my $tail = pop @lines;

		#print "\nChild $$ is working with a gene block... node range: $head $tail\n"  ;

		##### annotations   ###
		my %node_vec; #### big ....
		my %anno; # {gene}{node} = [direct * $SD , start,end]
		my %align_inf; # {read} : ?????{V} = [direct,star,end,type} {S} = AS  {L} = frag length
		my %gene_inf; #while samples have this gene;=
		my %destiny ;
		my %uni_reads;
		my %vecf; # for path_vec

		#
		# foreach island, all reads mapping and annotation were processed, then extending to upstream boundary.
		#
		my %highcov ; 
		foreach my $lr ( @lines  ){
			my @ar = split /\t/, $lr;
			
			### determing if it is high repeat reads
			my $high_rep = determine_repeat($ar[4],$ar[5],$ar[7],\%highcov);

			if($high_rep){
				print "highcov read $lr\n" if $debug ;
				next;
			}

			my ($readm,$sample) = @ar[6,7] ;
			my ($read,$mapn) = $readm=~ /(.+):(.+)$/;
		
			## calcualte the real multi idx.
			my ($t) = $mapn =~ /^(\d+)/ ;
			my $r = eval ($mapn);
			$mapn = "$t-$r";

			my @genes = split /,/, $ar[-1];
			my $path = $ar[5];

			if($debug){
				print "\nstep1:NewLine:\t>>@ar\n" ;
				print "\tread:$readm=$read  +   , $sample\n";
			}

			# ===> STEP1: GENE Annotation from bed
			#                  {L}  = length of gene
			# $anno{gene name} {N}{node_id} = [ direction, start, end]
			#                  {H}{chromose_linear} = [ start, end, direction]
			#
			#  if a node is within a annotion, then the gap of reads coverage will be ignored.
			#

			print "step1: gene annotation parse\n" if $debug ;
			for my $gene (@genes){
				if ( ! defined $anno{$gene} ){
					print "\tfound new genes, parse...\n" if $debug;
					my %nodes;
					my %haps;

					my ($g,$n,$p) = split /\|/, $gene;
					my $ps = $p -1 ;
					my $iter = $tabix_i -> query("$g:$ps-$p");
					my $gl;
					my $gene_found = 0;
					while (my $aline = $iter -> next) {
						my @ar = split /\t/, $aline ;
						if($ar[3] eq $gene){
							print "\tannotation extract $aline\n" if $debug;
							map{ $nodes{$_} = 1 } split /;/, $ar[5];
							map{ $haps{$_} = 1 } split /;/, $ar[4];
							$gl = $ar[6];
							$gene_found = 1;
						}elsif($ar[3] ne $gene and $gene_found){
							last;
						}
					}

					$anno{$gene}{L} = $gl;
					#print "N:@nodes\n" if $debug;
					foreach my $n ( sort keys %nodes ){
						my($id,$l,$s,$e,$t) =  split /,/, $n ; #=~ /(.)(\d+):(\d+),(\d+),(\d+)/;
						my $d;
						if ( $id =~ /^(=|\+|-)(\d+)/ ){
							if($1 eq "+"){
								$d = 1;
							}elsif($1 eq "-"){
								$d = -1;
							}else{
								$d = 0; # no direction
							}
							$id = $2;
							if($s == 0 and $e == 0){
								$vecf{$id}{F} = 1;
							}else{
								reg_vec($id,$s,$l-$e-1, $l,\%vecf);
							}
						}else{
							die  "Bad node id $id in $gene, skip this node\n" if $debug;
						}
						#$d = 0 if ($SD == 0); # adjust direction
						$anno{$gene}{N}{$id} = [$d, $s, $l - $e - 1];
					}
					foreach my $h ( sort keys %haps ){
						my($chr,$start,$end,$dir) = $h =~ /^(.+):(\d+)-(\d+)#(-|\+)/;
						$anno{$gene}{H}{$chr} = [$start,$end,$dir];
					}
				}
			}

			## ===> STEP 2 READ Alignmentes
			print "step2: get_vecs for this read: $read  $path\n" if $debug ;

			#my($mline, $uni_reads , $node_vec, $align_inf) = @_;
			record_align($lr,\%uni_reads,\%node_vec,\%align_inf,0);


			## ===> STEP 3 OVERLAP calculation
			print "step3: overlap calculation .. $read##$mapn##$sample\n" if $debug ;
			foreach my $gene (@genes){
				print "\tinitall cal_over subprocessss $read##$mapn##$sample\n" if $debug;
				my ($nrange) = cal_overs($gene, \%anno, $node_vec{R}{"$read##$mapn##$sample"} ,\%gene_inf ,  $align_inf{"$read##$mapn##$sample"}{L} ); # if have overlap, then update %overs and return 1. else do nothing, return 0
				if( $nrange) {
					$destiny{ "$read##$mapn##$sample"  }{1}  = 1 if $debug ;
					print "\tstep3 -> have overlap with $read##$mapn##$sample\n" if $debug;
					$gene_inf{$gene}{S}{$sample} ++;
					$gene_inf{$gene}{R}{"$read##$mapn##$sample"}  = 1;
					$gene_inf{$gene}{B} = $nrange; # save ranges for this read
					#$gene_inf{$gene}{L} = ($gene_inf{$gene}{L} // $anno{$gene}{L})   +  $elen;
				}else{
					if($debug){
						$destiny{ "$read##$mapn##$sample"  }{1}  = 0  ;
						print "\tstep3 -> no overlap with $read##$mapn##$sample\n" ;
					}
				}

			}
		}

        #####################################################################
		#### all known information has been processed, then finished them.###
        ####################################################################

		my @anno_nodes = keys %vecf;
		my @genes_with_hits;
		my @hits_on_genes;
		foreach my $g (keys %gene_inf){
			if(defined $gene_inf{$g}{R}){
				my %rs = %{ $gene_inf{$g}{R} };
				if( keys %rs ){
					push @genes_with_hits, $g;
					push @hits_on_genes, keys %rs;
				}
			}
		}

		##############################################################
		######## beging to expand into non-overlapped reads  #########
		## information will be writen into node_vec and algin_inf 
		unless(@genes_with_hits){
			print "no genes with hits, skip this block\n" if $debug;
			next;
		}else{
			## considering only uniquly mapped????
			print "continue expanding this block with genes @genes_with_hits\n" if $debug;
		}
		# expand mismated reads
		my $win = "G:$head-$tail";
		#print "\t->retriving downstream mismatch: $win\n" ;
		my $iter_i = $tabix_m -> query($win);
		my $mnode_last = $head;
		while(my $mline = $iter_i -> next){
			my @ar = split /\t/, $mline;
			
			my $high_rep = determine_repeat($ar[4],$ar[5],$ar[7],\%highcov);
			if($high_rep){
				next;
			}

			path_vec($ar[5],\%vecf);

			if( $mnode_last <=   $ar[1] ){
				last;
			}else{
				my @vks = sort { $a  <=> $b } (grep { $_  <=  $ar[1] } keys %{ $vecf{B}} );
				if(@vks){
					$mnode_last = $vks[0];
					last;
				}
			}
			$mnode_last = $ar[2];
			record_align($mline,\%uni_reads,\%node_vec,\%align_inf,0);
		}

		# backward retrieve mismated reads  <---------+++++++++++++++++++++++++---------|
		#print "\t->retrieving backtrack upstream mismatch ..\n" ;
		my $loop = 1;
		my %vec;
		while(1){
			my $ls = $head - $loop * 10 + 1;
			my $le = $head - ($loop - 1) * 10;
			if($le <= 0){
				last;
			}

			$ls  = $ls < 0 ? 0 : $ls;

			my $node_last = $le;
			my $uiter = $tabix_m -> query("G:$ls-$le");
			my @ulines = ($uiter -> next );
			print "\tupstream window G:$ls-$le\n" if $debug;
			my $last;
			foreach my $uline( reverse @ulines ){
				my @ar = split /\t/, $uline;
				
				my $high_rep = determine_repeat($ar[4],$ar[5],$ar[7],\%highcov);
				if($high_rep ){
					next;
				}

				path_vec($ar[5],\%vecf);

				if( $node_last >  $ar[2] ){
					$last = 1;
					last;
				}else{
					my @vks = sort {$b <=> $a} (grep { $_ >  $ar[2] } keys %{$vec{B}} );
					if(@vks){
						$last = 1;
						$node_last = $vks[0];
						last;
					}
				}
				$node_last = $ar[1] + 1;
				record_align($uline, \%uni_reads,\%node_vec,\%align_inf,0);
			}

			if($last){
				last;
			}else{
				$loop ++;
			}
		}

		#print "Child $$ is expanding  overlaps in block ... node range: $head $tail \n";

		my $collector = expand_overlap(\%node_vec ,\%anno,\%align_inf, \%gene_inf,\%destiny);
		my $collectorgenes = scalar keys %$collector;
		#print "Child $$ is clusting genes $collectorgenes in block .... ... node range: $head $tail\n";
		cluster_output($collector , \%anno, \%align_inf, $fh);
		########### summary  #####
		if($debug){
			foreach my $r (keys %destiny){
				my @steps;
				foreach my $step (1..3){
					if(exists $destiny{$r}{$step}){
						push @steps, "$step:$destiny{$r}{$step}";
					}else{
						push @steps, "$step:NA"; # if not exists, then set to 0
					}
				}
				print $fhd "DESTINATION:$r\t", join("\t", @steps), "\n";
			}
		}
		############ end of summary
	}
	$fhd -> flush if $debug ;
	$fh -> flush;
	if($use_fork){
		$pm -> finish(0);
	}
}

########################
# subfunctuinos

sub determine_repeat{
	my($len, $path,$sam,$highcov) = @_;
	
	#print "highcovB $path\n";
	if ($path =~ /;/) {
		my @parts = split /;/, $path;
		$path = join ";", $parts[0], $parts[-1] if @parts > 2;
	}
	
	#print "highcovA $path\n";

	if(exists $$highcov{$path} ){
		if(exists $$highcov{$path}{$len.$sam} and $$highcov{$path}{$len.$sam} > 5 ){
			return 1;
		}
	}else{
		undef %$highcov ;
	}
	$$highcov{$path}{$len.$sam} ++;
	return 0;
}


sub path_vec{
	my($p,$vec) = @_;
	my @ns = split /;/, $p;
	foreach my $n ( @ns ){
		my @nis = split /,/, $n;
		my $nid = abs($nis[0]);

		#next if $nis[1] < $hit_clus;
		next if ($nis[2] + $nis[3] == $nis[1] );
		next if $$vec{$nid}{F};

		if($nis[2] == 0 and $nis[3] == 0){
			$$vec{$nid}{F} = 1;
			delete $$vec{$nid}{B};
			next;
		}

		unless(defined $$vec{ $nid }{B} ){
			$$vec{ $nid }{B} = Bit::Vector->new($nis[1]);
		}

		my $fs = $nis[2] ;
		my $fe = $nis[1] - 1 - $nis[3];
		$$vec{ $nid }{B} -> Interval_Fill( $fs  , $fe );
		if( $$vec{$nid}{B} -> is_full() ) {
			$$vec{$nid}{F} = 1;
			delete $$vec{$nid}{B};
		}
	}
}
sub reg_vec{
	my($n,$s,$e,$l,$vec) = @_;
	
	if($$vec{$n}{F}){
		return ;
	}

	unless(defined $$vec{$n}{B}){
		$$vec{$n}{B} = Bit::Vector -> new($l);
	}
	$$vec{$n}{B} -> Interval_Fill($s,$e);
	if( $$vec{$n}{B} -> is_full()){
		$$vec{$n}{F} = 1;
		delete $$vec{$n}{B};
	}
}

### save reads to node_vec and align_inf
sub record_align{
	####
	# requred inpt:
	#  > alignment line
	#  > %uni_reads
	#  > %node_vec;
	#  > $align_inf
	#
	# for each line:
	# 1, check if a existing similar read exists?
	#	- if same path seen before. then %uni_reads are used to extract the original unique readID.
	#	- based on the UniqueReadID, te node vec informatio can be extracted then;
	#		
	#		%node_vec{R}{readID##Mapn##Sample} = \%nv;
	#                {N}{$n}{readID##Mapn##Sample} = $$node_vec{N}{$n}{$uread};
	#    Or:
	#		creat from scratch
	# 2, save node_vec
	# 3, save align_inf
	#	{read##Mapn##Sample}{S} = $as   : alignment score
	#	                    {L} = $frag : Length of frag

	my($mline, $uni_reads , $node_vec, $align_inf,$uniflag) = @_;

	my @ar = split /\t/,$mline;
	my $as = $ar[3];
	my $frag = $ar[4];
	my $path = $ar[5];

	my ($readm,$sample) = @ar[6,7] ;
 	my ($read,$mapn) = $readm=~ /(.+):(.+)$/;

	my ($t) = $mapn =~ /^(\d+)/ ;
	if($uniflg){
		if($t > 1){
			return;
		}
	}

	my $r = eval ($mapn);
	$mapn = "$t-$r";

	print "\tget_vecs for this read: $read  $path\n" if $debug;
	# can make a univeral graph, reduce memory .
	if( $$uni_reads{$path} ){
		my $uread = $$uni_reads{$path};

		if($debug){
			print "\tHave this $path before from $uread, reuse ...\n";
		}
		#print "WHYSHY @alls \n" unless $node_vec{R}{$uread};
		my $nv = $$node_vec{R}{$uread};
		$$node_vec{R}{"$read##$mapn##$sample"} = $nv;
	 	foreach my $n (keys %$nv ){
		 	print "\t\tnode $n linked $$node_vec{N}{$n}{$uread}\n" if $debug;
		 	$$node_vec{N}{$n}{"$read##$mapn##$sample"} = $$node_vec{N}{$n}{$uread};
		}
	}else{
		print "\tParse new path $path  ...\n" if $debug;
		get_vecs ( $path, "$read##$mapn##$sample", $node_vec );
		$$uni_reads{$path} = "$read##$mapn##$sample" ;
		print "\tnew path recoreded .$path saved $read##$mapn##$sample\n" if $debug;
	}

	$$align_inf{"$read##$mapn##$sample"}{S} = $as;
	$$align_inf{"$read##$mapn##$sample"}{L} = $frag;

}

sub genes_similar {
    my ($gene1, $gene2, $anno, $range, $share) = @_;
    $range //= 1000;  # default range 1000bp

    # Get chromosomes shared by both genes
	#my @share_haps;
	my $ct = 0;
	for my $chr (keys %{ $anno->{$gene1}{H} }) {
        next unless exists $anno->{$gene2}{H}{$chr};

        my $haps1 = $anno->{$gene1}{H}{$chr};
        my $haps2 = $anno->{$gene2}{H}{$chr};
		my ($start1,$end1,$dir1) = @$haps1;
		my ($start2,$end2,$dir2) = @$haps2;
        # Check same direction
		next unless $dir1 eq $dir2;
        my ($dis,$totl) =  intervals_within_range($start1, $end1, $start2, $end2);

		if($dis < $range){
			push @$share, $chr;  # true, similar
			#push @hdis,  $dis;
			$ct ++;
		}
    }
    return ( $ct ); # no matching intervals found
}

sub intervals_within_range {
    my ($s1, $e1, $s2, $e2) = @_;
    # Calculate minimal distance between intervals (overlap is 0 distance)
    if ($e1 < $s2) {
        return ( $s2 - $e1, $e2 - $s1);
    } elsif ($e2 < $s1) {
        return ( $s1 - $e2, $e1 - $s2);
    } else {
		my @ar = sort ($s1,$s2,$e1,$e2);
		my $tl = $ar[-1] - $ar[0];
        return (0, $tl) ;  # intervals overlap
    }
}

sub median {
    my @data = @_;
    return undef unless @data;  # handle empty array

    # sort numerically
    @data = sort { $a <=> $b } @data;

    my $count = @data;
    if($count == 1){
		return $data[0];
	}elsif ($count % 2) {
        # odd number of elements → return the middle one
        return $data[ int($count / 2) ];
    } else {
        # even number of elements → average the two middle ones
        return ($data[$count/2 - 1] + $data[$count/2]) / 2;
    }
}


#####

sub get_vecs{ # for a read, get all nodes and their vectors
	my ( $path , $read, $node_vec ) = @_;
	my @nodes = split /;/, $path;
	#my @nodes = $path =~ /(<|>)(\d+)/g ;
	#@nodes = map{ [split /,/]  } @nodes;
	for( my $i = 0; $i < @nodes; $i ++ ){
		my($id,$l,$s,$e,$t) =  split /,/, $nodes[$i]  ;
		########## overall linkage with other nodes;
		my $d = $id > 0	 ? 1 : -1; # direction
		#$d = $d * $SD; # adjust direction
		$id = abs($id);
		if($s + $e >= $l ){
			print "Bad node, $id  [$s, $e] for vec of size " . $l . "from $read\n" if $debug;
			next;
		}
		## transcorm $t to score
		my $score = 0 ;
		if ($t ne "" ){
			if($t > 0 ){
				$score = 2;
			}elsif($t == 0 ){
				$score = 1;
			}
		}
		my $inf  = [ $d, $s, $l - $e - 1, $score * ($l - $s - $e)];
		$$node_vec{N}{$id}{$read} = $inf;  # save vector for this read
		$$node_vec{R}{$read}{$id} = $inf;  # save vector for this read
	}
}

###############
sub cal_overs  {  # calculate overlaps for a gene, return ranges if have enough overlaps
	my($gene, $anno, $vecs, $gene_inf, $fragl) = @_;
	my $len_m = 0;
	my $len_t = 0;
	print "\n\t\tcal_overs:subfunction: cal_overs for gene: $gene\n" if $debug;
	my %ranges  = ();
	%ranges = %{ dclone ($gene_inf -> {$gene} -> {B})  } if $gene_inf->{$gene}{B}; # get ranges if have

	my %update_materials;
	foreach my $id ( keys %$vecs ){
		my $rvec = $$vecs{$id} ;
		if( $$anno{$gene}{N}{ $id } ){
			print "\t\tcal_overs:gene: $gene, id: $id, vec: @$rvec\n" if $debug;

			my $avec = $$anno{$gene}{N}{$id};
			my $olen = 0;
			if( $$rvec[0]  *  $$avec[0] * $SD >= 0 ){ # if the direction is different, then skip
				my $start = $$rvec[1] > $$avec[1] ? $$rvec[1] : $$avec[1];
				my $end = $$rvec[2] < $$avec[2] ? $$rvec[2] : $$avec[2];
				$olen = ($end - $start + 1) > 0 ? ($end - $start + 1) : 0; # overlap length

				if($debug){
					print "\t\tcal_overs:id: $id. olen $olen  $gene @$avec @$rvec\n";
				}
				$len_m +=  $olen ;
			}else{
				print "\t\tcal_overs:direction mismatch, skip $id, $gene, $$rvec[0] * $$avec[0] * $SD < 0\n" if $debug;
			}
			$len_t += $olen;
		}else{
			print "\t\tcal_overs:no anno for $gene at node: $id, @$rvec\n" if $debug;
		}
		# even thouth the node is outside the gene, we still need to save it
		print "\t\tcal_overs:update gene ranges for $gene, $id, @$rvec\n" if $debug;
		$update_materials{$id} = $rvec; # save for later updati
	}

	#### FILTER, VV find overlapped reads with gene
	print "\t\tcal_overs: len_m: $len_m, len_t: $len_t\n\n" if $debug;

	if( $len_m > $OLEN  and $len_m / $len_t > $orat ){ ## adjustahble FILTER
	
		#my $ext_len = $fragl - $len_m;

		foreach my $id  ( keys %update_materials ){
			my $rvec = $update_materials{$id};
			update_gene_ranges(\%ranges, $id, $rvec ); ## update gene ranges
		}
		if($debug){
			# inspect structure of ranges
			foreach my $n ( keys %ranges ){
				foreach my $d ( keys %{$ranges{$n}} ){
					my @r = @{$ranges{$n}{$d}{R}};
					my $score = $ranges{$n}{$d}{S} // 0;
					foreach my $r (@r){
						my ($start, $end) = @$r;
						print "\t\tcal_overs:structure: $gene, id: $n, d: $d, range: [$start, $end], score: $score\n" ;
					}
				}
			}
		}
		return (\%ranges); #, $ext_len) ; # return ranges for this gene
	}else{
		return 0;
	}
}

sub update_gene_ranges{
	my ($ranges, $id, $init_r) = @_;
	my( $d,$start,$end,$score) =  @$init_r;

	print "\n\t\tsubfunction: update_gene_range $id $d  $start $end \n" if $debug ;
	if(! exists $$ranges{$id} or ! exists $$ranges{$id}{$d} or ! exists $$ranges{$id}{$d}{R} ){
		print "\t\tupdate_gene_ranges:not exists, create new ranges for $id, $d\n" if $debug;
		$ranges -> {$id}{$d}{R} = [[$start, $end]]; # save the range
		$ranges -> {$id}{$d}{S} = $score; # save the score
	}else{
		my @merged;
		#my $overlap_len = 0;
		$$ranges{$id}{$d}{S} += $score; # Sum the scoresi
		my $merge_tag = 0;
		# Iterate through existing ranges and merge if necessary
		for ( my $i = 0; $i < @{$ranges->{$id}{$d}{R}}; $i++ ) {
			my ($s, $e) = @{ $ranges->{$id}{$d}{R}[$i] };
			print "\t\tupdate_gene_ranges:check range: $s, $e\n" if $debug;
			if ($start <= $e && $end >= $s) { # means have overlap????
				my @pos = sort {$a <=> $b} ($start, $end, $s, $e);
				#$overlap_len += $pos[2] - $pos[1] + 1; # Calculate overlap length
				# Merge it
				my $startn = $pos[0];
				my $endn   = $pos[3];
				### look ahead in the array
				while(1){
					if($i + 1 < @{$ranges->{$id}{$d}{R}} && $end >= $ranges->{$id}{$d}{R}[$i + 1][0]) {
						($s, $e) = @{$ranges->{$id}{$d}{R}[ $i + 1]};
						$end = max($end, $e);
						print "\t\tupdate_gene_ranges:merge with next range: $s, $e\n" if $debug;
						$i++;
					}else{
						print "\t\tupdate_gene_ranges:no more overlaps with next range\n" if $debug;
						last; # No more overlaps, break the loop
					}
				}
				push @merged, [$startn, $endn];  # Add merged range
				print "\t\tupdate_gene_ranges:update merged range: [$startn, $endn]\n" if $debug;
				$merge_tag = 1; # Set merge tag
			} else {
				print "\t\tupdate_gene_ranges:no overlap, keep range: $s, $e\n" if $debug;
				push @merged, [$s, $e];  # Keep unchanged
			}
		}
		unless ($merge_tag) { # If no merge happened, just add the new range
			print "\t\tupdate_gene_ranges:no merge happened, add new range: [$start, $end]\n" if $debug;
			push @merged, [$start, $end];
		}
		# Sort merged ranges by start
		@merged = sort { $a->[0] <=> $b->[0] } @merged;
		$ranges->{$id}{$d}{R} = \@merged; # Update the ranges for this id and direction
	}
	print "\n" if $debug ;
}


##################

### sumary
# test genes overlap from the same block
sub expand_overlap{
	# expand_overlap(\%node_vec ,\%anno,\%align_inf, \%gene_inf);
	my ($node_vec, $anno, $align_inf, $gene_inf, $destiny) = @_;

	################################
	## 
	## Iterating each gene
	#	if, gene only exists in < 2 samples, then skip this gene
	#	
	#	1, for reads and nodes  associated/overlapped  with genes from the previous collection, delete them as keys in %node_vec
	#		save the reads and genes relattion into %collector;
	#  
	#   2, get the overhang margins , which are saved in %gene_inf{$gene}{B}. 
	#
	#   3, enter a endless loop: and iterate node ids from overhang region
	#
	#		4, baed on node, use node_vec to extract all read contains this node.
	#		5, then iterat all these reads: and calculate the sum  up all the score.
	#
	#	6, ranging all associated reads with overhang nodes by score from high to low
	#
	
	## iterate all genes in %gene_inf, collect associated reads
	my %overhangs;
	my %collector; # collect associated reads
    print "\n Entered expanding overlap model...\n" if $debug;
	foreach my $gene ( sort  keys %$gene_inf ){
		print "\texpand_overlap for gene: $gene\n" if $debug;
		my %oh;
		my @samples = keys %{ $$gene_inf{$gene}{S} };

		if( @samples  < $asam ){ ### adjustahble FILTER  <<+++++
			print "\texpand_overlap:genes appeared in less than 2 samples.@samples skip this gene.... next one\n" if $debug;
			if($debug){
				foreach my $r (keys %{ $$gene_inf{$gene}{R} } ){
					$$destiny{$r}{2} += 0 ;
				}
			}
			next;
		}else{
			foreach my $r ( keys %{ $$gene_inf{$gene}{R} } ){ # collect reads for this gene
				$$destiny{$r}{2} += 1 if $debug ;
				$collector{$gene}{$r} = 1; # save this read
				print "\texpand_overlap 1:collecting $gene => $r   \n" if $debug;
				# delete nodevec have this read and associated nodes
				map { delete $$node_vec{N}{$_}{$r} } keys %{$$node_vec{R}{$r}}; # delete this read from node_vec
				delete $$node_vec{R}{$r}; # delete this read from node_vec
			}
			%oh  = %{$$gene_inf{$gene}{B}};
		}

		# iteratate all overhang nodes. %oh is the initial nodes for this gene;
		while(1){
			print "\texpand_overlap 2: expanding overhang nodes for gene $gene : ", join "\t" ,  keys %oh, "\n" if $debug;
			my %ranges = %{ dclone( \%oh ) }; # get ranges if have
			%oh = ();

			## iterating reads defined in %node_vec, select reads that have most of overlapped nodes in %oh
			my %candidate_reads ;
			foreach my $n (keys %ranges){ # iteration ONE: nodes from %oh
				my $node_keep = 0;
				foreach my $read ( keys %{$$node_vec{N}{$n}} ){ # iteration TWO: reads from %node_vec
					if ( $collector{$gene}{$read}  ){ # if this read is not in %node_vec{R}, then skip
						print "already in collector $gene: $read $gene\n" if $debug;
						map { delete $$node_vec{N}{$_}{$read} } keys %{$$node_vec{R}{$read}}; # delete this read from node_vec
						delete $$node_vec{R}{$read}; # delete this read from node_ve
						next;
					}
					$candidate_reads{$read} +=  $$node_vec{N}{$n}{$read}[3]; # sum up the scores for this read
					$node_keep  = 1;
				}
				unless ($node_keep){
					delete $ranges{$n}; # delete this node if no reads extending
					print "delete node $n, no new  reads appeared\n" if $debug;
				}
			}

			my $update_boo  = 0;
			foreach my $read ( sort { $candidate_reads{$b} <=> $candidate_reads{$a} } keys %candidate_reads ){ # iteration THREE: reads from %node_vec
				print "test read: $read, score: $candidate_reads{$read} expanding from gene $gene \n" if $debug;
				my $len_m  = 0;
				my $len_t = 0;
				my %ranges_r  = %{ dclone( \ %ranges) };  # get ranges if have

				foreach my $id ( sort keys %{$$node_vec{R}{$read}} ){ # iterate all nodes for this read
					my $rvec = $$node_vec{R}{$read}{$id};
					print "checking id $id in read $read, vec: @$rvec\n" if $debug;
					if( $$gene_inf{$gene}{B}{ $id } ){
					# choose direction have the highest score
						my $gvec_d = $$gene_inf{$gene}{B}{$id};
						my @dires = sort { $$gvec_d{$b}{S} <=> $$gvec_d{$a}{S} } keys %{$gvec_d};
						if($debug){
							my @score = map { $$gvec_d{$_}{S} } @dires;
							print "directions for $id: ", join("\t", @dires), "\n";
							print "scores for $id: ", join("\t", @score), "\n";
						}
						my $gd = $dires[0]; # get the best direction
						my $gvec = $$gvec_d{$gd}{R};

						foreach my $gr ( @{ $gvec } ){ # iterate all ranges for this gene
							# calculate overlap length
							my $start = $$rvec[1] > $gr->[0] ? $$rvec[1] : $gr->[0];
							my $end = $$rvec[2] < $gr->[1] ? $$rvec[2] : $gr->[1];
							my $olen = ($end - $start + 1) > 0 ? ($end - $start + 1) : 0; # overlap length
							if( $$rvec[0] *   $gd  * $SD >= 0 ){ # if the direction is different, then skip
                                print "kep $id in $read , directio match: $gd <=> $$rvec[0], olen: $olen\n" if $debug;
								$len_m +=  $olen ;
							}else{
                                print "skip $id in $read direction mismatch: $gd <=> $$rvec[0] \n" if $debug;
							}
							$len_t += $olen;
						}
					}else{
						print "no anno for $gene at node: $id, @$rvec\n" if $debug;
					}
					my @init_region = @{ $rvec };
					print "extending update gene ranges for $gene, $id, @$rvec\n" if $debug;
					update_gene_ranges(\%ranges_r, $id, \@init_region ); ## update gene ranges
				}
				print "growth lenthm : $len_m, lentht: $len_t\n" if $debug;
				if( $len_m > $OLEN and $len_m / $len_t > $orat ){ ## adjustahble FILTER p
					
					#my $ext_len =  $$align_inf{$read}{L}  - $len_m;
					#my $$gene_inf{$gene}{L}  += $ext_len; 

					$$destiny{$read}{3}  += 1  if $debug ;
					print "growth to read $read\n" if $debug;
					#$$gene_inf{$gene}{B} = \%oh;
					# inspect structure of ranges
					if($debug){
						foreach my $n ( sort keys %ranges ){
							foreach my $d ( keys %{$ranges{$n}} ){
								my @r = @{$ranges{$n}{$d}{R}};
								my $score = $ranges{$n}{$d}{S} // 0;
								foreach my $r (@r){
									my ($start, $end) = @$r;
									print "Bstructure: $gene, id: $n, d: $d, range: [$start, $end], score: $score\n" ;
								}
							}
						}
					}
					%ranges = %ranges_r ; # update ranges for this gene
					# inspect structure of ranges
					if($debug){
						foreach my $n ( sort keys %ranges ){
							foreach my $d ( keys %{$ranges{$n}} ){
								my @r = @{$ranges{$n}{$d}{R}};
								my $score = $ranges{$n}{$d}{S} // 0;
								foreach my $r (@r){
									my ($start, $end) = @$r;
									print "Astructure: $gene, id: $n, d: $d, range: [$start, $end], score: $score\n" if $debug;
								}
							}
						}
					}

					$collector{$gene}{$read} = 2; # save this read
					print "\texpand_overlap: rescue read $read for gene $gene\n" if $debug;
					$update_boo = 1;
					map { delete $$node_vec{N}{$_}{$read } } keys %{$$node_vec{R}{$read}}; # delete this read from node_vec
					delete $$node_vec{R}{$read}; # delete this read from node_vec
					#return \%ranges; # return ranges for this gene
				}else{
					if($debug){
						$$destiny{$read}{3}  += 0 ;
						print "trim $read . discard\n" ;
					}
						#delete $$node_vec{R}{$read}; # delete this read from node_vec
				}
			}
			if($update_boo){
				print "update overhang nodes for gene $gene\n" if $debug;
				%oh = %ranges; # update overhang nodes
			}else{
				print "no more reads to collect for gene $gene\n" if $debug;
				last; # no more reads to collect, break
			}
		}
		$overhangs{$gene} = \%oh;
	}
	return \%collector;
}


sub cluster_output{
	
	################
	##
	##   to culster annotations, have two types of requirements:
	##
	##     1 > the annotations should be in a range specified in some haplotypes. this can create the possible merge from reference.
	##     2 > in sequencing data, a minimum number of reads and  in a minimum number of samples are required to demonstrated the real situation in data.
	#


	my($collector , $anno,  $align_inf,  $fh )  = @_;
	#my @genes_sorted  =  sort { (split /\|/, $a)[-1] <=> (split /\|/, $b)[-1] } keys %collector ; 
	my @genes_sorted = 	keys %$collector;
	##############################
	####  gene merge disabled ####
	##############################
	my @clusters; ### keep clusters

	if($dist_gmer < 0){
		# Non-clustering mode: each gene becomes its own “cluster” with all associated values
		foreach my $gene (@genes_sorted) {
			my @haps = @{$$anno{$gene}{H}};
			push @clusters, { genes => [$gene], shared_values => \@haps };
		}
	}else{
		## link genes that share the same reads and samples in %collector
		my %gene_connections;
		for (my $i = 0; $i < $#genes_sorted; $i ++){
			my $gene1 = $genes_sorted[$i];
			
			
			for (my $j = $i+1 ; $j < @genes_sorted; $j ++){
				my $gene2 = $genes_sorted[$j];

				######################################
				####  FILTER based on annotation #####
				######################################
				my @shares; 
				my $sim_c  = genes_similar($gene1, $gene2, $anno, $dist_gmer, \@shares );  ### adjustable
				
				# sim_c the count of haplotyps that support the distance of two genes are within 1000 bp
				#print "gene $gene1, gene $gene2 have $sim_c haps supports,haps %shares\n" ;
				## two genes are neighbors in < 1 haplotyps
				# $sim_c are counts of haptypes connecting two genes.
				# $shares = reference to an array that have the corresponding haps.  
				
				if($sim_c < $hps_c){
					next;
				}
			
				############################################
				### #### ===>>> filter based on reads  #####
				############################################
				
				my $same_reads = 0; # count of common reads overap both genes
				my %reads_source; # count of sources for common reads  coverage both gene
				foreach my $read ( keys %{$$collector{$gene1}} ){
					if(exists $$collector{$gene2}{$read}){
						#print "shared read $read\n";
						my($sample) = $read =~ /.+##(.+)$/;
						$reads_source{$sample} ++;
						$same_reads ++;
					}else{
					}
				}
				my $rsc = scalar keys %reads_source;
				if($same_reads < $nread_gmer  or  $rsc  <  $nsamp_gmer  ){ ### adjustahble FILTER
					#print "reads source not sufficent, discard connection.. $same_reads $nread_gmer $rsc $nsamp_gmer  \n";
					next;
				}
				#print "connect: gene $gene1 and $gene2 have overpaed $sim_c haps; share $same_reads reads in $rsc samples\n" ;
				$gene_connections{$gene2}{$gene1}{H} = \@shares ;
				$gene_connections{$gene1}{$gene2}{H} = \@shares ;
				$gene_connections{$gene1}{$gene2}{R} = $same_reads ;
				$gene_connections{$gene2}{$gene1}{R} = $same_reads ;
			}
		}

		### total similarity for ranking
		my %total_similarity;
		foreach my $gene (keys %gene_connections) {
			$total_similarity{$gene} = 0;
			foreach my $neighbor (keys %{$gene_connections{$gene}}) {
				my $sharec =  $gene_connections{$gene}{$neighbor}{R};
				#print "$gene = $sharec\n";
				#print "sum total similarity $gene $sharec\n";
				$total_similarity{$gene} += $sharec ;
			}
		}
		my @gs_sorted = sort { $total_similarity{$b} <=> $total_similarity{$a} } keys %total_similarity;
		
		my %visited;
		foreach my $gene ( @gs_sorted ) {
			next if $visited{$gene};

			# Start with “all values” of the first gene
			my %cluster_shared = map {$_ => 1} keys %{$$anno{$gene}{H}};
			#my %cluster_shared = map { $_ => 1 } (map { @{$gene_connections{$gene}{$_} || []} } keys %{$gene_connections{$gene}});
			#print "start with $gene $total_similarity{$gene} \n" ;	
			my @cluster = ($gene);
			$visited{$gene} = 1;
			
			while(1) {
				# Sort neighbors by number of shared values with current gene descending
				my %neighbors;
				for my $g (@cluster) {
					foreach my $n ( keys %{ $gene_connections{$g} } ){
						next if $visited{$n};
						
						my $ts = $total_similarity{$n};
						#print "TS:$n\n" unless defined $ts;
						$neighbors{$n} += $ts;
					}
				}

				last unless %neighbors;

				my $break;
				foreach my $neighbor ( sort { $neighbors{$b} <=> $neighbors{$a} } keys %neighbors ) {
					next if $visited{$neighbor};
					#print "\tfecth nei: $neighbor $neighbors{$neighbor}  $neighbors{$neighbor}\n";
					# Compute minimum shared values with all genes in the cluster
					#my $min_shared = min(map { scalar(@{$gene_connections{$neighbor}{$_} || []}) } @cluster);
						
					my %neighbor_values; 
					foreach my $member (@cluster) { 
						#print "\ttesting $neighbor and $member\n";
						foreach my $v (@{ $gene_connections{$neighbor}{$member}{H}  || [] }) { 
							#print "$v...\n";
							$neighbor_values{$v} ++ ; 
						} 
					} # Intersect with current cluster shared values 
					
					foreach my $v (keys %neighbor_values){
						delete $neighbor_values{$v}  unless $neighbor_values{$v}  ==  @cluster;
					}
					
					if(%neighbor_values){
						#print "\tcollect...\n";
						%cluster_shared = %neighbor_values;

						my @shares = keys %cluster_shared;
						my $sharesc = scalar @shares;
						if($debug){
							print "iter nei share $sharesc = @shares\n";
							print "include $neighbor... $sharesc  @shares \n";
						}
						delete $total_similarity{$neighbor};
						push @cluster, $neighbor;
						$visited{$neighbor} = 1;
					}else{
						#print "\tfound a empty neight ... break\n" ;
						$break = 1;
						last;
					}
				}
				if($break){
					last;
				}
			}
			
			if($debug){
				print "CCC:@cluster\n";
				print "SHAR:", keys %cluster_shared ,  "\n";
			}
			push @clusters, {genes =>  \@cluster, shared_values => [keys %cluster_shared] };
		}
	}


	for my $i (0 .. $#clusters) {
		my @genes_cluster = sort @{ $clusters[$i]{genes} };

		#### obtain median length of reference annotation
		#print "\nstart cluster $i = calculate median lenth .....\n"  ; 
		my @haplens;
		foreach my $v ( @{ $clusters[$i]{shared_values} } ){
			my @pos;
			foreach my $g ( @genes_cluster ){
				#print "clus$i = chrom:$v gene: $g pos1:$$anno{$g}{H}{$v}[0]  pos2:$$anno{$g}{H}{$v}[1]\ \n" if $debug ;
				push @pos, $$anno{$g}{H}{$v}[0], $$anno{$g}{H}{$v}[1];
			}
			@pos = sort {$a <=> $b} @pos;
			print "total sorted pos @pos\n" if $debug ;
			push @haplens, $pos[-1] - $pos[0];
		}
		my $hapcts = scalar @haplens;
		my $haplens_median = median(@haplens);
		#print "get haplen $haplens_median from @haplens\n" ;	
		###########################################
		###### find from the biggest clusters:  ###
		my $genes_str = join ":", @genes_cluster;

		#print "Cluster $i: ", join( ", ", @genes_cluster ), "\n";
		
		### retrieve reads from ....
		my %reads;
		foreach my $gene ( @genes_cluster ){
			# collect samples and reads
			foreach my $r (keys %{$$collector{$gene}}) {
				#print "collect read $r\n";
				push @{$reads{$r} }, "$gene:$$collector{$gene}{$r}";
			}
		}
		foreach my $read ( keys %reads ){
			my ($readid ,$mapn,$sample) = split /##/, $read;
			my $as = $$align_inf{$read}{S} // 0;
			my $flen = $$align_inf{$read}{L} // 0;
			my $genes = join ";", @{ $reads{$read} };
			print $fh "$genes_str\t$haplens_median\t$hapcts\t$readid\t$as\t$flen\t$mapn\t$genes\t$sample\n";
			#print "$genes_str\t$haplens_median\t$hapcts\t$readid\t$as\t$flen\t$mapn\t$genes\t$sample\n";
		}
	}
}


# DFS subroutine
sub dfs {
    my ($gene, $gene_connections, $cluster_ref, $visited ) = @_;
    return if $$visited{$gene};
    $$visited{$gene} = 1;
    push @$cluster_ref, $gene;
    for my $neighbor ( sort { @{$$gene_connections{$gene}{$b}} <=> @{$$gene_connections{$gene}{$a}} } keys %{ $$gene_connections{$gene} } ) {
        dfs($neighbor, $gene_connections, $cluster_ref, $visited ) unless $$visited{$neighbor};
    }
}

sub uniq {
	my %seen;
    return grep { !$seen{$_}++ } @_;
}

sub getfh {
    my ($tempfiles,$fh_assign) =   @_;

    my $fhi;
    while(! defined $fhi){
        print "waiting file handle...\n";
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

sub array_range {
	my ($min, $max) = (undef, undef);
	for my $val (@_) {
		$min = $val if !defined($min) || $val < $min;
		$max = $val if !defined($max) || $val > $max;
	}
	return ($min, $max);
}

sub first_index_gt {
    my ($val, $arr) = @_;
    my ($low, $high) = (0, $#$arr);
    while ($low <= $high) {
        my $mid = ($low + $high) >> 1;  # faster int division by 2
        if ($arr->[$mid] > $val) {
            $high = $mid - 1;
        } else {
            $low = $mid + 1;
        }
    }
    return $low <= $#$arr ? $low : undef;
}



=head
				if( $uni_reads{$path} and $uni_reads{$path} ne "$read##$mapn##$sample"  ){
					my $uread = $uni_reads{$path};
					print "\tgenome overlap same as $uread with path $path\n" if $debug ;
					if($gene_inf{$gene}{R}{$uread}){
						$gene_inf{$gene}{R}{"$read##$mapn##$sample"} = 1;
						$gene_inf{$gene}{S}{$sample} ++;
					}else{
						$destiny{ "$read##$mapn##$sample"  }{1}  = 0 if $debug  ;
					}
				}else{
					print "\tinitall cal_over subprocessss $read##$mapn##$sample\n" if $debug;
					my $nrange = cal_overs($gene, \%anno, $node_vec{R}{"$read##$mapn##$sample"} ,\%gene_inf ); # if have overlap, then update %overs and return 1. else do nothing, return 0
					if( $nrange) {
						$destiny{ "$read##$mapn##$sample"  }{1}  = 1 if $debug ;
						print "\tstep3 -> have overlap with $read##$mapn##$sample\n" if $debug;
						$gene_inf{$gene}{S}{$sample} ++;
						$gene_inf{$gene}{R}{"$read##$mapn##$sample"}  = 1;
						$gene_inf{$gene}{B} = $nrange; # save ranges for this read
					}else{
						if($debug){
							$destiny{ "$read##$mapn##$sample"  }{1}  = 0  ;
							print "\tstep3 -> no overlap with $read##$mapn##$sample\n" ;
						}
					}
				}



if($gap == 0){
	my $win = "G:$node_last" . "-".  $end ;
	
	#if($debug){
	if(1){
		print "\na new gene @genes come up.. $gc  \n" ;
		print "\tmight a gap\n" ;
		print "\ttest jointing of two annotation.. $node_last <=> $ar[1] \n" ;
		print "\tchecking the non-annotation-overlapped reads tabix query $win\n";
	}

	### =====> second condtions from mis-overlapped reads
	my $iter = $tabix_m -> query( $win ); # query missed alignments
	
	my $miss_l = $node_last;
	while ( my $aline = $iter -> next ) {
		#print "miss query: $aline\n" if $debug;
		my @lar = split /\t/,  $aline;

		if($lar[1] >= $miss_l){
			print "found a gap skip node $miss_l \n" if $debug;
			$tail = $miss_l;
			$gap = 1;
			last;
		}
		unless($lar[5]){
			print "empty path @lar\n";
		}
		path_vec($lar[5],\%vec);

		my @vks =  sort { $a <=> $b } ( grep { $_ <= $lar[1] and $_ >= $miss_l } keys %{$vec{B}} ); # only test node that less than current stat, since grep is used here.

		if(@vks){
			print "found a gap within node @vks\n" if $debug;
			$gap = 1;
			$miss_l = $vks[0];
			last;
		}else{
			#print "no gap cotinue..\n" if $debug;
		}
		$miss_l = $lar[2];
	}
	print "no-annotation reads checking finished...\n";
}





=cut
