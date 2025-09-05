#!/usr/bin/perl
use warnings; use strict;
use List::Util qw(min max);
use Storable  qw(dclone) ;
use Parallel::ForkManager;
use Bit::Vector;
use Getopt::Long;
use File::Temp qw(tempfile);
use IO::Handle;
use File::Basename;
use File::Path qw(make_path);

####
my $input_f;
my $sampleinfo ;
my $batch_size = 10000;
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
$| = 1;

GetOptions(
    "i|input=s"  => \$input_f,
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

Optional arguments:
    -b, --batch <number>       Batch size for processing (default: 10000)
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
}

warn "use fork: $use_fork\n";
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
	open my $fh, "> $pre.$ar[0].assign.tsv" or die "Cannot open 02_quant/$ar[0].assin.tsv: $!";	
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
#########################################

my @tempfiles;
my %fh_assign;
for my $i (0 .. $num_processes  ) {                                                            
	my ($fh, $filename) = tempfile("childassign_$i-XXXX", UNLINK => 0, DIR => ".");
    my %fls = ( fh => $fh, file => $filename);
	if($debug){
		my ($fhd,$filenamed) = tempfile("childdestiny_$i-XXXX", UNLINK => 0, DIR => ".") ;
		$fls{fhd} = $fhd;
		$fls{filed} = $filenamed;
	}
    push @tempfiles, \%fls ;
}

################## 
##################
srand(1);
my $pm = Parallel::ForkManager->new($num_processes);
$pm->run_on_finish(sub {
    my ($pid, $exit_code, $ident, $exit_signal, $core_dump) = @_;
    if ($exit_signal) {
        print "Child $pid died from signal $exit_signal\n";
    }
    print "file handle $ident released by $pid..\n";
    delete $fh_assign{$ident};
});

#############################
##### main process ##########
##### read input flow==######
##############################
my $count = 0;
my $count_b = 0;
my $lgp = 0;
my @lines;
my $node_last;
my %gene_seen;

my $map_r;
my $bed_r;

while(<$in_fh>){
	chomp;
	#if($count_b % $batch_size == 0){
	#	print  "Batch LINE $count_b $count last gene line : $lgp\n";
	#}
	$count ++;
	$count_b ++;	
	my $line = $_ ;
	my @ar = split /\t/, $line;
	my @map = @ar[0..5];
	my $rid = join "###", @ar[6..7];
	my @bed = @ar[10..$#ar -1];

	my ($gene ) = $bed[3] ;
	## determine when to start child
	my $continue = 0;
	if( $count_b > $batch_size){	
		if ( $lgp  == 0  ){
			#print "nothing in this iteration...skip\n";
			print "block no intersection .. at $count $count_b. skip\n";
			@lines = ();
		}elsif($node_last  and $map[1] > $node_last  and ! exists $gene_seen{$gene} ){	
			print "collect a block at $count $count_b. \n";
			process_lines(\@lines);
		}elsif( $lgp > 0  and $count - $lgp > $batch_size){
			print "long tail at $count $count_b  ..\n";
			process_lines(\@lines);
		#}elsif( $count_b > $batch_size * 3){
			#process_lines(\@lines);
		}else{
			$continue  = 1;
		}
	}else{
		$continue = 1;
	}

	unless( $continue ){
		%gene_seen = ();
		$count_b = 0;
		$lgp = 0;
	}
	$node_last = $map[2];
		
	my @pack = ($rid);
	if($map_r and $map[5] eq  $$map_r[5]){
		push @pack, $map_r;
	}else{
		push @pack, \@map;
		$map_r = dclone(\@map);
	}
	if($bed_r and $$bed_r[3] eq $bed[3]){
		push @pack, $bed_r;
	}else{
		push @pack, \@bed;
		$bed_r = dclone(\@bed);
	}
	push @lines, \@pack;
	if($gene ne "."){
		$gene_seen{$gene} = 1;
		$lgp = $count;
	}
}

#############################################
# Process any remaining lines after the loop
if(@lines){
	process_lines( \@lines ); 
}

if($use_fork){
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

################
sub process_lines {
	my($lr) = @_;
	my $lines = dclone($lr);
	@$lr = ();

	### initiate working horse subprocess:
	my $fhi;	
	if($use_fork){
		($fhi ) = getfh(\@tempfiles,\%fh_assign) ;	
		my $pid = $pm->start($fhi);
		
		if($pid){
			return;
		}
		my $ppid = getppid();
		print "\t---CHILD $$ parent $ppid : get filehandle $fhi with Count $count $count_b\n";
	}else{
		$fhi = 0; 
		print "SINGLE PROCESS: Processing with fhi $fhi...\n";
		print "Line: $count $count_b in line to process\n" if $debug;
	}


	my $fh = $tempfiles[$fhi]->{fh};
	my $fhd = $tempfiles[$fhi]->{fhd} if $debug ;

	##### annotations   ###	
	my %node_vec; #### big ....
	my %anno; # {gene}{node} = [direct * $SD , start,end]
	my %align_inf; # {read} : {V} = [direct,star,end,type} {S} = AS  {L} = frag length
	my %gene_inf; #while samples have this gene;

	my %destiny ;

	my $node_last ; # the last node form record
	my %uni_reads;
	foreach my $lr ( @$lines ){		
		#my @ar = split /\t/, $lr;
		my ($readm,$sample) = split /###/, $$lr[0];
		my ($read,$mapn) = $readm=~ /(.+):(\d+)$/;

		my @map = @{$$lr[1]};
		my @bed = @{$$lr[2]};
		if($debug){
			print "read:$readm=$read+$mapn, $sample\n";	
			print "\nNewLine:\n>>BED: @bed\n" ;
			print ">>MAP: @map\nvvvvvvvvvvvv\n";
		}

		my $gene = $bed[3];	
		##########################
		#### **** VVVVVV **** ####
		### test if a break point detected. then summarize count.
		if( $node_last and $map[1] > $node_last   and !$anno{$gene} ){
			print "\nfound a disconnection... $map[1] $node_last\n" if $debug ;
		
			## check_overlap
			print "expanding genes:", scalar keys %gene_inf, "\n" if $debug;
			expand_overlap(\%node_vec ,\%anno,\%align_inf, \%gene_inf, $fh,\%destiny);
		
			undef %anno;
			undef %node_vec;
			undef %gene_inf;
			undef %align_inf;
		}else{
			if($debug){
				print "continue collect this read\n";
			}
		}
		$node_last = $map[2];
		########################################
		# get informaton from reads and anno ###
		# ######################################
		# ===> GENE Annotation from bed
		
		if($gene eq "."){
			1;
		}elsif(! defined  $anno{$gene} ){
			my @nodes = split /;/, $bed[-2];
			my $gl = $bed[-4];
			my @haps = split /;/, $bed[-1];

			$anno{$gene}{L} = $gl;
			#print "N:@nodes\n" if $debug;
			foreach my $n (@nodes){
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
				}else{
					die  "Bad node id $id in $gene, skip this node\n" if $debug;
				}
				#$d = 0 if ($SD == 0); # adjust direction
				print "\tannotation collect: $gene $id $d $s $e\n" if $debug;	
				$anno{$gene}{N}{$id} = [$d, $s, $l - $e - 1]; 
			}
			foreach my $h ( @haps ){
				my($chr,$start,$end,$dir) = $h =~ /^(.+):(\d+)-(\d+)#(-|\+)/;
				$anno{$gene}{H}{$chr} = [$start,$end,$dir];
			}
		}
		###########################
		## ===> READ Alignmentes
		my $as = $map[3];
		my $frag = $map[4];
		my $path = $map[5];
		#my ($chkn,$chkt) = @map[8,9];

		### node hash
		# can make a univeral graph, reduce memory . 
		print "\t>>get_vecs for this read: $read  $path\n" if $debug;
		if($uni_reads{$path}){
			my $uread = $uni_reads{$path};
			print "\tHave this $path before from $uread, reuse ...\n" if $debug;
			#print "WHYSHY @alls \n" unless $node_vec{R}{$uread};
			my $nv = $node_vec{R}{$uread};
			print "\t\treads saveed $nv\n" if $debug;
			$node_vec{R}{"$read##$mapn##$sample"} = $nv;
			foreach my $n (keys %$nv ){
				print "\t\tnode $n linked $node_vec{N}{$n}{$uread}\n" if $debug;
				$node_vec{N}{$n}{"$read##$mapn##$sample"} = $node_vec{N}{$n}{$uread};
			}
		}else{
			print "\tParse new path $path  ...\n" if $debug;
			get_vecs ( $path, "$read##$mapn##$sample", \%node_vec );
		}

		$align_inf{"$read##$mapn##$sample"}{S} = $as;
		$align_inf{"$read##$mapn##$sample"}{L} = $frag;
	
		if($gene ne "."){
			print "\t>>cal_overs for gene: $gene and $read $sample \n" if $debug;	
			if( $uni_reads{$path}  ){
				my $uread = $uni_reads{$path};
				print "\t\tgenome overlap same as $uread with path $path\n" if $debug ;
				if($gene_inf{$gene}{R}{$uread}){
					$gene_inf{$gene}{R}{"$read##$mapn##$sample"} = 1;
					$gene_inf{$gene}{S}{$sample} ++;
				}else{
					$destiny{ "$read##$mapn##$sample"  }{1}  = 0 if $debug  ;
				}
			}else{
				print "\t\tinitall cal_over subprocessss $read##$mapn##$sample\n" if $debug;
				my $nrange = cal_overs($gene, \%anno, $node_vec{R}{"$read##$mapn##$sample"} ,\%gene_inf ); # if have overlap, then update %overs and return 1. else do nothing, return 0
				if( $nrange) {
					$destiny{ "$read##$mapn##$sample"  }{1}  = 1 if $debug ;
					print " $gene have overlap with $read##$mapn##$sample\n" if $debug;
					$gene_inf{$gene}{S}{$sample} ++;
					$gene_inf{$gene}{R}{"$read##$mapn##$sample"}  = 1;
					$gene_inf{$gene}{B} = $nrange; # save ranges for this read
				}else{
					if($debug){
						$destiny{ "$read##$mapn##$sample"  }{1}  = 0  ;
						print " $gene no overlap with $read##$mapn##$sample\n" ;
					}
				}
			}
		}else{
			$destiny{ "$read##$mapn##$sample" }{1} = -1  if $debug ;
		}
		$uni_reads{$path} = "$read##$mapn##$sample" ;
	}
		
	# analys the last records
	expand_overlap(\%node_vec ,\%anno,\%align_inf, \%gene_inf,$fh, \%destiny);
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
	$fhd -> flush if $debug ;
	$fh -> flush;
	if($use_fork){
		$pm -> finish();
	}
}

########################
# subfunctuinos  

sub genes_similar {
    my ($gene1, $gene2, $anno, $range) = @_;
    $range //= 1000;  # default range 1000bp

    # Get chromosomes shared by both genes
	my $ct = 0;
    for my $chr (keys %{ $anno->{$gene1}{H} }) {
        next unless exists $anno->{$gene2}{H}{$chr};

        my $haps1 = $anno->{$gene1}{H}{$chr};
        my $haps2 = $anno->{$gene2}{H}{$chr};
		my ($start1,$end1,$dir1) = @$haps1;
		my ($start2,$end2,$dir2) = @$haps2;
        # Check same direction
		next unless $dir1 eq $dir2;
                # Check distance within range (overlap or nearby)
        if (intervals_within_range($start1, $end1, $start2, $end2, $range)) {
            $ct ++ ;  # true, similar
        }
    }
    return $ct; # no matching intervals found
}

sub intervals_within_range {
    my ($s1, $e1, $s2, $e2, $range) = @_;
    # Calculate minimal distance between intervals (overlap is 0 distance)
    if ($e1 < $s2) {
        return ($s2 - $e1) <= $range;
    } elsif ($e2 < $s1) {
        return ($s1 - $e2) <= $range;
    } else {
        return 1;  # intervals overlap
    }
}


#####

sub get_vecs{ # for a read, get all nodes and their vectors	
	my ( $path , $read, $node_vec ) = @_;
	my @nodes = split /;/, $path;
	#my @nodes = $path =~ /(<|>)(\d+)/g ;
	#@nodes = map{ [split /,/]  } @nodes;
	for( my $i = 0; $i < @nodes; $i ++ ){
		my($id,$l,$s,$e,$t) =  split /,/, $nodes[$i]  ; # this $id have direction symbol 	
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
	my($gene, $anno, $vecs, $gene_inf) = @_;
	my $len_m = 0;
	my $len_t = 0;
	print "\n\tsubfunction: cal_overs for gene: $gene\n" if $debug;
	my %ranges  = ();
	%ranges = %{ dclone ($gene_inf -> {$gene} -> {B})  } if $gene_inf->{$gene}{B}; # get ranges if have
	foreach my $id ( keys %$vecs ){
		my $rvec = $$vecs{$id} ;
		if( $$anno{$gene}{N}{ $id } ){	
			print "\tgene: $gene, id: $id, vec: @$rvec\n" if $debug;
			my $avec = $$anno{$gene}{N}{$id};
			my $start = $$rvec[1] > $$avec[1] ? $$rvec[1] : $$avec[1];
			my $end = $$rvec[2] < $$avec[2] ? $$rvec[2] : $$avec[2];
			my $olen = ($end - $start + 1) > 0 ? ($end - $start + 1) : 0; # overlap length
			 
			if( $$rvec[0]  *  $$avec[0] * $SD >= 0 ){ # if the direction is different, then skip
				if($debug){
					print "\tid: $id. olen $olen  $gene @$avec @$rvec\n";
				}	
				$len_m +=  $olen ;
			}else{
				print "\tdirection mismatch, skip $id, $gene, $$rvec[0] * $$avec[0] * $SD < 0\n" if $debug;

			}
			$len_t += $olen;
		}else{
			print "\tno anno for $gene at node: $id, @$rvec\n" if $debug;
		}
		# even thouth the node is outside the gene, we still need to save it
		print "\tupdate gene ranges for $gene, $id, @$rvec\n" if $debug;	
		update_gene_ranges(\%ranges, $id, $rvec ); ## update gene ranges
	}
	#### FILTER, VV find overlapped reads with gene
	print "gene: $gene, len_m: $len_m, len_t: $len_t\n\n" if $debug;
	
	if( $len_m > $OLEN  and $len_m / $len_t > $orat ){ ## adjustahble FILTER 
		# inspect structure of ranges
		foreach my $n ( keys %ranges ){
			foreach my $d ( keys %{$ranges{$n}} ){
				my @r = @{$ranges{$n}{$d}{R}};
				my $score = $ranges{$n}{$d}{S} // 0;
				foreach my $r (@r){
					my ($start, $end) = @$r;
					print "structure: $gene, id: $n, d: $d, range: [$start, $end], score: $score\n" if $debug;
				}
			}
		}
		return \%ranges; # return ranges for this gene
	}else{
		return 0;
	}
}

sub update_gene_ranges{
	my ($ranges, $id, $init_r) = @_;
	my( $d,$start,$end,$score) =  @$init_r;
	
	print "\nsubfunction: update_gene_range $id $d  $start $end \n" if $debug ;
	if(! exists $$ranges{$id} or ! exists $$ranges{$id}{$d} or ! exists $$ranges{$id}{$d}{R} ){
		print "not exists, create new ranges for $id, $d\n" if $debug;
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
			print "check range: $s, $e\n" if $debug;
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
						print "merge with next range: $s, $e\n" if $debug;
						$i++;
					}else{
						print "no more overlaps with next range\n" if $debug;
						last; # No more overlaps, break the loop
					}
				}
				push @merged, [$startn, $endn];  # Add merged range
				print "update merged range: [$startn, $endn]\n" if $debug;
				$merge_tag = 1; # Set merge tag
			} else {
				print "no overlap, keep range: $s, $e\n" if $debug;
				push @merged, [$s, $e];  # Keep unchanged
			}
		}
		unless ($merge_tag) { # If no merge happened, just add the new range
			print "no merge happened, add new range: [$start, $end]\n" if $debug;
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
sub expand_overlap{
	# expand_overlap(\%node_vec ,\%anno,\%align_inf, \%gene_inf);
	my ($node_vec, $anno, $align_inf, $gene_inf, $fh, $destiny) = @_;

	my %collector;
	## iterate all genes in %gene_inf
	foreach my $gene ( sort  keys %$gene_inf ){
		print "expanding gene: $gene\n" if $debug;
		my %overhang;
		my @samples = keys %{ $$gene_inf{$gene}{S} };


		if( @samples  < $asam ){ ### adjustahble FILTER  <<+++++
			print "genes appeared in less than 2 samples.@samples skip this gene...\n" if $debug;
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
				print "collecting first round read: $gene $r from sample \n" if $debug;
				# delete nodevec have this read and associated nodes
				map { delete $$node_vec{N}{$_}{$r} } keys %{$$node_vec{R}{$r}}; # delete this read from node_vec
				delete $$node_vec{R}{$r}; # delete this read from node_vec
			}
			%overhang  = %{$$gene_inf{$gene}{B}};
		}

		# iteratate all overhang nodes. %overhang is the initial nodes for this gene;
		while(1){
			print "expanding overhang nodes for gene $gene : ", join "\t" ,  keys %overhang, "\n" if $debug;
			my %ranges = %{ dclone( \%overhang) }; # get ranges if have
			%overhang = ();

			## iterating reads defined in %node_vec, select reads that have most of overlapped nodes in %overhang
			my %candidate_reads ;
			foreach my $n (keys %ranges){ # iteration ONE: nodes from %overhang
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
			
			print "candidate reads: ", join "\t",  keys  %candidate_reads, "\n" if $debug;
			
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
					$$destiny{$read}{3}  += 1  if $debug ;
					print "growth to read $read\n" if $debug;
					#$$gene_inf{$gene}{B} = \%overhang;
					# inspect structure of ranges
					if($debug){
						foreach my $n ( sort keys %ranges ){
							foreach my $d ( keys %{$ranges{$n}} ){
								my @r = @{$ranges{$n}{$d}{R}};
								my $score = $ranges{$n}{$d}{S} // 0;
								foreach my $r (@r){
									my ($start, $end) = @$r;
									print "Bstructure: $gene, id: $n, d: $d, range: [$start, $end], score: $score\n" if $debug;
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

					$collector{$gene}{$read} = 1; # save this read
					print "rescue read $read for gene $gene\n" if $debug;
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
				%overhang = %ranges; # update overhang nodes
			}else{
				print "no more reads to collect for gene $gene\n" if $debug;
				last; # no more reads to collect, break
			}
		}
	}

	## link genes that share the same reads and samples in %collector
	my %gene_connections;
	foreach my $gene1 ( keys %collector ){
		foreach my $gene2 ( keys %collector ){
			if($gene1 eq $gene2){
				$gene_connections{$gene1}{$gene2} = 1; # self connection
				next; # skip self
			}else{
			### test gene1 and gene2 at the same strand
				
				my $sim_c;
				if($dist_gmer < 0 ) {
					$sim_c = 0;
				}else{
					$sim_c = genes_similar($gene1,$gene2,$anno, $dist_gmer);  ### adjustable
				}

				if($sim_c < 1){
					next;
				}
			}
			my $same_reads = 0;
			my %reads_source;
			foreach my $read ( keys %{$collector{$gene1}} ){
				if(exists $collector{$gene2}{$read}){
					my($sample) = $read =~ /##(.+)$/;
					$reads_source{$sample} = 1;
					$same_reads ++;
				}
			}
			# $nread_gmer = 2;
			# $nsamp_gmer = 2;
			if($same_reads > $nread_gmer  and keys %reads_source > $nsamp_gmer  ){ ### adjustahble FILTER
				print "gene $gene1 and $gene2 share $same_reads reads in ", scalar keys %reads_source, " samples\n" if $debug;
				$gene_connections{$gene1}{$gene2} = 1;
				$gene_connections{$gene2}{$gene1} = 1;
			}
		}
	}

	my %visited;
	my @clusters;
	
	# Iterate through all genes
	for my $gene (sort keys %gene_connections) {
		next if $visited{$gene};
		my @cluster;
		dfs($gene, \%gene_connections, \@cluster, \%visited); # Perform DFS to find all connected genes
		push @clusters, \@cluster if @cluster;
	}

	# Output
	for my $i (0 .. $#clusters) {
		my @genes_cluster = sort @{$clusters[$i]};
		my $genes_str = join ":", @genes_cluster;
		
		if($debug){
			print "Cluster $i: ", join( ", ", @genes_cluster ), "\n";
		}	
		my %reads;	
		foreach my $gene ( @genes_cluster ){
			# collect samples and reads
			foreach my $r (keys %{$collector{$gene}}) {
				push @{$reads{$r} }, "$gene:$collector{$gene}{$r}";
			}
		}
		foreach my $read ( keys %reads ){
			my ($readid ,$mapn,$sample) = split /##/, $read;
			my $as = $$align_inf{$read}{S} // 0;
			my $flen = $$align_inf{$read}{L} // 0;
			my $genes = join ";", @{$reads{$read}};
			print $fh "$genes_str\t$readid\t$as\t$flen\t$mapn\t$genes\t$sample\n";
		}
	}
}
# DFS subroutine
sub dfs {
    my ($gene, $gene_connections, $cluster_ref, $visited ) = @_;
    return if $$visited{$gene};
    $$visited{$gene} = 1;
    push @$cluster_ref, $gene;
    for my $neighbor (keys %{ $$gene_connections{$gene} }) {
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
