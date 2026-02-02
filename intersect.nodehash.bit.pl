#!/usr/bin/perl
use warnings; use strict;
use List::Util qw(min max shuffle);
use Storable  qw(dclone) ;
use Bit::Vector::Overload;
use Parallel::ForkManager;
use Getopt::Long;
use File::Temp qw( tempfile tempdir );
use IO::Handle;
use File::Basename;
use File::Copy qw(copy);
use File::Path qw( make_path );
use Fcntl;
use Bio::DB::HTS::Tabix;

$| = 1;


print "Running: perl $0 @ARGV\n";
Getopt::Long::Configure("no_ignore_case");
####
my $input_gf;
my $input_ff;
my $infor_gf;
my $sam_f;

my $directory ;
my $pre ;
my $gbatch = 100 ; # one batch have at least 10 genes
my $lbatch = 1000;
my $samc = 10;
my $num_processes = 4; # number of processes to run in parallel
my $strand_type = "r";
my $debug = 0;
my $help;
my $estep = 10;

GetOptions(
    "g|genebody=s"  => \$input_gf,
    "f|flank=s"   => \$input_ff, 
	"i|infor=s"   => \$infor_gf,
    "s|sampleinfo=s" => \$sam_f,
	"D|Directory=s" => \$directory,
	"b|lbatch=i" => \$lbatch,
    "B|gbatch=i" => \$gbatch,
	"p|prefix=s" => \$pre, 
	"d|debug"     => \$debug,
    "n|num_processes=i" => \$num_processes, # Allow setting number of processes
	"t|strand=s" => \$strand_type, # strand type: f, r, or n
	"S|samc=i" => \$samc,
	"h|help" => \$help
) or die "Invalid options!\n";

# If --help is called or no arguments were passed
if ($help ) {
    print_help();
    exit;
}

###
# Extract directory part from prefix
my $pdir = "$directory/03_assign";
unless (-d $pdir) {
    make_path($pdir) or die "Failed to create directory $pdir: $!";
}
######## prepra required input files ######
##### determine input format  ############
unless( $input_gf){
	$input_gf  = "$directory/02_complete/all.genebody.intersect.tsv";
}
unless( $input_ff){
	$input_ff = "$directory/02_complete/all.flank.intersect.tsv";
}

my $in_fh;
if($input_gf =~ /\.gz$/) {
	open $in_fh, "zcat $input_gf |" or die $!;
} else {
	open $in_fh, "<", $input_gf or die "Cannot open $input_gf: $!";
}

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


###############################
## creat sample file handles ##
my %samples_fhs;
my @samids;
my %fragl;
open IN, $sam_f  or die "Cannot open $sam_f: $!";
while(<IN>){
	chomp;
	next if /^#/;
	my @ar = split /\t/;
	open my $fh, "> $pdir/$ar[0].assign" or die "Cannot open $pdir.$ar[0].assign : $!";
	push @samids, $ar[0];
	$fragl{$ar[0]} = $ar[1];
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
my %mis_files;

for my $i (0 .. $num_processes  ) {
	my $dir = tempdir("childassign_$i-XXXX", DIR => ".", CLEANUP => 1);

	my ($fh, $filename) = tempfile("data_$i-XXXX", DIR => $dir, UNLINK => 0);
    push @tmps, $filename;


	### collect
	my %fls = ( fh => $fh, file => $filename);
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


#############################
##### MAIN PROCESS ##########
##### read input flow==######
##############################
my $count = 0;
my $count_b = 0;
my %lines;

while( <$in_fh> ){
	my $line = $_;
	chomp ( $line );

	$count ++;
	$count_b ++;

	my @ar = split /\t/, $line;

	my @genes  = split /,/, $ar[-2] ;
	next unless ($samples_fhs{$ar[-8]});

	if ( !grep { exists $lines{$_} } @genes ) {
		print "No elements\n" if $debug ;
		my $gc = scalar keys %lines;
		if($gc  > $gbatch   or $count_b > $lbatch ){
			print ">>>>>process gene count $gc and lines $count_b\n";
			process_lines(\%lines);
			$count_b = 0;
		}
	}

	foreach my $g (@genes){
		push @{$lines{$g}{L}}, $line;	
		#### determing if it is a read with branch/skip.
		if($ar[-5] ne "."){
			my @brs = split /;/, $ar[-5];
			foreach my $b (@brs) {
				my($loc1,$loc2) = split /,/, $b;
				push @{ $lines{$g}{B} }, [$loc1,$loc2];
			}
		}
	}

}

#############################################
# Process any remaining lines after the loop
if(%lines){
	process_lines( \%lines );
}

if($use_fork){
	print "Parent is waiting childs\n";
	$pm->wait_all_children; # Wait for all child processes to finish
	print "All child processes finished.\n";
}

######## CLEANUP #######
### LAST combine outputs ###########
foreach my $tf (@tempfiles) {
    open my $in, '<', $tf->{file} or die "Can't open $tf->{file}: $!";
    print "--- Contents of $tf->{file} ---\n";
    while (<$in>) {
		chomp;
		my @ar = split /\t/;
		print "=:@ar\n" if $debug;
		my $sample = $ar[-3];
		unless(exists $samples_fhs{$sample}){
			next;
		}
		my $ofh = $samples_fhs{$sample};
		die "$sample $tf->{file} $_ \n" unless $ofh;
		#warn "$ofh" ,  join("\t", @ar[0..$#ar]), "\n";
		print $ofh join("\t", @ar[0..$#ar]), "\n"; # write to sample file
    }
	close $in;
	#unlink $tf->{file}; # delete temp file
}

#################################################
################ SUB FUNCTIONS ##################
#################################################

sub process_lines {
	
	my($ls) = @_;
	my $lines= dclone($ls);
	undef %$ls; 

	### initiate working horse subprocess: single or parallel:
	my $fhi;
	if($use_fork){
		($fhi ) = getfh(\@tempfiles,\%fh_assign) ;
		my $pid = $pm->start($fhi);

		if($pid){
			return;
		}
		my $ppid = getppid();
		print "\t---CHILD $$ parent $ppid : get filehandle $fhi with Count $count  count_b: $count_b\n";
	}else{
		$fhi = 0;
		print "\t---SINGLE  $$ : get filehandle $fhi with Count $count  count_b: $count_b\n";
	}

	my $fh = $tempfiles[$fhi]->{fh};

	my $tabix_i = Bio::DB::HTS::Tabix -> new(filename => $infor_gf) ;	
	my $tabix_f = Bio::DB::HTS::Tabix -> new(filename => "$input_ff");

	# iteration each island:
	my @gs_pre = sort keys %$lines;	
	foreach my $g  ( @gs_pre ){
		# get annotation
		# ===>  GENE Annotation from bed
		print "\n==============\ncollect gene path $g \n" if $debug ;	
		my ($gvec, $glen) = get_gene_vec( $g, $tabix_i ); # extract genes infors and save them to %anno and  also update %vecf;
		
		#################################################
		## accumulate read vec from overlapped reads#####
		##################################################
		my %visited;
		my %line_depo ;

		my @lines = @{$$lines{$g}{L}};
		get_read_vec(\@lines, \%line_depo, \%visited);	
		### reads extends
		extends(\%line_depo, $gvec, $g, $glen, $fh, 0 );	
	

		##########################
		### extends to flanks  ###
	    ##########################

		# extract lines;
		my $flines = get_flank_lines($g,$tabix_f);
		get_read_vec($flines, \%line_depo, \%visited);
		### reads extends
		extends (\%line_depo, $gvec, $g, $glen, $fh, 1 );	
	
	}

	my $time  =  scalar localtime();
	print "subprocess finished $time @gs_pre ..\n";
	$fh -> flush;
	if($use_fork){
		$pm -> finish(0);
	}
}

#                                   {V} = %vecf 
# %line_depo {k1:$nodeid}{$readm+id}{I} = @read_inf
#
# %vecf {$nodeid}{d : $direction, v = bit::vec, c = coverge}
#
# @read_inf = ($read, $score, $fraglen,multi-index(read), sample, ref)
#
########################
# subfunctuinos
# requirements:
# 1, can return how many samples are used to extends a boundary
# 2, return extended lines, save for the later use
# 3, 
sub extends {
	my ($line_depo, $gvec, $g,$glen,$fh, $sta) = @_;

	my $debug = 0;
	print "\t|\n\t---> extends module...\n" if $debug ;	
	my $loop = 0;

	while(1){
		####### 
		# variables used for Depth calculation ....
		#my $inct = 0;
		#my %ncov ;
		
		## calculate overlapes and save all reads to %lrcd;
		my %lrcd;	
		foreach my $n (keys %$gvec ){
			if( defined $$line_depo{$n} ){
				foreach my $read (keys %{$$line_depo{$n}}){
					unless( exists $lrcd{$read} ){
						$lrcd{$read} =  $$line_depo{$n}{$read};
						$lrcd{$read}{L} = 0;
					}
					######### 
					my $rvec = $$line_depo{$n}{$read}{V};
					my $gd = $$gvec{$n}{D};
					my $rd = $$rvec{$n}{d};
					if( $gd * $rd * $SD >= 0 ){		
						my $interl =  ( $$gvec{$n}{V} &  $$rvec{$n}{v} ) -> Norm() ;
						$lrcd{$read}{L} += $interl;
					}
				}
			}
		}

		#######################
		#   Filtering by overlap len and collect number of samples;
		my %sams;	
		foreach my $read (keys %lrcd ){
			my $ol = $lrcd{$read}{L};
			print "test reads: $read $ol\n" if $debug ;		
			if( $ol > 30 ){
				#my @read_inf = ($sample,$read,$ar[3],$ar[4],$mapn);
				$sams{ $lrcd{$read}{I}[0] } ++;
			}else{
				delete $lrcd{$read};
			}
		}
		
		my $last = 0;
		my @ss = ();
		my $sss = 0;
		if(%sams){
			@ss = keys %sams; # number of samples supporting extends
			map {$sss += $sams{$_}} @ss; # total number of reads support extends
		}
		if($sss < $samc*2 or @ss < $samc){  #### FILTER 
			print "\t==>Escape all loop, don't meet requires: $sss @ss...\n" if $debug;
			$last ++;
		}else{
			print "\t===>keep encountered reads, continue loop...\n" if $debug ;
		}
		# summarize Filters: 1 + 2
		if( $last ){
			last;
		}
		######## Quality passed then, collect incorpate informations
		my %news;
		foreach my $read ( keys %lrcd ){
			my $infs = join "\t", @{$lrcd{$read}{I}};
			print $fh "$g\t$glen\t$infs\t$sta\n";
			
			my @rinfo = @{$lrcd{$read}{I}};
			my $rvec = $lrcd{$read}{V};
			my $ovl = $lrcd{$read}{L};

			print "\t\tUpdate gene vecs  $read - $ovl : @rinfo...\n" if ($debug);
			foreach my $n ( keys %$rvec ){		
				delete $$line_depo{$n}{$read};

				my $rd = $$rvec{$n}{d};
				if( $$gvec{$n} ){
					my $gd = $$gvec{$n}{D};
					if($gd * $rd * $SD >= 0){
						my $mvec  =  $$gvec{$n}{V} |  $$rvec{$n}{v} ;
						
						$$gvec{$n}{V} = $mvec;
						$$gvec{$n}{C}  +=   $$rvec{$n}{c};
					}
				}else{
					$news{$n}{$rd}{C} +=  $$rvec{$n}{c};
					$news{$n}{$rd}{V} = defined $news{$n}{$rd}{V} ? ( $news{$n}{$rd}{V} | $$rvec{$n}{v} ) : $$rvec{$n}{v} ; # | is bitwise or operator
				}
			}
		}
			
		### update directions for new nodes
		# ATTENTION:::
		#  should i update directions infors for new encouted nodes? or it will be better to summarize by the overall reads not just focus on reads extended here.
		my $newnode_c = keys %news;
		print "\t NEW nodes count $newnode_c\n";
		foreach my $n ( sort keys %news ){
			my $d;
			my $v;
			my $c ; 
			if(exists $news{$n}{-1} and exists $news{$n}{1} ){
				my $rs = $news{$n}{-1}{C};
				my $fs = $news{$n}{1}{C};
				if($rs > 2 * $fs ){
					$d = -1;
					$c = $rs;
				}elsif($fs > 2 *  $rs){
					$d = 1;
					$c = $fs;
				}else{
					$d = 0;
					$c = $rs + $fs;
				}
				$v = $d? $news{$n}{$d}{V} : ($news{$n}{1}{V} | $news{$n}{-1 }{V} );
			}else{
				$d = exists $news{$n}{1}  ? 1 : -1;
				$v = $news{$n}{$d}{V} ;
				$c = $news{$n}{$d}{C} ;
			}
			print "\t||\n\t---->summarise new nodes:$n direct $d ,vector:", $v -> to_Bin() , "\n"  if $debug ;
			$$gvec{$n}{D} = $d;
			$$gvec{$n}{V} = $v;
			$$gvec{$n}{C} = $c;
		}	
		$loop ++;
	}
}

sub occu_len {
	my ($ndep) = @_;
	
	my $len;
	foreach my $n (keys %$ndep){
		$len += $$ndep{$n} -> Norm();
	}
	return $len;
}

sub get_read_vec{
	my ($lines,$line_depo, $visited ) = @_;	
	foreach my $lr (@$lines){
	
		my @ar = split /\t/, $lr;
		my ($readm,$sample,$as, $readl,$ref) = @ar[6,7,3,4,11] ;
		
		if ($$visited{$readm}){
			next;
		}else{
			$$visited{$readm} = 1;
		}

		## calcualte the real multi idx.
		my ($read,$mapn) = $readm=~ /(.+):(.+)$/;
		my ($t) = $mapn =~ /^(\d+)/ ;
		my $r = eval ($mapn);
		$mapn = "$t-$r";

		## ===> STEP 1 READ Alignmentes
		my $path = $ar[5];
		#my($mline, $uni_reads , $node_vec, $align_inf) = @_;
		my @read_inf = ($read,$as,$readl,$mapn,$sample,$ref);
		my @nodes = split /;/, $path;
		my %vecf;
		foreach my $n (@nodes){
			my($id,$d,$l,$s,$e,$t) =  split /,/, $n ; #=~ /(.)(\d+):(\d+),(\d+),(\d+)/;
			
			if($t eq ""){
				$t = 0;
			}
			my $vec =  Bit::Vector -> new($l);
			
			my $ie = $l - $e - 1;
			if($ie >= $s){
				$vec -> Interval_Fill($s, $ie );
			}
			
			$vecf{$id}{d} = $d;
			$vecf{$id}{v} = $vec ;
			$vecf{$id}{c} = $t > 0 ? ($l - $e - $s) / $l : 0 ;
			
			$$line_depo{$id}{$readm}{V} = \%vecf;
			$$line_depo{$id}{$readm}{I} = \@read_inf;  # $ar[3] score $ar[4] fragment length
		}
	}
}


#                                   {V} = %vecf 
# %line_depo {k1:$nodeid}{$readm+id}{I} = @read_inf
#
# %vecf {$nodeid}{d : $direction, v = bit::vec, c = coverge}
#
# @read_inf = ($read, $score, $fraglen,multi-index(read), sample, ref)
#



################
sub get_flank_lines{
	my ($g, $tabix) = @_;

	my @flines ;
	my ($t,$n,$p) = split /\|/, $g;

	my $iter = $tabix -> query("$t:$p-$p") ;

	my $line_number = 0;
	while (1) {
		my $line =  $iter -> next ;   # wrap in eval to catch errors
		last unless defined $line;
		$line_number++;
		# process your line here
		# print "$line\n";
		my @ar = split /\t/, $line ;
		if($ar[3] eq $g){
			#print "\tflank reads: $aline\n" if $debug;
			push @flines, join "\t", @ar[4..$#ar];
		}
	}
	die  "ERROR:gene flank : $t:$p-$p\n" unless $line_number ;
	return(\@flines);
}
sub get_gene_vec{
	my ($gene, $tabix ) = @_;

	my %vecf;
	my $glen;

	my ($g,$n,$p) = split /\|/, $gene;

	my $iter = $tabix -> query("$g:$p-$p");

	my %nodes;
	
	my $line_number = 0;
	while(1){
		my $line =  $iter -> next ;
		last unless defined $line;
		$line_number ++;
		my @ar = split /\t/, $line ;
		if($ar[3] eq $gene){
			#print "\tannotation extract $aline\n" if $debug;
			map{ $nodes{$_} = 1 } split /;/, $ar[5];
			$glen = $ar[-2];
		}
	}
	die "ERROR:gene annotation $g:$p-$p\n"  unless $line_number;

	foreach my $n ( sort keys %nodes ){
		my($id,$l,$s,$e) =  split /,/, $n ; #=~ /(.)(\d+):(\d+),(\d+),(\d+)/;
		
		if ( $id =~ /^(=|\+|-)(\d+)/ ){
			my $d;
			if($1 eq "+"){
				$d = 1;
			}elsif($1 eq "-"){
				$d = -1;
			}else{
				$d = 0; # no direction
			}
			$id = $2;
			my $vec =  Bit::Vector -> new($l);
			$vec -> Interval_Fill($s, $l-$e - 1);
			
			$vecf{$id}{V} = $vec;	
			$vecf{$id}{C} = 0;
			$vecf{$id}{D} = $d;
		}else{
			die  "Bad node id $id in $gene, skip this node\n" if $debug;
		}
	}
	return (\%vecf, $glen);
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
{
			############################
			###   backwards   ##########
			############################
			my $uloop = 1;	
			print "backwards...\n" if $debug ;
			while(1){
				my $ls = $head - $uloop * $estep + 1;
				my $le = $head - ($uloop - 1) * $estep;
				print "\t===Window $uloop: $ls - $le\n" if $debug ;
				if($le <= 0){
					last;
				}
				$ls  = $ls < 0 ? 0 : $ls;
				
				for my $i ($ls .. $le) {
					$headhash{$i} = 1;
				}

				my @ulines;
				foreach my $t_m(@tabix_m){
					my $ui  = $t_m -> query("G:$ls-$le");
					my @ul;
					while(my $l =  $ui -> next ){
						push @ul, $l;
					}
					push @ulines, @ul;
				}
				get_read_vec(\@ulines, \%line_depo, \%visited);

				# $sc means number of sampes encoutered when extending....
				my  ($bsam_r ,$eboo  ) = extends (\%line_depo, $gvec, $g, $glen, $fh, -1 * $uloop );
				my $sc = scalar keys %$bsam_r;
				if(  $sc  < $samc or ! $eboo ){	###### <<<<< FILTER
					# because no  lines are used and merge ontoo gvec, so, stop  ulooping...
					print "\tuloop $uloop end reads: $eboo samples: $sc \n" if $debug ;	
					last;
				}else{
					$hct += $eboo;
					print "\tuloop $uloop continue reads: $eboo samples: $sc \n" if $debug ;
					$uloop ++;
				}
			}
			####################
			##   Forward   #####
			####################
			my $floop = 1;
			print "forwards looking...\n" if $debug ;
			while(1){
				my $le = $head + $floop * $estep  - 1;
				my $ls = $head + ($floop - 1) * $estep;

				for my $i ($ls .. $le) {
					$headhash{$i} = 1;
				}
				
				print "\t===Window $floop $ls $le\n" if $debug ;
				my @flines;
				foreach my $t_m(@tabix_m){
					my $fi  = $t_m -> query("G:$ls-$le");
					my @fl;
					while(my $l = $fi -> next ){
						push @fl, $l;
					}
					push @flines, @fl;
				}
				get_read_vec(\@flines, \%line_depo, \%visited );

				my ($fsam_r , $eboo ) = extends (\%line_depo, $gvec, $g, $glen, $fh, $floop);
				my $sc = scalar keys %$fsam_r;
				if( $sc < $samc  or !$eboo ){  #### <<< FILTER
					# because no  lines are used and merge ontoo gvec, so, stop  looping...
					print "\tfloop $floop end reads: $eboo samples: $sc \n" if $debug ;	
					last;
				}else{
					print "\tfloop $floop continue reads: $eboo samples: $sc \n" if $debug ;	
					$hct += $eboo;
					$floop ++;
				}
			}

		foreach my $lr ( sort keys %lines  ){
			my @ar = split /\t/, $lr;
			my ($readm,$sample,$readl) = @ar[6,7,3] ;
			
			next if ($$visited{$readm});
			## calcualte the real multi idx.
			my ($read,$mapn) = $readm=~ /(.+):(.+)$/;
			my ($t) = $mapn =~ /^(\d+)/ ;
			my $r = eval ($mapn);
			$mapn = "$t-$r";

			## ===> STEP 1 READ Alignmentes
			my $path = $ar[5];
			print "\textends:step1: get_vecs for this read: $read  $path\n" if $debug ;
			#my($mline, $uni_reads , $node_vec, $align_inf) = @_;
			my $rvec = get_read_vec($path );

			## ===> STEP 2 OVERLAP calculation
			print "\textends:step2: overlap calculation .. $read##$mapn##$sample  \n" if $debug ;
			my $len_overlap = 0 ;
			foreach my $n ( keys  %$rvec ){
				if($$gvec{$n}){
					my $gd = $$gvec{$n}{D};
					my $rd = $$rvec{$n}{d};
					if( $gd * $rd * $SD >= 0 ){		
						my $interl =  ( $$gvec{$n}{V} &  $$rvec{$n}{v} ) -> Norm() ;
						$len_overlap += $interl;
					}
				}
			}
			
			### extending overhangs collect this reads
			if ( $len_overlap >  30 ) { ## overhang short than 30.  <<<< FILTER
				print "\t===> overlap long ($len_overlap) reads $lr..\n" if $debug ;
				$sams{$sample} ++;
				#print  "OUT:$g\t$glen\t$read\t$ar[3]\t$ar[4]\t$t-$r\t$sample\t$sta\n" if $debug ; 
				#print $fh "$g\t$glen\t$read\t$ar[3]\t$ar[4]\t$t-$r\t$sample\t$sta\n";
				#delete %lines{$lr};
				$outs{$readm}{S} = $sample;
				$outs{$readm}{L} = $lr;
				$outs{$readm}{O} =  "$read\t$ar[3]\t$ar[4]\t$t-$r\t$sample\t$sta";
				$outs{$readm}{V} = $rvec;
			}else{
				print "\t===> overlap short $len_overlap..\n";
			}
		}
	
		############### SUM up  FILTERS   ###############	
		my $last = 0;
		my @ss = ();
		my $sss = 0;
		if(%sams){
			@ss = keys %sams; # number of samples supporting extends
			map {$sss += $sams{$_}} @ss; # total number of reads support extends
		}
		print "\tEnd Loop $loop sumup:\n" if $debug;
		if($sss < $samc*2 or @ss < $samc){  #### FILTER 
			print "\t\tescape all loop, don't meet requires: $sss @ss...\n" if $debug;
			$last ++;
		}else{
			print "\t\tkeep encountered reads, continue loop...\n" if $debug ;
		}
		### summarize Filters: 1 + 2
		if( $last ){
			last;
		}

		foreach my $rm (sort  keys %outs){
			$$visited{$rm} = 1;	
			$samples_seen{ $outs{$rm}{S} } = 1;
			delete $lines{ $outs{$rm}{L} };
			print $fh "$g\t$glen\t$outs{$rm}{O}\n";
		
			### update gene
			my $rvec = $outs{$rm}{V};
			my $ohl = 0; # lengh of overhang
			print "\tUpdate gene vecs and collect new nodes from read $rm ...\n" if $debug ;
			foreach my $n ( keys %$rvec ){		
				my $rd = $$rvec{$n}{d};
				if( $$gvec{$n} ){
					my $gd = $$gvec{$n}{D};
					if($gd * $rd * $SD >= 0){
						my $mvec  =  $$gvec{$n}{V} |  $$rvec{$n}{v} ;
						
						##########################
						#  update boarder nodes  #
						if($$bnode{$n} and $mvec -> is_full()){
							delete $$bnode{$n};
						}

						$$gvec{$n}{V} = $mvec;
						$$gvec{$n}{C}  +=   $$rvec{$n}{c};
					}
				}else{
					$news{$n}{$rd}{C} +=  $$rvec{$n}{c};
					$news{$n}{$rd}{V} = defined $news{$n}{$rd}{V} ? ( $news{$n}{$rd}{V} | $$rvec{$n}{v} ) : $$rvec{$n}{v} ; # | is bitwise or operator
					$ohl += $$rvec{$n}{v} -> Norm2();
				}
			}
		}
		### update directions for new nodes
		# ATTENTION:::
		#  should i update directions infors for new encouted nodes? or it will be better to summarize by the overall reads not just focus on reads extended here.
		foreach my $n ( sort keys %news ){
			print "\textends:determine direc for new node $n\n" if $debug ;	
			my $d;
			my $v;
			my $c ; 
			if(exists $news{$n}{-1} and exists $news{$n}{1} ){
				my $rs = $news{$n}{-1}{C};
				my $fs = $news{$n}{1}{C};
				if($rs > 2 * $fs ){
					$d = -1;
					$c = $rs;
				}elsif($fs > 2 *  $rs){
					$d = 1;
					$c = $fs;
				}else{
					$d = 0;
					$c = $rs + $fs;
				}
				$v = $d? $news{$n}{$d}{V} : ($news{$n}{1}{V} | $news{$n}{-1 }{V} );
			}else{
				$d = exists $news{$n}{1}  ? 1 : -1;
				$v = $news{$n}{$d}{V} ;
				$c = $news{$n}{$d}{C} ;
			}
			print "\textends:new node $n direct $d ,", $v, "\n"  if $debug ;
			$$gvec{$n}{D} = $d;
			$$gvec{$n}{V} = $v;
			$$gvec{$n}{C} = $c;

			## encouter new boarder nodes
			if(! $v -> is_full()){
				$$bnode{$n} = 1;
			}
		}
}
