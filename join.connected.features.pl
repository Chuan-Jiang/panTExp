#!/usr/bin/perl
use warnings; use strict;
use Bio::DB::HTS::Tabix;
use Getopt::Long;
use Parallel::ForkManager;
use Fcntl ':flock';
use List::Util qw(sum);


$| = 1;

my $sam_f;
my $directory;
my $info_f;
my $cpu = 4;
my $hpsc = 10; # required haps to support connection
my $samc = 10; # (scalar keys %assign_fs ) / 4;		$sam_t = $sam_t > 2 ? $sam_t : 2;
my $ninfo_f;


GetOptions(
	"sam=s" => \$sam_f,
	"dir=s" => \$directory,
	"info=s" => \$info_f,
	"cpu=i" => \$cpu,
	"hpsc=i" => \$hpsc,
	"samc=i" => \$samc,
	"node=s" => \$ninfo_f
) or die "Error in command line arguments\n";

####
my $dist_gmer = 2000;
my $bsize = 1000;
my $debug = 0;

####
my $pm = Parallel::ForkManager -> new($cpu);
# Parent process
$pm->run_on_finish(
    sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core, $data) = @_;
        if ($data->{stderr}) {
            print "Child $pid error: $data->{stderr}";
        }
    }
);


#######################
my $linkf = "$directory/03_assign/feature.links";
open my $fh, '>', $linkf;
close $fh;
#open OUT , ">$directory/03_assign/feature.links" or die $!;

#### collect assign files

my %assign_fs;
open IN, $sam_f  or die "Cannot open $sam_f: $!";
while(<IN>){
	chomp;
	next if /^#/;
	my @ar = split /\t/;
	my $f = "$directory/03_assign/$ar[0].assign";
	$assign_fs{$ar[0]} =  $f;
}
close IN;


###########################################################
## determing whether to join in each sequencing sample  ###
###########################################################
my %gene_connections;
my %gene_covs;
foreach my $sam ( keys %assign_fs ){
	my $af = $assign_fs{$sam};
	
	print "reading $sam : $af\n";
	open AF, $af or die $!;
	
	### read in each file and assign read in as two hashes
	my %muls_r;
	while(<AF>){
		chomp;
		my @ar = split /\t/;
		$muls_r{"$ar[2]$ar[5]"}{$ar[0]} = $ar[-1]; # $ar[-1]  in which loop, the reads collected
		if($ar[-1] == 0){
			$gene_covs{ $ar[0] }{$sam} += 1 / $ar[1] ; # only count overlapped reads;
		}
	}
	close AF;


	#################################################
	### test and summarize connection supportings ###\
	
	print "\tgene connections\n";
	
	## sample process step 1 : iteratie and collect any connection supporing reads
	foreach my $r (sort keys %muls_r){
		if(scalar keys %{$muls_r{$r}} > 1){
			my @fs = sort keys %{$muls_r{$r}};
			for(my $i = 0; $i < @fs; $i ++){
				for (my $j = $i + 1; $j < @fs; $j ++){
					my $gene1 = $fs[$i];
					my $gene2 = $fs[$j];
					$gene_connections{$gene1}{$gene2}{$sam} ++;
					$gene_connections{$gene2}{$gene1}{$sam} ++;
				}
			}
		}
	}
	print "sam $sam processed\n";
}

#==============

my @gs_sorted = sort cmp_gene_covs keys %gene_covs;   

my $total = scalar @gs_sorted;
print "\ntotal genes to cluster $total \n" ;


#++++++++++++++++++++++++++++++++++
###################################################
#### collect batches and issue child process  #####
my @block;
my %visited;
my $bs = 0;
foreach my $gene ( @gs_sorted ) {
	next if $visited{$gene};
	# Start with “all values” of the first gene
	
	## iterate to extract all genes connected with $gene
	my %gene_connections_local  = get_all_connected($gene, \%gene_connections);	
	my @gbatch = sort keys %gene_connections_local;
	my %gene_covs_local = %gene_covs{ @gbatch };

	#print "current block batch @gbatch\n";
	if(@gbatch < 2){
		next;
	}

	$bs += @gbatch;
	push @block , [\%gene_connections_local, \%gene_covs_local];
	unless(@block < $bsize){
		$total -= $bs;
		my $bls = scalar @block;
		print "batch_size:$bs remain:$total blocksize $bls\n";
		process(@block);
		undef @block ;
		$bs = 0;
	}
	map{$visited{$_} = 1} @gbatch;
}
process (@block);


$pm->wait_all_children;
print "scriot finished\n";



#################################################

sub  process {
	my @block = @_;
	
	my @out;
	my $tabix_i = Bio::DB::HTS::Tabix -> new(filename => $info_f) ;

	my %gene_connections;
	my %gene_covs;

	my $pid = $pm -> start();
	if($pid){
		return;
	}else{
		my $ppid = getppid();
		print "\t---Child $$ parent $ppid started \n";
	}

	#### iterate blocks
	foreach my $bl (@block){
		my ($g_c_l, $e_s_l ) = @$bl;
		
		%gene_connections = %$g_c_l;
		%gene_covs = %$e_s_l;
		
		my @batch = sort keys %gene_covs;
		print "\n>>>\nBlock: $g_c_l $e_s_l genes: @batch\n" if $debug ;

		#####################################
		#### collection gene annotations ####
		# at last obtain a annotation hash %anno;
		my %gposs ; # gene position , that is gene ids
		for my $item (@batch) {
			my @ar = split /\|/, $item;
			push @{ $gposs{$ar[0]} }, $ar[-1];
		}
		my %tranges = cluster_ranges(%gposs);
		
		my %anno;
		print "Step1=> tab annotation for genes\n" if $debug ;
		tab_anno($tabix_i, \%tranges, \%anno);

		##############################
		##### calculate clusters #####
		my @clusters; ## calculate clusters;
		my @gs_sorted = sort cmp_gene_covs  @batch;
		
		print "Step2=> sort genes in batch according to gene coverage  @gs_sorted\n" if $debug ;
		
		print "Step3=>  Clustering and filtering teration init...\n" if $debug;
		my %visited;
		foreach my $gene (@gs_sorted){	
			next if $visited{$gene};
			$visited{$gene} = 1;
			my %cluster;
			my @gclu = ($gene);	
			print "\n\t==>query elements from cluster: @gclu\n"  if $debug ;	
			while(1) {
				# Sort neighbors by number of shared values with current gene descending
				my %news;
			
				
				foreach my $g (@gclu){
					my @neighbors = sort cmp_gene_covs keys %{ $gene_connections{ $g } } ; 
					print "\tNeighbors: @neighbors\n" if $debug ;	
					
					foreach my $nei ( @neighbors ) {
						my %shares;

						##### Filter baed on haplotypes
						my ($sct, $idis)  = genes_similar_haps($gene, $nei, \%anno, $dist_gmer, \%shares ); 
						
						if($sct  < $hpsc){  # if two genes have less than $hpsc haplotype supporting, then skip this pair.
							print "\t------>skip = $g and $nei have few haps support $sct\n"  if $debug ;
							goto DEL;
							#next; # <<<<= FILTER ===== shares
						}	
						#### Filters in sequenced samples;
						my $samples = 0;
						foreach my $sam (sort keys %{$gene_connections{$g}{$nei} } ){
							my $g1c = $gene_covs{$g}{$sam} // 0  ;
							my $g2c = $gene_covs{$nei}{$sam} // 0 ;
							my $inc = $gene_connections{$g}{$nei}{$sam}/$idis;

							my $g1cl = log2($g1c + 1);
							my $g2cl = log2($g2c + 1);
							my $incl = log2($inc + 1);
							my $diff = abs($g1cl - $g2cl);

							if($diff < 2 ){
								if(abs($incl - $g1cl * $g2cl) < 2){
									$samples ++;
								}else{
									print "\t\t$sam inter region few support $inc $g1c $g2c ($g + $nei)\n" if $debug ;
								}
							}else{
								print "\t\t$sam not cover two genes $g1c $g2c $inc ($g + $nei)\n" if $debug ;
							}
						}

						if($samples <  $samc){  # <<<<<==== FILTER based on population
							print "\t------>skip = $g and $nei have few sampls $samples\n" if $debug;
							goto DEL;
						}

						#######
						print "\t -----> kepp $g and $nei connected...\n" if $debug ;
						merge_gene_hash( \%cluster, \%shares );
						$news{$nei} = 1;
						$visited{$nei} = 1;
						
						DEL:
						delete $gene_connections{$g}{$nei};
						delete $gene_connections{$nei}{$g};

					}
				}
				if(%news){
					@gclu = sort keys %news;
					print "\t\t include new query elements: @gclu\n" if $debug ;
				}else{	
					last;
				}
			}

			my @ginc = sort keys %cluster;
			if( @ginc ){
				push @clusters, \%cluster;
				$total  = $total - scalar @ginc;
				print "clustering end with gene $gene: @ginc. Left: $total\n" if $debug  ;	
			}else{
				print "\t Empty cluster for gene $gene\n" if $debug ;
			}
		}

		################################
		#enter cluster length  estimate module
		if($ninfo_f){
			my $tabix_n = Bio::DB::HTS::Tabix -> new(filename => $ninfo_f) ;
			cluster_len_node(\@clusters,\%anno,$tabix_n, \@out);
		}else{
			cluster_len(\@clusters,\@out);
		}

		print "length esimated...\n" if $debug ;
	}

	open my $fh, ">>", $linkf or die $!;
	flock($fh, LOCK_EX) or die "Cannot lock: $!";   # exclusive lock
	map { print $fh "$_\n" } @out; 
	flock($fh, LOCK_UN);   # release lock


	$pm -> finish();
}


##########################
##########################
#
sub cluster_len_node {
	my ($cs, $anno, $tabix_n, $out) = @_;
	
	my @clusters = @$cs;

	print "\n>>>>estimate length of each cluster $#clusters...\n" if $debug ;
	for my $ci ( 0 .. $#clusters ){
		## estimate length of each cluster;
		my %clus = %{$clusters[$ci]};
		my @genes = sort { scalar keys %{$clus{$b}} <=> scalar keys %{$clus{$a}} } keys %clus;
	}

}



sub cluster_len {
	my ($cs,$out) = @_;
	my @clusters = @$cs;

	print "\n>>>>estimate length of each cluster $#clusters...\n" if $debug ;
	for my $ci ( 0 .. $#clusters ){
		## estimate length of each cluster;
		my %clus = %{$clusters[$ci]};
		my @genes = sort { scalar keys %{$clus{$b}} <=> scalar keys %{$clus{$a}} } keys %clus;
		my $clus_name = join ":", @genes;	
		
		####### use the first genes as iniciated...
		my $fg = shift @genes;
		my %universe_regions; # record the transform universe vregions # {gene} = [0,median] 
		my %usedchr; # recored the used chromosome regionsi # usedchr {hap} = 1
		my @flens;

		#### save the frontier gene to universe region and usedchr
		foreach my $chr ( keys %{$clus{$fg}} ){
			my @cors = @{$clus{$fg}{$chr}};
			push @flens, $cors[1] - $cors[0];
			$usedchr{$chr} = 1;
		}
		print "calcualte length for fortier $fg @flens\n" if $debug ;
		my $fgl = median(@flens);
		$universe_regions{$fg} = [0,$fgl,1];
		my @fchr = sort keys %usedchr;
		print "\nFrontier gene $fg length: $fgl. CLuster: @genes Haps:@fchr\n" if $debug ;
		
		### incorporate more genes
		my %seen;
		while( @genes ){
			
			my $ngi;
			my @haps;
			
			print "begin compare frontier to genes @genes\n" if $debug ;
			for (my $i = 0; $i < @genes; $i ++){
				my $g = $genes[$i];
				
				my @ucs = sort keys %usedchr;
				my @tcs = sort keys %{$clus{$g}};
				
				my @in_haps = grep { exists $usedchr{$_} } sort keys %{$clus{$g}};
				
				my $haps_ct = scalar @in_haps;
				print "\n\tcompare to $i $g share $haps_ct haps \n" if $debug ;
				if( $haps_ct   > scalar @haps ){
					$ngi = $i;
					@haps  =  @in_haps ;
				}
			}
			print "\tnext gene index $ngi @haps\n" if $debug ;

			my $ng = $genes[$ngi];
			splice(@genes, $ngi, 1);


			## update used haps
			foreach my $chr ( sort keys %{$clus{$ng}} ){
				if(! exists $usedchr{$chr} ){
					print "\t==include $chr\n";
					$usedchr{$chr} = 1;
				}
			}

			print "next gene $ng index is $ngi, share:\n" if $debug ;

			my @ngs_p; # information in populaiton -P
			my @nge_p;
			my @ngd_p;

			foreach my $h (@haps){
				my ($ngs, $nge, $ngd) = @{$clus{$ng}{$h}};
				print "\tnext gene: $ng at $h range : $ngs - $nge $ngd\n" if $debug ;
				foreach my $fg (keys %universe_regions ){  # gene have fixed on universe corrdi
					next unless ( exists $clus{$fg}{$h} );
					my ($fgs, $fge, $fgd) = @{$clus{$fg}{$h}};
					my ($fgs_r, $fge_r, $fgd_r ) = @{ $universe_regions{$fg} };

					print "\tpre gene: $fg  at $h range : $fgs - $fge  $fgd universe $fgs_r - $fge_r  $fgd_r \n" if $debug ;
					my $sshift = $fgs - $ngs;
					my $eshift = $fge - $nge;

					my $ngs_r = $fgs_r + ($fgd * $fgd_r) * $sshift;
					my $nge_r = $fge_r  + ($fgd * $fgd_r) *  $eshift;
					($ngs_r, $nge_r) = sort {$a <=> $b} ($ngs_r,$nge_r);
					print "\ttransform next gene to uiver region: $ngs_r - $nge_r\n" if $debug ;
					push @ngs_p, $ngs_r; #
					push @nge_p, $nge_r; #
					push @ngd_p , $fgd * $fgd_r * $ngd;
				}
			}
			my $ns = median(@ngs_p);
			my $ne = median(@nge_p);
			my $nd = median(@ngd_p);
			
			print "\testimated range $ng $ns $ne  $nd \n" if $debug ;
			$universe_regions{$ng} = [ $ns, $ne , $nd ];
		}
		
		print "length esimated for cluster $clus_name\n" if $debug ;
		my @uni_range ;
		foreach my $g (sort keys %universe_regions ){
			my @coor = @{$universe_regions{$g}};
			unless (@uni_range){
				@uni_range = @coor[0,1];
			}else{
				if($coor[0] < $uni_range[0]){
					$uni_range[0] = $coor[0];
				}elsif($coor[1] > $uni_range[1]){
					$uni_range[1] = $coor[1];
				}
			}
		}
		my $len = $uni_range[1] - $uni_range[0] + 1;
		push @$out, "$clus_name\t$len";
		#print "$clus_name\t$len\n";
	}

}


sub log2 {
    my ($x) = @_;
    return log($x) / log(2);
}


sub cluster_ranges {
    my %tnums = @_;
	my $debug = 0;
	my %tranges;
	foreach my $type (keys %tnums ){
		my @nums = @{$tnums{$type}};
		@nums = sort { $a <=> $b } @nums;
		print "\tCluster:$type:@nums\n\t\t" if $debug ;

		my %includs;
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
				my $rg = "$type:" . ($start-1)."-$end";
				my %inc = %includs;
				$tranges{$type}{$rg} = \%inc;
				my @ins = keys %includs;
				print "RANGE:$type $rg @ins\t" if $debug ;
				# start a new range
				$start = $n;
				$end   = $n;
				%includs = ($n , 1);
			}
		}

		# push last range
		my $rg = "$type:" . ($start-1) . "-$end";
		$tranges{$type}{$rg} = \%includs;
		my @ins = keys %includs;
		print "RANGE:$type $rg @ins\t" if $debug ;
	}
	print "\n" if $debug ;
	return %tranges;
}

sub get_all_connected {
    my ($start, $graph ) = @_;
    
	my %visited;
    my %sgraph;

	my @stack = ($start);

    while (@stack) {
        my $g = pop @stack;

        next if $visited{$g};
        $visited{$g} = 1;

        # neighbors of this gene
        foreach my $neighbor (keys %{ $graph->{$g} }) {
			unless ($visited{$neighbor}){
				$sgraph{$g}{$neighbor} = $$graph{$g}{$neighbor};
				$sgraph{$neighbor}{$g} = $$graph{$g}{$neighbor};
				push @stack, $neighbor ;
			}
        }
    }

    return  %sgraph;
}


sub median {
	my @data = @_;
    return undef unless @data;  # handle empty array
    # sort numerically
    @data = sort { $a <=> $b } @data;

    my $count = @data;
    if($count == 1){
        return $data[0];
    }else {
        # odd number of elements → return the middle one
        my $idx = int($count / 2);
		return $data[ $idx ];
    }
}

sub tab_anno{
	my($tabix, $tranges, $anno) = @_;
	my $debug = 0;	
	foreach my $t ( keys %$tranges ){
		my @ranges = keys %{$$tranges{$t}};
		print "\n\tAnno $t @ranges \n" if $debug ;
		foreach my $r ( @ranges ){
			my %aims = %{$$tranges{$t}{$r}};
			my $iter = $tabix -> query($r);
			
			my @des = keys %aims;
			print "query $r..@des\n" if $debug ;
			
			while (my $aline = $iter -> next) {
				my @ar = split /\t/, $aline ;
				next unless $aims{$ar[2]};
				my $gene = $ar[3];
				
				### --------------------------
                ### Parse hap strings (ar[4])
                ### --------------------------
                for my $h (split /;/, $ar[4]) {
                    my($chr,$start,$end,$dir) = 
                        $h =~ /^(.+):(\d+)-(\d+)#(-|\+)/;

                    next unless defined $chr;

                    $dir = $dir eq "+" ? 1 : -1;

                    # Directly store into anno
                    $$anno{$gene}{C}{$chr} = [$start, $end, $dir];

                    print "\t\t$gene $chr $start $end $dir\n"
                        if $debug;
                }



				my @ns = split /;/, $ar[5];
				foreach my $n (@ns){
					my ($id,$l,$h,$t) = split /,/, $n;
					$id = abs ($id); #$id > 1 ? 1 : -1;
					$$anno{$gene}{N}{$id} = [$h,$t];
				}
			}
		}
	}
}

#######
sub genes_similar_haps  {
    my ($gene1, $gene2, $anno, $range, $shares) = @_;
    $range //= 1000;  # default range 1000bp
	my $debug = 0;

	print "\n\tgenes similar module $gene1 $gene2\n" if $debug;
    # Get chromosomes shared by both genes
	#my @share_haps;
	my @haplens;
	my $cts = 0;
	
	my @chr1 = keys %{ $anno->{$gene1}{C} };

	my @ils;	
	for my $chr (keys %{ $anno->{$gene1}{C} }) {		
		unless (exists $anno->{$gene2}{C}{$chr}){
			print "\t->$gene2 don't have $chr\n" if $debug ;
			next;
		}

        my $haps1 = $anno->{$gene1}{C}{$chr};
        my $haps2 = $anno->{$gene2}{C}{$chr};
		my ($start1,$end1,$dir1) = @$haps1;
		my $gl1 = $end1 - $start1 + 1;
		my ($start2,$end2,$dir2) = @$haps2;
		my $gl2 = $end2 - $start2 + 1;
		
		### graph reference supporint informations 
		# Check same direction
		unless ($dir1 eq $dir2) {
			print "\t-> direction not match..\n" if $debug ;
			next;
		}
        my ($dis,$totl) =  intervals_within_range($start1, $end1, $start2, $end2);
		push @ils, $dis;

		if($dis > $range or ( $dis/$gl1 > 10 and  $dis/$gl2 > 10) ){
			print "long distance , skip... $gene1 $gene2 $dis $gl1 $gl2\n" if $debug ;
			next;
		}

		print "\t$gene1 $gene2 distance $dis at chr $chr : @$haps1 = @$haps2\n" if $debug ;
		$cts ++;
		$$shares{$gene1}{$chr} = [$start1,$end1,$dir1];
		$$shares{$gene2}{$chr} = [$start2,$end2,$dir2];
    }

	my $ilm = median (@ils);	
	$ilm = $ilm ? $ilm : 0;
	print "\tend of similar module $cts , $ilm ...\n\n" if $debug ;
	return ($cts, $ilm);
}

sub intervals_within_range{
	my ($s1, $e1, $s2, $e2) = @_;
    # Calculate minimal distance between intervals (overlap is 0 distance)
    if ($e1 < $s2) {
        return ( $s2 - $e1, $e2 - $s1);
    } elsif ($e2 < $s1) {
        return ( $s1 - $e2, $e1 - $s2 );
    } else {
		my @ar = sort ($s1,$s2,$e1,$e2);
		my $tl = $ar[-1] - $ar[0];
		return (1, $tl) ;  # intervals overlap
    }
}

sub cmp_gene_covs {
    my $sum_a = sum(values %{ $gene_covs{$a} // {} }) // 0;
    my $sum_b = sum(values %{ $gene_covs{$b} // {} }) // 0;
    return $sum_b <=> $sum_a;
}

sub merge_gene_hash {
    my ($main, $addon) = @_;

    foreach my $gene (keys %$addon) {

        # If this gene not in main hash, just copy it
        if (!exists $main->{$gene}) {
            $main->{$gene} = $addon->{$gene};
            next;
        }

        # Otherwise merge chr-level info
        foreach my $chr (keys %{ $addon->{$gene} }) {
            $main->{$gene}{$chr} = $addon->{$gene}{$chr};  
            # This overwrites per-chr entry if same chr exists
        }
    }
}






=head

	#### FILTERS in each sample
	# 1, both genes have matched coverage in gene body region
	# 2, the coverge can not have larger  difference
	# 3, overlapped reads in the middle of two genes , can bot be too few.
	#if($g1c == 0 or $g2c ==0){
	#	print "\t\t ----zero cover -1\n" if $debug ;
	#	$decision = -1;
	#}elsif($g1c / $g2c > 3 or $g2c / $g1c > 3){  #### Filter FILTER
	#	print "\t\t ----mismatch coverage. -1 \n" if $debug ;
	#	$decision = -1;
	#}elsif($ovc < $g1c / 2 or $ovc < $g2c / 2){
	#	print "\t\t ----internal too few overlap... -1\n" if $debug ;
	#	$decision = -1;    ###### <<<<<<<<<< FILTER  Filter
	#}else{
	#	print "\t\t ----connect ... + 1\n" if $debug ;
	#	$decision = 1
	#}
		

##################################################
#  Filtering based on connections in population  #
##################################################

print "calculate the weight of each edge.. and filter by population \n";
my %edge_score;
foreach my $g1 (keys %gene_connections){	
	foreach my $g2 (keys %{$gene_connections{$g1}}){
	
		my ($total, $samples) = (0,0);

        foreach my $sam (keys %{$gene_connections{$g1}{$g2}}){
            my $evi = $gene_connections{$g1}{$g2}{$sam};
            $total += $evi;
            $samples++ if $evi > 0 ; ### count the number samples have this connection
        }
		
		my $sam_t =  (scalar keys %assign_fs ) / 4  > 2 ?  (scalar keys %assign_fs ) / 4  : 2;
        if( $samples >= $sam_t  and $total > 0 ) {  ########## FILTER <<<<< 
            $edge_score{$g1}{T} += $total;
            $edge_score{$g1}{S} += $samples;
        }else{
			delete $gene_connections{$g1}{$g2};
			delete $gene_connections{$g2}{$g1};
		}
    }
}
=cut 
