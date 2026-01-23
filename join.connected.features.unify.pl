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
my $debug = 0;


GetOptions(
	"sam=s" => \$sam_f,
	"dir=s" => \$directory,
	"info=s" => \$info_f,
	"cpu=i" => \$cpu,
	"hpsc=i" => \$hpsc,
	"samc=i" => \$samc,
	"node=s" => \$ninfo_f,
	"debug" => \$debug
) or die "Error in command line arguments\n";

####
my $dist_gmer = 2000;
my $bsize = 1000;

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
my $linkf = "$directory/03_assign/feature.u.links";
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
			my $wei ;
			if($ar[-4] eq "0-0"){
				$wei = 1;
			}else{
				my ($t,$w) = split /-/, $ar[-4];
				$wei = 1/$t;
			}
			$gene_covs{ $ar[0] }{$sam} += $wei / $ar[1] ; # only count overlapped reads;
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

my $total = scalar keys %gene_connections;
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
	my $tabix_n = Bio::DB::HTS::Tabix -> new(filename => $ninfo_f) ;

	#### iterate blocks
	# each bllck is a group of inter-connected genes
	#
	foreach my $bl (@block){
		my ($g_c_l, $e_s_l ) = @$bl;
	
		## these are the two main input hashes to this gene block
		%gene_connections = %$g_c_l;
		%gene_covs = %$e_s_l;
		
		my @batch = sort keys %gene_covs;
		print "\n>>>\nBlock: $g_c_l $e_s_l genes: @batch\n" if $debug ;

		#####################################
		#### collection gene annotations ####
		# at last obtain a annotation hash %anno;
		
		#my %tranges = cluster_ranges(%gposs);
		my %tranges = genes_ranges(@batch);
		
		my %anno;
		print "Step1=> tab annotation for genes\n" if $debug ;
		tab_anno($tabix_i, \%tranges, \%anno);

		##############################
		##### calculate clusters #####
		my @clusters; ## calculate clusters;
		my @gs_sorted = sort cmp_gene_covs  @batch;
		
		print "Step2=> Sort genes in batch according to gene coverage  @gs_sorted\n" if $debug ;
		
		print "Step3=> Clustering and filtering teration init...\n" if $debug;
		my %visited;
		foreach my $gene (@gs_sorted){	
			next if $visited{$gene};
			$visited{$gene} = 1;
			
			#######
			my %cluster;	
			my @gclu = ($gene);	
			
			my %sam_share; # shareed samples, each cluster of genes should share in the same linear componet
			if(defined $gene_covs{$gene}){
				%sam_share = %{$gene_covs{$gene}};
			}else{
				next;
			}

			print "\n\t==>query elements from cluster: @gclu\n"  if $debug ;	
			while(1) {
				# =>>>>First, extract all neighbors. 	
				
				my %neighbors;
				## iteration to find the most connected genes
				# 
				# 1, C : common sampes that alerady shared by gene cluster
				# 2, S : for new gene, all shared samples with any of gene in cluster
				# 3, R : total  number of reads coverting the new gene the gene cluster
				#
				foreach my $g (@gclu){ 
					# iterate each genes that are already clustered.
					next unless ( exists $gene_connections{$g} );
					foreach my $n (keys %{ $gene_connections{$g}}){ ## $n means next connected gene
						next if $visited{$n};
						my $sc = 0;
						foreach my $s ( keys %{$gene_connections{$g}{$n}}){
							if( exists $sam_share{$s} ){
								print "sample $s exists in shared ...\n" if $debug;
								$sc ++;
								$neighbors{$n}{C} ++ ; # total number of supproitn samples (From sequecing samples)( if the final value equal to @gclu,
								                        # this gene share samples with all other clus. 
							}
							$neighbors{$n}{R} += $gene_connections{$g}{$n}{$s}; # total number of supportingt reads (From seuqncing samples)
							$neighbors{$n}{S}{$s} = 1;   # record all share samples for this gene $n.
						}
						unless($sc){
							print "next candidate gene have 0 share sampes, skip  $n\n" if $debug ;
							delete $neighbors{$n};
						}
					}
				}
				
				####################################
				#
				# =>>>>>>> Filteration of neighbors
				#
				sub cmpCR {
					my ($neighbors, $a, $b) = @_;
					return
						  $neighbors->{$b}{C}  <=> $neighbors->{$a}{C} 
					  ||  $neighbors->{$b}{R}  <=> $neighbors->{$a}{R} ;
				}
				my $growth;
				foreach my $nei ( sort  { cmpCR(\%neighbors, $a, $b) }  keys %neighbors ) {
					
					# Filter 1
					if( $neighbors{$nei}{C}   < 1){ # if 0 sequenced samples connect this neighbor, then skip it
						# this one should not happen.
						print "not enough supporint samples... $neighbors{$nei}{C} \n" if $debug ;
						next;  ########### <<<<<<<<<<<< FILTER 
					}

					# Filter 2
					foreach my $g (@gclu){
						my %shares ;  #  shareed haplotype components .   ###########

						##### Filter baed on haplotypes return counts of haps and median length
						my ($sct, $idis)  = genes_similar_haps($gene, $nei, \%anno, $dist_gmer, \%shares ); 	
						if($sct  < $hpsc){  # if two genes have less than $hpsc haplotype supporting, then skip this pair.
							print "\t------>skip = $g and $nei have few haps support $sct\n"  if $debug ;
							goto DEL;
							#next; # <<<<= FILTER ===== shares
						}	
						#### Filters based on sequenced samples;
						my $supp_sams = 0; # number of supporting samples;
						my $supp_sams_rs = 0; # total depth of both genes from supporing samples
						my $supp_sams_rt = 0; # total depth of both genes from all samples;
						foreach my $sam (sort keys %{$gene_connections{$g}{$nei} } ){
							my $g1c = $gene_covs{$g}{$sam} // 0  ;
							my $g2c = $gene_covs{$nei}{$sam} // 0 ;

							my $g1cl = log2($g1c + 1);
							my $g2cl = log2($g2c + 1);
							
							my $diff = abs($g1cl - $g2cl);
							print "connectivity: $g and $nei $sam  $g1c  $g2c  logfc:$diff \n" if $debug ;
							if($g1cl and $g2cl  and $diff < 1 ){ # two genes have comparable expression level, logFC smaller than 1
								$supp_sams ++;
								$supp_sams_rs += ($g1cl + $g2cl); # total number of coverage from support samples;	
							}
							$supp_sams_rt += ($g1cl + $g2cl);  # tatal number of coverage of all samples;
						}

						if( $supp_sams <  $samc  or  $supp_sams_rs / $supp_sams_rt < 0.5 ){  # <<<<<==== FILTER based on population
							print "\t------>skip = $g and $nei have few sampls $supp_sams $supp_sams_rs $supp_sams_rt\n" if $debug ;
							goto DEL;
						}
						#######
						print "\t -----> kepp $g and $nei connected...\n" if $debug ;
						merge_gene_hash( \%cluster, \%shares );
						
						### update shared samples
						%sam_share = %{$neighbors{$nei}{S}};	
						$growth ++;
						$visited{$nei} = 1;
						
						DEL:
						delete $gene_connections{$g}{$nei};
						delete $gene_connections{$nei}{$g};
					}
				}
				if( $growth ){
					@gclu = sort keys %cluster ;
					print "\t\t include new query elements: @gclu\n" if $debug ;
				}else{	
					last;
				}
			}

			my @ginc = sort keys %cluster;
			if( @ginc ){
				push @clusters, [ \%cluster,  \%sam_share ] ;
				$total  = $total - scalar @ginc;
				print "clustering end with gene $gene: @ginc. Left: $total\n" if $debug  ;	
			}else{
				print "\t Empty cluster for gene $gene\n" if $debug ;
			}
		}

		################################
		#enter cluster length  estimate module
		if($ninfo_f){
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

sub cluster_len_node {
	my ($cs, $anno, $tabix_n, $out) = @_;
	
	my @clusters = @$cs;

	print "\n>>>>estimate length of each cluster $#clusters...\n" if $debug ;
	for my $ci ( 0 .. $#clusters ){	
		## estimate length of each cluster;
		my %clus = %{$clusters[$ci][0]};
		my %samp = %{$clusters[$ci][1]};
	
		## find middle genes within a cluster;
		my @genes = sort { 
							my ($a_num) = $a =~ /\|(\d+)$/;
							my ($b_num) = $b =~ /\|(\d+)$/;
							$a_num <=> $b_num;   # numerical sort	
						} keys %clus;


		my %checked ;
		my %genes_exp = map{$_ => 1} @genes;
		my %genes_avoid; 

		foreach my $g (@genes){	
			unless(%checked){
				%checked = %{$clus{$g}};
				next;
			}	
			foreach my $q_chr (keys %{$clus{$g}} ){
				unless($checked{$q_chr}){
					$checked{$q_chr} = $clus{$g}{$q_chr};
					next;
				}
				
				my ( $c_start,$c_end ) = @{$checked{$q_chr}};
				my ( $g_start,$g_end ) = @{$clus{$g}{$q_chr}}[0,1];

				my ($r_start,$r_end ) = (sort {$a <=> $b} ($c_start,$c_end, $g_start,$g_end) ) [0, -1];

				my($q_start,$q_end) = (0, 0);
				if($g_end < $c_start ){
					$q_start = $g_end;
					$q_end = $c_start;
				}elsif($g_start > $c_end){
					$q_start = $c_end;
					$q_end = $g_start;
				}
				
				if($q_end - $q_start > 50){ ############ <- FILTER   
					foreach my $agene (keys %$anno) {
						# skip if chromosome not present
						next unless exists $$anno{$agene}{C}{$q_chr};
						next if $genes_exp{$agene};
						next if $genes_avoid{$agene};

						my ($start, $end, $dir) = @{$$anno{$agene}{C}{$q_chr}};

						# check overlap
						if ($start <= $q_end and  $end >= $q_start) {
							print "found a inter missed gene $agene for group : <@genes>..\n";
							$genes_exp{$agene} = 1;
							$r_start = $start if ($start < $r_start);
							$r_end = $end if ($end > $r_end);
						}elsif($start - $q_end > 5000 or $q_start - $end > 5000 ){
							$genes_avoid{$agene} = 1;
						}
					}
				}
				$checked{$q_chr} = [$r_start,$r_end];
			}
		}
		## expand finished
		###############
		@genes  = sort { scalar keys %{$clus{$b}} <=> scalar keys %{$clus{$a}} } keys %genes_exp;
		my $clus_name = join ":", sort {
											# extract the number after the last |
											my ($a_num) = $a =~ /\|(\d+)$/;
											my ($b_num) = $b =~ /\|(\d+)$/;
											$a_num <=> $b_num;   # numerical sort
										} @genes;
		my $sams = join ":", sort { ($a =~ /(\d+)/)[0] <=> ($b =~ /(\d+)/)[0] }  keys %samp;
		
		### collect all nodes
		my %nodes;
		for my $g (@genes) {
			next unless exists $$anno{$g}{N};
			@nodes{ keys %{ $$anno{$g}{N} } } = values %{ $$anno{$g}{N} };
		}
		
		my %nranges = nodes_ranges(keys %nodes);
		
		my %coors;
		foreach my $r ( keys %nranges ){
			my $iter = $tabix_n -> query($r);
				
			print "query $r..\n" if $debug ;
			while (my $aline = $iter -> next) {
				my @ar = split /\t/, $aline;
				if($nodes{$ar[2]}){
					my($h,$t) = @{$nodes{$ar[2]}};
					my @haps = split /,/, $ar[3];
					foreach my $hap (@haps){
						my($chr,$s,$e) = $hap =~ /(.*):(\d+)-(\d+)/;
						if($e > $s){
							push @{$coors{$chr}}, $s + $h, $e - $t;
						}else{
							push @{$coors{$chr}}, $s - $h, $e + $t;
						}
					}
				}
			}
		}
		my @lens;
		foreach my $chr (keys %coors){
			my @pos = sort {$a <=> $b} @{$coors{$chr}};
			push @lens, $pos[-1] - $pos[0];
		}

		my $len = median(@lens);

		push @$out, "$clus_name\t$len\t$sams\tNode";

	}

}



sub cluster_len {
	my ($cs,$out) = @_;
	my @clusters = @$cs;

	print "\n>>>>estimate length of each cluster $#clusters...\n" if $debug ;
	for my $ci ( 0 .. $#clusters ){
		## estimate length of each cluster;
		my %clus = %{$clusters[$ci][0]};
		my %samp = %{$clusters[$ci][1]};
		
		my @genes = sort { scalar keys %{$clus{$b}} <=> scalar keys %{$clus{$a}} } keys %clus;
		my $clus_name = join ":", sort {
											# extract the number after the last |
											my ($a_num) = $a =~ /\|(\d+)$/;
											my ($b_num) = $b =~ /\|(\d+)$/;
											$a_num <=> $b_num;   # numerical sort
										} @genes;
		my $sams = join ":", sort { ($a =~ /(\d+)/)[0] <=> ($b =~ /(\d+)/)[0] }  keys %samp;

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
		push @$out, "$clus_name\t$len\t$sams\tGene";
		#print "$clus_name\t$len\n";
	}
}


sub log2 {
    my ($x) = @_;
    return log($x) / log(2);
}


sub genes_ranges {
   
	## return the range of gene ids. which is used to extract gene annotatio from *infor.tsv.gz
	
	my @batch = @_;

	my %gposs ; # gene position , that is gene ids
	for my $item (@batch) {
		my @ar = split /\|/, $item;
		push @{ $gposs{$ar[0]} }, $ar[-1];
	}

	my $debug = 0;
	my %tranges;
	foreach my $type (keys %gposs  ){
		my @nums = @{$gposs{$type}};
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

sub nodes_ranges {
	
	my $debug = 0;
	my @nums = @_;

	@nums = sort { $a <=> $b } @nums;
	
	my %ranges;
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
	#
	# extracted information will be save to hash ref: $anno
	#   {C} = location on each linear component genome
	#   {N} = coordications on each occupied nodes
	#
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

		if($dis > $range or ( $dis/$gl1 > 5 and  $dis/$gl2 > 5 ) ){
			print "long distance , skip... $gene1 $gene2 $dis $gl1 $gl2\n" if $debug ;
			next;
		}
		push @ils, $dis;

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
