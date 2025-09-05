#!/usr/bin/perl
use warnings; use strict;
use Storable;


my $debug = 0;
my $strand_type = "";
my @strand;
if($strand_type eq "f"){
    @strand = (0,1,-1);
}elsif($strand_type eq "r") {
    @strand = (0,-1,1);
}else{
    @strand = (0,0,0);
}   

my $OVERLAP = 0;
my $COVERAGE = 1;


#\%inside, \%overhang, \%outside
#my $hrs = retrieve('anchors.dat');


#my %n2r = %{retrieve('n2r.dat')}; #%{$$hrs[0]};
my %r2g = %{retrieve('r2g.dat')};#%{$$hrs[1]};
my %f2r = %{retrieve('f2r.dat')}; #%{$$hrs[2]};
#################



foreach my $gene (sort keys %f2r){
	unless( exists $f2r{$gene} ){
		print "Gene $gene have encountered before.next...\n" if $debug;
		next;
	}
	
	print "\nInitial checking $gene\n" if $debug;	
	## for each genes, extract all mapped reads , and then extract suspending nodes from these readsi
	## and collect them into %nodes_ini		


	my %genes_asso;
	my %reads_asso;

	my @gs = ($gene);
	while(@gs){	
		foreach my $g (@gs){
			if($genes_asso{$g}){
				print "\t$g have checked\n" if $debug;
				next;
			}else{
				$genes_asso{$g} = 1;
			}
		}

		my @gs_for = @gs;
		@gs = ();

		foreach my $g (@gs_for){
			print ">>>walking $g \n" if $debug;
			my ($ghr) = gene_walk($g,\%reads_asso);	## keep nodes expanded outside of gene	
			@gs = sort keys %{$ghr};
			print "<<<Found Neibhour genes\t@gs\n" if $debug;
		}
	}

	## export genes and reads
	my @ids = keys %genes_asso;
	my $ids_s = join ":", @ids;

	my @rds = keys %reads_asso;
	foreach my $r (@rds){
		my ($rr,$ri) = split /##/,$r;
		my $as = $reads_asso{$r}{1} + $reads_asso{$r}{2};
		print "$ids_s\t$rr\t$ri\t$as\n";
	}
}

sub node_2_gene{

}


sub gene_walk{
	my($g,$reads) = @_;

	print ">>>>>>>In sub-fucntion gene_2_node $g \n" if $debug;	
	
	unless($f2r{$g}){
		print "\tgene $g have seen next\n" if $debug;
		return;
	}
			
	#my %nodes; ### nodes associated with genes
	my %genes ; #### genes, need to be return
	
	##################

	my @rs_g = sort keys %{$f2r{$g}};
	print "\t>>>$g @rs_g\n" if $debug ;
	delete $f2r{$g};
	foreach my $r ( @rs_g ){
		print "\tGN:checking $r\n" if $debug;
		if ($$reads{$r}){
			print "\t reads seen, next..\n " if $debug;
			next;
		}

		if(exists $r2g{$r}{1} and exists $r2g{$r}{2}){
			print "\tGN:Pair mapped \n" if $debug;	
			my $olen = 0;
			$olen += $r2g{$r}{2}{AN}{$g} if $r2g{$r}{2}{AN}{$g} ;
			$olen += $r2g{$r}{1}{AN}{$g} if $r2g{$r}{1}{AN}{$g} ;
			unless($olen >= $OVERLAP){
				next;
			}else{
				print "\ttotal len: $olen\n" if $debug;
			}

			foreach my $which ( 1,2 ){
				## collect other genes
				foreach my $og (keys %{$r2g{$r}{$which}{AN}}){
					if($og eq $g or $og eq "N" ){
						next;
					}else{
						print "\tfound new $og\n" if $debug;
						my $oolen = 0;
						$oolen += $r2g{$r}{2}{AN}{$og} if $r2g{$r}{2}{AN}{$og} ;
						$oolen += $r2g{$r}{1}{AN}{$og} if $r2g{$r}{1}{AN}{$og} ;
						if($oolen >= $OVERLAP){
							$genes{$og} ++;
							print "\tkeep\tover $oolen\n" if $debug;
						}else{
							print "\tskip\tover $oolen\n" if $debug;
						}

					}
				}

				my $rd = $strand[$which] ; # read direcition
				print "\tGN:ALN: $which " if $debug ;
				$$reads{$r}{$which} = $r2g{$r}{$which}{AS};
				
				if($debug){
					my @ans = keys %{$r2g{$r}{$which}{AN}};
					print "\tANNOS:@ans\n" if $debug;
				}
		
				#### find overhang reads
				### overlap length
=head
				unless($r2g{$r}{$which}{NS}){
					print "\tGN:no overhang \n";
					next;
				}else{
					my %sus = %{$r2g{$r}{$which}{NS}} ;
					print "\tGN:suspend \n";
					if($r2g{$r}{$which}{AN}{$g}){ ### overlap might be 0.
						delete $r2g{$r}{$which}{AN}{$g};
					}
					###########
					my @nodes = keys %sus;
					foreach my $n (@nodes){
						my $nd = $sus{$n}[0];	
						my($h,$t) = ($sus{$n}[2],$sus{$n}[3]);
						push @{$nodes{$n}{$rd*$nd}{H}}, $h;
						push @{$nodes{$n}{$rd*$nd}{T}}, $t;	
						$nodes{$n}{L} = $sus{$n}[4];
					}
					print "\t\tGN:NODES: @nodes\n" if $debug;
				}
=cut

			}
		}else{
			print "\t\tGN:No paired $r, skip this read ..\n" if $debug;
			next;
		}
	}
	return (\%genes);
}

=head
	# walk to another genes/nodes based on nodes
	my %nodes_w;
	foreach my $n ( sort keys %nodes){
		my @reads_have_n;
		if($n2r{$n}){
			@reads_have_n = unique(@{$n2r{$n}});
			delete $n2r{$n};
		}else{
			print "have seen node $n\n" if $debug;
			next;
		}
		# reads associated with nodes
		
		foreach my $r(@reads_have_n){
			if($$reads{$r}){
				print "\thave seen read $r $n \n" if $debug;
				next;
			}
			print "\tMAIN: New read $r\n" if $debug;
			if(exists  $r2g{$r}{1} and exists  $r2g{$r}{2}){
				for my $w ( 1, 2){
					if($r2g{$r}{$w}{AN}){ # step into another annotation
						my @gs = keys %{$r2g{$r}{$w}{AN}};
						print "\tMAIN:$r $w overlapes with @gs\n" if $debug;
						foreach my $g (@gs){
							if($g eq "N"){
								my %sus_r = %{$r2g{$r}{$w}{NS}} ;
								my %nodes_r;
								my $rd = $strand[$w];
                 foreach my $n (keys %sus_r){
                     print "\t>$n\n" if $debug;
                     my $nd = $sus_r{$n}[0];
                     my($h,$t) = ($sus_r{$n}[2],$sus_r{$n}[3]);

                     push @{$nodes_r{$n}{$rd*$nd}{H}}, $h;
                     push @{$nodes_r{$n}{$rd*$nd}{T}}, $t;
                     $nodes_r{$n}{L} = $sus_r{$n}[4];

                 }
                 print "\tMAIN:extract nodes from $g \n" if $debug;
                 mergehash(\%nodes_ini , \%nodes_r );
             }else{
                 print "\tMAIN:second expand gene $g\n" if $debug;
                 my $node_g = gene_nodes($g,\%seen,\%genes_asso);
                 mergehash(\%nodes_ini , $node_g);
             }
         }
     } 
=cut


sub mergehash{
	my ($hr1,$hr2) = @_;
	foreach my $n (keys %$hr2){
		if(defined $$hr1{$n}){
			foreach my $d  ( keys %{$$hr2{$n}} ){
				if(defined $$hr1{$n}{$d}){
					push @{$$hr1{$n}{$d}{H}}, @{$$hr2{$n}{$d}{H}};
					push @{$$hr1{$n}{$d}{T}}, @{$$hr2{$n}{$d}{T}};
				}else{
					$$hr1{$n}{$d} = $$hr2{$n}{$d};
				}
			}
		}else{
			$$hr1{$n} = $$hr2{$n};
		}
	}
}

sub unique {
	my @input = @_;
	my %seen;
	my @unique = grep { !$seen{$_}++ } @input;
	return @unique;
}
