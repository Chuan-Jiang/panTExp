#!/usr/bin/perl
use warnings; use strict;
use Storable;

my $debug = 1;
my $strand_type = "";
my @strand;
if($strand_type eq "f"){
    @strand = (1,-1);
}elsif($strand_type eq "r") {
    @strand = (-1,1);
}else{
    @strand = (0,0);
}

#\%inside, \%overhang, \%outside
my $hrs = retrieve('anchors.dat');

my %n2r = %{$$hrs[0]};
my %r2g = %{$$hrs[1]};
my %f2r = %{$$hrs[2]};

my %genes_cluster;
foreach my $g_ini (sort keys %f2r){
    print "\nFor checking $g_ini\n" if $debug;
    ## and collect them into %nodes_ini
    print ">>>colecting suspending nodes for $g_ini\n" if $debug;

	my @gs = ($g_ini);
	while(1){
		
		my @gs_for = @gs;
		@gs = ();
		foreach my $g ( @gs_for ){
			foreach my $r (keys %{$f2r{$g}}){
				print "\n\tchecking $r\n" if $debug;
				my($rn,$mapn) = split /##/, $r;
				if(defined $r2g{$r}{0} and defined $r2g{$r}{1}){
					my $as = $r2g{$r}{0}{AS} + $r2g{$r}{1}{AS};
					my ($rn,$rm) = split /##/, $r;
					$genes_cluster{$g}{$rn}{$rm} = $as;
				
					for (my $i = 0; ; $i++){
						next if $i == $mapn;

						if(defined $r2g{"$rn


    extract_reads($gene,\%gene_reads);
}


sub extract_reads{
    my($gene,$hr) = @_;

    foreach my $r ( keys %{$f2r{$gene}} ){
        print "\n\tchecking $r\n" if $debug;
        my($rn,$mapn) = split /##/, $r;
        if(defined $r2g{$r}{0} and defined $r2g{$r}{1}){
            my $as = $r2g{$r}{0}{AS} + $r2g{$r}{1}{AS};
			my ($rn,$rm) = split /##/, $r;
			$$hr{$gene}{$rn}{$rm} = $as;
		}
	}

}



