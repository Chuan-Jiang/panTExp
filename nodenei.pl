#!/usr/bin/perl
use warnings; use strict;
use Storable;
use Parallel::ForkManager;

my $max_processes = 4;
my $pm = Parallel::ForkManager->new($max_processes);
#use DBM::Deep;

# Tie the hash to a disk-based hash
#my $hash = DBM::Deep->new('hash_data.db');

# Use it like a normal hash
#$hash->{key1} = 'value1';
#$hash->{key2} = 'value2';

my $prefix = shift @ARGV;

my %nodelen;
open NEI, ">$prefix.nodenei.tsv" or die $!;
#my $nodenei = DBM::Deep -> new(file => "$prefix.nodenei.db",pack_size => 'large');

my %neis;
while(<>){
	chomp;
	if(/^S/){
		if(%neis){
			print "importing new ...\n";
			foreach my $n ( sort { $a <=> $b } keys %neis){

				my $pn;
				if( defined $neis{$n}{-1}){
					$pn = join ";", (sort {$a <=> $b}  keys %{$neis{$n}{-1}} );
				}else{
					$pn = "NA";
				}
				my $nn; 
				if( defined  $neis{$n}{1} ){
					$nn = join ";", ( sort {$a <=> $b} keys %{$neis{$n}{1}} );
				}else{
					$nn = "NA";
				}
			
				print NEI "G\t$n\t$n\t$pn\t$nn\n";
			}
			%neis  = ();
		}

		my @ar = split /\t/;
		my $l = length($ar[2]);
		#print "$ar[0]\t$ar[1]\t$ar[1]\t$l\n";
		$nodelen{$ar[1]} = $l if $l > 1;
	}elsif(/^W/){
		my @ar = split /\t/;
		print "@ar[0..5]\n";
		my @nodes = $ar[6] =~ /(<|>)(\d+)/g;
		for(my $i = 1; $i < @nodes; $i += 2){
			my $n = $nodes[$i];
			my $d = $nodes[$i - 1] eq ">" ? 1 : -1;
			#print  "Node $n\n";
				
			if($i >=3 ){
				my $pn = $nodes[$i - 2];
				my $pd = $nodes[$i - 3] eq ">" ? 1  : -1;
				$pd = $pd * $d;	
				$neis{$n}{ -1 * $d }{$pn*$pd} = 1;
			}
			if($i < @nodes - 2){
				my $nn = $nodes[$i+2];
				my $nd = $nodes[$i+1] eq ">" ? 1  : -1;
				$nd = $nd * $d;
				$neis{$n}{ 1 * $d }{$nn*$nd} = 1;
				
			}
		}
	}
}

store \%nodelen, "$prefix.nodelen.dat";
