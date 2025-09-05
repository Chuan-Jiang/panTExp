#!/usr/bin/perl
use warnings; use strict;
use DBM::Deep;

my $prefix = shift @ARGV;

my $nodenei = DBM::Deep -> new(file => "$prefix.nodenei.db", pack_size => "large");


my %neis;
while(<>){
	chomp;
	if(/^S/){
		if(%neis){
			print "importing new ...\n";
			$nodenei -> import(\%neis);
			undef %neis;
		}
	}elsif(/^W\t(\S+)/){
		my $hap = $1;
		while(/(\d+)(<|>)(\d+)(<|>)(\d+)/g){
			print "$&\n";
			my $pn = $1;
			my $nn = $5;
			my $cn = $3;

			my $d = $2 eq ">" ? 1 : -1;

			if($nn != $cn + $d){
				$neis{$cn * $d} = $nn;
			}
			if($pn !=  $cn - $d){
				$neis{$cn * $d * -1} = $pn;
			}
		}
	}
}
		
