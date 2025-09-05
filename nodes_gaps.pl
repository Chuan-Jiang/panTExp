#!/usr/bin/perl
use warnings; use strict;


while(<>){
	chomp;
	if(/^W\t(\S+)/){
		my $hap = $1;
		while(/(\d+)(<|>)(\d+)/g){
			pos() -= length($3);
			#print "GGGG $&\n";
			my $f = $1;
			my $s = $3;
			my $d = abs($f - $s);
			if($d > 20){
				print "$hap\t$d\t$f\t$s\n";
			}
		}
	}
}
		
