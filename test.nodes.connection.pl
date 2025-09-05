#!/usr/bin/perl
use warnings; use strict;


while(<>){
	chomp;
	if(/^W\t(\S+)/){
		my $hap = $1;
		if($_ =~ /(72243280.{0,8000}73694032)/){
			print "$hap\t$1\n";
		}elsif($ _ =~ /(73694032.{0,8000}72243280)/){
			print "$hap\t$1\n";
		}
	}
}
		
