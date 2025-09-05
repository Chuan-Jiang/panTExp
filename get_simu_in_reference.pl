#!/usr/bin/perl
use warnings; use strict;

while(<>){
	chomp;
	my @ar = split /\t/;
	my $per = ($ar[6] - $ar[5])/$ar[4];
	#if(($ar[6] - $ar[5])/$ar[4] > 0.9){
		#$ha{$ar[0]} = "$ar[1]:$ar[8]-$ar[9]";
		print "$ar[0]\t$ar[1]:$ar[8]-$ar[9]\t$per\t$ar[2]\n";
	#}
}


