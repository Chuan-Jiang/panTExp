#!/usr/bin/perl
use warnings; use strict;
use Storable;

my %nodelen;
while(<TAB>){
	chomp;
	my @ar = split /\t/;
	$nodelen{$ar[1]} = $ar[3] if $ar[3] > 1;
}

#store \%nodelen, "nodelen.dat";

