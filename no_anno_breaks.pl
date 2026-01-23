#!/usr/bin/perl
use warnings; use strict;

my $brks_file = shift @ARGV;

open OUT, ">$brks_file" or die $!;

my $last_node ;
while(<>){
	chomp;
	my @ar = split /\t/;
	my $curt_node = $ar[1] + 1;
	if($last_node and $curt_node  > $last_node){
		print OUT "$last_node\t$curt_node\n";
		#print  "$last_node\t$curt_node\n";
	}
	$last_node = $ar[2];
}

