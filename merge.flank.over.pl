#!/usr/bin/perl
use warnings; use strict;


while(<STDIN>){
	chomp;
	my @ar = split /\t/;
	my $gene = pop @ar;

	my @gar = split /\|/, $gene;
	my ($gs,$ge) = ($gar[2] - 1, $gar[2]);

	print join ("\t",$gar[0],$gs,$ge,$gene, @ar), "\n";
}

