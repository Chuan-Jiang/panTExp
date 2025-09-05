#!/usr/bin/perl
use warnings; use strict;

my $prefix = shift;
open H, "| bgzip > $prefix.Overlap.tsv.gz" or die $!;
open N, "| bgzip > $prefix.NoOverlap.tsv.gz" or die $!;

while(<>){
	chomp;
	my @ar = split /\t/;

	unless ($ar[8] eq "G"){
		print N "$_\n";
	}else{
		print H "$_\n";
	}
}

