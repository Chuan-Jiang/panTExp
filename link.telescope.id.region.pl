#!/usr/bin/perl
use strict;
use warnings;

my($report,$bed) = @ARGV;

my %cts;
open REP, $report or die "Cannot open $report: $!";
while(<REP>){
	chomp;
	if(/^#/ || /^transcript/ || /^__no_feature/ ){
		next;
	}

	my @fields = split(/\t/);
	$cts{$fields[0]} = $fields[2];
}
close REP;

open BED, $bed or die "Cannot open $bed: $!";
while(<BED>){
	chomp;

	my @fields = split(/\t/);
	my ($gid) = $fields[8] =~ /gene_id "(.+?)";/;
	if($cts{$gid}){
		print join ("\t", $fields[0], $fields[3] --, $fields[4], $gid, $cts{$gid}), "\n";
	}
}
