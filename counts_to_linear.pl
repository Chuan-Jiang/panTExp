#!/usr/bin/perl
use warnings; use strict;

my %tars;
my %over;

while(<>){
	chomp;
	my @ar = split /\t/;
	my $pid = $ar[4];
	my $sid = $ar[9];

	my $ct = $ar[5];
	next unless $ct;

	my $ol = $ar[12];

	$tars{$sid}{$pid} = $ct;
	$over{$pid}{O}{$sid} = $ol;
	$over{$pid}{C} = $ct;
	#	print "$pid ==> $ct\n";
}

foreach my $id (sort keys %tars){
	my $c = 0;
	# to check if the match id  is the longest match for it
	foreach my $id2 ( keys %{$tars{$id}} ){
		my @h1 = sort { $over{$id2}{O}{$b} <=> $over{$id2}{O}{$a} } keys %{ $over{$id2}{O} }; 
		unless (@h1){
			delete $over{$id2};
			next;
		}
		if($id eq $h1[0]) {
			$c += $tars{$id}{$id2};
			delete $over{$id2}; 
		}
	}
	print "$id\t$c\n" if $c > 0 ;
}

foreach my $id (keys %over ){
	print "__$id\t$over{$id}{C}\n";
}

