#!/usr/bin/perl
use strict;
use warnings;

my %targ;
my %query;
while(<>){
	chomp;
	my @ar = split /\t/;
	if($ar[-1]){
		$targ{$ar[8]}{$ar[3]}  = $ar[4]  if $ar[-1] ;
		if (!defined $query{$ar[3]}{B}{$ar[8]}) {
			$query{$ar[3]}{B}{$ar[8]} = $ar[-1];
		} elsif ($ar[-1] > $query{$ar[3]}{B}{$ar[8]}) {
			$query{$ar[3]}{B}{$ar[8]} = $ar[-1];
		}
	}
	
	$query{$ar[3]}{C} = $ar[4];
	#print "save $ar[3]  $ar[4]\n";

}
foreach my $k ( sort keys %targ ){
	my $v =  $targ{$k};
	my $ct = 0;
	my @asso;
	foreach my $k1 ( keys %$v ){
		if(  exists $query{$k1} and $query{$k1}{B}{$k} ){
			push @asso, $k1;
			$ct += $query{$k1}{C};
			delete $query{$k1};
			#	print "$k1 deleted\n";
		}
	}
	my $assos = join ";", @asso;
	print "$k\t$ct\t$assos\n" if $ct ;
}


foreach my $id (keys %query ){
	print "__$id\t$query{$id}{C}\n";
}
