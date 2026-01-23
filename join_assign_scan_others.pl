#!/usr/bin/perl
use warnings; use strict;
use Bit::Vector;
use File::Basename;

# kfile : complete and keep file
# rfile : raw assign result
# lfile : feature links


my ($kfile, $rfile, $lfile, $pre ) = @ARGV;

###
open OUT, ">$pre.join.assign" or die $! ;
open OTH, ">$pre.others" or die $!;


my %links;
open LINK, $lfile or die $!;
while(<LINK>){
	chomp;
	my ($id,$len) = split /\t/;
	my @ids = split /:/, $id;
	my $name = $ids[0].".J";

	foreach my $d (@ids){
		$links{$d} = "$name\t$len";
	}
}
close LINK;


my %remains;
my %seens;

open IN, $rfile or die $!;
while(<IN>){
	chomp;
	my @ar = split /\t/;
		
	my $r =  $ar[2];

	next if ( $seens{ "$r:$ar[-4]" } );

	####
	if ( $links{$ar[0]} ){
		@ar[0,1] = split /\t/, $links{$ar[0]};
	}

	print OUT (join "\t", @ar), "\n";
	$seens{ "$r:$ar[-4]" } = 1;


	#####################
	
	my @muls = split /-/, $ar[-4];
	
	if($muls[0] < 2){
		next;
	}else{
		#print "@muls\n";
	}

	unless ( defined $remains{ $r } ){
		#print "creat a vector $muls[0] + 1 \n";
		my $vec = Bit::Vector -> new( $muls[0] + 1 );
		$vec -> Bit_On(0);
		$remains{ $r } = $vec;

	}
	
	if ($muls[1]  < 0 || $muls[1]  >= $remains{$r}->Size()) {
		my $ts = $remains{$r} -> Size();
		die "Bit_On index out of range: @ar $muls[0] $muls[1] $ts\n";
	}

	$remains{ $r }  -> Bit_On($muls[1]);

	if( $remains{$r} -> is_full ){
		delete $remains{$r};
	}
}
close IN;

open KIN, $kfile or die $!;

while(<KIN>){
	chomp;
	my @ar = split /\t/;
	my($read,$mul) = $ar[6] =~ /(.+):(.+)/;
	
	$mul = eval( $mul );
	if(exists $remains{$read}){
		if(! $remains{$read} -> bit_test($mul) ){
			print OTH "$_\n";
		}
	}
mecfs_pan/02_complete/CF1.complete.keep.tsvmecfs_pan/02_complete/CF1.complete.keep.tsv}
