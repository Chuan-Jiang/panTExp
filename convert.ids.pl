#!/usr/bin/perl
use warnings; use strict;


my ($mat_f,$ref,$que_f,$q_col) = @ARGV;


open MAT, $mat_f or die $!;

my $index  ;
if($ref =~ /(.+):(\d+)/){
	my $ref = $1;
	$index = $2;
}else{
	chomp(my $fl = <MAT>);
	my @hds = split /\t/, $fl;

	for ( my $i = 0; $i < @hds ; $i++ ) {
		if ($hds[$i] eq $ref) {
			$index = $i;
			last;
		}
	}
}
##############
unless(defined $index ){
	warn "can not fiond..\n";
	exit;
}else{
	warn "$ref with index $index\n";
}


my %pos;
warn "matrix reading start..\n";
while(<MAT>){
	chomp;
	my @ar = split /\t/;
	
	$pos{$ar[0]} = $ar[$index];
}
close MAT;
warn "matrix reading finished.\n";

my %re;
my $lc = 1;

open QUE, "$que_f" or die $!;
while(<QUE>){
	chomp;
	my $line = $_;
	my @qar = split /\t/, $line;
	my $id = $qar[$q_col];
	my @qs = split /:/, $id;
	my @regions;
	foreach my $q (@qs){
		if($pos{$q} ne "NA"){
			my @s = split /;/, $pos{$q};
			push @regions, @s
		}
	}
	#	print "\nLine:$lc>>> @qs => @regions\n";
	next unless (@regions > 0);

	my $results = merge_regions(\@regions, 1000);
	foreach my $r ( @$results ){
		print "$r->{chr}\t$r->{start}\t$r->{end}\t$qar[0]\t$qar[1]\n"; #$r->{strand}\n";
	}
	$lc ++;
}


sub merge_regions {
    my ($regions_ref, $max_gap, $consider_strand) = @_;
    $max_gap //= 1000;  # Default gap
	$consider_strand //= 0;      # 默认不考虑方向
    
	my @parsed = map {
        /^(\w+):(\d+)-(\d+):([+-])$/ or die "Bad format: $_";
        { chr => $1, start => $2, end => $3, strand => $consider_strand ? $4 : "." }
    } @$regions_ref;

    # Sort by chromosome, strand, and start
    @parsed = sort {
        $a->{chr} cmp $b->{chr} ||
		($consider_strand ? $a->{strand} cmp $b->{strand} : 0) || #$a->{strand} cmp $b->{strand} ||
        $a->{start} <=> $b->{start}
    } @parsed;

    my @merged;
    my $current = shift @parsed;

    foreach my $r (@parsed) {
        if (
            $r->{chr} eq $current->{chr} &&
			(!$consider_strand || $r->{strand} eq $current->{strand}) &&  #$r->{strand} eq $current->{strand} &&
            $r->{start} - $current->{end} <= $max_gap
        ) {
            $current->{end} = $r->{end} if $r->{end} > $current->{end};
        } else {
            push @merged, $current;
            $current = $r;
        }
    }
    push @merged, $current;

    return \@merged;
}

