#!/usr/bin/perl
use warnings; use strict;
use File::Path qw(make_path);
use File::Basename;

my ($spacer_f,$prefix, $batch) =  @ARGV;

$batch ||= 1000000;

my $dir = dirname($prefix);
unless(-d $dir){
	make_path($dir) or die "Cannot create path $dir: $!";
}


my %snds;
open my $fh, "$spacer_f" or die "Cannot open $spacer_f: $!";
while(<$fh>){
	chomp;
	$snds{$_} = 1;
}
close $fh;


my $split = 1;
open my $out , "> $prefix.$split.tsv" or die "Cannot open $prefix.$split.tsv: $!";


my $last_line = <STDIN>;
chomp $last_line;
my $last_node = (split /\t/,$last_line)[2];
my $line_ct = 1;
while(<STDIN>){
	chomp;
	$line_ct ++;
	my @ar = split /\t/;
	my $curr_node = $ar[1] + 1;
	if($curr_node > $last_node and $line_ct > $batch  ){
		my @hits = grep { $_ > $last_node && $_ < $curr_node } keys %snds;
		if(@hits  ){
			print $out $last_line, "\n";
			print $out $_, "\n";
			$split++;
			close $out;
			$line_ct = 0;
			warn "create a new split $split \n";
			open $out , "> $prefix.$split.tsv" or die "Cannot open $prefix.$split.tsv: $!";
		}	
	}else{
		print $out $last_line, "\n";
	}
	$last_line = $_;
	$last_node = $ar[2];
}

print $out "$last_line\n";



