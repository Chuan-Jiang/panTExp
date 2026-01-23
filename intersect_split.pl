#!/usr/bin/perl
use strict;
use warnings;
use Bit::Vector;
use List::Util qw(min max);
use Getopt::Long;

my %opt;
GetOptions(
    "overlap|o=s"     => \$opt{overlap},
    "noverlap|n=s"  => \$opt{noverlap},
    "experiment|e=s"  => \$opt{experiment},   # e.g., f or r
) or die "Error in command line arguments\n";

# Check required args
die <<"USAGE" unless defined $opt{overlap} ;

Usage:
  $0 --overlap FILE --noverlap FILE  --experiment f|r

Example:
  $0 --overlap overlap.tsv --noverlap non.tsv --experiment f
USAGE

my $overlap_file    = $opt{overlap};
my $noverlap_file = $opt{noverlap} // ""; ;
my $experiment      = $opt{experiment} // '';   # optional

#############################################
#############################################


my $overlen = 30;

my $SD = 0;
if($experiment eq "f"){
	$SD = 1;
}elsif($experiment eq "r"){
	$SD = -1;
}

open my $fh_overlap, ">", $overlap_file or die "Cannot open $overlap_file: $!";

open my $fh_non,    ">", $noverlap_file or die "Cannot open $noverlap_file: $!" if $noverlap_file ;

# Disable autoflush for speed

my %highcov;
my %gvec;

while (<STDIN>) {
	#	print;
    chomp;
	my @fields = split /\t/, $_ ;

	## remove redundant
	my $high = determine_repeat($fields[4],$fields[5],$fields[6],\%highcov);
	if($high){
		#print "highcov: $_\n";
		next;
	}

	my $lout = join "\t", @fields[0..$#fields-5];
    my $anno = $fields[-2];
	my $apath = $fields[-1];

    if ($anno eq "." ) {
        print $fh_non  "$lout\tM\t.\n" if $noverlap_file ;
		undef %gvec;
		#print "MISS $lout\n";
    } else {

		unless ( exists $gvec{$anno}){
			$gvec{$anno} = get_gene_vec ( $apath );
		}

		### test strand direction
		my $score = 0 ;
		#my @rns = $fields[5] =~ /(?:^|;)([^,]+)/g;
		my @rns = split /;/, $fields[5];
		foreach my $r (@rns) {
			my ($rid, $rd, $rl, $rs, $re, $rt)  = split /,/, $r;
			#print "$r = $rid $rd $rl $rs $re $rt...\n";
			if( exists $gvec{$anno}{$rid} ){
				my $dd = $SD * $rd ;
				$dd = 0 if exists $gvec{$anno}{$rid}{0};

				if( exists $gvec{$anno}{$rid}{$dd}   ){
					my @gc = @{$gvec{$anno}{$rid}{$dd}};	
					my @rc = ( $rs , $rl - $re - 1 );
					my $overlap = 0;
					if ($rc[1] >= $gc[0] and $gc[1] >= $rc[0] ){
						$overlap = min( $gc[1], $rc[1] ) - max( $gc[0], $rc[0] ) + 1;
					}
					#print "$anno : $rid :  $rd $dd $overlap\n";
					$score += $overlap;
				}
			}
		}

		#print "total score $score\n";
		if($score <  $overlen ){
			print $fh_non  "$lout\tA\t$anno\n" if $noverlap_file ;
			# direction mismatch, discard this annoation
		}else{
			print $fh_overlap "$lout\tO\t$anno\n" ;
		}
	}
}

close $fh_overlap;
close $fh_non if $noverlap_file ;

sub get_gene_vec{
	my ( $path ) = @_;
	my %vecf;
	
	my @nodes = split /;/, $path;
	foreach my $n ( @nodes ){
		my($id,$d, $l,$s,$e) =  split /,/, $n ; #=~ /(.)(\d+):(\d+),(\d+),(\d+)/;

		$vecf{$id}{$d} = [$s,$l - $e - 1];
	}

	return (\%vecf);
}

sub determine_repeat{
    my($len, $path,$read,$highcov) = @_;

    #print "highcovB $path\n";
    if ($path =~ /;/) {
        my @parts = split /;/, $path;
        $path = join ";", $parts[0], $parts[-1] if @parts > 2;
    }

    #print "highcovA $path\n";

    if(exists $$highcov{"$path=$len"} ){
        if(scalar keys %{ $$highcov{"$path=$len"} } > 3){
			return 1;
		}
	}else{
        undef %$highcov ;
	}

	$$highcov{"$path=$len"}{$read} ++;
		
    return 0;
}
