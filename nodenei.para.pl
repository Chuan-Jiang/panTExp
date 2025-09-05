#!/usr/bin/perl
use warnings; use strict;
use Storable;
use Parallel::ForkManager;
use IO::Handle;
use File::Temp qw(tempfile);

my ($prefix,$sample,$hap)  = @ARGV;

if($sample){
	$prefix = "$prefix.$sample.$hap";
}

#open NEI, ">$prefix.nodenei.tsv" or die $!;

my $num_processes = 5;
my $pm = Parallel::ForkManager->new($num_processes);

##### create temp files

my @tempfiles;
my %fh_assign;
for my $i (0 .. $num_processes ) {
    my ($fhk, $filenamek) = tempfile("childnei_$i-XXXX", UNLINK => 1, DIR => ".");
    push @tempfiles, { fhk => $fhk, filek => $filenamek};
}


###### clean subprocess 
$pm->run_on_finish(sub {
    my ($pid, $exit_code, $ident, $exit_signal, $core_dump) = @_;
    if ($exit_signal) {
        print "Child $pid died from signal $exit_signal\n";
    }
    print "file handle $ident released by $pid..\n";
    delete $fh_assign{$ident};
});


my %nodelen;

my %neis;
my @ws;
while(<STDIN>){
	chomp;
	if(/^S\t/){
		if(@ws){
			wline(\@ws);
		}
		my @ar = split /\t/;
		my $l = length($ar[2]);
		#die "$_" unless $l;
		next unless $l;
		#print "$ar[0]\t$ar[1]\t$ar[1]\t$l\n";
		$nodelen{$ar[1]} = $l if $l > 1;
	}elsif(/^W\t/){
		push @ws, $_;
	}
}

wline(\@ws) if (@ws);

store \%nodelen, "$prefix.nodelen.dat" ;

$pm->wait_all_children;

#### summarize all outputs froom child
my @files = map { $_ -> {filek} } @tempfiles;
#system("cat @files  | sort -k 2,2n -k 3,3n -k 4,4 --parallel=10 --buffer-size=30G  > $output");             
system("cat @files |  sort -k 2,2n --parallel=$num_processes | bgzip -@ $num_processes - > $prefix.nodenei.gz  ");
system("tabix -p bed $prefix.nodenei.gz  ");


########### sub funciton 
sub wline {
	my $wr = shift @_;
	my @wsi = @{$wr};
	@{$wr} = ();
	
	warn "yes.. start a new child\n";
	my ($fhi, $fhk) = getfh(\@tempfiles,\%fh_assign) ;
	$pm -> start($fhi) and return;
	warn  "child is working..\n";
	my %neis;
	foreach my $l (@wsi){
		my @ar = split /\t/, $l;
		if($sample){
			unless($sample eq $ar[1] and $hap eq $ar[2]){
				next;
			}
		}
		#print "@ar[0..5]\n";
		my @nodes = $ar[6] =~ /(<|>)(\d+)/g;
		for(my $i = 1; $i < @nodes; $i += 2){
			my $n = $nodes[$i];
			my $d = $nodes[$i - 1] eq ">" ? 1 : -1;
			#print  "Node $n\n";
			### current , up and down are all in the same DNA string define by current node
			#. not the reverse complete string
			if($i >=3 ){
				my $pn = $nodes[$i - 2];
				my $pd = $nodes[$i - 3] eq ">" ? 1  : -1;
				$pd = $pd * $d;	
				## the second key indicate its up(-1)/down(1) stream nodei relative the to current node
				$neis{$n}{ -1 * $d }{$pn*$pd} = 1;
			}
			if($i < @nodes - 2){
				my $nn = $nodes[$i+2];
				my $nd = $nodes[$i+1] eq ">" ? 1  : -1;
				$nd = $nd * $d ;
				$neis{$n}{ 1 * $d }{$nn*$nd} = 1;	
			}
		}
	}

	#print "importing new ...\n";
	my @outs;
	foreach my $n (  keys %neis){
		my $pn;
		if( defined $neis{$n}{-1}){
			my @sn = sort { abs($a) <=> abs($b) }  keys %{$neis{$n}{-1}};
			#@sn = map {$nodelen{abs($_)} ? "$_:$nodelen{abs($_)}" : "$_:1"} @sn;
			$pn = join ";", @sn; 
		}else{
			$pn = "NA";
		}
		my $nn; 
		if( defined  $neis{$n}{1} ){
			my @sn = sort { abs($a) <=> abs($b) } keys %{$neis{$n}{1}};
			#@sn = map {$nodelen{abs($_)} ? "$_:$nodelen{abs($_)}" : "$_:1"} @sn;
			$nn = join ";",  @sn;
		}else{
			$nn = "NA";
		}
		my $len = $nodelen{$n} ? $nodelen{$n} : 1;
		#print "G\t$n\t$n\t$pn\t$nn\t$len\n";  
		my $n1 = $n - 1;
		#push @outs, "G\t$n1\t$n\t$pn\t$nn\t$len";
		print $fhk "G\t$n1\t$n\t$pn\t$nn\t$len\n"; 
	}
	$fhk -> flush;
	$pm -> finish;
}
sub getfh {
    my ($tempfiles,$fh_assign) =   @_;

    my $fhi;
    while(! defined $fhi){
        for(my $i = 0 ;$i < @$tempfiles; $i ++){
            unless (  defined $$fh_assign{$i} ){
                $fhi = $i;
                $$fh_assign{$i} = 1;
                last;
            }
        }
    }
    my $fhk;
    $fhk = $$tempfiles[$fhi] -> {fhk};
    return ($fhi, $fhk);
}

