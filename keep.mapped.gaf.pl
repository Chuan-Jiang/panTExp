#!/usr/bin/perl
use warnings; use strict;

while(<>){
        chomp;
        my @ar = split /\t/;
        if($ar[5] =~ /<|>/){
                print "$_\n";
        }
}

