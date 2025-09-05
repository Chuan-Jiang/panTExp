#!/usr/bin/perl
use warnings; use strict;

while(<>){
	s/GRCh38#0#//g;
	if(/CHM13/){
		next;
	}elsif(/AS:i:\d+/){
		print ;
	}elsif(/^@/){
		print;
	}
}

