#!/usr/bin/perl

use strict;
use warnings;

if ($#ARGV < 0){
	print "Usage: $0 <fof of gzipped fastq files>\n";
	exit 1;

}

my $fof = $ARGV[0];
my $indices;

open(FOF, $fof) || die "Can't open $fof...\n";
while(<FOF>){
	chomp;
	readFastq($_);
}
close FOF;

foreach my $i (sort {$indices->{$b} <=> $indices->{$a}} keys %$indices) {
	my $mult = $indices->{$i};
	print "$i,$mult\n";
}


sub readFastq {
	my $readFile = shift;
	my $cmd = "gunzip -c $readFile|";
	open(READS, $cmd) || die "Can't run $cmd...\n";
	while(<READS>){
		if(/^\@\S+\s+BX:Z:(\S+)/){
			my $ind = $1;
			$indices->{$ind}++;
		}
	}

	close READS;

}
