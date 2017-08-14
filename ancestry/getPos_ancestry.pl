#!/usr/bin/env perl

use strict;
use warnings;

use IO::File;

my $raw_file = $ARGV[0];

unless (-e $raw_file) {
    die "Could not locate raw file: $raw_file";
}

my $fh;
if ($raw_file =~ m/zip$/) {
    $fh = IO::File->new("zcat $ARGV[0]|");
} else {
    $fh = IO::File->new($ARGV[0]);
}

my $ref_path = $ARGV[1];


while (my $line = $fh->getline) {
    chomp $line;
    if (substr($line, 0, 1) eq '#') {
        next;
    }
    # ancestry mod: skip header
    if (substr($line, 0, 4) eq 'rsid') {
        next;
    }
    # ancestry mod: modify header
    my ($rsid, $chr, $pos, $a1, $a2) = split /\t/, $line;

    #change mitochondrial MT->M, 23->X, 24->Y, skip chr25 & chr26
    if ( ($chr eq '25') || ($chr eq '26') ) {
        next;
    }
    elsif ($chr eq 'MT') {
        $chr = 'M';
    }
    elsif ($chr eq '23') {
        $chr = 'X';
    }
    elsif ($chr eq '24') {
        $chr = 'Y';
    }
    # else {
    #     $chr = 'chr'.$chr
    # }

    my $faidx_string = "samtools faidx $ref_path chr$chr:$pos-$pos";

    my @result = `$faidx_string`;

    my $ref = $result[1];
    chomp $ref;
    print "$chr\t$pos\t$rsid\t$ref\n";
}
$fh->close;
