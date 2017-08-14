#!/usr/bin/env perl

use strict;
use warnings;
use IO::File;

#input
my $unphased_path = $ARGV[0];
my $ref_alt_path = $ARGV[1];

# # extract filename from unphased_path
# my @a = split('/', $unphased_path);
# my $b = $a[-1];
# my $filename = substr($b, 0, -4);

#open file
my $unphaed_fh = IO::File->new($unphased_path);
my $ref_alt_fh = IO::File->new("gunzip -c $ref_alt_path |");

#load 1000GP Phase3 ref_alt into memory
my %ref_alt = ();
while (<$ref_alt_fh>) {
	chomp $_;
	my ($chr, $pos, $ref, $alt) = split /\t/, $_;
	$ref_alt{$chr}{$pos} = {
		ref => $ref,
		alt => $alt,
	};
}
close $ref_alt_fh;

#read chrN.unphased.vcf line by line to check ref match and replace alt if necessary
my $totalC = 0;
my $missingC = 0;
my $aberrantC = 0;
my $replacedC = 0;
my $passC = 0;
while (<$unphaed_fh>) {
	chomp $_;
	if (substr($_, 0, 1) eq '#') {
		print "$_\n";
        next;
    }

    $totalC += 1;
	my ($chr, $pos, $rsid, $ref, $alt, $qual, $filter, $info, $format, $gt) = split /\t/, $_;
	my $REF;
	my $ALT;
	if ( exists( $ref_alt{$chr}{$pos} ) ) {
		$REF = $ref_alt{$chr}{$pos}{ref};
		$ALT = $ref_alt{$chr}{$pos}{alt};
	}
	else {
		$missingC += 1;
		# print "$chr $pos is not in ref_alt reference\n";
		print "$_\n";
		next;
	}

	if ($ref ne $REF and $ref ne $ALT) {
		$aberrantC += 1;
		# print "$chr $pos is aberrant\n";
		print "$_\n";
		next;
	}
	else {
		if ($alt eq ".") {
			if ($ref eq $REF){
				$alt = $ALT;
			}
			else {
				$alt = $REF;
			}
			$replacedC += 1;
			print "$chr\t$pos\t$rsid\t$ref\t$alt\t$qual\t$filter\t$info\t$format\t$gt\n";
		}
		else {
			$passC += 1;
			print "$_\n";
		}
	}
}

# print "total: $totalC\tmissing: $missingC\taberrant: $aberrantC\treplaced: $replacedC\tpass: $passC\n";