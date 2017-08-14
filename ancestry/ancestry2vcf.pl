#!/usr/bin/env perl

use strict;
use warnings;

use IO::File;
use POSIX qw/strftime/;

# need 3 arguments as input
usage() unless @ARGV >= 3;

my $PASS = "PASS";
my $skip_count = 0;
my $date = strftime('%Y%m%d',localtime);

# arg0->raw_path, arg1->output_dir, arg2->version 
my $raw_path = $ARGV[0];
my $output_dir = $ARGV[1];
my $version = $ARGV[2];
# my $hg = $ARGV[3];

# extract filename from raw_path
my @a = split('/', $raw_path);
my $b = $a[-1];
my $filename = substr($b, 0, -8);

# build output_path, ref_path, missing_ref_path
my $output_path = "$output_dir/$filename.vcf"; #print "$output_path\n";
my $ref_path = "/home/haoyuan/Downloads/reference/ancestry_v".$version."_hg19_ref.txt.gz"; #problem: cannot use ~/Documents/...
my $missing_ref_path = "$output_dir/$filename.sites_not_in_reference.txt";

# validate read files: raw_path and ref_path
missing($raw_path) unless -s $raw_path;
missing($ref_path) unless -s $ref_path;

#open files
my $fh = ($raw_path =~ m/zip$/) ? IO::File->new("gunzip -c $raw_path|") : IO::File->new($raw_path);
my $output_fh = IO::File->new(">$output_path");
my $ref_fh = IO::File->new("gunzip -c $ref_path|");
my $missing_ref_fh = -1;

#print the header for the VCF
print $output_fh "##fileformat=VCFv4.1\n";
print $output_fh "##fileDate=$date\n";
print $output_fh "##source=23andme2vcf.pl https://github.com/arrogantrobot/23andme2vcf\n";
print $output_fh "##reference=file://$ref_path\n";
print $output_fh "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print $output_fh "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGENOTYPE\n";

# load reference into memory, thank you @MattBrauer for this chunk of code! -- issue #4
my %ref = ();
while(<$ref_fh>) {
	chomp $_;
	my ($chr,$pos,$rsid,$ref) = split /\t/, $_;
	$ref{$chr}{$pos} = { 
		rsid => $rsid,
		ref => $ref,
	};
}
close $ref_fh;

#skip the meta-info and the header of the ancestry raw data file
my $line = $fh->getline;
while($line =~ m/^#/) {
	$line = $fh->getline;
}

#process ancestry raw data line by line
while(my $line = $fh->getline) {
	chomp $line;

	#read in a line of the ancestry data
	my ($rsid, $chr, $pos, $a1, $a2) = split /\t/, $line;
	if (not $rsid) { 
		$rsid = ".";
	}
	chomp $a2;
	#skip current line if the call was "--"
	if ($a1 eq '0') {
		next;
	}

	#skip insertions and deletions
	if (($a1 eq "D") || ($a1 eq "I")) { 
		next;
	}

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

	#get the reference base from ancestry ref 
	my $ref;
	if (exists($ref{$chr}{$pos})) {
		$ref = $ref{$chr}{$pos}{ref};
	} else {
		missing_sites($rsid,$chr,$pos);
		next;
	}

	if ($ref eq $PASS) {
		next;
	}

	#get the genotype
	my ($alt,$genotype) = getAltAndGenotype($ref, $a1, $a2);

	#capitalize ref & alt alleles
	$ref = uc($ref);
	$alt = uc($alt);

	#output a line of VCF data
	print $output_fh "$chr\t$pos\t$rsid\t$ref\t$alt\t.\t.\t.\tGT\t$genotype\n";
}

$fh->close;
$output_fh->close;
if ($missing_ref_fh != -1) {
	$missing_ref_fh->close;
}

skips();

#determine genotype
sub getAltAndGenotype {
	my $ref = shift;
	my $a = shift;
	my $c = shift;
	my $b = substr($c,0,1);

	if ($b !~ m/[A,C,G,T,N,a,c,g,t,n]/) {
		$b = undef;
	}
	my $alt;
	my $lc_ref = lc($ref);
	my $lc_a = lc($a);

	my $genotype;

	#determine which of the alleles are alts, if any
	if ($a && not $b) {
		if ($lc_a eq $lc_ref) {
			$alt = ".";
			$genotype = "0";
		} 
		else {
			$alt = $a;
			$genotype = "1";
		}
	} 
	else {
		my $lc_b = lc($b);
		if ($lc_a ne $lc_b) {
			if ($lc_ref eq $lc_a) {
				$alt = $b;
				$genotype = "0/1";
			} 
			elsif ($lc_ref eq $lc_b) {
				$alt = $a;
				$genotype = "0/1";
			} 
			else {
				$alt = "$a,$b";
				$genotype = "1/2";
			}
		} 
		else {
			if ($lc_a eq $lc_ref) {
				$alt = ".";
				$genotype = "0/0";
			} else {
				$alt = "$a";
				$genotype = "1/1";
			}
		}
	}
	return $alt, $genotype;
}

sub usage {
	print "usage: ./ancestry2vcf /path/to/ancestry/raw_data /path/to/output_dir version\n";
	exit(1);
}

sub missing {
	my $path = shift;
	print "Could not locate a file at: $path\n";
	usage();
}

sub missing_sites {
	my $rsid = shift;
	my $chr = shift;
	my $pos = shift;
	if ($missing_ref_fh == -1) {
		$missing_ref_fh = IO::File->new(">$missing_ref_path");
	}
	print $missing_ref_fh "$rsid\t$chr\t$pos\n";
	$skip_count++;
}

sub skips {
	if ($skip_count) {
		print "$skip_count sites were not included.\n"
	}
}