#!/usr/bin/perl

use strict;
use Getopt::Long;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--infile filename]
		The name of the file to read. 
  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $infile = undef;
my $help;

GetOptions(
    "help" => \$help,
    "infile=s" => \$infile,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

open (my $IN, '<', $infile) or die "FATAL: Can't open file: $infile for reading.\n$!\n";

printf "##gff-version 3\n";

my $gene_id;
my $transcript_id;
my $cds_counter = 0;

while (<$IN>) {
	chomp; 
	my $line = $_; 
	
	next if ($line =~ /^#.*/ );
	
	my @elements = split("\t", $line);
	
	my $feature = @elements[2];
	
	if ($feature eq "gene") {
		$gene_id = @elements[-1];
		@elements[-1] = "ID=$gene_id";
		printf join("\t",@elements) . "\n";
	} elsif ($feature eq "mRNA" || $feature eq "transcript") {
		$cds_counter = 0;
		@elements[2] = "mRNA";
		$transcript_id = @elements[-1];
		@elements[-1] = "ID=$transcript_id;Parent=$gene_id";
		printf join("\t",@elements) . "\n";
	} elsif ($feature eq "CDS" || $feature eq "cds") {
		$cds_counter += 1;
		@elements[-1] = "ID=$transcript_id.CDS-$cds_counter;Parent=$transcript_id";
		my $cds_line = join("\t",@elements);
		@elements[-1] = "ID=$transcript_id.EXON-$cds_counter;Parent=$transcript_id";
		@elements[2] = "exon";
		@elements[7] = ".";
		my $exon_line = join("\t",@elements);
		printf $exon_line . "\n";
		printf $cds_line . "\n";
	
	}
	
}

close $IN;
