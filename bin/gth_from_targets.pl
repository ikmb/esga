#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--targets  filename]
		a list of scaffolds<->proteins
		
  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $targets = undef;
my $options = "";
my $model = undef;
my $help;

GetOptions(
    "help" => \$help,
    "targets=s" => \$targets,
    "options=s" => \$options,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

# Read targets file to get list of scaffolds

open (my $TARGETS, '<', $targets) or die "FATAL: Can't open file: $targets for reading.\n$!\n";

my @jobs;
my @indexing;

while (<$TARGETS>) {

	my $line = $_;
	chomp($line);

	my ($prot,$chr) = split("\t",$line);

	my $prot_file = $prot . ".fa" ;
	my $chr_file = $chr . ".fa" ;
	my $gth_out = $chr . "." . $prot . ".gth.out" ;

	my $command = "gth -genomic genome_db/$chr_file -protein protein_db/$prot_file $options -gff3out -skipalignmentout -skipindexcheck -o $gth_out" ;
	push(@jobs,$command);
	my $command = "gth -genomic genome_db/$chr_file -protein protein_db/$prot_file -gff3out -skipalignmentout -o $gth_out -createindicesonly" ;
	push(@jobs,$command);
	
}

foreach my $c (@jobs) {
	printf $c . "\n";
}

close($TARGETS);

