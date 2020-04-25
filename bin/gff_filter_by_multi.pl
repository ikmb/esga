#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Data::Dumper;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--infile filename]
		The name of the file to read. 
    [--min_multi]
		The minimum coverage of this hint

  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $infile = undef;
my $min_cov = 5;

my $help;

GetOptions(
    "help" => \$help,
    "min_cov=i" => \$min_cov,
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

while (<$IN>) {

	my $line = $_;
	chomp $line;

	my ($seq,$source,$feature,$start,$stop,$phase,$strand,$score,$info) = split("\t",$line);	

	my %entry = { "seq" => $seq, "source" => $source, "feature" => $feature, "start" => $start, "stop" => $stop, "phase" => $phase, "strand" => $strand, "score" => $score } ;

	my %attribs;
	
	my @fields = split(";",$info);

	foreach my $f (@fields) {
		my ($key,$value) = split("=",$f);
		$attribs{$key} = $value;
	}

	if (defined $attribs{'mult'}) {	
		my $cov = $attribs{'mult'};

		if ($cov >= $min_cov) {
			printf $line . "\n";
		}
	}

}

