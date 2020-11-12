#!/usr/bin/env perl

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
my $pri = 3;

my $help;

GetOptions(
    "help" => \$help,
    "infile=s" => \$infile,
    "pri=i" => \$pri,
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

	printf $seq . "\t" . "repeatmasker" . "\t" . "nonexonpart" . "\t" .  $start .  "\t" . $stop . "\t" . $score . "\t" . $strand . "\t.\t" . "src=RM;pri=$pri\n";

	

}

close($IN);
