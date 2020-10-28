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
my $support = 1.0;
my $help;

GetOptions(
    "help" => \$help,
    "infile=s" => \$infile,
    "support=f" => \$support,
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

my @bucket;

foreach my $line (<$IN>) {

	chomp($line);

	
	if ($line =~ /.+ transcript supported by hints.*/) {
		my $this_support = (split " ", $line)[-1] ;
		if ($this_support >= $support) {
			foreach my $stored (@bucket) {
				printf $stored . "\n";
			}
			printf $line . "\n";
		}
		@bucket = ();
	} else {
		push (@bucket, $line);

	}
}

close($IN);

