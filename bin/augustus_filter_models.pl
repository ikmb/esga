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
};

my $infile = undef;
my $support = 1.0;
my $help;

GetOptions(
    "help" => \$help,
    "infile=s" => \$infile,
    "support=f" => \$support);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

my $good_models = $infile . ".good.gff";
my $bad_models = $infile . ".bad.gff";

open (my $IN, '<', $infile) or die "FATAL: Can't open file: $infile for reading.\n$!\n";
open (my $GOOD, '>', $good_models) or die $! ;
open (my $BAD, '>', $bad_models) or die $!;

my @bucket;

foreach my $line (<$IN>) {

	chomp($line);

	
	if ($line =~ /.+ transcript supported by hints.*/) {
		my $this_support = (split " ", $line)[-1] ;
		if ($this_support >= $support) {
			foreach my $stored (@bucket) {
				print $GOOD $stored . "\n";
			}
			print $GOOD $line . "\n";
		} else {
			foreach my $stored (@bucket) {
                                print $BAD $stored . "\n";
                        }
                        print $BAD $line . "\n";
		}
		@bucket = ();
	} else {
		push (@bucket, $line);

	}
}

close($IN);
close($GOOD);
close($BAD);
