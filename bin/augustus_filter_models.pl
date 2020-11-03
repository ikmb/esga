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
my $p = 0;
my $w = 0;
my $e = 0;
my $t = 0;
my $stats = 0;
my $valid = 0;
my $this_support = 0;

foreach my $line (<$IN>) {

	chomp($line);

	if ($line =~  /^# start gene.*/) {

		if ($p > 0 ) {
                        $valid = 1;
                } else {
                        $valid = 0;
                }

                if ($this_support >= $support && $valid == 1) {
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

		$p = 0;
                $t = 0;
                $w = 0;
                $e = 0;
                @bucket = ();

	}
	
	if ($line =~ /.+ transcript supported by hints.*/) {
                push (@bucket, $line);

		$stats = 1;
                $this_support = (split " ", $line)[-1] ;

	} elsif ($line =~ /.* incompatible hint groups.*/) {

	} elsif ( $line =~ /^#\h+[P,W,T,E]:\h+[0-9]*.*/ ) {

                push (@bucket, $line);

		my @elements = split /\h+/, $line ;
		my $source = @elements[1];
		$source =~ s/:// ;
		my $count = @elements[2];
		if ($source eq "P" && $count > 0) {
			$p = 1;
		} elsif ($source eq "E") {
			$e = 1;
		} elsif ($source eq "T") {
			$t = 1;
		}	

	} else {
		push (@bucket, $line);

	}
}

close($IN);
close($GOOD);
close($BAD);
