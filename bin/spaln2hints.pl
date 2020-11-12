#!/usr/bin/env perl
# Converts Exonerate output to GFF3 hints
use strict;
use Getopt::Long;


my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--infile filename]
		The name of the file to read. 
    [--pri integer]
		Priority of the resulting hints (default: 5)
  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $infile = undef;
my $pri = 5;
my $help;
my $src = "P";
my $method = "protein2genome";

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


my $GeneID;

while (<$IN>) {
	
	my $line = $_;
	chomp $line;

	# skip comment lines
	next if ($line =~ m/^#.*/ );

	my ($Chrom,$met,$feature,$start,$end,$score,$strand,$frame,$comment) = split(/\t+/,$line);
	my %attributes = ();
	my @attribs = split ";", $comment;
	foreach my $a (@attribs) {
		my ($key,$value) = (split "=", $a);
		$attributes{$key} = $value;
	}
	
	if ($feature eq "gene") {
		$GeneID = $attributes{'ID'};
		my $multi = $attributes{'isoforms'};
		printf $Chrom."\t".$method."\tgenicpart\t" . ($start-20) . "\t" . ($end+20) . "\t1000\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;mult=$multi;pri=$pri\n";
	} elsif ($feature eq "mRNA") {
		$GeneID = $attributes{'ID'};
		if ($strand eq "+") {
			printf $Chrom."\t".$method."\tstart\t". ($start-20) . "\t" . ($start+20) . "\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;pri=$pri\n";
			printf $Chrom."\t".$method."\tstop\t" . ($end-20) . "\t" . ($end+20) ."\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;pri=$pri\n";
		} else {
			printf $Chrom."\t".$method."\tstop\t". ($end-20) . "\t" . ($end+20) . "\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;pri=$pri\n";
                        printf $Chrom."\t".$method."\tstart\t" . ($start-20) . "\t" . ($start+20) ."\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;pri=$pri\n";
		}
	} elsif ($feature eq "cds") {
		$GeneID = $attributes{'Parent'};
		printf $Chrom."\t".$method."\tCDSpart\t".$start."\t".$end."\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;pri=$pri\n";
	} else {
		next;
	}
}

close $IN;
