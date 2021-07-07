#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--fasta filename]
		The name of the assembly fasta
    [--list  filename]
		a list of scaffolds<->proteins
		
  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $fasta = undef;
my $list = undef;
my $model = undef;
my $help;

GetOptions(
    "help" => \$help,
    "fasta=s" => \$fasta,
    "list=s" => \$list,
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

open (my $TARGETS, '<', $list) or die "FATAL: Can't open file: $list for reading.\n$!\n";

my @accs;

while (<$TARGETS>) {

	my $line = $_;
	chomp($line);

	my $acc = $line;

	push(@accs,$acc);
	
}

my @uaccs = uniq(@accs);

foreach my $a (@uaccs) {
	my $ua = $a . ".fa";
	my $ec = "samtools faidx $fasta $a > $ua";
	printf $ec . "\n";
}

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

