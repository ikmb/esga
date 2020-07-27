#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--targets filename]
		The name of the targets file to read. 
    [--proteins filename]
		The name of the protein fasta file

  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $targets = undef;
my $proteins = undef;

my $help;

GetOptions(
    "help" => \$help,
    "proteins=s" => \$proteins,
    "targets=s" => \$targets,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

open (my $IN, '<', $targets) or die "FATAL: Can't open file: $targets for reading.\n$!\n";

my %target_list;
my @protein_accs;

while (<$IN>) {

	my $line = $_;
	chomp $line;
	my $acc = (split /\t/, $line)[0];
	$target_list{$acc} = 1;
}

close($IN);

open (my $PROT, '<', $proteins) or die "FATAL: Can't open file: $proteins for readind.\n";

while (<$PROT>) {

	my $line = $_;
	chomp $line;

	if ($line =~ /^>.*/) {
		my $def_line = $line;
		$def_line =~  s/^>// ;
		my $acc = (split " ",$def_line)[0];
		push(@protein_accs,$acc);
	}
}

close($PROT);

foreach my $prot (@protein_accs) {

	if (!defined $target_list{$prot}) {
		printf $prot . "\n";
	}
}

