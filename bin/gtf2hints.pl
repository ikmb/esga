#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--gtf filename]
		The name of the file to read. 
    [--source name]

    [--pri value]

  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $gtf = undef;
my $source = "T";
my $pri = "4";
my $help;

GetOptions(
    "help" => \$help,
    "gtf=s" => \$gtf,
    "source=s" => \$source,
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

open (my $IN, '<', $gtf) or die "FATAL: Can't open file: $gtf for reading.\n$!\n";

my $is_first_cds;
my $is_last_cds;

while (<$IN>) {

        my $line = $_;
        chomp $line;

        my ($seq,$src,$feature,$start,$stop,$score,$strand,$phase,$info) = split("\t",$line);

	
	my $hint_type = "" ;

	if ($feature eq "transcript") {
		$is_first_cds = 0;
		$is_last_cds = 0;
	} elsif ($feature eq "CDS") {
		$hint_type = "CDSpart";
	} elsif ($feature eq "exon") {
		$hint_type = "exon";
	} else {
		next;
	}
	
	my %attribs;

        my @fields = split(";",$info);

        foreach my $f (@fields) {
                my ($key,$value) = split(" ",$f);
                $attribs{$key} = $value;
        }

	my $group = $attribs{"transcript_id"};

	$group =~ s/\"//g ;
	
	printf $seq . "\t" . "transmapped" . "\t" . $hint_type . "\t" . $start . "\t" . $stop . "\t" . "." . "\t" . $strand . "\t" . "." . "\t" . "group=$group;source=$source;pri=$pri\n";

	
}

close($IN);

