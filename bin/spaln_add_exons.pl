#!/usr/bin/perl

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
my $help;

GetOptions(
    "help" => \$help,
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

printf "##gff-version 3\n";

my $gene_id;
my $transcript_id;
my $cds_counter = 0;

while (<$IN>) {
	chomp; 
	my $line = $_; 
	
	next if ($line =~ /^#.*/ );
	
	my @elements = split("\t", $line);
	
	my $feature = @elements[2];

	my @attributes = split ";" , @elements[-1] ;
	my %attribs;
	
	foreach my $a (@attributes) {
		my ($key,$value) = split "=", $a;
		$attribs{$key} = $value ;
	}	
		
	if ($feature eq "gene") {
		printf $line . "\n" ;
		$gene_id = $attribs{"ID"};
	} elsif ($feature eq "mRNA" || $feature eq "transcript" ) {
		$cds_counter = 0;
		$transcript_id = $attribs{"ID"};
		$line =~ s/transcript/mRNA/ ;
		printf $line . "\n";
	} elsif ($feature eq "CDS" || $feature eq "cds") {
		$cds_counter += 1;
		my $cds_line = $line;

		$attribs{"ID"} = $attribs{"ID"} . "-E" ;
		my $string = "";
		foreach my $k (keys %attribs) {
			my $val = $attribs{$k};
			$string .= $k . "=" . $val . ";" ;
		}
		@elements[2] = "exon" ; 
		@elements[-1] = $string ; 

		printf join( "\t" , @elements ) . "\n";
		printf $cds_line . "\n";
		
	}
	
}

close $IN;
