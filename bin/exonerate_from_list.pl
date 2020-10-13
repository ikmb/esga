#!/usr/bin/env perl
# Run exonerate from a list of accession numbers

use strict;
use Getopt::Long;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
	[--genome filename]
		Name of the CDBtools genome assembly index
	[--db filename]
		Name of the CDBtools protein fasta index
	[--accessions filename]
		List of protein IDs to use for alignment
	[--max_intron_size]
		Maximum length of intron to consider for spliced alignments
		
  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $accessions = undef;
my $db = undef;
my $max_intron_size = undef;
my $genome = undef;
my $help;

GetOptions(
    "help" => \$help,
    "accessions=s" => \$accessions,
    "genome=s" => \$genome,
    "max_intron_size=i" => \$max_intron_size,
    "db=s" => \$db,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if (!defined $max_intron_size){
	die "Must provide a maximum intron size (--max_intron_size)\n" ;
}
if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

print STDERR "Preparing exonerate jobs!\n" ;

open (my $IN, '<', $accessions) or die "FATAL: Can't open file: $accessions for reading.\n$!\n";

while (<$IN>) {
	
	chomp; 
	my $line = $_; 
	my $acc = $line;
	my $query_fa = $acc . ".fasta" ;
	# Run exonerate on these data
	my $cmd_run = "exonerate --model protein2genome --softmasktarget --percent 70 --bestn 2 --minintron 20 --maxintron $max_intron_size  --showalignment false --showtargetgff true $query_fa $genome > $acc.exonerate.align\n";
	
	printf($cmd_run);

}

close($IN);

print STDERR "Finished building exonerate jobs!\n";

