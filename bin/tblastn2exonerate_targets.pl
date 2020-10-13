#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Data::Dumper;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--infile filename]
		The name of the file to read. 
	[--max_intron_size]
		Maximum length of intron for spliced alignments
	[--min_bit]
		Minimum bitscore for match to be considered valid (default: 50)
	[--min_id]
		Miniumum id% for match to be considered valid (default: 0.6)
	[--length_percent ]
		Minimum length percentage for the match to be considered valid (default: 0.6)
  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $infile = undef;
my %targets;
my $min_bit = 25.0; # bit score of the match must be at least this high
my $min_id = 60; # match must be at least this similar
my $length_percent = 0.7; # At least this much of the query has to align
my $max_intron_size = undef;

my $help;

GetOptions(
    "help" => \$help,
    "infile=s" => \$infile,
    "max_intron_size=i" => \$max_intron_size,
    "min_bit=i" => \$min_bit,
    "min_id=i" => \$min_id,
    "length_percent=f" => \$length_percent,
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

open (my $IN, '<', $infile) or die "FATAL: Can't open file: $infile for reading.\n$!\n";

my %bucket;

my $protein_id = undef;

while(<$IN>) {
	
	chomp; 
	my $line = $_; 

	my ($qseqid, $sseqid, $sstart, $send, $slen, $pident, $qlen, $qstart, $qend, $length, $mismatch, $gapopen, $evalue, $bitscore) = split("\t", $line);

	if (defined $protein_id && $protein_id ne $qseqid) {
		# process cluster
		process_protein_hits(%bucket);
		%bucket = ();
	}

	my %entry = (
		"query_id" => $qseqid,
		"target_id" => $sseqid,
		"target_start" => $sstart,
		"target_end" => $send,
		"target_length" => $slen,
		"query_length" => $qlen,
		"query_start" => $qstart,
		"query_end" => $qend,
		"aln_length" => $length,
		"mismatch" => $mismatch,
		"gapopen" => $gapopen,
		"evalue" => $evalue,
		"bitscore" => $bitscore
	);
	# If this match is in reverse direction
	if ($entry{"target_end"} < $entry{"target_start"} ) {
		my $new_start = $entry{"target_end"};
		my $new_end = $entry{"target_start"};
		$entry{"target_start"} = $new_start;
		$entry{"target_end"} = $new_end;
	}
	if ($entry{"query_end"} < $entry{"query_start"}){
		 my $new_start = $entry{"query_end"};
                my $new_end = $entry{"query_start"};
                $entry{"query_start"} = $new_start;
                $entry{"query_end"} = $new_end;

	}

	$protein_id = $qseqid;

        # add to our inventory ( a hash of hashes using the query and target id as keys)
        if (exists $bucket{$sseqid}) {
	        push( @{ $bucket{$sseqid} }, \%entry );
        } else {
        	$bucket{$sseqid} = [ \%entry ];
        }

}

close(IN);

# process final cluster
process_protein_hits(%bucket);


###################################
############## SUB PROCESS HERE ###
###################################

sub process_protein_hits {

	my %data = @_;

        foreach my $target ( keys %data ) {

                # Sort the query/target specific matches based on their starting position
                my $matches = $bucket{$target};

                my @sorted_matches = sort { $a->{"target_start"} <=> $b->{"target_start"} } @{$matches};
                my @sorted_query_matches = sort { $a->{"query_start"} <=> $b->{"query_start"} } @{$matches};

                # Sorted by query (protein) positions to determine how much of the query sequence is aligned in this scaffold
                my $first_query_entry = @sorted_query_matches[0];
                my $last_query_entry = @sorted_query_matches[-1];

                my $query_length = $first_query_entry->{"query_length"};
                my $query_covered = merge_ranges(@sorted_query_matches);
		
		my $sum = 0 ;
		foreach my $range (@$query_covered) {
			$sum += ($range->{"query_end"}-$range->{"query_start"});
		}
		
                my $fraction = $sum/$query_length;

                # Sorted by target positions to determine the genomic boundaries for subsequent exonerate alignments
                my $first_entry = @sorted_matches[0];
                my $last_entry = @sorted_matches[-1];

                #printf STDERR "${query}\t${target}\tQL:${query_length}\tML:${match_length}\tFR:${fraction}\n";

                # is this overall a valid match?
                if ($first_entry->{"bitscore"} >= $min_bit && $last_entry->{"bitscore"} >= $min_bit && $fraction >= $length_percent ) {

                        my $this_start = $first_entry->{"target_start"} ;
                        my $this_end = $last_entry->{"target_end"} ;

                        my $target_start;
                        my $target_end;
			my $target_length;

                        # Define the final genomic coordinates adding twice the maximum intron size as flanking regions (or start/end of scaffold, whichever comes first)
                        my $target_start = ($this_start <=  $max_intron_size*2) ? 1 : ($this_start-$max_intron_size*2);
                        my $target_end = ($this_end+$max_intron_size*2) >= $first_entry->{"target_length"} ? $first_entry->{"target_length"} : ($this_end+$max_intron_size*2);
			my $target_length = $first_entry->{"target_length"} ;
                        #printf $first_entry->{"query_id"} . "\t" . $first_entry->{"target_id"} . "\t" . 1 . "\t" . $target_length . " " . $fraction . "\n";

                         printf $first_entry->{"query_id"} . "\t" . $first_entry->{"target_id"} . "\t" . $target_start . "\t" . $target_end . "\n";

                }

        }

}

sub merge_ranges {
	
	my @matches =  @_;

	my @ranges;
	my $first_range = shift @matches;
	my %this_range = ( "query_start" => $first_range->{"query_start"} , "query_end" => $first_range->{"query_end"} );

	foreach my $match (@matches) {

		if ($match->{"query_start"} < $this_range{"query_end"} && $match->{"query_end"} > $this_range{"query_end"} ) {
			%this_range = ( "query_start" => $this_range{"query_start"} , "query_end" => $match->{"query_end"} ) ;
		} elsif ( $match->{"query_end"} < $this_range{"query_end"} && $match->{"query_start"} > $this_range{"query_start"} ) {
			# do nothing
		} elsif ( $match->{"query_start"} == $this_range{"query_start"} && $match->{"query_end"} == $this_range{"query_end"} ) {
			# do nothing
		} elsif ($match->{"query_start"} > $this_range{"query_end"}) {
			my %store = ( "query_start" => $this_range{"query_start"} , "query_end" => $this_range{"query_end"} ) ;
			push @ranges, \%store;
			%this_range = ( "query_start" => $match->{"query_start"} , "query_end" => $match->{"query_end"} ) ;
		}
	
	}

	push @ranges, \%this_range;
	return \@ranges;

}
