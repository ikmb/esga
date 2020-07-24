#!/usr/bin/env perl

use POSIX;
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
my $min_bit = 25; # bit score of the match must be at least this high
my $min_id = 70; # match must be at least this similar
my $length_percent = 0.6; # At least this much of the query has to align
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

while(<$IN>) {
	
	chomp; 
	my $line = $_; 
	my ($qseqid, $sseqid, $sstart, $send, $slen, $pident, $qlen, $qstart, $qend, $length, $mismatch, $gapopen, $evalue, $bitscore) = split("\t", $line);
	# ENSP00000332130.4       22      38025861        38026157        50818468        97.98   145     39      137     99      2       0       4e-64   190

	next if ($pident < $min_id);

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
		"pident" => $pident,
		"bitscore" => floor($bitscore)
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

	# add to our inventory ( a hash of hashes using the query and target id as keys)
	if ( exists $bucket{$qseqid} ) {
		if (exists $bucket{$qseqid}{$sseqid}) {
			push( @{ $bucket{$qseqid}{$sseqid} }, \%entry ); 
		} else {
			$bucket{$qseqid}{$sseqid} = [ \%entry ];
		}
	} else {
		$bucket{$qseqid}{$sseqid} = [ \%entry ] ;
	}

}

# All BLAST entries are grouped by query and target sequence, now stitch into clusters based on max_intron length

# Iterate over each protein to build clusters
foreach my $query ( keys %bucket ) {

	my $data = $bucket{$query};

	# every protein - chromosome/scaffold relationship
	# Our target are the outer most mapping coordinates +/- twice the max intron size
	foreach my $target ( keys %$data ) {
	
		# Sort the query/target specific matches based on their starting position
		my $matches = $bucket{$query}{$target};

		my @sorted_matches = sort { $a->{"target_start"} <=> $b->{"target_start"} } @{$matches};
		my @sorted_query_matches = sort { $a->{"query_start"} <=> $b->{"query_start"} } @{$matches};
	
		my $previous_m = @sorted_matches[0];
		my @linked_ms = ();
		foreach my $m (@sorted_matches) {
			my $distance = $m->{"target_start"}-$previous_m->{"target_end"} ;
			if ($distance > (2*$max_intron_size) ) {
				make_region(\@linked_ms);
				@linked_ms = ($m);
			} else {
	                        push(@linked_ms,$m);
			}
			$previous_m = $m;
		}

		make_region(\@linked_ms);

	}

}

sub make_region {

	my $matches = shift @_ ;
	#print Dumper($matches);

	my @sorted_query_matches = sort { $a->{"query_start"} <=> $b->{"query_start"} } @{$matches};

	my $first_entry = @{$matches}[0];
	my $last_entry = @{$matches}[-1];

	my $first_query_entry = @sorted_query_matches[0];

	my $qsum = 1;

	my @starts = ();
	my @ends = ();

	foreach my $qm (@sorted_query_matches) {
		push @starts, $qm->{"query_start"} ;
		push @ends, $qm->{"query_end"} ;		
	}

	# Condense the alignment intervals to a set of non-overlapping ranges
	for (1..$#starts) {
 		# extra check on array bounds, since we edit in-place
    		last unless $_ < @starts;
		# don't need to collapse if no overlap with previous end
		next unless $starts[$_] <= $ends[$_-1];
		# delete this start and the previous end
		splice(@starts,$_,1);
		splice(@ends,$_-1,1);
		# rerun this loop for the same value of $_ since it was deleted
		redo;
	}

	# Add up the total number of aminoacids aligned in this region
	foreach my $s (@starts) {
		my $e = shift @ends;
		$qsum += ($e-$s);
	}

	# How much of the total protein was aligned in this potential target region
	my $fraction = ($qsum/$first_entry->{"query_length"});
	warn $fraction . "\n";

	# This is a simplistic set of criteria to include a set of matches into the final list of targets...
        if ( $fraction > $length_percent ) {

		my $this_start = $first_entry->{"target_start"} ;
        	my $this_end = $last_entry->{"target_end"} ;
	        my $target_start;
        	my $target_end;

	        # Define the final genomic coordinates adding twice the maximum intron size as flanking regions (or start/end of scaffold, whichever comes first)
        	my $target_start = ($this_start <=  $max_intron_size*2) ? 1 : ($this_start-$max_intron_size*2);
	        my $target_end = ($this_end+$max_intron_size*2) >= $first_entry->{"target_length"} ? $first_entry->{"target_length"} : ($this_end+$max_intron_size*2);
		#printf "Stats " . $first_entry->{"target_id"} . " " . $first_entry->{"bitscore"} . " " .  $last_entry->{"bitscore"} . " " . $fraction . "\n";
        	printf $first_entry->{"query_id"} . "\t" . $first_entry->{"target_id"} . "\t" . $target_start . "\t" . $target_end . "\n";

	} 

}
