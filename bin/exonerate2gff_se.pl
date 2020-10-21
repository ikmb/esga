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
    [--source string]
		A valid source for processing (est, protein or trinity)
    [--pri integer]
		Priority of the resulting hints (default: 5)
  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $infile = undef;
my $source = undef;
my $GeneID = undef;
my $pri = 5;
my $help;

my %source_keys = (
	"est" => {
		"src" => "E",
		"pri" => 3,
		"source" => "est2genome"
	},
	"protein" => {
		"src" => "P",
		"pri" => 5,
		"source" => "protein2genome"
		
	},
	"trinity" => {
		"src" => "T",
		"pri" => 3,
		"source" => "trinity2genome"
	}
) ;
my $help;

GetOptions(
    "help" => \$help,
    "infile=s" => \$infile,
    "pri=i" => \$pri,
    "source=s" => \$source,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

# make sure the source is valid
my @sources = keys %source_keys;

unless ( grep( /^$source$/, @sources ) ) {
  exit 1, "Did not provide a valid source (${join(',', @{ keys %source_keys })})\n";
}

# open the exonerate report 
open (my $IN, '<', $infile) or die "FATAL: Can't open file: $infile for reading.\n$!\n";

my $src = $source_keys{$source}{"src"};
my $method = $source_keys{$source}{"source"};

my $exon_count = 0;
my @bucket;

while (<$IN>) {
	
	my $line = $_;
	chomp $line;

	# skip comment lines
	next if ($line =~ m/^#.*/ );

	my ($Chrom,$met,$feature,$start,$end,$score,$strand,$frame,$comment) = split(/\t+/, $line);

	if ($feature eq "gene") {

		if (scalar @bucket > 0) {
		
			// for est data, we only accept spliced hints / multi-exon genes
			if ($method eq "est2genome" || $method eq "trinity2genome" ) { 
				if ($exon_count > 1) {
					foreach my $e (@bucket) {
						printf $e . "\n";
					}
				}
			else {
                        	foreach my $e (@bucket) {
	                        	printf $e . "\n";
				}
			}
		}

		@bucket = ();		
		$exon_count = 0;

		($GeneID) =($comment =~/gene_id\s\w+\s;\ssequence\s(\S+)\s;\s/);		
		# if this is a protein alignment, we use the bounds as start and stop hints
		if ($method eq "protein2genome") {
			printf STDERR $strand . "\n";
			if ($strand eq "+") {
				my $shint = $Chrom."\t".$method."\tstart\t". ($start-20) . "\t" . ($start+20) . "\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;pri=$pr";
				my $ehint =  $Chrom."\t".$method."\tstop\t" . ($end-20) . "\t" . ($end+20) ."\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;pri=$pri";
				push @bucket, $shint;
				push @bucket, $ehint;
				
			} else {
				my $shint = $Chrom."\t".$method."\tstop\t". ($start-20) . "\t" . ($start+20) . "\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;pri=$pr";
                                my $ehint $Chrom."\t".$method."\tstart\t" . ($end-20) . "\t" . ($end+20) ."\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;pri=$pri";
				push @bucket, $shint;
                                push @bucket, $ehint;
			}
			my $ghint = $Chrom."\t".$method."\tgenicpart\t" . ($start-20) . "\t" . ($end+20) . "\t" . $score . "\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;pri=$pri";	
			push @bucket, $ghint;
		}
	} elsif ($feature eq "exon") {
		$exon_count += 1;
		my $exon_hint = $Chrom."\t".$method."\texonpart\t".$start."\t".$end."\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;pri=$pri";
		push @bucket, $exon_hint;
	} elsif ($feature eq "cds") {
		my $cds_hint = $Chrom."\t".$method."\tCDSpart\t".$start."\t".$end."\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;pri=$pri";
		push @bucket, $cds_hint;
	} elsif ($feature eq "intron") {
		my $intron_hint = "";
                # For transcriptome data, we assume that intron positions are accurate - otherwise merely use intronpart hints
                if ($method  eq "est2genome" || $method eq "trinity2genome") {
                        $intron_hint = $Chrom."\t".$method."\tintron\t".$start."\t".$end."\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;pri=$pri";
                } else {
                        $intron_hint = $Chrom."\t".$method."\tintronpart\t".$start."\t".$end."\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;pri=$pri";
                }
                push @bucket, $intron_hint;
	} elsif ($feature eq "utr5" || $feature eq "utr3") {
		my $utr_hint = $Chrom."\t".$method."\tUTRpart\t".$start."\t".$end."\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;pri=$pri";
		push @bucket, $utr_hint;
	} else {
		next;
	}
}

close $IN;
