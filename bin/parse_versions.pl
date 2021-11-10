#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Cwd qw(getcwd);

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $help;

GetOptions(
    "help" => \$help,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

my $header = qq(
id: 'software_versions'
section_name: 'Software Versions'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |\n  <dl class="dl-horizontal">
);

printf $header . "\n";

my $directory = getcwd;

opendir (DIR, $directory) or die $!;

my $version = undef;
my $tool = undef;

while (my $file = readdir(DIR)) {

	next unless ($file =~ /^v_.*\.txt$/);
	open (my $IN, '<', $file) or die "FATAL: Can't open file: $file for reading.\n$!\n";
	
	chomp(my @lines = <$IN>);	
	
	if ($file =~ /v_gatk.*/) {
		my $line = @lines[0];
		$tool = "GATK4" ;
		$version = (split " ", $line)[-1];
		
	} elsif ($file =~ /^v_nextflow\.txt$/ )  {
		my $line = @lines[0];
		$tool = "Nextflow";
		$version = $line;
	} elsif ($file =~ /^v_ikmb_esga\.txt$/) {
		my $line = @lines[0];
                $tool = "ESGA Pipeline";
		$version = $line;
	} elsif ($file =~ /^v_picard\.txt$/) {
		my $line = @lines[0];
                $tool = "Picard";
                $version = (split " ", $line)[-1];
	} elsif ($file =~ /^v_bwa\.txt$/) {
		my $line = @lines[2];
                $tool = "BWA";
                $version = (split " ", $line)[-1];
	} elsif ($file =~ /^v_trinity\.txt$/) {
		my $line = @lines[0];
		$tool = "Trinity";
		$version = (split " ", $line)[-1];
	} elsif ($file =~ /^v_fastp\.txt$/) {
		my $line = @lines[0];
		$tool = "fastP";
		$version = (split " ", $line)[-1];
	} elsif ($file =~ /^v_augustus\.txt$/) {
                my $line = @lines[0];
                $tool = "AUGUSTUS";
                $version = (split " ", $line)[1];
	} elsif ($file =~ /^v_minimap2\.txt$/) {
                my $line = @lines[0];
                $tool = "Minimap2";
                $version = (split " ", $line)[0];
	 } elsif ($file =~ /^v_gffread\.txt$/) {
                my $line = @lines[0];
                $tool = "GFFread";
                $version = (split " ", $line)[0];
	} elsif ($file =~ /^v_spaln\.txt$/) {
                my $line = @lines[2];
                $tool = "Spaln";
                $version = (split " ", $line)[3];

	} else {
		my $line = @lines[0];
		my @elements = (split " ",$line);
		$tool = @elements[0];
		$tool =~ s/\,//;
		$version = @elements[-1];
	}
	
	#my $entry = qq(
	my $entry = "<dt>$tool</dt><dd><samp>$version</samp></dd>" ;
	#);
	printf "    $entry\n";
	
	close($IN);
	
}

close(DIR);

printf "  </dl>\n";


